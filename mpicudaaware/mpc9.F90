PROGRAM MPICUDAAWARE

USE CUDAFOR
USE OPENACC
USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV, MPL_BARRIER, MPL_INIT, MPL_END, MPL_WAIT, &
                     & JP_NON_BLOCKING_STANDARD

IMPLICIT NONE

INCLUDE 'mpif.h'

REAL, ALLOCATABLE         :: C_SEND (:), C_RECV (:)
REAL, ALLOCATABLE, DEVICE :: D_SEND (:), D_RECV (:)
REAL, ALLOCATABLE         :: H_SEND (:), H_RECV (:)
REAL, ALLOCATABLE         :: H_SEND2 (:),H_RECV2(:)

INTEGER :: IREQ_RECV, IREQ_SEND
INTEGER :: ISTATUS (MPI_STATUS_SIZE)
INTEGER :: ISIZE, IERROR
INTEGER :: IRANK, IRANKP, IRANKN
INTEGER, PARAMETER :: taille=1024*1024*64
INTEGER     :: num_carte,compteur,repetition

integer(kind=8)  :: count_val1, count_val2
integer(kind=8)  :: count_rate
integer(kind=8)  :: count_max
real             :: taux

logical, parameter :: detail = .false.

CALL MPL_INIT ()

CALL MPI_COMM_RANK (MPI_COMM_WORLD, IRANK, IERROR)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, ISIZE, IERROR)

call acc_set_device_num(modulo(IRANK,4),ACC_DEVICE_NVIDIA)
call acc_init(ACC_DEVICE_NVIDIA)
num_carte = acc_get_device_num(ACC_DEVICE_NVIDIA)
print *, "========================="
print *, " test of transfer between ",isize," ranks ; size (real 4) : ",taille
print *, "Rank ",IRANK ," on card number ",num_carte

ALLOCATE (C_SEND (taille))
ALLOCATE (C_RECV (taille))
ALLOCATE (D_RECV (taille))
ALLOCATE (D_SEND (taille))
ALLOCATE (H_RECV (taille))
ALLOCATE (H_SEND (taille))
ALLOCATE (H_RECV2 (taille))
ALLOCATE (H_SEND2 (taille))

CALL RANDOM_NUMBER (C_SEND)
CALL RANDOM_NUMBER (H_SEND)

do compteur=1,taille
  H_SEND(compteur) = H_SEND(compteur)+irank
  H_RECV(compteur) = -1.0
enddo

do compteur=1,taille
  h_send2(compteur) = h_send(compteur)+3*irank
  H_RECV2(compteur) = -1.0
enddo

!$acc data copy(h_send2,h_recv2)

print *," "
print *,"====================="

if (detail) then
  print *,"rank : ",irank," ; data used : "
  print *, "c_send", C_SEND(1:10)
  print *, "h_send", H_SEND(1:10)
endif

do repetition=1,40

do compteur=1,taille
  C_SEND(compteur) = C_SEND(compteur)+irank
  C_RECV(compteur) = -1.0
enddo
D_SEND=C_SEND
D_RECV=C_RECV

IRANKP = MODULO (IRANK-1, ISIZE)
IRANKN = MODULO (IRANK+1, ISIZE)

!$acc wait
call system_clock(count_val1,count_rate,count_max)
taux=real(count_rate)

print *," "
print *,"=========================="
PRINT *, " DEVICE DATA "

CALL MPL_RECV (D_RECV, KSOURCE=IRANKP+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_RECV, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)

CALL MPL_BARRIER ()

CALL MPL_SEND (D_SEND, KDEST=IRANKN+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_SEND, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)

CALL MPL_WAIT (IREQ_RECV)
CALL MPL_WAIT (IREQ_SEND)

!$acc wait
call system_clock(count_val2)
print *,"rank : ",irank," ; time spent, case device : ",(count_val2-count_val1)/taux
print *," "

C_SEND=D_SEND
C_RECV=D_RECV
if (detail) then
  PRINT *, "rank : ",IRANK, " ; data received ==> ", C_RECV (1:10)
  PRINT *, "rank : ",IRANK, " ; data sent  ==> ", C_SEND(1:10)
endif

do compteur=1,taille
  H_SEND(compteur) = H_SEND(compteur)+irank
  H_RECV(compteur) = -1.0
enddo

!$acc wait
call system_clock(count_val1)

print *," "
print *,"====================="
PRINT *, " HOST DATA "

CALL MPL_RECV (H_RECV, KSOURCE=IRANKP+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_RECV, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)

CALL MPL_BARRIER ()

CALL MPL_SEND (H_SEND, KDEST=IRANKN+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_SEND, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)

CALL MPL_WAIT (IREQ_RECV)
CALL MPL_WAIT (IREQ_SEND)

!$acc wait
call system_clock(count_val2)
print *,"rank : ",irank," ; time spent, case host data ",(count_val2-count_val1)/taux
print *," "
if (detail) then
  PRINT *, "rank : ",IRANK, " ; data received  ==> ", H_RECV (1:10)
  PRINT *, "rank : ",IRANK, " ; data sent  ==> ", H_SEND(1:10)
endif

!$acc kernels
do compteur=1,taille
  H_SEND2(compteur) = H_SEND2(compteur)+irank
  H_RECV2(compteur) = -1.0
enddo
!$acc end kernels

!$acc wait
call system_clock(count_val1)

print *," "
print *,"=================================="
PRINT *, " HOST DATA (USE_DEVICE) "

!$ACC HOST_DATA USE_DEVICE (H_RECV2)
CALL MPL_RECV (H_RECV2, KSOURCE=IRANKP+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_RECV, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)
!$ACC END HOST_DATA

CALL MPL_BARRIER ()

!$ACC HOST_DATA USE_DEVICE (H_SEND2)
CALL MPL_SEND (H_SEND2, KDEST=IRANKN+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_SEND, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)
!$ACC END HOST_DATA

CALL MPL_WAIT (IREQ_RECV)
CALL MPL_WAIT (IREQ_SEND)

!$acc wait
call system_clock(count_val2)
print *,"rank : ",irank," ; time spent, case host_data use device ",(count_val2-count_val1)/taux
print *," "

if (detail) then
  !$acc update host(h_recv2,h_send2)
  PRINT *,"rank : ", IRANK, " ; received  ==> ", H_RECV2 (1:10)
  PRINT *,"rank : ", IRANK, " ; sent ==> ", H_SEND2(1:10)
endif

enddo !!repetition

!$acc end data

DEALLOCATE (C_SEND)
DEALLOCATE (C_RECV)
DEALLOCATE (D_SEND)
DEALLOCATE (D_RECV)
DEALLOCATE (H_SEND)
DEALLOCATE (H_RECV)
DEALLOCATE (H_RECV2)
DEALLOCATE (H_SEND2)

PRINT *, "FINALIZE"

CALL MPL_END
 
END PROGRAM MPICUDAAWARE
