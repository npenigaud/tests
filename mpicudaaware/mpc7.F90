PROGRAM MPICUDAAWARE

USE CUDAFOR
USE OPENACC
USE YOMHOOK,    ONLY : LHOOK, DR_HOOK, JPHOOK
USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV, MPL_BARRIER, MPL_INIT, MPL_END, MPL_WAIT, &
                     & JP_NON_BLOCKING_STANDARD

IMPLICIT NONE

INCLUDE 'mpif.h'

REAL, ALLOCATABLE         :: C_SEND (:), C_RECV (:)
REAL, ALLOCATABLE, DEVICE :: D_SEND (:), D_RECV (:)
REAL, ALLOCATABLE         :: H_SEND (:), H_RECV (:)
REAL, ALLOCATABLE         :: dummy1(:), dummy2(:)

INTEGER :: IREQ_RECV, IREQ_SEND
INTEGER :: ISTATUS (MPI_STATUS_SIZE)
INTEGER :: ISIZE, IERROR
INTEGER :: IRANK, IRANKP, IRANKN
INTEGER, PARAMETER :: taille=1024*1024
INTEGER     :: num_carte,compteur
REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

integer(kind=8)  :: count_val1, count_val2
integer(kind=8)  :: count_rate
integer(kind=8)  :: count_max
real             :: taux

CALL MPL_INIT ()

CALL MPI_COMM_RANK (MPI_COMM_WORLD, IRANK, IERROR)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, ISIZE, IERROR)

call acc_set_device_num(IRANK,ACC_DEVICE_NVIDIA)
call acc_init(ACC_DEVICE_NVIDIA)
num_carte = acc_get_device_num(ACC_DEVICE_NVIDIA)
print *, "Rang ",IRANK ," sur la carte numéro ",num_carte

ALLOCATE (C_SEND (taille))
ALLOCATE (C_RECV (taille))
ALLOCATE (D_RECV (taille))
ALLOCATE (D_SEND (taille))
ALLOCATE (H_RECV (taille))
ALLOCATE (H_SEND (taille))
allocate (dummy1(taille),dummy2(taille))

CALL RANDOM_NUMBER (C_SEND)
CALL RANDOM_NUMBER (H_SEND)
call random_number (dummy1)
call random_number (dummy2)
!$acc data copy(dummy1,dummy2)

print *," "
print *,"====================="
print *,"données utilisées"
print *, "c_send", C_SEND(1:10)
print *, "h_send", H_SEND(1:10)

do compteur=1,taille
  C_SEND(compteur) = C_SEND(compteur)+irank
  C_RECV(compteur) = -1.0
enddo
D_SEND=C_SEND
D_RECV=C_RECV

IRANKP = MODULO (IRANK-1, ISIZE)
IRANKN = MODULO (IRANK+1, ISIZE)

!$acc kernels
do compteur=1,taille
  dummy2(compteur)=dummy2(compteur)+1.01*dummy1(compteur)
enddo
!$acc end kernels

!$acc wait
call system_clock(count_val1,count_rate,count_max)
taux=real(count_rate)

CALL DR_HOOK('device',0,ZHOOK_HANDLE) 

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

C_SEND=D_SEND
C_RECV=D_RECV
PRINT *, IRANK, " reçu ==> ", C_RECV (1:10)
PRINT *, IRANK, " envoyé ==> ", C_SEND(1:10)

CALL DR_HOOK('device',1,ZHOOK_HANDLE)
!$acc wait
call system_clock(count_val2)
print *,irank," temps device : ",(count_val2-count_val1)/taux

!$acc kernels
do compteur=1,taille
  dummy2(compteur)=dummy2(compteur)+1.01*dummy1(compteur)
enddo
!$acc end kernels


do compteur=1,taille
  H_SEND(compteur) = H_SEND(compteur)+irank
  H_RECV(compteur) = -1.0
enddo

!$acc kernels
do compteur=1,taille
  dummy2(compteur)=dummy2(compteur)+1.01*dummy1(compteur)
enddo
!$acc end kernels

!$acc wait
call system_clock(count_val1)
CALL DR_HOOK('host data',0,ZHOOK_HANDLE)

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

PRINT *, IRANK, " reçu  ==> ", H_RECV (1:10)
PRINT *, IRANK, " envoyé ==> ", H_SEND(1:10)

CALL DR_HOOK('host data',1,ZHOOK_HANDLE)
!$acc wait
call system_clock(count_val2)
print *,irank," durée host data ",(count_val2-count_val1)/taux

!$acc kernels
do compteur=1,taille
  dummy2(compteur)=dummy2(compteur)+1.01*dummy1(compteur)
enddo
!$acc end kernels


do compteur=1,taille
  h_send(compteur) = h_send(compteur)+irank
  H_RECV(compteur) = -1.0
enddo

!$ACC DATA COPY (H_RECV, H_SEND)

!$acc kernels
do compteur=1,taille
  dummy2(compteur)=dummy2(compteur)+1.01*dummy1(compteur)
enddo
!$acc end kernels

!$acc wait
call system_clock(count_val1)
print *," "
print *,"======================="
PRINT *, "HOST DATA UPDATE"

CALL MPL_RECV (H_RECV, KSOURCE=IRANKP+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_RECV, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)

CALL MPL_BARRIER ()

!$acc update host(h_send)
CALL MPL_SEND (H_SEND, KDEST=IRANKN+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_SEND, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)

CALL MPL_WAIT (IREQ_RECV)
CALL MPL_WAIT (IREQ_SEND)
CALL MPL_BARRIER()
!$acc update device(h_recv)

print *, irank, " H_RECV ==> ", H_RECV (1:10)
PRINT *, irank, "H_SEND ==> ", H_SEND (1:10)

!$acc wait
call system_clock(count_val2)
print *,irank," durée update host ",(count_val2-count_val1)/taux

!$acc kernels
do compteur=1,taille
  dummy2(compteur)=dummy2(compteur)+1.01*dummy1(compteur)
enddo
!$acc end kernels

PRINT *, IRANK, " H_RECV ==> ", H_RECV (1)
PRINT *,irank, " H_SEND ==> ", H_SEND (1)

do compteur=1,taille
  h_send(compteur) = h_send(compteur)+irank
  H_RECV(compteur) = -1
enddo

!$acc kernels
do compteur=1,taille
  dummy2(compteur)=dummy2(compteur)+1.01*dummy1(compteur)
enddo
!$acc end kernels

!$acc wait
call system_clock(count_val1)

print *," "
print *,"=================================="
PRINT *, " HOST DATA (USE_DEVICE) "

!$ACC HOST_DATA USE_DEVICE (H_RECV)
CALL MPL_RECV (H_RECV, KSOURCE=IRANKP+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_RECV, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)
!$ACC END HOST_DATA

CALL MPL_BARRIER ()

!$ACC HOST_DATA USE_DEVICE (H_SEND)
CALL MPL_SEND (H_SEND, KDEST=IRANKN+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_SEND, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)
!$ACC END HOST_DATA

CALL MPL_WAIT (IREQ_RECV)
CALL MPL_WAIT (IREQ_SEND)

!$acc wait
call system_clock(count_val2)
print *,irank," durée host_data use device ",(count_val2-count_val1)/taux

!$acc kernels
do compteur=1,taille
  dummy2(compteur)=dummy2(compteur)+1.01*dummy1(compteur)
enddo
!$acc end kernels


!$ACC END DATA
!$acc end data

PRINT *, IRANK, " reçu  ==> ", H_RECV (1:10)
PRINT *, IRANK, " envoyé ==> ", H_SEND(1:10)
!!call system_clock(count_val2)
!!print *,irank," durée host_data use device ",(count_val2-count_val1)/taux

DEALLOCATE (C_SEND)
DEALLOCATE (C_RECV)
DEALLOCATE (D_SEND)
DEALLOCATE (D_RECV)
DEALLOCATE (H_SEND)
DEALLOCATE (H_RECV)

PRINT *, "FINALIZE"

CALL MPL_END
 
END PROGRAM MPICUDAAWARE
