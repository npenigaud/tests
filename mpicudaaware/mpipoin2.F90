PROGRAM MPICUDAAWARE

USE CUDAFOR
USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV, MPL_BARRIER, MPL_INIT, MPL_END, MPL_WAIT, &
                     & JP_NON_BLOCKING_STANDARD

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, ALLOCATABLE :: H_SEND (:,:), H_RECV (:,:)

INTEGER :: IREQ_RECV, IREQ_SEND
INTEGER :: ISTATUS (MPI_STATUS_SIZE)
INTEGER :: ISIZE, IERROR
INTEGER :: IRANK, IRANKP, IRANKN

type tchamp2d
  integer,pointer :: tA(:)=> null(), tB(:,:)=>NULL ()
  character(len=16) :: nom
end type

type tliste
  type(tchamp2d) :: sc1(3)
  integer       :: nombre=0
end type

integer, allocatable, target :: tabl1(:),tabl3(:)
integer, allocatable, target :: tabl2(:,:),tabl4(:,:)
type(tliste)              :: liste
type(tchamp2d)            :: champ
integer                   :: c1,c2,numero

CALL MPL_INIT ()

CALL MPI_COMM_RANK (MPI_COMM_WORLD, IRANK, IERROR)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, ISIZE, IERROR)

ALLOCATE (H_RECV (2,2))
ALLOCATE (H_SEND (2,2))

allocate(tabl1(4))
do c1=1,4
  tabl1(c1)=irank+1
enddo
allocate(tabl2(2,2))
do c1=1,2
  do c2=1,2
    tabl2(c1,c2)=irank+1
  enddo
enddo
!$acc data copy(tabl1) copy(tabl2)

allocate(tabl3(4))
do c1=1,4
  tabl3(c1)=-1
enddo
allocate(tabl4(2,2))
do c1=1,2
  do c2=1,2
    tabl4(c1,c2)=-1
  enddo
enddo
!$acc data copy(tabl3) copy(tabl4)

IRANKP = MODULO (IRANK-1, ISIZE)
IRANKN = MODULO (IRANK+1, ISIZE)

!$acc enter data copyin(liste)
call add2df(liste,tabl1,tabl2,"envoi")
call add2df(liste,tabl3,tabl4,"reception")

!$acc data create(h_recv,h_send)

!$acc kernels
do c1=1,2
  h_send(c1,2)=liste%sc1(1)%tB(c1,2)
  h_recv(c1,1)=liste%sc1(1)%tB(c1,1)
enddo
print *,"h_send 1 carte avant MPI",h_send(1,1)
print *,"h_send 2 carte avant MPI",h_send(1,2)
print *,"h_recv 1 carte avant MPI ",h_recv(1,1)
print *,"h_recv 2 carte avant MPI",h_recv(1,2)
!$acc end kernels

print *,"h_recv host 2 avant MPI avant update",h_recv(1,2)
print *,"h_send host 2 avant MPI avant update",h_send(1,2)
print *,"liste%sc1(2)%tB(1,2) avant MPI avant update",liste%sc1(2)%tB(1,2)
!$acc update host(h_recv,h_send,liste%sc1(2)%tB)
print *,"liste%sc1(2)%tB(1,2) avant MPI après update",liste%sc1(2)%tB(1,2)
print *,"h_recv host 2 avant MPI après update",h_recv(1,2)
print *,"h_send host 2 avant MPI après update",h_send(1,2)

!$ACC HOST_DATA USE_DEVICE (H_RECV(:,2))
CALL MPL_RECV (H_RECV(:,2), KSOURCE=IRANKP+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_RECV, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)
!$ACC END HOST_DATA

CALL MPL_BARRIER ()

!$ACC HOST_DATA USE_DEVICE (H_SEND(:,2))
CALL MPL_SEND (H_SEND(:,2), KDEST=IRANKN+1, KTAG=1001, KCOMM=MPI_COMM_WORLD, KREQUEST=IREQ_SEND, &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD)
!$ACC END HOST_DATA

CALL MPL_WAIT (IREQ_RECV)
CALL MPL_WAIT (IREQ_SEND)

print *,"h_recv host 2 après MPI avant update",h_recv(1,2)
print *,"liste%sc1(2)%tB(1,2) après MPI avant update",liste%sc1(2)%tB(1,2)
!$acc update host(h_recv,liste%sc1(2)%tB)
print *,"h_recv host 2 après MPI après update",h_recv(1,2)
print *,"liste%sc1(2)%tB(1,2) après MPI après update",liste%sc1(2)%tB(1,2)

!$acc kernels
do c1=1,2
  do c2=1,2
    liste%sc1(2)%tB(c1,c2)=h_recv(c1,c2)
  enddo
enddo
print *,"liste sc1(1) pspg ==> ",liste%sc1(1)%tB(1,2)
print *,"liste sc1(2) pspg ==> ",liste%sc1(2)%tB(1,2)
print *,"tabl2 1 sur carte ==> ",tabl2(1,1)
print *,"tabl2 2 sur carte ==> ",tabl2(1,2)

print *,"tabl4 1 sur carte ==> ",tabl4(1,1)
print *,"tabl2 2 sur carte ==> ",tabl2(1,2)
print *,"tabl4 1 sur carte ==> ",tabl4(1,1)

print *,"tabl4 2 sur carte ==> ",tabl4(1,2)
print *,"h_recv carte 2",h_recv(1,2)
print *,"h_send carte 2",h_send(1,2)
!$acc end kernels

!$ACC END DATA
!$acc end data
!$acc end data

PRINT *, IRANK, "tabl2 ==> ", tabl2(:,1),tabl2(:,2)
print *, irank, "tabl4 ==> ", tabl4(:,1),tabl4(:,2)

DEALLOCATE (H_SEND)
DEALLOCATE (H_RECV)

PRINT *, "FINALIZE"

CALL MPL_END
 
contains

SUBROUTINE ADD2DF (YDLIST, PSP, PSPG,cdname)

TYPE (tliste)                                 :: YDLIST
integer                ,OPTIONAL ,TARGET :: PSP (:)
integer, optional, target                ::  PSPG (:,:)
CHARACTER(LEN=*)  ,INTENT(IN)                    :: CDNAME

YDLIST%nombre = YDLIST%nombre + 1
!$acc update device(YDLIST%nombre)

ASSOCIATE (sc1 => YDLIST%sc1 (YDLIST%nombre))
  sc1%nom = CDNAME
  sc1%tA   => PSP
  sc1%tB  => PSPG
!$acc enter data copyin(sc1) attach(sc1%tA,sc1%tB)
!$acc update device(sc1%nom)

END ASSOCIATE

END SUBROUTINE


END PROGRAM MPICUDAAWARE
