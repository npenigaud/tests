PROGRAM MPICUDAAWARE

USE CUDAFOR
USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV, MPL_BARRIER, MPL_INIT, MPL_END, MPL_WAIT, &
                     & JP_NON_BLOCKING_STANDARD

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, ALLOCATABLE :: H_SEND (:), H_RECV (:)

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

ALLOCATE (H_RECV (4))
ALLOCATE (H_SEND (4))

allocate(tabl1(4))
do c1=1,4
  tabl1(c1)=irank
enddo
allocate(tabl2(2,2))
!$acc data copyin(tabl1) copyout(tabl2)

allocate(tabl3(4))
do c1=1,4
  tabl3(c1)=-1
enddo
allocate(tabl4(2,2))
!!$acc data copyin(tabl3) copyout(tabl4)
!$acc data copy(tabl3) copyout(tabl4)

IRANKP = MODULO (IRANK-1, ISIZE)
IRANKN = MODULO (IRANK+1, ISIZE)


!$acc data create(h_recv,h_send)

!$acc kernels
h_send(1:4)=tabl1(1:4)
!$acc end kernels

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

!$acc kernels
tabl3(1:4)=h_recv(1:4)
!$acc end kernels

!$ACC END DATA
!$acc end data
!$acc end data

PRINT *, IRANK, "tabl1 ==> ", tabl1(:)
print *, irank, "tabl3 ==> ", tabl3(:)

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
