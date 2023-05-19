program pointeur

implicit none

type tchamp2d
  real,pointer :: tA(:)=> null(), tB(:,:)=>NULL ()
  character(len=16) :: nom
end type

type tliste
  type(tchamp2d) :: sc1(3)
  integer       :: nombre=0
end type

real, allocatable, target :: tabl1(:),tabl3(:)
real, allocatable, target :: tabl2(:,:),tabl4(:,:)
type(tliste)              :: liste
type(tchamp2d)            :: champ
integer                   :: c1,c2,numero

allocate(tabl1(4))
tabl1=(/1,2,3,4/)
allocate(tabl2(2,2))
!$acc data copyin(tabl1) copyout(tabl2)

allocate(tabl3(4))
tabl3=(/9.5,8.5,7.5,6.5/)
allocate(tabl4(2,2))
!$acc data copyin(tabl3) copyout(tabl4)


!$acc kernels
tabl1(3)=2.1
tabl1(2)=3.1
!$acc end kernels

!$acc enter data copyin(liste)
call add2df(liste,tabl1,tabl2,"t1")
call add2df(liste,tabl3,tabl4,"t2")

!$acc kernels
do numero=1,liste%nombre
print *,"traitement",liste%sc1(numero)%nom
do c1=1,2
  do c2=1,2
    liste%sc1(numero)%tB(c1,c2)=liste%sc1(numero)%tA(2*(c1-1)+c2)
  enddo
enddo
enddo
print *,liste%nombre
!$acc end kernels
!$acc end data
!$acc end data
!!$acc exit data delete(liste%sc1(1))
!!$acc exit data delete(liste)

print *, tabl2
print *,tabl4

contains

SUBROUTINE ADD2DF (YDLIST, PSP, PSPG,cdname)

TYPE (tliste)                                 :: YDLIST
REAL                ,OPTIONAL ,TARGET :: PSP (:)
real, optional, target                ::  PSPG (:,:)
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

end program 





