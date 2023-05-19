program pointeur

implicit none

type tchamp2d
  real,pointer :: tA(:)=> null(), tB(:,:)=>NULL ()
end type

type tliste
  type(tchamp2d) :: sc1(3)
  integer       :: nombre=0
end type

real, allocatable, target :: tabl1(:)
real, allocatable, target :: tabl2(:,:)
type(tliste)              :: liste
type(tchamp2d)            :: champ
integer                   :: c1,c2

allocate(tabl1(4))
tabl1=(/1,2,3,4/)
allocate(tabl2(2,2))
!$acc data copyin(tabl1) copyout(tabl2)

!$acc kernels
tabl1(3)=2.1
tabl1(2)=3.1
!$acc end kernels

!$acc enter data copyin(liste)
call add2df(liste,tabl1,tabl2)

!liste%nombre=liste%nombre+1
!!$acc update device(liste%nombre)
!!liste%sc1(liste%nombre)%tA=>tabl1
!!liste%sc1(liste%nombre)%tB=>tabl2

!!$acc enter data copyin(liste%sc1(1)) attach(liste%sc1(1)%tA,liste%sc1(1)%tB)


!$acc kernels
do c1=1,2
  do c2=1,2
    liste%sc1(1)%tB(c1,c2)=liste%sc1(1)%tA(2*(c1-1)+c2)
  enddo
enddo
print *,liste%nombre
!$acc end kernels
!$acc end data
!!$acc exit data delete(liste%sc1(1))
!!$acc exit data delete(liste)

print *, tabl2

contains

SUBROUTINE ADD2DF (YDLIST, PSP, PSPG)!!,cdname,ldbcast)

TYPE (tliste)                                 :: YDLIST
REAL                ,OPTIONAL ,TARGET :: PSP (:)
real, optional, target                ::  PSPG (:,:)
!!CHARACTER(LEN=*)  ,INTENT(IN)                    :: CDNAME
!!LOGICAL           ,INTENT(IN)  ,OPTIONAL         :: LDBCAST

!!!$acc data present(psp,pspg)

YDLIST%nombre = YDLIST%nombre + 1
!$acc update device(YDLIST%nombre)

ASSOCIATE (sc1 => YDLIST%sc1 (YDLIST%nombre))
  !!sc1%CNAME = CDNAME
  sc1%tA   => PSP
  sc1%tB  => PSPG
!!$acc enter data copyin(ydlist%sc1(ydlist%nombre)) attach(ydlist%sc1(ydlist%nombre)%tA,ydlist%sc1(ydlist%nombre)%tB)
!$acc enter data copyin(sc1) attach(sc1%tA,sc1%tB)

  !!IF (PRESENT (LDBCAST)) YL2D%LBCAST = LDBCAST
END ASSOCIATE

!!!$acc end data
END SUBROUTINE

end program 





