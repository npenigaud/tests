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

liste%nombre=liste%nombre+1
liste%sc1(liste%nombre)%tA=>tabl1
liste%sc1(liste%nombre)%tB=>tabl2

!$acc data copyin(liste,liste%sc1(1)) attach(liste%sc1(1)%tA,liste%sc1(1)%tB)
!$acc kernels
do c1=1,2
  do c2=1,2
    liste%sc1(1)%tB(c1,c2)=liste%sc1(1)%tA(2*(c1-1)+c2)
  enddo
enddo
print *,liste%nombre
!$acc end kernels
!$acc end data
!$acc end data

print *, tabl2

contains

SUBROUTINE ADD2DF (YDLIST, PSP, PSPG)!!,cdname,ldbcast)

TYPE (tliste)                                 :: YDLIST
REAL                ,OPTIONAL ,TARGET :: PSP (:)
real, optional, target                ::  PSPG (:,:)
!!CHARACTER(LEN=*)  ,INTENT(IN)                    :: CDNAME
!!LOGICAL           ,INTENT(IN)  ,OPTIONAL         :: LDBCAST

IF (.NOT. PRESENT (PSP)) GOTO 999

IF (PRESENT(PSP) .NEQV. PRESENT (PSPG)) THEN
  print *,"erreur, les deux tableaux doivent être présents"
ENDIF

YDLIST%nombre = YDLIST%nombre + 1

IF (YDLIST%nombre > SIZE (YDLIST%sc1)) THEN
  print *,"la liste est pleine"
ENDIF

ASSOCIATE (sc1 => YDLIST%sc1 (YDLIST%nombre))
  !!sc1%CNAME = CDNAME
  sc1%tA   => PSP
  sc1%tB  => PSPG
  !!IF (PRESENT (LDBCAST)) YL2D%LBCAST = LDBCAST
END ASSOCIATE

999 CONTINUE

END SUBROUTINE

end program 





