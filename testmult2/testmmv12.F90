program testmm
        use cublas
        implicit none
        integer, parameter :: longueur=32
        double precision :: valeur
        double precision, dimension(longueur,longueur) :: a,inter2
        double precision, dimension(longueur,longueur)::inter
        double precision, dimension(longueur)::inter4
        double precision,dimension (longueur,100000) :: v1,v2b, v2
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer :: i,j,m,comptmult,nbcolonnes,indice,compteurboucle
        integer(kind=8)::nb_periodes_initial1,nb_periodes_initial2
        integer(kind=8)::nb_periodes_final1,nb_periodes_final2
        integer(kind=8)::nb_periodes_max,nb_periodes_sec,nb_periodes,temps_passe
        integer(kind=8)::compteur

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)

        !$acc data copyin(a,v1) copyout(v2,v2b) create(inter,inter2)

        do compteur=1,20

        call system_clock(COUNT_RATE=nb_periodes_sec,COUNT_MAX=nb_periodes_max)
        call system_clock(COUNT=nb_periodes_initial1)


        !$acc host_data use_device(a,v1,v2b)
        call cublasDgemm('N','N',longueur,100000,longueur,alpha,a,longueur,v1,longueur,beta,v2b,longueur)
        !$acc end host_data
        !$acc wait

        call system_clock(COUNT=nb_periodes_final1)
        nb_periodes=nb_periodes_final1-nb_periodes_initial1
        if (nb_periodes_final1<nb_periodes_initial1) nb_periodes=nb_periodes+nb_periodes_max
        temps_passe=REAL(nb_periodes)/nb_periodes_sec
        print *, "version1, duree (s) : ",temps_passe

        call system_clock(COUNT=nb_periodes_initial2)

        !$acc parallel private(valeur,i,j,m,nbcolonnes,compteurboucle,comptmult,indice,inter,inter2) vector_length(32)
        !$acc cache(inter2(longueur,longueur),inter(longueur,32))
        !$acc loop gang 
        do i = 1, 100000,32
          !$acc loop vector private(inter4)
          do j=1,32 
          !$acc loop seq
          do compteurboucle=1,32
            inter(j,compteurboucle)=v1(j,i+compteurboucle-1)
            inter2(j,compteurboucle)=a(j,compteurboucle)
          enddo

          nbcolonnes=min(100000-i+1,32)
          !$acc loop seq
          do m=1,nbcolonnes
             indice=i+m-1
             valeur=0.0
             !$acc loop seq
             do comptmult=1,32
                valeur=valeur+inter2(j,comptmult)*inter(comptmult,m)
             enddo
             inter4(m)=valeur
          enddo          
          v2(j,i:i+nbcolonnes-1)=inter4(1:nbcolonnes)
          enddo
        enddo
        !$acc end parallel
       

        call system_clock(COUNT=nb_periodes_final2)
        nb_periodes=nb_periodes_final2-nb_periodes_initial2
        if (nb_periodes_final2<nb_periodes_initial2) nb_periodes=nb_periodes+nb_periodes_max
        temps_passe=REAL(nb_periodes)/nb_periodes_sec
        print *, "version2, duree (s) : ",temps_passe

        enddo

        !$acc end data

        print *, SUM((v2(:,:)-v2b(:,:))*(v2(:,:)-v2b(:,:)))
        print *, v2(:,1)
        print *,v2b(:,1)
endprogram
subroutine tm4sprou(indiceSpec,indiceLev)
        use cublas
        use cutensorex
        implicit none
integer(kind=4),intent(in) :: indiceSpec
integer(kind=4),intent(in) :: indiceLev
#define dynamique 0
        !!integer(kind=4), parameter  :: indiceSpec = 250500 
        !!integer(kind=4), parameter  :: indiceLev  = 105     
#if dynamique
        real(kind=4), allocatable  :: a(:,:),at(:,:),b(:,:),bt(:,:)
        real(kind=4), allocatable  :: v1(:,:),v2(:,:)
        real(kind=4), allocatable  :: v3(:,:),v4(:,:),v5(:,:),v6(:,:)
#else
        real(kind=4), dimension(indiceLev,indiceLev) ::a,at,b,bt
        real(kind=4), dimension(indiceSpec,indiceLev)::v1,v2
        real(kind=4), dimension(indiceSpec,indiceLev)::v3,v4,v5,v6
#endif
        real(kind=4), dimension (indiceSpec)     :: rlapdi,tech1
        integer(kind=4) , dimension (indiceSpec)     :: nvalue
        real(kind=4), dimension (indiceLev)         :: sivp
        real(kind=4), dimension (indiceLev) :: inter1,inter2
        real(kind=4), parameter :: alpha=1.0
        real(kind=4), parameter :: beta=0.0
        real(kind=4), parameter :: zdbt=0.69123d0
        real(kind=4), parameter :: const2=0.231d0
        real(kind=4), parameter :: const3=0.4107d0
        integer(kind=4)             :: indice
        real(kind=4)            :: intermediaire
        integer(kind=4) :: jsp,jlev,compteur,nblignes

#if dynamique
        print *,"allocation dynamique"
        allocate(a(indiceLev,indiceLev))
        allocate(at(indiceLev,indiceLev))
        allocate(b(indiceLev,indiceLev))
        allocate(bt(indiceLev,indiceLev))
        allocate(v1(indiceSpec,indiceLev))
        allocate(v2(indiceSpec,indiceLev))
        allocate(v3(indiceSpec,indiceLev))
        allocate(v4(indiceSpec,indiceLev))
        allocate(v5(indiceSpec,indiceLev))
        allocate(v6(indiceSpec,indiceLev))
#else
        print *,"allocation statique"
#endif

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(b)
        CALL RANDOM_NUMBER(v1)
        CALL RANDOM_NUMBER(v3)
        CALL RANDOM_NUMBER(rlapdi)
        CALL RANDOM_NUMBER(tech1)
        do jsp=1,indiceSpec
          nvalue(jsp)=1+floor(indiceSpec*tech1(jsp))
        enddo
        CALL RANDOM_NUMBER(sivp)
        at=transpose(a/50.0d0)
        bt=transpose(b/50.0d0)

        !$acc data copy(at,bt,v1,v3,nvalue,sivp,rlapdi) copyout(v6) create(inter1,inter2,v2,v4,v5)
  
        do compteur=1,20
        !$acc parallel private(nblignes,indice,intermediaire) vector_length(32) default(none)
        !$acc loop gang vector private(inter1,inter2)
          do jsp=1,indiceSpec
            indice=nvalue(jsp)
            intermediaire=zdbt*rlapdi(indice)
            do jlev=1,indiceLev
              v1(jsp,jlev)=v3(jsp,jlev)-intermediaire*v1(jsp,jlev)
            enddo
            inter1(1:indiceLev)=v1(jsp,1:indiceLev)
            inter2 = matmul(at,inter1)
            v2(jsp,1:indiceLev)=inter2(1:indiceLev)
            do jlev=1,indiceLev
              intermediaire=sivp(jlev)*const2
              v4(jsp,jlev)=v2(jsp,jlev)/(intermediaire*rlapdi(nvalue(jsp)))
              !!v4b(jsp,jlev)=v2b(jsp,jlev)*(intermediaire*rlapdi(nvalue(jsp)))
            enddo
            inter1(1:indiceLev)=v4(jsp,1:indiceLev)
            inter2=matmul(bt,inter1)
            v5(jsp,1:indiceLev)=inter2(1:indiceLev)
            do jlev=1,indiceLev
              v6(jsp,jlev)=v5(jsp,jlev)*const3*const3
            enddo
        enddo
        !$acc end parallel

        enddo

        !$acc end data
    
        print *,v6(1,:)   
  
#if dynamique
deallocate(a,b,at,bt,v1,v2,v3,v4,v5,v6)
#endif

end subroutine
