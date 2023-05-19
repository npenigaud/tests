program testmm
        use cublas
        implicit none
        integer, parameter :: longueur1=16
        integer, parameter :: longueur2=30
        integer, parameter :: tuile=32
        integer, parameter :: nbvect=100000
        double precision, dimension(tuile,tuile) :: inter,inter2
        double precision, dimension(tuile) :: inter4
        double precision, allocatable :: a(:,:)
        double precision, allocatable ::v1(:,:),v2(:,:),v2b(:,:)
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer :: i,j,m,comptmult,nbcolonnes,compteurboucle
        integer::compteur,longueur

        allocate(v1(longueur2,nbvect),v2(longueur1,nbvect),v2b(longueur1,nbvect))
        allocate(a(longueur1,longueur2))

        longueur=max(longueur1,longueur2)

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)

        !$acc data copyin(a,v1) copyout(v2,v2b) create(inter,inter2,inter4)

        do compteur=1,20

        !$acc host_data use_device(a,v1,v2b)
        call cublasDgemm('N','N',longueur1,nbvect,longueur2,alpha,a,longueur1,v1,longueur2,beta,v2b,longueur1)
        !$acc end host_data
        !$acc wait

        !$acc parallel private(i,j,m,nbcolonnes,compteurboucle,comptmult,inter,inter2) vector_length(32)
        !$acc cache(inter2(tuile,tuile),inter(tuile,tuile))
        !$acc loop gang 
        do i = 1, nbvect,longueur
          !$acc loop vector private(inter4)
          do j=1,longueur
            !!load data to cache memory
            !$acc loop seq
            do compteurboucle=1,longueur2
              if (j.le.longueur2) inter(j,compteurboucle)=v1(j,i+compteurboucle-1)
              if (j.le.longueur1) inter2(j,compteurboucle)=a(j,compteurboucle)
            enddo

            !!matrix multiplication
            nbcolonnes=min(nbvect-i+1,longueur)
            if (j.le.longueur1) then
              !$acc loop seq
              do m=1,nbcolonnes
                inter4(m)=0.0
                !$acc loop seq
                do comptmult=1,longueur2
                  inter4(m)=inter4(m)+inter2(j,comptmult)*inter(comptmult,m)
                enddo
              enddo          
              v2(j,i:i+nbcolonnes-1)=inter4(1:nbcolonnes)
            endif
          enddo
        enddo
        !$acc end parallel
        !$acc wait

        enddo

        !$acc end data

        print *, SUM((v2(:,:)-v2b(:,:))*(v2(:,:)-v2b(:,:)))
        do compteur=1,32
          print *,"colonne ",compteur," ; matrice cache : "
          print *, v2(:,compteur)
          print *,"=== matrice DGEMM : "
          print *,v2b(:,compteur)
        enddo

        deallocate(a)
        deallocate(v1,v2,v2b)
endprogram
