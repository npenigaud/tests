program testmm
        use cublas
        implicit none
        integer, parameter :: longueur=32
        integer, parameter :: nbvect=100000
        double precision, dimension(longueur,longueur) :: inter,inter2
        double precision, dimension(longueur) :: inter4
        double precision, allocatable :: a(:,:)
        double precision, allocatable ::v1(:,:),v2(:,:),v2b(:,:)
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer :: i,j,m,comptmult,nbcolonnes,compteurboucle
        integer(kind=8)::compteur

        allocate(v1(longueur,nbvect),v2(longueur,nbvect),v2b(longueur,nbvect))
        allocate(a(longueur,longueur))

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)

        !$acc data copyin(a,v1) copyout(v2,v2b) create(inter,inter2,inter4)

        do compteur=1,20

        !$acc host_data use_device(a,v1,v2b)
        call cublasDgemm('N','N',longueur,nbvect,longueur,alpha,a,longueur,v1,longueur,beta,v2b,longueur)
        !$acc end host_data
        !$acc wait

        !$acc parallel private(i,j,m,nbcolonnes,compteurboucle,comptmult,inter,inter2) vector_length(32)
        !$acc cache(inter2(longueur,longueur),inter(longueur,longueur))
        !$acc loop gang 
        do i = 1, nbvect,32
          !$acc loop vector private(inter4)
          do j=1,32 
            !!load data to cache memory
            !$acc loop seq
            do compteurboucle=1,32
              inter(j,compteurboucle)=v1(j,i+compteurboucle-1)
              inter2(j,compteurboucle)=a(j,compteurboucle)
            enddo

            !!matrix multiplication
            nbcolonnes=min(nbvect-i+1,32)
            !$acc loop seq
            do m=1,nbcolonnes
              inter4(m)=0.0
              !$acc loop seq
              do comptmult=1,32
                inter4(m)=inter4(m)+inter2(j,comptmult)*inter(comptmult,m)
              enddo
            enddo          
            v2(j,i:i+nbcolonnes-1)=inter4(1:nbcolonnes)
          enddo
        enddo
        !$acc end parallel
        !$acc wait

        enddo

        !$acc end data

        print *, SUM((v2(:,:)-v2b(:,:))*(v2(:,:)-v2b(:,:)))
        print *, v2(:,1)
        print *,v2b(:,1)

        deallocate(a)
        deallocate(v1,v2,v2b)
endprogram
