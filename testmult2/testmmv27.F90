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
        integer :: i,j,m,comptmult,nblignes,compteurboucle,indice
        integer(kind=8)::compteur

        allocate(v1(nbvect,longueur),v2(nbvect,longueur),v2b(nbvect,longueur))
        allocate(a(longueur,longueur))

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)

        !$acc data copyin(a,v1) copyout(v2,v2b) create(inter,inter2,inter4)

        do compteur=1,20

        !$acc host_data use_device(a,v1,v2b)
        call cublasDgemm('N','T',nbvect,longueur,longueur,alpha,v1,nbvect,a,longueur,beta,v2b,nbvect)
        !$acc end host_data
        !$acc wait

        !$acc parallel private(indice,i,j,m,nblignes,compteurboucle,comptmult,inter,inter2) vector_length(32)
        !$acc cache(inter2(longueur,longueur),inter(longueur,longueur))
        !$acc loop gang 
        do i = 1, nbvect,32
          !$acc loop vector private(inter4)
          do j=1,32 
            !!load data to cache memory
            indice=i+j-1
            !$acc loop seq
            do compteurboucle=1,32
              if (indice .le. nbvect) inter(j,compteurboucle)=v1(indice,compteurboucle)
              inter2(j,compteurboucle)=a(j,compteurboucle)
            enddo

            !!matrix multiplication
            if (i+j-1 .le. nbvect) then
              !$acc loop seq
              do m=1,longueur
                 inter4(m)=0.0
                 !$acc loop seq
                 do comptmult=1,32
                   inter4(m)=inter4(m)+inter(j,comptmult)*inter2(m,comptmult)
                 enddo
              enddo
              v2(i+j-1,:)=inter4(:)
            endif
          enddo
        enddo
        !$acc end parallel
        !$acc wait

        enddo

        !$acc end data

        print *, SUM((v2(:,:)-v2b(:,:))*(v2(:,:)-v2b(:,:)))
        print *, v2(1,:)
        print *,v2b(1,:)

        deallocate(a)
        deallocate(v1,v2,v2b)
endprogram
