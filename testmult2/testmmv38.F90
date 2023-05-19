program testmm
        use cublas
        implicit none
        integer, parameter :: longueur=32
        integer, parameter :: nbvect=100000
        double precision, dimension(longueur,longueur) :: inter2
        double precision, dimension(2*longueur,longueur)::inter
        double precision, dimension(longueur) :: inter4,inter5
        double precision, allocatable :: a(:,:)
        double precision, allocatable ::v1(:,:),v2(:,:),v2b(:,:)
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer ::i,j,m,comptmult,nblignes,compteurboucle
        integer::indice,indice2,jp32
        integer(kind=8)::compteur
        logical::bool1,bool2

        allocate(v1(nbvect,longueur),v2(nbvect,longueur),v2b(nbvect,longueur))
        allocate(a(longueur,longueur))

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)

        !$acc data copyin(a,v1) copyout(v2,v2b) create(inter,inter2,inter4,inter5)
        do compteur=1,20

        !$acc host_data use_device(a,v1,v2b)
        call cublasDgemm('N','T',nbvect,longueur,longueur,alpha,v1,nbvect,a,longueur,beta,v2b,nbvect)
        !$acc end host_data
        !$acc wait

        !$acc parallel private(bool1,bool2,indice,indice2,i,j,jp32,m,nblignes,compteurboucle,comptmult,inter,inter2) vector_length(32)
        !$acc cache(inter2(longueur,longueur),inter(longueur,longueur))
        !$acc loop gang 
        do i = 1, nbvect,32*2
          !$acc loop vector private(inter4,inter5)
          do j=1,32 
            !!load data to cache memory
            indice=i+j-1
            indice2=indice+32
            jp32=j+32
            bool1=indice.le.nbvect
            bool2=indice2.le.nbvect
            !$acc loop seq
            do compteurboucle=1,32
              if (bool1) inter(j,compteurboucle)=v1(indice,compteurboucle)
              if (bool2) inter(jp32,compteurboucle)=v1(indice2,compteurboucle)
              if (j.le.32) inter2(j,compteurboucle)=a(j,compteurboucle)
            enddo

            !!matrix multiplication
            if (bool1) then
              !$acc loop seq
              do m=1,longueur
                 inter4(m)=0.0
                 if (bool2) inter5(m)=0.0
                 !$acc loop seq
                 do comptmult=1,32
                   inter4(m)=inter4(m)+inter(j,comptmult)*inter2(m,comptmult)
                   if (bool2) inter5(m)=inter5(m)+inter(jp32,comptmult)*inter2(m,comptmult)
                 enddo
              enddo
              v2(indice,:)=inter4(:)
              if (bool2) v2(indice2,:)=inter5(:)
            endif
            !!matrix multiplication
!            if (bool2) then 
!              !$acc loop seq
!              do m=1,longueur
!                 inter4(m)=0.0
!                 !$acc loop seq
!                 do comptmult=1,32
!                   inter4(m)=inter4(m)+inter(jp32,comptmult)*inter2(m,comptmult)
!                 enddo
!              enddo
!              v2(indice2,:)=inter4(:)
!            endif
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
