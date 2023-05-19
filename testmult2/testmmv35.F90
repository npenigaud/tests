program testmm
        use cublas
        implicit none
        integer, parameter :: longueur=100
        integer, parameter :: tuile=32
        integer, parameter :: nbvect=100000
        double precision, dimension (tuile,tuile) :: inter,inter2
        double precision, dimension (tuile,int(longueur/tuile)+1) :: inter4
        double precision, allocatable :: a(:,:)
        double precision, allocatable ::v1(:,:),v2(:,:),v2b(:,:)
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer ::i,j,m,comptmult,nbcolonnes,compteurboucle,bloc_x,bloc_y
        integer(kind=8)::compteur

        double precision,dimension(tuile) ::valeur

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

        !$acc parallel private(i,j,m,nbcolonnes,compteurboucle,comptmult,inter,inter2,bloc_x,bloc_y) vector_length(32)
        !$acc cache(inter2(tuile,tuile),inter(tuile,tuile))
        !$acc loop gang 
        do i = 1, nbvect,32 
          !!initialize values
          nbcolonnes=min(nbvect-i+1,32)
          !$acc loop vector private(inter4,valeur)
          do j=1,32 
            do m=1,nbcolonnes !!v2 tuile
              do bloc_y=1,int(longueur/tuile)+1
   !!v2             inter4(m,bloc_y)=0.0
                  if (j+tuile*(bloc_y-1) .le. longueur) v2(j+tuile*(bloc_y-1),i+m-1)=0.0
              enddo
            enddo
            do bloc_x=1,int(longueur/tuile)+1 !!attention si ne tombe pas juste

              !!load vector data to cache memory
              !$acc loop seq
              do compteurboucle=1,tuile
                if ((tuile*(bloc_x-1)+j .le. longueur) .and. (i+compteurboucle-1 .le. nbvect)) then
                  inter(j,compteurboucle)=v1(j+tuile*(bloc_x-1),i+compteurboucle-1)
                else
                  inter(j,compteurboucle)=0.0
                endif
              enddo

              do bloc_y=1,int(longueur/tuile)+1 !!attention si ne tombe pas juste
                !!load parameter matrix data to cache memory
                !$acc loop seq
                do compteurboucle=1,tuile
                  if ((tuile*(bloc_y-1)+j .le.longueur).AND.(tuile*(bloc_x-1)+compteurboucle .le. longueur)) then
                    inter2(j,compteurboucle)=a(j+tuile*(bloc_y-1),compteurboucle+tuile*(bloc_x-1))
                  else
                    inter2(j,compteurboucle)=0.0
                  endif
                enddo

                !!matrix multiplication
    !!v2            valeur(1:nbcolonnes)=inter4(1:nbcolonnes,bloc_y)
                  do compteurboucle=1,nbcolonnes
                    if (j+tuile*(bloc_y-1) .le. longueur) valeur(compteurboucle)=v2(j+tuile*(bloc_y-1),i-1+compteurboucle)
                  enddo

                !$acc loop seq
                do m=1,nbcolonnes
                  !$acc loop seq
                  do comptmult=1,32
     !!v1               inter4(m,bloc_y)=inter4(m,bloc_y)+inter2(j,comptmult)*inter(comptmult,m)
                      valeur(m)=valeur(m)+inter2(j,comptmult)*inter(comptmult,m)
                  enddo
                enddo    !!matrix multiplication
     !!v2           inter4(1:nbcolonnes,bloc_y)=valeur(1:nbcolonnes)
                  do compteurboucle=1,nbcolonnes
                    if (j+tuile*(bloc_y-1) .le. longueur) v2(j+tuile*(bloc_y-1),i-1+compteurboucle)=valeur(compteurboucle)
                  enddo
              enddo      !!bloc_y
            enddo        !!bloc_x
           
      !! v2     !!move result to v2
            !!do bloc_x=1,int(longueur/tuile)+1
            !!  if (j+tuile*(bloc_x-1) .le. longueur) v2(j+tuile*(bloc_x-1),i:i+nbcolonnes-1)=inter4(1:nbcolonnes,bloc_x)
            !!enddo

          enddo !!thread j
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
