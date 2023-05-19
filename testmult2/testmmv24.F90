program testmm
        use cublas
        implicit none
        integer, parameter :: longueur=100
        integer, parameter :: tuile=8
        integer, parameter :: nbvect=100000
        double precision, dimension (tuile,longueur) :: inter2
        double precision, dimension (longueur,tuile)::inter
        double precision, dimension (tuile,int(longueur/tuile)+1) :: inter4
        double precision, allocatable :: a(:,:)
        double precision, allocatable ::v1(:,:),v2(:,:),v2b(:,:)
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer ::i,j,m,comptmult,nbcolonnes,compteurboucle,bloc_x,bloc_y
        integer(kind=8)::compteur

        double precision,dimension(tuile) ::valeur
        double precision::pivot

        allocate(v1(longueur,nbvect),v2(longueur,nbvect),v2b(longueur,nbvect))
        allocate(a(longueur,longueur))

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)

        !$acc data copyin(a,v1) copyout(v2,v2b) create(inter,inter2,valeur)

        do compteur=1,20

        !$acc host_data use_device(a,v1,v2b)
        call cublasDgemm('N','N',longueur,nbvect,longueur,alpha,a,longueur,v1,longueur,beta,v2b,longueur)
        !$acc end host_data
        !$acc wait

        !$acc parallel private(i,j,m,nbcolonnes,compteurboucle,comptmult,inter,inter2,bloc_x,bloc_y) vector_length(32)
        !$acc cache(inter2(tuile,longueur),inter(longueur,tuile))
        !$acc loop gang 
        do i = 1, nbvect,8 
          !!remaining number of columns to compute
          nbcolonnes=min(nbvect-i+1,8)
          !$acc loop vector private(valeur,pivot)
          do j=1,32 
            !!initialize v2
           !! do m=1,nbcolonnes !!v2 tuile
           !!   do bloc_y=1,int(longueur/32)+1
           !!       if (j+32*(bloc_y-1) .le. longueur) v2(j+32*(bloc_y-1),i+m-1)=0.0
           !!   enddo
           !! enddo
            do bloc_x=1,int(longueur/32)+1 !!attention si ne tombe pas juste

              !!load vector data to cache memory
              !$acc loop seq
              do compteurboucle=1,tuile
                if ((32*(bloc_x-1)+j .le. longueur) .and. (compteurboucle .le.nbcolonnes)) then
                  inter(j,compteurboucle)=v1(j+tuile*(bloc_x-1),i+compteurboucle-1)
                else
                  inter(j,compteurboucle)=0.0
                endif
              enddo
            enddo

            do bloc_y=1,int(longueur/tuile)+1 !!attention si ne tombe pas juste
              !!load parameter matrix data to cache memory
              !$acc loop seq
              do compteurboucle=1,32
                if ((tuile*(bloc_y-1)+mod(j,8) .le.longueur).AND.(compteurboucle +(j/8)*32 .le. longueur)) then
                  inter2(j,compteurboucle)=a(mod(j,8)+tuile*(bloc_y-1),compteurboucle+(j/8)*32)
                else
                  inter2(j,compteurboucle)=0.0
                endif
              enddo

              !!for now only 8 threads used
              if ((j .le. 8) .and. (j+8*bloc_y .le. longueur)) then
                !!initialize valeur
         !!       do compteurboucle=1,nbcolonnes
         !!         valeur(compteurboucle)=0.0
         !!       enddo
                
                !!matrix multiplication
                !$acc loop seq
                do m=1,nbcolonnes
                  pivot=0.0
                  !$acc loop seq
                  do comptmult=1,longueur
                    !!valeur(m)=valeur(m)+inter2(j,comptmult)*inter(comptmult,m)
                    pivot=pivot+inter2(j,comptmult)*inter(comptmult,m)
                  enddo
                  valeur(m)=pivot
                enddo    !!matrix multiplication
   
                !!transfer intermediate array valeur to v2 in memory
                do compteurboucle=1,nbcolonnes
                   v2(j+tuile*(bloc_y-1),i-1+compteurboucle)=valeur(compteurboucle)
                enddo
              endif      
            enddo !!bloc_y      
           
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
