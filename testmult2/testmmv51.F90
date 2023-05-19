program testmm
        use cublas
        implicit none
        integer, parameter :: longueur=32
        integer, parameter :: nbvect=100000
        double precision, dimension(longueur,longueur) :: a,inter,inter2
        double precision,dimension (longueur,nbvect) :: v1,v2b,v2,v3,v4,v5,v6,v7,v8
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer :: i,j,nbcolonnes,k
        integer::compteur

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)
        CALL RANDOM_NUMBER(v3)
        CALL RANDOM_NUMBER(v4)
        CALL RANDOM_NUMBER(V5)
        CALL RANDOM_NUMBER(v6)
        CALL RANDOM_NUMBER(V7)
        CALL RANDOM_NUMBER(v8)

        !$acc data copyin(a,v1,v3,v4,v5,v6,v7,v8) copyout(v2,v2b)

        do compteur=1,20

        !$acc host_data use_device(a,v1,v2)
        call cublasDgemm('N','N',longueur,nbvect,longueur,alpha,a,longueur,v1,longueur,beta,v2b,longueur)
        !$acc end host_data
        !$acc wait

        !$acc parallel private(i,j,nbcolonnes,k) vector_length(32)
        !$acc loop gang
        do i = 1, nbvect,32
          nbcolonnes=min(32,nbvect-i+1)
          !$acc loop vector
          do j=1,32
                if (j .le. nbcolonnes) v2(:, i+j-1) = matmul(a, v1(:,i+j-1))
          enddo

          !$acc loop vector
          do j=1,32
            do k=1,32
              v3(j,i+k-1)=1.0/(10000.0+v3(j,i+k-1))
              v5(j,i+k-1)=sqrt(abs(v5(j,i+k-1)))
              v3(j,i+k-1)=v3(j,i+k-1)+v5(j,i+k-1)
              v7(j,i+k-1)=v3(j,i+k-1)*v7(j,i+k-1)+2.0*v5(j,i+k-1)
              v3(j,i+k-1)=sqrt(abs(v3(j,i+k-1)*v7(j,i+k-1)))
            enddo
          enddo
        enddo
        !$acc end parallel

        !$acc parallel private(i,j,nbcolonnes,k) vector_length(32)
        !$acc loop gang
        do i=1,nbvect,32
          nbcolonnes=min(32,nbvect-i+1)
          !$acc loop vector
          do j=1,32
            do k=1,32
              v4(j,i+k-1)=1.0/(10000.0+v4(k,i+k-1))
              v6(j,i+k-1)=sqrt(abs(v6(j,i+k-1)))
              v4(j,i+k-1)=v4(j,i+k-1)+v6(j,i+k-1)
              v8(j,i+k-1)=v4(j,i+k-1)*v8(j,i+k-1)+2.0*v5(j,i+k-1)
              v4(j,i+k-1)=sqrt(abs(v4(j,i+k-1)*v8(j,i+k-1)))
            enddo
          enddo
        enddo
        !$acc end parallel

        !$acc parallel private(i,j,nbcolonnes,k) vector_length(32)
        !$acc loop gang
        do i = 1, nbvect,32
          nbcolonnes=min(32,nbvect-i+1)
          !$acc loop vector
          do j=1,32
                if (j .le. nbcolonnes) v2(:, i+j-1) = matmul(a, v1(:,i+j-1))
          enddo
        enddo
        !$acc end parallel
        
        end do
        !$acc end data

        print *, SUM((v2(:,:)-v2b(:,:))*(v2(:,:)-v2b(:,:)))
endprogram
