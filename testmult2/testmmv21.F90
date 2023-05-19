program testmm
        use cublas
        use cutensorex
        implicit none
        double precision, dimension(100,100) :: a,at
        double precision,dimension (100000,100) :: v1,v2b, v2
        double precision, dimension (100) ::inter1,inter2
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer :: i,j,compteur,nblignes

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)
        at=transpose(a)

        !$acc data copyin(a,at,v1) copyout(v2,v2b) create(inter1,inter2)
  
        do compteur=1,20

        !$acc host_data use_device(a,v1,v2b)
        call cublasDgemm('N','N',100000,100,100,alpha,v1,100000,a,100,beta,v2b,100000)
        !$acc end host_data
        !$acc wait

        !$acc parallel vector_length(32)
        !$acc loop gang independent private(nblignes)
        !!do i = 1, 100000,32
          do i=0,100000-1,32
           !!nblignes=min(32,100000-i+1)
           nblignes=min(32,100000-i)
           !$acc loop vector independent private(inter1,inter2)
           do j=1,nblignes
                !!inter1(1:100)=v1(i+j-1,1:100)
                inter1(1:100)=v1(i+j,1:100)
                inter2 = matmul(at,inter1)
                !!v2(i+j-1,1:100)=inter2(1:100)
                v2(i+j,1:100)=inter2(1:100)
           enddo
        enddo
        !$acc end parallel

        enddo

        !$acc end data

        print *, SUM((v2(:,:)-v2b(:,:))*(v2(:,:)-v2b(:,:)))
        print *,"===="
        print *,"v2 : "
        print *,v2(1,:)
        print *,"v2b : "
        print *,v2b(1,:)
endprogram
program testmm
        use cublas
        use cutensorex
        implicit none
        double precision, dimension(100,100) :: a,at,b,bt
        double precision,dimension (200000,100) :: v1,v1b,v2b, v2
        double precision,dimension (200000,100) :: v3,v3b,v4,v4b,v5,v5b,v6,v6b
        double precision,dimension (200000)     ::param1,param2
        double precision,dimension(100)         ::param3
        double precision, dimension (100) ::inter1,inter2
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        double precision, parameter ::zdbt=6.9123d0
        double precision            ::in,intermediaire
        integer :: i,j,k,compteur,nblignes

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(b)
        CALL RANDOM_NUMBER(v1)
        CALL RANDOM_NUMBER(v3)
        CALL RANDOM_NUMBER(param1)
        CALL RANDOM_NUMBER(param2)
        CALL RANDOM_NUMBER(param3)
        at=transpose(a)
        bt=transpose(b)
        v1b(:,:)=v1(:,:)
        v3b(:,:)=v3(:,:)

        !$acc data copy(a,at,b,bt,v1,v1b,v3,v3b,param1,param2,param3) copyout(v6,v6b) create(inter1,inter2,v2,v2b,v4,v4b,v5,v5b)
  
        do compteur=1,20
        !$acc kernels
        do i=1,200000
          in=param1(i)
          intermediaire=zdbt*in*param2(i)
          do k=1,100
            v1(i,k)=v3(i,k)-intermediaire*v1(i,k)
          enddo
        enddo
        !$acc end kernels

        !$acc host_data use_device(a,v1,v2)
        call cublasDgemm('N','N',200000,100,100,alpha,v1,200000,a,100,beta,v2,200000)
        !$acc end host_data
        !$acc wait

        !$acc kernels
        do i=1,200000
          do k=1,100
             intermediaire=param3(k)*2.31d0
             v4(i,k)=v2(i,k)/(intermediaire*param1(i))
          enddo
        enddo
        !$acc end kernels

        !$acc host_data use_device(a,v4,v5)
        call cublasDgemm('N','N',200000,100,100,alpha,v4,200000,b,100,beta,v5,200000)
        !$acc end host_data
        !$acc wait
        
        !$acc kernels
        do i=1,200000
          do k=1,100
            v6(i,k)=v5(i,k)*4.107d0*4.107d0
          enddo
        enddo
        !$acc end kernels

        !$acc parallel vector_length(32)
        !$acc loop gang independent private(nblignes,in,intermediaire)
          do i=0,200000-1,32
           nblignes=min(32,200000-i)
           !$acc loop vector independent private(inter1,inter2)
           do j=1,nblignes
                in=param1(i+j)
                intermediaire=zdbt*in*param2(i+j)
                do k=1,100
                  v1b(i+j,k)=v3b(i+j,k)-intermediaire*v1b(i+j,k)
                enddo
                inter1(1:100)=v1b(i+j,1:100)
                inter2 = matmul(at,inter1)
                v2b(i+j,1:100)=inter2(1:100)
                do k=1,100
                  intermediaire=param3(k)*2.31d0
                  v4b(i+j,k)=v2b(i+j,k)/(intermediaire*param1(i+j))
                enddo
                inter1(1:100)=v4b(i+j,1:100)
                inter2=matmul(bt,inter1)
                v5b(i+j,1:100)=inter2(1:100)
                do k=1,100
                  v6b(i+j,k)=v5b(i+j,k)*4.107d0*4.107d0
                enddo
           enddo
        enddo
        !$acc end parallel

        enddo

        !$acc end data

        print *, SUM((v6(:,:)-v6b(:,:))*(v6(:,:)-v6b(:,:)))
        print *,"===="
        print *,"v6 : "
        print *,v6(1,:)
        print *,"v6b : "
        print *,v6b(1,:)
endprogram
