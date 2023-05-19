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
        !$acc loop gang vector private(inter1,inter2)
        do i=1,100000
             !!inter1(1:100)=v1(i+j-1,1:100)
             inter1(1:100)=v1(i,1:100)
             inter2 = matmul(at,inter1)
             !!v2(i+j-1,1:100)=inter2(1:100)
             v2(i,1:100)=inter2(1:100)
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
