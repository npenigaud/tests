program testmm
        use cutensorex
        implicit none
        real(kind=4), dimension(100,100) :: a,at
        real(kind=4),dimension (100000,100) :: v1,v2b, v2
        real(kind=4), dimension (100) ::inter1,inter2
        real(kind=4), parameter ::alpha=1.0
        real(kind=4), parameter ::beta=0.0
        integer :: i,j,compteur,nblignes

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)
        at=transpose(a)

        !$acc data copyin(a,at,v1) copyout(v2,v2b) create(inter1,inter2)
  
        do compteur=1,20

        !$acc kernels
        !$acc loop independent private(inter1,inter2)
        do i = 1, 100000
                inter1(1:100)=v1(i,1:100)
                inter2 = matmul(at,inter1)
                v2(i,1:100)=inter2(1:100)
                !!v2(i,1:100)=matmul(v1(i,1:100),a)
                !!!!inter2=matmul(v1(i,1:100),a)
                !!!!v2(i,1:100)=inter2(1:100)
        enddo
        !$acc end kernels

        enddo

        !$acc end data

        print *, SUM((v2(:,:)-v2b(:,:))*(v2(:,:)-v2b(:,:)))
        print *,"===="
        print *,"v2 : "
        print *,v2(1,:)
        print *,"v2b : "
        print *,v2b(1,:)
endprogram
