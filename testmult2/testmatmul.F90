program testmm
        implicit none
        double precision :: a(3, 3)
        double precision :: v1(3, 2), v2(3, 2)
        integer :: i

        a = reshape([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], [3, 3])
        v1 = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [3, 2])

        print *, v1

        !$acc kernels copyin(a, v1) copyout(v2)
        !$acc loop independent
        do i = 1, 2
                v2(:, i) = matmul(a, v1(:, i))
        enddo
        !$acc end kernels

        print *, v2
endprogram
