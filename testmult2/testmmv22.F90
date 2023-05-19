program testmm
        use cublas
        implicit none
        integer, parameter :: longueur=32
        double precision, dimension(longueur,longueur) :: a,inter,inter2
        double precision,dimension (longueur,100000) :: v1,v2b, v2
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer :: i,indice,j
        integer(kind=8)::nb_periodes_initial1,nb_periodes_initial2
        integer(kind=8)::nb_periodes_final1,nb_periodes_final2
        integer(kind=8)::nb_periodes_max,nb_periodes_sec,nb_periodes,temps_passe
        integer(kind=8)::compteur

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)

        !$acc data copyin(a,v1) copyout(v2,v2b) create(inter,inter2)

        do compteur=1,20

        call system_clock(COUNT_RATE=nb_periodes_sec,COUNT_MAX=nb_periodes_max)
        call system_clock(COUNT=nb_periodes_initial1)


        !$acc host_data use_device(a,v1,v2b)
        call cublasDgemm('N','N',longueur,100000,longueur,alpha,a,longueur,v1,longueur,beta,v2b,longueur)
        !$acc end host_data
        !$acc wait

        call system_clock(COUNT=nb_periodes_final1)
        nb_periodes=nb_periodes_final1-nb_periodes_initial1
        if (nb_periodes_final1<nb_periodes_initial1) nb_periodes=nb_periodes+nb_periodes_max
        temps_passe=REAL(nb_periodes)/nb_periodes_sec
        print *, "version1, duree (s) : ",temps_passe

        call system_clock(COUNT=nb_periodes_initial2)

        !$acc parallel private(i,j,indice,inter,inter2) vector_length(32)
        !$acc cache(inter2(longueur,longueur),inter(longueur,32))
        !$acc loop gang 
        do i = 1, 100000,32
          !!!!$acc cache(inter(longueur,32))
          !$acc loop vector
          do j=1,32
                indice=i+j-1
                inter(j,:)=v1(j,i:i+32)
                inter2(j,:)=a(j,:)
                if (indice .le.100000) v2(:, indice) = matmul(inter2, inter(:,j))
          enddo
        enddo
        !$acc end parallel
       

        call system_clock(COUNT=nb_periodes_final2)
        nb_periodes=nb_periodes_final2-nb_periodes_initial2
        if (nb_periodes_final2<nb_periodes_initial2) nb_periodes=nb_periodes+nb_periodes_max
        temps_passe=REAL(nb_periodes)/nb_periodes_sec
        print *, "version2, duree (s) : ",temps_passe

        enddo

        !$acc end data

        print *, SUM((v2(:,:)-v2b(:,:))*(v2(:,:)-v2b(:,:)))
        print *, v2(:,1)
        print *,v2b(:,1)
endprogram
