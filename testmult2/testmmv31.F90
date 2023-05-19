program testmm
        use cublas
        implicit none
        double precision ::coeff
        double precision, dimension(100,100) :: a,b
        double precision,dimension (100,100000) :: v1,v2b, v2
        double precision, parameter ::alpha=1.0
        double precision, parameter ::beta=0.0
        integer :: i,j,m,n,k,ligne,colonne,tailleligne,taillecolonne
        integer(kind=8)::nb_periodes_initial1,nb_periodes_initial2
        integer(kind=8)::nb_periodes_final1,nb_periodes_final2
        integer(kind=8)::nb_periodes_max,nb_periodes_sec,nb_periodes,temps_passe
        integer(kind=8)::compteur

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(v1)

        b(:,:)=a(:,:)
        v2(:,:)=0.0

        !$acc data copyin(a,v1) copyout(v2,v2b)

        do compteur=1,20

        call system_clock(COUNT_RATE=nb_periodes_sec,COUNT_MAX=nb_periodes_max)
        call system_clock(COUNT=nb_periodes_initial1)


        !$acc host_data use_device(a,v1,v2)
        call cublasDgemm('N','N',100,100000,100,alpha,a,100,v1,100,beta,v2b,100)
        !$acc end host_data
        !$acc wait

        call system_clock(COUNT=nb_periodes_final1)
        nb_periodes=nb_periodes_final1-nb_periodes_initial1
        if (nb_periodes_final1<nb_periodes_initial1) nb_periodes=nb_periodes+nb_periodes_max
        temps_passe=REAL(nb_periodes)/nb_periodes_sec
        print *, "version1, duree (s) : ",temps_passe

        call system_clock(COUNT=nb_periodes_initial2)

        !$acc parallel vector_length(32) private(tailleligne,taillecolonne,i,j,k,m,n,ligne,colonne,coeff,a)
        !$acc loop gang collapse(2) 
        do m=1,100000,64
          do n=1,100,128 
            !!!$acc cache(a(1:100,1:100))
            !$acc loop seq
            do j=1,128
              ligne=n+j-1
              !$acc loop seq
              do k=1,100
                coeff=a(ligne,k)
                !$acc loop vector
                do i=1,64 
                  colonne=m+i-1
                  if ((colonne.LE.100000).AND.(ligne .LE. 100)) v2(ligne,colonne)=v2(ligne,colonne)+a(ligne,k)*v1(k,colonne)
                enddo
              enddo
            enddo
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
endprogram
