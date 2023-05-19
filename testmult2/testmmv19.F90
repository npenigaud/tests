program testmm
        use cublas
        use cutensorex
        implicit none
#define dynamique 1
        integer(kind=4), parameter  :: indiceSpec = 250500   !!250500
        integer(kind=4), parameter  :: indiceLev  = 105      !!105
#if dynamique
        double precision, allocatable  :: a(:,:),at(:,:),b(:,:),bt(:,:)
        double precision, allocatable  :: v1(:,:),v1b(:,:),v2b(:,:), v2(:,:)
        double precision, allocatable  :: v3(:,:),v3b(:,:),v4(:,:),v4b(:,:),v5(:,:),v5b(:,:),v6(:,:),v6b(:,:)
#else
        double precision, dimension(indiceLev,indiceLev) ::a,at,b,bt
        double precision, dimension(indiceSpec,indiceLev)::v1,v1b,v2,v2b
        double precision, dimension(indiceSpec,indiceLev)::v3,v3b,v4,v4b,v5,v5b,v6,v6b
#endif
        double precision, dimension (indiceSpec)     :: rlapdi,tech1
        integer(kind=4) , dimension (indiceSpec)     :: nvalue
        double precision, dimension (indiceLev)         :: sivp
        double precision, dimension (indiceLev) :: inter1,inter2
        double precision, parameter :: alpha=1.0
        double precision, parameter :: beta=0.0
        double precision, parameter :: zdbt=6.9123d0
        integer(kind=4)             :: in
        double precision            :: intermediaire
        integer(kind=4) :: jsp,jlev,compteur,nblignes

#if dynamique
        print *,"allocation dynamique"
        allocate(a(indiceLev,indiceLev))
        allocate(at(indiceLev,indiceLev))
        allocate(b(indiceLev,indiceLev))
        allocate(bt(indiceLev,indiceLev))
        allocate(v1(indiceSpec,indiceLev))
        allocate(v2(indiceSpec,indiceLev))
        allocate(v3(indiceSpec,indiceLev))
        allocate(v4(indiceSpec,indiceLev))
        allocate(v5(indiceSpec,indiceLev))
        allocate(v6(indiceSpec,indiceLev))
        allocate(v1b(indiceSpec,indiceLev))
        allocate(v2b(indiceSpec,indiceLev))
        allocate(v3b(indiceSpec,indiceLev))
        allocate(v4b(indiceSpec,indiceLev))
        allocate(v5b(indiceSpec,indiceLev))
        allocate(v6b(indiceSpec,indiceLev))
#else
        print *,"allocation statique"
#endif

        CALL RANDOM_NUMBER(a)
        CALL RANDOM_NUMBER(b)
        CALL RANDOM_NUMBER(v1)
        CALL RANDOM_NUMBER(v3)
        CALL RANDOM_NUMBER(rlapdi)
        CALL RANDOM_NUMBER(tech1)
        do jsp=1,indiceSpec
          nvalue(jsp)=1+floor(indiceSpec*tech1(jsp))
        enddo
        CALL RANDOM_NUMBER(sivp)
        at=transpose(a)
        bt=transpose(b)
        v1b(:,:)=v1(:,:)
        v3b(:,:)=v3(:,:)

        !!do compteur=1,20
        !$acc data copy(a,at,b,bt,v1,v1b,v3,v3b,nvalue,sivp,rlapdi) copyout(v6,v6b) create(inter1,inter2,v2,v2b,v4,v4b,v5,v5b)
  
        do compteur=1,20
        !$acc kernels
        do jsp=1,indiceSpec
          in=nvalue(jsp)
          intermediaire=zdbt*rlapdi(in)
          do jlev=1,indiceLev
            v1(jsp,jlev)=v3(jsp,jlev)-intermediaire*v1(jsp,jlev)
          enddo
        enddo
        !$acc end kernels

        !$acc host_data use_device(a,v1,v2)
        call cublasDgemm('N','N',indiceSpec,indiceLev,indiceLev,alpha,v1,indiceSpec,a,indiceLev,beta,v2,indiceSpec)
        !$acc end host_data
        !$acc wait

        !$acc kernels
        do jsp=1,indiceSpec
          do jlev=1,indiceLev
             intermediaire=sivp(jlev)*2.31d0
             v4(jsp,jlev)=v2(jsp,jlev)/(intermediaire*rlapdi(nvalue(jsp)))
             !!v4(jsp,jlev)=v2(jsp,jlev)*(intermediaire*rlapdi(nvalue(jsp)))

          enddo
        enddo
        !$acc end kernels

        !$acc host_data use_device(a,v4,v5)
        call cublasDgemm('N','N',indiceSpec,indiceLev,indiceLev,alpha,v4,indiceSpec,b,indiceLev,beta,v5,IndiceSpec)
        !$acc end host_data
        !$acc wait
        
        !$acc kernels
        do jsp=1,indiceSpec
          do jlev=1,indiceLev
            v6(jsp,jlev)=v5(jsp,jlev)*4.107d0*4.107d0
          enddo
        enddo
        !$acc end kernels

        !$acc parallel private(nblignes,in,intermediaire) vector_length(32) default(none)
        !$acc loop gang vector private(inter1,inter2)
          do jsp=1,indiceSpec
            in=nvalue(jsp)
            intermediaire=zdbt*rlapdi(in)
            do jlev=1,indiceLev
              v1b(jsp,jlev)=v3b(jsp,jlev)-intermediaire*v1b(jsp,jlev)
            enddo
            inter1(1:indiceLev)=v1b(jsp,1:indiceLev)
            inter2 = matmul(at,inter1)
            v2b(jsp,1:indiceLev)=inter2(1:indiceLev)
            do jlev=1,indiceLev
              intermediaire=sivp(jlev)*2.31d0
              v4b(jsp,jlev)=v2b(jsp,jlev)/(intermediaire*rlapdi(nvalue(jsp)))
              !!v4b(jsp,jlev)=v2b(jsp,jlev)*(intermediaire*rlapdi(nvalue(jsp)))
            enddo
            inter1(1:indiceLev)=v4b(jsp,1:indiceLev)
            inter2=matmul(bt,inter1)
            v5b(jsp,1:indiceLev)=inter2(1:indiceLev)
            do jlev=1,indiceLev
              v6b(jsp,jlev)=v5b(jsp,jlev)*4.107d0*4.107d0
            enddo
        enddo
        !$acc end parallel

        enddo

        !$acc end data

        !!enddo

        print *, SUM((v6(:,:)-v6b(:,:))*(v6(:,:)-v6b(:,:)))
        print *,maxval(abs(v6(:,:)-v6b(:,:)))
        print *,"==== ligne 1"
        print *,"v6 : "
        print *,v6(1,:)
        print *,"v6b : "
        print *,v6b(1,:)
        print *,"==== différence lignes 1 à 8"
        do jsp=1,8
          print *,v6(jsp,:)-v6b(jsp,:) 
          print *," / "
        enddo
endprogram
