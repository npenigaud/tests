program testmm
        use cublas
        use cutensorex
        implicit none
#define dynamique 1
        integer(kind=4), parameter  :: indiceSpec = 250500 
        integer(kind=4), parameter  :: indiceLev  = 105     
#if dynamique
        double precision, allocatable  :: a(:,:),at(:,:),b(:,:),bt(:,:)
        double precision, allocatable  :: v1(:,:),v2(:,:)
        double precision, allocatable  :: v3(:,:),v4(:,:),v5(:,:),v6(:,:)
#else
        double precision, dimension(indiceLev,indiceLev) ::a,at,b,bt
        double precision, dimension(indiceSpec,indiceLev)::v1,v2
        double precision, dimension(indiceSpec,indiceLev)::v3,v4,v5,v6
#endif
        double precision, dimension (indiceSpec)     :: rlapdi,tech1
        integer(kind=4) , dimension (indiceSpec)     :: nvalue
        double precision, dimension (indiceLev)         :: sivp
        double precision, dimension (indiceLev) :: inter1,inter2
        double precision, parameter :: alpha=1.0
        double precision, parameter :: beta=0.0
        double precision, parameter :: zdbt=6.9123d0
        double precision, parameter :: const2=2.31d0
        double precision, parameter :: const3=4.107d0
        integer(kind=4)             :: indice
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

        !$acc data copy(at,bt,v1,v3,nvalue,sivp,rlapdi) copyout(v6) create(inter1,inter2,v2,v4,v5)
  
        do compteur=1,20
        !$acc parallel private(nblignes,indice,intermediaire) vector_length(32) default(none)
        !$acc loop gang vector private(inter1,inter2)
          do jsp=1,indiceSpec
            indice=nvalue(jsp)
            intermediaire=zdbt*rlapdi(indice)
            !$acc loop seq
            do jlev=1,indiceLev
              v1(jsp,jlev)=v3(jsp,jlev)-intermediaire*v1(jsp,jlev)
            enddo

            !$acc loop seq
            do jlev=1,indiceLev
              inter1(jlev)=v1(jsp,jlev)
            enddo
            inter2 = matmul(at,inter1)
            !$acc loop seq
            do jlev=1,indiceLev
              v2(jsp,jlev)=inter2(jlev)
            enddo 

            !$acc loop seq
            do jlev=1,indiceLev
              intermediaire=sivp(jlev)*const2
              v4(jsp,jlev)=v2(jsp,jlev)/(intermediaire*rlapdi(nvalue(jsp)))
              !!v4b(jsp,jlev)=v2b(jsp,jlev)*(intermediaire*rlapdi(nvalue(jsp)))
            enddo

            !$acc loop seq
            do jlev=1,indiceLev
              inter1(jlev)=v4(jsp,jlev)
            enddo
            inter2=matmul(bt,inter1)
            !$acc loop seq
            do jlev=1,indiceLev
              v5(jsp,jlev)=inter2(jlev)
            enddo

            !$acc loop seq
            do jlev=1,indiceLev
              v6(jsp,jlev)=v5(jsp,jlev)*const3*const3
            enddo
        enddo
        !$acc end parallel

        enddo

        !$acc end data
        
        print *,v6(1,:)
#if dynamique
        deallocate(a,at,b,bt,v1,v2,v3,v4,v5,v6)
#endif
endprogram
