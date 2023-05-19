program testmm
        use cublas
        use cutensorex
        implicit none
#define dynamique 1
        integer(kind=4), parameter  :: indiceSpec = 250500 
        integer(kind=4), parameter  :: indiceLev  = 105     
#if dynamique
        real(kind=4), allocatable  :: a(:,:),at(:,:),b(:,:),bt(:,:)
        real(kind=4), allocatable  :: v1(:,:),v2(:,:)
        real(kind=4), allocatable  :: v3(:,:),v4(:,:),v5(:,:),v6(:,:)
#else
        real(kind=4), dimension(indiceLev,indiceLev) ::a,at,b,bt
        real(kind=4), dimension(indiceSpec,indiceLev)::v1,v2
        real(kind=4), dimension(indiceSpec,indiceLev)::v3,v4,v5,v6
#endif
        real(kind=4), dimension (indiceSpec)     :: rlapdi,tech1
        integer(kind=4) , dimension (indiceSpec)     :: nvalue
        real(kind=4), dimension (indiceLev)         :: sivp
        real(kind=4), dimension (indiceLev) :: inter1,inter2
        real(kind=4), parameter :: alpha=1.0
        real(kind=4), parameter :: beta=0.0
        real(kind=4), parameter :: zdbt=0.69123d0
        real(kind=4), parameter :: const2=0.231d0
        real(kind=4), parameter :: const3=0.4107d0
        integer(kind=4)             :: indice
        real(kind=4)            :: intermediaire
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
        at=transpose(a/50.0d0)
        bt=transpose(b/50.0d0)

        !$acc data copy(at,bt,v1,v3,nvalue,sivp,rlapdi) copyout(v6) create(inter1,inter2,v2,v4,v5)
  
        do compteur=1,20
        !$acc parallel private(nblignes,indice,intermediaire) vector_length(32) default(none)
        !$acc loop gang vector private(inter1,inter2)
          do jsp=1,indiceSpec
            indice=nvalue(jsp)
            intermediaire=zdbt*rlapdi(indice)
            do jlev=1,indiceLev
              v1(jsp,jlev)=v3(jsp,jlev)-intermediaire*v1(jsp,jlev)
            enddo
            inter1(1:indiceLev)=v1(jsp,1:indiceLev)
            inter2 = matmul(at,inter1)
            v2(jsp,1:indiceLev)=inter2(1:indiceLev)
            do jlev=1,indiceLev
              intermediaire=sivp(jlev)*const2
              v4(jsp,jlev)=v2(jsp,jlev)/(intermediaire*rlapdi(nvalue(jsp)))
              !!v4b(jsp,jlev)=v2b(jsp,jlev)*(intermediaire*rlapdi(nvalue(jsp)))
            enddo
            inter1(1:indiceLev)=v4(jsp,1:indiceLev)
            inter2=matmul(bt,inter1)
            v5(jsp,1:indiceLev)=inter2(1:indiceLev)
            do jlev=1,indiceLev
              v6(jsp,jlev)=v5(jsp,jlev)*const3*const3
            enddo
        enddo
        !$acc end parallel

        enddo

        !$acc end data
    
        print *,v6(1,:)   
  
#if dynamique
deallocate(a,b,at,bt,v1,v2,v3,v4,v5,v6)
#endif

endprogram
