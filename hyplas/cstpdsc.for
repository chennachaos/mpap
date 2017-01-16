CDOC BEGIN_SUBROUTINE CSTPDS
CDOC Spatial tangent modulus for planar double slip single crystal model
CDOC
CDOC This routine computes the (elastic or elasto-plastic) spatial
CDOC tangent modulus, a, for the
CDOC planar double slip single crystal elasto-plastic model.
CDOC The tangent modulus computed here is fully consistent with the
CDOC exponential map-based elastic predictor/return mapping integration
CDOC algorithm implemented in subroutine \smparm{SUPDSC}.
CDOC This model is restricted to the plane strain case.
CDOC The elastic behaviour of the single crystal is assumed isotropic
CDOC (regularised) Neo-Hookean.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AMATX  <  Matrix of components of the spatial
CDOC C                          tangent (fourth order) modulus.
CDOC DOUBLE_PRECISION DGAM   >  Array of incremental plastic
CDOC C                          multipliers last determined in
CDOC C                          \sparm{SUPDSC}. All multipliers are set
CDOC C                          to zero for the first iteration of each
CDOC C                          load increment.
CDOC LOGICAL          EPFLAG >  Elasto-plastic flag.
CDOC C                          If \smparm{.FALSE.}, \smparm{AMATX}
CDOC C                          returns as the elastic modulus. If
CDOC C                          \sparm{.TRUE.}, \smparm{AMATX} returns
CDOC C                          as the elasto-plastic modulus.
CDOC DOUBLE_PRECISION FINCR  >  Last incremental deformation gradient.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC LOGICAL          LALGVA >  Array of logical algorithmic flags.
CDOC INTEGER          NTYPE  >  Stress state type flag. The present
CDOC C                          implementation is compatible only with
CDOC C                          plane strain (\smparm{NTYPE=2}).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA >  Array of current real state variables
CDOC C                          other than the stress tensor components.
CDOC DOUBLE_PRECISION RSTAVN >  Array of real state variables other than
CDOC C                          the stress tensor components at the
CDOC C                          beginning of the curent load increment.
CDOC DOUBLE_PRECISION STRES  >  Array of current Cauchy stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, November 1998: Initial coding (Miehe's approxim.)
CDOC
CDOC E.de Souza Neto, December 1998: Exact elasto-plastic tangent coded
CDOC
      SUBROUTINE CSTPDS
     1(   AMATX      ,DGAM       ,EPFLAG     ,FINCR      ,IPROPS     ,
     2    LALGVA     ,NTYPE      ,RPROPS     ,RSTAVA     ,RSTAVN     ,
     3    STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   IPHARD=6   ,IPHVAR=5   ,MADIM=5    ,NDIM=2     ,NGDIM=4    ,
     2    NIPROP=3   ,NLALGV=6   ,NRALGV=4   ,NRSTAV=5   ,NSTRE=4    ,
     3    NSYST=4    )
C Arguments
      LOGICAL
     1    EPFLAG     ,LALGVA(NLALGV)
      DIMENSION
     1    AMATX(MADIM,MADIM) ,DGAM(NRALGV)       ,FINCR(3,3)         ,
     2    IPROPS(NIPROP)     ,RPROPS(*)          ,RSTAVA(NRSTAV)     ,
     3    RSTAVN(NRSTAV)     ,STRES(NSTRE)
C Local arrays and variables
      LOGICAL
     1    S1ACT  ,S2ACT  ,S3ACT  ,S4ACT, NOCONV
      DIMENSION
     1    APMATX(NGDIM,NGDIM),AUXMTX(NGDIM,NGDIM),AUX2ND(NDIM,NDIM)  ,
     2    AUX4TH(NDIM,NDIM,NDIM,NDIM)            ,BEISO(3,3)         ,
     3    BMATX(NDIM,NDIM,NSYST)                 ,DELKRO(NDIM,NDIM)  ,
     4    DEREXP(3,3,3,3)    ,DEVPRJ(NGDIM,NGDIM),DEVSTR(NGDIM)      ,
     5    DS0M0(NDIM,NDIM)   ,DUMATX(NDIM,NDIM,NDIM,NDIM)            ,
     6    FE(2,2)            ,FEFETR(NGDIM,NGDIM),FEISO(2,2)         ,
     7    FEN(2,2)           ,FETISO(2,2)        ,FETRL(2,2)         ,
     8    FOIDS(NGDIM,NGDIM) ,FPILOG(3,3)        ,GINV(NSYST,NSYST)  ,
     9    GMATX(NSYST,NSYST) ,IACSET(NSYST)      ,
     O    SM0MS0(NDIM,NDIM,NSYST)                ,SOID(NGDIM)        ,
     1    S0M0(NDIM,NDIM,NSYST)                  ,
     2    UMATX(NDIM,NDIM,NDIM,NDIM)             ,VECM(NDIM,NSYST)   ,
     3    VECM0(NDIM,NSYST)  ,VECS(NDIM,NSYST)   ,VECS0(NDIM,NSYST)  ,
     4    VMATX(NGDIM,NGDIM)
C... Kroenecker delta
      DATA
     1    DELKRO(1,1),DELKRO(1,2)/
     2    1.0D0      ,0.D0       /
     3    DELKRO(2,1),DELKRO(2,2)/
     4    0.0D0      ,1.D0       /
C... fourth order (symmetric subspace) identity components stored in
C    matrix form using G matrix ordering (11,21,12,22)
      DATA
     1    FOIDS(1,1),FOIDS(1,2),FOIDS(1,3),FOIDS(1,4)/
     2    1.0D0     ,0.0D0     ,0.0D0     ,0.0D0     /
     3    FOIDS(2,1),FOIDS(2,2),FOIDS(2,3),FOIDS(2,4)/
     4    0.0D0     ,0.5D0     ,0.5D0     ,0.0D0     /
     5    FOIDS(3,1),FOIDS(3,2),FOIDS(3,3),FOIDS(3,4)/
     6    0.0D0     ,0.5D0     ,0.5D0     ,0.0D0     /
     7    FOIDS(4,1),FOIDS(4,2),FOIDS(4,3),FOIDS(4,4)/
     8    0.0D0     ,0.0D0     ,0.0D0     ,1.0D0     /
C... second order identity components in stored in vector form using G
C    matrix ordering
      DATA
     1    SOID(1)  ,SOID(2)  ,SOID(3)  ,SOID(4)  /
     2    1.0D0    ,0.0D0    ,0.0D0    ,1.D0     /
      DATA
     1    R0   ,R1   ,R2   ,R3   /
     2    0.0D0,1.0D0,2.0D0,3.0D0/
C***********************************************************************
C COMPUTATION OF THE CONSISTENT SPATIAL TANGENT MODULUS 'a' FOR
C THE PLANAR DOULBLE SLIP SINGLE CRYSTAL ELASTO-PLASTIC MODEL.
C MODEL VALID FOR PLANE STRAIN ONLY.
C***********************************************************************
C Stop program if not plane strain
      IF(NTYPE.NE.2)CALL ERRPRT('EI0035')
C Retrieve some state and algorithmic variables
C ---------------------------------------------
C... current hardening internal variable
      HRVAR=RSTAVA(IPHVAR)
C... current elastic deformation gradient
      FE(1,1)=RSTAVA(1)
      FE(2,1)=RSTAVA(2)
      FE(1,2)=RSTAVA(3)
      FE(2,2)=RSTAVA(4)
C... current active slip-systems logical flags
      S1ACT=LALGVA(3)
      S2ACT=LALGVA(4)
      S3ACT=LALGVA(5)
      S4ACT=LALGVA(6)
C... elastic deformation gradient at the beginning of the current load
C    increment
      FEN(1,1)=RSTAVN(1)
      FEN(2,1)=RSTAVN(2)
      FEN(1,2)=RSTAVN(3)
      FEN(2,2)=RSTAVN(4)
C Retrieve material properties
C ----------------------------
C... neo-Hookean constants
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C... initial system orientation
      THETA=RPROPS(4)
C... relative angle between systems
      BETA=RPROPS(5)
C... number of sampling points on hardening curve
      NHARD=IPROPS(3)
C Set some constants
C ------------------
      R1D3=R1/R3
      R2D3=R2*R1D3
C Assemble deviatoric projection tensor (use G matrix ordering)
      DO 20 I=1,NGDIM
        DO 10 J=1,NGDIM
          DEVPRJ(I,J)=FOIDS(I,J)-R1D3*SOID(I)*SOID(J)
   10   CONTINUE
   20 CONTINUE
C Get current Cauchy deviatoric stress and Cauchy hydrostatic pressure
      P=R1D3*(STRES(1)+STRES(2)+STRES(4))
C... use G matrix component ordering to store in-plane deviatoric
C    Cauchy stress components
      DEVSTR(1)=STRES(1)-P
      DEVSTR(2)=STRES(3)
      DEVSTR(3)=STRES(3)
      DEVSTR(4)=STRES(2)-P
C Compute isochoric component of Fe and Be
      DETFE=FE(1,1)*FE(2,2)-FE(1,2)*FE(2,1)
      FACTOR=DETFE**(-R1D3)
      FEISO(1,1)=FACTOR*FE(1,1)
      FEISO(1,2)=FACTOR*FE(1,2)
      FEISO(2,1)=FACTOR*FE(2,1)
      FEISO(2,2)=FACTOR*FE(2,2)
      CALL RVZERO(BEISO,9)
      DO 50 I=1,2
        DO 40 J=1,2
          DO 30 K=1,2
            BEISO(I,J)=BEISO(I,J)+FEISO(I,K)*FEISO(J,K)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
      BEISO(3,3)=FACTOR*FACTOR
C Trace of isochoric component of Be
      TRBISO=BEISO(1,1)+BEISO(2,2)+BEISO(3,3)
C
C
C Compute ELASTIC tangent modulus
C ===============================
C
      GFAC=R2D3*GMODU*TRBISO/DETFE
      BULFAC=BULK/DETFE
      R2P=R2*P
C... assemble tensorially compact part
      DO 70 I=1,NGDIM
        DO 60 J=1,NGDIM
          AMATX(I,J)=BULFAC*SOID(I)*SOID(J)-R2P*FOIDS(I,J)+
     1               GFAC*DEVPRJ(I,J)-
     2               R2D3*(DEVSTR(I)*SOID(J)+SOID(I)*DEVSTR(J))
   60   CONTINUE
   70 CONTINUE
C... add non-compact part: delta_ik sigma_jl
      AMATX(1,1)=AMATX(1,1)+STRES(1)
      AMATX(3,1)=AMATX(3,1)+STRES(3)
      AMATX(2,2)=AMATX(2,2)+STRES(1)
      AMATX(4,2)=AMATX(4,2)+STRES(3)
      AMATX(1,3)=AMATX(1,3)+STRES(3)
      AMATX(3,3)=AMATX(3,3)+STRES(2)
      AMATX(2,4)=AMATX(2,4)+STRES(3)
      AMATX(4,4)=AMATX(4,4)+STRES(2)
C
C
      IF(EPFLAG)THEN
C
C Compute and add algorithm-consistent PLASTIC contribution to spatial
C tangent modulus
C ====================================================================
C
C Compute individual terms needed to assemble the plastic contribution
C --------------------------------------------------------------------
C
C Last elastic trial deformation gradient
        CALL RVZERO(FETRL,4)
        DO 100 I=1,2
          DO 90 J=1,2
            DO 80 K=1,2
              FETRL(I,J)=FETRL(I,J)+FINCR(I,K)*FEN(K,J)
   80       CONTINUE
   90     CONTINUE
  100   CONTINUE
C... isochoric component
        DETFET=FETRL(1,1)*FETRL(2,2)-FETRL(1,2)*FETRL(2,1)
        VOLFAC=DETFET**(-R1D3)
        FETISO(1,1)=VOLFAC*FETRL(1,1)
        FETISO(2,1)=VOLFAC*FETRL(2,1)
        FETISO(1,2)=VOLFAC*FETRL(1,2)
        FETISO(2,2)=VOLFAC*FETRL(2,2)
C Assemble relevant fourth order tensor
        DO 140 I=1,NDIM
          DO 130 J=1,NDIM
            DO 120 K=1,NDIM
              DO 110 L=1,NDIM
                AUX4TH(I,J,K,L)=FEISO(I,L)*FETISO(J,K)
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
C... rearrange in matrix form with G matrix component ordering
        CALL ARRGO2(AUX4TH,FEFETR)
C
C Compute exact exponential map derivative
C ----------------------------------------
C
C retrieve information on current set of active slip-systems and set up
C corresponding initial slip-system vectors
        IACSYS=0
        IF(S1ACT)THEN
C... system 1:
          IACSYS=IACSYS+1
          IACSET(IACSYS)=1
          VECS0(1,1)=COS(THETA)
          VECS0(2,1)=SIN(THETA)
          VECM0(1,1)=-SIN(THETA)
          VECM0(2,1)=COS(THETA)
        ENDIF
        IF(S2ACT)THEN
C... system 2:
          IACSYS=IACSYS+1
          IACSET(IACSYS)=2
          VECS0(1,2)=COS(THETA+BETA)
          VECS0(2,2)=SIN(THETA+BETA)
          VECM0(1,2)=-SIN(THETA+BETA)
          VECM0(2,2)=COS(THETA+BETA)
        ENDIF
        IF(S3ACT)THEN
          IACSYS=IACSYS+1
          IACSET(IACSYS)=3
C... system 3:
          VECS0(1,3)=-COS(THETA)
          VECS0(2,3)=-SIN(THETA)
          VECM0(1,3)=-SIN(THETA)
          VECM0(2,3)=COS(THETA)
        ENDIF
        IF(S4ACT)THEN
C... system 4:
          IACSYS=IACSYS+1
          IACSET(IACSYS)=4
          VECS0(1,4)=-COS(THETA+BETA)
          VECS0(2,4)=-SIN(THETA+BETA)
          VECM0(1,4)=-SIN(THETA+BETA)
          VECM0(2,4)=COS(THETA+BETA)
        ENDIF
C... number of currently active systems
        NACSYS=IACSYS
C Compute current elastic push forward of the slip-system vectors
C of the active systems
        CALL RVZERO(VECS,NDIM*NSYST)
        CALL RVZERO(VECM,NDIM*NSYST)
        DO 170 I=1,NDIM
          DO 160 J=1,NDIM
            DO 150 II=1,NACSYS
              ISYST=IACSET(II)
              VECS(I,ISYST)=VECS(I,ISYST)+FEISO(I,J)*VECS0(J,ISYST)
              VECM(I,ISYST)=VECM(I,ISYST)+FEISO(I,J)*VECM0(J,ISYST)
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
C Current slope of Taylor hardening curve
        HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS(IPHARD))
C Compute logarithm of inverse of incremental plastic deformation
C gradient by summing up contributions from each active slip system
        CALL RVZERO(FPILOG,9)
        DO 200 II=1,NACSYS
          ISYST=IACSET(II)
          DO 190 I=1,2
            DO 180 J=1,2
              FPILOG(I,J)=FPILOG(I,J)-
     1                    DGAM(ISYST)*VECS0(I,ISYST)*VECM0(J,ISYST)
  180       CONTINUE
  190     CONTINUE
  200   CONTINUE
C... and the corresponding exact derivative of the exponential map
        CALL DEXPMP(   DEREXP      ,NOCONV     ,FPILOG     )
C
C Compute jacobian of non-linear system of return mapping equations
C -----------------------------------------------------------------
C compute some preliminary matrices
        CALL RVZERO(BMATX,NDIM*NDIM*NSYST)
        DO 300 II=1,NACSYS
          ISYST=IACSET(II)
          DO 220 I=1,2
            DO 210 J=1,2
              S0M0(I,J,ISYST)=VECS0(I,ISYST)*VECM0(J,ISYST)
              SM0MS0(I,J,ISYST)=VECS(I,ISYST)*VECM0(J,ISYST)+
     1                          VECM(I,ISYST)*VECS0(J,ISYST)
  210       CONTINUE
  220     CONTINUE
          CALL RVZERO(DS0M0,NDIM*NDIM)
          DO 260 I=1,2
            DO 250 J=1,2
              DO 240 K=1,2
                DO 230 L=1,2
                  DS0M0(I,J)=DS0M0(I,J)+
     1                       DEREXP(I,J,K,L)*S0M0(K,L,ISYST)
  230           CONTINUE
  240         CONTINUE
  250       CONTINUE
  260     CONTINUE
          DO 290 I=1,2
            DO 280 J=1,2
              DO 270 K=1,2
                BMATX(I,J,ISYST)=BMATX(I,J,ISYST)+
     1                           FETISO(I,K)*DS0M0(K,J)
  270         CONTINUE
  280       CONTINUE
  290     CONTINUE
  300   CONTINUE
C Assemble exact jacobian of non-linear system
        DO 320 II=1,NACSYS
          ISYST=IACSET(II)
          DO 310 JJ=1,NACSYS
            JSYST=IACSET(JJ)
            GMATX(II,JJ)=GMODU*
     1         SCAPRD(SM0MS0(1,1,ISYST),BMATX(1,1,JSYST),NDIM*NDIM)+
     2         HSLOPE
  310     CONTINUE
  320   CONTINUE
C Invert jacobian: Note that for the double slip model only one or two
C systems may be active
        IF(NACSYS.EQ.1)THEN
          IF(GMATX(1,1).EQ.R0)CALL ERRPRT('EE0006')
          GINV(1,1)=R1/GMATX(1,1)
        ELSEIF(NACSYS.EQ.2)THEN
          DETG=GMATX(1,1)*GMATX(2,2)-GMATX(1,2)*GMATX(2,1)
          IF(DETG.EQ.R0)CALL ERRPRT('EE0006')
          DETGIN=R1/DETG
          GINV(1,1)=GMATX(2,2)*DETGIN
          GINV(2,2)=GMATX(1,1)*DETGIN
          GINV(1,2)=-GMATX(1,2)*DETGIN
          GINV(2,1)=-GMATX(2,1)*DETGIN
        ENDIF
C Compute U matrix
C ----------------
        CALL RVZERO(UMATX,NDIM*NDIM*NDIM*NDIM)
        DO 380 I=1,NDIM
          DO 370 J=1,NDIM
            DO 360 K=1,NDIM
              DO 350 L=1,NDIM
                DO 340 II=1,NACSYS
                  ISYST=IACSET(II)
                  DO 330 JJ=1,NACSYS
                    JSYST=IACSET(JJ)
                    UMATX(I,J,K,L)=UMATX(I,J,K,L)+S0M0(I,J,ISYST)*
     1                             GINV(II,JJ)*SM0MS0(K,L,JSYST)
  330             CONTINUE
  340           CONTINUE
  350         CONTINUE
  360       CONTINUE
  370     CONTINUE
  380   CONTINUE
C Compute product [D:U]
        CALL RVZERO(DUMATX,NDIM*NDIM*NDIM*NDIM)
        DO 440 I=1,NDIM
          DO 430 J=1,NDIM
            DO 420 K=1,NDIM
              DO 410 L=1,NDIM
                DO 400 M=1,NDIM
                  DO 390 N=1,NDIM
                    DUMATX(I,J,K,L)=DUMATX(I,J,K,L)+
     1                              DEREXP(I,J,M,N)*UMATX(M,N,K,L)
  390             CONTINUE
  400           CONTINUE
  410         CONTINUE
  420       CONTINUE
  430     CONTINUE
  440   CONTINUE
C... and the contribution to a^p involving the product D:U
        CALL RVZERO(AUX4TH,NDIM*NDIM*NDIM*NDIM)
        DO 490 I=1,NDIM
          DO 480 J=1,NDIM
            DO 470 K=1,NDIM
              DO 460 L=1,NDIM
                DO 450 M=1,NDIM
                  AUX4TH(I,J,K,L)=AUX4TH(I,J,K,L)+
     1                            DUMATX(I,J,K,M)*FEISO(L,M)
  450           CONTINUE
  460         CONTINUE
  470       CONTINUE
  480     CONTINUE
  490   CONTINUE
        CALL RVZERO(AUX2ND,NDIM*NDIM)
        DO 530 I=1,NDIM
          DO 520 J=1,NDIM
            DO 510 K=1,NDIM
              DO 500 L=1,NDIM
                AUX2ND(I,J)=AUX2ND(I,J)+DUMATX(I,J,K,L)*FEISO(K,L)
  500         CONTINUE
  510       CONTINUE
  520     CONTINUE
  530   CONTINUE
        DO 570 I=1,NDIM
          DO 560 J=1,NDIM
            DO 550 K=1,NDIM
              DO 540 L=1,NDIM
                AUX4TH(I,J,K,L)=AUX4TH(I,J,K,L)-
     1                          R1D3*AUX2ND(I,J)*DELKRO(K,L)
  540         CONTINUE
  550       CONTINUE
  560     CONTINUE
  570   CONTINUE
C... rearrange in matrix form
        CALL ARRGO2(AUX4TH,VMATX)
        CALL RVZERO(AUXMTX,NGDIM*NGDIM)
        DO 600 I=1,NGDIM
          DO 590 J=1,NGDIM
            DO 580 K=1,NGDIM
              AUXMTX(I,J)=AUXMTX(I,J)+FEFETR(I,K)*VMATX(K,J)
  580       CONTINUE
  590     CONTINUE
  600   CONTINUE
C Compute plastic contribution
        CALL RVZERO(APMATX,NGDIM*NGDIM)
        AUX=R2*GMODU*GMODU/DETFE
        DO 630 I=1,NGDIM
          DO 620 J=1,NGDIM
            DO 610 K=1,NGDIM
              APMATX(I,J)=APMATX(I,J)-AUX*DEVPRJ(I,K)*AUXMTX(K,J)
  610       CONTINUE
  620     CONTINUE
  630   CONTINUE
C Add plastic contribution to spatial tangent modulus
        DO 650 I=1,NGDIM
          DO 640 J=1,NGDIM
            AMATX(I,J)=AMATX(I,J)+APMATX(I,J)
  640     CONTINUE
  650   CONTINUE
C
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE CSTPDS
