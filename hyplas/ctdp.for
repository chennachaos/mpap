CDOC BEGIN_SUBROUTINE CTDP
CDOC Consistent tangent matrix for the Drucker-Prager model
CDOC
CDOC The tangent matrix computed in this routine is consistent with
CDOC fully implicit elastic predictor/return mapping algorithm for the
CDOC Drucker-Prager elasto-plastic material model with piece-wise linear
CDOC isotropic hardening carried out in subroutine \smparm{SUDP}.
CDOC It contains only the plane strain and axisymmetric implementations.
CDOC It returns either the elastic tangent or the elasto-plastic
CDOC consistent tangent matrix depending on the input value of
CDOC the logical argument \smparm{EPFLAG}.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   >  Array of incremental plastic
CDOC C                          multipliers obtained in routine
CDOC C                          \smparm{SUDP}.
CDOC DOUBLE_PRECISION DMATX  <  Consistent tangent matrix.
CDOC LOGICAL          EPFLAG >  Elasto-plastic flag.
CDOC C                          If \smparm{.FALSE.}, \smparm{DMATX}
CDOC C                          returns as the elastic matrix. If
CDOC C                          \sparm{.TRUE.}, \smparm{DMATX} returns
CDOC C                          as the elasto-plastic tangent consistent
CDOC C                          with the return mapping algorithm
CDOC C                          implemented in routine \smparm{SUDP}.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDDP}.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC LOGICAL          LALGVA >  Array of logical algorithmic flags.
CDOC C                          See list of arguments of \smparm{SUDP}.
CDOC INTEGER          NTYPE  >  Stress state type.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          Same as in the argument list of
CDOC C                          subroutine \smparm{SUDP}.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Output of \smparm{SUDP}.
CDOC DOUBLE_PRECISION STRAN  >  Array of elastic trial (engineering)
CDOC C                          strain components last used as input of
CDOC C                          subroutine \smparm{SUDP}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto and P.H.Saksono, June 1996: Initial coding
CDOC
      SUBROUTINE CTDP
     1(   DGAM       ,DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,
     2    NTYPE      ,RPROPS     ,RSTAVA     ,STRAN      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=7  ,MSTRE=4)
      LOGICAL EPFLAG, LALGVA(3), TWOVEC
      DIMENSION
     1    DGAM(2)            ,DMATX(MSTRE,MSTRE),IPROPS(*)           ,
     2    RPROPS(*)          ,RSTAVA(MSTRE+1)   ,STRAN(MSTRE)
      DIMENSION
     1    EETD(MSTRE)        ,FOID(MSTRE,MSTRE)  ,SOID(MSTRE)        ,
     2    UNIDEV(MSTRE)
      DATA
     1    FOID(1,1),FOID(1,2),FOID(1,3),FOID(1,4)/
     2    1.0D0    ,0.0D0    ,0.0D0    ,0.0D0    /
     3    FOID(2,1),FOID(2,2),FOID(2,3),FOID(2,4)/
     4    0.0D0    ,1.0D0    ,0.0D0    ,0.0D0    /
     5    FOID(3,1),FOID(3,2),FOID(3,3),FOID(3,4)/
     6    0.0D0    ,0.0D0    ,0.5D0    ,0.0D0    /
     7    FOID(4,1),FOID(4,2),FOID(4,3),FOID(4,4)/
     8    0.0D0    ,0.0D0    ,0.0D0    ,1.0D0    /
      DATA
     1    SOID(1)  ,SOID(2)  ,SOID(3)  ,SOID(4)  /
     2    1.0D0    ,1.0D0    ,0.0D0    ,1.0D0    /
      DATA
     1    R0   ,R1   ,RP5  ,R2   ,R3   ,R9   /
     2    0.0D0,1.0D0,0.5D0,2.0D0,3.0D0,9.0D0/
C***********************************************************************
C COMPUTATION OF CONSISTENT TANGENT MODULUS FOR DRUCKER-PRAGER TYPE
C ELASTO-PLASTIC MATERIAL WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE AND
C PIECE-WISE LINEAR ISOTROPIC HARDENING
C***********************************************************************
      IF(NTYPE.EQ.2)THEN
        NSTRE=3
      ELSEIF(NTYPE.EQ.3)THEN
        NSTRE=4
      ELSE
        CALL ERRPRT('EI0017')
      ENDIF
C Accumulated plastic strain, DGAMA, DGAMB and two-vector algorithm flag
      EPBAR=RSTAVA(MSTRE+1)
      DGAMA=DGAM(1)
      DGAMB=DGAM(2)
      TWOVEC=LALGVA(3)
C Set some material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      ETA=RPROPS(4)
      XI=RPROPS(5)
      ETABAR=RPROPS(6)
      NHARD=IPROPS(3)
C and some constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R1D3=R1/R3
      R2D3=R2/R3
      R2D9=R2/R9
      ROOT2=SQRT(R2)
C
      IF(EPFLAG)THEN
C Compute elastoplastic consistent tangent
C ========================================
C Elastic trial deviatoric (physical) strain
        EEVD3=(STRAN(1)+STRAN(2)+STRAN(4))*R1D3
        EETD(1)=STRAN(1)-EEVD3
        EETD(2)=STRAN(2)-EEVD3
        EETD(3)=STRAN(3)*RP5
        EETD(4)=STRAN(4)-EEVD3
        ETDNOR=SQRT(EETD(1)*EETD(1)+EETD(2)*EETD(2)+
     1           R2*EETD(3)*EETD(3)+EETD(4)*EETD(4))
C Unit deviatoric flow vector
        IF(ETDNOR.NE.R0)THEN
          EDNINV=R1/ETDNOR
        ELSE
          EDNINV=R0
        ENDIF
        DO 10 I=1,NSTRE
          UNIDEV(I)=EETD(I)*EDNINV
   10   CONTINUE
C Hardening slope
        HSLOPE=DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
        IF(TWOVEC)THEN
C Elastoplastic tangent consistent with two-vector return (apex)
C --------------------------------------------------------------
          BDGAMA=ETABAR*DGAMA
          DGAMA2=DGAMA*DGAMA
          DEP=SQRT(R1D3*DGAMA2+R2D9*(BDGAMA+DGAMB)**2)
          AUX=BULK*ETA+XI*HSLOPE*R2D9/DEP*(BDGAMA+DGAMB)
          AFACT=BULK*(R1-BULK*ETA/AUX)
          BFACT=-BULK*(ROOT2*ETABAR-
     1          (BULK*ETA*ROOT2*ETABAR+XI*HSLOPE*R2D3/DEP*(ETDNOR+
     2          ROOT2/R3*ETABAR*(BDGAMA+DGAMB)))/AUX)
          DO 30 I=1,NSTRE
            DO 20 J=1,NSTRE
              DMATX(I,J)=SOID(I)*(AFACT*SOID(J)+BFACT*UNIDEV(J))
   20       CONTINUE       
   30     CONTINUE
        ELSE
C Elastoplastic tangent consistent with one-vector return (cone wall)
C -------------------------------------------------------------------
          AUX=R1/(GMODU+BULK*ETA*ETABAR+
     1            XI*HSLOPE*SQRT(R1D3+R2D9*ETABAR*ETABAR))
          AFACT=R2G*(R1-DGAMA/(ROOT2*ETDNOR))
          AFACD3=AFACT*R1D3
          BFACT=R2G*(DGAMA/(ROOT2*ETDNOR)-GMODU*AUX)
          CFACT=-ROOT2*GMODU*BULK*AUX
          DFACT=BULK*(R1-BULK*ETA*ETABAR*AUX)
          DO 50 I=1,NSTRE
            DO 40 J=1,NSTRE
              DMATX(I,J)=AFACT*FOID(I,J)+BFACT*UNIDEV(I)*UNIDEV(J)+
     1                   CFACT*(ETA*UNIDEV(I)*SOID(J)+
     2                           ETABAR*SOID(I)*UNIDEV(J))+
     3                   (DFACT-AFACD3)*SOID(I)*SOID(J)
   40       CONTINUE       
   50     CONTINUE
        ENDIF
      ELSE
C Compute elasticity matrix
C =========================
        FACTOR=BULK-R2G*R1D3
        DO 70 I=1,NSTRE
          DO 60 J=I,NSTRE
            DMATX(I,J)=R2G*FOID(I,J)+FACTOR*SOID(I)*SOID(J)
   60     CONTINUE       
   70   CONTINUE
        DO 90 J=1,NSTRE-1
          DO 80 I=J+1,NSTRE
            DMATX(I,J)=DMATX(J,I)
   80     CONTINUE
   90   CONTINUE
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE CTDP
