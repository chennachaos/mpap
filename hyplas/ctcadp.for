CDOC BEGIN_SUBROUTINE CTCADP
CDOC Consistent tangent matrix for the Capped Drucker-Prager model
CDOC
CDOC The tangent matrix computed in this routine is consistent with
CDOC fully implicit elastic predictor/return mapping algorithm for the
CDOC Capped Drucker-Prager elasto-plastic model with piece-wise linear
CDOC isotropic hardening carried out in subroutine \smparm{SUCADP}.
CDOC This routine contains only the plane strain and the axisymmetric
CDOC implementations of the model.
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
CDOC C                          implemented in routine \smparm{SUCADP}.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDCADP}.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC LOGICAL          LALGVA >  Array of logical algorithmic flags.
CDOC C                          See list of arguments of
CDOC C                          \smparm{SUCADP}.
CDOC INTEGER          NTYPE  >  Stress state type.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          Same as in the argument list of
CDOC C                          subroutine \smparm{SUCADP}.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Output of \smparm{SUCADP}.
CDOC DOUBLE_PRECISION STRAN  >  Array of elastic trial (engineering)
CDOC C                          strain components. Same as in the input
CDOC C                          \smparm{SUCADP}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto,  October 1997: Initial coding
CDOC
      SUBROUTINE CTCADP
     1(   DGAM       ,DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,
     2    NTYPE      ,RPROPS     ,RSTAVA     ,STRAN      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=8  ,MSTRE=4)
      LOGICAL CAP, CAPAPX, CONAPX, CONE, EPFLAG, INTSCT, LALGVA(7)
      DIMENSION
     1    DGAM(2)            ,DMATX(MSTRE,MSTRE),IPROPS(*)           ,
     2    RPROPS(*)          ,RSTAVA(MSTRE+1)   ,STRAN(*)
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
     1    R0   ,R1   ,RP5  ,R2   ,R3   /
     2    0.0D0,1.0D0,0.5D0,2.0D0,3.0D0/
C***********************************************************************
C COMPUTATION OF CONSISTENT TANGENT MODULUS FOR CAPPED DRUCKER-PRAGER
C ELASTO-PLASTIC MODEL WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE AND
C PIECE-WISE LINEAR ISOTROPIC HARDENING.
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C***********************************************************************
      IF(NTYPE.EQ.2)THEN
        NSTRE=3
      ELSEIF(NTYPE.EQ.3)THEN
        NSTRE=4
      ELSE
        CALL ERRPRT('EI0014')
      ENDIF
C Accumulated plastic strain, DGAMA, DGAMB and logical algorithm flags
      EPBAR=RSTAVA(MSTRE+1)
      DGAMA=DGAM(1)
      DGAMB=DGAM(2)
      CONE  =LALGVA(3)
      CONAPX=LALGVA(4)
      CAP   =LALGVA(5)
      CAPAPX=LALGVA(6)
      INTSCT=LALGVA(7)
C Set some material properties
      YOUNG =RPROPS(2)
      POISS =RPROPS(3)
      ETA   =RPROPS(4)
      ETABAR=RPROPS(5)
      ZETA  =RPROPS(7)
      NHARD =IPROPS(3)
C and some constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R1D3=R1/R3
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
        IF(CONAPX)THEN
C Elastoplastic tangent consistent with return to D-P CONE APEX
C -------------------------------------------------------------
          DO 30 I=1,NSTRE
            DO 20 J=1,NSTRE
              DMATX(I,J)=R0
   20       CONTINUE       
   30     CONTINUE
        ELSEIF(CONE)THEN
C Elastoplastic tangent consistent with one-vector return to D-P CONE
C -------------------------------------------------------------------
          AUX=R1/(GMODU+BULK*ETA*ETABAR)
          AFACT=R2G*(R1-DGAMA/(ROOT2*ETDNOR))
          AFACD3=AFACT*R1D3
          BFACT=R2G*(DGAMA/(ROOT2*ETDNOR)-GMODU*AUX)
          CFACT=-ROOT2*GMODU*BULK*AUX
          DFACT=BULK*(R1-BULK*ETA*ETABAR*AUX)
          DO 50 I=1,NSTRE
            DO 40 J=1,NSTRE
              DMATX(I,J)=AFACT*FOID(I,J)+BFACT*UNIDEV(I)*UNIDEV(J)+
     1                   CFACT*(ETA*UNIDEV(I)*SOID(J)+
     2                          ETABAR*SOID(I)*UNIDEV(J))+
     3                   (DFACT-AFACD3)*SOID(I)*SOID(J)
   40       CONTINUE       
   50     CONTINUE
        ELSEIF(CAP)THEN
C Elastoplastic tangent consistent with one-vector return to CAP
C --------------------------------------------------------------
          AUX=R1/(GMODU+ZETA*ZETA*(BULK+HSLOPE))
          AFACT=R2G*(R1-DGAMB/(ROOT2*ETDNOR))
          AFACD3=AFACT*R1D3
          BFACT=R2G*(DGAMB/(ROOT2*ETDNOR)-GMODU*AUX)
          CFACT=ROOT2*GMODU*BULK*AUX*ZETA
          DFACT=BULK*(R1-BULK*ZETA*ZETA*AUX)
          DO 70 I=1,NSTRE
            DO 60 J=1,NSTRE
              DMATX(I,J)=AFACT*FOID(I,J)+BFACT*UNIDEV(I)*UNIDEV(J)+
     1                   CFACT*(UNIDEV(I)*SOID(J)+SOID(I)*UNIDEV(J))+
     2                   (DFACT-AFACD3)*SOID(I)*SOID(J)
   60       CONTINUE       
   70     CONTINUE
        ELSEIF(CAPAPX)THEN
C Elastoplastic tangent consistent with two-vector return to CAP APEX
C -------------------------------------------------------------------
          AFACT=BULK*HSLOPE/(BULK+HSLOPE)
          DO 90 I=1,NSTRE
            DO 80 J=1,NSTRE
              DMATX(I,J)=AFACT*SOID(I)*SOID(J)
   80       CONTINUE       
   90     CONTINUE
        ELSEIF(INTSCT)THEN
C Elastoplastic tangent consistent with return to INTERSECTION CONE/CAP
C ---------------------------------------------------------------------
          ZETA2=ZETA*ZETA
          AUX1=GMODU+BULK*ETABAR*ETA
          AUX2=GMODU-BULK*ETA*ZETA
          A=ROOT2*GMODU*(R1-(GMODU+ZETA2*(BULK+HSLOPE))/AUX2)
          B=-BULK*(ZETA+ETA*(GMODU+ZETA2*(BULK+HSLOPE))/AUX2)
          C=GMODU*(R1-AUX1/AUX2)-ZETA*(BULK+HSLOPE)*
     1                                (ETABAR+ZETA*AUX1/AUX2)
          AFACT=R2G*(R1-(DGAMA+DGAMB)/(ROOT2*ETDNOR))
          AFACD3=AFACT*R1D3
          BFACT=ROOT2*GMODU*((DGAMA+DGAMB)/ETDNOR-ROOT2*GMODU/AUX2+
     1          A/C*(AUX1/AUX2-R1))
          CFACT=ROOT2*GMODU*(B/C*(AUX1/AUX2-R1)-BULK*ETA/AUX2)
          DFACT=BULK*(ROOT2*GMODU*ZETA/AUX2-A/C*(ETABAR+ZETA*AUX1/AUX2))
          EFACT=BULK*(R1+BULK*ZETA*ETA/AUX2-B/C*(ETABAR+ZETA*AUX1/AUX2))
          DO 110 I=1,NSTRE
            DO 100 J=1,NSTRE
              DMATX(I,J)=AFACT*FOID(I,J)+BFACT*UNIDEV(I)*UNIDEV(J)+
     1                   CFACT*UNIDEV(I)*SOID(J)+
     2                   DFACT*SOID(I)*UNIDEV(J)+
     3                   (EFACT-AFACD3)*SOID(I)*SOID(J)
  100       CONTINUE       
  110     CONTINUE
        ENDIF
      ELSE
C Compute elasticity matrix
C =========================
        FACTOR=BULK-R2G*R1D3
        DO 130 I=1,NSTRE
          DO 120 J=I,NSTRE
            DMATX(I,J)=R2G*FOID(I,J)+FACTOR*SOID(I)*SOID(J)
  120     CONTINUE       
  130   CONTINUE
        DO 150 J=1,NSTRE-1
          DO 140 I=J+1,NSTRE
            DMATX(I,J)=DMATX(J,I)
  140     CONTINUE
  150   CONTINUE
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE CTCADP


