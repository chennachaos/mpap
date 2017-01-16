CDOC BEGIN_SUBROUTINE SUCADP
CDOC State update procedure for the Capped Drucker-Prager model.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the
CDOC Capped Drucker-Prager elasto-plastic material model with
CDOC piece-wise linear isotropic hardening. This routine contains only
CDOC the plane strain and axisymmetric implementations of the model.
CDOC For the present model,
CDOC hardening is considered by having the consolidation pressure
CDOC (location of the cap apex on the hydrostatic line) as a
CDOC user-defined piece-wise linear function
CDOC of the volumetric plastic strain.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   <  Array of incremental plastic
CDOC C                          multipliers.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element stored in this array used here.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDCADP}.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic yielding
CDOC C                          flag, \smparm{IFPLAS}; the return
CDOC C                          algorithm failure flag, \smparm{SUFAIL};
CDOC INTEGER          NTYPE  >  Stress state type.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, and the plastic
CDOC C                          properties: \smparm{ETA},
CDOC C                          \smparm{SQRJ20}, \smparm{ETABAR}
CDOC C                          and the pairs
CDOC C                          ``alpha-consolidation pressure''
CDOC C                          defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          \smparm{RDCADP}.
CDOC DOUBLE_PRECISION RSTAVA <> Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Previous converged values on entry,
CDOC C                          updated values on exit.
CDOC C                          The state variables stored in
CDOC C                          this array are the (engineering)
CDOC C                          elastic strain components and the
CDOC C                          accumulated plastic strain.
CDOC DOUBLE_PRECISION STRAN  >  Array of elastic trial (engineering)
CDOC C                          strain components.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, October 1997: Initial coding
CDOC
      SUBROUTINE SUCADP
     1(   DGAM       ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=8  ,MSTRE=4)
      LOGICAL  CAP, CAPAPX, CONAPX, CONE, IFPLAS, INTSCT, LALGVA(7),
     1         SUFAIL
      DIMENSION
     1    DGAM(2)            ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)    ,STRAN(*)           ,STRES(*)
      DIMENSION
     1    STRIAL(MSTRE)
      DATA
     1    R0   ,RP5  ,R1   ,R2   ,R3   ,TOL   / 
     2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,1.D-08/
      DATA MXITER/ 50 /
C***********************************************************************
C STRESS UPDATE PROCEDURE FOR THE CAPPED DRUCKER-PRAGER ELASTO-PLASTIC
C MATERIAL WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE AND PIECE_WISE
C LINEAR ISOTROPIC HARDENING:
C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM.
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C***********************************************************************
C Stop program if neither plane strain nor axisymmetric
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0015')
C Initialize some algorithmic and internal variables
      DGAMA=R0
      DGAMB=R0
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
      CONE  =.FALSE.
      CONAPX=.FALSE.
      CAP   =.FALSE.
      CAPAPX=.FALSE.
      INTSCT=.FALSE.
      ALPHAN=RSTAVA(MSTRE+1)
      ALPHA=ALPHAN
C Set some material properties
      YOUNG =RPROPS(2)
      POISS =RPROPS(3)
      ETA   =RPROPS(4)
      ETABAR=RPROPS(5)
      SQRJ20=RPROPS(6)
      ZETA  =RPROPS(7)
      NHARD =IPROPS(3)
C and some constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R1D3=R1/R3
C Compute elastic trial state
C ---------------------------
C Elastic trial volumetric strain and pressure stress
      EETV=STRAN(1)+STRAN(2)+STRAN(4)
      PT=BULK*EETV
C Elastic trial deviatoric stress
      EEVD3=EETV*R1D3
      STRIAL(1)=R2G*(STRAN(1)-EEVD3)
      STRIAL(2)=R2G*(STRAN(2)-EEVD3)
      STRIAL(4)=R2G*(STRAN(4)-EEVD3)
C shear component
      STRIAL(3)=R2G*(STRAN(3)*RP5)
C Compute elastic trial stress J2 invariant and cohesion
      VARJ2T=STRIAL(3)*STRIAL(3)+RP5*(STRIAL(1)*STRIAL(1)+
     1       STRIAL(2)*STRIAL(2)+STRIAL(4)*STRIAL(4))
C Check for plastic consistency
C -----------------------------
      SQRJ2T=SQRT(VARJ2T)
c compute value of standard D-P yield function
      AUXA=ETA*PT-SQRJ20
      PHIAT=SQRJ2T+AUXA
c compute value of CAP yield function
      PC=PLFUN(ALPHA,NHARD,RPROPS(IPHARD))
      AUXB=-ZETA*(PT+PC)
      PHIBT=SQRJ2T+AUXB
      PHIANO=PHIAT
      IF(AUXA.NE.R0)PHIANO=PHIANO/ABS(AUXA)
      PHIBNO=PHIBT
      IF(AUXB.NE.R0)PHIBNO=PHIBNO/ABS(AUXB)
      IF(PHIANO.GT.TOL.OR.PHIBNO.GT.TOL)THEN
C Plastic step: Use return mapping
C ================================
        IFPLAS=.TRUE.
C
        IF(PHIANO.GT.TOL)THEN
C Apply (closed-form) one-vector return to smooth D-P CONE wall
C -------------------------------------------------------------
C Compute incremental plastic multiplier
          DGAMA=PHIAT/(GMODU+BULK*ETABAR*ETA)
C Update hydrostatic pressure
          DEPV=ETABAR*DGAMA
          P=PT-BULK*DEPV
C Update hardening internal variable
          ALPHA=ALPHAN-DEPV
C Evaluate cap yield function at the returned state
          SQRJ2=SQRJ2T-GMODU*DGAMA
          PC=PLFUN(ALPHA,NHARD,RPROPS(IPHARD))
C Check validity of present return
          AUXB=-ZETA*(P+PC)
          PHIB=SQRJ2+AUXB
          PHIBNO=PHIB
          IF(AUXB.NE.R0)PHIBNO=PHIBNO/ABS(AUXB)
          IF(SQRJ2.GE.R0.AND.PHIBNO.LE.TOL)THEN
            CONE=.TRUE.
            FACTOR=R1-GMODU*DGAMA/SQRJ2T
          ELSEIF(SQRJ2.LT.R0)THEN
C Apply (closed-form) two-vector return to D-P CONE APEX
C ------------------------------------------------------
C Compute hydrostatic pressure at apex
            P=SQRJ20/ETA
C Update hardening internal variable
            ALPHA=ALPHAN+(P-PT)/BULK
            CONAPX=.TRUE.
            FACTOR=R0
          ENDIF
        ENDIF
C
        IF(PHIBNO.GT.TOL)THEN
C Apply one-vector return to smooth portion of the CAP 
C ----------------------------------------------------
          PHIB=PHIBT
          ALPHA=ALPHAN
          DGAMB=R0
C Begin Newton-Raphson iterations for DGAMB
          DO 10 NRITER=1,MXITER
            HSLOPE=DPLFUN(ALPHA,NHARD,RPROPS(IPHARD))
            DENOM=-GMODU-ZETA*ZETA*(BULK+HSLOPE)
c new guess for DGAMB
            DGAMB=DGAMB-PHIB/DENOM
c update some variables and check convergence
            SQRJ2=SQRJ2T-GMODU*DGAMB
            ALPHA=ALPHAN+ZETA*DGAMB
            PC=PLFUN(ALPHA,NHARD,RPROPS(IPHARD))
            P=PT+BULK*ZETA*DGAMB
            PHIB=SQRJ2-ZETA*(P+PC)
            AUXB=-ZETA*(P+PC)
            PHIBNO=ABS(PHIB)
            IF(AUXB.NE.R0)PHIBNO=PHIBNO/ABS(AUXB)
            IF(PHIBNO.LE.TOL)THEN
              GOTO 20
            ENDIF
   10     CONTINUE
C failure of stress update procedure
          SUFAIL=.TRUE.
          CALL ERRPRT('WE0005')
          GOTO 999
   20     CONTINUE
C Check validity
          AUXA=ETA*P-SQRJ20
          PHIANO=SQRJ2+AUXA
          IF(AUXA.NE.R0)PHIANO=PHIANO/ABS(AUXA)
          IF(PHIANO.LE.TOL.AND.SQRJ2.GE.R0)THEN
            CAP=.TRUE.
C What if SQRJ2T=0?????
            FACTOR=R1-GMODU*DGAMB/SQRJ2T
          ELSEIF(SQRJ2.LT.R0)THEN
C Apply two-vector return to CAP APEX
C -----------------------------------
            DGAMA=R0
            DGAMB=SQRJ2T/GMODU
            ALPHA=ALPHAN+ZETA*DGAMB
            PC=PLFUN(ALPHA,NHARD,RPROPS(IPHARD))
            P=PT+BULK*ZETA*DGAMB
            RES=P+PC
C Begin Newton-Raphson iterations for DGAMA
            DO 25 NRITER=1,MXITER
              HSLOPE=DPLFUN(ALPHA,NHARD,RPROPS(IPHARD))
              DENOM=BULK+HSLOPE
C compute new guess for DGAMA
              DGAMA=DGAMA-RES/DENOM 
C update hardening internal variable and P
              DEPV=-(ZETA*DGAMB+DGAMA)
              P=PT-BULK*DEPV
              ALPHA=ALPHAN-DEPV
              PC=PLFUN(ALPHA,NHARD,RPROPS(IPHARD))
c check convergence
              RES=P+PC
              RESNOR=ABS(RES)
              IF(PC.NE.R0)RESNOR=RESNOR/ABS(PC)
              IF(RESNOR.LE.TOL)THEN
                CAPAPX=.TRUE.
                FACTOR=R0
                GOTO 26
              ENDIF
   25       CONTINUE
C return to cap apex failed to converge
            SUFAIL=.TRUE.
            CALL ERRPRT('WE0006')
            GOTO 999
   26       CONTINUE
          ELSEIF(PHIANO.GT.TOL)THEN
C Apply two-vector return to INTERSECTION of D-P cone and cap
C -----------------------------------------------------------
            DGAMA=R0
            DGAMB=PHIAT/(GMODU-ETA*ZETA*BULK)
            ALPHA=ALPHAN+ZETA*DGAMB
            PC=PLFUN(ALPHA,NHARD,RPROPS(IPHARD))
            P=PT-BULK*(-ZETA*DGAMB)
            SQRJ2=SQRJ2T-GMODU*DGAMB
            PHIB=SQRJ2-ZETA*(P+PC)
            CONSTA=-GMODU+ZETA*ETABAR*BULK+
     1             (GMODU+ZETA*ZETA*BULK)*(GMODU+ETA*ETABAR*BULK)/
     2             (GMODU-ETA*ZETA*BULK)
            CONSTB=ZETA*(ZETA*(GMODU+ETA*ETABAR*BULK)/
     1             (GMODU-ETA*ZETA*BULK)-ETABAR)
c Begin Newton-Raphson iterations for DGAMA
            DO 27 NRITER=1,MXITER
              HSLOPE=DPLFUN(ALPHA,NHARD,RPROPS(IPHARD))
              DENOM=CONSTA+HSLOPE*CONSTB
c new guess for DGAMA
              DGAMA=DGAMA-PHIB/DENOM
c check convergence
              DGAMB=(PHIAT-DGAMA*(GMODU+ETA*ETABAR*BULK))/
     1              (GMODU-ETA*ZETA*BULK)
              DEPV=ETABAR*DGAMA-ZETA*DGAMB
              P=PT-BULK*DEPV
              ALPHA=ALPHAN-DEPV
              PC=PLFUN(ALPHA,NHARD,RPROPS(IPHARD))
              SQRJ2=SQRJ2T-GMODU*(DGAMA+DGAMB)
              PHIB=SQRJ2-ZETA*(P+PC)
              AUXB=-ZETA*(P+PC)
              PHIBNO=ABS(PHIB)
              IF(AUXB.NE.R0)PHIBNO=PHIBNO/ABS(AUXB)
              IF(PHIBNO.LE.TOL)THEN
                INTSCT=.TRUE.
C What if SQRJ2T=0?????
                FACTOR=R1-GMODU*(DGAMA+DGAMB)/SQRJ2T
                GOTO 30
              ENDIF
   27       CONTINUE
C return to cap apex failed to converge
            SUFAIL=.TRUE.
            CALL ERRPRT('WE0007')
            GOTO 999
   30       CONTINUE
          ENDIF 
        ENDIF
C Update converged stress components and other state variables
C ------------------------------------------------------------
C update hardening internal variable
        RSTAVA(MSTRE+1)=ALPHA
C update stress components
        STRES(1)=FACTOR*STRIAL(1)+P
        STRES(2)=FACTOR*STRIAL(2)+P
        STRES(3)=FACTOR*STRIAL(3)
        STRES(4)=FACTOR*STRIAL(4)+P
C compute converged elastic (engineering) strain components
        FACTOR=FACTOR/R2G
        EEVD3=P/(BULK*R3)
        RSTAVA(1)=FACTOR*STRIAL(1)+EEVD3
        RSTAVA(2)=FACTOR*STRIAL(2)+EEVD3
        RSTAVA(3)=FACTOR*STRIAL(3)*R2
        RSTAVA(4)=FACTOR*STRIAL(4)+EEVD3
      ELSE
C Elastic step: update stress using linear elastic law
C ====================================================
        STRES(1)=STRIAL(1)+PT
        STRES(2)=STRIAL(2)+PT
        STRES(3)=STRIAL(3)
        STRES(4)=STRIAL(4)+PT
C elastic engineering strain
        RSTAVA(1)=STRAN(1)
        RSTAVA(2)=STRAN(2)
        RSTAVA(3)=STRAN(3)
        RSTAVA(4)=STRAN(4)
      ENDIF
  999 CONTINUE
C Update some algorithmic variables before exit
C =============================================
      LALGVA(1)=IFPLAS
      LALGVA(2)=SUFAIL
      LALGVA(3)=CONE
      LALGVA(4)=CONAPX
      LALGVA(5)=CAP
      LALGVA(6)=CAPAPX
      LALGVA(7)=INTSCT
      DGAM(1)=DGAMA
      DGAM(2)=DGAMB
      RETURN
      END
CDOC END_SUBROUTINE SUCADP
