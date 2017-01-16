CDOC BEGIN_SUBROUTINE SUDP
CDOC State update procedure for the Drucker-Prager material model.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the
CDOC Drucker elasto-plastic material model with piece-wise linear
CDOC isotropic hardening. This routine contains only the plane strain
CDOC and axisymmetric implementations of the model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   <  Array of incremental plastic
CDOC C                          multipliers.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element stored in this array used here.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDDP}.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic yielding
CDOC C                          flag, \smparm{IFPLAS}; the return
CDOC C                          algorithm failure flag, \smparm{SUFAIL};
CDOC C                          the two-vector or apex return flag,
CDOC C                          \smparm{TWOVEC}.
CDOC C                          The plastic yielding flag is set
CDOC C                          to \smparm{.TRUE.} if plastic yielding
CDOC C                          has occurred and to \smparm{.FALSE.} if
CDOC C                          the step is elastic. The algorithm
CDOC C                          failure flag is set to {\tt .FALSE.} if
CDOC C                          the state update algorithm has been
CDOC C                          successful and to \smparm{.TRUE.} if the
CDOC C                          return mapping algorithm has failed
CDOC C                          to converge.
CDOC C                          \smparm{TWOVEC} is set to
CDOC C                          \smparm{.TRUE.}
CDOC C                          if the selected return mapping is the
CDOC C                          two-vector return to the apex. It is
CDOC C                          set to \smparm{.FALSE.} otherwise.
CDOC INTEGER          NTYPE  >  Stress state type.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, and the plastic
CDOC C                          properties: \smparm{ETA},
CDOC C                          \smparm{XI}, \smparm{ETABAR}
CDOC C                          and the pairs
CDOC C                          ``accumulated plastic strain-cohesion''
CDOC C                          defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          \smparm{RDDP}.
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
CDOC E.de Souza Neto and P.H.Saksono, June 1996: Initial coding
CDOC
      SUBROUTINE SUDP
     1(   DGAM       ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=7  ,MSTRE=4)
      LOGICAL IFPLAS, LALGVA(3), SUFAIL, TWOVEC
      DIMENSION
     1    DGAM(2)            ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)    ,STRAN(MSTRE)       ,STRES(MSTRE)
      DIMENSION
     1    STRIAL(MSTRE)
      DATA
     1    R0   ,RP5  ,R1   ,R2   ,R3   ,R9   ,TOL   / 
     2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,9.0D0,1.D-08/
      DATA MAXRT / 50 /
C***********************************************************************
C STRESS UPDATE PROCEDURE FOR DRUCKER PRAGER TYPE ELASTO-PLASTIC
C MATERIAL WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE AND PIECE_WISE
C LINEAR ISOTROPIC HARDENING:
C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM (Boxes 8.8-10)
C***********************************************************************
C Stops program if neither plane strain nor axisymmetric
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0016')
C Initialize some algorithmic and internal variables
      DGAMA=R0
      DGAMB=R0
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
      EPBARN=RSTAVA(MSTRE+1)
      EPBAR=EPBARN
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
      R2D9=R2/R9
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
      COHE=PLFUN(EPBARN,NHARD,RPROPS(IPHARD))
C Check for plastic consistency
C -----------------------------
      SQRJ2T=SQRT(VARJ2T)
      APTMC=ETA*PT-XI*COHE
      PHI=SQRJ2T+APTMC
      RES=PHI
      IF(APTMC.NE.R0)RES=RES/ABS(APTMC)
      IF(RES.GT.TOL)THEN
C Plastic step: Use return mapping
C ================================
        IFPLAS=.TRUE.
        TWOVEC=.FALSE.
C Apply one-vector return mapping first (return to smooth cone wall)
C ------------------------------------------------------------------
        AUX=SQRT(R1D3+R2D9*ETABAR*ETABAR)
        DO 20 IPTER1=1,MAXRT
C Compute residual derivative
          DENOM=-GMODU-BULK*ETABAR*ETA-
     1           XI*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))*AUX
C Compute Newton-Raphson increment and update variable DGAMA
          DDGAMA=-PHI/DENOM
          DGAMA=DGAMA+DDGAMA
C Compute new residual
          EPBAR=EPBARN+AUX*DGAMA
          COHE=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          SQRJ2=SQRJ2T-GMODU*DGAMA
          P=PT-BULK*ETABAR*DGAMA
          APMC=ETA*P-XI*COHE
          PHI=SQRJ2+APMC
C Check convergence
          RESNOR=ABS(PHI)
          IF(APMC.NE.R0)RESNOR=RESNOR/ABS(APMC)
          IF(RESNOR.LE.TOL)THEN
C Check validity of one-vector return
            IF(SQRJ2.GE.R0)THEN
C results are valid, update stress components and other variables
              GOTO 50
            ELSE
C one-vector return not valid - go to two vector procedure
              GOTO 30
            ENDIF
          ENDIF
   20   CONTINUE
C failure of stress update procedure
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0002')
        GOTO 999
   30   CONTINUE
C Apply two-vector return mapping (return to APEX)
C -----------------------------------------------------
        TWOVEC=.TRUE.
C Solve for DGAMA first
        DGAMA=SQRJ2T/GMODU
        IF(SQRJ2T.EQ.R0)ROO2D9=SQRT(R2D9)
C and iterate for DGAMB
        BDGAMA=ETABAR*DGAMA
        PT2VEC=PT-BULK*BDGAMA
        DGAMA2=DGAMA*DGAMA
        DEP=SQRT(R1D3*DGAMA2+R2D9*BDGAMA**2)
        EPBAR=EPBARN+DEP
        COHE=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
        RES=ETA*PT2VEC-XI*COHE
        DO 40 IPTER2=1,MAXRT
          IF(SQRJ2T.EQ.R0)THEN
            DENOM=-BULK*ETA-
     1             XI*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))*ROO2D9
          ELSE
            DENOM=-BULK*ETA-XI*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))*
     1             R2D9*(BDGAMA+DGAMB)/DEP
          ENDIF
C Compute Newton-Raphson increment and update variable DGAMB
          DDGAMB=-RES/DENOM
          DGAMB=DGAMB+DDGAMB
C Compute new residual
          DEP=SQRT(R1D3*DGAMA2+R2D9*(BDGAMA+DGAMB)**2)
          EPBAR=EPBARN+DEP
          COHE=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          P=PT2VEC-BULK*DGAMB
          RES=ETA*P-XI*COHE
C Check convergence
          RESNOR=ABS(RES)
          IF(COHE.NE.R0)RESNOR=RESNOR/ABS(XI*COHE)
          IF(RESNOR.LE.TOL)THEN
C update stress components and other variables
            GOTO 50
          ENDIF 
   40   CONTINUE
C failure of stress update procedure
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0002')
        GOTO 999
C Update converged stress components and other state variables
C ------------------------------------------------------------
   50   CONTINUE
C update EPBAR
        RSTAVA(MSTRE+1)=EPBAR
C update stress components
        IF(SQRJ2T.EQ.R0)THEN
          FACTOR=R0
        ELSE
          FACTOR=R1-GMODU*DGAMA/SQRJ2T
        ENDIF
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
      LALGVA(3)=TWOVEC
      DGAM(1)=DGAMA
      DGAM(2)=DGAMB
      RETURN
      END
CDOC END_SUBROUTINE SUDP
