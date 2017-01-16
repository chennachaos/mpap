CDOC BEGIN_SUBROUTINE SUDPPN
CDOC State update procedure for the Drucker-Prager model. Plane stress.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the
CDOC Drucker-Prager
CDOC elasto-plastic material model with general non-linear (piece-wise
CDOC linear) isotropic hardening under plane stress condition.
CDOC The algorithm used here is based on the nested iteration approach
CDOC for enforcement of the plane stress constraint at the Gauss point
CDOC level.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION RALGVA <  Array of real algorithmic variables.
CDOC C                          For the present plane stress
CDOC C                          implementation, it contains the
CDOC C                          incremental plastic multipliers
CDOC C                          and the elastic trial thickness strain
CDOC C                          obtained as the solution of the plane
CDOC C                          stress enforcement loop.
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
CDOC C                          state update procedure has failed
CDOC C                          to converge.
CDOC C                          \smparm{TWOVEC} is set to
CDOC C                          \smparm{.TRUE.}
CDOC C                          if the selected return mapping is the
CDOC C                          two-vector return to the apex. It is
CDOC C                          set to \smparm{.FALSE.} otherwise.
CDOC INTEGER          NTYPE  >  Stress state type. Present routine is
CDOC C                          compatible only with \smparm{NTYPE=1}
CDOC C                          (plane stress).
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
CDOC DOUBLE_PRECISION STRAN  <> Array of elastic trial
CDOC C                          (engineering) strain components.
CDOC C                          Its first three components are the
CDOC C                          in-plane elastic trial components which
CDOC C                          are not updated in the present routine.
CDOC C                          Its fourth component - the thickness
CDOC C                          elastic trial strain - is determined
CDOC C                          here as the solution to the plane
CDOC C                          stress enforcement N-R loop.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, January 1999: Initial coding
CDOC
      SUBROUTINE SUDPPN
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=7  ,MSTRE=4)
C Arguments
      LOGICAL LALGVA(3)
      DIMENSION
     1    RALGVA(3)          ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)    ,STRAN(MSTRE)       ,STRES(MSTRE)
C Local arrays and variables
      LOGICAL EPFLAG ,IFPLAS ,SUFAIL
      DIMENSION
     1    DMATX(MSTRE,MSTRE) ,RSTAUX(MSTRE+1)
      DATA
     1    R0    ,TOL   / 
     2    0.0D0 ,1.D-08/
      DATA MXITER / 20 /
C***********************************************************************
C STATE UPDATE PROCEDURE FOR THE DRUCKER-PRAGER ELASTO-PLASTIC MODEL
C WITH NON-LINEAR (PIECEWISE LINEAR) ISOTROPIC HARDENING IN PLANE
C STRESS:
C               ---  NESTED ITERATION APPROACH  ---
C***********************************************************************
C Stop program if not plane stress
      IF(NTYPE.NE.1)CALL ERRPRT('EI0038')
C Initialise the state update failure flag
      SUFAIL=.FALSE.
C Set some material properties
      NHARD=IPROPS(3)
C
C Begin Newton-Raphson iteration loop for plane stress enforcement
C ----------------------------------------------------------------
C
C Set initial guess for elastic trial thickness strain. Use previously
C converged elastic thickness strain.
      E33TRL=RSTAVA(4)
C Start N-R loop
      DO 20 ITER=1,MXITER
C Set state variables to values at beginning of increment
        DO 10 I=1,MSTRE+1
          RSTAUX(I)=RSTAVA(I)
   10   CONTINUE
C Use axisymmetric integration algorithm to compute stresses, etc.
        STRAN(4)=E33TRL
        CALL SUDP
     1(   RALGVA     ,IPROPS     ,LALGVA     ,3          ,RPROPS     ,
     2    RSTAUX     ,STRAN      ,STRES      )
        SUFAIL=LALGVA(2)
        IF(SUFAIL)THEN
C... emergency exit in case of failure of the state apdate procedure
          GOTO 999
        ENDIF
        IFPLAS=LALGVA(1)
C Check plane stress convergence
        EPBAR=RSTAVA(MSTRE+1)
        COHE=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
        RES=ABS(STRES(4))
C...use normalised out-of-plane stress
        IF(COHE.NE.R0)RES=RES/ABS(COHE)
        IF(RES.LE.TOL)THEN
C...and break N-R loop in case of convergence
          GOTO 30
        ENDIF
C Compute axisymmetric consistent tangent components
        EPFLAG=IFPLAS
        CALL CTDP
     1(   RALGVA     ,DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,
     2    3          ,RPROPS     ,RSTAUX     ,STRAN      )
C Apply Newton-Raphson correction to normal elastic trial strain
        D22=DMATX(4,4)
        E33TRL=E33TRL-STRES(4)/D22
   20 CONTINUE
C Emergency exit in case of failure of the plane stress enforcement loop
      SUFAIL=.TRUE.
      LALGVA(2)=SUFAIL
      CALL ERRPRT('WE0016')
      GOTO 999
   30 CONTINUE
C Set state variables to current updated values
      DO 40 I=1,MSTRE+1
        RSTAVA(I)=RSTAUX(I)
   40 CONTINUE
C Store the converged elastic trial thickness strain in the array of
C real algorithmic variables
      RALGVA(3)=E33TRL
C
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE SUDPPN
