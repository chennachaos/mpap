CDOC BEGIN_SUBROUTINE SWEL
CDOC Initialise/switch state variables for the elastic material model
CDOC
CDOC This initialises and switches state variables (between current and
CDOC previous values) for the elastic material model (Hencky material
CDOC in large strain analysis).
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          MODE   >  Initialisation/Switching mode.
CDOC INTEGER          NLARGE >  Large strain flag.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RSTAVC <> Array of real state variables at Gauss
CDOC C                          point. Current values.
CDOC DOUBLE_PRECISION RSTAVL <> Array of real state variables at Gauss
CDOC C                          point. Last converged (equilibrium)
CDOC C                          values.
CDOC DOUBLE_PRECISION STRESC <> Array of stress (Cauchy in large strain)
CDOC C                          components. Current values.
CDOC DOUBLE_PRECISION STRESL <> Array of stress (Cauchy in large strain)
CDOC C                          components. Last converged (equilibrium)
CDOC C                          values.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, July      1999: Initial coding
CDOC
      SUBROUTINE SWEL
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Arguments
      DIMENSION
     1    RSTAVC(*)          ,RSTAVL(*)          ,STRESC(*)          ,
     2    STRESL(*)
C Local numerical constants
      DATA R0   ,R1   /
     1     0.0D0,1.0D0/
C***********************************************************************
C INITIALISE/SWITCH DATA FOR THE ELASTIC MATERIAL MODEL
C
C    MODE=0 -> Initialises the relevant data.
C
C    MODE=1 -> Assigns current values of the state variables to
C              converged solution (when the current iteration
C              satisfies the convergence criterion).
C
C    MODE=2 -> Assigns the last converged solution to current state
C              variables values (when a new iteration is required by
C              the iterative process).
C
C    MODE=3 -> Assigns the last converged solution to current state
C              variables values (when increment cutting is required).
C***********************************************************************
C
      IF(NTYPE.EQ.1.OR.NTYPE.EQ.2.OR.NTYPE.EQ.3)NSTRE=4
      IF(NTYPE.EQ.5)NSTRE=6
C
      IF(MODE.EQ.0)THEN
C Initialisation mode
C ===================
        CALL RVZERO(STRESC,NSTRE)
        IF(NLARGE.EQ.1)THEN
C Large strain analysis
C ---------------------
C RSTAVA stores the left Cauchy-Green strain tensor components
          IF(NTYPE.EQ.1.OR.NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
            RSTAVC(1)=R1
            RSTAVC(2)=R1
            RSTAVC(3)=R0
            RSTAVC(4)=R1
          ENDIF
        ELSE
C Small strain analysis
C ---------------------
C RSTAVA stores the infinitesimal egineering strain tensor components
          CALL RVZERO(RSTAVC,NSTRE)
        ENDIF
      ELSE
C Switching mode
C ==============
        IF(MODE.EQ.1)THEN
          DO 10 ISTRE=1,NSTRE
            STRESL(ISTRE)=STRESC(ISTRE)
            RSTAVL(ISTRE)=RSTAVC(ISTRE)
   10     CONTINUE
        ELSEIF(MODE.EQ.2.OR.MODE.EQ.3)THEN
          DO 20 ISTRE=1,NSTRE
            STRESC(ISTRE)=STRESL(ISTRE)
            RSTAVC(ISTRE)=RSTAVL(ISTRE)
   20     CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE SWEL
