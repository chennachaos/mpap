CDOC BEGIN_SUBROUTINE SUDAMA
CDOC State update procedure for Lemaitre's ductile damage model.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for Lemaitre's
CDOC ductile damage elasto-plastic model with general non-linear
CDOC (piece-wise linear) isotropic hardening under plane strain and
CDOC axisymmetric conditions.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAMA  <  Incremental plastic multiplier.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDDAMA}.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic yielding
CDOC C                          flag and the return algorithm failure
CDOC C                          flag. The plastic yielding flag is set
CDOC C                          to \smparm{.TRUE.} if plastic yielding
CDOC C                          has occurred and to \smparm{.FALSE.} if
CDOC C                          the step is elastic. The algorithm
CDOC C                          failure flag is set to {\tt .FALSE.} if
CDOC C                          the state update algorithm has been
CDOC C                          successful and to \smparm{.TRUE.} if the
CDOC C                          return mapping algorithm has failed
CDOC C                          to converge.
CDOC INTEGER          NTYPE  >  Stress state type. Present routine is
CDOC C                          compatible only with \smparm{NTYPE=2}
CDOC C                          (plane strain) and \smparm{NTYPE=3}
CDOC C                          (axisymmetric condition).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, the two damage
CDOC C                          evolution constants and the pairs
CDOC C                          ``hardening variable strain-uniaxial
CDOC C                          yield stress'' defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          \smparm{RDDAMA}.
CDOC DOUBLE_PRECISION RSTAVA <> Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Previous converged values on entry,
CDOC C                          updated values on exit.
CDOC C                          The state variables stored in
CDOC C                          this array are the (engineering)
CDOC C                          elastic strain components, the
CDOC C                          hardening variable and the
CDOC C                          damage internal variable.
CDOC DOUBLE_PRECISION STRAN  >  Array of elastic trial (engineering)
CDOC C                          strain components.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, Jan 2000: Initial coding
CDOC
      SUBROUTINE SUDAMA
     1(   DGAMA      ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=6  ,MSTRE=4)
      LOGICAL IFPLAS, LALGVA(2), SUFAIL
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)          ,RSTAVA(MSTRE+2)    ,
     2    STRAN(MSTRE)       ,STRES(MSTRE)
      DIMENSION
     1    EET(MSTRE)
      DATA
     1    R0   ,RP5  ,R1   ,R2   ,R3   ,R6   ,SMALL ,TOL   / 
     2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,6.0D0,1.D-20,1.D-08/
      DATA MXITER / 25 /
C***********************************************************************
C STATE UPDATE PROCEDURE FOR LEMAITRE'S DUCTILE DAMAGE MATERIAL MODEL
C WITH NON-LINEAR (PIECEWISE LINEAR) ISOTROPIC HARDENING:
C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM.
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C***********************************************************************
C Stop program if neither plane strain nor axisymmetric state
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0051')
C Initialise some algorithmic and internal variables
      DGAMA=R0
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
C Retrieve hardening and damage internal variables
      HVARN=RSTAVA(MSTRE+1)
      DAMAGN=RSTAVA(MSTRE+2)
C... integrity
      OMEGAN=R1-DAMAGN
C Retrieve some material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      DAMEXP=RPROPS(4)
      DAMDEN=RPROPS(5)
      NHARD=IPROPS(3)
C Shear and bulk moduli and other necessary constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2BULK=R2*BULK
      R2G=R2*GMODU
      R3G=R3*GMODU
      R6G=R6*GMODU
C Elastic predictor: Compute elastic trial state
C ==============================================
C Volumetric strain and (undamaged) pressure stress
      EEV=STRAN(1)+STRAN(2)+STRAN(4)
      PTRIAL=BULK*EEV
C Elastic trial deviatoric strain
      EEVD3=EEV/R3
      EET(1)=STRAN(1)-EEVD3
      EET(2)=STRAN(2)-EEVD3
      EET(4)=STRAN(4)-EEVD3
C Convert engineering shear component into physical component
      EET(3)=STRAN(3)/R2
C Compute trial (undamaged) von Mises effective stress and uniaxial
C yield stress
      VARJ2T=R2G*R2G*(EET(3)*EET(3)+RP5*(EET(1)*EET(1)+
     1                     EET(2)*EET(2)+EET(4)*EET(4)))
      QTRIAL=SQRT(R3*VARJ2T)
      SIGMAY=PLFUN(HVARN,NHARD,RPROPS(IPHARD))
C Check for plastic admissibility
C ===============================
      PHI=QTRIAL-SIGMAY
      IF(PHI/SIGMAY.GT.TOL)THEN
C Plastic step: Apply return mapping - use Newton-Raphson algorithm
C               to solve the return mapping equation for DGAMA
C =================================================================
C Reset plastic flag
        IFPLAS=.TRUE.
C Initial guess for DGAMA: Use perfectly plastic solution with frozen
C yield surface at the beginning of the load increment
        DGAMA=OMEGAN*PHI/R3G
C Initialise hardening variable
        HVAR=HVARN+DGAMA
C Start N-R iterations
C --------------------
        PTRIA2=PTRIAL**2
        DO 10 NRITER=1,MXITER
C yield stress
          SIGMAY=PLFUN(HVAR,NHARD,RPROPS(IPHARD))
C integrity
          OMEGA=R3G*DGAMA/(QTRIAL-SIGMAY)
          IF(ABS(OMEGA).LT.SMALL)THEN
C singular residual: Return map fails. Reset failure flag, break N-R
C                    loop and exit
            SUFAIL=.TRUE.
           CALL ERRPRT('WE0018')
           GOTO 999
         ENDIF
C stress triaxiality and damage energy release rate
          SIGMA2=SIGMAY**2
          Y=-SIGMA2/R6G-PTRIA2/R2BULK
C Compute residual function
          RES=OMEGA-OMEGAN+DGAMA/OMEGA*(-Y/DAMDEN)**DAMEXP
C Check for convergence
*************************************
**          write(*,*)'res = ',res
*************************************
          IF(ABS(RES).LE.TOL)THEN
C... update hardening and damage variables
**        write(*,*)'SU  domega= ',domega
**        write(*,*)'SU  dy    = ',dy  
**        write(*,*)'SU  dres  = ',dres
            DAMAGE=R1-OMEGA
            IF(DAMAGE.GT.R1)THEN
C... check if converged damage variable is physically acceptable
              SUFAIL=.TRUE.
              CALL ERRPRT('WE0019')
            ENDIF
            RSTAVA(MSTRE+1)=HVAR
            RSTAVA(MSTRE+2)=DAMAGE
C... update stress components
            P=OMEGA*PTRIAL
            Q=SIGMAY*OMEGA
            FACTOR=R2G*Q/QTRIAL
            STRES(1)=FACTOR*EET(1)+P
            STRES(2)=FACTOR*EET(2)+P
            STRES(3)=FACTOR*EET(3)
            STRES(4)=FACTOR*EET(4)+P
C... compute and store converged elastic (engineering) strain components
            FACTOR=R1-R3G*DGAMA/(OMEGA*QTRIAL)
            RSTAVA(1)=FACTOR*EET(1)+EEVD3
            RSTAVA(2)=FACTOR*EET(2)+EEVD3
            RSTAVA(3)=FACTOR*EET(3)*R2
            RSTAVA(4)=FACTOR*EET(4)+EEVD3
            GOTO 999
          ENDIF
C Compute derivative of residual function
C... slope of hardening function
          HSLOPE=DPLFUN(HVAR,NHARD,RPROPS(IPHARD))
C... derivative of OMEGA and Y
          OMEGA2=OMEGA**2
          DOMEGA=(R3G+OMEGA*HSLOPE)/(QTRIAL-SIGMAY)
          DY=-HSLOPE*SIGMAY/R3G
C... residual derivative
          AUX=-Y/DAMDEN
          DRES=AUX**DAMEXP/OMEGA+DOMEGA-(DGAMA*AUX**DAMEXP*DOMEGA)/
     1         OMEGA2-(DAMEXP*DGAMA*AUX**(DAMEXP-R1)*DY)/(DAMDEN*OMEGA)
C Apply N-R correction to DGAMA
          DDGAMA=-RES/DRES
          DGAMA=DGAMA+DDGAMA
C... update hardening variable
          HVAR=HVARN+DGAMA
   10   CONTINUE
C N-R loop failed to converge: Reset failure flag and issue warning
C                              message before exiting
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0018')
      ELSE
C Elastic step: Update stress using damaged elastic law
C =====================================================
        FACTOR=R2G*OMEGAN
        P=OMEGAN*PTRIAL
        STRES(1)=FACTOR*EET(1)+P
        STRES(2)=FACTOR*EET(2)+P
        STRES(3)=FACTOR*EET(3)
        STRES(4)=FACTOR*EET(4)+P
C elastic engineering strain
        RSTAVA(1)=STRAN(1)
        RSTAVA(2)=STRAN(2)
        RSTAVA(3)=STRAN(3)
        RSTAVA(4)=STRAN(4)
      ENDIF
  999 CONTINUE
C Update some algorithmic variables before exit
      LALGVA(1)=IFPLAS
      LALGVA(2)=SUFAIL
      RETURN
      END
CDOC END_SUBROUTINE SUDAMA

