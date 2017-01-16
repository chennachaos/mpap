CDOC BEGIN_SUBROUTINE SUVMMX
CDOC State update for the von Mises model with mixed hardening.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the von Mises
CDOC elasto-plastic material model with general non-linear (piece-wise
CDOC linear) mixed isotropic/kinematic hardening under plane strain and
CDOC axisymmetric conditions.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAMA  <  Incremental plastic multiplier.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curves is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines
CDOC C                          \smparm{MATIRD} and \smparm{RDVMMX}.
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
CDOC C                          Poisson's ratio, the pairs
CDOC C                          ``accumulated plastic strain-isotropic
CDOC C                          hard. stress'' defining the isotropic
CDOC C                          hardening curve and the pairs
CDOC C                          ``accumulated plastic strain-kinematic
CDOC C                          hardening stress'' defining the
CDOC C                          kinematic hardening curve.
CDOC C                          This array is set in routine
CDOC C                          \smparm{RDVMMX}.
CDOC DOUBLE_PRECISION RSTAVA <> Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Previous converged values on entry,
CDOC C                          updated values on exit.
CDOC C                          The state variables stored in
CDOC C                          this array are the (engineering)
CDOC C                          elastic strain components, the
CDOC C                          accumulated plastic strain and the
CDOC C                          backstress tensor components.
CDOC DOUBLE_PRECISION STRAN  >  Array of elastic trial (engineering)
CDOC C                          strain components.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1999: Initial coding
CDOC
      SUBROUTINE SUVMMX
     1(   DGAMA      ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=4  ,MSTRE=4)
C Arguments
      LOGICAL LALGVA(2)
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)          ,RSTAVA(2*MSTRE+1)  ,
     2    STRAN(MSTRE)       ,STRES(MSTRE)
C Local arrays and variables
      LOGICAL IFPLAS, SUFAIL
      DIMENSION
     1    BACSTN(MSTRE)      ,BACSTR(MSTRE)      ,EET(MSTRE)         ,
     2    ETATRL(MSTRE)      ,FLOVEC(MSTRE)      ,STRIAL(MSTRE)
      DATA
     1    R0   ,R1   ,R2   ,R3   ,TOL   / 
     2    0.0D0,1.0D0,2.0D0,3.0D0,1.D-08/
      DATA MXITER / 50 /
C***********************************************************************
C STATE UPDATE PROCEDURE FOR THE VON MISES ELASTO-PLASTIC MATERIAL MODEL
C WITH NON-LINEAR (PIECEWISE LINEAR) MIXED ISOTROPIC/KINEMATIC
C HARDENING:
C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM (Box 7.5).
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS ONLY.
C***********************************************************************
C Stop program if neither plane strain nor axisymmetric state
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0047')
C Retrieve some previous converged state variables
C... accumulated plastic strain (hardening internal variable)
      EPBARN=RSTAVA(5)
C... backstress tensor components
      BACSTN(1)=RSTAVA(6)
      BACSTN(2)=RSTAVA(7)
      BACSTN(3)=RSTAVA(8)
      BACSTN(4)=RSTAVA(9)
C Initialise some algorithmic and internal variables
      DGAMA=R0
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
C Retrieve some material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      NHARD=IPROPS(3)
C Set pointers to isotropic and kinematic hardening curves
      IPIHAR=IPHARD
      IPKHAR=IPHARD+2*NHARD
C Shear and bulk moduli and other necessary constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R3G=R3*GMODU
C Elastic predictor: Compute elastic trial state [item(i) Box 7.5]
C ----------------------------------------------------------------
C Volumetric strain and pressure stress
      EEV=STRAN(1)+STRAN(2)+STRAN(4)
      P=BULK*EEV
C Elastic trial deviatoric strain
      EEVD3=EEV/R3
      EET(1)=STRAN(1)-EEVD3
      EET(2)=STRAN(2)-EEVD3
      EET(4)=STRAN(4)-EEVD3
C Convert engineering shear component into physical component
      EET(3)=STRAN(3)/R2
C Compute trial deviatoric stress and trial relative stress
      DO 10 ISTRE=1,MSTRE
        STRIAL(ISTRE)=R2G*EET(ISTRE)
        ETATRL(ISTRE)=STRIAL(ISTRE)-BACSTN(ISTRE)
   10 CONTINUE
C Compute trial effective relative stress
      QBARTR=SQRT(R3/R2*(ETATRL(1)**2+ETATRL(2)**2+R2*ETATRL(3)**2+
     1            ETATRL(4)**2))
C and radius of von Mises cylinder
      SIGMAY=PLFUN(EPBARN,NHARD,RPROPS(IPIHAR))
C Check for plastic admissibility [item (ii) Box 7.5]
C ---------------------------------------------------
      PHI=QBARTR-SIGMAY
      IF(PHI/SIGMAY.GT.TOL)THEN
C Plastic step: Apply return mapping - use Newton-Raphson algorithm to
C               solve the return mapping equation [item (iii) Box 7.5]
C --------------------------------------------------------------------
        IFPLAS=.TRUE.
        EPBAR=EPBARN
        BETBAN=PLFUN(EPBARN,NHARD,RPROPS(IPKHAR))
        DO 40 NRITER=1,MXITER
C Compute residual derivative
          HISLOP=DPLFUN(EPBAR,NHARD,RPROPS(IPIHAR))
          HKSLOP=DPLFUN(EPBAR,NHARD,RPROPS(IPKHAR))
          DENOM=-R3G-HISLOP-HKSLOP
C Compute Newton-Raphson increment and update variable DGAMA
          DDGAMA=-PHI/DENOM
          DGAMA=DGAMA+DDGAMA
C Compute new residual
          EPBAR=EPBAR+DDGAMA
          SIGMAY=PLFUN(EPBAR,NHARD,RPROPS(IPIHAR))
          BETBAR=PLFUN(EPBAR,NHARD,RPROPS(IPKHAR))
          PHI=QBARTR-R3G*DGAMA-BETBAR+BETBAN-SIGMAY
C Check convergence
          RESNOR=ABS(PHI/SIGMAY)
          IF(RESNOR.LE.TOL)THEN
C update stress components
            ETANOR=SQRT(ETATRL(1)**2+ETATRL(2)**2+R2*ETATRL(3)**2+
     1                  ETATRL(4)**2)
            FACTOR=SQRT(R3/R2)/ETANOR
            DO 20 ISTRE=1,MSTRE
              FLOVEC(ISTRE)=FACTOR*ETATRL(ISTRE)
   20       CONTINUE
            FACTOR=R2G*DGAMA
ccc            write(*,*)'strial:',strial
ccc            write(*,*)'flovec:',flovec
ccc            write(*,*)'factor:',factor
            STRES(1)=STRIAL(1)-FACTOR*FLOVEC(1)+P
            STRES(2)=STRIAL(2)-FACTOR*FLOVEC(2)+P
            STRES(3)=STRIAL(3)-FACTOR*FLOVEC(3)
            STRES(4)=STRIAL(4)-FACTOR*FLOVEC(4)+P
C compute and store converged elastic (engineering) strain components
            RSTAVA(1)=EET(1)-DGAMA*FLOVEC(1)+EEVD3
            RSTAVA(2)=EET(2)-DGAMA*FLOVEC(2)+EEVD3
            RSTAVA(3)=R2*(EET(3)-DGAMA*FLOVEC(3))
            RSTAVA(4)=EET(4)-DGAMA*FLOVEC(4)+EEVD3
C store updated accumulated plastic strain
            RSTAVA(5)=EPBAR
C compute and store updated backstress tensor components
            FACTOR=R2/R3*(BETBAR-BETBAN)
            DO 30 ISTRE=1,MSTRE
              BACSTR(ISTRE)=BACSTN(ISTRE)+FACTOR*FLOVEC(ISTRE)
   30       CONTINUE
            RSTAVA(6)=BACSTR(1)
            RSTAVA(7)=BACSTR(2)
            RSTAVA(8)=BACSTR(3)
            RSTAVA(9)=BACSTR(4)
            GOTO 999
          ENDIF
   40   CONTINUE
C reset failure flag and print warning message if the algorithm fails
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0017')
      ELSE
C Elastic step: Update stress using linear elastic law
C ----------------------------------------------------
        STRES(1)=R2G*EET(1)+P
        STRES(2)=R2G*EET(2)+P
        STRES(3)=R2G*EET(3)
        STRES(4)=R2G*EET(4)+P
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
CDOC END_SUBROUTINE SUVMMX