CDOC BEGIN_SUBROUTINE SUOGD
CDOC State update procedure for the Ogden hyperelastic material model.
CDOC
CDOC This routine updates the Cauchy stress and other state variables
CDOC for the regularised (compressible) Ogden hyperelastic material
CDOC model under plane stress, plane strain and axisymmetric states.
CDOC Under plane stress, the model is exactly incompressible rather
CDOC than regularised. Under plane strain and axisymmetric conditions,
CDOC the regularised (penalty type) approach is adopted to enforce
CDOC incompressibility. In the regularised approach, the bulk modulus
CDOC plays the role of the incompressibility penalty factor.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION B      >  Array of components of the Left
CDOC C                          Cauchy-Green strain tensor.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC INTEGER          NTYPE  >  Stress state type flag. This routine
CDOC C                          accepts \smaprm{NTYPE=1} (plane stress)
CDOC C                          \smaprm{NTYPE=2} (plane strain) and
CDOC C                          \smaprm{NTYPE=3} (axisymmetric).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA <  Array of real state variables other than
CDOC C                          stress tensor components. Use here to
CDOC C                          store the left Cauchy-Green strain
CDOC C                          tensor.
CDOC DOUBLE_PRECISION STRES  <  Array of Cauchy stress components.
CDOC DOUBLE_PRECISION THICK  <> Thickness of the current point (previous
CDOC C                          converged value on entry - current
CDOC C                          updated value on return).
CDOC C                          Used only in plane stress analysis.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      SUBROUTINE SUOGD
     1(   B          ,IPROPS     ,NTYPE      ,RPROPS     ,RSTAVA    ,
     2    STRES      ,THICK      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPOGDC=2)
      LOGICAL DUMMY
      PARAMETER
     1(   MCOMP=4    ,MSTRE=4    ,NDIM=2     )
      DIMENSION
     1    B(MCOMP)           ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE)      ,STRES(MSTRE)
      DIMENSION
     1    EIGPRJ(MCOMP,NDIM) ,EIGB(NDIM)         ,PSTRES(3)          ,
     2    PSTRTC(3)
      DATA  R1   ,R3   /
     1      1.0D0,3.0D0/
C***********************************************************************
C STRESS UPDATE PROCEDURE FOR OGDEN TYPE HYPERELASTIC MATERIAL MODEL.
C PLANE STRESS, PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C***********************************************************************
C Set Ogden material constants
C ============================
C Number of terms in Ogden's strain-energy function
      NOGTRM=IPROPS(3)
C Bulk modulus
      BULK=RPROPS(IPOGDC+NOGTRM*2)
C Compute principal stretches
C ===========================
C Perform spectral decomposition of the left Cauchy-Green tensor B 
      CALL SPDEC2
     1(   EIGPRJ     ,EIGB       ,DUMMY      ,B          )
C Compute in-plane principal stretches
      PSTRTC(1)=SQRT(EIGB(1))
      PSTRTC(2)=SQRT(EIGB(2))
C...and out-of-plane stretches
      IF(NTYPE.EQ.1)THEN
        PSTRTC(3)=R1/(PSTRTC(1)*PSTRTC(2))
      ELSEIF(NTYPE.EQ.2)THEN
        PSTRTC(3)=R1
      ELSEIF(NTYPE.EQ.3)THEN
        PSTRTC(3)=SQRT(B(4))
      ENDIF
C Compute principal Kirchhoff stresses
C ====================================
      CALL RVZERO(PSTRES,3)
      IF(NTYPE.EQ.1) THEN
C Plane stress: Exact incompressibility assumed
C ---------------------------------------------
        DO 10 I=1,NOGTRM
          CMU=RPROPS(IPOGDC-1+I*2-1)
          ALPHA=RPROPS(IPOGDC-1+I*2)
          PSTRES(1)=PSTRES(1)+CMU*(PSTRTC(1)**ALPHA-
     1              (PSTRTC(1)*PSTRTC(2))**(-ALPHA))
          PSTRES(2)=PSTRES(2)+CMU*(PSTRTC(2)**ALPHA-
     1              (PSTRTC(1)*PSTRTC(2))**(-ALPHA))
   10   CONTINUE
        DETF=R1
      ELSE IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
C Plane strain and axisymmetric: Regularised Ogden constitutive law
C -----------------------------------------------------------------
C Compute principal deviatoric Kirchhoff stresses
        R1D3=R1/R3
        DETF=PSTRTC(1)*PSTRTC(2)
        IF(NTYPE.EQ.3)DETF=DETF*PSTRTC(3)
        DO 20 I=1,NOGTRM
          CMU=RPROPS(IPOGDC-1+I*2-1)
          ALPHA=RPROPS(IPOGDC-1+I*2)
          FACTOR=R1D3*(PSTRTC(1)**ALPHA+PSTRTC(2)**ALPHA+
     1           PSTRTC(3)**ALPHA)
          FACVOL=DETF**(-ALPHA*R1D3)
          PSTRES(1)=PSTRES(1)+CMU*FACVOL*(PSTRTC(1)**ALPHA-FACTOR)
          PSTRES(2)=PSTRES(2)+CMU*FACVOL*(PSTRTC(2)**ALPHA-FACTOR)
          PSTRES(3)=PSTRES(3)+CMU*FACVOL*(PSTRTC(3)**ALPHA-FACTOR)
   20   CONTINUE
C Add hydrostatic Kirchhoff pressure (incompressibility penalty term)
        PRESS=BULK*LOG(DETF)
        DO 30 I=1,3
          PSTRES(I)=PSTRES(I)+PRESS
  30    CONTINUE
      ENDIF
C Assemble vector of Cauchy stress tensor components
C ==================================================
      CALL RVZERO(STRES,3)
      R1DDET=R1/DETF
      PSTRES(1)=PSTRES(1)*R1DDET
      PSTRES(2)=PSTRES(2)*R1DDET
      DO 50 ICOMP=1,3
        DO 40 IDIR=1,2
          STRES(ICOMP)=STRES(ICOMP)+PSTRES(IDIR)*EIGPRJ(ICOMP,IDIR)
   40   CONTINUE
   50 CONTINUE
      IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)STRES(4)=PSTRES(3)*R1DDET
C Update thickness (plane stress only) and state variables vector RSTAVA
C ======================================================================
      RSTAVA(1)=B(1)
      RSTAVA(2)=B(2)
      RSTAVA(3)=B(3)
      IF(NTYPE.EQ.1)THEN
        THICK=THICK*PSTRTC(3)/(SQRT(RSTAVA(4)))
        RSTAVA(4)=PSTRTC(3)*PSTRTC(3)
      ELSEIF(NTYPE.EQ.2)THEN
        RSTAVA(4)=R1
      ELSEIF(NTYPE.EQ.3)THEN
        RSTAVA(4)=B(4)
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SUOGD