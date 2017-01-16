CDOC BEGIN_SUBROUTINE SUEL
CDOC State update procedure for the linear elastic material model.
CDOC
CDOC Given the total strain, this routine computes the corresponding
CDOC stress using the standard generalised Hooke's law for linear
CDOC elastic materials. This routine contains the plane stress, plane
CDOC strain and axisymmetric implementations of the linear elastic
CDOC model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NTYPE  >  Stress state type. Present routine is
CDOC C                          compatible with \smparm{NTYPE=1} (plane
CDOC C                          stress), \smparm{NTYPE=2}
CDOC C                          (plane strain) and \smparm{NTYPE=3}
CDOC C                          (axisymmetric condition).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA <  Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          For the linear elastic model, this array
CDOC C                          stores the (engineering) strain
CDOC C                          components.
CDOC DOUBLE_PRECISION STRAN  >  Array of current total (engineering)
CDOC C                          strain components.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, September 1996; Initial coding
CDOC
      SUBROUTINE SUEL
     1(   NTYPE      ,RPROPS     ,RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MSTRE=4)
      DIMENSION
     1    RPROPS(*)          ,RSTAVA(MSTRE)      ,STRAN(*)           ,
     2    STRES(*)    
      DIMENSION
     1    EED(MSTRE)
      DATA
     1    RP5  ,R2   ,R3   ,R4   / 
     2    0.5D0,2.0D0,3.0D0,4.0D0/
C***********************************************************************
C STATE UPDATE PROCEDURE FOR LINEAR ELASTIC MATERIAL
C***********************************************************************
C
C Set shear and bulk modulus
C
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C
C Decompose strain into deviatoric and volumetric components
C ----------------------------------------------------------
C

      R2G=R2*GMODU
      IF(NTYPE.EQ.1)THEN
C for plane stress
        R4G=R4*GMODU
        R4GD3=R4G/R3
        FACTOR=R2G/(BULK+R4GD3)
        EEV=(STRAN(1)+STRAN(2))*FACTOR
        EEVD3=EEV/R3
        EED(1)=STRAN(1)-EEVD3
        EED(2)=STRAN(2)-EEVD3
      ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
C for plane strain and axisymmetric cases
        EEV=STRAN(1)+STRAN(2)+STRAN(4)
        EEVD3=EEV/R3
        EED(1)=STRAN(1)-EEVD3
        EED(2)=STRAN(2)-EEVD3
        EED(4)=STRAN(4)-EEVD3
      ELSE
        CALL ERRPRT('EI0018')
      ENDIF
C Convert engineering shear component into physical component
      EED(3)=STRAN(3)*RP5
C
C Update stress using linear elastic law
C ---------------------------------------
C
C hydrostatic stress
      P=BULK*EEV
C stress tensor components
      STRES(1)=R2G*EED(1)+P
      STRES(2)=R2G*EED(2)+P
      STRES(3)=R2G*EED(3)
      IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)STRES(4)=R2G*EED(4)+P
C
C Store elastic engineering strain in RSTAVA
C ------------------------------------------
C
      RSTAVA(1)=STRAN(1)
      RSTAVA(2)=STRAN(2)
      RSTAVA(3)=STRAN(3)
      IF(NTYPE.EQ.1)THEN
        R3BULK=R3*BULK
        RSTAVA(4)=-(STRAN(1)+STRAN(2))*(R3BULK-R2G)/(R3BULK+R4G)
      ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        RSTAVA(4)=STRAN(4)
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SUEL