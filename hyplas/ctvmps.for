CDOC BEGIN_SUBROUTINE CTVMPS
CDOC Consistent tangent matrix for von Mises model in plane stress.
CDOC
CDOC This routine computes the tangent matrix consistent with the
CDOC plane stress-projected implicit elastic predictor/return mapping
CDOC algorithm for the von Mises model coded in subroutine
CDOC \smparm{SUVMPS}.
CDOC It returns either the elastic tangent or the elasto-plastic
CDOC consistent tangent matrix depending on the input value of
CDOC the logical argument \smparm{EPFLAG}.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAMA  >  Incremental plastic multiplier (outcome
CDOC C                          of \smparm{SUVMPS}.
CDOC DOUBLE_PRECISION DMATX  <  Plane stress consistent tangent matrix.
CDOC LOGICAL          EPFLAG >  Elasto-plastic flag.
CDOC C                          If \smparm{.FALSE.}, \smparm{DMATX}
CDOC C                          returns as the elastic matrix (the
CDOC C                          standard plane stress linear elasticity
CDOC C                          matrix). If
CDOC C                          \sparm{.TRUE.}, \smparm{DMATX} returns
CDOC C                          as the elasto-plastic tangent consistent
CDOC C                          with the plane stress-projected return
CDOC C                          mapping algorithm
CDOC C                          implemented in routine \smparm{SUVMPS}.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDVM}.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, and the pairs
CDOC C                          ``accumulated plastic strain-uniaxial
CDOC C                          yield stress'' defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          \smparm{RDVM}.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Current values (outcome of
CDOC C                          \smparm{SUVMPS}).
CDOC C                          The state variables stored in
CDOC C                          this array are the (engineering)
CDOC C                          elastic strain components and the
CDOC C                          accumulated plastic strain.
CDOC DOUBLE_PRECISION STRES  >  Array of current stress tensor
CDOC C                          components (outcome of \smparm{SUVMPS}).
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, Sept 1998: Initial coding
CDOC
      SUBROUTINE CTVMPS
     1(   DGAMA      ,DMATX      ,EPFLAG     ,IPROPS     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=4  ,MSTRE=4)
      LOGICAL EPFLAG
C Array arguments
      DIMENSION
     1    DMATX(MSTRE,MSTRE),IPROPS(*)           ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)   ,STRES(MSTRE)
C Local arrays
      DIMENSION
     1    FOID(MSTRE,MSTRE)  ,SOID(MSTRE)        ,VECN(3)
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
     1    RP5  ,R1   ,R2   ,R3   ,R4   /
     2    0.5D0,1.0D0,2.0D0,3.0D0,4.0D0/
C***********************************************************************
C COMPUTATION OF THE CONSISTENT TANGENT MODULUS FOR VON MISES TYPE
C ELASTO-PLASTIC MATERIAL WITH PIECE-WISE LINEAR ISOTROPIC HARDENING.
C PLANE STRESS IMPLEMENTATION ONLY.
C***********************************************************************
C Stops program if neither not plane stress
      IF(NTYPE.NE.1)CALL prgError(1,'ctvmps','plane stress only!')
C Current accumulated plastic strain
      EPBAR=RSTAVA(MSTRE+1)
C Set material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      NHARD=IPROPS(3)
C Shear and bulk moduli
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R1D3=R1/R3
      R2D3=R2*R1D3
      IF(EPFLAG)THEN
C Compute elastoplastic consistent tangent (Box 9.6)
C ==================================================
C Item (i):
C ---------
C Compute XI
        XI=R2D3*(STRES(1)*STRES(1)+STRES(2)*STRES(2)-STRES(1)*STRES(2))+
     1     R2*STRES(3)*STRES(3)
C Hardening slope
        HSLOPE=DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
C Matrix E components
        ESTAR1=R3*YOUNG/(R3*(R1-POISS)+YOUNG*DGAMA)
        ESTAR2=R2G/(R1+R2G*DGAMA)
        ESTAR3=GMODU/(R1+R2G*DGAMA)
        E11=RP5*(ESTAR1+ESTAR2)
        E22=E11
        E12=RP5*(ESTAR1-ESTAR2)
        E21=E12
        E33=ESTAR3
C Components of the matrix product EP
        EPSTA1=R1D3*ESTAR1
        EPSTA2=ESTAR2
        EPSTA3=EPSTA2
        EP11=RP5*(EPSTA1+EPSTA2)
        EP22=EP11
        EP12=RP5*(EPSTA1-EPSTA2)
        EP21=EP12
        EP33=EPSTA3
C Vector n
        VECN(1)=EP11*STRES(1)+EP12*STRES(2)
        VECN(2)=EP21*STRES(1)+EP22*STRES(2)
        VECN(3)=EP33*STRES(3)
C Scalar alpha
        DENOM1=STRES(1)*(R2D3*VECN(1)-R1D3*VECN(2))+
     1         STRES(2)*(R2D3*VECN(2)-R1D3*VECN(1))+
     2         STRES(3)*R2*VECN(3)
        DENOM2=R2*XI*HSLOPE/(R3-R2*HSLOPE*DGAMA)
        ALPHA=R1/(DENOM1+DENOM2)
C Item (ii): Assemble elasto-plastic tangent
C ------------------------------------------
        DMATX(1,1)=E11-ALPHA*VECN(1)*VECN(1)
        DMATX(1,2)=E12-ALPHA*VECN(1)*VECN(2)
        DMATX(1,3)=-ALPHA*VECN(1)*VECN(3)
        DMATX(2,1)=DMATX(1,2)
        DMATX(2,2)=E22-ALPHA*VECN(2)*VECN(2)
        DMATX(2,3)=-ALPHA*VECN(2)*VECN(3)
        DMATX(3,1)=DMATX(1,3)
        DMATX(3,2)=DMATX(2,3)
        DMATX(3,3)=E33-ALPHA*VECN(3)*VECN(3)
      ELSE
C Compute plane stress elasticity matrix
C ======================================
        NSTRE=3
        R4GD3=R4*GMODU/R3
        FACTOR=(BULK-R2G/R3)*(R2G/(BULK+R4GD3))
        DO 20 I=1,NSTRE
          DO 10 J=I,NSTRE
            DMATX(I,J)=R2G*FOID(I,J)+FACTOR*SOID(I)*SOID(J)
   10     CONTINUE
   20   CONTINUE
C lower triangle
        DO 40 J=1,NSTRE-1
          DO 30 I=J+1,NSTRE
            DMATX(I,J)=DMATX(J,I)
   30     CONTINUE
   40   CONTINUE
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE CTVMPS
