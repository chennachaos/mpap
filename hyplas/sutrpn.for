CDOC BEGIN_SUBROUTINE SUTRPN
CDOC State update procedure for the Tresca model in plane stress.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the
CDOC Tresca elasto-plastic material model with piece-wise linear
CDOC isotropic hardening under plane stress condition.
CDOC The algorithm used here is based on the nested iteration approach
CDOC for enforcement of the plane stress constraint at the Gauss point
CDOC level.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   <  Array of incremental plastic
CDOC C                          multipliers.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDTR}.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic yielding
CDOC C                          flag, \smparm{IFPLAS}; the return
CDOC C                          algorithm failure flag, \smparm{SUFAIL};
CDOC C                          the two-vector return flag,
CDOC C                          \smparm{TWOVEC} and
CDOC C                          the right corner return flag,
CDOC C                          \smparm{RIGHT}.
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
CDOC C                          if the selected return mapping is to
CDOC C                          a corner (right or left) and is set to
CDOC C                          \smparm{.FALSE.} otherwise.
CDOC C                          \smparm{RIGHT} is set to \smparm{.TRUE.}
CDOC C                          if the selected return is to the right
CDOC C                          corner and is set to \smparm{.FALSE.}
CDOC C                          otherwise.
CDOC INTEGER          NTYPE  >  Stress state type.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, and the plastic
CDOC C                          properties: the pairs
CDOC C                          ``accumulated plastic strain-cohesion''
CDOC C                          defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          \smparm{RDTR}.
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
CDOC E.de Souza Neto, Sept 1998: Initial coding
CDOC
      SUBROUTINE SUTRPN
     1(   DGAM       ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=4  ,MSTRE=4)
C Arguments
      LOGICAL
     1    LALGVA(4)
      DIMENSION
     1    DGAM(2)            ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)    ,STRAN(*)           ,STRES(*)
C Local arrays and variables
      LOGICAL EPFLAG, IFPLAS, SUFAIL
      DIMENSION
     1    DMATX(MSTRE,MSTRE) ,RSTAUX(MSTRE+1)
      DATA
     1    TOL   / 
     2    1.D-08/
      DATA MXITER / 20 /
C***********************************************************************
C STRESS UPDATE PROCEDURE FOR TRESCA TYPE ELASTO-PLASTIC MATERIAL WITH
C PIECE-WISE LINEAR ISOTROPIC HARDENING IN PLANE STRESS:
C               ---  NESTED ITERATION APPROACH  ---
C***********************************************************************
C Stops program if not plane stress
cccccccccccc      IF(NTYPE.NE.1)CALL ERRPRT('EI0029')
C
C Newton-Raphson iteration loop for plane stress enforcement
C
C Set initial guess for elastic trial thickness strain. Use previously
C converged elastic thickness strain.
      E33TRL=RSTAVA(4)
C Set axisymmetric state flag
      NTYPAX=3
C Start N-R loop
      DO 20 ITER=1,MXITER
C Set state variables to values at beginning of increment
        DO 10 I=1,MSTRE+1
          RSTAUX(I)=RSTAVA(I)
   10   CONTINUE
C Use axisymmetric integration algorithm to compute stresses
        STRAN(4)=E33TRL
ccc        stran(1)=-0.3d-4
ccc        stran(2)=0.00003d0
ccc        stran(3)=0.d0
ccc        stran(4)=0.d0
c        write(*,*)(RSTAUX(i),i=1,mstre+1)
c        write(*,*)stran(1),stran(2),stran(3),stran(4)
        CALL SUTR
     1(   DGAM       ,IPROPS     ,LALGVA     ,NTYPAX     ,RPROPS     ,
     2    RSTAUX     ,STRAN      ,STRES      )
        IFPLAS=LALGVA(1)
C Check plane stress convergence
cccccccccccccccc
c        write(*,*)stres(1),stres(2),stres(3)
c        if(ifplas)write(*,*)stres(4),' plast'
c        if(.not.ifplas)write(*,*)stres(4)
cccccccccccccccc
        IF(ABS(STRES(4)).LE.TOL)THEN
C...and break N-R loop in case of convergence
ccccccccc
c          write(*,*)'         OK'
          GOTO 30
        ENDIF
C Compute axisymmetric consistent tangent components
        EPFLAG=IFPLAS
        CALL CTTR
     1(   DGAM       ,DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,
     2    NTYPAX     ,RPROPS     ,RSTAUX     ,STRAN      ,STRES      )
C Apply Newton-Raphson correction to normal elastic trial strain
        D22=DMATX(4,4)
        E33TRL=E33TRL-STRES(4)/D22
   20 CONTINUE
ccccccccc
c      write(*,*)'***********'
      stop
ccccccccc
   30 CONTINUE
C Set state variables to current updated values
      DO 40 I=1,MSTRE+1
        RSTAVA(I)=RSTAUX(I)
   40 CONTINUE
ccccccccc
      RSTAVA(4)=E33TRL
ccccccccc
c        write(*,*)' This is the output ',stres(1),stres(2),stres(3)
      RETURN
      END
CDOC END_SUBROUTINE SUTRPN