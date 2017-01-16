CDOC BEGIN_SUBROUTINE CTDPPN
CDOC Consistent tangent for the Drucker-Prager model in plane stress.
CDOC
CDOC This routine computes the tangent matrix consistent with the
CDOC nested iteration algorithm for the Drucker-Prager model in
CDOC plane stress coded in subroutine \smparm{SUDPPN}.
CDOC It returns either the elastic tangent or the elasto-plastic
CDOC consistent tangent matrix depending on the input value of
CDOC the logical argument \smparm{EPFLAG}.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION RALGVA >  Array of real algorithmic variables.
CDOC C                          For the present plane stress
CDOC C                          implementation, it contains the
CDOC C                          incremental plastic multipliers obtained
CDOC C                          in routine \smparm{SUDPPN} and the
CDOC C                          elastic trial thickness strain obtained
CDOC C                          as the solution of the plane stress
CDOC C                          enforcement loop of \smparm{SUDPPN}.
CDOC C                          Note that for the first iteration of a
CDOC C                          load increment the incremental plastic
CDOC C                          multipliers must be set to zero and
CDOC C                          elastic trial thickness strain must be
CDOC C                          set equal to the elastic thickness
CDOC C                          strain of the previous equilibrium
CDOC C                          solution.
CDOC DOUBLE_PRECISION DMATX  <  Plane stress consistent tangent matrix.
CDOC LOGICAL          EPFLAG >  Elasto-plastic flag.
CDOC C                          If \smparm{.FALSE.}, \smparm{DMATX}
CDOC C                          returns as the elastic matrix (the
CDOC C                          standard plane stress linear elasticity
CDOC C                          matrix). If
CDOC C                          \sparm{.TRUE.}, \smparm{DMATX} returns
CDOC C                          as the elasto-plastic tangent consistent
CDOC C                          with the nested iteration algorithm
CDOC C                          implemented in routine \smparm{SUDPPN}.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDDP}.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC LOGICAL          LALGVA >  Array of logical algorithmic flags.
CDOC C                          See list of arguments of \smparm{SUDP}.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          Same as in the argument list of
CDOC C                          subroutines
CDOC C                          \smparm{SUDP}-\smparm{SUDPPN}.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Output of \smparm{SUDP}-\smparm{SUDPPN}.
CDOC DOUBLE_PRECISION STRES  >  Array of current stress tensor
CDOC C                          components (outcome of \smparm{SUDPPN}).
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, January 1999: Initial coding
CDOC
      SUBROUTINE CTDPPN
     1(   RALGVA     ,DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,
     2    NTYPE      ,RPROPS     ,RSTAVA     ,STRAN      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(  MSTRE=4  )
C Arguments
      LOGICAL EPFLAG, LALGVA(3)
      DIMENSION
     1    RALGVA(3)          ,DMATX(MSTRE,MSTRE),IPROPS(*)           ,
     2    RPROPS(*)          ,RSTAVA(MSTRE+1)   ,STRAN(MSTRE)
C Local arrays
      DIMENSION
     1    D12(3)             ,D21(3)
C***********************************************************************
C COMPUTATION OF THE CONSISTENT TANGENT MODULUS FOR THE DRUCKER-PRAGER
C ELASTO-PLASTIC MATERIAL WITH PIECE-WISE LINEAR ISOTROPIC HARDENING.
C PLANE STRESS IMPLEMENTATION ONLY.
C***********************************************************************
C Stops program if not plane stress
      IF(NTYPE.NE.1)CALL ERRPRT('EI0039')
C Retrieve the elastic trial THICKNESS STRAIN last determined in the
C plane stress enforcement loop of subroutine SUDPPN. The in-plane
C elastic trial components have already been stored in the first
C three components of STRAN before the present routine was called.
      STRAN(4)=RALGVA(3)
C Compute the axisymmetric consistent tangent matrix
      CALL CTDP
     1(   RALGVA     ,DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,
     2    3          ,RPROPS     ,RSTAVA     ,STRAN      )
C Decompose into submatrices
      D12(1)=DMATX(1,4)
      D12(2)=DMATX(2,4)
      D12(3)=DMATX(3,4)
      D21(1)=DMATX(4,1)
      D21(2)=DMATX(4,2)
      D21(3)=DMATX(4,3)
      D22=DMATX(4,4)
C Assemble plane stress consistent tangent matrix
      DO 20 I=1,3
        DO 10 J=1,3
          DMATX(I,J)=DMATX(I,J)-D12(I)*D21(J)/D22
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE CTDPPN