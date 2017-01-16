CDOC BEGIN_SUBROUTINE ARRGAX
CDOC Arrange a fourth order tensor in matrix form with G matrix ordering
CDOC
CDOC This routine re-arranges a given fourth order tensor, stored as a
CDOC 4-index array, in matrix form (2-index array) using G matrix
CDOC component ordering.
CDOC Implemented only for axisymmetric problems.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A4TH   >  Fourth order tensor stored as a 4-index
CDOC C                          array.
CDOC DOUBLE_PRECISION AMATX  <  2-index array containing the components
CDOC C                          of the given 4th order tensor stored
CDOC C                          using G matrix ordering.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, May  1999: Initial coding
CDOC
      SUBROUTINE ARRGAX
     1(   A4TH       ,AMATX      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=3     ,NGDIM=5    )
C Arguments
      DIMENSION
     1    A4TH(NDIM,NDIM,NDIM,NDIM)              ,AMATX(NGDIM,NGDIM)
C***********************************************************************
C RE-ARRANGES A FOURTH ORDER TENSOR, STORED AS A 4-INDEX ARRAY, IN
C MATRIX FORM (2-INDEX ARRAY) USING G MATRIX COMPONENT ORDERING
C C (11,21,12,22,33). FOR AXISYMMETRIC CASE ONLY.
C***********************************************************************
C
      AMATX(1,1)=A4TH(1,1,1,1) 
      AMATX(1,2)=A4TH(1,1,2,1) 
      AMATX(1,3)=A4TH(1,1,1,2)
      AMATX(1,4)=A4TH(1,1,2,2) 
      AMATX(1,5)=A4TH(1,1,3,3) 
C
      AMATX(2,1)=A4TH(2,1,1,1) 
      AMATX(2,2)=A4TH(2,1,2,1) 
      AMATX(2,3)=A4TH(2,1,1,2) 
      AMATX(2,4)=A4TH(2,1,2,2) 
      AMATX(2,5)=A4TH(2,1,3,3) 
C
      AMATX(3,1)=A4TH(1,2,1,1) 
      AMATX(3,2)=A4TH(1,2,2,1) 
      AMATX(3,3)=A4TH(1,2,1,2) 
      AMATX(3,4)=A4TH(1,2,2,2) 
      AMATX(3,5)=A4TH(1,2,3,3) 
C
      AMATX(4,1)=A4TH(2,2,1,1) 
      AMATX(4,2)=A4TH(2,2,2,1) 
      AMATX(4,3)=A4TH(2,2,1,2) 
      AMATX(4,4)=A4TH(2,2,2,2) 
      AMATX(4,5)=A4TH(2,2,3,3) 
C
      AMATX(5,1)=A4TH(3,3,1,1) 
      AMATX(5,2)=A4TH(3,3,2,1) 
      AMATX(5,3)=A4TH(3,3,1,2) 
      AMATX(5,4)=A4TH(3,3,2,2) 
      AMATX(5,5)=A4TH(3,3,3,3) 
C
      RETURN
      END
CDOC END_SUBROUTINE ARRGAX
CDOC BEGIN_SUBROUTINE ARRGO2
CDOC Arrange a fourth order tensor in matrix form with G matrix ordering
CDOC
CDOC This routine re-arranges a given fourth order tensor, stored as a
CDOC 4-index array, in matrix form (2-index array) using G matrix
CDOC component ordering.
CDOC Implemented only for 2-D tensors.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A4TH   >  Fourth order tensor stored as a 4-index
CDOC C                          array.
CDOC DOUBLE_PRECISION AMATX  <  2-index array containing the components
CDOC C                          of the given 4th order tensor stored
CDOC C                          using G matrix ordering.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, January  1999: Initial coding
CDOC
      SUBROUTINE ARRGO2
     1(   A4TH       ,AMATX      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=2     ,NGDIM=4    )
C Arguments
      DIMENSION
     1    A4TH(NDIM,NDIM,NDIM,NDIM)              ,AMATX(NGDIM,NGDIM)
C***********************************************************************
C RE-ARRANGES A FOURTH ORDER TENSOR, STORED AS A 4-INDEX ARRAY, IN
C MATRIX FORM (2-INDEX ARRAY) USING G MATRIX COMPONENT ORDERING
C C (11,21,12,22). FOR 2-D ONLY.
C***********************************************************************
C
      AMATX(1,1)=A4TH(1,1,1,1) 
      AMATX(1,2)=A4TH(1,1,2,1) 
      AMATX(1,3)=A4TH(1,1,1,2)
      AMATX(1,4)=A4TH(1,1,2,2) 
C
      AMATX(2,1)=A4TH(2,1,1,1) 
      AMATX(2,2)=A4TH(2,1,2,1) 
      AMATX(2,3)=A4TH(2,1,1,2) 
      AMATX(2,4)=A4TH(2,1,2,2) 
C
      AMATX(3,1)=A4TH(1,2,1,1) 
      AMATX(3,2)=A4TH(1,2,2,1) 
      AMATX(3,3)=A4TH(1,2,1,2) 
      AMATX(3,4)=A4TH(1,2,2,2) 
C
      AMATX(4,1)=A4TH(2,2,1,1) 
      AMATX(4,2)=A4TH(2,2,2,1) 
      AMATX(4,3)=A4TH(2,2,1,2) 
      AMATX(4,4)=A4TH(2,2,2,2) 
C
      RETURN
      END
CDOC END_SUBROUTINE ARRGO2
CDOC BEGIN_SUBROUTINE DEXPMP
CDOC Derivative of exponential map for general three-dimensional tensors
CDOC
CDOC This routine computes the derivative of the exponential of a
CDOC generally unsymmetric three-dimensional tensor.
CDOC It uses the series definition of the tensor exponential.
CDOC The exponential map itself is implemented in subroutine
CDOC \smparm{EXPMAP}.
CDOC 
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DEXPX  <  Derivative of the exponential map at X.
CDOC C                          This derivative is a fourth order
CDOC C                          tensor stored here as a 4-index array.
CDOC LOGICAL          NOCONV <  Logical convergence flag. Set to
CDOC C                          \smparm{.TRUE.} if the series fail to
CDOC C                          converge. Set to \smparm{.FALSE.}
CDOC C                          otherwise.
CDOC DOUBLE_PRECISION X      >  Tensor at which exponential derivative
CDOC C                          is to be computed.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, December 1998: Initial coding
CDOC
      SUBROUTINE DEXPMP
     1(   DEXPX      ,NOCONV     ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=3     ,NDIM2=9    ,NDIM4=81   ,MAXN=100   )
C Arguments
      LOGICAL  NOCONV
      DIMENSION
     1    DEXPX(NDIM,NDIM,NDIM,NDIM),X(NDIM,NDIM)
C Local arrays and variables
C...matrix of powers of X
      DIMENSION
     1    R1DFAC(MAXN)       ,XMATX(NDIM,NDIM,0:MAXN)
C...initialise identity matrix: X to the power 0
      DATA
     1    XMATX(1,1,0)  ,XMATX(1,2,0)  ,XMATX(1,3,0)  /
     2    1.D0          ,0.D0          ,0.D0          /
     3    XMATX(2,1,0)  ,XMATX(2,2,0)  ,XMATX(2,3,0)  /
     4    0.D0          ,1.D0          ,0.D0          /
     5    XMATX(3,1,0)  ,XMATX(3,2,0)  ,XMATX(3,3,0)  /
     6    0.D0          ,0.D0          ,1.D0          /
      DATA
     1    R0   ,RP5  ,R1   ,R2   ,TOL    ,OVER    ,UNDER   /
     2    0.0D0,0.5D0,1.0D0,2.0D0,1.0D-10,1.0D+100,1.0D-100/
C***********************************************************************
C COMPUTES THE DERIVATIVE OF THE EXPONENTIAL OF A (GENERALLY
C UNSYMMETRIC) 3-D TENSOR:
C     --  USES THE SERIES DEFINITION OF THE TENSOR EXPONENTIAL  --
C***********************************************************************
C Initialise convergence flag
      NOCONV=.FALSE.
C X to the power 1
      DO 20 I=1,NDIM
        DO 10 J=1,NDIM
          XMATX(I,J,1)=X(I,J)
   10   CONTINUE
   20 CONTINUE
C Zero remaining powers of X
      CALL RVZERO(XMATX(1,1,2),NDIM*NDIM*(MAXN-1))
C Compute X square
      DO 50 I=1,NDIM
        DO 40 J=1,NDIM
          DO 30 K=1,NDIM
            XMATX(I,J,2)=XMATX(I,J,2)+X(I,K)*X(K,J)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
C Compute principal invariants of X
      C1=X(1,1)+X(2,2)+X(3,3)
      C2=RP5*(C1*C1-(XMATX(1,1,2)+XMATX(2,2,2)+XMATX(3,3,2)))
      C3=X(1,1)*X(2,2)*X(3,3)+X(1,2)*X(2,3)*X(3,1)+
     1   X(1,3)*X(2,1)*X(3,2)-X(1,2)*X(2,1)*X(3,3)-
     2   X(1,1)*X(2,3)*X(3,2)-X(1,3)*X(2,2)*X(3,1)
C Compute X to the powers 3,4,...,NMAX using recursive formula
      R1DFAC(1)=R1
      R1DFAC(2)=RP5
      DO 80 N=3,MAXN
        R1DFAC(N)=R1DFAC(N-1)/DBLE(N)
        DO 70 I=1,NDIM
          DO 60 J=1,NDIM
            XMATX(I,J,N)=C1*XMATX(I,J,N-1)-C2*XMATX(I,J,N-2)+
     1                   C3*XMATX(I,J,N-3)
   60     CONTINUE
   70   CONTINUE
        XNNORM=SQRT(SCAPRD(XMATX(1,1,N),XMATX(1,1,N),NDIM2))
C...check number of terms required for series convergence
        IF(XNNORM.GT.OVER.OR.(XNNORM.LT.UNDER.AND.XNNORM.GT.R0)
     1                                  .OR.R1DFAC(N).LT.UNDER)THEN
C...numbers are to small or too big: Exit without computing derivative
          NOCONV=.TRUE.
          GOTO 999
        ELSEIF(XNNORM*R1DFAC(N).LT.TOL)THEN
C...series will converge with NMAX terms:
C   Carry on to derivative computation
          NMAX=N
          GOTO 90
        ENDIF
   80 CONTINUE
C...series will not converge for the currently prescribed tolerance
C   with the currently prescribed maximum number of terms MAXN:
C   Exit without computing derivative
      NOCONV=.TRUE.
      GOTO 999
   90 CONTINUE
C Compute the derivative of exponential map
      CALL RVZERO(DEXPX,NDIM4)
      DO 150 I=1,NDIM
        DO 140 J=1,NDIM
          DO 130 K=1,NDIM
            DO 120 L=1,NDIM
              DO 110 N=1,NMAX
                DO 100 M=1,N
                  DEXPX(I,J,K,L)=DEXPX(I,J,K,L)+
     1                           R1DFAC(N)*XMATX(I,K,M-1)*XMATX(L,J,N-M)
  100           CONTINUE
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
C
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE DEXPMP
CDOC BEGIN_SUBROUTINE EXPMAP
CDOC Exponential map for general three-dimensional tensors
CDOC
CDOC This routine computes the exponential of a generally unsymmetric
CDOC three-dimensional tensor. It uses the series definition of the
CDOC tensor exponential
CDOC 
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION EXPX   <  Tensor exponential of X.
CDOC LOGICAL          NOCONV <  Logical convergence flag. Set to
CDOC C                          \smparm{.TRUE.} if the series fail to
CDOC C                          converge. Set to \smparm{.FALSE.}
CDOC C                          otherwise.
CDOC DOUBLE_PRECISION X      >  Tensor whose exponential is to be
CDOC C                          computed.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, Sept 1998: Initial coding
CDOC
      SUBROUTINE EXPMAP
     1(   EXPX       ,NOCONV     ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=3     ,NDIM2=9    )
C Arguments
      LOGICAL  NOCONV
      DIMENSION
     1    EXPX(NDIM,NDIM)           ,X(NDIM,NDIM)
C Local arrays and variables
      DIMENSION
     1    XN(NDIM,NDIM)     ,XNM1(NDIM,NDIM)     ,XNM2(NDIM,NDIM)     ,
     2    XNM3(NDIM,NDIM)   ,X2(NDIM,NDIM)
      DATA
     1    R0   ,RP5  ,R1   ,R2   ,TOL    ,OVER    ,UNDER   /
     2    0.0D0,0.5D0,1.0D0,2.0D0,1.0D-10,1.0D+100,1.0D-100/
      DATA
     1    NMAX / 100 /
C***********************************************************************
C COMPUTES THE EXPONENTIAL OF A (GENERALLY UNSYMMETRIC) 3-D TENSOR:
C     --  USES THE SERIES DEFINITION OF THE TENSOR EXPONENTIAL  --
C
C Reference:
C C. Miehe, Exponential map algorithm for stress updates in
C anisotropic multiplicative elastoplasticity for single crystals,
C Int.J.Num.Meth.Engng., 39, pp 3367-3390, 1996.
C***********************************************************************
C Initialise series convergence flag
      NOCONV=.FALSE.
C Compute X square
      CALL RVZERO(X2,NDIM2)
      DO 30 I=1,NDIM
        DO 20 J=1,NDIM
          DO 10 K=1,NDIM
            X2(I,J)= X2(I,J)+X(I,K)*X(K,J)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C Compute principal invariants of X
      C1=X(1,1)+X(2,2)+X(3,3)
      C2=RP5*(C1*C1-(X2(1,1)+X2(2,2)+X2(3,3)))
      C3=X(1,1)*X(2,2)*X(3,3)+X(1,2)*X(2,3)*X(3,1)+
     1   X(1,3)*X(2,1)*X(3,2)-X(1,2)*X(2,1)*X(3,3)-
     2   X(1,1)*X(2,3)*X(3,2)-X(1,3)*X(2,2)*X(3,1)
C Start computation of exponential using its series definition
C ============================================================
      DO 50 I=1,NDIM
        DO 40 J=1,NDIM
          XNM1(I,J)=X2(I,J)
          XNM2(I,J)=X(I,J)
   40   CONTINUE
   50 CONTINUE
      XNM3(1,1)=R1
      XNM3(1,2)=R0
      XNM3(1,3)=R0
      XNM3(2,1)=R0
      XNM3(2,2)=R1
      XNM3(2,3)=R0
      XNM3(3,1)=R0
      XNM3(3,2)=R0
      XNM3(3,3)=R1
C Add first three terms of series
C -------------------------------
      DO 70 I=1,NDIM
        DO 60 J=1,NDIM
          EXPX(I,J)=RP5*XNM1(I,J)+XNM2(I,J)+XNM3(I,J)
   60   CONTINUE
   70 CONTINUE
C Add remaining terms (with X to the powers 3 to NMAX)
C ----------------------------------------------------
      FACTOR=R2
      DO 140 N=3,NMAX
C Use recursive formula to obtain X to the power N
        DO 90 I=1,NDIM
          DO 80 J=1,NDIM
            XN(I,J)=C1*XNM1(I,J)-C2*XNM2(I,J)+C3*XNM3(I,J)
   80     CONTINUE
   90   CONTINUE
C Update factorial
        FACTOR=DBLE(N)*FACTOR
        R1DFAC=R1/FACTOR
C Add Nth term of the series
        DO 110 I=1,NDIM
          DO 100 J=1,NDIM
            EXPX(I,J)=EXPX(I,J)+R1DFAC*XN(I,J)
  100     CONTINUE
  110   CONTINUE
C Check convergence of series
        XNNORM=SQRT(SCAPRD(XN(1,1),XN(1,1),NDIM2))
        IF(XNNORM.GT.OVER.OR.(XNNORM.LT.UNDER.AND.XNNORM.GT.R0)
     1                                     .OR.R1DFAC.LT.UNDER)THEN
C...first check possibility of overflow or underflow.
C...numbers are to small or too big: Break (unconverged) loop and exit
          NOCONV=.TRUE.
          GOTO 999
        ELSEIF(XNNORM*R1DFAC.LT.TOL)THEN
C...converged: Break series summation loop and exit with success
          GOTO 999
        ENDIF 
        DO 130 I=1,NDIM
          DO 120 J=1,NDIM
            XNM3(I,J)=XNM2(I,J)
            XNM2(I,J)=XNM1(I,J)
            XNM1(I,J)=XN(I,J)
  120     CONTINUE
  130   CONTINUE
  140 CONTINUE
C Re-set convergence flag if series did not converge
      NOCONV=.TRUE.
C
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE EXPMAP
      SUBROUTINE GAUSEL
     1(   A          ,B          ,N          ,SINGUL     )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C Arguments
      LOGICAL  SINGUL
      DIMENSION
     1    A(N,N)             ,B(N)
C Local definitions
      DATA
     1    SMALL  /
     2    1.0D-10/
C***********************************************************************
C Routine for solution of simultaneous linear algebraic equations:
C
C                                A x = b
C
C by the GAUSS ELIMINATION method WITH ROW SWAPPING.
C***********************************************************************
C Initialise singular matrix flag
      SINGUL=.FALSE.
C Get norm of matrix A
      ANORM=SQRT(SCAPRD(A,A,N*N))
C
C=============================
C Start loop over pivot rows |
C=============================
C
      DO 40 IPASS=1,N
C
        IF(ABS(A(IPASS,IPASS))/ANORM.LT.SMALL)THEN
C Current pivot is zero (too small): Perform row swapping
          CALL ROWSWP(A,ANORM,B,IPASS,N,SINGUL)
          IF(SINGUL)THEN
C Matrix is singular: exit without solving the system
            GOTO 999
          ENDIF
        ENDIF
C
C STEP 1.  Divide the entire pivot row by the pivot element to get a 1
C in the diagonal position of the pivot row
C
        PIVOT=A(IPASS,IPASS)
        DO 10 ICOL=1,N
          A(IPASS,ICOL)=A(IPASS,ICOL)/PIVOT
   10   CONTINUE
C... the same for the right hand side vector
        B(IPASS)=B(IPASS)/PIVOT
C
C STEP 2.  Replace each row other than the pivot row by that row plus a
C multiple of the pivot row to get a 0 in the pivot column
C
        DO 30 IROW=1,N
          IF(IROW.NE.IPASS)THEN
            FACTOR=A(IROW,IPASS)
            DO 20 ICOL=1,N
              A(IROW,ICOL)=A(IROW,ICOL)-FACTOR*A(IPASS,ICOL) 
   20       CONTINUE
C... the same for the right hand side vector
            B(IROW)=B(IROW)-FACTOR*B(IPASS)
          ENDIF
   30   CONTINUE
C
   40 CONTINUE
C
C==============================
C End of loop over pivot rows |
C==============================
C
C Now, "A" is the IDENTITY matrix and "B" is the solution
C vector.  Print solution vector.
C
  999 CONTINUE
      RETURN
      END




      SUBROUTINE ROWSWP
     1(   A          ,ANORM      ,B          ,I          ,N          ,
     2    SINGUL     )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C Arguments
      LOGICAL  SINGUL
      DIMENSION
     1    A(N,N)             ,B(N)
C Local definitions
      DATA
     1    SMALL  /
     2    1.0d-10/
C***********************************************************************
C If row "i" has a zero (sufficiently small) entry on the diagonal
C (column "i") then swap that row with the next row below row "i" with
C a non-zero (sufficiently large) entry in column "i"
C***********************************************************************
      DO 20 J=I,N
        IF(ABS(A(J,I))/ANORM.GE.SMALL)THEN
C Non-zero (sufficiently large) element found: Swap rows "i" and "j"
          DO 10 K=1,N
            ATMP=A(I,K)
            A(I,K)=A(J,K)
            A(J,K)=ATMP
   10     CONTINUE
          BTMP=B(I)
          B(I)=B(J)
          B(J)=BTMP
C Row swapping complete: break loop and exit successfully
          SINGUL=.FALSE.
          GOTO 30
        ENDIF
   20 CONTINUE
C Row swapping failed (no non-zero element found): matrix is singular
      SINGUL=.TRUE.
C
   30 CONTINUE
      RETURN
      END
C  @(#)   Module:<lubksb.f>   Version:1.4   Date:05/03/94
      SUBROUTINE LUBKSB
     1(   A          ,B          ,INDX       ,
     2    N          ,NP         )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1    A(NP,NP)   ,INDX(N)    ,B(N)
      DATA NAME/'LUBKSB'/
C***********************************************************************
C Routine to solve the set of N linear equations AX=B
C See NUMERICAL RECIPES p36.
C*ACRONYM
C LU_BacK_SuBustitutions
C*DESCRIPTION
C*HISTORY
C Name          Date         Comment
C G.C.Huang    Oct,92      initial coding
C*EXTERNAL
C Arrays
C A      - LU decomposed matrix
C=B      - Right hand side matrix as input and stored solutions as output
C INDX   - Permutation vector
C Variables
C N      - Size of the problem
C NP     - Physical size of A matrix
C (c) Copyright 1992, Rockfield Software Limited, Swansea, UK
C***********************************************************************
cccccccccccccccD     CALL SENTRY(NAME,MODEDB)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
   11     CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
   12 CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
   13     CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
   14 CONTINUE
cccccccccccccccD     CALL SEXIT(MODEDB)
      RETURN
      END
C  @(#)   Module:<ludcmp.f>   Version:1.6   Date:10/02/95
      SUBROUTINE LUDCMP
     1(   A          ,INDX       ,
     2    D          ,N          ,NP         ,ERROR      )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      PARAMETER
     1(   NMAX=50    )
      LOGICAL
     1    ERROR
      DIMENSION
     1    A(NP,NP)   ,INDX(N)    ,VV(NMAX)
C$DP,3
      DATA
     1    R0         ,R1         ,TINY       /
     2    0.0D0      ,1.0D0      ,1.0D-19    /
C$SP,3
C     DATA
C    1    R0         ,R1         ,TINY       /
C    2    0.0        ,1.0        ,1.0E-19    /
      DATA NAME/'LUDCMP'/
C***********************************************************************
C Routine to do LU decomposition
C See NUMERICAL RECIPES p35.
C*ACRONYM
C LU_DeCOmPosition
C*DESCRIPTION
C*HISTORY
C Name          Date         Comment
C G.C.Huang    Oct,92      initial coding
C*EXTERNAL
C Arrays
C=A      - LU decomposed matrix
C=INDX   - Permutation vector
C Variables
C=D      - Row interchange indicator
C N      - Size of the problem
C NP     - Physical size of A matrix
C ERROR  - Error flag
C (c) Copyright 1992, Rockfield Software Limited, Swansea, UK
C***********************************************************************
ccccccccccccccD     CALL SENTRY(NAME,MODEDB)
      ERROR=.FALSE.
      D=R1
      DO 12 I=1,N
        AAMAX=R0
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
   11   CONTINUE
        IF(AAMAX.EQ.R0)THEN
C Error. Singular matrix encountered.
cccccccccccccc          CALL WRTER('A0176E',NAME,0 )
C Flag the error
          ERROR=.TRUE.
          GOTO 999
        ENDIF
        VV(I)=R1/AAMAX
   12 CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                IF(ABS(A(I,K)).LT.TINY.AND.ABS(A(K,J)).LT.TINY)GOTO 13
C If A(I,K) and A(K,J) are very small, skip the following line. Otherwise
C we have problem here. (G.C.)
                SUM=SUM-A(I,K)*A(K,J)
   13         CONTINUE
              A(I,J)=SUM
            ENDIF
   14     CONTINUE
        ENDIF
        AAMAX=R0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              IF(ABS(A(I,K)).LT.TINY.AND.ABS(A(K,J)).LT.TINY)GOTO 15
C If A(I,K) and A(K,J) are very small, skip the following line. Otherwise
C we have problem here. (G.C.)
              SUM=SUM-A(I,K)*A(K,J)
   15       CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
   16   CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
   17     CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.TINY)A(J,J)=TINY
          DUM=R1/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
   18     CONTINUE
        ENDIF
   19 CONTINUE
      IF(A(N,N).EQ.R0)A(N,N)=TINY
 999  CONTINUE
ccccccccccccccD     CALL SEXIT(MODEDB)
      RETURN
      END




CDOC BEGIN_SUBROUTINE PODEC2
CDOC Polar decomposition of 2-D tensors
CDOC
CDOC This routine performs the right polar decomposition of 2-D
CDOC tensors:   F = R U, where R is a rotation (orthogonal tensor)
CDOC and U is a symmetric tensor.
CDOC 
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION F      >  2-D tensor to be decomposed.
CDOC C                          Dimension 2x2.
CDOC DOUBLE_PRECISION R      <  Rotation matrix resulting from the polar
CDOC C                          decomposition. Dimension 2x2.
CDOC DOUBLE_PRECISION U      <  Right symmetric tensor resulting from
CDOC C                          the polar decomposition. Dimension 2x2.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, November 1998: Initial coding
CDOC
      SUBROUTINE PODEC2
     1(   F          ,R          ,U          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=2     )
C Arguments
      DIMENSION
     1    F(NDIM,NDIM)       ,R(NDIM,NDIM)       ,U(NDIM,NDIM)
C Local variables and arrays
      LOGICAL DUMMY
      DIMENSION
     1    C(NDIM,NDIM)       ,CVEC(4)            ,EIGPRJ(4,NDIM)     ,
     2    EIGC(NDIM)         ,UM1(NDIM,NDIM)     ,UM1VEC(3)          ,
     3    UVEC(3)
      DATA
     1    R1   /
     2    1.0D0/
C***********************************************************************
C PERFORMS THE RIGHT POLAR DECOMPOSITION OF A 2-D TENSOR
C***********************************************************************
C                T
C Compute  C := F  F
C -------------------
      CALL RVZERO(C,NDIM*NDIM)
      DO 30 I=1,NDIM
        DO 20 J=1,NDIM
          DO 10 K=1,NDIM
            C(I,J)=C(I,J)+F(K,I)*F(K,J)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C Perform spectral decomposition of C
C -----------------------------------
      CVEC(1)=C(1,1)
      CVEC(2)=C(2,2)
      CVEC(3)=C(1,2)
      CALL SPDEC2(EIGPRJ     ,EIGC       ,DUMMY      ,CVEC       )
C
C                 1/2         -1
C Compute  U := (C)    and   U
C ------------------------------
C assemble in vector form
      CALL RVZERO(UVEC,3)
      CALL RVZERO(UM1VEC,3)
      DO 50 IDIM=1,NDIM
        UEIG=SQRT(EIGC(IDIM))
        UM1EIG=R1/UEIG
        DO 40 ICOMP=1,3
          UVEC(ICOMP)=UVEC(ICOMP)+UEIG*EIGPRJ(ICOMP,IDIM)
          UM1VEC(ICOMP)=UM1VEC(ICOMP)+UM1EIG*EIGPRJ(ICOMP,IDIM)
   40   CONTINUE
   50 CONTINUE
C and matrix form
      U(1,1)=UVEC(1)
      U(2,2)=UVEC(2)
      U(1,2)=UVEC(3)
      U(2,1)=UVEC(3)
      UM1(1,1)=UM1VEC(1)
      UM1(2,2)=UM1VEC(2)
      UM1(1,2)=UM1VEC(3)
      UM1(2,1)=UM1VEC(3)
C                           -1
C Compute rotation  R := F U
C ----------------------------
      CALL RVZERO(R,NDIM*NDIM)
      DO 80 I=1,NDIM
        DO 70 J=1,NDIM
          DO 60 K=1,NDIM
            R(I,J)=R(I,J)+F(I,K)*UM1(K,J)
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE PODEC2
