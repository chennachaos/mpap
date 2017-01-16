CDOC BEGIN_SUBROUTINE CTTR
CDOC Consistent tangent matrix for the Tresca model
CDOC
CDOC The tangent matrix computed in this routine is consistent with
CDOC fully implicit elastic predictor/return mapping algorithm for the
CDOC Tresca elasto-plastic material model with piece-wise linear
CDOC isotropic hardening carried out in subroutine \smparm{SUTR}.
CDOC This routine contains the plane strain and axisymmetric
CDOC implementations of the model.
CDOC It returns either the elastic tangent or the elasto-plastic
CDOC consistent tangent matrix depending on the input value of
CDOC the logical argument \smparm{EPFLAG}.
CDOC The present implementation is based on the use of closed formulae
CDOC for the derivatives of general isotropic
CDOC tensor functions of one tensor in which the eigenvalues of the
CDOC function are expressed as functions (the algorithmic incremental
CDOC principal stress-based constitutive relation in the present case)
CDOC of the eigenvalues of its tensor argument.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   >  Array of incremental plastic
CDOC C                          multipliers obtained in routine
CDOC C                          \smparm{SUTR}.
CDOC DOUBLE_PRECISION DMATX  <  Consistent tangent matrix.
CDOC LOGICAL          EPFLAG >  Elasto-plastic flag.
CDOC C                          If \smparm{.FALSE.}, \smparm{DMATX}
CDOC C                          returns as the elastic matrix. If
CDOC C                          \sparm{.TRUE.}, \smparm{DMATX} returns
CDOC C                          as the elasto-plastic tangent consistent
CDOC C                          with the return mapping algorithm
CDOC C                          implemented in routine \smparm{SUTR}.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDTR}.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC LOGICAL          LALGVA >  Array of logical algorithmic flags.
CDOC C                          See list of arguments of \smparm{SUTR}.
CDOC INTEGER          NTYPE  >  Stress state type.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          Same as in the argument list of
CDOC C                          subroutine \smparm{SUTR}.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other
CDOC C                           than the stress tensor components.
CDOC C                          Output of \smparm{SUTR}.
CDOC DOUBLE_PRECISION STRAN  >  Array of elastic trial (engineering)
CDOC C                          strain components. Same as in the input
CDOC C                          of \smparm{SUTR}.
CDOC DOUBLE_PRECISION STRES  >  Array of current stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, June 1996: Initial coding
CDOC
      SUBROUTINE CTTR
     1(   DGAM       ,DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,
     2    NTYPE      ,RPROPS     ,RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=4  ,MDIM=3,  MSTRE=4)
      LOGICAL EPFLAG, LALGVA(4), OUTOFP, RIGHT, REPEAT, TWOVEC
      DIMENSION
     1    DGAM(2)            ,DMATX(MSTRE,MSTRE) ,IPROPS(*)          ,
     2    RPROPS(*)          ,RSTAVA(MSTRE+1)    ,STRAN(*)           ,
     3    STRES(*)
      DIMENSION
     1    DPSTRS(MDIM,MDIM)  ,DPSTRE(MDIM,MDIM)  ,EIGPRJ(MSTRE,2)    ,
     2    FOID(MSTRE,MSTRE)  ,PSTRS(MDIM)        ,PSTRA(MDIM)        ,
     3    SOID(MSTRE)        ,STRAT(MSTRE)
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
     1    R0   ,RP5  ,R1   ,R2   ,R3   ,R4   / 
     2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,4.0D0/
C***********************************************************************
C COMPUTATION OF CONSISTENT TANGENT MODULUS FOR TRESCA TYPE
C ELASTO-PLASTIC MATERIAL WITH PIECE-WISE LINEAR ISOTROPIC HARDENING.
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C***********************************************************************
C Stops program if neither plane strain nor axisymmetric state
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0028')
C Current accumulated plastic strain
      EPBAR=RSTAVA(MSTRE+1)
C Set material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      NHARD=IPROPS(3)
C Set needed algorithmic variables
      DGAMA=DGAM(1)
      DGAMB=DGAM(2)
      TWOVEC=LALGVA(3)
      RIGHT=LALGVA(4)
C Set some constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R4G=R4*GMODU
      R1D3=R1/R3
      R2D3=R2*R1D3
      ROO4D3=SQRT(R4*R1D3)
      IF(EPFLAG)THEN
C Compute elastoplastic consistent tangent
C ----------------------------------------
C Spectral decomposition of the elastic trial strain
        STRAT(1)=STRAN(1)
        STRAT(2)=STRAN(2)
        STRAT(3)=STRAN(3)*RP5
        CALL SPDEC2(EIGPRJ,PSTRA,REPEAT,STRAT)
        PSTRA(3)=STRAN(4)
C and current total stress
        PSTRS(1)=STRES(1)*EIGPRJ(1,1)+STRES(2)*EIGPRJ(2,1)+
     1                             R2*STRES(3)*EIGPRJ(3,1)
        PSTRS(2)=STRES(1)*EIGPRJ(1,2)+STRES(2)*EIGPRJ(2,2)+
     1                             R2*STRES(3)*EIGPRJ(3,2)
        PSTRS(3)=STRES(4)
C Identify directions of maximum and minimum principal trial stresses
        II=1
        JJ=1
        PSTMAX=PSTRA(II)
        PSTMIN=PSTRA(JJ)
        DO 10 I=2,3
          IF(PSTRA(I).GT.PSTMAX)THEN
            II=I
            PSTMAX=PSTRA(II)
          ENDIF
          IF(PSTRA(I).LT.PSTMIN)THEN
            JJ=I
            PSTMIN=PSTRA(JJ)
          ENDIF
   10   CONTINUE
        IF(II.NE.1.AND.JJ.NE.1)MM=1
        IF(II.NE.2.AND.JJ.NE.2)MM=2
        IF(II.NE.3.AND.JJ.NE.3)MM=3
        IF(TWOVEC)THEN
C Tangent consistent with two-vector return algorithm
          IF(DGAMA.EQ.DGAMB)THEN
            DUNIDA=DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
            DUNIDB=DUNIDA
          ELSE
            DGABAR=SQRT(DGAMA*DGAMA+DGAMA*DGAMB+DGAMB*DGAMB)
            FACTOR=ROO4D3/(R2*DGABAR)*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
            DUNIDA=FACTOR*(R2*DGAMA+DGAMB)
            DUNIDB=FACTOR*(DGAMA+R2*DGAMB)
          ENDIF
          DAA=R4G+DUNIDA
          DAB=R2G+DUNIDB
          DBA=R2G+DUNIDA
          DBB=R4G+DUNIDB
          DET=DAA*DBB-DAB*DBA
          R2GDD=R2G/DET
          R4G2DD=R2G*R2GDD
          IF(RIGHT)THEN
C ...returned to right corner
            DPSTRS(II,II)=R2G*(R1-R2GDD*R4G)
            DPSTRS(II,MM)=R4G2DD*(DAA-DAB)
            DPSTRS(II,JJ)=R4G2DD*(DBB-DBA)
            DPSTRS(MM,II)=R4G2DD*R2G
            DPSTRS(MM,MM)=R2G*(R1-R2GDD*DAA)
            DPSTRS(MM,JJ)=R4G2DD*DBA
            DPSTRS(JJ,II)=R4G2DD*R2G
            DPSTRS(JJ,MM)=R4G2DD*DAB
            DPSTRS(JJ,JJ)=R2G*(R1-R2GDD*DBB)
          ELSE
C ...returned to left corner
            DPSTRS(II,II)=R2G*(R1-R2GDD*DBB)
            DPSTRS(II,MM)=R4G2DD*DAB
            DPSTRS(II,JJ)=R4G2DD*R2G
            DPSTRS(MM,II)=R4G2DD*DBA
            DPSTRS(MM,MM)=R2G*(R1-R2GDD*DAA)
            DPSTRS(MM,JJ)=R4G2DD*R2G
            DPSTRS(JJ,II)=R4G2DD*(DBB-DBA)
            DPSTRS(JJ,MM)=R4G2DD*(DAA-DAB)
            DPSTRS(JJ,JJ)=R2G*(R1-R2GDD*R4G)
          ENDIF
        ELSE
C Tangent consistent with one-vector return algorithm
          FACTOR=R2G/(R4G+ROO4D3*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD)))
          DPSTRS(II,II)=R2G*(R1-FACTOR)
          DPSTRS(II,MM)=R0
          DPSTRS(II,JJ)=R2G*FACTOR
          DPSTRS(MM,II)=DPSTRS(II,MM)
          DPSTRS(MM,MM)=R2G
          DPSTRS(MM,JJ)=R0
          DPSTRS(JJ,II)=DPSTRS(II,JJ)
          DPSTRS(JJ,MM)=DPSTRS(MM,JJ)
          DPSTRS(JJ,JJ)=DPSTRS(II,II)
        ENDIF
        DPSTRE(1,1)=+DPSTRS(1,1)*R2D3-DPSTRS(1,2)*R1D3-DPSTRS(1,3)*R1D3+
     1                                                             BULK
        DPSTRE(2,1)=+DPSTRS(2,1)*R2D3-DPSTRS(2,2)*R1D3-DPSTRS(2,3)*R1D3+
     1                                                             BULK
        DPSTRE(3,1)=+DPSTRS(3,1)*R2D3-DPSTRS(3,2)*R1D3-DPSTRS(3,3)*R1D3+
     1                                                             BULK
        DPSTRE(1,2)=-DPSTRS(1,1)*R1D3+DPSTRS(1,2)*R2D3-DPSTRS(1,3)*R1D3+
     1                                                             BULK
        DPSTRE(2,2)=-DPSTRS(2,1)*R1D3+DPSTRS(2,2)*R2D3-DPSTRS(2,3)*R1D3+
     1                                                             BULK
        DPSTRE(3,2)=-DPSTRS(3,1)*R1D3+DPSTRS(3,2)*R2D3-DPSTRS(3,3)*R1D3+
     1                                                             BULK
        DPSTRE(1,3)=-DPSTRS(1,1)*R1D3-DPSTRS(1,2)*R1D3+DPSTRS(1,3)*R2D3+
     1                                                             BULK
        DPSTRE(2,3)=-DPSTRS(2,1)*R1D3-DPSTRS(2,2)*R1D3+DPSTRS(2,3)*R2D3+
     1                                                             BULK
        DPSTRE(3,3)=-DPSTRS(3,1)*R1D3-DPSTRS(3,2)*R1D3+DPSTRS(3,3)*R2D3+
     1                                                             BULK
        IF(NTYPE.EQ.2)THEN
          OUTOFP=.FALSE.
        ELSEIF(NTYPE.EQ.3)THEN
          OUTOFP=.TRUE.
        ENDIF
        CALL DGISO2
     1(   DPSTRE     ,DMATX      ,EIGPRJ     ,PSTRA      ,PSTRS      ,
     2    OUTOFP     ,REPEAT     )
      ELSE
C Compute elasticity matrix
C -------------------------
        IF(NTYPE.EQ.2)THEN
          NSTRE=3
        ELSEIF(NTYPE.EQ.3)THEN
          NSTRE=4
        ENDIF
C
        FACTOR=BULK-R2G*R1D3
        DO 50 I=1,NSTRE
          DO 40 J=I,NSTRE
            DMATX(I,J)=R2G*FOID(I,J)+FACTOR*SOID(I)*SOID(J)
   40     CONTINUE
   50   CONTINUE
        DO 70 J=1,NSTRE-1
          DO 60 I=J+1,NSTRE
            DMATX(I,J)=DMATX(J,I)
   60     CONTINUE
   70   CONTINUE
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE CTTR