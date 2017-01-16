
CDOC BEGIN_SUBROUTINE BETRIA
CDOC Computes the left Cauchy-Green strain tensor
CDOC
CDOC Given the previous elastic (total) left Cauchy-Green strain tensor
CDOC and the incremental deformation gradient between the previous and
CDOC current configuration, this routine computes the current elastic
CDOC trial (total) left Cauchy-Green strain tensor. This routine
CDOC contains the plane strain, plane stress and axisymmetric
CDOC implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION BEN    >  Array of components of the previous
CDOC C                          elastic (or total) left Cauchy-Green
CDOC C                          strain tensor.
CDOC DOUBLE_PRECISION BETRL  <  Array of components of the current
CDOC C                          elastic trial (or total) left
CDOC C                          Cauchy-Green strain tensor.
CDOC DOUBLE_PRECISION FINCR  >  Incremental deformation gradient between
CDOC C                          the previous and current configuration.
CDOC INTEGER          NTYPE  >  Stress state type. Present routine is
CDOC C                          compatible with \smparm{NTYPE=1} (plane
CDOC C                          stress), \smparm{NTYPE=2} (plane strain)
CDOC C                          and \smparm{NTYPE=3} (axisymmetric
CDOC C                          condition).
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      SUBROUTINE BETRIA
     1(   BEN        ,BETRL      ,FINCR      ,NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    BEN(*)             ,BETRL(*)           ,FINCR(3,3)
      DIMENSION
     1    AUXM(2,2)          ,BENMTX(2,2)        ,BETRLM(2,2)
C***********************************************************************
C COMPUTES THE ("ELASTIC TRIAL" FOR HYPERELASTIC-BASED LARGE STRAIN
C ELASTO-PLASTIC MODELS) LEFT CAUCHY-GREEN STRAIN TENSOR:
C
C                     e trial             e   T
C                    B        :=   F     B   F
C                     n+1           incr  n   incr
C
C***********************************************************************
C Convert previously converged elastic left Cauchy-Green strain tensor
C from vector form to matrix form
      BENMTX(1,1)=BEN(1)
      BENMTX(2,1)=BEN(3)
      BENMTX(1,2)=BEN(3)
      BENMTX(2,2)=BEN(2)


C
C In-plane components of the elastic trial left Cauchy-Green tensor
C
      CALL RVZERO(AUXM,4)
      DO 30 I=1,2
        DO 20 J=1,2
          DO 10 K=1,2
            AUXM(I,J)=AUXM(I,J)+FINCR(I,K)*BENMTX(K,J)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      CALL RVZERO(BETRLM,4)
      DO 60 I=1,2
        DO 50 J=1,2
          DO 40 K=1,2
            BETRLM(I,J)=BETRLM(I,J)+AUXM(I,K)*FINCR(J,K)
   40     CONTINUE
   50   CONTINUE
   60 CONTINUE
C
C        e trial
C Store B        in vector form
C        n+1 
C
      BETRL(1)=BETRLM(1,1)
      BETRL(2)=BETRLM(2,2)
      BETRL(3)=BETRLM(1,2)
C out-of-plane component
      IF(NTYPE.EQ.2)THEN
        BETRL(4)=BEN(4)
      ELSEIF(NTYPE.EQ.3)THEN
        BETRL(4)=BEN(4)*FINCR(3,3)*FINCR(3,3)
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE BETRIA
CDOC BEGIN_SUBROUTINE CSTEP2
CDOC Consistent spatial tangent modulus for finite elasto-plastic models
CDOC
CDOC This routine computes the spatial tangent modulus, a, for
CDOC hyperelastic based (logarithmic strain-based) finite
CDOC elasto-plasticity models in 2-D: Plane strain, plane stress and
CDOC axisymmetric states.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AMATX  <  Matrix of components of the spatial
CDOC C                          tangent modulus, a.
CDOC DOUBLE_PRECISION BETRL  >  Array of components of the elastic trial
CDOC C                          left Cauchy-Green strain tensor.
CDOC DOUBLE_PRECISION DMATX  >  Array of components of the infinitesimal
CDOC C                          consistent tangent modulus.
CDOC DOUBLE_PRECISION STRES  >  Array of Cauchy stress tensor
CDOC C                          components
CDOC DOUBLE_PRECISION DETF   >  Determinant of the current deformation
CDOC C                          gradient.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      SUBROUTINE CSTEP2
     1(   AMATX      ,BETRL      ,DMATX      ,STRES      ,DETF       ,
     2    NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MADIM=5    ,MSTRE=4    )
      EXTERNAL  DDLGD2 ,DLGD2
      LOGICAL   OUTOFP
      DIMENSION
     1    AMATX(MADIM,MADIM)        ,BETRL(MSTRE)              ,
     2    DMATX(MSTRE,MSTRE)        ,STRES(MSTRE)
      DIMENSION
     1    AUXMTX(MADIM,MADIM)       ,BMTX(MADIM,MADIM)         ,
     2    DLGAUX(MSTRE,MSTRE)       ,DLGMTX(MADIM,MADIM)       ,
     3    DMATX2(MADIM,MADIM)       ,IG(MADIM)
      DATA
     1    R1   ,R2    /
     2    1.0D0,2.0D0 /
      DATA
     1    IG(1),IG(2),IG(3),IG(4),IG(5)  /
     2    1    ,3    ,3    ,2    ,4      /
C***********************************************************************
C COMPUTE THE CONSISTENT SPATIAL TANGENT MODULUS 'a' FOR LARGE STRAIN
C HYPERELASTIC-BASED ELASTOPLASTIC MATERIAL MODELS
C***********************************************************************
      IF(NTYPE.EQ.3)THEN
        OUTOFP=.TRUE.
        NADIM=5
      ELSEIF(NTYPE.EQ.1.OR.NTYPE.EQ.2)THEN
        OUTOFP=.FALSE.
        NADIM=4
      ELSE
        CALL ERRPRT('EI0020')
      ENDIF

C                           e trial   e trial
C Compute the derivative  dE       /dB
      CALL DISO2
     1(   DLGAUX     ,DDLGD2     ,DLGD2      ,OUTOFP     ,BETRL      )
      FACTOR=R1/DETF
      DO 20 I=1,MSTRE
        DO 10 J=1,MSTRE
          DLGAUX(I,J)=FACTOR*DLGAUX(I,J)
   10   CONTINUE
   20 CONTINUE
C                                      e trial   e trial
C Rearrange components of DMATX and  dE       /dB         into the
C ordering (11,21,12,22,33), compatible with the discrete gradient G.
      CALL RVZERO(DMATX2,MADIM*MADIM)
      CALL RVZERO(DLGMTX,MADIM*MADIM)
      DO 40 INEW=1,NADIM
        DO 30 JNEW=1,NADIM
          IOLD=IG(INEW)
          JOLD=IG(JNEW)
          DLGMTX(INEW,JNEW)=DLGAUX(IOLD,JOLD)
          DMATX2(INEW,JNEW)=DMATX(IOLD,JOLD)
   30   CONTINUE
   40 CONTINUE



C Compute remaining needed matrix [DELTA_ik BETRL_jl+DELTA_jl BETRL_il]
      CALL RVZERO(BMTX,MADIM*MADIM)
      BMTX(1,1)=R2*BETRL(1)
      BMTX(1,3)=R2*BETRL(3)
      BMTX(2,1)=BETRL(3)
      BMTX(2,2)=BETRL(1)
      BMTX(2,3)=BETRL(2)
      BMTX(2,4)=BETRL(3)
      BMTX(3,1)=BETRL(3)
      BMTX(3,2)=BETRL(1)
      BMTX(3,3)=BETRL(2)
      BMTX(3,4)=BETRL(3)
      BMTX(4,2)=R2*BETRL(3)
      BMTX(4,4)=R2*BETRL(2)
      IF(OUTOFP)BMTX(5,5)=R2*BETRL(4)


C Assemble the spatial tangent modulus a
C --------------------------------------
C compute the product  D:L:B
      CALL RVZERO(AUXMTX,MADIM*MADIM)
      DO 70 I=1,NADIM
        DO 60 J=1,NADIM
          DO 50 K=1,NADIM
            AUXMTX(I,J)=AUXMTX(I,J)+DMATX2(I,K)*DLGMTX(K,J)
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
      CALL RVZERO(AMATX,MADIM*MADIM)



      DO 100 I=1,NADIM
        DO 90 J=1,NADIM
          DO 80 K=1,NADIM
            AMATX(I,J)=AMATX(I,J)+AUXMTX(I,K)*BMTX(K,J)
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE



C subtract  [SIGMA_il DELTA_jk]
      AMATX(1,1)=AMATX(1,1)-STRES(1)
      AMATX(1,3)=AMATX(1,3)-STRES(3)
      AMATX(2,1)=AMATX(2,1)-STRES(3)
      AMATX(2,3)=AMATX(2,3)-STRES(2)
      AMATX(3,2)=AMATX(3,2)-STRES(1)
      AMATX(3,4)=AMATX(3,4)-STRES(3)
      AMATX(4,2)=AMATX(4,2)-STRES(3)
      AMATX(4,4)=AMATX(4,4)-STRES(2)
      IF(OUTOFP)AMATX(5,5)=AMATX(5,5)-STRES(4)


      RETURN
      END
CDOC END_SUBROUTINE CSTEP2
CDOC BEGIN_SUBROUTINE DEFGRA
CDOC Deformation gradient for 2D isoparametric finite element
CDOC
CDOC Given the element nodal displacements and the discrete gradient
CDOC operator, G-matrix, at a point, this routine computes the
CDOC corresponding deformation gradient at that point. This routine
CDOC contains the plane strain, plane stress and axisymmetric
CDOC implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION ELDISP >  Array of nodal displacements of the
CDOC C                          finite element.
CDOC DOUBLE_PRECISION F      <  Deformation gradient.
CDOC DOUBLE_PRECISION GMATX  >  Discrete (full) gradient operator,
CDOC C                          G-matrix, at the point of interest.
CDOC INTEGER          MDOFN  >  Dimensioning parameter: number of
CDOC C                          rows of array \smparm{ELDISP}.
CDOC INTEGER          MGDIM  >  Dimensioning parameter: number of
CDOC C                          rows of array \smparm{GMATX}.
CDOC INTEGER          NDOFN  >  Number of degrees of freedom per node.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      SUBROUTINE DEFGRA
     1(   ELDISP     ,F          ,GMATX      ,MDOFN      ,MGDIM      ,
     2    NDOFN      ,NTYPE      ,NNODE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=5    )
      DIMENSION
     1    ELDISP(MDOFN,*)    ,F(3,3)             ,GMATX(MGDIM,*)
      DIMENSION
     1    FVEC(MCOMP)
      DATA  R0   ,R1    / 0.0D0,1.0D0 /
C***********************************************************************
C COMPUTES THE DEFORMATION GRADIENT TENSOR ASSOCIATED WITH THE ELEMENT
C DISPLACEMENT 'ELDISP'
C***********************************************************************
C Set total number of deformation gradient components
      IF(NTYPE.EQ.1.OR.NTYPE.EQ.2)THEN
        NCOMP=4
      ELSEIF(NTYPE.EQ.3)THEN
        NCOMP=5
      ELSE
        CALL ERRPRT('EI0021')
      ENDIF
C Evaluate the deformation gradient stored in vector form
      CALL RVZERO(FVEC,NCOMP)
      DO 30 ICOMP=1,NCOMP
        IEVAB=0
        DO 20 INODE=1,NNODE
          DO 10 IDOFN=1,NDOFN
            IEVAB=IEVAB+1
            FVEC(ICOMP)=FVEC(ICOMP)+
     1                  GMATX(ICOMP,IEVAB)*ELDISP(IDOFN,INODE)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C Store the deformation gradient in matrix form
      F(1,1)=FVEC(1)+R1
      F(2,1)=FVEC(2)
      F(3,1)=R0
      F(1,2)=FVEC(3)
      F(2,2)=FVEC(4)+R1
      F(3,2)=R0
      F(1,3)=R0
      F(2,3)=R0
      IF(NTYPE.EQ.1)THEN
        F(3,3)=R0
      ELSEIF(NTYPE.EQ.2)THEN
        F(3,3)=R1
      ELSEIF(NTYPE.EQ.3)THEN
        F(3,3)=FVEC(5)+R1
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE DEFGRA
CDOC BEGIN_SUBROUTINE ERRPRT
CDOC Prints error message and can abort the program if requested
CDOC
CDOC This routine prints error/warning messages to the standard output
CDOC and results file and may abort the program depending on the entry
CDOC value of its argument. The character string passed as its argument
CDOC is an error code which will be searched for in the file ERROR.RUN,
CDOC assumed to be kept in the directory defined by the HYPLASHOME
CDOC environment variable. If the correponding character string (error
CDOC code) is found, the associated error/warning message is printed
CDOC in the standard output and results file.
CDOC The are 5 types of errors/warnings: Input data error, input data
CDOC warning, internal error, execution error and execution warning.
CDOC These are characterised, respectively, by error codes of the
CDOC types: 'ED????', 'WD????', 'EI????', 'EE????' and 'WE????'.
CDOC This routine aborts the program in case of input data error
CDOC (ED????), internal error ('EI????') or execution error (EE????).
CDOC The execution of the program is not interrupted in case of
CDOC warnings ('WE????' or 'WD????').
CDOC See file ERROR.RUN for more details.
CDOC 
CDOC BEGIN_PARAMETERS
CDOC CHARACTEER       ERRCOD <  Character string containing the error
CDOC C                          code.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, July 1996: Initial coding
CDOC
      SUBROUTINE ERRPRT(ERRCOD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL       AVAIL ,FOUND
      CHARACTER*6   ERRCOD
      CHARACTER*22  HEADER
      CHARACTER*72  INLINE
      CHARACTER*256 HYPLASHOME
      DIMENSION IWBEG(40), IWEND(40)
C***********************************************************************
C PRINT ERROR MESSAGE
C***********************************************************************
 1000 FORMAT(///' ',74('*')/' *',72X,'*'/' *',25X,A22,25X,'*'/' *',72X,
     1      '*'/' ',74('*')/' *',72X,'*'/' * Code:         ',A6,51X,'*')
 1010 FORMAT(' *',A72,'*')
 1020 FORMAT(' *',72X,'*'/' ',74('*'))
 1030 FORMAT(//'   ASSOCIATED MESSAGES WILL NOT BE PRINTED !',//,
     1  '   The above error code was not found in file ERROR.RUN.',//)
 1040 FORMAT(//'   ASSOCIATED MESSAGES WILL NOT BE PRINTED !',//,
     1   '   File ERROR.RUN does not exist in HYPLASHOME directory.',//)
C
      IF(ERRCOD(1:2).EQ.'ED')THEN
        HEADER='  INPUT  DATA  ERROR  '
      ELSEIF(ERRCOD(1:2).EQ.'WD')THEN
        HEADER=' INPUT  DATA  WARNING '
      ELSEIF(ERRCOD(1:2).EQ.'EI')THEN
        HEADER='   INTERNAL   ERROR   '
      ELSEIF(ERRCOD(1:2).EQ.'EE')THEN
        HEADER='    EXECUTION ERROR   '
      ELSEIF(ERRCOD(1:2).EQ.'WE')THEN
        HEADER='   EXECUTION WARNING  '
      ELSE
        HEADER='  UNKNOWN ERROR TYPE  '
      ENDIF
      WRITE(*,1000)HEADER,ERRCOD
      WRITE(16,1000)HEADER,ERRCOD
C
CC      CALL GETENV('HYPLASHOME',HYPLASHOME)
CC      LENGTH=INDEX(HYPLASHOME,' ')-1
CC      INQUIRE(FILE=HYPLASHOME(1:LENGTH)//'/ERROR.RUN',EXIST=AVAIL) 
CC      IF(.NOT.AVAIL)THEN
CC        WRITE(*,1040)
CC        WRITE(16,1040)
CC        GOTO 999
CC      ENDIF
CC      OPEN(23,FILE=HYPLASHOME(1:LENGTH)//'/ERROR.RUN',STATUS='OLD')
C      OPEN(23,FILE='ERROR.RUN',STATUS='OLD')
C
      CALL FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,ERRCOD     ,INLINE     ,
     2    23         ,NWRD       )
      IF(.NOT.FOUND)THEN
        WRITE(*,1030)
        WRITE(16,1030)
        GOTO 998
      ENDIF
      NLINES=INTNUM(INLINE(IWBEG(2):IWEND(2)))
      DO 10 I=1,NLINES
        READ(23,'(A72)')INLINE
        WRITE(*,1010)INLINE
        WRITE(16,1010)INLINE
   10 CONTINUE
      WRITE(*,1020)
      WRITE(16,1020)
  998 CONTINUE
      CLOSE(23,STATUS='KEEP')
  999 CONTINUE
C Aborts the program if not execution warning
      IF(ERRCOD(1:2).NE.'WE'.AND.ERRCOD(1:2).NE.'WD')CALL PEXIT
      RETURN
      END
CDOC END_SUBROUTINE ERRPRT
CDOC BEGIN_SUBROUTINE FCLOSE
CDOC Closes data file and results file
CDOC
      SUBROUTINE FCLOSE
C***********************************************************************
C CLOSES DATA AND RESULTS FILES
C***********************************************************************
      CLOSE(UNIT=15,STATUS='KEEP')
      CLOSE(UNIT=16,STATUS='KEEP')
      RETURN
      END
CDOC END_SUBROUTINE FCLOSE
CDOC BEGIN_SUBROUTINE GETBMX
CDOC Computes the discrete symmetric gradient operator for 2-D elements
CDOC 
CDOC This routine assembles the discrete symmetric gradient operator
CDOC (strain-displacement matrix in small strains), the B-matrix, for
CDOC isoparametric 2-D finite elements: Plane strain, plane stress and
CDOC axisymmetric cases.
CDOC 
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION BMATX  <  The discrete symmetric gradient
CDOC C                          operator, B-matrix.
CDOC DOUBLE_PRECISION CARTCO >  Cartesian coordinates of the point where
CDOC C                          the B-matrix is to be computed.
CDOC DOUBLE_PRECISION CARTD  >  Array of cartesian derivatives of the
CDOC C                          element shape functions at the point of
CDOC C                          interest.
CDOC INTEGER          NDIME  >  Dimensioning parameter: Number of rows
CDOC C                          of \smparm{CARTCO} and \smparm{CARTD}.
CDOC INTEGER          MBDIM  >  Dimensioning parameter: Number of rows
CDOC C                          of \smparm{BMATX}.
CDOC INTEGER          NAXIS  >  Axis of symmetry flag. Used only for the
CDOC C                          axisymmetric case.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION SHAPE  >  Array containing the value of the shape
CDOC C                          functions at the point of interest.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC 
      SUBROUTINE GETBMX
     1(   BMATX      ,CARTCO     ,CARTD      ,NDIME      ,MBDIM      ,  
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    BMATX(MBDIM,*)     ,CARTCO(NDIME)      ,CARTD(NDIME,*)     ,
     2    SHAPE(*)
      DATA R0/0.0D0/
C***********************************************************************
C EVALUATES THE DISCRETE SYMMETRIC GRADIENT OPERATOR 'B' (SMALL
C STRAIN-DISPLACEMENT MATRIX) FOR PLANE STRESS/STRAIN AND AXISYMMETRIC
C PROBLEMS
C***********************************************************************
C Plane strain/stress
C -------------------
      IY=0
      DO 10 INODE=1,NNODE
        IX=IY+1
        IY=IX+1
        BMATX(1,IX)=CARTD(1,INODE)
        BMATX(1,IY)=R0
        BMATX(2,IX)=R0
        BMATX(2,IY)=CARTD(2,INODE)
        BMATX(3,IX)=CARTD(2,INODE)
        BMATX(3,IY)=CARTD(1,INODE)
   10 CONTINUE
      IF(NTYPE.EQ.3)THEN
C Axisymmetric problem
C --------------------
        IY=0
        DO 20 INODE=1,NNODE
          IX=IY+1
          IY=IX+1
          IF(NAXIS.EQ.1)THEN
C Axisymmetric about Y axis
            BMATX(4,IX)=SHAPE(INODE)/CARTCO(NAXIS)
            BMATX(4,IY)=R0
          ELSE IF(NAXIS.EQ.2)THEN
C Axisymmetric about X axis
            BMATX(4,IX)=R0
            BMATX(4,IY)=SHAPE(INODE)/CARTCO(NAXIS)
          ENDIF
   20   CONTINUE
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE GETBMX
CDOC BEGIN_SUBROUTINE GETGCO
CDOC Gets coordinates of a point within an element by interpolation
CDOC
CDOC This routine computes the global cartesian coordinates of a point
CDOC within a finite element by interpolation of its nodal coordinates.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION CARTCO <  Cartesian coordinates of the point of
CDOC C                          interest.
CDOC DOUBLE_PRECISION ELCOD  >  Array of nodal coordinates of the
CDOC C                          element.
CDOC INTEGER          MDIME  >  Dimensioning parameter: Number of rows
CDOC C                          of \smparm{ELCOD}.
CDOC INTEGER          NDIME  >  Number of spatial dimensions.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC DOUBLE_PRECISION SHAPE  >  Array containing the value of the shape
CDOC C                          function at the point of interest.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto  August  1996
CDOC
      SUBROUTINE GETGCO
     1(   CARTCO     ,ELCOD      ,MDIME      ,NDIME      ,NNODE      ,
     2    SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    CARTCO(NDIME)      ,ELCOD(MDIME,NNODE) ,SHAPE(NNODE)
      DATA R0/0.0D0/
C***********************************************************************
C EVALUATES THE GLOBAL CARTESIAN COORDINATES OF A POINT WITHIN AN
C ELEMENT BY INTERPOLATION OF THE ELEMENT NODAL COORDINATES
C***********************************************************************
      DO 20 IDIME=1,NDIME
        CARTCO(IDIME)=R0
        DO 10 INODE=1,NNODE
          CARTCO(IDIME)=CARTCO(IDIME)+ELCOD(IDIME,INODE)*SHAPE(INODE)
   10   CONTINUE
   20 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE GETGCO
CDOC BEGIN_SUBROUTINE GETGMX
CDOC Computes the discrete (full) gradient operator for 2-D elements
CDOC
CDOC This routine assembles the discrete gradient operator, the G-matrix
CDOC for isoparametric 2-D finite elements: Plane strain, plane stress
CDOC and axisymmetric cases.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION CARTCO >  Cartesian coordinates of the point where
CDOC C                          the G-matrix is to be computed.
CDOC DOUBLE_PRECISION CARTD  >  Array of cartesian derivatives of the
CDOC C                          element shape functions at the point of
CDOC C                          interest.
CDOC DOUBLE_PRECISION GMATX  <  The discrete gradient operator,
CDOC C                          G-matrix.
CDOC INTEGER          MDIME  >  Dimensioning parameter: Number of rows
CDOC C                          of \smparm{CARTD}.
CDOC INTEGER          MGDIM  >  Dimensioning parameter: Number of rows
CDOC C                          of \smparm{GMATX}.
CDOC INTEGER          NAXIS  >  Axis of symmetry flag. Used only for the
CDOC C                          axisymmetric case.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION SHAPE  >  Array containing the value of the shape
CDOC C                          functions at the point of interest.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto  August  1996
CDOC
      SUBROUTINE GETGMX
     1(   CARTCO     ,CARTD      ,GMATX      ,MDIME      ,MGDIM      ,
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    CARTCO(MDIME)      ,CARTD(MDIME,*)     ,GMATX(MGDIM,*)     ,
     2    SHAPE(*)
      DATA R0/0.0D0/
C***********************************************************************
C EVALUATES THE DISCRETE (FULL) GRADIENT OPERATOR 'G' FOR PLANE
C STRESS/STRAIN AND AXISYMMETRIC PROBLEMS. COMPONENT ORDERING
C (11,21,12,22,33).
C***********************************************************************
C Plane strain/stress
C -------------------
      IY=0
      DO 10 INODE=1,NNODE
        IX=IY+1
        IY=IX+1
        GMATX(1,IX)=CARTD(1,INODE)
        GMATX(1,IY)=R0
        GMATX(2,IX)=R0
        GMATX(2,IY)=CARTD(1,INODE)
        GMATX(3,IX)=CARTD(2,INODE)
        GMATX(3,IY)=R0
        GMATX(4,IX)=R0
        GMATX(4,IY)=CARTD(2,INODE)
   10 CONTINUE
      IF(NTYPE.EQ.3)THEN
C Axisymmetric problem
C --------------------
        IY=0
        DO 20 INODE=1,NNODE
          IX=IY+1
          IY=IX+1
          IF(NAXIS.EQ.1)THEN
C Axisymmetric about Y axis
            GMATX(5,IX)=SHAPE(INODE)/CARTCO(NAXIS)
            GMATX(5,IY)=R0
          ELSE IF(NAXIS.EQ.2)THEN
C Axisymmetric about X axis
            GMATX(5,IX)=R0
            GMATX(5,IY)=SHAPE(INODE)/CARTCO(NAXIS)
          ENDIF
   20   CONTINUE
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE GETGMX
CDOC BEGIN_INTEGER_FUNCTION INTNUM
CDOC Converts a character string into an integer
CDOC
CDOC This function returns the integer corresponding to the number
CDOC contained in the character string passed as argument.
CDOC
CDOC BEGIN_PARAMETERS
CDOC CHARACTER        CHRSTR >  Character string containing a number.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto  July  1996
CDOC
      INTEGER FUNCTION INTNUM(CHRSTR)
      IMPLICIT NONE
      CHARACTER*(*) CHRSTR
      INTEGER I, IASCII, IEND, IPOWER, LEN, LENGTH, NUMBER
C***********************************************************************
C CONVERTS A NUMBER CONTAINED IN A CHARACTER STRING INTO AN INTEGER
C***********************************************************************
 1000 FORMAT(/15X,'ERROR: String of blank characters passed'/
     1        22X,'into integer conversion function INTMUN')
 1100 FORMAT(/15X,'ERROR: Invalid character in string ''',A,''' passed'/
     1        22X,'into integer conversion function INTMUN')
C
      LENGTH=LEN(CHRSTR)
      DO 10 I=LENGTH,1,-1
        IF(CHRSTR(I:I).NE.' ')THEN
          IEND=I
          GOTO 20
        ENDIF
   10 CONTINUE
      WRITE(*,1000)
      WRITE(16,1000)
      CALL PEXIT
   20 CONTINUE
      INTNUM=0
      IPOWER=0
      DO 30 I=IEND,1,-1
        IASCII=ICHAR(CHRSTR(I:I))
        IF(IASCII.GE.48.AND.IASCII.LE.57)THEN
          NUMBER=IASCII-48
          INTNUM=INTNUM+NUMBER*(10**IPOWER)  
          IPOWER=IPOWER+1
        ELSEIF(CHRSTR(I:I).EQ.' ')THEN
          GOTO 40
        ELSEIF(CHRSTR(I:I).EQ.'-'.OR.CHRSTR(I:I).EQ.'+')THEN
          IF(I.NE.IEND)THEN
            IF(CHRSTR(I:I).EQ.'-')INTNUM=-INTNUM
            GOTO 40
          ELSE
            WRITE(*,1100)CHRSTR(1:IEND)
            WRITE(16,1100)CHRSTR(1:IEND)
            CALL PEXIT
          ENDIF
        ELSE
          WRITE(*,1100)CHRSTR(1:IEND)
          WRITE(16,1100)CHRSTR(1:IEND)
          CALL PEXIT
        ENDIF
   30 CONTINUE
   40 CONTINUE
      RETURN
      END
CDOC END_INTEGER_FUNCTION INTNUM
CDOC BEGIN_SUBROUTINE INVF2
CDOC Inverts the deformation gradient for 2-D problems
CDOC
CDOC This routine inverts deformation gradient tensors for plane strain,
CDOC plane stress and axisymmetric problems.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION F      >  Deformation gradient.
CDOC DOUBLE_PRECISION FINV   <  Inverse of the deformation gradient.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      SUBROUTINE INVF2
     1(   F          ,FINV       ,NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    F(3,3)             ,FINV(3,3)
      DATA
     1   R0   ,R1   /
     2   0.0D0,1.0D0/
C***********************************************************************
C INVERT DEFORMATION GRADIENT TENSORS FOR PLANE STRESS/STRAIN AND
C AXISYMMETRIC PROBLEMS
C***********************************************************************
      DETFPL=F(1,1)*F(2,2)-F(1,2)*F(2,1)
      IF(DETFPL.EQ.R0)CALL ERRPRT('EE0001')
      IF(NTYPE.EQ.3.AND.F(3,3).EQ.R0)CALL ERRPRT('EE0001')
C
      DETFIN=R1/DETFPL
      FINV(1,1)=F(2,2)*DETFIN
      FINV(2,2)=F(1,1)*DETFIN
      FINV(1,2)=-F(1,2)*DETFIN
      FINV(2,1)=-F(2,1)*DETFIN
      IF(NTYPE.EQ.2)THEN
        FINV(3,3)=R1
      ELSEIF(NTYPE.EQ.3)THEN
        FINV(3,3)=R1/F(3,3)
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE INVF2
CDOC BEGIN_SUBROUTINE LISTRA
CDOC Computes the infinitesimal strain components in 2-D
CDOC
CDOC Given the nodal displacements of the element and the B-matrix
CDOC (discrete symmetric gradient operator) at a point in the element
CDOC domain, this routine computes the corresponding (engineering)
CDOC infinitesimal strain components at that point by performing
CDOC the standard operation: e = B u, where e is the array of
CDOC engineering strain components and u is the array of nodal
CDOC displacements. This routine contains the plane strain/stress and
CDOC axisymmetric implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION BMATX  >  The discrete symmetric gradient
CDOC C                          operator, B-matrix.
CDOC DOUBLE_PRECISION ELDISP >  Array containing the element nodal
CDOC C                          displacements.
CDOC INTEGER          MDOFN  >  Dimensioning parameter: Number of rows
CDOC C                          of \smparm{ELDISP}.
CDOC INTEGER          MBDIM  >  Dimensioning parameter: Number of rows
CDOC C                          of \smparm{BMATX}.
CDOC INTEGER          NDOFN  >  Number of degrees of freedom per node.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION STRAN  <  Array of engineering infinitesimal
CDOC C                          strain components.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, June 1996: Initial coding
CDOC
      SUBROUTINE LISTRA
     1(   BMATX      ,ELDISP     ,MDOFN      ,MBDIM      ,NDOFN      ,
     2    NNODE      ,NTYPE      ,STRAN      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    BMATX(MBDIM,*)     ,ELDISP(MDOFN,*)    ,STRAN(*)
C***********************************************************************
C COMPUTES THE SYMMETRIC GRADIENT (LINEAR STRAIN MEASURE) ASSOCIATED
C WITH THE ELEMENT DISPLACEMENT 'ELDISP' IN 2-D: PLANE STRAIN, PLANE
C STRESS AND AXISYMMETRIC PROBLEMS
C***********************************************************************
      IF(NTYPE.EQ.1)THEN
        NSTRE=3
        NBDIM=3
      ELSEIF(NTYPE.EQ.2)THEN
        NSTRE=4
        NBDIM=3
      ELSEIF(NTYPE.EQ.3)THEN
        NSTRE=4
        NBDIM=4
      ELSE
        CALL ERRPRT('EI0023')
      ENDIF
C
      CALL RVZERO(STRAN,NSTRE)
      DO 30 ISTRE=1,NBDIM
        IEVAB=0
        DO 20 INODE=1,NNODE
          DO 10 IDOFN=1,NDOFN
            IEVAB=IEVAB+1
            STRAN(ISTRE)=STRAN(ISTRE)+
     1                   BMATX(ISTRE,IEVAB)*ELDISP(IDOFN,INODE)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE LISTRA
CDOC BEGIN_SUBROUTINE LOGSTR
CDOC Logarithmic strain computation
CDOC
CDOC Given the left (right) Cauchy-Green strain tensor, this routine
CDOC computes the corresponding Eulerian (Lagrangian) logarithmic strain
CDOC tensor (engineering components).
CDOC Plane strain, plane stress and axisymmetric implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION B      >  Array of components of the Cauchy-Green
CDOC C                          strain tensor.
CDOC DOUBLE_PRECISION E      <  Array of (engineering) components of the
CDOC C                          logarithmic strain tensor.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding (as EETRIA)
CDOC
CDOC E.de Souza Neto, May    1998: Routine and some variables renamed
CDOC
      SUBROUTINE LOGSTR
     1(   B          ,E          ,NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL  DLGD2
      LOGICAL   OUTOFP
      DIMENSION
     1    B(*)               ,E(*)
      DATA  R2   /2.0D0/
C***********************************************************************
C COMPUTES THE LOGARITHMIC STRAIN TENSOR:
C
C                       E :=  1/2 ln[ B ]
C
C***********************************************************************
      IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        OUTOFP=.TRUE.
      ELSEIF(NTYPE.EQ.1)THEN
        OUTOFP=.FALSE.
      ELSE
        CALL ERRPRT('EI0022')
      ENDIF
C
C Use isotropic tensor function to compute the logarithmic (physical)
C strain components
C
      CALL ISO2
     1(   DLGD2      ,OUTOFP     ,B          ,E          )
C
C Convert physical components into engineering strain components 
C
      E(3)=R2*E(3)
C
      RETURN
      END
CDOC END_SUBROUTINE LOGSTR
CDOC BEGIN_INTEGER_FUNCTION NWORD
CDOC Returns the number of words contained in a character string
CDOC
CDOC The return value of this function is the number of words contained
CDOC in the character string passed in its argument list. The function
CDOC also sets the pointers to the beginning and end of each word. The
CDOC pointers are returned via argument list.
CDOC
CDOC BEGIN_PARAMETERS
CDOC CHARACTER        CHRSTR >  Character string.
CDOC INTEGER          IWBEG  <  Array of pointers to the beginning of
CDOC C                          the words contained in \smparm{CHRSTR}.
CDOC C                          For the \smparm{N}th word, beginning at
CDOC C                          \smparm{CHRSTR(I:I)}, the function sets
CDOC C                          \smparm{IWBEG(N)=I}.
CDOC INTEGER          IWEND  <  Array of pointers to the end of
CDOC C                          the words contained in \smparm{CHRSTR}.
CDOC C                          For the \smparm{N}th word, ending at
CDOC C                          \smparm{CHRSTR(I:I)}, the function sets
CDOC C                          \smparm{IWEND(N)=I}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, July 1996: Initial coding
CDOC
      INTEGER FUNCTION NWORD(CHRSTR,IWBEG,IWEND)
      IMPLICIT NONE
      CHARACTER*(*) CHRSTR
      INTEGER IWBEG(*), IWEND(*)
      LOGICAL OUT
      INTEGER I, LEN, LENGTH
C***********************************************************************
C FIND NUMBER OF WORDS CONTAINED IN A CHARACTER STRING AND SET POINTERS
C TO BEGINNING AND END OF EACH WORD
C***********************************************************************
      LENGTH=LEN(CHRSTR)
      NWORD=0
      OUT=.TRUE.
      DO 10 I=1,LENGTH
        IF(OUT)THEN
          IF(CHRSTR(I:I).NE.' ')THEN
            OUT=.FALSE.
            NWORD=NWORD+1
            IWBEG(NWORD)=I
          ENDIF
        ELSE
          IF(CHRSTR(I:I).EQ.' ')THEN
            OUT=.TRUE.
            IWEND(NWORD)=I-1
          ELSEIF(I.EQ.LENGTH)THEN
            IWEND(NWORD)=I
          ENDIF
        ENDIF
   10 CONTINUE
      RETURN
      END
CDOC END_INTEGER_FUNCTION NWORD
CDOC BEGIN_SUBROUTINE PEXIT
CDOC Aborts execution of HYPLAS
CDOC
CDOC This routine closes the open files and stops the execution of
CDOC HYPLAS, sending a message to the results file and to the standard
CDOC output. It is called in emergency situations when a irrecoverable
CDOC error occurs.
CDOC
      SUBROUTINE PEXIT
C Print message
      WRITE(*,'(///15X,A,///)')'Program HYPLAS aborted.'
      WRITE(16,'(///15X,A,///)')'Program HYPLAS aborted.'
C Close files
      CALL FCLOSE
C and exit program
      STOP ' '
      END
CDOC END_SUBROUTINE PEXIT

CDOC BEGIN_SUBROUTINE RTSR
CDOC Matrix product s.Rt S R
CDOC
CDOC This routine performs the matrix product s Rt S R, where s is a
CDOC scalar, R a rectangular real matrix and S a square real matrix.
CDOC Rt denotes the transpose of R.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AUXM   <  Auxiliary matrix used to store partial
CDOC C                          results of the calculation.
CDOC INTEGER          MODE   >  If set to 1, the argument \smparm{Q}
CDOC C                          returns the resulting matrix Rt S R.
CDOC C                          Otherwise, Rt S R is added to the input
CDOC C                          value of \smparm{Q}.
CDOC INTEGER          MROWQ  >  Dimensioning parameter: maximum
CDOC C                          dimension of the square matrix
CDOC C                          \smparm{Q}.
CDOC INTEGER          MROWR  >  Dimensioning parameter: maximum number
CDOC C                          of rows of \smparm{R} (same as the
CDOC C                          maximum dimension of square matrix
CDOC C                          \smparm{S}).
CDOC INTEGER          NCOLR  >  Number of columns of \smparm{R}.
CDOC INTEGER          NROWR  >  Number of rows of \smparm{R}.
CDOC DOUBLE_PRECISION Q      <> Matrix where results are stored.
CDOC DOUBLE_PRECISION R      >  Rectangular real matrix.
CDOC DOUBLE_PRECISION S      >  Square real matrix.
CDOC DOUBLE_PRECISION SCAL   >  Real scalar.
CDOC LOGICAL          UNSYM  >  Unsymmetry flag.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      SUBROUTINE RTSR
     1(   AUXM       ,MODE       ,MROWQ      ,MROWR      ,NCOLR      ,
     2    NROWR      ,Q          ,R          ,S          ,SCAL       ,
     3    UNSYM      )  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UNSYM
      DIMENSION
     1    AUXM(NCOLR,NROWR)  ,Q(MROWQ,MROWQ)     ,R(MROWR,NCOLR)     ,
     2    S(MROWR,MROWR)
      DATA  R0   /0.0D0/
C***********************************************************************
C PERFORMS THE MATRIX PRODUCTS
C                                 T
C                    Q := SCAL * R  S R        (IF MODE=1)
C OR
C                                     T
C                    Q := Q + SCAL * R  S R    (OTHERWISE)
C
C WHERE 'R' IS A REAL RECTANGULAR MATRIX, 'S' A REAL SQUARE MATRIX
C AND 'SCAL' A SCALAR.
C***********************************************************************
      CALL RVZERO(AUXM,NCOLR*NROWR)
      DO 30 I=1,NCOLR
        DO 20 K=1,NROWR
          IF(R(K,I).NE.R0)THEN
            DO 10 J=1,NROWR
              AUXM(I,J)=AUXM(I,J)+SCAL*R(K,I)*S(K,J)
   10       CONTINUE
          ENDIF
   20   CONTINUE
   30 CONTINUE
C
      IF(MODE.EQ.1)THEN
        DO 50 I=1,NCOLR
          DO 40 J=1,NCOLR
            Q(I,J)=R0
   40     CONTINUE
   50   CONTINUE
      ENDIF
C
      IF(UNSYM)THEN
C Construct the whole matrix Q at once
        DO 80 J=1,NCOLR
          DO 70 K=1,NROWR
            IF(R(K,J).NE.R0)THEN
              DO 60 I=1,NCOLR
                Q(I,J)=Q(I,J)+AUXM(I,K)*R(K,J)
   60         CONTINUE
            ENDIF
   70     CONTINUE
   80   CONTINUE
      ELSE
C Construct the lower triangle of Q first
        DO 110 J=1,NCOLR
          DO 100 K=1,NROWR
            IF(R(K,J).NE.R0)THEN
              DO 90 I=J,NCOLR
                Q(I,J)=Q(I,J)+AUXM(I,K)*R(K,J)
   90         CONTINUE
            ENDIF
  100     CONTINUE
  110   CONTINUE
C and then assemble the upper triangle
        DO 130 I=1,NCOLR
          DO 120 J=I+1,NCOLR
            Q(I,J)=Q(J,I)
  120     CONTINUE
  130   CONTINUE
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE RTSR
CDOC BEGIN_SUBROUTINE RTSX
CDOC Matrix product s.Rt S X
CDOC
CDOC This routine performs the matrix product s Rt S X, where s is a
CDOC scalar, R and X rectangular real matrices of identical dimensions
CDOC and S a square real matrix.  Rt denotes the transpose of R.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AUXM   <  Auxiliary matrix used to store partial
CDOC C                          results of the calculation.
CDOC INTEGER          MODE   >  If set to 1, the argument \smparm{Q}
CDOC C                          returns the resulting matrix Rt S R.
CDOC C                          Otherwise, Rt S R is added to the input
CDOC C                          value of \smparm{Q}.
CDOC INTEGER          MROWQ  >  Dimensioning parameter: maximum
CDOC C                          dimension of the square matrix
CDOC C                          \smparm{Q}.
CDOC INTEGER          MROWR  >  Dimensioning parameter: maximum number
CDOC C                          of rows of \smparm{R} (same as the
CDOC C                          maximum dimension of square matrix
CDOC C                          \smparm{S}).
CDOC INTEGER          NCOLR  >  Number of columns of \smparm{R}.
CDOC INTEGER          NROWR  >  Number of rows of \smparm{R}.
CDOC DOUBLE_PRECISION Q      <> Matrix where results are stored.
CDOC DOUBLE_PRECISION R      >  Rectangular real matrix.
CDOC DOUBLE_PRECISION S      >  Square real matrix.
CDOC DOUBLE_PRECISION X      >  Rectangular real matrix.
CDOC DOUBLE_PRECISION SCAL   >  Real scalar.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, September 1996: Initial coding
CDOC
CDOC E.de Souza Neto & F.M.A.Pires , April 2002:
CDOC                          Bug fix in skipping zero multiplication
CDOC
      SUBROUTINE RTSX
     1(   AUXM       ,MODE       ,MROWQ      ,MROWR      ,NCOLR      ,
     2    NROWR      ,Q          ,R          ,S          ,X          ,
     3    SCAL       )  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    AUXM(NCOLR,NROWR)  ,Q(MROWQ,MROWQ)     ,R(MROWR,NCOLR)     ,
     2    S(MROWR,MROWR)     ,X(MROWR,NCOLR)
      DATA  R0   /0.0D0/
C***********************************************************************
C PERFORMS THE MATRIX PRODUCTS
C                                  T
C                     Q := SCAL * R  S X        (IF MODE=1)
C OR
C                                      T
C                     Q := Q + SCAL * R  S X    (OTHERWISE)
C
C WHERE 'R' AND 'X' ARE REAL RECTANGULAR MATRICES OF IDENTICAL
C DIMENSIONS, 'S' A REAL SQUARE MATRIX AND 'SCAL' A SCALAR.
C***********************************************************************
      CALL RVZERO(AUXM,NCOLR*NROWR)
      DO 30 I=1,NCOLR
        DO 20 K=1,NROWR
          IF(R(K,I).NE.R0)THEN
            DO 10 J=1,NROWR
              AUXM(I,J)=AUXM(I,J)+SCAL*R(K,I)*S(K,J)
   10       CONTINUE
          ENDIF
   20   CONTINUE
   30 CONTINUE
C
      IF(MODE.EQ.1)THEN
        DO 50 I=1,NCOLR
          DO 40 J=1,NCOLR
            Q(I,J)=R0
   40     CONTINUE
   50   CONTINUE
      ENDIF
C
C Construct the matrix Q
      DO 80 J=1,NCOLR
        DO 70 K=1,NROWR
          IF(X(K,J).NE.R0)THEN
            DO 60 I=1,NCOLR
              Q(I,J)=Q(I,J)+AUXM(I,K)*X(K,J)
   60       CONTINUE
          ENDIF
   70   CONTINUE
   80 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE RTSX
CDOC BEGIN_SUBROUTINE RTV
CDOC Matrix-vectot product s.Rt V
CDOC
CDOC This routine performs the matrix-vector product s Rt V, where s is
CDOC a scalar, R a real rectangular matrix and v a real vector.
CDOC Rt denotes the transpose of R.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          MODE   >  If set to 1, the argument \smparm{V}
CDOC C                          returns the resulting vector s Rt V.
CDOC C                          Otherwise, s Rt v is added to the input
CDOC C                          value of \smparm{P}.
CDOC INTEGER          MROWR  >  Dimensioning parameter: maximum number
CDOC C                          of rows of \smparm{R}.
CDOC INTEGER          NCOLR  >  Number of columns of \smparm{R}.
CDOC INTEGER          NROWR  >  Number of rows of \smparm{R}.
CDOC DOUBLE_PRECISION P      <> Vector where results are stored.
CDOC DOUBLE_PRECISION R      >  Rectangular real matrix.
CDOC DOUBLE_PRECISION V      >  Real vector.
CDOC DOUBLE_PRECISION SCAL   >  Real scalar.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      SUBROUTINE RTV
     1(   MODE       ,MROWR      ,NCOLR      ,NROWR      ,P          ,
     2    R          ,V          ,SCAL       )  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    P(NCOLR)           ,R(MROWR,NCOLR)     ,V(NROWR)
      DATA  R0   /0.0D0/
C***********************************************************************
C PERFORMS THE PRODUCT
C                                  T
C                     P := SCAL * R  V          (IF MODE=1)
C OR
C                                      T
C                     P := P + SCAL * R  V      (OTHERWISE)
C
C WHERE 'R' IS A REAL RECTANGULAR MATRIX, 'V' A REAL VECTOR AND
C 'SCAL' A SCALAR.
C***********************************************************************
      IF(MODE.EQ.1)CALL RVZERO(P,NCOLR)
      DO 30 I=1,NCOLR
        DO 20 J=1,NROWR
          IF(R(J,I).NE.R0)THEN
            P(I)=P(I)+SCAL*R(J,I)*V(J)
          ENDIF
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE RTV
CDOC BEGIN_SUBROUTINE SETBE
CDOC Obtains the left Cauchy-Green strain tensor from the log. strain
CDOC
CDOC Given the Eulerian (Lagrangian) logarithmic strain tensor, this
CDOC routine computes the corresponding left (right) Cauchy-Green
CDOC strain tensor.
CDOC Plane strain, plane stress and axisymmetric implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION B      <> Array of engineering logarithmic strain
CDOC C                          components on entry. Returns as the
CDOC C                          array of components of the corresponding
CDOC C                          Cauchy-Green strain tensor.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      SUBROUTINE SETBE
     1(   BE         ,NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL  EXP2X
      LOGICAL   OUTOFP
      DIMENSION BE(*)
      DATA  RP5  /0.5D0/
C***********************************************************************
C COMPUTES THE CAUCHY-GREEN TENSOR AS A FUNCTION OF
C THE LOGARITHMIC STRAIN TENSOR:
C
C                      B    :=  exp[ 2 E   ]
C
C***********************************************************************
      IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        OUTOFP=.TRUE.
      ELSEIF(NTYPE.EQ.1)THEN
        OUTOFP=.FALSE.
      ELSE
        CALL ERRPRT('EI0024')
      ENDIF
C Convert engineering elastic strain components into physical components
      BE(3)=RP5*BE(3)
C Use isotropic tensor function to compute elastic Cauchy-Green tensor
      CALL ISO2
     1(   EXP2X      ,OUTOFP     ,BE         ,BE         )
C
      RETURN
      END
CDOC END_SUBROUTINE SETBE

CDOC BEGIN_DOUBLE_PRECISION_FUNCTION DDLGD2
CDOC ddlg2(x)=1/(2x)
CDOC
CDOC This is the derivative of the function defined in \smparm{DLGD2},
CDOC that relates the principal logarithmic stretches and the
CDOC eigenvalues of the Cauchy-Green strain tensor.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the function will be
CDOC C                          evaluated.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      DOUBLE PRECISION FUNCTION DDLGD2(X)
C
      DOUBLE PRECISION X
      DOUBLE PRECISION RP5
      DATA   RP5 /0.5D0/
C***********************************************************************
C DERIVATIVE OF THE SCALAR FUNCTION 'DLGD2' THAT RELATES PRINCIPAL
C LOGARITHMIC STRECTHES AND EIGENVALUES OF THE CAUCHY-GREEN TENSOR
C***********************************************************************
      DDLGD2=RP5/X
C
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION DDLGD2
CDOC BEGIN_SUBROUTINE DGISO2
CDOC Derivative of a general isotropic tensor function of one tensor
CDOC
CDOC This function computes the derivative, dY(X)/dX, of a general
CDOC isotropic tensor function of one tensor, Y(X). This implementation
CDOC is restricted to 2-D with one possible out-of-plane component
CDOC (normally needed in axisymmetric problems). The tensor function
CDOC Y(X) is assumed to be defined as Y(X)= Sum[yi(x1,x2,x3) ei(x)ei],
CDOC where yi are the eigenvalues of the tensor Y and xi the eigenvalues
CDOC of the tensor X. ei are the egenvectors of X (which by definition
CDOC of Y(X), coincide with those of tensor Y) and "(x)" denotes the
CDOC tensor product.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DEIGY  >  Matrix containing the derivatives
CDOC C                          dyi/dxj of the eigenvalues of Y(X) with
CDOC C                          respect to the eigenvalues of X.
CDOC DOUBLE_PRECISION DYDX   <  Matrix of components of the derivative
CDOC C                          (fourth order tensor) dY/dX.
CDOC DOUBLE_PRECISION EIGPRJ >  Matrix with each column containing the
CDOC C                          components of one eigenprojection tensor
CDOC C                          ei(x)ei.
CDOC DOUBLE_PRECISION EIGX   >  Array of eigenvalues of X.
CDOC DOUBLE_PRECISION EIGY   >  Array of eigenvalues of Y.
CDOC LOGICAL          OUTOFP >  Out-of-plane component flag. If set to
CDOC C                          \smparm{.TRUE.} the out-of-plane
CDOC C                          component (normally required in
CDOC C                          axisymmetric problems) is computed.
CDOC C                          The out-of-plane component is not
CDOC C                          computed otherwise.
CDOC LOGICAL          REPEAT >  Repeated in-plane eigenvalues flag.
CDOC C                          If the in-plane eigenvalues are repeated
CDOC C                          this argument must be set to
CDOC C                          \smaparm{.TRUE.} on entry, so that the
CDOC C                          appropriate limit expression for the
CDOC C                          the derivative is employed.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, May 1996: Initial coding
CDOC
      SUBROUTINE DGISO2
     1(   DEIGY      ,DYDX       ,EIGPRJ     ,EIGX       ,EIGY       ,
     2    OUTOFP     ,REPEAT     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=4    ,MDIM=3     ,NDIM=2     )
C Arguments
      LOGICAL OUTOFP, REPEAT
      DIMENSION
     1    DEIGY(MDIM,MDIM)  ,DYDX(MCOMP,MCOMP) ,EIGPRJ(MCOMP,NDIM),
     2    EIGX(NDIM)        ,EIGY(NDIM)
C Local arrays
      DIMENSION
     1    EIGPR3(MCOMP)     ,FOID(MCOMP,MCOMP) ,SOPID(MCOMP)
      DATA
     1    FOID(1,1)     ,FOID(1,2)     ,FOID(1,3)     /
     2    1.0D0         ,0.0D0         ,0.0D0         /
     3    FOID(2,1)     ,FOID(2,2)     ,FOID(2,3)     /
     4    0.0D0         ,1.0D0         ,0.0D0         /
     5    FOID(3,1)     ,FOID(3,2)     ,FOID(3,3)     /
     6    0.0D0         ,0.0D0         ,0.5D0         /
      DATA
     1    SOPID(1)      ,SOPID(2)      ,SOPID(3)      ,SOPID(4)        /
     2    1.0D0         ,1.0D0         ,0.0D0         ,0.0D0           /
      DATA
     1    EIGPR3(1)     ,EIGPR3(2)     ,EIGPR3(3)     ,EIGPR3(4)       /
     2    0.0D0         ,0.0D0         ,0.0D0         ,1.0D0           /
C***********************************************************************
C COMPUTE THE DERIVATIVE OF A GENERAL ISOTROPIC TENSOR FUNCTION OF ONE
C TENSOR IN 2-D (WITH ONE POSSIBLE OUT-OF-PLANE COMPONENT)
C***********************************************************************
      CALL RVZERO(DYDX,MCOMP*MCOMP)
      IF(REPEAT)THEN
C Derivative dY/dX for repeated in-plane eigenvalues of X
C -------------------------------------------------------
C In-plane component
        DO 20 I=1,3
          DO 10 J=1,3
            DYDX(I,J)=(DEIGY(1,1)-DEIGY(1,2))*FOID(I,J)+
     1                 DEIGY(1,2)*SOPID(I)*SOPID(J)
   10     CONTINUE
   20   CONTINUE
        IF(OUTOFP)THEN
C out-of-plane components required
          DO 40 I=1,4
            DO 30 J=1,4
              IF(I.EQ.4.OR.J.EQ.4)DYDX(I,J)=
     1                    DEIGY(1,3)*SOPID(I)*EIGPR3(J)+
     2                    DEIGY(3,1)*EIGPR3(I)*SOPID(J)+
     3                    DEIGY(3,3)*EIGPR3(I)*EIGPR3(J)
   30       CONTINUE
   40     CONTINUE
        ENDIF
      ELSE
C Derivative dY/dX for distinct in-plane eigenvalues of X
C -------------------------------------------------------
C Assemble in-plane DYDX
        A1=(EIGY(1)-EIGY(2))/(EIGX(1)-EIGX(2))
        DO 70 I=1,3
          DO 60 J=1,3
            DYDX(I,J)=A1*(FOID(I,J)-EIGPRJ(I,1)*EIGPRJ(J,1)-
     1                EIGPRJ(I,2)*EIGPRJ(J,2))+
     2                DEIGY(1,1)*EIGPRJ(I,1)*EIGPRJ(J,1)+
     3                DEIGY(1,2)*EIGPRJ(I,1)*EIGPRJ(J,2)+
     4                DEIGY(2,1)*EIGPRJ(I,2)*EIGPRJ(J,1)+
     5                DEIGY(2,2)*EIGPRJ(I,2)*EIGPRJ(J,2)
   60     CONTINUE
   70   CONTINUE
        IF(OUTOFP) THEN
C out-of-plane components required
          DO 90 I=1,4
            DO 80 J=1,4
              IF(I.EQ.4.OR.J.EQ.4)DYDX(I,J)=
     1                DEIGY(1,3)*EIGPRJ(I,1)*EIGPR3(J)+
     2                DEIGY(2,3)*EIGPRJ(I,2)*EIGPR3(J)+
     3                DEIGY(3,1)*EIGPR3(I)*EIGPRJ(J,1)+
     4                DEIGY(3,2)*EIGPR3(I)*EIGPRJ(J,2)+
     5                DEIGY(3,3)*EIGPR3(I)*EIGPR3(J)
   80       CONTINUE
   90     CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE DGISO2
CDOC BEGIN_SUBROUTINE DISO2
CDOC Derivative of a class of isotropic tensor function of one tensor
CDOC
CDOC This function computes the derivative, dY(X)/dX, of a particular
CDOC class of isotropic tensor valued function of one tensor, Y(X).
CDOC This implementation is restricted to 2-D with one possible
CDOC out-of-plane component (normally needed in axisymmetric problems).
CDOC The class of tensor functions Y(X) is assumed to be defined as
CDOC Y(X)= Sum[y(xi) ei(x)ei], where the scalar function y(xi) defines
CDOC the eigenvalues of the
CDOC tensor Y and xi the eigenvalues of the tensor X. ei are the
CDOC egenvectors of X (which by definition of Y(X), coincide with
CDOC those of tensor Y) and "(x)" denotes the tensor product.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DYDX   <  Matrix of components of the derivative
CDOC C                          (fourth order tensor) dY/dX.
CDOC LOGICAL          OUTOFP >  Out-of-plane component flag. If set to
CDOC C                          \smparm{.TRUE.} the out-of-plane
CDOC C                          component (normally required in
CDOC C                          axisymmetric problems) is computed.
CDOC C                          The out-of-plane component is not
CDOC C                          computed otherwise.
CDOC DOUBLE_PRECISION X      >  Point at which the derivative is to be
CDOC C                          computed.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
cccccCDOC ???????          DFUNC  >  Symbolic name of the double precision
cccccCDOC C                          function defining dy(x)/dx.
cccccCDOC ???????          FUNC   >  Symbolic name of the double precision
cccccCDOC C                          function defining y(x).
      SUBROUTINE DISO2
     1(   DYDX       ,DFUNC      ,FUNC       ,OUTOFP     ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL
     1    DFUNC      ,FUNC
      PARAMETER
     1(   MCOMP=4    ,NDIM=2     )
C Arguments
      LOGICAL OUTOFP
      DIMENSION
     1    DYDX(MCOMP,MCOMP) ,X(*)
C Local variables and arrays
      LOGICAL REPEAT
      DIMENSION
     1    DEIGY(NDIM)       ,EIGPRJ(MCOMP,NDIM),EIGX(NDIM)        ,
     2    EIGY(NDIM)        ,FOID(MCOMP,MCOMP)
      DATA
     1    FOID(1,1)     ,FOID(1,2)     ,FOID(1,3)     /
     2    1.0D0         ,0.0D0         ,0.0D0         /
     3    FOID(2,1)     ,FOID(2,2)     ,FOID(2,3)     /
     4    0.0D0         ,1.0D0         ,0.0D0         /
     5    FOID(3,1)     ,FOID(3,2)     ,FOID(3,3)     /
     6    0.0D0         ,0.0D0         ,0.5D0         /
C***********************************************************************
C COMPUTE (AND STORE IN MATRIX FORM) THE DERIVATIVE dY/dX OF AN
C ISOTROPIC TENSOR FUNCTION OF THE TYPE:
C
C                        Y(X) = sum{ y(x_i) E_i }
C
C WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
C THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
C FUNCTION. THIS ROUTINE IS RESTRICTED TO 2-D TENSORS WITH ONE
C POSSIBLE (TRANSVERSAL) OUT-OF-PLANE COMPONENT.
C***********************************************************************
C Spectral decomposition of X
      CALL SPDEC2
     1(   EIGPRJ     ,EIGX       ,REPEAT     ,X          )
C In-plane eigenvalues of Y (and derivatives)
      DO 10 IDIR=1,2
        EIGY(IDIR)=FUNC(EIGX(IDIR))
        DEIGY(IDIR)=DFUNC(EIGX(IDIR))
   10 CONTINUE
C
C In-plane components of dY/dX
C ----------------------------
      CALL RVZERO(DYDX,MCOMP*MCOMP)
      IF(REPEAT)THEN
C for repeated in-plane eigenvalues of X
        DO 20 I=1,3
          DYDX(I,I)=DEIGY(1)*FOID(I,I)
   20   CONTINUE
      ELSE
C for distinct in-plane eigenvalues of X
        A1=(EIGY(1)-EIGY(2))/(EIGX(1)-EIGX(2))
        DO 40 I=1,3
          DO 30 J=I,3
            DYDX(I,J)=A1*(FOID(I,J)-EIGPRJ(I,1)*EIGPRJ(J,1)-
     1                EIGPRJ(I,2)*EIGPRJ(J,2))+
     2                DEIGY(1)*EIGPRJ(I,1)*EIGPRJ(J,1)+
     3                DEIGY(2)*EIGPRJ(I,2)*EIGPRJ(J,2)
            IF(I.NE.J)DYDX(J,I)=DYDX(I,J)
   30     CONTINUE
   40   CONTINUE
      ENDIF
C
C Out-of-plane component required
C -------------------------------
      IF(OUTOFP)DYDX(4,4)=DFUNC(X(4))
C
      RETURN
      END
CDOC END_SUBROUTINE DISO2
CDOC BEGIN_DOUBLE_PRECISION_FUNCTION DLGD2
CDOC dlg2(x)=log(x)/2
CDOC
CDOC This function relates the principal logarithmic stretches and the
CDOC eigenvalues of the Cauchy-Green strain tensor.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the function will be
CDOC C                          evaluated.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996; Initial coding
CDOC
      DOUBLE PRECISION FUNCTION DLGD2(X)
C
      DOUBLE PRECISION X
      DOUBLE PRECISION RP5
      DATA   RP5 /0.5D0/
C***********************************************************************
C SCALAR FUNCTION THAT RELATES PRINCIPAL LOGARITHMIC STRECTHES AND
C EIGENVALUES OF THE CAUCHY-GREEN TENSOR
C***********************************************************************
      DLGD2=RP5*DLOG(X)
C
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION DLGD2
CDOC BEGIN_DOUBLE_PRECISION_FUNCTION DPLFUN
CDOC Returns the derivative of piece-wise linear scalar function
CDOC
CDOC This procedure returns the derivative of the piece-wise linear
CDOC scalar function of procedure \smparm{PLFUN}.
CDOC The piece-wise linear function F(X) is defined by a set
CDOC of \smparm{NPOINT} pairs (X,F(X)) passed in the matrix argument
CDOC \smparm{XFX}.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the derivative will be
CDOC C                          evaluated.
CDOC INTEGER          NPOINT >  Number of points defining the piece-wise
CDOC C                          linear function.
CDOC DOUBLE_PRECISION XFX    >  Matrix (dimension 2*\smparm{NPOINT})
CDOC C                          containing the pairs (x,f(x)) which
CDOC C                          define the piece-wise linear function.
CDOC C                          evaluated. Each column of \smparm{XFX}
CDOC C                          contains a pair (xi,f(xi)). The pairs
CDOC C                          supplied in \smparm{XFX}
CDOC C                          must be ordered such that the x's are
CDOC C                          monotonically increasing. That is, the
CDOC C                          x [XFX(1,i+1)] of a column i+1 must be
CDOC C                          greater than XFX(1,i) (x of column i).
CDOC C                          If \smparm{X}<\smparm{XFX(1,1)} the
CDOC C                          piece-wise linear function is assumed
CDOC C                          constant equal to \smparm{XFX(1,1)}.
CDOC C                          If \smparm{X}>\smparm{XFX(1,NPOINT)} the
CDOC C                          piece-wise linear function is assumed
CDOC C                          constant equal to
CDOC C                          \smparm{XFX(1,NPOINT)}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1992: Initial coding
CDOC
      DOUBLE PRECISION FUNCTION DPLFUN(X, NPOINT, XFX)
C
      INTEGER NPOINT, I
      DOUBLE PRECISION X, XFX(2,NPOINT), R0
      DATA R0 / 0.0D0 /
C***********************************************************************
C DERIVATIVE OF THE PIECEWISE LINEAR FUNCTION DEFINED BY A SET OF
C NPOINT PAIRS {X,F(X)} STORED IN THE MATRIX XFX (DIM. 2*NPOINT).
C***********************************************************************
      DO 100 I=1,NPOINT 
        IF (X.GE.XFX(1,I)) THEN
          GOTO 100
        ELSE
          IF (I.EQ.1) THEN
C           -- x < x1   --> f(x)=f(x1) --> df(x)/dx=0 --- 
            DPLFUN=R0
            GOTO 999
          ELSE
C           -- x(i-1) <= x < x(i) ---
            DPLFUN=(XFX(2,I)-XFX(2,I-1))/
     1             (XFX(1,I)-XFX(1,I-1))
            GOTO 999
          ENDIF
        ENDIF
 100  CONTINUE
C     ---- x >= x(npoint) --> f(x) = f(x(npoint)) --> df/dx=0 ---
      DPLFUN=R0 
 999  CONTINUE
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION DPLFUN
CDOC BEGIN_DOUBLE_PRECISION_FUNCTION EXP2X
CDOC f(x)=exp(2x)
CDOC
CDOC This function relates the eigenvalues of the Cauchy-Green strain
CDOC tensors to the principal logarithmic stretches.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the function will be
CDOC C                          evaluated.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
      DOUBLE PRECISION FUNCTION EXP2X(X)
C
      DOUBLE PRECISION X
      DOUBLE PRECISION R2
      DATA   R2  /2.0D0/
C***********************************************************************
C SCALAR FUNCTION THAT RELATES EIGENVALUES OF THE CAUCHY-GREEN
C TENSOR TO THE PRINCIPAL LOGARITHMIC STRECTHES
C***********************************************************************
      EXP2X=DEXP(R2*X)
C
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION EXP2X
CDOC BEGIN_SUBROUTINE INVMT3
CDOC Inverts a 3x3 double precision matrix
CDOC
CDOC This routine inverts a generally unsymmetric 3x3 double precision
CDOC matrix.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION S      >  Matrix to be einverted
CDOC DOUBLE_PRECISION SINV   <  Inverse matrix
CDOC DOUBLE_PRECISION DETS   <  Determinant of matrix S
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, July 2001: Initial coding
CDOC
      SUBROUTINE INVMT3
     1(   S          ,SINV       ,DETS       )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    S(3,3)             ,SINV(3,3)
      DATA
     1   R0   ,R1   /
     2   0.0D0,1.0D0/
C***********************************************************************
C INVERT A REAL 3x3 MATRIX
C***********************************************************************
      DETS=S(1,1)*S(2,2)*S(3,3)+S(1,2)*S(2,3)*S(3,1)+
     1     S(1,3)*S(2,1)*S(3,2)-S(1,2)*S(2,1)*S(3,3)-
     2     S(1,1)*S(2,3)*S(3,2)-S(1,3)*S(2,2)*S(3,1)
      IF(DETS.EQ.R0)CALL ERRPRT('EE????')
C
      DETSIN=R1/DETS
      SINV(1,1)=+DETSIN*(S(2,2)*S(3,3)-S(2,3)*S(3,2))
      SINV(2,1)=-DETSIN*(S(2,1)*S(3,3)-S(2,3)*S(3,1))
      SINV(3,1)=+DETSIN*(S(2,1)*S(3,2)-S(2,2)*S(3,1))
      SINV(1,2)=-DETSIN*(S(1,2)*S(3,3)-S(1,3)*S(3,2))
      SINV(2,2)=+DETSIN*(S(1,1)*S(3,3)-S(1,3)*S(3,1))
      SINV(3,2)=-DETSIN*(S(1,1)*S(3,2)-S(1,2)*S(3,1))
      SINV(1,3)=+DETSIN*(S(1,2)*S(2,3)-S(1,3)*S(2,2))
      SINV(2,3)=-DETSIN*(S(1,1)*S(2,3)-S(1,3)*S(2,1))
      SINV(3,3)=+DETSIN*(S(1,1)*S(2,2)-S(1,2)*S(2,1))
C
      RETURN
      END
CDOC END_SUBROUTINE INVMT3
CDOC BEGIN_SUBROUTINE ISO2
CDOC Computes the value isotropic tensor functions of one tensor.
CDOC
CDOC This subroutine evaluates isotropic tensor functions Y(X), of one
CDOC tensor belongin to the class described below.
CDOC This implementation is restricted to 2-D with one possible
CDOC out-of-plane component (normally needed in axisymmetric problems).
CDOC The class of tensor functions Y(X) is assumed to be defined as
CDOC Y(X)= Sum[y(xi) ei(x)ei], where the scalar function y(xi) defines
CDOC the eigenvalues of the
CDOC tensor Y and xi the eigenvalues of the tensor X. ei are the
CDOC egenvectors of X (which by definition of Y(X), coincide with
CDOC those of tensor Y) and "(x)" denotes the tensor product.
CDOC
CDOC BEGIN_PARAMETERS
CDOC LOGICAL          OUTOFP >  Out-of-plane component flag. If set to
CDOC C                          \smparm{.TRUE.} the out-of-plane
CDOC C                          component (normally required in
CDOC C                          axisymmetric problems) is computed.
CDOC C                          The out-of-plane component is not
CDOC C                          computed otherwise.
CDOC DOUBLE_PRECISION X      >  Array of components of the tensor at
CDOC C                          which the function is to be evaluated.
CDOC DOUBLE_PRECISION Y      <  Array of components of the tensor
CDOC C                          function at \smparm{X}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1996: Initial coding
CDOC
cccccCDOC ???????          FUNC   >  Symbolic name of the double
cccccCDOC C                          precision function defining y(x).
      SUBROUTINE ISO2
     1(   FUNC       ,OUTOFP     ,X          ,Y          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FUNC
      PARAMETER
     1(   MCOMP=4    ,NDIM=2     )
      LOGICAL OUTOFP ,REPEAT
      DIMENSION
     1    X(*)                      ,Y(*)
      DIMENSION
     1    EIGPRJ(MCOMP,NDIM)        ,EIGX(NDIM)                ,
     1    EIGY(NDIM)
C***********************************************************************
C COMPUTE THE TENSOR Y (STORED IN VECTOR FORM) AS AN ISOTROPIC
C FUNCTION OF THE TYPE:
C
C                     Y(X) = sum{ y(x_i) E_i }
C
C WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
C THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
C FUNCTION. THIS ROUTINE IS RESTRICTED TO 2-D TENSORS WITH ONE
C POSSIBLE (TRANSVERSAL) OUT-OF-PLANE COMPONENT.
C***********************************************************************
C Performs the spectral decomposition of X
      CALL SPDEC2
     1(   EIGPRJ     ,EIGX       ,REPEAT     ,X          )
C Computes the in-plane eigenvalues of Y
      DO 10 IDIR=1,2
        EIGY(IDIR)=FUNC(EIGX(IDIR))
   10 CONTINUE
C Assembles in-plane component of Y (in vector form)
      CALL RVZERO(Y,3)
      DO 30 ICOMP=1,3
        DO 20 IDIR=1,2
          Y(ICOMP)=Y(ICOMP)+EIGY(IDIR)*EIGPRJ(ICOMP,IDIR)
   20   CONTINUE
   30 CONTINUE
C Out-of-plane component required
      IF(OUTOFP)Y(4)=FUNC(X(4))
C
      RETURN
      END
CDOC END_SUBROUTINE ISO2
CDOC BEGIN_SUBROUTINE IVZERO
CDOC Zero an integer array
CDOC
CDOC This routine initialises to zero the \smparm{N} components of the
CDOC integer array argument \smparm{IV}.
CDOC 
CDOC BEGIN_PARAMETERS 
CDOC INTEGER          IV     <  Zeroed integer array.
CDOC INTEGER          N      >  Dimension of \smparm{IV}.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE IVZERO
     1(   IV         ,N          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IV(N)
C***********************************************************************
C INITIALISES TO ZERO AN INTEGER ARRAY OF DIMENSION N
C***********************************************************************
      DO 10 I=1,N
        IV(I)=0
   10 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE IVZERO
CDOC BEGIN_SUBROUTINE JACOB
CDOC Jacobi procedure for spectral decomposition of a symmetric matrix
CDOC
CDOC This routine uses the Jacobi iterative procedure for the spectral
CDOC decomposition (decomposition into eigenvalues and eigenvectors)
CDOC of a symmetric matrix.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A      <> Matrix to be decomposed.
CDOC DOUBLE_PRECISION D      <  Array containing eigenvalues of
CDOC C                          \smparm{A}.
CDOC DOUBLE_PRECISION V      <  Matrix containing one eigenvector of
CDOC C                          \smparm{A} in each column.
CDOC INTEGER          N      >  Dimension of \smparm{A}.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE JACOB(A,D,V,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MJITER=50,NMAX=100)
      DIMENSION
     1    A(N,N)     ,D(N)      ,V(N,N)
      DIMENSION
     1    B(NMAX)    ,Z(NMAX)
      DATA R0   ,RP2  ,RP5  ,R1   ,R100   /
     1     0.0D0,0.2D0,0.5D0,1.0D0,100.0D0/
      DATA TOLER  /
     1     1.0D-12/
C***********************************************************************
C JACOBI ITERATIVE PROCEDURE FOR SPECTRAL DECOMPOSITION OF A
C N-DIMENSIONAL SYMMETRIC MATRIX
C***********************************************************************
      IF(N.GT.NMAX)THEN
        CALL ERRPRT('EI0025')
      ENDIF
      DO 20 IP=1,N
        DO 10 IQ=1,N
          V(IP,IQ)=R0
   10   CONTINUE
        V(IP,IP)=R1
   20 CONTINUE
      DO 30 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=R0
   30 CONTINUE
      DO 130 I=1,MJITER
        SM=R0
        DO 50 IP=1,N-1
          DO 40 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
   40     CONTINUE
   50   CONTINUE
        IF(SM.LT.TOLER)GOTO 999
        IF(I.LT.4)THEN
          TRESH=RP2*SM/DBLE(N**2)
        ELSE
          TRESH=R0
        ENDIF
        DO 110 IP=1,N-1
          DO 100 IQ=IP+1,N
            G=R100*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     1         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=R0
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=RP5*H/A(IP,IQ)
                T=R1/(ABS(THETA)+SQRT(R1+THETA**2))
                IF(THETA.LT.R0)T=-T
              ENDIF
              C=R1/SQRT(R1+T**2)
              S=T*C
              TAU=S/(R1+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=R0
              DO 60 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
   60         CONTINUE
              DO 70 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
   70         CONTINUE
              DO 80 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
   80         CONTINUE
              DO 90 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
   90         CONTINUE
            ENDIF
  100     CONTINUE
  110   CONTINUE
        DO 120 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=R0
  120   CONTINUE
  130 CONTINUE
      CALL ERRPRT('EE0005')
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE JACOB
CDOC BEGIN_DOUBLE_PRECISION_FUNCTION PLFUN
CDOC Returns the value of a piece-wise linear scalar function
CDOC
CDOC This function returns the value of 
CDOC a piece-wise linear scalar function F(X) defined by a set
CDOC of \smparm{NPOINT} pairs (X,F(X)) passed in the matrix argument
CDOC \smparm{XFX}.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the function will be
CDOC C                          evaluated.
CDOC INTEGER          NPOINT >  Number of points defining the piece-wise
CDOC C                          linear function.
CDOC DOUBLE_PRECISION XFX    >  Matrix (dimension 2*\smparm{NPOINT})
CDOC C                          containing the pairs (x,f(x)) which
CDOC C                          define the piece-wise linear function.
CDOC C                          evaluated. Each column of \smparm{XFX}
CDOC C                          contains a pair (xi,f(xi)). The pairs
CDOC C                          supplied in \smparm{XFX}
CDOC C                          must be ordered such that the x's are
CDOC C                          monotonically increasing. That is, the
CDOC C                          x [XFX(1,i+1)] of a column i+1 must be
CDOC C                          greater than XFX(1,i) (x of column i).
CDOC C                          If \smparm{X}<\smparm{XFX(1,1)} the
CDOC C                          piece-wise linear function is assumed
CDOC C                          constant equal to \smparm{XFX(1,1)}.
CDOC C                          If \smparm{X}>\smparm{XFX(1,NPOINT)} the
CDOC C                          piece-wise linear function is assumed
CDOC C                          constant equal to
CDOC C                          \smparm{XFX(1,NPOINT)}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, August 1992: Initial coding
CDOC
      DOUBLE PRECISION FUNCTION PLFUN(X, NPOINT, XFX)
C
      INTEGER NPOINT, I
      DOUBLE PRECISION X, XFX(2,*)
C***********************************************************************
C PIECEWISE LINEAR FUNCTION DEFINED BY A SET OF NPOINT PAIRS
C {X,F(X)} STORED IN THE MATRIX XFX (DIM. 2*NPOINT).
C***********************************************************************
      DO 100 I=1,NPOINT 
        IF (X.GE.XFX(1,I)) THEN
          GOTO 100
        ELSE  
          IF (I.EQ.1) THEN
C           -- x < x1 --> f(x)=f(x1) --- 
            PLFUN=XFX(2,1)
            GOTO 999
          ELSE
C           -- x(i-1) <= x < x(i) ---
            PLFUN=XFX(2,I-1)+(X-XFX(1,I-1))*
     1                 (XFX(2,I)-XFX(2,I-1))/
     2                 (XFX(1,I)-XFX(1,I-1))
            GOTO 999
          ENDIF
        ENDIF
 100  CONTINUE
C     ----  x >= x(npoint) --> f(x) = f(x(npoint))  ---
 999  CONTINUE
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION PLFUN
CDOC BEGIN_SUBROUTINE QUATEX
CDOC Quaternion/Rotation extraction procedure.
CDOC
CDOC This routine either extracts the unit quaternion associated with
CDOC a given rotation matrix or extracts the rotation matrix associated
CDOC with a given unit quaternion.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          MODE   >  Mode flag. If \sparm{MODE}=0, this
CDOC C                          subroutine computes the unit quaternion
CDOC C                          associated with the given rotation
CDOC C                          matrix \smparm{R} and returns it stored
CDOC C                          in \smparm{Q}. If \sparm{MODE}=1, then
CDOC C                          the rotation matrix associated with the
CDOC C                          given unit quaternion \smparm{Q} is
CDOC C                          computed and returned stored in
CDOC C                          \smparm{R}.
CDOC DOUBLE_PRECISION Q      <> Unit quaternion.
CDOC DOUBLE_PRECISION R      <> Rotation matrix.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, February 1998: Initial coding
CDOC
      SUBROUTINE QUATEX( MODE   ,Q      ,R      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION
     1    Q(0:3)             ,R(3,3)
      DATA
     1    RP25  ,RP5  ,R1   ,R2   /
     2    0.25D0,0.5D0,1.0D0,2.0D0/
C***********************************************************************
C PERFORMS QUATERNION/ROTATION EXTRACTION.
C
C References:
C
C R.A.Spurrier, Comment on ``singularity-free extraction of a quaternion
C from a direction-cosine matrix''. J.Spacecraft, 15, pp.255, 1978.
C
C J.C.SImo and L.Vu-Quoc, A three-dimensional finite strain rod model.
C Part II: Computational aspects. Comp.Meth.Appl.Mech.Engng., 58,
C pp. 79-116, 1988.
C***********************************************************************
      IF(MODE.EQ.0)THEN
C extracts unit quaternion associated with rotation matrix R
        TRACE=R(1,1)+R(2,2)+R(3,3)
        RMAX=TRACE
        IF(R(1,1).GT.RMAX)THEN
          I=1
          J=2
          K=3
          RMAX=R(1,1)
        ENDIF
        IF(R(2,2).GT.RMAX)THEN
          I=2
          J=3
          K=1
          RMAX=R(2,2)
        ENDIF
        IF(R(3,3).GT.RMAX)THEN
          I=3
          J=1
          K=2
          RMAX=R(3,3)
        ENDIF
        IF(RMAX.EQ.TRACE)THEN
          Q(0)=RP5*SQRT(R1+TRACE)
          Q(1)=RP25*(R(3,2)-R(2,3))/Q(0)
          Q(2)=RP25*(R(1,3)-R(3,1))/Q(0)
          Q(3)=RP25*(R(2,1)-R(1,2))/Q(0)
        ELSE
          Q(I)=SQRT(RP5*R(I,I)+RP25*(R1-TRACE))
          Q(0)=RP25*(R(K,J)-R(J,K))/Q(I)
          Q(J)=RP25*(R(J,I)+R(I,J))/Q(I)
          Q(K)=RP25*(R(K,I)+R(I,K))/Q(I)
        ENDIF
      ELSEIF(MODE.EQ.1)THEN
C computes the rotation matrix R associated with the unit quaternion Q
        R(1,1)=R2*(Q(0)*Q(0)+Q(1)*Q(1)-RP5)
        R(1,2)=R2*(Q(1)*Q(2)-Q(3)*Q(0))
        R(1,3)=R2*(Q(1)*Q(3)+Q(2)*Q(0))
        R(2,1)=R2*(Q(2)*Q(1)+Q(3)*Q(0))
        R(2,2)=R2*(Q(0)*Q(0)+Q(2)*Q(2)-RP5)
        R(2,3)=R2*(Q(2)*Q(3)-Q(1)*Q(0))
        R(3,1)=R2*(Q(3)*Q(1)-Q(2)*Q(0))
        R(3,2)=R2*(Q(3)*Q(2)+Q(1)*Q(0))
        R(3,3)=R2*(Q(0)*Q(0)+Q(3)*Q(3)-RP5)
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE QUATEX
CDOC BEGIN_SUBROUTINE RVSCAL
CDOC Multiplies a double precision vector by a scalar
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION V      <> Double precision vector.
CDOC INTEGER          N      >  Dimension of \smparm{V}.
CDOC DOUBLE_PRECISION SCAL   <> Scalar by which \smparm{V} will be
CDOC C                          multiplied.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE RVSCAL
     1(   V          ,N          ,SCAL       )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N)
C***********************************************************************
C MULTIPLIES THE DOUBLE PRECISION VECTOR 'V', OF DIMENSION 'N',
C BY THE SCALAR 'SCAL'
C***********************************************************************
      DO 10 I=1,N
        V(I)=SCAL*V(I)
   10 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RVSCAL
CDOC BEGIN_SUBROUTINE RVSUB
CDOC Subtracts two doulble precision vectors
CDOC
CDOC This function subtracts two double precision vectors passed as
CDOC arguments and stores the result in another vector \smparm{U=V-W}.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION U      >  Double precision vector.
CDOC DOUBLE_PRECISION V      >  Double precision vector.
CDOC DOUBLE_PRECISION W      <  Double precision vector with the result
CDOC C                          \smparm{V-W}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, September 1996: Initial coding
CDOC
      SUBROUTINE RVSUB
     1(   U          ,V          ,W          ,N          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    U(N)               ,V(N)               ,W(N)
C***********************************************************************
C SUBTRACTS THE VECTOR 'W' FROM THE VECTOR 'V' AND STORE THE RESULT
C IN 'U'. U ,V AND W ARE DOUBLE PRECISION VECTORS OF DIMENSION N.
C***********************************************************************
      DO 10 I=1,N
        U(I)=V(I)-W(I)
   10 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RVSUB
CDOC BEGIN_SUBROUTINE RVZERO
CDOC Zero a double precision array
CDOC
CDOC This routine initialises to zero the \smparm{N} components of the 
CDOC double precision array argument \smparm{V}.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION V      <  Zeroed double precision array.
CDOC INTEGER          N      >  Dimension of \smparm{V}.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE RVZERO
     1(   V          ,N          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N)
      DATA R0/0.0D0/
C***********************************************************************
C INITIALISES TO ZERO A DOUBLE PRECISION ARRAY OF DIMENSION N
C***********************************************************************
      DO 10 I=1,N
        V(I)=R0
   10 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RVZERO
CDOC BEGIN_DOUBLE_PRECISION_FUNCTION SCAPRD
CDOC Scalar product of double precision vectors
CDOC
CDOC This function returns the scalar product between its two double
CDOC precision vector arguments \smparm{U}.\smparm{V}.
CDOC
CDOC BEGIN_PARAMETERS 
CDOC DOUBLE_PRECISION U      >  Array of components of a double
CDOC C                          precision vector.
CDOC DOUBLE_PRECISION V      >  Array of components of a double
CDOC C                          precision vector.
CDOC INTEGER          N      >  Dimension of \smparm{U} and \smparm{V}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, May 1996: Initial coding
CDOC
      DOUBLE PRECISION FUNCTION SCAPRD(U  ,V  ,N  ) 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(N), V(N)
      DATA  R0  / 0.0D0 /
C***********************************************************************
C SCALAR PRODUCT OF DOUBLE PRECISION VECTORS U AND V OF DIMENSION N
C***********************************************************************
      SCAPRD=R0
      DO 10 I=1,N
        SCAPRD=SCAPRD+U(I)*V(I)
  10  CONTINUE
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION SCAPRD
CDOC BEGIN_SUBROUTINE SOLQUA
CDOC Solves a quadratic equation:  a x**2 + b x + c = 0
CDOC
CDOC Given the coeficients a, b and c, this routine computes the real
CDOC roots of the associated quadratic equation, a x**2 + b x + c = 0.
CDOC The return values of the arguments \smparm{ROOT1} and
CDOC \smparm{ROOT2} (the roots) are set only if real roots exist.
CDOC If the equation admits only one real solution, the return value
CDOC of the logical argument \smparm{ONEROO} is set to \smparm{.TRUE.}
CDOC (set to \smparm{.FALSE.} otherwise).
CDOC If the equation admits two real roots, the return value
CDOC of the logical argument \smparm{TWOROO} is set to \smparm{.TRUE.}
CDOC (set to \smparm{.FALSE.} otherwise).
CDOC Consequently, if the roots are complex or no roots/infinite
CDOC number of roots exist (the equation is ill-defined), the return
CDOC value of the logical arguments \smparm{ONEROO} and
CDOC \smparm{TWOROO} is set to \smparm{.FALSE.}.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A      >  Coefficient of the quadratic term.
CDOC DOUBLE_PRECISION B      >  Coefficient of the linear term.
CDOC DOUBLE_PRECISION C      >  Coefficient of the constant term.
CDOC LOGICAL          ONEROO <  Logical flag. Set to \smparm{.TRUE.} if
CDOC C                          there is only one (real) root.
CDOC C                          Set to \smparm{.FALSE.} otherwise.
CDOC LOGICAL          TWOROO <  Logical flag. Set to \smparm{.TRUE.} if
CDOC C                          there are two distinct (real) root.
CDOC C                          Set to \smparm{.FALSE.} otherwise.
CDOC DOUBLE_PRECISION ROOT1  <  One of the roots (the only root if there
CDOC C                          is only one real root - not set if there
CDOC C                          are no real roots).
CDOC DOUBLE_PRECISION ROOT2  <  The other root (not set if there is only
CDOC C                          one real root or no real roots).
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, June 1998: Initial coding
CDOC
      SUBROUTINE SOLQUA
     1(   A          ,B          ,C          ,ONEROO     ,TWOROO     ,
     2    ROOT1      ,ROOT2      )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL  ONEROO ,TWOROO
      DATA
     1    R0   ,R1   ,R2   ,R4   ,SMALL /
     2    0.D0 ,1.0D0,2.0D0,4.0D0,1.D-12/
C***********************************************************************
C FINDS THE REAL ROOTS OF A QUADRATIC EQUATION:  A X**2 + B X + C = 0.
C
C REFERENCE:
C W.H.Press, S.A.Teukolsky, W.T.Vetterling and B.P.Flannery. Numerical
C recipes in FORTRAN. The art of scientific computing. 2nd Ed.,
C Cambridge Univ. Press, 1992. (Section 5.6)
C***********************************************************************
C Initialises logical flags
      ONEROO=.FALSE.
      TWOROO=.FALSE.
      IF(A.NE.R0)THEN
C The equation is non-linear in fact
C ----------------------------------
        IF(B.NE.R0)THEN
          SIGNB=B/ABS(B)
        ELSE
          SIGNB=R1
        ENDIF
        B2=B*B
        R4AC=R4*A*C
        SQUAR=B2-R4AC
        IF(SQUAR.GT.R0)THEN
C there are two distinct real roots: uses formula which minimises
C round-off errors when the coefficients A and/or C are small
          TWOROO=.TRUE.
          SQUAR=SQRT(SQUAR)
          Q=-(B+SIGNB*SQUAR)/R2
          ROOT1=Q/A
          ROOT2=C/Q
        ELSEIF(SQUAR.EQ.R0.OR.
     1         (SQUAR/DMAX1(B2,ABS(R4AC))+SMALL).GE.R0)THEN
C there is only one root
          ONEROO=.TRUE.
          ROOT1=-B/(R2*A)
        ENDIF
      ELSE
C The equation is linear
C ----------------------
        IF(B.NE.R0)THEN
C and well defined -> (only) one root exists
          ONEROO=.TRUE.
          ROOT1=-C/B
        ENDIF
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE SOLQUA
CDOC BEGIN_SUBROUTINE SPDEC2
CDOC Closed form spectral decomposition of 2-D symmetric tensors
CDOC
CDOC This routine performs the spectral decomposition of 2-D symmetric
CDOC tensors in closed form. The tensor is passed as argument (stored in
CDOC vector form).
CDOC 
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION EIGPRJ <  Matrix with one eigenprojection tensor
CDOC C                          of \smparm{X} stored in each column.
CDOC DOUBLE_PRECISION EIGX   <  Array containing the eigenvalues of
CDOC C                          \smparm{X}.
CDOC LOGICAL          REPEAT <  Repeated eigenvalues flag. Set to
CDOC C                          \smparm{.TRUE.} if the eigenvalues of
CDOC C                          \smparm{X} are repeated (within a small
CDOC C                          tolerance).
CDOC DOUBLE_PRECISION X      >  Array containing the components of a
CDOC C                          symmetric tensor.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, May 1996: Initial coding
CDOC
      SUBROUTINE SPDEC2
     1(   EIGPRJ     ,EIGX       ,REPEAT     ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=4    ,NDIM=2     )
      LOGICAL REPEAT
      DIMENSION
     1    EIGPRJ(MCOMP,NDIM)        ,EIGX(NDIM)                ,
     2    X(MCOMP)
      DIMENSION
     1    AUXMTX(NDIM,NDIM)         ,EIGVEC(NDIM,NDIM)
      DATA
     1    R0   ,RP5  ,R1   ,R4   ,SMALL  /
     2    0.0D0,0.5D0,1.0D0,4.0D0,1.D-5  /
C***********************************************************************
C PERFORMS THE CLOSED FORM SPECTRAL DECOMPOSITION OF A
C SYMMETRIC 2-D TENSOR STORED IN VECTOR FORM
C***********************************************************************
      REPEAT=.FALSE.
C Compute eigenvalues of X
C ------------------------
      TRX=X(1)+X(2)
      B=SQRT((X(1)-X(2))**2+R4*X(3)*X(3))
      EIGX(1)=RP5*(TRX+B)
      EIGX(2)=RP5*(TRX-B)
C Compute eigenprojection tensors
C -------------------------------
      DIFFER=ABS(EIGX(1)-EIGX(2))
      AMXEIG=DMAX1(ABS(EIGX(1)),ABS(EIGX(2)))
      IF(AMXEIG.NE.R0)DIFFER=DIFFER/AMXEIG
      IF(DIFFER.LT.SMALL)THEN
        REPEAT=.TRUE.
C for repeated (or nearly repeated) eigenvalues, re-compute eigenvalues
C and compute eigenvectors using the iterative procedure. In such cases,
C the closed formula for the eigenvectors is singular (or dominated by
C round-off errors)
        AUXMTX(1,1)=X(1)
        AUXMTX(2,2)=X(2)
        AUXMTX(1,2)=X(3)
        AUXMTX(2,1)=AUXMTX(1,2)
        CALL JACOB(AUXMTX,EIGX,EIGVEC,2)
        DO 10 IDIR=1,2
          EIGPRJ(1,IDIR)=EIGVEC(1,IDIR)*EIGVEC(1,IDIR)
          EIGPRJ(2,IDIR)=EIGVEC(2,IDIR)*EIGVEC(2,IDIR)
          EIGPRJ(3,IDIR)=EIGVEC(1,IDIR)*EIGVEC(2,IDIR)
          EIGPRJ(4,IDIR)=R0
 10     CONTINUE
      ELSE
C Use closed formula to compute eigenprojection tensors
        DO 20 IDIR=1,2
          B=EIGX(IDIR)-TRX
          C=R1/(EIGX(IDIR)+B)
          EIGPRJ(1,IDIR)=C*(X(1)+B)
          EIGPRJ(2,IDIR)=C*(X(2)+B)
          EIGPRJ(3,IDIR)=C*X(3)
          EIGPRJ(4,IDIR)=R0
 20     CONTINUE
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE SPDEC2
CDOC BEGIN_SUBROUTINE VECPRD
CDOC Vector product of two double precision vectors.
CDOC
CDOC This routine computes the vector product of two double precision
CDOC vectors passed in its argument list and store the resulting vector
CDOC in another vector also passed through its argument list.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION U      >  Array containing first vector.
CDOC DOUBLE_PRECISION V      >  Array containing second vector.
CDOC DOUBLE_PRECISION W      <  Array containing the vector product
CDOC C                          \smparm{U x V}.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, April 1998: Initial coding
CDOC
      SUBROUTINE VECPRD(U  ,V  ,W  ) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(3), V(3) ,W(3)
C***********************************************************************
C VECTOR PRODUCT OF DOUBLE PRECISION VECTORS U AND V:
C
C                            W = U x V
C***********************************************************************
      W(1)=U(2)*V(3)-U(3)*V(2)
      W(2)=U(3)*V(1)-U(1)*V(3)
      W(3)=U(1)*V(2)-U(2)*V(1)
      RETURN
      END
CDOC END_SUBROUTINE VECPRD
CDOC BEGIN_SUBROUTINE FNDKEY
CDOC Finds and reads a line containing a specified keyword from a file.
CDOC
CDOC This routine searches for a given keyword positioned at the
CDOC beginning of a line in a file.
CDOC If the given keyword is found then the corresponding line is read
CDOC and returned together with the number of words in the line and two
CDOC integer arrays containing the position of the beginning and end of
CDOC each word.
CDOC
CDOC BEGIN_PARAMETERS
CDOC LOGICAL          FOUND  <  Logical flag. Its return value is set to
CDOC C                          \smparm{.TRUE.} if the specified keyword
CDOC C                          is found at the beginning of a line in
CDOC C                          the specified file and set to
CDOC C                          \smparm{.FALSE.} otherwise.
CDOC INTEGER          IWBEG  <  Integer array containg the position of
CDOC C                          the beginning of each word in the line
CDOC C                          containing the specified keyword.
CDOC C                          \smparm{IWBEG(K)} is the position of the
CDOC C                          beginning of the \smparm{K}th word in
CDOC C                          the line. 
CDOC INTEGER          IWEND  <  Integer array containg the position of
CDOC C                          the end of each word in the line
CDOC C                          containing the specified keyword.
CDOC C                          \smparm{IWEND(K)} is the position of the
CDOC C                          end of the \smparm{K}th word in the
CDOC C                          line.
CDOC CHARACTER        KEYWRD >  Keyword.
CDOC CHARACTER        INLINE <  Contents of the line which contains the
CDOC C                          specified  keyword at the its beginning.
CDOC INTEGER          NFILE  >  Unit identifier of the file to be
CDOC C                          searched.
CDOC INTEGER          NWRD   <  Total number of words in the line
CDOC C                          containing the specified keyword.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, July 1996: Initial coding
CDOC
      SUBROUTINE FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,KEYWRD     ,INLINE     ,
     2    NFILE      ,NWRD       )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FOUND
      CHARACTER*80 INLINE
      CHARACTER*(*) KEYWRD
      DIMENSION
     1    IWBEG(40), IWEND(40)
C***********************************************************************
C FIND AND READ A LINE CONTAINING A SPECIFIED KEYWORD FROM A FILE
C***********************************************************************
 1000 FORMAT(A80)
C
      FOUND=.TRUE.
      IEND=0
   10 READ(NFILE,1000,END=20,ERR=20)INLINE
      NWRD=NWORD(INLINE,IWBEG,IWEND)
      IF(NWRD.NE.0)THEN
        IF(INLINE(IWBEG(1):IWEND(1)).EQ.KEYWRD)THEN
          GOTO 999
        ENDIF
      ENDIF
      GOTO 10
   20 IF(IEND.EQ.0)THEN
        IEND=1
        REWIND NFILE
        GOTO 10
      ELSE
        FOUND=.FALSE.
      ENDIF
  999 RETURN
      END
CDOC END_SUBROUTINE FNDKEY
