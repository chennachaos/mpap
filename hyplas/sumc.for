CDOC BEGIN_SUBROUTINE SUMC
CDOC State update procedure for the Mohr-Coulomb type material model.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the
CDOC Mohr-Coulomb elasto-plastic material model with piece-wise linear
CDOC isotropic hardening. It contains the plane strain and axisymmetric
CDOC implementations of the model. The essential return mapping is
CDOC carried out here in principal stress space.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   <  Array of incremental plastic
CDOC C                          multipliers.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines
CDOC C                          \smparm{INDATA} and \smparm{RDMC}.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic yielding
CDOC C                          flag, \smparm{IFPLAS}; the return
CDOC C                          algorithm failure flag, \smparm{SUFAIL};
CDOC C                          the edge return flag, \smparm{EDGE};
CDOC C                          the right edge return flag,
CDOC C                          \smparm{RIGHT}, and; the apex return
CDOC C                          flag, \smparm{APEX}.
CDOC C                          The plastic yielding flag is set
CDOC C                          to \smparm{.TRUE.} if plastic yielding
CDOC C                          has occurred and to \smparm{.FALSE.} if
CDOC C                          the step is elastic. The algorithm
CDOC C                          failure flag is set to {\tt .FALSE.} if
CDOC C                          the state update algorithm has been
CDOC C                          successful and to \smparm{.TRUE.} if the
CDOC C                          return mapping algorithm has failed
CDOC C                          to converge.
CDOC C                          \smparm{EDGE} is set to \smparm{.TRUE.}
CDOC C                          if the selected return mapping is to
CDOC C                          an edge (right or left) and is set to
CDOC C                          \smparm{.FALSE.} otherwise.
CDOC C                          \smparm{RIGHT} is set to \smparm{.TRUE.}
CDOC C                          if the selected return is to the right
CDOC C                          edge and is set to \smparm{.FALSE.}
CDOC C                          otherwise. \smparm{APEX} is set to
CDOC C                          \smparm{.TRUE.} if the selected return
CDOC C                          is to the apex and is \smparm{.FALSE.}
CDOC C                          otherwise.
CDOC INTEGER          NTYPE  >  Stress state type. This routine is
CDOC C                          compatible only with plane strain
CDOC C                          \smparm{NTYPE=2} and axisymmetric
CDOC C                          \smparm{NTYPE=3} states.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, and the plastic
CDOC C                          properties: Sine and cosine of
CDOC C                          frictional angle, cosine of the
CDOC C                          dilatancy angle and the pairs
CDOC C                          ``accumulated plastic strain-cohesion''
CDOC C                          defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          \smparm{RDMC}.
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
CDOC E.de Souza Neto and P.H.Saksono, July 1996: Initial coding
CDOC
      SUBROUTINE SUMC
     1(   DGAM       ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=7  ,MSTRE=4)
C Arguments
      LOGICAL
     1    LALGVA(5) 
      DIMENSION
     1    DGAM(2)            ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)    ,STRAN(MSTRE)       ,STRES(MSTRE)
C Local variables and arrays
      LOGICAL
     1    APEX, DUMMY, EDGE, IFPLAS, RIGHT, SUFAIL
      DIMENSION
     1    EIGPRJ(MSTRE,2)    ,PSTRS(3)           ,STREST(3)
      DATA
     1    R0   ,R1   ,R2   ,R3   ,R4   ,SMALL ,TOL   / 
     2    0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,1.D-06,1.D-10/
      DATA MXITER / 50 /
C***********************************************************************
C STATE UPDATE PROCEDURE FOR MOHR-COULOMB TYPE ELASTO-PLASTIC MATERIAL
C WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE AND PIECE-WISE LINEAR
C ISOTROPIC HARDENING:
C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM (BOXES 8.4-7).
C PLANE STRAIN AND AXISYMMETRIC IMPLMENTATIONS.
C***********************************************************************
C Stops program if neither plane strain nor plane stress state
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0027')
C Initialize some algorithmic and internal variables
      DGAMA=R0
      DGAMB=R0
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
      EDGE=.FALSE.
      APEX=.FALSE.
      EPBARN=RSTAVA(MSTRE+1)
      EPBAR=EPBARN
C Set some material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      SINPHI=RPROPS(4)
      COSPHI=RPROPS(5)
      SINPSI=RPROPS(6)
      NHARD=IPROPS(3)
C Set some constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R4G=R4*GMODU
      R2BULK=R2*BULK
      R2CPHI=R2*COSPHI
      R1D3=R1/R3
      ROO2D3=SQRT(R2*R1D3)
      ROO4D3=SQRT(R4*R1D3)
C Compute elastic trial state
C ---------------------------
C Elastic trial volumetric strain and pressure stress
      EETV=STRAN(1)+STRAN(2)+STRAN(4)
      PT=BULK*EETV
C Spectral decomposition of the elastic trial stress
      EETVD3=EETV*R1D3
      STREST(1)=R2G*(STRAN(1)-EETVD3)+PT
      STREST(2)=R2G*(STRAN(2)-EETVD3)+PT
      STREST(3)=GMODU*STRAN(3)
      CALL SPDEC2(EIGPRJ,PSTRS,DUMMY,STREST)
      PSTRS(3)=R2G*(STRAN(4)-EETVD3)+PT
C Identify maximum (PSTRS1) and minimum (PSTRS3) principal stresses
      II=1
      JJ=1
      PSTRS1=PSTRS(II)
      PSTRS3=PSTRS(JJ)
      DO 10 I=2,3
        IF(PSTRS(I).GT.PSTRS1)THEN
          II=I
          PSTRS1=PSTRS(II)
        ENDIF
        IF(PSTRS(I).LT.PSTRS3)THEN
          JJ=I
          PSTRS3=PSTRS(JJ)
        ENDIF
   10 CONTINUE
      IF(II.NE.1.AND.JJ.NE.1)MM=1
      IF(II.NE.2.AND.JJ.NE.2)MM=2
      IF(II.NE.3.AND.JJ.NE.3)MM=3
      PSTRS2=PSTRS(MM)
C Compute trial yield function and check for plastic consistency
C --------------------------------------------------------------
      COHE=PLFUN(EPBARN,NHARD,RPROPS(IPHARD))
      SMCT=PSTRS1-PSTRS3+(PSTRS1+PSTRS3)*SINPHI
      PHIA=SMCT-R2CPHI*COHE
      RES=PHIA
      IF(COHE.NE.R0)RES=RES/ABS(COHE)
      IF(RES.GT.TOL)THEN
C Plastic step: Apply return mapping
C ==================================
        IFPLAS=.TRUE.
C identify possible edge return: either right or left of main plane
        SCAPRD=PSTRS1*(R1-SINPSI)+PSTRS2*(-R2)+PSTRS3*(R1+SINPSI)
        IF(SCAPRD.GE.R0)THEN
          RIGHT=.TRUE.
        ELSE
          RIGHT=.FALSE.
        ENDIF
C Apply one-vector return mapping first (return to MAIN PLANE)
C ------------------------------------------------------------
        SPHSPS=SINPHI*SINPSI
        CONSTA=R4G*(R1+R1D3*SPHSPS)+R4*BULK*SPHSPS
        DEPFAC=ROO4D3*SQRT(R1+SINPSI**2)
C Start Newton-Raphson iterations for DGAMA
        DO 20 NRITER=1,MXITER
C Compute residual derivative
          DENOM=-CONSTA-R2CPHI*DEPFAC*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
C Compute Newton-Raphson increment and update variable DGAMA
          DDGAMA=-PHIA/DENOM
          DGAMA=DGAMA+DDGAMA
C Compute new residual
          EPBAR=EPBARN+DEPFAC*DGAMA
          COHE=PLFUN(EPBAR,NHARD,RPROPS(IPHARD)) 
          PHIA=SMCT-CONSTA*DGAMA-R2CPHI*COHE
C Check convergence
          RESNOR=ABS(PHIA)
          IF(SMCT.NE.R0)RESNOR=RESNOR/ABS(SMCT)
          IF(RESNOR.LE.TOL)THEN
C Check validity of 1-vector return (check sextant of converged stress)
            S1=PSTRS1-(R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI)*DGAMA
            S2=PSTRS2+(R4G*R1D3-R2BULK)*SINPSI*DGAMA
            S3=PSTRS3+(R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI)*DGAMA
            DELTA=DMAX1(ABS(S1),ABS(S2),ABS(S3))*SMALL
            IF(S1+DELTA.GE.S2.AND.S2+DELTA.GE.S3)THEN
C converged stress is in the same sextant as trial stress -> 1-vector
C return is valid.
              P=(S1+S2+S3)*R1D3
              GOTO 70
            ELSE
C converged stress is not in the same sextant -> 1-vector result is
C not valid. Go to two-vector return map to edge
              GOTO 30
            ENDIF
          ENDIF
   20   CONTINUE
C failure of stress update procedure
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0003')
        GOTO 999
   30   CONTINUE
C Apply two-vector return mapping to appropriate EDGE
C ---------------------------------------------------
        DGAMA=R0
        DGBAR=R1
        EPBAR=EPBARN
        COHE=PLFUN(EPBARN,NHARD,RPROPS(IPHARD))
        SMCTA=PSTRS1-PSTRS3+(PSTRS1+PSTRS3)*SINPHI
        IF(RIGHT)THEN
          SMCTB=PSTRS1-PSTRS2+(PSTRS1+PSTRS2)*SINPHI
        ELSE
          SMCTB=PSTRS2-PSTRS3+(PSTRS2+PSTRS3)*SINPHI
        ENDIF
        PHIA=SMCTA-R2CPHI*COHE
        PHIB=SMCTB-R2CPHI*COHE
        IF(RIGHT)THEN
          CONSTB=R2G*(R1+SINPHI+SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS
          CONSTG=(R1-SINPSI)**2
          CONSTH=(R1+SINPSI)**2
        ELSE
          CONSTB=R2G*(R1-SINPHI-SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS
          CONSTG=(R1+SINPSI)**2
          CONSTH=(R1-SINPSI)**2
        ENDIF
C Start Newton-Raphson iterations for DGAMA and DGAMB
        DO 40 NRITER=1,MXITER
C Compute residual derivative matrix
          FACTA=R2CPHI*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))*ROO2D3/DGBAR*
     1          (DGAMA*CONSTG+(DGAMA+DGAMB)*CONSTH)
          FACTB=R2CPHI*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))*ROO2D3/DGBAR*
     1          (DGAMB*CONSTG+(DGAMA+DGAMB)*CONSTH)
          DRVAA=-CONSTA-FACTA
          DRVAB=-CONSTB-FACTB
          DRVBA=-CONSTB-FACTA
          DRVBB=-CONSTA-FACTB
C Compute Newton-Raphson increment and update variables DGAMA and DGAMB
          R1DDET=R1/(DRVAA*DRVBB-DRVAB*DRVBA)
          DDGAMA=(-DRVBB*PHIA+DRVAB*PHIB)*R1DDET
          DDGAMB=(DRVBA*PHIA-DRVAA*PHIB)*R1DDET
          DGAMA=DGAMA+DDGAMA
          DGAMB=DGAMB+DDGAMB
C Compute new residual
          DGBAR=SQRT(CONSTG*(DGAMA**2+DGAMB**2)+CONSTH*(DGAMA+DGAMB)**2)
          EPBAR=EPBARN+ROO2D3*DGBAR
          COHE=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          PHIA=SMCTA-CONSTA*DGAMA-CONSTB*DGAMB-R2CPHI*COHE
          PHIB=SMCTB-CONSTB*DGAMA-CONSTA*DGAMB-R2CPHI*COHE
C Check convergence
          RESNOR=(ABS(PHIA)+ABS(PHIB))
          FACTOR=(ABS(SMCTA)+ABS(SMCTB))
          IF(FACTOR.NE.R0)RESNOR=RESNOR/FACTOR
          IF(RESNOR.LE.TOL)THEN
C Check validity of 2-vector return to edge
            AUX1=R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI
            AUX2=(R4G*R1D3-R2BULK)*SINPSI
            AUX3=R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI
            IF(RIGHT)THEN
              S1=PSTRS1-AUX1*(DGAMA+DGAMB)
              S2=PSTRS2+AUX2*DGAMA+AUX3*DGAMB
              S3=PSTRS3+AUX3*DGAMA+AUX2*DGAMB
            ELSE
              S1=PSTRS1-AUX1*DGAMA+AUX2*DGAMB
              S2=PSTRS2+AUX2*DGAMA-AUX1*DGAMB
              S3=PSTRS3+AUX3*(DGAMA+DGAMB)
            ENDIF
            DELTA=DMAX1(ABS(S1),ABS(S2),ABS(S3))*SMALL
            IF(S1+DELTA.GE.S2.AND.S2+DELTA.GE.S3)THEN
C converged stress is in the same sextant as trial stress -> 2-vector
C return to edge is valid.
              EDGE=.TRUE.
              P=(S1+S2+S3)*R1D3
              GOTO 70
            ELSE
C converged stress is not in the same sextant -> 2-vector return to edge
C is not valid. Go to two-vector return map to APEX
              GOTO 50
            ENDIF
          ENDIF
   40   CONTINUE
C failure of stress update procedure
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0003')
        GOTO 999
   50   CONTINUE
C Apply two-vector return mapping to APEX
C ---------------------------------------
C First set DGAMA and initial guess for DGAMB
        DGAMB=R0
        DGAMA=R1/R2G
        STRLNR=SQRT((PSTRS1-PT)**2+(PSTRS2-PT)**2+(PSTRS3-PT)**2)
        ETDNOR=DGAMA*STRLNR
        IF(ETDNOR.EQ.R0)ROO1D3=SQRT(R1D3)
        DGBAR=ETDNOR
        EPBAR=EPBARN+ROO2D3*DGBAR
        COHE=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
        COTPHI=COSPHI/SINPHI
        RES=COTPHI*COHE-PT
C Newton-Raphson iterations for DGAMB
        DO 60 NRITER=1,MXITER
          IF(ETDNOR.EQ.R0)THEN
            DENOM=COTPHI*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))*
     1            ROO2D3*ROO1D3+BULK
          ELSE
            DENOM=COTPHI*DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))*
     1            ROO2D3/DGBAR*R1D3*DGAMB+BULK
          ENDIF
          DDGAMB=-RES/DENOM
          DGAMB=DGAMB+DDGAMB
          DGBAR=SQRT(ETDNOR**2+R1D3*DGAMB**2)
          EPBAR=EPBARN+ROO2D3*DGBAR
          COHE=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          P=PT-BULK*DGAMB
          RES=COTPHI*COHE-P
C check for convergence
          RESNOR=ABS(RES)
          IF(PT.NE.R0)RESNOR=RESNOR/ABS(PT)
          IF(RESNOR.LE.TOL)THEN
            APEX=.TRUE.
C update principal stresses
            S1=P
            S2=P
            S3=P
            GOTO 70
          ENDIF
   60   CONTINUE
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0003')
        GOTO 999
   70   CONTINUE
C update internal variable EPBAR  and stress components
C -----------------------------------------------------
        RSTAVA(MSTRE+1)=EPBAR
        PSTRS(II)=S1
        PSTRS(JJ)=S3
        PSTRS(MM)=S2
        STRES(1)=PSTRS(1)*EIGPRJ(1,1)+PSTRS(2)*EIGPRJ(1,2)
        STRES(2)=PSTRS(1)*EIGPRJ(2,1)+PSTRS(2)*EIGPRJ(2,2)
        STRES(3)=PSTRS(1)*EIGPRJ(3,1)+PSTRS(2)*EIGPRJ(3,2)
        STRES(4)=PSTRS(3)
C and elastic engineering strain
        EEVD3=P/BULK*R1D3
        RSTAVA(1)=(STRES(1)-P)/R2G+EEVD3
        RSTAVA(2)=(STRES(2)-P)/R2G+EEVD3
        RSTAVA(3)=STRES(3)/GMODU
        RSTAVA(4)=(STRES(4)-P)/R2G+EEVD3
      ELSE
C Elastic step: update stress using linear elastic law
C ====================================================
        STRES(1)=STREST(1)
        STRES(2)=STREST(2)
        STRES(3)=STREST(3)
        STRES(4)=PSTRS(3)
C elastic engineering strain
        RSTAVA(1)=STRAN(1)
        RSTAVA(2)=STRAN(2)
        RSTAVA(3)=STRAN(3)
        RSTAVA(4)=STRAN(4)
      ENDIF
  999 CONTINUE
C Update algorithmic variables before exit
C ========================================
      DGAM(1)=DGAMA
      DGAM(2)=DGAMB
      LALGVA(1)=IFPLAS
      LALGVA(2)=SUFAIL
      LALGVA(3)=EDGE
      LALGVA(4)=RIGHT
      LALGVA(5)=APEX
      RETURN
      END
CDOC END_SUBROUTINE SUMC