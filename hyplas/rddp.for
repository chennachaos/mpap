      SUBROUTINE RDDP
     1(   RPROPS     ,IFLAG      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    RPROPS(*)
      DATA R0   ,R1   ,R2   ,R3   ,R6   ,R9   ,R12   ,R90   ,R180   / 
     1     0.0D0,1.0D0,2.0D0,3.0D0,6.0D0,9.0D0,12.0D0,90.0D0,180.0D0/
C***********************************************************************
C ECHO MATERIAL PROPERTIES FOR DRUCKER-PRAGER TYPE
C ELASTO-PLASTIC MATERIAL WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW
C RULE AND NON-LINEAR ISOTROPIC HARDENING
C***********************************************************************
C Read and echo some of the real properties
      PHI=RPROPS(4)
      PSI=RPROPS(5)
C Check validity if some material properties
      IF(PHI.LT.R0.OR.PHI.GE.R90.OR.PSI.LT.R0.OR.PSI.GE.R90)THEN
        CALL ERRPRT('ED0115')
      ENDIF
C Echo selected approximation of Mohr-Coulomb criterion
C and set related material constants
      ROOT3=SQRT(R3)
      RADEG=ACOS(-R1)/R180
      PHIRAD=PHI*RADEG
      SINPHI=SIN(PHIRAD)
      COSPHI=COS(PHIRAD)
      TANPHI=TAN(PHIRAD)
      PSIRAD=PSI*RADEG
      SINPSI=SIN(PSIRAD)
      TANPSI=TAN(PSIRAD)
      IF(IFLAG.EQ.0)THEN
C Outer edge match with Mohr-Coulomb criterion
        DENOMA=ROOT3*(R3-SINPHI)
        DENOMB=ROOT3*(R3-SINPSI)
        ETA=R6*SINPHI/DENOMA
        XI=R6*COSPHI/DENOMA
        ETABAR=R6*SINPSI/DENOMB
      ELSEIF(IFLAG.EQ.1)THEN
C Inner edge match with Mohr-Coulomb criterion
        DENOMA=ROOT3*(R3+SINPHI)
        DENOMB=ROOT3*(R3+SINPSI)
        ETA=R6*SINPHI/DENOMA
        XI=R6*COSPHI/DENOMA
        ETABAR=R6*SINPSI/DENOMB
      ELSEIF(IFLAG.EQ.2)THEN
C Plane strain match with Mohr-Coulomb criterion
        DENOMA=SQRT(R9+R12*TANPHI**2)
        DENOMB=SQRT(R9+R12*TANPSI**2)
        ETA=R3*TANPHI/DENOMA
        XI=R3/DENOMA
        ETABAR=R3*TANPSI/DENOMB
      ELSEIF(IFLAG.EQ.3)THEN
C Match Mohr-Coulomb criterion in uniaxial compression and uniaxial
C tension
        ETA=R3*SINPHI/ROOT3
        XI=R2*COSPHI/ROOT3
        ETABAR=R3*SINPSI/ROOT3
      ELSEIF(IFLAG.EQ.4)THEN
C Match Mohr-Coulomb criterion in biaxial compression and biaxial
C tension
        ETA=R3*SINPHI/(R2*ROOT3)
        XI=R2*COSPHI/ROOT3
        ETABAR=R3*SINPSI/(R2*ROOT3)
      ELSE
        CALL ERRPRT('ED0116')
      ENDIF

      RPROPS(4)=ETA
      RPROPS(5)=XI
      RPROPS(6)=ETABAR
      RETURN
      END