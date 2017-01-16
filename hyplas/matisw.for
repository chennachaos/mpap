CDOC BEGIN_SUBROUTINE MATISW
CDOC Material interface for initialisation and switching state variables
CDOC
CDOC This routine calls the state/algorithmic variables
CDOC initialising/switching routines according to the material type.
CDOC Each material type has its own routine that initialises and
CDOC switches Gauss point data (between current and previous values).
CDOC The initialised/switched data comprises state and algorithmic
CDOC variables whose initialisation/switching rules depend on the
CDOC particular material type considered.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DETF   >  Determinant of the total deformation
CDOC C                          gradient at the current Gauss point.
CDOC C                          Used for large strain analysis only.
CDOC INTEGER          NLARGE >  Large strain flag. Large strain
CDOC C                          analysis if \smparm{NLARGE=1} and
CDOC C                          infinitesimal strain analysis otherwise.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC LOGICAL          SUFAIL <  State update failure flag. Return value
CDOC C                          set to \smparm{.FALSE.} if the state
CDOC C                          update procedure was performed
CDOC C                          successfully and set to \smparm{.TRUE.}
CDOC C                          otherwise.
CDOC DOUBLE_PRECISION EINCR  >  Array of incremental engineering strain
CDOC C                          components. Used in infinitesimal strain
CDOC C                          analysis only.
CDOC DOUBLE_PRECISION FINCR  >  Incremental deformation gradient. Used
CDOC C                          in large strain analysis only.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC LOGICAL          LALGVA <> Array of current logical algorithmic
CDOC C                          variables at the current Gauss point.
CDOC DOUBLE_PRECISION RALGVA <> Array of current real algorithmic
CDOC C                          variables at the current Gauss point.
CDOC C                          Previous converged (equilibrium) value
CDOC C                          on entry. Returns as updated value.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA <> Array of current real state variables
CDOC C                          at the current Gauss point. Previous
CDOC C                          converged (equilibrium) value on entry.
CDOC C                          Returns as updated value.
CDOC DOUBLE_PRECISION STRES  <> Array of current (Cauchy) stress
CDOC C                          components at the current Gauss point.
CDOC C                          Previous converged (equilibrium) value 
CDOC C                          on entry. Returns as updated value.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, July      1999: Initial coding
CDOC
      SUBROUTINE MATISW
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    IPROPS     ,LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,
     3    RPROPS     ,RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'hyplasmat.fi'
C Arguments
      LOGICAL   
     1     LALGVC            ,LALGVL
      DIMENSION
     1    IPROPS(*)          ,LALGVC(*)          ,LALGVL(*)          ,
     2    RALGVC(*)          ,RALGVL(*)          ,RPROPS(*)          ,
     3    RSTAVC(*)          ,RSTAVL(*)          ,STRESC(*)          ,
     4    STRESL(*)
C***********************************************************************
C MATERIAL INTERFACE FOR INITIALISATION/SWITCHING ROUTINE CALLS:
C ACCORDING TO THE MATERIAL TYPE, CALLS MATERIAL-SPECIFIC ROUTINE TO
C INITIALISE/SWITCH GAUSS POINT STATE AND ALGORITHMIC VARIABLES
C***********************************************************************
C First identify material type and class
C --------------------------------------
      MATTYP=IPROPS(2)
C	
C Then call material type-specific routines
C -----------------------------------------
      IF(MATTYP.EQ.ELASTC)THEN
C Elastic (Hencky material in large strains)
        CALL SWEL
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.TRESCA)THEN
C Tresca elasto-plastic
        CALL SWTR
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.VMISES)THEN
C von Mises elasto-plastic
        CALL SWVM
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.MOHCOU)THEN
C Mohr-Coulomb elasto-plastic
        CALL SWMC
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.DRUPRA)THEN
C Drucker-Prager elasto-plastic
        CALL SWDP
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.CAPDP)THEN
C Capped Drucker-Prager elasto-plastic
        CALL SWCADP
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.LEMDAM)THEN
C Lemaitre's ductile damage elasto-plastic model
        CALL SWDAMA
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.DAMELA)THEN
C Isotropically damaged isotropic elastic material with crack closure
C effects
        CALL SWDMEL
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.PDSCRY)THEN
C Planar double-slip single crystal
        CALL SWPDSC
     1(   MODE       ,
     2    LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.OGDEN)THEN
C Ogden hyperelasticity model
        CALL SWOGD
     1(   MODE       ,NTYPE      ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.VMMIXD)THEN
C von Mises with mixed isotropic/kinematic hardening (infinitesimal
C only)
        CALL SWVMMX
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    LALGVC     ,LALGVL     ,RALGVC     ,RALGVL     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      ELSE
C Error: Material type not recognised
        CALL ERRPRT('EI0046')
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE MATISW
