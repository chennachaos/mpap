
      subroutine von_Mises_elastoPlasticity_with_kinematic_hardening
     *            (matData,F,stre,cc,iv1,iv2,niv,finite,isw,err)

      implicit none

      integer          niv, finite, isw, err,
     *                 i, j, nh, ctrlInt(30), round

      double precision matData(*), F(3,*), stre(*), cc(6,*), iv1(*),
     *                 iv2(*),
     *                 C(6), S(6), d2SdC(6,6), fact, ctrlDbl(30),
     *                 det_u


c
c     3D MATERIAL SUBROUTINE
c
c     elasto-plasticity 
c     with isotropic and kinematic hardening 
c     --------------------------------------
c
c     isw
c      1   set nivGp and exit
c      2   initialise internal variables
c      3   full stress update (compute stress, internal variables and tangent tensor)
c

 
      return
     
      end



