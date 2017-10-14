
      subroutine matlib3D(matData,F,stre,cc,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)

      implicit none

      integer          ELM ! do not touch

      integer          matId, nivGp, finite, isw, gp, err,
     *                 matDim, round

      double precision matData(*), F, stre, cc, iv1(*), iv2(*), dt,
     *                 eps(6)

c c
c   isw
c    1   set nivGp and exit
c    2   initialise internal variables
c    3   full stress update (compute stress, internal variables and tangent tensor)
c
      if (matDim(matId).ne.3) 
     *  call prgError(1,'matlib3d','invalid matId!')

          ! 1 2 3 4 5 6  7 8   9 10 11 12,13,14,15 16 17 18,19,20
      goto (1,2,3,4,5,6,99,8,999,99,11,12,13,14,15,16,17,18,99,99),matId

99    call prgError(2,'matlib3d','invalid matId!')

999   continue ! dummy material
      if (isw.eq.3) call prgError(2,'matlib3d','dummy material called!')
      nivGp = 0
      return

c...............................................................................................

 1    continue ! simple small strain elasticity

      nivGp = 0

      if (finite.ge.1) call prgError(3,'matlib3d',
     *      'small strain elasticity: invalid value of finite!')

      if (isw.ne.3) return

      call small_strain_elasticity (matData,F,stre,cc)

      return

c..............................................................................................

 2     continue ! Neo Hooke 

       nivGp = 0

       if (isw.ne.3) return

       call Neo_Hooke_elasticity (matData,F,stre,cc,finite)

       return

c..............................................................................................

 3    continue ! Ogden

      nivGp = 0

      if (isw.ne.3) return

      call Ogden_elasticity (matData,
     *                       matData(3),
     *                       matData(round(matData(2))+3),
     *                       F,stre,cc,finite)

      return

c..............................................................................................

 4    continue ! von Mises elasto-plasticity with isotropic hardening

      if (finite.eq.0) call prgError(1,"matlib3d",
     *   "no small strain version for von_Mises_elastoPlasticity yet!")

      if (isw.eq.1) then

        nivGp = 7
        return

      else if (isw.eq.2) then

        call unit_s_p(iv1)
        iv1(7) = 0.d0
        call pmove(iv1,iv2,7)
        return

      else
      
        call von_Mises_elastoPlasticity
     *        (matData,F,stre,cc,iv2,iv2(7),finite,err)
 
        return

      endif

c..............................................................................................

 5    continue ! fast small strain von Mises elasto-plasticity with linear isotropic hardening

      if (finite.ne.0) 
     *  call prgError(2,"matlib3d","invalid value of finite!")

      if (isw.eq.1) then

        nivGp = 7
        return

      else if (isw.eq.2) then

        call pzero(iv1,7)
        call pzero(iv2,7)
        return

      else
      
        call smallStrainTensor(eps,F)

        call small_strain_von_Mises_elastoPlasticity
     *                (iv2,iv2(7),eps,stre,cc,matData,err)
 
        return

      endif

c..............................................................................................

 6    continue ! von Mises elasto-plasticity with kinematic hardening

      call von_Mises_elastoPlasticity_with_kinematic_hardening
     *        (matData,F,stre,cc,iv1,iv2,nivGp,finite,isw,err)

      return

c..............................................................................................

 8    continue ! multiscale material 

      nivGp = 0

      if (isw.ne.3) return

c      call multiscaleMaterial3D(matData, F,stre,cc,finite,gp,err,ELM)

      return

c..............................................................................................

11    continue ! Small/Finite Strains Von Mises Isotropic Elastoplasticity (Isotropic Hardening)
               !   (Deniz)

      nivGp = 7

      if (isw.eq.1) return

c      call sf_vonmises_isotropic_elastoplasticity
c     * (matData,F,stre,cc,iv1,iv2,finite,isw, err)

    
      return
c..............................................................................................

12    continue ! Small Strains Hill Initially Anisotropic Elastoplasticity (No Hardening)
               !   (Deniz)

      nivGp = 7
      
      if (finite.eq.1) call prgError(2,'matlib3d',
     *      'small strain model: invalid value of finite!')

      if (isw.eq.1) return

c      call s_hill_anisotropic_elastoplasticity
c     * (matData,F,stre,cc,iv1,iv2,finite,isw, err)

      return

c..............................................................................................

13    continue ! testTangent

      nivGp = 0

      if (isw.ne.3) return

      call testTangent(matData,F,stre,cc,finite)

      return

c..............................................................................................

14    continue ! (E)lastic(V)isco(E)lastic(P)lastic deformations EVEP (Wulf's MSc)

      call evep (matData,F,stre,cc,iv1,iv2,dt,finite,nivGP,isw,err)

      return

c..............................................................................................

15    continue ! Hencky

      nivGp = 0

      if (isw.ne.3) return

      call Hencky_elasticity (matData,F,stre,cc,finite)

      return

c..............................................................................................

16    continue ! StVenantKirchhoff

      nivGp = 0

      if (isw.ne.3) return

      call StVenantKirchhoff_elasticity (matData,F,stre,cc,finite)

      return

c..............................................................................................

17    continue ! generic hyperelasticity

      nivGp = 0

      if (isw.ne.3) return

      call generic_hyperelasticity (matData,F,stre,cc,finite)

      return

c..............................................................................................

18    continue ! small strain anisotropic elasticity

      nivGp = 0

      if (finite.ne.0) call prgError(1,"matlib3d",
     *                    "this material works for small strains only")

      if (isw.ne.3) return

      call small_strain_anisotropic_elasticity(matData,F,stre,cc)

      return

c..............................................................................................

      end


