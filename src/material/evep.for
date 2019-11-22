
      subroutine evep(matData,F,sig,cc,iv1,iv2,dt,finite,nivGP,isw,err)
      implicit none

      double precision matData(*), F(3,*), sig(*), cc(6,*), iv1(*),
     *                 iv2(*), dt,
     *                 tol, K, detF, det_u, detFm1d3, r1ddetF,r2ddetF,
     *                 invF(3,3), betr(6), lametr(3), lame(3), fact,
     *                 H(3,3), lametr2(3), tau6(6), dtau(3,3), E(6,3),
     *                 dtaudbetr(6,6), lame2(3), press, be(6), tau(3),
     *                 r0, r1d3, r1d2, r1

      integer nivGP, isw, finite, err,
     *        i, nCmp, par(4), niv(4), i0, typ, j, j0, mxi, c,
     *        round, spectralDecomposition

      logical outp

      data r0, r1d2, r1 / 0.d0, 0.5d0, 1.d0 /

      r1d3 = r1 / 3.d0

c------------------------------------------------------------------------
c
c     nearly incompressible generalised material (see master thesis),
c
c     multiple parallel decomposition of F
c
c     formulation in principal directions
c
c------------------------------------------------------------------------

c
c   isw
c    1   set nivGp and exit
c    2   initialise internal variables
c    3   full stress update (compute stress, internal variables and tangent tensor)
c
      err = 0

      par(1) = 2
      par(2) = 5
      par(3) = 6
      par(4) = 3

      niv(1) = 0
      niv(2) = 7
      niv(3) = 7
      niv(4) = 6

      if (isw.eq.3) goto 300

c     CALCULATE nivGP

      nCmp = round(matData(2))

      nivGP = 0

      i0 = 6

      !write(*,*) nCmp

      do i=1, nCmp

        typ  = round(matData(i0))

        !write(*,*) "type = ", typ

        if (typ.lt.1.or.typ.gt.4) 
     *    call prgError(1,"evep","error in material data!")

        nivGP = nivGP + niv(typ)

        i0 = i0 + round(matData(i0+1)) * 2 + par(typ)

      enddo

      !write(*,*) nivGP

      !stop

      if (isw.eq.1) return

c     INITIALISE INTERNAL VARIABLES

      call pzero(iv1,nivGP)
      call pzero(iv2,nivGP)

      if (finite.eq.0) return

      i0 = 6
      j0 = 0

      do i=1, nCmp

        typ  = round(matData(i0))

        if (niv(typ).gt.0) then
          iv1(j0+1) = r1
          iv1(j0+2) = r1
          iv1(j0+3) = r1
          iv2(j0+1) = r1
          iv2(j0+2) = r1
          iv2(j0+3) = r1
        endif
        j0 = j0 + niv(typ)

        i0 = i0 + round(matData(i0+1)) * 2 + par(typ)

      enddo

      return

 300  continue

c     FULL STRESS UPDATE

      K    = matData(1)
      nCmp = round(matData(2))
      tol  = matData(3)
      mxi  = round(matData(4))
      outp = (round(matData(5)).ne.0)

      if (finite.eq.0) then

        !call evep_small_strain()

        return

      endif

      detF     = det_u(F)
      detFm1d3 = detF**(-r1d3)
      r1ddetF  = r1 / detF
      r2ddetF  = r1ddetF + r1ddetF

      call inverse_u(invF,F)

      call pzero(sig,6)
      call pzero(cc,36)

c...  start loop over rheological components

      i0 = 6
      j0 = 1

      do c=1, nCmp

        typ = round(matData(i0))

c...    trial elastic left Cauchy Green tensor,  be_trial = F Ci_n^{-1) F^T

        if (typ.ne.1) then
 
          call mult_u_us_p(H,F,iv2(j0))

          call mult_s_uuT_p(betr,H,F)

        else ! purely elastic ogden ("infinity") spring

          call mult_s_uuT_p(betr,F,F)

        endif

c...    spectral decomposition of be_trial

        if (spectralDecomposition(lametr2,E,betr).eq.0) then
          err = 1
          return
        endif

        lametr(1) = sqrt(lametr2(1))
        lametr(2) = sqrt(lametr2(2))
        lametr(3) = sqrt(lametr2(3))

c       initial values lame for local Newton iterations
        lame(1)   = lametr(1)
        lame(2)   = lametr(2)
        lame(3)   = lametr(3)

c...    principal Kirchhoff stresses and 
c       their derivatives with respect to lametr^2

        if     (typ.eq.1) then
c         elastic ogden spring

          call ogdenpd(matData(i0),tau,dtau,lame,detFm1d3)

        elseif (typ.eq.2) then
c         rate independent Prandtl element (elasto-plasticity)
                                         ! xi
          call rielplpd(matData(i0),tau,dtau,iv2(j0+6),lame,lametr,
     *                  detFm1d3,tol,mxi,c,outp,err)

        elseif (typ.eq.3) then
c         rate dependent Prandtl element (elasto-viscoplasticity)

          call elviplpd(matData(i0),tau,dtau,iv2(j0+6),lame,lametr,
     *                  detFm1d3,dt,tol,mxi,c,outp,err)
      
        else
c         Maxwell element (viscoelasticity)

          call vielpd  (matData(i0),tau,dtau,lame,lametr,
     *                  detFm1d3,dt,tol,mxi,c,outp,err)
        endif

        if (err.ne.0) return

        do j=1, 3
          fact = r1d2 / lametr2(j)
          do i=1, 3
            dtau(i,j) = dtau(i,j) * fact
          enddo
        enddo

c...    assemble tau and compute derivative d tau / d B

        call isotropicTensorFunction(tau6,dtaudbetr,
     *                               lametr2,tau,dtau,E,betr)

c...    transform (tau,dtau/dbetr) -> (sig,cc)

        call taub2sigccPart1(sig,cc,tau6,dtaudbetr,betr,detF)

c...    elastic left Cauchy Green tensor and Ci^{-1} for storage

        if (typ.ne.1) then

          lame2(1) = lame(1) * lame(1)
          lame2(2) = lame(2) * lame(2)
          lame2(3) = lame(3) * lame(3)

          call pzero(be,6)
          
          do i=1, 3
            do j=1, 6
              be(j) = be(j) + lame2(i) * E(j,i)
            enddo
          enddo

          call mult_u_us_p(H,invF,be)

          call mult_s_uuT_p(iv2(j0),H,invF)

        endif

        i0 = i0 + round(matData(i0+1)) * 2 + par(typ)

        j0 = j0 + niv(typ)

      enddo

c...  compute volumetric response (pressure)      

      press = r1d2 * K * (detF - r1ddetF)
      fact  = K * detF
      do i=1, 3
        sig(i)  = sig(i)  + press
        cc(i,1) = cc(i,1) + fact
        cc(i,2) = cc(i,2) + fact
        cc(i,3) = cc(i,3) + fact
      enddo

c...  add geometrical part to tangent tensor

      call taub2sigccPart2(sig,cc)

      return

      end












c=============================================================================

      subroutine ogdenpd(d,tau,dtau,lam,detFm1d3)
      implicit none
      integer          nogd, i
      double precision d(*), tau(*), dtau(3,*), lam(*), detFm1d3,
     *                 mu, alpha, fact, fact1, fact2, fact3, 
     *                 fact4, tau1, tau2, tau3,
     *                 r1d3, r2

      r2     = 2.d0
      r1d3   = 1.d0 / 3.d0

      nogd   = nint(d(2))

      call pzero(tau,3)
      call pzero(dtau,9)

      do i=1, nogd
        mu        = d(2+i)
        alpha     = d(2+nogd+i)
        !write(*,*) mu, alpha
        fact      = mu * detFm1d3**alpha
        fact1     = lam(1)**alpha
        fact2     = lam(2)**alpha
        fact3     = lam(3)**alpha
        fact4     = r1d3 * (fact1 + fact2 + fact3)
        tau1      = fact * (fact1 - fact4)
        tau2      = fact * (fact2 - fact4)
        tau3      = fact * (fact3 - fact4)
        tau(1)    = tau(1) + tau1
        tau(2)    = tau(2) + tau2
        tau(3)    = tau(3) + tau3
        fact4     = r1d3 * alpha
        tau1      = tau1 * fact4
        tau2      = tau2 * fact4
        tau3      = tau3 * fact4
        fact      = fact * r1d3 * alpha
        dtau(1,1) = dtau(1,1) + fact * fact1 * r2 - tau1
        dtau(1,2) = dtau(1,2) - fact * fact2      - tau1
        dtau(1,3) = dtau(1,3) - fact * fact3      - tau1
        dtau(2,1) = dtau(2,1) - fact * fact1      - tau2
        dtau(2,2) = dtau(2,2) + fact * fact2 * r2 - tau2
        dtau(2,3) = dtau(2,3) - fact * fact3      - tau2
        dtau(3,1) = dtau(3,1) - fact * fact1      - tau3
        dtau(3,2) = dtau(3,2) - fact * fact2      - tau3
        dtau(3,3) = dtau(3,3) + fact * fact3 * r2 - tau3
      enddo

      if (abs(tau(1)+tau(2)+tau(3)).gt.1.d-10)
     *  write(*,'(/9x,a,a/)') 'Warning in subroutine ogdenpd:',
     *                        ' |tau1 + tau2 + tau3| > 1.d-10!' 

      return
      end

c=============================================================================

      subroutine vielpd(d,tau,dtau,lame,lametr,detFm1d3,dt,
     *                  tol,mxi,c,outp,err)
      implicit none

      integer          mxi, c, err, 
     *                 i, j, iter
      double precision d(*), tau(*), dtau(3,*), lame(*), lametr(*), 
     *                 detFm1d3, dt, tol, 
     *                 lnletr(3), lnle(3), r(2), dtd2eta,
     *                 rnorm, rmax, k(3,3), invk(2,2), fact1, fact2,
     *                 eta, r1ddetk, h(3,3), detF, 
     *                 r1
      logical          outp

      r1        = 1.d0

      eta       = d(3+nint(d(2))*2)

      detF      = lametr(1) * lametr(2) * lametr(3)

      dtd2eta   = dt / (2.d0 * eta)

      lnletr(1) = log(lametr(1))
      lnletr(2) = log(lametr(2))
      lnletr(3) = log(lametr(3))

      lnle(1)   = lnletr(1)
      lnle(2)   = lnletr(2)
      lnle(3)   = lnletr(3)

      call ogdenpd(d,tau,dtau,lame,detFm1d3)

      r(1)      = lnletr(1) - lnle(1) - dtd2eta * tau(1)
      r(2)      = lnletr(2) - lnle(2) - dtd2eta * tau(2)

      rnorm = sqrt(r(1)*r(1)+r(2)*r(2))
      rmax  = max(abs(r(1)),abs(r(2)))

      if (outp) then
        write(*,'(/''  comp.'',i1,'' (Maxwell)   rnorm'')') c
        write(*,'(''                   '',g12.5)') rnorm
      endif

      iter = 0

      do while ((rnorm.gt.tol.or.rmax.gt.tol).and.iter.le.mxi) 

        iter = iter + 1

        k(1,1)    = r1 + dtd2eta * (dtau(1,1) - dtau(1,3))
        k(1,2)    =      dtd2eta * (dtau(1,2) - dtau(1,3))
        k(2,1)    =      dtd2eta * (dtau(2,1) - dtau(2,3))
        k(2,2)    = r1 + dtd2eta * (dtau(2,2) - dtau(2,3))

        r1ddetk   = r1 / (k(1,1) * k(2,2) - k(1,2) * k(2,1))

        invk(1,1) =   k(2,2) * r1ddetk
        invk(1,2) = - k(1,2) * r1ddetk
        invk(2,1) = - k(2,1) * r1ddetk
        invk(2,2) =   k(1,1) * r1ddetk

        lnle(1)   = lnle(1) + invk(1,1) * r(1) + invk(1,2) * r(2)
        lnle(2)   = lnle(2) + invk(2,1) * r(1) + invk(2,2) * r(2)

        lame(1)   = exp(lnle(1))
        lame(2)   = exp(lnle(2))
        lame(3)   = detF / (lame(1)*lame(2))

        call ogdenpd(d,tau,dtau,lame,detFm1d3)

        r(1)      = lnletr(1) - lnle(1) - dtd2eta * tau(1)
        r(2)      = lnletr(2) - lnle(2) - dtd2eta * tau(2)

        rnorm = sqrt(r(1)*r(1)+r(2)*r(2))
        rmax  = max(abs(r(1)),abs(r(2)))

        if (outp) write(*,'(''                   '',g12.5)') rnorm

      enddo

      if (iter.gt.mxi) then
        err = 1
        return
      endif

      k(1,1)    = r1 + dtd2eta * (dtau(1,1) - dtau(1,3))
      k(1,2)    =      dtd2eta * (dtau(1,2) - dtau(1,3))
      k(2,1)    =      dtd2eta * (dtau(2,1) - dtau(2,3))
      k(2,2)    = r1 + dtd2eta * (dtau(2,2) - dtau(2,3))

      r1ddetk   = r1 / (k(1,1) * k(2,2) - k(1,2) * k(2,1))

      invk(1,1) =   k(2,2) * r1ddetk
      invk(1,2) = - k(1,2) * r1ddetk
      invk(2,1) = - k(2,1) * r1ddetk
      invk(2,2) =   k(1,1) * r1ddetk
     
      fact1     = - dtd2eta * dtau(1,3)
      fact2     = - dtd2eta * dtau(2,3)

      h(1,1) = invk(1,1) * (r1 + fact1) + invk(1,2) *       fact2
      h(1,2) = invk(1,1) *       fact1  + invk(1,2) * (r1 + fact2)
      h(1,3) = invk(1,1) *       fact1  + invk(1,2) *       fact2
      h(2,1) = invk(2,1) * (r1 + fact1) + invk(2,2) *       fact2
      h(2,2) = invk(2,1) *       fact1  + invk(2,2) * (r1 + fact2)
      h(2,3) = invk(2,1) *       fact1  + invk(2,2) *       fact2
      h(3,1) = r1 - h(1,1) - h(2,1)
      h(3,2) = r1 - h(1,2) - h(2,2)
      h(3,3) = r1 - h(1,3) - h(2,3)

      call mult_u_uu_p(k,dtau,h)

      call set_u_p(dtau,k)

      return
      end

c=============================================================================

      subroutine rielplpd(d,tau,dtau,xi,lame,lametr,detFm1d3,
     *                    tol,mxi,c,outp,err)
      implicit none

      integer          mxi, c, err,
     *                 i, j, iter
      double precision d(*), tau(*), dtau(3,*), xi, lame(*), lametr(*), 
     *                 detFm1d3, tol, 
     *                 lnletr(3), lnle(3), r(2), ntau, r1dntau, 
     *                 r1dntau2, Dgdntau, Dg, dDg, v(2), w(2),
     *                 rnorm, rmax, k(2,2), invk(2,2), Phi,
     *                 r1ddetk, h(2,3), detF, n(2,2), invh(2,2),
     *                 m(2,3), fact, sigY, sigInf, HH, Dsig, mHdDsig, 
     *                 expFact,
     *                 r0, r1, r2d3, sqrt2d3
      logical          outp

      data r0, r1 / 0.d0, 1.d0 /

      r2d3      = 2.d0/3.d0
      sqrt2d3   = sqrt(r2d3)

      sigY      = d(3+nint(d(2))*2)
      HH        = d(4+nint(d(2))*2)
      sigInf    = d(5+nint(d(2))*2)

      Dsig      = sigInf - sigY
      if (Dsig .lt. 1.d-10) then
        HH = r0
        mHdDsig = r0
      else
        mHdDsig   = - HH / Dsig
      endif

      expFact   = exp(mHdDsig * xi)

      detF      = lametr(1) * lametr(2) * lametr(3)

      lnletr(1) = log(lametr(1))
      lnletr(2) = log(lametr(2))
      lnletr(3) = log(lametr(3))

      lnle(1)   = lnletr(1)
      lnle(2)   = lnletr(2)
      lnle(3)   = lnletr(3)

      call ogdenpd(d,tau,dtau,lame,detFm1d3)

      ntau      = sqrt(tau(1)*tau(1) + tau(2)*tau(2) + tau(3)*tau(3))

      Phi       = ntau - sqrt2d3 * (sigInf - Dsig*expFact)

      if (Phi.le.r0) return

      r1dntau   = r1 /  ntau
      r1dntau2  = r1dntau * r1dntau

      Dg        = r0

      rmax      = abs(Phi)
      rnorm     = rmax

      if (outp) then
        write(*,'(/''  comp.'',i1,'' (ri Prandtl)  rnorm'')') c
        write(*,'(''                     '',g12.5)') rnorm
      endif

      iter = 0

      invk(1,1) = r1
      invk(2,1) = r0
      invk(1,2) = r0
      invk(2,2) = r1

      do while ((rnorm.gt.tol.or.rmax.gt.tol).and.iter.le.mxi) 

        iter = iter + 1

        v(1)      =   tau(1) * (dtau(1,1)-dtau(1,3))
     *              + tau(2) * (dtau(2,1)-dtau(2,3))
     *              + tau(3) * (dtau(3,1)-dtau(3,3))
        v(2)      =   tau(1) * (dtau(1,2)-dtau(1,3))
     *              + tau(2) * (dtau(2,2)-dtau(2,3))
     *              + tau(3) * (dtau(3,2)-dtau(3,3))

        w(1)      =   v(1) * invk(1,1) + v(2) * invk(2,1)
        w(2)      =   v(1) * invk(1,2) + v(2) * invk(2,2)

        dDg       =  (Phi + r1dntau * (w(1)*r(1)+w(2)*r(2))) 
     *              / (r1dntau2 * (w(1)*tau(1)+w(2)*tau(2)) 
     *                   + r2d3 * HH * expFact)

        v(1)      = r(1) - tau(1) * r1dntau * dDg
        v(2)      = r(2) - tau(2) * r1dntau * dDg

        lnle(1)   = lnle(1) + invk(1,1) * v(1) + invk(1,2) * v(2)
        lnle(2)   = lnle(2) + invk(2,1) * v(1) + invk(2,2) * v(2)

        Dg        = Dg + dDg
        xi        = xi + sqrt2d3 * dDg

        expFact   = exp(mHdDsig * xi)

        lame(1)   = exp(lnle(1))
        lame(2)   = exp(lnle(2))
        lame(3)   = detF / (lame(1)*lame(2))

        call ogdenpd(d,tau,dtau,lame,detFm1d3)

        ntau      = sqrt(tau(1)*tau(1) + tau(2)*tau(2) + tau(3)*tau(3))
        r1dntau   = r1 /  ntau
        r1dntau2  = r1dntau * r1dntau

        Phi       = ntau - sqrt2d3 * (sigInf - Dsig*expFact)

        r(1)      = lnletr(1) - lnle(1) - Dg * tau(1) / ntau
        r(2)      = lnletr(2) - lnle(2) - Dg * tau(2) / ntau

        rnorm = sqrt(r(1)*r(1)+r(2)*r(2)+Phi*Phi)
        rmax  = max(abs(r(1)),max(abs(r(2)),Phi))

        if (outp) write(*,'(''                     '',g12.5)') rnorm

        Dgdntau   = Dg * r1dntau

        h(1,1)    = (r1 - tau(1) * tau(1) * r1dntau2) * Dgdntau
        h(1,2)    = (   - tau(1) * tau(2) * r1dntau2) * Dgdntau
        h(1,3)    = (   - tau(1) * tau(3) * r1dntau2) * Dgdntau
        h(2,1)    = (   - tau(2) * tau(1) * r1dntau2) * Dgdntau
        h(2,2)    = (r1 - tau(2) * tau(2) * r1dntau2) * Dgdntau
        h(2,3)    = (   - tau(2) * tau(3) * r1dntau2) * Dgdntau

        k(1,1)    =   h(1,1) * (dtau(1,1)-dtau(1,3)) 
     *              + h(1,2) * (dtau(2,1)-dtau(2,3))
     *              + h(1,3) * (dtau(3,1)-dtau(3,3)) + r1
        k(1,2)    =   h(1,1) * (dtau(1,2)-dtau(1,3)) 
     *              + h(1,2) * (dtau(2,2)-dtau(2,3))
     *              + h(1,3) * (dtau(3,2)-dtau(3,3))
        k(2,1)    =   h(2,1) * (dtau(1,1)-dtau(1,3)) 
     *              + h(2,2) * (dtau(2,1)-dtau(2,3))
     *              + h(2,3) * (dtau(3,1)-dtau(3,3))
        k(2,2)    =   h(2,1) * (dtau(1,2)-dtau(1,3)) 
     *              + h(2,2) * (dtau(2,2)-dtau(2,3))
     *              + h(2,3) * (dtau(3,2)-dtau(3,3)) + r1

        r1ddetk   = r1 / (k(1,1) * k(2,2) - k(1,2) * k(2,1))

        invk(1,1) =   k(2,2) * r1ddetk
        invk(1,2) = - k(1,2) * r1ddetk
        invk(2,1) = - k(2,1) * r1ddetk
        invk(2,2) =   k(1,1) * r1ddetk

      enddo

      if (iter.gt.mxi) then
        err = 1
        return
      endif

      k(1,1) = dtau(1,1) - dtau(1,3) 
      k(1,2) = dtau(1,2) - dtau(1,3) 
      k(2,1) = dtau(2,1) - dtau(2,3) 
      k(2,2) = dtau(2,2) - dtau(2,3) 
      
      r1ddetk   = r1 / (k(1,1) * k(2,2) - k(1,2) * k(2,1))

      m(1,1) = r1 + (k(2,2)*dtau(1,3) - k(1,2)*dtau(2,3)) * r1ddetk
      m(1,2) = m(1,1) - r1
      m(1,3) = m(1,2)
      m(2,1) = (- k(2,1)*dtau(1,3) + k(1,1)*dtau(2,3)) * r1ddetk
      m(2,2) = r1 + m(2,1)
      m(2,3) = m(2,1)

      h(1,1) = + k(2,2) * r1ddetk
     *         + r1dntau*Dg * (r1 - r1dntau2 * tau(1) * (tau(1)-tau(3)))
      h(1,2) = - k(1,2) * r1ddetk
     *         + r1dntau*Dg * (   - r1dntau2 * tau(1) * (tau(2)-tau(3)))
      h(2,1) = - k(2,1) * r1ddetk
     *         + r1dntau*Dg * (   - r1dntau2 * tau(2) * (tau(1)-tau(3)))
      h(2,2) = + k(1,1) * r1ddetk
     *         + r1dntau*Dg * (r1 - r1dntau2 * tau(2) * (tau(2)-tau(3)))

      r1ddetk = r1 / (h(1,1) * h(2,2) - h(1,2) * h(2,1))

      invh(1,1) =   h(2,2) * r1ddetk
      invh(1,2) = - h(1,2) * r1ddetk
      invh(2,1) = - h(2,1) * r1ddetk
      invh(2,2) =   h(1,1) * r1ddetk

      v(1) = (tau(1)-tau(3))*invh(1,1) + (tau(2)-tau(3))*invh(2,1)
      v(2) = (tau(1)-tau(3))*invh(1,2) + (tau(2)-tau(3))*invh(2,2)

      fact = r1dntau2 / (r1dntau2 * (v(1)*tau(1) + v(2)*tau(2))
     *                     + r2d3 * HH * expFact)

      n(1,1) = tau(1) * v(1) * fact
      n(1,2) = tau(1) * v(2) * fact
      n(2,1) = tau(2) * v(1) * fact
      n(2,2) = tau(2) * v(2) * fact

      h(1,1) = m(1,1) - n(1,1)*m(1,1)-n(1,2)*m(2,1)
      h(1,2) = m(1,2) - n(1,1)*m(1,2)-n(1,2)*m(2,2)
      h(1,3) = m(1,3) - n(1,1)*m(1,3)-n(1,2)*m(2,3)
      h(2,1) = m(2,1) - n(2,1)*m(1,1)-n(2,2)*m(2,1)
      h(2,2) = m(2,2) - n(2,1)*m(1,2)-n(2,2)*m(2,2)
      h(2,3) = m(2,3) - n(2,1)*m(1,3)-n(2,2)*m(2,3)

      dtau(1,1) = invh(1,1) * h(1,1) + invh(1,2) * h(2,1)
      dtau(1,2) = invh(1,1) * h(1,2) + invh(1,2) * h(2,2)
      dtau(1,3) = invh(1,1) * h(1,3) + invh(1,2) * h(2,3)
      dtau(2,1) = invh(2,1) * h(1,1) + invh(2,2) * h(2,1)
      dtau(2,2) = invh(2,1) * h(1,2) + invh(2,2) * h(2,2)
      dtau(2,3) = invh(2,1) * h(1,3) + invh(2,2) * h(2,3)
      dtau(3,1) = - dtau(1,1) - dtau(2,1)
      dtau(3,2) = - dtau(1,2) - dtau(2,2)
      dtau(3,3) = - dtau(1,3) - dtau(2,3)

      return
      end

c=============================================================================

      subroutine elviplpd(d,tau,dtau,xi,lame,lametr,detFm1d3,dt,
     *                    tol,mxi,c,outp,err)
      implicit none

      integer          mxi, c, err, 
     *                 i, j, iter
      double precision d(*), tau(*), dtau(3,*), xi, lame(*), lametr(*), 
     *                 detFm1d3, dt, tol, 
     *                 lnletr(3), lnle(3), r(2), ntau, r1dntau, 
     *                 r1dntau2, Dgdntau, Dg, dDg, v(2), w(2),
     *                 rnorm, rmax, k(2,2), invk(2,2), Phi,
     *                 r1ddetk, h(2,3), detF, n(2,2), invh(2,2),
     *                 m(2,3), fact, dtd2eta, eta, 
     *                 sigY, sigInf, HH, Dsig, mHdDsig, expFact,
     *                 r0, r1, r2d3, r3, sqrt2d3
      logical          outp

      data r0, r1 / 0.d0, 1.d0 /

      r2d3      = 2.d0/3.d0
      sqrt2d3   = sqrt(r2d3)

      sigY      = d(3+nint(d(2))*2)
      HH        = d(4+nint(d(2))*2)
      sigInf    = d(5+nint(d(2))*2)
      eta       = d(6+nint(d(2))*2)

      Dsig      = sigInf - sigY
      if (Dsig .lt. 1.d-10) then
        HH = r0
        mHdDsig = r0
      else
        mHdDsig   = - HH / Dsig
      endif

      expFact   = exp(mHdDsig * xi)

      dtd2eta   = dt / (eta + eta)

      detF      = lametr(1) * lametr(2) * lametr(3)

      lnletr(1) = log(lametr(1))
      lnletr(2) = log(lametr(2))
      lnletr(3) = log(lametr(3))

      lnle(1)   = lnletr(1)
      lnle(2)   = lnletr(2)
      lnle(3)   = lnletr(3)

      call ogdenpd(d,tau,dtau,lame,detFm1d3)

      ntau      = sqrt(tau(1)*tau(1) + tau(2)*tau(2) + tau(3)*tau(3))

      Phi       = ntau - sqrt2d3 * (sigInf - Dsig*expFact)

      if (Phi.le.r0) return

      r1dntau   = r1 /  ntau
      r1dntau2  = r1dntau * r1dntau

      Dg        = r0

      r3        = Phi * dtd2eta

      rmax      = abs(r3)
      rnorm     = rmax

      if (outp) then
        write(*,'(/''  comp.'',i1,'' (rd Prandtl)  rnorm'')') c
        write(*,'(''                     '',g12.5)') rnorm
      endif

      iter = 0

      invk(1,1) = r1
      invk(2,1) = r0
      invk(1,2) = r0
      invk(2,2) = r1

      do while ((rnorm.gt.tol.or.rmax.gt.tol).and.iter.le.mxi) 

        iter = iter + 1

        v(1)      =   tau(1) * (dtau(1,1)-dtau(1,3))
     *              + tau(2) * (dtau(2,1)-dtau(2,3))
     *              + tau(3) * (dtau(3,1)-dtau(3,3))
        v(2)      =   tau(1) * (dtau(1,2)-dtau(1,3))
     *              + tau(2) * (dtau(2,2)-dtau(2,3))
     *              + tau(3) * (dtau(3,2)-dtau(3,3))

        w(1)      =   v(1) * invk(1,1) + v(2) * invk(2,1)
        w(2)      =   v(1) * invk(1,2) + v(2) * invk(2,2)

        dDg       =  (r3 + dtd2eta * r1dntau * (w(1)*r(1)+w(2)*r(2))) 
     *              / (dtd2eta * r1dntau2 * (w(1)*tau(1)+w(2)*tau(2))
     *                    + r1 + dtd2eta * r2d3 * HH * expFact)

        v(1)      = r(1) - tau(1) * r1dntau * dDg
        v(2)      = r(2) - tau(2) * r1dntau * dDg

        lnle(1)   = lnle(1) + invk(1,1) * v(1) + invk(1,2) * v(2)
        lnle(2)   = lnle(2) + invk(2,1) * v(1) + invk(2,2) * v(2)

        Dg        = Dg + dDg
        xi        = xi + sqrt2d3 * dDg

        expFact   = exp(mHdDsig * xi)

        lame(1)   = exp(lnle(1))
        lame(2)   = exp(lnle(2))
        lame(3)   = detF / (lame(1)*lame(2))

        call ogdenpd(d,tau,dtau,lame,detFm1d3)

        ntau      = sqrt(tau(1)*tau(1) + tau(2)*tau(2) + tau(3)*tau(3))
        r1dntau   = r1 /  ntau
        r1dntau2  = r1dntau * r1dntau

        Phi       = ntau - sqrt2d3 * (sigInf - Dsig*expFact)

        r(1)      = lnletr(1) - lnle(1) - Dg * tau(1) / ntau
        r(2)      = lnletr(2) - lnle(2) - Dg * tau(2) / ntau
        r3        = Phi * dtd2eta - Dg

        rnorm = sqrt(r(1)*r(1)+r(2)*r(2)+r3*r3)
        rmax  = max(abs(r(1)),max(abs(r(2)),r3))

        if (outp) write(*,'(''                     '',g12.5)') rnorm

        Dgdntau   = Dg * r1dntau

        h(1,1)    = (r1 - tau(1) * tau(1) * r1dntau2) * Dgdntau
        h(1,2)    = (   - tau(1) * tau(2) * r1dntau2) * Dgdntau
        h(1,3)    = (   - tau(1) * tau(3) * r1dntau2) * Dgdntau
        h(2,1)    = (   - tau(2) * tau(1) * r1dntau2) * Dgdntau
        h(2,2)    = (r1 - tau(2) * tau(2) * r1dntau2) * Dgdntau
        h(2,3)    = (   - tau(2) * tau(3) * r1dntau2) * Dgdntau

        k(1,1)    =   h(1,1) * (dtau(1,1)-dtau(1,3)) 
     *              + h(1,2) * (dtau(2,1)-dtau(2,3))
     *              + h(1,3) * (dtau(3,1)-dtau(3,3)) + r1
        k(1,2)    =   h(1,1) * (dtau(1,2)-dtau(1,3)) 
     *              + h(1,2) * (dtau(2,2)-dtau(2,3))
     *              + h(1,3) * (dtau(3,2)-dtau(3,3))
        k(2,1)    =   h(2,1) * (dtau(1,1)-dtau(1,3)) 
     *              + h(2,2) * (dtau(2,1)-dtau(2,3))
     *              + h(2,3) * (dtau(3,1)-dtau(3,3))
        k(2,2)    =   h(2,1) * (dtau(1,2)-dtau(1,3)) 
     *              + h(2,2) * (dtau(2,2)-dtau(2,3))
     *              + h(2,3) * (dtau(3,2)-dtau(3,3)) + r1

        r1ddetk   = r1 / (k(1,1) * k(2,2) - k(1,2) * k(2,1))

        invk(1,1) =   k(2,2) * r1ddetk
        invk(1,2) = - k(1,2) * r1ddetk
        invk(2,1) = - k(2,1) * r1ddetk
        invk(2,2) =   k(1,1) * r1ddetk

      enddo

      if (iter.gt.mxi) then
        err = 1
        return
      endif

      k(1,1) = dtau(1,1) - dtau(1,3) 
      k(1,2) = dtau(1,2) - dtau(1,3) 
      k(2,1) = dtau(2,1) - dtau(2,3) 
      k(2,2) = dtau(2,2) - dtau(2,3) 
      
      r1ddetk   = r1 / (k(1,1) * k(2,2) - k(1,2) * k(2,1))

      m(1,1) = r1 + (k(2,2)*dtau(1,3) - k(1,2)*dtau(2,3)) * r1ddetk
      m(1,2) = m(1,1) - r1
      m(1,3) = m(1,2)
      m(2,1) = (- k(2,1)*dtau(1,3) + k(1,1)*dtau(2,3)) * r1ddetk
      m(2,2) = r1 + m(2,1)
      m(2,3) = m(2,1)

      h(1,1) = + k(2,2) * r1ddetk
     *         + r1dntau*Dg * (r1 - r1dntau2 * tau(1) * (tau(1)-tau(3)))
      h(1,2) = - k(1,2) * r1ddetk
     *         + r1dntau*Dg * (   - r1dntau2 * tau(1) * (tau(2)-tau(3)))
      h(2,1) = - k(2,1) * r1ddetk
     *         + r1dntau*Dg * (   - r1dntau2 * tau(2) * (tau(1)-tau(3)))
      h(2,2) = + k(1,1) * r1ddetk
     *         + r1dntau*Dg * (r1 - r1dntau2 * tau(2) * (tau(2)-tau(3)))

      r1ddetk = r1 / (h(1,1) * h(2,2) - h(1,2) * h(2,1))

      invh(1,1) =   h(2,2) * r1ddetk
      invh(1,2) = - h(1,2) * r1ddetk
      invh(2,1) = - h(2,1) * r1ddetk
      invh(2,2) =   h(1,1) * r1ddetk

      v(1) = (tau(1)-tau(3))*invh(1,1) + (tau(2)-tau(3))*invh(2,1)
      v(2) = (tau(1)-tau(3))*invh(1,2) + (tau(2)-tau(3))*invh(2,2)

      fact = r1dntau2 / ((r1dntau2 * (v(1)*tau(1) + v(2)*tau(2))
     *           + r2d3 * HH * expFact)*dtd2eta + r1) * dtd2eta

      n(1,1) = tau(1) * v(1) * fact
      n(1,2) = tau(1) * v(2) * fact
      n(2,1) = tau(2) * v(1) * fact
      n(2,2) = tau(2) * v(2) * fact

      h(1,1) = m(1,1) - n(1,1)*m(1,1)-n(1,2)*m(2,1)
      h(1,2) = m(1,2) - n(1,1)*m(1,2)-n(1,2)*m(2,2)
      h(1,3) = m(1,3) - n(1,1)*m(1,3)-n(1,2)*m(2,3)
      h(2,1) = m(2,1) - n(2,1)*m(1,1)-n(2,2)*m(2,1)
      h(2,2) = m(2,2) - n(2,1)*m(1,2)-n(2,2)*m(2,2)
      h(2,3) = m(2,3) - n(2,1)*m(1,3)-n(2,2)*m(2,3)

      dtau(1,1) = invh(1,1) * h(1,1) + invh(1,2) * h(2,1)
      dtau(1,2) = invh(1,1) * h(1,2) + invh(1,2) * h(2,2)
      dtau(1,3) = invh(1,1) * h(1,3) + invh(1,2) * h(2,3)
      dtau(2,1) = invh(2,1) * h(1,1) + invh(2,2) * h(2,1)
      dtau(2,2) = invh(2,1) * h(1,2) + invh(2,2) * h(2,2)
      dtau(2,3) = invh(2,1) * h(1,3) + invh(2,2) * h(2,3)
      dtau(3,1) = - dtau(1,1) - dtau(2,1)
      dtau(3,2) = - dtau(1,2) - dtau(2,2)
      dtau(3,3) = - dtau(1,3) - dtau(2,3)

      return
      end

c=============================================================================




