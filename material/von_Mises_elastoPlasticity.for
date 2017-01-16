      subroutine von_Mises_elastoPlasticity
     *            (matData,F,stre,cc,Cp,xi,finite,err)

      implicit none

      integer          finite, err,
     *                 i, j, P(7), elTyp, mxi, iter, nip,
     *                 round

      double precision matData(*), F(3,*), stre(*), cc(6,*), Cp(*), xi,
     *                 C(6), S(6), invC(6), diC(6,6), invCp(6), 
     *                 diCp(6,6), Cpn(6), invCpn(6), xin, Phi, 
     *                 dRdCp(6,6), dRdfp(6,6), dRdDg(6),
     *                 dSdC(6,6), dSdCp(6,6), dfpdDg(6), dfpdCp(6,6),
     *                 CS(3,3), gp(3,3), fp(6), n, m1dn, r2dn, Dg,
     *                 Rnorm, R(7), K(7,7), upd(7), dndCp(6),
     *                 h21(9), h22(6), h41(6,6), h42(7,7), h43(6,7), 
     *                 h44(7,6), sigY, sigInf, H, Dsig, mHdDsig,
     *                 tol, expFact, fact, fact2, bulk, Hlin,
     *                 r2d3, sr2d3, r0, r1d2, r1, r2,
     *                 norm2_s, norm2s_u, det_u

      logical          outp

      data r0, r1d2, r1, r2 / 0.d0, 0.5d0, 1.d0, 2.d0 /
      r2d3  = r2 / 3.d0
      sr2d3 = sqrt(r2d3)

c
c     3D MATERIAL SUBROUTINE
c
c     elasto-plasticity 
c     with isotropic hardening 
c     ----------------------------------

c...  small strains .........

      if (finite.eq.0)

     *  call prgError(1,"von_Mises_elastoPlasticity",
     *                  "no small strain version implemented!")


c...  finite strains .........  
c
c    bulk, elTyp, tol, mxi, nip, outp, mu,                        sigY, H, sigInfty, Hlin
c    bulk, elTyp, tol, mxi, nip, outp, nOgd, mu1, ..., alp1, ..., sigY, H, sigInfty, Hlin
c
      bulk      = matData(1)
      elTyp     = round(matData(2))
      tol       = matData(3)
      mxi       = round(matData(4))
      nip       = round(matData(5))
      outp      = (round(matData(6)).ne.0)

      i = 7
      if (elTyp.eq.4) i = i + round(matData(7))*2

      sigY      = matData(i+1)
      H         = matData(i+2)
      sigInf    = matData(i+3)
      Hlin      = matData(i+4)

      if (nip.lt.0)      nip =  4
      if (mxi.le.0)      mxi = 10
      if (tol.lt.1.d-14) tol = 1.d-14

      Dsig      = sigInf - sigY
      if (Dsig .lt. 1.d-10) then
        H = r0
        mHdDsig = r0
      else
        mHdDsig = - H / Dsig
      endif

      expFact   = exp(mHdDsig * xi)

c...  stress update  -> S, Cp, xi

      call mult_s_uTu_p(C,F,F)
      call inverse_s (invC,C)
      call diff_ATA_m(diC,invC)      

c     compute trial state

      call inverse_s (invCp,Cp)                  ! -> invCp
      call diff_ATA_m(diCp,invCp)                ! -> diCp
      call computeElasticResponse(S,dSdC,dSdCp,  !
     *         C,invC,diC,Cp,invCp,diCp,         !
     *         bulk,matData(7),elTyp)            ! -> S, dSdC, dSdCp    

      call mult_u_ss_p(CS,C,S)                   ! -> CS
      call set_u_p    (gp,CS)                    !
      call dev_u_r    (gp)                       ! -> gp
      n    = norm2s_u(gp)                        ! -> n
      m1dn = - r1 / n                            ! -> m1dn
      r2dn = - m1dn - m1dn                       ! -> r2dn
      call mult_s_us_f(fp,gp,Cp,r2dn)            ! -> fp

      Phi  = n - sr2d3 * (sigInf - Dsig*expFact + Hlin*xi) !

c      Phi  = n - sr2d3 * (sigY + Dsig*(1.0-expFact) + Hlin*xi) !


      !write(*,*) Phi

      Dg = r0

      if (Phi.gt.tol) then

        call pmove(Cp,    Cpn,    6   )
        call pmove(invCp, invCpn, 6   )
        xin = xi

        if (outp) write(*,'(/5x,g12.5)') abs(Phi)

        iter = 0
        do while ((Rnorm.gt.tol.or.Phi.gt.tol).and.iter.le.mxi)

          iter = iter + 1

c         partial derivatives dfp/dCp, dfp/dDg, dn/dCp, dn/dDg 
          call diff_norm2sATd_p(h21,C,CS,n)             ! now: h21 = d n / d S
          call mult_s_s4_p     (dndCp,h21,dSdCp)        ! now: dndCp = d n / d Cp
          call diff_UT_p       (h41,gp)                 ! now: h41 = d dev(CS)T / d T
          call diff_ATdB_p     (h42,C,Cp)               ! now: h42 = d dev(CS)Cp / d S
          call mult_4_44_ap    (h41,h42,dSdCp)          ! now: h41 = d dev(CS)Cp / d Cp
          call set_4_f         (dfpdCp,h41,r2dn)        ! 
          call mult_4_ss_af    (dfpdCp,fp,dndCp,m1dn)   ! now: dfpdCp = d fp / d Cp

c         R, dR/dCp, dR/dfp, dR/dDg                 
          call comp_R1111(R,dRdCp,dRdfp,dRdDg,                !
     *                    fp,Cp,invCp,diCp,Cpn,invCpn,Dg,nip) ! -> R
          call mult_4_44_ap(dRdCp,dRdfp,dfpdCp)               !

          do i=1, 6              
            do j=1, 6            
              K(i,j) = dRdCp(i,j)
            enddo                
            K(i,7) = dRdDg(i)
            K(7,i) = dndCp(i)    
          enddo
          K(7,7) = - r2d3 * (H * expFact + Hlin)

c         solve and update Cp, Dg
          R(7) = Phi
          call decompLR_matrix(K,P,7,0)
          call solve_matrix   (K,P,R,upd,7)
          Dg = Dg - upd(7)
          xi = xin + sr2d3 * Dg
          expFact = exp(mHdDsig * xi)
          do i=4, 6
            upd(i) = r1d2 * upd(i)
          enddo
          call set_s_am(Cp,upd)

c         update S, fp, Phi
          call inverse_s (invCp,Cp)                  ! -> invCp
          call diff_ATA_m(diCp,invCp)                ! -> diCp
          call computeElasticResponse(S,dSdC,dSdCp,  !
     *         C,invC,diC,Cp,invCp,diCp,             !
     *         bulk,matData(7),elTyp)                ! -> S, dSdC, dSdCp    

          call mult_u_ss_p(CS,C,S)                   ! -> CS
          call set_u_p    (gp,CS)                    !
          call dev_u_r    (gp)                       ! -> gp
          n    = norm2s_u(gp)                        ! -> n
          m1dn = - r1 / n                            ! -> m1dn
          r2dn = - m1dn - m1dn                       ! -> r2dn
          call mult_s_us_f(fp,gp,Cp,r2dn)            ! -> fp

          Phi  = n - sr2d3 * (sigInf - Dsig*expFact + Hlin*xi) ! -> Phi

c          Phi  = n - sr2d3 * (sigY + Dsig*(1.0-expFact) + Hlin*xi) !

          call comp_R1000(R,fp,Cp,invCp,Cpn,invCpn,Dg,nip)

          Rnorm = norm2_s(R)

          if (outp) write(*,'(5x,g12.5)') sqrt(Rnorm*Rnorm + Phi*Phi)

        enddo

        if (iter.gt.mxi) then
          write(*,'(/a/)')' strangely, it did not converge (Cp,Dg)!'
          err = 1
          return
        endif

      endif

c...  consistent tangent = dS/dC

      if (Dg.gt.tol) then

c       dR/dC, dPhi/dC
        call diff_norm2sATd_p (h21,C,CS,n)         ! now: h21 = d n / d S
        call mult_s_s4_p      (h22,h21,dSdC)       ! now: h22 = d n / d S * d S / d C
        call transpose_r      (CS)                 !
        call diff_norm2sATd_ap(h22,S,CS,n)         ! now: h22 = d n / d C
        call diff_ATdB_f      (h41,C,Cp,r2dn)      !
        call mult_4_44_p      (h42,h41,dSdC)       !
        call diff_TAdB_af     (h42,S,Cp,r2dn)      !
        call mult_4_ss_af     (h42,fp,h22,m1dn)    ! now: h42 = d fp / d C
        call mult_4_44_p      (h41,dRdfp,h42)      ! now: h41 = d R / d C
                                                   !
        call inverse_matrix(K,h42,P,7)             ! now: h42 = K^{-1}
                                                   !
        do i=1, 6
          do j=1, 6
            h43(i,j) = h42(i,j)
            h44(i,j) = h41(i,j) 
          enddo
          h43(i,7) = h42(i,7)  !now: h43 = d Cp / d(R,Phi)   -> 6x7-matrix 
          h44(7,i) = h22(i)    !now: h44 = d(R,Phi) / d C    -> 7x6-matrix  
        enddo

        call mult_matrix_m (h41, h43,  h44,6,7,6)
        call mult_matrix_ap(dSdC,dSdCp,h41,6,6,6)

      endif

c...  push forward of stre and tangent tensors

      call pushfwrd_s(stre,F,S)

      call pushfwrd_4(cc,F,dSdC)

      fact  = r1 / det_u(F)
      fact2 = fact + fact

      do i=1, 6
        stre(i) = stre(i) * fact
        do j=1, 6
          cc(i,j) = cc(i,j) * fact2
        enddo
      enddo
 
      return
     
 100  continue

      end


