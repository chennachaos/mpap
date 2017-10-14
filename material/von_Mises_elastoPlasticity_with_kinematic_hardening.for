
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

c...  small strains .........

      if (finite.eq.0)

     *  call prgError(1,
     *    "von_Mises_elastoPlasticity_with_kinematic_hardening",
     *    "small strain version implemented not implemented!")


c...  finite strains .........  
c
c   K, nHard, elTyp, tol, mxi, nip, outp
c   mu  (or)  nOgd, mu1, mu2, ..., alpha1, alpha2, ...
c   sigY, H, sigInf
c   muH, kH, nipH
c   muH, kH, nipH
c   muH, kH, nipH
c

      nh  = round(matData(2))

      niv = (nh+1) * 6 + 1

      !write(*,*) niv

      if (isw.eq.1) then

        return

      else if (isw.eq.2) then

        do i=1, nh+1
          call unit_s_p(iv2(1+(i-1)*6))
        enddo
        iv2(niv) = 0.d0
        call pmove(iv2,iv1,niv)
        return

      else if (isw.ne.3) then

        call prgError(2,
     *  "von_Mises_elastoPlasticity_with_kinematic_hardening",
     *  "invalid value of isw!")

      endif

      call mult_s_uTu_p(C,F,F)

               !    Cp     Cpi     xi
      call nlkh_lC(iv2(1),iv2(7),iv2(niv),C,S,d2SdC,matData,
     *             ctrlDbl,ctrlInt)

c...  push forward of stre and tangent tensors

      call pushfwrd_s(stre,F,S)

      call pushfwrd_4(cc,F,d2SdC)

      fact  = 1.d0 / det_u(F)

      do i=1, 6
        stre(i) = stre(i) * fact
        do j=1, 6
          cc(i,j) = cc(i,j) * fact
        enddo
      enddo
 
      return
     
      end






c*****************************************************************************
c*****************************************************************************


      subroutine nlkh_lC(Cp,Cpi,xi,C,S,d2SdC,matData,ctlr,ctli)
c
c           n on
c           l inear
c           k inematic
c           h ardening
c
c           l arge strains
c     model C
c
c=============================================================================
c
c      +--------------------------------------------------------------+
c      |                                                              |
c      |  NONLINEAR KINEMATIC HARDENING IN FINITE STRAINS (MODEL C)   |
c      |      WITH SEVERAL HARDENING COMPONENTS IN PARALLEL AND       |
c      |             WITH NONLINEAR ISOTROPIC HARDENING               |
c      |                                                              |
c      +--------------------------------------------------------------+
c      |           UPDATE OF STRESS AND INTERNAL VARIABLES            |
c      |             & COMPUTATION OF CONSISTENT TANGENT              |
c      +--------------------------------------------------------------+
c      | USING THE BACKWARD EULER EXPONENTIAL MAP FOR BOTH Cp AND Cpi,|
c      |         LOCAL SYSTEM OF EQUATIONS:  7x7 (+ nh * 6x6)         |
c      +--------------------------------------------------------------+
c      |                    Wulf Dettmer, May 2003                    |
c      |                    modified, October 2007                    |
c      +--------------------------------------------------------------+
c
c=============================================================================
c
c  INPUT:  Cp(6)      - plastic right Cauchy-Green tensor    at t_n
c          Cpi(6,nh)  - plastic inelastic right C.-G. tensor at t_n
c          xi         - accumulated plastic strain           at t_n
c          C(6)       - right Cauchy-Green tensor            at t_{n+1}
c          matData(*) - material parameters, control data
c
c=============================================================================
c
c  OUTPUT: Cp(6)      - plastic right Cauchy-Green tensor    at t_{n+1}
c          Cpi(6,nh)  - plastic inelastic right C.-G. tensor at t_{n+1}
c          xi         - accumulated plastic strain           at t_{n+1}
c          S(6)       - 2. Piola-Kirchhoff stress tensor     at t_{n+1}
c          d2SdC(6,6) - 4. order tangent tensor   2 dS/dC    at t_{n+1}
c          ctlr(nh+1) - double precision control variables
c          clti(nh+1) - integer control variables
c
c             ctlr - determinants of plastic strain tensors
c                         ctlr(1) - det Cp(t_{n+1})
c                         ctlr(2) - det Cpi_1(t_{n+1})
c                         ctlr(3) - det Cpi_2(t_{n+1})
c                         ctlr(.) - ... ..............
c
c             ctli - number of iterations required to achieve tol
c                         ctli(1) - Cp(t_{n+1})
c                         ctli(2) - Cpi_1(t_{n+1})  (max_iter)
c                         ctli(3) - Cpi_2(t_{n+1})  (max_iter)
c                         ctli(.) - ..............  ..........
c
c=============================================================================
c
      implicit none
      integer          ctli(*), 
     *                 i, j, i0, l, nh, P(7), maxh, mxi, 
     *                 iter1, iter2, nipCp, nipCpi, elTyp, round
      parameter        (maxh = 9)

      double precision Cp(6), Cpi(6,*), xi, C(6), S(6), d2SdC(6,6), 
     *                 matData(*), ctlr(*),

     *                 invC(6),        diC(6,6), 
     *                 invCp(6),       diCp(6,6), 
     *                 invCpi(6,maxh), diCpi(6,6,maxh), 
     *                 Cpn(6), invCpn(6), Cpin(6,maxh), invCpin(6,maxh), 
     *                 xin,

     *                 Rp(6), Rpi(6,maxh), Phi, 
   
     *                 dRdA(6,6), dRdf(6,6), dRdDg(6),

     *                 X(6,maxh), Xt(6), 

     *                 dSdC(6,6), dSdCp(6,6),
     *                 dXdCp(6,6,maxh), dXdCpi(6,6,maxh), dXtdCp(6,6), 
     *                 
     *                 CS(3,3), mCpXt(3,3), 
     *                 gp(3,3), fp(6), n, m1dn, r2dn, Dg, 
     *                 expFact, sigY, sigInf, Dsig, H, mHdDsig,

     *                 dfpidX(6,6), dCpidCp(6,6,maxh), dCpidDg(6,maxh),

     *                 R(7), K(7,7), upd(7), 

     *                 h21(9), h22(6), h23(6), h24(6), h25(6), 
     *                 h41(6,6), h42(7,7), h43(6,7), h44(7,6), 
     *                 fact,

     *                 tol, r2dk, r2d3, sr2d3, r0, r1, r1d2, r2,
     *                 norm2_s, normx_s, det_s

      logical          outR, outW, noconv

      data r0, r1d2, r1, r2 / 0.d0, 0.5d0, 1.d0, 2.d0 /
      r2d3  = r2 / 3.d0
      sr2d3 = sqrt(r2d3)

      nh      = round(matData(2))
      elTyp   = round(matData(3))
      tol     = matData(4)
      mxi     = round(matData(5))
      nipCp   = round(matData(6))
      i       = round(matData(7))
      outR    = (i.eq.1.or.i.eq.3)
      outW    = (i.eq.2.or.i.eq.3)
      i0 = 9                                        ! matData(i0,...) holds isotropic
      if (elTyp.eq.4) i0 = i0 + round(matData(8))*2 ! hardening properties
      sigY    = matData(i0)
      H       = matData(i0+1)
      sigInf  = matData(i0+2)
                                                    ! matData(i0+j+j+j,...) holds
                                                    ! hardening properties of comp. j
      if (mxi.le.0)      mxi = 10
      if (tol.lt.1.d-16) tol = 1.d-16

      Dsig      = sigInf - sigY
      if (Dsig .lt. 1.d-10) then
        H = r0
        mHdDsig = r0
      else
        mHdDsig = - H / Dsig
      endif

      sigInf = sigInf * sr2d3
      Dsig   = Dsig   * sr2d3
      H      = H      * r2d3

c*****************************************************************************
c     CHECK FOR SOME ERRORS IN INPUT DATA
c*****************************************************************************

      if (nh.gt.maxh) call prgError(1,"nlkh_lC","maxh < nh!")
      
      if (outW) then
        if (abs(det_s(Cp)-r1).gt.1.d-2) 
     *    write(*,'(/a/)') ' warning: |det(Cp)-1| > 0.01 !'
        do j=1, nh
          if (abs(det_s(Cpi(1,j))-r1).gt.1.d-2) 
     *      write(*,'(/a,i1,a/)')' warning: |det(Cpi_',j,')-1| > 0.01 !'
        enddo
      endif

c*****************************************************************************
c     STRESS UPDATE   S, Cp, Cpi, xi
c*****************************************************************************

      call inverse_s (invC,C)
      call diff_ATA_m(diC,invC)      

      expFact = exp(mHdDsig * xi)

c     compute trial state
      call comp_S_X(invCp,diCp,S,dSdC,dSdCp,invCpi,diCpi,X,dXdCp,
     *              dXdCpi,Xt,
     *              C,invC,diC,Cp,Cpi,matData,nh,i0)
      call comp_Phi(CS,mCpXt,gp,n,m1dn,r2dn,fp,Phi,
     *              C,Cp,S,Xt,expFact,sigInf,Dsig)

c...  if plastic deformation then ............................................

      Dg = r0

      if (Phi.gt.tol) then

        call pmove(Cp,    Cpn,    6   )
        call pmove(invCp, invCpn, 6   )
        call pmove(Cpi,   Cpin,   6*nh)
        call pmove(invCpi,invCpin,6*nh)
        xin = xi

        call pzeroi(ctli,nh+1)
        call pzero (Rp,6)
        call pzero (Rpi,6*nh)

        if (outR) then
          write(*,'(/a,10(a))')
     *       '    Rp        Phi    ',('    Rpi    ',j=1,nh)
          write(*,'(12(g10.3,1x))') 
     *      norm2_s(Rp), abs(Phi), (norm2_s(Rpi(1,j)),j=1,nh)
        endif

c...    begin Newton loop for (Cp,Dg) ........................................

        iter1 = 0
        do while (noconv(Phi,Rp,Rpi,tol,nh).and.iter1.le.mxi)
          iter1 = iter1 + 1

          call pzero(Xt,6)
          call pzero(dXtdCp,36)

c...      begin loop over hardening components ...............................

          do j=1, nh
       
            r2dk   = r2 / matData(i0+j+j+j+1)

            nipCpi = nint(matData(i0+j+j+j+2))

c           compute Rpi, dRpi/dCpi
            call comp_Rpi(dfpidX,Rpi,dRdA,dRdf,dRdDg,
     *                    Cp,Cpi(1,j),X(1,j),invCpi(1,j),diCpi(1,1,j),
     *                    Cpin(1,j),invCpin(1,j),dXdCpi(1,1,j),Dg,r2dk,
     *                    h21,h22,h41,nipCpi)

c...        begin Newton loop for (Cpi) ......................................

            iter2 = 0
            do while ((norm2_s(Rpi).gt.tol.or.normx_s(Rpi).gt.tol)
     *      .and.iter2.le.mxi)
              iter2 = iter2 + 1

c             solve and update Cpi
              call decompLR_matrix(dRdA,P,6,0)
              call solve_matrix   (dRdA,P,Rpi,upd,6)
              do i=1, 3
                Cpi(i  ,j) = Cpi(i  ,j) -        upd(i)
                Cpi(i+3,j) = Cpi(i+3,j) - r1d2 * upd(i+3)
              enddo
              call inverse_s (invCpi(1,j),Cpi(1,j))
              call diff_ATA_m(diCpi(1,1,j),invCpi(1,j))

c             compute new X
              call comp_X(X(1,j),dXdCp(1,1,j),dXdCpi(1,1,j),
     *                    Cp,invCp,diCp,Cpi(1,j),invCpi(1,j),
     *                    diCpi(1,1,j),matData(i0+j+j+j))

c             compute Rpi, dRpi/dCpi
              call comp_Rpi(dfpidX,Rpi,dRdA,dRdf,dRdDg,
     *                      Cp,Cpi(1,j),X(1,j),invCpi(1,j),diCpi(1,1,j),
     *                      Cpin(1,j),invCpin(1,j),dXdCpi(1,1,j),Dg,
     *                      r2dk,
     *                      h21,h22,h41,nipCpi)

              if (outR) then
                do i=0, j
                  write(*,'(10x,''.'',$)')
                enddo
                write(*,'(g10.3)') norm2_s(Rpi)
              endif
            enddo
            if (iter2.gt.mxi) then
              write(*,'(/a/)')' strangely, it did not converge (Cpi)!'
              stop
            endif
            ctli(1+j) = max(ctli(1+j),iter2)

c...        end Newton loop for (Cpi) ........................................

c           compute dCpi(j)/dCp
            call comp_dCpi(dCpidCp(1,1,j),dCpidDg(1,j),
     *                     X(1,j),Cpi(1,j),r2dk,dfpidX,dXdCp(1,1,j),
     *                     dRdA,dRdf,dRdDg,
     *                     h41,h42,P)

c           add X(j) to Xtotal
            call set_s_ap(Xt,X(1,j))
            call set_4_ap(dXtdCp,dXdCp(1,1,j))

          enddo

c...      end loop over hardening components .................................

c         gp, n, fp, Phi
          call comp_Phi(CS,mCpXt,gp,n,m1dn,r2dn,fp,Phi,
     *                  C,Cp,S,Xt,expFact,sigInf,Dsig)

c         K
          call comp_K(K,Rp,dRdf,
     *                C,Cp,Cpn,Dg,invCp,diCp,invCpn,Xt,mCpXt,CS,gp,n,
     *                r2dn,m1dn,fp,dSdCp,dXtdCp,dXdCpi,dCpidCp,dCpidDg,
     *                expFact,H,
     *                h41,h42,h43,h44,h21,h22,h23,h24,h25,
     *                nipCp,nh)

c         solve and update Cp, Dg, Cpi, xi
          call set_s_p(R,Rp)
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
          fact = - upd(7)
          do j=1, nh
            call mult_s_4s_am(Cpi(1,j),dCpidCp(1,1,j),upd)
            call set_s_af(Cpi(1,j),dCpidDg(1,j),fact)
          enddo

c         update invCp, invCpi, S, X, Xt, ...
          call comp_S_X(invCp,diCp,S,dSdC,dSdCp,invCpi,diCpi,X,dXdCp,
     *                  dXdCpi,Xt,
     *                  C,invC,diC,Cp,Cpi,matData,nh,i0)

c         Phi, Rp, Rpi
          call comp_Phi  (CS,mCpXt,gp,n,m1dn,r2dn,fp,Phi,
     *                    C,Cp,S,Xt,expFact,sigInf,Dsig)
          call comp_R1000(Rp,fp,Cp,invCp,Cpn,invCpn,Dg,nipCp)

          do j=1, nh
            r2dk   = r2 / matData(i0+j+j+j+1)
            nipCpi = nint(matData(i0+j+j+j+2))
            call mult_u_ss_f(h21,Cp,X(1,j),r2dk)
            call dev_u_r    (h21)              ! now: h21 = gpi = 2/k dev(Cp X)
            call mult_s_us_p(h22,h21,Cpi(1,j)) ! now: h22 = fpi = gpi Cpi
            call comp_R1000 (Rpi(1,j),h22,Cpi(1,j),invCpi(1,j),
     *                       Cpin(1,j),invCpin(1,j),Dg,nipCpi)
          enddo

          if (outR) write(*,'(12(g10.3,1x))') 
     *      norm2_s(Rp), abs(Phi), (norm2_s(Rpi(1,j)),j=1,nh)

        enddo
        if (iter1.gt.mxi) then
          write(*,'(/a/)')' strangely, it did not converge (Cp,Dg)!'
          stop
        endif
        ctli(1) = iter1

c...    end Newton loop for (Cp,Dg) ..........................................

      endif
        
c...  end elasto-plastic stress update .......................................

c     compute determinants of Cp, Cpi
      ctlr(1) = det_s(Cp)
      do j=1, nh
        ctlr(1+j) = det_s(Cpi(1,j))
      enddo

c     write(*,'(6(f8.5)) ') xi
c     write(*,'(6(f8.5)) ') (Cp (  j),j=1,6)
c     write(*,'(6(f8.5)/)') (Cpi(j,1),j=1,6)

c*****************************************************************************
c     CONSISTENT TANGENT     2 dS/dC 
c*****************************************************************************

      if (Dg.gt.tol) then

c       dRp/dC, dPhi/dC
        call transpose_r        (mCpXt)            !
        call diff_norm2sTApUd_p (h21,C,mCpXt,CS,n) !
        call mult_s_s4_p        (h22,h21,dSdC)     !
        call transpose_r        (CS)               !
        call diff_norm2sATpUd_ap(h22,S,mCpXt,CS,n) ! now: h22 = d n / d C
        call diff_ATdB_p        (h41,C,Cp)         !
        call mult_4_44_f        (h42,h41,dSdC,r2dn)!
        call diff_TAdB_af       (h42,S,Cp,r2dn)    !
        call mult_4_ss_af       (h42,fp,h22,m1dn)  ! now: h42 = d fp / d C
        call mult_4_44_p        (h41,dRdf,h42)     ! now: h41 = d Rp / d C
                                                   !
        call inverse_matrix(K,h42,P,7)             ! now: h42 = K^{-1}
                                                   !
        do i=1, 6
          do j=1, 6
            h43(i,j) = h42(i,j)
            h44(i,j) = h41(i,j) 
          enddo
          h43(i,7) = h42(i,7)  !now: h43 = d Cp / d(Rp,Phi)   -> 6x7-matrix 
          h44(7,i) = h22(i)    !now: h44 = d(Rp,Phi) / d C    -> 7x6-matrix  
        enddo

        call mult_matrix_m (h41, h43,  h44,6,7,6)
        call mult_matrix_ap(dSdC,dSdCp,h41,6,6,6)

      endif

      call set_4_f(d2SdC,dSdC,r2)

      return
      end

c=============================================================================

      subroutine comp_S_X(invCp,diCp,S,dSdC,dSdCp,invCpi,diCpi,X,dXdCp,
     *                    dXdCpi,Xt,
     *                    C,invC,diC,Cp,Cpi,matData,nh,i0)
c
c     input:  C, invC, diC, Cp, Cpi, d
c
c     output: invCp  = Cp^{-1} 
c             diCp   = d Cp^{-1} / d Cp
c             S      = S(C,Cp)
c             dSdC   = d S / d C
c             dSdCp  = d S / d Cp
c             invCpi = Cpi^{-1} 
c             diCpi  = d Cpi^{-1} / d Cpi
c             X      = X(Cp,Cpi)
c             dXdCp  = d X / d Cp
c             dXdCpi = d X / d Cpi
c             Xt     = sum X
c
      implicit none
      integer          j, nh, i0, elTyp
      double precision invCp, diCp, S, dSdC, dSdCp, invCpi(6,*), 
     *                 diCpi(6,6,*), X(6,*), dXdCp(6,6,*), 
     *                 dXdCpi(6,6,*), Xt,
     *                 C, invC, diC, Cp, Cpi(6,*), matData(*),
     *                 bulk
      bulk  = matData(1)
      elTyp = nint(matData(3))
      call inverse_s (invCp,Cp)                     ! -> invCp
      call diff_ATA_m(diCp,invCp)                   ! -> diCp
      call computeElasticResponse(S,dSdC,dSdCp,     !
     *                C,invC,diC,Cp,invCp,diCp,     !
     *                bulk,matData(8),elTyp)        ! -> S, dSdC, dSdCp    
      call pzero     (Xt,6)                         !
      do j=1, nh                                    !
        call inverse_s (invCpi(1,j),Cpi(1,j))       ! -> invCpi
        call diff_ATA_m(diCpi(1,1,j),invCpi(1,j))   ! -> diCpi
        call comp_X  (X(1,j),dXdCp(1,1,j),          !
     *                dXdCpi(1,1,j),                !
     *                Cp,invCp,diCp,Cpi(1,j),       !
     *                invCpi(1,j),diCpi(1,1,j),     !
     *                matData(i0+j+j+j))            ! -> X, dXdCp, dXdCpi
        call set_s_ap(Xt,X(1,j))                    ! -> Xt      
      enddo                                         !
      return
      end

c=============================================================================

      subroutine comp_Phi(CS,mCpXt,gp,n,m1dn,r2dn,fp,Phi,
     *                    C,Cp,S,Xt,expFact,sigInf,Dsig)
c
c     input:  C, Cp, S, Xt, xi, d, sr2d3 = sqrt(2/3)
c
c     output: CS    = C S
c             mCpXt = Cp Xt
c             gp    = dev(C S - Cp Xt)
c             n     = sqrt( gp : gp^T )
c             m1dn  = - 1 / n
c             r2dn  = 2 / n
c             fp    = 2 / n  gp Cp
c             Phi   = n - sqrt(2/3) sigY(xi)
c
      implicit none
      double precision CS, mCpXt, gp, n, m1dn, r2dn, fp, Phi,
     *                 C, Cp, S, Xt, expFact, sigInf, Dsig,
     *                 norm2s_u
c     gp, n, fp, Phi
      call mult_u_ss_p(CS,C,S)           ! -> CS
      call set_u_p    (gp,CS)            !
      call mult_u_ss_m(mCpXt,Cp,Xt)      ! -> mCpXt
      call set_u_ap   (gp,mCpXt)         !   
      call dev_u_r    (gp)               ! -> gp
      n    = norm2s_u(gp)                ! -> n
      m1dn = - 1.d0 / n                  ! -> m1dn
      r2dn = - m1dn - m1dn               ! -> r2dn
      call mult_s_us_f(fp,gp,Cp,r2dn)    ! -> fp
      Phi  = n - sigInf + Dsig * expFact ! -> Phi
      return
      end

c=============================================================================

      subroutine comp_Rpi(dfpidX,Rpi,dRdA,dRdf,dRdDg,
     *                    Cp,Cpi,X,invCpi,diCpi,Cpin,invCpin,dXdCpi,Dg,
     *                    r2dk,
     *                    gpi,fpi,h4,nip)
c
c     input:  Cp, Cpi, X, invCpi, diCpi, Cpin, invCpin, dXdCpi, Dg, 
c             r2dk   = 2 / k
c
c     output: dfpidX = d fpi / d X
c             Rpi    = R(Cpi,fpi,Dg)
c             dRdA   = d Rpi / d Cpi
c             dRdf   = d Rpi / d fpi
c             dRdDg  = d Rpi / d Dg
c
      implicit none
      integer          nip
      double precision dfpidX, Rpi, dRdA, dRdf, dRdDg,
     *                 Cp, Cpi, X, invCpi, diCpi, Cpin, invCpin, dXdCpi, 
     *                 Dg, r2dk,
     *                 gpi, fpi, h4
c     compute Rpi, dRpi/dCpi
      call mult_u_ss_f (gpi,Cp,X,r2dk)           !
      call dev_u_r     (gpi)                     ! gpi = 2 b/c dev(Cp X)
      call mult_s_us_p (fpi,gpi,Cpi)             ! fpi = 2 b/c dev(Cp X) Cpi
      call diff_ATdB_f (dfpidX,Cp,Cpi,r2dk)      ! -> dfpidX
      call mult_4_44_p (h4,dfpidX,dXdCpi)        ! 
      call diff_UT_ap  (h4,gpi)                  ! now: h4 = d fpi/d Cpi
      call comp_R1111  (Rpi,dRdA,dRdf,dRdDg,fpi, ! 
     *                  Cpi,invCpi,diCpi,        !
     *                  Cpin,invCpin,Dg,nip)     ! -> Rpi, dRdf, dRdDg
      call mult_4_44_ap(dRdA,dRdf,h4)            ! -> dRdA
      return
      end

c=============================================================================

      subroutine comp_dCpi(dCpidCp,dCpidDg,
     *                     X,Cpi,r2dk,dfpidX,dXdCp,dRdA,dRdf,dRdDg,
     *                     h4,h42,P)
c
c     input:  X, Cpi, r2dk = 2 / k
c             dfpidX = d fpi / d X
c             dXdCp  = d X / d Cp
c             dRdA   = d Rpi / d Cpi 
c             dRdf   = d Rpi / d fpi
c             dRdDg  = d Rpi / d Dg
c
c     output: dCpidCp = d Cpi / d Cp
c             dCpidDg = d Cpi / d Dg
c
      implicit none
      integer          P, i, j
      double precision dCpidCp, dCpidDg, 
     *                 X, Cpi, r2dk, dfpidX, dXdCp, dRdA, dRdf, dRdDg,
     *                 h4(6,6), h42,
     *                 r1d4, r1d2
      data r1d4, r1d2 / 0.25d0, 0.5d0 /
c     compute dCpi(j)/dCp
      call diff_TAdB_f    (h4,X,Cpi,r2dk)   !
      call mult_4_44_ap   (h4,dfpidX,dXdCp) ! now: h4  = d fpi / d Cp
      call mult_4_44_p    (h42,dRdf,h4)     ! now: h42 = d Rpi / d Cp
      call decompLR_matrix(dRdA,P,6,0)      !
      call inverse_matrix (dRdA,h4,P,6)     !
      do i=4, 6                             !
        do j=1, 3                           !
          h4(i  ,j) = h4(i  ,j) * r1d2      !
          h4(j  ,i) = h4(j  ,i) * r1d2      !
          h4(j+3,i) = h4(j+3,i) * r1d4      !
        enddo                               ! now: h4 = [d Rpi / d Cpi]^{-1}
      enddo                                 !
      call mult_4_44_m(dCpidCp,h4,h42)      ! -> dCpidCp 
c     compute dCpi(j)/dDgam                 !
      call mult_s_4s_m(dCpidDg,h4,dRdDg)    ! -> dCpidDg
      return
      end

c=============================================================================

      subroutine comp_K(K,Rp,dRdf,
     *                  C,Cp,Cpn,Dg,invCp,diCp,invCpn,Xt,mCpXt,CS,gp,n,
     *                  r2dn,m1dn,fp,dSdCp,dXtdCp,dXdCpi,dCpidCp,
     *                  dCpidDg,expFact,H,
     *                  dfpdCp,dRdA,h4,h42,dndCp,dfpdDg,dRdDg,hs,hs2,
     *                  nip,nh)
c
c     input:  C, Cp, Cpn, Dg
c             invCp   = Cp^{-1}
c             diCp    = d Cp^{-1} / d Cp
c             invCpn  = Cp_n^{-1}
c             Xt      = sum X
c             mCpXt   = - Cp Xt      NOTE: at output mCpXt is transposed!
c             CS      = C S
c             gp      = dev(C S - Cp Xt)
c             n       = sqrt( gp : gp^T )
c             r2dn    = 2 / n
c             m1dn    = - 1 / n
c             fp      = 2 / n gp Cp
c             dSdCp   = d S / d Cp
c             dXtdCp  = d Xt / d Cp
c             dXdCpi  = d X / d Cpi
c             dCpidCp = d Cpi / d Cp
c             dCpidDg = d Cpi / d Dg
c             xi
c             d
c
c     output: K
c             Rp      = R(Cp,fp,Dg)
c             dRdf    = d R(Cp,fp,Dg) / d fp
c
      implicit none
      integer          i, j, nh, nip
      double precision K(7,*), Rp, dRdf,
     *                 C, Cp, Cpn, Dg, invCp, diCp, invCpn, Xt, mCpXt, 
     *                 CS, gp, n, r2dn, m1dn, fp, dSdCp, dXtdCp,
     *                 dXdCpi(6,6,*), dCpidCp(6,6,*), dCpidDg(6,*), 
     *                 expFact, H,
     *                 dfpdCp, dRdA(6,*), h4, h42, 
     *                 dndCp(*), dfpdDg, dRdDg(*), hs, hs2,      fact,
     *                 dot_ss

c     partial derivatives dfp/dCp, dfp/dDg, dn/dCp, dn/dDg 
      call diff_norm2sATpUd_p(hs,C,mCpXt,CS,n)      !    
      call mult_s_s4_p (dndCp,hs,dSdCp)             !
      call diff_norm2sATpUd_m(hs,Cp,CS,mCpXt,n)     ! now: hs = d n / d Xt
      call mult_s_s4_ap(dndCp,hs,dXtdCp)            ! 
      call transpose_r (mCpXt)                      !
      call diff_norm2sTApUd_am(dndCp,Xt,CS,mCpXt,n) ! now: dndCp = d n / d Cp
      call diff_UT_p   (h4,gp)                      !
      call diff_TAdB_am(h4,Xt,Cp)                   ! 
      call diff_ATdB_p (h42,C,Cp)                   !
      call mult_4_44_ap(h4,h42,dSdCp)               !
      call diff_ATdB_p (h42,Cp,Cp)                  ! now: h42 = d(n/2 fp)/d Xt
      call mult_4_44_am(h4,h42,dXtdCp)              ! 
      call set_4_f     (dfpdCp,h4,r2dn)             ! 
      call mult_4_ss_af(dfpdCp,fp,dndCp,m1dn)       ! now: dfpdCp = d fp / d Cp
      call pzero       (dfpdDg,6)                   !
      K(7,7) = - H * expFact                        !

c     + contribution arising from Cpi(Cp,Dg), dfp/dCp, dfp/dDg, dn/dCp, dn/dDg 
      fact = - r2dn                                  !
      do j=1, nh                                     !
        call mult_s_s4_p (hs2,hs,dXdCpi(1,1,j))      ! now: hs2 = dn/dCpi
        call mult_4_ss_f (h4,fp,hs2,m1dn)            ! 
        call mult_4_44_af(h4,h42,dXdCpi(1,1,j),fact) ! h4 = dfp/dCpi
        call mult_4_44_ap(dfpdCp,h4,dCpidCp(1,1,j))  !
        call mult_s_4s_ap(dfpdDg,h4,dCpidDg(1,j))    !
        call mult_s_s4_ap(dndCp,hs2,dCpidCp(1,1,j))  !
        K(7,7) = K(7,7) +dot_ss(hs2,dCpidDg(1,j))    !
      enddo                                          !

c     Rp, dRp/dCp, dRp/dfp, dRp/dDg                 
      call comp_R1111  (Rp,dRdA,dRdf,dRdDg,fp,Cp,    ! 
     *                  invCp,diCp,Cpn,invCpn,Dg,nip)! -> Rp
      call mult_4_44_ap(dRdA,dRdf,dfpdCp)            !
      call mult_s_4s_ap(dRdDg,dRdf,dfpdDg)           !

c     set matrix K                                  
      do i=1, 6                                      !
        do j=1, 6                                    !
          K(i,j) = dRdA(i,j)                         !
        enddo                                        !
        K(i,7) = dRdDg(i)                            !
        K(7,i) = dndCp(i)                            !
      enddo                                          ! -> K
      return
      end

c=============================================================================

      logical function noconv(Phi,Rp,Rpi,tol,nh)
c
c     check whether maximum and Euclidian norms are less or equal to tol
c
      implicit none
      integer          nh, i
      double precision Phi, Rp(*), Rpi(*), tol, norm2_s
      if (abs(Phi).gt.tol) then
        noconv = .true.
        return
      endif
      i = 1
      do while (abs(Rp(i)).le.tol.and.i.le.6)
        i = i + 1
      enddo
      if (i.le.6) then
        noconv = .true.
        return
      endif
      i = 1
      do while (abs(Rpi(i)).le.tol.and.i.le.6*nh)
        i = i + 1
      enddo
      if (i.le.6*nh) then
        noconv = .true.
        return
      endif
      if (norm2_s(Rp).gt.tol) then
        noconv = .true.
        return
      endif
      do i=1, nh
        if (norm2_s(Rpi((i-1)*6+1)).gt.tol) then
          noconv = .true.
          return
        endif
      enddo
      noconv = .false.
      return
      end

c=============================================================================

      subroutine comp_X(X,dXdCp,dXdCpi,
     *                  Cp,invCp,diCp,Cpi,invCpi,diCpi,mu)
      implicit none
      double precision X, dXdCp, dXdCpi, Cp, invCp, diCp, 
     *                 Cpi, invCpi, diCpi, mu,
     *                 h4(6,6), fact, trp, traceAB_ss
     
c    Simple NEO-HOOKE type back stress
c
c    chi = mu Jpe^(-5/3) dev(Bpe)
c   -> S = mu Jpe^(-2/3) (invCpi - tr(Cp invCpi)/3 invCp)
c
c    we use assumption of inelastic incompressibility 
c    det(Cp) = det(Cpi) = 1 -> Jpe = 1
c     
c       here: C  represents Cp,
c             Cp represents Cpi 
c             S  represents X
c
      fact = - mu / 3.d0
      trp  = - fact * traceAB_ss(Cp,invCpi)

c     X
      call set_s_f (X,invCpi,mu)
      call set_s_af(X,invCp,trp)

c     dXdCp
      call set_4_f     (dXdCp,diCp,trp)
      call mult_4_ss_af(dXdCp,invCp,invCpi,fact)

c     dXdCpi
      call mult_4_ss_f (h4,invCp,Cp,fact)
      call unit_4_af   (h4,mu)
      call mult_4_44_p (dXdCpi,h4,diCpi)

      return
      end

c=============================================================================

