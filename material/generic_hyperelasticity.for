
      subroutine generic_hyperelasticity(matData,F,stre,cc,finite)
      implicit none

      integer          mxI
      parameter        (mxI = 10) ! number of invariants defined below

      integer          finite,
     *                 ptyp, j, jj, l, l0, nI, nFib

      double precision matData(*), F(3,*), stre(*), cc(6,*), 
     *                 S(6), C4(6,6), I(mxI), dI(6,mxI), ddI(6,6,mxI), 
     *                 dW(mxI), ddW(mxI,mxI), C(6), M(6,mxI), CM(3,3), 
     *                 CCM(3,3), invC(6), One(6), fact, fct1, fct2,

     *                 mu, K, lam, EA, E, GA, nu, alp, beta, gam, 
     *                 mum, muf, Km, Kf,

     *                 r0, r1d4, r1d2, r1, r3d2, r2, r5d2, r3, r4, r8,
     *                 r1d3, m2d3, m5d3, m8d3,

     *                 traceAB_ss, det_s

      data   r0,  r1d4,r1d2, r1, r3d2,  r2, r5d2,  r3,  r4,  r8
     *    / 0.d0,.25d0,.5d0,1.d0,1.5d0,2.d0,2.5d0,3.d0,4.d0,8.d0/

      r1d3 = r1 / 3.d0
      m2d3 = - r1d3 - r1d3
      m5d3 = m2d3 - r1
      m8d3 = m5d3 - r1

      nI   = mxI

      nFib = nint(matData(1))

      if (nFib+nFib+4.gt.mxI) 
     *  call prgError(1,"generic_hyperelasticity",
     *                  "increase mxI in the code!")

      l0 = 1

      do l=1, nFib

        M(1,l) = matData(1+l0) * matData(1+l0)
        M(2,l) = matData(2+l0) * matData(2+l0)
        M(3,l) = matData(3+l0) * matData(3+l0)
        M(4,l) = matData(1+l0) * matData(2+l0)
        M(5,l) = matData(2+l0) * matData(3+l0)
        M(6,l) = matData(3+l0) * matData(1+l0)

        fact = r1 / (M(1,l) + M(2,l) + M(3,l))

        do j=1, 6
          M(j,l) = M(j,l) * fact
        enddo

        l0 = l0 + 3

      enddo

      l0 = l0 + 1

      ptyp = nint(matData(l0))

      call unit_s_p(One)

      call mult_s_uTu_p(C,F,F)

      call inverse_s(invC,C)

c-------------------------------------------------------
c
c     tr(C) 
     
      I(1) = C(1) + C(2) + C(3)

      call unit_s_p(dI(1,1)) 

      call pzero(ddI(1,1,1),36)

c-------------------------------------------------------
c
c     0.5 * (tr(C)^2 - tr(C*C)) 

      I(2) = r1d2 * (I(1)*I(1) - traceAB_ss(C,C))

      call unit_s_f(dI(1,2),I(1))
      call set_s_ap(dI(1,2),C)

      call unit_4_p(ddI(1,1,2))
      call mult_4_ss_ap(ddI(1,1,2),One,One)
          
c-------------------------------------------------------
c
c     J = sqrt(det(C))

      I(3) = sqrt(det_s(C))

      fact = I(3) * 0.5d0
      call set_s_f(dI(1,3),invC,fact)

      fact = - fact
      call diff_ATA_f(ddI(1,1,3),invC,fact)
      fact = - fact * 0.5d0
      call mult_4_ss_af(ddI(1,1,3),invC,invC,fact)

c-------------------------------------------------------
c
c     tr(CC)

      I(4) = traceAB_ss(C,C)

      call set_s_f(dI(1,4),C,r2)

      call unit_4_f(ddI(1,1,4),r2)

c-------------------------------------------------------
c
      do l=1, nFib

        call mult_u_ss_p(CM,C,M(1,l))

        call mult_u_su_p(CCM,C,CM)

c       tr(CM)

        j = 3 + l + l

        I(j) = CM(1,1) + CM(2,2) + CM(3,3)

        call set_s_p(dI(1,j),M(1,l))

        call pzero(ddI(1,1,j),36)

c-------------------------------------------------------
c
c       tr(CCM)

        j = j + 1

        I(j) = CCM(1,1) + CCM(2,2) + CCM(3,3)

        call sympart_f(dI(1,j),CM,r2)

        call diff_symAT_f(ddI(1,1,j),M(1,l),r2)

      enddo

c-------------------------------------------------------

      call pzero(dW,nI)
      call pzero(ddW,nI*mxI)

c...  derivatives of free energy function

      goto (1,2,3,4,10,10,10,10,10), ptyp

c=========================================================================
 1    continue

c     NEARLY INCOMPRESSIBLE NEO-HOOKE  

c     Bonet & Wood, eqs. (5.49),(5.57)

c     P. Wriggers, Nichtlineare FEM, page 46

c     W = 1/2 mu ( J^(-2/3) tr(C) - 3 ) + 1/4 K (J^2-1) - K/2 lnJ

c       = 1/2 mu ( I3^(-2/3) I1   - 3 ) + 1/4 K (I3^2-1) - K/2 lnI3

      K  = matData(l0+1)
      mu = matData(l0+2)

      dW(1)    = r1d2 * mu * I(3)**m2d3 
      dW(3)    = - r1d3 * mu * I(3)**m5d3 * I(1)
     *             + K * r1d2 * (I(3) - r1/I(3))

      ddW(1,3) = - r1d3 * mu * I(3)**m5d3

      ddW(3,1) = ddW(1,3)
      ddW(3,3) = - r1d3 * m5d3 * mu * I(3)**m8d3 * I(1)
     *           + K * (r1d2 + r1d2/(I(3)*I(3)))

      nI = 3

      goto 999
c=========================================================================
 2    continue

c     TRANSVERSELY ISOTROPIC ST.VENANT-KIRCHHOFF 

c     Bonet & Burton, 1998, eqs. (46),(50)

c     W =     1/2 lam tr(E)^2 + mu E:E        (isotropic part with E = 1/2 (C-I))

c          +  (alp + beta (I1-3) + gam (I4-1)) (I4-1) - 1/2 alp (I5-1)    (anisotropic part)


c       = 1/4 lam (I1 - 3)^2 + 1/4 mu (I4- 2 I1 + 3)
c          +  (alp + beta (I1-3) + gam (I5-1)) (I5-1) - 1/2 alp (I6-1)

      EA       = matData(l0+1)
      E        = matData(l0+2)
      GA       = matData(l0+3)
      nu       = matData(l0+4)
      
      fct2     = EA / E
      fct1     = r1 - nu - r2*fct2*nu*nu

      lam      = E * (nu + fct2*nu*nu) / (fct1 * (r1 + nu))
      mu       = E / (r2 * (r1 + nu))
      alp      = mu - GA
      beta     = E * nu*nu * (r1 - fct2) / (r4 * fct1 * (r1 + nu))
      gam      = EA * (r1 - nu) / (r8 * fct1)
     *              - (lam + mu+mu) / r8 + alp*r1d2 - beta

      dW(1)    = r1d2 * (lam*(I(1)-r3) - mu)  +  (I(5)-r1) * beta
      dW(4)    = r1d4 * mu
      dW(5)    = alp + beta * (I(1)-r3) + gam * (I(5)-r1) * r2
      dW(6)    = - r1d2 * alp

      ddW(1,1) = r1d2 * lam
      ddW(1,5) = beta

      ddW(5,1) = beta
      ddW(5,5) = gam + gam

      nI = 6

      goto 999
c=========================================================================
 3    continue

c     TRANSVERSELY ISOTROPIC NEO-HOOKE (I)

c     Bonet & Burton, 1998, eqs. (35),(50)

c     W = 1/2 mu (I1 - 3) - mu ln(J) + 1/2 lam (J-1)^2                    (isotropic part)

c          +  [alp + beta (I1-3) + gam (I4-1)] (I4-1) - 1/2 alp (I5-1)    (anisotropic part)


c       = 1/2 mu (I1 - 3) - mu ln(I3) + 1/2 lam (I3-1)^2
c          +  [alp + beta (I1-3) + gam (I5-1)] (I5-1) - 1/2 alp (I6-1)

      EA       = matData(l0+1)
      E        = matData(l0+2)
      GA       = matData(l0+3)
      nu       = matData(l0+4)
      
      fct2     = EA / E
      fct1     = r1 - nu - r2*fct2*nu*nu

      lam      = E * (nu + fct2*nu*nu) / (fct1 * (r1 + nu))
      mu       = E / (r2 * (r1 + nu))
      alp      = mu - GA
      beta     = E * nu*nu * (r1 - fct2) / (r4 * fct1 * (r1 + nu))
      gam      = EA * (r1 - nu) / (r8 * fct1)
     *              - (lam + mu+mu) / r8 + alp*r1d2 - beta

      dW(1)    = r1d2 * mu   +  (I(5)-r1) * beta
      dW(3)    = - mu / I(3) + lam * (I(3) - r1)
      dW(5)    = alp + beta * (I(1)-r3) + gam * (I(5)-r1) * r2
      dW(6)    = - r1d2 * alp

      ddW(1,5) = beta

      ddW(3,3) = mu / (I(3)*I(3)) + lam

      ddW(5,1) = beta
      ddW(5,5) = gam + gam

      nI = 6

      goto 999
c=========================================================================
 4    continue

c     TRANSVERSELY ISOTROPIC NEO-HOOKE (II)

c     Bonet & Burton, 1998, eqs. (35),(50)

c     W = 1/2 mu (I1 - 3) - mu ln(J) + 1/2 lam (J-1)^2                    (isotropic part)

c          +  (alp + beta ln(J) + gam (I4-1)) (I4-1) - 1/2 alp (I5-1)    (anisotropic part)


c       = 1/2 mu (I1 - 3) - mu ln(I3) + 1/2 lam (I3-1)^2
c          +  (alp + beta ln(I3) + gam (I5-1)) (I5-1) - 1/2 alp (I6-1)

      EA       = matData(l0+1)
      E        = matData(l0+2)
      GA       = matData(l0+3)
      nu       = matData(l0+4)
      
      fct2     = EA / E
      fct1     = r1 - nu - r2*fct2*nu*nu

      lam      = E * (nu + fct2*nu*nu) / (fct1 * (r1 + nu))
      mu       = E / (r2 * (r1 + nu))
      alp      = mu - GA
      beta     = E * nu*nu * (r1 - fct2) / (r4 * fct1 * (r1 + nu))
      gam      = EA * (r1 - nu) / (r8 * fct1)
     *              - (lam + mu+mu) / r8 + alp*r1d2 - beta

      dW(1)    = r1d2 * mu
      dW(3)    = - mu/I(3) + lam*(I(3)-r1) + (I(5)-r1)*beta/I(3)
      dW(5)    = alp + beta * log(I(3)) + gam * (I(5)-r1) * r2
      dW(6)    = - r1d2 * alp

      ddW(3,3) = mu/(I(3)*I(3)) + lam - (I(5)-r1)*beta/(I(3)*I(3))
      ddW(3,5) = beta / I(3)

      ddW(5,3) = beta / I(3)
      ddW(5,5) = gam + gam

      nI = 6

      goto 999
c=========================================================================
 10   continue

      call prgError(1,'generic_hyperelasticity',
     *           'what shall we do with the drunken sailor?!')

999   continue

c...  second Piola Kirchhoff stress

      call pzero(S,6)

      do j=1, nI
        fact = dW(j) + dW(j)
        call set_s_af(S,dI(1,j),fact)
      enddo

c...  material tangent tensor 2 dS/dC

      call pzero(C4,36)

      do j=1, nI
        fact = r4 * dW(j)
        call set_4_af(C4,ddI(1,1,j),fact)
        do jj=1, nI
          fact = r4 * ddW(j,jj)
          call mult_4_ss_af(C4,dI(1,j),dI(1,jj),fact)
        enddo
      enddo

c...  push forward

      call pushfwrd_s(stre,F,S)

      call pushfwrd_4(cc,F,C4)

      fact = r1 / I(3)

      do l=1, 6
        stre(l) = stre(l) * fact
        do j=1, 6
          cc(l,j) = cc(l,j) * fact
        enddo
      enddo
 
      return

      end

