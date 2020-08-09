
      subroutine ogden_elasticity(matData,mu,alpha,F,sig,cc,finite)

      implicit none

      integer          finite,
     *                 nOgd, i, l, m, n, round,
     *                 spectralDecomposition
      double precision matData(*), mu(*), alpha(*), F(3,*), sig(*), 
     *                 cc(6,*),
     *                 matDataSmall(2), K, b(6), lam2(3), dtau(3,3), 
     *                 tau6(6), h4(6,6), J, lam(3), lampa(3), press, 
     *                 tau(3), fact, fact2, E(6,3),
     *                 r1, r1d2, r1d3, r1d6

      r1   = 1.d0
      r1d2 = 0.5d0
      r1d3 = r1 / 3.d0
      r1d6 = r1d3 * r1d2

      K    = matData(1)
      nOgd = round(matData(2))
 
      if (finite.lt.1) then
        matDataSmall(1) = 0.d0
        do i=1, nOgd
          matDataSmall(1) = matDataSmall(1) + mu(i) * alpha(i)
        enddo
        matDataSmall(1) = matDataSmall(1) * 0.5d0
        matDataSmall(2) = K
        call small_strain_elasticity(matDataSmall,F,sig,cc)
        return
      endif

      call calcB(b,F)

      if (spectralDecomposition(lam2,E,b).ne.1)

     *   call prgError(1,"ogden_elasticity","fatal error!")

      lam(1) = sqrt(lam2(1))
      lam(2) = sqrt(lam2(2))
      lam(3) = sqrt(lam2(3))

      J = lam(1) * lam(2) * lam(3)

      call pzero(tau,3)
      call pzero(dtau,9)

      do i=1, nOgd

        !write(*,*) mu(i), alpha(i)
 
        lampa(1) = lam(1)**alpha(i) 
        lampa(2) = lam(2)**alpha(i) 
        lampa(3) = lam(3)**alpha(i) 

        fact2    = mu(i) * J**(-alpha(i)*r1d3) * r1d6

        do l=1, 3

          m = l + 1
          if (m .eq. 4) m = 1
          n = m + 1
          if (n .eq. 4) n = 1

          fact  = lampa(l) + lampa(l) - lampa(m) - lampa(n)

          tau(l) = tau(l) + (fact2 + fact2) * fact

          dtau(l,l) = dtau(l,l) + fact2 * alpha(i) / lam2(l) *
     *                         (lampa(l) + lampa(l) - r1d3 * fact)

          dtau(l,m) = dtau(l,m) - fact2 * alpha(i) / lam2(m) * 
     *                         (lampa(m)            + r1d3 * fact)

          dtau(l,n) = dtau(l,n) - fact2 * alpha(i) / lam2(n) * 
     *                         (lampa(n)            + r1d3 * fact)

        enddo

      enddo

      call isotropicTensorFunction(tau6,h4,lam2,tau,dtau,E,b)

      call pzero(sig,6)
      call pzero(cc,36)

      call taub2sigccPart1(sig,cc,tau6,h4,b,J)

c...  compute volumetric response (pressure)

      press = r1d2 * K * (J - r1/J)
      fact  = K * J
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

