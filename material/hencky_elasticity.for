

      subroutine hencky_elasticity(matData,F,stre,cc,finite)

      implicit none

      integer          finite,
     *                 i, round
      double precision matData(*), F(3,*), stre(*), cc(6,*),
     *                 mu, K, Lamb, b(6), fact, h(6),
     *                 tau(6), dtaudb(6,6),
     *                 det_u, r1dJ, J, J2

      if (finite.eq.0) then
        call small_strain_elasticity(matData,F,stre,cc)
        return
      endif

      K    = matData(1)
      mu   = matData(2)

      call calcB(b,F)

      J    = det_u(F)

      J2   = J * J

      call log_s(tau,dtaudb,b,10)

      call set_s_f(tau,tau,mu)

      fact = K * .5d0 * (J2 - 1.d0) - mu / 1.5d0 * log(J)

      call unit_s_p(h)
      call set_s_af(tau,h,fact)

      call set_4_f(dtaudb,dtaudb,mu)

      call pzero(stre,6)
      call pzero(cc, 36)

      call taub2sigccPart1(stre,cc,tau,dtaudb,b,J)

      fact = K * J - mu / (J*1.5d0)

      do i=1, 3
        cc(i,1) = cc(i,1) + fact
        cc(i,2) = cc(i,2) + fact
        cc(i,3) = cc(i,3) + fact
      enddo

      call taub2sigccPart2(stre,cc)

      return

      end





