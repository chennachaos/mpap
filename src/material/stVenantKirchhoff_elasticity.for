

      subroutine stVenantKirchhoff_elasticity(matData,F,stre,cc,finite)

      implicit none

      integer          finite,
     *                 i, round
      double precision matData(*), F(3,*), stre(*), cc(6,*),
     *                 mu, K, Lamb, b(6), fact,
     *                 det_u, r1dJ, J

      if (finite.eq.0) then
        call small_strain_elasticity(matData,F,stre,cc)
        return
      endif

      K    = matData(1)
      mu   = matData(2)
      Lamb = K - (mu + mu) / 3.d0

      call calcB(b,F)

      J    = det_u(F)
      r1dJ = 1.d0 / J

      fact = r1dJ * (Lamb * 0.5d0 * (b(1)+b(2)+b(3)-3.d0) - mu)

      call set_s_f(stre,b,fact)

      fact = r1dJ * mu

      call mult_s_ss_af(stre,b,b,fact)

      fact = r1dJ * Lamb

      call mult_4_ss_f(cc,b,b,fact)

      fact = r1dJ * (mu + mu)

      cc(1,1) = cc(1,1) + fact * b(1)*b(1)
      cc(1,2) = cc(1,2) + fact * b(4)*b(4)
      cc(1,3) = cc(1,3) + fact * b(6)*b(6)
      cc(2,1) = cc(1,2)
      cc(2,2) = cc(2,2) + fact * b(2)*b(2)
      cc(2,3) = cc(2,3) + fact * b(5)*b(5)
      cc(3,1) = cc(1,3)
      cc(3,2) = cc(2,3)
      cc(3,3) = cc(3,3) + fact * b(3)*b(3)

      cc(1,4) = cc(1,4) + fact * b(1)*b(4)
      cc(1,5) = cc(1,5) + fact * b(4)*b(6)
      cc(1,6) = cc(1,6) + fact * b(1)*b(6)
      cc(2,4) = cc(2,4) + fact * b(2)*b(4)
      cc(2,5) = cc(2,5) + fact * b(2)*b(5)
      cc(2,6) = cc(2,6) + fact * b(4)*b(5)
      cc(3,4) = cc(3,4) + fact * b(5)*b(6)
      cc(3,5) = cc(3,5) + fact * b(3)*b(5)
      cc(3,6) = cc(3,6) + fact * b(3)*b(6)

      cc(4,1) = cc(1,4)      
      cc(5,1) = cc(1,5)
      cc(6,1) = cc(1,6)
      cc(4,2) = cc(2,4)
      cc(5,2) = cc(2,5)
      cc(6,2) = cc(2,6)
      cc(4,3) = cc(3,4)
      cc(5,3) = cc(3,5)
      cc(6,3) = cc(3,6)

      fact = fact * 0.5d0

      cc(4,4) = cc(4,4) + fact * (b(1)*b(2) + b(4)*b(4))
      cc(4,5) = cc(4,5) + fact * (b(2)*b(6) + b(4)*b(5))
      cc(4,6) = cc(4,6) + fact * (b(1)*b(5) + b(4)*b(6))
      cc(5,4) = cc(4,5)
      cc(5,5) = cc(5,5) + fact * (b(2)*b(3) + b(5)*b(5))
      cc(5,6) = cc(5,6) + fact * (b(3)*b(4) + b(5)*b(6))
      cc(6,4) = cc(4,6)
      cc(6,5) = cc(5,6)
      cc(6,6) = cc(6,6) + fact * (b(1)*b(3) + b(6)*b(6))

      return

      end





