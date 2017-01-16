

      subroutine small_strain_elastic_lamina2D(matData,F,stre,cc)

      implicit none

      integer          i, j

      double precision matData(*), F(2,*), stre(*), cc(4,*),
     *                 E1, E2, nu12, nu21, G12, alpha, eps(3), fact,
     *                 Q11, Q22, Q12, Q66, m, m2, m4, n, n2, n4, m2n2,
     *                 m3n, mn3, Qxx, Qyy, Qxy, Qxs, Qys, Qss,
     *                 r1d2, r1, r2, r4

      data r1d2, r1, r2, r4 / 0.5d0, 1.d0, 2.d0, 4.d0 /

      E1    = matData(1)
      E2    = matData(2)
      nu12  = matData(3)
      nu21  = E2 / E1 * nu12
      G12   = matData(5)
      alpha = matData(6) * acos(-r1) / 180.d0

c     strain tensor

      eps(1) = F(1,1) - r1
      eps(2) = F(2,2) - r1
      eps(3) = r1d2 * (F(1,2) + F(2,1))

c     tangent tensor

      fact = r1 / (r1 - nu12*nu21)

      Q11  = fact * E1
      Q22  = fact * E2
      Q12  = nu21 * Q11
      Q66  = G12

      m    = cos(alpha)
      n    = sin(alpha)

      m2   = m * m
      m4   = m2 * m2
      n2   = n * n
      n4   = n2 * n2
      m2n2 = m2 * n2
      mn3  = m * n * n2
      m3n  = m * n * m2

      Qxx  = m4 * Q11 + n4 * Q22 + (m2n2+m2n2) * Q12 + r4*m2n2 * Q66
      Qyy  = n4 * Q11 + m4 * Q22 + (m2n2+m2n2) * Q12 + r4*m2n2 * Q66
      Qxy  = m2n2 * (Q11 + Q22) + (m4 + n4) * Q12 - r4*m2n2 * Q66
      Qxs  = m3n * Q11 - mn3 * Q22 + (mn3-m3n) * (Q12 + Q66 + Q66)
      Qys  = mn3 * Q11 - m3n * Q22 + (m3n-mn3) * (Q12 + Q66 + Q66)
      Qss  = m2n2 * (Q11 + Q22 - Q12 - Q12) + (m4-m2n2-m2n2+n4) * Q66

      cc(1,1) = Qxx
      cc(1,2) = Qxy
      cc(1,3) = Qxs
      cc(2,1) = Qxy
      cc(2,2) = Qyy
      cc(2,3) = Qys
      cc(3,1) = Qxs
      cc(3,2) = Qys
      cc(3,3) = Qss

c     stress

      stre(1) = cc(1,1)*eps(1) + cc(1,2)*eps(2) + r2*cc(1,3)*eps(3)
      stre(2) = cc(2,1)*eps(1) + cc(2,2)*eps(2) + r2*cc(2,3)*eps(3)
      stre(3) = cc(3,1)*eps(1) + cc(3,2)*eps(2) + r2*cc(3,3)*eps(3)

      return

      end






