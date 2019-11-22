

      subroutine small_strain_anisotropic_elasticity(matData,F,stre,cc)

      implicit none

      integer          i,
     *                 ptyp, nFib, round

      double precision matData(*), F(3,*), stre(*), cc(6,*),
     *                 a, b, c, d, e, m1, m2, m3, M(6), eps(6), 
     *                 E1, E2, nu12, nu23, G12, fact, 
     *                 r1d2, r1, r4

      data r1d2, r1, r4 / 0.5d0, 1.d0, 4.d0 /

c     strain tensor

      eps(1) = F(1,1) - r1
      eps(2) = F(2,2) - r1
      eps(3) = F(3,3) - r1
      eps(4) = r1d2 * (F(1,2) + F(2,1))
      eps(5) = r1d2 * (F(2,3) + F(3,2))
      eps(6) = r1d2 * (F(3,1) + F(1,3))

c     tangent tensor

      ptyp   = round(matData(1))

      if (ptyp.eq.1) then

        a    = matData(2)
        b    = matData(3)
        nFib = round(matData(4))

        call ccIsotropicPart(cc,a,b)

        do i=1, nFib

          c  = matData( 5+(i-1)*6)
          d  = matData( 6+(i-1)*6)
          e  = matData( 7+(i-1)*6)
          m1 = matData( 8+(i-1)*6)
          m2 = matData( 9+(i-1)*6)
          m3 = matData(10+(i-1)*6)

          fact = 1.d0 / (m1*m1 + m2*m2 + m3*m3)

          M(1) = m1 * m1 * fact
          M(2) = m2 * m2 * fact
          M(3) = m3 * m3 * fact
          M(4) = m1 * m2 * fact
          M(5) = m2 * m3 * fact
          M(6) = m3 * m1 * fact

          call ccFibrePart(cc,c,d,e,M)

        enddo

      else if (ptyp.eq.2) then

        E1   = matData(2)
        E2   = matData(3)
        nu12 = matData(4)
        nu23 = matData(5)
        G12  = matData(6)
        m1   = matData(7)
        m2   = matData(8)
        m3   = matData(9)

        fact = r1 / ((E1+E1) * nu12 * nu12 + E2 * (nu23-r1))

        a    = - E2 * (E1 * nu12 * nu12 + E2 * nu23) / (r1+nu23) * fact
        b    = E2 / (r1 + nu23)
        c    = E2 * (E1*nu12*(nu12-nu23-r1)+E2*nu23) / (r1+nu23) * fact
        d    = E1-r4*G12-r1d2*fact*(E2-(E1+E1)*nu12)*(E2-(E1+E1)*nu12)
        e    = G12+G12 - b

        call prgError(2,"small_strain_anisotropic_elasticity",
     *                  "ptyp = 2 has not yet been implemented!")

        fact = 1.d0 / (m1*m1 + m2*m2 + m3*m3)

        M(1) = m1 * m1 * fact
        M(2) = m2 * m2 * fact
        M(3) = m3 * m3 * fact
        M(4) = m1 * m2 * fact
        M(5) = m2 * m3 * fact
        M(6) = m3 * m1 * fact

        call ccIsotropicPart(cc,a,b)

        call ccFibrePart(cc,c,d,e,M)

      else

        call prgError(1,"small_strain_anisotropic_elasticity",
     *                     "invalid ptyp!")

      endif

c     stress

      call mult_s_4s_p(stre,cc,eps)

      return

      end







      subroutine ccIsotropicPart(cc,a,b)
      implicit none
      double precision cc(6,*), a, b, fact, r0

      r0 = 0.d0

      fact = b * .5d0

      cc(1,1) = a + b
      cc(1,2) = a
      cc(1,3) = a
      cc(1,4) = r0
      cc(1,5) = r0
      cc(1,6) = r0
      cc(2,1) = a
      cc(2,2) = a + b
      cc(2,3) = a
      cc(2,4) = r0
      cc(2,5) = r0
      cc(2,6) = r0
      cc(3,1) = a
      cc(3,2) = a
      cc(3,3) = a + b
      cc(3,4) = r0
      cc(3,5) = r0
      cc(3,6) = r0
      cc(4,1) = r0
      cc(4,2) = r0
      cc(4,3) = r0
      cc(4,4) = fact
      cc(4,5) = r0
      cc(4,6) = r0
      cc(5,1) = r0
      cc(5,2) = r0
      cc(5,3) = r0
      cc(5,4) = r0
      cc(5,5) = fact
      cc(5,6) = r0
      cc(6,1) = r0
      cc(6,2) = r0
      cc(6,3) = r0
      cc(6,4) = r0
      cc(6,5) = r0
      cc(6,6) = fact

      return

      end







      subroutine ccFibrePart(cc,c,d,e,M)
      implicit none
      double precision cc(6,*), c, d, e, M(*), ee

      ee = 0.5d0 * e

      cc(1,1) = cc(1,1)+ c *(M(1)+M(1)) + d *M(1)*M(1) + e  *(M(1)+M(1))
      cc(1,2) = cc(1,2)+ c *(M(2)+M(1)) + d *M(1)*M(2)
      cc(1,3) = cc(1,3)+ c *(M(3)+M(1)) + d *M(1)*M(3)
      cc(1,4) = cc(1,4)+ c * M(4)       + d *M(1)*M(4) + e  * M(4)
      cc(1,5) = cc(1,5)+ c * M(5)       + d *M(1)*M(5)
      cc(1,6) = cc(1,6)+ c * M(6)       + d *M(1)*M(6) + e  * M(6)
      cc(2,1) = cc(2,1)+ c *(M(1)+M(2)) + d *M(2)*M(1)
      cc(2,2) = cc(2,2)+ c *(M(2)+M(2)) + d *M(2)*M(2) + e  *(M(2)+M(2))
      cc(2,3) = cc(2,3)+ c *(M(3)+M(2)) + d *M(2)*M(3)
      cc(2,4) = cc(2,4)+ c * M(4)       + d *M(2)*M(4) + e  * M(4)
      cc(2,5) = cc(2,5)+ c * M(5)       + d *M(2)*M(5) + e  * M(5)
      cc(2,6) = cc(2,6)+ c * M(6)       + d *M(2)*M(6)
      cc(3,1) = cc(3,1)+ c *(M(1)+M(3)) + d *M(3)*M(1)
      cc(3,2) = cc(3,2)+ c *(M(2)+M(3)) + d *M(3)*M(2)
      cc(3,3) = cc(3,3)+ c *(M(3)+M(3)) + d *M(3)*M(3) + e  *(M(3)+M(3))
      cc(3,4) = cc(3,4)+ c * M(4)       + d *M(3)*M(4)
      cc(3,5) = cc(3,5)+ c * M(5)       + d *M(3)*M(5) + e  * M(5)
      cc(3,6) = cc(3,6)+ c * M(6)       + d *M(3)*M(6) + e  * M(6)
      cc(4,1) = cc(4,1)+ c *      M(4)  + d *M(4)*M(1) + e  * M(4)
      cc(4,2) = cc(4,2)+ c *      M(4)  + d *M(4)*M(2) + e  * M(4)
      cc(4,3) = cc(4,3)+ c *      M(4)  + d *M(4)*M(3)
      cc(4,4) = cc(4,4)                 + d *M(4)*M(4) + ee *(M(1)+M(2))
      cc(4,5) = cc(4,5)                 + d *M(4)*M(5) + ee * M(6)
      cc(4,6) = cc(4,6)                 + d *M(4)*M(6) + ee * M(5)
      cc(5,1) = cc(5,1)+ c *      M(5)  + d *M(5)*M(1)
      cc(5,2) = cc(5,2)+ c *      M(5)  + d *M(5)*M(2) + e  * M(5)
      cc(5,3) = cc(5,3)+ c *      M(5)  + d *M(5)*M(3) + e  * M(5)
      cc(5,4) = cc(5,4)                 + d *M(5)*M(4) + ee * M(6)
      cc(5,5) = cc(5,5)                 + d *M(5)*M(5) + ee *(M(2)+M(3))
      cc(5,6) = cc(5,6)                 + d *M(5)*M(6) + ee * M(4)
      cc(6,1) = cc(6,1)+ c *      M(6)  + d *M(6)*M(1) + e  * M(6)
      cc(6,2) = cc(6,2)+ c *      M(6)  + d *M(6)*M(2)
      cc(6,3) = cc(6,3)+ c *      M(6)  + d *M(6)*M(3) + e  * M(6)
      cc(6,4) = cc(6,4)                 + d *M(6)*M(4) + ee * M(5)
      cc(6,5) = cc(6,5)                 + d *M(6)*M(5) + ee * M(4)
      cc(6,6) = cc(6,6)                 + d *M(6)*M(6) + ee *(M(3)+M(1))

      return

      end


