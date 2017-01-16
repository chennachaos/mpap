

      subroutine JacobianForQuad4AndItsDerivatives(detJ,dJ,ddJ,x,xi)

      implicit none

      double precision detJ, dJ(8), ddJ(8,8), x(2,4), xi(2),
     *                 dshp(2,4), Jmtx(2,2),
     *                 r1d4, r1

      integer          i

      data r1d4, r1 / 0.25d0, 1.d0 /

      ! shape function derivatives

      dshp(1,1) = - r1d4 * (r1 - xi(2))
      dshp(2,1) = - r1d4 * (r1 - xi(1))
      dshp(1,2) = + r1d4 * (r1 - xi(2))
      dshp(2,2) = - r1d4 * (r1 + xi(1))
      dshp(1,3) = + r1d4 * (r1 + xi(2))
      dshp(2,3) = + r1d4 * (r1 + xi(1))
      dshp(1,4) = - r1d4 * (r1 + xi(2))
      dshp(2,4) = + r1d4 * (r1 - xi(1))

      ! Jacobian matrix

      call pzero(Jmtx,4)

      do i=1, 4
        Jmtx(1,1) = Jmtx(1,1) + dshp(1,i) * x(1,i)
        Jmtx(1,2) = Jmtx(1,2) + dshp(2,i) * x(1,i)
        Jmtx(2,1) = Jmtx(2,1) + dshp(1,i) * x(2,i)
        Jmtx(2,2) = Jmtx(2,2) + dshp(2,i) * x(2,i)
      enddo

      ! Jacobian determinant

      detJ = Jmtx(1,1) * Jmtx(2,2) - Jmtx(1,2) * Jmtx(2,1)

      ! first derivatives of determinant

      do i=1, 4
        dJ(i+i-1) = dshp(1,i) * Jmtx(2,2) - dshp(2,i) * Jmtx(2,1)
        dJ(i+i  ) = dshp(2,i) * Jmtx(1,1) - dshp(1,i) * Jmtx(1,2)
      enddo

      ! second derivatives of determinant

      call pzero(ddJ,64)

      ddJ(1,4) = dshp(1,1) * dshp(2,2) - dshp(2,1) * dshp(1,2)
      ddJ(1,6) = dshp(1,1) * dshp(2,3) - dshp(2,1) * dshp(1,3)
      ddJ(1,8) = dshp(1,1) * dshp(2,4) - dshp(2,1) * dshp(1,4)

      ddJ(2,3) = - ddJ(1,4)
      ddJ(2,5) = - ddJ(1,6)
      ddJ(2,7) = - ddJ(1,8)

      ddJ(3,2) = dshp(1,2) * dshp(2,1) - dshp(2,2) * dshp(1,1)
      ddJ(3,6) = dshp(1,2) * dshp(2,3) - dshp(2,2) * dshp(1,3)
      ddJ(3,8) = dshp(1,2) * dshp(2,4) - dshp(2,2) * dshp(1,4)

      ddJ(4,1) = - ddJ(3,2)
      ddJ(4,5) = - ddJ(3,6)
      ddJ(4,7) = - ddJ(3,8)

      ddJ(5,2) = dshp(1,3) * dshp(2,1) - dshp(2,3) * dshp(1,1)
      ddJ(5,4) = dshp(1,3) * dshp(2,2) - dshp(2,3) * dshp(1,2)
      ddJ(5,8) = dshp(1,3) * dshp(2,4) - dshp(2,3) * dshp(1,4)

      ddJ(6,1) = - ddJ(5,2)
      ddJ(6,3) = - ddJ(5,4)
      ddJ(6,7) = - ddJ(5,8)

      ddJ(7,2) = dshp(1,4) * dshp(2,1) - dshp(2,4) * dshp(1,1)
      ddJ(7,4) = dshp(1,4) * dshp(2,2) - dshp(2,4) * dshp(1,2)
      ddJ(7,6) = dshp(1,4) * dshp(2,3) - dshp(2,4) * dshp(1,3)

      ddJ(8,1) = - ddJ(7,2)
      ddJ(8,3) = - ddJ(7,4)
      ddJ(8,5) = - ddJ(7,6)

      return

      end









      subroutine quadOptShape(x,p,s)

      implicit none

      double precision x(2,4), p(8), s(8,8),
     *                 J, dJ(8), ddJ(8,8), Jc, dJc(8), ddJc(8,8), 
     *                 xi(10), fact, fact1, fact2, Jc2, Jc3, dfact1dJ,
     *                 dfact1dJc, dfact2dJ, dfact2dJc,
     *                 r1, r2, r3

      integer i, k, l

      data r1, r2, r3 / 1.d0, 2.d0, 3.d0 /

      fact    = r1 / sqrt(3.d0)
      xi(1) = - fact
      xi(2) = - fact
      xi(3) = + fact
      xi(4) = - fact
      xi(5) = + fact
      xi(6) = + fact
      xi(7) = - fact
      xi(8) = + fact

      xi( 9) = 0.d0
      xi(10) = 0.d0


      call pzero(p,8)
      call pzero(s,64)

      call JacobianForQuad4AndItsDerivatives(Jc,dJc,ddJc,x,xi(9))

      Jc2 = Jc * Jc
      Jc3 = Jc * Jc2

      write(*,*) " Jc = ", Jc

      do l=1, 4

        call JacobianForQuad4AndItsDerivatives(J,dJ,ddJ,x,xi(l+l-1))

        write(*,'(a,f10.5,a,f10.5)') " J = ", J, " W = ", ((J-Jc)/Jc)**2

        fact1 = (J - Jc) / Jc2
        fact2 = - fact1 * J / Jc

        do i=1, 8
          p(i) = p(i) + fact1 * dJ(i) + fact2 * dJc(i)
        enddo

        fact = r2 * J / Jc3

        dfact1dJ  = r1 / Jc2
        dfact1dJc = dfact1dJ - fact
        dfact2dJ  = dfact1dJc
        dfact2dJc = r3 * J * J / (Jc2 * Jc2) - fact

        do i=1, 8
          do k = 1, 8
            s(i,k) = s(i,k) - fact1 * ddJ(i,k) - fact2 * ddJc(i,k)
     *                      - dfact1dJ  * dJ(i)  * dJ(k)
     *                      - dfact1dJc * dJ(i)  * dJc(k)
     *                      - dfact2dJ  * dJc(i) * dJ(k)
     *                      - dfact2dJc * dJc(i) * dJc(k)
          enddo
        enddo

      enddo

      return

      end
