
      integer function ale2D3nodedTriangleHyperElasticity
     *                    (d,xl,s,p,ndm,nst,ndat)
      implicit none

      integer          ndm, nst, ndat,
     *                 i, j, m, ii, jj, i0, j0
      double precision d(*), xl(ndm,*), s(nst,*), p(*),
     *                 x(2,5), x0(2,5), u(2,5), mu, bulk, shp(3,3), 
     *                 fact, r(2), k(2,2), area, detF, F(2,2), sig(4), 
     *                 cc(4,4), dmy
c
c     LINEAR TRIANGLES,
c
c     MOVE MESH BY SOLVING PSEUDO-HYPERELASTIC PROBLEM
c 
c
c     d(1) = mu
c     d(2) = K
c     d(3) = 0 -> x0 is reference config.
c            1 -> x1 is reference config.

      if (ndat.lt.3) 
     *  call prgerror(1,'ale2D3nodedTriangleHyperElasticity',
     *                  '3 parameters required!')

      ale2D3nodedTriangleHyperElasticity = 0

c...  current configuration
      do i=1, 3
        do j=1, ndm
          x(j,i) = xl(j,i) + xl(j,i+3)
        enddo
      enddo

c...  reference configuration and current displacement
      if (nint(d(3)).ne.0) then
        do i=1, 3
          do j=1, ndm
            x0(j,i) = xl(j,i)
            u (j,i) = xl(j,i+3)
          enddo
        enddo
      else
        do i=1, 3
          do j=1, ndm
            x0(j,i) = xl(j,i+6)
            u (j,i) = x(j,i) - x0(j,i)
          enddo
        enddo

      endif

      fact =  (x(1,1)-x(1,3))*(x(2,2)-x(2,3)) 
     *      - (x(1,2)-x(1,3))*(x(2,1)-x(2,3))

      if (fact.lt.0.d0) then
        ale2D3nodedTriangleHyperElasticity = -1
        return
      endif

      area = .5d0 * fact

c...  derivatives of shape functions
      do i=1, 3 
        j = i+1
        m = i+2
        if (j.gt.3) j = j - 3
        if (m.gt.3) m = m - 3
        fact     = 1.d0 / ( (x(2,i) - x(2,j)) * (x(1,j) - x(1,m))
     *                     -(x(2,j) - x(2,m)) * (x(1,i) - x(1,j)))
        shp(1,i) = - fact * (x(2,j) - x(2,m))
        shp(2,i) =   fact * (x(1,j) - x(1,m))
      enddo

c...  compute F, sig, c
      call f2d(shp,u,F,detF,2,3,2)
      call mat2d01(d,F,dmy,sig,cc,dmy,dmy,.true.,2,2)

c      do i=1, 4
c        write(*,'(6(1x,f6.4))') (cc(i,j),j=1,4)
c      enddo
c      write(*,*)


      do i=1, 3
        i0 = (i-1)*2

c...    residual
        r(1) = shp(1,i)*sig(1) + shp(2,i)*sig(3)
        r(2) = shp(1,i)*sig(3) + shp(2,i)*sig(2)
        do ii=1, 2
          p(i0+ii) = p(i0+ii) - r(ii) * area
        enddo

        do j=1, 3
          j0 = (j-1)*2

c...      material and geometrical stiffness
          fact  =   shp(1,i) * sig(1) * shp(1,j)
     *            + shp(1,i) * sig(3) * shp(2,j)
     *            + shp(2,i) * sig(3) * shp(1,j)
     *            + shp(2,i) * sig(2) * shp(2,j)
          k(1,1) =  shp(1,i) * cc(1,1) * shp(1,j)
     *            + shp(1,i) * cc(1,3) * shp(2,j)
     *            + shp(2,i) * cc(3,1) * shp(1,j)
     *            + shp(2,i) * cc(3,3) * shp(2,j) + fact
          k(1,2) =  shp(1,i) * cc(1,3) * shp(1,j)
     *            + shp(1,i) * cc(1,2) * shp(2,j)
     *            + shp(2,i) * cc(3,3) * shp(1,j)
     *            + shp(2,i) * cc(3,2) * shp(2,j)
          k(2,1) =  shp(1,i) * cc(3,1) * shp(1,j)
     *            + shp(1,i) * cc(3,3) * shp(2,j)
     *            + shp(2,i) * cc(2,1) * shp(1,j)
     *            + shp(2,i) * cc(2,3) * shp(2,j)
          k(2,2) =  shp(1,i) * cc(3,3) * shp(1,j)
     *            + shp(1,i) * cc(3,2) * shp(2,j)
     *            + shp(2,i) * cc(2,3) * shp(1,j)
     *            + shp(2,i) * cc(2,2) * shp(2,j) + fact
          do ii=1, 2
            do jj=1, 2
              s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * area
            enddo
          enddo

        enddo

      enddo

c      do i=1, 6
c        write(*,'(6(1x,f6.4),3x,g12.5)') (s(i,j),j=1,6), p(i)
c      enddo
c      write(*,*)

      return
      end









      subroutine mat2d01(d,F,F33,stre,cc,iv1,iv2,finite,ss,isw)
      implicit none

      integer          ss, isw, i, j

      double precision d(*), F(2,*), F33, stre(*), cc(4,*), iv1(*), 
     *                 iv2(*),
     *                 fact, fact1, fact2, fact3, mu, K, Lamb, eps(3),
     *                 detF, b(3),
     *                 r1d2, r1, r2

      logical          finite

c
c     2D MATERIAL SUBROUTINE
c
c     Neo-Hookean elasticity  (Wriggers, page 45)
c     -------------------------------------------
c
c     isw = 1 -> initialize internal variables 
c
c     isw = 2 -> compute Cauchy stress, update internal variables and
c                compute consistent spatial tangent tensor
c
c     ss  = 1 -> plane stress
c     ss  = 2 -> plane strain
c     ss  = 3 -> axisymmetric
c    (ss  = 4 -> plane stress, incompressible
c     ss  = 5 -> plane strain, incompressible
c     ss  = 6 -> axisymmetric, incompressible)
c
c     stre: 11  22  12  33
c
c     cc:  1111 1122 1112 1133
c          2211 2222 2212 2233 
c          1211 1222 1212 1233
c          3311 3322 3312 3333
c
      data r1d2, r1, r2 / .5d0, 1.d0, 2.d0 /

      if (isw.eq.1) return ! no internal variables

      mu     = d(1)
      K      = d(2)
      Lamb   = K - 2.d0/3.d0 * mu

      if (.not.finite) goto 100

c...  finite strains .........

      b(1) = F(1,1) * F(1,1) + F(1,2) * F(1,2)
      b(2) = F(2,2) * F(2,2) + F(2,1) * F(2,1)
      b(3) = F(2,1) * F(1,1) + F(1,2) * F(2,2)

      detF = F(1,1) * F(2,2) - F(1,2) * F(2,1)

      if (ss.le.1) then     ! plane stress

        F33    = sqrt((Lamb+r2*mu)/(Lamb*detF*detF+r2*mu))
        fact2  =Lamb*detF*F33*(r1-detF*detF*Lamb/(Lamb*detF*detF+r2*mu))
        detF   = detF * F33
        fact1  = mu / detF
        fact   = Lamb / (r2 * detF) * (detF*detF-r1) - mu / detF

        stre(1) = fact1 * b(1) + fact
        stre(2) = fact1 * b(2) + fact
        stre(3) = fact1 * b(3)

        call pzero(cc,16)
        fact = - fact * r2
        do i=1, 2
          do j=1, 2
            cc(i,j) = fact2
          enddo
          cc(i,i) = cc(i,i) + fact
        enddo
        cc(3,3) = r1d2 * fact

      elseif (ss.eq.2) then ! plane strain

        fact2  = detF * Lamb
        fact   = Lamb / (r2 * detF) * (detF*detF-r1) - mu / detF
        fact1  = mu / detF

        stre(1) = fact1 * b(1) + fact
        stre(2) = fact1 * b(2) + fact
        stre(3) = fact1 * b(3)

        call pzero(cc,16)

        fact = - fact * r2
        do i=1, 2
          do j=1, 2
            cc(i,j) = fact2
          enddo
          cc(i,i) = cc(i,i) + fact
        enddo
        cc(3,3) = r1d2 * fact

      elseif (ss.ge.3) then ! axisymmetric

        call prgerror(1,'mat2d01/large','option not yet available')

      endif

      return

 100  continue

c...  small strains .........

      eps(1) = F(1,1) - r1
      eps(2) = F(2,2) - r1
      eps(3) = r1d2 * (F(1,2) + F(2,1))

      if (ss.le.1) then     ! plane stress

        fact1 = Lamb * (r1 - Lamb/(Lamb+r2*mu)) + r2*mu
        fact2 = Lamb * (r1 - Lamb/(Lamb+r2*mu))
        fact3 = r2*mu

      elseif (ss.eq.2) then ! plane strain

        fact1 = Lamb + r2*mu
        fact2 = Lamb
        fact3 = r2*mu

      elseif (ss.ge.3) then ! axisymmetric

        call prgerror(1,'mat2d01/small','option not yet available')

      endif

      stre(1) = fact1 * eps(1) + fact2 * eps(2)
      stre(2) = fact2 * eps(1) + fact1 * eps(2)
      stre(3) = fact3 * eps(3)

      call pzero(cc,16)
      cc(1,1) = fact1
      cc(2,2) = fact1
      cc(1,2) = fact2
      cc(2,1) = cc(1,2)
      cc(3,3) = r1d2 * fact3 

      return

      end



