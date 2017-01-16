
      integer function element1DSolidALE(elmDat,timDat,xl,ul,s,p)

      implicit none

      integer          i, j, ngp, ii, jj, l, i0, j0,
     *                 nen, ndf, 
     *                 round

      double precision elmDat(*), timDat(*), 
     *                 xl(*), ul(*), s(6,*), p(*),
     *                 ua(2), x0a(2), xa(2), dua(2), dx0a(2), va(2),
     *                 shp(2,2), wgp(10), xi(10), r(3), k(3,3),
     *                 lam, lam2, sig, C, dvol, rho, dl, mu, f, diam2, 
     *                 h, u, du, x0, dx0, gu, gx0, uv, detJ, ddl(2),
     *                 ddvol(3,2), tau,
     *                 fact, p6, p7, p8, td1, td2, dshp(2,2), h2,
     *                 r1d2, r1, r2, r3

      data r1d2, r1, r2, r3 / 0.5d0, 1.d0, 2.d0, 3.d0 /

      ngp  = round(elmDat(1))

      nen = 2
      ndf = 3

      diam2 = elmDat(2) * elmDat(2)

      f    = elmDat(3)
      rho  = elmDat(4)

      mu   = elmDat(5)

      p6   = timDat(6)
      p7   = timDat(7)
      p8   = timDat(8)

      td1  = timDat(11)   ! d du / d u
      td2  = timDat(31)   ! d x0 / d dx0   and    d x / d v

      element1DSolidALE = 0

c                .      .
c...  compute u, u, x0, x0, x, v ...........................

      do i=1, nen
        ua  (i) = p7 * ul(-2+i*ndf) + (r1-p7) * ul( 4+i*ndf) 
        dx0a(i) = p8 * ul(-1+i*ndf) + (r1-p8) * ul( 5+i*ndf)
        va  (i) = p8 * ul(   i*ndf) + (r1-p8) * ul( 6+i*ndf)
        dua (i) = p8 * ul(10+i*ndf) + (r1-p8) * ul(16+i*ndf)
        x0a (i) = p7 * ul(11+i*ndf) + (r1-p7) * ul(17+i*ndf)
        xa  (i) = p7 * ul(12+i*ndf) + (r1-p7) * ul(18+i*ndf)
      enddo

      h = xa(2) - xa(1)  ! element length

      h2 = h * h

c...  get Gauss point locations and weights

      call comp_xigp1D(xi,wgp,ngp)

c...  ALE shape function derivatives

      dshp(1,1) = - p7 / h2 * td2
      dshp(1,2) = - dshp(1,1)

      dshp(2,1) = - dshp(1,1)
      dshp(2,2) = + dshp(1,1)

c...  start loop over Gauss points and get shape functions

      do l=1, ngp

        call comp_shp1D(shp,detJ,xi(l),xa,nen)

c        write(*,*) xi(l), wgp(l), detJ

        dl = detJ * wgp(l)

        if (ngp.eq.2) then
          ddl(1) = - .5d0 * p7 * td2 * wgp(l)
          ddl(2) = - ddl(1)
        else 
          ddl(1) = - .5d0 * p7 * td2 * wgp(l)
          ddl(2) = - ddl(1)
        endif

c...    compute u, u-v, x0, du/dt, dx0/dt, grad(u), grad0(d)->lam

        u   = ua(1) * shp(2,1) + ua(2) * shp(2,2)

        uv  = (ua(1)-va(1)) * shp(2,1) + (ua(2)-va(2)) * shp(2,2)

        x0  = x0a(1) * shp(2,1) + x0a(2) * shp(2,2)

        du  = dua(1) * shp(2,1) + dua(2) * shp(2,2)

        dx0 = dx0a(1) * shp(2,1) + dx0a(2) * shp(2,2)

        gu  = ua(1) * shp(1,1) + ua(2) * shp(1,2)

        gx0 = x0a(1) * shp(1,1) + x0a(2) * shp(1,2)
  
        lam = r1 / (x0a(1) * shp(1,1) + x0a(2) * shp(1,2))

!        if (l.eq.1) write(*,*) u!, uv, x0, du, dx0, gu, gx0, lam

c...    calculate stress

c       incompressible Neo-Hooke (see my PhD)

        lam2 = lam * lam

        dvol = dl * diam2 / lam

        ddvol(2,1) = dl * diam2 * shp(1,1) * p7 * td2
        ddvol(2,2) = dl * diam2 * shp(1,2) * p7 * td2

        ddvol(3,1) = diam2 * (ddl(1) / lam + dl * (dshp(1,1) * x0a(1) +
     *                                             dshp(2,1) * x0a(2)))
        ddvol(3,2) = diam2 * (ddl(2) / lam + dl * (dshp(1,2) * x0a(1) +
     *                                             dshp(2,2) * x0a(2)))

        sig = mu * (lam2*lam2 - r1) / lam2

        !write(*,*) lam, sig

c       d sig / d lam

        C = mu * (r2 + r2/(lam2*lam2)) * lam

c:::::: RESIDUAL VECTOR ::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !write(*,*) f * rho * dvol

        do i=1, nen
          i0 = (i-1) * ndf

c         standard Galerkin terms

          r(1) = shp(2,i) * rho * (du + gu*uv - f) + shp(1,i) * sig

          r(2) = shp(2,i) * (dx0 + gx0*uv)

          if (i.eq.1) then 
            r(3) = h
          else
            r(3) = -h
          endif

c         least squares stabilisation terms

          tau = 0.

          r(1) = r(1) + tau * (rho * (du + gu*uv - f))
     *        * (rho * (shp(2,i) * p8 * td1 + shp(1,i) * p7 * uv))


          tau = 0.

          r(2) = r(2) + tau * (dx0 + gx0*uv)
     *        * (shp(2,i) * p8 + shp(1,i) * p7 * td2 * uv)


c         assembly

          p(i0+1) = p(i0+1) - r(1) * dvol
          p(i0+2) = p(i0+2) - r(2) * dl
          p(i0+3) = p(i0+3) - r(3)

c:::::::: STIFFNESS MATRIX FOR u_{n+1} ::::::::::::::::::::::::::::::::::::::

          do j=1, 2
            j0 = (j-1) * ndf

c           time derivative term and convective acceleration

            k(1,1) = shp(2,i) * rho * (shp(2,j) * p8 * td1
     *               + (gu*shp(2,j) + shp(1,j)*uv) * p7)

            k(1,3) = shp(2,i) * rho * (uv
     *               * (dshp(1,j) * ua(1) + dshp(2,j) * ua(2))
     *               -  gu * shp(2,j) * p8)

c           stress

            k(1,2) = - shp(1,i) * C * lam2 * shp(1,j) * p7 * td2

            k(1,3) = k(1,3) + dshp(i,j) * sig - shp(1,i) * C
     *               * lam2*(x0a(1)*dshp(1,j)+x0a(2)*dshp(2,j))

c           convection of x0

            k(2,1) = shp(2,i) * gx0 * shp(2,j) * p7

            k(2,2) = shp(2,i) * (shp(2,j) * p8 + shp(1,j)*uv*p7*td2)
 
            k(2,3) = shp(2,i) * ((dshp(1,j)*x0a(1)+dshp(2,j)*x0a(2))
     *                  * uv - gx0 * shp(2,j) * p8)

c           equal spacing

           !k(3,1) = 0.

           !k(3,2) = 0.

            k(3,3) = p7 * td2
            if (i.eq.1) k(3,3) = - k(3,3)
            if (j.eq.2) k(3,3) = - k(3,3)

c           assembly

            s(i0+1,j0+1) = s(i0+1,j0+1) + k(1,1) * dvol
            s(i0+1,j0+2) = s(i0+1,j0+2) + k(1,2) * dvol
     *                                  + r(1) * ddvol(2,j)
            s(i0+1,j0+3) = s(i0+1,j0+3) + k(1,3) * dvol
     *                                  + r(1) * ddvol(3,j)

            s(i0+2,j0+1) = s(i0+2,j0+1) + k(2,1) * dl
            s(i0+2,j0+2) = s(i0+2,j0+2) + k(2,2) * dl
            s(i0+2,j0+3) = s(i0+2,j0+3) + k(2,3) * dl
     *                                  + r(2) * ddl(j)

           !s(i0+3,j0+1) = s(i0+3,j0+1) + k(3,1)
           !s(i0+3,j0+2) = s(i0+3,j0+2) + k(3,2)
            s(i0+3,j0+3) = s(i0+3,j0+3) + k(3,3)

          enddo

        enddo

      enddo

c     stop

      return

      end



