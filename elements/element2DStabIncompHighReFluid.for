
      integer function element2DStabIncompHighReFluid
     *                    (elmDat,timDat,xl,ul,dt,s,p,
     *                     ndf,ndm,nen,nelm,isw)

c  ***************************************************************************
c  *                                                                         *
c  *                                                                         *
c  *                                                                         *
c  *  ---------------------------------------------------------------------  *
c  *    parameters and constants                                             *
c  *                                                                         *
c  *    d(1) = 1  x  stabilise momentum equation + discontinuity capturing   *
c  *                 type 1                                                  *
c  *           2  x  stabilise momentum equation + discontinuity capturing   *
c  *                 type 2                                                  *
c  *                                                                         *
c  *    d(2)          not used                                               * 
c  *                                                                         * 
c  *    d(3)      x  beta1 for tau_velocity                                  * 
c  *    d(4)      x  beta2                                                   * 
c  *                                                                         * 
c  *    d(5)      x  beta1 for tau_pressure                                  * 
c  *    d(6)      x  beta2                                                   * 
c  *                                                                         * 
c  *    d(7)      x  not used                                                * 
c  *    d(8)      x  not used                                                * 
c  *                                                                         *
c  *    d(9)      x  not used                                                *
c  *    d(10)     x  not used                                                *
c  *                                                                         *
c  *    d(11)     x  mu                                                      *
c  *    d(12)     x  rho                                                     *
c  *    d(13)     x  f(1)   body force vector                                *
c  *    d(14)     x  f(2)                                                    *
c  *                                                                         * 
c  *    d(15)     x  not used                                                * 
c  *                                                                         * 
c  *  ---------------------------------------------------------------------  * 
c  *        ul(j,i),  j = 1, 2, 3=ndf,  i = 1, 2, 3=nen                      *
c  *                                                                         * 
c  *    j = 1. dof  -  x velocity                                            *
c  *        2. dof  -  y velocity                                            *
c  *        3. dof  -  pressure                                              *
c  *                                                                         * 
c  *    i = 1...3   -  current nodal values    u^k(t)                        * 
c  *        4...6   -  increment  u^k(t) - u(t_n)  .                         *
c  *        7...9   -  previous solution rate      u(t_n)   (if required)    *
c  *                                                                         *
c  *  ---------------------------------------------------------------------  * 
c  *        xl(j,i),  i = 1, 2=ndm,     j = 1, 2, 3=nen                      *
c  *                                                                         * 
c  *    j = 1. dof  -  x                                                     *
c  *        2. dof  -  y                                                     *
c  *                                                                         *
c  *    i = 1..3    -  nodal coordinates  x_n                                *
c  *        4..6    -  nodal coordinates  x_{n+1} - x_n                      *
c  *                                                                         *
c  *  ---------------------------------------------------------------------  * 
c  *   isw                                                                   *
c  *    1   check elements and compute volume                                *
c  *    2 x compute residual                                                 *
c  *    3 x compute residual + tangent for  u_{n+1}                          *
c  *    4 x compute residual + tangent for du_{n+1}                          *
c  *    5 x compute residual + tangent for  u_{n+1}, du_{n+1}(u_{n+1})       *
c  *    6 x not used                                                         *
c  *    7 x not used                                                         *
c  *    8 x not used                                                         *
c  *    9 x not used                                                         *
c  *   10 x not used                                                         *
c  *   11 x not used                                                         *
c  *   12 x not used                                                         *
c  *   13 x not used                                                         *
c  *   14 x not used                                                         *
c  *   15 x not used                                                         *
c  *   16 x not used                                                         *
c  *   17 x not used                                                         *
c  *   18 x not used                                                         *
c  *   19 x not used                                                         *
c  *   20 x not used                                                         *
c  *                                                                         *
c  ***************************************************************************

      implicit none

      integer          ndf, ndm, nen, nelm, isw
      double precision elmDat(*), timDat(*), dt, ul(ndf,*), xl(ndm,*), 
     *                 s(nen*ndf,*), p(*)

      double precision mu, rho, area2, h, nRC(3), z(3), beta(6),norm_u,
     *                 wa(2,3), ua(3,18), dua(3,18), u(3), uw(2), rs,
     *                 gu(2,2), f(2), xa(2,3), urdr, urdrr, rad, rad2, 
     *                 gp(2), du(2), trgu, rstab(2), rstabi(3), tauDC,
     *                 drdui(2,2), drduj(2,2), tau(3), dtauDC(2), rs2,
     *                 dtau(2,3), dtaudx(2,3,2), dtaudh(3), unit_s(2), 
     *                 sg(4,4), shp(3,3), r(3), k(3,3), unit_r(2),
     *                 dvol, xsj, dgu(2,2,2,3), dtrgu(2,3), dgp(2,2,3), 
     *                 darea2(2,3), darea2r(2,3), dhfact, Dshp(2,3,2,3), 
     *                 drdxj(2,2), drduidw(2), fact, hlp, rr(2,2),
     *                 r0, r1, r2, r4, r1d3, r1d2, r1d4, r1d6, r1d216,
     *                 pi, r2pi, p6, p7, p8, b17, ss(2,2), dbdu(3,2),
     *                 dadu(4,2),

     *                 dtemp, b(3), c(2,2), st, a(4)  ! Sony

      integer i, j, i0, j0, ii, jj, l, kk, lint, iord, ie, stab

      logical errck, stab1, stab2

      data r0, r1d4, r1d2, r1, r2, r4 
     *         / 0.d0, 0.25d0, 0.5d0, 1.d0, 2.d0, 4.d0 /

      r1d3  = 1.d0 /  3.d0
      r1d6  = r1d2 * r1d3
      r1d216= r1d6 * r1d6 * r1d6
      pi    = acos(-1.d0)
      r2pi  = r2 * pi
      
c      okflg = .true.

      p6    = timDat(6)
      p7    = timDat(7)
      p8    = timDat(8)

      b17   = timDat(9)

      area2 =  (xl(1,1)-xl(1,3))*(xl(2,2)-xl(2,3)) 
     *       - (xl(1,2)-xl(1,3))*(xl(2,1)-xl(2,3))

      if (area2.lt.1.d-12) then
c        okflg = .false.
        return
      endif

c...  set control paramters and element constants

      stab  = nint(elmDat(1))
      do i=1, 6
        beta(i) = elmDat(2+i)
      enddo
      mu    = elmDat(11)
      rho   = elmDat(12)
      f(1)  = elmDat(13)
      f(2)  = elmDat(14)
      stab1 = stab.eq.1
      stab2 = stab.eq.2

      element2DStabIncompHighReFluid = 0


c            1 2 3 4 5  6  7  8  9 10 11 12 13 14 15 16,17,18,19,20
      goto (14,2,2,2,2,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14),isw

      call prgerror(1,'elmt03fld','what ?')
      
c=============================================================================
 2    continue ! isw 2 -> compute residual
               !     3 -> compute residual + tangent for  u_{n+1}
               !     4 -> compute residual + tangent for du_{n+1}
               !     5 -> compute residual + full tangent for u_{n+1} 
c-----------------------------------------------------------------------------

c                .     
c...  compute u, u, x, w ...........................

      do i=1, 3
        do ii=1, 2
          ua (ii,i) = p7 * ul(ii,i  ) + (r1-p7) * ul(ii,i+3) 
          dua(ii,i) = p8 * ul(ii,i+6) + (r1-p8) * ul(ii,i+9) 
          xa (ii,i) = p6 * xl(ii,i  ) + (r1-p6) * xl(ii,i+3) 
          wa (ii,i) = p7 * xl(ii,i+6) + (r1-p7) * xl(ii,i+9) 
        enddo
        ua(3,i) = ul(3,i) ! pressure
      enddo

c...  characteristic element size
   
      area2 =  (xa(1,1)-xa(1,3))*(xa(2,2)-xa(2,3)) 
     *       - (xa(1,2)-xa(1,3))*(xa(2,1)-xa(2,3))

      h     = sqrt(r2/pi * area2)  ! diameter of circle

c...  evaluate tau and its derivatives at centroid

      u(1) = r1d3 * (ua(1,1)+ua(1,2)+ua(1,3)-wa(1,1)-wa(1,2)-wa(1,3))
      u(2) = r1d3 * (ua(2,1)+ua(2,2)+ua(2,3)-wa(2,1)-wa(2,2)-wa(2,3))
      call eval_tau (u,h,rho,mu,dt,nRC,z,beta,tau)
      call eval_dtau(u,h,rho,mu,dt,nRC,z,beta,dtau,r1d3)

c...  evaluate tauDC and its derivatives at centroid

      fact      = h / r2
      norm_u    = sqrt(u(1) * u(1) + u(2) * u(2))
      tauDC     = fact * norm_u
      if (norm_u.lt.1.0d-16) then
        dtauDC(1) = 0.0
        dtauDC(2) = 0.0
      else
        dtauDC(1) = r1d3 * fact * u(1) / norm_u
        dtauDC(2) = r1d3 * fact * u(2) / norm_u
      endif

c...  get Gauss point coordinates and weights

      l    = -3
      iord = 1
      call tint2d(l,lint,sg)

c...  start loop over Gauss points and get shape functions

      do l=1, lint

        call trishp(sg(1,l),xa,ndm,iord,xsj,shp)

        dvol = xsj * sg(4,l)

c...    compute u, p, grad(u), tr(grad(u)), du/dt

        call pzero(u,3)                                    
        call pzero(uw,2)                                    
        call pzero(gu,4)                                   
        call pzero(gp,2)
        call pzero(du,2)
        do i=1, 3
          do ii=1, 2                                              !
            u(ii)       = u(ii)  +  ua(ii,i)           * shp(3,i) ! u
            uw(ii)      = uw(ii) + (ua(ii,i)-wa(ii,i)) * shp(3,i) ! u - w
            gp(ii)      = gp(ii)    +  ua( 3,i) * shp(ii,i)       ! grad(p)
            du(ii)      = du(ii)    + dua(ii,i) * shp( 3,i)       ! du/dt
            do jj=1, 2                                            !
              gu(ii,jj) = gu(ii,jj) + ua(ii,i) * shp(jj,i)        ! grad(u)
            enddo                                                 !
          enddo                                                   !
          u(3)          = u(3)      + ua(3,i) * shp(3,i)          ! p
        enddo                                                     !
        trgu = gu(1,1) + gu(2,2)                                  ! tr(grad(u))
                                                                  !
c... compute Gauss point values for stabilization

        rstab(1) = rho*(du(1)+gu(1,1)*uw(1)+gu(1,2)*uw(2)-f(1))+gp(1)
        rstab(2) = rho*(du(2)+gu(2,1)*uw(1)+gu(2,2)*uw(2)-f(2))+gp(2)

        if (stab2) then
          call get_b(shp,u,gu,unit_r,unit_s,b)
*         write(*,*)(b(ii), ii = 1, 3)
        endif

c:::::: RESIDUAL VECTOR ::::::::::::::::::::::::::::::::::::::::::::::::::::::

        do i=1, 3
          i0 = (i-1) * ndf

c         time derivative term and convective acceleration, body load
          fact = shp(3,i) * rho
          r(1) = fact * (du(1) + gu(1,1)*uw(1) + gu(1,2)*uw(2) - f(1))
          r(2) = fact * (du(2) + gu(2,1)*uw(1) + gu(2,2)*uw(2) - f(2))

c         viscous part
          r(1) = r(1) + mu * (shp(1,i) * (gu(1,1) + gu(1,1)) +
     &                        shp(2,i) * (gu(1,2) + gu(2,1)))
          r(2) = r(2) + mu * (shp(1,i) * (gu(1,2) + gu(2,1)) +
     &                        shp(2,i) * (gu(2,2) + gu(2,2)))

c         pressure terms
          r(1) = r(1) - u(3) * shp(1,i)
          r(2) = r(2) - u(3) * shp(2,i)
          r(3) = trgu * shp(3,i)

c         stabilization of momentum equation

          call rupwind(shp,uw,gu,tau,rstab,drdui,rstabi,r,i)

c         discontinuity capturing term

          fact = tauDC * rho
          if (stab1) then
            r(1) = r(1) + fact * shp(1,i) * trgu
            r(2) = r(2) + fact * shp(2,i) * trgu
          elseif (stab2) then
            call get_a(b,gu,a)

*           write(*,*)(a(ii), ii = 1, 4)

            r(1) = r(1) + fact * (shp(1,i) * a(1) + shp(2,i) * a(3))
            r(2) = r(2) + fact * (shp(1,i) * a(2) + shp(2,i) * a(4))
          endif

c         assembly

          do ii=1, 3
            p(i0+ii) = p(i0+ii) - r(ii) * dvol
          enddo

c:::::::: STIFFNESS MATRIX FOR u_{n+1} ::::::::::::::::::::::::::::::::::::::

          if (isw.ne.3.and.isw.ne.5) goto 201

          do j=1, 3
            j0 = (j-1) * ndf

            fact = shp(1,j) * uw(1) + shp(2,j) * uw(2)
            drduj(1,1) = rho * (gu(1,1) * shp(3,j) + fact)
            drduj(1,2) = rho *  gu(1,2) * shp(3,j)
            drduj(2,1) = rho *  gu(2,1) * shp(3,j)
            drduj(2,2) = rho * (gu(2,2) * shp(3,j) + fact)

c           time derivative term and convective acceleration
            k(1,1) = shp(3,i) * drduj(1,1)
            k(1,2) = shp(3,i) * drduj(1,2)
            k(2,1) = shp(3,i) * drduj(2,1)
            k(2,2) = shp(3,i) * drduj(2,2)

c           viscous part   ! could be integrated with only one GP
            fact = shp(1,i)*shp(1,j) + shp(2,i)*shp(2,j)
            k(1,1) = k(1,1) + mu * (shp(1,j)*shp(1,i) + fact)
            k(1,2) = k(1,2) + mu *  shp(1,j)*shp(2,i)
            k(2,1) = k(2,1) + mu *  shp(2,j)*shp(1,i)
            k(2,2) = k(2,2) + mu * (shp(2,j)*shp(2,i) + fact)

c           pressure terms   ! could be integrated with only one GP  
            k(1,3) = - shp(3,j) * shp(1,i)
            k(2,3) = - shp(3,j) * shp(2,i)
            k(3,1) = shp(3,i) * shp(1,j)
            k(3,2) = shp(3,i) * shp(2,j)
            k(3,3) = r0

c           stabilization of momentum equation

            call kupwind(shp,tau,dtau,rstab,rstabi,drdui,drduj,k,i,j)

c           discontinuity capturing term

            if (stab1) then
              k(1,1) = k(1,1) +
     &                 rho*shp(1,i)*(dtauDC(1)*trgu+tauDC*shp(1,j))
              k(1,2) = k(1,2) +
     &                 rho*shp(1,i)*(dtauDC(2)*trgu+tauDC*shp(2,j))
              k(2,1) = k(2,1) +
     &                 rho*shp(2,i)*(dtauDC(1)*trgu+tauDC*shp(1,j))
              k(2,2) = k(2,2) +
     &                 rho*shp(2,i)*(dtauDC(2)*trgu+tauDC*shp(2,j))
            elseif (stab2) then

              call get_dbdu(unit_r,unit_s,dbdu)

              call get_dadu(shp,gu,b,dbdu,dadu)

              k(1,1) = k(1,1) +
     &                 rho*(shp(1,i)*(dtauDC(1)*a(1)+tauDC*dadu(1,1))+
     &                      shp(2,i)*(dtauDC(1)*a(3)+tauDC*dadu(3,1)))
              k(1,2) = k(1,2) +
     &                 rho*(shp(1,i)*(dtauDC(2)*a(1)+tauDC*dadu(1,2))+
     &                      shp(2,i)*(dtauDC(2)*a(3)+tauDC*dadu(3,2)))
              k(2,1) = k(2,1) +
     &                 rho*(shp(1,i)*(dtauDC(1)*a(2)+tauDC*dadu(2,1))+
     &                      shp(2,i)*(dtauDC(1)*a(4)+tauDC*dadu(4,1)))
              k(2,2) = k(2,2) +
     &                 rho*(shp(1,i)*(dtauDC(2)*a(2)+tauDC*dadu(2,1))+
     &                      shp(2,i)*(dtauDC(2)*a(4)+tauDC*dadu(4,2)))
            endif

c           assembly

            fact  = dvol * p7
            do ii=1, 3
              do jj=1, 2
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * fact
              enddo
              s(i0+ii,j0+ 3)   = s(i0+ii,j0+ 3) + k(ii, 3) * dvol
            enddo

          enddo

c:::::::: STIFFNESS MATRIX FOR du_{n+1} :::::::::::::::::::::::::::::::::::::

 201      if (isw.ne.4.and.isw.ne.5) goto 202

          k(1,2) = r0
          k(2,1) = r0
          k(3,1) = r0
          k(3,2) = r0

          do j=1, 3
            j0 = (j-1) * ndf

c           time derivative term and convective acceleration
            fact = shp(3,i) * rho * shp(3,j)
            k(1,1) = fact
            k(2,2) = fact

c           stabilization of momentum equation

            fact   = tau(1) * (shp(1,i)*uw(1) + shp(2,i)*uw(2)) *
     &               rho * shp(3,j)
            k(1,1) = k(1,1) + fact
            k(2,2) = k(2,2) + fact
            fact   = tau(2) * rho * shp(3,j)
            k(3,1) = fact * shp(1,i)
            k(3,2) = fact * shp(2,i)

c           assembly
            fact  = dvol * p8
            if (isw.eq.5) fact = fact * b17
            do ii=1, 3
              do jj=1, 2
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * fact
              enddo
            enddo

          enddo

 202      continue

        enddo

      enddo

      return
c==============================================================================
 14   continue ! isw 14-15 -> not used
c-----------------------------------------------------------------------------

      write(*,'(/9x,a/)') 'elmt03fld: what do you want ?'

      return

cc=============================================================================
      end
















