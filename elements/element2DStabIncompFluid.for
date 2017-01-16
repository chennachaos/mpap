
      integer function element2DStabIncompFluid
     *                    (elmDat,timDat,xl,ul,dt,s,p,
     *                     ndf,ndm,nen,nelm,isw)

c  ***************************************************************************
c  *                                                                         *
c  *    LINEAR TRIANGULAR 2D STABILIZED VELOCITY-PRESSURE FLUID ELEMENT      *
c  *                                                                         *
c  *                  INCOMPRESSIBLE NEWTONIAN FLUID                         *
c  *                                                                         *
c  *              LAGRANGIAN / ALE / EULERIAN  FORMULATION                   *
c  *                                                                         *
c  *         STEADY / UNSTEADY STATE (DISCRETE TIME STEPPING SCHEME)         *
c  *                                                                         *
c  *  ---------------------------------------------------------------------  *
c  *    paramter and constants                                               *
c  *                                                                         *
c  *    d(1) = 0   x  standard Galerkin formulation (does not make sense)    *
c  *           1   x  stabilize momentum equation                            *
c  *           2   x  stabilize continuity equation                          *
c  *           3   x  stabilize momentum and continuity equation             *
c  *                                                                         *
c  *    d(2)          not used                                               * 
c  *                                                                         * 
c  *    d(3)       x  beta1 for tau_velocity                                 * 
c  *    d(4)       x  beta2                                                  * 
c  *                                                                         * 
c  *    d(5)       x  beta1 for tau_pressure                                 * 
c  *    d(6)       x  beta2                                                  * 
c  *                                                                         * 
c  *    d(7)       x  beta1 for tau_continuity                               * 
c  *    d(8)       x  beta2                                                  * 
c  *                                                                         *
c  *    d(9)          not used                                               *
c  *                                                                         *
c  *    d(10)  0,2 x  2D                                                     *
c  *           1   x  axisymmetric                                           *
c  *                                                                         *
c  *    d(11)         mu                                                     *
c  *    d(12)         rho                                                    *
c  *    d(13)         f(1)   body force vector                               *
c  *    d(14)         f(2)                                                   *
c  *                                                                         * 
c  *    d(15)         not used                                               * 
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
c  *    6 x compute residual + tangent for  x_{n+1}                          *
c  *    7 x compute            tangent for dx_{n+1}                          *
c  *    8 x compute residual + tangent for  x_{n+1}, dx_{n+1}(x_{n+1}), Lagr *
c  *    9 x compute residual + tangent for  x_{n+1}, dx_{n+1}(x_{n+1}), ALE  *
c  *   10 x compute volume, kinetic and potential energy,                    *
c  *        kinetic + potential energy, vorticity and div(u)                 *
c  *   11 x project to nodes, vort(u)                                        *
c  *   12 x project to nodes, div(u)                                         *
c  *   13 x project to nodes, ||grad(u)||                                    *
c  *   14   not used                                                         *
c  *   15   not used                                                         *
c  *                                                                         *
c  *        for error indication we need the following stuff:                *
c  *   16 x evaluate element aspect ratio                                    *
c  *   17 x project to nodes element size                                    *
c  *   18 x project to nodes, flow complexity                                *
c  *   19   not used                                                         *
c  *   20   not used                                                         *
c  *                                                                         *
c  ***************************************************************************

      implicit none

      integer          ndf, ndm, nen, nelm, isw
      double precision elmDat(*), timDat(*), dt, ul(ndf,*), xl(ndm,*), 
     *                 s(nen*ndf,*), p(*)

      double precision mu, rho, area2, h, nRC(3), z(3), beta(6),
     *                 wa(2,3), ua(3,18), dua(3,18), u(3), uw(2), 
     *                 gu(2,2), f(2), xa(2,3), urdr, urdrr, rad, rad2, 
     *                 gp(2), du(2), trgu, rstab(2), rstabi(3), 
     *                 drdui(2,2), drduj(2,2), tau(3),  
     *                 dtau(2,3), dtaudx(2,3,2), dtaudh(3),  
     *                 sg(4,4), shp(3,3), r(3), k(3,3), 
     *                 dvol, xsj, dgu(2,2,2,3), dtrgu(2,3), dgp(2,2,3), 
     *                 darea2(2,3), darea2r(2,3), dhfact, Dshp(2,3,2,3), 
     *                 drdxj(2,2), drduidw(2), fact, hlp,
     *                 r0, r1, r2, r4, r1d3, r1d2, r1d4, r1d6, r1d216,
     *                 pi, r2pi, p6, p7, p8, b17, b33, bALE,

     *                 dtemp, b, c, st, a(3)  ! Sony

      integer i, j, i0, j0, ii, jj, l, kk, lint, iord, ie, stab

      logical errck, stabm, stabc, axsy

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
      b33   = timDat(21)  ! ??????????????
      bALE  = timDat(21)

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
      axsy  = nint(elmDat(10)).ne.0
      mu    = elmDat(11)
      rho   = elmDat(12)
      f(1)  = elmDat(13)
      f(2)  = elmDat(14)
      stabm = stab.ne.0.and.stab.ne.2
      stabc = stab.ge.2

      element2DStabIncompFluid = 0

c            1  2 3 4 5 6 7 8 9 10 11 12 13 14 15 16,17,18,19,20
      go to (14,2,2,2,2,6,7,6,6,14,11,11,11,14,14,14,14,14,14,14),isw

      call prgerror(1,'elmt03fld','what ?')
      
c=============================================================================
 2    continue ! isw 2 -> compute residual
               !     3 -> compute residual + tangent for  u_{n+1}
               !     4 -> compute residual + tangent for du_{n+1}
               !     5 -> compute residual + full tangent for u_{n+1} 
c-----------------------------------------------------------------------------

c               .     
c...  compute u, u, x, w ...........................

      do i=1, 3
        do ii=1, 2
          ua (ii,i) = p7 * ul(ii,i  ) + (r1-p7) * ul(ii,i+3) 
          dua(ii,i) = p8 * ul(ii,i+6) + (r1-p8) * ul(ii,i+9) 
          xa (ii,i) = p6 * xl(ii,i  ) + (r1-p6) * xl(ii,i+3) 
          wa (ii,i) = p7 * xl(ii,i+6) + (r1-p7) * xl(ii,i+9) 

          !write(*,*) wa(ii,i)

        enddo
        ua(3,i) = ul(3,i) ! pressure
      enddo

c...  characteristic element size

      area2 =  (xa(1,1)-xa(1,3))*(xa(2,2)-xa(2,3)) 
     *       - (xa(1,2)-xa(1,3))*(xa(2,1)-xa(2,3))

      if (area2 .lt. 0) then
        write(*,100) isw
        return
      endif

      h     = sqrt(r2/pi * area2)  ! diameter of circle

c...  evaluate tau and its derivatives at centroid

      u(1) = r1d3 * (ua(1,1)+ua(1,2)+ua(1,3)-wa(1,1)-wa(1,2)-wa(1,3))
      u(2) = r1d3 * (ua(2,1)+ua(2,2)+ua(2,3)-wa(2,1)-wa(2,2)-wa(2,3))
      call eval_tau (u,h,rho,mu,dt,nRC,z,beta,tau)
      call eval_dtau(u,h,rho,mu,dt,nRC,z,beta,dtau,r1d3)

c...  get Gauss point coordinates and weights

      l    = -3
      iord = 1
      call tint2d(l,lint,sg)

c...  start loop over Gauss points and get shape functions

      do l=1, lint

        call trishp(sg(1,l),xa,ndm,iord,xsj,shp)

        dvol = xsj * sg(4,l)

c...    compute radius r, u_r / r (if axisymmetric)

        if (axsy) then
          rad  = r0
          urdr = r0
          do i=1, 3
            rad  = rad  + shp(3,i) * xa(1,i) 
            urdr = urdr + shp(3,i) * ua(1,i)
          enddo
          rad2  = rad * rad
          urdr  = urdr / rad
          urdrr = urdr / rad
          dvol  = dvol * rad * r2pi
        endif

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

        if (stabm) then
          rstab(1) = rho*(du(1)+gu(1,1)*uw(1)+gu(1,2)*uw(2)-f(1))+gp(1)
          rstab(2) = rho*(du(2)+gu(2,1)*uw(1)+gu(2,2)*uw(2)-f(2))+gp(2)
        endif

c:::::: RESIDUAL VECTOR ::::::::::::::::::::::::::::::::::::::::::::::::::::::

        do i=1, 3
          i0 = (i-1) * ndf

c         time derivative term and convective acceleration, body load
          fact = shp(3,i) * rho
          r(1) = fact * (du(1) + gu(1,1)*uw(1) + gu(1,2)*uw(2) - f(1))
          r(2) = fact * (du(2) + gu(2,1)*uw(1) + gu(2,2)*uw(2) - f(2))

c         viscous part
          r(1) = r(1) + mu * (shp(1,i) * (gu(1,1) + gu(1,1))
     *                      + shp(2,i) * (gu(1,2) + gu(2,1)))
          r(2) = r(2) + mu * (shp(1,i) * (gu(1,2) + gu(2,1))
     *                      + shp(2,i) * (gu(2,2) + gu(2,2)))

c         pressure terms
          r(1) = r(1) - u(3) * shp(1,i)
          r(2) = r(2) - u(3) * shp(2,i)
          r(3) = trgu * shp(3,i)

c         additional terms for axisymmetric case 
          if (axsy) then
            r(1) = r(1) + mu * shp(3,i) * (urdrr + urdrr) ! viscous term
     *                  - u(3) * shp(3,i) / rad           ! pressure term
            r(3) = r(3) + urdr * shp(3,i)                 ! pressure term
          endif

c         stabilization of momentum equation
          if (stabm) call rupwind(shp,uw,gu,tau,rstab,drdui,rstabi,r,i)

c         stabilization of continuity equation
          if (stabc) then
            if (.not.axsy) then              ! 2D
              fact = tau(3) * rho * trgu
              r(1) = r(1) + fact * shp(1,i)
              r(2) = r(2) + fact * shp(2,i)
            else                             ! axisymmetric
              fact = tau(3) * rho * (trgu + urdr)
              r(1) = r(1) + fact * (shp(1,i) + shp(3,i) / rad)
              r(2) = r(2) + fact * shp(2,i)
            endif
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

c           additional terms for axisymmetric case 
            if (axsy) then
              k(1,1) = k(1,1) +mu*shp(3,i)*(shp(3,j)+shp(3,j))/rad2 ! viscosity 
              k(1,3) = k(1,3) - shp(3,i) * shp(3,j) / rad           ! pressure
              k(3,1) = k(3,1) + shp(3,i) * shp(3,j) / rad           ! pressure
            endif

c           stabilization of momentum equation
            if (stabm) 
     *        call kupwind(shp,tau,dtau,rstab,rstabi,drdui,drduj,k,i,j)

c           stabilization of continuity equation          
            if (stabc) then
              if (.not.axsy) then              ! 2D
                fact   = tau(3) * rho
                k(1,1) = k(1,1) + fact * shp(1,i) * shp(1,j)
                k(1,2) = k(1,2) + fact * shp(1,i) * shp(2,j)
                k(2,1) = k(2,1) + fact * shp(2,i) * shp(1,j)
                k(2,2) = k(2,2) + fact * shp(2,i) * shp(2,j)
                fact   = trgu * rho
                k(1,1) = k(1,1) + fact * shp(1,i) * dtau(1,3)
                k(1,2) = k(1,2) + fact * shp(1,i) * dtau(2,3)
                k(2,1) = k(2,1) + fact * shp(2,i) * dtau(1,3)
                k(2,2) = k(2,2) + fact * shp(2,i) * dtau(2,3)
              else                             ! axisymmetric
                fact   = tau(3) * rho
                k(1,1) = k(1,1) + fact * (shp(1,i)+shp(3,i)/rad)
     *                                 * (shp(1,j)+shp(3,j)/rad)
                k(1,2) = k(1,2) + fact * (shp(1,i)+shp(3,i)/rad)
     *                                 *  shp(2,j)
                k(2,1) = k(2,1) + fact *  shp(2,i)
     *                                 * (shp(1,j)+shp(3,j)/rad)
                k(2,2) = k(2,2) + fact *  shp(2,i) * shp(2,j)
                fact   = rho * (trgu + urdr)
                k(1,1) = k(1,1) + fact*(shp(1,i)+shp(3,i)/rad)*dtau(1,3)
                k(1,2) = k(1,2) + fact*(shp(1,i)+shp(3,i)/rad)*dtau(2,3)
                k(2,1) = k(2,1) + fact * shp(2,i) * dtau(1,3)
                k(2,2) = k(2,2) + fact * shp(2,i) * dtau(2,3)
              endif
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
            if (stabm) then
              fact   = tau(1) * (shp(1,i)*uw(1) + shp(2,i)*uw(2))
     *                    * rho * shp(3,j)
              k(1,1) = k(1,1) + fact
              k(2,2) = k(2,2) + fact
              fact   = tau(2) * rho * shp(3,j)
              k(3,1) = fact * shp(1,i)
              k(3,2) = fact * shp(2,i)
            endif

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
c==========================================================================
 6    continue ! isw 6 -> compute residual + tangent for  x_{n+1}
               !     8 -> compute residual + full tangent for x_{n+1} Lagr
               !     9 -> compute residual + full tangent for x_{n+1} ALE
c-----------------------------------------------------------------------------

c               .     
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
   
      area2  =  (xa(1,1)-xa(1,3))*(xa(2,2)-xa(2,3)) 
     *        - (xa(1,2)-xa(1,3))*(xa(2,1)-xa(2,3))

      if (area2 .lt. 0) then
        write(*,100) isw
        return
      endif

      h      = sqrt(r2/pi * area2)  ! diameter of circle

      dhfact = r1 / sqrt(pi * r2 * area2)

      darea2(1,1) =   (xa(2,2) - xa(2,3))
      darea2(2,1) = - (xa(1,2) - xa(1,3))
      darea2(1,2) = - (xa(2,1) - xa(2,3))
      darea2(2,2) =   (xa(1,1) - xa(1,3))
      darea2(1,3) = - (xa(2,2) - xa(2,1))
      darea2(2,3) = - (xa(1,1) - xa(1,2))

      call pmove(darea2,darea2r,6)

c...  evaluate tau and its derivatives at centroid

      u(1)  = r1d3 * (ua(1,1)+ua(1,2)+ua(1,3)-wa(1,1)-wa(1,2)-wa(1,3))
      u(2)  = r1d3 * (ua(2,1)+ua(2,2)+ua(2,3)-wa(2,1)-wa(2,2)-wa(2,3))
      call eval_tau   (u,h,rho,mu,dt,nRC,z,beta,tau)
      call eval_dtau  (u,h,rho,mu,dt,nRC,z,beta,dtau,r1d3)
      call eval_dtaudh(u,h,rho,mu,dt,nRC,z,beta,dtaudh)
      do ii=1, 2
        fact           = dtaudh(ii) * dhfact
        dtaudx(1,1,ii) = fact * darea2(1,1)
        dtaudx(2,1,ii) = fact * darea2(2,1)
        dtaudx(1,2,ii) = fact * darea2(1,2)
        dtaudx(2,2,ii) = fact * darea2(2,2)
        dtaudx(1,3,ii) = fact * darea2(1,3)
        dtaudx(2,3,ii) = fact * darea2(2,3)
      enddo

c...  get Gauss point coordinates and weights

      l    = -3
      iord = 1
      call tint2d(l,lint,sg)

c...  start loop over Gauss points and get shape functions

      do l=1, lint

        call shp_spaceGP(xa,shp,Dshp,sg(1,l))

        dvol = area2 * r1d6

c...    compute radius r, u_r / r (if axisymmetric)

        if (axsy) then
          rad  = r0
          urdr = r0
          do i=1, 3
            rad  = rad  + shp(3,i) * xa(1,i) 
            urdr = urdr + shp(3,i) * ua(1,i)
          enddo
          rad2  = rad * rad
          urdr  = urdr / rad
          urdrr = urdr / rad
          dvol  = dvol * rad * r2pi
          darea2r(1,1) = (darea2(1,1) * rad + area2 * shp(3,1))*r2pi
          darea2r(2,1) = (darea2(2,1) * rad                   )*r2pi
          darea2r(1,2) = (darea2(1,2) * rad + area2 * shp(3,2))*r2pi
          darea2r(2,2) = (darea2(2,2) * rad                   )*r2pi
          darea2r(1,3) = (darea2(1,3) * rad + area2 * shp(3,3))*r2pi
          darea2r(2,3) = (darea2(2,3) * rad                   )*r2pi
        endif

c...    compute u, p, grad(u), tr(grad(u)), du/dt and their derivatives with 
c       respect to xa 

        call pzero(u,3)
        call pzero(uw,2)                                    
        call pzero(gu,4)                                   
        call pzero(gp,2)
        call pzero(du,2)
        call pzero(dgp,12)
        call pzero(dgu,24)
        call pzero(dtrgu,6)
        do i=1, 3
          do ii=1, 2                                              !
            u(ii)       = u(ii)  +  ua(ii,i)           * shp(3,i) ! u
            uw(ii)      = uw(ii) + (ua(ii,i)-wa(ii,i)) * shp(3,i) ! u - w
            gp(ii)      = gp(ii)  +  ua( 3,i) * shp(ii,i)         ! grad(p)
            du(ii)      = du(ii)  + dua(ii,i) * shp(3,i)          ! du/dt
            do jj=1, 2                                            !
              gu(ii,jj) = gu(ii,jj) + ua(ii,i) * shp(jj,i)        ! grad(u)
            enddo                                                 !
          enddo                                                   !
          u(3)          = u(3)      + ua(3,i)  * shp(3,i)         ! p
                                                                  !
          do jj=1, 2
            do j=1, 3
              do ii=1, 2
                dgp(ii,jj,i) = dgp(ii,jj,i) + ua (3,j)*Dshp(ii,j,jj,i)
                do kk=1, 2
                  dgu(ii,kk,jj,i) = dgu(ii,kk,jj,i) 
     *                               + ua(ii,j) * Dshp(kk,j,jj,i)
                enddo ! kk
              enddo ! ii
            enddo ! j
          enddo ! jj
          dtrgu(1,i) = dgu(1,1,1,i) + dgu(2,2,1,i)
          dtrgu(2,i) = dgu(1,1,2,i) + dgu(2,2,2,i)
        enddo ! i                                                 !
        trgu = gu(1,1) + gu(2,2)                                  ! tr(grad(u))
                                                                  !
c...    compute Gauss point values for stabilization

        if (stabm) then
          rstab(1) = rho*(du(1)+gu(1,1)*uw(1)+gu(1,2)*uw(2)-f(1))+gp(1)
          rstab(2) = rho*(du(2)+gu(2,1)*uw(1)+gu(2,2)*uw(2)-f(2))+gp(2)
        endif

c:::::: RESIDUAL VECTOR ::::::::::::::::::::::::::::::::::::::::::::::::::::::

        do i=1, 3
          i0 = (i-1) * ndf

c         time derivative term and convective acceleration, body load
          fact = shp(3,i) * rho
          r(1) = fact * (du(1) + gu(1,1)*uw(1) + gu(1,2)*uw(2) - f(1))
          r(2) = fact * (du(2) + gu(2,1)*uw(1) + gu(2,2)*uw(2) - f(2))

c         viscous part
          r(1) = r(1) + mu * (shp(1,i) * (gu(1,1) + gu(1,1))
     *                      + shp(2,i) * (gu(1,2) + gu(2,1)))
          r(2) = r(2) + mu * (shp(1,i) * (gu(1,2) + gu(2,1))
     *                      + shp(2,i) * (gu(2,2) + gu(2,2)))

c         pressure terms
          r(1) = r(1) - u(3) * shp(1,i)
          r(2) = r(2) - u(3) * shp(2,i)
          r(3) = trgu * shp(3,i)

c         additional terms for axisymmetric case 
          if (axsy) then
            r(1) = r(1) + mu * shp(3,i) * (urdrr + urdrr) ! viscous term
     *                  - u(3) * shp(3,i) / rad           ! pressure term
            r(3) = r(3) + urdr * shp(3,i)                 ! pressure term
          endif

c         stabilization of momentum equation
          if (stabm) call rupwind(shp,uw,gu,tau,rstab,drdui,rstabi,r,i)

c         stabilization of continuity equation
          if (stabc) then
            if (.not.axsy) then              ! 2D
              fact = tau(3) * rho * trgu
              r(1) = r(1) + fact * shp(1,i)
              r(2) = r(2) + fact * shp(2,i)
            else                             ! axisymmetric
              fact = tau(3) * rho * (trgu + urdr)
              r(1) = r(1) + fact * (shp(1,i) + shp(3,i) / rad)
              r(2) = r(2) + fact * shp(2,i)
            endif
          endif

c         assembly
          do ii=1, 3
            p(i0+ii) = p(i0+ii) - r(ii) * dvol
          enddo

c:::::::: STIFFNESS MATRIX FOR x_{n+1} :::::::::::::::::::::::::::::::::::::::

          do j=1, 3
            j0 = (j-1) * ndm

            drdxj(1,1) = rho * (dgu(1,1,1,j)*uw(1) + dgu(1,2,1,j)*uw(2))
            drdxj(1,2) = rho * (dgu(1,1,2,j)*uw(1) + dgu(1,2,2,j)*uw(2))
            drdxj(2,1) = rho * (dgu(2,1,1,j)*uw(1) + dgu(2,2,1,j)*uw(2))
            drdxj(2,2) = rho * (dgu(2,1,2,j)*uw(1) + dgu(2,2,2,j)*uw(2))

c           time derivative term and convective acceleration
            k(1,1) = shp(3,i) * drdxj(1,1)
            k(1,2) = shp(3,i) * drdxj(1,2)
            k(2,1) = shp(3,i) * drdxj(2,1)
            k(2,2) = shp(3,i) * drdxj(2,2)

c           viscous part  ! could be integrated with only one GP
            k(1,1) = k(1,1) + mu * (Dshp(1,i,1,j) * (gu(1,1) + gu(1,1))
     *                            + Dshp(2,i,1,j) * (gu(1,2) + gu(2,1))
     *                       + shp(1,i) * (dgu(1,1,1,j) + dgu(1,1,1,j))
     *                       + shp(2,i) * (dgu(1,2,1,j) + dgu(2,1,1,j)))
            k(1,2) = k(1,2) + mu * (Dshp(1,i,2,j) * (gu(1,1) + gu(1,1))
     *                            + Dshp(2,i,2,j) * (gu(1,2) + gu(2,1))
     *                       + shp(1,i) * (dgu(1,1,2,j) + dgu(1,1,2,j))
     *                       + shp(2,i) * (dgu(1,2,2,j) + dgu(2,1,2,j)))
            k(2,1) = k(2,1) + mu * (Dshp(1,i,1,j) * (gu(2,1) + gu(1,2))
     *                            + Dshp(2,i,1,j) * (gu(2,2) + gu(2,2))
     *                       + shp(1,i) * (dgu(2,1,1,j) + dgu(1,2,1,j))
     *                       + shp(2,i) * (dgu(2,2,1,j) + dgu(2,2,1,j)))
            k(2,2) = k(2,2) + mu * (Dshp(1,i,2,j) * (gu(2,1) + gu(1,2))
     *                            + Dshp(2,i,2,j) * (gu(2,2) + gu(2,2))
     *                       + shp(1,i) * (dgu(2,1,2,j) + dgu(1,2,2,j))
     *                       + shp(2,i) * (dgu(2,2,2,j) + dgu(2,2,2,j)))

c           pressure terms  ! could be integrated with only one GP  
            k(1,1) = k(1,1) - u(3) * Dshp(1,i,1,j)
            k(1,2) = k(1,2) - u(3) * Dshp(1,i,2,j)
            k(2,1) = k(2,1) - u(3) * Dshp(2,i,1,j)
            k(2,2) = k(2,2) - u(3) * Dshp(2,i,2,j)
            k(3,1) = dtrgu(1,j) * shp(3,i)
            k(3,2) = dtrgu(2,j) * shp(3,i)

c           additional terms for axisymmetric case 
            if (axsy) then
              k(1,1) = k(1,1) - (mu * r4 * urdrr / rad - u(3) / rad2 )
     *                          * shp(3,i) * shp(3,j)
              k(3,1) = k(3,1) - urdrr * shp(3,j) * shp(3,i)  
            endif

c           stabilization of momentum equation
            if (stabm) then
              drdxj(1,1) = drdxj(1,1) + dgp(1,1,j)
              drdxj(1,2) = drdxj(1,2) + dgp(1,2,j)
              drdxj(2,1) = drdxj(2,1) + dgp(2,1,j)
              drdxj(2,2) = drdxj(2,2) + dgp(2,2,j)
              drduidw(1) = Dshp(1,i,1,j)*uw(1) + Dshp(2,i,1,j)*uw(2)
              drduidw(2) = Dshp(1,i,2,j)*uw(1) + Dshp(2,i,2,j)*uw(2)
              call kupwindw(shp,tau,dtaudx,rstab,rstabi,drdui,drduidw,
     *                      drdxj,Dshp,k,i,j)
            endif

c           stabilization of continuity equation          
            if (stabc) then
              if (.not.axsy) then              ! 2D
                fact   = tau(3) * rho * trgu
                k(1,1) = k(1,1) + fact * Dshp(1,i,1,j)
                k(1,2) = k(1,2) + fact * Dshp(1,i,2,j)
                k(2,1) = k(2,1) + fact * Dshp(2,i,1,j)
                k(2,2) = k(2,2) + fact * Dshp(2,i,2,j)
                fact   = tau(3) * rho
                k(1,1) = k(1,1) + fact * shp(1,i) * dtrgu(1,j)
                k(1,2) = k(1,2) + fact * shp(1,i) * dtrgu(2,j)
                k(2,1) = k(2,1) + fact * shp(2,i) * dtrgu(1,j)
                k(2,2) = k(2,2) + fact * shp(2,i) * dtrgu(2,j)
                fact   = dtaudh(3) * dhfact * rho * trgu
                k(1,1) = k(1,1) + fact * shp(1,i) * darea2(1,j)
                k(1,2) = k(1,2) + fact * shp(1,i) * darea2(2,j)
                k(2,1) = k(2,1) + fact * shp(2,i) * darea2(1,j)
                k(2,2) = k(2,2) + fact * shp(2,i) * darea2(2,j)
              else                             ! axisymmetric
                fact   = tau(3) * rho * (trgu + urdr)
                k(1,1) = k(1,1) + fact * 
     *                   (Dshp(1,i,1,j)-shp(3,i)/rad2**shp(3,j))
                k(1,2) = k(1,2) + fact * Dshp(1,i,2,j)
                k(2,1) = k(2,1) + fact * Dshp(2,i,1,j)
                k(2,2) = k(2,2) + fact * Dshp(2,i,2,j)
                fact   = tau(3) * rho
                k(1,1) = k(1,1) + fact * (shp(1,i)+shp(3,i)/rad) * 
     *                   (dtrgu(1,j)-urdrr*shp(3,j))
                k(1,2) = k(1,2) + fact * (shp(1,i)+shp(3,i)/rad) * 
     *                   dtrgu(2,j)
                k(2,1) = k(2,1) + fact * shp(2,i) * 
     *                   (dtrgu(1,j)-urdrr*shp(3,j))
                k(2,2) = k(2,2) + fact * shp(2,i) * dtrgu(2,j)
                fact   = dtaudh(3) * dhfact * rho * (trgu + urdr)
                k(1,1) = k(1,1) + fact * (shp(1,i)+shp(3,i)/rad) * 
     *                   darea2(1,j)
                k(1,2) = k(1,2) + fact * (shp(1,i)+shp(3,i)/rad) * 
     *                   darea2(2,j)
                k(2,1) = k(2,1) + fact * shp(2,i) * darea2(1,j)
                k(2,2) = k(2,2) + fact * shp(2,i) * darea2(2,j)
              endif
            endif

c           assembly
            do ii=1, 3
              do jj=1, 2
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + (k(ii,jj) * dvol
     *                     + r(ii) * darea2r(jj,j) * r1d6) * p6
              enddo
            enddo

          enddo

          if (isw.eq.6) goto 602

c:::::::: STIFFNESS MATRIX FOR dx_{n+1} ::::::::::::::::::::::::::::::::::::::

          do j=1, 3
            j0 = (j-1) * ndm

c           time derivative term and convective acceleration
            fact   = - shp(3,i) * rho * shp(3,j)
            k(1,1) = fact * gu(1,1)
            k(1,2) = fact * gu(1,2)
            k(2,1) = fact * gu(2,1)
            k(2,2) = fact * gu(2,2)

c           stabilization of momentum equation
            if (stabm) then
              fact   = (shp(1,i) * uw(1) + shp(2,i) * uw(2)) * rho
              k(1,1) = k(1,1) - tau(1) * shp(3,j) * 
     *                           (fact * gu(1,1) + rstab(1) * shp(1,i))
              k(1,2) = k(1,2) - tau(1) * shp(3,j) * 
     *                           (fact * gu(1,2) + rstab(1) * shp(2,i)) 
              k(2,1) = k(2,1)  - tau(1) * shp(3,j) * 
     *                           (fact * gu(2,1) + rstab(2) * shp(1,i))
              k(2,2) = k(2,2)  - tau(1) * shp(3,j) * 
     *                           (fact * gu(2,2) + rstab(2) * shp(2,i))
              fact   = - tau(2) * rho * shp(3,j)
              k(3,1) = fact * (shp(1,i) * gu(1,1) + shp(2,i) * gu(2,1))
              k(3,2) = fact * (shp(1,i) * gu(1,2) + shp(2,i) * gu(2,2))
              fact   = dtau(1,1) / tau(1)
              k(1,1) = k(1,1) - fact * rstabi(1)  
              k(2,1) = k(2,1) - fact * rstabi(2)  
              fact   = dtau(2,1) / tau(1)
              k(1,2) = k(1,2) - fact * rstabi(1)  
              k(2,2) = k(2,2) - fact * rstabi(2)
              k(3,1) = k(3,1) - dtau(1,2) / tau(2) * rstabi(3)
              k(3,2) = k(3,2) - dtau(2,2) / tau(2) * rstabi(3)
            endif

c           assembly
            if (isw.eq.8) then
              fact = dvol * p7 * b33
            else
              fact = dvol * p7 * bALE
            endif
            do ii=1, 3
              do jj=1, 2
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * fact
              enddo
            enddo

          enddo

 602      continue

        enddo

      enddo

      return

c==========================================================================
 7    continue ! isw 7 -> compute tangent for dx_{n+1}
c-----------------------------------------------------------------------------

c               .     
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

      if (area2 .lt. 0) then
        write(*,100) isw
        return
      endif

      h     = sqrt(r2/pi * area2)  ! diameter of circle

c...  evaluate tau and its derivatives at centroid

      u(1) = r1d3 * (ua(1,1)+ua(1,2)+ua(1,3)-wa(1,1)-wa(1,2)-wa(1,3))
      u(2) = r1d3 * (ua(2,1)+ua(2,2)+ua(2,3)-wa(2,1)-wa(2,2)-wa(2,3))
      call eval_tau (u,h,rho,mu,dt,nRC,z,beta,tau)
      call eval_dtau(u,h,rho,mu,dt,nRC,z,beta,dtau,r1d3)

c...  get Gauss point coordinates and weights

      l    = -3
      iord = 1
      call tint2d(l,lint,sg)

c...  start loop over Gauss points and get shape functions

      do l=1, lint

        call trishp(sg(1,l),xa,ndm,iord,xsj,shp)

        dvol = xsj * sg(4,l)

c...    compute radius r, u_r / r (if axisymmetric)

        if (axsy) then
          rad  = r0
          urdr = r0
          do i=1, 3
            rad  = rad  + shp(3,i) * xa(1,i) 
            urdr = urdr + shp(3,i) * ua(1,i)
          enddo
          rad2   = rad * rad
          urdr  = urdr / rad
          urdrr = urdr / rad
          dvol  = dvol * rad * r2pi
        endif

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

        if (stabm) then
          rstab(1) = rho*(du(1)+gu(1,1)*uw(1)+gu(1,2)*uw(2)-f(1))+gp(1)
          rstab(2) = rho*(du(2)+gu(2,1)*uw(1)+gu(2,2)*uw(2)-f(2))+gp(2)
        endif

        do i=1, 3
          i0 = (i-1) * ndf

c         stabilization of momentum equation
          if (stabm) call rupwind(shp,uw,gu,tau,rstab,drdui,rstabi,r,i)

c:::::::: STIFFNESS MATRIX FOR dx_{n+1} ::::::::::::::::::::::::::::::::::::::

          k(3,1) = r0
          k(3,2) = r0

          do j=1, 3
            j0 = (j-1) * ndm

c           time derivative term and convective acceleration
            fact   = - shp(3,i) * rho * shp(3,j)
            k(1,1) = fact * gu(1,1)
            k(1,2) = fact * gu(1,2)
            k(2,1) = fact * gu(2,1)
            k(2,2) = fact * gu(2,2)

c           stabilization of momentum equation
            if (stabm) then
              fact   = (shp(1,i) * uw(1) + shp(2,i) * uw(2)) * rho
              k(1,1) = k(1,1) - tau(1) * shp(3,j) * 
     *                           (fact * gu(1,1) + rstab(1) * shp(1,i))
              k(1,2) = k(1,2) - tau(1) * shp(3,j) * 
     *                           (fact * gu(1,2) + rstab(1) * shp(2,i)) 
              k(2,1) = k(2,1)  - tau(1) * shp(3,j) * 
     *                           (fact * gu(2,1) + rstab(2) * shp(1,i))
              k(2,2) = k(2,2)  - tau(1) * shp(3,j) * 
     *                           (fact * gu(2,2) + rstab(2) * shp(2,i))
              fact   = - tau(2) * rho * shp(3,j)
              k(3,1) = fact * (shp(1,i) * gu(1,1) + shp(2,i) * gu(2,1))
              k(3,2) = fact * (shp(1,i) * gu(1,2) + shp(2,i) * gu(2,2))
              fact   = dtau(1,1) / tau(1)
              k(1,1) = k(1,1) - fact * rstabi(1)  
              k(2,1) = k(2,1) - fact * rstabi(2)  
              fact   = dtau(2,1) / tau(1)
              k(1,2) = k(1,2) - fact * rstabi(1)  
              k(2,2) = k(2,2) - fact * rstabi(2)
              k(3,1) = k(3,1) - dtau(1,2) / tau(2) * rstabi(3)
              k(3,2) = k(3,2) - dtau(2,2) / tau(2) * rstabi(3)
            endif

c           assembly
            fact = dvol * p7
            do ii=1, 3
              do jj=1, 2
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * fact
              enddo
            enddo

          enddo

        enddo

      enddo

      return

cc==============================================================================
c 10   continue ! isw 10 -> compute volume, kinetic & potential energy
cc-----------------------------------------------------------------------------
c
cc     OUTPUT:    dp(1) = dp(1) + element volume
cc                dp(2) = dp(2) + element kinetic energy
cc                dp(3) = dp(3) + element potential energy
cc                dp(4) = dp(4) + element kinetic + potential energy
cc                dp(5) = dp(5) + vorticity * area
cc                dp(6) = dp(6) + div(u) * area
c
c      axsy = nint(d(10)).ne.0
c      rho  = d(12)
c
c      if (axsy) rad = r1d3 * (xl(1,1) + xl(1,2) + xl(1,3))
c
c      do i=1, 3
c        j = i + 1 - max(0,i-2)*3
c        l = j + 1 - max(0,j-2)*3
c        fact     = r1 / ( (xl(2,i)-xl(2,j)) * (xl(1,j)-xl(1,l))
c     *                   -(xl(2,j)-xl(2,l)) * (xl(1,i)-xl(1,j)))
c        shp(1,i) = - fact * (xl(2,j) - xl(2,l))
c        shp(2,i) =   fact * (xl(1,j) - xl(1,l))
c      enddo
c
c      area2 =  (xl(1,1)-xl(1,3))*(xl(2,2)-xl(2,3)) 
c     *       - (xl(1,2)-xl(1,3))*(xl(2,1)-xl(2,3))
c
c      if (area2.le.r0) then
c        write(*,*) '       element ',n
c        call error('elmt02fld',45,' ')
c      endif
c
cc     compute element volume (area)
c      fact = r1d2 * area2
c      if (axsy) fact = fact * rad * r2pi
c      dp(1) = dp(1) + fact
c
cc     compute kinetic energy
c      do i=1, 3
c        j = i + 1 - max(0,i-2)*3
c        l = j + 1 - max(0,j-2)*3
c        rad2 = r1
c        if (axsy) rad2 = r1d6 * (xl(1,i)+xl(1,j)+r4*xl(1,l)) / rad
c        hlp   = r1d216 * rho * fact * rad2 * 
c     *       ((ul(1,i)+ul(1,j)+r4*ul(1,l))*(ul(1,i)+ul(1,j)+r4*ul(1,l))
c     *       +(ul(2,i)+ul(2,j)+r4*ul(2,l))*(ul(2,i)+ul(2,j)+r4*ul(2,l)))
c        dp(2) = dp(2) + hlp
c        dp(4) = dp(4) + hlp
c      enddo
c
cc     compute potential energy with respect to y=0/ z=0
c      dp(3) = dp(3) + fact * rho * r1d3 * (xl(2,1)+xl(2,2)+xl(2,3))
c      dp(4) = dp(4) + fact * rho * r1d3 * (xl(2,1)+xl(2,2)+xl(2,3))
c
cc     compute vorticity  
c      do i=1, 3
c        dp(5) = dp(5) + fact * (ul(2,i) * shp(1,i) - ul(1,i) * shp(2,i))
c      enddo
c
cc     div(u) * area
c      do i=1, 3
c        dp(6) = dp(6) + (ul(1,i) * shp(1,i) + ul(2,i) * shp(2,i)) * fact
c      enddo
c      if (axsy)
c     *  dp(6) = dp(6) + ((ul(1,1)+ul(1,2)+ul(1,3))
c     *                  /(xl(1,1)+xl(1,2)+xl(1,3))) * fact
c
c      return
c==============================================================================
 11   continue ! isw 11 -> project to nodes vort(u)
               !     12 -> project to nodes div(u)
               !     13 -> project to nodes ||grad(u)||
c-----------------------------------------------------------------------------

      if (axsy) rad = r1d3 * (xl(1,1) + xl(1,2) + xl(1,3))

      do i=1, 3
        j = i + 1 - max(0,i-2)*3
        l = j + 1 - max(0,j-2)*3
        fact     = r1 / ( (xl(2,i)-xl(2,j)) * (xl(1,j)-xl(1,l))
     *                   -(xl(2,j)-xl(2,l)) * (xl(1,i)-xl(1,j)))
        shp(1,i) = - fact * (xl(2,j) - xl(2,l))
        shp(2,i) =   fact * (xl(1,j) - xl(1,l))
      enddo
      
      fact = r0

      if (isw.eq.11) then

c       vorticity  
        do i=1, 3
          fact = fact + (ul(2,i)*shp(1,i) - ul(1,i)*shp(2,i))
        enddo

      elseif (isw.eq.12) then

c       divergence
        do i=1, 3
          fact = fact + (ul(1,i)*shp(1,i) + ul(2,i)*shp(2,i))
        enddo
        if (axsy)
     *    fact = fact + ((ul(1,1)+ul(1,2)+ul(1,3))
     *                  /(xl(1,1)+xl(1,2)+xl(1,3)))

      else

c       ||grad(u)||
        gu(1,1) = ul(1,1)*shp(1,1) + ul(1,2)*shp(1,2) + ul(1,3)*shp(1,3)
        gu(1,2) = ul(1,1)*shp(2,1) + ul(1,2)*shp(2,2) + ul(1,3)*shp(2,3)
        gu(2,1) = ul(2,1)*shp(1,1) + ul(2,2)*shp(1,2) + ul(2,3)*shp(1,3)
        gu(2,2) = ul(2,1)*shp(2,1) + ul(2,2)*shp(2,2) + ul(2,3)*shp(2,3)
        fact =  sqrt(abs(gu(1,1)*gu(1,1) + gu(1,2)*gu(1,2) 
     *                 + gu(2,1)*gu(2,1) + gu(2,2)*gu(2,2)))
 
      endif

c...  project to nodes

      p(1) = fact
      p(2) = fact
      p(3) = fact

      return

c==============================================================================
 14   continue ! isw 14-15 -> not used
c-----------------------------------------------------------------------------

      write(*,'(/9x,a/)') 'elmt03fld: what do you want ?'

      return

c==============================================================================









cc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc
cc     The following is needed for error indication.
cc
cc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c
c
cc===============================================================================
c 16   continue ! isw 16 -> element aspect ratio
cc-------------------------------------------------------------------------------
c
c      st=r0
c      do i=1,3
c        j=i+1
c        if(j.gt.3) j=j-3
c        dtemp=r0
c        do l=1,ndm
c          dtemp=dtemp+(xl(l,j)-xl(l,i))*(xl(l,j)-xl(l,i))
c        enddo
c        a(i)=sqrt(dtemp)
c        st=st+a(i)
c      enddo
c      st=r1d2*st
cc
cc ... compute AREA2
cc
c      b=r1
c      c=r1
c      do i=1,3
c        b=b*(st-a(i))
c        c=c*a(i)
c      enddo
c      area2=st*b
cc
cc ... compute aspect ratio
cc
c      fact=r1d4*c*st/area2
cc
c
c      dp(1)=max(dp(1),fact)
c
c      return
c
cc=============================================================================
c 17   continue ! isw 17 -> project to nodes, element size
cc-----------------------------------------------------------------------------
c
c      area2 =  (xl(1,1)-xl(1,3))*(xl(2,2)-xl(2,3))
c     *       - (xl(1,2)-xl(1,3))*(xl(2,1)-xl(2,3))
c
c      fact  = sqrt(area2) * 1.07457d0
c
c      p(0+1) = fact
c      p(3+1) = fact
c      p(6+1) = fact
c
c      return
c
cc=============================================================================
c 18   continue ! isw 18 -> project to nodes, flow complexity
cc
cc                          flow complexity =  |vort(u)| + ||grad(u)||
cc
cc-----------------------------------------------------------------------------
c
c      fact = r1 / ((xl(1,1)-xl(1,3))*(xl(2,2)-xl(2,3))
c     *           - (xl(1,2)-xl(1,3))*(xl(2,1)-xl(2,3)))
c
c      shp(1,1) = fact * (xl(2,2)-xl(2,3))
c      shp(2,1) = fact * (xl(1,3)-xl(1,2))
c
c      shp(1,2) = fact * (xl(2,3)-xl(2,1))
c      shp(2,2) = fact * (xl(1,1)-xl(1,3))
c
c      shp(1,3) = fact * (xl(2,1)-xl(2,2))
c      shp(2,3) = fact * (xl(1,2)-xl(1,1))
c
ccc     vorticity
c      fact = 0.d0
cc      do i=1, 3
cc        fact = fact + (ul(2,i)*shp(1,i) - ul(1,i)*shp(2,i))
cc      enddo
cc      fact = abs(fact)
c
cc     velocity gradient 
c      gu(1,1) = ul(1,1)*shp(1,1) + ul(1,2)*shp(1,2) + ul(1,3)*shp(1,3)
c      gu(1,2) = ul(1,1)*shp(2,1) + ul(1,2)*shp(2,2) + ul(1,3)*shp(2,3)
c      gu(2,1) = ul(2,1)*shp(1,1) + ul(2,2)*shp(1,2) + ul(2,3)*shp(1,3)
c      gu(2,2) = ul(2,1)*shp(2,1) + ul(2,2)*shp(2,2) + ul(2,3)*shp(2,3)
c      fact    =  fact + sqrt(abs(gu(1,1)*gu(1,1) + gu(1,2)*gu(1,2) 
c     *                         + gu(2,1)*gu(2,1) + gu(2,2)*gu(2,2)))
c
c      p(0+1) = fact
c      p(3+1) = fact
c      p(6+1) = fact
c
c      return
c
cc=============================================================================
c 19   continue ! isw 19,20, ...  -> not used yet
cc-----------------------------------------------------------------------------
c
c      write(*,'(/9x,a/)') 'elmt03fld: what do you want ?'
c
c      return
c
cc=============================================================================

100   format('   element2DStabIncompFluid (isw =',I2,
     *       '): A < 0 => ignore this element!')

      end
















