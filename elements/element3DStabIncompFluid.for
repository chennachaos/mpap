
      integer function element3DStabIncompFluid
     *                    (elmDat,timDat,xl,ul,dt,s,p,
     *                     ndf,ndm,nen,nelm,isw)

c  ***************************************************************************
c  *                                                                         *
c  *              3D STABILIZED VELOCITY-PRESSURE FLUID ELEMENT              *
c  *                                                                         *
c  *                       TRI-LINEAR TETRAHEDRA                             *
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
c  *    d(10)         not used                                               *
c  *                                                                         *
c  *    d(11)         mu                                                     *
c  *    d(12)         rho                                                    *
c  *    d(13)         f(1)   body force vector                               *
c  *    d(14)         f(2)                                                   *
c  *    d(15)         f(3)                                               * 
c  *                                                                         * 
c  *  ---------------------------------------------------------------------  * 
c  *        ul(j,i),  j = 1, 2, 3, 4=ndf,  i = 1, 2, 3, 4=nen                *
c  *                                                                         * 
c  *    j = 1. dof  -  x velocity                                            *
c  *        2. dof  -  y velocity                                            *
c  *        3. dof  -  z velocity                                            *
c  *        4. dof  -  pressure                                              *
c  *                                                                         * 
c  *    i = 1...4   -  current nodal values    u^k(t)                        * 
c  *        5...8   -  increment  u^k(t) - u(t_n)  .                         *
c  *        9...12  -  previous solution rate      u(t_n)   (if required)    *
c  *                                                                         *
c  *  ---------------------------------------------------------------------  * 
c  *        xl(j,i),  i = 1, 2, 3=ndm,     j = 1, 2, 3, 4=nen                *
c  *                                                                         * 
c  *    j = 1. dof  -  x                                                     *
c  *        2. dof  -  y                                                     *
c  *        3. dof  -  z                                                     *
c  *                                                                         *
c  *    i = 1..4    -  nodal coordinates  x_n                                *
c  *        5..8    -  nodal coordinates  x_{n+1} - x_n                      *
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
c  *   10 x not used                                                         *
c  *   11 xproject vort(u) to nodes                                          *
c  *   12 xproject vort(v) to nodes                                          *
c  *   13 xproject vort(w) to nodes                                          *
c  ***************************************************************************

      implicit none

      integer          ndf, ndm, nen, nelm, isw, itau,
     *                 i, j, i0, j0, ii, jj, l, kk, stab,
     *                 ngp, nen2,nen3,err   
     
      double precision elmDat(*), timDat(*), dt, ul(ndf,*), xl(ndm,*), 
     *                 s(nen*ndf,*), p(*),
     *                 mu, rho, h, nRC(3), z(3), beta(6),
     *                 xa(3,4),wa(3,4),ua(4,4),dua(4,4),
     *                 u(4), uw(3), gu(3,3), f(3),  
     *                 gp(3), du(3),wgp(4),xi(3,4),
     *                 shp(3+1,4), r(4), 
     *                 trgu, rstab(3), rstabi(4), 
     *                 drdui(3,3), drduj(3,3), tau(3), 
     *                 dtau(3,3), dtaudx(3,4,3), dtaudh(3),  
     *                 k(4,4), dgp(3,3,4), 
     *                 dvol, dgu(3,3,3,4), dtrgu(3,4), 
     *                 dvolume(3,4),  Dshp(3,4,3,4),
     *                 drdxj(3,3), drduidw(3), fact, dhfact,
     *                 ddV(3,4,3,4),
     *                 r0, r1, r2, r4, r1d3, r1d2, r1d4, r1d6, r1d216,
     *                 pi, r2pi, p6, p7, p8, b17, b33, bALE,
     *                 r3, detJ, fact1, fact2

                    

      logical stabm, stabc
      
      data r0, r1d4, r1d2, r1, r2, r3, r4
     *         / 0.d0, 0.25d0, 0.5d0, 1.d0, 2.d0, 3.d0, 4.d0  /

      r1d3  = r1 /  r3
      r1d4  = r1 /  r4
      r1d6  = r1d2 * r1d3
      r1d216= r1d6 * r1d6 * r1d6
      pi    = acos(-1.d0)
      r2pi  = r2 * pi
      
      nen2=nen+nen
      nen3=nen2+nen
      
c      okflg = .true.

      p6    = timDat(6)
      p7    = timDat(7)
      p8    = timDat(8)

      b17   = timDat(9)
      b33   = timDat(21)  ! ??????????????
      bALE  = timDat(21)
      
      


c...  set control paramters and element constants

      stab  = nint(elmDat(1))
      do i=1, 6
        beta(i) = elmDat(2+i)
      enddo
      mu    = elmDat(11)
      rho   = elmDat(12)
      do i=1, ndm
      f(i)  = elmDat(12+i)
        enddo
      stabm = stab.ne.0.and.stab.ne.2
      stabc = stab.ge.2

      element3DStabIncompFluid = 0

c...  compute u, u, x, w ...........................

      do i=1, nen
        do ii=1, ndm
          ua (ii,i) = p7 * ul(ii,i     ) + (r1-p7) * ul(ii,i+nen ) 
          dua(ii,i) = p8 * ul(ii,i+nen2) + (r1-p8) * ul(ii,i+nen3) 
          xa (ii,i) = p6 * xl(ii,i     ) + (r1-p6) * xl(ii,i+nen ) 
          wa (ii,i) = p7 * xl(ii,i+nen2) + (r1-p7) * xl(ii,i+nen3) 
        enddo
        ua(ndf,i) = ul(ndf,i) ! pressure
      enddo
      
c            1  2 3 4 5 6 7 8 9 10 11 12 13 14 15 
      go to (14,2,2,2,2,6,7,6,6,14,11,11,11,14,14),isw

      call prgerror(1,'element3DStabIncompFluid','what ?')
      
c=============================================================================
 2    continue ! isw 2 -> compute residual
               !     3 -> compute residual + tangent for  u_{n+1}
               !     4 -> compute residual + tangent for du_{n+1}
               !     5 -> compute residual + full tangent for u_{n+1} 
c-----------------------------------------------------------------------------

c...  characteristic element size
c...  get shape functions at centroid   

      if(nen.eq.8) then
      call compxigp3D(xi,wgp,1)
      ngp=8
      else if(nen.eq.4) then
      xi(1,1) = r1d4
      xi(2,1) = r1d4
      xi(3,1) = r1d4
      wgp(1)  = r1d6
      ngp=4
      else
      call prgerror(3,'element3DStabIncompFluid','only linear elements')
      endif
      
      call compshp3D(shp,detJ,xi(1,1),xa,nen)

      h=r2*(r3*detJ*wgp(1)/(r4*pi))**r1d3 ! diameter of sphere
     
c...  evaluate tau and its derivatives at centroid

       do ii=1, ndm  
       u(ii)  = r0
       do i=1, nen                                        
        u(ii)   =  u(ii)+ (ua(ii,i)-wa(ii,i)) * shp(ndm+1,i)! u - w
        enddo
        enddo
     
      call eval_tau3d (u,h,rho,mu,dt,nRC,z,beta,tau,ndm)
      call eval_dtau3d(u,h,rho,mu,dt,nRC,z,beta,dtau,r1d4,ndm)

      
c...  Compute Gauss points coordinates and weights

      call compxigp3D(xi,wgp,ngp)

c...  loop over Gauss points
      do l=1, ngp

c...    compute shape functions  

      call compshp3D(shp,detJ,xi(1,l),xa,nen)
        
c...    compute volume for current Gauss point 

      dvol = wgp(l)  * detJ

c...    compute u, p, grad(u), tr(grad(u)), du/dt

        call pzero(u,ndf)                                    
        call pzero(uw,ndm)                                    
        call pzero(gu,ndm*ndm)                                   
        call pzero(gp,ndm)
        call pzero(du,ndm)
              
      do i=1, nen
       do ii=1, ndm                                     
        u(ii)       = u(ii)+ ua(ii,i)          * shp(ndm+1,i) ! u
        uw(ii)      = uw(ii)+(ua(ii,i)-wa(ii,i))* shp(ndm+1,i) ! u - w
        gp(ii)      = gp(ii)+ ua(ndf,i)         * shp(ii   ,i) ! grad(p)
        du(ii)      = du(ii)+dua(ii,i)          * shp(ndm+1,i) ! du/dt
         do jj=1, ndm                                            
          gu(ii,jj) = gu(ii,jj)+ ua(ii,i)      * shp(jj   ,i) ! grad(u)
         enddo                                                 
       enddo                                                 
       u(ndf)       =  u(ndf)+ua(ndf,i)         * shp(ndm+1,i) ! p
      enddo  
      
      trgu=0
      do ii=1, ndm                                            
          trgu = trgu+gu(ii,ii)                            !tr(grad(u))
      enddo 
                                                 

                                                                 
c... compute Gauss point values for stabilization

      if (stabm) then
       do ii=1, ndm
       rstab(ii) = rho*(du(ii)-f(ii))+gp(ii)
        do jj=1, ndm
        rstab(ii) = rstab(ii)+rho*gu(ii,jj)*uw(jj)
        enddo
       enddo
      endif
      

c:::::: RESIDUAL VECTOR ::::::::::::::::::::::::::::::::::::::::::::::::::::::

        do i=1, nen
          i0 = (i-1) * ndf

c         time derivative term and convective acceleration, body load
      fact = shp(ndm+1,i) * rho
      do ii=1, ndm
      r(ii) = fact *(du(ii)-f(ii))
       do jj=1, ndm
       r(ii) = r(ii)+fact *gu(ii,jj)*uw(jj)
       enddo
      enddo
      
c         viscous part
          r(1) = r(1) + mu * (shp(1,i) * (gu(1,1) + gu(1,1))
     *                      + shp(2,i) * (gu(1,2) + gu(2,1))
     *                      + shp(3,i) * (gu(1,3) + gu(3,1)))
          r(2) = r(2) + mu * (shp(1,i) * (gu(2,1) + gu(1,2))
     *                      + shp(2,i) * (gu(2,2) + gu(2,2))
     *                      + shp(3,i) * (gu(2,3) + gu(3,2)))
          r(3) = r(3) + mu * (shp(1,i) * (gu(3,1) + gu(1,3))
     *                      + shp(2,i) * (gu(3,2) + gu(2,3))
     *                      + shp(3,i) * (gu(3,3) + gu(3,3)))

c         pressure terms
      do ii=1, ndm
      r(ii) = r(ii) - u(ndf) * shp(ii,i)
      enddo
      r(ndf) = trgu * shp(ndm+1,i)


c         stabilization of momentum equation
          if (stabm) 
     *      call rupwind3d(shp,uw,gu,tau,rstab,drdui,rstabi,r,i,ndm,ndf)

c         stabilization of continuity equation
          if (stabc) then
              fact = tau(3) * rho * trgu
              r(1) = r(1) + fact * shp(1,i)
              r(2) = r(2) + fact * shp(2,i)
              r(3) = r(3) + fact * shp(3,i)
          endif

c         assembly
          do ii=1, ndf
            p(i0+ii) = p(i0+ii) - r(ii) * dvol
          enddo

c:::::::: STIFFNESS MATRIX FOR u_{n+1} ::::::::::::::::::::::::::::::::::::::

          if (isw.ne.3.and.isw.ne.5) goto 201

          do j=1, nen
            j0 = (j-1) * ndf

            fact = shp(1,j)*uw(1) + shp(2,j)*uw(2) + shp(3,j)*uw(3)
            drduj(1,1) = rho * (gu(1,1) * shp(ndm+1,j) + fact)
            drduj(1,2) = rho *  gu(1,2) * shp(ndm+1,j)
            drduj(1,3) = rho *  gu(1,3) * shp(ndm+1,j)
            drduj(2,1) = rho *  gu(2,1) * shp(ndm+1,j)
            drduj(2,2) = rho * (gu(2,2) * shp(ndm+1,j) + fact)
            drduj(2,3) = rho *  gu(2,3) * shp(ndm+1,j)
            drduj(3,1) = rho *  gu(3,1) * shp(ndm+1,j)
            drduj(3,2) = rho *  gu(3,2) * shp(ndm+1,j) 
            drduj(3,3) = rho * (gu(3,3) * shp(ndm+1,j) + fact)
            
c           stiffness terms that require single pt integration      
c           viscous part   
            fact = shp(1,i)*shp(1,j)+shp(2,i)*shp(2,j)+shp(3,i)*shp(3,j)
            k(1,1) = mu * (shp(1,j)*shp(1,i) + fact)
            k(1,2) = mu *  shp(1,j)*shp(2,i)
            k(1,3) = mu *  shp(1,j)*shp(3,i)
            k(2,1) = mu *  shp(2,j)*shp(1,i)
            k(2,2) = mu * (shp(2,j)*shp(2,i) + fact)
            k(2,3) = mu *  shp(2,j)*shp(3,i)
            k(3,1) = mu *  shp(3,j)*shp(1,i)
            k(3,2) = mu *  shp(3,j)*shp(2,i) 
            k(3,3) = mu * (shp(3,j)*shp(3,i) + fact)

c           pressure terms    
            k(1,4) = - shp(ndm+1,j) * shp(1,i)
            k(2,4) = - shp(ndm+1,j) * shp(2,i)
            k(3,4) = - shp(ndm+1,j) * shp(3,i)
            k(4,1) = shp(ndm+1,i) * shp(1,j)
            k(4,2) = shp(ndm+1,i) * shp(2,j)
            k(4,3) = shp(ndm+1,i) * shp(3,j)
            k(4,4) = r0

c           time derivative term and convective acceleration
            k(1,1) = k(1,1) + shp(ndm+1,i) * drduj(1,1)
            k(1,2) = k(1,2) + shp(ndm+1,i) * drduj(1,2)
            k(1,3) = k(1,3) + shp(ndm+1,i) * drduj(1,3)
            k(2,1) = k(2,1) + shp(ndm+1,i) * drduj(2,1)
            k(2,2) = k(2,2) + shp(ndm+1,i) * drduj(2,2)
            k(2,3) = k(2,3) + shp(ndm+1,i) * drduj(2,3)
            k(3,1) = k(3,1) + shp(ndm+1,i) * drduj(3,1)
            k(3,2) = k(3,2) + shp(ndm+1,i) * drduj(3,2)
            k(3,3) = k(3,3) + shp(ndm+1,i) * drduj(3,3)

c           stabilization of momentum equation
            if (stabm)  call kupwind3d(shp,tau,dtau,rstab,rstabi
     *                         ,drdui,drduj,k,i,j, ndm,ndf)

c           stabilization of continuity equation 
                   
            if (stabc) then
                fact   = tau(3) * rho
                k(1,1) = k(1,1) + fact * shp(1,i) * shp(1,j)
                k(1,2) = k(1,2) + fact * shp(1,i) * shp(2,j)
                k(1,3) = k(1,3) + fact * shp(1,i) * shp(3,j)
                k(2,1) = k(2,1) + fact * shp(2,i) * shp(1,j)
                k(2,2) = k(2,2) + fact * shp(2,i) * shp(2,j)
                k(2,3) = k(2,3) + fact * shp(2,i) * shp(3,j)
                k(3,1) = k(3,1) + fact * shp(3,i) * shp(1,j)
                k(3,2) = k(3,2) + fact * shp(3,i) * shp(2,j)
                k(3,3) = k(3,3) + fact * shp(3,i) * shp(3,j)
                fact   = trgu * rho
                k(1,1) = k(1,1) + fact * shp(1,i) * dtau(1,3)
                k(1,2) = k(1,2) + fact * shp(1,i) * dtau(2,3)
                k(1,3) = k(1,3) + fact * shp(1,i) * dtau(3,3)
                k(2,1) = k(2,1) + fact * shp(2,i) * dtau(1,3)
                k(2,2) = k(2,2) + fact * shp(2,i) * dtau(2,3)
                k(2,3) = k(2,3) + fact * shp(2,i) * dtau(3,3)
                k(3,1) = k(3,1) + fact * shp(3,i) * dtau(1,3)
                k(3,2) = k(3,2) + fact * shp(3,i) * dtau(2,3)
                k(3,3) = k(3,3) + fact * shp(3,i) * dtau(3,3)
            endif

c           assembly
            fact  = dvol * p7
            do ii=1, ndf
              do jj=1, ndf-1
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * fact
              enddo
              s(i0+ii,j0+ ndf)   = s(i0+ii,j0+ ndf) + k(ii, ndf) * dvol
            enddo

          enddo
          

c:::::::: STIFFNESS MATRIX FOR du_{n+1} :::::::::::::::::::::::::::::::::::::
!
 201      if (isw.ne.4.and.isw.ne.5) goto 202

          k(1,2) = r0
          k(1,3) = r0
          k(2,1) = r0
          k(2,3) = r0
          k(3,1) = r0
          k(3,2) = r0

          do j=1, nen
            j0 = (j-1) * ndf

c           time derivative term and convective acceleration
            fact = shp(ndm+1,i) * rho * shp(ndm+1,j)
            k(1,1) = fact
            k(2,2) = fact
            k(3,3) = fact

c           stabilization of momentum equation
            if (stabm) then
              fact   = tau(1) * (shp(1,i)*uw(1) 
     *         + shp(2,i)*uw(2) + shp(3,i)*uw(3))
     *                    * rho * shp(ndm+1,j)
              k(1,1) = k(1,1) + fact
              k(2,2) = k(2,2) + fact
              k(3,3) = k(3,3) + fact
              fact   = tau(2) * rho * shp(ndm+1,j)
              k(ndf,1) = fact * shp(1,i)
              k(ndf,2) = fact * shp(2,i)
              k(ndf,3) = fact * shp(3,i)
            endif

c           assembly
            fact  = dvol * p8
            if (isw.eq.5) fact = fact * b17
            do ii=1, ndf
              do jj=1, ndf-1
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * fact
              enddo
            enddo

          enddo

 202      continue

        enddo !nen

      enddo ! gauss

      return
c==========================================================================
 6    continue ! isw 6 -> compute residual + tangent for  x_{n+1}
               !     8 -> compute residual + full tangent for x_{n+1} Lagr
               !     9 -> compute residual + full tangent for x_{n+1} ALE
c-----------------------------------------------------------------------------

c...  characteristic element size

      if(nen.ne.4) then
      call prgerror(3,'element3DStabIncompFluid','only tetra4 for ale')
      endif
      
      xi(1,1) = r1d4
      xi(2,1) = r1d4
      xi(3,1) = r1d4

      call shp_spaceGP3d(xa,shp,Dshp,xi(1,1))
      
      call comp_tV(xa,dvol,dvolume,ddV,err)

            
      h=r2*(r3*dvol/(r4*pi))**r1d3 ! diameter of sphere with same volume
      dhfact = (r1/r2/pi)*(h*r1d2)**(-r2)
      

c...  evaluate tau and its derivatives at centroid

      
      
       do ii=1, ndm  
       u(ii)  = r0
       do i=1, nen                                        
        u(ii)   =  u(ii)+ (ua(ii,i)-wa(ii,i)) * shp(ndm+1,i)! u - w
        enddo
        enddo
     
      call eval_tau3d (u,h,rho,mu,dt,nRC,z,beta,tau,ndm)
      call eval_dtau3d(u,h,rho,mu,dt,nRC,z,beta,dtau,r1d4,ndm)
      call eval_dtaudh(u,h,rho,mu,dt,nRC,z,beta,dtaudh)
      
      do itau=1, 2 !3 is not used
        fact = dtaudh(itau) * dhfact
        do i=1, ndm
          do j=1, nen 
            dtaudx(i,j,itau) = fact * dvolume(i,j)
          enddo
        enddo
      enddo
     
      dvol = dvol * r1d4

      do i=1, nen
        do ii=1, ndm
          dvolume(ii,i) = dvolume(ii,i) * r1d4
        enddo
      enddo

c...  Compute Gauss points coordinates and weights
      ngp = 4
      call compxigp3D(xi,wgp,ngp)

c...  loop over Gauss points
      do l=1, ngp

c...    compute shape functions  

      call shp_spaceGP3d(xa,shp,Dshp,xi(1,l))

        
c       respect to xa 

        call pzero(u,ndf)                                    
        call pzero(uw,ndm)                                    
        call pzero(gu,ndm*ndm)                                   
        call pzero(gp,ndm)
        call pzero(du,ndm)
        call pzero(dgp,ndm*ndm*nen)
        call pzero(dgu,ndm*ndm*ndm*nen)
        call pzero(dtrgu,ndm*nen)
              
      do i=1, nen
       do ii=1, ndm                                     
        u(ii)       = u(ii)+ ua(ii,i)           * shp(ndm+1,i) ! u
        uw(ii)      = uw(ii)+(ua(ii,i)-wa(ii,i))* shp(ndm+1,i) ! u - w
        gp(ii)      = gp(ii)+ ua(ndf,i)         * shp(ii   ,i) ! grad(p)
        du(ii)      = du(ii)+dua(ii,i)          * shp(ndm+1,i) ! du/dt
         do jj=1, ndm                                            
          gu(ii,jj) = gu(ii,jj)+ ua(ii,i)       * shp(jj   ,i) ! grad(u)
         enddo                                                 
       enddo                                                 
       u(ndf)       =  u(ndf)+ua(ndf,i)         * shp(ndm+1,i) ! p
      
        do jj=1, ndm
            do j=1, nen
              do ii=1, ndm
              dgp(ii,jj,i) = dgp(ii,jj,i) + ua (ndm+1,j)*Dshp(ii,j,jj,i)
                do kk=1, ndm
                  dgu(ii,kk,jj,i) = dgu(ii,kk,jj,i) 
     *                               + ua(ii,j) * Dshp(kk,j,jj,i)
                enddo ! kk
              enddo ! ii
            enddo ! j
            
            do kk=1, ndm
            dtrgu(jj,i) = dtrgu(jj,i)+dgu(kk,kk,jj,i)
            enddo
          enddo ! jj
          

        enddo ! i 
      
      trgu=0
      do ii=1, ndm                                            
          trgu = trgu+gu(ii,ii)                            !tr(grad(u))
      enddo 
      
 
                                                 

                                                                 
c... compute Gauss point values for stabilization

      if (stabm) then
       do ii=1, ndm
       rstab(ii) = rho*(du(ii)-f(ii))+gp(ii)
        do jj=1, ndm
        rstab(ii) = rstab(ii)+rho*gu(ii,jj)*uw(jj)
        enddo
       enddo
      endif
      

c:::::: RESIDUAL VECTOR ::::::::::::::::::::::::::::::::::::::::::::::::::::::

        do i=1, nen
          i0 = (i-1) * ndf

c         time derivative term and convective acceleration, body load
      fact = shp(ndm+1,i) * rho
      do ii=1, ndm
      r(ii) = fact *(du(ii)-f(ii))
       do jj=1, ndm
       r(ii) = r(ii)+fact *gu(ii,jj)*uw(jj)
       enddo
      enddo
      
c         viscous part
      do ii=1, ndm
       do jj=1, ndm
       r(ii) = r(ii)+mu * shp(jj,i) * (gu(ii,jj) + gu(jj,ii))
       enddo
      enddo


c         pressure terms
      do ii=1, ndm
      r(ii) = r(ii) - u(ndf) * shp(ii,i)
      enddo
      r(ndf) = trgu * shp(ndm+1,i)


c         stabilization of momentum equation
          if (stabm) 
     *      call rupwind3d(shp,uw,gu,tau,rstab,drdui,rstabi,r,i,ndm)

c         stabilization of continuity equation
          if (stabc) then
              fact = tau(3) * rho * trgu
               do ii=1, ndm
               r(ii) = r(ii) + fact * shp(ii,i)
               enddo

          endif

c         assembly
          do ii=1, ndf
            p(i0+ii) = p(i0+ii) - r(ii) * dvol
          enddo

c:::::::: STIFFNESS MATRIX FOR x_{n+1} ::::::::::::::::::::::::::::::::::::::


          do j=1, nen
            j0 = (j-1) * ndm
            
          call pzero(drdxj,ndm*ndm)
          do ii=1, ndm
            do jj=1, ndm
              do kk=1, ndm
              drdxj(ii,jj) = drdxj(ii,jj) + rho * dgu(ii,kk,jj,j)*uw(kk)
              enddo
             enddo
           enddo

c           time derivative term and convective acceleration
          do ii=1, ndm
            do jj=1, ndm
                k(ii,jj) = shp(ndm+1,i) * drdxj(ii,jj)
             enddo
           enddo


c           stiffness terms that require single pt integration      
c           viscous part  
          do ii=1, ndm
            do jj=1, ndm
            fact1 = r0
            fact2 = r0
              do kk=1, ndm
          fact1 = fact1 + Dshp(kk,i,jj,j) * (gu(ii,kk) + gu(kk,ii))
          fact2 = fact2 + shp(kk,i) *(dgu(ii,kk,jj,j) + dgu(kk,ii,jj,j))
              enddo
              k(ii,jj) = k(ii,jj) + mu * (fact1 + fact2)
             enddo
           enddo
           
c           pressure terms   
          do ii=1, ndm
            do jj=1, ndm
                k(ii,jj) =  k(ii,jj)- u(ndm+1) * Dshp(ii,i,jj,j)
             enddo
             k(ndf,ii) = dtrgu(ii,j) * shp(ndm+1,i)
           enddo 
             k(ndf,ndf) = r0 


c           stabilization of momentum equation
            if (stabm) then
             do ii=1, ndm
            do jj=1, ndm
                drdxj(ii,jj) = drdxj(ii,jj) + dgp(ii,jj,j)
             enddo
           enddo
           
            call pzero(drduidw,ndm)
            do ii=1, ndm
             do jj=1, ndm
             drduidw(ii) = drduidw(ii)+Dshp(jj,i,ii,j)*uw(jj)
             enddo
           enddo
              

              call kupwindw3d(shp,tau,dtaudx,rstab,rstabi,drdui,drduidw,
     *                      drdxj,Dshp,k,i,j,ndm,nen)
              endif
              
              
c           stabilization of continuity equation          
            if (stabc) then
            
            fact   = tau(3) * rho * trgu
            do ii=1, ndm
            do jj=1, ndm
                k(ii,jj) =  k(ii,jj) + fact * Dshp(ii,i,jj,j)
             enddo
           enddo 

            fact   = tau(3) * rho
            do ii=1, ndm
            do jj=1, ndm
                k(ii,jj) =  k(ii,jj) + fact * shp(ii,i) * dtrgu(jj,j)
             enddo
           enddo 

 
            fact   = dtaudh(3) * dhfact * rho * trgu
            do ii=1, ndm
            do jj=1, ndm
                k(ii,jj) =  k(ii,jj) + fact * shp(ii,i) * dvolume(jj,j)
             enddo
           enddo 
            endif

c           assembly
            do ii=1, ndf
              do jj=1, ndm
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + (k(ii,jj) * dvol
     *                     + r(ii) * dvolume(jj,j)) * p6

              enddo
            enddo

          enddo
          

c:::::::: STIFFNESS MATRIX FOR dx_{n+1} :::::::::::::::::::::::::::::::::::::
!
         if (isw.eq.6) goto 602

          do j=1, nen
            j0 = (j-1) * ndm

c           time derivative term and convective acceleration
            fact = - shp(ndm+1,i) * rho * shp(ndm+1,j)
            
            do ii=1, ndm
            do jj=1, ndm
                k(ii,jj) =  fact * gu(ii,jj)
             enddo
           enddo 
           


c           stabilization of momentum equation
            if (stabm) then
            
            fact   = r0
              do jj=1, ndm
                fact   = fact + shp(jj,i) * uw(jj) * rho
             enddo
     
            do ii=1, ndm
            do jj=1, ndm
              k(ii,jj) = k(ii,jj) - tau(1) * shp(ndm+1,j) * 
     *                 (fact * gu(ii,jj) + rstab(ii) * shp(jj,i))
             enddo
           enddo 
           
           
              fact   = - tau(2) * rho * shp(ndm+1,j)
              do ii=1, ndm
              k(ndm+1,ii) = r0
              do jj=1, ndm
              k(ndm+1,ii) = k(ndm+1,ii)+ 
     *                      fact * shp(jj,i) * gu(jj,ii)
             enddo
             enddo 
             
              

              do ii=1, ndm
              fact   = dtau(ii,1) / tau(1)
              do jj=1, ndm
              k(jj,ii) = k(jj,ii)- 
     *                      fact * rstabi(jj)
             enddo
             enddo 


              do ii=1, ndm
              k(ndm+1,ii) = k(ndm+1,ii)- 
     *                      dtau(ii,2) / tau(2) * rstabi(ndm+1)
             enddo

 
            endif

c           assembly
            if (isw.eq.8) then
              fact = dvol * p7 * b33
            else
              fact = dvol * p7 * bALE
            endif
            
           
            do ii=1, ndf
              do jj=1, ndf-1
                s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * fact
              enddo
            enddo

          enddo

 602      continue

        enddo !nen

      enddo ! gauss

      return

c==========================================================================
 7    continue ! isw 7 -> compute tangent for dx_{n+1}
c-----------------------------------------------------------------------------
     
c...  characteristic element size
c...  get shape functions at centroid   

      if(nen.eq.8) then
      call compxigp3D(xi,wgp,1)
      ngp=8
      else if(nen.eq.4) then
      xi(1,1) = r1d4
      xi(2,1) = r1d4
      xi(3,1) = r1d4
      wgp(1)  = r1d6
      ngp=4
      else
      call prgerror(3,'element3DStabIncompFluid','only linear elements')
      endif
      
      call compshp3D(shp,detJ,xi(1,1),xa,nen)

      h=r2*(r3*detJ*wgp(1)/(r4*pi))**r1d3 ! diameter of sphere
     
c...  evaluate tau and its derivatives at centroid

       do ii=1, ndm  
       u(ii)  = r0
       do i=1, nen                                        
        u(ii)   =  u(ii)+ (ua(ii,i)-wa(ii,i)) * shp(ndm+1,i)! u - w
        enddo
        enddo
     
      call eval_tau3d (u,h,rho,mu,dt,nRC,z,beta,tau,ndm)
      call eval_dtau3d(u,h,rho,mu,dt,nRC,z,beta,dtau,r1d4,ndm)

      
c...  Compute Gauss points coordinates and weights

      call compxigp3D(xi,wgp,ngp)

c...  loop over Gauss points
      do l=1, ngp

c...    compute shape functions  

      call compshp3D(shp,detJ,xi(1,l),xa,nen)
        
c...    compute volume for current Gauss point 

      dvol = wgp(l)  * detJ

c...    compute u, p, grad(u), tr(grad(u)), du/dt

        call pzero(u,ndf)                                    
        call pzero(uw,ndm)                                    
        call pzero(gu,ndm*ndm)                                   
        call pzero(gp,ndm)
        call pzero(du,ndm)
              
      do i=1, nen
       do ii=1, ndm                                     
        u(ii)       = u(ii)+ ua(ii,i)          * shp(ndm+1,i) ! u
        uw(ii)      = uw(ii)+(ua(ii,i)-wa(ii,i))* shp(ndm+1,i) ! u - w
        gp(ii)      = gp(ii)+ ua(ndf,i)         * shp(ii   ,i) ! grad(p)
        du(ii)      = du(ii)+dua(ii,i)          * shp(ndm+1,i) ! du/dt
         do jj=1, ndm                                            
          gu(ii,jj) = gu(ii,jj)+ ua(ii,i)      * shp(jj   ,i) ! grad(u)
         enddo                                                 
       enddo                                                 
       u(ndf)       =  u(ndf)+ua(ndf,i)         * shp(ndm+1,i) ! p
      enddo  
      
      trgu=0
      do ii=1, ndm                                            
          trgu = trgu+gu(ii,ii)                            !tr(grad(u))
      enddo                                  ! tr(grad(u))
                                                                  !
c... compute Gauss point values for stabilization

      if (stabm) then
       do ii=1, ndm
       rstab(ii) = rho*(du(ii)-f(ii))+gp(ii)
        do jj=1, ndm
        rstab(ii) = rstab(ii)+rho*gu(ii,jj)*uw(jj)
        enddo
       enddo
      endif
      
        do i=1, nen
          i0 = (i-1) * ndf

c         stabilization of momentum equation
          if (stabm) 
     *     call rupwind3d(shp,uw,gu,tau,rstab,drdui,rstabi,r,i,ndm,ndf)

c:::::::: STIFFNESS MATRIX FOR dx_{n+1} ::::::::::::::::::::::::::::::::::::::

          do ii=1, ndm
          k(ndf,ii) = r0
          enddo

          do j=1, nen
            j0 = (j-1) * ndm
            
c           time derivative term and convective acceleration
            fact = - shp(ndm+1,i) * rho * shp(ndm+1,j)
            
            do ii=1, ndm
            do jj=1, ndm
                k(ii,jj) =  fact * gu(ii,jj)
             enddo
           enddo 

c           stabilization of momentum equation
            if (stabm) then
            
            fact   = r0
              do jj=1, ndm
                fact   = fact + shp(jj,i) * uw(jj) * rho
             enddo
     
            do ii=1, ndm
            do jj=1, ndm
              k(ii,jj) = k(ii,jj) - tau(1) * shp(ndm+1,j) * 
     *                 (fact * gu(ii,jj) + rstab(ii) * shp(jj,i))
             enddo
           enddo 
           
           
              fact   = - tau(2) * rho * shp(ndm+1,j)
              do ii=1, ndm
              k(ndm+1,ii) = r0
              do jj=1, ndm
              k(ndm+1,ii) = k(ndm+1,ii)+ 
     *                      fact * shp(jj,i) * gu(jj,ii)
             enddo
             enddo 
             
              

              do ii=1, ndm
              fact   = dtau(ii,1) / tau(1)
              do jj=1, ndm
              k(jj,ii) = k(jj,ii)- 
     *                      fact * rstabi(jj)
             enddo
             enddo 


              do ii=1, ndm
              k(ndm+1,ii) = k(ndm+1,ii)- 
     *                      dtau(ii,2) / tau(2) * rstabi(ndm+1)
             enddo

 
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


c==============================================================================
 11   continue ! isw 11 -> project to nodes vort(u)
               !     12 -> project to nodes div(u)
               !     13 -> project to nodes ||grad(u)||
c-----------------------------------------------------------------------------


c     get shape functions at centroid      
      xi(1,1) = r1d4
      xi(2,1) = r1d4
      xi(3,1) = r1d4
      wgp(1)  = r1d6
      
      call compshp3D(shp,detJ,xi(1,1),xa,nen)
    
      fact = r0
      
      if (isw.eq.11) then
      ii=3
      jj=2
      else if (isw.eq.12) then
      ii=1
      jj=3
      else
      ii=2
      jj=1
      endif

        do i=1, nen
          fact = fact + ul(ii,i)*shp(jj,i) - ul(jj,i)*shp(ii,i)
        enddo
        
c...  project to nodes

        do i=1, nen
          p(i) = fact
        enddo

      return
c==============================================================================
 14   continue 
       call prgwarning(1,'element3DStabIncompFluid','what ?')
 
 
      return

      end
















