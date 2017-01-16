c
c  subroutine rupwind
c  subroutine kupwind
c  subroutine kupwindw
c
c  subroutine eval_tau
c  subroutine eval_dtau
c  subroutine eval_dtauLD
c
c  function GQ
c
c  subroutine shp_space
c  subroutine shp_spaceGP
c
c  subroutine comp_hydrostatic   (d,ul,xl,s,p,ndf,ndm,nst)
c  subroutine comp_hydrostaticLD6(d,ul,xl,s,p,ndf,ndm,nst)
c  subroutine comp_hydrostaticLD5(d,ul,xl,s,p,ndf,ndm,nst)
c
c
c
c============================================================================
      subroutine rupwind(shp,u,gu,tau,rstab,drdui,rstabi,r,i)
c----------------------------------------------------------------------------
c  output: drdui       -  weighting term
c          rstabi      -  stabilization term
c----------------------------------------------------------------------------
      implicit none
      integer i
      double precision shp(3,3), u(3), gu(2,2), tau(*), 
     *                 rstab(2), drdui(2,2), rstabi(3), r(3), 
     *                 fact, r0
      data r0 / 0.d0 /

      fact = (shp(1,i) * u(1) + shp(2,i) * u(2))
      drdui(1,1) = fact
      drdui(1,2) = r0
      drdui(2,1) = r0
      drdui(2,2) = fact

      rstabi(1) = tau(1) * drdui(1,1) * rstab(1)
      rstabi(2) = tau(1) * drdui(2,2) * rstab(2)
      rstabi(3) = tau(2) * (shp(1,i)*rstab(1) + shp(2,i)*rstab(2))
 
      r(1) = r(1) + rstabi(1)
      r(2) = r(2) + rstabi(2)
      r(3) = r(3) + rstabi(3)

      return
      end

c============================================================================
      subroutine kupwind(shp,tau,dtau,rstab,rstabi,drdui,drduj,k,i,j)
c----------------------------------------------------------------------------
c  output: k   -  stabilization stiffness matrix for nodes(i,j)
c----------------------------------------------------------------------------
      implicit none
      integer i, j
      double precision shp(3,3), tau(*), dtau(2,*), 
     *                 drdui(2,2), drduj(2,2), rstab(2), 
     *                 rstabi(3), k(3,3), fact

      fact   = dtau(1,1) / tau(1)
      k(1,1) = k(1,1) + fact * rstabi(1)  
      k(2,1) = k(2,1) + fact * rstabi(2)  
      fact   = dtau(2,1) / tau(1)
      k(1,2) = k(1,2) + fact * rstabi(1)  
      k(2,2) = k(2,2) + fact * rstabi(2)

      fact   = tau(1) * shp(3,j) 
      k(1,1) = k(1,1) + fact * shp(1,i) * rstab(1)   
      k(1,2) = k(1,2) + fact * shp(2,i) * rstab(1)
      k(2,1) = k(2,1) + fact * shp(1,i) * rstab(2)
      k(2,2) = k(2,2) + fact * shp(2,i) * rstab(2)

      k(1,1) = k(1,1) + tau(1) * drdui(1,1) * drduj(1,1)
      k(1,2) = k(1,2) + tau(1) * drdui(1,1) * drduj(1,2)
      k(1,3) = k(1,3) + tau(1) * drdui(1,1) * shp(1,j) 
      k(2,1) = k(2,1) + tau(1) * drdui(2,2) * drduj(2,1)
      k(2,2) = k(2,2) + tau(1) * drdui(2,2) * drduj(2,2)
      k(2,3) = k(2,3) + tau(1) * drdui(2,2) * shp(2,j) 

      k(3,1) = k(3,1) + dtau(1,2) / tau(2) * rstabi(3)
      k(3,2) = k(3,2) + dtau(2,2) / tau(2) * rstabi(3)

      k(3,1) = k(3,1) + tau(2) * (shp(1,i)*drduj(1,1) 
     *                                       + shp(2,i)*drduj(2,1))
      k(3,2) = k(3,2) + tau(2) * (shp(1,i)*drduj(1,2) 
     *                                       + shp(2,i)*drduj(2,2))
      k(3,3) = tau(2) * (shp(1,i)*shp(1,j) + shp(2,i)*shp(2,j))

      return
      end

c============================================================================
      subroutine kupwindw(shp,tau,dtaudw,rstab,rstabi,
     *                    drdui,drduidw,drdxj,Dshp,k,i,j)
c----------------------------------------------------------------------------
c  output: k   -  stabilization stiffness matrix for nodes(i,j)
c----------------------------------------------------------------------------
      implicit none
      integer i, j
      double precision shp(3,3), tau(*), dtaudw(2,3,2), 
     *                 drdui(2,2), drdxj(2,2), rstab(2), 
     *                 rstabi(3), k(3,3), drduidw(2),
     *                 Dshp(2,3,2,3), fact

      fact   = dtaudw(1,j,1) / tau(1)       
      k(1,1) = k(1,1) + fact * rstabi(1) 
      k(2,1) = k(2,1) + fact * rstabi(2)  
      fact   = dtaudw(2,j,1) / tau(1)        
      k(1,2) = k(1,2) + fact * rstabi(1)
      k(2,2) = k(2,2) + fact * rstabi(2) 

      k(1,1) = k(1,1) + tau(1) * drduidw(1) * rstab(1) 
      k(1,2) = k(1,2) + tau(1) * drduidw(2) * rstab(1)
      k(2,1) = k(2,1) + tau(1) * drduidw(1) * rstab(2)
      k(2,2) = k(2,2) + tau(1) * drduidw(2) * rstab(2)

      k(1,1) = k(1,1) + tau(1) * drdxj(1,1) * drdui(1,1)
      k(1,2) = k(1,2) + tau(1) * drdxj(1,2) * drdui(1,1)
      k(2,1) = k(2,1) + tau(1) * drdui(2,2) * drdxj(2,1)
      k(2,2) = k(2,2) + tau(1) * drdui(2,2) * drdxj(2,2)

      k(3,1) = k(3,1) + dtaudw(1,j,2) / tau(2) * rstabi(3)
      k(3,2) = k(3,2) + dtaudw(2,j,2) / tau(2) * rstabi(3)

      k(3,1) = k(3,1) + tau(2) * ( shp(1,i)*drdxj(1,1) 
     *                          +  shp(2,i)*drdxj(2,1)
     *                        + Dshp(1,i,1,j)*rstab(1)
     *                        + Dshp(2,i,1,j)*rstab(2))
      k(3,2) = k(3,2) + tau(2) * (shp(1,i)*drdxj(1,2) 
     *                          + shp(2,i)*drdxj(2,2)
     *                        + Dshp(1,i,2,j)*rstab(1)
     *                        + Dshp(2,i,2,j)*rstab(2))

      return
      end

c============================================================================
      subroutine eval_tau(u,h,rho,mu,dt,nRC,z,beta,tau)
c----------------------------------------------------------------------------
c  stabilisation parameters tau_u, tau_p and tau_c
c
c                         beta1i
c  zi(Re) = ---------------------------------- ,  i=1,2,3
c            sqrt(1 + (beta1i/(beta2i Re))^2)
c
c             h
c  tau_u = ------- z1(Re)
c           2 |u|
c
c               h
c  tau_p = ----------- z2(Re)
c           2 |u| rho
c
c  tau_c = h |u| z3(Re)
c
c----------------------------------------------------------------------------

      implicit none
      integer          i
      double precision u(2), h, rho, mu, dt, nRC(*), z(*), 
     *                 beta(2,*), tau(*), 
     *                 Re, nu, fact

      nu       = sqrt(u(1)*u(1) + u(2)*u(2))   !   |u|
      Re       = nu * h * rho / (2.d0 * mu)    !  Re(u)  Reynolds no.
 !    C        = nu * dt / h                   !  C(u)   Courant no.

      nRC(1)   = nu
      nRC(2)   = Re
 !    nRC(3)   = C

      if (Re.lt.1.d-16) then  
        fact   = h * h / (4.d0*mu)
        tau(1) = fact * beta(2,1) * rho    !  tau_u
        tau(2) = fact * beta(2,2)          !  tau_p
        tau(3) = h * nu * Re * beta(2,3)   !  tau_c
      else
        do i=1, 3
          z(i) = beta(1,i) / sqrt(1.d0 + (beta(1,i)/(beta(2,i)*Re))**2)  
        enddo
        fact   = h / (2.d0 * nu)
        tau(1) = fact * z(1)               !  tau_u
        tau(2) = fact * z(2) / rho         !  tau_p
        tau(3) = h * nu * z(3)             !  tau_c
      endif

      return
      end

c============================================================================
      subroutine eval_dtaudh(u,h,rho,mu,dt,nRC,z,beta,dtaudh)
c----------------------------------------------------------------------------
c  derivatives of eval_tau with respect to h
c----------------------------------------------------------------------------
      implicit none
      integer          i
      double precision uw(2), u(2), h, rho, mu, dt, nRC(*), z(*),   
     *                 beta(2,*), dtaudh(*),  
     *                 Re, nu, dz(3), dRe, fact

      nu          = nRC(1)
      Re          = nRC(2)
      dRe         = Re / h

      if (Re.lt.1.d-16) then                           
        dtaudh(1) = 0.5d0 * h / mu * beta(2,1) * rho       
        dtaudh(2) = 0.5d0 * h / mu * beta(2,2)
        dtaudh(3) = nu * beta(2,3) * (Re + h * dRe)
      else
        do i=1, 3
          dz(i)   = (z(i)/Re)**3 / (beta(2,i)*beta(2,i))
        enddo
        fact      = h / (2.d0 * nu) * dRe
        dtaudh(1) =  fact * dz(1) + z(1) / (2.d0*nu)
        dtaudh(2) = (fact * dz(2) + z(2) / (2.d0*nu)) / rho
        dtaudh(3) =  h * nu * dz(3) * dRe + nu * z(3)
      endif

      return
      end

c============================================================================
      subroutine eval_dtau(u,h,rho,mu,dt,nRC,z,beta,dtau,mult)
c----------------------------------------------------------------------------
c  derivatives of eval_tau 
c----------------------------------------------------------------------------
      implicit none
      integer          i
      double precision u(*), h, rho, mu, dt, nRC(*), z(*),  
     *                 beta(2,*), mult, dtau(2,*),  
     *                 nu, Re, dRe, dnu(2), dz(3), fact, fact1, dfact

      nu          = nRC(1)
      Re          = nRC(2)

c      if (nu.lt.1.d-10) then
c        call pzero(dtau,6)
c        return
c      else
c        dRe         = Re / nu
c        dnu(1)      = u(1) / nu
c        dnu(2)      = u(2) / nu
c      endif

      if (Re.lt.1.d-16) then                           
        call pzero(dtau,4)
c        fact      = h * beta(2,3) * (Re + nu * dRe) * mult
c        dtau(1,3) = fact * dnu(1)
c        dtau(2,3) = fact * dnu(2)
      else
 
        dRe         = Re / nu
        dnu(1)      = u(1) / nu
        dnu(2)      = u(2) / nu

        do i=1, 3
          dz(i)   = (z(i)/Re)**3 / (beta(2,i)*beta(2,i))
        enddo
        fact      = h / (2.d0 * nu)
        dfact     = - h / (2.d0 * nu * nu)
        fact1     = (dfact * z(1) + fact * dz(1) * dRe) * mult
        dtau(1,1) = fact1 * dnu(1)
        dtau(2,1) = fact1 * dnu(2)
        fact1     = (dfact * z(2) + fact * dz(2) * dRe) * mult / rho
        dtau(1,2) = fact1 * dnu(1)
        dtau(2,2) = fact1 * dnu(2)
        fact1     = (h * z(3) + h * nu * dz(3) * dRe) * mult
        dtau(1,3) = fact1 * dnu(1)
        dtau(2,3) = fact1 * dnu(2)
      endif

      return
      end

c============================================================================
      subroutine eval_dtauLD(u,h,rho,mu,nu,dt,Re,C,zu,zp,beta1,beta2,
     *                       dtauu,dtaup,multt,multtn)
c----------------------------------------------------------------------------
c  derivatives of eval_tau for LD
c----------------------------------------------------------------------------
      implicit none
      double precision u(2), h, rho, mu, nu, dt, C, Re, zu, zp,  
     *                 beta1, beta2, multt, multtn, dtauu(4), dtaup(4), 
     *                 dnu(2), dz(2), fact, fact1, fact2, hlp1, hlp2, 
     *                 r0, r2, r45

      data r0, r2, r45   / 0.d0, 2.d0, 4.5d0 /

      if (Re.lt.1.d-16) then                  
        dtauu(1) = r0        
        dtaup(1) = r0
        dtauu(2) = r0        
        dtaup(2) = r0
      else

        dnu(1) = u(1) / nu
        dnu(2) = u(2) / nu

c... tau u .............................

        fact     = h*(r45*rho/(Re*Re*Re*mu))
        fact1    = zu*zu*zu*fact
        fact2    = h / (r2*nu)
        dz(1)    = fact1 * dnu(1)
        dz(2)    = fact1 * dnu(2)
        hlp1     = fact2 * (dz(1) - dnu(1) * zu / nu)
        hlp2     = fact2 * (dz(2) - dnu(2) * zu / nu)
        dtauu(1) = hlp1 * multt
        dtauu(2) = hlp2 * multt 
        dtauu(3) = hlp1 * multtn
        dtauu(4) = hlp2 * multtn 

c... tau p .............................

        fact1    = zp*zp*zp*fact/(beta1*beta1)
        fact2    = fact2 / rho
        dz(1)    = fact1 * dnu(1)
        dz(2)    = fact1 * dnu(2)
        hlp1     = fact2 * (dz(1) - dnu(1) * zp / nu)
        hlp2     = fact2 * (dz(2) - dnu(2) * zp / nu)
        dtaup(1) = hlp1 * multt
        dtaup(2) = hlp2 * multt 
        dtaup(3) = hlp1 * multtn
        dtaup(4) = hlp2 * multtn 

      endif
   
      return
      end

c============================================================================
      double precision function GQ(l,i,lint)
c----------------------------------------------------------------------------
c  GAUSS INTEGRATION 
c----------------------------------------------------------------------------
      implicit none
      integer          l, i, lint
      double precision fact

      if (lint.eq.1) then
        GQ = 0.5d0
        if (i.eq.4) GQ = 1.0d0
      elseif (lint.eq.2) then
        fact = 1.d0
        if (i.eq.3) fact = - fact
        if (l.eq.1) then
          GQ = (3.d0 - fact * sqrt(3.d0)) / 6.d0
        elseif (l.eq.2) then
          GQ = (3.d0 + fact * sqrt(3.d0)) / 6.d0
        endif
        if (i.eq.4) GQ = 0.5d0
      elseif (lint.eq.3) then
        fact = 1.d0
        if (i.eq.3) fact = - fact
        if (l.eq.1) then
          GQ = (5.d0 - fact * sqrt(15.d0)) * 0.1d0
        elseif (l.eq.2) then
          GQ = 0.5d0
        elseif (l.eq.3) then
          GQ = (5.d0 + fact * sqrt(15.d0)) * 0.1d0
        endif
        if (i.eq.4) then
          if (l.eq.2) then
            GQ = 4.d0 / 9.d0
          else
            GQ = 5.d0 / 18.d0
          endif
        endif
      endif
      return
      end

c============================================================================
      subroutine shp_space(xa,shp,Dshp)
c----------------------------------------------------------------------------
c  shp   SPATIAL SHAPE FUNCTIONS 
c
c  Dshp  DERIVATIVE OF SHAPE FUNCTIONS WITH RESPECT TO NODAL COORDINATES
c
c     evaluated in centroid of triangle
c----------------------------------------------------------------------------
      implicit none
      integer          i, j, k
      double precision xa(2,3), shp(3,3), Dshp(2,3,2,3),
     *                 fact, fact1, c(3), Dfact(2,3), Dc(3,2,3),
     *                 x1, x2, r1, r1d3
      r1   = 1.d0
      r1d3 = r1 / 3.d0
      
      do i=1, 3   !  loop over nodes
        j = i+1
        k = i+2
        if (j.gt.3) j = j - 3
        if (k.gt.3) k = k - 3

        x1         = r1d3 * (xa(1,1) + xa(1,2) + xa(1,3))
        x2         = r1d3 * (xa(2,1) + xa(2,2) + xa(2,3))

        fact       = r1 / (  (xa(2,i)-xa(2,j)) * (xa(1,j)-xa(1,k))
     *                      -(xa(2,j)-xa(2,k)) * (xa(1,i)-xa(1,j)) )

        c(1)       = - fact * (xa(2,j)         - xa(2,k)        )
        c(2)       =   fact * (xa(1,j)         - xa(1,k)        )
c        c(3)       =   fact * (xa(2,j)*xa(1,k) - xa(2,k)*xa(1,j))

        shp(1,i)   = c(1)
        shp(2,i)   = c(2)
        shp(3,i)   = r1d3   ! = c(1)*x1 + c(2)*x2 + c(3)
            
        fact1      = - fact * fact

        Dfact(1,i) = fact1 * (-xa(2,j)+xa(2,k))
        Dfact(1,j) = fact1 * ( xa(2,i)-xa(2,j) + xa(2,j)-xa(2,k))
        Dfact(1,k) = fact1 * (-xa(2,i)+xa(2,j))

        Dfact(2,i) = fact1 * ( xa(1,j)-xa(1,k))
        Dfact(2,j) = fact1 * (-xa(1,j)+xa(1,k) - xa(1,i)+xa(1,j)) 
        Dfact(2,k) = fact1 * ( xa(1,i)-xa(1,j)) 

        fact1      = - (xa(2,j) - xa(2,k))
        Dc(1,1,i)  = Dfact(1,i) * fact1
        Dc(1,1,j)  = Dfact(1,j) * fact1
        Dc(1,1,k)  = Dfact(1,k) * fact1

        fact1      = xa(1,j) - xa(1,k)
        Dc(2,1,i)  = Dfact(1,i) * fact1
        Dc(2,1,j)  = Dfact(1,j) * fact1 + fact
        Dc(2,1,k)  = Dfact(1,k) * fact1 - fact

c        fact1      = xa(2,j)*xa(1,k) - xa(2,k)*xa(1,j)
c        Dc(3,1,i)  = Dfact(1,i) * fact1
c        Dc(3,1,j)  = Dfact(1,j) * fact1 - fact*xa(2,k)
c        Dc(3,1,k)  = Dfact(1,k) * fact1 + fact*xa(2,j)

        fact1      = - (xa(2,j) - xa(2,k))
        Dc(1,2,i)  = Dfact(2,i) * fact1
        Dc(1,2,j)  = Dfact(2,j) * fact1 - fact
        Dc(1,2,k)  = Dfact(2,k) * fact1 + fact

        fact1      = xa(1,j) - xa(1,k)
        Dc(2,2,i)  = Dfact(2,i) * fact1
        Dc(2,2,j)  = Dfact(2,j) * fact1
        Dc(2,2,k)  = Dfact(2,k) * fact1

c        fact1      = xa(2,j)*xa(1,k) - xa(2,k)*xa(1,j)
c        Dc(3,2,i)  = Dfact(2,i) * fact1 
c        Dc(3,2,j)  = Dfact(2,j) * fact1 + fact*xa(1,k)
c        Dc(3,2,k)  = Dfact(2,k) * fact1 - fact*xa(1,j)

        Dshp(1,i,1,i) = Dc(1,1,i)
        Dshp(2,i,1,i) = Dc(2,1,i)
c        Dshp(3,i,1,i) = Dc(1,1,i)*x1+Dc(2,1,i)*x2+Dc(3,1,i)+c(1)*r1d3

        Dshp(1,i,1,j) = Dc(1,1,j)
        Dshp(2,i,1,j) = Dc(2,1,j)
c        Dshp(3,i,1,j) = Dc(1,1,j)*x1+Dc(2,1,j)*x2+Dc(3,1,j)+c(1)*r1d3

        Dshp(1,i,1,k) = Dc(1,1,k)
        Dshp(2,i,1,k) = Dc(2,1,k)
c        Dshp(3,i,1,k) = Dc(1,1,k)*x1+Dc(2,1,k)*x2+Dc(3,1,k)+c(1)*r1d3

        Dshp(1,i,2,i) = Dc(1,2,i)
        Dshp(2,i,2,i) = Dc(2,2,i)
c        Dshp(3,i,2,i) = Dc(1,2,i)*x1+Dc(2,2,i)*x2+Dc(3,2,i)+c(2)*r1d3

        Dshp(1,i,2,j) = Dc(1,2,j)
        Dshp(2,i,2,j) = Dc(2,2,j)
c        Dshp(3,i,2,j) = Dc(1,2,j)*x1+Dc(2,2,j)*x2+Dc(3,2,j)+c(2)*r1d3

        Dshp(1,i,2,k) = Dc(1,2,k)
        Dshp(2,i,2,k) = Dc(2,2,k)
c        Dshp(3,i,2,k) = Dc(1,2,k)*x1+Dc(2,2,k)*x2+Dc(3,2,k)+c(2)*r1d3

      enddo
      return
      end

c============================================================================
      subroutine shp_spaceGP(x,shp,Dshp,sg)
c----------------------------------------------------------------------------
c  shp   SPATIAL SHAPE FUNCTIONS 
c
c  Dshp  DERIVATIVE OF SHAPE FUNCTIONS WITH RESPECT TO NODAL COORDINATES
c
c----------------------------------------------------------------------------
      implicit none
      integer          i, j, k
      double precision x(2,3), shp(3,3), Dshp(2,3,2,3), sg(3), fact
      
      fact = 1.d0 / (  x(1,2)*x(2,1) - x(1,3)*x(2,1) - x(1,1)*x(2,2) 
     *               + x(1,3)*x(2,2) + x(1,1)*x(2,3) - x(1,2)*x(2,3))

      shp(1,1)   = fact * (x(2,3) - x(2,2))
      shp(2,1)   = fact * (x(1,2) - x(1,3))
      shp(3,1)   = sg(1)

      shp(1,2)   = fact * (x(2,1) - x(2,3))
      shp(2,2)   = fact * (x(1,3) - x(1,1))
      shp(3,2)   = sg(2)

      shp(1,3)   = fact * (x(2,2) - x(2,1))
      shp(2,3)   = fact * (x(1,1) - x(1,2))
      shp(3,3)   = sg(3)
            
      fact =    x(1,3)*(x(2,2)-x(2,1))
     *        + x(1,2)*(x(2,1)-x(2,3))
     *        + x(1,1)*(x(2,3)-x(2,2))
      fact = 1.d0 / (fact * fact)

      Dshp(1,1,1,1) = fact * (x(2,2)-x(2,3)) * (x(2,3)-x(2,2))
      Dshp(1,1,2,1) = fact * (x(2,2)-x(2,3)) * (x(1,2)-x(1,3))
      Dshp(1,1,1,2) = fact * (x(2,1)-x(2,3)) * (x(2,2)-x(2,3))
      Dshp(1,1,2,2) = fact * (x(2,1)-x(2,3)) * (x(1,3)-x(1,2))
      Dshp(1,1,1,3) = fact * (x(2,2)-x(2,1)) * (x(2,2)-x(2,3))
      Dshp(1,1,2,3) = fact * (x(2,1)-x(2,2)) * (x(1,2)-x(1,3))

      Dshp(2,1,1,1) = fact * (x(1,2)-x(1,3)) * (x(2,2)-x(2,3))
      Dshp(2,1,2,1) = fact * (x(1,2)-x(1,3)) * (x(1,3)-x(1,2))
      Dshp(2,1,1,2) = fact * (x(1,3)-x(1,1)) * (x(2,2)-x(2,3))
      Dshp(2,1,2,2) = fact * (x(1,1)-x(1,3)) * (x(1,2)-x(1,3))
      Dshp(2,1,1,3) = fact * (x(1,1)-x(1,2)) * (x(2,2)-x(2,3))
      Dshp(2,1,2,3) = fact * (x(1,2)-x(1,1)) * (x(1,2)-x(1,3))

      Dshp(1,2,1,1) = fact * (x(2,1)-x(2,3)) * (x(2,2)-x(2,3))
      Dshp(1,2,2,1) = fact * (x(2,2)-x(2,3)) * (x(1,3)-x(1,1))
      Dshp(1,2,1,2) = fact * (x(2,1)-x(2,3)) * (x(2,3)-x(2,1))
      Dshp(1,2,2,2) = fact * (x(2,1)-x(2,3)) * (x(1,1)-x(1,3))
      Dshp(1,2,1,3) = fact * (x(2,1)-x(2,2)) * (x(2,1)-x(2,3))
      Dshp(1,2,2,3) = fact * (x(2,1)-x(2,2)) * (x(1,3)-x(1,1))

      Dshp(2,2,1,1) = fact * (x(1,3)-x(1,2)) * (x(2,1)-x(2,3))
      Dshp(2,2,2,1) = fact * (x(1,1)-x(1,3)) * (x(1,2)-x(1,3))
      Dshp(2,2,1,2) = fact * (x(1,1)-x(1,3)) * (x(2,1)-x(2,3))
      Dshp(2,2,2,2) = fact * (x(1,1)-x(1,3)) * (x(1,3)-x(1,1))
      Dshp(2,2,1,3) = fact * (x(1,2)-x(1,1)) * (x(2,1)-x(2,3))
      Dshp(2,2,2,3) = fact * (x(1,1)-x(1,2)) * (x(1,1)-x(1,3))

      Dshp(1,3,1,1) = fact * (x(2,2)-x(2,1)) * (x(2,2)-x(2,3))
      Dshp(1,3,2,1) = fact * (x(2,2)-x(2,3)) * (x(1,1)-x(1,2))
      Dshp(1,3,1,2) = fact * (x(2,1)-x(2,2)) * (x(2,1)-x(2,3))
      Dshp(1,3,2,2) = fact * (x(2,1)-x(2,3)) * (x(1,2)-x(1,1))
      Dshp(1,3,1,3) = fact * (x(2,1)-x(2,2)) * (x(2,2)-x(2,1))
      Dshp(1,3,2,3) = fact * (x(2,1)-x(2,2)) * (x(1,1)-x(1,2))

      Dshp(2,3,1,1) = fact * (x(1,2)-x(1,3)) * (x(2,1)-x(2,2))
      Dshp(2,3,2,1) = fact * (x(1,2)-x(1,1)) * (x(1,2)-x(1,3))
      Dshp(2,3,1,2) = fact * (x(1,3)-x(1,1)) * (x(2,1)-x(2,2))
      Dshp(2,3,2,2) = fact * (x(1,1)-x(1,2)) * (x(1,1)-x(1,3))
      Dshp(2,3,1,3) = fact * (x(1,1)-x(1,2)) * (x(2,1)-x(2,2))
      Dshp(2,3,2,3) = fact * (x(1,1)-x(1,2)) * (x(1,2)-x(1,1))

      return
      end

c=============================================================================
      subroutine comp_hydrostatic(d,ul,xl,s,p,ndf,ndm,nst)
c----------------------------------------------------------------------------
c  
c     COMPUTE HYDROSTATIC PRESSURE DISTRIBUTION
c
c----------------------------------------------------------------------------
      implicit none
      integer          ndf, ndm, nst, i, j, l, j0, i0, ii
      double precision d(*), ul(ndf,*), xl(ndm,*), s(nst,*), p(*),
     *                 rho, f(2), fact, gp(2), dvol, shp(3,3), k(3,3), 
     *                 r(3), r1d3, r1, r0, r1d2
      data r0, r1d2, r1 / 0.d0, 0.5d0, 1.d0 /
      r1d3 = r1 / 3.d0

      rho   = d(12)
      f(1)  = d(13)
      f(2)  = d(14)
      dvol = r1d2 * ((xl(1,1)-xl(1,3))*(xl(2,2)-xl(2,3)) 
     *             - (xl(1,2)-xl(1,3))*(xl(2,1)-xl(2,3)))
      if (int(d(7)+.1).eq.1) 
     *  dvol = dvol * (xl(1,1)+xl(1,2)+xl(1,3))
      do i=1, 3  
        j = i+1
        l = i+2
        if (j.gt.3) j = j - 3
        if (l.gt.3) l = l - 3
        fact     = r1 / (  (xl(2,i)-xl(2,j)) * (xl(1,j)-xl(1,l))
     *                    -(xl(2,j)-xl(2,l)) * (xl(1,i)-xl(1,j)) )
        shp(1,i) = - fact * (xl(2,j) - xl(2,l))
        shp(2,i) =   fact * (xl(1,j) - xl(1,l))
        shp(3,i) =   r1d3
      enddo
      call pzero(gp,2)
      do i=1, 3
        do ii=1, 2
          gp(ii) = gp(ii) + ul(3,i) * shp(ii,i) 
        enddo
      enddo
      do i=1, 3
        i0 = (i-1)*3
        r(1) = ul(1,i)
        r(2) = ul(2,i)
        r(3) = shp(1,i)*(gp(1)-rho*f(1)) + shp(2,i)*(gp(2)-rho*f(2)) 
        do ii=1, 3
          p(i0+ii) = p(i0+ii) - r(ii) * dvol
        enddo
        do j=1, 3
          j0 = (j-1)*3
          if (i.eq.j) then
            k(1,1) = r1
            k(2,2) = r1
          else
            k(1,1) = r0
            k(2,2) = r0
          endif
          k(3,3) = shp(1,i)*shp(1,j) + shp(2,i)*shp(2,j)
          do ii=1, 3
            s(i0+ii,j0+ii) = s(i0+ii,j0+ii) + k(ii,ii) * dvol
          enddo
        enddo
      enddo
      return 
      end

c=============================================================================
      subroutine comp_hydrostaticLD6(d,ul,xl,s,p,ndf,ndm,nst)
c----------------------------------------------------------------------------
c  
c     COMPUTE HYDROSTATIC PRESSURE DISTRIBUTION
c
c----------------------------------------------------------------------------
      implicit none
      integer          ndf, ndm, nst, i, j, l, j0, i0, ii
      double precision d(*), ul(ndf,*), xl(ndm,*), s(nst,*), p(*),
     *                 rho, f(2), fact, gp(2), dvol, shp(3,3), k(6,6), 
     *                 r(6), r1d3, r1, r0, r1d2
      data r0, r1d2, r1 / 0.d0, 0.5d0, 1.d0 /
      r1d3 = r1 / 3.d0

      rho   = d(12)
      f(1)  = d(13)
      f(2)  = d(14)
      dvol = r1d2 * ((xl(1,1)-xl(1,3))*(xl(2,2)-xl(2,3)) 
     *             - (xl(1,2)-xl(1,3))*(xl(2,1)-xl(2,3)))
      if (int(d(7)+.1).eq.1) 
     *  dvol = dvol * (xl(1,1)+xl(1,2)+xl(1,3))
      do i=1, 3  
        j = i+1
        l = i+2
        if (j.gt.3) j = j - 3
        if (l.gt.3) l = l - 3
        fact     = r1 / (  (xl(2,i)-xl(2,j)) * (xl(1,j)-xl(1,l))
     *                    -(xl(2,j)-xl(2,l)) * (xl(1,i)-xl(1,j)) )
        shp(1,i) = - fact * (xl(2,j) - xl(2,l))
        shp(2,i) =   fact * (xl(1,j) - xl(1,l))
        shp(3,i) =   r1d3
      enddo
      call pzero(gp,2)
      do i=1, 3
        do ii=1, 2
          gp(ii) = gp(ii) + ul(3,i) * shp(ii,i) 
        enddo
      enddo
      do i=1, 3
        i0 = (i-1)*6
        r(1) = ul(1,i)
        r(2) = ul(2,i)
        r(3) = shp(1,i)*(gp(1)-rho*f(1)) + shp(2,i)*(gp(2)-rho*f(2)) 
        r(4) = ul(4,i)
        r(5) = ul(5,i)
        r(6) = shp(1,i)*(gp(1)-rho*f(1)) + shp(2,i)*(gp(2)-rho*f(2)) 
        do ii=1, 6
          p(i0+ii) = p(i0+ii) - r(ii) * dvol
        enddo
        do j=1, 3
          j0 = (j-1)*6
          if (i.eq.j) then
            k(1,1) = r1
            k(2,2) = r1
            k(4,4) = r1
            k(5,5) = r1
          else
            k(1,1) = r0
            k(2,2) = r0
            k(4,4) = r0
            k(5,5) = r0
          endif
          k(3,3) = shp(1,i)*shp(1,j) + shp(2,i)*shp(2,j)
          k(6,6) = shp(1,i)*shp(1,j) + shp(2,i)*shp(2,j)
          do ii=1, 6
            s(i0+ii,j0+ii) = s(i0+ii,j0+ii) + k(ii,ii) * dvol
          enddo
        enddo
      enddo
      return 
      end

c=============================================================================
      subroutine comp_hydrostaticLD5(d,ul,xl,s,p,ndf,ndm,nst)
c----------------------------------------------------------------------------
c  
c     COMPUTE HYDROSTATIC PRESSURE DISTRIBUTION
c
c----------------------------------------------------------------------------
      implicit none
      integer          ndf, ndm, nst, i, j, l, j0, i0, ii
      double precision d(*), ul(ndf,*), xl(ndm,*), s(nst,*), p(*),
     *                 rho, f(2), fact, gp(2), dvol, shp(3,3), k(5,5), 
     *                 r(5), r1d3, r1, r0, r1d2
      data r0, r1d2, r1 / 0.d0, 0.5d0, 1.d0 /
      r1d3 = r1 / 3.d0

      rho   = d(12)
      f(1)  = d(13)
      f(2)  = d(14)
      dvol = r1d2 * ((xl(1,1)-xl(1,3))*(xl(2,2)-xl(2,3)) 
     *             - (xl(1,2)-xl(1,3))*(xl(2,1)-xl(2,3)))
      if (int(d(7)+.1).eq.1) 
     *  dvol = dvol * (xl(1,1)+xl(1,2)+xl(1,3))
      do i=1, 3  
        j = i+1
        l = i+2
        if (j.gt.3) j = j - 3
        if (l.gt.3) l = l - 3
        fact     = r1 / (  (xl(2,i)-xl(2,j)) * (xl(1,j)-xl(1,l))
     *                    -(xl(2,j)-xl(2,l)) * (xl(1,i)-xl(1,j)) )
        shp(1,i) = - fact * (xl(2,j) - xl(2,l))
        shp(2,i) =   fact * (xl(1,j) - xl(1,l))
        shp(3,i) =   r1d3
      enddo
      call pzero(gp,2)
      do i=1, 3
        do ii=1, 2
          gp(ii) = gp(ii) + ul(3,i) * shp(ii,i) 
        enddo
      enddo
      do i=1, 3
        i0 = (i-1)*5
        r(1) = ul(1,i)
        r(2) = ul(2,i)
        r(3) = shp(1,i)*(gp(1)-rho*f(1)) + shp(2,i)*(gp(2)-rho*f(2)) 
        r(4) = ul(4,i)
        r(5) = ul(5,i)
        do ii=1, 5
          p(i0+ii) = p(i0+ii) - r(ii) * dvol
        enddo
        do j=1, 3
          j0 = (j-1)*5
          if (i.eq.j) then
            k(1,1) = r1
            k(2,2) = r1
            k(4,4) = r1
            k(5,5) = r1
          else
            k(1,1) = r0
            k(2,2) = r0
            k(4,4) = r0
            k(5,5) = r0
          endif
          k(3,3) = shp(1,i)*shp(1,j) + shp(2,i)*shp(2,j)
          do ii=1, 5
            s(i0+ii,j0+ii) = s(i0+ii,j0+ii) + k(ii,ii) * dvol
          enddo
        enddo
      enddo
      return 
      end

c=============================================================================





