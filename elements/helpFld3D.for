c
c  subroutine rupwind3d
c  subroutine kupwind3d
c  subroutine kupwindw3d
c
c  subroutine eval_tau3d
c  subroutine eval_dtau3d
c  subroutine shp_spaceGP3d
c

c============================================================================
      subroutine rupwind3d(shp,u,gu,tau,rstab,drdui,rstabi,r,i,ndm)
c----------------------------------------------------------------------------
c  output: drdui       -  weighting term
c          rstabi      -  stabilization term
c----------------------------------------------------------------------------
      implicit none
      integer i, ndm, j
      double precision shp(ndm+1,*),u(ndm+1),gu(ndm,ndm), tau(*), 
     *                 rstab(ndm), drdui(ndm,ndm), rstabi(ndm+1),  
     *                 r(ndm+1),fact, r0
      data r0 / 0.d0 /

      fact = r0
      do j=1, ndm
      fact = fact + shp(j,i) * u(j)                 
      enddo
      
      call pzero(drdui,ndm*ndm)
      do j=1, ndm
      drdui(j,j) = fact
      enddo
            
            
      do j=1, ndm
      rstabi(j) = tau(1) * drdui(j,j) * rstab(j)
      enddo
      
      
      fact = r0
      do j=1, ndm
      fact = fact + shp(j,i) * rstab(j)                 
      enddo
      rstabi(ndm+1) = tau(2) * fact
     
      
      do j=1, ndm+1
      r(j) = r(j) + rstabi(j)
      enddo
 

      return
      end

c============================================================================
      subroutine kupwind3d(shp,tau,dtau,rstab,rstabi,
     *                   drdui,drduj,k,i,j,ndm)
c----------------------------------------------------------------------------
c  output: k   -  stabilization stiffness matrix for nodes(i,j)
c----------------------------------------------------------------------------
      implicit none
      integer i, j, ndm, ii, jj
      double precision shp(ndm+1,*), tau(*), dtau(ndm,3), 
     *                 drdui(ndm,ndm), drduj(ndm,ndm), rstab(ndm), 
     *                 rstabi(ndm+1), k(ndm+1,ndm+1), fact, r0
      
      data r0 / 0.d0 /


      do ii=1, ndm
      fact   = dtau(ii,1) / tau(1)
      do jj=1, ndm
      k(jj,ii) = k(jj,ii) + fact * rstabi(jj)  
      enddo
      enddo
      

      fact   = tau(1) * shp(ndm+1,j) 
      
      do ii=1, ndm
      do jj=1, ndm
       k(ii,jj) = k(ii,jj) + fact * shp(jj,i) * rstab(ii)
     *                     + tau(1) * drdui(ii,ii) * drduj(ii,jj)

      enddo
      k(ii,ndm+1) = k(ii,ndm+1) + tau(1) * drdui(ii,ii) * shp(ii,j)   
      k(ndm+1,ii) = k(ndm+1,ii) + dtau(ii,2) / tau(2) * rstabi(ndm+1)
      enddo
      
      k(ndm+1,ndm+1) = r0
      do ii=1, ndm
      do jj=1, ndm
      k(ndm+1,ii) = k(ndm+1,ii) + tau(2) * shp(jj,i)*drduj(jj,ii) 
      enddo
      k(ndm+1,ndm+1) = k(ndm+1,ndm+1) +tau(2) * shp(ii,i)*shp(ii,j)
      enddo




      return
      end



c============================================================================
      subroutine kupwindw3d(shp,tau,dtaudw,rstab,rstabi,
     *                    drdui,drduidw,drdxj,Dshp,k,i,j,ndm,nen)
c----------------------------------------------------------------------------
c  output: k   -  stabilization stiffness matrix for nodes(i,j)
c----------------------------------------------------------------------------
      implicit none
      integer i, j, ii, jj, ndm, nen
      double precision shp(ndm+1,nen), tau(*), dtaudw(ndm,nen,3), 
     *                 drdui(ndm,ndm), drdxj(ndm,ndm), rstab(ndm), 
     *                 rstabi(ndm+1), k(ndm+1,ndm+1), drduidw(ndm),
     *                 Dshp(ndm,nen,ndm,nen), fact


      do ii=1, ndm
         fact   = dtaudw(ii,j,1) / tau(1)
       do jj=1, ndm 
        k(jj,ii) = k(jj,ii) + fact * rstabi(jj) 
        enddo 
       enddo 

      do ii=1, ndm
       do jj=1, ndm 
        k(ii,jj) = k(ii,jj) + tau(1) * (drduidw(jj) * rstab(ii)+
     *                                 drdxj(ii,jj) * drdui(ii,ii))
       enddo 
      enddo 
      


      do jj=1, ndm 
      k(ndm+1,jj) = k(ndm+1,jj) + dtaudw(jj,j,2) 
     *   / tau(2) * rstabi(ndm+1)
      enddo  
      
      
      do ii=1, ndm 
      do jj=1, ndm 
      k(ndm+1,ii) = k(ndm+1,ii) + tau(2) * ( shp(jj,i)*drdxj(jj,ii)+ 
     *                                 Dshp(jj,i,ii,j)*rstab(jj))
      enddo  
      enddo        



      return
      end
      
      
c============================================================================
      subroutine eval_tau3d(u,h,rho,mu,dt,nRC,z,beta,tau,ndm)
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
      integer          i, ndm
      double precision u(ndm), h, rho, mu, dt, nRC(*), z(*), 
     *                 beta(2,*), tau(*), 
     *                 Re, nu, fact

      nu       =u(1)*u(1)
      do i=2, ndm
      nu       = nu +u(i)*u(i)                 
      enddo
      nu =sqrt(nu)                             !   |u|
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
      subroutine eval_dtau3d(u,h,rho,mu,dt,nRC,z,beta,dtau,mult,ndm)
c----------------------------------------------------------------------------
c  derivatives of eval_tau 
c----------------------------------------------------------------------------
      implicit none
      integer          i,ndm,j
      double precision u(*), h, rho, mu, dt, nRC(*), z(*),  
     *                 beta(2,3), mult, dtau(ndm,3),  
     *                 nu, Re, dRe, dnu(ndm), dz(3), fact, fact1, dfact

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

      !j=0

      !if (j.eq.0) then  !
      if (Re.lt.1.d-16) then     
        call pzero(dtau,ndm*3)
c        fact      = h * beta(2,3) * (Re + nu * dRe) * mult
c        dtau(1,3) = fact * dnu(1)
c        dtau(2,3) = fact * dnu(2)
      else
 
        dRe         = Re / nu
        
        do i=1, ndm
        dnu(i)      = u(i) / nu                 
        enddo
      
        do i=1, 3
          dz(i)   = (z(i)/Re)**3 / (beta(2,i)*beta(2,i))
        enddo
        fact      = h / (2.d0 * nu)
        dfact     = - h / (2.d0 * nu * nu)
        fact1     = (dfact * z(1) + fact * dz(1) * dRe) * mult
        do i=1, ndm
        dtau(i,1) = fact1 * dnu(i)
        enddo
        
        fact1     = (dfact * z(2) + fact * dz(2) * dRe) * mult / rho
        do i=1, ndm
        dtau(i,2) = fact1 * dnu(i)
        enddo
        
        fact1     = (h * z(3) + h * nu * dz(3) * dRe) * mult
        do i=1, ndm
        dtau(i,3) = fact1 * dnu(i)
        enddo
        
        
      endif

      return
      end


      
c============================================================================
      subroutine shp_spaceGP3D(x,shp,Dshp,sg)
c----------------------------------------------------------------------------
c  shp   SPATIAL SHAPE FUNCTIONS 
c
c  Dshp  DERIVATIVE OF SHAPE FUNCTIONS WITH RESPECT TO NODAL COORDINATES
c
c----------------------------------------------------------------------------
      implicit none
      integer          i, j, k,l, ndm, nen
      double precision x(3,4), shp(4,4), Dshp(3,4,3,4), sg(3), fact
     
      ndm=3
      nen=4
      
     
      fact =  1.d0 /(x(1,2)*x(2,3)*x(3,1) - x(1,2)*x(2,4)*x(3,1) - 
     *  x(1,1)*x(2,3)*x(3,2) + x(1,1)*x(2,4)*x(3,2) - 
     *  x(1,2)*x(2,1)*x(3,3) + x(1,1)*x(2,2)*x(3,3) - 
     *  x(1,1)*x(2,4)*x(3,3) + x(1,2)*x(2,4)*x(3,3) + 
     *  x(1,4)*(x(2,2)*x(3,1) - x(2,3)*x(3,1) - x(2,1)*x(3,2) + 
     *  x(2,3)*x(3,2) + x(2,1)*x(3,3) - x(2,2)*x(3,3)) + 
     *  x(1,2)*x(2,1)*x(3,4) - x(1,1)*x(2,2)*x(3,4) + 
     *  x(1,1)*x(2,3)*x(3,4) - x(1,2)*x(2,3)*x(3,4) + 
     *  x(1,3)*(x(2,4)*x(3,1) + x(2,1)*x(3,2) - 
     *  x(2,4)*x(3,2) - x(2,1)*x(3,4) + x(2,2)*(-x(3,1) + x(3,4))))
     

       shp(1,1)   = fact * (x(2,4)*(x(3,2) - x(3,3)) + x(2,2)*(x(3,3) - 
     *  x(3,4)) + x(2,3)*(-x(3,2) + x(3,4)))
       shp(2,1)   = fact * (x(1,4)*(-x(3,2) + x(3,3)) + x(1,3)*(x(3,2) -
     *  x(3,4)) + x(1,2)*(-x(3,3) + x(3,4)))
       shp(3,1)   = fact * (x(1,4)*(x(2,2) - x(2,3)) + x(1,2)*(x(2,3) - 
     *  x(2,4)) + x(1,3)*(-x(2,2) + x(2,4)))
       shp(4,1)   = 1.d0 - sg(1) - sg(2) - sg(3)
 
       shp(1,2)   = fact * (x(2,4)*(-x(3,1) + x(3,3)) + x(2,3)*(x(3,1) -
     *   x(3,4)) + x(2,1)*(-x(3,3) + x(3,4)))
       shp(2,2)   = fact * (x(1,4)*(x(3,1) - x(3,3)) + x(1,1)*(x(3,3) - 
     *  x(3,4)) + x(1,3)*(-x(3,1) + x(3,4)))
       shp(3,2)   = fact * (x(1,4)*(-x(2,1) + x(2,3)) + x(1,3)*(x(2,1) -
     *   x(2,4)) + x(1,1)*(-x(2,3) + x(2,4)))
       shp(4,2)   = sg(1)

         
       shp(1,3)   = fact * (x(2,4)*(x(3,1) - x(3,2)) + x(2,1)*(x(3,2) - 
     *  x(3,4)) + x(2,2)*(-x(3,1) + x(3,4)))
       shp(2,3)   = fact * (x(1,4)*(-x(3,1) + x(3,2)) + x(1,2)*(x(3,1) -
     *   x(3,4)) + x(1,1)*(-x(3,2) + x(3,4)))
       shp(3,3)   = fact * (x(1,4)*(x(2,1) - x(2,2)) + x(1,1)*(x(2,2) - 
     *  x(2,4)) + x(1,2)*(-x(2,1) + x(2,4)))
       shp(4,3)   = sg(2)

       shp(1,4)   = fact * (x(2,3)*(-x(3,1) + x(3,2)) + x(2,2)*(x(3,1) -
     *   x(3,3)) + x(2,1)*(-x(3,2) + x(3,3)))
       shp(2,4)   = fact * (x(1,3)*(x(3,1) - x(3,2)) + x(1,1)*(x(3,2) - 
     *  x(3,3)) + x(1,2)*(-x(3,1) + x(3,3)))
       shp(3,4)   = fact * (x(1,3)*(-x(2,1) + x(2,2)) + x(1,2)*(x(2,1) -
     *   x(2,3)) + x(1,1)*(-x(2,2) + x(2,3)))
       shp(4,4)   = sg(3)
            

      
      do i=1, ndm
        do j=1, nen
          do k=1, ndm
            do l=1, nen
           Dshp(i,j,k,l)= - shp(i,l) * shp(k,j)
         enddo
        enddo
       enddo
      enddo

      return
      end

