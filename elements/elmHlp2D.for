c
c  subroutine compxigp2D(xi,wgp,ngp)
c  subroutine compshp2D(shp,detJ,xi,x,nen)
c
c  subroutine f2d(shp,ul,F,detF,ndf,nel,isw)
c
c=============================================================================
c
c     CALCULATION OF THE GAUSS-POINT COORDINATES AND WEIGHTS FOR  
c                     2D  ELEMENTS
c
      subroutine compxigp2D(xi,wgp,ngp)
c-----------------------------------------------------------------------------
      implicit none

      integer          ngp, i

      double precision xi(2,*), wgp(*), fact, 
     *                 r0, r1d2, r1, r3,
     *                 r1d3
      
      data    r0,  r1d2,   r1,  r3  
     *    / 0.d0, .5d0, 1.d0, 3.d0 /
     
      r1d3 = r1   / r3

      if (ngp.eq.1) then
        xi(1,1) = 0.d0
        xi(2,1) = 0.d0
        wgp(1)  = 4.d0
      elseif (ngp.eq.3) then
        xi(1,1) = r1d2
        xi(2,1) = r1d2
        xi(1,2) = r0
        xi(2,2) = r1d2
        xi(1,3) = r1d2
        xi(2,3) = r0
        do i=1, 3
          wgp(i) = r1d3*r1d2
        enddo
      elseif (ngp.eq.4) then
        fact    = 1.d0 / sqrt(3.d0)
        xi(1,1) = - fact
        xi(2,1) = - fact
        xi(1,2) = + fact
        xi(2,2) = - fact
        xi(1,3) = + fact
        xi(2,3) = + fact
        xi(1,4) = - fact
        xi(2,4) = + fact
        do i=1, 4
          wgp(i) = 1.d0
        enddo
      elseif (ngp.eq.9) then
        fact    = sqrt(.6d0)
        xi(1,1) = - fact
        xi(2,1) = - fact
        xi(1,2) = r0
        xi(2,2) = - fact
        xi(1,3) = + fact
        xi(2,3) = - fact
        xi(1,4) = - fact
        xi(2,4) = r0
        xi(1,5) = r0
        xi(2,5) = r0
        xi(1,6) = + fact
        xi(2,6) = r0
        xi(1,7) = - fact
        xi(2,7) = + fact
        xi(1,8) = r0
        xi(2,8) = + fact
        xi(1,9) = + fact
        xi(2,9) = + fact
        fact    = 25.d0 / 81.d0
        wgp(1)  = fact
        wgp(3)  = fact
        wgp(7)  = fact
        wgp(9)  = fact
        fact    = 40.d0 / 81.d0
        wgp(2)  = fact
        wgp(4)  = fact
        wgp(6)  = fact
        wgp(8)  = fact
        fact    = 64.d0 / 81.d0
        wgp(5)  = fact
      else

      call prgError
     * (1,'compxigp2D(elmlib2D)','choose appropriate ngp!')
     
      endif

      return

      end

c=============================================================================
c
c     CALCULATION OF THE SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR 
c         2D TRIANGULAR  ELEMENTS WITH 3 OR 6 NODES AND FOR
c           2D QUADRILATERAL ELEMENTS WITH 4, 8 OR 9 NODES
c
      subroutine compshp2D(shp,detJ,xi,x,nen)
c-----------------------------------------------------------------------------
      implicit none

      integer          nen, i, j, k

      double precision shp(3,*), detJ, xi(*), x(2,*),
     *                 dshp(2,9), Jmtx(2,2), fact, iJmtx(2,2), xi2(2),
     *                 r0, r1d4, r1d2, r1, r2, r4

      data    r0,   r1d4, r1d2,  r1,   r2 ,   r4 
     *     / 0.d0, .25d0, .5d0, 1.d0, 2.d0, 4.d0 /
      
      xi2(1) = xi(1) * xi(1)
      xi2(2) = xi(2) * xi(2)
c
c        nen = 3                
c       	
c      
c      2
c      |\
c      | \			                 
c      |  \			             
c      |   \			            
c      |    \				     
c      |     \
c      |      \		            
c      |       \  		        
c      3--------1	         
c
c
c        nen = 6                
c      	
c      
c      2
c      |\
c      | \			                 
c      |  \			             
c      |   \			            
c      5    4				     
c      |     \
c      |      \		            
c      |       \  		        
c      3---6----1	         
c
c
c        nen = 4                nen = 8                nen = 9
c
c           ^ xi(2)                ^ xi(2)                ^ xi(2)         
c           |                      |                      |               
c      4----+----3            4----7----3            4----7----3          
c      |    |    |            |    |    |            |    |    |          
c      |    |    |  xi(1)     |    |    |  xi(1)     |    |    |          
c    --+----+----+--->      --8----+----6--->      --8----9----6---> xi(1)
c      |    |    |            |    |    |            |    |    |          
c      |    |    |            |    |    |            |    |    |          
c      1----+----2            1----5----2            1----5----2          
c           |                      |                      |               
c

      if (nen.eq.3) then

        shp(3,1)  = xi(1)
        shp(3,2)  = xi(2)
        shp(3,3)  = r1-xi(1)-xi(2)

        dshp(1,1) =   r1
        dshp(2,1) =   r0
        dshp(1,2) =   r0
        dshp(2,2) =   r1
        dshp(1,3) = - r1
        dshp(2,3) = - r1
        
      elseif (nen.eq.6) then
      
        fact=r1-xi(1)-xi(2)
        
        shp(3,1)  = xi(1) * ( r2*xi(1) - r1)
        shp(3,2)  = xi(2) * ( r2*xi(2) - r1)
        shp(3,3)  = fact  * ( r2*fact  - r1)
        shp(3,4)  = r4 * xi(1) * xi(2)
        shp(3,5)  = r4 * xi(2) * fact
        shp(3,6)  = r4 * fact * xi(1)

        dshp(1,1) =  r4 * xi(1)-r1
        dshp(2,1) =  r0
        dshp(1,2) =  r0
        dshp(2,2) =  r4 * xi(2)-r1
        dshp(1,3) =  r1 - r4 * fact 
        dshp(2,3) =  r1 - r4 * fact  
        dshp(1,4) =  r4 * xi(2)
        dshp(2,4) =  r4 * xi(1)
        dshp(1,5) = -r4 * xi(2)
        dshp(2,5) =  r4 * (fact-xi(2))
        dshp(1,6) =  r4 * (fact-xi(1))
        dshp(2,6) = -r4 * xi(1)
        
      elseif (nen.eq.4) then

        shp(3,1)  = r1d4 * (r1 - xi(1)) * (r1 - xi(2))
        shp(3,2)  = r1d4 * (r1 + xi(1)) * (r1 - xi(2))
        shp(3,3)  = r1d4 * (r1 + xi(1)) * (r1 + xi(2))
        shp(3,4)  = r1d4 * (r1 - xi(1)) * (r1 + xi(2))

        dshp(1,1) = - r1d4 * (r1 - xi(2))
        dshp(2,1) = - r1d4 * (r1 - xi(1))
        dshp(1,2) = + r1d4 * (r1 - xi(2))
        dshp(2,2) = - r1d4 * (r1 + xi(1))
        dshp(1,3) = + r1d4 * (r1 + xi(2))
        dshp(2,3) = + r1d4 * (r1 + xi(1))
        dshp(1,4) = - r1d4 * (r1 + xi(2))
        dshp(2,4) = + r1d4 * (r1 - xi(1))

      elseif (nen.eq.8) then

        shp(3,5)  = r1d2 * ((r1 - xi2(1)) * (r1 - xi(2)) )
        shp(3,6)  = r1d2 * ((r1 + xi(1))  * (r1 - xi2(2)))
        shp(3,7)  = r1d2 * ((r1 - xi2(1)) * (r1 + xi(2)) )
        shp(3,8)  = r1d2 * ((r1 - xi(1))  * (r1 - xi2(2)))

        shp(3,1)  = r1d4*((r1-xi(1))*(r1-xi(2))-r2*(shp(3,5)+shp(3,8)))
        shp(3,2)  = r1d4*((r1+xi(1))*(r1-xi(2))-r2*(shp(3,5)+shp(3,6)))
        shp(3,3)  = r1d4*((r1+xi(1))*(r1+xi(2))-r2*(shp(3,6)+shp(3,7)))
        shp(3,4)  = r1d4*((r1-xi(1))*(r1+xi(2))-r2*(shp(3,7)+shp(3,8)))

        dshp(1,5) = r1d2 * (- r2 * xi(1) * (r1 - xi(2)))
        dshp(2,5) = r1d2 * (- (r1 - xi2(1))            )
        dshp(1,6) = r1d2 * (  (r1 - xi2(2))            )
        dshp(2,6) = r1d2 * (- r2 * xi(2) * (r1 + xi(1)))
        dshp(1,7) = r1d2 * (- r2 * xi(1) * (r1 + xi(2)))
        dshp(2,7) = r1d2 * (  (r1 - xi2(1))            )
        dshp(1,8) = r1d2 * (- (r1 - xi2(2))            )
        dshp(2,8) = r1d2 * (- r2 * xi(2) * (r1 - xi(1)))

        dshp(1,1) = r1d4*(-r1+xi(2)-r2*(dshp(1,8)+dshp(1,5)))
        dshp(2,1) = r1d4*(-r1+xi(1)-r2*(dshp(2,8)+dshp(2,5)))
        dshp(1,2) = r1d4*( r1-xi(2)-r2*(dshp(1,5)+dshp(1,6)))
        dshp(2,2) = r1d4*(-r1-xi(1)-r2*(dshp(2,5)+dshp(2,6)))
        dshp(1,3) = r1d4*( r1+xi(2)-r2*(dshp(1,6)+dshp(1,7)))
        dshp(2,3) = r1d4*( r1+xi(1)-r2*(dshp(2,6)+dshp(2,7)))
        dshp(1,4) = r1d4*(-r1-xi(2)-r2*(dshp(1,7)+dshp(1,8)))
        dshp(2,4) = r1d4*( r1-xi(1)-r2*(dshp(2,7)+dshp(2,8)))

      elseif (nen.eq.9) then

        shp(3,9) = (r1 - xi2(1)) * (r1 - xi2(2))

        shp(3,5)  = r1d2 * ((r1 - xi2(1)) * (r1 - xi(2))  - shp(3,9))
        shp(3,6)  = r1d2 * ((r1 + xi(1))  * (r1 - xi2(2)) - shp(3,9))
        shp(3,7)  = r1d2 * ((r1 - xi2(1)) * (r1 + xi(2))  - shp(3,9))
        shp(3,8)  = r1d2 * ((r1 - xi(1))  * (r1 - xi2(2)) - shp(3,9))

        shp(3,1)  = r1d4*((r1-xi(1))*(r1-xi(2))-r2*(shp(3,5)+shp(3,8))
     *                                                       -shp(3,9))
        shp(3,2)  = r1d4*((r1+xi(1))*(r1-xi(2))-r2*(shp(3,5)+shp(3,6))
     *                                                       -shp(3,9))
        shp(3,3)  = r1d4*((r1+xi(1))*(r1+xi(2))-r2*(shp(3,6)+shp(3,7))
     *                                                       -shp(3,9))
        shp(3,4)  = r1d4*((r1-xi(1))*(r1+xi(2))-r2*(shp(3,7)+shp(3,8))
     *                                                       -shp(3,9))

        dshp(1,9) = - r2 * xi(1) * (r1 - xi2(2))
        dshp(2,9) = - r2 * xi(2) * (r1 - xi2(1))

        dshp(1,5) = r1d2 * (- r2 * xi(1) * (r1 - xi(2)) - dshp(1,9))
        dshp(2,5) = r1d2 * (- (r1 - xi2(1))             - dshp(2,9))
        dshp(1,6) = r1d2 * (  (r1 - xi2(2))             - dshp(1,9))
        dshp(2,6) = r1d2 * (- r2 * xi(2) * (r1 + xi(1)) - dshp(2,9))
        dshp(1,7) = r1d2 * (- r2 * xi(1) * (r1 + xi(2)) - dshp(1,9))
        dshp(2,7) = r1d2 * (  (r1 - xi2(1))             - dshp(2,9))
        dshp(1,8) = r1d2 * (- (r1 - xi2(2))             - dshp(1,9))
        dshp(2,8) = r1d2 * (- r2 * xi(2) * (r1 - xi(1)) - dshp(2,9))

        dshp(1,1) = r1d4*(-r1+xi(2)-r2*(dshp(1,8)+dshp(1,5))-dshp(1,9))
        dshp(2,1) = r1d4*(-r1+xi(1)-r2*(dshp(2,8)+dshp(2,5))-dshp(2,9))
        dshp(1,2) = r1d4*( r1-xi(2)-r2*(dshp(1,5)+dshp(1,6))-dshp(1,9))
        dshp(2,2) = r1d4*(-r1-xi(1)-r2*(dshp(2,5)+dshp(2,6))-dshp(2,9))
        dshp(1,3) = r1d4*( r1+xi(2)-r2*(dshp(1,6)+dshp(1,7))-dshp(1,9))
        dshp(2,3) = r1d4*( r1+xi(1)-r2*(dshp(2,6)+dshp(2,7))-dshp(2,9))
        dshp(1,4) = r1d4*(-r1-xi(2)-r2*(dshp(1,7)+dshp(1,8))-dshp(1,9))
        dshp(2,4) = r1d4*( r1-xi(1)-r2*(dshp(2,7)+dshp(2,8))-dshp(2,9))

      else

        call prgError
     * (1,'compshpQ2D(elmlib2D)','choose appropriate nen!')

      endif

      do i=1, 2
        do j=1, 2
          Jmtx(i,j) = r0
          do k=1, nen
            Jmtx(i,j) = Jmtx(i,j) + dshp(j,k) * x(i,k)
          enddo
        enddo
      enddo

      detJ = Jmtx(1,1) * Jmtx(2,2) - Jmtx(2,1) * Jmtx(1,2)

      fact = r1 / detJ

      iJmtx(1,1) = + fact * Jmtx(2,2)
      iJmtx(2,1) = - fact * Jmtx(2,1)
      iJmtx(1,2) = - fact * Jmtx(1,2)
      iJmtx(2,2) = + fact * Jmtx(1,1)

      do j=1, 2
        do k=1, nen
          shp(j,k) = dshp(1,k) * iJmtx(1,j) + dshp(2,k) * iJmtx(2,j)
        enddo
      enddo

      return

      end

c=============================================================================
c
c     CALCULATION OF THE DEFORMATION GRADIENT  F  FOR 2D PROBLEMS
c
      subroutine f2d(shp,ul,F,detF,ndf,nel,isw)
c-----------------------------------------------------------------------------
      implicit none

      integer          ndf, nel, isw, i, j, k
      double precision shp(3,*), ul(ndf,*), F(2,*), Finv(2,2), detF,
     *                 r0, r1

      data r0, r1 / 0.d0, 1.d0 /
c
c     isw1 = 1 ==> reference configuration, known X_A: F = 1 + Grad u
c     isw1 = 2 ==> current configuration, known x_i: F^-1 = 1 - grad u
c
      goto(1,2), isw
c
1     continue
      do i=1, 2
        do j=1, 2
          F(i,j) = r0
          do k=1, nel
            F(i,j) = F(i,j) + ul(i,k)*shp(j,k)
          enddo
        enddo
        F(i,i) = F(i,i) + r1
      enddo
      detF = F(1,1)*F(2,2) - F(1,2)*F(2,1)
      return

2     continue

c.... calculation of F^-1
      do i=1, 2
        do j=1, 2
          Finv(i,j) = r0
          do k=1, nel
            Finv(i,j) = Finv(i,j) - ul(i,k)*shp(j,k)
          enddo
        enddo
        Finv(i,i) = Finv(i,i) + r1
      enddo

c.... compute det(F)
      detF = r1 / (Finv(1,1)*Finv(2,2) - Finv(1,2)*Finv(2,1))

c.... calculation of F
      F(2,2) =   Finv(1,1) * detF
      F(1,1) =   Finv(2,2) * detF
      F(1,2) = - Finv(1,2) * detF
      F(2,1) = - Finv(2,1) * detF

      return
      end

