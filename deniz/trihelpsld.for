c  subroutine compxigpT(xi,wgp,ngp)
c  subroutine  areaTrianglefor(area, x1,x2,x3)
c  subroutine compshpT(shp,detJ,xi,x,nen)
c  subroutine compshpbQ(shp,Jmtx,xi,x,nen)
c  subroutine compshpbT(shp,Jmtx,xi,x,nen)




c=============================================================================
c
c     CALCULATION OF THE GAUSS-POINT COORDINATES AND WEIGHTS FOR  
c                     2D QUADRILATERAL ELEMENTS
c
      subroutine compxigpT(xi,wgp,ngp)
c-----------------------------------------------------------------------------
      implicit none

      integer          ngp, i

      double precision xi(2,*), wgp(*), fact,
     *                 r0, r1d2, r1, r2, r3,
     *                 r1d3, r1d6, r2d3

      data    r0,  r1d2,   r1,   r2,  r3  
     *    / 0.d0, .5d0, 1.d0, 2.d0, 3.d0 /
      

      r1d3 = r1   / r3
      r1d6 = r1d2 * r1d3
      r2d3 = r2   / r3


      if (ngp.eq.1) then
        xi(1,1) = r1d3
        xi(2,1) = r1d3
        wgp(1)  = r1
      elseif (ngp.eq.3) then
        xi(1,1) = r2d3
        xi(2,1) = r1d6
        xi(1,2) = r1d6
        xi(2,2) = r2d3
        xi(1,3) = r1d6
        xi(2,3) = r1d6
        do i=1, 3
          wgp(i) = r1d3
        enddo
      elseif (ngp.eq.4) then
        xi(1,1) = r1d3
        xi(2,1) = r1d3
        xi(1,2) = .6d0
        xi(2,2) = .2d0
        xi(1,3) = .2d0
        xi(2,3) = .2d0
        xi(1,4) = .2d0
        xi(2,4) = .6d0

        wgp(1) = -0.5625d0
        fact = r1 / 1.92d0
        do i=2, 4
          wgp(i) = fact
        enddo

      else

        call prgError(1,'compxigpT',
     *                  'can''t handle this number of Gauss points !')
      endif


      return

      end

c=============================================================================
c
c     CALCULATION OF AREA OF A TRIANGLE
c
      subroutine  areaTrianglefor(area, x1,x2,x3)
c------------------------------------------------------------------------

      double precision area, x1(2), x2(2) ,x3(2)
      
      integer k
      
       area=0.5d0 * (x1(1)*(x2(2)-x3(2))
     *		 +x2(1)*(x3(2)-x1(2))
     *		 +x3(1)*(x1(2)-x2(2)))
      
      return
      
      end
      
c=============================================================================
c
c     CALCULATION OF THE SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR 
c           2D QUADRILATERAL ELEMENTS WITH 4, 8 OR 9 NODES
c
      subroutine compshpT(shp,detJ,xi,x,nen)
c-----------------------------------------------------------------------------
      implicit none

      integer          nen, i, j, k

      double precision shp(3,*), detJ, xi(*), x(2,*),tarea(6),xc(2),
     *                 dshp(2,6), Jmtx(2,2), fact, iJmtx(2,2), xi2(2),
     *                 r0, r1d4, r1d2, r1, r2, r3,r4

      data    r0,   r1d4, r1d2,  r1,   r2,    r3,    r4 
     *     / 0.d0, .25d0, .5d0, 1.d0, 2.d0 , 3.d0 , 4.d0/
     
      
c
c        nen = 3                
c       	
c      
c      3
c      |\
c	 | \			                 
c      |  \			             
c      |   \			            
c      |    \				     
c      |     \
c	   |      \		            
c      |       \  		        
c      1--------2	         
c
c
c        nen = 6                
c      	
c      
c      3
c      |\
c	 | \			                 
c      |  \			             
c      |   \			            
c      6    5				     
c      |     \
c	 |      \		            
c      |       \  		        
c      1---4----2	         
c

      if (nen.eq.3) then

        shp(3,1)  = r1-xi(1)-xi(2)
        shp(3,2)  = xi(1)
        shp(3,3)  = xi(2)


        dshp(1,1) = - r1
        dshp(2,1) = - r1
        dshp(1,2) = + r1
        dshp(2,2) =   r0
        dshp(1,3) =   r0
        dshp(2,3) = + r1

        
      else if (nen.eq.6) then

        shp(3,1)  = r2 * (xi(1)+xi(2)-r1) * (xi(1)+xi(2)-r1d2)
        shp(3,2)  = xi(1) * (r2*xi(1)-r1)
        shp(3,3)  = xi(2) * (r2*xi(2)-r1)
        shp(3,4)  = -r4 * xi(1) * (xi(1)+xi(2)-r1)
        shp(3,5)  = +r4 * xi(1) * xi(2)
        shp(3,6)  = -r4 * xi(2) * (xi(1)+xi(2)-r1)

        dshp(1,1) = +r4 * (xi(1)+xi(2))-r3
        dshp(2,1) = +r4 * (xi(1)+xi(2))-r3
        dshp(1,2) = +r4 * xi(1)-r1
        dshp(2,2) =  r0
        dshp(1,3) =  r0
        dshp(2,3) = +r4 * xi(2)-r1
        dshp(1,4) = -r4 * (r2*xi(1)+xi(2)-r1)
        dshp(2,4) = -r4 * xi(1)
        dshp(1,5) = +r4 * xi(2)
        dshp(2,5) = +r4 * xi(1)
        dshp(1,6) = -r4 * xi(2)
        dshp(2,6) = -r4 * (r2*xi(2)+xi(1)-r1)


      else

       call prgError(1,'compshpT','choose nen = 3 or 6 !')

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
c     CALCULATION OF THE SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR 
c           2D QUADRILATERAL ELEMENTS WITH 4, 8 OR 9 NODES
c
      subroutine compshpbQ(shp,Jmtx,xi,x,nen)
c-----------------------------------------------------------------------------
      implicit none

      integer          nen, i, j, k

      double precision shp(*), detJ, xi(*), x(2,*),
     *                 dshp(2,9), Jmtx(2,2), fact, iJmtx(2,2), xi2(2),
     *                 r0, r1d4, r1d2, r1, r2

      data    r0,   r1d4, r1d2,  r1,   r2  
     *     / 0.d0, .25d0, .5d0, 1.d0, 2.d0 /
      
      xi2(1) = xi(1) * xi(1)
      xi2(2) = xi(2) * xi(2)

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

      if (nen.eq.4) then

        shp(1)  = r1d4 * (r1 - xi(1)) * (r1 - xi(2))
        shp(2)  = r1d4 * (r1 + xi(1)) * (r1 - xi(2))
        shp(3)  = r1d4 * (r1 + xi(1)) * (r1 + xi(2))
        shp(4)  = r1d4 * (r1 - xi(1)) * (r1 + xi(2))

        dshp(1,1) = - r1d4 * (r1 - xi(2))
        dshp(2,1) = - r1d4 * (r1 - xi(1))
        dshp(1,2) = + r1d4 * (r1 - xi(2))
        dshp(2,2) = - r1d4 * (r1 + xi(1))
        dshp(1,3) = + r1d4 * (r1 + xi(2))
        dshp(2,3) = + r1d4 * (r1 + xi(1))
        dshp(1,4) = - r1d4 * (r1 + xi(2))
        dshp(2,4) = + r1d4 * (r1 - xi(1))

      elseif (nen.eq.8) then

        shp(5)  = r1d2 * ((r1 - xi2(1)) * (r1 - xi(2)) )
        shp(6)  = r1d2 * ((r1 + xi(1))  * (r1 - xi2(2)))
        shp(7)  = r1d2 * ((r1 - xi2(1)) * (r1 + xi(2)) )
        shp(8)  = r1d2 * ((r1 - xi(1))  * (r1 - xi2(2)))

        shp(1)  = r1d4*((r1-xi(1))*(r1-xi(2))-r2*(shp(5)+shp(8)))
        shp(2)  = r1d4*((r1+xi(1))*(r1-xi(2))-r2*(shp(5)+shp(6)))
        shp(3)  = r1d4*((r1+xi(1))*(r1+xi(2))-r2*(shp(6)+shp(7)))
        shp(4)  = r1d4*((r1-xi(1))*(r1+xi(2))-r2*(shp(7)+shp(8)))

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

        shp(9) = (r1 - xi2(1)) * (r1 - xi2(2))

        shp(5)  = r1d2 * ((r1 - xi2(1)) * (r1 - xi(2))  - shp(9))
        shp(6)  = r1d2 * ((r1 + xi(1))  * (r1 - xi2(2)) - shp(9))
        shp(7)  = r1d2 * ((r1 - xi2(1)) * (r1 + xi(2))  - shp(9))
        shp(8)  = r1d2 * ((r1 - xi(1))  * (r1 - xi2(2)) - shp(9))

        shp(1)  = r1d4*((r1-xi(1))*(r1-xi(2))-r2*(shp(5)+shp(8))
     *                                                       -shp(9))
        shp(2)  = r1d4*((r1+xi(1))*(r1-xi(2))-r2*(shp(5)+shp(6))
     *                                                       -shp(9))
        shp(3)  = r1d4*((r1+xi(1))*(r1+xi(2))-r2*(shp(6)+shp(7))
     *                                                       -shp(9))
        shp(4)  = r1d4*((r1-xi(1))*(r1+xi(2))-r2*(shp(7)+shp(8))
     *                                                       -shp(9))

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

        call prgError(1,'compshpbQ)','choose nen = 4, 8 or 9 !')

      endif

      do i=1, 2
        do j=1, 2
          Jmtx(i,j) = r0
          do k=1, nen
            Jmtx(i,j) = Jmtx(i,j) + dshp(j,k) * x(i,k)
          enddo
        enddo
      enddo

      return

      end 
      

c=============================================================================
c
c     CALCULATION OF THE SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR 
c           2D QUADRILATERAL ELEMENTS WITH 4, 8 OR 9 NODES
c
      subroutine compshpbT(shp,Jmtx,xi,x,nen)
c-----------------------------------------------------------------------------
      implicit none

      integer          nen, i, j, k
       
      double precision shp(*), detJ, xi(*), x(2,*),tarea(6),xc(2),
     *                 dshp(2,6), Jmtx(2,2), fact, iJmtx(2,2), xi2(2),
     *                 r0, r1d4, r1d2, r1, r2, r3,r4

      data    r0,   r1d4, r1d2,  r1,   r2,    r3,    r4 
     *     / 0.d0, .25d0, .5d0, 1.d0, 2.d0 , 3.d0 , 4.d0/
     
      
c
c        nen = 3                
c       	
c      
c      3
c      |\
c	 | \			                 
c      |  \			             
c      |   \			            
c      |    \				     
c      |     \
c	   |      \		            
c      |       \  		        
c      1--------2	         
c
c
c        nen = 6                
c      	
c      
c      3
c      |\
c	 | \			                 
c      |  \			             
c      |   \			            
c      6    5				     
c      |     \
c	 |      \		            
c      |       \  		        
c      1---4----2	         
c

      if (nen.eq.3) then

        shp(1)  = r1-xi(1)-xi(2)
        shp(2)  = xi(1)
        shp(3)  = xi(2)


        dshp(1,1) = - r1
        dshp(2,1) = - r1
        dshp(1,2) = + r1
        dshp(2,2) =   r0
        dshp(1,3) =   r0
        dshp(2,3) = + r1

        
      else if (nen.eq.6) then

        shp(1)  = r2 * (xi(1)+xi(2)-r1) * (xi(1)+xi(2)-r1d2)
        shp(2)  = xi(1) * (r2*xi(1)-r1)
        shp(3)  = xi(2) * (r2*xi(2)-r1)
        shp(4)  = -r4 * xi(1) * (xi(1)+xi(2)-r1)
        shp(5)  = +r4 * xi(1) * xi(2)
        shp(6)  = -r4 * xi(2) * (xi(1)+xi(2)-r1)

        dshp(1,1) = r4 * xi(2) + xi(1)
        dshp(2,1) = r4 * xi(1) + xi(2)
        dshp(1,2) = r4 * xi(1)-r1
        dshp(2,2) =  r0
        dshp(1,3) =  r0
        dshp(2,3) = r4 * xi(2)-r1
        dshp(1,4) = -r4 * (r2*xi(1)+xi(2)-r1)
        dshp(2,4) = -r4 * xi(1)
        dshp(1,5) = +r4 * xi(2)
        dshp(2,5) = +r4 * xi(1)
        dshp(1,6) = -r4 * xi(2)
        dshp(2,6) = -r4 * (r2*xi(2)+xi(1)-r1)


      else

       call prgError(1,'compshpbT','choose nen = 3 or 6 !')

      endif

      do i=1, 2
        do j=1, 2
          Jmtx(i,j) = r0
          do k=1, nen
            Jmtx(i,j) = Jmtx(i,j) + dshp(j,k) * x(i,k)
          enddo
        enddo
      enddo

     
      return

      end