
c  subroutine compxigp3D(xi,wgp,ngp)
c  subroutine compshp3D(shp,detJ,xi,x,nen)
c  subroutine f3D(shp,ul,F,detF,ndf,nel,isw)
c
c=============================================================================
c
c     CALCULATION OF THE GAUSS-POINT COORDINATES AND WEIGHTS FOR  
c                     3D QUADRILATERAL ELEMENTS
c
      subroutine compxigp3D(xi,wgp,ngp)
c-----------------------------------------------------------------------------
      implicit none

      integer          ngp, i, j

      double precision xi(3,*), wgp(*), fact, 
     *                 fact1, xp(3,27)

      data  (xp(i,1),  i=1,3) / -1.d0, -1.d0, -1.d0 /
      data  (xp(i,2),  i=1,3) /  1.d0, -1.d0, -1.d0 /
      data  (xp(i,3),  i=1,3) /  1.d0,  1.d0, -1.d0 /
      data  (xp(i,4),  i=1,3) / -1.d0,  1.d0, -1.d0 /
      data  (xp(i,5),  i=1,3) / -1.d0, -1.d0,  1.d0 /
      data  (xp(i,6),  i=1,3) /  1.d0, -1.d0,  1.d0 /
      data  (xp(i,7),  i=1,3) /  1.d0,  1.d0,  1.d0 /
      data  (xp(i,8),  i=1,3) / -1.d0,  1.d0,  1.d0 /

      if (ngp.eq.1) then

        xi(1,1) = 0.d0
        xi(2,1) = 0.d0
        xi(3,1) = 0.d0
        wgp(1)  = 8.d0
        
      else if (ngp.eq.4) then

        fact=0.5854102
        fact1=0.1381966
      
        xi(1,1) = fact
        xi(2,1) = fact1
        xi(3,1) = fact1
        
        xi(1,2) = fact1
        xi(2,2) = fact
        xi(3,2) = fact1
        
        xi(1,3) = fact1
        xi(2,3) = fact1
        xi(3,3) = fact
        
        xi(1,4) = fact1
        xi(2,4) = fact1
        xi(3,4) = fact1
        
        do i=1, 4
          wgp(i) = .25d0/6.d0
        enddo
          
      else if (ngp.eq.8) then
          
        fact    = 1.d0 / sqrt(3.d0)
        do i=1, ngp
          do j=1,3
            xi(j,i)=xp(j,i)*fact
          enddo
          wgp(i) = 1.d0
        enddo
        
      else if (ngp.eq.27) then
       
        data  (xp(i,9),  i=1,3) / 0.d0, -1.d0, -1.d0 /
        data  (xp(i,11), i=1,3) / 0.d0,  1.d0, -1.d0 /
        data  (xp(i,13), i=1,3) / 0.d0, -1.d0,  1.d0 /
        data  (xp(i,15), i=1,3) / 0.d0,  1.d0,  1.d0 /
        
        data  (xp(i,10), i=1,3) /  1.d0, 0.d0, -1.d0 /
        data  (xp(i,12), i=1,3) / -1.d0, 0.d0, -1.d0 /
        data  (xp(i,14), i=1,3) /  1.d0, 0.d0,  1.d0 /
        data  (xp(i,16), i=1,3) / -1.d0, 0.d0,  1.d0 /
        
        data  (xp(i,17), i=1,3) / -1.d0, -1.d0, 0.d0 /
        data  (xp(i,18), i=1,3) /  1.d0, -1.d0, 0.d0 /
        data  (xp(i,19), i=1,3) /  1.d0,  1.d0, 0.d0 /
        data  (xp(i,20), i=1,3) / -1.d0,  1.d0, 0.d0 /
          
        data  (xp(i,21), i=1,3) /  1.d0,  0.d0,  0.d0 /
        data  (xp(i,22), i=1,3) / -1.d0,  0.d0,  0.d0 /
        data  (xp(i,23), i=1,3) /  0.d0,  1.d0,  0.d0 /
        data  (xp(i,24), i=1,3) /  0.d0, -1.d0,  0.d0 /
        data  (xp(i,25), i=1,3) /  0.d0,  0.d0,  1.d0 /
        data  (xp(i,26), i=1,3) /  0.d0,  0.d0, -1.d0 /
        data  (xp(i,27), i=1,3) /  0.d0,  0.d0,  0.d0 /

        fact    = sqrt(6.d-1)
        fact1=125.d0/729.d0
        do i=1, 8
          wgp(i) = fact1
          do j=1,3
            xi(j,i)=xp(j,i)*fact
          enddo
        enddo
          
        fact1=200.d0/729.d0
        do i=9, 20
          wgp(i) = fact1
          do j=1,3
            xi(j,i)=xp(j,i)*fact
          enddo
        enddo
          
        fact1=320.d0/729.d0
        do i=21, 26
          wgp(i) = fact1
          do j=1,3
            xi(j,i)=xp(j,i)*fact
          enddo
        enddo
          
        wgp(27) = 512.d0/729.d0
        do j=1,3
          xi(j,27)=xp(j,27)*fact
        enddo

      else

        call prgError
     *     (1,'compxigp3D(helpsld3D)','choose appropriate ngp!')
               
      endif

      return

      end

c=============================================================================
c
c     CALCULATION OF THE SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR 
c                            3D ELEMENTS 
c
      subroutine compshp3D(shp,detJ,xi,x,nen)
c-----------------------------------------------------------------------------
      implicit none

      integer          nen, i, j, k

      double precision shp(4,*), detJ, xi(*), x(3,*),
     *                 dshp(3,20), Jmtx(3,3), fact, iJmtx(3,3), 
     *                 xp(3,20),det_u,
     *                 r0, r1d8, r1d4,r1d2, r1, r2, r4

      data    r0 ,   r1d8,  r1d4,  r1d2,   r1,   r2,   r4  
     *     / 0.d0, .125d0, .25d0,  .5d0, 1.d0, 2.d0 , 4.d0/
     
!         nen = 4     
! 
!              1            
!              /\ .                  
!             /  \   .            
!            /    \     .     
!           /      \       .4        
!          /        \     .     
!         /          \  .
!        3------------2               
!         nen = 10 :
!        1--5--2
!        2--6--3 
!        3--7--1
!        1--8--4
!        3--10--4
!        2--9--4

      if (nen.eq.4) then

        shp(4,1)  = r1-xi(1)-xi(2)-xi(3)
        shp(4,2)  = xi(1)
        shp(4,3)  = xi(2)
        shp(4,4)  = xi(3)

        dshp(1,1) = - r1
        dshp(2,1) = - r1
        dshp(3,1) = - r1
        dshp(1,2) =   r1
        dshp(2,2) =   r0
        dshp(3,2) =   r0
        dshp(1,3) =   r0
        dshp(2,3) =   r1
        dshp(3,3) =   r0
        dshp(1,4) =   r0
        dshp(2,4) =   r0
        dshp(3,4) =   r1

      else if (nen.eq.10) then
      
        fact  = r1-xi(1)-xi(2)-xi(3)

        shp(4,1)   = fact  * ( r2*fact  - r1)     
        shp(4,2)   = xi(1) * ( r2*xi(1) - r1)
        shp(4,3)   = xi(2) * ( r2*xi(2) - r1)
        shp(4,4)   = xi(3) * ( r2*xi(3) - r1)
        shp(4,5)   = r4 * xi(1) * fact
        shp(4,6)   = r4 * xi(1) * xi(2)
        shp(4,7)   = r4 * xi(2) * fact
        shp(4,8)   = r4 * xi(3) * fact
        shp(4,9)   = r4 * xi(3) * xi(1)
        shp(4,10)  = r4 * xi(2) * xi(3)
        
        dshp(1,1) =  r1 - r4 * fact 
        dshp(2,1) =  r1 - r4 * fact
        dshp(3,1) =  r1 - r4 * fact
        dshp(1,2) =  r4 * xi(1)-r1
        dshp(2,2) =  r0
        dshp(3,2) =  r0
        dshp(1,3) =  r0
        dshp(2,3) =  r4 * xi(2)-r1
        dshp(3,3) =  r0
        dshp(1,4) =  r0
        dshp(2,4) =  r0
        dshp(3,4) =  r4 * xi(2)-r1

        dshp(1,5) =  r4 * (fact-xi(1))
        dshp(2,5) = -r4 * xi(1)
        dshp(3,5) = -r4 * xi(1)
        dshp(1,6) =  r4 * xi(2)
        dshp(2,6) =  r4 * xi(1) 
        dshp(3,6) =  r0       
        dshp(1,7) = -r4 * xi(2)
        dshp(2,7) =  r4 * (fact-xi(2))
        dshp(3,7) = -r4 * xi(2)
        dshp(1,8) = -r4 * xi(3) 
        dshp(2,8) = -r4 * xi(3)
        dshp(3,8) =  r4 * (fact-xi(3))
        dshp(1,9) =  r4 * xi(3)
        dshp(2,9) =  r0 
        dshp(3,9) =  r4 * xi(1)  
        dshp(1,10) =  r0
        dshp(2,10) =  r4 * xi(3) 
        dshp(3,10) =  r4 * xi(2)    
c       
c        nen = 8      
c
c                         
c          8----+----7           
c         .|        .|
c       .  |^ xi(3). |
c      5---|+----6   |         
c      |   ||    |   |         
c      |   ||    |  xi(1)     
c    --+---|+----+--->     
c      |   4-----|---3           
c      |  . |    | .         
c      |.   |    |.           
c      1----+----2             
c           |                    
c

      elseif ((nen.eq.8).or.(nen.eq.20)) then
      
        data  (xp(i,1),  i=1,3) / -1.d0, -1.d0, -1.d0 /
        data  (xp(i,2),  i=1,3) /  1.d0, -1.d0, -1.d0 /
        data  (xp(i,3),  i=1,3) /  1.d0,  1.d0, -1.d0 /
        data  (xp(i,4),  i=1,3) / -1.d0,  1.d0, -1.d0 /
        data  (xp(i,5),  i=1,3) / -1.d0, -1.d0,  1.d0 /
        data  (xp(i,6),  i=1,3) /  1.d0, -1.d0,  1.d0 /
        data  (xp(i,7),  i=1,3) /  1.d0,  1.d0,  1.d0 /
        data  (xp(i,8),  i=1,3) / -1.d0,  1.d0,  1.d0 /
        
        do i=1, 8
        
          shp(4,i)  = r1d8 * (r1 + xp(1,i)*xi(1)) * 
     *                (r1 + xp(2,i)*xi(2)) * (r1 + xp(3,i)*xi(3))
        
          dshp(1,i)  = r1d8 * xp(1,i)* 
     *                (r1 + xp(2,i)*xi(2)) * (r1 + xp(3,i)*xi(3))
        
          dshp(2,i)  = r1d8 * (r1 + xp(1,i)*xi(1)) * 
     *                xp(2,i) * (r1 + xp(3,i)*xi(3))
        
          dshp(3,i)  = r1d8 * (r1 + xp(1,i)*xi(1)) * 
     *                (r1 + xp(2,i)*xi(2)) * xp(3,i)

        enddo
   
   
c          nen = 20      
c
c               ^ xi(3)
c              
c            8---15----7           
c           .|   .    .|
c          16|  .    14|
c         .  | .    .  |
c        5---|13---6   |         
c        |  20|....|..19         
c        |   ||    |   | xi(1)     
c      --+---|+----+---|-->     
c       17   4---11|---3           
c        |  . |  .18  .     
c        | 12.|....|.10    
c        |.   |.   |.           
c        1----9----2             
c             |                    
c
         
   
        if (nen.eq.20) then
 
          data  (xp(i,9),  i=1,3) / 0.d0, -1.d0, -1.d0 /
          data  (xp(i,11), i=1,3) / 0.d0,  1.d0, -1.d0 /
          data  (xp(i,13), i=1,3) / 0.d0, -1.d0,  1.d0 /
          data  (xp(i,15), i=1,3) / 0.d0,  1.d0,  1.d0 /
          
          data  (xp(i,10), i=1,3) /  1.d0, 0.d0, -1.d0 /
          data  (xp(i,12), i=1,3) / -1.d0, 0.d0, -1.d0 /
          data  (xp(i,14), i=1,3) /  1.d0, 0.d0,  1.d0 /
          data  (xp(i,16), i=1,3) / -1.d0, 0.d0,  1.d0 /
          
          data  (xp(i,17), i=1,3) / -1.d0, -1.d0, 0.d0 /
          data  (xp(i,18), i=1,3) /  1.d0, -1.d0, 0.d0 /
          data  (xp(i,19), i=1,3) /  1.d0,  1.d0, 0.d0 /
          data  (xp(i,20), i=1,3) / -1.d0,  1.d0, 0.d0 /
           
          do i=1, 8
          
            do j=1, 3
              dshp(j,i)  = dshp(j,i) * ( xp(1,i)*xi(1)
     *                      + xp(2,i)*xi(2) + xp(3,i)*xi(3)- r2)
     *                      + shp(4,i) * xp(j,i)
            enddo
              
            shp(4,i)  = shp(4,i) * ( xp(1,i)*xi(1)
     *                      + xp(2,i)*xi(2) + xp(3,i)*xi(3)- r2) 
   
     
            if(i.lt.5) then

              j = 7 + 2*i
   
              shp(4,j)  = r1d4 * (r1 - xi(1)*xi(1)) *
     *                   (r1 + xp(2,j)*xi(2))*(r1 + xp(3,j)*xi(3))
              dshp(1,j)  = -r1d2 * xi(1) *
     *                    (r1 + xp(2,j)*xi(2))*(r1 + xp(3,j)*xi(3))
              dshp(2,j)  = r1d4 * (r1 - xi(1)*xi(1)) *
     *                    xp(2,j)*(r1 + xp(3,j)*xi(3))
              dshp(3,j)  = r1d4 * (r1 - xi(1)*xi(1)) *
     *                    (r1 + xp(2,j)*xi(2))*xp(3,j)
   
              shp(4,(j+1)) = r1d4 * (r1 - xi(2)*xi(2)) * 
     *                 (r1 + xp(1,(j+1))*xi(1))*(r1 + xp(3,(j+1))*xi(3))
              dshp(1,(j+1)) = r1d4 * (r1 - xi(2)*xi(2)) * 
     *                  xp(1,(j+1))*(r1 + xp(3,(j+1))*xi(3))
              dshp(2,(j+1)) = -r1d2 * xi(2) * 
     *               (r1 + xp(1,(j+1))*xi(1))*(r1 + xp(3,(j+1))*xi(3))
              dshp(3,(j+1)) = r1d4 * (r1 - xi(2)*xi(2)) * 
     *                       (r1 + xp(1,(j+1))*xi(1))*xp(3,(j+1))
   
              shp(4,(i+16)) = r1d4 * (r1 - xi(3)*xi(3)) * 
     *               (r1 + xp(1,(i+16))*xi(1))*(r1 + xp(2,(i+16))*xi(2))
              dshp(1,(i+16)) = r1d4 * (r1 - xi(3)*xi(3)) * 
     *                        xp(1,(i+16))*(r1 + xp(2,(i+16))*xi(2))
              dshp(2,(i+16)) = r1d4 * (r1 - xi(3)*xi(3)) * 
     *                        (r1 + xp(1,(i+16))*xi(1))*xp(2,(i+16))
              dshp(3,(i+16)) = -r1d2 * xi(3) * 
     *               (r1 + xp(1,(i+16))*xi(1))*(r1 + xp(2,(i+16))*xi(2))
   
            endif
   
          enddo  
 
        endif !if 20

      else   !if 8 or 20

        call prgError
     *         (1,'compshp3D(helpsld3D)','choose appropriate nen!')

      endif

      do i=1, 3
        do j=1, 3
          Jmtx(i,j) = r0
          do k=1, nen
            Jmtx(i,j) = Jmtx(i,j) + dshp(j,k) * x(i,k)
          enddo
        enddo
      enddo

      detJ = det_u(Jmtx)

      iJmtx(1,1)=(Jmtx(2,2)*Jmtx(3,3)-Jmtx(3,2)*Jmtx(2,3))/detJ
      iJmtx(1,2)=(Jmtx(1,3)*Jmtx(3,2)-Jmtx(3,3)*Jmtx(1,2))/detJ
      iJmtx(1,3)=(Jmtx(1,2)*Jmtx(2,3)-Jmtx(2,2)*Jmtx(1,3))/detJ
      iJmtx(2,1)=(Jmtx(2,3)*Jmtx(3,1)-Jmtx(3,3)*Jmtx(2,1))/detJ
      iJmtx(2,2)=(Jmtx(1,1)*Jmtx(3,3)-Jmtx(3,1)*Jmtx(1,3))/detJ
      iJmtx(2,3)=(Jmtx(1,3)*Jmtx(2,1)-Jmtx(2,3)*Jmtx(1,1))/detJ
      iJmtx(3,1)=(Jmtx(2,1)*Jmtx(3,2)-Jmtx(3,1)*Jmtx(2,2))/detJ
      iJmtx(3,2)=(Jmtx(1,2)*Jmtx(3,1)-Jmtx(3,2)*Jmtx(1,1))/detJ
      iJmtx(3,3)=(Jmtx(1,1)*Jmtx(2,2)-Jmtx(2,1)*Jmtx(1,2))/detJ 
     


      do k=1, nen
          shp(1,k) = dshp(1,k) * iJmtx(1,1) + 
     *               dshp(2,k) * iJmtx(2,1) +
     *               dshp(3,k) * iJmtx(3,1)
          shp(2,k) = dshp(1,k) * iJmtx(1,2) + 
     *               dshp(2,k) * iJmtx(2,2) +
     *               dshp(3,k) * iJmtx(3,2)
          shp(3,k) = dshp(1,k) * iJmtx(1,3) + 
     *               dshp(2,k) * iJmtx(2,3) +
     *               dshp(3,k) * iJmtx(3,3)
      enddo


      return

      end

c=============================================================================
c
c     CALCULATION OF THE DEFORMATION GRADIENT  F  FOR 3D PROBLEMS
c
      subroutine f3D(shp,ul,F,detF,ndf,nel,isw)
c-----------------------------------------------------------------------------
      implicit none

      integer          ndf, nel, isw, i, j, k
      double precision shp(4,*), ul(ndf,*), F(3,*), Finv(3,3), detF,
     *                 r0, r1, det_u

      integer P(3)

      data r0, r1 / 0.d0, 1.d0 /
c
c     isw1 = 1 ==> reference configuration, known X_A: F = 1 + Grad u
c     isw1 = 2 ==> current configuration, known x_i: F^-1 = 1 - grad u
c
      goto(1,2), isw
c
1     continue
      do i=1, 3
        do j=1, 3
          F(i,j) = r0
          do k=1, nel
            F(i,j) = F(i,j) + ul(i,k)*shp(j,k)
          enddo
        enddo
        F(i,i) = F(i,i) + r1
      enddo
      detF = det_u(F)
      return

2     continue

c.... calculation of F^-1
      do i=1, 3
        do j=1, 3
          Finv(i,j) = r0
          do k=1, nel
            Finv(i,j) = Finv(i,j) - ul(i,k)*shp(j,k)
          enddo
        enddo
        Finv(i,i) = Finv(i,i) + r1
      enddo


c.... calculation of F

      call decompLR_matrix( Finv,P,3,0)
      call inverse_matrix( Finv,F,P,3)

c.... compute det(F)

      detF = det_u(F)

      return
      end

