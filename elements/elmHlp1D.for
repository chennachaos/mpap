
c
c  subroutine comp_xigp1D(xi,wgp,ngp)
c  subroutine comp_shp1D(shp,detJ,xi,x,nen)
c
c
c=============================================================================
c
c   CALCULATION OF THE GAUSS-POINT COORDINATES AND WEIGHTS FOR 1D ELEMENTS
c
      subroutine comp_xigp1D(xi,wgp,ngp)
c-----------------------------------------------------------------------------
      implicit none

      integer          ngp, i

      double precision xi(*), wgp(*), fact, r0

      r0 = 0.d0

      if (ngp.eq.1) then
        xi(1)  = 0.d0
        wgp(1) = 2.d0
      elseif (ngp.eq.2) then
        fact   = 1.d0 / sqrt(3.d0)
        xi(1)  = - fact
        xi(2)  = + fact
        wgp(1) = 1.d0
        wgp(2) = 1.d0
      elseif (ngp.eq.3) then
        fact    = sqrt(.6d0)
        xi(1)  = - fact
        xi(2)  = r0
        xi(3)  = + fact
        fact   = 5.d0 / 9.d0
        wgp(1) = fact
        wgp(2) = fact * 1.6d0
        wgp(3) = fact
      else

        call prgError(1,'comp_xigp1D',
     *                  'can''t handle this number of Gauss points !')
      endif

      return

      end

c=============================================================================
c
c     CALCULATION OF THE SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR 
c                    1D ELEMENTS WITH 2 OR 3 NODES
c
      subroutine comp_shp1D(shp,detJ,xi,x,nen)
c-----------------------------------------------------------------------------
      implicit none

      integer          nen, k

      double precision shp(2,*), detJ, xi, x(*),
     *                 dshp(9), fact, xi2,
     *                 r1d2, r1, r2

      data   r1d2,  r1,   r2  
     *     / .5d0, 1.d0, 2.d0 /
      
      xi2 = xi * xi

      if (nen.eq.2) then

c       1-------2         

        shp(2,1)  = r1d2 * (r1 - xi)
        shp(2,2)  = r1d2 * (r1 + xi)

        dshp(1) = - r1d2
        dshp(2) = + r1d2

      elseif (nen.eq.3) then

c       1---3---2         

        shp(2,1)  = r1d2 * (xi2 - xi)
        shp(2,2)  = r1d2 * (xi2 + xi)
        shp(2,3)  = r1 - xi2

        dshp(1) = xi - r1d2
        dshp(2) = xi + r1d2
        dshp(3) = - r2 * xi

      else

        call prgError(1,'comp_shp1D(helpsld)','choose nen = 2 or 3 !')

      endif

      detJ = 0.d0
      do k=1, nen
        detJ = detJ + dshp(k) * x(k)
      enddo

      fact = r1 / detJ

      do k=1, nen
        shp(1,k) = dshp(k) * fact
      enddo

      return

      end

