
      integer function ale2D4nodedQuadrilateralAspectRatio
     *                    (d,xl,s,p,ndm,nst,ndat)
      implicit none

      integer          ndm, nst, ndat,
     *                 i, j
      double precision d(*), xl(ndm,*), s(nst,*), p(*),
     *                 xt(2,6), st(6,6), pt(6)

c
c     LINEAR QUADRIALTERALS,
c

      ale2D4nodedQuadrilateralAspectRatio = 0

      call pzero(p,8)
      call pzero(s,64)

c  triangle 1 - 2 - 3

      xt(1,1) = xl(1,1)
      xt(2,1) = xl(2,1)
      xt(1,2) = xl(1,2)
      xt(2,2) = xl(2,2)
      xt(1,3) = xl(1,3)
      xt(2,3) = xl(2,3)

      xt(1,4) = xl(1,5)
      xt(2,4) = xl(2,5)
      xt(1,5) = xl(1,6)
      xt(2,5) = xl(2,6)
      xt(1,6) = xl(1,7)
      xt(2,6) = xl(2,7)

      call pzero(pt,6)
      call pzero(st,36)

      call ale2D3nodedTriangleAspectRatio(d,xt,st,pt,2,6,ndat)

      do i=1, 6
        p(i) = pt(i)
        do j=1, 6
          s(i,j) = st(i,j)
        enddo
      enddo

c  triangle 2 - 3 - 4

      xt(1,1) = xl(1,2)
      xt(2,1) = xl(2,2)
      xt(1,2) = xl(1,3)
      xt(2,2) = xl(2,3)
      xt(1,3) = xl(1,4)
      xt(2,3) = xl(2,4)

      xt(1,4) = xl(1,6)
      xt(2,4) = xl(2,6)
      xt(1,5) = xl(1,7)
      xt(2,5) = xl(2,7)
      xt(1,6) = xl(1,8)
      xt(2,6) = xl(2,8)

      call pzero(pt,6)
      call pzero(st,36)

      call ale2D3nodedTriangleAspectRatio(d,xt,st,pt,2,6,ndat)

      do i=1, 6
        p(i+2) = p(i+2) + pt(i)
        do j=1, 6
          s(i+2,j+2) = s(i+2,j+2) + st(i,j)
        enddo
      enddo

c  triangle 3 - 4 - 1

      xt(1,1) = xl(1,3)
      xt(2,1) = xl(2,3)
      xt(1,2) = xl(1,4)
      xt(2,2) = xl(2,4)
      xt(1,3) = xl(1,1)
      xt(2,3) = xl(2,1)

      xt(1,4) = xl(1,7)
      xt(2,4) = xl(2,7)
      xt(1,5) = xl(1,8)
      xt(2,5) = xl(2,8)
      xt(1,6) = xl(1,5)
      xt(2,6) = xl(2,5)

      call pzero(pt,6)
      call pzero(st,36)

      call ale2D3nodedTriangleAspectRatio(d,xt,st,pt,2,6,ndat)

      do i=1, 4
        p(i+4) = p(i+4) + pt(i)
        do j=1, 4
          s(i+4,j+4) = s(i+4,j+4) + st(i,j)
        enddo
        do j=1, 2
          s(i+4,j) = s(i+4,j) + st(i,j+4)
          s(j,i+4) = s(j,i+4) + st(j+4,i)
        enddo
      enddo

      do i=1, 2
        p(i) = p(i) + pt(i+4)
        do j=1, 2
          s(i,j) = s(i,j) + st(i+4,j+4)
        enddo
      enddo 

c  triangle 4 - 1 - 2

      xt(1,1) = xl(1,4)
      xt(2,1) = xl(2,4)
      xt(1,2) = xl(1,1)
      xt(2,2) = xl(2,1)
      xt(1,3) = xl(1,2)
      xt(2,3) = xl(2,2)

      xt(1,4) = xl(1,8)
      xt(2,4) = xl(2,8)
      xt(1,5) = xl(1,5)
      xt(2,5) = xl(2,5)
      xt(1,6) = xl(1,6)
      xt(2,6) = xl(2,6)

      call pzero(pt,6)
      call pzero(st,36)

      call ale2D3nodedTriangleAspectRatio(d,xt,st,pt,2,6,ndat)

      do i=1, 4
        p(i) = p(i) + pt(i+2)
        do j=1, 4
          s(i,j) = s(i,j) + st(i+2,j+2)
        enddo
        do j=1, 2
          s(i,j+6) = s(i,j+6) + st(i+2,j)
          s(j+6,i) = s(j+6,i) + st(j,i+2)
        enddo
      enddo

      do i=1, 2
        p(i+6) = p(i+6) + pt(i)
        do j=1, 2
          s(i+6,j+6) = s(i+6,j+6) + st(i,j)
        enddo
      enddo 

      return

      end


