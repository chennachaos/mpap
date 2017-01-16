
      integer function ale2D3nodedTriangleLinearElasticity
     *                    (d,xl,s,p,ndm,nst,ndat)
      implicit none
 
      integer          ndm, nst, ndat,
     *                 i, j, m, ii, jj, i0, j0
      double precision d(*), xl(ndm,*), s(nst,*), p(*),
     *                 x(2,3), x0(2,3), u(2,3), mu, bulk, shp(2,3), 
     *                 fact, gu(2,2), trgu, r(2), k(2,2), area,
     *                 r0, r1d3, r1d2, r2d3, r1
c
c     LINEAR TRIANGLES,
c
c     MOVE MESH BY SOLVING LINEAR PSEUDO-ELASTIC PROBLEM
c 

      data r0, r1d2, r1 / 0.d0, 0.5d0, 1.d0 /

      r1d3 = r1 / 3.d0
      r2d3 = r1d3 * 2.d0

      mu   = d(1)
      bulk = d(2)

      if (ndat.lt.3) 
     *  call prgerror(1,'ale2D3nodedTriangleLinearElasticity',
     *                  '2 parameters required!')

      ale2D3nodedTriangleLinearElasticity = 0

      do i=1, 3
        do j=1, ndm
          x(j,i) = xl(j,i) + xl(j,i+3)
        enddo
      enddo
      if (nint(d(3)).ne.0) then
        do i=1, 3
          do j=1, ndm
            x0(j,i) = xl(j,i)
            u (j,i) = xl(j,i+3)
          enddo
        enddo
      else
        do i=1, 3
          do j=1, ndm
            x0(j,i) = xl(j,i+6)
            u (j,i) = x(j,i) - x0(j,i) 
          enddo
        enddo
      endif

c     new area
      area = r1d2 * ((x(1,1)-x(1,3))*(x(2,2)-x(2,3)) 
     *             - (x(1,2)-x(1,3))*(x(2,1)-x(2,3)))

      if (area.lt.1.d-16) then
        ale2D3nodedTriangleLinearElasticity = -1
        return
      endif

c     old area
      area = r1d2 * ((x0(1,1)-x0(1,3))*(x0(2,2)-x0(2,3)) 
     *             - (x0(1,2)-x0(1,3))*(x0(2,1)-x0(2,3)))

      do i=1, 3   !  loop over nodes
        j = i+1
        m = i+2
        if (j.gt.3) j = j - 3
        if (m.gt.3) m = m - 3
        fact     = r1 / ( (x0(2,i)-x0(2,j)) * (x0(1,j)-x0(1,m))
     *                   -(x0(2,j)-x0(2,m)) * (x0(1,i)-x0(1,j)))
        shp(1,i) = - fact * (x0(2,j) - x0(2,m))
        shp(2,i) =   fact * (x0(1,j) - x0(1,m))
      enddo
      
      call pzero(gu,4)                                   
      trgu = r0                                            !
      do i=1, 3                                            !
        do ii=1, 2                                         !
          trgu        = trgu      + u(ii,i) * shp(ii,i) ! tr(grad(u))
          do jj=1, 2                                       !
            gu(ii,jj) = gu(ii,jj) + u(ii,i) * shp(jj,i) ! grad(u)
          enddo                                            !
        enddo                                              !
      enddo                                                !

      do i=1, 3   !   nodes 
        i0 = (i-1)*2

        fact = r2d3 * trgu
        r(1) = mu * (shp(1,i) * (gu(1,1) + gu(1,1) - fact)
     *             + shp(2,i) * (gu(1,2) + gu(2,1)))
        r(2) = mu * (shp(1,i) * (gu(1,2) + gu(2,1))
     *             + shp(2,i) * (gu(2,2) + gu(2,2) - fact))

        fact = bulk * trgu * r1d3
        r(1) = r(1) + fact * shp(1,i) 
        r(2) = r(2) + fact * shp(2,i)
        
        do ii=1, 2
          p(i0+ii) = p(i0+ii) - r(ii) * area
        enddo

        do j=1, 3
          j0 = (j-1)*2

          fact = shp(1,i)*shp(1,j) + shp(2,i)*shp(2,j)
          k(1,1) = mu * (r1d3*shp(1,j)*shp(1,i) + fact)
          k(1,2) = mu * (shp(1,j)*shp(2,i) - r2d3*shp(1,i)*shp(2,j))
          k(2,1) = mu * (shp(2,j)*shp(1,i) - r2d3*shp(2,i)*shp(1,j))
          k(2,2) = mu * (r1d3*shp(2,j)*shp(2,i) + fact)

          fact = bulk * r1d3
          k(1,1) = k(1,1) + fact * shp(1,i) * shp(1,j)
          k(1,2) = k(1,2) + fact * shp(1,i) * shp(2,j)
          k(2,1) = k(2,1) + fact * shp(2,i) * shp(1,j)
          k(2,2) = k(2,2) + fact * shp(2,i) * shp(2,j)

          do ii=1, 2
            do jj=1, 2
              s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj) * area
            enddo
          enddo
        enddo

      enddo

      return
      end











