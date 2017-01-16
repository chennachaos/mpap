
      integer function ale2D3nodedTriangleCellCentroid
     *                    (d,xl,s,p,ndm,nst,ndat)
      implicit none

      integer          ndm, nst, ndat,
     *                 i, j, ii, i0, j0, jj
      double precision d(*), xl(ndm,*), s(nst,*), p(*), 
     *                 x(2,3), A, dA(2,3), xc(2), Dx(2,3),
     *                 dDx(2,3,2,3),
     *                 r1d3, r2d3

      r1d3 = 1.d0 / 3.d0
      r2d3 = r1d3 + r1d3

c
c     UPDATE FOR TRIANGLES, BASED ON CELL CENTROID
c
      ale2D3nodedTriangleCellCentroid = 0

      do i=1, 3
        do j=1, ndm
          x(j,i) = xl(j,i) + xl(j,i+3)
        enddo
      enddo

      A =  (x(1,1)-x(1,3))*(x(2,2)-x(2,3)) 
     *   - (x(1,2)-x(1,3))*(x(2,1)-x(2,3)) 

      if (A.lt.1.d-20) then 
        ale2D3nodedTriangleCellCentroid = -1
        return
      endif

      dA(1,1) = x(2,2) - x(2,3)
      dA(2,1) = x(1,3) - x(1,2)
      dA(1,2) = x(2,3) - x(2,1)
      dA(2,2) = x(1,1) - x(1,3)
      dA(1,3) = x(2,1) - x(2,2)
      dA(2,3) = x(1,2) - x(1,1)

      xc(1) = r1d3 * (x(1,1) + x(1,2) + x(1,3))
      xc(2) = r1d3 * (x(2,1) + x(2,2) + x(2,3))

      Dx(1,1) = x(1,1) - xc(1)
      Dx(2,1) = x(2,1) - xc(2)
      Dx(1,2) = x(1,2) - xc(1)
      Dx(2,2) = x(2,2) - xc(2)
      Dx(1,3) = x(1,3) - xc(1)
      Dx(2,3) = x(2,3) - xc(2)

      call pzero(dDx,36)

      dDx(1,1,1,1) =   r2d3
      dDx(1,1,1,2) = - r1d3
      dDx(1,1,1,3) = - r1d3

      dDx(2,1,2,1) =   r2d3
      dDx(2,1,2,2) = - r1d3
      dDx(2,1,2,3) = - r1d3

      dDx(1,2,1,1) = - r1d3
      dDx(1,2,1,2) =   r2d3
      dDx(1,2,1,3) = - r1d3

      dDx(2,2,2,1) = - r1d3
      dDx(2,2,2,2) =   r2d3
      dDx(2,2,2,3) = - r1d3

      dDx(1,3,1,1) = - r1d3
      dDx(1,3,1,2) = - r1d3
      dDx(1,3,1,3) =   r2d3

      dDx(2,3,2,1) = - r1d3
      dDx(2,3,2,2) = - r1d3
      dDx(2,3,2,3) =   r2d3

      do i=1, 3
        i0 = (i-1)*2

        do ii=1, 2
          p(i0+ii) = p(i0+ii) - A * Dx(ii,i)

          do j=1, 3
            j0 = (j-1)*2

            do jj=1, 2

              s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + Dx(ii,i) * dA(jj,j)
     *                                        + A * dDx(ii,i,jj,j)
 
            enddo

          enddo

        enddo

      enddo

      !write(*,'(6(1x,g12.5)/)') (p(j),j=1,6)
      !do i=1, 6
      !  write(*,'(6(1x,g12.5))') (s(i,j),j=1,6)
      !enddo
      !call presskey(1)

      return
      end




