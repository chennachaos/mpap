
      integer function ale2D3nodedTriangleAspectRatio
     *                    (d,xl,s,p,ndm,nst,ndat)
      implicit none

      integer          ndm, nst, ndat,
     *                 ABC(3,3,3), i, j, ii, i0, j0, jj, pow
      double precision d(*), xl(ndm,*), s(nst,*), p(*), 
     *                 x(2,3), r(2), k(2,2), s2(3), ds2(3,2,3), ri2(4), 
     *                 ro2(4), dri2(3,3), dro2(3,3), dWds(3), 
     *                 ddWds(3,3), fact, r2

      data r2 / 2.d0 /

      data ABC 
     *  /  1,  0,  1, -1,  0,  0,  0,  0, -1,  
     *    -1,  0,  0,  1,  1,  0,  0, -1,  0,
     *     0,  0, -1,  0, -1,  0,  0,  1,  1  /
c
c     LINEAR TRIANGLES,
c
c     OPTIMISE MESH BY MINIMISING   W = Sum (ro/ri)^2
c 
      
      ale2D3nodedTriangleAspectRatio = 0

      pow = 2

      do i=1, 3
        do j=1, ndm
          x(j,i) = xl(j,i) + xl(j,i+3)
        enddo
      enddo

      if ((x(1,1)-x(1,3))*(x(2,2)-x(2,3)) 
     *  - (x(1,2)-x(1,3))*(x(2,1)-x(2,3)).lt.1.d-20) then 
        ale2D3nodedTriangleAspectRatio = -1
        return
      endif

      call pzero(ds2,18)

      do i=1, 3  
        j = i + 1
        if (j.gt.3) j = j - 3
        s2(i) =  (x(1,i)-x(1,j)) * (x(1,i)-x(1,j))
     *         + (x(2,i)-x(2,j)) * (x(2,i)-x(2,j))
        ds2(i,1,i) =   r2 * (x(1,i)-x(1,j))
        ds2(i,2,i) =   r2 * (x(2,i)-x(2,j))
        ds2(i,1,j) = - r2 * (x(1,i)-x(1,j))
        ds2(i,2,j) = - r2 * (x(2,i)-x(2,j))
      enddo

      if (sqrt(min(s2(1),min(s2(2),s2(3)))).lt.1.d-12) then
        ale2D3nodedTriangleAspectRatio = -1
        return
      endif

      call comp_rin2 (s2,ri2,dri2)
      call comp_rout2(s2,ro2,dro2)

      call comp_dW(ri2,ro2,dri2,dro2,pow,dWds,ddWds)

      do i=1, 3      !   nodes 
        i0 = (i-1)*2
        call pzero(r,2)

        do ii=1, 3
          r(1) = r(1) + dWds(ii) * ds2(ii,1,i)
          r(2) = r(2) + dWds(ii) * ds2(ii,2,i)
        enddo

        do ii=1, 2
          p(i0+ii) = p(i0+ii) - r(ii)
        enddo

        do j=1, 3
          j0 = (j-1)*2
          call pzero(k,4)

          do ii=1, 3 
            do jj=1, 3
              fact = ddWds(ii,jj) * ds2(ii,1,i)
              k(1,1) = k(1,1) + fact * ds2(jj,1,j)
              k(1,2) = k(1,2) + fact * ds2(jj,2,j)
              fact = ddWds(ii,jj) * ds2(ii,2,i) 
              k(2,1) = k(2,1) + fact * ds2(jj,1,j) 
              k(2,2) = k(2,2) + fact * ds2(jj,2,j) 
            enddo
            fact = r2 * dWds(ii) * ABC(ii,i,j)
            k(1,1) = k(1,1) + fact
            k(2,2) = k(2,2) + fact
          enddo

          do ii=1, 2
            do jj=1, 2
              s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj)
            enddo
          enddo
        enddo

      enddo

      return

      end





c===========================================================================
      subroutine comp_rin2(s2,ri2,dri2)
c---------------------------------------------------------------------------
      implicit none
      integer          i, j, k
      double precision s2(3), ri2(4), dri2(3,3), s(3), 
     *                 h1, h2, h3, h4, dh1(3), dh2(3), dh3(3), dh4(3),
     *                 ddh1ii, ddh2ii, ddh2ij, ddh2ik, ddh3ii, ddh3ij, 
     *                 ddh3ik, ddh4ii, r1, r2, r4, r1d2, r1d4, r3d2, 
     *                 r3d4

      data r1, r2, r4, r1d2, r1d4, r3d2, r3d4 
     *        / 1.d0, 2.d0, 4.d0, 0.5d0, 0.25d0, 1.5d0, 0.75d0 /

      do i=1, 3
        s(i) = sqrt(s2(i))
      enddo

      h1 = s2(1)*s(1) + s2(2)*s(2) + s2(3)*s(3)
      h2 = r2 * s(1) * s(2) * s(3)
      h3 =      s2(1)*s(2) + s2(2)*s(1)
     *        + s2(2)*s(3) + s2(3)*s(2)
     *        + s2(3)*s(1) + s2(1)*s(3)
      h4 = r4 * (s(1)+s(2)+s(3))

      ri2(4) = (h3 - h1 - h2) / h4 

      do i=1, 3
        j = i + 1
        k = i + 2
        if (j.gt.3) j = j - 3
        if (k.gt.3) k = k - 3

        dh1(i) = r3d2 * s(i)
        dh2(i) = s(j) * s(k) / s(i)
        dh3(i) = s(j) + s(k) + r1d2 * (s2(j) + s2(k)) / s(i)
        dh4(i) = r2 / s(i)

        ri2(i) = (dh3(i) - dh1(i) - dh2(i) - ri2(4) * dh4(i)) / h4
      enddo

      do i=1, 3
        j = i + 1
        k = i + 2
        if (j.gt.3) j = j - 3
        if (k.gt.3) k = k - 3

        ddh1ii = r3d4 / s(i)
        !ddh1ij = 0
        !ddh1ik = 0

        ddh2ii = - r1d2 * dh2(i) / s2(i) 
        ddh2ij =   r1d2 * dh2(i) / s2(j)
        ddh2ik =   r1d2 * dh2(i) / s2(k)

        ddh3ii = - r1d4 * (s2(j) + s2(k)) / (s(i) * s2(i))
        ddh3ij = r1d2 / s(j) + r1d2 / s(i)
        ddh3ik = r1d2 / s(k) + r1d2 / s(i)

        ddh4ii = - r1 / (s(i) * s2(i))
        !ddh4ij = 0
        !ddh4ik = 0

        dri2(i,i) = (ddh3ii - ddh1ii   - ddh2ii - ri2(4)*ddh4ii 
     *               - ri2(i) * dh4(i) - ri2(i) * dh4(i)) / h4
        dri2(i,j) = (ddh3ij - ddh2ij ! - ddh1ij - ri2(4)*ddh4ij 
     *               - ri2(j) * dh4(i) - ri2(i) * dh4(j)) / h4
        dri2(i,k) = (ddh3ik - ddh2ik ! - ddh1ik - ri2(4)*ddh4ik 
     *               - ri2(k) * dh4(i) - ri2(i) * dh4(k)) / h4

      enddo

      return
      end


c===========================================================================
      subroutine comp_rout2(s2,ro2,dro2)
c---------------------------------------------------------------------------
      implicit none
      integer          i, j, k
      double precision s2(3), ro2(4), dro2(3,3), h1, h2, h3, h4, 
     *                 dh1, dh2, dh3, dh4(3), r2, fact

      data r2 / 2.d0 /

      h1 = s2(1)*s2(2)*s2(3)
      h2 = s2(1)*s2(2) + s2(2)*s2(3) + s2(3)*s2(1)
      h3 = s2(1)*s2(1) + s2(2)*s2(2) + s2(3)*s2(3)
      h4 = r2 * h2 - h3

      fact = 1.d0 / h4

      ro2(4) = h1 * fact

      do i=1, 3
        j = i + 1
        k = i + 2
        if (j.gt.3) j = j - 3
        if (k.gt.3) k = k - 3
        dh1    = s2(j) * s2(k)
        dh2    = s2(j) + s2(k)
        dh3    = r2 * s2(i)
        dh4(i) = r2 * dh2 - dh3

        ro2(i) = (dh1 - ro2(4) * dh4(i)) * fact

      enddo

      do i=1, 3
        j = i + 1
        k = i + 2
        if (j.gt.3) j = j - 3
        if (k.gt.3) k = k - 3

        !ddh1ii = 0
        !ddh1ij = s2(k)
        !ddh1ik = s2(j)

        !ddh4ii = - 2
        !ddh4ij =  2
        !ddh4ik =  2
                            
        dro2(i,i) = (     -ro2(i)*dh4(i)+r2*ro2(4)-ro2(i)*dh4(i))*fact
        dro2(i,j) = (s2(k)-ro2(j)*dh4(i)-r2*ro2(4)-ro2(i)*dh4(j))*fact
        dro2(i,k) = (s2(j)-ro2(k)*dh4(i)-r2*ro2(4)-ro2(i)*dh4(k))*fact
      enddo
      return
      end


c===========================================================================
      subroutine comp_dW(ri2,ro2,dri2,dro2,pow,dWds,ddWds)
c---------------------------------------------------------------------------
      implicit none
      integer          i, j, pow
      double precision ri2(4), ro2(4), dri2(3,3), dro2(3,3), 
     *                 dWds(3), ddWds(3,3), fact, fact1, fact2, fact3,
     *                 fact4, fact5, fact6, W,
     *                 r1, r1d2

      data r1, r1d2 / 1.d0, 0.5d0 /

      fact6 = - r1 / ri2(4)
      fact  = - ro2(4) * fact6
      fact1 = sqrt(fact)
      W     = fact1
      do i=2, pow
        W = W * fact1
      enddo

      fact1 = W / ro2(4)
      fact2 = W * fact6

      fact3 = (float(pow) * r1d2 - r1) / W
      fact4 = W * fact6
      fact5 = W / ro2(4)

      do i=1, 3
        dWds(i) = fact1 * ro2(i) + fact2 * ri2(i)
        do j=1, i
          ddWds(i,j) =  fact3 * dWds(i) * dWds(j)
     *                + fact4 * dri2(i,j) + fact5 * dro2(i,j)
     *                + fact6 * (dWds(i) * ri2(j) + dWds(j) * ri2(i))
        enddo
      enddo
      ddWds(1,2) = ddWds(2,1)
      ddWds(1,3) = ddWds(3,1)
      ddWds(2,3) = ddWds(3,2)

      return
      end

c===========================================================================





