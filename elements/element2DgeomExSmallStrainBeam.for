
      integer function element2dgeomexsmallstrainbeam
     *           (elmDat,matDat,timDat,xl,ul,iv1,iv2,s,p,
     *            ndf,ndm,nen,nelm,nmat,niv,matId,isw,ELM)


c  ***************************************************************************
c  *                                                                         *
c  *  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  *
c  *  O                                                                   O  *
c  *  O            2D GEOMETRICALLY EXACT LINEAR BEAM ELEMENT             O  *
c  *  O                                                                   O  *
c  *  O                  WITH SMALL STRAIN ELASTICITY                     O  *
c  *  O                                                                   O  *
c  *  O          taken from Zienkiewicz & Taylor, fifth edition           O  *
c  *  O                                                                   O  *
c  *  OoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO  *
c  *                                                                         *
c  *                                                                         *
c  *                                                                         *
c  ***************************************************************************

      implicit none

      integer          ELM ! do not touch

      integer          ndf, ndm, nen, nelm, nmat, niv, matId, isw, 
     *                 l, i, ii, j, jj, ngp,
     *                 nen2, nen3, nen4, nen5, m1, m2, m3

      double precision elmDat(*), timDat(*), matDat(*), 
     *                 ul(ndf,*), xl(ndm,*), 
     *                 iv1(*), iv2(*),
     *                 s(nen*ndf,*), p(*),
     *
     *                 b(2), p6, p7, p8, ddfact,
     *                 dvol, detF, ux, uz, bt,
     *                 dux, duz, dbt, sbt, cbt, detJ, lth0,
     *                 xi(10), wgp(10), shp(2,10), uxn(10), 
     *                 uzn(10), cth0, sth0, E, K, G, Arho,
     *                 NF, SF, BM, EA, EI, GA, EAdv, EIdv, GAdv,
     *                 x(2,10), xn(2,10), x0(2,10), u(3,10), ddu(3,10),
     *                 xb0(5), Bi(3,3), Bj(3,3), hlp(2,2),
     *                 fact, fact1, fact2, fact3, fact4,
     *                 r0, r1d8, r1d4, r1d3, r1d2, r1, r3d2, r2, r4,
     *                 r1d6

      data   r0,   r1d8,   r1d4, r1d2,  r1,  r3d2,   r2,   r4 
     *    / 0.d0, .125d0, .25d0, .5d0, 1.d0, 1.5d0, 2.d0, 4.d0 /
      r1d6    = r1 / 6.d0
      r1d3    = r1 / 3.d0

      if (nelm.lt.6) 
     *  call prgerror(1,'element_2d_geom_ex_small_strain_beam',
     *                  'error in element data!')

      nen2   = nen + nen
      nen3   = nen + nen2
      nen4   = nen + nen3
      nen5   = nen + nen4

      ngp    = 1
      b(1)   = elmDat(1)
      b(2)   = elmDat(2)
      Arho   = elmDat(3)

      EA     = elmDat(4)
      EI     = elmDat(5)
      GA     = elmDat(6)

      p6     = timDat(6)
      p7     = timDat(7)
      p8     = timDat(8)
      ddfact = timDat(9) * p8

      
      element2dgeomexsmallstrainbeam = 0



c=============================================================================
 2    continue ! isw 2 -> int.forces
               ! isw 3 -> int.forces + inertia
               !     4 -> int.forces              & tangent for u_{n+1}
               !     5 -> int.forces + inertia    & tangent for u_{n+1}
               !     6 -> int.forces + inertia    & tangent for u_{n+1}, 
               !                                     du(u_{n+1}), ddu(u_{n+1})
c-----------------------------------------------------------------------------

c...  nodal values

      do i=1, nen
        do j=1, ndf
          u  (j,i) = p6 * ul(j,i     ) + (r1-p6) * ul(j,i+nen )
!         du (j,i) = p7 * ul(j,i+nen2) + (r1-p7) * ul(j,i+nen3)
          ddu(j,i) = p8 * ul(j,i+nen4) + (r1-p8) * ul(j,i+nen5)
        enddo
        do j=1, ndm
          x  (j,i) = p6 * xl(j,i     ) + (r1-p6) * xl(j,i+nen )
          xn (j,i) = xl(j,i+nen)
          x0 (j,i) = xl(j,i+nen2)
          if (abs(x(j,i)-x0(j,i)-u(j,i)).gt.1.d-10) then
            write(*,'(/9x,2(f14.8))') x(j,i), x0(j,i)+u(j,i)
            call prgerror(2,'element_2d_geom_ex_small_strain_beam','A')
          endif
        enddo
      enddo

c...  rotate nodal displacements and compute nodal positions on element axis

      lth0   = sqrt((x0(1,1)-x0(1,2))*(x0(1,1)-x0(1,2))
     *             +(x0(2,1)-x0(2,2))*(x0(2,1)-x0(2,2)))
      xb0(1) = r0
      xb0(2) = lth0
      sth0   = (x0(2,2) - x0(2,1)) / lth0
      cth0   = (x0(1,2) - x0(1,1)) / lth0
      uxn(1) = + u(1,1)*cth0 + u(2,1)*sth0
      uzn(1) = - u(1,1)*sth0 + u(2,1)*cth0
      uxn(2) = + u(1,2)*cth0 + u(2,2)*sth0
      uzn(2) = - u(1,2)*sth0 + u(2,2)*cth0

c...  loop over Gauss points (length)

      call comp_xigp1D(xi,wgp,ngp)

      do l=1, ngp

c...    compute shape functions

        call comp_shp1D(shp,detJ,xi(l),xb0,nen)

c...    compute ux, uz, beta in current Gauss point

        ux  = r0
        uz  = r0
        bt  = r0
        dux = r0
        duz = r0
        dbt = r0
        do i=1, nen
          ux  = ux  + uxn(i) * shp(2,i)
          uz  = uz  + uzn(i) * shp(2,i)
          bt  = bt  + u(3,i) * shp(2,i)
          dux = dux + uxn(i) * shp(1,i)
          duz = duz + uzn(i) * shp(1,i)
          dbt = dbt + u(3,i) * shp(1,i)
        enddo
        sbt = sin(bt)
        cbt = cos(bt)

c...    compute average normal strain, shear strain and curvature

        fact = (r1+dux)*cbt - duz*sbt

        E = dux + r1d2 * (dux*dux + duz*duz)
        G = (r1+dux)*sbt + duz*cbt
        K = dbt * fact

c...    compute material response (elastic)

        NF = EA * E  ! normal force
        SF = GA * G  ! shear force
        BM = EI * K  ! bending moment

c...    multiply with volume element

        dvol  = wgp(l) * detJ
        fact1 = dvol * p6
        NF    = NF * dvol
        SF    = SF * dvol
        BM    = BM * dvol
        EAdv  = EA * fact1
        GAdv  = GA * fact1
        EIdv  = EI * fact1

        fact1 = (+ SF * cbt - BM * dbt * sbt) * p6
        fact2 = (- SF * sbt - BM * dbt * cbt) * p6

c...    calculate residual p and stiffness s

        ii = 0
        do i=1, nen

c         compute B-matrix for node i 

          Bi(1,1) = (r1+dux) * shp(1,i)
          Bi(1,2) = duz * shp(1,i)
          Bi(1,3) = r0
          Bi(2,1) = sbt * shp(1,i)
          Bi(2,2) = cbt * shp(1,i)
          Bi(2,3) = fact * shp(2,i)
          Bi(3,1) = dbt*cbt * shp(1,i)
          Bi(3,2) = - dbt*sbt * shp(1,i)
          Bi(3,3) = fact * shp(1,i) - G*dbt * shp(2,i)

c         internal forces Bi^T S

          p(ii+1) = p(ii+1) - Bi(1,1)*NF - Bi(2,1)*SF - Bi(3,1)*BM
          p(ii+2) = p(ii+2) - Bi(1,2)*NF - Bi(2,2)*SF - Bi(3,2)*BM
          p(ii+3) = p(ii+3) - Bi(1,3)*NF - Bi(2,3)*SF - Bi(3,3)*BM

          if (isw.gt.3) then

c           calculate tangent stiffness s

c           compute C Bi

            do m1=1, 3
              Bi(1,m1) = Bi(1,m1) * EAdv
              Bi(2,m1) = Bi(2,m1) * GAdv
              Bi(3,m1) = Bi(3,m1) * EIdv
            enddo

            jj = 0
            do j=1, nen

c             compute B-matrix for node j

              Bj(1,1) = (r1+dux) * shp(1,j)
              Bj(1,2) = duz * shp(1,j)
              Bj(1,3) = r0
              Bj(2,1) = sbt * shp(1,j)
              Bj(2,2) = cbt * shp(1,j)
              Bj(2,3) = fact * shp(2,j)
              Bj(3,1) = dbt*cbt * shp(1,j)
              Bj(3,2) = - dbt*sbt * shp(1,j)
              Bj(3,3) = fact * shp(1,j) - G*dbt * shp(2,j)

c             compute material part of stiffness matrix Bi^T C Bj

              do m1=1, 3
                do m2=1, 3
                  do m3=1, 3
                    s(ii+m1,jj+m2) = s(ii+m1,jj+m2)
     *                                 + Bi(m3,m1) * Bj(m3,m2)
                  enddo
                enddo
              enddo

c             compute geometrical part of stiffness matrix

              s(ii+1,jj+1) = s(ii+1,jj+1) + shp(1,i)*NF*shp(1,j) * p6
              s(ii+2,jj+2) = s(ii+2,jj+2) + shp(1,i)*NF*shp(1,j) * p6

              fact3 =  shp(1,i)*BM*cbt*shp(1,j) * p6
              fact4 = -shp(1,i)*BM*sbt*shp(1,j) * p6

              s(ii+1,jj+3) = s(ii+1,jj+3) +fact3+shp(1,i)*fact1*shp(2,j)
              s(ii+2,jj+3) = s(ii+2,jj+3) +fact4+shp(1,i)*fact2*shp(2,j)
              s(ii+3,jj+1) = s(ii+3,jj+1) +fact3+shp(2,i)*fact1*shp(1,j)
              s(ii+3,jj+2) = s(ii+3,jj+2) +fact4+shp(2,i)*fact2*shp(1,j)

              s(ii+3,jj+3) = s(ii+3,jj+3) -  
     *                        (shp(2,i)*(SF*G+BM*dbt*fact)*shp(2,j)
     *                        +shp(1,i)*BM*G*shp(2,j)
     *                        +shp(2,i)*BM*G*shp(1,j)) * p6
              jj = jj + ndf
            enddo

          endif

          ii = ii + ndf
        enddo

      enddo

c...  rotate residual and stiffness matrix

      ii = 0
      do i=1, nen

        hlp(1,1) = p(ii+1)
        p(ii+1)  = cth0 * p(ii+1)  - sth0 * p(ii+2)
        p(ii+2)  = sth0 * hlp(1,1) + cth0 * p(ii+2)

        jj = 0
        do j=1, nen
          
          hlp(1,1) = cth0 * s(ii+1,jj+1) - sth0 * s(ii+2,jj+1)
          hlp(2,1) = sth0 * s(ii+1,jj+1) + cth0 * s(ii+2,jj+1)
          hlp(1,2) = cth0 * s(ii+1,jj+2) - sth0 * s(ii+2,jj+2)
          hlp(2,2) = sth0 * s(ii+1,jj+2) + cth0 * s(ii+2,jj+2)

          s(ii+1,jj+1) = hlp(1,1) * cth0 - hlp(1,2) * sth0
          s(ii+1,jj+2) = hlp(1,1) * sth0 + hlp(1,2) * cth0
          s(ii+2,jj+1) = hlp(2,1) * cth0 - hlp(2,2) * sth0
          s(ii+2,jj+2) = hlp(2,1) * sth0 + hlp(2,2) * cth0

          hlp(1,1) = s(ii+3,jj+1)
          s(ii+3,jj+1) = s(ii+3,jj+1) * cth0 - s(ii+3,jj+2) * sth0
          s(ii+3,jj+2) = hlp(1,1)     * sth0 + s(ii+3,jj+2) * cth0

          hlp(1,1) = s(ii+1,jj+3)
          s(ii+1,jj+3) = cth0 * s(ii+1,jj+3) - sth0 * s(ii+2,jj+3)
          s(ii+2,jj+3) = sth0 * hlp(1,1)     + cth0 * s(ii+2,jj+3)

          jj = jj + ndf
        enddo

        ii = ii + ndf
      enddo

c...  body forces

      fact = Arho * lth0 * r1d2

      p(1) = p(1) + fact * b(1)
      p(2) = p(2) + fact * b(2)

      p(4) = p(4) + fact * b(1)
      p(5) = p(5) + fact * b(2)

c...  inertia

      if (isw.ne.2.and.isw.ne.4) then

        fact  = Arho * lth0 * r1d6
        fact2 = fact + fact

        p(1) = p(1) - fact2 * ddu(1,1) - fact * ddu(1,2)
        p(2) = p(2) - fact2 * ddu(2,1) - fact * ddu(2,2)

        p(4) = p(4) - fact2 * ddu(1,2) - fact * ddu(1,1)
        p(5) = p(5) - fact2 * ddu(2,2) - fact * ddu(2,1)

      endif

      if (isw.eq.6) then
       
        fact   = Arho * lth0 * r1d6 * ddfact
        fact2  = fact + fact
 
        s(1,1) = s(1,1) + fact2
        s(2,2) = s(2,2) + fact2

        s(4,4) = s(4,4) + fact2
        s(5,5) = s(5,5) + fact2

        s(1,4) = s(1,4) + fact
        s(2,5) = s(2,5) + fact

        s(4,1) = s(4,1) + fact
        s(5,2) = s(5,2) + fact

      endif

      return

cc==========================================================================
c 7    continue ! tangent for ddu_{n+1} (consistent mass matrix)
cc--------------------------------------------------------------------------
c
cc...  nodal values
c
c      do i=1, nen
c        do j=1, ndm
c          x0 (j,i) = xl(j,i+nen2)
c        enddo
c      enddo
c
c      lth0   = sqrt((x0(1,1)-x0(1,2))*(x0(1,1)-x0(1,2))
c     *             +(x0(2,1)-x0(2,2))*(x0(2,1)-x0(2,2)))
c      xb0(1) = r0
c      xb0(2) = lth0
c
cc     loop over Gauss points (length)
c
c      call comp_xigp1D(xi,wgp,ngp)
c
c      do l=1, ngp
c
cc       compute shape functions and volume element
c
c        call comp_shp1D(shp,detJ,xi(l),xb0,nen)
c
c        dvol = wgp(l) * detJ
c
c        fact4 = dvol * (dvol*dvol*Arho*r1d3+Irho)
c
c        ii = 0
c        do i=1, nen
c          fact1 = shp(2,i) * dvol  * Arho * p8
c          fact2 = shp(2,i) * fact4 * p8
c          jj = 0
c          do j=1, nen
c            fact = fact1 * shp(2,j)
c            s(ii+1,jj+1) = s(ii+1,jj+1) + fact
c            s(ii+2,jj+2) = s(ii+2,jj+2) + fact
c            jj = jj + ndf
c          enddo
c          ii = ii + ndf
c        enddo
c
c      enddo

      return

c=================================================================================

      end


      
                 
