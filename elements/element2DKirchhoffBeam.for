
      integer function element2dkirchhoffbeam
     *           (elmDat,matDat,timDat,xl,ul,iv1,iv2,s,p,
     *            ndf,ndm,nen,nelm,nmat,niv,matId,isw,ELM)

      implicit none

      integer          ELM ! do not touch

      integer          ndf, ndm, nen, nelm, nmat, niv, matId, isw, 
     *                 l, i, ii, j, jj, ngp,
     *                 nen2, nen3, nen4, nen5, m1, m2, m3

      double precision elmDat(*), timDat(*), matDat(*), 
     *                 ul(ndf,*), xl(ndm,*), 
     *                 iv1(*), iv2(*), s(nen*ndf,*), p(*),
     *                 b(2), p6, p7, p8, ddfact,
     *                 lth0, ux(2), uz(2), cth0, sth0, Arho,
     *                 EA, EI, x0(2,2), u(3,2), ddu(3,2), hlp(2,2),
     *                 fact, fact2, fact3, 
     *                 r0, r1d2, r1, r4, r6, r12, r1d6

      data   r0,   r1d2,  r1,   r4,   r6,   r12 
     *    / 0.d0, 0.5d0, 1.d0, 4.d0, 6.d0, 12.d0 /
      r1d6 = r1 / 6.d0

      if (nelm.lt.6) 
     *  call prgerror(1,'element2dkirchhoffbeam',
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

      p6     = timDat(6)
      p7     = timDat(7)
      p8     = timDat(8)
      ddfact = timDat(9) * p8

      
      element2dkirchhoffbeam = 0



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
         !du (j,i) = p7 * ul(j,i+nen2) + (r1-p7) * ul(j,i+nen3)
          ddu(j,i) = p8 * ul(j,i+nen4) + (r1-p8) * ul(j,i+nen5)
        enddo
        do j=1, ndm
         !x  (j,i) = p6 * xl(j,i     ) + (r1-p6) * xl(j,i+nen )
         !xn (j,i) = xl(j,i+nen)
          x0 (j,i) = xl(j,i+nen2)
         !if (abs(x(j,i)-x0(j,i)-u(j,i)).gt.1.d-10) then
         !  write(*,'(/9x,2(f14.8))') x(j,i), x0(j,i)+u(j,i)
         !  call prgerror(2,'element2dkirchhoffbeam','A')
         !endif
        enddo
      enddo

c...  rotate nodal displacements and compute nodal positions on element axis

      lth0  = sqrt((x0(1,1)-x0(1,2))*(x0(1,1)-x0(1,2))
     *            +(x0(2,1)-x0(2,2))*(x0(2,1)-x0(2,2)))
      sth0  = (x0(2,2) - x0(2,1)) / lth0
      cth0  = (x0(1,2) - x0(1,1)) / lth0
      ux(1) = + u(1,1)*cth0 + u(2,1)*sth0
      uz(1) = - u(1,1)*sth0 + u(2,1)*cth0
      ux(2) = + u(1,2)*cth0 + u(2,2)*sth0
      uz(2) = - u(1,2)*sth0 + u(2,2)*cth0

c...  residual and stiffness matrix with respect to axial direction

      s(1,1) = EA / lth0
      s(1,4) = - s(1,1)
      s(4,1) = - s(1,1)
      s(4,4) = s(1,1)

      fact   = EI / lth0
      fact2  = fact / lth0
      fact3  = fact2 / lth0

      s(2,2) = r12 * fact3
      s(2,3) =  r6 * fact2
      s(3,2) = s(2,3)
      s(3,3) = r4 * fact

      s(5,5) = s(2,2)
      s(5,6) = - s(2,3)
      s(6,5) = - s(3,2)
      s(6,6) = s(3,3)

      s(2,5) = - s(2,2)
      s(2,6) = s(2,3)
      s(3,5) = - s(2,3)
      s(3,6) = fact + fact

      s(5,2) = s(2,5)
      s(6,2) = s(2,6)
      s(5,3) = s(3,5)
      s(6,3) = s(3,6)

      p(1) = -s(1,1)*ux(1) -s(1,4)*ux(2)
      p(2) = -s(2,2)*uz(1) -s(2,5)*uz(2) -s(2,3)*u(3,1) -s(2,6)*u(3,2)
      p(3) = -s(3,2)*uz(1) -s(3,5)*uz(2) -s(3,3)*u(3,1) -s(3,6)*u(3,2)

      p(4) = -s(4,1)*ux(1) -s(4,4)*ux(2)
      p(5) = -s(5,2)*uz(1) -s(5,5)*uz(2) -s(5,3)*u(3,1) -s(5,6)*u(3,2)
      p(6) = -s(6,2)*uz(1) -s(6,5)*uz(2) -s(6,3)*u(3,1) -s(6,6)*u(3,2)

      do i=1, 2
        do j=1, 2
          s(1+i,1+j) = s(1+i,1+j) * p6
          s(1+i,4+j) = s(1+i,4+j) * p6
          s(4+i,4+j) = s(4+i,4+j) * p6
          s(4+i,1+j) = s(4+i,1+j) * p6
        enddo
      enddo
      s(1,1) = s(1,1) * p6
      s(1,4) = s(1,4) * p6
      s(4,4) = s(4,4) * p6
      s(4,1) = s(4,1) * p6

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

c==========================================================================

      end


      
                 
