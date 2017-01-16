
      integer function element2DTruss
     *           (elmDat,matDat,timDat,xl,ul,iv1,iv2,s,p,dpress,
     *            ndf,ndm,nen,nelm,nmat,niv,matId,isw,ELM)



c    ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c    O                                                                   O
c    O                 2D 2 NODED LINEAR TRUSS ELEMENT                   O
c    O                                                                   O
c    O                     SMALL / LARGE  STRAINS                        O
c    O                                                                   O
c    O                   ALSO AXISYMMETRIC MEMBRANE                      O
c    O                                                                   O
c    OoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO

      implicit none

      integer          ELM ! do not touch

      integer          ndf, ndm, nen, nelm, nmat, niv, matId, isw,
     *                 i, l, ii, j, jj, ngp, nivgp,  
     *                 nen2, nen3, nen4, nen5, finInt, sss

      double precision elmDat(*), matDat(*), timDat(*), 
     *                 ul(ndf,*), xl(ndm,*), iv1(*), iv2(*),
     *                 s(ndf*nen,*), p(*), dpress,
     *                 x(2,2), x0(2,2), xn(2,2), u(2,2), du(2,2), b(2),
     *                 ddu(2,2), stre(3), cc(2,2), dlam(3), lam(3),  
     *                 rho0, A0, lam0,
     *                 lth, lth0, rad, rad0, dvol, dvol0,
     *                 xi(10), wgp(10), shp(2), 
     *                 factx, facty, factxx, factxy, factyy,
     *                 fact1x, fact1y, fact2x, fact2y, factl, factr,
     *                 factl0, factr0,
     *                 fctxx, fctxy, fctyy, fctdxx, fctdxy, fctdyy, 
     *                 fct1x, fct1y, fct2x, fct2y, fct,
     *                 p6, p7, p8, ddfact, 
     *                 r1d2, r1, r2, pi, r2pi, pid3

      logical          finite, axsy, rope

      data r1d2, r1, r2 / .5d0, 1.d0, 2.d0 /

      pi     = acos(-r1)
      r2pi   = pi + pi
      pid3   = pi / 3.d0

      if (nen.ne.2) call prgError(1,"element2DTruss","nen <> 2!")

      nen2   = nen + nen
      nen3   = nen + nen2
      nen4   = nen + nen3
      nen5   = nen + nen4

      ngp    = min(3,max(1,nint(elmDat(1))))
      finInt = nint(elmDat(2))
      finite = (finInt.ge.1)
      sss    = nint(elmDat(3))
      axsy   = (sss.eq.3.or.sss.eq.6)
      A0     = elmDat(4)       ! A0 is thickness for axisymmetry
      b(1)   = elmDat(5)
      b(2)   = elmDat(6)
      rho0   = elmDat(7)
      lam0   = elmDat(8) + r1

      rope   = (elmDat(11).eq.1)
      if (rope.and.axsy)
     *     call prgError(2,"element2DTruss","rope and axsy impossible!")

      p6     = timDat(6)
      p7     = timDat(7)
      p8     = timDat(8)
      ddfact = timDat(9) * p8

      if (lam0.lt.1.d-16) 
     *  call prgError(1,'element2DTruss','invalid initial strain!')

      if (abs(lam0-r1).gt.1.d-16.and.axsy) 
     *  call prgError(1,'element2DTruss',
     *         'initial strain not applicable for axisymmetric case !')


      element2DTruss = 0


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
        do j=1, ndm
          x  (j,i) = p6 * xl(j,i     ) + (r1-p6) * xl(j,i+nen )
          u  (j,i) = p6 * ul(j,i     ) + (r1-p6) * ul(j,i+nen )
!         du (j,i) = p7 * ul(j,i+nen2) + (r1-p7) * ul(j,i+nen3)
          ddu(j,i) = p8 * ul(j,i+nen4) + (r1-p8) * ul(j,i+nen5)
          xn (j,i) = xl(j,i+nen)
          x0 (j,i) = xl(j,i+nen2)
          if (abs(x(j,i)-x0(j,i)-u(j,i)).gt.1.d-10) then
            write(*,'(/9x,2(f14.8))') x(j,i), x0(j,i)+u(j,i)
            call prgError(3,'element2DTruss','A')
          endif
        enddo
      enddo

c...  prestressed element

      fct     = r1 / lam0
      x0(1,2) = x0(1,1) + fct * (x0(1,2)-x0(1,1))
      x0(2,2) = x0(2,1) + fct * (x0(2,2)-x0(2,1))

c...  initial and current element length

      lth0    = sqrt( (x0(1,2)-x0(1,1))*(x0(1,2)-x0(1,1))
     *               +(x0(2,2)-x0(2,1))*(x0(2,2)-x0(2,1)))
      if (finite) then
        lth   = sqrt( (x (1,2)-x (1,1))*(x (1,2)-x (1,1))
     *               +(x (2,2)-x (2,1))*(x (2,2)-x (2,1)))
        factx = (x (1,2) - x (1,1)) / lth
        facty = (x (2,2) - x (2,1)) / lth
      else
        lth   =     ( (x (1,2)-x (1,1))*(x0(1,2)-x0(1,1))
     *               +(x (2,2)-x (2,1))*(x0(2,2)-x0(2,1))) / lth0
        factx = (x0(1,2) - x0(1,1)) / lth0
        facty = (x0(2,2) - x0(2,1)) / lth0
      endif
      lam(1)  = lth / lth0
      factxx  = factx * factx
      factxy  = factx * facty
      factyy  = facty * facty

c...  external pressure

      if (axsy) then

        fct1x = x(1,1) + x(1,1) + x(1,2)
        fct2x = x(1,1) + x(1,2) + x(1,2)

        fct  = dpress * pid3
        p(1) = (x(2,1)-x(2,2)) * fct1x * fct
        p(2) = (x(1,2)-x(1,1)) * fct1x * fct
        p(3) = (x(2,1)-x(2,2)) * fct2x * fct
        p(4) = (x(1,2)-x(1,1)) * fct2x * fct

        fct    = - fct * p6 

        s(1,1) = + (x(2,1)-x(2,2)) * (fct + fct)
        s(1,2) = + fct1x * fct 
        s(1,3) = + (x(2,1)-x(2,2)) * fct 
        s(1,4) = - fct1x * fct

        s(2,1) = - (fct1x - (x(1,2)-x(1,1))*r2) * fct
        s(2,3) = + (fct1x +  x(1,2)-x(1,1))     * fct

        s(3,1) = + (x(2,1)-x(2,2)) * fct
        s(3,2) = + fct2x * fct
        s(3,3) = + (x(2,1)-x(2,2)) * (fct + fct)
        s(3,4) = - fct2x * fct

        s(4,1) = - (fct2x -  x(1,2)+x(1,1))     * fct
        s(4,3) = + (fct2x + (x(1,2)-x(1,1))*r2) * fct

      else

        fct  = dpress * r1d2
        p(1) = (x(2,1) - x(2,2)) * fct
        p(2) = (x(1,2) - x(1,1)) * fct
        p(3) = (x(2,1) - x(2,2)) * fct
        p(4) = (x(1,2) - x(1,1)) * fct

        fct    = - fct * p6 
        s(1,2) = + fct 
        s(1,4) = - fct
        s(2,1) = - fct
        s(2,3) = + fct
        s(3,2) = + fct
        s(3,4) = - fct
        s(4,1) = - fct
        s(4,3) = + fct

      endif

c...  body forces

      if (axsy) then

        fct1x = x0(1,1) + x0(1,1) + x0(1,2)
        fct2x = x0(1,1) + x0(1,2) + x0(1,2)

        fct  = rho0 * A0 * lth0 * pid3
        p(1) = p(1) + fct * fct1x * b(1)
        p(2) = p(2) + fct * fct1x * b(2)
        p(3) = p(3) + fct * fct2x * b(1)
        p(4) = p(4) + fct * fct2x * b(2)

      else

        fct  = rho0 * A0 * lth0 * r1d2 
        p(1) = p(1) + fct * b(1)
        p(2) = p(2) + fct * b(2)
        p(3) = p(3) + fct * b(1)
        p(4) = p(4) + fct * b(2)

        !write(*,*) fct * b(1)

      endif

      if (isw.ne.2.and.isw.ne.4) then

c...    inertia

        if (axsy) then

          p(1) = p(1) - fct * fct1x * ddu(1,1)
          p(2) = p(2) - fct * fct1x * ddu(2,1)
          p(3) = p(3) - fct * fct2x * ddu(1,2)
          p(4) = p(4) - fct * fct2x * ddu(2,2)

          if (isw.eq.6) then
            fct = fct * ddfact
            s(1,1) = s(1,1) + fct * fct1x
            s(2,2) = s(2,2) + fct * fct1x
            s(3,3) = s(3,3) + fct * fct2x
            s(4,4) = s(4,4) + fct * fct2x
          endif

        else

          p(1) = p(1) - fct * ddu(1,1)
          p(2) = p(2) - fct * ddu(2,1)
          p(3) = p(3) - fct * ddu(1,2)
          p(4) = p(4) - fct * ddu(2,2)

          if (isw.eq.6) then
            fct = fct * ddfact
            do i=1, 4
              s(i,i) = s(i,i) + fct
            enddo
          endif

        endif

      endif

c...  loop over Gauss points x

      call comp_xigp1D(xi,wgp,ngp)

      do l=1, ngp

c...    compute shape functions & radius for axisymmetry

        if (axsy) then
          shp(1) = r1d2 * (r1 - xi(l))
          shp(2) = r1d2 * (r1 + xi(l))
          rad0   = shp(1) * x0(1,1) + shp(2) * x0(1,2)
          rad    = shp(1) * x (1,1) + shp(2) * x (1,2)
          lam(3) = rad / rad0
          fact1x = shp(1) * factx
          fact1y = shp(1) * facty
          fact2x = shp(2) * factx
          fact2y = shp(2) * facty
        endif

c...    compute material response

        call matlib1D(matDat,lam,stre,cc,dlam,iv1,iv2,timDat(5),
     *                matId,nivGp,finInt,sss,3,l,ELM)

        !write(*,*) lam(1), stre(1)

        if (rope.and.stre(1).lt.0.d0) then
          stre(1) = 0.d0
          cc(1,1) = 0.d0
        endif

c...    calculate residual

c       internal forces

        factl0 = A0 * wgp(l) * r1d2
        if (axsy) factl0 = factl0 * rad0 * r2pi
        if (finite) then 
          factl = factl0 * lam(2) * lam(3)
        else
          factl = factl0
        endif

        fct  = factl * stre(1)
        p(1) = p(1) + fct * factx
        p(2) = p(2) + fct * facty
        p(3) = p(3) - fct * factx
        p(4) = p(4) - fct * facty

        if (axsy) then

          factr0 = A0 * wgp(l) * lth0 * pi
          if (finite) then
            factr = factr0 * lam(1) * lam(2)
          else
            factr = factr0
          endif

          fct  = factr * stre(2)
          p(1) = p(1) - fct * shp(1)
          p(3) = p(3) - fct * shp(2)

        endif

        if (isw.gt.3) then

c...      calculate tangent stiffness

c         internal forces
  
c         material stiffness
   
          fct    = factl / lth0 * cc(1,1) * p6
          fctxx  = fct * factxx
          fctxy  = fct * factxy
          fctyy  = fct * factyy

          s(1,1) = s(1,1) + fctxx
          s(1,2) = s(1,2) + fctxy
          s(1,3) = s(1,3) - fctxx
          s(1,4) = s(1,4) - fctxy

          s(2,1) = s(2,1) + fctxy
          s(2,2) = s(2,2) + fctyy
          s(2,3) = s(2,3) - fctxy
          s(2,4) = s(2,4) - fctyy

          s(3,1) = s(3,1) - fctxx
          s(3,2) = s(3,2) - fctxy
          s(3,3) = s(3,3) + fctxx
          s(3,4) = s(3,4) + fctxy

          s(4,1) = s(4,1) - fctxy
          s(4,2) = s(4,2) - fctyy
          s(4,3) = s(4,3) + fctxy
          s(4,4) = s(4,4) + fctyy

          if (axsy) then

            fct    = factl / rad0 * cc(1,2) * p6
            fct1x  = fct * fact1x
            fct1y  = fct * fact1y
            fct2x  = fct * fact2x
            fct2y  = fct * fact2y

            s(1,1) = s(1,1) - fct1x
            s(1,3) = s(1,3) - fct2x

            s(2,1) = s(2,1) - fct1y
            s(2,3) = s(2,3) - fct2y

            s(3,1) = s(3,1) + fct1x
            s(3,3) = s(3,3) + fct2x

            s(4,1) = s(4,1) + fct1y
            s(4,3) = s(4,3) + fct2y


            fct    = factr / lth0 * cc(2,1) * p6 
            fct1x  = fct * fact1x
            fct1y  = fct * fact1y
            fct2x  = fct * fact2x
            fct2y  = fct * fact2y

            s(1,1) = s(1,1) - fct1x
            s(1,2) = s(1,2) - fct1y
            s(1,3) = s(1,3) + fct1x
            s(1,4) = s(1,4) + fct1y

            s(3,1) = s(3,1) - fct2x
            s(3,2) = s(3,2) - fct2y
            s(3,3) = s(3,3) + fct2x
            s(3,4) = s(3,4) + fct2y


            fct    = factr / rad0 * cc(2,2) * p6 

            s(1,1) = s(1,1) + fct * shp(1) * shp(1)
            s(1,3) = s(1,3) + fct * shp(1) * shp(2)

            s(3,1) = s(3,1) + fct * shp(2) * shp(1)
            s(3,3) = s(3,3) + fct * shp(2) * shp(2)

          endif

c         geometrical stiffness
   
          if (finite) then

c           derivative with respect to lambda_1
            if (axsy) then
              fct  = lam(3) * dlam(1) * factl0 * stre(1) / lth0 * p6
            else
              fct  = (dlam(2) * lam(3) + lam(2) * dlam(3))
     *              * factl0 * stre(1) / lth0 * p6
            endif
            fctxx  = fct * factxx
            fctxy  = fct * factxy
            fctyy  = fct * factyy

c           derivative with respect to xfact and yfact
            fct    = stre(1) * factl / lth * p6
            fctdxx = (- r1 + factxx) * fct
            fctdxy = (     + factxy) * fct
            fctdyy = (- r1 + factyy) * fct

            s(1,1) = s(1,1) + fctxx - fctdxx
            s(1,2) = s(1,2) + fctxy - fctdxy
            s(1,3) = s(1,3) - fctxx + fctdxx
            s(1,4) = s(1,4) - fctxy + fctdxy
                                   
            s(2,1) = s(2,1) + fctxy - fctdxy
            s(2,2) = s(2,2) + fctyy - fctdyy
            s(2,3) = s(2,3) - fctxy + fctdxy
            s(2,4) = s(2,4) - fctyy + fctdyy
                                   
            s(3,1) = s(3,1) - fctxx + fctdxx
            s(3,2) = s(3,2) - fctxy + fctdxy
            s(3,3) = s(3,3) + fctxx - fctdxx
            s(3,4) = s(3,4) + fctxy - fctdxy
                                   
            s(4,1) = s(4,1) - fctxy + fctdxy
            s(4,2) = s(4,2) - fctyy + fctdyy
            s(4,3) = s(4,3) + fctxy - fctdxy
            s(4,4) = s(4,4) + fctyy - fctdyy

            if (axsy) then

              fct    = (lam(2) + lam(3) * dlam(3))
     *                * factl0 * stre(1) / rad0 * p6

              s(1,1) = s(1,1) - fct * fact1x
              s(1,3) = s(1,3) - fct * fact2x
                                           
              s(2,1) = s(2,1) - fct * fact1y
              s(2,3) = s(2,3) - fct * fact2y
                                           
              s(3,1) = s(3,1) + fct * fact1x
              s(3,3) = s(3,3) + fct * fact2x
                                           
              s(4,1) = s(4,1) + fct * fact1y
              s(4,3) = s(4,3) + fct * fact2y


              fct    =  ( lam(1) * dlam(1) +  lam(2))
     *                * factr0 * stre(2) / lth0 * p6

              s(1,1) = s(1,1) - fct * fact1x
              s(1,2) = s(1,2) - fct * fact1y
              s(1,3) = s(1,3) + fct * fact1x
              s(1,4) = s(1,4) + fct * fact1y
                                           
              s(3,1) = s(3,1) - fct * fact2x
              s(3,2) = s(3,2) - fct * fact2y
              s(3,3) = s(3,3) + fct * fact2x
              s(3,4) = s(3,4) + fct * fact2y


              fct    = lam(1) * dlam(3) * factr0 * stre(2) / rad0 * p6

              s(1,1) = s(1,1) + fct * shp(1) * shp(1)
              s(1,3) = s(1,3) + fct * shp(1) * shp(2)

              s(3,1) = s(3,1) + fct * shp(2) * shp(1)
              s(3,3) = s(3,3) + fct * shp(2) * shp(2)

            endif

          endif

        endif

      enddo

      return

c=============================================================================

      end



