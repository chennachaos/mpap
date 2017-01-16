  
      integer function element2dFBarSolid
     *           (elmDat,matDat,timDat,xl,ul,iv1,iv2,s,p,
     *            ndf,ndm,nen,nelm,nmat,nivGP,matId,isw,err,ELM)

      implicit none

      integer          ELM ! do not touch

      integer          ndf, ndm, nen, nelm, nmat, nivGP, matId,
     *                 i, j, ngp, sss, finiteInt, isw, err,
     *                 nen2, nen3, nen4, nen5, ii, jj, l,
     *                 round

      double precision elmDat(*), matDat(*), timDat(*), 
     *                 xl(ndm,*), ul(ndf,*), 
     *                 iv1(*), iv2(*), s(ndf*nen,*), p(*),
     *                 x(2,4), xn(2,4), x0(2,4), u(2,4), ddu(2,4),
     *                 shp(3,4), wgp(9), xi(2,9),
     *                 F(4), stre(4), cc(4,4), b(2), dvol0, rho0,
     *                 detJ, detF, dvol, thick, F33, r1dF33, r1drad0,
     *                 Jc, trc,cch(4),shpc(3,4),JcdJ,
     *                 bc(2,4), fact, fact1, fact2, fact3, 
     *                 fact4, ddfact, rad0, rad, p6, p7, p8,
     *                 r1d2, r0, r1, r2, r2pi, pi

      logical          finite

      data r0, r1d2, r1, r2 / 0.d0, .5d0, 1.d0, 2.d0 /

            if (nen.ne.4) call prgError(2,'element2d4nodedfbarsolid',
     *      '4 nodes only!')
     
      pi        = acos(-r1)
      r2pi      = pi * r2

      nen2      = nen + nen
      nen3      = nen + nen2
      nen4      = nen + nen3
      nen5      = nen + nen4

      ngp       = round(elmDat(1))
      finiteInt = round(elmDat(2))
      sss       = round(elmDat(3))
      thick     = round(elmDat(4))

      finite    = (finiteInt.eq.1)
      
    
      if (sss.eq.3) call prgError(2,'element2d4nodedfbarsolid(f)',
     *      'axisymmetry not working!')

      if (ngp.ne.4) call prgError(2,'element2d4nodedfbarsolid(f)',
     *      'use four Gauss points: ngp = 4 !')

      b(1)      = elmDat(5)
      b(2)      = elmDat(6)
      rho0      = elmDat(7)

      p6        = timDat(6)
      p7        = timDat(7)
      p8        = timDat(8)
      ddfact    = timDat(9) * p8 * rho0

      element2DFBarSolid = 0
      
c=============================================================================
      continue ! isw 2 -> int.forces
               ! isw 3 -> int.forces + inertia
               !     4 -> int.forces              & tangent for u_{n+1}
               !     5 -> int.forces + inertia    & tangent for u_{n+1}
               !     6 -> int.forces + inertia    & tangent for u_{n+1}, 
               !                                     du(u_{n+1}), ddu(u_{n+1})
               !     7 -> store the stresses in p, leave s untouched
c-----------------------------------------------------------------------------

c...  nodal values

      do i=1, nen
        do j=1, ndm
          x  (j,i) = p6 * xl(j,i     ) + (r1-p6) * xl(j,i+nen )
          xn (j,i) = xl(j,i+nen)
          x0 (j,i) = xl(j,i+nen2)
          u  (j,i) = p6 * ul(j,i     ) + (r1-p6) * ul(j,i+nen )
          ddu(j,i) = p8 * ul(j,i+nen4) + (r1-p8) * ul(j,i+nen5)
          if (abs(x(j,i)-x0(j,i)-u(j,i)).gt.1.d-10) then
            write(*,'(/9x,2(f14.8))') x(j,i), x0(j,i)+u(j,i)
            call prgError(3,'element_2d_489noded_solid','u + x0 != x')
          endif 
        enddo
      enddo
      
      
c...  compute determinant detF (tr(F)) in element centre 

        call compxigp2D(xi,wgp,1)
        if (finite) then
          call compshp2D(shpc,detJ,xi,x,nen)
          call F2d(shpc,u,F,detF,ndf,nen,2)
          Jc = detF
        else
          call compshp2D(shpc,detJ,xi,x0,nen)
          call F2d(shpc,u,F,detF,ndf,nen,1)
          trc = F(1) + F(4)
        endif


c...  loop over Gauss points

      call compxigp2D(xi,wgp,ngp)
      
      do l=1, ngp 
      

c...    compute shape functions on undeformed configuration 
c       (for small strain formulation, also for inertia and body forces)

        call compshp2D(shp,detJ,xi(1,l),x0,nen)


c...    compute volume for current Gauss point in undeformed configuration
c       (for inertia and body forces)

        dvol0 = detJ * wgp(l) * thick

c...    compute shape funtions and deformation gradient F

        if (finite) then
          call compshp2D(shp,detJ,xi(1,l),x,nen)
          call F2d(shp,u,F,detF,ndf,nen,2)
        else
          call F2d(shp,u,F,detF,ndf,nen,1)
        endif

        if     (sss.le.1.or.sss.eq.4) then ! plane stress
          F33 = r1 / sqrt(detF)
        elseif (sss.eq.2.or.sss.eq.5) then ! plane strain
          F33 = r1
        endif
        r1dF33 = r1 / F33

c...    replace detF (tr(F)) with detF (tr(F)) from element centre

          if (finite) then
            JcdJ = Jc / detF
            fact = sqrt(JcdJ)
            do i=1, 4
              F(i) = F(i) * fact
            enddo
          else
            fact = r1d2 * (trc - F(1) - F(4))
            F(1) = F(1) + fact
            F(4) = F(4) + fact
          endif
          

c...    compute material response

        call matlib2d(matDat,F,F33,stre,cc,
     *                iv1((l-1)*nivGP+1),iv2((l-1)*nivGP+1),timDat(5),
     *                matId,nivGP,finiteInt,sss,3,err,l,ELM)
    
        if (err.ne.0) return

        if (isw.eq.7) then

          do i=1, 4
            p((l-1)*4+i) = stre(i)
          enddo   

          goto 7

        endif
        
c...    multiply stress and tangent tensor with volume element

        dvol = wgp(l) * detJ * thick
        if (finite) dvol = dvol * F33
        do i=1, 4
          stre(i) = stre(i) * dvol
          do j=1, 4
            cc(i,j) = cc(i,j) * dvol
          enddo
        enddo

c...    calculate residual p

c       internal forces

        ii = 0

        do i=1, nen
          p(ii+1) = p(ii+1) - shp(1,i)*stre(1) - shp(2,i)*stre(3)
          p(ii+2) = p(ii+2) - shp(1,i)*stre(3) - shp(2,i)*stre(2)
          ii = ii + ndf
        enddo

c       body forces

        fact  = rho0 * dvol0
        fact1 = b(1) * fact
        fact2 = b(2) * fact
        ii = 0
        do i=1, nen
          p(ii+1) = p(ii+1) + shp(3,i) * fact1
          p(ii+2) = p(ii+2) + shp(3,i) * fact2
          ii = ii + ndf
        enddo

c       inertia

        if (isw.ne.2.and.isw.ne.4) then
          fact1 = 0.d0
          fact2 = 0.d0
          do i=1, nen
            fact1 = fact1 + shp(3,i) * ddu(1,i)
            fact2 = fact2 + shp(3,i) * ddu(2,i)
          enddo
          fact  = rho0 * dvol0
          fact1 = fact1 * fact
          fact2 = fact2 * fact
          ii = 0
          do i=1, nen
            p(ii+1) = p(ii+1) - shp(3,i) * fact1
            p(ii+2) = p(ii+2) - shp(3,i) * fact2
            ii = ii + ndf
          enddo
        endif

        if (isw.gt.3) then
        
c...      calculate tangent stiffness s

c         internal forces

c         part 1. -- geometrical matrix  (if geometry nonlinear)           

          if (finite) then
            ii = 0
            do i=1, nen
              fact1 = (shp(1,i) * stre(1) + shp(2,i) * stre(3)) * p6
              fact2 = (shp(1,i) * stre(3) + shp(2,i) * stre(2)) * p6
              jj  = 0
              do j=1, nen
                fact         = fact1 * shp(1,j) + fact2 * shp(2,j)
                s(ii+1,jj+1) = s(ii+1,jj+1) + fact
                s(ii+2,jj+2) = s(ii+2,jj+2) + fact
                jj = jj + ndm
              enddo
              ii = ii + ndm
            enddo
          endif

c         part 2. -- material part (not necessarily symmetric!!)

c         correct tangent cc for fbar formulation
            cch(1)  = r1d2 * (cc(1,1) + cc(1,2))
            cch(2)  = r1d2 * (cc(2,1) + cc(2,2))
            cch(3)  = r1d2 * (cc(3,1) + cc(3,2))
            cc(1,1) = cc(1,1) - cch(1)
            cc(2,1) = cc(2,1) - cch(2)
            cc(3,1) = cc(3,1) - cch(3)
            cc(1,2) = cc(1,2) - cch(1) 
            cc(2,2) = cc(2,2) - cch(2)
            cc(3,2) = cc(3,2) - cch(3)


          ii = 0
          do i=1, nen
            bc(1,1) = (shp(1,i) * cc(1,1) + shp(2,i) * cc(3,1)) * p6
            bc(1,2) = (shp(1,i) * cc(1,2) + shp(2,i) * cc(3,2)) * p6
            bc(1,3) = (shp(1,i) * cc(1,3) + shp(2,i) * cc(3,3)) * p6
            bc(2,1) = (shp(2,i) * cc(2,1) + shp(1,i) * cc(3,1)) * p6
            bc(2,2) = (shp(2,i) * cc(2,2) + shp(1,i) * cc(3,2)) * p6
            bc(2,3) = (shp(2,i) * cc(2,3) + shp(1,i) * cc(3,3)) * p6

            jj = 0
            do j=1, nen
              s(ii+1,jj+1) = s(ii+1,jj+1) + bc(1,1) * shp(1,j) 
     *                                    + bc(1,3) * shp(2,j) 
              s(ii+2,jj+1) = s(ii+2,jj+1) + bc(2,1) * shp(1,j) 
     *                                    + bc(2,3) * shp(2,j) 
              s(ii+1,jj+2) = s(ii+1,jj+2) + bc(1,2) * shp(2,j) 
     *                                    + bc(1,3) * shp(1,j) 
              s(ii+2,jj+2) = s(ii+2,jj+2) + bc(2,2) * shp(2,j) 
     *                                    + bc(2,3) * shp(1,j) 

              jj = jj + ndm
            enddo

c...        derivative with respect to shp(centroid)
            
            fact1 = (shp(1,i) * cch(1) + shp(2,i) * cch(3)) * p6
            fact2 = (shp(2,i) * cch(2) + shp(1,i) * cch(3)) * p6
            jj = 0
            do j=1, nen
              s(ii+1,jj+1) = s(ii+1,jj+1) + fact1 * shpc(1,j) 
              s(ii+2,jj+1) = s(ii+2,jj+1) + fact2 * shpc(1,j) 
              s(ii+1,jj+2) = s(ii+1,jj+2) + fact1 * shpc(2,j) 
              s(ii+2,jj+2) = s(ii+2,jj+2) + fact2 * shpc(2,j) 
              jj = jj + ndm
            enddo

            ii = ii + ndm
          enddo
          


          if (isw.eq.6) then

c           inertia

            ii = 0
            do i=1, nen
              fact = shp(3,i) * dvol0 * ddfact
              jj = 0
              do j=1, nen
                fact1 = fact * shp(3,j)
                s(ii+1,jj+1) = s(ii+1,jj+1) + fact1
                s(ii+2,jj+2) = s(ii+2,jj+2) + fact1
                jj = jj + ndf
              enddo

              ii = ii + ndf
            enddo

          endif
          
        endif
        


 7      continue

      enddo


 8      continue

      return


      end



