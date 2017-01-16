
      integer function element3dsolid
     *           (elmDat,matDat,timDat,xl,ul,iv1,iv2,s,p,
     *            ndf,ndm,nen,nelm,nmat,nivGP,matId,isw,err,ELM)

      implicit none

      integer          ELM ! do not touch

      integer          ndf, ndm, nen, nelm, nmat, nivGP, matId, isw,
     *                 err,
     *                 i, j, ngp, finInt,
     *                 nen2, nen3, nen4, nen5, ii, jj, l,
     *                 round, k
 
           double precision elmDat(*), matDat(*), timDat(*), 
     *                 xl(ndm,*), ul(ndf,*), 
     *                 iv1(*), iv2(*), s(ndf*nen,*), p(*),
     *                 x(3,20), xn(3,20), x0(3,20), u(3,20), ddu(3,20),
     *                 shp(4,20), wgp(27), xi(3,27),
     *                 F(3,3), stre(6), cc(6,6), b(3), dvol0, rho0,
     *                 detJ, detF, dvol,  
     *                 bc(3,6), fact, fact1, fact2, fact3, 
     *                 ddfact, p6, p7, p8,
     *                 r1

      logical          finite

      r1 = 1.d0

      nen2   = nen + nen
      nen3   = nen + nen2
      nen4   = nen + nen3
      nen5   = nen + nen4
      
      ngp    = round(elmDat(1))
      finInt = round(elmDat(2))
      finite = (finInt.ge.1)
      
      rho0   = elmDat(3)
      b(1)   = elmDat(4)
      b(2)   = elmDat(5)
      b(3)   = elmDat(6)

      
      p6     = timDat(6)
      p7     = timDat(7)
      p8     = timDat(8)
      ddfact = timDat(9) * p8 * rho0
      
      element3DSolid = 0
      

     
      if ((nen.ne.4).and.(nen.ne.10).and.(nen.ne.8).and.
     *(nen.ne.20)) call prgError
     * (1,'element3DSolid','nen does not exist in 3D')

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
            call prgWarning(2,'element3DSolid','u + x0 != x')
            err=1
            return
          endif 
        enddo
      enddo

c...  Compute Gauss points coordinates and weights

        if (nen.eq.4) then
            ngp=1
            xi(1,1) = .25d0
            xi(2,1) = .25d0
            xi(3,1) = .25d0
            wgp(1)  = .25d0/6.d0
        else
            if (nen.eq.10) ngp=4
            call compxigp3D(xi,wgp,ngp)
        endif

c...  loop over Gauss points
      do l=1, ngp

c...    compute shape functions on undeformed configuration 
c       (for small strain formulation, also for inertia and body forces)

        call compshp3D(shp,detJ,xi(1,l),x0,nen)

c...    compute volume for current Gauss point in undeformed configuration
c       (for inertia and body forces)

        dvol0 = wgp(l)  * detJ

c...    compute shape funtions and deformation gradient F

        if (finite) then
          call compshp3D(shp,detJ,xi(1,l),x,nen)
          call F3d(shp,u,F,detF,ndf,nen,2)
        else
          call F3d(shp,u,F,detF,ndf,nen,1)
        endif
	
c...    compute material response

        call matlib3d(matDat,F,stre,cc,
     *                iv1((l-1)*nivGP+1),iv2((l-1)*nivGP+1),timDat(5),
     *                matId,nivGP,finInt,3,err,l,ELM)

  

	if (err.ne.0) return

        if (isw.eq.7) then

          do i=1, 6
            p((l-1)*6+i) = stre(i)
          enddo   

          goto 7

        endif
c...    multiply stress and tangent tensor with volume element

        dvol = wgp(l) * detJ

        do i=1, 6
          stre(i) = stre(i) * dvol
          do j=1, 6
            cc(i,j) = cc(i,j) * dvol
          enddo
        enddo

c...    calculate residual: p

        
c       internal forces

        ii = 0
        do i=1, nen
      p(ii+1)=p(ii+1)-shp(1,i)*stre(1)-shp(2,i)*stre(4)-shp(3,i)*stre(6)
      p(ii+2)=p(ii+2)-shp(1,i)*stre(4)-shp(2,i)*stre(2)-shp(3,i)*stre(5)
      p(ii+3)=p(ii+3)-shp(1,i)*stre(6)-shp(2,i)*stre(5)-shp(3,i)*stre(3)
          ii = ii + ndf
        enddo
 	
c       body forces

        fact  = rho0 * dvol0
        fact1 = b(1) * fact
        fact2 = b(2) * fact
        fact3 = b(3) * fact
        ii = 0
        do i=1, nen
          p(ii+1) = p(ii+1) + shp(4,i) * fact1
          p(ii+2) = p(ii+2) + shp(4,i) * fact2
          p(ii+3) = p(ii+3) + shp(4,i) * fact3
          ii = ii + ndf
        enddo

c       inertia  

        if (isw.ne.2.and.isw.ne.4) then
          fact1 = 0.d0
          fact2 = 0.d0
          fact3 = 0.d0
          do i=1, nen
            fact1 = fact1 + shp(4,i) * ddu(1,i)
            fact2 = fact2 + shp(4,i) * ddu(2,i)
            fact3 = fact3 + shp(4,i) * ddu(3,i)
          enddo
          fact  = rho0 * dvol0
          fact1 = fact1 * fact
          fact2 = fact2 * fact
          fact3 = fact3 * fact
          ii = 0
          do i=1, nen
            p(ii+1) = p(ii+1) - shp(4,i) * fact1
            p(ii+2) = p(ii+2) - shp(4,i) * fact2
            p(ii+3) = p(ii+3) - shp(4,i) * fact3
            ii = ii + ndf
          enddo
        endif

        if (isw.gt.3) then

c...      calculate tangent stiffness: s

c         internal forces

c         part 1. -- geometrical matrix  (if geometry nonlinear)

          if (finite) then
            ii = 0
            do i=1, nen
      fact1 = (shp(1,i)*stre(1)+shp(2,i)*stre(4)+shp(3,i)*stre(6)) * p6
      fact2 = (shp(1,i)*stre(4)+shp(2,i)*stre(2)+shp(3,i)*stre(5)) * p6
      fact3 = (shp(1,i)*stre(6)+shp(2,i)*stre(5)+shp(3,i)*stre(3)) * p6
              jj  = 0
              do j=1, nen
           fact = fact1 * shp(1,j) + fact2 * shp(2,j) + fact3 * shp(3,j)
                s(ii+1,jj+1) = s(ii+1,jj+1) + fact
                s(ii+2,jj+2) = s(ii+2,jj+2) + fact
                s(ii+3,jj+3) = s(ii+3,jj+3) + fact
                jj = jj + ndm
              enddo
              ii = ii + ndm
            enddo
          endif


c         part 2. -- material part (not necessarily symmetric!!)

          ii = 0
          do i=1, nen
          
          do k=1, 6
       bc(1,k)=(shp(1,i)*cc(1,k)+shp(2,i)*cc(4,k)+shp(3,i)*cc(6,k))*p6
       bc(2,k)=(shp(1,i)*cc(4,k)+shp(2,i)*cc(2,k)+shp(3,i)*cc(5,k))*p6
       bc(3,k)=(shp(1,i)*cc(6,k)+shp(2,i)*cc(5,k)+shp(3,i)*cc(3,k))*p6
          enddo

            jj = 0
            do j=1, nen
            
            do k=1, 3
              s(ii+k,jj+1) = s(ii+k,jj+1)
     *           + shp(1,j)*bc(k,1)+ shp(2,j)*bc(k,4)+ shp(3,j)*bc(k,6) 
              s(ii+k,jj+2) = s(ii+k,jj+2)
     *           + shp(1,j)*bc(k,4)+ shp(2,j)*bc(k,2)+ shp(3,j)*bc(k,5)
              s(ii+k,jj+3) = s(ii+k,jj+3)
     *           + shp(1,j)*bc(k,6)+ shp(2,j)*bc(k,5)+ shp(3,j)*bc(k,3) 
            enddo
 

              jj = jj + ndm
            enddo
            ii = ii + ndm
          enddo



          if (isw.eq.6) then

c           inertia

            ii = 0
            do i=1, nen
              fact = shp(4,i) * dvol0 * ddfact
              jj = 0
              do j=1, nen
                fact1 = fact * shp(4,j)
                s(ii+1,jj+1) = s(ii+1,jj+1) + fact1
                s(ii+2,jj+2) = s(ii+2,jj+2) + fact1
                s(ii+3,jj+3) = s(ii+3,jj+3) + fact1
                jj = jj + ndf
              enddo
              ii = ii + ndf
            enddo

          endif
        endif
7      continue
      enddo
      


      end

