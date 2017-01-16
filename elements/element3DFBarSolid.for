
      integer function element3dFBarSolid
     *           (elmDat,matDat,timDat,xl,ul,iv1,iv2,s,p,
     *            ndf,ndm,nen,nelm,nmat,nivGP,matId,isw,err,GPDOM)

      implicit none

      integer          GPDOM

      integer          ndf, ndm, nen, nelm, nmat, nivGP, matId, isw,
     *                 err,
     *                 i, j, ngp, finInt,
     *                 nen2, nen3, nen4, nen5, ii, jj, l,
     *                 round, k
 
           double precision elmDat(*), matDat(*), timDat(*), 
     *                 xl(ndm,*), ul(ndf,*), 
     *                 iv1(*), iv2(*), s(ndf*nen,*), p(*),
     *                 x(3,8), xn(3,8), x0(3,8), u(3,8), ddu(3,8),
     *                 shp(4,8), wgp(27), xi(3,27),
     *                 F(9), stre(7), cc(6,6), b(3), dvol0, rho0,
     *                 detJ, detF, dvol,
     *                 bc(3,6), fact, fact1, fact2, fact3, 
     *                 ddfact, p6, p7, p8,  
     *                 r1, r2, r1d3, r2d3, r3,r1d2,
     *                 trc, cch(6), shpc(4,8), Jc, JcdJ,
     *                 pres, sdev(6), dot_ss

      logical          finite
      
      data r1, r2, r3/ 1.d0, 2.d0 , 3.d0 /

      if (nen.ne.8) call prgError(1,'element3d8nodedfbarsolid',
     *      '8 nodes only!')
     
      r1d3=r1/r3
      r2d3=r2/r3
      r1d2 = 0.5d0
      
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
      
      element3DFBarSolid = 0
      
      
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
            call prgWarning(3,'element3D820nodedSolid','u + x0 != x')
            err=1
            return
          endif 
        enddo
      enddo
      

c...  compute determinant detF (tr(F)) in element centre 

        call compxigp3D(xi,wgp,1)
        if (finite) then
          call compshp3D(shpc,detJ,xi,x,nen)
          call F3d(shpc,u,F,detF,ndf,nen,2)
          Jc = detF
        else
          call compshp3D(shpc,detJ,xi,x0,nen)
          call F3d(shpc,u,F,detF,ndf,nen,1)
          trc = F(1) + F(5) + F(9)
        endif
       
      
       call compxigp3D(xi,wgp,ngp)

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
        
     
c...    replace detF with detF(bar) from element centre
            
            if (finite) then
            JcdJ = Jc / detF
            do i=1, 9
            F(i) = F(i) * JcdJ**r1d3
            enddo
            else
            fact = r1d3 * (trc - F(1) - F(5)- F(9))
            F(1) = F(1) + fact
            F(5) = F(5) + fact
            F(9) = F(9) + fact
            endif
            



c...    compute material response

        call matlib3d(matDat,F,stre,cc,
     *                iv1((l-1)*nivGP+1),iv2((l-1)*nivGP+1),timDat(5),
     *                matId,nivGP,finInt,3,err,l-1,GPDOM)

        if (err.ne.0) return


        if (isw.eq.7) then
        
           ! von mises stress
          
            pres=(stre(1)+stre(2)+stre(3))/3.0

            sdev(1)=stre(1)-pres
            sdev(2)=stre(2)-pres
            sdev(3)=stre(3)-pres
            sdev(4)=stre(4)
            sdev(5)=stre(5)
            sdev(6)=stre(6)

            stre(7)=(dot_ss(sdev,sdev)*1.5)**r1d2
  
          do i=1, 7
            p((l-1)*7+i) = stre(i)
          enddo   

          goto 7
          
        endif
        
c...    multiply stress and tangent tensor with volume element

        dvol = wgp(l) * detJ
        
         do i=1, 6
          stre(i) = stre(i) * dvol 
          if (finite) stre(i) = stre(i) * JcdJ**r1d3

          do j=1, 6
            cc(i,j) = cc(i,j) * dvol 
            if (finite) cc(i,j) = cc(i,j) * JcdJ**r1d3
          enddo
          
             cch(i)  =  r1d3 * (cc(i,1) + cc(i,2) + cc(i,3))
             
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
      fact1 = (shp(1,i)*stre(1)+shp(2,i)*stre(4)+shp(3,i)*stre(6))* p6
      fact2 = (shp(1,i)*stre(4)+shp(2,i)*stre(2)+shp(3,i)*stre(5))* p6 
      fact3 = (shp(1,i)*stre(6)+shp(2,i)*stre(5)+shp(3,i)*stre(3))* p6
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
          
          
            do j=1, nen
            jj = (j-1)*ndm
            do k=1, 3
              s(ii+k,jj+1) = s(ii+k,jj+1)
     *           + shp(1,j)*bc(k,1)+ shp(2,j)*bc(k,4)+ shp(3,j)*bc(k,6) 
              s(ii+k,jj+2) = s(ii+k,jj+2)
     *           + shp(1,j)*bc(k,4)+ shp(2,j)*bc(k,2)+ shp(3,j)*bc(k,5)
              s(ii+k,jj+3) = s(ii+k,jj+3)
     *           + shp(1,j)*bc(k,6)+ shp(2,j)*bc(k,5)+ shp(3,j)*bc(k,3) 
            enddo
            enddo
            
c...        derivative with respect to shp(centroid)

            fact1 = (shp(1,i) * cch(1) + shp(2,i) * cch(4)
     *             + shp(3,i) * cch(6)) * p6 
            fact2 = (shp(1,i) * cch(4) + shp(2,i) * cch(2)
     *             + shp(3,i) * cch(5)) * p6 
            fact3 = (shp(1,i) * cch(6) + shp(2,i) * cch(5)
     *             + shp(3,i) * cch(3)) * p6 


              do j=1, nen
         jj = (j-1)*ndm
         s(ii+1,jj+1) = s(ii+1,jj+1) + fact1 * (shpc(1,j)-shp(1,j))
         s(ii+2,jj+1) = s(ii+2,jj+1) + fact2 * (shpc(1,j)-shp(1,j))
         s(ii+3,jj+1) = s(ii+3,jj+1) + fact3 * (shpc(1,j)-shp(1,j))
          
         s(ii+1,jj+2) = s(ii+1,jj+2) + fact1 * (shpc(2,j)-shp(2,j))
         s(ii+2,jj+2) = s(ii+2,jj+2) + fact2 * (shpc(2,j)-shp(2,j))
         s(ii+3,jj+2) = s(ii+3,jj+2) + fact3 * (shpc(2,j)-shp(2,j))
          
         s(ii+1,jj+3) = s(ii+1,jj+3) + fact1 * (shpc(3,j)-shp(3,j))
         s(ii+2,jj+3) = s(ii+2,jj+3) + fact2 * (shpc(3,j)-shp(3,j))
         s(ii+3,jj+3) = s(ii+3,jj+3) + fact3 * (shpc(3,j)-shp(3,j))
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

