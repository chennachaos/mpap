
      integer function ale3D4nodedTetrahedronHyperElasticity
     *           (matDat,xl,s,p,matId,finInt)

      implicit none

      integer          ndf, ndm, nen, nelm, nmat, nivGP, matId, 
     *                 err,
     *                 i, j, ngp, finInt,
     *                 nen2, nen3, nen4, nen5, ii, jj, l, k
 
      double precision matDat(*),  xl(3,12), s(12,*), p(*),
     *                 x(3,4), x0(3,4), u(3,4), shp(4,4), wgp(27), 
     *                 xi(3,27), F(3,3), stre(6), cc(6,6), dvol0, 
     *                 detJ, detF, dvol,bc(3,6),   
     *                 fact, fact1, fact2, fact3, r1

      logical          finite
      r1 = 1.d0
      
      nen=4
      ndf=3
      ndm=3

      nen2   = nen + nen
      nen3   = nen + nen2
      nen4   = nen + nen3
      nen5   = nen + nen4
      
      finite = (finInt.ge.1)

      
      ale3D4nodedTetrahedronHyperElasticity= 0


c...  nodal values

      do i=1, 4
        do j=1, ndm
          x(j,i) = xl(j,i) + xl(j,i+4)
        enddo
      enddo
      if (nint(matDat(3)).eq.0) then
        do i=1, 4
          do j=1, ndm
            x0(j,i) = xl(j,i)
            u (j,i) = xl(j,i+4)
          enddo
        enddo
      else
        do i=1, 4
          do j=1, ndm
            x0(j,i) = xl(j,i+8)
            u (j,i) = x(j,i) - x0(j,i) 
          enddo
        enddo
      endif

c...  Compute Gauss points coordinates and weights


            ngp=1
            xi(1,1) = .25d0
            xi(2,1) = .25d0
            xi(3,1) = .25d0
            wgp(1)  = r1/6.d0
            l=1


c...    compute shape functions on undeformed configuration 
c       (for small strain formulation, also for inertia and body forces)

        call compshp3D(shp,detJ,xi,x0,nen)

c...    compute volume for current Gauss point in undeformed configuration
c       (for inertia and body forces)

        dvol0 = wgp(l)  * detJ

c...    compute shape funtions and deformation gradient F

        if (finite) then
          call compshp3D(shp,detJ,xi,x,nen)
          call F3d(shp,u,F,detF,ndf,nen,2)
        else
          call F3d(shp,u,F,detF,ndf,nen,1)
        endif
	
c...    compute material response

        call matlib3d(matDat,F,stre,cc,
     *                0.,0.,0.,
     *                matId,nivGP,finInt,3,err,l,0)

  

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
 	

c...      calculate tangent stiffness: s

c         internal forces

c         part 1. -- geometrical matrix  (if geometry nonlinear)

          if (finite) then
            ii = 0
            do i=1, nen
      fact1 = (shp(1,i)*stre(1)+shp(2,i)*stre(4)+shp(3,i)*stre(6)) 
      fact2 = (shp(1,i)*stre(4)+shp(2,i)*stre(2)+shp(3,i)*stre(5)) 
      fact3 = (shp(1,i)*stre(6)+shp(2,i)*stre(5)+shp(3,i)*stre(3)) 
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
       bc(1,k)=(shp(1,i)*cc(1,k)+shp(2,i)*cc(4,k)+shp(3,i)*cc(6,k))
       bc(2,k)=(shp(1,i)*cc(4,k)+shp(2,i)*cc(2,k)+shp(3,i)*cc(5,k))
       bc(3,k)=(shp(1,i)*cc(6,k)+shp(2,i)*cc(5,k)+shp(3,i)*cc(3,k))
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

      end

