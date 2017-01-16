
      subroutine small_strain_von_Mises_elastoPlasticity
     *              (epsp,xi,eps,stre,cc,matdata,err)

      implicit none

      integer          err,
     *                 i, j

      double precision epsp(*), xi, eps(*), stre(*), cc(6,*), 
     *                 matData(*),
     *                 Dg, nstrialD, nsD, Phitrial, treps,
     *                 mu, K, sY, H,
     *                 fact1, fact2, fact3,
     *                 r0, r1d2, r1, r2, r3, r1d3, r2d3, sq2d3

      data r0, r1d2, r1, r2, r3 / 0.d0, 0.5d0, 1.d0, 2.d0, 3.d0 /

      r1d3 = r1 / r3
      r2d3 = r1d3 + r1d3

      sq2d3 = sqrt(r2d3)

      K  = matData(1)
      mu = matData(2)
      sY = matData(3)
      H  = matData(4)

      !write(*,*) K, mu, sY, H

      if (abs(epsp(1)+epsp(2)+epsp(3)).gt.1.d-10) 
     *  write(*,'(9x,a/)')'WARNING! tr(epsp) <> 0 !'

c...  update for deviatoric stress and internal variables

      !write(*,'(6(1x,g12.5))') (epsp(i),i=1,6)
      !write(*,'(6(1x,g12.5))') (eps(i),i=1,6)

      fact1 = mu + mu
      treps = eps(1) + eps(2) + eps(3)
      fact2 = r1d3 * treps

      do i=1, 3
        stre(i)   = fact1 * (eps(i)  -epsp(i)  - fact2)
        stre(3+i) = fact1 * (eps(3+i)-epsp(3+i))
      enddo 
      
      nstrialD = r0
      do i=1, 3
        nstrialD = nstrialD + stre(i)  *stre(i)
     *                 + r2 * stre(3+i)*stre(3+i)
      enddo
      nstrialD = sqrt(nstrialD)

      Phitrial = nstrialD - sq2d3 * (sY + xi*H)

      if (Phitrial.gt.r0) then

        !write(*,*) Phitrial

        Dg = Phitrial / (r2d3 * H + fact1)

        !write(*,*) 'Dg = ', Dg

        nsD   = nstrialD - fact1 * Dg 
        fact3 = nsD / nstrialD         ! beta
        do i=1, 6
          stre(i) = fact3 * stre(i)
        enddo

        fact2 = Dg / nsD
        do i=1, 6
          epsp(i) = epsp(i) + fact2 * stre(i)
        enddo

        xi = xi + sq2d3 * Dg

        fact2 = fact1 * fact3   ! 2 mu beta

      else

        fact2 = fact1

      endif

c...  computation of consistent tantgent tensor
 
      fact3 = K - r1d3 * fact2
      do i=1, 3
        do j=1, 3
          cc(i,j)     = fact3
          cc(3+i,j)   = r0
          cc(i,3+j)   = r0
          cc(3+i,3+j) = r0
        enddo
        cc(i,i)     = cc(i,i) + fact2
        cc(3+i,3+i) = fact2 * r1d2
      enddo

      if (Phitrial.gt.r0) then

        fact3 = (fact1/(r1+H/(mu+fact1)) - fact1+fact2) / (nsD*nsD) ! 2 mu gam
        do i=1, 6
          do j=1, 6
            cc(i,j) = cc(i,j) - fact3 * stre(i) * stre(j)
          enddo
        enddo

      endif

c...  volumetric stress

      fact3 = K * treps
      do i=1, 3
        stre(i) = stre(i) + fact3
      enddo

      return

      end

