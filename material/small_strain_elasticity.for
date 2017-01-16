

      subroutine small_strain_elasticity(matData,F,stre,cc)

      implicit none

      integer          i, j

      double precision matData(*), F(3,*), stre(*), cc(6,*),
     *                 K, mu, lamb, eps(6), fact, fact2,
     *                 r1d2, r2d3, r1, r2

      data r1d2, r1, r2 / 0.5d0, 1.d0, 2.d0 /

      r2d3 = r2 / 3.d0

      K      = matData(1)
      mu     = matData(2)
      lamb   = K - r2d3 * mu

c     strain tensor

      eps(1) = F(1,1) - r1
      eps(2) = F(2,2) - r1
      eps(3) = F(3,3) - r1
      eps(4) = r1d2 * (F(1,2) + F(2,1))
      eps(5) = r1d2 * (F(2,3) + F(3,2))
      eps(6) = r1d2 * (F(3,1) + F(1,3))

c      do i=1, 6
c        write(*,'(A, F20.18)') 'eps=', eps(i)
c      enddo

c     stress

      fact   = lamb * (eps(1) + eps(2) + eps(3))
      fact2  = r2 * mu
      
c      write(*,'(A, F20.18)') 'fact=', fact
c      write(*,'(A, F12.10)') 'fact=', fact2

      stre(1) = fact + fact2 * eps(1)
      stre(2) = fact + fact2 * eps(2)
      stre(3) = fact + fact2 * eps(3)
      stre(4) =        fact2 * eps(4)
      stre(5) =        fact2 * eps(5)
      stre(6) =        fact2 * eps(6)

c      do i=1, 6
c        write(*,'(A, F20.18)') 'stre=', stre(i)
c      enddo

c     tangent tensor

      call pzero(cc,36)

      do i=1, 3
        do j=1, 3
          cc(i,j) = lamb
        enddo
        cc(i,i)     = cc(i,i) + r2 * mu
        cc(i+3,i+3) = mu
      enddo

      return

      end






