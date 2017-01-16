c
c         +-----------------------------------------+
c         | MATHEMATISCHE SUBROUTINES UND FUNCTIONS |
c         |                                         |
c         |          ZUR LINEAREN ALGEBRA           |
c         |         ======================          |
c         |                                         |
c         |          Wulf Dettmer, 2000 -           |
c         +-----------------------------------------+
c
c-----------------------------------------------------------
c
c  integer function decompLR_matrix(A,P,n,isw)
c
c  subroutine solve_matrix(A,P,b,x,n)
c
c  subroutine inverse_matrix(A,invA,P,n)
c
c  double precision function det_matrixLR(A,n)
c
c  subroutine mult_matrix_p (S,A,B,l,m,n)
c  subroutine mult_matrix_m (S,A,B,l,m,n)
c  subroutine mult_matrix_f (S,A,B,l,m,n,fact)
c  subroutine mult_matrix_ap(S,A,B,l,m,n)
c  subroutine mult_matrix_am(S,A,B,l,m,n)
c  subroutine mult_matrix_af(S,A,B,l,m,n,fact)
c
c  subroutine specmtx(lam,x,A,K,M,tol,mxi,n)
c
c*****************************************************************************
c*****************************************************************************

c=============================================================================
c    LR-Zerlegung einer quadratischen Matrix A  (Gauss-Elimination)
c
c    isw - switch:
c     -1  pivot line i <- A(i,i)      (kein Pivoting, keine Warnungen)
c      0  pivot line i <- A(i,i)                    (kein Pivoting)
c      1  pivot line i <- max[abs{A(i,j)},i=j,n]    (Spaltenpivot-Wahl)
c      2  pivot line i <- max[abs{A(i,j)}/(sum[abs{A(i,l)},l=j+1,n]),i=j,n]
c
c     P A = L R       P  Permutationsmatrix (vertauscht Zeilen)
c                     L  untere Dreiecksmatrix (alle Diagonalelemente = 1)
c                     R  obere Dreiecksmatrix
c
c    input:  Matrix A
c    output: Integervektor P
c            Matrix A überschrieben mit L und R
c-----------------------------------------------------------------------------
      integer function decompLR_matrix(A,P,n,isw)
      implicit none
      integer i, j, k, l, n, P(n), piv, r, s, isw
      double precision A(n,n), hlp, alpha, beta, r0
      r0 = 0.d0

      do i=1, n 
        P(i) = i
      enddo

      r = 0

      do i=1, n-1 

        if (isw.eq.1) goto 10
        if (isw.eq.2) goto 20
        goto 40

 10     piv = i
        do j=i+1, n
          if (abs(A(j,i)).gt.abs(A(piv,i))) piv = j
        enddo
        goto 30

 20     piv = i
        beta = r0
        do k=i, n
          beta = beta + abs(A(piv,k))
        enddo 
        do j=i+1, n
          alpha = r0
          do k=i, n
            alpha = alpha + abs(A(j,k))
          enddo
          if (abs(A(j,i)/alpha).gt.abs(A(piv,i)/beta)) then
            piv = j
            beta = r0
            do k=i, n
              beta = beta + abs(A(piv,k))
            enddo 
          endif
        enddo

 30     if (piv.ne.i) then

c          do j=piv-1, i, -1
c            do l=1, n
c              hlp = A(j,l)
c              A(j,l) = A(j+1,l)
c              A(j+1,l) = hlp
c            enddo
c            l      = P(j)
c            P(j)   = P(j+1)
c            P(j+1) = l
c          enddo

          do j=1, n
            hlp = A(i,j)
            A(i,j) = A(piv,j)
            A(piv,j) = hlp
          enddo
          l = P(i)
          P(i) = P(piv)
          P(piv) = l

        endif

 40     if (abs(A(i,i)).lt.1.d-15) then

          do j=i+1, n
            A(i+1,j) = A(i+1,j) + A(i,j)
          enddo

          r = r + 1

        else

          do j=i+1, n
            A(j,i) = A(j,i) / A(i,i)
          enddo

          do j=i+1, n
            do l=i+1, n
              A(j,l) = A(j,l) - A(i,l)*A(j,i) 
            enddo
          enddo
        endif

      enddo

      if (abs(A(n,n)).lt.1.d-15) then

        r = r + 1

      endif

      if (r.gt.0) then
        if (isw.gt.0) then
          write(*,2000) r
        elseif (isw.eq.0) then
          write(*,1000) r
        endif
      endif

      decompLR_matrix = r

      return 
1000  format(/,'  ERROR in decompLR_matrix:  reduced diagonal is zero',
     *           ' in ',i3,' equations')
2000  format(/,'  ERROR in decompLR_matrix:  matrix is singular, ',
     *           'rank reduced by ',i3)
      end
c
c=============================================================================
c    Lösung des linearen Gleichungssystems  A x = b
c
c       input:  Matrix A in LR-Zerlegung
c               zur LR-Zerlegung zugehörige Zeilenpermutation P
c               Vektor der rechten Seite b
c       output: Lösungsvektor x
c-----------------------------------------------------------------------------
      subroutine solve_matrix(A,P,b,x,n)
      implicit none

      integer n, P(n), i, j
      double precision A(n,n), b(n), x(n)

      do i=1, n
        x(i) = b(P(i))
      enddo

      do i=1, n
        b(i) = x(i) 
        do j=1, i-1
          b(i) = b(i) - b(j)*A(i,j)
        enddo
      enddo

      do i=n, 1, -1
        x(i) = b(i)
        do j=n, i+1, -1
          x(i) = x(i) - x(j)*A(i,j)
        enddo
        x(i) = x(i) / A(i,i)
      enddo

      return
      end
c
c=============================================================================
c   Berechnung der Inverse einer Matrix A
c
c    input:  Matrix A in LR-Zerlegung
c            zur LR-Zerlegung zugehörige Zeilenpermutation P
c    output: Inverse Matrix invA
c-----------------------------------------------------------------------------
      subroutine inverse_matrix(A,invA,P,n)
      implicit none
      integer i, n, P(n)
      double precision A(n,n), invA(n,n), e(1000), r1
      r1 = 1.d0

      do i=1, n
        call pzero(e,n)
        e(i) = r1
        call solve_matrix(A,P,e,invA(1,i),n)
      enddo

      return
      end
c
c=============================================================================
c   det(A)             n x n Matrix A in LR-Zerlegung
c-----------------------------------------------------------------------------
      double precision function det_matrixLR(A,n)
      implicit none
      integer n, i
      double precision A(n,n), x
      x = 1.d0
      do i=1, n
        x = x * A(i,i)
      enddo
      det_matrixLR = x
      return
      end
c
c=============================================================================
c   Multiplikation zweier Matritzen     S  =  A   B
c                                      lxn   lxm mxn 
c-----------------------------------------------------------------------------
      subroutine mult_matrix_p(S,A,B,l,m,n)
      implicit none
      integer i, j, k, l, m, n
      double precision A(l,m), B(m,n), S(l,n), r0
      r0 = 0.d0
      do i=1, l 
        do j=1, n
          S(i,j) = r0
          do k=1, m
            S(i,j) = S(i,j) + A(i,k)*B(k,j)
          enddo 
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_matrix_m(S,A,B,l,m,n)
      implicit none
      integer i, j, k, l, m, n
      double precision A(l,m), B(m,n), S(l,n), r0
      r0 = 0.d0
      do i=1, l 
        do j=1, n
          S(i,j) = r0
          do k=1, m
            S(i,j) = S(i,j) - A(i,k)*B(k,j)
          enddo 
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_matrix_f(S,A,B,l,m,n,fact)
      implicit none
      integer i, j, k, l, m, n
      double precision A(l,m), B(m,n), S(l,n), r0, fact
      r0 = 0.d0
      do i=1, l 
        do j=1, n
          S(i,j) = 0.d0
          do k=1, m
            S(i,j) = S(i,j) + fact * A(i,k)*B(k,j)
          enddo 
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_matrix_ap(S,A,B,l,m,n)
      implicit none
      integer i, j, k, l, m, n
      double precision A(l,m), B(m,n), S(l,n)
      do i=1, l 
        do j=1, n
          do k=1, m
            S(i,j) = S(i,j) + A(i,k)*B(k,j)
          enddo 
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_matrix_am(S,A,B,l,m,n)
      implicit none
      integer i, j, k, l, m, n
      double precision A(l,m), B(m,n), S(l,n)
      do i=1, l 
        do j=1, n
          do k=1, m
            S(i,j) = S(i,j) - A(i,k)*B(k,j)
          enddo 
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_matrix_af(S,A,B,l,m,n,fact)
      implicit none
      integer i, j, k, l, m, n
      double precision A(l,m), B(m,n), S(l,n), fact
      do i=1, l 
        do j=1, n
          do k=1, m
            S(i,j) = S(i,j) + fact * A(i,k)*B(k,j)
          enddo 
        enddo
      enddo
      return
      end
c
c=============================================================================
c     compute Eigenvalues of a symmetric matrix with Jacobi method 
c-----------------------------------------------------------------------------
      subroutine specmtx(lam,x,A,K,M,tol,mxi,n)
      implicit none
      integer          mxi, n, i, j, l, iter
      double precision lam(*), A(n,*), x(n,*), K(n,*), M(n,*), tol,
     *                 r, hlp, Kii, Kjj, xx, yy, alp, gam,
     *                 r0, r1d2, r1

      logical NaN

      data r0, r1d2, r1 / 0.d0, .5d0, 1.d0 /
      
      i = n * n
      call pmove(A,K,i)
      call pzero(M,i)
      call pzero(x,i)
      do i=1, n
        M(i,i) = r1
        x(i,i) = r1
      enddo
        
      r = r0
      do i=1, n-1
        do j=i+1, n
          r = r + K(i,j)*K(i,j) + K(j,i)*K(j,i)
        enddo
      enddo
      r = sqrt(r)

      iter = 0

      do while (r.gt.tol.and.iter.lt.mxi)

        iter = iter + 1

        do i=1, n-1

          do j=i+1, n

            if (abs(K(i,j)).le.0.1d0*tol) goto 999


c            Kii = K(i,i) * M(i,j) - K(i,j) * M(i,i)
c            Kjj = K(j,j) * M(i,j) - K(i,j) * M(j,j)

            Kii = - K(i,j) * M(i,i)
            Kjj = - K(i,j) * M(j,j)

            yy  = r1d2 * (K(i,i) * M(j,j) - K(j,j) * M(i,i))
            xx  = yy + sign(sqrt(yy*yy + Kii * Kjj),yy)
            gam = - Kii / xx
            alp =   Kjj / xx

            do l=1, n
              hlp    = x(i,l) + x(j,l) * gam
              x(j,l) = x(i,l) * alp + x(j,l)
              x(i,l) = hlp
            enddo

            do l=1, n
              hlp    = K(l,i) + K(l,j) * gam
              K(l,j) = K(l,i) * alp + K(l,j)
              K(l,i) = hlp
            enddo

            do l=1, n
              hlp    = K(i,l) + K(j,l) * gam
              K(j,l) = K(i,l) * alp + K(j,l)
              K(i,l) = hlp
            enddo

c            do l=1, n
c              hlp    = M(l,i) + M(l,j) * gam
c              M(l,j) = M(l,i) * alp + M(l,j)
c              M(l,i) = hlp
c            enddo

            M(i,j) = M(i,i) * alp

            M(j,i) = M(j,j) * gam

c            do l=1, n
c              hlp    = M(i,l) + M(j,l) * gam
c              M(j,l) = M(i,l) * alp + M(j,l)
c              M(i,l) = hlp
c            enddo

            hlp    = M(i,i) + M(j,i) * gam
            M(j,i) = M(i,i) * alp + M(j,i)
            M(i,i) = hlp

            hlp    = M(i,j) + M(j,j) * gam
            M(j,j) = M(i,j) * alp + M(j,j)
            M(i,j) = hlp

 999        continue

          enddo

        enddo

        r = r0
        do i=1, n-1
          do j=i+1, n
            r = r + K(i,j)*K(i,j) + K(j,i)*K(j,i)
c    *            + M(i,j)*M(i,j) + M(j,i)*M(j,i)
          enddo
        enddo
        r = sqrt(r)

      enddo

      if (r.gt.tol) then

        write(*,'(/9x,a,g12.5/)') 'Warning! Jacobi iteration: r = ', r

      else

        do i=1, n
          lam(i) = K(i,i) / M(i,i)
        enddo
        
        do i=1, n
          hlp = r0
          do j=1, n
            hlp = hlp + x(j,i)*x(j,i)
          enddo
          hlp = r1 / sqrt(hlp)
          do j=1, n
            x(j,i) = x(j,i) * hlp
          enddo
        enddo

      endif


      do i=1, n
        if (NaN(lam(i))) then
          write(*,*) 'Error! specmtx! Jacobi Iteration!'
          write(*,*) 'lam ',i,lam(i)
          write(*,*) K(i,i), M(i,i)
          do l=1, n
            write(*,'(10(1x,g12.5))') (A(l,j),j=1,n)
          enddo
          stop
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------------





