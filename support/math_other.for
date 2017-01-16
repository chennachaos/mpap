c
c
c   double precision function myatan(y,x)
c   double precision function dmyatandx(y,x)
c   double precision function dmyatandy(y,x)
c
c   double precision function polyarea(x,n)
c   subroutine polyprop(xs,x,n)
c
c   integer function round(x)
c   integer function trunc(x)
c
c   integer function Hs(x)
c
c   double precision function LP(i,m,s)
c   double precision function dLP(i,m,s)
c   double precision function ddLP(i,m,s)
c
c   double precision function dot(a,b,n)
c
c   logical function NaN(x)
c
c   double precision function value(ad,al,au,jpl,jpu,i,j)
c
c   subroutine eliminate_column(A,l,m,n)
c   subroutine eliminate_row(A,l,m,n)
c
c
c
c
c
c
c
c
c
c=======================================================
c (A.9)  myatan(y,x)   takes into account the signs of x and y 
c-------------------------------------------------------
      double precision function myatan(y,x)
      implicit none
      double precision x, y, r0, pi, phi
      r0 = 0.d0
      pi = acos(-1.d0)
      phi = atan(y/x)
      if (x.lt.r0) then
        myatan = pi + phi
      elseif ((x.gt.r0).and.(y.lt.r0)) then
        myatan = 2.d0 * pi + phi
      else
        myatan = phi
      endif
      return
      end
c
c=======================================================
c (A.10)  dmyatandx(y,x)   takes into account the signs of x and y 
c-------------------------------------------------------
      double precision function dmyatandx(y,x)
      implicit none
      double precision x, y
      dmyatandx = - y / (x*x + y*y)
      return
      end
c
c=======================================================
c (A.11)  dmyatandy(y,x)   takes into account the signs of x and y 
c-------------------------------------------------------
      double precision function dmyatandy(y,x)
      implicit none
      double precision x, y
      dmyatandy = x / (x*x + y*y)
      return
      end
c
c=======================================================
c polyarea(x,n)
c-------------------------------------------------------
      double precision function polyarea(x,n)
      implicit none
      integer          n, i
      double precision x(2,*), A, x0(2)
      x0(1) = 0.d0
      x0(2) = 0.d0
      A     = 0.d0
      do i=1, n-1
        A = A + ((x(1,i  )-x0(1))*(x(2,i+1)-x0(2)) 
     *         - (x(1,i+1)-x0(1))*(x(2,i  )-x0(2))) * 0.5d0
      enddo
      polyarea = A + ((x(1,n)-x0(1))*(x(2,1)-x0(2)) 
     *              - (x(1,1)-x0(1))*(x(2,n)-x0(2))) * 0.5d0
      return
      end
c
c=======================================================
c polyprop(prop,x,n)
c-------------------------------------------------------
      subroutine polyprop(prop,x,n)
      implicit none
      integer          n, i
      double precision prop(*), x(2,*), x0(2), Dx(2), Ax(2), DA, A
      x0(1) = 0.d0
      x0(2) = 0.d0
      A     = 0.d0
      Ax(1) = 0.d0
      Ax(2) = 0.d0
      do i=1, n-1
        Dx(1) = ((x(1,i) + x(1,i+1)) - 2.d0*x0(1)) / 3.d0
        Dx(2) = ((x(2,i) + x(2,i+1)) - 2.d0*x0(2)) / 3.d0
        DA    = ((x(1,i  )-x0(1))*(x(2,i+1)-x0(2)) 
     *         - (x(1,i+1)-x0(1))*(x(2,i  )-x0(2))) * 0.5d0
        A     = A     + DA
        Ax(1) = Ax(1) + DA * Dx(1)
        Ax(2) = Ax(2) + DA * Dx(2)
      enddo
      Dx(1) = ((x(1,n) + x(1,1)) - 2.d0*x0(1)) / 3.d0
      Dx(2) = ((x(2,n) + x(2,1)) - 2.d0*x0(2)) / 3.d0
      DA    = ((x(1,n)-x0(1))*(x(2,1)-x0(2)) 
     *       - (x(1,1)-x0(1))*(x(2,n)-x0(2))) * 0.5d0
      A     = A     + DA
      Ax(1) = Ax(1) + DA * Dx(1)
      Ax(2) = Ax(2) + DA * Dx(2)
      prop(1) = Ax(1) / A  ! centre of gravity, x-coordinate 
      prop(2) = Ax(2) / A  ! centre of gravity, y-coordinate 
      prop(3) = A          ! area

      Ax(1) = 0.d0
      Ax(2) = 0.d0
      do i=1, n-1
        Dx(1) = ((x(1,i) + x(1,i+1)) - 2.d0*prop(1)) / 3.d0
        Dx(2) = ((x(2,i) + x(2,i+1)) - 2.d0*prop(2)) / 3.d0
        DA    = ((x(1,i  )-prop(1))*(x(2,i+1)-prop(2)) 
     *         - (x(1,i+1)-prop(1))*(x(2,i  )-prop(2))) * 0.5d0
        Ax(1) = Ax(1) + DA * Dx(1) * Dx(1)
        Ax(2) = Ax(2) + DA * Dx(2) * Dx(2)
      enddo
      Dx(1) = ((x(1,n) + x(1,1)) - 2.d0*prop(1)) / 3.d0
      Dx(2) = ((x(2,n) + x(2,1)) - 2.d0*prop(2)) / 3.d0
      DA    = ((x(1,n)-prop(1))*(x(2,1)-prop(2)) 
     *       - (x(1,1)-prop(1))*(x(2,n)-prop(2))) * 0.5d0
      A     = A     + DA
      Ax(1) = Ax(1) + DA * Dx(1) * Dx(1)
      Ax(2) = Ax(2) + DA * Dx(2) * Dx(2)
      prop(4) = Ax(1)         ! I_x 
      prop(5) = Ax(2)         ! I_y 
      prop(6) = Ax(1) + Ax(2) ! I_p
      return
      end
c
c=======================================================
c (A.12) round(x)    
c-------------------------------------------------------
      integer function round(x)
      implicit none
      double precision x
      integer          h
      h     = int(2*abs(x))
      h     = h + mod(h,2)
      round = sign(.5d0,x) * h 
      return
      end
c=======================================================
c (A.12) trunc(x)    
c-------------------------------------------------------
      integer function trunc(x)
      implicit none
      double precision x
      integer          h
      h = int(x)
      if (x.lt.h) then
        trunc = h - 1
      else
        trunc = h
      endif 
      return
      end
cc=======================================================
c (A.13) Hs(x)    
c-------------------------------------------------------
      integer function Hs(x)
      implicit none
      double precision x
      if (x.gt.1.d-15) then
        Hs = 1
      elseif (x.lt.-1.d-15) then
        Hs = -1
      else
        Hs = 0
      endif      
      return
      end
c
c=======================================================
c (A.15) LP(i,m,s)   
c
c    evaluate Lagrangian Polynomial i
c 
c    0 <= x <= 1,  
c    m equally distributed interpolation points
c
c    example: m = 3 -> x1=0.0, x2=0.5, x3=1.0
c-------------------------------------------------------
      double precision function LP(i,m,s)
      implicit none
      integer          i, m, j
      double precision ds, s, si, sj, N 
      ds = 1.d0 / real(m-1)
      N  = 1.d0
      si = (i-1) * ds
      do j=1, m
        sj = (j-1) * ds
        if (j.ne.i) N = N * (s-sj) / (si-sj)
      enddo
      LP = N
      return
      end
c
c=======================================================
c (A.16) dLP(i,m,s)   
c
c    evaluate first derivative of Lagrangian Polynomial i
c 
c    0 <= x <= 1,  
c    m equally distributed interpolation points
c
c    example: m = 3 -> x1=0.0, x2=0.5, x3=1.0
c-------------------------------------------------------
      double precision function dLP(i,m,s)
      implicit none
      integer          i, m, j, k
      double precision ds, s, si, sj, sk, dN, h
      save
      ds = 1.d0 / real(m-1)
      dN = 0.d0
      si = (i-1) * ds
      do k=1, m
        if (k.ne.i) then
          sk = (k-1) * ds
          h = 1.d0 / (si-sk)
          do j=1, m
            sj = (j-1) * ds
            if (j.ne.i.and.j.ne.k) h = h *  (s-sj) / (si-sj)       
          enddo
          dN = dN + h
        endif
      enddo
      dLP = dN
      return
      end
c
c=======================================================
c (A.17) ddLP(i,m,s)   
c
c    evaluate second derivative of Lagrangian Polynomial i
c 
c    0 <= x <= 1,  
c    m equally distributed interpolation points
c
c    example: m = 3 -> x1=0.0, x2=0.5, x3=1.0
c-------------------------------------------------------
      double precision function ddLP(i,m,s)
      implicit none
      integer          i, m, j, k, l
      double precision ds, s, si, sj, sk, sl, ddN, h
      save
      ds = 1.d0 / real(m-1)
      ddN = 0.d0
      si = (i-1) * ds
      do l=1, m
        if (l.ne.i) then
          sl = (l-1) * ds
          do k=1, m
            if (k.ne.i.and.k.ne.l) then
              sk = (k-1) * ds
              h = 1.d0 / ((si-sk) * (si-sl))
              do j=1, m
                sj = (j-1) * ds
                if (j.ne.i.and.j.ne.k.and.j.ne.l) h = h * (s-sj)/(si-sj)
              enddo
              ddN = ddN + h
            endif
          enddo
        endif
      enddo
      ddLP = ddN
      return
      end
c
c=======================================================
c (A.18) dot(a,b,n)    
c
c    dot product of two vectors of dimension n 
c
c-------------------------------------------------------
      double precision function dot(a,b,n)
      implicit none
      integer n, i
      double precision a(*), b(*), d
      save
      d = 0.d0
      do i=1, n
        d = d + a(i) * b(i)
      enddo
      dot = d
      return
      end
c
c=======================================================
c (A.18) nan(x)    
c
c     check whether x is a well defined number or not
c
c-------------------------------------------------------
      logical function NaN(x)
      implicit none
      double precision x
      character strg*3
      write(strg,'(f3.1)') x
c      write(*,*) strg
      NaN = .false.
      if (strg(1:1).eq.'N'.or.strg(1:1).eq.'n'.or.strg(1:1).eq.'?'.or.
     *    strg(1:1).eq.' ') 
     *  NaN = .true.
      return
      end
c
c========================================================================
c
c (D.19)   return value of coefficient A[i,j] in profile matrix
c
c------------------------------------------------------------------------
      double precision function value(ad,al,au,jpl,jpu,i,j)
      implicit none
      integer i, j, jpl(*), jpu(*), n
      double precision ad(*), al(*), au(*)
      if (i.gt.j) then
        if (jpl(i)-jpl(i-1).ge.i-j) then
          value = al(jpl(i)-i+j+1)
        else
          value = 0.d0
        endif
      elseif (i.eq.j) then
        value = ad(i)
      else
        if (jpu(j)-jpu(j-1).ge.j-i) then
          value = au(jpu(j)-j+i+1)
        else
          value = 0.d0
        endif
      endif
      return
      end
c
c=======================================================
c (D.20)   Eliminiere die l-te Spalte aus einer Matrix
c           Matrixdimension:   m x n  ->  m x n-1
c-------------------------------------------------------
      subroutine eliminate_column(A,l,m,n)
      implicit none
      integer i, l, m, n
      double precision A(m*n)
      do i=(l-1)*m+1, l*m
        A(i) = A(i+m)
      enddo
      return
      end
c
c=======================================================
c (D.12)   Eliminiere die l-te Zeile aus einer Matrix
c           Matrixdimension:   m x n  ->  m-1 x n
c-------------------------------------------------------
      subroutine eliminate_row(A,l,m,n)
      implicit none
      integer i, j, l, m, n
      double precision A(m*n)
      do i=1, n
        do j=1, m-1
          A((i-1)*(m-1)+l + j-1) = A((i-1)*(m-1)+l + j-1+i)
        enddo
      enddo
      return
      end
c
c=======================================================

