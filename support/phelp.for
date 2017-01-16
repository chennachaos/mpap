
      subroutine pzeroi(a,n)
      implicit none
      integer n, a(*), i
      do i=1, n
        a(i) = 0
      enddo
      return
      end


      subroutine pzero(a,n)
      implicit none
      integer          n, i
      double precision a(*), r0
      r0 = 0.d0
      do i=1, n
        a(i) = r0
      enddo
      return
      end


      subroutine pmovei(a,b,n)
      implicit none
      integer a(*), b(*), n, i
      do i=1, n
        b(i) = a(i)
      enddo
      return
      end


      subroutine pmove(a,b,n)
      implicit none
      integer          n, i
      double precision a(*), b(*)
      do i=1, n
        b(i) = a(i)
      enddo
      return
      end


      subroutine paddi(a,b,n)
      implicit none
      integer a(*), b(*), n, i
      do i=1, n
        b(i) = b(i) + a(i)
      enddo
      return
      end


      subroutine padd(a,b,n)
      implicit none
      integer          n, i
      double precision a(*), b(*)
      do i=1, n
        b(i) = b(i) + a(i)
      enddo
      return
      end


      subroutine psubi(a,b,n)
      implicit none
      integer a(*), b(*), n, i
      do i=1, n
        b(i) = b(i) - a(i)
      enddo
      return
      end


      subroutine psub(a,b,n)
      implicit none
      integer          n, i
      double precision a(*), b(*)
      do i=1, n
        b(i) = b(i) - a(i)
      enddo
      return
      end








