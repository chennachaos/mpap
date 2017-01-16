

      subroutine myma38ad(N, NE, JOB, TRANSA_, LVALUE, LINDEX, VALUE,
     *                    INDEX, KEEP, CNTL, ICNTL, INFO, RINFO)

      implicit none

      integer          N, NE, JOB, TRANSA_, LVALUE, LINDEX, INDEX, KEEP,
     *                 ICNTL, INFO

      double precision VALUE, CNTL, RINFO

      logical          TRANSA

      TRANSA = .true.
      if (TRANSA_ .eq. 0) TRANSA = .false.

      call ma38ad(N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     *            INDEX, KEEP, CNTL, ICNTL, INFO, RINFO)

      return

      end


