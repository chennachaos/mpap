

      subroutine myma38cd(N, JOB, TRANSC_, LVALUE, LINDEX, VALUE, INDEX,
     *                    KEEP, B, X, W, CNTL, ICNTL, INFO, RINFO)

      implicit none

      integer          N, JOB, TRANSC_, LVALUE, LINDEX, INDEX, KEEP,
     *                 ICNTL, INFO

      double precision VALUE, B, X, W, CNTL, RINFO

      logical          TRANSC

      TRANSC = .true.
      if (TRANSC_ .eq. 0) TRANSC = .false.

      call ma38cd(N, JOB, TRANSC, LVALUE, LINDEX, VALUE, INDEX, 
     *            KEEP, B, X, W, CNTL, ICNTL, INFO, RINFO)

      return

      end


