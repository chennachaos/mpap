
      subroutine testTangent(matData,F,stre,cc,finite)

      IMPLICIT NONE

      integer          finite

      double precision matData(*), F(3,*), stre(*), cc(6,*)

      logical spectralDecomposition


!      call unit_s_p(stre)
!  
!      do i=1, 6
!        c(i) = dlamb2(j,i)
!      enddo
!
!      call mult_4_ss_p(cc,stre,c)
!
!      call set_s_f(stre,stre,lamb2(j))


c      call sq_s(stre,cc,b)
c      !call set_4_f(cc,cc,1.d0)

      return

      end





