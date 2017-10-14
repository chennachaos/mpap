c
c         +-----------------------------------------+
c         | MATHEMATISCHE SUBROUTINES UND FUNCTIONS |
c         | mathematical subroutines  and functions |
c         |                                         |
c         |           ZUR TENSORANALYSIS            |
c         |           for tensor analysis           |
c         |          ====================           |
c         |                                         |
c         |       Wulf G. Dettmer, 2000 - today     |
c         +-----------------------------------------+
c
c-----------------------------------------------------------
c  OPERATIONEN, DIE EINEN SKALAR LIEFERN
c  operations which render a scalar
c
c  double precision function det_s(A)
c  double precision function det_u(A)
c
c  double precision function traceAB_ss(A,B)
c  double precision function traceAB_su(A,B)
c  double precision function traceAB_uu(A,B)
c
c  double precision function dot_ss(A,B)
c  double precision function dot_su(A,B)
c  double precision function dot_uu(A,B)
c
c  double precision function norm2_u (A)
c  double precision function norm2_s (A)
c  double precision function norm2s_u(A)
c  double precision function normx_s (A)
c  double precision function normx_u (A)
c      
c-----------------------------------------------------------
c  OPERATIONEN, DIE EINEN ZWEISTUFIGEN TENSOR LIEFERN
c  operations that render a second order tensor
c
c  subroutine unit_s_p (S)
c  subroutine unit_s_m (S)
c  subroutine unit_s_f (S,fact)
c  subroutine unit_s_ap(S)
c  subroutine unit_s_am(S)
c  subroutine unit_s_af(S,fact)
c
c  subroutine unit_u_p (S)
c  subroutine unit_u_m (S)
c  subroutine unit_u_f (S)
c  subroutine unit_u_ap(S)
c  subroutine unit_u_am(S)
c  subroutine unit_u_af(S)
c
c  subroutine set_s_p (S,A)
c  subroutine set_s_m (S,A)
c  subroutine set_s_f (S,A,fact)
c  subroutine set_s_ap(S,A)
c  subroutine set_s_am(S,A)
c  subroutine set_s_af(S,A,fact)
c
c  subroutine set_u_p (S,A)
c  subroutine set_u_m (S,A)
c  subroutine set_u_f (S,A,fact)
c  subroutine set_u_ap(S,A)
c  subroutine set_u_am(S,A)
c  subroutine set_u_af(S,A,fact)
c
c  subroutine set_u_s_p (S,A)
c  subroutine set_u_s_m (S,A)
c  subroutine set_u_s_f (S,A,fact)
c  subroutine set_u_s_ap(S,A)
c  subroutine set_u_s_am(S,A)
c  subroutine set_u_s_af(S,A,fact)
c
c  subroutine sympart_p (S,A)
c  subroutine sympart_m (S,A)
c  subroutine sympart_f (S,A,fact)
c  subroutine sympart_ap(S,A)
c  subroutine sympart_am(S,A)
c  subroutine sympart_af(S,A,fact)
c
c  subroutine sympart_AU_p (S,A,B)
c  subroutine sympart_AU_m (S,A,B)
c  subroutine sympart_AU_f (S,A,B,fact)
c  subroutine sympart_AU_ap(S,A,B)
c  subroutine sympart_AU_am(S,A,B)
c  subroutine sympart_AU_af(S,A,B,fact)
c
c  subroutine sympart_UA_p (S,A,B)
c  subroutine sympart_UA_m (S,A,B)
c  subroutine sympart_UA_f (S,A,B,fact)
c  subroutine sympart_UA_ap(S,A,B)
c  subroutine sympart_UA_am(S,A,B)
c  subroutine sympart_UA_af(S,A,B,fact)
c
c  subroutine mult_s_ss_p (S,A,B)
c  subroutine mult_s_ss_m (S,A,B)
c  subroutine mult_s_ss_f (S,A,B,fact)
c  subroutine mult_s_ss_ap(S,A,B)
c  subroutine mult_s_ss_am(S,A,B)
c  subroutine mult_s_ss_af(S,A,B)
c
c  subroutine mult_s_su_p (S,A,B)
c  subroutine mult_s_su_m (S,A,B)
c  subroutine mult_s_su_f (S,A,B,fact)
c  subroutine mult_s_su_ap(S,A,B)
c  subroutine mult_s_su_am(S,A,B)
c  subroutine mult_s_su_af(S,A,B,fact)
c    
c  subroutine mult_s_us_p (S,A,B)
c  subroutine mult_s_us_m (S,A,B)
c  subroutine mult_s_us_f (S,A,B,fact)
c  subroutine mult_s_us_ap(S,A,B)
c  subroutine mult_s_us_am(S,A,B)
c  subroutine mult_s_us_af(S,A,B,fact)
c
c  subroutine mult_s_uu_p (S,A,B)
c  subroutine mult_s_uu_m (S,A,B)
c  subroutine mult_s_uu_f (S,A,B,fact)
c  subroutine mult_s_uu_ap(S,A,B)
c  subroutine mult_s_uu_am(S,A,B)
c  subroutine mult_s_uu_af(S,A,B,fact)
c
c  subroutine mult_s_uuT_p (S,A,B)
c  subroutine mult_s_uuT_m (S,A,B)
c  subroutine mult_s_uuT_f (S,A,B,fact)
c  subroutine mult_s_uuT_ap(S,A,B)
c  subroutine mult_s_uuT_am(S,A,B)
c  subroutine mult_s_uuT_af(S,A,B,fact)
c
c  subroutine mult_s_uTu_p (S,A,B)
c  subroutine mult_s_uTu_m (S,A,B)
c  subroutine mult_s_uTu_f (S,A,B,fact)
c  subroutine mult_s_uTu_ap(S,A,B)
c  subroutine mult_s_uTu_am(S,A,B)
c  subroutine mult_s_uTu_af(S,A,B,fact)
c
c  subroutine mult_sd_ss_p (S,A,B)
c  subroutine mult_sd_ss_m (S,A,B)
c  subroutine mult_sd_ss_f (S,A,B,fact)
c  subroutine mult_sd_ss_ap(S,A,B)
c  subroutine mult_sd_ss_am(S,A,B)
c  subroutine mult_sd_ss_af(S,A,B,fact)
c
c  subroutine mult_sd_su_p (S,A,B)
c  subroutine mult_sd_su_m (S,A,B)
c  subroutine mult_sd_su_f (S,A,B,fact)
c  subroutine mult_sd_su_ap(S,A,B)
c  subroutine mult_sd_su_am(S,A,B)
c  subroutine mult_sd_su_af(S,A,B,fact)
c    
c  subroutine mult_sd_us_p (S,A,B)
c  subroutine mult_sd_us_m (S,A,B)
c  subroutine mult_sd_us_f (S,A,B,fact)
c  subroutine mult_sd_us_ap(S,A,B)
c  subroutine mult_sd_us_am(S,A,B)
c  subroutine mult_sd_us_af(S,A,B,fact)
c
c  subroutine mult_sd_uu_p (S,A,B)
c  subroutine mult_sd_uu_m (S,A,B)
c  subroutine mult_sd_uu_f (S,A,B,fact)
c  subroutine mult_sd_uu_ap(S,A,B)
c  subroutine mult_sd_uu_am(S,A,B)
c  subroutine mult_sd_uu_af(S,A,B,fact)
c
c  subroutine mult_u_ss_p (S,A,B)
c  subroutine mult_u_ss_m (S,A,B)
c  subroutine mult_u_ss_f (S,A,B,fact)
c  subroutine mult_u_ss_ap(S,A,B)
c  subroutine mult_u_ss_am(S,A,B)
c  subroutine mult_u_ss_af(S,A,B)
c
c  subroutine mult_u_su_p (S,A,B)
c  subroutine mult_u_su_m (S,A,B)
c  subroutine mult_u_su_f (S,A,B,fact)
c  subroutine mult_u_su_ap(S,A,B)
c  subroutine mult_u_su_am(S,A,B)
c  subroutine mult_u_su_af(S,A,B,fact)
c    
c  subroutine mult_u_us_p (S,A,B)
c  subroutine mult_u_us_m (S,A,B)
c  subroutine mult_u_us_f (S,A,B,fact)
c  subroutine mult_u_us_ap(S,A,B)
c  subroutine mult_u_us_am(S,A,B)
c  subroutine mult_u_us_af(S,A,B,fact)
c
c  subroutine mult_u_uu_p (S,A,B)
c  subroutine mult_u_uu_m (S,A,B)
c  subroutine mult_u_uu_f (S,A,B,fact)
c  subroutine mult_u_uu_ap(S,A,B)
c  subroutine mult_u_uu_am(S,A,B)
c  subroutine mult_u_uu_af(S,A,B,fact)
c
c  subroutine mult_s_s4_p (S,A,X)
c  subroutine mult_s_s4_m (S,A,X)
c  subroutine mult_s_s4_f (S,A,X,fact)
c  subroutine mult_s_s4_ap(S,A,X)
c  subroutine mult_s_s4_am(S,A,X)
c  subroutine mult_s_s4_af(S,A,X,fact)
c
c  subroutine mult_s_4s_p (S,X,A)
c  subroutine mult_s_4s_m (S,X,A)
c  subroutine mult_s_4s_f (S,X,A,fact)
c  subroutine mult_s_4s_ap(S,X,A)
c  subroutine mult_s_4s_am(S,X,A)
c  subroutine mult_s_4s_af(S,X,A,fact)
c
c  subroutine power_s(S,A,n)
c  subroutine power_u(S,A,n)
c
c  subroutine inverse_s(S,A)
c  subroutine inverse_u(S,A)
c
c  subroutine transpose(S,A)
c  subroutine transpose_r(S)
c
c  subroutine dev_s_r(S)
c  subroutine dev_u_r(S)
c
c  subroutine diff_norm2sATd_p (S,A,AT,n2s)
c  subroutine diff_norm2sATd_m (S,A,AT,n2s)
c  subroutine diff_norm2sATd_f (S,A,AT,n2s,fact)
c  subroutine diff_norm2sATd_ap(S,A,AT,n2s)
c  subroutine diff_norm2sATd_am(S,A,AT,n2s)
c  subroutine diff_norm2sATd_af(S,A,AT,n2s,fact)
c
c  subroutine diff_norm2sATpUd_p (S,A,B,AT,n2s)
c  subroutine diff_norm2sATpUd_m (S,A,B,AT,n2s)
c  subroutine diff_norm2sATpUd_f (S,A,B,AT,n2s,fact)
c  subroutine diff_norm2sATpUd_ap(S,A,B,AT,n2s)
c  subroutine diff_norm2sATpUd_am(S,A,B,AT,n2s)
c  subroutine diff_norm2sATpUd_af(S,A,B,AT,n2s,fact)
c
c  subroutine diff_norm2sTApUd_p (S,A,B,AT,n2s)
c  subroutine diff_norm2sTApUd_m (S,A,B,AT,n2s)
c  subroutine diff_norm2sTApUd_f (S,A,B,AT,n2s,fact)
c  subroutine diff_norm2sTApUd_ap(S,A,B,AT,n2s)
c  subroutine diff_norm2sTApUd_am(S,A,B,AT,n2s)
c  subroutine diff_norm2sTApUd_af(S,A,B,AT,n2s,fact)
c
c--------------------------------------------------------------------
c  OPERATIONEN, DIE EINEN VIERSTUFIGEN TENSOR (6x6) LIEFERN
c  operations that render a fourth order tensor (6x6 matrix)
c
c  subroutine set_4_p (X,A)
c  subroutine set_4_m (X,A)
c  subroutine set_4_f (X,A,fact)
c  subroutine set_4_ap(X,A)
c  subroutine set_4_am(X,A)
c  subroutine set_4_af(X,A,fact)
c
c  subroutine unit_4_p (X)
c  subroutine unit_4_m (X)
c  subroutine unit_4_f (X,fact)
c  subroutine unit_4_ap(X)
c  subroutine unit_4_am(X)
c  subroutine unit_4_af(X,fact)
c
c  subroutine mult_4_ss_p (X,A,B)
c  subroutine mult_4_ss_m (X,A,B)
c  subroutine mult_4_ss_f (X,A,B,fact)
c  subroutine mult_4_ss_ap(X,A,B)
c  subroutine mult_4_ss_am(X,A,B)
c  subroutine mult_4_ss_af(X,A,B,fact)
c
c  subroutine mult_4_44_p (X,A,B)
c  subroutine mult_4_44_m (X,A,B)
c  subroutine mult_4_44_f (X,A,B,fact)
c  subroutine mult_4_44_ap(X,A,B)
c  subroutine mult_4_44_am(X,A,B)
c  subroutine mult_4_44_af(X,A,B,fact)
c
c  subroutine mult_4_4s_p (X,A,B)
c  subroutine mult_4_4s_m (X,A,B)
c  subroutine mult_4_4s_f (X,A,B,fact)
c  subroutine mult_4_4s_ap(X,A,B)
c  subroutine mult_4_4s_am(X,A,B)
c  subroutine mult_4_4s_af(X,A,B,fact)
c
c  subroutine diff_AT_p (X,A)
c  subroutine diff_AT_m (X,A)
c  subroutine diff_AT_f (X,A,fact)
c  subroutine diff_AT_ap(X,A)
c  subroutine diff_AT_am(X,A)
c  subroutine diff_AT_af(X,A,fact)
c
c  subroutine diff_UT_p (X,A)
c  subroutine diff_UT_m (X,A)
c  subroutine diff_UT_f (X,A,fact)
c  subroutine diff_UT_ap(X,A)
c  subroutine diff_UT_am(X,A)
c  subroutine diff_UT_af(X,A,fact)
c
c  subroutine diff_TU_p (X,A)
c  subroutine diff_TU_m (X,A)
c  subroutine diff_TU_f (X,A,fact)
c  subroutine diff_TU_ap(X,A)
c  subroutine diff_TU_am(X,A)
c  subroutine diff_TU_af(X,A,fact)
c
c  subroutine diff_TAT_p (X,AT)
c  subroutine diff_TAT_m (X,AT)
c  subroutine diff_TAT_f (X,AT,fact)
c  subroutine diff_TAT_ap(X,AT)
c  subroutine diff_TAT_am(X,AT)
c  subroutine diff_TAT_af(X,AT,fact)
c
c  subroutine diff_ATA_p (X,A)
c  subroutine diff_ATA_m (X,A)
c  subroutine diff_ATA_f (X,A,fact)
c  subroutine diff_ATA_ap(X,A)
c  subroutine diff_ATA_am(X,A)
c  subroutine diff_ATA_af(X,A,fact)
c
c  subroutine diff_UTUt_p (X,A)
c  subroutine diff_UTUt_m (X,A)
c  subroutine diff_UTUt_f (X,A,fact)
c  subroutine diff_UTUt_ap(X,A)
c  subroutine diff_UTUt_am(X,A)
c  subroutine diff_UTUt_af(X,A,fact)
c
c  subroutine diff_ATB_p (X,A,B)
c  subroutine diff_ATB_m (X,A,B)
c  subroutine diff_ATB_f (X,A,B,fact)
c  subroutine diff_ATB_ap(X,A,B)
c  subroutine diff_ATB_am(X,A,B)
c  subroutine diff_ATB_af(X,A,B,fact)
c
c  subroutine diff_UTV_p (X,A,B)
c  subroutine diff_UTV_m (X,A,B)
c  subroutine diff_UTV_f (X,A,B,fact)
c  subroutine diff_UTV_ap(X,A,B)
c  subroutine diff_UTV_am(X,A,B)
c  subroutine diff_UTV_af(X,A,B,fact)
c
c  subroutine diff_symAT_p (X,A)
c  subroutine diff_symAT_m (X,A)
c  subroutine diff_symAT_f (X,A,fact)
c  subroutine diff_symAT_ap(X,A)
c  subroutine diff_symAT_am(X,A)
c  subroutine diff_symAT_af(X,A,fact)
c
c  subroutine diff_symUT_p (X,A)
c  subroutine diff_symUT_m (X,A)
c  subroutine diff_symUT_f (X,A,fact)
c  subroutine diff_symUT_ap(X,A)
c  subroutine diff_symUT_am(X,A)
c  subroutine diff_symUT_af(X,A,fact)
c
c  subroutine diff_symATB_p (X,A,B)
c  subroutine diff_symATB_m (X,A,B)
c  subroutine diff_symATB_f (X,A,B,fact)
c  subroutine diff_symATB_ap(X,A,B)
c  subroutine diff_symATB_am(X,A,B)
c  subroutine diff_symATB_af(X,A,B,fact)
c
c  subroutine diff_UTA_p (X,A,B)
c  subroutine diff_UTA_m (X,A,B)
c  subroutine diff_UTA_f (X,A,B,fact)
c  subroutine diff_UTA_ap(X,A,B)
c  subroutine diff_UTA_am(X,A,B)
c  subroutine diff_UTA_af(X,A,B,fact)
c
c  subroutine diff_ATdB_p (X,A,B)
c  subroutine diff_ATdB_m (X,A,B)
c  subroutine diff_ATdB_f (X,A,B,fact)
c  subroutine diff_ATdB_ap(X,A,B)
c  subroutine diff_ATdB_am(X,A,B)
c  subroutine diff_ATdB_af(X,A,B,fact)
c
c  subroutine diff_TAdB_p (X,A,B)
c  subroutine diff_TAdB_m (X,A,B)
c  subroutine diff_TAdB_f (X,A,B,fact)
c  subroutine diff_TAdB_ap(X,A,B)
c  subroutine diff_TAdB_am(X,A,B)
c  subroutine diff_TAdB_af(X,A,B,fact)
c
c--------------------------------------------------------------------
c  OPERATIONEN, DIE EINEN VIERSTUFIGEN TENSOR (6x3x3) LIEFERN
c  operations that render a fourth order tensor (6x3x3) matrix
c
c  subroutine diff_TATt_p (X,TA)
c
c
c
c--------------------------------------------------------------------
c  SPECTRAL DECOMPOSITION OF SECOND ORDER TENSOR
c
c  integer function spectralDecomposition(x,E,A)
c
c--------------------------------------------------------------------
c  ISOTROPIC TENSOR FUNCTIONS AND THEIR DERIVATIVES
c
c  subroutine exp_s(S,X,A,n)
c
c  subroutine sq_s(S,X,A)
c
c  subroutine log_s(S,X,A,n)
c
c  subroutine sqrt_s(S,X,A)
c
c  subroutine isotropicTensorFunction(S,C,x,y,dydx,E,A)
c
c
c*****************************************************************************
c*****************************************************************************

c=============================================================================
c Determinante eines Tensors 2.Stufe
c-----------------------------------------------------------------------------
      double precision function det_s(A) 
      implicit none
      double precision A(6)
      det_s =  A(1)*A(2)*A(3) 
     *       + A(4)*A(5)*A(6)
     *       + A(6)*A(4)*A(5)
     *       - A(1)*A(5)*A(5)
     *       - A(4)*A(4)*A(3)
     *       - A(6)*A(2)*A(6)
      return
      end
c-----------------------------------------------------------------------------
      double precision function det_u(A) 
      implicit none
      double precision A(3,3)
      det_u =  A(1,1)*A(2,2)*A(3,3) 
     *       + A(1,2)*A(2,3)*A(3,1)
     *       + A(1,3)*A(2,1)*A(3,2)
     *       - A(1,3)*A(2,2)*A(3,1)
     *       - A(1,1)*A(2,3)*A(3,2)
     *       - A(1,2)*A(2,1)*A(3,3)
      return
      end
c
c=============================================================================
c Spur des Produktes zweier Tensoren     s = tr(AB)
c-----------------------------------------------------------------------------
      double precision function traceAB_ss(A,B)
      implicit none
      double precision A(6), B(6)
      traceAB_ss =  A(1)*B(1) + A(4)*B(4) + A(6)*B(6)
     *            + A(4)*B(4) + A(2)*B(2) + A(5)*B(5)
     *            + A(6)*B(6) + A(5)*B(5) + A(3)*B(3)
      return
      end
c-----------------------------------------------------------------------------
      double precision function traceAB_su(A,B)
      implicit none
      double precision A(6), B(3,3)
      traceAB_su =  A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
     *            + A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
     *            + A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      double precision function traceAB_uu(A,B)
      implicit none
      double precision A(3,3), B(3,3)
      traceAB_uu =  A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
     *            + A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
     *            + A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)
      return
      end
c
c=============================================================================
c Skalarprodukt zweier Tensoren      s = A : B   
c-----------------------------------------------------------------------------
      double precision function dot_ss(A,B)
      implicit none 
      integer i
      double precision A(6), B(6), s
      s = 0.0d0
      do i=1, 3 
        s = s + A(i)*B(i) + 2.d0 * A(i+3)*B(i+3)
      enddo
      dot_ss = s
      return
      end
c-----------------------------------------------------------------------------
      double precision function dot_su(A,B)
      implicit none 
      integer i
      double precision A(6), B(3,3), s
      s = 0.0d0
      do i=1, 3 
        s = s + A(i)*B(i,i)
      enddo
      s = s + A(4) * (B(1,2) + B(2,1))
      s = s + A(5) * (B(2,3) + B(3,2))
      s = s + A(6) * (B(3,1) + B(1,3))
      dot_su = s
      return
      end
c-----------------------------------------------------------------------------
      double precision function dot_uu(A,B)
      implicit none 
      integer i, j
      double precision A(3,3), B(3,3), s
      s = 0.0d0
      do i=1, 3 
        do j=1, 3
          s = s + A(i,j)*B(i,j)
        enddo
      enddo
      dot_uu = s
      return
      end
c
c=============================================================================
c   Euclidian norm of a tensor
c-----------------------------------------------------------------------------
      double precision function norm2_u(A)
      implicit none
      integer          i
      double precision A(9), h
      h = 0.d0
      do i=1, 9
        h = h + A(i) * A(i)
      enddo
      norm2_u = sqrt(h)
      return
      end
c-----------------------------------------------------------------------------
      double precision function norm2_s(A)
      implicit none
      double precision A(6)
      norm2_s = sqrt(A(1)*A(1) + A(2)*A(2) + A(3)*A(3)
     *     + 2.d0 * (A(4)*A(4) + A(5)*A(5) + A(6)*A(6)))
      return
      end
c-----------------------------------------------------------------------------
      double precision function norm2s_u(A)
      implicit none
      double precision A(3,3)
      norm2s_u = sqrt( A(1,1)*A(1,1) + A(2,2)*A(2,2) + A(3,3)*A(3,3) 
     *       + 2.d0 * (A(1,2)*A(2,1) + A(1,3)*A(3,1) + A(3,2)*A(2,3)))
      return
      end
c
c=============================================================================
c   maximum norm of a tensor
c-----------------------------------------------------------------------------
      double precision function normx_s(A)
      implicit none
      integer          i
      double precision A(6), h
      h = abs(A(1))
      do i=2, 6
        h = max(h,abs(A(i)))
      enddo
      normx_s = h
      return
      end
c-----------------------------------------------------------------------------
      double precision function normx_u(A)
      implicit none
      integer          i
      double precision A(9), h
      h = abs(A(1))
      do i=2, 9
        h = max(h,abs(A(i)))
      enddo
      normx_u = h
      return
      end
c
c=============================================================================
c*****************************************************************************
c*****************************************************************************

c=============================================================================
c  Einheitstensor
c-----------------------------------------------------------------------------
      subroutine unit_s_p(S)
      implicit none
      double precision S(6)
      S(1) = 1.d0
      S(2) = S(1)
      S(3) = S(1)
      S(4) = 0.d0
      S(5) = S(4)
      S(6) = S(4)
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_s_m(S)
      implicit none
      double precision S(6)
      S(1) = - 1.d0
      S(2) = S(1)
      S(3) = S(1)
      S(4) = 0.d0
      S(5) = S(4)
      S(6) = S(4)
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_s_f(S,fact)
      implicit none
      double precision S(6), fact
      S(1) = fact
      S(2) = fact
      S(3) = fact
      S(4) = 0.d0
      S(5) = S(4)
      S(6) = S(4)
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_s_ap(S)
      implicit none
      double precision S(6), r1
      r1 = 1.d0
      S(1) = S(1) + r1
      S(2) = S(2) + r1
      S(3) = S(3) + r1
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_s_am(S)
      implicit none
      double precision S(6), r1
      r1 = 1.d0
      S(1) = S(1) - r1
      S(2) = S(2) - r1
      S(3) = S(3) - r1
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_s_af(S,fact)
      implicit none
      double precision S(6), fact
      S(1) = S(1) + fact
      S(2) = S(2) + fact
      S(3) = S(3) + fact
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_u_p(S)
      implicit none
      double precision S(3,3), r1
      r1 = 1.d0
      call pzero(S,9)
      S(1,1) = r1
      S(2,2) = r1
      S(3,3) = r1
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_u_m(S)
      implicit none
      double precision S(3,3), r1
      r1 = - 1.d0
      call pzero(S,9)
      S(1,1) = r1
      S(2,2) = r1
      S(3,3) = r1
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_u_f(S,fact)
      implicit none
      double precision S(3,3), fact
      call pzero(S,9)
      S(1,1) = fact
      S(2,2) = fact
      S(3,3) = fact
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_u_ap(S)
      implicit none
      double precision S(3,3), r1
      r1 = 1.d0
      S(1,1) = S(1,1) + r1
      S(2,2) = S(2,2) + r1
      S(3,3) = S(3,3) + r1
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_u_am(S)
      implicit none
      double precision S(3,3), r1
      r1 = 1.d0
      S(1,1) = S(1,1) - r1
      S(2,2) = S(2,2) - r1
      S(3,3) = S(3,3) - r1
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_u_af(S,fact)
      implicit none
      double precision S(3,3), fact
      S(1,1) = S(1,1) + fact
      S(2,2) = S(2,2) + fact
      S(3,3) = S(3,3) + fact
      return
      end
c
c=============================================================================
c  Zuweisung    S := A
c-----------------------------------------------------------------------------
      subroutine set_s_p(S,A)
      implicit none
      integer i
      double precision S(6), A(6)
      do i=1, 6
        S(i) = A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_s_m(S,A)
      implicit none
      integer i
      double precision S(6), A(6)
      do i=1, 6
        S(i) = - A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_s_f(S,A,fact)
      implicit none
      integer i
      double precision S(6), A(6), fact
      do i=1, 6
        S(i) = fact * A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_s_ap(S,A)
      implicit none
      integer i
      double precision S(6), A(6)
      do i=1, 6
        S(i) = S(i) + A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_s_am(S,A)
      implicit none
      integer i
      double precision S(6), A(6)
      do i=1, 6
        S(i) = S(i) - A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_s_af(S,A,fact)
      implicit none
      integer i
      double precision S(6), A(6), fact
      do i=1, 6
        S(i) = S(i) + fact * A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_p(S,A)
      implicit none
      integer i
      double precision S(9), A(9)
      do i=1, 9
        S(i) = A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_m(S,A)
      implicit none
      integer i
      double precision S(9), A(9)
      do i=1, 9
        S(i) = - A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_f(S,A,fact)
      implicit none
      integer i
      double precision S(9), A(9), fact
      do i=1, 9
        S(i) = fact * A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_ap(S,A)
      implicit none
      integer i
      double precision S(9), A(9)
      do i=1, 9
        S(i) = S(i) + A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_am(S,A)
      implicit none
      integer i
      double precision S(9), A(9)
      do i=1, 9
        S(i) = S(i) - A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_af(S,A,fact)
      implicit none
      integer i
      double precision S(9), A(9), fact
      do i=1, 9
        S(i) = S(i) + fact * A(i)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_s_p(S,A)
      implicit none
      integer i
      double precision S(3,*), A(*)
      S(1,1) = A(1)
      S(2,1) = A(4)
      S(3,1) = A(6)
      S(1,2) = A(4)
      S(2,2) = A(2)
      S(3,2) = A(5)
      S(1,3) = A(6)
      S(2,3) = A(5)
      S(3,3) = A(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_s_m(S,A)
      implicit none
      integer i
      double precision S(3,*), A(*)
      S(1,1) = - A(1)
      S(2,1) = - A(4)
      S(3,1) = - A(6)
      S(1,2) = - A(4)
      S(2,2) = - A(2)
      S(3,2) = - A(5)
      S(1,3) = - A(6)
      S(2,3) = - A(5)
      S(3,3) = - A(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_s_f(S,A,fact)
      implicit none
      integer i
      double precision S(3,*), A(*), fact
      S(1,1) = fact * A(1)
      S(2,1) = fact * A(4)
      S(3,1) = fact * A(6)
      S(1,2) = fact * A(4)
      S(2,2) = fact * A(2)
      S(3,2) = fact * A(5)
      S(1,3) = fact * A(6)
      S(2,3) = fact * A(5)
      S(3,3) = fact * A(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_s_ap(S,A)
      implicit none
      integer i
      double precision S(3,*), A(*)
      S(1,1) = S(1,1) + A(1)
      S(2,1) = S(2,1) + A(4)
      S(3,1) = S(3,1) + A(6)
      S(1,2) = S(1,2) + A(4)
      S(2,2) = S(2,2) + A(2)
      S(3,2) = S(3,2) + A(5)
      S(1,3) = S(1,3) + A(6)
      S(2,3) = S(2,3) + A(5)
      S(3,3) = S(3,3) + A(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_s_am(S,A)
      implicit none
      integer i
      double precision S(3,*), A(*)
      S(1,1) = S(1,1) - A(1)
      S(2,1) = S(2,1) - A(4)
      S(3,1) = S(3,1) - A(6)
      S(1,2) = S(1,2) - A(4)
      S(2,2) = S(2,2) - A(2)
      S(3,2) = S(3,2) - A(5)
      S(1,3) = S(1,3) - A(6)
      S(2,3) = S(2,3) - A(5)
      S(3,3) = S(3,3) - A(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_u_s_af(S,A,fact)
      implicit none
      integer i
      double precision S(3,*), A(*), fact
      S(1,1) = S(1,1) + fact * A(1)
      S(2,1) = S(2,1) + fact * A(4)
      S(3,1) = S(3,1) + fact * A(6)
      S(1,2) = S(1,2) + fact * A(4)
      S(2,2) = S(2,2) + fact * A(2)
      S(3,2) = S(3,2) + fact * A(5)
      S(1,3) = S(1,3) + fact * A(6)
      S(2,3) = S(2,3) + fact * A(5)
      S(3,3) = S(3,3) + fact * A(3)
      return
      end
c
c=============================================================================
c  Symmetrischer Teil eines Tensors 2.Stufe   S = (A + A^t) / 2
c-----------------------------------------------------------------------------
      subroutine sympart_p(S,A)
      implicit none
      double precision A(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = A(1,1)
      S(2) = A(2,2)
      S(3) = A(3,3)
      S(4) = r1d2 * (A(1,2) + A(2,1))
      S(5) = r1d2 * (A(2,3) + A(3,2))
      S(6) = r1d2 * (A(1,3) + A(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_m(S,A)
      implicit none
      double precision A(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = - A(1,1)
      S(2) = - A(2,2)
      S(3) = - A(3,3)
      S(4) = - r1d2 * (A(1,2) + A(2,1))
      S(5) = - r1d2 * (A(2,3) + A(3,2))
      S(6) = - r1d2 * (A(1,3) + A(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_f(S,A,fact)
      implicit none
      double precision A(3,3), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = fact * A(1,1)
      S(2) = fact * A(2,2)
      S(3) = fact * A(3,3)
      S(4) = factd2 * (A(1,2) + A(2,1))
      S(5) = factd2 * (A(2,3) + A(3,2))
      S(6) = factd2 * (A(1,3) + A(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_ap(S,A)
      implicit none
      double precision A(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) + A(1,1)
      S(2) = S(2) + A(2,2)
      S(3) = S(3) + A(3,3)
      S(4) = S(4) + r1d2 * (A(1,2) + A(2,1))
      S(5) = S(5) + r1d2 * (A(2,3) + A(3,2))
      S(6) = S(6) + r1d2 * (A(1,3) + A(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_am(S,A)
      implicit none
      double precision A(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) - A(1,1)
      S(2) = S(2) - A(2,2)
      S(3) = S(3) - A(3,3)
      S(4) = S(4) - r1d2 * (A(1,2) + A(2,1))
      S(5) = S(5) - r1d2 * (A(2,3) + A(3,2))
      S(6) = S(6) - r1d2 * (A(1,3) + A(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_af(S,A,fact)
      implicit none
      double precision A(3,3), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = S(1) + fact * A(1,1)
      S(2) = S(2) + fact * A(2,2)
      S(3) = S(3) + fact * A(3,3)
      S(4) = S(4) + factd2 * (A(1,2) + A(2,1))
      S(5) = S(5) + factd2 * (A(2,3) + A(3,2))
      S(6) = S(6) + factd2 * (A(1,3) + A(3,1))
      return
      end
c
c=============================================================================
c    S = sym(A B)   A,S sym., B unsym. 2.Stufe
c-----------------------------------------------------------------------------
      subroutine sympart_AU_p(S,A,B)
      implicit none
      double precision A(6), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
      S(2) = A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
      S(3) = A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      S(4) = r1d2 * (A(1)*B(1,2)+A(4)*(B(2,2)+B(1,1))+A(6)*B(3,2)
     *                          +A(2)*B(2,1)+A(5)*B(3,1))
      S(5) = r1d2 * (A(4)*B(1,3)+A(2)*B(2,3)+A(5)*(B(3,3)+B(2,2))
     *             + A(6)*B(1,2)            +A(3)*B(3,2))
      S(6) = r1d2 * (A(1)*B(1,3)+A(4)*B(2,3)+A(6)*(B(3,3)+B(1,1))
     *                          +A(5)*B(2,1)+A(3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_AU_m(S,A,B)
      implicit none
      double precision A(6), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = - A(1)*B(1,1) - A(4)*B(2,1) - A(6)*B(3,1)
      S(2) = - A(4)*B(1,2) - A(2)*B(2,2) - A(5)*B(3,2)
      S(3) = - A(6)*B(1,3) - A(5)*B(2,3) - A(3)*B(3,3)
      S(4) = - r1d2 * (A(1)*B(1,2)+A(4)*(B(2,2)+B(1,1))+A(6)*B(3,2)
     *                          +A(2)*B(2,1)+A(5)*B(3,1))
      S(5) = - r1d2 * (A(4)*B(1,3)+A(2)*B(2,3)+A(5)*(B(3,3)+B(2,2))
     *             + A(6)*B(1,2)            +A(3)*B(3,2))
      S(6) = - r1d2 * (A(1)*B(1,3)+A(4)*B(2,3)+A(6)*(B(3,3)+B(1,1))
     *                          +A(5)*B(2,1)+A(3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_AU_f(S,A,B,fact)
      implicit none
      double precision A(6), B(3,3), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = fact * (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(2) = fact * (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(3) = fact * (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      S(4) = factd2 * (A(1)*B(1,2)+A(4)*(B(2,2)+B(1,1))+A(6)*B(3,2)
     *                          +A(2)*B(2,1)+A(5)*B(3,1))
      S(5) = factd2 * (A(4)*B(1,3)+A(2)*B(2,3)+A(5)*(B(3,3)+B(2,2))
     *             + A(6)*B(1,2)            +A(3)*B(3,2))
      S(6) = factd2 * (A(1)*B(1,3)+A(4)*B(2,3)+A(6)*(B(3,3)+B(1,1))
     *                          +A(5)*B(2,1)+A(3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_AU_ap(S,A,B)
      implicit none
      double precision A(6), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) + A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
      S(2) = S(2) + A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
      S(3) = S(3) + A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      S(4) = S(4) + r1d2 * (A(1)*B(1,2)+A(4)*(B(2,2)+B(1,1))+A(6)*B(3,2)
     *                          +A(2)*B(2,1)+A(5)*B(3,1))
      S(5) = S(5) + r1d2 * (A(4)*B(1,3)+A(2)*B(2,3)+A(5)*(B(3,3)+B(2,2))
     *             + A(6)*B(1,2)            +A(3)*B(3,2))
      S(6) = S(6) + r1d2 * (A(1)*B(1,3)+A(4)*B(2,3)+A(6)*(B(3,3)+B(1,1))
     *                          +A(5)*B(2,1)+A(3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_AU_am(S,A,B)
      implicit none
      double precision A(6), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) - A(1)*B(1,1) - A(4)*B(2,1) - A(6)*B(3,1)
      S(2) = S(2) - A(4)*B(1,2) - A(2)*B(2,2) - A(5)*B(3,2)
      S(3) = S(3) - A(6)*B(1,3) - A(5)*B(2,3) - A(3)*B(3,3)
      S(4) = S(4) - r1d2 * (A(1)*B(1,2)+A(4)*(B(2,2)+B(1,1))+A(6)*B(3,2)
     *                          +A(2)*B(2,1)+A(5)*B(3,1))
      S(5) = S(5) - r1d2 * (A(4)*B(1,3)+A(2)*B(2,3)+A(5)*(B(3,3)+B(2,2))
     *             + A(6)*B(1,2)            +A(3)*B(3,2))
      S(6) = S(6) - r1d2 * (A(1)*B(1,3)+A(4)*B(2,3)+A(6)*(B(3,3)+B(1,1))
     *                          +A(5)*B(2,1)+A(3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_AU_af(S,A,B,fact)
      implicit none
      double precision A(6), B(3,3), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = S(1) + fact * (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(2) = S(2) + fact * (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(3) = S(3) + fact * (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      S(4) = S(4) + factd2 * (A(1)*B(1,2)+A(4)*(B(2,2)+B(1,1))
     *                         +A(6)*B(3,2)+A(2)*B(2,1)+A(5)*B(3,1))
      S(5) = S(5) + factd2 * (A(4)*B(1,3)+A(2)*B(2,3)
     *              +A(5)*(B(3,3)+B(2,2)) + A(6)*B(1,2)+A(3)*B(3,2))
      S(6) = S(6) + factd2 * (A(1)*B(1,3)+A(4)*B(2,3)
     *                +A(6)*(B(3,3)+B(1,1))+A(5)*B(2,1)+A(3)*B(3,1))
      return
      end
c
c=============================================================================
c  S = sym(A B)   B,S sym., A unsym. 2.Stufe
c-----------------------------------------------------------------------------
      subroutine sympart_UA_p(S,A,B)
      implicit none
      double precision B(6), A(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6)
      S(2) = A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5)
      S(3) = A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3)
      S(4) = r1d2 * ((A(1,1)+A(2,2))*B(4)+A(1,2)*B(2)+A(1,3)*B(5) 
     *              +A(2,1)*B(1)+A(2,3)*B(6))
      S(5) = r1d2 * (A(2,1)*B(6)+(A(2,2)+A(3,3))*B(5)+A(2,3)*B(3) 
     *              +A(3,1)*B(4)+A(3,2)*B(2))
      S(6) = r1d2 * ((A(1,1)+A(3,3))*B(6)+A(1,2)*B(5)+A(1,3)*B(3) 
     *              +A(3,1)*B(1)+A(3,2)*B(4))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_UA_m(S,A,B)
      implicit none
      double precision B(6), A(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = - A(1,1)*B(1) - A(1,2)*B(4) - A(1,3)*B(6)
      S(2) = - A(2,1)*B(4) - A(2,2)*B(2) - A(2,3)*B(5)
      S(3) = - A(3,1)*B(6) - A(3,2)*B(5) - A(3,3)*B(3)
      S(4) = - r1d2 * ((A(1,1)+A(2,2))*B(4)+A(1,2)*B(2)+A(1,3)*B(5)
     *              +A(2,1)*B(1)+A(2,3)*B(6))
      S(5) = - r1d2 * (A(2,1)*B(6)+(A(2,2)+A(3,3))*B(5)+A(2,3)*B(3)
     *              +A(3,1)*B(4)+A(3,2)*B(2))
      S(6) = - r1d2 * ((A(1,1)+A(3,3))*B(6)+A(1,2)*B(5)+A(1,3)*B(3)
     *              +A(3,1)*B(1)+A(3,2)*B(4))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_UA_f(S,A,B,fact)
      implicit none
      double precision B(6), A(3,3), S(6), factd2, fact
      factd2 = fact * 0.5d0
      S(1) = fact * (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(2) = fact * (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(3) = fact * (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      S(4) = factd2 * ((A(1,1)+A(2,2))*B(4)+A(1,2)*B(2)+A(1,3)*B(5) 
     *              +A(2,1)*B(1)+A(2,3)*B(6))
      S(5) = factd2 * (A(2,1)*B(6)+(A(2,2)+A(3,3))*B(5)+A(2,3)*B(3) 
     *              +A(3,1)*B(4)+A(3,2)*B(2))
      S(6) = factd2 * ((A(1,1)+A(3,3))*B(6)+A(1,2)*B(5)+A(1,3)*B(3) 
     *              +A(3,1)*B(1)+A(3,2)*B(4))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_UA_ap(S,A,B)
      implicit none
      double precision B(6), A(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) + A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6)
      S(2) = S(2) + A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5)
      S(3) = S(3) + A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3)
      S(4) = S(4) + r1d2 * ((A(1,1)+A(2,2))*B(4)+A(1,2)*B(2)+A(1,3)*B(5) 
     *              +A(2,1)*B(1)+A(2,3)*B(6))
      S(5) = S(5) + r1d2 * (A(2,1)*B(6)+(A(2,2)+A(3,3))*B(5)+A(2,3)*B(3) 
     *              +A(3,1)*B(4)+A(3,2)*B(2))
      S(6) = S(6) + r1d2 * ((A(1,1)+A(3,3))*B(6)+A(1,2)*B(5)+A(1,3)*B(3) 
     *              +A(3,1)*B(1)+A(3,2)*B(4))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_UA_am(S,A,B)
      implicit none
      double precision B(6), A(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) - A(1,1)*B(1) - A(1,2)*B(4) - A(1,3)*B(6)
      S(2) = S(2) - A(2,1)*B(4) - A(2,2)*B(2) - A(2,3)*B(5)
      S(3) = S(3) - A(3,1)*B(6) - A(3,2)*B(5) - A(3,3)*B(3)
      S(4) = S(4) - r1d2 * ((A(1,1)+A(2,2))*B(4)+A(1,2)*B(2)+A(1,3)*B(5)
     *              +A(2,1)*B(1)+A(2,3)*B(6))
      S(5) = S(5) - r1d2 * (A(2,1)*B(6)+(A(2,2)+A(3,3))*B(5)+A(2,3)*B(3)
     *              +A(3,1)*B(4)+A(3,2)*B(2))
      S(6) = S(6) - r1d2 * ((A(1,1)+A(3,3))*B(6)+A(1,2)*B(5)+A(1,3)*B(3)
     *              +A(3,1)*B(1)+A(3,2)*B(4))
      return
      end
c-----------------------------------------------------------------------------
      subroutine sympart_UA_af(S,A,B,fact)
      implicit none
      double precision B(6), A(3,3), S(6), factd2, fact
      factd2 = fact * 0.5d0
      S(1) = S(1) + fact * (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(2) = S(2) + fact * (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(3) = S(3) + fact * (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      S(4) = S(4) + factd2 * ((A(1,1)+A(2,2))*B(4)+A(1,2)*B(2) 
     *              +A(1,3)*B(5)+A(2,1)*B(1)+A(2,3)*B(6))
      S(5) = S(5) + factd2 * (A(2,1)*B(6)+(A(2,2)+A(3,3))*B(5) 
     *              +A(2,3)*B(3)+A(3,1)*B(4)+A(3,2)*B(2))
      S(6) = S(6) + factd2 * ((A(1,1)+A(3,3))*B(6)+A(1,2)*B(5) 
     *              +A(1,3)*B(3)+A(3,1)*B(1)+A(3,2)*B(4))
      return
      end
c
c=============================================================================
c  Multiplikation von Tensoren       sym = sym x sym
c-----------------------------------------------------------------------------
      subroutine mult_s_ss_p(S,A,B)
      implicit none
      double precision A(6), B(6), S(6)
      S(1) = A(1)*B(1) + A(4)*B(4) + A(6)*B(6)
      S(2) = A(4)*B(4) + A(2)*B(2) + A(5)*B(5)
      S(3) = A(6)*B(6) + A(5)*B(5) + A(3)*B(3)
      S(4) = A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
      S(5) = A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
      S(6) = A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_ss_m(S,A,B)
      implicit none
      double precision A(6), B(6), S(6)
      S(1) = - (A(1)*B(1) + A(4)*B(4) + A(6)*B(6))
      S(2) = - (A(4)*B(4) + A(2)*B(2) + A(5)*B(5))
      S(3) = - (A(6)*B(6) + A(5)*B(5) + A(3)*B(3))
      S(4) = - (A(4)*B(1) + A(2)*B(4) + A(5)*B(6))
      S(5) = - (A(6)*B(4) + A(5)*B(2) + A(3)*B(5))
      S(6) = - (A(6)*B(1) + A(5)*B(4) + A(3)*B(6))
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_ss_f(S,A,B,fact)
      implicit none
      double precision A(6), B(6), S(6), fact
      S(1) = fact * (A(1)*B(1) + A(4)*B(4) + A(6)*B(6))
      S(2) = fact * (A(4)*B(4) + A(2)*B(2) + A(5)*B(5))
      S(3) = fact * (A(6)*B(6) + A(5)*B(5) + A(3)*B(3))
      S(4) = fact * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6))
      S(5) = fact * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5))
      S(6) = fact * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6))
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_ss_ap(S,A,B)
      implicit none
      double precision A(6), B(6), S(6)
      S(1) = S(1) + A(1)*B(1) + A(4)*B(4) + A(6)*B(6)
      S(2) = S(2) + A(4)*B(4) + A(2)*B(2) + A(5)*B(5)
      S(3) = S(3) + A(6)*B(6) + A(5)*B(5) + A(3)*B(3)
      S(4) = S(4) + A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
      S(5) = S(5) + A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
      S(6) = S(6) + A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_ss_am(S,A,B)
      implicit none
      double precision A(6), B(6), S(6)
      S(1) = S(1) - A(1)*B(1) - A(4)*B(4) - A(6)*B(6)
      S(2) = S(2) - A(4)*B(4) - A(2)*B(2) - A(5)*B(5)
      S(3) = S(3) - A(6)*B(6) - A(5)*B(5) - A(3)*B(3)
      S(4) = S(4) - A(4)*B(1) - A(2)*B(4) - A(5)*B(6)
      S(5) = S(5) - A(6)*B(4) - A(5)*B(2) - A(3)*B(5)
      S(6) = S(6) - A(6)*B(1) - A(5)*B(4) - A(3)*B(6)
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_ss_af(S,A,B,fact)
      implicit none
      double precision A(6), B(6), S(6), fact
      S(1) = S(1) + fact * (A(1)*B(1) + A(4)*B(4) + A(6)*B(6))
      S(2) = S(2) + fact * (A(4)*B(4) + A(2)*B(2) + A(5)*B(5))
      S(3) = S(3) + fact * (A(6)*B(6) + A(5)*B(5) + A(3)*B(3))
      S(4) = S(4) + fact * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6))
      S(5) = S(5) + fact * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5))
      S(6) = S(6) + fact * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6))
      return 
      end
c
c============================================================================
c  Multiplikation von Tensoren  sym = sym x unsym
c----------------------------------------------------------------------------
      subroutine mult_s_su_p(S,A,B)
      implicit none 
      double precision A(6), B(3,3), S(6)
      S(1) = A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
      S(2) = A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
      S(3) = A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      S(4) = A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
      S(5) = A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
      S(6) = A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)     
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_s_su_m(S,A,B)
      implicit none 
      double precision A(6), B(3,3), S(6)
      S(1) = - (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(2) = - (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(3) = - (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      S(4) = - (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1))
      S(5) = - (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2))
      S(6) = - (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1))     
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_s_su_f(S,A,B,fact)
      implicit none 
      double precision A(6), B(3,3), S(6), fact
      S(1) = fact * (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(2) = fact * (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(3) = fact * (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      S(4) = fact * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1))
      S(5) = fact * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2))
      S(6) = fact * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1))     
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_s_su_ap(S,A,B)
      implicit none 
      double precision A(6), B(3,3), S(6)
      S(1) = S(1) + A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
      S(2) = S(2) + A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
      S(3) = S(3) + A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      S(4) = S(4) + A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
      S(5) = S(5) + A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
      S(6) = S(6) + A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)     
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_s_su_am(S,A,B)
      implicit none 
      double precision A(6), B(3,3), S(6)
      S(1) = S(1) - A(1)*B(1,1) - A(4)*B(2,1) - A(6)*B(3,1)
      S(2) = S(2) - A(4)*B(1,2) - A(2)*B(2,2) - A(5)*B(3,2)
      S(3) = S(3) - A(6)*B(1,3) - A(5)*B(2,3) - A(3)*B(3,3)
      S(4) = S(4) - A(4)*B(1,1) - A(2)*B(2,1) - A(5)*B(3,1)
      S(5) = S(5) - A(6)*B(1,2) - A(5)*B(2,2) - A(3)*B(3,2)
      S(6) = S(6) - A(6)*B(1,1) - A(5)*B(2,1) - A(3)*B(3,1)     
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_s_su_af(S,A,B,fact)
      implicit none 
      double precision A(6), B(3,3), S(6), fact
      S(1) = S(1) + fact * (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(2) = S(2) + fact * (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(3) = S(3) + fact * (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      S(4) = S(4) + fact * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1))
      S(5) = S(5) + fact * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2))
      S(6) = S(6) + fact * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1))     
      return
      end
c
c=============================================================================
c  Multiplikation von Tensoren     sym = unsym x sym
c-----------------------------------------------------------------------------
      subroutine mult_s_us_p(S,A,B)
      implicit none
      double precision A(3,3), B(6), S(6)
      S(1) = A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6)
      S(2) = A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5)
      S(3) = A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3)
      S(4) = A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
      S(5) = A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
      S(6) = A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)     
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_us_m(S,A,B)
      implicit none
      double precision A(3,3), B(6), S(6)
      S(1) = - (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(2) = - (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(3) = - (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      S(4) = - (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5))
      S(5) = - (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3))
      S(6) = - (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3))     
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_us_f(S,A,B,fact)
      implicit none
      double precision A(3,3), B(6), S(6), fact
      S(1) = fact * (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(2) = fact * (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(3) = fact * (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      S(4) = fact * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5))
      S(5) = fact * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3))
      S(6) = fact * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3))     
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_us_ap(S,A,B)
      implicit none
      double precision A(3,3), B(6), S(6)
      S(1) = S(1) + A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6)
      S(2) = S(2) + A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5)
      S(3) = S(3) + A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3)
      S(4) = S(4) + A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
      S(5) = S(5) + A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
      S(6) = S(6) + A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)     
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_us_am(S,A,B)
      implicit none
      double precision A(3,3), B(6), S(6), fact
      S(1) = S(1) - A(1,1)*B(1) - A(1,2)*B(4) - A(1,3)*B(6)
      S(2) = S(2) - A(2,1)*B(4) - A(2,2)*B(2) - A(2,3)*B(5)
      S(3) = S(3) - A(3,1)*B(6) - A(3,2)*B(5) - A(3,3)*B(3)
      S(4) = S(4) - A(1,1)*B(4) - A(1,2)*B(2) - A(1,3)*B(5)
      S(5) = S(5) - A(2,1)*B(6) - A(2,2)*B(5) - A(2,3)*B(3)
      S(6) = S(6) - A(1,1)*B(6) - A(1,2)*B(5) - A(1,3)*B(3)     
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_us_af(S,A,B,fact)
      implicit none
      double precision A(3,3), B(6), S(6), fact
      S(1) = S(1) + fact * (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(2) = S(2) + fact * (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(3) = S(3) + fact * (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      S(4) = S(4) + fact * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5))
      S(5) = S(5) + fact * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3))
      S(6) = S(6) + fact * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3))     
      return
      end
c
c=============================================================================
c   Multiplikation von Tensoren    sym = unsym x unsym
c-----------------------------------------------------------------------------
      subroutine mult_s_uu_p(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
      S(2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
      S(3) = A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)
      S(4) = A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2)
      S(5) = A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)
      S(6) = A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uu_m(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = - (A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1))
      S(2) = - (A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2))
      S(3) = - (A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3))
      S(4) = - (A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2))
      S(5) = - (A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3))
      S(6) = - (A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uu_f(S,A,B,fact)
      implicit none
      double precision A(3,3), B(3,3), S(6), fact
      S(1) = fact * (A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1))
      S(2) = fact * (A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2))
      S(3) = fact * (A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3))
      S(4) = fact * (A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2))
      S(5) = fact * (A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3))
      S(6) = fact * (A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uu_ap(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = S(1) + A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
      S(2) = S(2) + A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
      S(3) = S(3) + A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)
      S(4) = S(4) + A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2)
      S(5) = S(5) + A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)
      S(6) = S(6) + A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uu_am(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = S(1) - A(1,1)*B(1,1) - A(1,2)*B(2,1) - A(1,3)*B(3,1)
      S(2) = S(2) - A(2,1)*B(1,2) - A(2,2)*B(2,2) - A(2,3)*B(3,2)
      S(3) = S(3) - A(3,1)*B(1,3) - A(3,2)*B(2,3) - A(3,3)*B(3,3)
      S(4) = S(4) - A(1,1)*B(1,2) - A(1,2)*B(2,2) - A(1,3)*B(3,2)
      S(5) = S(5) - A(2,1)*B(1,3) - A(2,2)*B(2,3) - A(2,3)*B(3,3)
      S(6) = S(6) - A(1,1)*B(1,3) - A(1,2)*B(2,3) - A(1,3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uu_af(S,A,B,fact)
      implicit none
      double precision A(3,3), B(3,3), S(6), fact
      S(1) = S(1) + fact * (A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1))
      S(2) = S(2) + fact * (A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2))
      S(3) = S(3) + fact * (A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3))
      S(4) = S(4) + fact * (A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2))
      S(5) = S(5) + fact * (A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3))
      S(6) = S(6) + fact * (A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3))
      return
      end
c
c=============================================================================
c   Multiplikation von Tensoren    sym = unsym x unsym^T
c-----------------------------------------------------------------------------
      subroutine mult_s_uuT_p(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = A(1,1)*B(1,1) + A(1,2)*B(1,2) + A(1,3)*B(1,3)
      S(2) = A(2,1)*B(2,1) + A(2,2)*B(2,2) + A(2,3)*B(2,3)
      S(3) = A(3,1)*B(3,1) + A(3,2)*B(3,2) + A(3,3)*B(3,3)
      S(4) = A(1,1)*B(2,1) + A(1,2)*B(2,2) + A(1,3)*B(2,3)
      S(5) = A(2,1)*B(3,1) + A(2,2)*B(3,2) + A(2,3)*B(3,3)
      S(6) = A(1,1)*B(3,1) + A(1,2)*B(3,2) + A(1,3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uuT_m(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = - A(1,1)*B(1,1) - A(1,2)*B(1,2) - A(1,3)*B(1,3)
      S(2) = - A(2,1)*B(2,1) - A(2,2)*B(2,2) - A(2,3)*B(2,3)
      S(3) = - A(3,1)*B(3,1) - A(3,2)*B(3,2) - A(3,3)*B(3,3)
      S(4) = - A(1,1)*B(2,1) - A(1,2)*B(2,2) - A(1,3)*B(2,3)
      S(5) = - A(2,1)*B(3,1) - A(2,2)*B(3,2) - A(2,3)*B(3,3)
      S(6) = - A(1,1)*B(3,1) - A(1,2)*B(3,2) - A(1,3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uuT_f(S,A,B,fact)
      implicit none
      double precision A(3,3), B(3,3), S(6), fact
      S(1) = fact * (A(1,1)*B(1,1) + A(1,2)*B(1,2) + A(1,3)*B(1,3))
      S(2) = fact * (A(2,1)*B(2,1) + A(2,2)*B(2,2) + A(2,3)*B(2,3))
      S(3) = fact * (A(3,1)*B(3,1) + A(3,2)*B(3,2) + A(3,3)*B(3,3))
      S(4) = fact * (A(1,1)*B(2,1) + A(1,2)*B(2,2) + A(1,3)*B(2,3))
      S(5) = fact * (A(2,1)*B(3,1) + A(2,2)*B(3,2) + A(2,3)*B(3,3))
      S(6) = fact * (A(1,1)*B(3,1) + A(1,2)*B(3,2) + A(1,3)*B(3,3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uuT_ap(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = S(1) + A(1,1)*B(1,1) + A(1,2)*B(1,2) + A(1,3)*B(1,3)
      S(2) = S(2) + A(2,1)*B(2,1) + A(2,2)*B(2,2) + A(2,3)*B(2,3)
      S(3) = S(3) + A(3,1)*B(3,1) + A(3,2)*B(3,2) + A(3,3)*B(3,3)
      S(4) = S(4) + A(1,1)*B(2,1) + A(1,2)*B(2,2) + A(1,3)*B(2,3)
      S(5) = S(5) + A(2,1)*B(3,1) + A(2,2)*B(3,2) + A(2,3)*B(3,3)
      S(6) = S(6) + A(1,1)*B(3,1) + A(1,2)*B(3,2) + A(1,3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uuT_am(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = S(1) - A(1,1)*B(1,1) - A(1,2)*B(1,2) - A(1,3)*B(1,3)
      S(2) = S(2) - A(2,1)*B(2,1) - A(2,2)*B(2,2) - A(2,3)*B(2,3)
      S(3) = S(3) - A(3,1)*B(3,1) - A(3,2)*B(3,2) - A(3,3)*B(3,3)
      S(4) = S(4) - A(1,1)*B(2,1) - A(1,2)*B(2,2) - A(1,3)*B(2,3)
      S(5) = S(5) - A(2,1)*B(3,1) - A(2,2)*B(3,2) - A(2,3)*B(3,3)
      S(6) = S(6) - A(1,1)*B(3,1) - A(1,2)*B(3,2) - A(1,3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uuT_af(S,A,B,fact)
      implicit none
      double precision A(3,3), B(3,3), S(6), fact
      S(1) = S(1) + fact*(A(1,1)*B(1,1) + A(1,2)*B(1,2) + A(1,3)*B(1,3))
      S(2) = S(2) + fact*(A(2,1)*B(2,1) + A(2,2)*B(2,2) + A(2,3)*B(2,3))
      S(3) = S(3) + fact*(A(3,1)*B(3,1) + A(3,2)*B(3,2) + A(3,3)*B(3,3))
      S(4) = S(4) + fact*(A(1,1)*B(2,1) + A(1,2)*B(2,2) + A(1,3)*B(2,3))
      S(5) = S(5) + fact*(A(2,1)*B(3,1) + A(2,2)*B(3,2) + A(2,3)*B(3,3))
      S(6) = S(6) + fact*(A(1,1)*B(3,1) + A(1,2)*B(3,2) + A(1,3)*B(3,3))
      return
      end
c
c=============================================================================
c   Multiplikation von Tensoren    sym = unsym^T x unsym
c-----------------------------------------------------------------------------
      subroutine mult_s_uTu_p(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = A(1,1)*B(1,1) + A(2,1)*B(2,1) + A(3,1)*B(3,1)
      S(2) = A(1,2)*B(1,2) + A(2,2)*B(2,2) + A(3,2)*B(3,2)
      S(3) = A(1,3)*B(1,3) + A(2,3)*B(2,3) + A(3,3)*B(3,3)
      S(4) = A(1,1)*B(1,2) + A(2,1)*B(2,2) + A(3,1)*B(3,2)
      S(5) = A(1,2)*B(1,3) + A(2,2)*B(2,3) + A(3,2)*B(3,3)
      S(6) = A(1,1)*B(1,3) + A(2,1)*B(2,3) + A(3,1)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uTu_m(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = - A(1,1)*B(1,1) - A(2,1)*B(2,1) - A(3,1)*B(3,1)
      S(2) = - A(1,2)*B(1,2) - A(2,2)*B(2,2) - A(3,2)*B(3,2)
      S(3) = - A(1,3)*B(1,3) - A(2,3)*B(2,3) - A(3,3)*B(3,3)
      S(4) = - A(1,1)*B(1,2) - A(2,1)*B(2,2) - A(3,1)*B(3,2)
      S(5) = - A(1,2)*B(1,3) - A(2,2)*B(2,3) - A(3,2)*B(3,3)
      S(6) = - A(1,1)*B(1,3) - A(2,1)*B(2,3) - A(3,1)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uTu_f(S,A,B,fact)
      implicit none
      double precision A(3,3), B(3,3), S(6), fact
      S(1) = fact * (A(1,1)*B(1,1) + A(2,1)*B(2,1) + A(3,1)*B(3,1))
      S(2) = fact * (A(1,2)*B(1,2) + A(2,2)*B(2,2) + A(3,2)*B(3,2))
      S(3) = fact * (A(1,3)*B(1,3) + A(2,3)*B(2,3) + A(3,3)*B(3,3))
      S(4) = fact * (A(1,1)*B(1,2) + A(2,1)*B(2,2) + A(3,1)*B(3,2))
      S(5) = fact * (A(1,2)*B(1,3) + A(2,2)*B(2,3) + A(3,2)*B(3,3))
      S(6) = fact * (A(1,1)*B(1,3) + A(2,1)*B(2,3) + A(3,1)*B(3,3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uTu_ap(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = S(1) + A(1,1)*B(1,1) + A(2,1)*B(2,1) + A(3,1)*B(3,1)
      S(2) = S(2) + A(1,2)*B(1,2) + A(2,2)*B(2,2) + A(3,2)*B(3,2)
      S(3) = S(3) + A(1,3)*B(1,3) + A(2,3)*B(2,3) + A(3,3)*B(3,3)
      S(4) = S(4) + A(1,1)*B(1,2) + A(2,1)*B(2,2) + A(3,1)*B(3,2)
      S(5) = S(5) + A(1,2)*B(1,3) + A(2,2)*B(2,3) + A(3,2)*B(3,3)
      S(6) = S(6) + A(1,1)*B(1,3) + A(2,1)*B(2,3) + A(3,1)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uTu_am(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6)
      S(1) = S(1) - A(1,1)*B(1,1) - A(2,1)*B(2,1) - A(3,1)*B(3,1)
      S(2) = S(2) - A(1,2)*B(1,2) - A(2,2)*B(2,2) - A(3,2)*B(3,2)
      S(3) = S(3) - A(1,3)*B(1,3) - A(2,3)*B(2,3) - A(3,3)*B(3,3)
      S(4) = S(4) - A(1,1)*B(1,2) - A(2,1)*B(2,2) - A(3,1)*B(3,2)
      S(5) = S(5) - A(1,2)*B(1,3) - A(2,2)*B(2,3) - A(3,2)*B(3,3)
      S(6) = S(6) - A(1,1)*B(1,3) - A(2,1)*B(2,3) - A(3,1)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_uTu_af(S,A,B,fact)
      implicit none
      double precision A(3,3), B(3,3), S(6), fact
      S(1) = S(1) + fact*(A(1,1)*B(1,1) + A(2,1)*B(2,1) + A(3,1)*B(3,1))
      S(2) = S(2) + fact*(A(1,2)*B(1,2) + A(2,2)*B(2,2) + A(3,2)*B(3,2))
      S(3) = S(3) + fact*(A(1,3)*B(1,3) + A(2,3)*B(2,3) + A(3,3)*B(3,3))
      S(4) = S(4) + fact*(A(1,1)*B(1,2) + A(2,1)*B(2,2) + A(3,1)*B(3,2))
      S(5) = S(5) + fact*(A(1,2)*B(1,3) + A(2,2)*B(2,3) + A(3,2)*B(3,3))
      S(6) = S(6) + fact*(A(1,1)*B(1,3) + A(2,1)*B(2,3) + A(3,1)*B(3,3))
      return
      end
c
c============================================================================
c  Multiplikation von Tensoren  sym = sym x sym   mit Mittelung
c   (Mittelung kann notwendig sein in Verbindung mit diff-Routinen)
c----------------------------------------------------------------------------
      subroutine mult_sd_ss_p(S,A,B)
      implicit none 
      double precision A(6), B(6), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = A(1)*B(1) + A(4)*B(4) + A(6)*B(6)
      S(2) = A(4)*B(4) + A(2)*B(2) + A(5)*B(5)
      S(3) = A(6)*B(6) + A(5)*B(5) + A(3)*B(3)
      S(4) = r1d2 * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
     *             + A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(5) = r1d2 * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
     *             + A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(6) = r1d2 * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
     *             + A(1)*B(6) + A(4)*B(5) + A(6)*B(3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_ss_m(S,A,B)
      implicit none 
      double precision A(6), B(6), S(6), r1d2
      r1d2 = - 0.5d0
      S(1) = - A(1)*B(1) - A(4)*B(4) - A(6)*B(6)
      S(2) = - A(4)*B(4) - A(2)*B(2) - A(5)*B(5)
      S(3) = - A(6)*B(6) - A(5)*B(5) - A(3)*B(3)
      S(4) = r1d2 * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
     *             + A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(5) = r1d2 * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
     *             + A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(6) = r1d2 * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
     *             + A(1)*B(6) + A(4)*B(5) + A(6)*B(3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_ss_f(S,A,B,fact)
      implicit none 
      double precision A(6), B(6), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = fact   * (A(1)*B(1) + A(4)*B(4) + A(6)*B(6))
      S(2) = fact   * (A(4)*B(4) + A(2)*B(2) + A(5)*B(5))
      S(3) = fact   * (A(6)*B(6) + A(5)*B(5) + A(3)*B(3))
      S(4) = factd2 * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
     *               + A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(5) = factd2 * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
     *               + A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(6) = factd2 * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
     *               + A(1)*B(6) + A(4)*B(5) + A(6)*B(3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_ss_ap(S,A,B)
      implicit none 
      double precision A(6), B(6), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) + A(1)*B(1) + A(4)*B(4) + A(6)*B(6)
      S(2) = S(2) + A(4)*B(4) + A(2)*B(2) + A(5)*B(5)
      S(3) = S(3) + A(6)*B(6) + A(5)*B(5) + A(3)*B(3)
      S(4) = S(4) + r1d2 * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
     *                    + A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(5) = S(5) + r1d2 * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
     *                    + A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(6) = S(6) + r1d2 * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
     *                    + A(1)*B(6) + A(4)*B(5) + A(6)*B(3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_ss_am(S,A,B)
      implicit none 
      double precision A(6), B(6), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) - A(1)*B(1) - A(4)*B(4) - A(6)*B(6)
      S(2) = S(2) - A(4)*B(4) - A(2)*B(2) - A(5)*B(5)
      S(3) = S(3) - A(6)*B(6) - A(5)*B(5) - A(3)*B(3)
      S(4) = S(4) - r1d2 * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
     *                    + A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(5) = S(5) - r1d2 * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
     *                    + A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(6) = S(6) - r1d2 * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
     *                    + A(1)*B(6) + A(4)*B(5) + A(6)*B(3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_ss_af(S,A,B,fact)
      implicit none 
      double precision A(6), B(6), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = S(1) + fact   * (A(1)*B(1) + A(4)*B(4) + A(6)*B(6))
      S(2) = S(2) + fact   * (A(4)*B(4) + A(2)*B(2) + A(5)*B(5))
      S(3) = S(3) + fact   * (A(6)*B(6) + A(5)*B(5) + A(3)*B(3))
      S(4) = S(4) + factd2 * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
     *                      + A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(5) = S(5) + factd2 * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
     *                      + A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(6) = S(6) + factd2 * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
     *                      + A(1)*B(6) + A(4)*B(5) + A(6)*B(3))
      return
      end
c
c============================================================================
c  Multiplikation von Tensoren  sym = sym x unsym   mit Mittelung
c   (Mittelung kann notwendig sein in Verbindung mit diff-Routinen)
c----------------------------------------------------------------------------
      subroutine mult_sd_su_p(S,A,B)
      implicit none 
      double precision A(6), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
      S(2) = A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
      S(3) = A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      S(4) = r1d2 * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
     *             + A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(5) = r1d2 * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
     *             + A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(6) = r1d2 * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)     
     *             + A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_su_m(S,A,B)
      implicit none 
      double precision A(6), B(3,3), S(6), r1d2
      r1d2 = - 0.5d0
      S(1) = - A(1)*B(1,1) - A(4)*B(2,1) - A(6)*B(3,1)
      S(2) = - A(4)*B(1,2) - A(2)*B(2,2) - A(5)*B(3,2)
      S(3) = - A(6)*B(1,3) - A(5)*B(2,3) - A(3)*B(3,3)
      S(4) = r1d2 * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
     *             + A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(5) = r1d2 * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
     *             + A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(6) = r1d2 * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)     
     *             + A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_su_f(S,A,B,fact)
      implicit none 
      double precision A(6), B(3,3), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = fact   * (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(2) = fact   * (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(3) = fact   * (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      S(4) = factd2 * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
     *               + A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(5) = factd2 * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
     *               + A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(6) = factd2 * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)     
     *               + A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_su_ap(S,A,B)
      implicit none 
      double precision A(6), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) + A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
      S(2) = S(2) + A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
      S(3) = S(3) + A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      S(4) = S(4) + r1d2 * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
     *                    + A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(5) = S(5) + r1d2 * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
     *                    + A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(6) = S(6) + r1d2 * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)     
     *                    + A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_su_am(S,A,B)
      implicit none 
      double precision A(6), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) - A(1)*B(1,1) - A(4)*B(2,1) - A(6)*B(3,1)
      S(2) = S(2) - A(4)*B(1,2) - A(2)*B(2,2) - A(5)*B(3,2)
      S(3) = S(3) - A(6)*B(1,3) - A(5)*B(2,3) - A(3)*B(3,3)
      S(4) = S(4) - r1d2 * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
     *                    + A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(5) = S(5) - r1d2 * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
     *                    + A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(6) = S(6) - r1d2 * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)     
     *                    + A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      return
      end
c----------------------------------------------------------------------------
      subroutine mult_sd_su_af(S,A,B,fact)
      implicit none 
      double precision A(6), B(3,3), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = S(1) + fact   * (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(2) = S(2) + fact   * (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(3) = S(3) + fact   * (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      S(4) = S(4) + factd2 * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
     *                      + A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(5) = S(5) + factd2 * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
     *                      + A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(6) = S(6) + factd2 * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)     
     *                      + A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      return
      end
c
c=============================================================================
c  Multiplikation von Tensoren     sym = unsym x sym    mit Mittelung
c   (Mittelung kann notwendig sein in Verbindung mit diff-Routinen)
c-----------------------------------------------------------------------------
      subroutine mult_sd_us_p(S,A,B)
      implicit none
      double precision S(6), A(3,3), B(6), r1d2
      r1d2 = 0.5d0
      S(1) = A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6)
      S(2) = A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5)
      S(3) = A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3)
      S(4) = r1d2 * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
     *             + A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(5) = r1d2 * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
     *             + A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(6) = r1d2 * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)
     *             + A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_us_m(S,A,B)
      implicit none
      double precision S(6), A(3,3), B(6), r1d2
      r1d2 = - 0.5d0
      S(1) = - A(1,1)*B(1) - A(1,2)*B(4) - A(1,3)*B(6)
      S(2) = - A(2,1)*B(4) - A(2,2)*B(2) - A(2,3)*B(5)
      S(3) = - A(3,1)*B(6) - A(3,2)*B(5) - A(3,3)*B(3)
      S(4) = r1d2 * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
     *             + A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(5) = r1d2 * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
     *             + A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(6) = r1d2 * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)
     *             + A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_us_f(S,A,B,fact)
      implicit none
      double precision S(6), A(3,3), B(6), fact, factd2
      factd2 = 0.5d0
      S(1) = fact   * (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(2) = fact   * (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(3) = fact   * (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      S(4) = factd2 * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
     *               + A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(5) = factd2 * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
     *               + A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(6) = factd2 * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)
     *               + A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_us_ap(S,A,B)
      implicit none
      double precision S(6), A(3,3), B(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) + A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6)
      S(2) = S(2) + A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5)
      S(3) = S(3) + A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3)
      S(4) = S(4) + r1d2 * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
     *                    + A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(5) = S(5) + r1d2 * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
     *                    + A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(6) = S(6) + r1d2 * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)
     *                    + A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_us_am(S,A,B)
      implicit none
      double precision S(6), A(3,3), B(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) - A(1,1)*B(1) - A(1,2)*B(4) - A(1,3)*B(6)
      S(2) = S(2) - A(2,1)*B(4) - A(2,2)*B(2) - A(2,3)*B(5)
      S(3) = S(3) - A(3,1)*B(6) - A(3,2)*B(5) - A(3,3)*B(3)
      S(4) = S(4) - r1d2 * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
     *                    + A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(5) = S(5) - r1d2 * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
     *                    + A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(6) = S(6) - r1d2 * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)
     *                    + A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_us_af(S,A,B,fact)
      implicit none
      double precision S(6), A(3,3), B(6), fact, factd2
      factd2 = 0.5d0
      S(1) = S(1) + fact   * (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(2) = S(2) + fact   * (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(3) = S(3) + fact   * (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      S(4) = S(4) + factd2 * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
     *                      + A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(5) = S(5) + factd2 * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
     *                      + A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(6) = S(6) + factd2 * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)
     *                      + A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      return
      end
c
c=============================================================================
c  Multiplikation von Tensoren     sym = unsym x unsym    mit Mittelung
c   (Mittelung kann notwendig sein in Verbindung mit diff-Routinen)
c-----------------------------------------------------------------------------
      subroutine mult_sd_uu_p(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
      S(2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
      S(3) = A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)
      S(4) = r1d2 * (A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2)
     *             + A(2,1)*B(1,1) + A(2,2)*B(2,1) + A(2,3)*B(3,1))
      S(5) = r1d2 * (A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)
     *             + A(3,1)*B(1,2) + A(3,2)*B(2,2) + A(3,3)*B(3,2))
      S(6) = r1d2 * (A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)
     *             + A(3,1)*B(1,1) + A(3,2)*B(2,1) + A(3,3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_uu_m(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6), r1d2
      r1d2 = - 0.5d0
      S(1) = - A(1,1)*B(1,1) - A(1,2)*B(2,1) - A(1,3)*B(3,1)
      S(2) = - A(2,1)*B(1,2) - A(2,2)*B(2,2) - A(2,3)*B(3,2)
      S(3) = - A(3,1)*B(1,3) - A(3,2)*B(2,3) - A(3,3)*B(3,3)
      S(4) = r1d2 * (A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2)
     *             + A(2,1)*B(1,1) + A(2,2)*B(2,1) + A(2,3)*B(3,1))
      S(5) = r1d2 * (A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)
     *             + A(3,1)*B(1,2) + A(3,2)*B(2,2) + A(3,3)*B(3,2))
      S(6) = r1d2 * (A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)
     *             + A(3,1)*B(1,1) + A(3,2)*B(2,1) + A(3,3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_uu_f(S,A,B,fact)
      implicit none
      double precision A(3,3), B(3,3), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = fact   * (A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1))
      S(2) = fact   * (A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2))
      S(3) = fact   * (A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3))
      S(4) = factd2 * (A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2)
     *               + A(2,1)*B(1,1) + A(2,2)*B(2,1) + A(2,3)*B(3,1))
      S(5) = factd2 * (A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)
     *               + A(3,1)*B(1,2) + A(3,2)*B(2,2) + A(3,3)*B(3,2))
      S(6) = factd2 * (A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)
     *               + A(3,1)*B(1,1) + A(3,2)*B(2,1) + A(3,3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_uu_ap(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) + A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
      S(2) = S(2) + A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
      S(3) = S(3) + A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)
      S(4) = S(4) + r1d2 * (A(1,1)*B(1,2) + A(1,2)*B(2,2)+A(1,3)*B(3,2)
     *                    + A(2,1)*B(1,1) + A(2,2)*B(2,1)+A(2,3)*B(3,1))
      S(5) = S(5) + r1d2 * (A(2,1)*B(1,3) + A(2,2)*B(2,3)+A(2,3)*B(3,3)
     *                    + A(3,1)*B(1,2) + A(3,2)*B(2,2)+A(3,3)*B(3,2))
      S(6) = S(6) + r1d2 * (A(1,1)*B(1,3) + A(1,2)*B(2,3)+A(1,3)*B(3,3)
     *                    + A(3,1)*B(1,1) + A(3,2)*B(2,1)+A(3,3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_uu_am(S,A,B)
      implicit none
      double precision A(3,3), B(3,3), S(6), r1d2
      r1d2 = 0.5d0
      S(1) = S(1) - A(1,1)*B(1,1) - A(1,2)*B(2,1) - A(1,3)*B(3,1)
      S(2) = S(2) - A(2,1)*B(1,2) - A(2,2)*B(2,2) - A(2,3)*B(3,2)
      S(3) = S(3) - A(3,1)*B(1,3) - A(3,2)*B(2,3) - A(3,3)*B(3,3)
      S(4) = S(4) - r1d2 * (A(1,1)*B(1,2) + A(1,2)*B(2,2)+A(1,3)*B(3,2)
     *                    + A(2,1)*B(1,1) + A(2,2)*B(2,1)+A(2,3)*B(3,1))
      S(5) = S(5) - r1d2 * (A(2,1)*B(1,3) + A(2,2)*B(2,3)+A(2,3)*B(3,3)
     *                    + A(3,1)*B(1,2) + A(3,2)*B(2,2)+A(3,3)*B(3,2))
      S(6) = S(6) - r1d2 * (A(1,1)*B(1,3) + A(1,2)*B(2,3)+A(1,3)*B(3,3)
     *                    + A(3,1)*B(1,1) + A(3,2)*B(2,1)+A(3,3)*B(3,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_sd_uu_af(S,A,B,fact)
      implicit none
      double precision A(3,3), B(3,3), S(6), fact, factd2
      factd2 = fact * 0.5d0
      S(1) = S(1) + fact   * (A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1))
      S(2) = S(2) + fact   * (A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2))
      S(3) = S(3) + fact   * (A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3))
      S(4) = S(4) + factd2 * (A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
     *                      + A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1))
      S(5) = S(5) + factd2 * (A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
     *                      + A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2))
      S(6) = S(6) + factd2 * (A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
     *                      + A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1))
      return
      end
c
c=============================================================================
c  Multiplikation von Tensoren   unsym = sym x sym
c-----------------------------------------------------------------------------
      subroutine mult_u_ss_p(S,A,B)
      implicit none
      double precision A(6), B(6), S(3,3)
      S(1,1) = A(1)*B(1) + A(4)*B(4) + A(6)*B(6)
      S(1,2) = A(1)*B(4) + A(4)*B(2) + A(6)*B(5)
      S(1,3) = A(1)*B(6) + A(4)*B(5) + A(6)*B(3)      
      S(2,1) = A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
      S(2,2) = A(4)*B(4) + A(2)*B(2) + A(5)*B(5)
      S(2,3) = A(4)*B(6) + A(2)*B(5) + A(5)*B(3)
      S(3,1) = A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
      S(3,2) = A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
      S(3,3) = A(6)*B(6) + A(5)*B(5) + A(3)*B(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_ss_m(S,A,B)
      implicit none
      double precision A(6), B(6), S(3,3)
      S(1,1) = - (A(1)*B(1) + A(4)*B(4) + A(6)*B(6))
      S(1,2) = - (A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(1,3) = - (A(1)*B(6) + A(4)*B(5) + A(6)*B(3))      
      S(2,1) = - (A(4)*B(1) + A(2)*B(4) + A(5)*B(6))
      S(2,2) = - (A(4)*B(4) + A(2)*B(2) + A(5)*B(5))
      S(2,3) = - (A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(3,1) = - (A(6)*B(1) + A(5)*B(4) + A(3)*B(6))
      S(3,2) = - (A(6)*B(4) + A(5)*B(2) + A(3)*B(5))
      S(3,3) = - (A(6)*B(6) + A(5)*B(5) + A(3)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_ss_f(S,A,B,fact)
      implicit none
      double precision A(6), B(6), S(3,3), fact
      S(1,1) = fact * (A(1)*B(1) + A(4)*B(4) + A(6)*B(6))
      S(1,2) = fact * (A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(1,3) = fact * (A(1)*B(6) + A(4)*B(5) + A(6)*B(3))      
      S(2,1) = fact * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6))
      S(2,2) = fact * (A(4)*B(4) + A(2)*B(2) + A(5)*B(5))
      S(2,3) = fact * (A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(3,1) = fact * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6))
      S(3,2) = fact * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5))
      S(3,3) = fact * (A(6)*B(6) + A(5)*B(5) + A(3)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_ss_ap(S,A,B)
      implicit none
      double precision A(6), B(6), S(3,3)
      S(1,1) = S(1,1) + A(1)*B(1) + A(4)*B(4) + A(6)*B(6)
      S(1,2) = S(1,2) + A(1)*B(4) + A(4)*B(2) + A(6)*B(5)
      S(1,3) = S(1,3) + A(1)*B(6) + A(4)*B(5) + A(6)*B(3)      
      S(2,1) = S(2,1) + A(4)*B(1) + A(2)*B(4) + A(5)*B(6)
      S(2,2) = S(2,2) + A(4)*B(4) + A(2)*B(2) + A(5)*B(5)
      S(2,3) = S(2,3) + A(4)*B(6) + A(2)*B(5) + A(5)*B(3)
      S(3,1) = S(3,1) + A(6)*B(1) + A(5)*B(4) + A(3)*B(6)
      S(3,2) = S(3,2) + A(6)*B(4) + A(5)*B(2) + A(3)*B(5)
      S(3,3) = S(3,3) + A(6)*B(6) + A(5)*B(5) + A(3)*B(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_ss_am(S,A,B)
      implicit none
      double precision A(6), B(6), S(3,3)
      S(1,1) = S(1,1) - A(1)*B(1) - A(4)*B(4) - A(6)*B(6)
      S(1,2) = S(1,2) - A(1)*B(4) - A(4)*B(2) - A(6)*B(5)
      S(1,3) = S(1,3) - A(1)*B(6) - A(4)*B(5) - A(6)*B(3)      
      S(2,1) = S(2,1) - A(4)*B(1) - A(2)*B(4) - A(5)*B(6)
      S(2,2) = S(2,2) - A(4)*B(4) - A(2)*B(2) - A(5)*B(5)
      S(2,3) = S(2,3) - A(4)*B(6) - A(2)*B(5) - A(5)*B(3)
      S(3,1) = S(3,1) - A(6)*B(1) - A(5)*B(4) - A(3)*B(6)
      S(3,2) = S(3,2) - A(6)*B(4) - A(5)*B(2) - A(3)*B(5)
      S(3,3) = S(3,3) - A(6)*B(6) - A(5)*B(5) - A(3)*B(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_ss_af(S,A,B,fact)
      implicit none
      double precision A(6), B(6), S(3,3), fact
      S(1,1) = S(1,1) + fact * (A(1)*B(1) + A(4)*B(4) + A(6)*B(6))
      S(1,2) = S(1,2) + fact * (A(1)*B(4) + A(4)*B(2) + A(6)*B(5))
      S(1,3) = S(1,3) + fact * (A(1)*B(6) + A(4)*B(5) + A(6)*B(3))      
      S(2,1) = S(2,1) + fact * (A(4)*B(1) + A(2)*B(4) + A(5)*B(6))
      S(2,2) = S(2,2) + fact * (A(4)*B(4) + A(2)*B(2) + A(5)*B(5))
      S(2,3) = S(2,3) + fact * (A(4)*B(6) + A(2)*B(5) + A(5)*B(3))
      S(3,1) = S(3,1) + fact * (A(6)*B(1) + A(5)*B(4) + A(3)*B(6))
      S(3,2) = S(3,2) + fact * (A(6)*B(4) + A(5)*B(2) + A(3)*B(5))
      S(3,3) = S(3,3) + fact * (A(6)*B(6) + A(5)*B(5) + A(3)*B(3))
      return
      end
c
c=============================================================================
c  Multiplikation von Tensoren     unsym = sym x unsym 
c-----------------------------------------------------------------------------
      subroutine mult_u_su_p(S,A,B)
      implicit none
      double precision A(6), B(3,3), S(3,3)
      S(1,1) = A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
      S(1,2) = A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2)
      S(1,3) = A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3)
      S(2,1) = A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
      S(2,2) = A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
      S(2,3) = A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3)
      S(3,1) = A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)
      S(3,2) = A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
      S(3,3) = A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_su_m(S,A,B)
      implicit none
      double precision A(6), B(3,3), S(3,3)
      S(1,1) = - (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(1,2) = - (A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(1,3) = - (A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      S(2,1) = - (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1))
      S(2,2) = - (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(2,3) = - (A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(3,1) = - (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1))
      S(3,2) = - (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2))
      S(3,3) = - (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_su_f(S,A,B,fact)
      implicit none
      double precision A(6), B(3,3), S(3,3), fact
      S(1,1) = fact * (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(1,2) = fact * (A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(1,3) = fact * (A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      S(2,1) = fact * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1))
      S(2,2) = fact * (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(2,3) = fact * (A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(3,1) = fact * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1))
      S(3,2) = fact * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2))
      S(3,3) = fact * (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_su_ap(S,A,B)
      implicit none
      double precision A(6), B(3,3), S(3,3)
      S(1,1) = S(1,1) + A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1)
      S(1,2) = S(1,2) + A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2)
      S(1,3) = S(1,3) + A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3)
      S(2,1) = S(2,1) + A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1)
      S(2,2) = S(2,2) + A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2)
      S(2,3) = S(2,3) + A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3)
      S(3,1) = S(3,1) + A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1)
      S(3,2) = S(3,2) + A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2)
      S(3,3) = S(3,3) + A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_su_am(S,A,B)
      implicit none
      double precision A(6), B(3,3), S(3,3)
      S(1,1) = S(1,1) - A(1)*B(1,1) - A(4)*B(2,1) - A(6)*B(3,1)
      S(1,2) = S(1,2) - A(1)*B(1,2) - A(4)*B(2,2) - A(6)*B(3,2)
      S(1,3) = S(1,3) - A(1)*B(1,3) - A(4)*B(2,3) - A(6)*B(3,3)
      S(2,1) = S(2,1) - A(4)*B(1,1) - A(2)*B(2,1) - A(5)*B(3,1)
      S(2,2) = S(2,2) - A(4)*B(1,2) - A(2)*B(2,2) - A(5)*B(3,2)
      S(2,3) = S(2,3) - A(4)*B(1,3) - A(2)*B(2,3) - A(5)*B(3,3)
      S(3,1) = S(3,1) - A(6)*B(1,1) - A(5)*B(2,1) - A(3)*B(3,1)
      S(3,2) = S(3,2) - A(6)*B(1,2) - A(5)*B(2,2) - A(3)*B(3,2)
      S(3,3) = S(3,3) - A(6)*B(1,3) - A(5)*B(2,3) - A(3)*B(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_su_af(S,A,B,fact)
      implicit none
      double precision A(6), B(3,3), S(3,3), fact
      S(1,1) = S(1,1) + fact * (A(1)*B(1,1) + A(4)*B(2,1) + A(6)*B(3,1))
      S(1,2) = S(1,2) + fact * (A(1)*B(1,2) + A(4)*B(2,2) + A(6)*B(3,2))
      S(1,3) = S(1,3) + fact * (A(1)*B(1,3) + A(4)*B(2,3) + A(6)*B(3,3))
      S(2,1) = S(2,1) + fact * (A(4)*B(1,1) + A(2)*B(2,1) + A(5)*B(3,1))
      S(2,2) = S(2,2) + fact * (A(4)*B(1,2) + A(2)*B(2,2) + A(5)*B(3,2))
      S(2,3) = S(2,3) + fact * (A(4)*B(1,3) + A(2)*B(2,3) + A(5)*B(3,3))
      S(3,1) = S(3,1) + fact * (A(6)*B(1,1) + A(5)*B(2,1) + A(3)*B(3,1))
      S(3,2) = S(3,2) + fact * (A(6)*B(1,2) + A(5)*B(2,2) + A(3)*B(3,2))
      S(3,3) = S(3,3) + fact * (A(6)*B(1,3) + A(5)*B(2,3) + A(3)*B(3,3))
      return
      end
c
c=============================================================================
c  Multiplikation von Tensoren     unsym = unsym x sym
c-----------------------------------------------------------------------------
      subroutine mult_u_us_p(S,A,B)
      implicit none
      double precision A(3,3), B(6), S(3,3)
      S(1,1) = A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6)
      S(1,2) = A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
      S(1,3) = A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)     
      S(2,1) = A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6)
      S(2,2) = A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5)
      S(2,3) = A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
      S(3,1) = A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6)
      S(3,2) = A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5)
      S(3,3) = A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_us_m(S,A,B)
      implicit none
      double precision A(3,3), B(6), S(3,3)
      S(1,1) = - (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(1,2) = - (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5))
      S(1,3) = - (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3))     
      S(2,1) = - (A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(2,2) = - (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(2,3) = - (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3))
      S(3,1) = - (A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      S(3,2) = - (A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(3,3) = - (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_us_f(S,A,B,fact)
      implicit none
      double precision A(3,3), B(6), S(3,3), fact
      S(1,1) = fact * (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(1,2) = fact * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5))
      S(1,3) = fact * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3))     
      S(2,1) = fact * (A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(2,2) = fact * (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(2,3) = fact * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3))
      S(3,1) = fact * (A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      S(3,2) = fact * (A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(3,3) = fact * (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_us_ap(S,A,B)
      implicit none
      double precision A(3,3), B(6), S(3,3)
      S(1,1) = S(1,1) + A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6)
      S(1,2) = S(1,2) + A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5)
      S(1,3) = S(1,3) + A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3)     
      S(2,1) = S(2,1) + A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6)
      S(2,2) = S(2,2) + A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5)
      S(2,3) = S(2,3) + A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3)
      S(3,1) = S(3,1) + A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6)
      S(3,2) = S(3,2) + A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5)
      S(3,3) = S(3,3) + A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_us_am(S,A,B)
      implicit none
      double precision A(3,3), B(6), S(3,3)
      S(1,1) = S(1,1) - A(1,1)*B(1) - A(1,2)*B(4) - A(1,3)*B(6)
      S(1,2) = S(1,2) - A(1,1)*B(4) - A(1,2)*B(2) - A(1,3)*B(5)
      S(1,3) = S(1,3) - A(1,1)*B(6) - A(1,2)*B(5) - A(1,3)*B(3)     
      S(2,1) = S(2,1) - A(2,1)*B(1) - A(2,2)*B(4) - A(2,3)*B(6)
      S(2,2) = S(2,2) - A(2,1)*B(4) - A(2,2)*B(2) - A(2,3)*B(5)
      S(2,3) = S(2,3) - A(2,1)*B(6) - A(2,2)*B(5) - A(2,3)*B(3)
      S(3,1) = S(3,1) - A(3,1)*B(1) - A(3,2)*B(4) - A(3,3)*B(6)
      S(3,2) = S(3,2) - A(3,1)*B(4) - A(3,2)*B(2) - A(3,3)*B(5)
      S(3,3) = S(3,3) - A(3,1)*B(6) - A(3,2)*B(5) - A(3,3)*B(3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_us_af(S,A,B,fact)
      implicit none
      double precision A(3,3), B(6), S(3,3), fact
      S(1,1) = S(1,1) + fact * (A(1,1)*B(1) + A(1,2)*B(4) + A(1,3)*B(6))
      S(1,2) = S(1,2) + fact * (A(1,1)*B(4) + A(1,2)*B(2) + A(1,3)*B(5))
      S(1,3) = S(1,3) + fact * (A(1,1)*B(6) + A(1,2)*B(5) + A(1,3)*B(3))     
      S(2,1) = S(2,1) + fact * (A(2,1)*B(1) + A(2,2)*B(4) + A(2,3)*B(6))
      S(2,2) = S(2,2) + fact * (A(2,1)*B(4) + A(2,2)*B(2) + A(2,3)*B(5))
      S(2,3) = S(2,3) + fact * (A(2,1)*B(6) + A(2,2)*B(5) + A(2,3)*B(3))
      S(3,1) = S(3,1) + fact * (A(3,1)*B(1) + A(3,2)*B(4) + A(3,3)*B(6))
      S(3,2) = S(3,2) + fact * (A(3,1)*B(4) + A(3,2)*B(2) + A(3,3)*B(5))
      S(3,3) = S(3,3) + fact * (A(3,1)*B(6) + A(3,2)*B(5) + A(3,3)*B(3))
      return
      end
c
c=============================================================================
c  Multiplikation von Tensoren     unsym = unsym x unsym
c-----------------------------------------------------------------------------
      subroutine mult_u_uu_p(S,A,B)
      implicit none
      integer          i, j, k
      double precision A(3,3), B(3,3), S(3,3), r0
      r0 = 0.d0
      do i=1, 3
        do j=1, 3
          S(i,j) = r0
          do k=1, 3
            S(i,j) = S(i,j) + A(i,k) * B(k,j) 
          enddo
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_uu_m(S,A,B)
      implicit none
      integer          i, j, k
      double precision A(3,3), B(3,3), S(3,3), r0
      r0 = 0.d0
      do i=1, 3
        do j=1, 3
          S(i,j) = r0
          do k=1, 3
            S(i,j) = S(i,j) - A(i,k) * B(k,j) 
          enddo
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_uu_f(S,A,B,fact)
      implicit none
      integer          i, j, k
      double precision A(3,3), B(3,3), S(3,3), fact, r0
      r0 = 0.d0
      do i=1, 3
        do j=1, 3
          S(i,j) = r0
          do k=1, 3
            S(i,j) = S(i,j) + A(i,k) * B(k,j) 
          enddo
          S(i,j) = fact * S(i,j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_uu_ap(S,A,B)
      implicit none
      integer          i, j, k
      double precision A(3,3), B(3,3), S(3,3)
      do i=1, 3
        do j=1, 3
          do k=1, 3
            S(i,j) = S(i,j) + A(i,k) * B(k,j) 
          enddo
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_uu_am(S,A,B)
      implicit none
      integer          i, j, k
      double precision A(3,3), B(3,3), S(3,3)
      do i=1, 3
        do j=1, 3
          do k=1, 3
            S(i,j) = S(i,j) - A(i,k) * B(k,j) 
          enddo
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_u_uu_af(S,A,B,fact)
      implicit none
      integer          i, j, k
      double precision A(3,3), B(3,3), S(3,3), H(3,3), fact, r0
      r0 = 0.d0
      do i=1, 3
        do j=1, 3
          H(i,j) = r0
          do k=1, 3
            H(i,j) = H(i,j) + A(i,k) * B(k,j) 
          enddo
          S(i,j) = S(i,j) + fact * H(i,j)
        enddo
      enddo
      return
      end
c
c=============================================================================
c  double contraction product eines symmetrischen zweistufigen Tensors 
c  mit einem Tensor 4.Stufe -> symmetrischen Tensor 2.Stufe
c-----------------------------------------------------------------------------
      subroutine mult_s_s4_p(S,A,X)
      implicit none
      integer i
      double precision A(6), X(6,6), S(6), r2
      r2 = 2.d0
      do i=1, 6
        S(i)   =          A(1)*X(1,i) + A(2)*X(2,i) + A(3)*X(3,i)
     *            + r2 * (A(4)*X(4,i) + A(5)*X(5,i) + A(6)*X(6,i))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_s4_m(S,A,X)
      implicit none
      integer i
      double precision A(6), X(6,6), S(6), r2
      r2 = 2.d0
      do i=1, 6
        S(i)   =        - A(1)*X(1,i) - A(2)*X(2,i) - A(3)*X(3,i)
     *            - r2 * (A(4)*X(4,i) + A(5)*X(5,i) + A(6)*X(6,i))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_s4_f(S,A,X,fact)
      implicit none
      integer i
      double precision A(6), X(6,6), S(6), fact, r2
      r2 = 2.d0
      do i=1, 6
        S(i)   = fact * (A(1)*X(1,i) + A(2)*X(2,i) + A(3)*X(3,i)
     *           + r2 * (A(4)*X(4,i) + A(5)*X(5,i) + A(6)*X(6,i)))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_s4_ap(S,A,X)
      implicit none
      integer i
      double precision A(6), X(6,6), S(6), r2
      r2 = 2.d0
      do i=1, 6
        S(i)   = S(i) +   A(1)*X(1,i) + A(2)*X(2,i) + A(3)*X(3,i)
     *            + r2 * (A(4)*X(4,i) + A(5)*X(5,i) + A(6)*X(6,i))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_s4_am(S,A,X)
      implicit none
      integer i
      double precision A(6), X(6,6), S(6), r2
      r2 = 2.d0
      do i=1, 6
        S(i)   = S(i) -   A(1)*X(1,i) - A(2)*X(2,i) - A(3)*X(3,i)
     *            - r2 * (A(4)*X(4,i) + A(5)*X(5,i) + A(6)*X(6,i))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_s4_af(S,A,X,fact)
      implicit none
      integer i
      double precision A(6), X(6,6), S(6), fact, r2
      r2 = 2.d0
      do i=1, 6
        S(i) = S(i) + fact* (A(1)*X(1,i) + A(2)*X(2,i) + A(3)*X(3,i)
     *               + r2 * (A(4)*X(4,i) + A(5)*X(5,i) + A(6)*X(6,i)))
      enddo
      return
      end
c
c=============================================================================
c  double contraction product eines Tensors 4.Stufe mit einem symmetrischen 
c  zweistufigen Tensor -> symmetrischen Tensor 2.Stufe
c-----------------------------------------------------------------------------
      subroutine mult_s_4s_p(S,X,A)
      implicit none
      integer i
      double precision X(6,6), A(6), S(6), r2
      r2 = 2.d0
      do i=1, 6
        S(i) =       X(i,1)*A(1) + X(i,2)*A(2) + X(i,3)*A(3)
     *       + r2 * (X(i,4)*A(4) + X(i,5)*A(5) + X(i,6)*A(6))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_4s_m(S,X,A)
      implicit none
      integer i
      double precision X(6,6), A(6), S(6), r2
      r2 = 2.d0
      do i=1, 6
        S(i) =     - X(i,1)*A(1) - X(i,2)*A(2) - X(i,3)*A(3)
     *       - r2 * (X(i,4)*A(4) + X(i,5)*A(5) + X(i,6)*A(6))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_4s_f(S,X,A,fact)
      implicit none
      integer i
      double precision X(6,6), A(6), S(6), fact, r2
      r2 = 2.d0
      do i=1, 6
        S(i) =  fact * ( X(i,1)*A(1) + X(i,2)*A(2) + X(i,3)*A(3)
     *           + r2 * (X(i,4)*A(4) + X(i,5)*A(5) + X(i,6)*A(6)))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_4s_ap(S,X,A)
      implicit none
      integer i
      double precision X(6,6), A(6), S(6), r2
      r2 = 2.d0
      do i=1, 6
        S(i) = S(i) + X(i,1)*A(1) + X(i,2)*A(2) + X(i,3)*A(3)
     *        + r2 * (X(i,4)*A(4) + X(i,5)*A(5) + X(i,6)*A(6))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_4s_am(S,X,A)
      implicit none
      integer i
      double precision X(6,6), A(6), S(6), r2
      r2 = 2.d0
      do i=1, 6
        S(i) = S(i) - X(i,1)*A(1) - X(i,2)*A(2) - X(i,3)*A(3)
     *        - r2 * (X(i,4)*A(4) + X(i,5)*A(5) + X(i,6)*A(6))
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_s_4s_af(S,X,A,fact)
      implicit none
      integer i
      double precision X(6,6), A(6), S(6), fact, r2
      r2 = 2.d0
      do i=1, 6
        S(i) = S(i) + fact* ( X(i,1)*A(1) + X(i,2)*A(2) + X(i,3)*A(3)
     *                + r2 * (X(i,4)*A(4) + X(i,5)*A(5) + X(i,6)*A(6)))
      enddo
      return
      end
c
c=============================================================================
c  Tensor hoch n            S = A^n
c-----------------------------------------------------------------------------
      subroutine power_s(S,A,n)
      implicit none 
      integer i, n
      double precision A(6), S(6), H(6)
      call set_s_p(S,A)
      do i=2, n
        call mult_s_ss_p(H,S,A)
        call set_s_p(S,H)
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine power_u(S,A,n)
      implicit none 
      integer i, n
      double precision A(3,3), S(3,3), H(3,3)
      call set_u_p(S,A)
      do i=2, n
        call mult_u_uu_p(H,S,A)
        call set_u_p(S,H)
      enddo
      return
      end
c
c=============================================================================
c   Inverse eines Tensors 2.Stufe        S = A^-1   
c-----------------------------------------------------------------------------
      subroutine inverse_s(S,A)
      implicit none
      integer i
      double precision A(6), S(6), detA
      S(1) = A(2)*A(3) - A(5)*A(5)
      S(4) = A(6)*A(5) - A(4)*A(3)
      S(6) = A(4)*A(5) - A(6)*A(2)
      S(2) = A(1)*A(3) - A(6)*A(6)
      S(5) = A(6)*A(4) - A(1)*A(5)
      S(3) = A(1)*A(2) - A(4)*A(4)
      detA = A(1)*S(1) + A(4)*S(4) + A(6)*S(6)
      do i=1, 6
        S(i) = S(i) / detA
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine inverse_u(S,A)
      implicit none
      integer i, j
      double precision A(3,3), S(3,3), detA
      S(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
      S(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
      S(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
      S(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
      S(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
      S(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
      S(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
      S(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
      S(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      detA   = A(1,1)*S(1,1) + A(1,2)*S(2,1) + A(1,3)*S(3,1)
      do i=1, 3
        do j=1, 3
          S(i,j) = S(i,j) / detA
        enddo
      enddo
      return
      end
c
c=============================================================================
c   Transponierte eines Tensors   S = A^T
c-----------------------------------------------------------------------------
      subroutine transpose(S,A)
      implicit none
      double precision A(3,3), S(3,3)
      S(1,1) = A(1,1)
      S(1,2) = A(2,1)
      S(1,3) = A(3,1)
      S(2,1) = A(1,2)
      S(2,2) = A(2,2)
      S(2,3) = A(3,2)
      S(3,1) = A(1,3)
      S(3,2) = A(2,3)
      S(3,3) = A(3,3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine transpose_r(S)
      implicit none
      double precision S(3,3), hlp
      hlp    = S(1,2)
      S(1,2) = S(2,1)
      S(2,1) = hlp
      hlp    = S(1,3)
      S(1,3) = S(3,1)
      S(3,1) = hlp
      hlp    = S(2,3)
      S(2,3) = S(3,2)
      S(3,2) = hlp
      return
      end
c
c=============================================================================
c   Deviator eines Tensors   S = A - tr(A)/3 I
c-----------------------------------------------------------------------------
      subroutine dev_u_r(S)
      implicit none
      double precision S(3,3), fact
      fact = (S(1,1) + S(2,2) + S(3,3)) / 3.d0
      S(1,1) = S(1,1) - fact
      S(2,2) = S(2,2) - fact
      S(3,3) = S(3,3) - fact
      return
      end
c-----------------------------------------------------------------------------
      subroutine dev_s_r(S)
      implicit none
      double precision S(6), fact
      fact = (S(1) + S(2) + S(3)) / 3.d0
      S(1) = S(1) - fact
      S(2) = S(2) - fact
      S(3) = S(3) - fact
      return
      end
c
cc=============================================================================
c  Differentiation
c
c         D sqrt{dev(AT):[dev(AT)]^T}           mit A,T sym., AT unsym. 2.Stufe
c    S = ---------------------------------      
c                    D T                    -> sym. Tensor 2.Stufe
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATd_p(S,A,AT,n2s)
      implicit none
      double precision S(6), A(6), AT(3,3), n2s, fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      fact = 1.d0 / n2s
      call set_s_f(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATd_m(S,A,AT,n2s)
      implicit none
      double precision S(6), A(6), AT(3,3), n2s, fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      fact = - 1.d0 / n2s
      call set_s_f(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATd_f(S,A,AT,n2s,fact)
      implicit none
      double precision S(6), A(6), AT(3,3), n2s, fact, fact2, hs(6)
      fact2 = - (AT(1,1)+AT(2,2)+AT(3,3)) / 3.d0
      call set_s_f(hs,A,fact2)
      call mult_s_us_ap(hs,AT,A)
      fact2 = fact / n2s
      call set_s_f(S,hs,fact2)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATd_ap(S,A,AT,n2s)
      implicit none
      double precision S(6), A(6), AT(3,3), n2s, fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      fact = 1.d0 / n2s
      call set_s_af(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATd_am(S,A,AT,n2s)
      implicit none
      double precision S(6), A(6), AT(3,3), n2s, fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      fact = - 1.d0 / n2s
      call set_s_af(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATd_af(S,A,AT,n2s,fact)
      implicit none
      double precision S(6), A(6), AT(3,3), n2s, fact, fact2, hs(6)
      fact2 = - (AT(1,1)+AT(2,2)+AT(3,3)) / 3.d0
      call set_s_f(hs,A,fact2)
      call mult_s_us_ap(hs,AT,A)
      fact2 = fact / n2s
      call set_s_af(S,hs,fact2)
      return
      end
c
cc=============================================================================
c  Differentiation
c
c         D sqrt{dev(AT+B):[dev(AT+B)]^T}     mit A,T sym., B unsym. 2.Stufe
c    S = ---------------------------------      
c                      D T                    -> sym. Tensor 2.Stufe
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATpUd_p(S,A,B,AT,n2s)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, 
     *                 fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      call sympart_UA_ap(hs,B,A)
      fact = 1.d0 / n2s
      call set_s_f(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATpUd_m(S,A,B,AT,n2s)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, 
     *                 fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      call sympart_UA_ap(hs,B,A)
      fact = - 1.d0 / n2s
      call set_s_f(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATpUd_f(S,A,B,AT,n2s,fact)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, fact,
     *                 fact2, hs(6)
      fact2 = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact2)
      call mult_s_us_ap(hs,AT,A)
      call sympart_UA_ap(hs,B,A)
      fact2 = fact / n2s
      call set_s_f(S,hs,fact2)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATpUd_ap(S,A,B,AT,n2s)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, 
     *                 fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      call sympart_UA_ap(hs,B,A)
      fact = 1.d0 / n2s
      call set_s_af(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATpUd_am(S,A,B,AT,n2s)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, 
     *                 fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      call sympart_UA_ap(hs,B,A)
      fact = - 1.d0 / n2s
      call set_s_af(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sATpUd_af(S,A,B,AT,n2s,fact)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, fact,
     *                 fact2, hs(6)
      fact2 = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact2)
      call mult_s_us_ap(hs,AT,A)
      call sympart_UA_ap(hs,B,A)
      fact2 = fact / n2s
      call set_s_af(S,hs,fact2)
      return
      end
c
c=============================================================================
c  Differentiation
c
c         D sqrt{dev(TA+B):[dev(TA+B)]^T}     mit A,T sym., B unsym. 2.Stufe
c    S = ---------------------------------      
c                      D T                    -> sym. Tensor 2.Stufe
c-----------------------------------------------------------------------------
      subroutine diff_norm2sTApUd_p(S,A,B,AT,n2s)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, 
     *                 fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      call sympart_AU_ap(hs,A,B)
      fact = 1.d0 / n2s
      call set_s_f(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sTApUd_m(S,A,B,AT,n2s)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, 
     *                 fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      call sympart_AU_ap(hs,A,B)
      fact = - 1.d0 / n2s
      call set_s_f(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sTApUd_f(S,A,B,AT,n2s,fact)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, fact,
     *                 fact2, hs(6), hu(3,3)
      fact2 = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact2)
      call mult_s_us_ap(hs,AT,A)
      call sympart_AU_ap(hs,A,B)
      fact2 = fact / n2s
      call set_s_f(S,hs,fact2)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sTApUd_ap(S,A,B,AT,n2s)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, 
     *                 fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      call sympart_AU_ap(hs,A,B)
      fact = 1.d0 / n2s
      call set_s_af(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sTApUd_am(S,A,B,AT,n2s)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, 
     *                 fact, hs(6)
      fact = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact)
      call mult_s_us_ap(hs,AT,A)
      call sympart_AU_ap(hs,A,B)
      fact = - 1.d0 / n2s
      call set_s_af(S,hs,fact)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_norm2sTApUd_af(S,A,B,AT,n2s,fact)
      implicit none
      double precision S(6), A(6), B(3,3), AT(3,3), n2s, fact,
     *                 fact2, hs(6)
      fact2 = - (AT(1,1)+AT(2,2)+AT(3,3) + B(1,1)+B(2,2)+B(3,3)) / 3.d0
      call set_s_f(hs,A,fact2)
      call mult_s_us_ap(hs,AT,A)
      call sympart_AU_ap(hs,A,B)
      fact2 = fact / n2s
      call set_s_af(S,hs,fact2)
      return
      end
c
c=============================================================================
c*****************************************************************************
c*****************************************************************************

c=============================================================================
c  Zuweisung   X := A
c-----------------------------------------------------------------------------
      subroutine set_4_p(X,A)
      implicit none
      integer i, j
      double precision X(6,6), A(6,6)
      do i=1, 6
        do j=1, 6
          X(i,j) = A(i,j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_4_m(X,A)
      implicit none
      integer i, j
      double precision A(6,6), X(6,6)
      do i=1, 6
        do j=1, 6
          X(i,j) = - A(i,j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_4_f(X,A,fact)
      implicit none
      integer i, j
      double precision A(6,6), X(6,6), fact
      do i=1, 6
        do j=1, 6
          X(i,j) = fact * A(i,j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_4_ap(X,A)
      implicit none
      integer i, j
      double precision A(6,6), X(6,6)
      do i=1, 6
        do j=1, 6
          X(i,j) = X(i,j) + A(i,j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_4_am(X,A)
      implicit none
      integer i, j
      double precision A(6,6), X(6,6)
      do i=1, 6
        do j=1, 6
          X(i,j) = X(i,j) - A(i,j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine set_4_af(X,A,fact)
      implicit none
      integer i, j
      double precision A(6,6), X(6,6), fact
      do i=1, 6
        do j=1, 6
          X(i,j) = X(i,j) + fact * A(i,j)
        enddo
      enddo
      return
      end
c
c=============================================================================
c  Einheitstensor 4. Stufe  X_ijkl = delta_ik delta_jl
c-----------------------------------------------------------------------------
      subroutine unit_4_p(X)
      implicit none
      double precision X(6,6)
      call pzero(X,36)
      X(1,1) = 1.d0
      X(2,2) = X(1,1)
      X(3,3) = X(1,1)
      X(4,4) = 0.5d0
      X(5,5) = X(4,4)
      X(6,6) = X(4,4)
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_4_m(X)
      implicit none
      double precision X(6,6)
      call pzero(X,36)
      X(1,1) = - 1.d0
      X(2,2) = X(1,1)
      X(3,3) = X(1,1)
      X(4,4) = - 0.5d0
      X(5,5) = X(4,4)
      X(6,6) = X(4,4)
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_4_f(X,fact)
      implicit none
      double precision X(6,6), fact
      call pzero(X,36)
      X(1,1) = fact
      X(2,2) = fact
      X(3,3) = fact
      X(4,4) = 0.5d0 * fact
      X(5,5) = X(4,4)
      X(6,6) = X(4,4)
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_4_ap(X)
      implicit none
      double precision X(6,6), r1, r1d2
      data r1d2, r1 / 0.5d0, 1.d0 /
      X(1,1) = X(1,1) + r1
      X(2,2) = X(2,2) + r1
      X(3,3) = X(3,3) + r1
      X(4,4) = X(4,4) + r1d2
      X(5,5) = X(5,5) + r1d2
      X(6,6) = X(6,6) + r1d2
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_4_am(X)
      implicit none
      double precision X(6,6), r1, r1d2
      data r1d2, r1 / 0.5d0, 1.d0 /
      X(1,1) = X(1,1) - r1
      X(2,2) = X(2,2) - r1
      X(3,3) = X(3,3) - r1
      X(4,4) = X(4,4) - r1d2
      X(5,5) = X(5,5) - r1d2
      X(6,6) = X(6,6) - r1d2
      return
      end
c-----------------------------------------------------------------------------
      subroutine unit_4_af(X,fact)
      implicit none
      double precision X(6,6), fact, factd2
      factd2 = fact * 0.5d0
      X(1,1) = X(1,1) + fact
      X(2,2) = X(2,2) + fact
      X(3,3) = X(3,3) + fact
      X(4,4) = X(4,4) + factd2
      X(5,5) = X(5,5) + factd2
      X(6,6) = X(6,6) + factd2
      return
      end
c
c=============================================================================
c  Dyadisches Produkt     X = A o B     (-> vierstufiger Tensor)
c-----------------------------------------------------------------------------
      subroutine mult_4_ss_p(X,A,B)
      implicit none 
      integer i,j
      double precision A(6), B(6), X(6,6)
      do i=1, 6 
        do j=1, 6
          X(i,j) = A(i) * B(j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_ss_m(X,A,B)
      implicit none 
      integer i,j
      double precision A(6), B(6), X(6,6)
      do i=1, 6 
        do j=1, 6
          X(i,j) = - A(i) * B(j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_ss_f(X,A,B,fact)
      implicit none 
      integer i,j
      double precision A(6), B(6), X(6,6), fact
      do i=1, 6 
        do j=1, 6
          X(i,j) = fact * A(i) * B(j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_ss_ap(X,A,B)
      implicit none 
      integer i,j
      double precision A(6), B(6), X(6,6)
      do i=1, 6 
        do j=1, 6
          X(i,j) = X(i,j) + A(i) * B(j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_ss_am(X,A,B)
      implicit none 
      integer i,j
      double precision A(6), B(6), X(6,6)
      do i=1, 6 
        do j=1, 6
          X(i,j) = X(i,j) - A(i) * B(j)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_ss_af(X,A,B,fact)
      implicit none 
      integer i,j
      double precision A(6), B(6), X(6,6), fact
      do i=1, 6 
        do j=1, 6
          X(i,j) = X(i,j) + fact * A(i) * B(j)
        enddo
      enddo
      return
      end
c
c=============================================================================
c  double contraction product zweier vierstufiger Tensoren    X = A : B
c-----------------------------------------------------------------------------
      subroutine mult_4_44_p(X,A,B)
      implicit none 
      integer i, j
      double precision A(6,6), B(6,6), X(6,6), r2
      r2 = 2.d0
      do i=1, 6
        do j=1, 6
          X(i,j) =   A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) 
     *       + r2 * (A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j)) 
        enddo
      enddo 
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_44_m(X,A,B)
      implicit none 
      integer i, j
      double precision A(6,6), B(6,6), X(6,6), r2
      r2 = 2.d0
      do i=1, 6
        do j=1, 6
          X(i,j) = - A(i,1)*B(1,j) - A(i,2)*B(2,j) - A(i,3)*B(3,j) 
     *       - r2 * (A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j)) 
        enddo
      enddo 
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_44_f(X,A,B,fact)
      implicit none 
      integer i, j
      double precision A(6,6), B(6,6), X(6,6), fact, r2
      r2 = 2.d0
      do i=1, 6
        do j=1, 6
          X(i,j) = fact * 
     *             ( A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) 
     *       + r2 * (A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j))) 
        enddo
      enddo 
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_44_ap(X,A,B)
      implicit none 
      integer i, j
      double precision A(6,6), B(6,6), X(6,6), r2
      r2 = 2.d0
      do i=1, 6
        do j=1, 6
          X(i,j) = X(i,j) + 
     *                A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) 
     *        + r2 * (A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j)) 
        enddo
      enddo 
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_44_am(X,A,B)
      implicit none 
      integer i, j
      double precision A(6,6), B(6,6), X(6,6), r2
      r2 = 2.d0
      do i=1, 6
        do j=1, 6
          X(i,j) = X(i,j) 
     *              - A(i,1)*B(1,j) - A(i,2)*B(2,j) - A(i,3)*B(3,j) 
     *        - r2 * (A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j)) 
        enddo
      enddo 
      return 
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_44_af(X,A,B,fact)
      implicit none 
      integer i, j
      double precision A(6,6), B(6,6), X(6,6), fact, r2
      r2 = 2.d0
      do i=1, 6
        do j=1, 6
          X(i,j) = X(i,j) + fact * 
     *             ( A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) 
     *       + r2 * (A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j))) 
        enddo
      enddo 
      return 
      end
c
c=============================================================================
c  single contraction product eines vierstufigen mit einem zweistufingen 
c  Tensor     X = A B
c    beachte Symmetrieeigenschaften: geht nur, wenn das Ergebnis sich 
c    wieder als 6x6 Matrix speichern laesst!
c-----------------------------------------------------------------------------
      subroutine mult_4_4s_p(X,A,B)
      implicit none
      double precision X(6,*), A(6,*), B(*), r1d2
      r1d2 = 0.5d0
      X(1,1) = A(1,1)*B(1) + A(1,4)*B(4) + A(1,6)*B(6)
      X(2,1) = A(2,1)*B(1) + A(2,4)*B(4) + A(2,6)*B(6)
      X(3,1) = A(3,1)*B(1) + A(3,4)*B(4) + A(3,6)*B(6)
      X(4,1) = A(4,1)*B(1) + A(4,4)*B(4) + A(4,6)*B(6)
      X(5,1) = A(5,1)*B(1) + A(5,4)*B(4) + A(5,6)*B(6)
      X(6,1) = A(6,1)*B(1) + A(6,4)*B(4) + A(6,6)*B(6)
      X(1,2) = A(1,4)*B(4) + A(1,2)*B(2) + A(1,5)*B(5)
      X(2,2) = A(2,4)*B(4) + A(2,2)*B(2) + A(2,5)*B(5)
      X(3,2) = A(3,4)*B(4) + A(3,2)*B(2) + A(3,5)*B(5)
      X(4,2) = A(4,4)*B(4) + A(4,2)*B(2) + A(4,5)*B(5)
      X(5,2) = A(5,4)*B(4) + A(5,2)*B(2) + A(5,5)*B(5)
      X(6,2) = A(6,4)*B(4) + A(6,2)*B(2) + A(6,5)*B(5)
      X(1,3) = A(1,6)*B(6) + A(1,5)*B(5) + A(1,3)*B(3)
      X(2,3) = A(2,6)*B(6) + A(2,5)*B(5) + A(2,3)*B(3)
      X(3,3) = A(3,6)*B(6) + A(3,5)*B(5) + A(3,3)*B(3)
      X(4,3) = A(4,6)*B(6) + A(4,5)*B(5) + A(4,3)*B(3)
      X(5,3) = A(5,6)*B(6) + A(5,5)*B(5) + A(5,3)*B(3)
      X(6,3) = A(6,6)*B(6) + A(6,5)*B(5) + A(6,3)*B(3)
      X(1,4) = r1d2 * (A(1,1)*B(4) + A(1,4)*B(2) + A(1,6)*B(5)
     *               + A(1,4)*B(1) + A(1,2)*B(4) + A(1,5)*B(6))
      X(2,4) = r1d2 * (A(2,1)*B(4) + A(2,4)*B(2) + A(2,6)*B(5)
     *               + A(2,4)*B(1) + A(2,2)*B(4) + A(2,5)*B(6))
      X(3,4) = r1d2 * (A(3,1)*B(4) + A(3,4)*B(2) + A(3,6)*B(5)
     *               + A(3,4)*B(1) + A(3,2)*B(4) + A(3,5)*B(6))
      X(4,4) = r1d2 * (A(4,1)*B(4) + A(4,4)*B(2) + A(4,6)*B(5)
     *               + A(4,4)*B(1) + A(4,2)*B(4) + A(4,5)*B(6))
      X(5,4) = r1d2 * (A(5,1)*B(4) + A(5,4)*B(2) + A(5,6)*B(5)
     *               + A(5,4)*B(1) + A(5,2)*B(4) + A(5,5)*B(6))
      X(6,4) = r1d2 * (A(6,1)*B(4) + A(6,4)*B(2) + A(6,6)*B(5)
     *               + A(6,4)*B(1) + A(6,2)*B(4) + A(6,5)*B(6))
      X(1,5) = r1d2 * (A(1,4)*B(6) + A(1,2)*B(5) + A(1,5)*B(3)
     *               + A(1,6)*B(4) + A(1,5)*B(2) + A(1,3)*B(5))
      X(2,5) = r1d2 * (A(2,4)*B(6) + A(2,2)*B(5) + A(2,5)*B(3)
     *               + A(2,6)*B(4) + A(2,5)*B(2) + A(2,3)*B(5))
      X(3,5) = r1d2 * (A(3,4)*B(6) + A(3,2)*B(5) + A(3,5)*B(3)
     *               + A(3,6)*B(4) + A(3,5)*B(2) + A(3,3)*B(5))
      X(4,5) = r1d2 * (A(4,4)*B(6) + A(4,2)*B(5) + A(4,5)*B(3)
     *               + A(4,6)*B(4) + A(4,5)*B(2) + A(4,3)*B(5))
      X(5,5) = r1d2 * (A(5,4)*B(6) + A(5,2)*B(5) + A(5,5)*B(3)
     *               + A(5,6)*B(4) + A(5,5)*B(2) + A(5,3)*B(5))
      X(6,5) = r1d2 * (A(6,4)*B(6) + A(6,2)*B(5) + A(6,5)*B(3)
     *               + A(6,6)*B(4) + A(6,5)*B(2) + A(6,3)*B(5))
      X(1,6) = r1d2 * (A(1,6)*B(1) + A(1,5)*B(4) + A(1,3)*B(6)
     *               + A(1,1)*B(6) + A(1,4)*B(5) + A(1,6)*B(3))
      X(2,6) = r1d2 * (A(2,6)*B(1) + A(2,5)*B(4) + A(2,3)*B(6)
     *               + A(2,1)*B(6) + A(2,4)*B(5) + A(2,6)*B(3))
      X(3,6) = r1d2 * (A(3,6)*B(1) + A(3,5)*B(4) + A(3,3)*B(6)
     *               + A(3,1)*B(6) + A(3,4)*B(5) + A(3,6)*B(3))
      X(4,6) = r1d2 * (A(4,6)*B(1) + A(4,5)*B(4) + A(4,3)*B(6)
     *               + A(4,1)*B(6) + A(4,4)*B(5) + A(4,6)*B(3))
      X(5,6) = r1d2 * (A(5,6)*B(1) + A(5,5)*B(4) + A(5,3)*B(6)
     *               + A(5,1)*B(6) + A(5,4)*B(5) + A(5,6)*B(3))
      X(6,6) = r1d2 * (A(6,6)*B(1) + A(6,5)*B(4) + A(6,3)*B(6)
     *               + A(6,1)*B(6) + A(6,4)*B(5) + A(6,6)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_4s_m(X,A,B)
      implicit none
      double precision X(6,*), A(6,*), B(*), m1d2
      m1d2 = - 0.5d0
      X(1,1) = - A(1,1)*B(1) - A(1,4)*B(4) - A(1,6)*B(6)
      X(2,1) = - A(2,1)*B(1) - A(2,4)*B(4) - A(2,6)*B(6)
      X(3,1) = - A(3,1)*B(1) - A(3,4)*B(4) - A(3,6)*B(6)
      X(4,1) = - A(4,1)*B(1) - A(4,4)*B(4) - A(4,6)*B(6)
      X(5,1) = - A(5,1)*B(1) - A(5,4)*B(4) - A(5,6)*B(6)
      X(6,1) = - A(6,1)*B(1) - A(6,4)*B(4) - A(6,6)*B(6)
      X(1,2) = - A(1,4)*B(4) - A(1,2)*B(2) - A(1,5)*B(5)
      X(2,2) = - A(2,4)*B(4) - A(2,2)*B(2) - A(2,5)*B(5)
      X(3,2) = - A(3,4)*B(4) - A(3,2)*B(2) - A(3,5)*B(5)
      X(4,2) = - A(4,4)*B(4) - A(4,2)*B(2) - A(4,5)*B(5)
      X(5,2) = - A(5,4)*B(4) - A(5,2)*B(2) - A(5,5)*B(5)
      X(6,2) = - A(6,4)*B(4) - A(6,2)*B(2) - A(6,5)*B(5)
      X(1,3) = - A(1,6)*B(6) - A(1,5)*B(5) - A(1,3)*B(3)
      X(2,3) = - A(2,6)*B(6) - A(2,5)*B(5) - A(2,3)*B(3)
      X(3,3) = - A(3,6)*B(6) - A(3,5)*B(5) - A(3,3)*B(3)
      X(4,3) = - A(4,6)*B(6) - A(4,5)*B(5) - A(4,3)*B(3)
      X(5,3) = - A(5,6)*B(6) - A(5,5)*B(5) - A(5,3)*B(3)
      X(6,3) = - A(6,6)*B(6) - A(6,5)*B(5) - A(6,3)*B(3)
      X(1,4) = m1d2 * (A(1,1)*B(4) + A(1,4)*B(2) + A(1,6)*B(5)
     *               + A(1,4)*B(1) + A(1,2)*B(4) + A(1,5)*B(6))
      X(2,4) = m1d2 * (A(2,1)*B(4) + A(2,4)*B(2) + A(2,6)*B(5)
     *               + A(2,4)*B(1) + A(2,2)*B(4) + A(2,5)*B(6))
      X(3,4) = m1d2 * (A(3,1)*B(4) + A(3,4)*B(2) + A(3,6)*B(5)
     *               + A(3,4)*B(1) + A(3,2)*B(4) + A(3,5)*B(6))
      X(4,4) = m1d2 * (A(4,1)*B(4) + A(4,4)*B(2) + A(4,6)*B(5)
     *               + A(4,4)*B(1) + A(4,2)*B(4) + A(4,5)*B(6))
      X(5,4) = m1d2 * (A(5,1)*B(4) + A(5,4)*B(2) + A(5,6)*B(5)
     *               + A(5,4)*B(1) + A(5,2)*B(4) + A(5,5)*B(6))
      X(6,4) = m1d2 * (A(6,1)*B(4) + A(6,4)*B(2) + A(6,6)*B(5)
     *               + A(6,4)*B(1) + A(6,2)*B(4) + A(6,5)*B(6))
      X(1,5) = m1d2 * (A(1,4)*B(6) + A(1,2)*B(5) + A(1,5)*B(3)
     *               + A(1,6)*B(4) + A(1,5)*B(2) + A(1,3)*B(5))
      X(2,5) = m1d2 * (A(2,4)*B(6) + A(2,2)*B(5) + A(2,5)*B(3)
     *               + A(2,6)*B(4) + A(2,5)*B(2) + A(2,3)*B(5))
      X(3,5) = m1d2 * (A(3,4)*B(6) + A(3,2)*B(5) + A(3,5)*B(3)
     *               + A(3,6)*B(4) + A(3,5)*B(2) + A(3,3)*B(5))
      X(4,5) = m1d2 * (A(4,4)*B(6) + A(4,2)*B(5) + A(4,5)*B(3)
     *               + A(4,6)*B(4) + A(4,5)*B(2) + A(4,3)*B(5))
      X(5,5) = m1d2 * (A(5,4)*B(6) + A(5,2)*B(5) + A(5,5)*B(3)
     *               + A(5,6)*B(4) + A(5,5)*B(2) + A(5,3)*B(5))
      X(6,5) = m1d2 * (A(6,4)*B(6) + A(6,2)*B(5) + A(6,5)*B(3)
     *               + A(6,6)*B(4) + A(6,5)*B(2) + A(6,3)*B(5))
      X(1,6) = m1d2 * (A(1,6)*B(1) + A(1,5)*B(4) + A(1,3)*B(6)
     *               + A(1,1)*B(6) + A(1,4)*B(5) + A(1,6)*B(3))
      X(2,6) = m1d2 * (A(2,6)*B(1) + A(2,5)*B(4) + A(2,3)*B(6)
     *               + A(2,1)*B(6) + A(2,4)*B(5) + A(2,6)*B(3))
      X(3,6) = m1d2 * (A(3,6)*B(1) + A(3,5)*B(4) + A(3,3)*B(6)
     *               + A(3,1)*B(6) + A(3,4)*B(5) + A(3,6)*B(3))
      X(4,6) = m1d2 * (A(4,6)*B(1) + A(4,5)*B(4) + A(4,3)*B(6)
     *               + A(4,1)*B(6) + A(4,4)*B(5) + A(4,6)*B(3))
      X(5,6) = m1d2 * (A(5,6)*B(1) + A(5,5)*B(4) + A(5,3)*B(6)
     *               + A(5,1)*B(6) + A(5,4)*B(5) + A(5,6)*B(3))
      X(6,6) = m1d2 * (A(6,6)*B(1) + A(6,5)*B(4) + A(6,3)*B(6)
     *               + A(6,1)*B(6) + A(6,4)*B(5) + A(6,6)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_4s_f(X,A,B,fact)
      implicit none
      double precision X(6,*), A(6,*), B(*), fact, factd2
      factd2 = 0.5d0 * fact
      X(1,1) = fact * (A(1,1)*B(1) + A(1,4)*B(4) + A(1,6)*B(6))
      X(2,1) = fact * (A(2,1)*B(1) + A(2,4)*B(4) + A(2,6)*B(6))
      X(3,1) = fact * (A(3,1)*B(1) + A(3,4)*B(4) + A(3,6)*B(6))
      X(4,1) = fact * (A(4,1)*B(1) + A(4,4)*B(4) + A(4,6)*B(6))
      X(5,1) = fact * (A(5,1)*B(1) + A(5,4)*B(4) + A(5,6)*B(6))
      X(6,1) = fact * (A(6,1)*B(1) + A(6,4)*B(4) + A(6,6)*B(6))
      X(1,2) = fact * (A(1,4)*B(4) + A(1,2)*B(2) + A(1,5)*B(5))
      X(2,2) = fact * (A(2,4)*B(4) + A(2,2)*B(2) + A(2,5)*B(5))
      X(3,2) = fact * (A(3,4)*B(4) + A(3,2)*B(2) + A(3,5)*B(5))
      X(4,2) = fact * (A(4,4)*B(4) + A(4,2)*B(2) + A(4,5)*B(5))
      X(5,2) = fact * (A(5,4)*B(4) + A(5,2)*B(2) + A(5,5)*B(5))
      X(6,2) = fact * (A(6,4)*B(4) + A(6,2)*B(2) + A(6,5)*B(5))
      X(1,3) = fact * (A(1,6)*B(6) + A(1,5)*B(5) + A(1,3)*B(3))
      X(2,3) = fact * (A(2,6)*B(6) + A(2,5)*B(5) + A(2,3)*B(3))
      X(3,3) = fact * (A(3,6)*B(6) + A(3,5)*B(5) + A(3,3)*B(3))
      X(4,3) = fact * (A(4,6)*B(6) + A(4,5)*B(5) + A(4,3)*B(3))
      X(5,3) = fact * (A(5,6)*B(6) + A(5,5)*B(5) + A(5,3)*B(3))
      X(6,3) = fact * (A(6,6)*B(6) + A(6,5)*B(5) + A(6,3)*B(3))
      X(1,4) = factd2 * (A(1,1)*B(4) + A(1,4)*B(2) + A(1,6)*B(5)
     *                 + A(1,4)*B(1) + A(1,2)*B(4) + A(1,5)*B(6))
      X(2,4) = factd2 * (A(2,1)*B(4) + A(2,4)*B(2) + A(2,6)*B(5)
     *                 + A(2,4)*B(1) + A(2,2)*B(4) + A(2,5)*B(6))
      X(3,4) = factd2 * (A(3,1)*B(4) + A(3,4)*B(2) + A(3,6)*B(5)
     *                 + A(3,4)*B(1) + A(3,2)*B(4) + A(3,5)*B(6))
      X(4,4) = factd2 * (A(4,1)*B(4) + A(4,4)*B(2) + A(4,6)*B(5)
     *                 + A(4,4)*B(1) + A(4,2)*B(4) + A(4,5)*B(6))
      X(5,4) = factd2 * (A(5,1)*B(4) + A(5,4)*B(2) + A(5,6)*B(5)
     *                 + A(5,4)*B(1) + A(5,2)*B(4) + A(5,5)*B(6))
      X(6,4) = factd2 * (A(6,1)*B(4) + A(6,4)*B(2) + A(6,6)*B(5)
     *                 + A(6,4)*B(1) + A(6,2)*B(4) + A(6,5)*B(6))
      X(1,5) = factd2 * (A(1,4)*B(6) + A(1,2)*B(5) + A(1,5)*B(3)
     *                 + A(1,6)*B(4) + A(1,5)*B(2) + A(1,3)*B(5))
      X(2,5) = factd2 * (A(2,4)*B(6) + A(2,2)*B(5) + A(2,5)*B(3)
     *                 + A(2,6)*B(4) + A(2,5)*B(2) + A(2,3)*B(5))
      X(3,5) = factd2 * (A(3,4)*B(6) + A(3,2)*B(5) + A(3,5)*B(3)
     *                 + A(3,6)*B(4) + A(3,5)*B(2) + A(3,3)*B(5))
      X(4,5) = factd2 * (A(4,4)*B(6) + A(4,2)*B(5) + A(4,5)*B(3)
     *                 + A(4,6)*B(4) + A(4,5)*B(2) + A(4,3)*B(5))
      X(5,5) = factd2 * (A(5,4)*B(6) + A(5,2)*B(5) + A(5,5)*B(3)
     *                 + A(5,6)*B(4) + A(5,5)*B(2) + A(5,3)*B(5))
      X(6,5) = factd2 * (A(6,4)*B(6) + A(6,2)*B(5) + A(6,5)*B(3)
     *                 + A(6,6)*B(4) + A(6,5)*B(2) + A(6,3)*B(5))
      X(1,6) = factd2 * (A(1,6)*B(1) + A(1,5)*B(4) + A(1,3)*B(6)
     *                 + A(1,1)*B(6) + A(1,4)*B(5) + A(1,6)*B(3))
      X(2,6) = factd2 * (A(2,6)*B(1) + A(2,5)*B(4) + A(2,3)*B(6)
     *                 + A(2,1)*B(6) + A(2,4)*B(5) + A(2,6)*B(3))
      X(3,6) = factd2 * (A(3,6)*B(1) + A(3,5)*B(4) + A(3,3)*B(6)
     *                 + A(3,1)*B(6) + A(3,4)*B(5) + A(3,6)*B(3))
      X(4,6) = factd2 * (A(4,6)*B(1) + A(4,5)*B(4) + A(4,3)*B(6)
     *                 + A(4,1)*B(6) + A(4,4)*B(5) + A(4,6)*B(3))
      X(5,6) = factd2 * (A(5,6)*B(1) + A(5,5)*B(4) + A(5,3)*B(6)
     *                 + A(5,1)*B(6) + A(5,4)*B(5) + A(5,6)*B(3))
      X(6,6) = factd2 * (A(6,6)*B(1) + A(6,5)*B(4) + A(6,3)*B(6)
     *                 + A(6,1)*B(6) + A(6,4)*B(5) + A(6,6)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_4s_ap(X,A,B)
      implicit none
      double precision X(6,*), A(6,*), B(*), r1d2
      r1d2 = 0.5d0
      X(1,1) = X(1,1) + A(1,1)*B(1) + A(1,4)*B(4) + A(1,6)*B(6)
      X(2,1) = X(2,1) + A(2,1)*B(1) + A(2,4)*B(4) + A(2,6)*B(6)
      X(3,1) = X(3,1) + A(3,1)*B(1) + A(3,4)*B(4) + A(3,6)*B(6)
      X(4,1) = X(4,1) + A(4,1)*B(1) + A(4,4)*B(4) + A(4,6)*B(6)
      X(5,1) = X(5,1) + A(5,1)*B(1) + A(5,4)*B(4) + A(5,6)*B(6)
      X(6,1) = X(6,1) + A(6,1)*B(1) + A(6,4)*B(4) + A(6,6)*B(6)
      X(1,2) = X(1,2) + A(1,4)*B(4) + A(1,2)*B(2) + A(1,5)*B(5)
      X(2,2) = X(2,2) + A(2,4)*B(4) + A(2,2)*B(2) + A(2,5)*B(5)
      X(3,2) = X(3,2) + A(3,4)*B(4) + A(3,2)*B(2) + A(3,5)*B(5)
      X(4,2) = X(4,2) + A(4,4)*B(4) + A(4,2)*B(2) + A(4,5)*B(5)
      X(5,2) = X(5,2) + A(5,4)*B(4) + A(5,2)*B(2) + A(5,5)*B(5)
      X(6,2) = X(6,2) + A(6,4)*B(4) + A(6,2)*B(2) + A(6,5)*B(5)
      X(1,3) = X(1,3) + A(1,6)*B(6) + A(1,5)*B(5) + A(1,3)*B(3)
      X(2,3) = X(2,3) + A(2,6)*B(6) + A(2,5)*B(5) + A(2,3)*B(3)
      X(3,3) = X(3,3) + A(3,6)*B(6) + A(3,5)*B(5) + A(3,3)*B(3)
      X(4,3) = X(4,3) + A(4,6)*B(6) + A(4,5)*B(5) + A(4,3)*B(3)
      X(5,3) = X(5,3) + A(5,6)*B(6) + A(5,5)*B(5) + A(5,3)*B(3)
      X(6,3) = X(6,3) + A(6,6)*B(6) + A(6,5)*B(5) + A(6,3)*B(3)
      X(1,4) = X(1,4) + r1d2 * (A(1,1)*B(4) + A(1,4)*B(2) + A(1,6)*B(5)
     *                        + A(1,4)*B(1) + A(1,2)*B(4) + A(1,5)*B(6))
      X(2,4) = X(2,4) + r1d2 * (A(2,1)*B(4) + A(2,4)*B(2) + A(2,6)*B(5)
     *                        + A(2,4)*B(1) + A(2,2)*B(4) + A(2,5)*B(6))
      X(3,4) = X(3,4) + r1d2 * (A(3,1)*B(4) + A(3,4)*B(2) + A(3,6)*B(5)
     *                        + A(3,4)*B(1) + A(3,2)*B(4) + A(3,5)*B(6))
      X(4,4) = X(4,4) + r1d2 * (A(4,1)*B(4) + A(4,4)*B(2) + A(4,6)*B(5)
     *                        + A(4,4)*B(1) + A(4,2)*B(4) + A(4,5)*B(6))
      X(5,4) = X(5,4) + r1d2 * (A(5,1)*B(4) + A(5,4)*B(2) + A(5,6)*B(5)
     *                        + A(5,4)*B(1) + A(5,2)*B(4) + A(5,5)*B(6))
      X(6,4) = X(6,4) + r1d2 * (A(6,1)*B(4) + A(6,4)*B(2) + A(6,6)*B(5)
     *                        + A(6,4)*B(1) + A(6,2)*B(4) + A(6,5)*B(6))
      X(1,5) = X(1,5) + r1d2 * (A(1,4)*B(6) + A(1,2)*B(5) + A(1,5)*B(3)
     *                        + A(1,6)*B(4) + A(1,5)*B(2) + A(1,3)*B(5))
      X(2,5) = X(2,5) + r1d2 * (A(2,4)*B(6) + A(2,2)*B(5) + A(2,5)*B(3)
     *                        + A(2,6)*B(4) + A(2,5)*B(2) + A(2,3)*B(5))
      X(3,5) = X(3,5) + r1d2 * (A(3,4)*B(6) + A(3,2)*B(5) + A(3,5)*B(3)
     *                        + A(3,6)*B(4) + A(3,5)*B(2) + A(3,3)*B(5))
      X(4,5) = X(4,5) + r1d2 * (A(4,4)*B(6) + A(4,2)*B(5) + A(4,5)*B(3)
     *                        + A(4,6)*B(4) + A(4,5)*B(2) + A(4,3)*B(5))
      X(5,5) = X(5,5) + r1d2 * (A(5,4)*B(6) + A(5,2)*B(5) + A(5,5)*B(3)
     *                        + A(5,6)*B(4) + A(5,5)*B(2) + A(5,3)*B(5))
      X(6,5) = X(6,5) + r1d2 * (A(6,4)*B(6) + A(6,2)*B(5) + A(6,5)*B(3)
     *                        + A(6,6)*B(4) + A(6,5)*B(2) + A(6,3)*B(5))
      X(1,6) = X(1,6) + r1d2 * (A(1,6)*B(1) + A(1,5)*B(4) + A(1,3)*B(6)
     *                        + A(1,1)*B(6) + A(1,4)*B(5) + A(1,6)*B(3))
      X(2,6) = X(2,6) + r1d2 * (A(2,6)*B(1) + A(2,5)*B(4) + A(2,3)*B(6)
     *                        + A(2,1)*B(6) + A(2,4)*B(5) + A(2,6)*B(3))
      X(3,6) = X(3,6) + r1d2 * (A(3,6)*B(1) + A(3,5)*B(4) + A(3,3)*B(6)
     *                        + A(3,1)*B(6) + A(3,4)*B(5) + A(3,6)*B(3))
      X(4,6) = X(4,6) + r1d2 * (A(4,6)*B(1) + A(4,5)*B(4) + A(4,3)*B(6)
     *                        + A(4,1)*B(6) + A(4,4)*B(5) + A(4,6)*B(3))
      X(5,6) = X(5,6) + r1d2 * (A(5,6)*B(1) + A(5,5)*B(4) + A(5,3)*B(6)
     *                        + A(5,1)*B(6) + A(5,4)*B(5) + A(5,6)*B(3))
      X(6,6) = X(6,6) + r1d2 * (A(6,6)*B(1) + A(6,5)*B(4) + A(6,3)*B(6)
     *                        + A(6,1)*B(6) + A(6,4)*B(5) + A(6,6)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_4s_am(X,A,B)
      implicit none
      double precision X(6,*), A(6,*), B(*), r1d2
      r1d2 = 0.5d0
      X(1,1) = X(1,1) - A(1,1)*B(1) - A(1,4)*B(4) - A(1,6)*B(6)
      X(2,1) = X(2,1) - A(2,1)*B(1) - A(2,4)*B(4) - A(2,6)*B(6)
      X(3,1) = X(3,1) - A(3,1)*B(1) - A(3,4)*B(4) - A(3,6)*B(6)
      X(4,1) = X(4,1) - A(4,1)*B(1) - A(4,4)*B(4) - A(4,6)*B(6)
      X(5,1) = X(5,1) - A(5,1)*B(1) - A(5,4)*B(4) - A(5,6)*B(6)
      X(6,1) = X(6,1) - A(6,1)*B(1) - A(6,4)*B(4) - A(6,6)*B(6)
      X(1,2) = X(1,2) - A(1,4)*B(4) - A(1,2)*B(2) - A(1,5)*B(5)
      X(2,2) = X(2,2) - A(2,4)*B(4) - A(2,2)*B(2) - A(2,5)*B(5)
      X(3,2) = X(3,2) - A(3,4)*B(4) - A(3,2)*B(2) - A(3,5)*B(5)
      X(4,2) = X(4,2) - A(4,4)*B(4) - A(4,2)*B(2) - A(4,5)*B(5)
      X(5,2) = X(5,2) - A(5,4)*B(4) - A(5,2)*B(2) - A(5,5)*B(5)
      X(6,2) = X(6,2) - A(6,4)*B(4) - A(6,2)*B(2) - A(6,5)*B(5)
      X(1,3) = X(1,3) - A(1,6)*B(6) - A(1,5)*B(5) - A(1,3)*B(3)
      X(2,3) = X(2,3) - A(2,6)*B(6) - A(2,5)*B(5) - A(2,3)*B(3)
      X(3,3) = X(3,3) - A(3,6)*B(6) - A(3,5)*B(5) - A(3,3)*B(3)
      X(4,3) = X(4,3) - A(4,6)*B(6) - A(4,5)*B(5) - A(4,3)*B(3)
      X(5,3) = X(5,3) - A(5,6)*B(6) - A(5,5)*B(5) - A(5,3)*B(3)
      X(6,3) = X(6,3) - A(6,6)*B(6) - A(6,5)*B(5) - A(6,3)*B(3)
      X(1,4) = X(1,4) - r1d2 * (A(1,1)*B(4) + A(1,4)*B(2) + A(1,6)*B(5)
     *                        + A(1,4)*B(1) + A(1,2)*B(4) + A(1,5)*B(6))
      X(2,4) = X(2,4) - r1d2 * (A(2,1)*B(4) + A(2,4)*B(2) + A(2,6)*B(5)
     *                        + A(2,4)*B(1) + A(2,2)*B(4) + A(2,5)*B(6))
      X(3,4) = X(3,4) - r1d2 * (A(3,1)*B(4) + A(3,4)*B(2) + A(3,6)*B(5)
     *                        + A(3,4)*B(1) + A(3,2)*B(4) + A(3,5)*B(6))
      X(4,4) = X(4,4) - r1d2 * (A(4,1)*B(4) + A(4,4)*B(2) + A(4,6)*B(5)
     *                        + A(4,4)*B(1) + A(4,2)*B(4) + A(4,5)*B(6))
      X(5,4) = X(5,4) - r1d2 * (A(5,1)*B(4) + A(5,4)*B(2) + A(5,6)*B(5)
     *                        + A(5,4)*B(1) + A(5,2)*B(4) + A(5,5)*B(6))
      X(6,4) = X(6,4) - r1d2 * (A(6,1)*B(4) + A(6,4)*B(2) + A(6,6)*B(5)
     *                        + A(6,4)*B(1) + A(6,2)*B(4) + A(6,5)*B(6))
      X(1,5) = X(1,5) - r1d2 * (A(1,4)*B(6) + A(1,2)*B(5) + A(1,5)*B(3)
     *                        + A(1,6)*B(4) + A(1,5)*B(2) + A(1,3)*B(5))
      X(2,5) = X(2,5) - r1d2 * (A(2,4)*B(6) + A(2,2)*B(5) + A(2,5)*B(3)
     *                        + A(2,6)*B(4) + A(2,5)*B(2) + A(2,3)*B(5))
      X(3,5) = X(3,5) - r1d2 * (A(3,4)*B(6) + A(3,2)*B(5) + A(3,5)*B(3)
     *                        + A(3,6)*B(4) + A(3,5)*B(2) + A(3,3)*B(5))
      X(4,5) = X(4,5) - r1d2 * (A(4,4)*B(6) + A(4,2)*B(5) + A(4,5)*B(3)
     *                        + A(4,6)*B(4) + A(4,5)*B(2) + A(4,3)*B(5))
      X(5,5) = X(5,5) - r1d2 * (A(5,4)*B(6) + A(5,2)*B(5) + A(5,5)*B(3)
     *                        + A(5,6)*B(4) + A(5,5)*B(2) + A(5,3)*B(5))
      X(6,5) = X(6,5) - r1d2 * (A(6,4)*B(6) + A(6,2)*B(5) + A(6,5)*B(3)
     *                        + A(6,6)*B(4) + A(6,5)*B(2) + A(6,3)*B(5))
      X(1,6) = X(1,6) - r1d2 * (A(1,6)*B(1) + A(1,5)*B(4) + A(1,3)*B(6)
     *                        + A(1,1)*B(6) + A(1,4)*B(5) + A(1,6)*B(3))
      X(2,6) = X(2,6) - r1d2 * (A(2,6)*B(1) + A(2,5)*B(4) + A(2,3)*B(6)
     *                        + A(2,1)*B(6) + A(2,4)*B(5) + A(2,6)*B(3))
      X(3,6) = X(3,6) - r1d2 * (A(3,6)*B(1) + A(3,5)*B(4) + A(3,3)*B(6)
     *                        + A(3,1)*B(6) + A(3,4)*B(5) + A(3,6)*B(3))
      X(4,6) = X(4,6) - r1d2 * (A(4,6)*B(1) + A(4,5)*B(4) + A(4,3)*B(6)
     *                        + A(4,1)*B(6) + A(4,4)*B(5) + A(4,6)*B(3))
      X(5,6) = X(5,6) - r1d2 * (A(5,6)*B(1) + A(5,5)*B(4) + A(5,3)*B(6)
     *                        + A(5,1)*B(6) + A(5,4)*B(5) + A(5,6)*B(3))
      X(6,6) = X(6,6) - r1d2 * (A(6,6)*B(1) + A(6,5)*B(4) + A(6,3)*B(6)
     *                        + A(6,1)*B(6) + A(6,4)*B(5) + A(6,6)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine mult_4_4s_af(X,A,B,fact)
      implicit none
      double precision X(6,*), A(6,*), B(*), fact, factd2
      factd2 = 0.5d0 * fact
      X(1,1) = X(1,1) + fact * (A(1,1)*B(1) + A(1,4)*B(4) + A(1,6)*B(6))
      X(2,1) = X(2,1) + fact * (A(2,1)*B(1) + A(2,4)*B(4) + A(2,6)*B(6))
      X(3,1) = X(3,1) + fact * (A(3,1)*B(1) + A(3,4)*B(4) + A(3,6)*B(6))
      X(4,1) = X(4,1) + fact * (A(4,1)*B(1) + A(4,4)*B(4) + A(4,6)*B(6))
      X(5,1) = X(5,1) + fact * (A(5,1)*B(1) + A(5,4)*B(4) + A(5,6)*B(6))
      X(6,1) = X(6,1) + fact * (A(6,1)*B(1) + A(6,4)*B(4) + A(6,6)*B(6))
      X(1,2) = X(1,2) + fact * (A(1,4)*B(4) + A(1,2)*B(2) + A(1,5)*B(5))
      X(2,2) = X(2,2) + fact * (A(2,4)*B(4) + A(2,2)*B(2) + A(2,5)*B(5))
      X(3,2) = X(3,2) + fact * (A(3,4)*B(4) + A(3,2)*B(2) + A(3,5)*B(5))
      X(4,2) = X(4,2) + fact * (A(4,4)*B(4) + A(4,2)*B(2) + A(4,5)*B(5))
      X(5,2) = X(5,2) + fact * (A(5,4)*B(4) + A(5,2)*B(2) + A(5,5)*B(5))
      X(6,2) = X(6,2) + fact * (A(6,4)*B(4) + A(6,2)*B(2) + A(6,5)*B(5))
      X(1,3) = X(1,3) + fact * (A(1,6)*B(6) + A(1,5)*B(5) + A(1,3)*B(3))
      X(2,3) = X(2,3) + fact * (A(2,6)*B(6) + A(2,5)*B(5) + A(2,3)*B(3))
      X(3,3) = X(3,3) + fact * (A(3,6)*B(6) + A(3,5)*B(5) + A(3,3)*B(3))
      X(4,3) = X(4,3) + fact * (A(4,6)*B(6) + A(4,5)*B(5) + A(4,3)*B(3))
      X(5,3) = X(5,3) + fact * (A(5,6)*B(6) + A(5,5)*B(5) + A(5,3)*B(3))
      X(6,3) = X(6,3) + fact * (A(6,6)*B(6) + A(6,5)*B(5) + A(6,3)*B(3))
      X(1,4) = X(1,4) + factd2*(A(1,1)*B(4) + A(1,4)*B(2) + A(1,6)*B(5)
     *                        + A(1,4)*B(1) + A(1,2)*B(4) + A(1,5)*B(6))
      X(2,4) = X(2,4) + factd2*(A(2,1)*B(4) + A(2,4)*B(2) + A(2,6)*B(5)
     *                        + A(2,4)*B(1) + A(2,2)*B(4) + A(2,5)*B(6))
      X(3,4) = X(3,4) + factd2*(A(3,1)*B(4) + A(3,4)*B(2) + A(3,6)*B(5)
     *                        + A(3,4)*B(1) + A(3,2)*B(4) + A(3,5)*B(6))
      X(4,4) = X(4,4) + factd2*(A(4,1)*B(4) + A(4,4)*B(2) + A(4,6)*B(5)
     *                        + A(4,4)*B(1) + A(4,2)*B(4) + A(4,5)*B(6))
      X(5,4) = X(5,4) + factd2*(A(5,1)*B(4) + A(5,4)*B(2) + A(5,6)*B(5)
     *                        + A(5,4)*B(1) + A(5,2)*B(4) + A(5,5)*B(6))
      X(6,4) = X(6,4) + factd2*(A(6,1)*B(4) + A(6,4)*B(2) + A(6,6)*B(5)
     *                        + A(6,4)*B(1) + A(6,2)*B(4) + A(6,5)*B(6))
      X(1,5) = X(1,5) + factd2*(A(1,4)*B(6) + A(1,2)*B(5) + A(1,5)*B(3)
     *                        + A(1,6)*B(4) + A(1,5)*B(2) + A(1,3)*B(5))
      X(2,5) = X(2,5) + factd2*(A(2,4)*B(6) + A(2,2)*B(5) + A(2,5)*B(3)
     *                        + A(2,6)*B(4) + A(2,5)*B(2) + A(2,3)*B(5))
      X(3,5) = X(3,5) + factd2*(A(3,4)*B(6) + A(3,2)*B(5) + A(3,5)*B(3)
     *                        + A(3,6)*B(4) + A(3,5)*B(2) + A(3,3)*B(5))
      X(4,5) = X(4,5) + factd2*(A(4,4)*B(6) + A(4,2)*B(5) + A(4,5)*B(3)
     *                        + A(4,6)*B(4) + A(4,5)*B(2) + A(4,3)*B(5))
      X(5,5) = X(5,5) + factd2*(A(5,4)*B(6) + A(5,2)*B(5) + A(5,5)*B(3)
     *                        + A(5,6)*B(4) + A(5,5)*B(2) + A(5,3)*B(5))
      X(6,5) = X(6,5) + factd2*(A(6,4)*B(6) + A(6,2)*B(5) + A(6,5)*B(3)
     *                        + A(6,6)*B(4) + A(6,5)*B(2) + A(6,3)*B(5))
      X(1,6) = X(1,6) + factd2*(A(1,6)*B(1) + A(1,5)*B(4) + A(1,3)*B(6)
     *                        + A(1,1)*B(6) + A(1,4)*B(5) + A(1,6)*B(3))
      X(2,6) = X(2,6) + factd2*(A(2,6)*B(1) + A(2,5)*B(4) + A(2,3)*B(6)
     *                        + A(2,1)*B(6) + A(2,4)*B(5) + A(2,6)*B(3))
      X(3,6) = X(3,6) + factd2*(A(3,6)*B(1) + A(3,5)*B(4) + A(3,3)*B(6)
     *                        + A(3,1)*B(6) + A(3,4)*B(5) + A(3,6)*B(3))
      X(4,6) = X(4,6) + factd2*(A(4,6)*B(1) + A(4,5)*B(4) + A(4,3)*B(6)
     *                        + A(4,1)*B(6) + A(4,4)*B(5) + A(4,6)*B(3))
      X(5,6) = X(5,6) + factd2*(A(5,6)*B(1) + A(5,5)*B(4) + A(5,3)*B(6)
     *                        + A(5,1)*B(6) + A(5,4)*B(5) + A(5,6)*B(3))
      X(6,6) = X(6,6) + factd2*(A(6,6)*B(1) + A(6,5)*B(4) + A(6,3)*B(6)
     *                        + A(6,1)*B(6) + A(6,4)*B(5) + A(6,6)*B(3))
      return
      end
c
c=============================================================================
c  Differentiation
c                       D AT      mit  A,T,AT  sym. 2.Stufe
c                 X =  ------       
c                       D T      -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_AT_p(X,A)
      implicit none
      double precision A(6), X(6,6), r0, r1d4, r1d2
      r0   = 0.d0
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = A(1)
      X(2,2) = A(2)
      X(3,3) = A(3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,4) = r1d2 * A(4)
      X(4,1) = X(1,4)
      X(1,5) = r0
      X(5,1) = r0
      X(1,6) = r1d2 * A(6)
      X(6,1) = X(1,6)
      X(2,4) = r1d2 * A(4)
      X(4,2) = X(2,4)
      X(2,5) = r1d2 * A(5)
      X(5,2) = X(2,5)
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(3,5) = r1d2 * A(5)
      X(5,3) = X(3,5)
      X(3,6) = r1d2 * A(6)
      X(6,3) = X(3,6)

      X(4,4) = r1d4 * (A(1) + A(2))
      X(5,5) = r1d4 * (A(2) + A(3))
      X(6,6) = r1d4 * (A(3) + A(1))

      X(4,5) = r1d4 * A(6)
      X(5,4) = X(4,5)
      X(4,6) = r1d4 * A(5)
      X(6,4) = X(4,6)
      X(5,6) = r1d4 * A(4)
      X(6,5) = X(5,6)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_AT_m(X,A)
      implicit none
      double precision A(6), X(6,6), r0, r1d4, r1d2
      r0   = 0.d0
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = - A(1)
      X(2,2) = - A(2)
      X(3,3) = - A(3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,4) = - r1d2 * A(4)
      X(4,1) = X(1,4)
      X(1,5) = r0
      X(5,1) = r0
      X(1,6) = - r1d2 * A(6)
      X(6,1) = X(1,6)
      X(2,4) = - r1d2 * A(4)
      X(4,2) = X(2,4)
      X(2,5) = - r1d2 * A(5)
      X(5,2) = X(2,5)
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(3,5) = - r1d2 * A(5)
      X(5,3) = X(3,5)
      X(3,6) = - r1d2 * A(6)
      X(6,3) = X(3,6)

      X(4,4) = - r1d4 * (A(1) + A(2))
      X(5,5) = - r1d4 * (A(2) + A(3))
      X(6,6) = - r1d4 * (A(3) + A(1))

      X(4,5) = - r1d4 * A(6)
      X(5,4) = X(4,5)
      X(4,6) = - r1d4 * A(5)
      X(6,4) = X(4,6)
      X(5,6) = - r1d4 * A(4)
      X(6,5) = X(5,6)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_AT_f(X,A,fact)
      implicit none
      double precision A(6), X(6,6), r0, fact, factd2, factd4
      r0   = 0.d0
      factd4 = fact * 0.25d0
      factd2 = fact * 0.5d0

      X(1,1) = fact * A(1)
      X(2,2) = fact * A(2)
      X(3,3) = fact * A(3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,4) = factd2 * A(4)
      X(4,1) = X(1,4)
      X(1,5) = r0
      X(5,1) = r0
      X(1,6) = factd2 * A(6)
      X(6,1) = X(1,6)
      X(2,4) = factd2 * A(4)
      X(4,2) = X(2,4)
      X(2,5) = factd2 * A(5)
      X(5,2) = X(2,5)
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(3,5) = factd2 * A(5)
      X(5,3) = X(3,5)
      X(3,6) = factd2 * A(6)
      X(6,3) = X(3,6)

      X(4,4) = factd4 * (A(1) + A(2))
      X(5,5) = factd4 * (A(2) + A(3))
      X(6,6) = factd4 * (A(3) + A(1))

      X(4,5) = factd4 * A(6)
      X(5,4) = X(4,5)
      X(4,6) = factd4 * A(5)
      X(6,4) = X(4,6)
      X(5,6) = factd4 * A(4)
      X(6,5) = X(5,6)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_AT_ap(X,A)
      implicit none
      double precision A(6), X(6,6), r1d4, r1d2, z
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = X(1,1) + A(1)
      X(2,2) = X(2,2) + A(2)
      X(3,3) = X(3,3) + A(3)

      z = r1d2 * A(4)
      X(1,4) = X(1,4) + z
      X(4,1) = X(4,1) + z
      z = r1d2 * A(6)
      X(1,6) = X(1,6) + z
      X(6,1) = X(6,1) + z
      z = r1d2 * A(4)
      X(2,4) = X(2,4) + z
      X(4,2) = X(4,2) + z
      z = r1d2 * A(5)
      X(2,5) = X(2,5) + z
      X(5,2) = X(5,2) + z
      z = r1d2 * A(5)
      X(3,5) = X(3,5) + z
      X(5,3) = X(5,3) + z
      z = r1d2 * A(6)
      X(3,6) = X(3,6) + z
      X(6,3) = X(6,3) + z

      X(4,4) = X(4,4) + r1d4 * (A(1) + A(2))
      X(5,5) = X(5,5) + r1d4 * (A(2) + A(3))
      X(6,6) = X(6,6) + r1d4 * (A(3) + A(1))

      z = r1d4 * A(6)
      X(4,5) = X(4,5) + z
      X(5,4) = X(5,4) + z
      z = r1d4 * A(5)
      X(4,6) = X(4,6) + z
      X(6,4) = X(6,4) + z
      z = r1d4 * A(4)
      X(5,6) = X(5,6) + z
      X(6,5) = X(6,5) + z

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_AT_am(X,A)
      implicit none
      double precision A(6), X(6,6), r1d4, r1d2, z
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = X(1,1) - A(1)
      X(2,2) = X(2,2) - A(2)
      X(3,3) = X(3,3) - A(3)

      z = r1d2 * A(4)
      X(1,4) = X(1,4) - z
      X(4,1) = X(4,1) - z
      z = r1d2 * A(6)
      X(1,6) = X(1,6) - z
      X(6,1) = X(6,1) - z
      z = r1d2 * A(4)
      X(2,4) = X(2,4) - z
      X(4,2) = X(4,2) - z
      z = r1d2 * A(5)
      X(2,5) = X(2,5) - z
      X(5,2) = X(5,2) - z
      z = r1d2 * A(5)
      X(3,5) = X(3,5) - z
      X(5,3) = X(5,3) - z
      z = r1d2 * A(6)
      X(3,6) = X(3,6) - z
      X(6,3) = X(6,3) - z

      X(4,4) = X(4,4) - r1d4 * (A(1) + A(2))
      X(5,5) = X(5,5) - r1d4 * (A(2) + A(3))
      X(6,6) = X(6,6) - r1d4 * (A(3) + A(1))

      z = r1d4 * A(6)
      X(4,5) = X(4,5) - z
      X(5,4) = X(5,4) - z
      z = r1d4 * A(5)
      X(4,6) = X(4,6) - z
      X(6,4) = X(6,4) - z
      z = r1d4 * A(4)
      X(5,6) = X(5,6) - z
      X(6,5) = X(6,5) - z

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_AT_af(X,A,fact)
      implicit none
      double precision A(6), X(6,6), fact, factd4, factd2, z
      factd4 = fact * 0.25d0
      factd2 = fact * 0.5d0

      X(1,1) = X(1,1) + fact * A(1)
      X(2,2) = X(2,2) + fact * A(2)
      X(3,3) = X(3,3) + fact * A(3)

      z = factd2 * A(4)
      X(1,4) = X(1,4) + z
      X(4,1) = X(4,1) + z
      z = factd2 * A(6)
      X(1,6) = X(1,6) + z
      X(6,1) = X(6,1) + z
      z = factd2 * A(4)
      X(2,4) = X(2,4) + z
      X(4,2) = X(4,2) + z
      z = factd2 * A(5)
      X(2,5) = X(2,5) + z
      X(5,2) = X(5,2) + z
      z = factd2 * A(5)
      X(3,5) = X(3,5) + z
      X(5,3) = X(5,3) + z
      z = factd2 * A(6)
      X(3,6) = X(3,6) + z
      X(6,3) = X(6,3) + z

      X(4,4) = X(4,4) + factd4 * (A(1) + A(2))
      X(5,5) = X(5,5) + factd4 * (A(2) + A(3))
      X(6,6) = X(6,6) + factd4 * (A(3) + A(1))

      z = factd4 * A(6)
      X(4,5) = X(4,5) + z
      X(5,4) = X(5,4) + z
      z = factd4 * A(5)
      X(4,6) = X(4,6) + z
      X(6,4) = X(6,4) + z
      z = factd4 * A(4)
      X(5,6) = X(5,6) + z
      X(6,5) = X(6,5) + z

      return
      end
c
c=============================================================================
c  Differentiation
c                       D AT      mit  AT,T  sym., A unsym. 2.Stufe
c                 X =  ------       
c                       D T       -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_UT_p(X,A)
      implicit none
      double precision A(3,3), X(6,6), r0, r1d4, r1d2
      r0   = 0.d0
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = A(1,1)
      X(2,2) = A(2,2)
      X(3,3) = A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,5) = r0
      X(5,1) = r0
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(4,1) = r1d2 * A(2,1)
      X(2,4) = X(4,1)
      X(1,4) = r1d2 * A(1,2)
      X(4,2) = X(1,4)
      X(5,2) = r1d2 * A(3,2)
      X(3,5) = X(5,2)
      X(2,5) = r1d2 * A(2,3)
      X(5,3) = X(2,5)
      X(6,1) = r1d2 * A(3,1)
      X(3,6) = X(6,1)
      X(1,6) = r1d2 * A(1,3)
      X(6,3) = X(1,6)

      X(4,4) = r1d4 * (A(1,1) + A(2,2))
      X(5,5) = r1d4 * (A(2,2) + A(3,3))
      X(6,6) = r1d4 * (A(3,3) + A(1,1))

      X(4,5) = r1d4 * A(1,3)
      X(5,4) = r1d4 * A(3,1)
      X(4,6) = r1d4 * A(2,3)
      X(6,4) = r1d4 * A(3,2)
      X(5,6) = r1d4 * A(2,1)
      X(6,5) = r1d4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UT_m(X,A)
      implicit none
      double precision A(3,3), X(6,6), r0, r1d4, r1d2
      r0   = 0.d0
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = - A(1,1)
      X(2,2) = - A(2,2)
      X(3,3) = - A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,5) = r0
      X(5,1) = r0
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(4,1) = - r1d2 * A(2,1)
      X(2,4) = X(4,1)
      X(1,4) = - r1d2 * A(1,2)
      X(4,2) = X(1,4)
      X(5,2) = - r1d2 * A(3,2)
      X(3,5) = X(5,2)
      X(2,5) = - r1d2 * A(2,3)
      X(5,3) = X(2,5)
      X(6,1) = - r1d2 * A(3,1)
      X(3,6) = X(6,1)
      X(1,6) = - r1d2 * A(1,3)
      X(6,3) = X(1,6)

      X(4,4) = - r1d4 * (A(1,1) + A(2,2))
      X(5,5) = - r1d4 * (A(2,2) + A(3,3))
      X(6,6) = - r1d4 * (A(3,3) + A(1,1))

      X(4,5) = - r1d4 * A(1,3)
      X(5,4) = - r1d4 * A(3,1)
      X(4,6) = - r1d4 * A(2,3)
      X(6,4) = - r1d4 * A(3,2)
      X(5,6) = - r1d4 * A(2,1)
      X(6,5) = - r1d4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UT_f(X,A,fact)
      implicit none
      double precision A(3,3), X(6,6), r0, fact,  factd4, factd2
      r0     = 0.d0
      factd4 = fact * 0.25d0
      factd2 = fact * 0.5d0

      X(1,1) = fact * A(1,1)
      X(2,2) = fact * A(2,2)
      X(3,3) = fact * A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,5) = r0
      X(5,1) = r0
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(4,1) = factd2 * A(2,1)
      X(2,4) = X(4,1)
      X(1,4) = factd2 * A(1,2)
      X(4,2) = X(1,4)
      X(5,2) = factd2 * A(3,2)
      X(3,5) = X(5,2)
      X(2,5) = factd2 * A(2,3)
      X(5,3) = X(2,5)
      X(6,1) = factd2 * A(3,1)
      X(3,6) = X(6,1)
      X(1,6) = factd2 * A(1,3)
      X(6,3) = X(1,6)

      X(4,4) = factd4 * (A(1,1) + A(2,2))
      X(5,5) = factd4 * (A(2,2) + A(3,3))
      X(6,6) = factd4 * (A(3,3) + A(1,1))

      X(4,5) = factd4 * A(1,3)
      X(5,4) = factd4 * A(3,1)
      X(4,6) = factd4 * A(2,3)
      X(6,4) = factd4 * A(3,2)
      X(5,6) = factd4 * A(2,1)
      X(6,5) = factd4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UT_ap(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d4, r1d2, z
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = X(1,1) + A(1,1)
      X(2,2) = X(2,2) + A(2,2)
      X(3,3) = X(3,3) + A(3,3)

      z = r1d2 * A(2,1)
      X(4,1) = X(4,1) + z
      X(2,4) = X(2,4) + z
      z = r1d2 * A(1,2)
      X(1,4) = X(1,4) + z
      X(4,2) = X(4,2) + z
      z = r1d2 * A(3,2)
      X(5,2) = X(5,2) + z
      X(3,5) = X(3,5) + z
      z = r1d2 * A(2,3)
      X(2,5) = X(2,5) + z
      X(5,3) = X(5,3) + z
      z = r1d2 * A(3,1)
      X(6,1) = X(6,1) + z
      X(3,6) = X(3,6) + z
      z = r1d2 * A(1,3)
      X(1,6) = X(1,6) + z
      X(6,3) = X(6,3) + z

      X(4,4) = X(4,4) + r1d4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) + r1d4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) + r1d4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) + r1d4 * A(1,3)
      X(5,4) = X(5,4) + r1d4 * A(3,1)
      X(4,6) = X(4,6) + r1d4 * A(2,3)
      X(6,4) = X(6,4) + r1d4 * A(3,2)
      X(5,6) = X(5,6) + r1d4 * A(2,1)
      X(6,5) = X(6,5) + r1d4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UT_am(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d4, r1d2, z
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = X(1,1) - A(1,1)
      X(2,2) = X(2,2) - A(2,2)
      X(3,3) = X(3,3) - A(3,3)

      z = r1d2 * A(2,1)
      X(4,1) = X(4,1) - z
      X(2,4) = X(2,4) - z
      z = r1d2 * A(1,2)
      X(1,4) = X(1,4) - z
      X(4,2) = X(4,2) - z
      z = r1d2 * A(3,2)
      X(5,2) = X(5,2) - z
      X(3,5) = X(3,5) - z
      z = r1d2 * A(2,3)
      X(2,5) = X(2,5) - z
      X(5,3) = X(5,3) - z
      z = r1d2 * A(3,1)
      X(6,1) = X(6,1) - z
      X(3,6) = X(3,6) - z
      z = r1d2 * A(1,3)
      X(1,6) = X(1,6) - z
      X(6,3) = X(6,3) - z

      X(4,4) = X(4,4) - r1d4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) - r1d4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) - r1d4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) - r1d4 * A(1,3)
      X(5,4) = X(5,4) - r1d4 * A(3,1)
      X(4,6) = X(4,6) - r1d4 * A(2,3)
      X(6,4) = X(6,4) - r1d4 * A(3,2)
      X(5,6) = X(5,6) - r1d4 * A(2,1)
      X(6,5) = X(6,5) - r1d4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UT_af(X,A,fact)
      implicit none
      double precision A(3,3), X(6,6), fact, factd4, factd2, z
      factd4 = fact * 0.25d0
      factd2 = fact * 0.5d0

      X(1,1) = X(1,1) + fact * A(1,1)
      X(2,2) = X(2,2) + fact * A(2,2)
      X(3,3) = X(3,3) + fact * A(3,3)

      z = factd2 * A(2,1)
      X(4,1) = X(4,1) + z
      X(2,4) = X(2,4) + z
      z = factd2 * A(1,2)
      X(1,4) = X(1,4) + z
      X(4,2) = X(4,2) + z
      z = factd2 * A(3,2)
      X(5,2) = X(5,2) + z
      X(3,5) = X(3,5) + z
      z = factd2 * A(2,3)
      X(2,5) = X(2,5) + z
      X(5,3) = X(5,3) + z
      z = factd2 * A(3,1)
      X(6,1) = X(6,1) + z
      X(3,6) = X(3,6) + z
      z = factd2 * A(1,3)
      X(1,6) = X(1,6) + z
      X(6,3) = X(6,3) + z

      X(4,4) = X(4,4) + factd4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) + factd4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) + factd4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) + factd4 * A(1,3)
      X(5,4) = X(5,4) + factd4 * A(3,1)
      X(4,6) = X(4,6) + factd4 * A(2,3)
      X(6,4) = X(6,4) + factd4 * A(3,2)
      X(5,6) = X(5,6) + factd4 * A(2,1)
      X(6,5) = X(6,5) + factd4 * A(1,2)

      return
      end
c
c=============================================================================
c  Differentiation
c                        D TA      mit  AT,T  sym., A unsym. 2.Stufe
c                  X =  ------       
c                        D T       -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_TU_p(X,A)
      implicit none
      double precision A(3,3), X(6,6), r0, r1d4, r1d2
      r0   = 0.d0
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = A(1,1)
      X(2,2) = A(2,2)
      X(3,3) = A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,5) = r0
      X(5,1) = r0
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(4,1) = r1d2 * A(1,2)
      X(2,4) = X(4,1)
      X(1,4) = r1d2 * A(2,1)
      X(4,2) = X(1,4)
      X(5,2) = r1d2 * A(2,3)
      X(3,5) = X(5,2)
      X(2,5) = r1d2 * A(3,2)
      X(5,3) = X(2,5)
      X(6,1) = r1d2 * A(1,3)
      X(3,6) = X(6,1)
      X(1,6) = r1d2 * A(3,1)
      X(6,3) = X(1,6)

      X(4,4) = r1d4 * (A(1,1) + A(2,2))
      X(5,5) = r1d4 * (A(2,2) + A(3,3))
      X(6,6) = r1d4 * (A(3,3) + A(1,1))

      X(4,5) = r1d4 * A(3,1)
      X(5,4) = r1d4 * A(1,3)
      X(4,6) = r1d4 * A(3,2)
      X(6,4) = r1d4 * A(2,3)
      X(5,6) = r1d4 * A(1,2)
      X(6,5) = r1d4 * A(2,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TU_m(X,A)
      implicit none
      double precision A(3,3), X(6,6), r0, r1d4, r1d2
      r0   = 0.d0
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = - A(1,1)
      X(2,2) = - A(2,2)
      X(3,3) = - A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,5) = r0
      X(5,1) = r0
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(4,1) = - r1d2 * A(1,2)
      X(2,4) = X(4,1)
      X(1,4) = - r1d2 * A(2,1)
      X(4,2) = X(1,4)
      X(5,2) = - r1d2 * A(2,3)
      X(3,5) = X(5,2)
      X(2,5) = - r1d2 * A(3,2)
      X(5,3) = X(2,5)
      X(6,1) = - r1d2 * A(1,3)
      X(3,6) = X(6,1)
      X(1,6) = - r1d2 * A(3,1)
      X(6,3) = X(1,6)

      X(4,4) = - r1d4 * (A(1,1) + A(2,2))
      X(5,5) = - r1d4 * (A(2,2) + A(3,3))
      X(6,6) = - r1d4 * (A(3,3) + A(1,1))

      X(4,5) = - r1d4 * A(3,1)
      X(5,4) = - r1d4 * A(1,3)
      X(4,6) = - r1d4 * A(3,2)
      X(6,4) = - r1d4 * A(2,3)
      X(5,6) = - r1d4 * A(1,2)
      X(6,5) = - r1d4 * A(2,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TU_f(X,A,fact)
      implicit none
      double precision A(3,3), X(6,6), r0, fact,  factd4, factd2
      r0     = 0.d0
      factd4 = fact * 0.25d0
      factd2 = fact * 0.5d0

      X(1,1) = fact * A(1,1)
      X(2,2) = fact * A(2,2)
      X(3,3) = fact * A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0 
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,5) = r0
      X(5,1) = r0
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(4,1) = factd2 * A(1,2)
      X(2,4) = X(4,1)
      X(1,4) = factd2 * A(2,1)
      X(4,2) = X(1,4)
      X(5,2) = factd2 * A(2,3)
      X(3,5) = X(5,2)
      X(2,5) = factd2 * A(3,2)
      X(5,3) = X(2,5)
      X(6,1) = factd2 * A(1,3)
      X(3,6) = X(6,1)
      X(1,6) = factd2 * A(3,1)
      X(6,3) = X(1,6)

      X(4,4) = factd4 * (A(1,1) + A(2,2))
      X(5,5) = factd4 * (A(2,2) + A(3,3))
      X(6,6) = factd4 * (A(3,3) + A(1,1))

      X(4,5) = factd4 * A(3,1)
      X(5,4) = factd4 * A(1,3)
      X(4,6) = factd4 * A(3,2)
      X(6,4) = factd4 * A(2,3)
      X(5,6) = factd4 * A(1,2)
      X(6,5) = factd4 * A(2,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TU_ap(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d4, r1d2, z
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = X(1,1) + A(1,1)
      X(2,2) = X(2,2) + A(2,2)
      X(3,3) = X(3,3) + A(3,3)

      z = r1d2 * A(1,2)
      X(4,1) = X(4,1) + z
      X(2,4) = X(2,4) + z
      z = r1d2 * A(2,1)
      X(1,4) = X(1,4) + z
      X(4,2) = X(4,2) + z
      z = r1d2 * A(2,3)
      X(5,2) = X(5,2) + z
      X(3,5) = X(3,5) + z
      z = r1d2 * A(3,2)
      X(2,5) = X(2,5) + z
      X(5,3) = X(5,3) + z
      z = r1d2 * A(1,3)
      X(6,1) = X(6,1) + z
      X(3,6) = X(3,6) + z
      z = r1d2 * A(3,1)
      X(1,6) = X(1,6) + z
      X(6,3) = X(6,3) + z

      X(4,4) = X(4,4) + r1d4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) + r1d4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) + r1d4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) + r1d4 * A(3,1)
      X(5,4) = X(5,4) + r1d4 * A(1,3)
      X(4,6) = X(4,6) + r1d4 * A(3,2)
      X(6,4) = X(6,4) + r1d4 * A(2,3)
      X(5,6) = X(5,6) + r1d4 * A(1,2)
      X(6,5) = X(6,5) + r1d4 * A(2,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TU_am(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d4, r1d2, z
      r1d4 = 0.25d0
      r1d2 = 0.5d0

      X(1,1) = X(1,1) - A(1,1)
      X(2,2) = X(2,2) - A(2,2)
      X(3,3) = X(3,3) - A(3,3)

      z = r1d2 * A(1,2)
      X(4,1) = X(4,1) - z
      X(2,4) = X(2,4) - z
      z = r1d2 * A(2,1)
      X(1,4) = X(1,4) - z
      X(4,2) = X(4,2) - z
      z = r1d2 * A(2,3)
      X(5,2) = X(5,2) - z
      X(3,5) = X(3,5) - z
      z = r1d2 * A(3,2)
      X(2,5) = X(2,5) - z
      X(5,3) = X(5,3) - z
      z = r1d2 * A(1,3)
      X(6,1) = X(6,1) - z
      X(3,6) = X(3,6) - z
      z = r1d2 * A(3,1)
      X(1,6) = X(1,6) - z
      X(6,3) = X(6,3) - z

      X(4,4) = X(4,4) - r1d4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) - r1d4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) - r1d4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) - r1d4 * A(3,1)
      X(5,4) = X(5,4) - r1d4 * A(1,3)
      X(4,6) = X(4,6) - r1d4 * A(3,2)
      X(6,4) = X(6,4) - r1d4 * A(2,3)
      X(5,6) = X(5,6) - r1d4 * A(1,2)
      X(6,5) = X(6,5) - r1d4 * A(2,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TU_af(X,A,fact)
      implicit none
      double precision A(3,3), X(6,6), fact, factd4, factd2, z
      factd4 = fact * 0.25d0
      factd2 = fact * 0.5d0

      X(1,1) = X(1,1) + fact * A(1,1)
      X(2,2) = X(2,2) + fact * A(2,2)
      X(3,3) = X(3,3) + fact * A(3,3)

      z = factd2 * A(1,2)
      X(4,1) = X(4,1) + z
      X(2,4) = X(2,4) + z
      z = factd2 * A(2,1)
      X(1,4) = X(1,4) + z
      X(4,2) = X(4,2) + z
      z = factd2 * A(2,3)
      X(5,2) = X(5,2) + z
      X(3,5) = X(3,5) + z
      z = factd2 * A(3,2)
      X(2,5) = X(2,5) + z
      X(5,3) = X(5,3) + z
      z = factd2 * A(1,3)
      X(6,1) = X(6,1) + z
      X(3,6) = X(3,6) + z
      z = factd2 * A(3,1)
      X(1,6) = X(1,6) + z
      X(6,3) = X(6,3) + z

      X(4,4) = X(4,4) + factd4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) + factd4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) + factd4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) + factd4 * A(3,1)
      X(5,4) = X(5,4) + factd4 * A(1,3)
      X(4,6) = X(4,6) + factd4 * A(3,2)
      X(6,4) = X(6,4) + factd4 * A(2,3)
      X(5,6) = X(5,6) + factd4 * A(1,2)
      X(6,5) = X(6,5) + factd4 * A(2,1)

      return
      end
c
c=============================================================================
c  Differentiation
c                      D TAT      mit  A,T,TAT  sym.,  AT unsym.  2.Stufe
c                X =  -------          
c                      D T       -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_TAT_p(X,AT)
      implicit none
      double precision X(6,6), AT(3,3), r0, r2, r1d2
      data r0, r1d2, r2 / 0.d0, 0.5d0, 2.d0 /
      X(1,1) = r2 * AT(1,1)
      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0
      X(2,2) = r2 * AT(2,2)
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0
      X(3,3) = r2 * AT(3,3)
      X(4,1) =      AT(1,2)
      X(4,2) =      AT(2,1)
      X(4,3) = r0
      X(5,1) = r0
      X(5,2) =      AT(2,3)
      X(5,3) =      AT(3,2)
      X(6,1) =      AT(1,3)
      X(6,2) = r0
      X(6,3) =      AT(3,1)
      X(1,4) =      AT(2,1)
      X(2,4) =      AT(1,2)
      X(3,4) = r0
      X(1,5) = r0
      X(2,5) =      AT(3,2)
      X(3,5) =      AT(2,3)
      X(1,6) =      AT(3,1)
      X(2,6) = r0
      X(3,6) =      AT(1,3)
      X(4,4) = r1d2 * (AT(1,1)+AT(2,2))
      X(4,5) = r1d2 *  AT(3,1)
      X(4,6) = r1d2 *  AT(3,2)
      X(5,4) = r1d2 *  AT(1,3)
      X(5,5) = r1d2 * (AT(2,2)+AT(3,3))
      X(5,6) = r1d2 *  AT(1,2)
      X(6,4) = r1d2 *  AT(2,3)
      X(6,5) = r1d2 *  AT(2,1)
      X(6,6) = r1d2 * (AT(3,3)+AT(1,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAT_m(X,AT)
      implicit none
      double precision X(6,6), AT(3,3), r0, r2, r1d2
      data r0, r1d2, r2 / 0.d0, 0.5d0, 2.d0 /
      X(1,1) = - r2 * AT(1,1)
      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0
      X(2,2) = - r2 * AT(2,2)
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0
      X(3,3) = - r2 * AT(3,3)
      X(4,1) = -      AT(1,2)
      X(4,2) = -      AT(2,1)
      X(4,3) = r0
      X(5,1) = r0
      X(5,2) = -      AT(2,3)
      X(5,3) = -      AT(3,2)
      X(6,1) = -      AT(1,3)
      X(6,2) = r0
      X(6,3) = -      AT(3,1)
      X(1,4) = -      AT(2,1)
      X(2,4) = -      AT(1,2)
      X(3,4) = r0
      X(1,5) = r0
      X(2,5) = -      AT(3,2)
      X(3,5) = -      AT(2,3)
      X(1,6) = -      AT(3,1)
      X(2,6) = r0
      X(3,6) = -      AT(1,3)
      X(4,4) = - r1d2 * (AT(1,1)+AT(2,2))
      X(4,5) = - r1d2 *  AT(3,1)
      X(4,6) = - r1d2 *  AT(3,2)
      X(5,4) = - r1d2 *  AT(1,3)
      X(5,5) = - r1d2 * (AT(2,2)+AT(3,3))
      X(5,6) = - r1d2 *  AT(1,2)
      X(6,4) = - r1d2 *  AT(2,3)
      X(6,5) = - r1d2 *  AT(2,1)
      X(6,6) = - r1d2 * (AT(3,3)+AT(1,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAT_f(X,AT,fact)
      implicit none
      double precision X(6,6), AT(3,3), r0, fact, factm2, factd2
      r0 = 0.d0
      factm2 = fact * 2.0d0
      factd2 = fact * 0.5d0
      X(1,1) = factm2 *  AT(1,1)
      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0
      X(2,2) = factm2 *  AT(2,2)
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0
      X(3,3) = factm2 *  AT(3,3)
      X(4,1) =           AT(1,2)
      X(4,2) =           AT(2,1)
      X(4,3) = r0
      X(5,1) = r0
      X(5,2) =           AT(2,3)
      X(5,3) =           AT(3,2)
      X(6,1) =           AT(1,3)
      X(6,2) = r0
      X(6,3) =           AT(3,1)
      X(1,4) =           AT(2,1)
      X(2,4) =           AT(1,2)
      X(3,4) = r0
      X(1,5) = r0
      X(2,5) =           AT(3,2)
      X(3,5) =           AT(2,3)
      X(1,6) =           AT(3,1)
      X(2,6) = r0
      X(3,6) =           AT(1,3)
      X(4,4) = factd2 * (AT(1,1)+AT(2,2))
      X(4,5) = factd2 *  AT(3,1)
      X(4,6) = factd2 *  AT(3,2)
      X(5,4) = factd2 *  AT(1,3)
      X(5,5) = factd2 * (AT(2,2)+AT(3,3))
      X(5,6) = factd2 *  AT(1,2)
      X(6,4) = factd2 *  AT(2,3)
      X(6,5) = factd2 *  AT(2,1)
      X(6,6) = factd2 * (AT(3,3)+AT(1,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAT_ap(X,AT)
      implicit none
      double precision X(6,6), AT(3,3), r2, r1d2
      data r1d2, r2 / 0.5d0, 2.d0 /
      X(1,1) = X(1,1) + r2 * AT(1,1)
      X(2,2) = X(2,2) + r2 * AT(2,2)
      X(3,3) = X(3,3) + r2 * AT(3,3)
      X(4,1) = X(4,1) +      AT(1,2)
      X(4,2) = X(4,2) +      AT(2,1)
      X(5,2) = X(5,2) +      AT(2,3)
      X(5,3) = X(5,3) +      AT(3,2)
      X(6,3) = X(6,3) +      AT(3,1)
      X(6,1) = X(6,1) +      AT(1,3)
      X(1,4) = X(1,4) +      AT(2,1)
      X(2,4) = X(2,4) +      AT(1,2)
      X(2,5) = X(2,5) +      AT(3,2)
      X(3,5) = X(3,5) +      AT(2,3)
      X(3,6) = X(3,6) +      AT(1,3)
      X(1,6) = X(1,6) +      AT(3,1)
      X(4,4) = X(4,4) + r1d2 * (AT(1,1)+AT(2,2))
      X(4,5) = X(4,5) + r1d2 *  AT(3,1)
      X(4,6) = X(4,6) + r1d2 *  AT(3,2)
      X(5,4) = X(5,4) + r1d2 *  AT(1,3)
      X(5,5) = X(5,5) + r1d2 * (AT(2,2)+AT(3,3))
      X(5,6) = X(5,6) + r1d2 *  AT(1,2)
      X(6,4) = X(6,4) + r1d2 *  AT(2,3)
      X(6,5) = X(6,5) + r1d2 *  AT(2,1)
      X(6,6) = X(6,6) + r1d2 * (AT(3,3)+AT(1,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAT_am(X,AT)
      implicit none
      double precision X(6,6), AT(3,3), r2, r1d2
      data r1d2, r2 / 0.5d0, 2.d0 /
      X(1,1) = X(1,1) - r2 * AT(1,1)
      X(2,2) = X(2,2) - r2 * AT(2,2)
      X(3,3) = X(3,3) - r2 * AT(3,3)
      X(4,1) = X(4,1) -      AT(1,2)
      X(4,2) = X(4,2) -      AT(2,1)
      X(5,2) = X(5,2) -      AT(2,3)
      X(5,3) = X(5,3) -      AT(3,2)
      X(6,3) = X(6,3) -      AT(3,1)
      X(6,1) = X(6,1) -      AT(1,3)
      X(1,4) = X(1,4) -      AT(2,1)
      X(2,4) = X(2,4) -      AT(1,2)
      X(2,5) = X(2,5) -      AT(3,2)
      X(3,5) = X(3,5) -      AT(2,3)
      X(3,6) = X(3,6) -      AT(1,3)
      X(1,6) = X(1,6) -      AT(3,1)
      X(4,4) = X(4,4) - r1d2 * (AT(1,1)+AT(2,2))
      X(4,5) = X(4,5) - r1d2 *  AT(3,1)
      X(4,6) = X(4,6) - r1d2 *  AT(3,2)
      X(5,4) = X(5,4) - r1d2 *  AT(1,3)
      X(5,5) = X(5,5) - r1d2 * (AT(2,2)+AT(3,3))
      X(5,6) = X(5,6) - r1d2 *  AT(1,2)
      X(6,4) = X(6,4) - r1d2 *  AT(2,3)
      X(6,5) = X(6,5) - r1d2 *  AT(2,1)
      X(6,6) = X(6,6) - r1d2 * (AT(3,3)+AT(1,1))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAT_af(X,AT,fact)
      implicit none
      double precision X(6,6), AT(3,3), fact, factm2, factd2
      factm2 = 2.0d0 * fact
      factd2 = 0.5d0 * fact
      X(1,1) = X(1,1) + factm2 * AT(1,1)
      X(2,2) = X(2,2) + factm2 * AT(2,2)
      X(3,3) = X(3,3) + factm2 * AT(3,3)
      X(4,1) = X(4,1) + fact   * AT(1,2)
      X(4,2) = X(4,2) + fact   * AT(2,1)
      X(5,2) = X(5,2) + fact   * AT(2,3)
      X(5,3) = X(5,3) + fact   * AT(3,2)
      X(6,3) = X(6,3) + fact   * AT(3,1)
      X(6,1) = X(6,1) + fact   * AT(1,3)
      X(1,4) = X(1,4) + fact   * AT(2,1)
      X(2,4) = X(2,4) + fact   * AT(1,2)
      X(2,5) = X(2,5) + fact   * AT(3,2)
      X(3,5) = X(3,5) + fact   * AT(2,3)
      X(3,6) = X(3,6) + fact   * AT(1,3)
      X(1,6) = X(1,6) + fact   * AT(3,1)
      X(4,4) = X(4,4) + factd2 *(AT(1,1)+AT(2,2))
      X(4,5) = X(4,5) + factd2 * AT(3,1)
      X(4,6) = X(4,6) + factd2 * AT(3,2)
      X(5,4) = X(5,4) + factd2 * AT(1,3)
      X(5,5) = X(5,5) + factd2 *(AT(2,2)+AT(3,3))
      X(5,6) = X(5,6) + factd2 * AT(1,2)
      X(6,4) = X(6,4) + factd2 * AT(2,3)
      X(6,5) = X(6,5) + factd2 * AT(2,1)
      X(6,6) = X(6,6) + factd2 *(AT(3,3)+AT(1,1))
      return
      end
c
c=============================================================================
c  Differentiation
c                      D ATA      mit  A,T,ATA  sym. 2.Stufe
c                X =  -------       
c                       D T       -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_ATA_p(X,A)
      implicit none
      double precision A(6), X(6,6), r1d2
      r1d2 = 0.5d0

      X(1,1) = A(1)*A(1)
      X(2,2) = A(2)*A(2)
      X(3,3) = A(3)*A(3)

      X(1,2) = A(4)*A(4)
      X(2,1) = X(1,2)
      X(2,3) = A(5)*A(5)
      X(3,2) = X(2,3)
      X(1,3) = A(6)*A(6)
      X(3,1) = X(1,3)

      X(4,4) = r1d2 * (A(4)*A(4) + A(1)*A(2))
      X(5,5) = r1d2 * (A(5)*A(5) + A(2)*A(3))
      X(6,6) = r1d2 * (A(6)*A(6) + A(1)*A(3))

      X(4,5) = r1d2 * (A(4)*A(5) + A(2)*A(6))
      X(5,4) = X(4,5)
      X(5,6) = r1d2 * (A(5)*A(6) + A(4)*A(3))
      X(6,5) = X(5,6)
      X(4,6) = r1d2 * (A(4)*A(6) + A(1)*A(5))
      X(6,4) = X(4,6)

      X(4,1) = A(1)*A(4)
      X(1,4) = X(4,1)
      X(4,2) = A(2)*A(4)
      X(2,4) = X(4,2)
      X(4,3) = A(5)*A(6)
      X(3,4) = X(4,3)
      X(5,1) = A(4)*A(6)
      X(1,5) = X(5,1)
      X(5,2) = A(2)*A(5)
      X(2,5) = X(5,2)
      X(5,3) = A(3)*A(5)
      X(3,5) = X(5,3)
      X(6,1) = A(1)*A(6)
      X(1,6) = X(6,1)
      X(6,2) = A(4)*A(5)
      X(2,6) = X(6,2)
      X(6,3) = A(3)*A(6)
      X(3,6) = X(6,3)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATA_m(X,A)
      implicit none
      double precision A(6), X(6,6), r1d2
      r1d2 = 0.5d0

      X(1,1) = - A(1)*A(1)
      X(2,2) = - A(2)*A(2)
      X(3,3) = - A(3)*A(3)

      X(1,2) = - A(4)*A(4)
      X(2,1) =   X(1,2)
      X(2,3) = - A(5)*A(5)
      X(3,2) =   X(2,3)
      X(1,3) = - A(6)*A(6)
      X(3,1) =   X(1,3)

      X(4,4) = - r1d2 * (A(4)*A(4) + A(1)*A(2))
      X(5,5) = - r1d2 * (A(5)*A(5) + A(2)*A(3))
      X(6,6) = - r1d2 * (A(6)*A(6) + A(1)*A(3))

      X(4,5) = - r1d2 * (A(4)*A(5) + A(2)*A(6))
      X(5,4) =   X(4,5)
      X(5,6) = - r1d2 * (A(5)*A(6) + A(4)*A(3))
      X(6,5) =   X(5,6)
      X(4,6) = - r1d2 * (A(4)*A(6) + A(1)*A(5))
      X(6,4) =   X(4,6)

      X(4,1) = - A(1)*A(4)
      X(1,4) =   X(4,1)
      X(4,2) = - A(2)*A(4)
      X(2,4) =   X(4,2)
      X(4,3) = - A(5)*A(6)
      X(3,4) =   X(4,3)
      X(5,1) = - A(4)*A(6)
      X(1,5) =   X(5,1)
      X(5,2) = - A(2)*A(5)
      X(2,5) =   X(5,2)
      X(5,3) = - A(3)*A(5)
      X(3,5) =   X(5,3)
      X(6,1) = - A(1)*A(6)
      X(1,6) =   X(6,1)
      X(6,2) = - A(4)*A(5)
      X(2,6) =   X(6,2)
      X(6,3) = - A(3)*A(6)
      X(3,6) =   X(6,3)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATA_f(X,A,fact)
      implicit none
      double precision A(6), X(6,6), fact, factd2
      factd2 = fact / 2.d0

      X(1,1) = fact * A(1)*A(1)
      X(2,2) = fact * A(2)*A(2)
      X(3,3) = fact * A(3)*A(3)

      X(1,2) = fact * A(4)*A(4)
      X(2,1) = X(1,2)
      X(2,3) = fact * A(5)*A(5)
      X(3,2) = X(2,3)
      X(1,3) = fact * A(6)*A(6)
      X(3,1) = X(1,3)

      X(4,4) = factd2 * (A(4)*A(4) + A(1)*A(2))
      X(5,5) = factd2 * (A(5)*A(5) + A(2)*A(3))
      X(6,6) = factd2 * (A(6)*A(6) + A(1)*A(3))

      X(4,5) = factd2 * (A(4)*A(5) + A(2)*A(6))
      X(5,4) = X(4,5)
      X(5,6) = factd2 * (A(5)*A(6) + A(4)*A(3))
      X(6,5) = X(5,6)
      X(4,6) = factd2 * (A(4)*A(6) + A(1)*A(5))
      X(6,4) = X(4,6)

      X(4,1) = fact * A(1)*A(4)
      X(1,4) = X(4,1)
      X(4,2) = fact * A(2)*A(4)
      X(2,4) = X(4,2)
      X(4,3) = fact * A(5)*A(6)
      X(3,4) = X(4,3)
      X(5,1) = fact * A(4)*A(6)
      X(1,5) = X(5,1)
      X(5,2) = fact * A(2)*A(5)
      X(2,5) = X(5,2)
      X(5,3) = fact * A(3)*A(5)
      X(3,5) = X(5,3)
      X(6,1) = fact * A(1)*A(6)
      X(1,6) = X(6,1)
      X(6,2) = fact * A(4)*A(5)
      X(2,6) = X(6,2)
      X(6,3) = fact * A(3)*A(6)
      X(3,6) = X(6,3)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATA_ap(X,A)
      implicit none
      double precision A(6), X(6,6), r1d2, z
      r1d2 = 0.5d0

      X(1,1) = X(1,1) + A(1)*A(1)
      X(2,2) = X(2,2) + A(2)*A(2)
      X(3,3) = X(3,3) + A(3)*A(3)

      z = A(4)*A(4)
      X(1,2) = X(1,2) + z
      X(2,1) = X(2,1) + z
      z = A(5)*A(5)
      X(2,3) = X(2,3) + z
      X(3,2) = X(3,2) + z
      z = A(6)*A(6)
      X(1,3) = X(1,3) + z
      X(3,1) = X(3,1) + z

      X(4,4) = X(4,4) + r1d2 * (A(4)*A(4) + A(1)*A(2))
      X(5,5) = X(5,5) + r1d2 * (A(5)*A(5) + A(2)*A(3))
      X(6,6) = X(6,6) + r1d2 * (A(6)*A(6) + A(1)*A(3))

      z = r1d2 * (A(4)*A(5) + A(2)*A(6))
      X(4,5) = X(4,5) + z
      X(5,4) = X(5,4) + z
      z = r1d2 * (A(5)*A(6) + A(4)*A(3))
      X(5,6) = X(5,6) + z
      X(6,5) = X(6,5) + z
      z = r1d2 * (A(4)*A(6) + A(1)*A(5))
      X(4,6) = X(4,6) + z
      X(6,4) = X(6,4) + z

      z = A(1)*A(4)
      X(4,1) = X(4,1) + z
      X(1,4) = X(1,4) + z
      z = A(2)*A(4)
      X(4,2) = X(4,2) + z
      X(2,4) = X(2,4) + z
      z = A(5)*A(6)
      X(4,3) = X(4,3) + z
      X(3,4) = X(3,4) + z
      z = A(4)*A(6)
      X(5,1) = X(5,1) + z
      X(1,5) = X(1,5) + z
      z = A(2)*A(5)
      X(5,2) = X(5,2) + z
      X(2,5) = X(2,5) + z
      z = A(3)*A(5)
      X(5,3) = X(5,3) + z
      X(3,5) = X(3,5) + z
      z = A(1)*A(6)
      X(6,1) = X(6,1) + z
      X(1,6) = X(1,6) + z
      z = A(4)*A(5)
      X(6,2) = X(6,2) + z
      X(2,6) = X(2,6) + z
      z = A(3)*A(6)
      X(6,3) = X(6,3) + z
      X(3,6) = X(3,6) + z

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATA_am(X,A)
      implicit none
      double precision A(6), X(6,6), r1d2, z
      r1d2 = 0.5d0

      X(1,1) = X(1,1) - A(1)*A(1)
      X(2,2) = X(2,2) - A(2)*A(2)
      X(3,3) = X(3,3) - A(3)*A(3)

      z = A(4)*A(4)
      X(1,2) = X(1,2) - z
      X(2,1) = X(2,1) - z
      z = A(5)*A(5)
      X(2,3) = X(2,3) - z
      X(3,2) = X(3,2) - z
      z = A(6)*A(6)
      X(1,3) = X(1,3) - z
      X(3,1) = X(3,1) - z

      X(4,4) = X(4,4) - r1d2 * (A(4)*A(4) + A(1)*A(2))
      X(5,5) = X(5,5) - r1d2 * (A(5)*A(5) + A(2)*A(3))
      X(6,6) = X(6,6) - r1d2 * (A(6)*A(6) + A(1)*A(3))

      z = r1d2 * (A(4)*A(5) + A(2)*A(6))
      X(4,5) = X(4,5) - z
      X(5,4) = X(5,4) - z
      z = r1d2 * (A(5)*A(6) + A(4)*A(3))
      X(5,6) = X(5,6) - z
      X(6,5) = X(6,5) - z
      z = r1d2 * (A(4)*A(6) + A(1)*A(5))
      X(4,6) = X(4,6) - z
      X(6,4) = X(6,4) - z

      z = A(1)*A(4)
      X(4,1) = X(4,1) - z
      X(1,4) = X(1,4) - z
      z = A(2)*A(4)
      X(4,2) = X(4,2) - z
      X(2,4) = X(2,4) - z
      z = A(5)*A(6)
      X(4,3) = X(4,3) - z
      X(3,4) = X(3,4) - z
      z = A(4)*A(6)
      X(5,1) = X(5,1) - z
      X(1,5) = X(1,5) - z
      z = A(2)*A(5)
      X(5,2) = X(5,2) - z
      X(2,5) = X(2,5) - z
      z = A(3)*A(5)
      X(5,3) = X(5,3) - z
      X(3,5) = X(3,5) - z
      z = A(1)*A(6)
      X(6,1) = X(6,1) - z
      X(1,6) = X(1,6) - z
      z = A(4)*A(5)
      X(6,2) = X(6,2) - z
      X(2,6) = X(2,6) - z
      z = A(3)*A(6)
      X(6,3) = X(6,3) - z
      X(3,6) = X(3,6) - z

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATA_af(X,A,fact)
      implicit none
      double precision A(6), X(6,6), fact, factd2, z
      factd2 = fact / 2.d0

      X(1,1) = X(1,1) + fact * A(1)*A(1)
      X(2,2) = X(2,2) + fact * A(2)*A(2)
      X(3,3) = X(3,3) + fact * A(3)*A(3)

      z = fact * A(4)*A(4)
      X(1,2) = X(1,2) + z
      X(2,1) = X(2,1) + z
      z = fact * A(5)*A(5)
      X(2,3) = X(2,3) + z
      X(3,2) = X(3,2) + z
      z = fact * A(6)*A(6)
      X(1,3) = X(1,3) + z
      X(3,1) = X(3,1) + z

      X(4,4) = X(4,4) + factd2 * (A(4)*A(4) + A(1)*A(2))
      X(5,5) = X(5,5) + factd2 * (A(5)*A(5) + A(2)*A(3))
      X(6,6) = X(6,6) + factd2 * (A(6)*A(6) + A(1)*A(3))

      z = factd2 * (A(4)*A(5) + A(2)*A(6))
      X(4,5) = X(4,5) + z
      X(5,4) = X(5,4) + z
      z = factd2 * (A(5)*A(6) + A(4)*A(3))
      X(5,6) = X(5,6) + z
      X(6,5) = X(6,5) + z
      z = factd2 * (A(4)*A(6) + A(1)*A(5))
      X(4,6) = X(4,6) + z
      X(6,4) = X(6,4) + z

      z = fact * A(1)*A(4)
      X(4,1) = X(4,1) + z
      X(1,4) = X(1,4) + z
      z = fact * A(2)*A(4)
      X(4,2) = X(4,2) + z
      X(2,4) = X(2,4) + z
      z = fact * A(5)*A(6)
      X(4,3) = X(4,3) + z
      X(3,4) = X(3,4) + z
      z = fact * A(4)*A(6)
      X(5,1) = X(5,1) + z
      X(1,5) = X(1,5) + z
      z = fact * A(2)*A(5)
      X(5,2) = X(5,2) + z
      X(2,5) = X(2,5) + z
      z = fact * A(3)*A(5)
      X(5,3) = X(5,3) + z
      X(3,5) = X(3,5) + z
      z = fact * A(1)*A(6)
      X(6,1) = X(6,1) + z
      X(1,6) = X(1,6) + z
      z = fact * A(4)*A(5)
      X(6,2) = X(6,2) + z
      X(2,6) = X(2,6) + z
      z = fact * A(3)*A(6)
      X(6,3) = X(6,3) + z
      X(3,6) = X(3,6) + z

      return
      end
c
c=============================================================================
c   Differentiation
c                     D ATA^t       mit  T,ATA^t sym., A unsym. 2.Stufe
c                X = ---------       
c                       D T         -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_UTUt_p(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d2
      r1d2 = 0.5d0

      X(1,1) = A(1,1)*A(1,1)
      X(2,2) = A(2,2)*A(2,2)
      X(3,3) = A(3,3)*A(3,3)

      X(1,2) = A(1,2)*A(1,2)
      X(2,1) = A(2,1)*A(2,1)
      X(2,3) = A(2,3)*A(2,3)
      X(3,2) = A(3,2)*A(3,2)
      X(1,3) = A(1,3)*A(1,3)
      X(3,1) = A(3,1)*A(3,1)

      X(4,4) = r1d2 * (A(1,2)*A(2,1) + A(1,1)*A(2,2))
      X(5,5) = r1d2 * (A(2,3)*A(3,2) + A(2,2)*A(3,3))
      X(6,6) = r1d2 * (A(1,3)*A(3,1) + A(1,1)*A(3,3))

      X(4,5) = r1d2 * (A(1,2)*A(2,3) + A(1,3)*A(2,2))
      X(5,4) = r1d2 * (A(2,1)*A(3,2) + A(3,1)*A(2,2))
      X(5,6) = r1d2 * (A(2,3)*A(3,1) + A(2,1)*A(3,3))
      X(6,5) = r1d2 * (A(3,2)*A(1,3) + A(1,2)*A(3,3))
      X(4,6) = r1d2 * (A(1,3)*A(2,1) + A(2,3)*A(1,1))
      X(6,4) = r1d2 * (A(3,1)*A(1,2) + A(3,2)*A(1,1))

      X(4,1) = A(1,1)*A(2,1)
      X(1,4) = A(1,1)*A(1,2)
      X(4,2) = A(2,2)*A(1,2)
      X(2,4) = A(2,2)*A(2,1)
      X(4,3) = A(1,3)*A(2,3)
      X(3,4) = A(3,1)*A(3,2)
      X(5,1) = A(2,1)*A(3,1)
      X(1,5) = A(1,2)*A(1,3)
      X(5,2) = A(2,2)*A(3,2)
      X(2,5) = A(2,2)*A(2,3)
      X(5,3) = A(3,3)*A(2,3)
      X(3,5) = A(3,3)*A(3,2)
      X(6,1) = A(1,1)*A(3,1)
      X(1,6) = A(1,1)*A(1,3)
      X(6,2) = A(1,2)*A(3,2)
      X(2,6) = A(2,1)*A(2,3)
      X(6,3) = A(3,3)*A(1,3)
      X(3,6) = A(3,3)*A(3,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTUt_m(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d2
      r1d2 = 0.5d0

      X(1,1) = - A(1,1)*A(1,1)
      X(2,2) = - A(2,2)*A(2,2)
      X(3,3) = - A(3,3)*A(3,3)

      X(1,2) = - A(1,2)*A(1,2)
      X(2,1) = - A(2,1)*A(2,1)
      X(2,3) = - A(2,3)*A(2,3)
      X(3,2) = - A(3,2)*A(3,2)
      X(1,3) = - A(1,3)*A(1,3)
      X(3,1) = - A(3,1)*A(3,1)

      X(4,4) = - r1d2 * (A(1,2)*A(2,1) + A(1,1)*A(2,2))
      X(5,5) = - r1d2 * (A(2,3)*A(3,2) + A(2,2)*A(3,3))
      X(6,6) = - r1d2 * (A(1,3)*A(3,1) + A(1,1)*A(3,3))

      X(4,5) = - r1d2 * (A(1,2)*A(2,3) + A(1,3)*A(2,2))
      X(5,4) = - r1d2 * (A(2,1)*A(3,2) + A(3,1)*A(2,2))
      X(5,6) = - r1d2 * (A(2,3)*A(3,1) + A(2,1)*A(3,3))
      X(6,5) = - r1d2 * (A(3,2)*A(1,3) + A(1,2)*A(3,3))
      X(4,6) = - r1d2 * (A(1,3)*A(2,1) + A(2,3)*A(1,1))
      X(6,4) = - r1d2 * (A(3,1)*A(1,2) + A(3,2)*A(1,1))

      X(4,1) = - A(1,1)*A(2,1)
      X(1,4) = - A(1,1)*A(1,2)
      X(4,2) = - A(2,2)*A(1,2)
      X(2,4) = - A(2,2)*A(2,1)
      X(4,3) = - A(1,3)*A(2,3)
      X(3,4) = - A(3,1)*A(3,2)
      X(5,1) = - A(2,1)*A(3,1)
      X(1,5) = - A(1,2)*A(1,3)
      X(5,2) = - A(2,2)*A(3,2)
      X(2,5) = - A(2,2)*A(2,3)
      X(5,3) = - A(3,3)*A(2,3)
      X(3,5) = - A(3,3)*A(3,2)
      X(6,1) = - A(1,1)*A(3,1)
      X(1,6) = - A(1,1)*A(1,3)
      X(6,2) = - A(1,2)*A(3,2)
      X(2,6) = - A(2,1)*A(2,3)
      X(6,3) = - A(3,3)*A(1,3)
      X(3,6) = - A(3,3)*A(3,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTUt_f(X,A,fact)
      implicit none
      double precision A(3,3), X(6,6), fact, factd2
      factd2 = fact / 2.d0

      X(1,1) = fact * A(1,1)*A(1,1)
      X(2,2) = fact * A(2,2)*A(2,2)
      X(3,3) = fact * A(3,3)*A(3,3)

      X(1,2) = fact * A(1,2)*A(1,2)
      X(2,1) = fact * A(2,1)*A(2,1)
      X(2,3) = fact * A(2,3)*A(2,3)
      X(3,2) = fact * A(3,2)*A(3,2)
      X(1,3) = fact * A(1,3)*A(1,3)
      X(3,1) = fact * A(3,1)*A(3,1)

      X(4,4) = factd2 * (A(1,2)*A(2,1) + A(1,1)*A(2,2))
      X(5,5) = factd2 * (A(2,3)*A(3,2) + A(2,2)*A(3,3))
      X(6,6) = factd2 * (A(1,3)*A(3,1) + A(1,1)*A(3,3))

      X(4,5) = factd2 * (A(1,2)*A(2,3) + A(1,3)*A(2,2))
      X(5,4) = factd2 * (A(2,1)*A(3,2) + A(3,1)*A(2,2))
      X(5,6) = factd2 * (A(2,3)*A(3,1) + A(2,1)*A(3,3))
      X(6,5) = factd2 * (A(3,2)*A(1,3) + A(1,2)*A(3,3))
      X(4,6) = factd2 * (A(1,3)*A(2,1) + A(2,3)*A(1,1))
      X(6,4) = factd2 * (A(3,1)*A(1,2) + A(3,2)*A(1,1))

      X(4,1) = fact * A(1,1)*A(2,1)
      X(1,4) = fact * A(1,1)*A(1,2)
      X(4,2) = fact * A(2,2)*A(1,2)
      X(2,4) = fact * A(2,2)*A(2,1)
      X(4,3) = fact * A(1,3)*A(2,3)
      X(3,4) = fact * A(3,1)*A(3,2)
      X(5,1) = fact * A(2,1)*A(3,1)
      X(1,5) = fact * A(1,2)*A(1,3)
      X(5,2) = fact * A(2,2)*A(3,2)
      X(2,5) = fact * A(2,2)*A(2,3)
      X(5,3) = fact * A(3,3)*A(2,3)
      X(3,5) = fact * A(3,3)*A(3,2)
      X(6,1) = fact * A(1,1)*A(3,1)
      X(1,6) = fact * A(1,1)*A(1,3)
      X(6,2) = fact * A(1,2)*A(3,2)
      X(2,6) = fact * A(2,1)*A(2,3)
      X(6,3) = fact * A(3,3)*A(1,3)
      X(3,6) = fact * A(3,3)*A(3,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTUt_ap(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d2
      r1d2 = 0.5d0

      X(1,1) = X(1,1) + A(1,1)*A(1,1)
      X(2,2) = X(2,2) + A(2,2)*A(2,2)
      X(3,3) = X(3,3) + A(3,3)*A(3,3)

      X(1,2) = X(1,2) + A(1,2)*A(1,2)
      X(2,1) = X(2,1) + A(2,1)*A(2,1)
      X(2,3) = X(2,3) + A(2,3)*A(2,3)
      X(3,2) = X(3,2) + A(3,2)*A(3,2)
      X(1,3) = X(1,3) + A(1,3)*A(1,3)
      X(3,1) = X(3,1) + A(3,1)*A(3,1)

      X(4,4) = X(4,4) + r1d2 * (A(1,2)*A(2,1) + A(1,1)*A(2,2))
      X(5,5) = X(5,5) + r1d2 * (A(2,3)*A(3,2) + A(2,2)*A(3,3))
      X(6,6) = X(6,6) + r1d2 * (A(1,3)*A(3,1) + A(1,1)*A(3,3))

      X(4,5) = X(4,5) + r1d2 * (A(1,2)*A(2,3) + A(1,3)*A(2,2))
      X(5,4) = X(5,4) + r1d2 * (A(2,1)*A(3,2) + A(3,1)*A(2,2))
      X(5,6) = X(5,6) + r1d2 * (A(2,3)*A(3,1) + A(2,1)*A(3,3))
      X(6,5) = X(6,5) + r1d2 * (A(3,2)*A(1,3) + A(1,2)*A(3,3))
      X(4,6) = X(4,6) + r1d2 * (A(1,3)*A(2,1) + A(2,3)*A(1,1))
      X(6,4) = X(6,4) + r1d2 * (A(3,1)*A(1,2) + A(3,2)*A(1,1))

      X(4,1) = X(4,1) + A(1,1)*A(2,1)
      X(1,4) = X(1,4) + A(1,1)*A(1,2)
      X(4,2) = X(4,2) + A(2,2)*A(1,2)
      X(2,4) = X(2,4) + A(2,2)*A(2,1)
      X(4,3) = X(4,3) + A(1,3)*A(2,3)
      X(3,4) = X(3,4) + A(3,1)*A(3,2)
      X(5,1) = X(5,1) + A(2,1)*A(3,1)
      X(1,5) = X(1,5) + A(1,2)*A(1,3)
      X(5,2) = X(5,2) + A(2,2)*A(3,2)
      X(2,5) = X(2,5) + A(2,2)*A(2,3)
      X(5,3) = X(5,3) + A(3,3)*A(2,3)
      X(3,5) = X(3,5) + A(3,3)*A(3,2)
      X(6,1) = X(6,1) + A(1,1)*A(3,1)
      X(1,6) = X(1,6) + A(1,1)*A(1,3)
      X(6,2) = X(6,2) + A(1,2)*A(3,2)
      X(2,6) = X(2,6) + A(2,1)*A(2,3)
      X(6,3) = X(6,3) + A(3,3)*A(1,3)
      X(3,6) = X(3,6) + A(3,3)*A(3,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTUt_am(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d2
      r1d2 = 0.5d0

      X(1,1) = X(1,1) - A(1,1)*A(1,1)
      X(2,2) = X(2,2) - A(2,2)*A(2,2)
      X(3,3) = X(3,3) - A(3,3)*A(3,3)

      X(1,2) = X(1,2) - A(1,2)*A(1,2)
      X(2,1) = X(2,1) - A(2,1)*A(2,1)
      X(2,3) = X(2,3) - A(2,3)*A(2,3)
      X(3,2) = X(3,2) - A(3,2)*A(3,2)
      X(1,3) = X(1,3) - A(1,3)*A(1,3)
      X(3,1) = X(3,1) - A(3,1)*A(3,1)

      X(4,4) = X(4,4) - r1d2 * (A(1,2)*A(2,1) + A(1,1)*A(2,2))
      X(5,5) = X(5,5) - r1d2 * (A(2,3)*A(3,2) + A(2,2)*A(3,3))
      X(6,6) = X(6,6) - r1d2 * (A(1,3)*A(3,1) + A(1,1)*A(3,3))

      X(4,5) = X(4,5) - r1d2 * (A(1,2)*A(2,3) + A(1,3)*A(2,2))
      X(5,4) = X(5,4) - r1d2 * (A(2,1)*A(3,2) + A(3,1)*A(2,2))
      X(5,6) = X(5,6) - r1d2 * (A(2,3)*A(3,1) + A(2,1)*A(3,3))
      X(6,5) = X(6,5) - r1d2 * (A(3,2)*A(1,3) + A(1,2)*A(3,3))
      X(4,6) = X(4,6) - r1d2 * (A(1,3)*A(2,1) + A(2,3)*A(1,1))
      X(6,4) = X(6,4) - r1d2 * (A(3,1)*A(1,2) + A(3,2)*A(1,1))

      X(4,1) = X(4,1) - A(1,1)*A(2,1)
      X(1,4) = X(1,4) - A(1,1)*A(1,2)
      X(4,2) = X(4,2) - A(2,2)*A(1,2)
      X(2,4) = X(2,4) - A(2,2)*A(2,1)
      X(4,3) = X(4,3) - A(1,3)*A(2,3)
      X(3,4) = X(3,4) - A(3,1)*A(3,2)
      X(5,1) = X(5,1) - A(2,1)*A(3,1)
      X(1,5) = X(1,5) - A(1,2)*A(1,3)
      X(5,2) = X(5,2) - A(2,2)*A(3,2)
      X(2,5) = X(2,5) - A(2,2)*A(2,3)
      X(5,3) = X(5,3) - A(3,3)*A(2,3)
      X(3,5) = X(3,5) - A(3,3)*A(3,2)
      X(6,1) = X(6,1) - A(1,1)*A(3,1)
      X(1,6) = X(1,6) - A(1,1)*A(1,3)
      X(6,2) = X(6,2) - A(1,2)*A(3,2)
      X(2,6) = X(2,6) - A(2,1)*A(2,3)
      X(6,3) = X(6,3) - A(3,3)*A(1,3)
      X(3,6) = X(3,6) - A(3,3)*A(3,1)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTUt_af(X,A,fact)
      implicit none
      double precision A(3,3), X(6,6), fact, factd2
      factd2 = fact / 2.d0

      X(1,1) = X(1,1) + fact * A(1,1)*A(1,1)
      X(2,2) = X(2,2) + fact * A(2,2)*A(2,2)
      X(3,3) = X(3,3) + fact * A(3,3)*A(3,3)

      X(1,2) = X(1,2) + fact * A(1,2)*A(1,2)
      X(2,1) = X(2,1) + fact * A(2,1)*A(2,1)
      X(2,3) = X(2,3) + fact * A(2,3)*A(2,3)
      X(3,2) = X(3,2) + fact * A(3,2)*A(3,2)
      X(1,3) = X(1,3) + fact * A(1,3)*A(1,3)
      X(3,1) = X(3,1) + fact * A(3,1)*A(3,1)

      X(4,4) = X(4,4) + factd2 * (A(1,2)*A(2,1) + A(1,1)*A(2,2))
      X(5,5) = X(5,5) + factd2 * (A(2,3)*A(3,2) + A(2,2)*A(3,3))
      X(6,6) = X(6,6) + factd2 * (A(1,3)*A(3,1) + A(1,1)*A(3,3))

      X(4,5) = X(4,5) + factd2 * (A(1,2)*A(2,3) + A(1,3)*A(2,2))
      X(5,4) = X(5,4) + factd2 * (A(2,1)*A(3,2) + A(3,1)*A(2,2))
      X(5,6) = X(5,6) + factd2 * (A(2,3)*A(3,1) + A(2,1)*A(3,3))
      X(6,5) = X(6,5) + factd2 * (A(3,2)*A(1,3) + A(1,2)*A(3,3))
      X(4,6) = X(4,6) + factd2 * (A(1,3)*A(2,1) + A(2,3)*A(1,1))
      X(6,4) = X(6,4) + factd2 * (A(3,1)*A(1,2) + A(3,2)*A(1,1))

      X(4,1) = X(4,1) + fact * A(1,1)*A(2,1)
      X(1,4) = X(1,4) + fact * A(1,1)*A(1,2)
      X(4,2) = X(4,2) + fact * A(2,2)*A(1,2)
      X(2,4) = X(2,4) + fact * A(2,2)*A(2,1)
      X(4,3) = X(4,3) + fact * A(1,3)*A(2,3)
      X(3,4) = X(3,4) + fact * A(3,1)*A(3,2)
      X(5,1) = X(5,1) + fact * A(2,1)*A(3,1)
      X(1,5) = X(1,5) + fact * A(1,2)*A(1,3)
      X(5,2) = X(5,2) + fact * A(2,2)*A(3,2)
      X(2,5) = X(2,5) + fact * A(2,2)*A(2,3)
      X(5,3) = X(5,3) + fact * A(3,3)*A(2,3)
      X(3,5) = X(3,5) + fact * A(3,3)*A(3,2)
      X(6,1) = X(6,1) + fact * A(1,1)*A(3,1)
      X(1,6) = X(1,6) + fact * A(1,1)*A(1,3)
      X(6,2) = X(6,2) + fact * A(1,2)*A(3,2)
      X(2,6) = X(2,6) + fact * A(2,1)*A(2,3)
      X(6,3) = X(6,3) + fact * A(3,3)*A(1,3)
      X(3,6) = X(3,6) + fact * A(3,3)*A(3,1)

      return
      end
c
c=============================================================================
c   Differentiation
c                        D ATB        mit  A,B,T,ATB sym. 2.Stufe
c                   X = -------       
c                         D T         -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_ATB_p(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d2, r1d4
      r1d2 = 0.50d0
      r1d4 = 0.25d0

      X(1,1) = A(1)*B(1)
      X(2,2) = A(2)*B(2)
      X(3,3) = A(3)*B(3)

      X(1,2) = A(4)*B(4)
      X(2,1) = X(1,2)
      X(2,3) = A(5)*B(5)
      X(3,2) = X(2,3)
      X(1,3) = A(6)*B(6)
      X(3,1) = X(1,3)

      X(4,4) = r1d4 * (A(1)*B(2) + A(4)*B(4) + A(4)*B(4) + A(2)*B(1))
      X(4,5) = r1d4 * (A(4)*B(5) + A(6)*B(2) + A(2)*B(6) + A(5)*B(4))
      X(4,6) = r1d4 * (A(1)*B(5) + A(6)*B(4) + A(4)*B(6) + A(5)*B(1))
      X(5,4) = r1d4 * (A(4)*B(5) + A(2)*B(6) + A(6)*B(2) + A(5)*B(4))
      X(5,5) = r1d4 * (A(2)*B(3) + A(5)*B(5) + A(5)*B(5) + A(3)*B(2))
      X(5,6) = r1d4 * (A(4)*B(3) + A(5)*B(6) + A(6)*B(5) + A(3)*B(4))
      X(6,4) = r1d4 * (A(1)*B(5) + A(4)*B(6) + A(6)*B(4) + A(5)*B(1))
      X(6,5) = r1d4 * (A(4)*B(3) + A(6)*B(5) + A(5)*B(6) + A(3)*B(4))
      X(6,6) = r1d4 * (A(1)*B(3) + A(6)*B(6) + A(6)*B(6) + A(3)*B(1))

      X(4,1) = r1d2 * (A(1)*B(4) + A(4)*B(1))
      X(1,4) = X(4,1)
      X(4,2) = r1d2 * (A(4)*B(2) + A(2)*B(4))
      X(2,4) = X(4,2)
      X(4,3) = r1d2 * (A(6)*B(5) + A(5)*B(6))
      X(3,4) = X(4,3)
      X(5,1) = r1d2 * (A(4)*B(6) + A(6)*B(4))
      X(1,5) = X(5,1)
      X(5,2) = r1d2 * (A(2)*B(5) + A(5)*B(2))
      X(2,5) = X(5,2)
      X(5,3) = r1d2 * (A(5)*B(3) + A(3)*B(5))
      X(3,5) = X(5,3)
      X(6,1) = r1d2 * (A(1)*B(6) + A(6)*B(1))
      X(1,6) = X(6,1)
      X(6,2) = r1d2 * (A(4)*B(5) + A(5)*B(4))
      X(2,6) = X(6,2)
      X(6,3) = r1d2 * (A(6)*B(3) + A(3)*B(6))
      X(3,6) = X(6,3)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATB_m(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d2, r1d4
      r1d2 = 0.50d0
      r1d4 = 0.25d0

      X(1,1) = - A(1)*B(1)
      X(2,2) = - A(2)*B(2)
      X(3,3) = - A(3)*B(3)

      X(1,2) = - A(4)*B(4)
      X(2,1) = X(1,2)
      X(2,3) = - A(5)*B(5)
      X(3,2) = X(2,3)
      X(1,3) = - A(6)*B(6)
      X(3,1) = X(1,3)

      X(4,4) = - r1d4 * (A(1)*B(2) + A(4)*B(4) + A(4)*B(4) + A(2)*B(1))
      X(4,5) = - r1d4 * (A(4)*B(5) + A(6)*B(2) + A(2)*B(6) + A(5)*B(4))
      X(4,6) = - r1d4 * (A(1)*B(5) + A(6)*B(4) + A(4)*B(6) + A(5)*B(1))
      X(5,4) = - r1d4 * (A(4)*B(5) + A(2)*B(6) + A(6)*B(2) + A(5)*B(4))
      X(5,5) = - r1d4 * (A(2)*B(3) + A(5)*B(5) + A(5)*B(5) + A(3)*B(2))
      X(5,6) = - r1d4 * (A(4)*B(3) + A(5)*B(6) + A(6)*B(5) + A(3)*B(4))
      X(6,4) = - r1d4 * (A(1)*B(5) + A(4)*B(6) + A(6)*B(4) + A(5)*B(1))
      X(6,5) = - r1d4 * (A(4)*B(3) + A(6)*B(5) + A(5)*B(6) + A(3)*B(4))
      X(6,6) = - r1d4 * (A(1)*B(3) + A(6)*B(6) + A(6)*B(6) + A(3)*B(1))

      X(4,1) = - r1d2 * (A(1)*B(4) + A(4)*B(1))
      X(1,4) = X(4,1)
      X(4,2) = - r1d2 * (A(4)*B(2) + A(2)*B(4))
      X(2,4) = X(4,2)
      X(4,3) = - r1d2 * (A(6)*B(5) + A(5)*B(6))
      X(3,4) = X(4,3)
      X(5,1) = - r1d2 * (A(4)*B(6) + A(6)*B(4))
      X(1,5) = X(5,1)
      X(5,2) = - r1d2 * (A(2)*B(5) + A(5)*B(2))
      X(2,5) = X(5,2)
      X(5,3) = - r1d2 * (A(5)*B(3) + A(3)*B(5))
      X(3,5) = X(5,3)
      X(6,1) = - r1d2 * (A(1)*B(6) + A(6)*B(1))
      X(1,6) = X(6,1)
      X(6,2) = - r1d2 * (A(4)*B(5) + A(5)*B(4))
      X(2,6) = X(6,2)
      X(6,3) = - r1d2 * (A(6)*B(3) + A(3)*B(6))
      X(3,6) = X(6,3)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATB_f(X,A,B,fact)
      implicit none
      double precision A(6), B(6), X(6,6), fact, factd2, factd4
      factd2 = fact / 2.d0
      factd4 = fact / 4.d0

      X(1,1) = fact * A(1)*B(1)
      X(2,2) = fact * A(2)*B(2)
      X(3,3) = fact * A(3)*B(3)

      X(1,2) = fact * A(4)*B(4)
      X(2,1) = X(1,2)
      X(2,3) = fact * A(5)*B(5)
      X(3,2) = X(2,3)
      X(1,3) = fact * A(6)*B(6)
      X(3,1) = X(1,3)

      X(4,4) = factd4 * (A(1)*B(2) + A(4)*B(4) + A(4)*B(4) + A(2)*B(1))
      X(4,5) = factd4 * (A(4)*B(5) + A(6)*B(2) + A(2)*B(6) + A(5)*B(4))
      X(4,6) = factd4 * (A(1)*B(5) + A(6)*B(4) + A(4)*B(6) + A(5)*B(1))
      X(5,4) = factd4 * (A(4)*B(5) + A(2)*B(6) + A(6)*B(2) + A(5)*B(4))
      X(5,5) = factd4 * (A(2)*B(3) + A(5)*B(5) + A(5)*B(5) + A(3)*B(2))
      X(5,6) = factd4 * (A(4)*B(3) + A(5)*B(6) + A(6)*B(5) + A(3)*B(4))
      X(6,4) = factd4 * (A(1)*B(5) + A(4)*B(6) + A(6)*B(4) + A(5)*B(1))
      X(6,5) = factd4 * (A(4)*B(3) + A(6)*B(5) + A(5)*B(6) + A(3)*B(4))
      X(6,6) = factd4 * (A(1)*B(3) + A(6)*B(6) + A(6)*B(6) + A(3)*B(1))

      X(4,1) = factd2 * (A(1)*B(4) + A(4)*B(1))
      X(1,4) = X(4,1)
      X(4,2) = factd2 * (A(4)*B(2) + A(2)*B(4))
      X(2,4) = X(4,2)
      X(4,3) = factd2 * (A(6)*B(5) + A(5)*B(6))
      X(3,4) = X(4,3)
      X(5,1) = factd2 * (A(4)*B(6) + A(6)*B(4))
      X(1,5) = X(5,1)
      X(5,2) = factd2 * (A(2)*B(5) + A(5)*B(2))
      X(2,5) = X(5,2)
      X(5,3) = factd2 * (A(5)*B(3) + A(3)*B(5))
      X(3,5) = X(5,3)
      X(6,1) = factd2 * (A(1)*B(6) + A(6)*B(1))
      X(1,6) = X(6,1)
      X(6,2) = factd2 * (A(4)*B(5) + A(5)*B(4))
      X(2,6) = X(6,2)
      X(6,3) = factd2 * (A(6)*B(3) + A(3)*B(6))
      X(3,6) = X(6,3)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATB_ap(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d2, r1d4, z
      r1d2 = 0.50d0
      r1d4 = 0.25d0

      X(1,1) = X(1,1) + A(1)*B(1)
      X(2,2) = X(2,2) + A(2)*B(2)
      X(3,3) = X(3,3) + A(3)*B(3)

      z = A(4)*B(4)
      X(1,2) = X(1,2) + z
      X(2,1) = X(2,1) + z
      z = A(5)*B(5)
      X(2,3) = X(2,3) + z
      X(3,2) = X(3,2) + z
      z = A(6)*B(6)
      X(1,3) = X(1,3) + z
      X(3,1) = X(3,1) + z

      X(4,4) = X(4,4) + r1d4 * (A(1)*B(2) + A(4)*B(4) 
     *                        + A(4)*B(4) + A(2)*B(1))
      X(4,5) = X(4,5) + r1d4 * (A(4)*B(5) + A(6)*B(2)
     *                        + A(2)*B(6) + A(5)*B(4))
      X(4,6) = X(4,6) + r1d4 * (A(1)*B(5) + A(6)*B(4)
     *                        + A(4)*B(6) + A(5)*B(1))
      X(5,4) = X(5,4) + r1d4 * (A(4)*B(5) + A(2)*B(6)
     *                        + A(6)*B(2) + A(5)*B(4))
      X(5,5) = X(5,5) + r1d4 * (A(2)*B(3) + A(5)*B(5)
     *                        + A(5)*B(5) + A(3)*B(2))
      X(5,6) = X(5,6) + r1d4 * (A(4)*B(3) + A(5)*B(6)
     *                        + A(6)*B(5) + A(3)*B(4))
      X(6,4) = X(6,4) + r1d4 * (A(1)*B(5) + A(4)*B(6)
     *                        + A(6)*B(4) + A(5)*B(1))
      X(6,5) = X(6,5) + r1d4 * (A(4)*B(3) + A(6)*B(5)
     *                        + A(5)*B(6) + A(3)*B(4))
      X(6,6) = X(6,6) + r1d4 * (A(1)*B(3) + A(6)*B(6)
     *                        + A(6)*B(6) + A(3)*B(1))

      z = r1d2 * (A(1)*B(4) + A(4)*B(1))
      X(4,1) = X(4,1) + z
      X(1,4) = X(1,4) + z
      z = r1d2 * (A(4)*B(2) + A(2)*B(4))
      X(4,2) = X(4,2) + z
      X(2,4) = X(2,4) + z
      z = r1d2 * (A(6)*B(5) + A(5)*B(6))
      X(4,3) = X(4,3) + z
      X(3,4) = X(3,4) + z
      z = r1d2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(5,1) + z
      X(1,5) = X(1,5) + z
      z = r1d2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(5,2) + z
      X(2,5) = X(2,5) + z
      z = r1d2 * (A(5)*B(3) + A(3)*B(5))
      X(5,3) = X(5,3) + z
      X(3,5) = X(3,5) + z
      z = r1d2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(6,1) + z
      X(1,6) = X(1,6) + z
      z = r1d2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(6,2) + z
      X(2,6) = X(2,6) + z
      z = r1d2 * (A(6)*B(3) + A(3)*B(6))
      X(6,3) = X(6,3) + z
      X(3,6) = X(3,6) + z

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATB_am(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d2, r1d4, z
      r1d2 = 0.50d0
      r1d4 = 0.25d0

      X(1,1) = X(1,1) - A(1)*B(1)
      X(2,2) = X(2,2) - A(2)*B(2)
      X(3,3) = X(3,3) - A(3)*B(3)

      z = A(4)*B(4)
      X(1,2) = X(1,2) - z
      X(2,1) = X(2,1) - z
      z = A(5)*B(5)
      X(2,3) = X(2,3) - z
      X(3,2) = X(3,2) - z
      z = A(6)*B(6)
      X(1,3) = X(1,3) - z
      X(3,1) = X(3,1) - z

      X(4,4) = X(4,4) - r1d4 * (A(1)*B(2) + A(4)*B(4) 
     *                        + A(4)*B(4) + A(2)*B(1))
      X(4,5) = X(4,5) - r1d4 * (A(4)*B(5) + A(6)*B(2)
     *                        + A(2)*B(6) + A(5)*B(4))
      X(4,6) = X(4,6) - r1d4 * (A(1)*B(5) + A(6)*B(4)
     *                        + A(4)*B(6) + A(5)*B(1))
      X(5,4) = X(5,4) - r1d4 * (A(4)*B(5) + A(2)*B(6)
     *                        + A(6)*B(2) + A(5)*B(4))
      X(5,5) = X(5,5) - r1d4 * (A(2)*B(3) + A(5)*B(5)
     *                        + A(5)*B(5) + A(3)*B(2))
      X(5,6) = X(5,6) - r1d4 * (A(4)*B(3) + A(5)*B(6)
     *                        + A(6)*B(5) + A(3)*B(4))
      X(6,4) = X(6,4) - r1d4 * (A(1)*B(5) + A(4)*B(6)
     *                        + A(6)*B(4) + A(5)*B(1))
      X(6,5) = X(6,5) - r1d4 * (A(4)*B(3) + A(6)*B(5)
     *                        + A(5)*B(6) + A(3)*B(4))
      X(6,6) = X(6,6) - r1d4 * (A(1)*B(3) + A(6)*B(6)
     *                        + A(6)*B(6) + A(3)*B(1))

      z = r1d2 * (A(1)*B(4) + A(4)*B(1))
      X(4,1) = X(4,1) - z
      X(1,4) = X(1,4) - z
      z = r1d2 * (A(4)*B(2) + A(2)*B(4))
      X(4,2) = X(4,2) - z
      X(2,4) = X(2,4) - z
      z = r1d2 * (A(6)*B(5) + A(5)*B(6))
      X(4,3) = X(4,3) - z
      X(3,4) = X(3,4) - z
      z = r1d2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(5,1) - z
      X(1,5) = X(1,5) - z
      z = r1d2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(5,2) - z
      X(2,5) = X(2,5) - z
      z = r1d2 * (A(5)*B(3) + A(3)*B(5))
      X(5,3) = X(5,3) - z
      X(3,5) = X(3,5) - z
      z = r1d2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(6,1) - z
      X(1,6) = X(1,6) - z
      z = r1d2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(6,2) - z
      X(2,6) = X(2,6) - z
      z = r1d2 * (A(6)*B(3) + A(3)*B(6))
      X(6,3) = X(6,3) - z
      X(3,6) = X(3,6) - z

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATB_af(X,A,B,fact)
      implicit none
      double precision A(6), B(6), X(6,6), fact, factd2, factd4, z
      factd2 = fact / 2.d0
      factd4 = fact / 4.d0

      X(1,1) = X(1,1) + fact * A(1)*B(1)
      X(2,2) = X(2,2) + fact * A(2)*B(2)
      X(3,3) = X(3,3) + fact * A(3)*B(3)

      z = fact * A(4)*B(4)
      X(1,2) = X(1,2) + z
      X(2,1) = X(2,1) + z
      z = fact * A(5)*B(5)
      X(2,3) = X(2,3) + z
      X(3,2) = X(3,2) + z
      z = fact * A(6)*B(6)
      X(1,3) = X(1,3) + z
      X(3,1) = X(3,1) + z

      X(4,4) = X(4,4) + factd4 * (A(1)*B(2) + A(4)*B(4) 
     *                          + A(4)*B(4) + A(2)*B(1))
      X(4,5) = X(4,5) + factd4 * (A(4)*B(5) + A(6)*B(2)
     *                          + A(2)*B(6) + A(5)*B(4))
      X(4,6) = X(4,6) + factd4 * (A(1)*B(5) + A(6)*B(4)
     *                          + A(4)*B(6) + A(5)*B(1))
      X(5,4) = X(5,4) + factd4 * (A(4)*B(5) + A(2)*B(6)
     *                          + A(6)*B(2) + A(5)*B(4))
      X(5,5) = X(5,5) + factd4 * (A(2)*B(3) + A(5)*B(5)
     *                          + A(5)*B(5) + A(3)*B(2))
      X(5,6) = X(5,6) + factd4 * (A(4)*B(3) + A(5)*B(6)
     *                          + A(6)*B(5) + A(3)*B(4))
      X(6,4) = X(6,4) + factd4 * (A(1)*B(5) + A(4)*B(6)
     *                          + A(6)*B(4) + A(5)*B(1))
      X(6,5) = X(6,5) + factd4 * (A(4)*B(3) + A(6)*B(5)
     *                          + A(5)*B(6) + A(3)*B(4))
      X(6,6) = X(6,6) + factd4 * (A(1)*B(3) + A(6)*B(6)
     *                          + A(6)*B(6) + A(3)*B(1))

      z = factd2 * (A(1)*B(4) + A(4)*B(1))
      X(4,1) = X(4,1) + z
      X(1,4) = X(1,4) + z
      z = factd2 * (A(4)*B(2) + A(2)*B(4))
      X(4,2) = X(4,2) + z
      X(2,4) = X(2,4) + z
      z = factd2 * (A(6)*B(5) + A(5)*B(6))
      X(4,3) = X(4,3) + z
      X(3,4) = X(3,4) + z
      z = factd2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(5,1) + z
      X(1,5) = X(1,5) + z
      z = factd2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(5,2) + z
      X(2,5) = X(2,5) + z
      z = factd2 * (A(5)*B(3) + A(3)*B(5))
      X(5,3) = X(5,3) + z
      X(3,5) = X(3,5) + z
      z = factd2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(6,1) + z
      X(1,6) = X(1,6) + z
      z = factd2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(6,2) + z
      X(2,6) = X(2,6) + z
      z = factd2 * (A(6)*B(3) + A(3)*B(6))
      X(6,3) = X(6,3) + z
      X(3,6) = X(3,6) + z

      return
      end
c
c=============================================================================
c   Differentiation
c                        D ATB      mit  A,B unsym., T,ATB sym. 2.Stufe
c                   X = -------       
c                         D T       -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_UTV_p(X,A,B)
      implicit none
      integer i, j
      double precision A(3,3), B(3,3), X(6,6), r1d2, r1d4
      r1d2 = 0.50d0
      r1d4 = 0.25d0

      do i=1, 3
        do j=1, 3
          X(i,j) = A(i,j)*B(j,i)
        enddo
      enddo

      X(4,4) = r1d4 * (A(1,1)*B(2,2) + A(1,2)*B(1,2)
     *               + A(2,1)*B(2,1) + A(2,2)*B(1,1))
      X(4,5) = r1d4 * (A(1,2)*B(3,2) + A(1,3)*B(2,2)
     *               + A(2,2)*B(3,1) + A(2,3)*B(2,1))
      X(4,6) = r1d4 * (A(1,1)*B(3,2) + A(1,3)*B(1,2)
     *               + A(2,1)*B(3,1) + A(2,3)*B(1,1))
      X(5,4) = r1d4 * (A(2,1)*B(2,3) + A(2,2)*B(1,3)
     *               + A(3,1)*B(2,2) + A(3,2)*B(1,2))
      X(5,5) = r1d4 * (A(2,2)*B(3,3) + A(2,3)*B(2,3)
     *               + A(3,2)*B(3,2) + A(3,3)*B(2,2))
      X(5,6) = r1d4 * (A(2,1)*B(3,3) + A(2,3)*B(1,3)
     *               + A(3,1)*B(3,2) + A(3,3)*B(1,2))
      X(6,4) = r1d4 * (A(1,1)*B(2,3) + A(1,2)*B(1,3)
     *               + A(3,1)*B(2,1) + A(3,2)*B(1,1))
      X(6,5) = r1d4 * (A(1,2)*B(3,3) + A(1,3)*B(2,3)
     *               + A(3,2)*B(3,1) + A(3,3)*B(2,1))
      X(6,6) = r1d4 * (A(1,1)*B(3,3) + A(1,3)*B(1,3)
     *               + A(3,1)*B(3,1) + A(3,3)*B(1,1))

      X(4,1) = r1d2 * (A(1,1)*B(1,2) + A(2,1)*B(1,1))
      X(1,4) = r1d2 * (A(1,1)*B(2,1) + A(1,2)*B(1,1))
      X(4,2) = r1d2 * (A(1,2)*B(2,2) + A(2,2)*B(2,1))
      X(2,4) = r1d2 * (A(2,1)*B(2,2) + A(2,2)*B(1,2))
      X(4,3) = r1d2 * (A(1,3)*B(3,2) + A(2,3)*B(3,1))
      X(3,4) = r1d2 * (A(3,1)*B(2,3) + A(3,2)*B(1,3))
      X(5,1) = r1d2 * (A(2,1)*B(1,3) + A(3,1)*B(1,2))
      X(1,5) = r1d2 * (A(1,2)*B(3,1) + A(1,3)*B(2,1))
      X(5,2) = r1d2 * (A(2,2)*B(2,3) + A(3,2)*B(2,2))
      X(2,5) = r1d2 * (A(2,2)*B(3,2) + A(2,3)*B(2,2))
      X(5,3) = r1d2 * (A(2,3)*B(3,3) + A(3,3)*B(3,2))
      X(3,5) = r1d2 * (A(3,2)*B(3,3) + A(3,3)*B(2,3))
      X(6,1) = r1d2 * (A(1,1)*B(1,3) + A(3,1)*B(1,1))
      X(1,6) = r1d2 * (A(1,1)*B(3,1) + A(1,3)*B(1,1))
      X(6,2) = r1d2 * (A(1,2)*B(2,3) + A(3,2)*B(2,1))
      X(2,6) = r1d2 * (A(2,1)*B(3,2) + A(2,3)*B(1,2))
      X(6,3) = r1d2 * (A(1,3)*B(3,3) + A(3,3)*B(3,1))
      X(3,6) = r1d2 * (A(3,1)*B(3,3) + A(3,3)*B(1,3))

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTV_m(X,A,B)
      implicit none
      integer i, j
      double precision A(3,3), B(3,3), X(6,6), r1d2, r1d4
      r1d2 = 0.50d0
      r1d4 = 0.25d0

      do i=1, 3
        do j=1, 3
          X(i,j) = - A(i,j)*B(j,i)
        enddo
      enddo

      X(4,4) = - r1d4 * (A(1,1)*B(2,2) + A(1,2)*B(1,2)
     *                 + A(2,1)*B(2,1) + A(2,2)*B(1,1))
      X(4,5) = - r1d4 * (A(1,2)*B(3,2) + A(1,3)*B(2,2)
     *                 + A(2,2)*B(3,1) + A(2,3)*B(2,1))
      X(4,6) = - r1d4 * (A(1,1)*B(3,2) + A(1,3)*B(1,2)
     *                 + A(2,1)*B(3,1) + A(2,3)*B(1,1))
      X(5,4) = - r1d4 * (A(2,1)*B(2,3) + A(2,2)*B(1,3)
     *                 + A(3,1)*B(2,2) + A(3,2)*B(1,2))
      X(5,5) = - r1d4 * (A(2,2)*B(3,3) + A(2,3)*B(2,3)
     *                 + A(3,2)*B(3,2) + A(3,3)*B(2,2))
      X(5,6) = - r1d4 * (A(2,1)*B(3,3) + A(2,3)*B(1,3)
     *                 + A(3,1)*B(3,2) + A(3,3)*B(1,2))
      X(6,4) = - r1d4 * (A(1,1)*B(2,3) + A(1,2)*B(1,3)
     *                 + A(3,1)*B(2,1) + A(3,2)*B(1,1))
      X(6,5) = - r1d4 * (A(1,2)*B(3,3) + A(1,3)*B(2,3)
     *                 + A(3,2)*B(3,1) + A(3,3)*B(2,1))
      X(6,6) = - r1d4 * (A(1,1)*B(3,3) + A(1,3)*B(1,3)
     *                 + A(3,1)*B(3,1) + A(3,3)*B(1,1))

      X(4,1) = - r1d2 * (A(1,1)*B(1,2) + A(2,1)*B(1,1))
      X(1,4) = - r1d2 * (A(1,1)*B(2,1) + A(1,2)*B(1,1))
      X(4,2) = - r1d2 * (A(1,2)*B(2,2) + A(2,2)*B(2,1))
      X(2,4) = - r1d2 * (A(2,1)*B(2,2) + A(2,2)*B(1,2))
      X(4,3) = - r1d2 * (A(1,3)*B(3,2) + A(2,3)*B(3,1))
      X(3,4) = - r1d2 * (A(3,1)*B(2,3) + A(3,2)*B(1,3))
      X(5,1) = - r1d2 * (A(2,1)*B(1,3) + A(3,1)*B(1,2))
      X(1,5) = - r1d2 * (A(1,2)*B(3,1) + A(1,3)*B(2,1))
      X(5,2) = - r1d2 * (A(2,2)*B(2,3) + A(3,2)*B(2,2))
      X(2,5) = - r1d2 * (A(2,2)*B(3,2) + A(2,3)*B(2,2))
      X(5,3) = - r1d2 * (A(2,3)*B(3,3) + A(3,3)*B(3,2))
      X(3,5) = - r1d2 * (A(3,2)*B(3,3) + A(3,3)*B(2,3))
      X(6,1) = - r1d2 * (A(1,1)*B(1,3) + A(3,1)*B(1,1))
      X(1,6) = - r1d2 * (A(1,1)*B(3,1) + A(1,3)*B(1,1))
      X(6,2) = - r1d2 * (A(1,2)*B(2,3) + A(3,2)*B(2,1))
      X(2,6) = - r1d2 * (A(2,1)*B(3,2) + A(2,3)*B(1,2))
      X(6,3) = - r1d2 * (A(1,3)*B(3,3) + A(3,3)*B(3,1))
      X(3,6) = - r1d2 * (A(3,1)*B(3,3) + A(3,3)*B(1,3))

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTV_f(X,A,B,fact)
      implicit none
      integer i, j
      double precision A(3,3), B(3,3), X(6,6), fact, factd2, factd4
      factd2 = fact / 2.d0
      factd4 = fact / 4.d0

      do i=1, 3
        do j=1, 3
          X(i,j) = fact * A(i,j)*B(j,i)
        enddo
      enddo

      X(4,4) = factd4 * (A(1,1)*B(2,2) + A(1,2)*B(1,2)
     *                 + A(2,1)*B(2,1) + A(2,2)*B(1,1))
      X(4,5) = factd4 * (A(1,2)*B(3,2) + A(1,3)*B(2,2)
     *                 + A(2,2)*B(3,1) + A(2,3)*B(2,1))
      X(4,6) = factd4 * (A(1,1)*B(3,2) + A(1,3)*B(1,2)
     *                 + A(2,1)*B(3,1) + A(2,3)*B(1,1))
      X(5,4) = factd4 * (A(2,1)*B(2,3) + A(2,2)*B(1,3)
     *                 + A(3,1)*B(2,2) + A(3,2)*B(1,2))
      X(5,5) = factd4 * (A(2,2)*B(3,3) + A(2,3)*B(2,3)
     *                 + A(3,2)*B(3,2) + A(3,3)*B(2,2))
      X(5,6) = factd4 * (A(2,1)*B(3,3) + A(2,3)*B(1,3)
     *                 + A(3,1)*B(3,2) + A(3,3)*B(1,2))
      X(6,4) = factd4 * (A(1,1)*B(2,3) + A(1,2)*B(1,3)
     *                 + A(3,1)*B(2,1) + A(3,2)*B(1,1))
      X(6,5) = factd4 * (A(1,2)*B(3,3) + A(1,3)*B(2,3)
     *                 + A(3,2)*B(3,1) + A(3,3)*B(2,1))
      X(6,6) = factd4 * (A(1,1)*B(3,3) + A(1,3)*B(1,3)
     *                 + A(3,1)*B(3,1) + A(3,3)*B(1,1))

      X(4,1) = factd2 * (A(1,1)*B(1,2) + A(2,1)*B(1,1))
      X(1,4) = factd2 * (A(1,1)*B(2,1) + A(1,2)*B(1,1))
      X(4,2) = factd2 * (A(1,2)*B(2,2) + A(2,2)*B(2,1))
      X(2,4) = factd2 * (A(2,1)*B(2,2) + A(2,2)*B(1,2))
      X(4,3) = factd2 * (A(1,3)*B(3,2) + A(2,3)*B(3,1))
      X(3,4) = factd2 * (A(3,1)*B(2,3) + A(3,2)*B(1,3))
      X(5,1) = factd2 * (A(2,1)*B(1,3) + A(3,1)*B(1,2))
      X(1,5) = factd2 * (A(1,2)*B(3,1) + A(1,3)*B(2,1))
      X(5,2) = factd2 * (A(2,2)*B(2,3) + A(3,2)*B(2,2))
      X(2,5) = factd2 * (A(2,2)*B(3,2) + A(2,3)*B(2,2))
      X(5,3) = factd2 * (A(2,3)*B(3,3) + A(3,3)*B(3,2))
      X(3,5) = factd2 * (A(3,2)*B(3,3) + A(3,3)*B(2,3))
      X(6,1) = factd2 * (A(1,1)*B(1,3) + A(3,1)*B(1,1))
      X(1,6) = factd2 * (A(1,1)*B(3,1) + A(1,3)*B(1,1))
      X(6,2) = factd2 * (A(1,2)*B(2,3) + A(3,2)*B(2,1))
      X(2,6) = factd2 * (A(2,1)*B(3,2) + A(2,3)*B(1,2))
      X(6,3) = factd2 * (A(1,3)*B(3,3) + A(3,3)*B(3,1))
      X(3,6) = factd2 * (A(3,1)*B(3,3) + A(3,3)*B(1,3))

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTV_ap(X,A,B)
      implicit none
      integer i, j
      double precision A(3,3), B(3,3), X(6,6), r1d2, r1d4
      r1d2 = 0.50d0
      r1d4 = 0.25d0

      do i=1, 3
        do j=1, 3
          X(i,j) = X(i,j) + A(i,j)*B(j,i)
        enddo
      enddo

      X(4,4) = X(4,4) + r1d4 * (A(1,1)*B(2,2) + A(1,2)*B(1,2)
     *                        + A(2,1)*B(2,1) + A(2,2)*B(1,1))
      X(4,5) = X(4,5) + r1d4 * (A(1,2)*B(3,2) + A(1,3)*B(2,2)
     *                        + A(2,2)*B(3,1) + A(2,3)*B(2,1))
      X(4,6) = X(4,6) + r1d4 * (A(1,1)*B(3,2) + A(1,3)*B(1,2)
     *                        + A(2,1)*B(3,1) + A(2,3)*B(1,1))
      X(5,4) = X(5,4) + r1d4 * (A(2,1)*B(2,3) + A(2,2)*B(1,3)
     *                        + A(3,1)*B(2,2) + A(3,2)*B(1,2))
      X(5,5) = X(5,5) + r1d4 * (A(2,2)*B(3,3) + A(2,3)*B(2,3)
     *                        + A(3,2)*B(3,2) + A(3,3)*B(2,2))
      X(5,6) = X(5,6) + r1d4 * (A(2,1)*B(3,3) + A(2,3)*B(1,3)
     *                        + A(3,1)*B(3,2) + A(3,3)*B(1,2))
      X(6,4) = X(6,4) + r1d4 * (A(1,1)*B(2,3) + A(1,2)*B(1,3)
     *                        + A(3,1)*B(2,1) + A(3,2)*B(1,1))
      X(6,5) = X(6,5) + r1d4 * (A(1,2)*B(3,3) + A(1,3)*B(2,3)
     *                        + A(3,2)*B(3,1) + A(3,3)*B(2,1))
      X(6,6) = X(6,6) + r1d4 * (A(1,1)*B(3,3) + A(1,3)*B(1,3)
     *                        + A(3,1)*B(3,1) + A(3,3)*B(1,1))

      X(4,1) = X(4,1) + r1d2 * (A(1,1)*B(1,2) + A(2,1)*B(1,1))
      X(1,4) = X(1,4) + r1d2 * (A(1,1)*B(2,1) + A(1,2)*B(1,1))
      X(4,2) = X(4,2) + r1d2 * (A(1,2)*B(2,2) + A(2,2)*B(2,1))
      X(2,4) = X(2,4) + r1d2 * (A(2,1)*B(2,2) + A(2,2)*B(1,2))
      X(4,3) = X(4,3) + r1d2 * (A(1,3)*B(3,2) + A(2,3)*B(3,1))
      X(3,4) = X(3,4) + r1d2 * (A(3,1)*B(2,3) + A(3,2)*B(1,3))
      X(5,1) = X(5,1) + r1d2 * (A(2,1)*B(1,3) + A(3,1)*B(1,2))
      X(1,5) = X(1,5) + r1d2 * (A(1,2)*B(3,1) + A(1,3)*B(2,1))
      X(5,2) = X(5,2) + r1d2 * (A(2,2)*B(2,3) + A(3,2)*B(2,2))
      X(2,5) = X(2,5) + r1d2 * (A(2,2)*B(3,2) + A(2,3)*B(2,2))
      X(5,3) = X(5,3) + r1d2 * (A(2,3)*B(3,3) + A(3,3)*B(3,2))
      X(3,5) = X(3,5) + r1d2 * (A(3,2)*B(3,3) + A(3,3)*B(2,3))
      X(6,1) = X(6,1) + r1d2 * (A(1,1)*B(1,3) + A(3,1)*B(1,1))
      X(1,6) = X(1,6) + r1d2 * (A(1,1)*B(3,1) + A(1,3)*B(1,1))
      X(6,2) = X(6,2) + r1d2 * (A(1,2)*B(2,3) + A(3,2)*B(2,1))
      X(2,6) = X(2,6) + r1d2 * (A(2,1)*B(3,2) + A(2,3)*B(1,2))
      X(6,3) = X(6,3) + r1d2 * (A(1,3)*B(3,3) + A(3,3)*B(3,1))
      X(3,6) = X(3,6) + r1d2 * (A(3,1)*B(3,3) + A(3,3)*B(1,3))

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTV_am(X,A,B)
      implicit none
      integer i, j
      double precision A(3,3), B(3,3), X(6,6), r1d2, r1d4
      r1d2 = 0.50d0
      r1d4 = 0.25d0

      do i=1, 3
        do j=1, 3
          X(i,j) = X(i,j) - A(i,j)*B(j,i)
        enddo
      enddo

      X(4,4) = X(4,4) - r1d4 * (A(1,1)*B(2,2) + A(1,2)*B(1,2)
     *                        + A(2,1)*B(2,1) + A(2,2)*B(1,1))
      X(4,5) = X(4,5) - r1d4 * (A(1,2)*B(3,2) + A(1,3)*B(2,2)
     *                        + A(2,2)*B(3,1) + A(2,3)*B(2,1))
      X(4,6) = X(4,6) - r1d4 * (A(1,1)*B(3,2) + A(1,3)*B(1,2)
     *                        + A(2,1)*B(3,1) + A(2,3)*B(1,1))
      X(5,4) = X(5,4) - r1d4 * (A(2,1)*B(2,3) + A(2,2)*B(1,3)
     *                        + A(3,1)*B(2,2) + A(3,2)*B(1,2))
      X(5,5) = X(5,5) - r1d4 * (A(2,2)*B(3,3) + A(2,3)*B(2,3)
     *                        + A(3,2)*B(3,2) + A(3,3)*B(2,2))
      X(5,6) = X(5,6) - r1d4 * (A(2,1)*B(3,3) + A(2,3)*B(1,3)
     *                        + A(3,1)*B(3,2) + A(3,3)*B(1,2))
      X(6,4) = X(6,4) - r1d4 * (A(1,1)*B(2,3) + A(1,2)*B(1,3)
     *                        + A(3,1)*B(2,1) + A(3,2)*B(1,1))
      X(6,5) = X(6,5) - r1d4 * (A(1,2)*B(3,3) + A(1,3)*B(2,3)
     *                        + A(3,2)*B(3,1) + A(3,3)*B(2,1))
      X(6,6) = X(6,6) - r1d4 * (A(1,1)*B(3,3) + A(1,3)*B(1,3)
     *                        + A(3,1)*B(3,1) + A(3,3)*B(1,1))

      X(4,1) = X(4,1) - r1d2 * (A(1,1)*B(1,2) + A(2,1)*B(1,1))
      X(1,4) = X(1,4) - r1d2 * (A(1,1)*B(2,1) + A(1,2)*B(1,1))
      X(4,2) = X(4,2) - r1d2 * (A(1,2)*B(2,2) + A(2,2)*B(2,1))
      X(2,4) = X(2,4) - r1d2 * (A(2,1)*B(2,2) + A(2,2)*B(1,2))
      X(4,3) = X(4,3) - r1d2 * (A(1,3)*B(3,2) + A(2,3)*B(3,1))
      X(3,4) = X(3,4) - r1d2 * (A(3,1)*B(2,3) + A(3,2)*B(1,3))
      X(5,1) = X(5,1) - r1d2 * (A(2,1)*B(1,3) + A(3,1)*B(1,2))
      X(1,5) = X(1,5) - r1d2 * (A(1,2)*B(3,1) + A(1,3)*B(2,1))
      X(5,2) = X(5,2) - r1d2 * (A(2,2)*B(2,3) + A(3,2)*B(2,2))
      X(2,5) = X(2,5) - r1d2 * (A(2,2)*B(3,2) + A(2,3)*B(2,2))
      X(5,3) = X(5,3) - r1d2 * (A(2,3)*B(3,3) + A(3,3)*B(3,2))
      X(3,5) = X(3,5) - r1d2 * (A(3,2)*B(3,3) + A(3,3)*B(2,3))
      X(6,1) = X(6,1) - r1d2 * (A(1,1)*B(1,3) + A(3,1)*B(1,1))
      X(1,6) = X(1,6) - r1d2 * (A(1,1)*B(3,1) + A(1,3)*B(1,1))
      X(6,2) = X(6,2) - r1d2 * (A(1,2)*B(2,3) + A(3,2)*B(2,1))
      X(2,6) = X(2,6) - r1d2 * (A(2,1)*B(3,2) + A(2,3)*B(1,2))
      X(6,3) = X(6,3) - r1d2 * (A(1,3)*B(3,3) + A(3,3)*B(3,1))
      X(3,6) = X(3,6) - r1d2 * (A(3,1)*B(3,3) + A(3,3)*B(1,3))

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTV_af(X,A,B,fact)
      implicit none
      integer i, j
      double precision A(3,3), B(3,3), X(6,6), fact, factd2, factd4
      factd2 = fact / 2.d0
      factd4 = fact / 4.d0

      do i=1, 3
        do j=1, 3
          X(i,j) = X(i,j) + fact * A(i,j)*B(j,i)
        enddo
      enddo

      X(4,4) = X(4,4) + factd4 * (A(1,1)*B(2,2) + A(1,2)*B(1,2)
     *                          + A(2,1)*B(2,1) + A(2,2)*B(1,1))
      X(4,5) = X(4,5) + factd4 * (A(1,2)*B(3,2) + A(1,3)*B(2,2)
     *                          + A(2,2)*B(3,1) + A(2,3)*B(2,1))
      X(4,6) = X(4,6) + factd4 * (A(1,1)*B(3,2) + A(1,3)*B(1,2)
     *                          + A(2,1)*B(3,1) + A(2,3)*B(1,1))
      X(5,4) = X(5,4) + factd4 * (A(2,1)*B(2,3) + A(2,2)*B(1,3)
     *                          + A(3,1)*B(2,2) + A(3,2)*B(1,2))
      X(5,5) = X(5,5) + factd4 * (A(2,2)*B(3,3) + A(2,3)*B(2,3)
     *                          + A(3,2)*B(3,2) + A(3,3)*B(2,2))
      X(5,6) = X(5,6) + factd4 * (A(2,1)*B(3,3) + A(2,3)*B(1,3)
     *                          + A(3,1)*B(3,2) + A(3,3)*B(1,2))
      X(6,4) = X(6,4) + factd4 * (A(1,1)*B(2,3) + A(1,2)*B(1,3)
     *                          + A(3,1)*B(2,1) + A(3,2)*B(1,1))
      X(6,5) = X(6,5) + factd4 * (A(1,2)*B(3,3) + A(1,3)*B(2,3)
     *                          + A(3,2)*B(3,1) + A(3,3)*B(2,1))
      X(6,6) = X(6,6) + factd4 * (A(1,1)*B(3,3) + A(1,3)*B(1,3)
     *                          + A(3,1)*B(3,1) + A(3,3)*B(1,1))

      X(4,1) = X(4,1) + factd2 * (A(1,1)*B(1,2) + A(2,1)*B(1,1))
      X(1,4) = X(1,4) + factd2 * (A(1,1)*B(2,1) + A(1,2)*B(1,1))
      X(4,2) = X(4,2) + factd2 * (A(1,2)*B(2,2) + A(2,2)*B(2,1))
      X(2,4) = X(2,4) + factd2 * (A(2,1)*B(2,2) + A(2,2)*B(1,2))
      X(4,3) = X(4,3) + factd2 * (A(1,3)*B(3,2) + A(2,3)*B(3,1))
      X(3,4) = X(3,4) + factd2 * (A(3,1)*B(2,3) + A(3,2)*B(1,3))
      X(5,1) = X(5,1) + factd2 * (A(2,1)*B(1,3) + A(3,1)*B(1,2))
      X(1,5) = X(1,5) + factd2 * (A(1,2)*B(3,1) + A(1,3)*B(2,1))
      X(5,2) = X(5,2) + factd2 * (A(2,2)*B(2,3) + A(3,2)*B(2,2))
      X(2,5) = X(2,5) + factd2 * (A(2,2)*B(3,2) + A(2,3)*B(2,2))
      X(5,3) = X(5,3) + factd2 * (A(2,3)*B(3,3) + A(3,3)*B(3,2))
      X(3,5) = X(3,5) + factd2 * (A(3,2)*B(3,3) + A(3,3)*B(2,3))
      X(6,1) = X(6,1) + factd2 * (A(1,1)*B(1,3) + A(3,1)*B(1,1))
      X(1,6) = X(1,6) + factd2 * (A(1,1)*B(3,1) + A(1,3)*B(1,1))
      X(6,2) = X(6,2) + factd2 * (A(1,2)*B(2,3) + A(3,2)*B(2,1))
      X(2,6) = X(2,6) + factd2 * (A(2,1)*B(3,2) + A(2,3)*B(1,2))
      X(6,3) = X(6,3) + factd2 * (A(1,3)*B(3,3) + A(3,3)*B(3,1))
      X(3,6) = X(3,6) + factd2 * (A(3,1)*B(3,3) + A(3,3)*B(1,3))

      return
      end
c
c=============================================================================
c  Differentiation
c                    D sym(AT)       mit  A, T sym. 2.Stufe
c               X = -----------      
c                      D T           -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_symAT_p(X,A)
      implicit none
      integer          i, j
      double precision A(6), X(6,6), r0, r1d2, r1d4
      r0   = 0.d0
      r1d2 = 0.5d0
      r1d4 = 0.25d0

      X(1,1) = A(1)
      X(2,2) = A(2)
      X(3,3) = A(3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,3) = r0

      X(1,4) = r1d2 * A(4)
      X(1,5) = r0
      X(1,6) = r1d2 * A(6)
      X(2,4) = r1d2 * A(4)
      X(2,5) = r1d2 * A(5)
      X(2,6) = r0
      X(3,4) = r0
      X(3,5) = r1d2 * A(5)
      X(3,6) = r1d2 * A(6)

      X(4,4) = r1d4 * (A(1) + A(2))
      X(5,5) = r1d4 * (A(2) + A(3))
      X(6,6) = r1d4 * (A(3) + A(1))

      X(4,5) = r1d4 * A(6)
      X(4,6) = r1d4 * A(5)
      X(5,6) = r1d4 * A(4)

      do i=1, 6
        do j=i+1, 6
          X(j,i) = X(i,j)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symAT_m(X,A)
      implicit none
      integer          i, j
      double precision A(6), X(6,6), r0, r1d2, r1d4
      r0   = 0.d0
      r1d2 = 0.5d0
      r1d4 = 0.25d0

      X(1,1) = - A(1)
      X(2,2) = - A(2)
      X(3,3) = - A(3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,3) = r0

      X(1,4) = - r1d2 * A(4)
      X(1,5) = r0
      X(1,6) = - r1d2 * A(6)
      X(2,4) = - r1d2 * A(4)
      X(2,5) = - r1d2 * A(5)
      X(2,6) = r0
      X(3,4) = r0
      X(3,5) = - r1d2 * A(5)
      X(3,6) = - r1d2 * A(6)

      X(4,4) = - r1d4 * (A(1) + A(2))
      X(5,5) = - r1d4 * (A(2) + A(3))
      X(6,6) = - r1d4 * (A(3) + A(1))

      X(4,5) = - r1d4 * A(6)
      X(4,6) = - r1d4 * A(5)
      X(5,6) = - r1d4 * A(4)

      do i=1, 6
        do j=i+1, 6
          X(j,i) = X(i,j)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symAT_f(X,A,fact)
      implicit none
      integer          i, j
      double precision A(6), X(6,6), r0, fact, factd2, factd4
      r0   = 0.d0
      factd2 = fact * 0.50d0
      factd4 = fact * 0.25d0

      X(1,1) = fact * A(1)
      X(2,2) = fact * A(2)
      X(3,3) = fact * A(3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,3) = r0

      X(1,4) = factd2 * A(4)
      X(1,5) = r0
      X(1,6) = factd2 * A(6)
      X(2,4) = factd2 * A(4)
      X(2,5) = factd2 * A(5)
      X(2,6) = r0
      X(3,4) = r0
      X(3,5) = factd2 * A(5)
      X(3,6) = factd2 * A(6)

      X(4,4) = factd4 * (A(1) + A(2))
      X(5,5) = factd4 * (A(2) + A(3))
      X(6,6) = factd4 * (A(3) + A(1))

      X(4,5) = factd4 * A(6)
      X(4,6) = factd4 * A(5)
      X(5,6) = factd4 * A(4)

      do i=1, 6
        do j=i+1, 6
          X(j,i) = X(i,j)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symAT_ap(X,A)
      implicit none
      double precision A(6), X(6,6), r1d2, r1d4
      r1d2 = 0.5d0
      r1d4 = 0.25d0

      X(1,1) = X(1,1) + A(1)
      X(2,2) = X(2,2) + A(2)
      X(3,3) = X(3,3) + A(3)

      X(1,4) = X(1,4) + r1d2 * A(4)
      X(4,1) = X(4,1) + r1d2 * A(4)
      X(1,6) = X(1,6) + r1d2 * A(6)
      X(6,1) = X(6,1) + r1d2 * A(6)
      X(2,4) = X(2,4) + r1d2 * A(4)
      X(4,2) = X(4,2) + r1d2 * A(4)
      X(2,5) = X(2,5) + r1d2 * A(5)
      X(5,2) = X(5,2) + r1d2 * A(5)
      X(3,5) = X(3,5) + r1d2 * A(5)
      X(5,3) = X(5,3) + r1d2 * A(5)
      X(3,6) = X(3,6) + r1d2 * A(6)
      X(6,3) = X(6,3) + r1d2 * A(6)

      X(4,4) = X(4,4) + r1d4 * (A(1) + A(2))
      X(5,5) = X(5,5) + r1d4 * (A(2) + A(3))
      X(6,6) = X(6,6) + r1d4 * (A(3) + A(1))

      X(4,5) = X(4,5) + r1d4 * A(6)
      X(5,4) = X(5,4) + r1d4 * A(6)
      X(4,6) = X(4,6) + r1d4 * A(5)
      X(6,4) = X(6,4) + r1d4 * A(5)
      X(5,6) = X(5,6) + r1d4 * A(4)
      X(6,5) = X(6,5) + r1d4 * A(4)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symAT_am(X,A)
      implicit none
      double precision A(6), X(6,6), r1d2, r1d4
      r1d2 = 0.5d0
      r1d4 = 0.25d0

      X(1,1) = X(1,1) - A(1)
      X(2,2) = X(2,2) - A(2)
      X(3,3) = X(3,3) - A(3)

      X(1,4) = X(1,4) - r1d2 * A(4)
      X(4,1) = X(4,1) - r1d2 * A(4)
      X(1,6) = X(1,6) - r1d2 * A(6)
      X(6,1) = X(6,1) - r1d2 * A(6)
      X(2,4) = X(2,4) - r1d2 * A(4)
      X(4,2) = X(4,2) - r1d2 * A(4)
      X(2,5) = X(2,5) - r1d2 * A(5)
      X(5,2) = X(5,2) - r1d2 * A(5)
      X(3,5) = X(3,5) - r1d2 * A(5)
      X(5,3) = X(5,3) - r1d2 * A(5)
      X(3,6) = X(3,6) - r1d2 * A(6)
      X(6,3) = X(6,3) - r1d2 * A(6)

      X(4,4) = X(4,4) - r1d4 * (A(1) + A(2))
      X(5,5) = X(5,5) - r1d4 * (A(2) + A(3))
      X(6,6) = X(6,6) - r1d4 * (A(3) + A(1))

      X(4,5) = X(4,5) - r1d4 * A(6)
      X(5,4) = X(5,4) - r1d4 * A(6)
      X(4,6) = X(4,6) - r1d4 * A(5)
      X(6,4) = X(6,4) - r1d4 * A(5)
      X(5,6) = X(5,6) - r1d4 * A(4)
      X(6,5) = X(6,5) - r1d4 * A(4)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symAT_af(X,A,fact)
      implicit none
      double precision A(6), X(6,6), fact, factd2, factd4
      factd2 = fact * 0.5d0
      factd4 = fact * 0.25d0

      X(1,1) = X(1,1) + fact * A(1)
      X(2,2) = X(2,2) + fact * A(2)
      X(3,3) = X(3,3) + fact * A(3)

      X(1,4) = X(1,4) + factd2 * A(4)
      X(4,1) = X(4,1) + factd2 * A(4)
      X(1,6) = X(1,6) + factd2 * A(6)
      X(6,1) = X(6,1) + factd2 * A(6)
      X(2,4) = X(2,4) + factd2 * A(4)
      X(4,2) = X(4,2) + factd2 * A(4)
      X(2,5) = X(2,5) + factd2 * A(5)
      X(5,2) = X(5,2) + factd2 * A(5)
      X(3,5) = X(3,5) + factd2 * A(5)
      X(5,3) = X(5,3) + factd2 * A(5)
      X(3,6) = X(3,6) + factd2 * A(6)
      X(6,3) = X(6,3) + factd2 * A(6)

      X(4,4) = X(4,4) + factd4 * (A(1) + A(2))
      X(5,5) = X(5,5) + factd4 * (A(2) + A(3))
      X(6,6) = X(6,6) + factd4 * (A(3) + A(1))

      X(4,5) = X(4,5) + factd4 * A(6)
      X(5,4) = X(5,4) + factd4 * A(6)
      X(4,6) = X(4,6) + factd4 * A(5)
      X(6,4) = X(6,4) + factd4 * A(5)
      X(5,6) = X(5,6) + factd4 * A(4)
      X(6,5) = X(6,5) + factd4 * A(4)

      return
      end
c
c=============================================================================
c  Differentiation
c                    D sym(AT)       mit  T sym., A unsym. 2.Stufe
c               X = -----------      
c                      D T           -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_symUT_p(X,A)
      implicit none
      double precision A(3,3), X(6,6), r0, r1d2, r1d4
      r0   = 0.d0
      r1d2 = 0.5d0
      r1d4 = 0.25d0

      X(1,1) = A(1,1)
      X(2,2) = A(2,2)
      X(3,3) = A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,4) = r1d2 * A(1,2)
      X(4,1) = r1d2 * A(2,1)
      X(1,5) = r0
      X(5,1) = r0
      X(1,6) = r1d2 * A(1,3)
      X(6,1) = r1d2 * A(3,1)
      X(2,4) = r1d2 * A(2,1)
      X(4,2) = r1d2 * A(1,2)
      X(2,5) = r1d2 * A(2,3)
      X(5,2) = r1d2 * A(3,2)
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(3,5) = r1d2 * A(3,2)
      X(5,3) = r1d2 * A(2,3)
      X(3,6) = r1d2 * A(3,1)
      X(6,3) = r1d2 * A(1,3)

      X(4,4) = r1d4 * (A(1,1) + A(2,2))
      X(5,5) = r1d4 * (A(2,2) + A(3,3))
      X(6,6) = r1d4 * (A(3,3) + A(1,1))

      X(4,5) = r1d4 * A(1,3)
      X(5,4) = r1d4 * A(3,1)
      X(4,6) = r1d4 * A(2,3)
      X(6,4) = r1d4 * A(3,2)
      X(5,6) = r1d4 * A(2,1)
      X(6,5) = r1d4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symUT_m(X,A)
      implicit none
      double precision A(3,3), X(6,6), r0, r1d2, r1d4
      r0   = 0.d0
      r1d2 = 0.5d0
      r1d4 = 0.25d0

      X(1,1) = - A(1,1)
      X(2,2) = - A(2,2)
      X(3,3) = - A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,4) = - r1d2 * A(1,2)
      X(4,1) = - r1d2 * A(2,1)
      X(1,5) = r0
      X(5,1) = r0
      X(1,6) = - r1d2 * A(1,3)
      X(6,1) = - r1d2 * A(3,1)
      X(2,4) = - r1d2 * A(2,1)
      X(4,2) = - r1d2 * A(1,2)
      X(2,5) = - r1d2 * A(2,3)
      X(5,2) = - r1d2 * A(3,2)
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(3,5) = - r1d2 * A(3,2)
      X(5,3) = - r1d2 * A(2,3)
      X(3,6) = - r1d2 * A(3,1)
      X(6,3) = - r1d2 * A(1,3)

      X(4,4) = - r1d4 * (A(1,1) + A(2,2))
      X(5,5) = - r1d4 * (A(2,2) + A(3,3))
      X(6,6) = - r1d4 * (A(3,3) + A(1,1))

      X(4,5) = - r1d4 * A(1,3)
      X(5,4) = - r1d4 * A(3,1)
      X(4,6) = - r1d4 * A(2,3)
      X(6,4) = - r1d4 * A(3,2)
      X(5,6) = - r1d4 * A(2,1)
      X(6,5) = - r1d4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symUT_f(X,A,fact)
      implicit none
      double precision A(3,3), X(6,6), r0, fact, factd2, factd4
      r0   = 0.d0
      factd2 = fact * 0.50d0
      factd4 = fact * 0.25d0

      X(1,1) = fact * A(1,1)
      X(2,2) = fact * A(2,2)
      X(3,3) = fact * A(3,3)

      X(1,2) = r0
      X(1,3) = r0
      X(2,1) = r0
      X(2,3) = r0
      X(3,1) = r0
      X(3,2) = r0

      X(1,4) = factd2 * A(1,2)
      X(4,1) = factd2 * A(2,1)
      X(1,5) = r0
      X(5,1) = r0
      X(1,6) = factd2 * A(1,3)
      X(6,1) = factd2 * A(3,1)
      X(2,4) = factd2 * A(2,1)
      X(4,2) = factd2 * A(1,2)
      X(2,5) = factd2 * A(2,3)
      X(5,2) = factd2 * A(3,2)
      X(2,6) = r0
      X(6,2) = r0
      X(3,4) = r0
      X(4,3) = r0
      X(3,5) = factd2 * A(3,2)
      X(5,3) = factd2 * A(2,3)
      X(3,6) = factd2 * A(3,1)
      X(6,3) = factd2 * A(1,3)

      X(4,4) = factd4 * (A(1,1) + A(2,2))
      X(5,5) = factd4 * (A(2,2) + A(3,3))
      X(6,6) = factd4 * (A(3,3) + A(1,1))

      X(4,5) = factd4 * A(1,3)
      X(5,4) = factd4 * A(3,1)
      X(4,6) = factd4 * A(2,3)
      X(6,4) = factd4 * A(3,2)
      X(5,6) = factd4 * A(2,1)
      X(6,5) = factd4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symUT_ap(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d2, r1d4
      r1d2 = 0.5d0
      r1d4 = 0.25d0

      X(1,1) = X(1,1) + A(1,1)
      X(2,2) = X(2,2) + A(2,2)
      X(3,3) = X(3,3) + A(3,3)

      X(1,4) = X(1,4) + r1d2 * A(1,2)
      X(4,1) = X(4,1) + r1d2 * A(2,1)
      X(1,6) = X(1,6) + r1d2 * A(1,3)
      X(6,1) = X(6,1) + r1d2 * A(3,1)
      X(2,4) = X(2,4) + r1d2 * A(2,1)
      X(4,2) = X(4,2) + r1d2 * A(1,2)
      X(2,5) = X(2,5) + r1d2 * A(2,3)
      X(5,2) = X(5,2) + r1d2 * A(3,2)
      X(3,5) = X(3,5) + r1d2 * A(3,2)
      X(5,3) = X(5,3) + r1d2 * A(2,3)
      X(3,6) = X(3,6) + r1d2 * A(3,1)
      X(6,3) = X(6,3) + r1d2 * A(1,3)

      X(4,4) = X(4,4) + r1d4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) + r1d4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) + r1d4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) + r1d4 * A(1,3)
      X(5,4) = X(5,4) + r1d4 * A(3,1)
      X(4,6) = X(4,6) + r1d4 * A(2,3)
      X(6,4) = X(6,4) + r1d4 * A(3,2)
      X(5,6) = X(5,6) + r1d4 * A(2,1)
      X(6,5) = X(6,5) + r1d4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symUT_am(X,A)
      implicit none
      double precision A(3,3), X(6,6), r1d2, r1d4
      r1d2 = 0.5d0
      r1d4 = 0.25d0

      X(1,1) = X(1,1) - A(1,1)
      X(2,2) = X(2,2) - A(2,2)
      X(3,3) = X(3,3) - A(3,3)

      X(1,4) = X(1,4) - r1d2 * A(1,2)
      X(4,1) = X(4,1) - r1d2 * A(2,1)
      X(1,6) = X(1,6) - r1d2 * A(1,3)
      X(6,1) = X(6,1) - r1d2 * A(3,1)
      X(2,4) = X(2,4) - r1d2 * A(2,1)
      X(4,2) = X(4,2) - r1d2 * A(1,2)
      X(2,5) = X(2,5) - r1d2 * A(2,3)
      X(5,2) = X(5,2) - r1d2 * A(3,2)
      X(3,5) = X(3,5) - r1d2 * A(3,2)
      X(5,3) = X(5,3) - r1d2 * A(2,3)
      X(3,6) = X(3,6) - r1d2 * A(3,1)
      X(6,3) = X(6,3) - r1d2 * A(1,3)

      X(4,4) = X(4,4) - r1d4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) - r1d4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) - r1d4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) - r1d4 * A(1,3)
      X(5,4) = X(5,4) - r1d4 * A(3,1)
      X(4,6) = X(4,6) - r1d4 * A(2,3)
      X(6,4) = X(6,4) - r1d4 * A(3,2)
      X(5,6) = X(5,6) - r1d4 * A(2,1)
      X(6,5) = X(6,5) - r1d4 * A(1,2)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symUT_af(X,A,fact)
      implicit none
      double precision A(3,3), X(6,6), fact, factd2, factd4
      factd2 = fact * 0.5d0
      factd4 = fact * 0.25d0

      X(1,1) = X(1,1) + fact * A(1,1)
      X(2,2) = X(2,2) + fact * A(2,2)
      X(3,3) = X(3,3) + fact * A(3,3)

      X(1,4) = X(1,4) + factd2 * A(1,2)
      X(4,1) = X(4,1) + factd2 * A(2,1)
      X(1,6) = X(1,6) + factd2 * A(1,3)
      X(6,1) = X(6,1) + factd2 * A(3,1)
      X(2,4) = X(2,4) + factd2 * A(2,1)
      X(4,2) = X(4,2) + factd2 * A(1,2)
      X(2,5) = X(2,5) + factd2 * A(2,3)
      X(5,2) = X(5,2) + factd2 * A(3,2)
      X(3,5) = X(3,5) + factd2 * A(3,2)
      X(5,3) = X(5,3) + factd2 * A(2,3)
      X(3,6) = X(3,6) + factd2 * A(3,1)
      X(6,3) = X(6,3) + factd2 * A(1,3)

      X(4,4) = X(4,4) + factd4 * (A(1,1) + A(2,2))
      X(5,5) = X(5,5) + factd4 * (A(2,2) + A(3,3))
      X(6,6) = X(6,6) + factd4 * (A(3,3) + A(1,1))

      X(4,5) = X(4,5) + factd4 * A(1,3)
      X(5,4) = X(5,4) + factd4 * A(3,1)
      X(4,6) = X(4,6) + factd4 * A(2,3)
      X(6,4) = X(6,4) + factd4 * A(3,2)
      X(5,6) = X(5,6) + factd4 * A(2,1)
      X(6,5) = X(6,5) + factd4 * A(1,2)

      return
      end
c
c=============================================================================
c  Differentiation
c                     D sym(ATB)       mit A,T,B sym. 2.Stufe
c                X = ------------      
c                        D T           -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_symATB_p(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d4, r1d2, r2
      r1d4 = 0.25d0
      r1d2 = 0.5d0
      r2 = 2.d0

      X(1,1) = A(1)*B(1)
      X(2,2) = A(2)*B(2)
      X(3,3) = A(3)*B(3)

      X(1,2) = A(4)*B(4)
      X(1,3) = A(6)*B(6)
      X(2,1) = X(1,2) 
      X(2,3) = A(5)*B(5)
      X(3,1) = X(1,3)
      X(3,2) = X(2,3)

      X(1,4) = r1d2 * (A(4)*B(1) + A(1)*B(4))
      X(4,1) = X(1,4)
      X(1,5) = r1d2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(1,5)
      X(1,6) = r1d2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(1,6)
      X(2,4) = r1d2 * (A(2)*B(4) + A(4)*B(2))
      X(4,2) = X(2,4)
      X(2,5) = r1d2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(2,5)
      X(2,6) = r1d2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(2,6)
      X(3,4) = r1d2 * (A(5)*B(6) + A(6)*B(5))
      X(4,3) = X(3,4)
      X(3,5) = r1d2 * (A(3)*B(5) + A(5)*B(3))
      X(5,3) = X(3,5)
      X(3,6) = r1d2 * (A(3)*B(6) + A(6)*B(3))
      X(6,3) = X(3,6)

      X(4,4) = r1d4 * (A(1)*B(2) + r2 * A(4)*B(4) + A(2)*B(1))
      X(5,5) = r1d4 * (A(2)*B(3) + r2 * A(5)*B(5) + A(3)*B(2))
      X(6,6) = r1d4 * (A(3)*B(1) + r2 * A(6)*B(6) + A(1)*B(3))

      X(4,5) = r1d4 * (A(4)*B(5) + A(5)*B(4) + A(2)*B(6) + A(6)*B(2))
      X(5,4) = X(4,5)
      X(4,6) = r1d4 * (A(4)*B(6) + A(6)*B(4) + A(1)*B(5) + A(5)*B(1))
      X(6,4) = X(4,6)
      X(5,6) = r1d4 * (A(5)*B(6) + A(6)*B(5) + A(3)*B(4) + A(4)*B(3))
      X(6,5) = X(5,6)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symATB_m(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d4, r1d2, r2
      r1d4 = 0.25d0
      r1d2 = 0.5d0
      r2 = 2.d0

      X(1,1) = - A(1)*B(1)
      X(2,2) = - A(2)*B(2)
      X(3,3) = - A(3)*B(3)

      X(1,2) = - A(4)*B(4)
      X(1,3) = - A(6)*B(6)
      X(2,1) = X(1,2) 
      X(2,3) = - A(5)*B(5)
      X(3,1) = X(1,3)
      X(3,2) = X(2,3)

      X(1,4) = - r1d2 * (A(4)*B(1) + A(1)*B(4))
      X(4,1) = X(1,4)
      X(1,5) = - r1d2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(1,5)
      X(1,6) = - r1d2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(1,6)
      X(2,4) = - r1d2 * (A(2)*B(4) + A(4)*B(2))
      X(4,2) = X(2,4)
      X(2,5) = - r1d2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(2,5)
      X(2,6) = - r1d2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(2,6)
      X(3,4) = - r1d2 * (A(5)*B(6) + A(6)*B(5))
      X(4,3) = X(3,4)
      X(3,5) = - r1d2 * (A(3)*B(5) + A(5)*B(3))
      X(5,3) = X(3,5)
      X(3,6) = - r1d2 * (A(3)*B(6) + A(6)*B(3))
      X(6,3) = X(3,6)

      X(4,4) = - r1d4 * (A(1)*B(2) + r2 * A(4)*B(4) + A(2)*B(1))
      X(5,5) = - r1d4 * (A(2)*B(3) + r2 * A(5)*B(5) + A(3)*B(2))
      X(6,6) = - r1d4 * (A(3)*B(1) + r2 * A(6)*B(6) + A(1)*B(3))

      X(4,5) = - r1d4 * (A(4)*B(5) + A(5)*B(4) + A(2)*B(6) + A(6)*B(2))
      X(5,4) = X(4,5)
      X(4,6) = - r1d4 * (A(4)*B(6) + A(6)*B(4) + A(1)*B(5) + A(5)*B(1))
      X(6,4) = X(4,6)
      X(5,6) = - r1d4 * (A(5)*B(6) + A(6)*B(5) + A(3)*B(4) + A(4)*B(3))
      X(6,5) = X(5,6)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symATB_f(X,A,B,fact)
      implicit none
      double precision A(6), B(6), X(6,6), fact, factd2, factd4, r2
      r2     = 2.d0
      factd2 = fact / r2
      factd4 = fact / 4.d0

      X(1,1) = fact * A(1)*B(1)
      X(2,2) = fact * A(2)*B(2)
      X(3,3) = fact * A(3)*B(3)

      X(1,2) = fact * A(4)*B(4)
      X(1,3) = fact * A(6)*B(6)
      X(2,1) = X(1,2) 
      X(2,3) = fact * A(5)*B(5)
      X(3,1) = X(1,3)
      X(3,2) = X(2,3)

      X(1,4) = factd2 * (A(4)*B(1) + A(1)*B(4))
      X(4,1) = X(1,4)
      X(1,5) = factd2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(1,5)
      X(1,6) = factd2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(1,6)
      X(2,4) = factd2 * (A(2)*B(4) + A(4)*B(2))
      X(4,2) = X(2,4)
      X(2,5) = factd2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(2,5)
      X(2,6) = factd2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(2,6)
      X(3,4) = factd2 * (A(5)*B(6) + A(6)*B(5))
      X(4,3) = X(3,4)
      X(3,5) = factd2 * (A(3)*B(5) + A(5)*B(3))
      X(5,3) = X(3,5)
      X(3,6) = factd2 * (A(3)*B(6) + A(6)*B(3))
      X(6,3) = X(3,6)

      X(4,4) = factd4 * (A(1)*B(2) + r2 * A(4)*B(4) + A(2)*B(1))
      X(5,5) = factd4 * (A(2)*B(3) + r2 * A(5)*B(5) + A(3)*B(2))
      X(6,6) = factd4 * (A(3)*B(1) + r2 * A(6)*B(6) + A(1)*B(3))

      X(4,5) = factd4 * (A(4)*B(5) + A(5)*B(4) + A(2)*B(6) + A(6)*B(2))
      X(5,4) = X(4,5)
      X(4,6) = factd4 * (A(4)*B(6) + A(6)*B(4) + A(1)*B(5) + A(5)*B(1))
      X(6,4) = X(4,6)
      X(5,6) = factd4 * (A(5)*B(6) + A(6)*B(5) + A(3)*B(4) + A(4)*B(3))
      X(6,5) = X(5,6)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symATB_ap(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d4, r1d2, r2
      r1d4 = 0.25d0
      r1d2 = 0.5d0
      r2 = 2.d0

      X(1,1) = A(1)*B(1)
      X(2,2) = A(2)*B(2)
      X(3,3) = A(3)*B(3)

      X(1,2) = A(4)*B(4)
      X(1,3) = A(6)*B(6)
      X(2,1) = X(1,2) 
      X(2,3) = A(5)*B(5)
      X(3,1) = X(1,3)
      X(3,2) = X(2,3)

      X(1,4) = r1d2 * (A(4)*B(1) + A(1)*B(4))
      X(4,1) = X(1,4)
      X(1,5) = r1d2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(1,5)
      X(1,6) = r1d2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(1,6)
      X(2,4) = r1d2 * (A(2)*B(4) + A(4)*B(2))
      X(4,2) = X(2,4)
      X(2,5) = r1d2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(2,5)
      X(2,6) = r1d2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(2,6)
      X(3,4) = r1d2 * (A(5)*B(6) + A(6)*B(5))
      X(4,3) = X(3,4)
      X(3,5) = r1d2 * (A(3)*B(5) + A(5)*B(3))
      X(5,3) = X(3,5)
      X(3,6) = r1d2 * (A(3)*B(6) + A(6)*B(3))
      X(6,3) = X(3,6)

      X(4,4) = r1d4 * (A(1)*B(2) + r2 * A(4)*B(4) + A(2)*B(1))
      X(5,5) = r1d4 * (A(2)*B(3) + r2 * A(5)*B(5) + A(3)*B(2))
      X(6,6) = r1d4 * (A(3)*B(1) + r2 * A(6)*B(6) + A(1)*B(3))

      X(4,5) = r1d4 * (A(4)*B(5) + A(5)*B(4) + A(2)*B(6) + A(6)*B(2))
      X(5,4) = X(4,5)
      X(4,6) = r1d4 * (A(4)*B(6) + A(6)*B(4) + A(1)*B(5) + A(5)*B(1))
      X(6,4) = X(4,6)
      X(5,6) = r1d4 * (A(5)*B(6) + A(6)*B(5) + A(3)*B(4) + A(4)*B(3))
      X(6,5) = X(5,6)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symATB_am(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d4, r1d2, r2
      r1d4 = 0.25d0
      r1d2 = 0.5d0
      r2 = 2.d0

      X(1,1) = - A(1)*B(1)
      X(2,2) = - A(2)*B(2)
      X(3,3) = - A(3)*B(3)

      X(1,2) = - A(4)*B(4)
      X(1,3) = - A(6)*B(6)
      X(2,1) = X(1,2) 
      X(2,3) = - A(5)*B(5)
      X(3,1) = X(1,3)
      X(3,2) = X(2,3)

      X(1,4) = - r1d2 * (A(4)*B(1) + A(1)*B(4))
      X(4,1) = X(1,4)
      X(1,5) = - r1d2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(1,5)
      X(1,6) = - r1d2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(1,6)
      X(2,4) = - r1d2 * (A(2)*B(4) + A(4)*B(2))
      X(4,2) = X(2,4)
      X(2,5) = - r1d2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(2,5)
      X(2,6) = - r1d2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(2,6)
      X(3,4) = - r1d2 * (A(5)*B(6) + A(6)*B(5))
      X(4,3) = X(3,4)
      X(3,5) = - r1d2 * (A(3)*B(5) + A(5)*B(3))
      X(5,3) = X(3,5)
      X(3,6) = - r1d2 * (A(3)*B(6) + A(6)*B(3))
      X(6,3) = X(3,6)

      X(4,4) = - r1d4 * (A(1)*B(2) + r2 * A(4)*B(4) + A(2)*B(1))
      X(5,5) = - r1d4 * (A(2)*B(3) + r2 * A(5)*B(5) + A(3)*B(2))
      X(6,6) = - r1d4 * (A(3)*B(1) + r2 * A(6)*B(6) + A(1)*B(3))

      X(4,5) = - r1d4 * (A(4)*B(5) + A(5)*B(4) + A(2)*B(6) + A(6)*B(2))
      X(5,4) = X(4,5)
      X(4,6) = - r1d4 * (A(4)*B(6) + A(6)*B(4) + A(1)*B(5) + A(5)*B(1))
      X(6,4) = X(4,6)
      X(5,6) = - r1d4 * (A(5)*B(6) + A(6)*B(5) + A(3)*B(4) + A(4)*B(3))
      X(6,5) = X(5,6)

      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_symATB_af(X,A,B,fact)
      implicit none
      double precision A(6), B(6), X(6,6), fact, factd2, factd4, r2
      r2     = 2.d0
      factd2 = fact / r2
      factd4 = fact / 4.d0

      X(1,1) = fact * A(1)*B(1)
      X(2,2) = fact * A(2)*B(2)
      X(3,3) = fact * A(3)*B(3)

      X(1,2) = fact * A(4)*B(4)
      X(1,3) = fact * A(6)*B(6)
      X(2,1) = X(1,2) 
      X(2,3) = fact * A(5)*B(5)
      X(3,1) = X(1,3)
      X(3,2) = X(2,3)

      X(1,4) = factd2 * (A(4)*B(1) + A(1)*B(4))
      X(4,1) = X(1,4)
      X(1,5) = factd2 * (A(4)*B(6) + A(6)*B(4))
      X(5,1) = X(1,5)
      X(1,6) = factd2 * (A(1)*B(6) + A(6)*B(1))
      X(6,1) = X(1,6)
      X(2,4) = factd2 * (A(2)*B(4) + A(4)*B(2))
      X(4,2) = X(2,4)
      X(2,5) = factd2 * (A(2)*B(5) + A(5)*B(2))
      X(5,2) = X(2,5)
      X(2,6) = factd2 * (A(4)*B(5) + A(5)*B(4))
      X(6,2) = X(2,6)
      X(3,4) = factd2 * (A(5)*B(6) + A(6)*B(5))
      X(4,3) = X(3,4)
      X(3,5) = factd2 * (A(3)*B(5) + A(5)*B(3))
      X(5,3) = X(3,5)
      X(3,6) = factd2 * (A(3)*B(6) + A(6)*B(3))
      X(6,3) = X(3,6)

      X(4,4) = factd4 * (A(1)*B(2) + r2 * A(4)*B(4) + A(2)*B(1))
      X(5,5) = factd4 * (A(2)*B(3) + r2 * A(5)*B(5) + A(3)*B(2))
      X(6,6) = factd4 * (A(3)*B(1) + r2 * A(6)*B(6) + A(1)*B(3))

      X(4,5) = factd4 * (A(4)*B(5) + A(5)*B(4) + A(2)*B(6) + A(6)*B(2))
      X(5,4) = X(4,5)
      X(4,6) = factd4 * (A(4)*B(6) + A(6)*B(4) + A(1)*B(5) + A(5)*B(1))
      X(6,4) = X(4,6)
      X(5,6) = factd4 * (A(5)*B(6) + A(6)*B(5) + A(3)*B(4) + A(4)*B(3))
      X(6,5) = X(5,6)

      return
      end
c
c=============================================================================
c  Differentiation
c                     D UTA       mit A,T,UTA sym., U unsym. 2.Stufe
c                X = -------      
c                      D T        -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_UTA_p(X,A,B)
      implicit none
      double precision A(3,3), B(6), X(6,6), r1d2, r1d4
      data r1d4, r1d2 / 0.25d0, 0.5d0 /
      X(1,1) = A(1,1) * B(1)
      X(1,2) = A(1,2) * B(4)
      X(1,3) = A(1,3) * B(6)
      X(2,1) = A(2,1) * B(4)
      X(2,2) = A(2,2) * B(2)
      X(2,3) = A(2,3) * B(5)
      X(3,1) = A(3,1) * B(6)
      X(3,2) = A(3,2) * B(5)
      X(3,3) = A(3,3) * B(3)

      X(1,4) = r1d2 * (A(1,1)*B(4) + A(1,2)*B(1))
      X(2,4) = r1d2 * (A(2,1)*B(2) + A(2,2)*B(4))
      X(3,4) = r1d2 * (A(3,1)*B(5) + A(3,2)*B(6))
      X(1,5) = r1d2 * (A(1,2)*B(6) + A(1,3)*B(4))
      X(2,5) = r1d2 * (A(2,2)*B(5) + A(2,3)*B(2))
      X(3,5) = r1d2 * (A(3,2)*B(3) + A(3,3)*B(5))
      X(1,6) = r1d2 * (A(1,3)*B(1) + A(1,1)*B(6))
      X(2,6) = r1d2 * (A(2,3)*B(4) + A(2,1)*B(5))
      X(3,6) = r1d2 * (A(3,3)*B(6) + A(3,1)*B(3))
      
      X(4,1) = r1d2 * (A(1,1)*B(4) + A(2,1)*B(1))
      X(4,2) = r1d2 * (A(1,2)*B(2) + A(2,2)*B(4))
      X(4,3) = r1d2 * (A(1,3)*B(5) + A(2,3)*B(6))
      X(5,1) = r1d2 * (A(2,1)*B(6) + A(3,1)*B(4))
      X(5,2) = r1d2 * (A(2,2)*B(5) + A(3,2)*B(2))
      X(5,3) = r1d2 * (A(2,3)*B(3) + A(3,3)*B(5))
      X(6,1) = r1d2 * (A(3,1)*B(1) + A(1,1)*B(6))
      X(6,2) = r1d2 * (A(3,2)*B(4) + A(1,2)*B(5))
      X(6,3) = r1d2 * (A(3,3)*B(6) + A(1,3)*B(3))

      X(4,4) = r1d4 * (A(1,1)*B(2)+ (A(1,2)+A(2,1))*B(4)  +A(2,2)*B(1))
      X(4,5) = r1d4 * (A(1,2)*B(5)+A(1,3)*B(2)+A(2,2)*B(6)+A(2,3)*B(4))
      X(4,6) = r1d4 * (A(1,3)*B(4)+A(1,1)*B(5)+A(2,3)*B(1)+A(2,1)*B(6))
      X(5,4) = r1d4 * (A(2,1)*B(5)+A(2,2)*B(6)+A(3,1)*B(2)+A(3,2)*B(4))
      X(5,5) = r1d4 * (A(2,2)*B(3)+ (A(2,3)+A(3,2))*B(5)  +A(3,3)*B(2))
      X(5,6) = r1d4 * (A(2,3)*B(6)+A(2,1)*B(3)+A(3,3)*B(4)+A(3,1)*B(5))
      X(6,4) = r1d4 * (A(3,1)*B(4)+A(3,2)*B(1)+A(1,1)*B(5)+A(1,2)*B(6))
      X(6,5) = r1d4 * (A(3,2)*B(6)+A(3,3)*B(4)+A(1,2)*B(3)+A(1,3)*B(5))
      X(6,6) = r1d4 * (A(3,3)*B(1)+ (A(3,1)+A(1,3))*B(6)  +A(1,1)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTA_m(X,A,B)
      implicit none
      double precision A(3,3), B(6), X(6,6), r1d2, r1d4
      data r1d4, r1d2 / -0.25d0, -0.5d0 /
      X(1,1) = - A(1,1) * B(1)
      X(1,2) = - A(1,2) * B(4)
      X(1,3) = - A(1,3) * B(6)
      X(2,1) = - A(2,1) * B(4)
      X(2,2) = - A(2,2) * B(2)
      X(2,3) = - A(2,3) * B(5)
      X(3,1) = - A(3,1) * B(6)
      X(3,2) = - A(3,2) * B(5)
      X(3,3) = - A(3,3) * B(3)
      
      X(4,1) = r1d2 * (A(1,1)*B(4) + A(2,1)*B(1))
      X(4,2) = r1d2 * (A(1,2)*B(2) + A(2,2)*B(4))
      X(4,3) = r1d2 * (A(1,3)*B(5) + A(2,3)*B(6))
      X(5,1) = r1d2 * (A(2,1)*B(6) + A(3,1)*B(4))
      X(5,2) = r1d2 * (A(2,2)*B(5) + A(3,2)*B(2))
      X(5,3) = r1d2 * (A(2,3)*B(3) + A(3,3)*B(5))
      X(6,1) = r1d2 * (A(3,1)*B(1) + A(1,1)*B(6))
      X(6,2) = r1d2 * (A(3,2)*B(4) + A(1,2)*B(5))
      X(6,3) = r1d2 * (A(3,3)*B(6) + A(1,3)*B(3))

      X(1,4) = r1d2 * (A(1,1)*B(4) + A(1,2)*B(1))
      X(2,4) = r1d2 * (A(2,1)*B(2) + A(2,2)*B(4))
      X(3,4) = r1d2 * (A(3,1)*B(5) + A(3,2)*B(6))
      X(1,5) = r1d2 * (A(1,2)*B(6) + A(1,3)*B(4))
      X(2,5) = r1d2 * (A(2,2)*B(5) + A(2,3)*B(2))
      X(3,5) = r1d2 * (A(3,2)*B(3) + A(3,3)*B(5))
      X(1,6) = r1d2 * (A(1,3)*B(1) + A(1,1)*B(6))
      X(2,6) = r1d2 * (A(2,3)*B(4) + A(2,1)*B(5))
      X(3,6) = r1d2 * (A(3,3)*B(6) + A(3,1)*B(3))

      X(4,4) = r1d4 * (A(1,1)*B(2)+ (A(1,2)+A(2,1))*B(4)  +A(2,2)*B(1))
      X(4,5) = r1d4 * (A(1,2)*B(5)+A(1,3)*B(2)+A(2,2)*B(6)+A(2,3)*B(4))
      X(4,6) = r1d4 * (A(1,3)*B(4)+A(1,1)*B(5)+A(2,3)*B(1)+A(2,1)*B(6))
      X(5,4) = r1d4 * (A(2,1)*B(5)+A(2,2)*B(6)+A(3,1)*B(2)+A(3,2)*B(4))
      X(5,5) = r1d4 * (A(2,2)*B(3)+ (A(2,3)+A(3,2))*B(5)  +A(3,3)*B(2))
      X(5,6) = r1d4 * (A(2,3)*B(6)+A(2,1)*B(3)+A(3,3)*B(4)+A(3,1)*B(5))
      X(6,4) = r1d4 * (A(3,1)*B(4)+A(3,2)*B(1)+A(1,1)*B(5)+A(1,2)*B(6))
      X(6,5) = r1d4 * (A(3,2)*B(6)+A(3,3)*B(4)+A(1,2)*B(3)+A(1,3)*B(5))
      X(6,6) = r1d4 * (A(3,3)*B(1)+ (A(3,1)+A(1,3))*B(6)  +A(1,1)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTA_f(X,A,B,fact)
      implicit none
      double precision A(3,3), B(6), X(6,6), fact, factd2, factd4
      factd2 = fact * 0.50d0
      factd4 = fact * 0.25d0
      X(1,1) = fact * A(1,1) * B(1)
      X(1,2) = fact * A(1,2) * B(4)
      X(1,3) = fact * A(1,3) * B(6)
      X(2,1) = fact * A(2,1) * B(4)
      X(2,2) = fact * A(2,2) * B(2)
      X(2,3) = fact * A(2,3) * B(5)
      X(3,1) = fact * A(3,1) * B(6)
      X(3,2) = fact * A(3,2) * B(5)
      X(3,3) = fact * A(3,3) * B(3)

      X(1,4) = factd2 * (A(1,1)*B(4) + A(1,2)*B(1))
      X(2,4) = factd2 * (A(2,1)*B(2) + A(2,2)*B(4))
      X(3,4) = factd2 * (A(3,1)*B(5) + A(3,2)*B(6))
      X(1,5) = factd2 * (A(1,2)*B(6) + A(1,3)*B(4))
      X(2,5) = factd2 * (A(2,2)*B(5) + A(2,3)*B(2))
      X(3,5) = factd2 * (A(3,2)*B(3) + A(3,3)*B(5))
      X(1,6) = factd2 * (A(1,3)*B(1) + A(1,1)*B(6))
      X(2,6) = factd2 * (A(2,3)*B(4) + A(2,1)*B(5))
      X(3,6) = factd2 * (A(3,3)*B(6) + A(3,1)*B(3))
      
      X(4,1) = factd2 * (A(1,1)*B(4) + A(2,1)*B(1))
      X(4,2) = factd2 * (A(1,2)*B(2) + A(2,2)*B(4))
      X(4,3) = factd2 * (A(1,3)*B(5) + A(2,3)*B(6))
      X(5,1) = factd2 * (A(2,1)*B(6) + A(3,1)*B(4))
      X(5,2) = factd2 * (A(2,2)*B(5) + A(3,2)*B(2))
      X(5,3) = factd2 * (A(2,3)*B(3) + A(3,3)*B(5))
      X(6,1) = factd2 * (A(3,1)*B(1) + A(1,1)*B(6))
      X(6,2) = factd2 * (A(3,2)*B(4) + A(1,2)*B(5))
      X(6,3) = factd2 * (A(3,3)*B(6) + A(1,3)*B(3))

      X(4,4) = factd4 *(A(1,1)*B(2)+ (A(1,2)+A(2,1))*B(4)  +A(2,2)*B(1))
      X(4,5) = factd4 *(A(1,2)*B(5)+A(1,3)*B(2)+A(2,2)*B(6)+A(2,3)*B(4))
      X(4,6) = factd4 *(A(1,3)*B(4)+A(1,1)*B(5)+A(2,3)*B(1)+A(2,1)*B(6))
      X(5,4) = factd4 *(A(2,1)*B(5)+A(2,2)*B(6)+A(3,1)*B(2)+A(3,2)*B(4))
      X(5,5) = factd4 *(A(2,2)*B(3)+ (A(2,3)+A(3,2))*B(5)  +A(3,3)*B(2))
      X(5,6) = factd4 *(A(2,3)*B(6)+A(2,1)*B(3)+A(3,3)*B(4)+A(3,1)*B(5))
      X(6,4) = factd4 *(A(3,1)*B(4)+A(3,2)*B(1)+A(1,1)*B(5)+A(1,2)*B(6))
      X(6,5) = factd4 *(A(3,2)*B(6)+A(3,3)*B(4)+A(1,2)*B(3)+A(1,3)*B(5))
      X(6,6) = factd4 *(A(3,3)*B(1)+ (A(3,1)+A(1,3))*B(6)  +A(1,1)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTA_ap(X,A,B)
      implicit none
      double precision A(3,3), B(6), X(6,6), r1d2, r1d4
      data r1d4, r1d2 / 0.25d0, 0.5d0 /
      X(1,1) = X(1,1) + A(1,1) * B(1)
      X(1,2) = X(1,2) + A(1,2) * B(4)
      X(1,3) = X(1,3) + A(1,3) * B(6)
      X(2,1) = X(2,1) + A(2,1) * B(4)
      X(2,2) = X(2,2) + A(2,2) * B(2)
      X(2,3) = X(2,3) + A(2,3) * B(5)
      X(3,1) = X(3,1) + A(3,1) * B(6)
      X(3,2) = X(3,2) + A(3,2) * B(5)
      X(3,3) = X(3,3) + A(3,3) * B(3)
      
      X(4,1) = X(4,1) + r1d2 * (A(1,1)*B(4) + A(2,1)*B(1))
      X(4,2) = X(4,2) + r1d2 * (A(1,2)*B(2) + A(2,2)*B(4))
      X(4,3) = X(4,3) + r1d2 * (A(1,3)*B(5) + A(2,3)*B(6))
      X(5,1) = X(5,1) + r1d2 * (A(2,1)*B(6) + A(3,1)*B(4))
      X(5,2) = X(5,2) + r1d2 * (A(2,2)*B(5) + A(3,2)*B(2))
      X(5,3) = X(5,3) + r1d2 * (A(2,3)*B(3) + A(3,3)*B(5))
      X(6,1) = X(6,1) + r1d2 * (A(3,1)*B(1) + A(1,1)*B(6))
      X(6,2) = X(6,2) + r1d2 * (A(3,2)*B(4) + A(1,2)*B(5))
      X(6,3) = X(6,3) + r1d2 * (A(3,3)*B(6) + A(1,3)*B(3))

      X(1,4) = X(1,4) + r1d2 * (A(1,1)*B(4) + A(1,2)*B(1))
      X(2,4) = X(2,4) + r1d2 * (A(2,1)*B(2) + A(2,2)*B(4))
      X(3,4) = X(3,4) + r1d2 * (A(3,1)*B(5) + A(3,2)*B(6))
      X(1,5) = X(1,5) + r1d2 * (A(1,2)*B(6) + A(1,3)*B(4))
      X(2,5) = X(2,5) + r1d2 * (A(2,2)*B(5) + A(2,3)*B(2))
      X(3,5) = X(3,5) + r1d2 * (A(3,2)*B(3) + A(3,3)*B(5))
      X(1,6) = X(1,6) + r1d2 * (A(1,3)*B(1) + A(1,1)*B(6))
      X(2,6) = X(2,6) + r1d2 * (A(2,3)*B(4) + A(2,1)*B(5))
      X(3,6) = X(3,6) + r1d2 * (A(3,3)*B(6) + A(3,1)*B(3))

      X(4,4) = X(4,4)
     *    + r1d4 * (A(1,1)*B(2)+ (A(1,2)+A(2,1))*B(4)  +A(2,2)*B(1))
      X(4,5) = X(4,5)
     *    + r1d4 * (A(1,2)*B(5)+A(1,3)*B(2)+A(2,2)*B(6)+A(2,3)*B(4))
      X(4,6) = X(4,6)
     *    + r1d4 * (A(1,3)*B(4)+A(1,1)*B(5)+A(2,3)*B(1)+A(2,1)*B(6))
      X(5,4) = X(5,4)
     *    + r1d4 * (A(2,1)*B(5)+A(2,2)*B(6)+A(3,1)*B(2)+A(3,2)*B(4))
      X(5,5) = X(5,5)
     *    + r1d4 * (A(2,2)*B(3)+ (A(2,3)+A(3,2))*B(5)  +A(3,3)*B(2))
      X(5,6) = X(5,6)
     *    + r1d4 * (A(2,3)*B(6)+A(2,1)*B(3)+A(3,3)*B(4)+A(3,1)*B(5))
      X(6,4) = X(6,4)
     *    + r1d4 * (A(3,1)*B(4)+A(3,2)*B(1)+A(1,1)*B(5)+A(1,2)*B(6))
      X(6,5) = X(6,5)
     *    + r1d4 * (A(3,2)*B(6)+A(3,3)*B(4)+A(1,2)*B(3)+A(1,3)*B(5))
      X(6,6) = X(6,6)
     *    + r1d4 * (A(3,3)*B(1)+ (A(3,1)+A(1,3))*B(6)  +A(1,1)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTA_am(X,A,B)
      implicit none
      double precision A(3,3), B(6), X(6,6), r1d2, r1d4
      data r1d4, r1d2 / 0.25d0, 0.5d0 /
      X(1,1) = X(1,1) - A(1,1) * B(1)
      X(1,2) = X(1,2) - A(1,2) * B(4)
      X(1,3) = X(1,3) - A(1,3) * B(6)
      X(2,1) = X(2,1) - A(2,1) * B(4)
      X(2,2) = X(2,2) - A(2,2) * B(2)
      X(2,3) = X(2,3) - A(2,3) * B(5)
      X(3,1) = X(3,1) - A(3,1) * B(6)
      X(3,2) = X(3,2) - A(3,2) * B(5)
      X(3,3) = X(3,3) - A(3,3) * B(3)
      
      X(4,1) = X(4,1) - r1d2 * (A(1,1)*B(4) + A(2,1)*B(1))
      X(4,2) = X(4,2) - r1d2 * (A(1,2)*B(2) + A(2,2)*B(4))
      X(4,3) = X(4,3) - r1d2 * (A(1,3)*B(5) + A(2,3)*B(6))
      X(5,1) = X(5,1) - r1d2 * (A(2,1)*B(6) + A(3,1)*B(4))
      X(5,2) = X(5,2) - r1d2 * (A(2,2)*B(5) + A(3,2)*B(2))
      X(5,3) = X(5,3) - r1d2 * (A(2,3)*B(3) + A(3,3)*B(5))
      X(6,1) = X(6,1) - r1d2 * (A(3,1)*B(1) + A(1,1)*B(6))
      X(6,2) = X(6,2) - r1d2 * (A(3,2)*B(4) + A(1,2)*B(5))
      X(6,3) = X(6,3) - r1d2 * (A(3,3)*B(6) + A(1,3)*B(3))

      X(1,4) = X(1,4) - r1d2 * (A(1,1)*B(4) + A(1,2)*B(1))
      X(2,4) = X(2,4) - r1d2 * (A(2,1)*B(2) + A(2,2)*B(4))
      X(3,4) = X(3,4) - r1d2 * (A(3,1)*B(5) + A(3,2)*B(6))
      X(1,5) = X(1,5) - r1d2 * (A(1,2)*B(6) + A(1,3)*B(4))
      X(2,5) = X(2,5) - r1d2 * (A(2,2)*B(5) + A(2,3)*B(2))
      X(3,5) = X(3,5) - r1d2 * (A(3,2)*B(3) + A(3,3)*B(5))
      X(1,6) = X(1,6) - r1d2 * (A(1,3)*B(1) + A(1,1)*B(6))
      X(2,6) = X(2,6) - r1d2 * (A(2,3)*B(4) + A(2,1)*B(5))
      X(3,6) = X(3,6) - r1d2 * (A(3,3)*B(6) + A(3,1)*B(3))

      X(4,4) = X(4,4)
     *    - r1d4 * (A(1,1)*B(2)+ (A(1,2)+A(2,1))*B(4)  +A(2,2)*B(1))
      X(4,5) = X(4,5)
     *    - r1d4 * (A(1,2)*B(5)+A(1,3)*B(2)+A(2,2)*B(6)+A(2,3)*B(4))
      X(4,6) = X(4,6)
     *    - r1d4 * (A(1,3)*B(4)+A(1,1)*B(5)+A(2,3)*B(1)+A(2,1)*B(6))
      X(5,4) = X(5,4)
     *    - r1d4 * (A(2,1)*B(5)+A(2,2)*B(6)+A(3,1)*B(2)+A(3,2)*B(4))
      X(5,5) = X(5,5)
     *    - r1d4 * (A(2,2)*B(3)+ (A(2,3)+A(3,2))*B(5)  +A(3,3)*B(2))
      X(5,6) = X(5,6)
     *    - r1d4 * (A(2,3)*B(6)+A(2,1)*B(3)+A(3,3)*B(4)+A(3,1)*B(5))
      X(6,4) = X(6,4)
     *    - r1d4 * (A(3,1)*B(4)+A(3,2)*B(1)+A(1,1)*B(5)+A(1,2)*B(6))
      X(6,5) = X(6,5)
     *    - r1d4 * (A(3,2)*B(6)+A(3,3)*B(4)+A(1,2)*B(3)+A(1,3)*B(5))
      X(6,6) = X(6,6)
     *    - r1d4 * (A(3,3)*B(1)+ (A(3,1)+A(1,3))*B(6)  +A(1,1)*B(3))
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_UTA_af(X,A,B,fact)
      implicit none
      double precision A(3,3), B(6), X(6,6), fact, factd2, factd4
      factd2 = fact * 0.50d0
      factd4 = fact * 0.25d0
      X(1,1) = X(1,1) + fact * A(1,1) * B(1)
      X(1,2) = X(1,2) + fact * A(1,2) * B(4)
      X(1,3) = X(1,3) + fact * A(1,3) * B(6)
      X(2,1) = X(2,1) + fact * A(2,1) * B(4)
      X(2,2) = X(2,2) + fact * A(2,2) * B(2)
      X(2,3) = X(2,3) + fact * A(2,3) * B(5)
      X(3,1) = X(3,1) + fact * A(3,1) * B(6)
      X(3,2) = X(3,2) + fact * A(3,2) * B(5)
      X(3,3) = X(3,3) + fact * A(3,3) * B(3)

      X(1,4) = X(1,4) + factd2 * (A(1,1)*B(4) + A(1,2)*B(1))
      X(2,4) = X(2,4) + factd2 * (A(2,1)*B(2) + A(2,2)*B(4))
      X(3,4) = X(3,4) + factd2 * (A(3,1)*B(5) + A(3,2)*B(6))
      X(1,5) = X(1,5) + factd2 * (A(1,2)*B(6) + A(1,3)*B(4))
      X(2,5) = X(2,5) + factd2 * (A(2,2)*B(5) + A(2,3)*B(2))
      X(3,5) = X(3,5) + factd2 * (A(3,2)*B(3) + A(3,3)*B(5))
      X(1,6) = X(1,6) + factd2 * (A(1,3)*B(1) + A(1,1)*B(6))
      X(2,6) = X(2,6) + factd2 * (A(2,3)*B(4) + A(2,1)*B(5))
      X(3,6) = X(3,6) + factd2 * (A(3,3)*B(6) + A(3,1)*B(3))
      
      X(4,1) = X(4,1) + factd2 * (A(1,1)*B(4) + A(2,1)*B(1))
      X(4,2) = X(4,2) + factd2 * (A(1,2)*B(2) + A(2,2)*B(4))
      X(4,3) = X(4,3) + factd2 * (A(1,3)*B(5) + A(2,3)*B(6))
      X(5,1) = X(5,1) + factd2 * (A(2,1)*B(6) + A(3,1)*B(4))
      X(5,2) = X(5,2) + factd2 * (A(2,2)*B(5) + A(3,2)*B(2))
      X(5,3) = X(5,3) + factd2 * (A(2,3)*B(3) + A(3,3)*B(5))
      X(6,1) = X(6,1) + factd2 * (A(3,1)*B(1) + A(1,1)*B(6))
      X(6,2) = X(6,2) + factd2 * (A(3,2)*B(4) + A(1,2)*B(5))
      X(6,3) = X(6,3) + factd2 * (A(3,3)*B(6) + A(1,3)*B(3))

      X(4,4) = X(4,4)
     *      + factd4 * (A(1,1)*B(2)+ (A(1,2)+A(2,1))*B(4)  +A(2,2)*B(1))
      X(4,5) = X(4,5)
     *      + factd4 * (A(1,2)*B(5)+A(1,3)*B(2)+A(2,2)*B(6)+A(2,3)*B(4))
      X(4,6) = X(4,6)
     *      + factd4 * (A(1,3)*B(4)+A(1,1)*B(5)+A(2,3)*B(1)+A(2,1)*B(6))
      X(5,4) = X(5,4)
     *      + factd4 * (A(2,1)*B(5)+A(2,2)*B(6)+A(3,1)*B(2)+A(3,2)*B(4))
      X(5,5) = X(5,5)
     *      + factd4 * (A(2,2)*B(3)+ (A(2,3)+A(3,2))*B(5)  +A(3,3)*B(2))
      X(5,6) = X(5,6)
     *      + factd4 * (A(2,3)*B(6)+A(2,1)*B(3)+A(3,3)*B(4)+A(3,1)*B(5))
      X(6,4) = X(6,4)
     *      + factd4 * (A(3,1)*B(4)+A(3,2)*B(1)+A(1,1)*B(5)+A(1,2)*B(6))
      X(6,5) = X(6,5)
     *      + factd4 * (A(3,2)*B(6)+A(3,3)*B(4)+A(1,2)*B(3)+A(1,3)*B(5))
      X(6,6) = X(6,6)
     *      + factd4 * (A(3,3)*B(1)+ (A(3,1)+A(1,3))*B(6)  +A(1,1)*B(3))
      return
      end
c
c=============================================================================
c  Differentiation
c                     D dev(AT)B       mit A,T,B,dev(AT)B  sym. 2.Stufe
c                X = ------------      
c                        D T           -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_ATdB_p(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d3
      r1d3 = - 1.d0 / 3.d0
      call diff_ATB_p(X,A,B)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATdB_m(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d3
      r1d3 = 1.d0 / 3.d0
      call diff_ATB_m(X,A,B)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATdB_f(X,A,B,fact)
      implicit none
      double precision A(6), B(6), X(6,6), fact, r1d3
      r1d3 = - fact / 3.d0
      call diff_ATB_f(X,A,B,fact)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATdB_ap(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d3
      r1d3 = - 1.d0 / 3.d0
      call diff_ATB_ap(X,A,B)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATdB_am(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), r1d3
      r1d3 = 1.d0 / 3.d0
      call diff_ATB_am(X,A,B)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_ATdB_af(X,A,B,fact)
      implicit none
      double precision A(6), B(6), X(6,6), fact, r1d3
      r1d3 = - fact / 3.d0
      call diff_ATB_af(X,A,B,fact)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c
c=============================================================================
c  Differentiation
c                     D dev(TA)B       mit A,T,B,dev(TA)B  sym. 2.Stufe
c                X = ------------      
c                        D T           -> Tensor 4.Stufe -> 6x6 Matrix
c-----------------------------------------------------------------------------
      subroutine diff_TAdB_p(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), hu(3,3), r1d3
      r1d3 = - 1.d0 / 3.d0
      call mult_u_ss_p(hu,A,B)
      call diff_TU_p(X,hu)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAdB_m(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), hu(3,3), r1d3
      r1d3 = 1.d0 / 3.d0
      call mult_u_ss_p(hu,A,B)
      call diff_TU_m(X,hu)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAdB_f(X,A,B,fact)
      implicit none
      double precision A(6), B(6), X(6,6), hu(3,3), fact, r1d3
      r1d3 = - fact / 3.d0
      call mult_u_ss_f(hu,A,B,fact)
      call diff_TU_p(X,hu)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAdB_ap(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), hu(3,3), r1d3
      r1d3 = - 1.d0 / 3.d0
      call mult_u_ss_p(hu,A,B)
      call diff_TU_ap(X,hu)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAdB_am(X,A,B)
      implicit none
      double precision A(6), B(6), X(6,6), hu(3,3), r1d3
      r1d3 = 1.d0 / 3.d0
      call mult_u_ss_p(hu,A,B)
      call diff_TU_am(X,hu)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c-----------------------------------------------------------------------------
      subroutine diff_TAdB_af(X,A,B,fact)
      implicit none
      double precision A(6), B(6), X(6,6), hu(3,3), fact, r1d3
      r1d3 = - fact / 3.d0
      call mult_u_ss_f(hu,A,B,fact)
      call diff_TU_ap(X,hu)
      call mult_4_ss_af(X,B,A,r1d3)
      return
      end
c=============================================================================

c*****************************************************************************
c*****************************************************************************

c=============================================================================
c  differentiation
c                     D TATt       where  A sym. 2.order, T unsym. 2.order
c                X = ------------      
c                       D T           -> 4.order tensor -> 6x3x3 matrix
c-----------------------------------------------------------------------------
      subroutine diff_TATt_p(X,TA)
      implicit none
      double precision X(6,3,3), TA(3,3), r0

      r0 = 0.d0

      X(1,1,1) = TA(1,1) + TA(1,1)
      X(1,1,2) = TA(1,2) + TA(1,2)
      X(1,1,3) = TA(1,3) + TA(1,3)
      X(1,2,1) = r0
      X(1,2,2) = r0
      X(1,2,3) = r0
      X(1,3,1) = r0
      X(1,3,2) = r0
      X(1,3,3) = r0

      X(2,1,1) = r0
      X(2,1,2) = r0
      X(2,1,3) = r0
      X(2,2,1) = TA(2,1) + TA(2,1)
      X(2,2,2) = TA(2,2) + TA(2,2)
      X(2,2,3) = TA(2,3) + TA(2,3)
      X(2,3,1) = r0
      X(2,3,2) = r0
      X(2,3,3) = r0

      X(3,1,1) = r0
      X(3,1,2) = r0
      X(3,1,3) = r0
      X(3,2,1) = r0
      X(3,2,2) = r0
      X(3,2,3) = r0
      X(3,3,1) = TA(3,1) + TA(3,1)
      X(3,3,2) = TA(3,2) + TA(3,2)
      X(3,3,3) = TA(3,3) + TA(3,3)

      X(4,1,1) = TA(2,1)
      X(4,1,2) = TA(2,2)
      X(4,1,3) = TA(2,3)
      X(4,2,1) = TA(1,1)
      X(4,2,2) = TA(1,2)
      X(4,2,3) = TA(1,3)
      X(4,3,1) = r0
      X(4,3,2) = r0
      X(4,3,3) = r0

      X(5,1,1) = r0
      X(5,1,2) = r0
      X(5,1,3) = r0
      X(5,2,1) = TA(3,1)
      X(5,2,2) = TA(3,2)
      X(5,2,3) = TA(3,3)
      X(5,3,1) = TA(2,1)
      X(5,3,2) = TA(2,2)
      X(5,3,3) = TA(2,3)

      X(6,1,1) = TA(3,1)
      X(6,1,2) = TA(3,2)
      X(6,1,3) = TA(3,3)
      X(6,2,1) = r0
      X(6,2,2) = r0
      X(6,2,3) = r0
      X(6,3,1) = TA(1,1)
      X(6,3,2) = TA(1,2)
      X(6,3,3) = TA(1,3)

      return
      end

c=============================================================================

c*****************************************************************************
c*****************************************************************************

c=============================================================================
c
c   spectral decomposition of a symmetric second order tensor
c
c                            
c   A = Sum x(i) * E(i)               A, E(i)  sym., 2.order
c
c   input: A(6)   tensor A
c
c   output:  x(3)  eigenvalues of A
c           E(6,3) eigenprojections of A
c
c-----------------------------------------------------------------------------
      integer function spectralDecomposition(x,E,A)
      implicit none
      double precision x(3), E(6,3), A(6),
     *                 H(3,3), V(3,3), L(6), Vn(3,3),
     *                 Lmx, fact, fact1, s, c,
     *                 r1d2, r1dsqrt2

      integer i, j, k, ij, iter

      r1d2     = .5d0
      r1dsqrt2 = 1.d0 / sqrt(2.d0)

      call unit_u_p(V)

      call set_s_p(L,A)

      iter = 0

 100  continue

      iter = iter + 1

      if (iter.gt.30) goto 300

      Lmx = abs(L(4))
      i   = 1
      if (Lmx < abs(L(5))) then
        Lmx = abs(L(5))
        i = 2
      endif
      if (Lmx < abs(L(6))) then
        Lmx = abs(L(6))
        i = 3
      endif

      if (Lmx .lt. 1.d-12) goto 200

      call pmove(V,Vn,9)

      j  = i + 1
      if (j.eq.4) j = 1
      ij = 3 + i
 
      fact  = r1d2 * abs(L(i) - L(j))

      fact1 = sqrt(L(ij)*L(ij) + fact*fact)

      fact = sqrt((fact + fact1) / fact1)

      c = r1dsqrt2 * fact

      s = L(ij) * r1dsqrt2 / (fact1 * fact)

      if (L(j).gt.L(i)) s = -s

      k = j + 1
      if (k.eq.4) k = 1
 
      V(i,i) = c * Vn(i,i) + s * Vn(j,i)
      V(i,j) = c * Vn(i,j) + s * Vn(j,j)
      V(i,k) = c * Vn(i,k) + s * Vn(j,k)
      V(j,i) = c * Vn(j,i) - s * Vn(i,i)
      V(j,j) = c * Vn(j,j) - s * Vn(i,j)
      V(j,k) = c * Vn(j,k) - s * Vn(i,k)
 
      call mult_u_us_p(H,V,A)
 
      call mult_s_uut_p(L,H,V)

      goto 100

 200  continue

      spectralDecomposition = 1

      x(1) = L(1)
      x(2) = L(2)
      x(3) = L(3)

      do j=1, 3
        E(1,j) = V(j,1) * V(j,1)
        E(2,j) = V(j,2) * V(j,2)
        E(3,j) = V(j,3) * V(j,3)
        E(4,j) = V(j,1) * V(j,2)
        E(5,j) = V(j,2) * V(j,3)
        E(6,j) = V(j,3) * V(j,1)
      enddo

      return

 300  continue

      call prgWarning(1,'spectralDecomposition','no convergence!')

      spectralDecomposition = 0

      return

      end

c=============================================================================

c*****************************************************************************
c*****************************************************************************

c=============================================================================
c
c   isotropic tensor function: EXPONENTIAL FUNCTION
c
c                      d exp(A)
c   S = exp(A),   X = ----------          A, S sym., 2.order, X 4.order
c                        d A
c-----------------------------------------------------------------------------
      subroutine exp_s(S,X,A,n)
      implicit none
      integer          n, i, j
      double precision A(6), S(6), X(6,6), 
     *                 H2(6), H22(6), fact, H4(6,6),
     *                 H42(6,6), dTAdT(6,6)

      call unit_s_p (S)
      call set_s_ap (S,A)

      call set_s_p  (H2,A)
      call unit_4_p (H4)

      call diff_AT_p(dTAdT,A) 

      call unit_4_p (X)

      fact = 1.d0

      do i=2, n

        fact = fact / dble(i)

        call mult_s_ss_p(H22,H2,A)
        call set_s_af   (S,H22,fact)

        call mult_4_44_p(H42,dTAdT,H4)
        call diff_AT_ap (H42,H2)
        call set_4_af   (X,H42,fact)

        call set_s_p    (H2,H22)
        call set_4_p    (H4,H42)

      enddo
      
      return
      end

c=============================================================================
c
c   isotropic tensor function: SQUARE
c
c                   d A A
c   S = A A,   X = -------          A, S sym., 2.order, X 4.order
c                    d A
c-----------------------------------------------------------------------------
      subroutine sq_s(S,X,A)
      implicit none
      double precision A(6), S(6), X(6,6), r2

      r2 = 2.d0

      call mult_s_ss_p(S,A,A)

      call diff_AT_f  (X,A,r2)

      return
      end

c=============================================================================
c
c   isotropic tensor function: LOGARITHM
c
c                      d log(A)
c   S = log(A),   X = ----------          A, S sym., 2.order, X 4.order
c                        d A
c-----------------------------------------------------------------------------
      subroutine log_s(S,X,A,n)
      implicit none
      integer          n, P(6), i, j, mxi
      double precision A(6), S(6), X(6,6), 
     *                 R(6), dRdS(6,6), upd(6), Rnorm, Rmax, tol, 
     *                 r0, r1d2

      r0   = 0.d0
      r1d2 = 0.5d0

      tol = 1.d-12

      mxi = 30
      i   =  0

      call set_s_p  (S,A)

 1    continue

      i = i + 1

      call exp_s    (R,dRdS,S,n)
      call set_s_am (R,A)

      Rnorm = r0
      Rmax  = r0
      do j=1, 6
        Rnorm = Rnorm + R(j)*R(j)
        Rmax  = max(Rmax,abs(R(j)))
      enddo

      call decompLR_matrix(dRdS,P,6,0)

      if (sqrt(Rnorm).le.tol.and.Rmax.le.tol) goto 2

      call solve_matrix   (dRdS,P,R,upd,6)
      do j=4, 6
        upd(j) = r1d2 * upd(j)
      enddo
      call set_s_am(S,upd)

      if (i.lt.mxi) goto 1

      write(*,'(/9x,a/)') 'log_s: surprisingly it did not converge!'

 2    continue

      call inverse_matrix(dRdS,X,P,6)
      do i=1, 6
        do j=4, 6
          X(i,j) = X(i,j) * r1d2
          X(j,i) = X(j,i) * r1d2
        enddo
      enddo

      return
      end

c=============================================================================
c
c   isotropic tensor function: SQUARE ROOT
c
c                       d sqrt(A)
c   S = sqrt(A),   X = ----------          A, S sym., 2.order, X 4.order
c                          d A
c-----------------------------------------------------------------------------
      subroutine sqrt_s(S,X,A)
      implicit none
      integer          P(6), i, j, mxi
      double precision A(6), S(6), X(6,6), 
     *                 R(6), dRdS(6,6), upd(6), Rnorm, Rmax, tol, 
     *                 r0, r1d2

      r0   = 0.d0
      r1d2 = 0.5d0

      tol = 1.d-12

      mxi = 30
      i   =  0

      call set_s_p  (S,A)

 1    continue

      i = i + 1

      call sq_s     (R,dRdS,S)
      call set_s_am (R,A)

      Rnorm = r0
      Rmax  = r0
      do j=1, 6
        Rnorm = Rnorm + R(j)*R(j)
        Rmax  = max(Rmax,abs(R(j)))
      enddo

      call decompLR_matrix(dRdS,P,6,0)

      if (sqrt(Rnorm).le.tol.and.Rmax.le.tol) goto 2

      call solve_matrix   (dRdS,P,R,upd,6)
      do j=4, 6
        upd(j) = r1d2 * upd(j)
      enddo
      call set_s_am(S,upd)

      if (i.lt.mxi) goto 1

      write(*,'(/9x,a/)') 'sqrt_s: surprisingly it did not converge!'

 2    continue

      call inverse_matrix(dRdS,X,P,6)
      do i=1, 6
        do j=4, 6
          X(i,j) = X(i,j) * r1d2
          X(j,i) = X(j,i) * r1d2
        enddo
      enddo

      return
      end

c=============================================================================
c
c   generic isotropic tensor function
c
c                       d S
c   S = Func(A),   C = -----      A, S  sym., 2.order, C  4.order
c                       d A
c
c   input: x(3)    eigenvalues of A
c          y(3)    functions of x1,x2,x3, eigenvalues of S
c       dydx(3,3)  derivatives of y with respect to x
c         E(6,3)   eigenprojections of A
c          A(6)    tensor A
c
c   output:  S(6)   tensor S
c           C(6,6)  tensor C
c
c-----------------------------------------------------------------------------
      subroutine isotropicTensorFunction(S,C,x,y,dydx,E,A)

      implicit none

      integer          p, i, j, k, l, m

      double precision S(*), C(6,*), x(*), y(*), dydx(3,*), E(6,*),A(*),
     *                 fact1, fact2, fact3, fact4, fact5, fact6, fact7, 
     *                 fact, dA2(6,6), xa, xc, ya, yc, dxac, dxac2, 
     *                 dxac3, yac, daa, dcc, dac, dca, dcb, dacs, daacc,
     *                 r0, r1d2, r1, r2, tol
 
      logical          flg12, flg23, flg31

      data  r0,  r1d2,   r1,   r2 
     *   / 0.d0, 0.5d0, 1.d0, 2.d0 /


c:::: tensor function :::::::::::::::::::::::

      call pzero(S,6)

      do l=1, 3
        do i=1, 6
          S(i) = S(i) + y(l) * E(i,l)
        enddo
      enddo

c:::: derivative ::::::::::::::::::::::::::::

      tol = 1.d-6

      flg12 = (abs(x(1)-x(2))/(abs(x(1))+abs(x(2))).lt.tol)
      flg23 = (abs(x(2)-x(3))/(abs(x(2))+abs(x(3))).lt.tol)
      flg31 = (abs(x(3)-x(1))/(abs(x(3))+abs(x(1))).lt.tol)

      if (flg12) then                       !
        if (flg23) then                     !  p =  0 -> three distinct
          p = -1  ! x1==x2, x2==x3, x3==x1  !            eigenvalues
        else                                ! 
          p = 3   ! x1==x2, x2<>x3, x3<>x1  !  p >  0 -> two identical 
        endif                               !            eigenvalues, 
      else                                  !            p points to the 
        if (flg23) then                     !            distinct value
          p = 1   ! x1<>x2, x2==x3, x3<>x1  !
        else                                !  p = -1 -> three identical 
          if (flg31) then                   !            eigenvalues
            p = 2 ! x1<>x2, x2<>x3, x3==x1  !
          else                              !
            p = 0 ! x1<>x2, x2<>x3, x3<>x1  !
          endif                             !
        endif                               !
      endif                                 !

      !if (p.ne.0) write(*,*) p

      call pzero(C,36)

      if (p.lt.0) then  ! three identical eigenvalues

        fact1 = dydx(1,1) - dydx(1,2)
        fact2 = fact1 * r1d2
        fact3 = dydx(1,2)

        do i=1, 3
          do j=1, 3
            C(i,j) = fact3
          enddo
          C(i  ,i  ) = C(i,i) + fact1
          C(i+3,i+3) = fact2
        enddo

        return

      endif

      dA2(1,1) = A(1) + A(1)
      dA2(2,2) = A(2) + A(2)
      dA2(3,3) = A(3) + A(3)
      dA2(1,2) = r0
      dA2(1,3) = r0
      dA2(2,1) = r0 
      dA2(2,3) = r0
      dA2(3,1) = r0
      dA2(3,2) = r0
      dA2(1,4) = A(4)
      dA2(4,1) = dA2(1,4)
      dA2(1,5) = r0
      dA2(5,1) = r0
      dA2(1,6) = A(6)
      dA2(6,1) = dA2(1,6)
      dA2(2,4) = A(4)
      dA2(4,2) = dA2(2,4)
      dA2(2,5) = A(5)
      dA2(5,2) = dA2(2,5)
      dA2(2,6) = r0
      dA2(6,2) = r0
      dA2(3,4) = r0
      dA2(4,3) = r0
      dA2(3,5) = A(5)
      dA2(5,3) = dA2(3,5)
      dA2(3,6) = A(6)
      dA2(6,3) = dA2(3,6)
      dA2(4,4) = r1d2 * (A(1) + A(2))
      dA2(5,5) = r1d2 * (A(2) + A(3))
      dA2(6,6) = r1d2 * (A(3) + A(1))
      dA2(4,5) = r1d2 * A(6)
      dA2(5,4) = dA2(4,5)
      dA2(4,6) = r1d2 * A(5)
      dA2(6,4) = dA2(4,6)
      dA2(5,6) = r1d2 * A(4)
      dA2(6,5) = dA2(5,6)

      if (p.eq.0) then ! three different eigenvalues

        do k=1, 3
          l = mod(k,3) + 1
          m = mod(l,3) + 1
          fact  = (x(k)-x(l))*(x(k)-x(m))
          if (abs(fact).lt.1.d-12) fact = sign(1.d-12,fact)
          fact1 = y(k) / fact
          fact2 = x(l) + x(m)
          fact3 = x(k)-x(l) + x(k)-x(m)
          fact4 = x(l) - x(m)
          fact5 = fact1 * fact2
          fact6 = r1d2 * fact5
          do i=1, 3
            C(i  ,i  ) = C(i  ,i  ) - fact5
            C(i+3,i+3) = C(i+3,i+3) - fact6
          enddo
          do i=1, 6
            do j=1, 6
              C(i,j) = C(i,j)
     *                + fact1 * (dA2(i,j) 
     *                          - fact3 *  E(i,k)*E(j,k)
     *                          - fact4 * (E(i,l)*E(j,l)-E(i,m)*E(j,m)))
            enddo
          enddo
        enddo
        do k=1, 3
          do l=1, 3
            do i=1, 6
              do j=1, 6
                C(i,j) = C(i,j) + dydx(k,l) * E(i,k)*E(j,l)
              enddo
            enddo
          enddo
        enddo

        return
  
      endif

      ! two identical eigenvalues, x(k) <> x(l) == x(m)

      k     = p
      l     = mod(k,3) + 1
      m     = mod(l,3) + 1

      xa    = x(k)
      xc    = x(m)
      ya    = y(k)
      yc    = y(m)
      daa   = dydx(k,k)
      dcc   = dydx(m,m)
      dac   = dydx(k,m)
      dca   = dydx(m,k)
      dcb   = dydx(m,l)

      dxac  = r1 / (xa - xc)
      dxac2 = dxac * dxac
      dxac3 = dxac * dxac2
      yac   = ya - yc
      dacs  = dac + dca
      daacc = daa + dcc

      fact  = r2*xc*yac*dxac3 + xc*dxac2*(dacs-daacc)
      fact1 = yac*dxac2 + dxac*(dcb-dcc)
      fact2 = r2*xc*yac*dxac2 + (xa+xc)*dxac*(dcb-dcc)
      fact3 = r2*yac*dxac3 + dxac2*(dacs-daacc)
      fact4 = fact + dxac*(dac-dcb)
      fact5 = fact + dxac*(dca-dcb)
      fact6 = r2*xc*xc*yac*dxac3 + xa*xc*dxac2*dacs 
     *        - xc*xc*dxac2*daacc - (xa+xc)*dxac*dcb
      fact7 = r1d2 * fact2

      do i=1, 3
        C(i  ,i  ) = C(i  ,i  ) - fact2
        C(i+3,i+3) = C(i+3,i+3) - fact7
        do j=1, 3
          C(i,j) = C(i,j) - fact6
        enddo
      enddo

      do i=1, 6
        do j=1, 6
          C(i,j) = C(i,j) + fact1 * dA2(i,j) - fact3 * A(i)*A(j)
        enddo
        do j=1, 3
          C(i,j) = C(i,j) + fact4 * A(i)
          C(j,i) = C(j,i) + fact5 * A(i)
        enddo
      enddo

      return

      end






