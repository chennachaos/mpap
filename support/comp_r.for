c
c
c
c    comp_R1111 -> R, dRdA, dRdf, dRdDg
c    comp_R1000 -> R
c    comp_R0111 ->    dRdA, dRdf, dRdDg
c    comp_R1100 -> R, dRdA
c    comp_R0100 ->    dRdA 
c
c
c==============================================================================

      subroutine comp_R1111(R,dRdA,dRdf,dRdDg,f,A,invA,diA,An,invAn,Dg,
     *                      nip)
c
c     compute R, dRdA, dRdf, dRdDg
c
c     nip = 0 -> simple backward Euler
c         > 0 -> backward Euler exponential map with nip series terms
c
      implicit none
      integer          nip, i
      double precision R(6), dRdf(6,6), dRdA(6,6), dRdDg(6), 
     *                 f(6), A(6), invA(6), diA(6,6), An(6), invAn(6), 
     *                 Dg, dRdiA(6,6), hs(6), hs2(6), hu(3,3), hu2(3,3), 
     *                 h4(6,6), h42(6,6), h4df(6,6), h4diA(6,6), 
     *                 fact, fact2

      if (nip.eq.0) then

c       simple backward Euler
c    ------------------------------
        fact = - Dg
        call set_s_p(R,A)        ! R
        call set_s_am(R,An)      !
        call set_s_af(R,f,fact)  !
        call unit_4_f(dRdf,fact) ! dRdf
        call unit_4_p(dRdA)      ! dRdA
        call set_s_m(dRdDg,f)    ! dRdDg

      else

c       backward Euler exponential map
c    -------------------------------------
c       for i=1

        call set_s_p(R,A)
        call mult_u_ss_p(hu,invAn,A)
        call mult_s_su_am(R,A,hu)
        call set_s_p(hs,f)
        call set_s_af(R,hs,Dg)
        
        call unit_4_p(dRdA)
        call diff_TAT_am(dRdA,hu)

        call unit_4_f(dRdf,Dg)         

        call set_s_p(dRdDg,f)
        
c       for i=2,3,..nip
        if (nip.gt.1) then

          call mult_u_ss_p(hu,invA,f)
          call diff_TU_p(h4,hu)
          call pzero(dRdiA,36)
          call set_4_p(h4df,h4)

c         for i=2

          fact  = Dg * Dg * 0.5d0

c         dR/df
          call mult_u_ss_p(hu2,hs,invA)
          call diff_UT_ap(h4df,hu2)
          call set_4_af(dRdf,h4df,fact)

c         dR/dinvA
          call diff_ATB_p(h4diA,hs,f)
          call set_4_af(dRdiA,h4diA,fact)

c         R
          call mult_s_su_p(hs2,hs,hu)
          call set_s_p(hs,hs2)
          call set_s_af(R,hs,fact)

c         dR/dDg
          call set_s_af(dRdDg,hs,Dg)
          fact2 = fact

c         for i=3,4...nip
          do i=3, nip

            fact  = fact * Dg / real(i)

c           dR/df
            call mult_4_44_p(h42,h4,h4df)
            call set_4_p(h4df,h42)
            call mult_u_ss_p(hu2,hs,invA)
            call diff_UT_ap(h4df,hu2)
            call set_4_af(dRdf,h4df,fact)

c           dR/dinvA
            call mult_4_44_p(h42,h4,h4diA)
            call set_4_p(h4diA,h42)
            call diff_ATB_ap(h4diA,hs,f)
            call set_4_af(dRdiA,h4diA,fact)

c           R
            call mult_s_su_p(hs2,hs,hu)
            call set_s_p(hs,hs2)
            call set_s_af(R,hs,fact)

c           dR/dDg
            call set_s_af(dRdDg,hs,fact2)
            fact2 = fact2 * Dg / real(i)

          enddo

c         dR/dA = dR/dinvA : dinvA/dA
          call mult_4_44_ap(dRdA,dRdiA,diA)

        endif

      endif

      return

      end


c=============================================================================

      subroutine comp_R1000(R,f,A,invA,An,invAn,Dg,nip)
c
c     compute R
c
c     nip = 0 -> simple backward Euler
c         > 0 -> backward Euler exponential map with nip series terms
c
      implicit none
      integer          nip, i
      double precision R(6), 
     *                 f(6), A(6), invA(6), An(6), invAn(6), Dg,
     *                 hs(6), hs2(6), hu(3,3), fact

      if (nip.eq.0) then

c       simple backward Euler
c    ------------------------------
        fact = - Dg
        call set_s_p(R,A)        ! R
        call set_s_am(R,An)      !
        call set_s_af(R,f,fact)  !

      else

c       backward Euler exponential map
c    -------------------------------------
c       for i=1

        call set_s_p(R,A)
        call mult_u_ss_p(hu,invAn,A)
        call mult_s_su_am(R,A,hu)
        call set_s_p(hs,f)
        call set_s_af(R,hs,Dg)
        
c       for i=2,3,..nip
        if (nip.gt.1) then

          call mult_u_ss_p(hu,invA,f)

          fact  = Dg

          do i=2, nip

            fact  = fact * Dg / real(i)

c           R
            call mult_s_su_p(hs2,hs,hu)
            call set_s_p(hs,hs2)
            call set_s_af(R,hs,fact)

          enddo

        endif

      endif

      return

      end



c=============================================================================



      subroutine comp_R0111(dRdA,dRdf,dRdDg,f,A,invA,diA,An,invAn,Dg,
     *                      nip)
c
c     compute dRdA, dRdf, dRdDg
c
c     nip = 0 -> simple backward Euler
c         > 0 -> backward Euler exponential map with nip series terms
c
      implicit none
      integer          nip, i
      double precision dRdf(6,6), dRdA(6,6), dRdDg(6), 
     *                 f(6), A(6), invA(6), diA(6,6), An(6), invAn(6), 
     *                 Dg, dRdiA(6,6), hs(6), hs2(6), hu(3,3), hu2(3,3), 
     *                 h4(6,6), h42(6,6), h4df(6,6), h4diA(6,6), 
     *                 fact, fact2

      if (nip.eq.0) then

c       simple backward Euler
c    ------------------------------
        fact = - Dg
        call unit_4_f(dRdf,fact) ! dRdf
        call unit_4_p(dRdA)      ! dRdA
        call set_s_m(dRdDg,f)    ! dRdDg

      else

c       backward Euler exponential map
c    -------------------------------------
c       for i=1

        call mult_u_ss_p(hu,invAn,A)
        call set_s_p(hs,f)
        
        call unit_4_p(dRdA)
        call diff_TAT_am(dRdA,hu)

        call unit_4_f(dRdf,Dg)         

        call set_s_p(dRdDg,f)
        
c       for i=2,3,..nip
        if (nip.gt.1) then

          call mult_u_ss_p(hu,invA,f)
          call diff_TU_p(h4,hu)
          call pzero(dRdiA,36)
          call set_4_p(h4df,h4)

c         for i=2

          fact  = Dg * Dg * 0.5d0

c         dR/df
          call mult_u_ss_p(hu2,hs,invA)
          call diff_UT_ap(h4df,hu2)
          call set_4_af(dRdf,h4df,fact)

c         dR/dinvA
          call diff_ATB_p(h4diA,hs,f)
          call set_4_af(dRdiA,h4diA,fact)

          call mult_s_su_p(hs2,hs,hu)
          call set_s_p(hs,hs2)

c         dR/dDg
          call set_s_af(dRdDg,hs,Dg)
          fact2 = fact

c         for i=3,4...nip
          do i=3, nip

            fact  = fact * Dg / real(i)

c           dR/df
            call mult_4_44_p(h42,h4,h4df)
            call set_4_p(h4df,h42)
            call mult_u_ss_p(hu2,hs,invA)
            call diff_UT_ap(h4df,hu2)
            call set_4_af(dRdf,h4df,fact)

c           dR/dinvA
            call mult_4_44_p(h42,h4,h4diA)
            call set_4_p(h4diA,h42)
            call diff_ATB_ap(h4diA,hs,f)
            call set_4_af(dRdiA,h4diA,fact)

            call mult_s_su_p(hs2,hs,hu)
            call set_s_p(hs,hs2)

c           dR/dDg
            call set_s_af(dRdDg,hs,fact2)
            fact2 = fact2 * Dg / real(i)

          enddo

c         dR/dA = dR/dinvA : dinvA/dA
          call mult_4_44_ap(dRdA,dRdiA,diA)

        endif

      endif

      return

      end


c=============================================================================


      subroutine comp_R1100(R,dRdA,f,A,invA,diA,An,invAn,Dg,nip)
c
c     compute R, dRdA
c
c     nip = 0 -> simple backward Euler
c         > 0 -> backward Euler exponential map with nip series terms
c
      implicit none
      integer          nip, i
      double precision R(6), dRdA(6,6), 
     *                 f(6), A(6), invA(6), diA(6,6), An(6), invAn(6), 
     *                 Dg, dRdiA(6,6), hs(6), hs2(6), hu(3,3), hu2(3,3), 
     *                 h4(6,6), h42(6,6), h4diA(6,6), 
     *                 fact

      if (nip.eq.0) then

c       simple backward Euler
c    ------------------------------
        fact = - Dg
        call set_s_p(R,A)        ! R
        call set_s_am(R,An)      !
        call set_s_af(R,f,fact)  !
        call unit_4_p(dRdA)      ! dRdA

      else

c       backward Euler exponential map
c    -------------------------------------
c       for i=1

        call set_s_p(R,A)
        call mult_u_ss_p(hu,invAn,A)
        call mult_s_su_am(R,A,hu)
        call set_s_p(hs,f)
        call set_s_af(R,hs,Dg)
        
        call unit_4_p(dRdA)
        call diff_TAT_am(dRdA,hu)

c       for i=2,3,..nip
        if (nip.gt.1) then

          call mult_u_ss_p(hu,invA,f)
          call diff_TU_p(h4,hu)
          call pzero(dRdiA,36)

c         for i=2

          fact  = Dg * Dg * 0.5d0

c         dR/dinvA
          call diff_ATB_p(h4diA,hs,f)
          call set_4_af(dRdiA,h4diA,fact)

c         R
          call mult_s_su_p(hs2,hs,hu)
          call set_s_p(hs,hs2)
          call set_s_af(R,hs,fact)

c         for i=3,4...nip
          do i=3, nip

            fact  = fact * Dg / real(i)

c           dR/dinvA
            call mult_4_44_p(h42,h4,h4diA)
            call set_4_p(h4diA,h42)
            call diff_ATB_ap(h4diA,hs,f)
            call set_4_af(dRdiA,h4diA,fact)

c           R
            call mult_s_su_p(hs2,hs,hu)
            call set_s_p(hs,hs2)
            call set_s_af(R,hs,fact)

          enddo

c         dR/dA = dR/dinvA : dinvA/dA
          call mult_4_44_ap(dRdA,dRdiA,diA)

        endif

      endif

      return

      end


c=============================================================================

      subroutine comp_R0100(dRdA,f,A,invA,diA,An,invAn,Dg,nip)
c
c     compute R, dRdA
c
c     nip = 0 -> simple backward Euler
c         > 0 -> backward Euler exponential map with nip series terms
c
      implicit none
      integer          nip, i
      double precision dRdA(6,6), 
     *                 f(6), A(6), invA(6), diA(6,6), An(6), invAn(6), 
     *                 Dg, dRdiA(6,6), hs(6), hs2(6), hu(3,3), hu2(3,3), 
     *                 h4(6,6), h42(6,6), h4diA(6,6), 
     *                 fact

      if (nip.eq.0) then

c       simple backward Euler
c    ------------------------------
        call unit_4_p(dRdA)      ! dRdA

      else

c       backward Euler exponential map
c    -------------------------------------
c       for i=1

        call mult_u_ss_p(hu,invAn,A)
        call set_s_p(hs,f)
        
        call unit_4_p(dRdA)
        call diff_TAT_am(dRdA,hu)

c       for i=2,3,..nip
        if (nip.gt.1) then

          call mult_u_ss_p(hu,invA,f)
          call diff_TU_p(h4,hu)
          call pzero(dRdiA,36)

c         for i=2

          fact  = Dg * Dg * 0.5d0

c         dR/dinvA
          call diff_ATB_p(h4diA,hs,f)
          call set_4_af(dRdiA,h4diA,fact)

          call mult_s_su_p(hs2,hs,hu)
          call set_s_p(hs,hs2)

c         for i=3,4...nip
          do i=3, nip

            fact  = fact * Dg / real(i)

c           dR/dinvA
            call mult_4_44_p(h42,h4,h4diA)
            call set_4_p(h4diA,h42)
            call diff_ATB_ap(h4diA,hs,f)
            call set_4_af(dRdiA,h4diA,fact)

            call mult_s_su_p(hs2,hs,hu)
            call set_s_p(hs,hs2)

          enddo

c         dR/dA = dR/dinvA : dinvA/dA
          call mult_4_44_ap(dRdA,dRdiA,diA)

        endif

      endif

      return

      end


c=============================================================================
