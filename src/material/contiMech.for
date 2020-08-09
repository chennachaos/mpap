c
c  subroutine calcC(C,F)
c  subroutine calcB(B,F)
c
c  subroutine pushfwrd_s(B,F,C)
c  subroutine pushfwrd2D_s(B,F,F33,C)
c
c  subroutine pushfwrd_4(cc,F,C)
c  subroutine pushfwrd2D_4(cc,F,F33,C)
c
c  subroutine taub2sigccPart1(stre,cc,b,dtaudb,detF)
c  subroutine taub2sigccPart2(stre,cc)
c
c  subroutine computeElasticResponse(S,dSdC,dSdCp,C,invC,diCdC,Cp,invCp,diCpdCp,d)
c
c  subroutine smallStrainTensor(eps,F)
c
c
c
c

c=============================================================================
c
c     CALCULATE THE RIGHT CAUCHY-GREEN TENSOR
c
      subroutine calcC(C,F)
c-----------------------------------------------------------------------------
      implicit none
      double precision C(*), F(3,*)

      C(1) = F(1,1)*F(1,1) + F(2,1)*F(2,1) + F(3,1)*F(3,1)
      C(2) = F(1,2)*F(1,2) + F(2,2)*F(2,2) + F(3,2)*F(3,2)
      C(3) = F(1,3)*F(1,3) + F(2,3)*F(2,3) + F(3,3)*F(3,3)
      C(4) = F(1,1)*F(1,2) + F(2,1)*F(2,2) + F(3,1)*F(3,2)
      C(5) = F(1,2)*F(1,3) + F(2,2)*F(2,3) + F(3,2)*F(3,3)
      C(6) = F(1,3)*F(1,1) + F(2,3)*F(2,1) + F(3,3)*F(3,1)

      return
      end

c=============================================================================
c
c     CALCULATE THE LEFT CAUCHY-GREEN TENSOR
c
      subroutine calcB(b,F)
c-----------------------------------------------------------------------------
      implicit none
      double precision b(*), F(3,*)

      b(1) = F(1,1)*F(1,1) + F(1,2)*F(1,2) + F(1,3)*F(1,3)
      b(2) = F(2,1)*F(2,1) + F(2,2)*F(2,2) + F(2,3)*F(2,3)
      b(3) = F(3,1)*F(3,1) + F(3,2)*F(3,2) + F(3,3)*F(3,3)
      b(4) = F(1,1)*F(2,1) + F(1,2)*F(2,2) + F(1,3)*F(2,3)
      b(5) = F(2,1)*F(3,1) + F(2,2)*F(3,2) + F(2,3)*F(3,3)
      b(6) = F(3,1)*F(1,1) + F(3,2)*F(1,2) + F(3,3)*F(1,3)

      return
      end

c=============================================================================
c
c     PUSHFORWARD OF SYMMETRIC SECOND ORDER TENSOR
c
      subroutine pushfwrd_s(B,F,C)
c-----------------------------------------------------------------------------
      implicit none
      double precision B(*), F(3,*), C(*), H(3,3)

      call mult_u_us_p (H,F,C)

      call mult_s_uuT_p(B,H,F)

      return
      end

c=============================================================================
c
c     PUSHFORWARD OF SYMMETRIC SECOND ORDER TENSOR IN 2D
c
      subroutine pushfwrd2D_s(B,F,F33,C)
c-----------------------------------------------------------------------------
      implicit none
      double precision B(*), F(2,*), F33, C(*), H(2,2)

      H(1,1) = F(1,1)*C(1) + F(1,2)*C(3)
      H(1,2) = F(1,1)*C(3) + F(1,2)*C(2)
      H(2,1) = F(2,1)*C(1) + F(2,2)*C(3)
      H(2,2) = F(2,1)*C(3) + F(2,2)*C(2)
      B(1)   = H(1,1)*F(1,1) + H(1,2)*F(1,2)
      B(2)   = H(2,1)*F(2,1) + H(2,2)*F(2,2)
      B(3)   = H(1,1)*F(2,1) + H(1,2)*F(2,2)
      B(4)   = C(4) * F33 * F33
      return
      end

c=============================================================================
c
c     PUSHFORWARD OF FOURTH ORDER TENSOR (format: 6x6 matrix)
c
      subroutine pushfwrd_4(cc,F,C)
c-----------------------------------------------------------------------------
      implicit none
      double precision cc(6,*), F(3,*), C(6,*), fact, H(6,6)

      integer I, J, K, L, ii, jj, kk, ll, 
     *        indx1(2,6), indx2(3,3), A, B, rw, cl

      data indx1 / 1, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3, 1 /

      data indx2 / 1, 4, 6, 
     *             4, 2, 5,
     *             6, 5, 3 /
c
c  this subroutine performs 756 multiplications !!!!!!  avoid using it !!!
c
      call pzero(H,36)
      do rw=1, 6
        ii = indx1(1,rw) 
        jj = indx1(2,rw)
        do I=1, 3
          do J=1, 3
            A = indx2(I,J)
            fact = F(ii,I) * F(jj,J) ! x 54
            do cl=1, 6
              H(rw,cl) = H(rw,cl) + fact * C(A,cl) ! x 324
            enddo
          enddo
        enddo
      enddo

      call pzero(cc,36)
      do cl=1, 6
        kk = indx1(1,cl) 
        ll = indx1(2,cl)
        do K=1, 3
          do L=1, 3
            B = indx2(K,L)
            fact = F(kk,K) * F(ll,L) ! x 54
            do rw=1, 6
              cc(rw,cl) = cc(rw,cl) + fact * H(rw,B) ! x 324
            enddo
          enddo
        enddo
      enddo

      return
      end

c=============================================================================
c
c     PUSHFORWARD OF FOURTH ORDER TENSOR IN 2D (format: 4x4 matrix)
c
      subroutine pushfwrd2D_4(cc,F,F33,C)
c-----------------------------------------------------------------------------
      implicit none
      double precision cc(4,*), F(2,*), F33, C(4,*), 
     *                 fact, H(3,3), fact1, F332

      integer          I, J, K, L, ii, jj, kk, ll, 
     *                 indx1(2,3), indx2(2,2), A, B, rw, cl

      data indx1 / 1, 1, 2, 2, 1, 2 /

      data indx2 / 1, 3,
     *             3, 2 /

c
c  this subroutine performs 159 multiplications !!!!!!  avoid using it !!!
c
      call pzero(H,9)
      do rw=1, 3
        ii = indx1(1,rw) 
        jj = indx1(2,rw)
        do I=1, 2
          do J=1, 2
            A = indx2(I,J)
            fact = F(ii,I) * F(jj,J) ! x 12
            do cl=1, 3
              H(rw,cl) = H(rw,cl) + fact * C(A,cl) ! x 36
            enddo
          enddo
        enddo
      enddo

      call pzero(cc,16)
      do cl=1, 3
        kk = indx1(1,cl) 
        ll = indx1(2,cl)
        do K=1, 2
          do L=1, 2
            B = indx2(K,L)
            fact = F(kk,K) * F(ll,L) ! x 12
            do rw=1, 3
              cc(rw,cl) = cc(rw,cl) + fact * H(rw,B) ! x 36
            enddo
          enddo
        enddo
      enddo

      F332 = F33 * F33 ! x 1

      do rw=1, 3
        ii = indx1(1,rw) 
        jj = indx1(2,rw)
        do I=1, 2
          fact1 = F332 * F(ii,I) ! x 6
          do J=1, 2
            A = indx2(I,J)
            cc(rw,4) = cc(rw,4) + fact1 * F(jj,J) * C(A,4) ! x 12 + 12
          enddo
        enddo
      enddo

      do cl=1, 3
        kk = indx1(1,cl) 
        ll = indx1(2,cl)
        do K=1, 2
          fact1 = F332 * F(kk,K) ! x 6
          do L=1, 2
            B = indx2(K,L)
c           if (cl.eq.3) write(*,'(10(i3))') cl, kk, K, ll, L, B
            cc(4,cl) = cc(4,cl) + fact1 * F(ll,L) * C(4,B) ! x 12 + 12
          enddo
        enddo
      enddo

      cc(4,4) = C(4,4) * F332 * F332 ! x 2

      return
      end
c=============================================================================
c
c     CALCULATION OF CAUCHY STRESS AND SPATIAL TANGENT
c
c     input:  tau    - Kirchhoff stress
c             b      - left Cauchy-Green tensor b = F Ft
c             dtaudb - derivative of tau with respect to b
c             detF   - determinant of F
c
c     output: sig    - Cauchy stress sigma
c             cc     - spatial tangent tensor !!without geometrical part!!
c
      subroutine taub2sigccPart1(sig,cc,tau,dtaudb,b,detF)
c-----------------------------------------------------------------------------
      implicit none

      double precision sig(6), cc(6,6), tau(6), dtaudb(6,6), b(6), detF,
     *                 fact, r1d2

      r1d2 = .5d0

      fact = 1.d0 / detF

      call set_s_af(sig,tau,fact)

      fact = fact + fact

      call mult_4_4s_af(cc,dtaudb,b,fact)

      return

      end

c=============================================================================
c
c     CALCULATION OF CAUCHY STRESS AND SPATIAL TANGENT
c
c     input:  sig - Cauchy stress
c             cc  - spatial tangent without geometrical part
c
c     output: cc  - complete spatial tangent tensor
c
      subroutine taub2sigccPart2(sig,cc)
c-----------------------------------------------------------------------------
      implicit none

      double precision sig(6), cc(6,6),
     *                 r1d2

      r1d2 = .5d0

      cc(1,1) = cc(1,1) - sig(1) - sig(1)
      cc(1,4) = cc(1,4) - sig(4)
      cc(1,6) = cc(1,6) - sig(6)
      cc(2,2) = cc(2,2) - sig(2) - sig(2)
      cc(2,4) = cc(2,4) - sig(4)
      cc(2,5) = cc(2,5) - sig(5)
      cc(3,3) = cc(3,3) - sig(3) - sig(3)
      cc(3,5) = cc(3,5) - sig(5)
      cc(3,6) = cc(3,6) - sig(6)
      cc(4,1) = cc(4,1) - sig(4)
      cc(4,2) = cc(4,2) - sig(4)
      cc(4,4) = cc(4,4) -(sig(1)+sig(2)) * r1d2
      cc(4,5) = cc(4,5) - sig(6) * r1d2
      cc(4,6) = cc(4,6) - sig(5) * r1d2
      cc(5,2) = cc(5,2) - sig(5)
      cc(5,3) = cc(5,3) - sig(5)
      cc(5,4) = cc(5,4) - sig(6) * r1d2
      cc(5,5) = cc(5,5) -(sig(2)+sig(3)) * r1d2
      cc(5,6) = cc(5,6) - sig(4) * r1d2
      cc(6,1) = cc(6,1) - sig(6)
      cc(6,3) = cc(6,3) - sig(6)
      cc(6,4) = cc(6,4) - sig(5) * r1d2
      cc(6,5) = cc(6,5) - sig(4) * r1d2
      cc(6,6) = cc(6,6) -(sig(3)+sig(1)) * r1d2

      return

      end

c=============================================================================
c
c     COMPUTATION OF ELASTIC MATERIAL RESPONSE
c     BASED ON MULTIPLICATIVE SPLIT OF DEFORMATION GRADIENT 
c     INTO ELASTIC AND INELASTIC PARTS
c
c     F = Fe Fp
c
      subroutine computeElasticResponse(S,dSdC,dSdCp,
     *                                  C,invC,diCdC,Cp,invCp,diCpdCp,
     *                                  K,matData,typ)
c
c     input: tensors C, invC, dinvC/dC, Cp, invCp, dinvCp/dCp
c            bulk modulus K, shear stiffness in matData, elasticity typ
c
c     output: S, dSdC, dSdCp
c
c-----------------------------------------------------------------------------
      implicit none
      integer          typ,
     *                 i, l, round
      double precision S(6), dSdC(6,6), dSdCp(6,6), C(6), invC(6),
     *                 diCdC(6,6), Cp(6), invCp(6), diCpdCp(6,6),
     *                 K, matData(*),
     *                 h4(6,6), h2(9), fact, srC(6), isrC(6), dsrC(6,6),
     *                 disrC(6,6), h22(9), h23(9), h42(6,6), h43(6,6),
     *                 mu, Lamb, J2, J, muJm2d3, p, pJ, trp,
     *                 det_s, traceAB_ss,
     *                 r1d3, r1d2, r2d3, r5d6, r1, r4d3

      data r1d2, r1 / 0.5d0, 1.d0 /
      r1d3 = 1.d0 / 3.d0
      r2d3 = r1d3 + r1d3
      r4d3 = r2d3 + r2d3
      r5d6 = 2.5d0 * r1d3

      goto (10,20,30,40), typ

c-----------------------------------------------------------------------------
c     (1) Neo-Hooke
c
c    sig = mu Je^(-5/3) dev(Be) + K (Je-1)
c
c   -> S = mu Je^(-2/3) (invCp - tr(C invCp)/3 invC) + K Je(Je-1) invC
c
c    Je = sqrt(det(C)/det(Cp))
c
c   Wriggers, page 46  XXXXXXXXXXXXXXXXXXXXX???????XXXXXXXXX
c
c-----------------------------------------------------------------------------
 10   continue

      mu      = matData(1)

      J2      = det_s(C) / det_s(Cp)
      muJm2d3 = J2**(-r1d3) * mu
      pJ      = r1d2 * K * (J2 - r1)
      trp     = - muJm2d3 * r1d3 * traceAB_ss(C,invCp) + pJ

c     S
      call set_s_f     (S,invCp,muJm2d3)
      call set_s_af    (S,invC, trp)

c      do i=1, 3
c        S(i) = S(i) + muJm2d3
c      enddo
      

c     dSdC
      fact = - r1d3
      call set_s_f     (h2,S,fact) 
      fact = r2d3 * K * J2 - K / 6.d0
      call set_s_af  (h2,invC,fact)
      call mult_4_ss_p (dSdC,h2,invC)
      fact = - r1d3 * muJm2d3
      call mult_4_ss_af(dSdC,invC,invCp,fact)
      call set_4_af    (dSdC,diCdC,trp)

c     dSdCp
      call mult_4_ss_m (dSdCp,h2,invCp)
      call mult_4_ss_f (h4,invC,C,fact)
      call unit_4_af   (h4,muJm2d3)
      call mult_4_44_ap(dSdCp,h4,diCpdCp)

      return

c-----------------------------------------------------------------------------
c     (2) St. Venant-Kirchhoff
c-----------------------------------------------------------------------------
 20   continue

      mu      = matData(1)
      Lamb    = K - r2d3 * mu

      call mult_u_ss_p (h2,invCp,C)
      call mult_s_us_f (S,h2,invCp,mu)
      fact = - mu
      call set_s_af    (S,invCp,fact)
      fact = Lamb * 0.5d0 * (h2(1)+h2(5)+h2(9)-3.d0)
      call set_s_af    (S,invC,fact)

      call diff_ATA_f  (dSdC,invCp,mu)
      call set_4_af    (dSdC,diCdC,fact)
      fact = Lamb * 0.5d0
      call mult_4_ss_af(dSdC,invC,invCp,fact)        

      call diff_UT_f   (h4,h2,mu)
      call transpose_r (h2)
      call diff_TU_af  (h4,h2,mu)
      fact = Lamb * 0.5d0
      call mult_4_ss_af(h4,invC,C,fact)
      call mult_4_44_p (dSdCp,h4,diCpdCp)
      fact = - mu
      call set_4_af    (dSdCp,diCpdCp,fact)

      return

c-----------------------------------------------------------------------------
c     (3) Hencky
c-----------------------------------------------------------------------------
 30   continue

c      write(*,*) "   Warning! this Hencky stuff may be dodgy!"

      mu = matData(1)

      call sqrt_s      (srC,dsrC,C)    ! -> srC, dsrC
      call inverse_s   (isrC,srC)      ! -> isrC
      call diff_ATA_m  (h4,isrC)       !
      call mult_4_44_p (disrC,h4,dsrC) ! -> disrC
      call inverse_s   (invCp,Cp)      ! -> iCp

      !call mult_s_ss_p (h2,isrC,isrC)
      !call set_s_am    (h2,invC)
      
      call mult_u_ss_p (h2,invCp,srC)  !
      call mult_s_su_p (h22,srC,h2)    !
      call log_s       (h23,h4,h22,10) ! -> h23 = log(srC invCp srC)
                                       ! -> h4  = d log(") d "
      call diff_TAT_p  (h42,h2)
      call mult_4_44_p (h43,h4,h42)    
      call mult_4_44_p (h42,h43,dsrC)  ! -> h42 = d log(srC invCp srC) / d C

      call mult_u_ss_p (h2,h23,isrC)   ! -> h2  = log(srC invCp srC) isrC

      call mult_s_su_f (S,isrC,h2,mu)

      call diff_TAT_p  (h43,h2)
      call mult_4_44_f (dSdC,h43,disrC,mu)
      call diff_ATA_p  (h43,isrC)      ! -> h43 = d isrC T isrC / d T
      call mult_4_44_af(dSdC,h43,h42,mu)

      fact = - mu / 3.d0 * (h23(1)+h23(2)+h23(3))
     *        + K / 2.d0 * (det_s(C) - 1.d0)
      call set_s_af(S,invC,fact)

      call set_4_af    (dSdC,diCdC,fact)
      call unit_s_p    (h23)
      fact = - mu / 3.d0
      call mult_s_s4_f (h22,h23,h42,fact)
      fact =    K / 2.d0 * det_s(C)
      call set_s_af    (h22,invC,fact)
      call mult_4_ss_ap(dSdC,invC,h22)

      call set_4_f     (h42,h43,mu)
      call unit_s_p    (h2)
      fact = - mu / 3.d0
      call mult_4_ss_af(h42,invC,h2,fact)
      call mult_4_44_p (h43,h42,h4)
      call diff_ATA_p  (h42,srC)
      call mult_4_44_p (h4,h43,h42)
      call mult_4_44_p (dSdCp,h4,diCpdCp)

      return

c-----------------------------------------------------------------------------
c     (4) Ogden
c-----------------------------------------------------------------------------
 40   continue

      call prgError(1,"computeElasticResponse(contiMech.for)",
     *                "this has not yet been implemented!")

      return

      end
c=============================================================================
c
c     CALCULATION OF SMALL STRAIN TENSOR FROM DEFORMATION GRADIENT
c
      subroutine smallStrainTensor(eps,F)
c-----------------------------------------------------------------------------
      implicit none
      double precision eps(*), F(3,*), r1d2, r1

      data r1d2, r1 / 0.5d0, 1.d0 /

      eps(1) = F(1,1) - r1
      eps(2) = F(2,2) - r1
      eps(3) = F(3,3) - r1

      eps(4) = r1d2 * (F(1,2)+F(2,1))
      eps(5) = r1d2 * (F(2,3)+F(3,2))
      eps(6) = r1d2 * (F(3,1)+F(1,3))

      return

      end
c=============================================================================


