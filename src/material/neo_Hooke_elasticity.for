

      subroutine neo_Hooke_elasticity(matData,F,stre,cc,finite)

      implicit none

      integer          finite,
     *                 NHtype, config, j, i, round
      double precision matData(*), F(3,*), stre(*), cc(6,*),
     *                 mu, K, Lamb, b(6), detF, C(6), invC(6),
     *                 fact, fact1, fact2, detC, detFm2d3,
     *                 det_u, det_s,
     *                 r1d2, r1, r2, r1d3, r2d3

      data r1d2, r1, r2 / .5d0, 1.d0, 2.d0 /
      r1d3 = r1 / 3.d0
      r2d3 = r1d3 + r1d3

      NHtype = max(1,round(matData(3)))
      config = 1    ! 1 -> current configuration
                    ! 2 -> initial configuration

      if (finite.lt.1) then
        call small_strain_elasticity(matData,F,stre,cc)
        return
      endif

      K      = matData(1)
      mu     = matData(2)
      Lamb   = K - r2d3 * mu

      if (NHtype.gt.5) 
     *  call prgerror(1,'neo_Hooke_elasticity','NHtype > 4 !')

      goto (1,2,3,4,5), NHtype

c========================================================================
 1    continue 
c
c     Psi = Lamb/4 (J^2-1) - Lamb/2 lnJ - mu lnJ + mu/2 (trb-3)
c
c     Wriggers, Nichtlin.FEM, page 45
c
      goto (11,12),config
c------------------------------------------------------------------------
 11   continue ! F -> sig, cc

      call mult_s_uuT_p(b,F,F)

      detF  = det_u(F)
      fact  = Lamb / (r2 * detF) * (detF*detF-r1) - mu / detF
      fact2 = mu / detF

c     Cauchy stress tensor sig

      stre(1) = fact + fact2 * b(1)
      stre(2) = fact + fact2 * b(2)
      stre(3) = fact + fact2 * b(3)
      stre(4) =        fact2 * b(4)
      stre(5) =        fact2 * b(5)
      stre(6) =        fact2 * b(6)

c     tangent tensor

      fact  = - fact
      fact1 =   fact * r2
      fact2 = Lamb * detF

      call pzero(cc,36)

      do i=1, 3
        do j=1, 3
          cc(i,j) = fact2
        enddo
        cc(i,i)     = cc(i,i) + fact1
        cc(i+3,i+3) =           fact
      enddo

      return
c------------------------------------------------------------------------

 12   continue ! F -> S, CC

      call mult_s_uTu_p(C,F,F)
      call inverse_s(invC,C)

      detC  = det_s(C)
      fact  = Lamb * r1d2 * (detC-r1) - mu

c     2. Piola-Kirchhoff stress tensor S

      Stre(1) = fact * invC(1) + mu
      Stre(2) = fact * invC(2) + mu
      Stre(3) = fact * invC(3) + mu
      Stre(4) = fact * invC(4)
      Stre(5) = fact * invC(5)
      Stre(6) = fact * invC(6)

c     tangent tensor

      fact = - fact - fact
      call diff_ATA_f(CC,invC,fact)
      fact = Lamb * detC
      call mult_4_ss_af(CC,invC,invC,fact)

      return
c========================================================================
 2    continue 
c
c     Psi = K/4 (J^2-1) - K/2 lnJ + mu/2 (J^(-2/3)trb-3)
c
c     Wriggers, Nichtlin.FEM, page 46
c
c     identical to Ogden with n      = 1, 
c                             alph_1 = 2, 
c                             mu_1   = mu
      goto (21,22),config

c------------------------------------------------------------------------
 21   continue ! F -> sig, cc

      call mult_s_uuT_p(b,F,F)

      detF     = det_u(F)
      detFm2d3 = detF**(-r2d3)
      fact     = K / (r2 * detF) * (detF*detF-r1)
     *            - mu / detF * detFm2d3 * r1d3 * (b(1)+b(2)+b(3))
      fact2    = mu / detF * detFm2d3

c     Cauchy stress tensor sig

      stre(1) = fact + fact2 * b(1)
      stre(2) = fact + fact2 * b(2)
      stre(3) = fact + fact2 * b(3)
      stre(4) =        fact2 * b(4)
      stre(5) =        fact2 * b(5)
      stre(6) =        fact2 * b(6)

c     tangent tensor

      fact  = - fact
      fact1 =   fact * r2
      fact2 = K * detF
     *          + r2d3 * r1d3 * mu / detF * detFm2d3 * (b(1)+b(2)+b(3))

      call pzero(cc,36)

      do i=1, 3
        do j=1, 3
          cc(i,j) = fact2
        enddo
        cc(i,i)     = cc(i,i) + fact1
        cc(i+3,i+3) =           fact
      enddo

      fact = - r2d3 * mu / detF * detFm2d3

      do i=1, 3
        do j=1, 6
          cc(i,j) = cc(i,j) + fact * b(j)
          cc(j,i) = cc(j,i) + fact * b(j)
        enddo
      enddo

      return
c------------------------------------------------------------------------
 22   continue ! F -> S, CC

      call mult_s_uTu_p(C,F,F)
      call inverse_s(invC,C)

      detC     = det_s(C)
      detFm2d3 = detC**(-r1d3)
      fact     = K * r1d2 * (detC-r1)
     *            - mu * detFm2d3 * r1d3 * (C(1)+C(2)+C(3))
      fact2    = mu * detFm2d3

c     2. Piola-Kirchhoff stress tensor S

      Stre(1) = fact * invC(1) + fact2
      Stre(2) = fact * invC(2) + fact2
      Stre(3) = fact * invC(3) + fact2
      Stre(4) = fact * invC(4)
      Stre(5) = fact * invC(5)
      Stre(6) = fact * invC(6)

c     tangent tensor

      fact = - fact - fact
      call diff_ATA_f(CC,invC,fact)
      fact = K * detC + r2d3 * r1d3 * detFm2d3*mu*(C(1)+C(2)+C(3))
      call mult_4_ss_af(CC,invC,invC,fact)
      fact = - r2d3 * mu * detFm2d3
      do i=1, 3
        do j=1, 6
          CC(i,j) = CC(i,j) + fact * invC(j)
          CC(j,i) = CC(j,i) + fact * invC(j)
        enddo
      enddo

      return
c========================================================================
 3    continue 
c
c     Psi = Lamb/2 (lnJ)^2 - mu lnJ + mu/2 (trb-3)
c
c     Reese, Habil, page 183
c
      goto (31,32),config

c------------------------------------------------------------------------
 31   continue ! F -> sig, cc

      call mult_s_uuT_p(b,F,F)

      detF  = det_u(F)
      fact  = Lamb / detF * log(detF) - mu / detF
      fact2 = mu / detF

c     Cauchy stress tensor sig

      stre(1) = fact + fact2 * b(1)
      stre(2) = fact + fact2 * b(2)
      stre(3) = fact + fact2 * b(3)
      stre(4) =        fact2 * b(4)
      stre(5) =        fact2 * b(5)
      stre(6) =        fact2 * b(6)
      
c      do i=1, 6
c        write(*,'(A, F8.3)') 'stre=', stre(i)
c      enddo

c     tangent tensor

      fact  = - fact
      fact1 =   fact * r2
      fact2 = Lamb / detF

      call pzero(cc,36)

      do i=1, 3
        do j=1, 3
          cc(i,j) = fact2
        enddo
        cc(i,i)     = cc(i,i) + fact1
        cc(i+3,i+3) =           fact
      enddo

      return
c------------------------------------------------------------------------
 32   continue ! F -> S, CC

      call mult_s_uTu_p(C,F,F)
      call inverse_s(invC,C)

      detF  = det_u(F)
      fact  = Lamb * log(detF) - mu

c     2. Piola-Kirchhoff stress tensor S

      Stre(1) = fact * invC(1) + mu
      Stre(2) = fact * invC(2) + mu
      Stre(3) = fact * invC(3) + mu
      Stre(4) = fact * invC(4)
      Stre(5) = fact * invC(5)
      Stre(6) = fact * invC(6)

c     tangent tensor

      fact = - fact - fact
      call diff_ATA_f(CC,invC,fact)
      call mult_4_ss_af(CC,invC,invC,Lamb)

      return
c========================================================================
 4    continue 
c
c     Psi = Lamb/4 (J^2-1) - Lamb/2 lnJ - mu lnJ + mu/2 (J^(-2/3)trb-3)

      goto (41,42),config

c------------------------------------------------------------------------
 41   continue ! F -> sig, cc

      call mult_s_uuT_p(b,F,F)

      detF     = det_u(F)
      detFm2d3 = detF**(-r2d3)
      fact     = Lamb / (r2 * detF) * (detF*detF-r1) - mu / detF 
     *               - mu / detF * detFm2d3 * r1d3 * (b(1)+b(2)+b(3))
      fact2    = mu / detF * detFm2d3

c     Cauchy stress tensor sig

      stre(1) = fact + fact2 * b(1)
      stre(2) = fact + fact2 * b(2)
      stre(3) = fact + fact2 * b(3)
      stre(4) =        fact2 * b(4)
      stre(5) =        fact2 * b(5)
      stre(6) =        fact2 * b(6)
      
c     tangent tensor

      fact  = - fact
      fact1 =   fact * r2
      fact2 = Lamb * detF 
     *         + r2d3 * r1d3 * mu / detF * detFm2d3 * (b(1)+b(2)+b(3))

      call pzero(cc,36)

      do i=1, 3
        do j=1, 3
          cc(i,j) = fact2
        enddo
        cc(i,i)     = cc(i,i) + fact1
        cc(i+3,i+3) =           fact
      enddo

      fact = - r2d3 * mu / detF * detFm2d3

      do i=1, 3
        do j=1, 6
          cc(i,j) = cc(i,j) + fact * b(j)
          cc(j,i) = cc(j,i) + fact * b(j)
        enddo
      enddo

      return
c------------------------------------------------------------------------
 42   continue ! F -> S, CC

      call mult_s_uTu_p(C,F,F)
      call inverse_s(invC,C)

      detC     = det_s(C)
      detFm2d3 = detC**(-r1d3)
      fact     = Lamb * r1d2 * (detC-r1) - mu
     *            - mu * detFm2d3 * r1d3 * (C(1)+C(2)+C(3))
      fact2    = mu * detFm2d3

c     2. Piola-Kirchhoff stress tensor S

      Stre(1) = fact * invC(1) + fact2
      Stre(2) = fact * invC(2) + fact2
      Stre(3) = fact * invC(3) + fact2
      Stre(4) = fact * invC(4)
      Stre(5) = fact * invC(5)
      Stre(6) = fact * invC(6)

c     tangent tensor

      fact = - fact - fact
      call diff_ATA_f(CC,invC,fact)
      fact = Lamb * detC + r2d3 * r1d3 * detFm2d3*mu*(C(1)+C(2)+C(3))
      call mult_4_ss_af(CC,invC,invC,fact)
      fact = - r2d3 * mu * detFm2d3
      do i=1, 3
        do j=1, 6
          CC(i,j) = CC(i,j) + fact * invC(j)
          CC(j,i) = CC(j,i) + fact * invC(j)
        enddo
      enddo

      return
c========================================================================

 5    continue 
c

c     Psi = Lamb/2 (lnJ)^2 - mu lnJ + mu/2 (J^(-2/3)trb-3)

      goto (51,52),config

c------------------------------------------------------------------------
 51   continue ! F -> sig, cc

      call mult_s_uuT_p(b,F,F)

      detF     = det_u(F)
      detFm2d3 = detF**(-r2d3)

      fact  = (Lamb * log(detF) - mu 
     *                 - mu * detFm2d3 * r1d3 * (b(1)+b(2)+b(3)) )/ detF

      fact2    = mu * detFm2d3 / detF

c          write(*,*) fact
c          write(*,*) fact2
c          write(*,*) b(1)

c     Cauchy stress tensor sig

      stre(1) = fact + fact2 * b(1)
      stre(2) = fact + fact2 * b(2)
      stre(3) = fact + fact2 * b(3)
      stre(4) =        fact2 * b(4)
      stre(5) =        fact2 * b(5)
      stre(6) =        fact2 * b(6)
      
c     tangent tensor

      fact  = - fact
      fact1 =   fact * r2
      fact2 = ( Lamb + r2d3 * r1d3 * mu * detFm2d3 * (b(1)+b(2)+b(3)) ) / detF

      call pzero(cc,36)

      do i=1, 3
        do j=1, 3
          cc(i,j) = fact2
        enddo
        cc(i,i)     = cc(i,i) + fact1
        cc(i+3,i+3) =           fact
      enddo

      fact = - r2d3 * mu * detFm2d3 / detF

      do i=1, 3
        do j=1, 6
          cc(i,j) = cc(i,j) + fact * b(j)
          cc(j,i) = cc(j,i) + fact * b(j)
        enddo
      enddo

      return
c------------------------------------------------------------------------
 52   continue ! F -> S, CC

      call mult_s_uTu_p(C,F,F)
      call inverse_s(invC,C)

      detF     = det_u(F)
      detC     = det_s(C)
      detFm2d3 = detF**(-r2d3)

      fact  = Lamb * log(detF) - mu - mu*detFm2d3*r1d3*(C(1)+C(2)+C(3))

      fact2 = mu * detFm2d3

c     2. Piola-Kirchhoff stress tensor S

      Stre(1) = fact * invC(1) + fact2
      Stre(2) = fact * invC(2) + fact2
      Stre(3) = fact * invC(3) + fact2
      Stre(4) = fact * invC(4)
      Stre(5) = fact * invC(5)
      Stre(6) = fact * invC(6)

c     tangent tensor

      fact = - fact - fact
      call diff_ATA_f(CC,invC,fact)
      fact = Lamb + r2d3 * r1d3 * detFm2d3*mu*(C(1)+C(2)+C(3))
      call mult_4_ss_af(CC,invC,invC,fact)
      fact = - r2d3 * mu * detFm2d3
      do i=1, 3
        do j=1, 6
          CC(i,j) = CC(i,j) + fact * invC(j)
          CC(j,i) = CC(j,i) + fact * invC(j)
        enddo
      enddo

      return
c========================================================================
      end





