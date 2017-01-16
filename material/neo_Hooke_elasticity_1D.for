
      subroutine neo_hooke_elasticity_1D(matDat,lam,stre,cc,dlam,
     *                                   finite,sss)
      implicit none

      integer          finite, sss

      double precision matDat(*), lam(*), stre(*), cc(2,*), dlam(*),  
     *                 E, nu, mu, K, Lamb, eps11, eps33, J, dJ, stre0,
     *                 lam0, 
     *                 fct, fct1, fct2, fct3, dJ1, dJ3

c
c     1D (TRUSS) MATERIAL SUBROUTINE
c
c     Neo-Hookean elasticity  (Wriggers, page 45)
c     -------------------------------------------
c
c     ss  = 1 -> uniaxial stress (plane stress), stre_22 = stre_33 = 0
c     ss  = 2 -> plane strain (membrane), stre_22 = 0, strain_33 = 0
c     ss  = 3 -> axisymmetric (membrane), stre_22 = 0
c    (ss  = 4 -> like 1, but incompressible
c     ss  = 5 -> like 2, but incompressible
c     ss  = 6 -> like 3, but incompressible)
c
c     stre: 11  33
c
c     cc:  1111  1133
c          3311  3333
c



      K      = matDat(1)
      mu     = matDat(2)
      Lamb   = K - 2.d0/3.d0 * mu

      if (finite.lt.1) goto 100

c...  finite strains .........

      if (sss.le.1) then

c       uniaxial stress

        fct     = 1.d0 / (lam(1) * lam(1))
        fct2    = mu / Lamb * fct
        fct3    = fct2 * fct2
        lam(2)  = sqrt( sqrt( fct3 + 2.d0*fct2 + fct) - fct2)
        lam(3)  = lam(2)
        dlam(2) = .5d0 / lam(2)
     *              * (.5d0 / sqrt( fct3 + 2.d0*fct2 + fct)
     *                  * (- fct3 - fct2 - .5d0 * fct) * 4.d0 / lam(1)
     *                     + 2.d0 * fct2 / lam(1))
        dlam(3) = dlam(2)
        J       = lam(1) * lam(2) * lam(3)
        dJ      = lam(2) * lam(3) + 2.d0 * lam(1) * lam(2) * dlam(2)
        stre(1) = .5d0 * Lamb / J * ( J    * J     - 1.d0)
     *                   + mu / J * (lam(1)*lam(1) - 1.d0)
        cc(1,1) = .5d0*Lamb * ((-1.d0/(J*J)*dJ)*(J*J-1.d0) + 2.d0*dJ)
     *                 + mu * ((-1.d0/(J*J)*dJ)*(lam(1)*lam(1)-1.d0)
     *                               + 2.d0/J*lam(1))

      elseif (sss.eq.2) then

c       plane strain
           
        lam(2)  = sqrt((Lamb+2.d0*mu)/(Lamb*lam(1)*lam(1)+2.d0*mu))
        lam(3)  = 1.d0
        dlam(2) = - lam(2)/(Lamb*lam(1)*lam(1)+2.d0*mu)*Lamb*lam(1)
        dlam(3) = 0.d0
        J       = lam(1) * lam(2)
        dJ      = lam(2) + lam(1) * dlam(2)
        stre(1) = .5d0 * Lamb / J * ( J    * J     - 1.d0)
     *                   + mu / J * (lam(1)*lam(1) - 1.d0)
        cc(1,1) = .5d0*Lamb * ((-1.d0/(J*J)*dJ)*(J*J-1.d0) + 2.d0*dJ)
     *                 + mu * ((-1.d0/(J*J)*dJ)*(lam(1)*lam(1)-1.d0)
     *                               + 2.d0/J*lam(1))

      elseif (sss.eq.3) then

c       axisymmetric

        fct     = Lamb * lam(1)*lam(1) * lam(3)*lam(3) + 2.d0*mu
        lam(2)  = sqrt( (Lamb+2.d0*mu) / fct)
        J       = lam(1) * lam(2) * lam(3)
        dlam(1) = - Lamb * J * lam(3) / fct     ! d lam2 / d lam1
        dlam(3) = dlam(1) / lam(3) * lam(1)     ! d lam2 / d lam3
        dJ1     = lam(2) * lam(3) + lam(1) * dlam(1) * lam(3)
        dJ3     = lam(1) * lam(2) + lam(1) * dlam(3) * lam(3)
        stre(1) = .5d0 * Lamb / J * ( J    * J     - 1.d0)
     *                   + mu / J * (lam(1)*lam(1) - 1.d0)

        stre(2) = .5d0 * Lamb / J * ( J    * J     - 1.d0)
     *                   + mu / J * (lam(3)*lam(3) - 1.d0)
        cc(1,1) = .5d0*Lamb * ((-1.d0/(J*J)*dJ1)*(J*J-1.d0) + 2.d0*dJ1)
     *                 + mu * ((-1.d0/(J*J)*dJ1)*(lam(1)*lam(1)-1.d0)
     *                               + 2.d0/J*lam(1))
        cc(1,2) = .5d0*Lamb * ((-1.d0/(J*J)*dJ3)*(J*J-1.d0) + 2.d0*dJ3)
     *                 + mu * ((-1.d0/(J*J)*dJ3)*(lam(1)*lam(1)-1.d0))
        cc(2,1) = .5d0*Lamb * ((-1.d0/(J*J)*dJ1)*(J*J-1.d0) + 2.d0*dJ1)
     *                 + mu * ((-1.d0/(J*J)*dJ1)*(lam(3)*lam(3)-1.d0))
        cc(2,2) = .5d0*Lamb * ((-1.d0/(J*J)*dJ3)*(J*J-1.d0) + 2.d0*dJ3)
     *                 + mu * ((-1.d0/(J*J)*dJ3)*(lam(3)*lam(3)-1.d0)
     *                               + 2.d0/J*lam(3))

      elseif (sss.eq.4) then

c       uniaxial stress, incompressible

        lam(2)  = 1.d0 / sqrt(lam(1))
        lam(3)  = lam(2)
        dlam(2) = - .5d0 / lam(1) * lam(2)
        dlam(3) = dlam(2)
        fct     = lam(1) * lam(1)
        stre(1) = mu * (fct*fct - 1.d0) / fct
        cc(1,1) = mu*(4.d0*lam(1)-2.d0*(fct*fct-1.d0)/(fct*lam(1)))

      elseif (sss.eq.5) then

c       plane strain, incompressible

        lam(2)  = 1.d0 / lam(1)
        lam(3)  = 1.d0
        dlam(2) = - 1.d0 / lam(1) * lam(2)
        dlam(3) = 0.d0
        fct     = lam(1) * lam(1)
        stre(1) = mu * (fct*fct - 1.d0) / fct 
        cc(1,1) = mu*(4.d0*lam(1)-2.d0*(fct*fct-1.d0)/(fct*lam(1)))

      elseif (sss.eq.6) then

c       axisymmetric, incompressible
           
        lam(2)  = 1.d0 / (lam(1) * lam(3))
        dlam(1) = - lam(2) / lam(1)   ! d lam2 / d lam1
        dlam(3) = - lam(2) / lam(3)   ! d lam2 / d lam3
        fct     = lam(1) * lam(1) * lam(3) * lam(3)
        stre(1) = mu * (fct * lam(1) * lam(1) - 1.d0) / fct 
        stre(2) = mu * (fct * lam(3) * lam(3) - 1.d0) / fct 
        cc(1,1) = mu * (2.d0*lam(1) + 2.d0/(fct*lam(1)))
        cc(1,2) = mu * (            + 2.d0/(fct*lam(3)))
        cc(2,1) = mu * (            + 2.d0/(fct*lam(1)))
        cc(2,2) = mu * (2.d0*lam(3) + 2.d0/(fct*lam(3)))

      else

        call prgError(1,'neo_Hooke_elasticity_1D','what do you want?')

      endif

      return

 100  continue

c...  small strains .........

      eps11 = lam(1) - 1.d0
      eps33 = lam(3) - 1.d0
      E     = 9.d0 * K * mu / (3.d0 * K + mu)
      nu    = Lamb / (2.d0 * (Lamb+mu))

      if (sss.le.1) then

c       uniaxial stress

        stre(1) = E * eps11
        cc(1,1) = E

      elseif (sss.eq.2) then

c       plane strain
           
        stre(1) = E / ((1.d0+nu) * (1.d0-nu)) * eps11
        cc(1,1) = E / ((1.d0+nu) * (1.d0-nu))

      elseif (sss.eq.3) then

c       axisymmetric

        fct1 = 4.d0 * mu * (mu + Lamb) / (2.d0 * mu + Lamb)
        fct2 = 2.d0 * mu * Lamb / (2.d0 * mu + Lamb)
        stre(1) = fct1 * eps11 + fct2 * eps33
        stre(2) = fct1 * eps33 + fct2 * eps11
        cc(1,1) = fct1
        cc(1,2) = fct2
        cc(2,1) = fct2
        cc(2,2) = fct1

      elseif (sss.eq.4) then

c       uniaxial stress, incompressible

        stre(1) = 4.d0 * mu * eps11
        cc(1,1) = 4.d0 * mu

      elseif (sss.eq.5) then

c       plane strain, incompressible

        stre(1) = 4.d0 * mu * eps11
        cc(1,1) = 4.d0 * mu

      elseif (sss.eq.6) then

c       axisymmetric, incompressible
           
        fct2 = mu + mu
        fct1 = fct1 + fct1
        stre(1) = fct1 * eps11 + fct2 * eps33
        stre(2) = fct2 * eps11 + fct1 * eps33
        cc(1,1) = fct1
        cc(1,2) = fct2
        cc(2,1) = fct2
        cc(2,2) = fct1

      else

        call prgError(2,'neo_Hooke_elasticity_1D','what do you want?')

      endif

      return

      end


