
      subroutine matlib1d(matDat,lam,stre,cc,dlam,iv1,iv2,dt,
     *                    matId,nivGp,finite,sss,isw,err,gp,ELM)

      implicit none

      integer          ELM ! do not touch

      integer          matId, nivGp, finite, sss, isw, gp,
     *                 err,
     *                 iter, ndm, matDim

      double precision matDat(*), lam(*), stre(*), cc(2,*), dlam(*), 
     *                 iv1, iv2, dt,
     *                 F3D(3,3), stre3D(6), stre0, cc3D(6,6), K, dmy(10)

c
c   isw
c    1   set nivGp and exit
c    2   initialise internal variables
c    3   full stress update (compute stress, internal variables and tangent tensor)
c

      if (isw.lt.1 .or. isw.gt.3) 
     *  call prgError(1,'matlib1d','invalid isw!')

      ndm = matDim(matId)

      if (ndm.eq.2) call prgError(2,'matlib1d',
     * 'use 3D material! no interface to 2D materials implemented yet!')


      if (ndm.eq.3 .and. isw.eq.3) goto 500

c.... EVERYTHING FOR: (isw=1,2,3 with ndm=1)  and  (isw=1,2 with ndm=3) ........................

          ! 1 2 3 4 5 6  7 8 9 10 11 12
      goto (1,1,1,1,1,1,99,1,1,10, 1, 1), matId

 99   call prgError(3,'matlib1d','matId / ndm mismatch ?!')

c.........................................................................................

 1    continue ! 3D material subroutines

      call matlib3d(dmy,dmy,dmy,dmy,iv1,iv2,dt,
     *              matId,nivGp,finite,isw,err,gp,ELM)

      return

c.........................................................................................

 10   continue  ! 1D Neo Hooke elasticity

      nivGp = 0 

      if (isw.ne.3) return

      call Neo_Hooke_Elasticity_1D(matDat,lam,stre,cc,dlam,finite,sss)

      return

c.........................................................................................

 500  continue

c.... FULL STRESS UPDATE BASED ON 3D MATERIAL SUBROUTINE: (isw=3 with ndm=3) .....................

      if (sss.le.1) then

c...    uniaxial stress (plane stress, truss) ......

        call pzero(F3D,9)
        lam(2) = lam(1)**(-1.d0/3.d0)
        F3D(1,1) = lam(1)
        F3D(2,2) = lam(2)
        F3D(3,3) = lam(2)
        call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                matId,nivGp,finite,isw,err,gp,ELM)
        if (isw.ne.3) return

        if (finite.ge.1) then

c         finite strains           

          K      = (cc3D(2,2)+cc3D(2,3)) / lam(2)
          stre0  = stre3D(2)
          !write(*,'(2(g12.5))') lam(2), stre3D(2)
          iter = 0

          do while (abs(stre3D(2)/stre0).gt.1.d-10
     *                 .and.abs(stre3D(2)).gt.1.d-10.and.iter.le.50)

            iter   = iter + 1
            call pmove(iv1,iv2,nivGp)
            lam(2) = lam(2) - stre3D(2) / K
            F3D(2,2) = lam(2)
            F3D(3,3) = lam(2)
            call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)
            K      = (cc3D(2,2)+cc3D(2,3)) / lam(2)
            !write(*,'(2(g12.5))') lam(2), stre3D(2)
          enddo
          !write(*,*)
          if (iter.gt.50) call prgError(1,'matlib1d','no convergence!')
          stre(1) = stre3D(1)
          lam(3)  = lam(2)
          dlam(2) = - cc3D(2,1)/lam(1) / K 
          dlam(3) = dlam(2)
          cc(1,1) = (1.d0/lam(1) - 2.d0/lam(2)*dlam(2)) * stre(1)
     *         + cc3D(1,1)/lam(1) + (cc3D(1,2)+cc3D(1,3))/lam(2)*dlam(2)

        else

c         small strains

          K      = cc3D(2,2) + cc3D(2,3)
          stre0  = stre3D(2)
          !write(*,'(2(g12.5))') lam(2), stre3D(2)
          iter = 0
          do while (abs(stre3D(2)/stre0).gt.1.d-10
     *                 .and.abs(stre3D(2)).gt.1.d-10.and.iter.le.50)
            iter   = iter + 1
            call pmove(iv1,iv2,nivGp)
            lam(2) = lam(2) - stre3D(2) / K
            F3D(2,2) = lam(2)
            F3D(3,3) = lam(2)
            call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)
            K      = cc3D(2,2) + cc3D(2,3)
            !write(*,'(2(g12.5))') lam(2), stre3D(2)
          enddo
          !write(*,*)
          if (iter.gt.50) call prgError(2,'matlib1d','no convergence!')
          stre(1) = stre3D(1)
          lam(3)  = lam(2)
          dlam(2) = - cc3D(2,1) / K 
          dlam(3) = dlam(2)
          cc(1,1) = cc3D(1,1) + (cc3D(1,2)+cc3D(1,3))*dlam(2)

        endif

      elseif (sss.eq.2) then 

c...    plane strain (membrane) ....................

        call pzero(F3D,9)
        lam(2) = 1.d0 / sqrt(lam(1))
        lam(3) = 1.d0
        F3D(1,1) = lam(1)
        F3D(2,2) = lam(2)
        F3D(3,3) = lam(3)
        call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                matId,nivGp,finite,isw,err,gp,ELM)
        if (isw.ne.3) return

        if (finite.ge.1) then

c         finite strains

          K      = (cc3D(2,2)+stre3D(2)) / lam(2)
          stre0  = stre3D(2)
          !write(*,'(2(g12.5))') lam(2), stre3D(2)
          iter = 0
          do while (abs(stre3D(2)/stre0).gt.1.d-10
     *               .and.abs(stre3D(2)).gt.1.d-10.and.iter.le.50)
            iter   = iter + 1
            call pmove(iv1,iv2,nivGp)
            lam(2)   = lam(2) - stre3D(2) / K
            F3D(2,2) = lam(2)
            call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)
            K      = (cc3D(2,2)+stre3D(2)) / lam(2)
            !write(*,'(2(g12.5))') lam(2), stre3D(2)
          enddo
          !write(*,*)
          if (iter.gt.50) call prgError(3,'matlib1d','no convergence!')
          stre(1) = stre3D(1)
          dlam(2) = - cc3D(2,1)/lam(1) / K 
          dlam(3) = 0.d0
          cc(1,1) = (cc3D(1,1)+stre(1))/lam(1) 
     *             +(cc3D(1,2)-stre(1))/lam(2) * dlam(2)

        else

c         small strains

          K      = cc3D(2,2)
          stre0  = stre3D(2)
          !write(*,'(2(g12.5))') lam(2), stre3D(2)
          iter = 0
          do while (abs(stre3D(2)/stre0).gt.1.d-10
     *               .and.abs(stre3D(2)).gt.1.d-10.and.iter.le.50)
            iter   = iter + 1
            call pmove(iv1,iv2,nivGp)
            lam(2)   = lam(2) - stre3D(2) / K
            F3D(2,2) = lam(2)
            call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)
            K      = cc3D(2,2)
            !write(*,'(2(g12.5))') lam(2), stre3D(2)
          enddo
          !write(*,*)
          if (iter.gt.50) call prgError(4,'matlib1d','no convergence!')
          stre(1) = stre3D(1)
          dlam(2) = - cc3D(2,1) / K 
          dlam(3) = 0.d0
          cc(1,1) = cc3D(1,1) + cc3D(1,2) * dlam(2)
        endif

        stre(2) = stre3D(3)

      elseif (sss.eq.3) then 

c...    axisymmetric (membrane shell) ..............

        call pzero(F3D,9)
        lam(2) = 1.d0 / (lam(1)*lam(3))
        F3D(1,1) = lam(1)
        F3D(2,2) = lam(2)
        F3D(3,3) = lam(3)
        call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                matId,nivGp,finite,isw,err,gp,ELM)
        if (isw.ne.3) return

        if (finite.ge.1) then

c         finite strains

          K      = (cc3D(2,2)+stre3D(2)) / lam(2)
          stre0  = stre3D(2)
          !write(*,'(2(g12.5))') lam(2), stre3D(2)
          iter = 0
          do while (abs(stre3D(2)/stre0).gt.1.d-10
     *               .and.abs(stre3D(2)).gt.1.d-10.and.iter.le.50)
            iter   = iter + 1
            call pmove(iv1,iv2,nivGp)
            lam(2)   = lam(2) - stre3D(2) / K
            F3D(2,2) = lam(2)
            call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)
            K      = (cc3D(2,2)+stre3D(2)) / lam(2)
            !write(*,'(2(g12.5))') lam(2), stre3D(2)
          enddo
          !write(*,*)
          if (iter.gt.50) call prgError(5,'matlib1d','no convergence!')
          stre(1)    = stre3D(1)
          dlam(1) = - cc3D(2,1)/lam(1) / K  ! d lam2 / d lam1
          dlam(3) = - cc3D(2,3)/lam(3) / K  ! d lam2 / d lam3
          cc(1,1) = (cc3D(1,1)+stre3D(1))/lam(1) 
     *             +(cc3D(1,2)-stre3D(1))/lam(2) * dlam(1)
          cc(1,2) = (cc3D(1,3)-stre3D(1))/lam(3) 
     *             +(cc3D(1,2)-stre3D(1))/lam(2) * dlam(3)
          cc(2,1) = (cc3D(3,1)-stre3D(3))/lam(1) 
     *             +(cc3D(3,2)-stre3D(3))/lam(2) * dlam(1)
          cc(2,2) = (cc3D(3,3)+stre3D(3))/lam(3) 
     *             +(cc3D(3,2)-stre3D(3))/lam(2) * dlam(3)

        else

c         small strains

          K      = cc3D(2,2)
          stre0  = stre3D(2)
          !write(*,'(2(g12.5))') lam(2), stre3D(2)
          iter = 0
          do while (abs(stre3D(2)/stre0).gt.1.d-10
     *               .and.abs(stre3D(2)).gt.1.d-10.and.iter.le.50)
            iter   = iter + 1
            call pmove(iv1,iv2,nivGp)
            lam(2)   = lam(2) - stre3D(2) / K
            F3D(2,2) = lam(2)
            call matlib3d(matDat,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)
            K      = cc3D(2,2)
            !write(*,'(2(g12.5))') lam(2), stre3D(2)
          enddo
          !write(*,*)
          if (iter.gt.50) call prgError(6,'matlib1d','no convergence!')
          stre(1) = stre3D(1)
          dlam(1) = - cc3D(2,1) / K ! d lam2 / d lam1
          dlam(3) = - cc3D(2,3) / K ! d lam2 / d lam3
          cc(1,1) = cc3D(1,1) + cc3D(1,2) * dlam(1)
          cc(1,2) = cc3D(1,3) + cc3D(1,2) * dlam(3)
          cc(2,1) = cc3D(3,1) + cc3D(3,2) * dlam(1)
          cc(2,2) = cc3D(3,3) + cc3D(3,2) * dlam(3)
        endif

        stre(2) = stre3D(3)

      else 

        call prgError(3,'matlib1d',
     *                    'for 3D material choose sss = 1,2,3 !')

      endif

      return

      end



