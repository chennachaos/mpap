
      subroutine matlib2d(matData,F,F33,stre,cc,iv1,iv2,dt,
     *                    matId,nivGp,finite,sss,isw,err,gp,ELM)

      implicit none

      integer          ELM ! do not touch

      integer          matId, nivGp, finite, sss, isw, gp, err,
     *                 iter, ndm, matDim, i, j

      double precision matData, F(2,*), F33, stre(*), cc(4,*), iv1, 
     *                 iv2, dt,
     *                 F3D(3,3), stre3D(6), cc3D(6,6), stre0, K, C(4,4),
     *                 invF(2,2), invF33, fact, dmy(5), cc1(4,4), r0

c
c   isw
c    1   set nivGp and exit
c    2   initialise internal variables
c    3   full stress update (compute stress, internal variables and tangent tensor)
c

      data r0 / 0.d0 /

      if (isw.lt.1 .or. isw.gt.3) 
     *  call prgError(1,'matlib2d','invalid isw!')


      ndm = matDim(matId)

      if (ndm.lt.2) call prgError(3,'matlib2d','invalid matId!')

      !write(*,*) matId


      if (ndm.eq.3 .and. isw.eq.3) goto 500

c.... EVERYTHING FOR: (isw=1,2,3 with ndm=2)  and  (isw=1,2 with ndm=3) ........................

          ! 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
      goto (1,1,1,1,1,1,7,1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1,19,20,20,20,20,

          !24 25 26 27 28 29 30
     *     20,20,20,20,20,20,20), matId

      call prgError(4,'matlib2d','matId / ndm mismatch ?!')

c.........................................................................................

 1    continue ! 3D material subroutines

      call matlib3d(dmy,dmy,dmy,dmy,iv1,iv2,dt,
     *              matId,nivGp,finite,isw,err,gp,ELM)

      return

c.........................................................................................

 7    continue  ! multiscale material 2D: call only for isw = 3

      nivGp = 0 

      if (isw.ne.3) return
      
c      call multiscaleMaterial2D(matData,F,stre,cc,finite,gp,err,ELM)

      return

c.........................................................................................

 19   continue ! small strain elastic lamina in 2D

      nivGp = 0

      if (isw.eq.1) return

      if (sss   .gt.1) 
     *  call prgError(1,"matlib2d","this works only for plane stress!")
      if (finite.ne.0) 
     *  call prgError(2,"matlib2d","this works only for small strains!")

      call small_strain_elastic_lamina2D(matData,F,stre,cc)

      return

c.........................................................................................

 20   continue ! Hyplas material

      nivGp = 17

      if (isw.eq.1) return

c      call hyplasMaterial(matData,F,F33,stre,cc,iv1,iv2,
c     *                    matId,finite,sss,isw,err,ELM)
     
    

      return

c.........................................................................................

 500  continue

c.... FULL STRESS UPDATE BASED ON 3D MATERIAL SUBROUTINE: (isw=3 with ndm=3) .....................

      if (sss.le.1) then

c...    plane stress ..............

        call pzero(F3D,9)
        F3D(1,1) = F(1,1)
        F3D(2,1) = F(2,1)
        F3D(1,2) = F(1,2)
        F3D(2,2) = F(2,2)
        F3D(3,3) = F33

        call matlib3d(matData,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                matId,nivGp,finite,isw,err,gp,ELM)

        if (finite.ge.1) then

c         finite strains

          K     = cc3D(3,3) / F3D(3,3)
          stre0 = stre3D(3)
          !write(*,'(2(g12.5))') F3D(3,3), stre3D(3)
          iter = 0
          do while (abs(stre3D(3)/stre0).gt.1.d-14
     *             .and.abs(stre3D(3)).gt.1.d-10.and.iter.le.50)
            iter = iter + 1
            call pmove(iv1,iv2,nivGp)
            F3D(3,3) = F3D(3,3) - stre3D(3) / K
            call matlib3d(matData,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)
            K = cc3D(3,3) / F3D(3,3)
            !write(*,'(2(g12.5))') F3D(3,3), stre3D(3)
          enddo
          !write(*,*)
          F33 = F3D(3,3)
          fact = 1.d0 / (F(1,1)*F(2,2) - F(1,2)*F(2,1))
          invF(1,1) =   F(2,2) * fact
          invF(1,2) = - F(1,2) * fact
          invF(2,1) = - F(2,1) * fact
          invF(2,2) =   F(1,1) * fact
          invF33    =  1.d0 / F33
          cc(1,1)   = cc3D(1,1)
          cc(2,1)   = cc3D(2,1)
          cc(3,1)   = cc3D(4,1)
          cc(4,1)   = cc3D(3,1)
          cc(1,2)   = cc3D(1,2)
          cc(2,2)   = cc3D(2,2)
          cc(3,2)   = cc3D(4,2)
          cc(4,2)   = cc3D(3,2)
          cc(1,3)   = cc3D(1,4)
          cc(2,3)   = cc3D(2,4)
          cc(3,3)   = cc3D(4,4)
          cc(4,3)   = cc3D(3,4)
          cc(1,4)   = cc3D(1,3)
          cc(2,4)   = cc3D(2,3)
          cc(3,4)   = cc3D(4,3)
          cc(4,4)   = cc3D(3,3)
          call pushfwrd2D_4(C,invF,invF33,cc) ! in fact it's "pull-back"
          K = 1.d0 / C(4,4)
          C(1,1) = C(1,1) - K * C(1,4) * C(4,1)
          C(1,2) = C(1,2) - K * C(1,4) * C(4,2)
          C(1,3) = C(1,3) - K * C(1,4) * C(4,3)
          C(2,1) = C(2,1) - K * C(2,4) * C(4,1)
          C(2,2) = C(2,2) - K * C(2,4) * C(4,2)
          C(2,3) = C(2,3) - K * C(2,4) * C(4,3)
          C(3,1) = C(3,1) - K * C(3,4) * C(4,1)
          C(3,2) = C(3,2) - K * C(3,4) * C(4,2)
          C(3,3) = C(3,3) - K * C(3,4) * C(4,3)
          call pushfwrd2D_4(cc,F,F33,C)

        else

c         small strains

          K     = cc3D(3,3)
          stre0 = stre3D(3)
          !write(*,'(2(g12.5))') F3D(3,3), stre3D(3)
          iter = 0
          do while (abs(stre3D(3)/stre0).gt.1.d-14
     *             .and.abs(stre3D(3)).gt.1.d-10.and.iter.le.50)
            iter = iter + 1
            call pmove(iv1,iv2,nivGp)
            F3D(3,3) = F3D(3,3) - stre3D(3) / K
            call matlib3d(matData,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                    matId,nivGp,finite,isw,err,gp,ELM)
            K = cc3D(3,3)
            !write(*,'(5(g12.5))') F3D(3,3), stre3D(3), 
     *      !                      iv2(3), iv2(5), iv2(6)
          enddo
          !write(*,*)
          F33 = F3D(3,3)
          K   = 1.d0 / K
          cc(1,1) = cc3D(1,1) - K * cc3D(1,3) * cc3D(3,1)
          cc(1,2) = cc3D(1,2) - K * cc3D(1,3) * cc3D(3,2)
          cc(1,3) = cc3D(1,4) - K * cc3D(1,3) * cc3D(3,4)
          cc(1,4) = r0
          cc(2,1) = cc3D(2,1) - K * cc3D(2,3) * cc3D(3,1)
          cc(2,2) = cc3D(2,2) - K * cc3D(2,3) * cc3D(3,2)
          cc(2,3) = cc3D(2,4) - K * cc3D(2,3) * cc3D(3,4)
          cc(2,4) = r0
          cc(3,1) = cc3D(4,1) - K * cc3D(4,3) * cc3D(3,1)
          cc(3,2) = cc3D(4,2) - K * cc3D(4,3) * cc3D(3,2)
          cc(3,3) = cc3D(4,4) - K * cc3D(4,3) * cc3D(3,4)
          cc(3,4) = r0
          cc(4,1) = r0
          cc(4,2) = r0
          cc(4,3) = r0
          cc(4,4) = r0

        endif


        do i=1,4
           do j=1,4
            cc1(i,j) = cc(i,j)
           enddo
        enddo

        cc(1,1) = cc1(1,1)
        cc(2,1) = cc1(2,1)
        cc(3,1) = cc1(4,1)
        cc(4,1) = cc1(3,1)
        
        cc(1,2) = cc1(1,2)
        cc(2,2) = cc1(2,2)
        cc(3,2) = cc1(4,2)
        cc(4,2) = cc1(3,2)
        
        cc(1,3) = cc1(1,4)
        cc(2,3) = cc1(2,4)
        cc(3,3) = cc1(4,4)
        cc(4,3) = cc1(3,4)
        
        cc(1,4) = cc1(1,3)
        cc(2,4) = cc1(2,3)
        cc(3,4) = cc1(4,3)
        cc(4,4) = cc1(3,3)
        
        stre(1) = stre3D(1)
        stre(2) = stre3D(2)
        stre(3) = 0.d0
        stre(4) = stre3D(4)

      elseif (sss.eq.2) then 

c...    plane strain ..............

        call pzero(F3D,9)
        F3D(1,1) = F(1,1)
        F3D(2,1) = F(2,1)
        F3D(1,2) = F(1,2)
        F3D(2,2) = F(2,2)
        F3D(3,3) = F33
        call matlib3d(matData,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                matId,nivGp,finite,isw,err,gp,ELM)

        do i=1,4
           stre(i) = stre3D(i)
           do j=1,4
            cc(i,j) = cc3D(i,j)
           enddo
        enddo

c        stre(1) = stre3D(1)
c        stre(2) = stre3D(2)
c        stre(3) = stre3D(4)
c        stre(4) = stre3D(3)


      elseif (sss.eq.3) then 

c...    axisymmetric ..............

        call pzero(F3D,9)
        F3D(1,1) = F(1,1)
        F3D(2,1) = F(2,1)
        F3D(1,2) = F(1,2)
        F3D(2,2) = F(2,2)
        F3D(3,3) = F33

        call matlib3d(matData,F3D,stre3D,cc3D,iv1,iv2,dt,
     *                matId,nivGp,finite,isw,err,gp,ELM)


        do i=1,4
           stre(i) = stre3D(i)
           do j=1,4
            cc(i,j) = cc3D(i,j)
           enddo
        enddo

      else

        call prgError(2,'matlib2d','choose pssPsnAsy = 1,2,or 3 !')

      endif

      return

      end





