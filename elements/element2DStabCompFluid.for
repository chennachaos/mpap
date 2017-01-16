
      integer function element2DStabCompFluid
     *                    (elmDat,timDat,xl,ul,dt,s,p,
     *                     ndf,ndm,nen,nelm,isw)

c  ***************************************************************************
c  *                                                                         *
c  *                                                                         *
c  *                                                                         *
c  *  ---------------------------------------------------------------------  *
c  *    paramter and constants                                               *
c  *                                                                         *
c  *                                                                         * 
c  *  ---------------------------------------------------------------------  * 
c  *                                                                         *
c  *  ---------------------------------------------------------------------  * 
c  *                                                                         *
c  *  ---------------------------------------------------------------------  * 
c  *                                                                         *
c  ***************************************************************************

      implicit none

      integer          ndf, ndm, nen, nelm, isw
      double precision elmDat(*), timDat(*), dt, ul(ndf,*), xl(ndm,*), 
     *                 s(nen*ndf,*), p(*)

      double precision mu, rho, area2, h, nRC(3), z(3), beta(6),
     *                 wa(2,3), ua(3,18), dua(3,18), u(3), uw(2), 
     *                 gu(2,2), f(2), xa(2,3), urdr, urdrr, rad, rad2, 
     *                 gp(2), du(2), trgu, rstab(2), rstabi(3), 
     *                 drdui(2,2), drduj(2,2), tau(3),  
     *                 dtau(2,3), dtaudx(2,3,2), dtaudh(3),  
     *                 sg(4,4), shp(3,3), r(3), k(3,3), 
     *                 dvol, xsj, dgu(2,2,2,3), dtrgu(2,3), dgp(2,2,3), 
     *                 darea2(2,3), darea2r(2,3), dhfact, Dshp(2,3,2,3), 
     *                 drdxj(2,2), drduidw(2), fact, hlp,
     *                 r0, r1, r2, r4, r1d3, r1d2, r1d4, r1d6, r1d216,
     *                 pi, r2pi, p6, p7, p8, b17, b33, bALE,

     *                 dtemp, b, c, st, a(3)  ! Sony

      integer i, j, i0, j0, ii, jj, l, kk, lint, iord, ie, stab

      logical errck, stabm, stabc, axsy

      data r0, r1d4, r1d2, r1, r2, r4 
     *         / 0.d0, 0.25d0, 0.5d0, 1.d0, 2.d0, 4.d0 /

      r1d3  = 1.d0 /  3.d0
      r1d6  = r1d2 * r1d3
      r1d216= r1d6 * r1d6 * r1d6
      pi    = acos(-1.d0)
      r2pi  = r2 * pi
      

      p6    = timDat(6)
      p7    = timDat(7)
      p8    = timDat(8)

c...  set control paramters and element constants







      element2DStabCompFluid = 0


c            1 2 3 4 5  6  7  8  9 10 11 12 13 14 15 16,17,18,19,20
      go to(14,2,2,2,2,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14),isw

      call prgerror(1,'element3DStabIncompHighReFluid','what ?')
      
c=============================================================================
 2    continue ! isw 2 -> compute residual
               !     3 -> compute residual + tangent for  u_{n+1}
               !     4 -> compute residual + tangent for du_{n+1}
               !     5 -> compute residual + full tangent for u_{n+1} 
c-----------------------------------------------------------------------------

c               .     
c...  compute u, u, x, w ...........................

      do i=1, 3
        do ii=1, 2
          ua (ii,i) = p7 * ul(ii,i  ) + (r1-p7) * ul(ii,i+3) 
          dua(ii,i) = p8 * ul(ii,i+6) + (r1-p8) * ul(ii,i+9) 
          xa (ii,i) = p6 * xl(ii,i  ) + (r1-p6) * xl(ii,i+3) 
          wa (ii,i) = p7 * xl(ii,i+6) + (r1-p7) * xl(ii,i+9) 

          !write(*,*) wa(ii,i)

        enddo
        ua(3,i) = ul(3,i) ! pressure
      enddo

c...  characteristic element size
   
      area2 =  (xa(1,1)-xa(1,3))*(xa(2,2)-xa(2,3)) 
     *       - (xa(1,2)-xa(1,3))*(xa(2,1)-xa(2,3))

      h     = sqrt(r2/pi * area2)  ! diameter of circle

c...




















      return
c==========================================================================
 6    continue ! isw 6 -> compute residual + tangent for  x_{n+1}
               !     8 -> compute residual + full tangent for x_{n+1} Lagr
               !     9 -> compute residual + full tangent for x_{n+1} ALE
c-----------------------------------------------------------------------------


      return

c==========================================================================
 7    continue ! isw 7 -> compute tangent for dx_{n+1}
c-----------------------------------------------------------------------------


      return

c==============================================================================
 14   continue ! isw 14-15 -> not used
c-----------------------------------------------------------------------------

      write(*,'(/9x,a/)') 
     *   'element3DStabIncompHighReFluid: what do you want ?'

      return

c==============================================================================
      end
















