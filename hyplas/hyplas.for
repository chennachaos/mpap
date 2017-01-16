
      subroutine hyplas
     * (matData,F,F33,stre,cc,iv1,iv2, THKGP,
     *  NLARGE,NTYPE,isw,error,MATID,NPROPS,IITER)


      implicit none

      integer          isw, KUNLD, IITER, error,
     *                 i,j, NPROPS, IFLAG , MATID,
     *                 IPROPS(3), NTYPE, NLARGE, ii, index

      double precision matData(*), F(2,*), F33, stre(*), cc(4,*), 
     *                 iv1(*), iv2(*), EINCR(4), FINCR(3,3),
     *                 RPROPS(50), STRES(4), RSTAVA(10), RSTAVA2(10),
     *                 RALGVA(10), DMATX(4,4), AMATX(5,5), DETF, THKGP,
     *                 sigY, mu, K, H, fact,dummy(6),
     *                 r0, r1d2, r1

      logical          LALGVA(10), SUFAIL

      data r0, r1d2, r1 / 0.d0, 0.5d0, 1.d0 /

		!material properties

	error=0


	IF(MATID.GT.19) THEN
	IPROPS(1)=MATID
	ELSE 
	IPROPS(1)=0
	ENDIF

      IPROPS(2)=MATID

      if (isw.eq.2) then
	  call pzero(iv1,17)
        call pzero(iv2,17)

      IF(NLARGE.EQ.1)THEN
	  do i=13,15
        iv1(i) = r1
	  iv2(i) = r1  
	  enddo
	ENDIF


	CALL MATISW
     1(   0       ,NLARGE     ,NTYPE      ,
     2    IPROPS     ,LALGVA     ,LALGVA     ,RALGVA     ,RALGVA     ,
     3    RPROPS     ,RSTAVA     ,RSTAVA     ,dummy     ,dummy     )

	call pmove (RALGVA,iv1(1),2)
	call pmove (RALGVA,iv2(1),2)
	iv1(12)=RALGVA(3)
	iv2(12)=RALGVA(3)

	call pmove (RSTAVA,iv1(3) ,9)
	call pmove (RSTAVA,iv2(3) ,9)

       return
      endif



      IPROPS(3)=matData(1)

	do i=2,NPROPS
         RPROPS(i)=matData(i)
	enddo


	if(MATID.eq.5) then
	IFLAG=RPROPS(6)
	call RDDP (   RPROPS     ,IFLAG      )
	endif


      ! use of internal variable vectors iv1 and iv2
	! 1 - 2    : RALGVA     :2
	! 3 - 11   : RSTAV      :9 RSTAV(5) -> iv2(7)
      ! 12       : RALGVA(3)  :1  
	! 13 - 17  : F^-1(n)    :5 for large strains 
	! 13 - 16  : e(n)       :4 for small strains 
 

		 call pmove (iv1(1), RALGVA,2)
		 call pmove (iv1(3) ,RSTAVA,9)
		 call pmove (iv1(3) ,RSTAVA2,9)
           RALGVA(3)=iv1(12)

     	

	IF(NLARGE.EQ.1)THEN

c	STORE  F^-1(n)
	fact=r1/(F(1,1)*F(2,2)-F(1,2)*F(2,1))
	iv2(13) = F(2,2) * fact
      iv2(14) = F(1,1) * fact    
	iv2(15) = r1/F33
      iv2(16) = -F(1,2)  * fact 
	iv2(17) = -F(2,1)  * fact 

C	CALCULATE F(INCR)=F^-1(n) F(n+1)
	do i=1,2
	   FINCR(i,3)=0.d0
	   FINCR(3,i)=0.d0
	 enddo
	 	 FINCR(1,1)=iv1(13)*F(1,1) + iv1(17)*F(1,2)
		 FINCR(1,2)=iv1(16)*F(1,1) + iv1(14)*F(1,2)
	 	 FINCR(2,1)=iv1(13)*F(2,1) + iv1(17)*F(2,2)
		 FINCR(2,2)=iv1(16)*F(2,1) + iv1(14)*F(2,2)
		 FINCR(3,3)=iv1(15)*F33

	ELSE

c	STORE  E(n)
	fact=r1/(F(1,1)*F(2,2)-F(1,2)*F(2,1))
	iv2(13) = F(1,1)- r1
      iv2(14) = F(2,2) - r1    
	iv2(15) = F(1,2)+F(2,1)
      iv2(16) = F33 - r1 
C	CALCULATE E(INCR)= e(n+1) - e(n) 
           EINCR(1)  = F(1,1)-r1-iv1(13)
           EINCR(2)  = F(2,2)-r1-iv1(14)
           EINCR(3)  = F(1,2)+F(2,1)-iv1(15) !2e12!!
	     EINCR(4)  = F33-r1-iv1(16)

	ENDIF

      DETF=F33*(F(1,1)*F(2,2)-F(1,2)*F(2,1))

C ----------------------------------------------------------------------
	KUNLD=0 !UNLOADING

	call MATISU
     1(   DETF       ,NLARGE     ,NTYPE      ,SUFAIL     ,THKGP      ,
     3    EINCR      ,FINCR      ,IPROPS     ,LALGVA     ,RALGVA     ,
     4    RPROPS     ,RSTAVA     ,STRES      )

	if(SUFAIL) 	then
	error=-1
	return
	endif


      call MATICT
     1(   DETF       ,IITER      ,KUNLD      ,4          ,5      ,
     2    NLARGE     ,NTYPE      ,
     3    AMATX      ,DMATX      ,EINCR      ,FINCR      ,IPROPS  ,
     4    LALGVA     ,RALGVA     ,RPROPS     ,RSTAVA     ,RSTAVA2 ,
     5    STRES      )





	IF(NLARGE.EQ.1)THEN

	cc(1,1)=AMATX(1,1)-STRES(1)
	cc(1,2)=AMATX(1,4)
	cc(1,3)=r1d2*(AMATX(1,2)-STRES(3)+AMATX(1,3))
	cc(2,1)=AMATX(4,1)
	cc(2,2)=AMATX(4,4)-STRES(2)
	cc(2,3)=r1d2*(AMATX(4,2)+AMATX(4,3)-STRES(3))
	cc(3,1)=r1d2*(AMATX(2,1)-STRES(3)+AMATX(3,1))
	cc(3,2)=r1d2*(AMATX(2,4)+AMATX(3,4)-STRES(3))
	cc(3,3)=r1d2*r1d2*(AMATX(2,2)-STRES(2)+AMATX(2,3)
     * +  AMATX(3,2)+AMATX(3,3)-STRES(1) )

	do i=1,4
	cc(i,4)=0.d0
	cc(4,i)=0.d0
	enddo

	ELSE
      call pmove(DMATX(1,1),cc(1,1),16)
	do i=1,4
           cc(i,4)=0.d0
           cc(4,i)=0.d0
	enddo
	ENDIF


C      call pmove(STRES,stre,4)

        stre(1) = STRES(1)
        stre(2) = STRES(2)
        stre(3) = STRES(4)
        stre(4) = STRES(3)



	call pmove (RALGVA,iv2(1),2)	
	call pmove (RSTAVA,iv2(3),9)
      iv2(12)=RALGVA(3)

      return

      end
