	subroutine sf_vonmises_isotropic_elastoplasticity
     * (matData,F,stre,cc,iv1,iv2,finiteInt,isw,error)
c
c	Infinitesimal/Finite Strain Isotropic Elastoplasticity (Von Mises) 
c
c            matData(1)   - mu
c
c            matData(2)  - K
c
c            matData(3)   - Yield Stress
c
c            matData(4)   - H

c            matData(5)   - tolerance for local Newton iterations
c
c            matData(6)  - maximum number of iterations for Newton procedures
c
c            matData(7)  - hardening type, htype
c                            |  htype < = 1 -> linear hardening (use with H=0 for perfect plasticity)
c                            |      > 1 -> exponential hardening
c            matData(8) - integration parameter parameter nip
c                            |  nip < = 2 -> EM with 2 terms
c                            |      > 2 -> EM with nip
c                            |             terms
c
c     isw = 2 -> initialize internal variables 
c
c     isw = 3 -> compute stress, update internal variables and
c                compute consistent spatial tangent tensor
c
c
c=============================================================================

      implicit none
      double precision matData(*), F(3,*), stre(*), cc(6,*), iv1(*), 
     * iv2(*), S(6), dSdC(6,6), dSdCp(6,6),
     * fun(6), dfdCp(6,6), dfdC(6,6), 
     * Rg, dRgdCp(6),  dRgdC(6),  dRgdg, 
     * Rp(6), dRpdCp(6,6),  dRpdC(6,6),  dRpdg(6), 
     * C(6), invC(6),diCdC(6,6), r1d3,r1,
     * Cp(6), invCp(6), diCpdCp(6,6), Cpn(6), invCpn(6), 		
     * R(7), dR(7,7), dvar(7), det_u, Rnorm, 
     * tmp41(6,6), imtx77(7,7), imtx67(6,7), ccref(6,6), mtx76(7,6),
     * tol, fact,detF, alpha, alphan, g, level,
     * r2, r3, r1d2, r2d3, sqr2d3, det_s

	integer    i, j, iter, maxiter, isw, P(7), 
     * htype, finiteInt, error
	logical    finite

      data r1, r3, r1d2, r2 / 1.d0, 3.d0, 0.5d0, 2.d0 /
      
              finite = (finiteInt.ge.1)
              	error=0
              	

        if (isw.ne.3) then
c initialize internal variables

       	if (finite)then
        call unit_s_p(iv1)
        call unit_s_p(iv2)
        else
        call pzero(iv1(1),6)
        call pzero(iv2(1),6)
        endif
        
        iv1(7) = 0.d0
        iv2(7) = 0.d0

		return
        endif
        
c	initialize Cp_n, alpha_n
	do i=1,6
	Cpn(i)=iv2(i)
	enddo
	alphan = iv2(7)
c	Cpn=Cp_n, alpha=alpha_n
	do i=1,6
	Cp(i)=Cpn(i)
	enddo
	alpha=alphan
       
        

       	if (.not.finite)then
       	      call sstrain(matData,F,stre,cc, alphan, Cpn, alpha, Cp)
       	      goto 100
       	endif
   
c	set parameters
	r2d3=r2/r3
	sqr2d3=r2d3**r1d2

      tol     = matData(5)
      maxiter = nint(matData(6))
	htype= nint(matData(7))
	level=r2
	if(matData(8).gt.r2) level=float(nint(matData(8)))
	
c	maxiter=100
c	tol=1.d-6

c	compute C, invC, J from F
	call mult_s_uTu_p (C,F,F)
	call inverse_s(invC, C)
	detF=det_u(F)
	
c	d(C^-1)/dC=-(C^-1)(C)(C^-1)
	call diff_ATA_m(diCdC,invC)

c	compute Cpn^(-1), Cp^(-1)
	call inverse_s(invCpn, Cpn)
	call set_s_p(invCp,invCpn)

c	d(Cp^-1)/dCp=-(Cp^-1)(Cp)(Cp^-1)
	call diff_ATA_m(diCpdCp,invCp)
	g=0.d0

c	update stress
      call piola2(S,dSdC,dSdCp,
     * C,invC,diCdC,Cp,invCp,diCpdCp,matData)

	call set_4_f (ccref, dSdC,r2)


c	compute Rg, f and their derivatives wrt C and Cp and g
	call residg( Rg, dRgdCp, dRgdC, dRgdg, 
     * fun, dfdCp,  dfdC,
     * S, dSdC, dSdCp, C, Cp, alpha, matData, htype)
      
	if (Rg.gt.0.d0) then   !PLASTIC
	
c	start iterations!!!!!!!!!!!!!!!!!!!!!!!!!!!!ITERATION START
	do iter=1, maxiter
	
c
c	compute Rp and its derivatives wrt C, Cp and g
	call residp(Rp, dRpdCp, dRpdC, dRpdg, 
     * fun, dfdCp,  dfdC, diCpdCp, 
     * Cp, invCp, invCpn, g, iter, level)
c
c	Construct and solve the matrix system for Cp, g:
c	dR.dvar=R
	do i=1,6
		do j=1,6
		dR(i,j)=dRpdCp(i,j)
		enddo
	enddo
	do i=1,6
		dR(i,7)=dRpdg(i)
	enddo
	do i=1,6
		dR(7,i)=dRgdCp(i)
	enddo
		dR(7,7)=dRgdg
	do i=1,6
		R(i)=-Rp(i)
	enddo	
		R(7)=-Rg

c	check convergence on normR
	Rnorm=0.d0
	do i=1,7
	Rnorm=Rnorm+R(i)**r2
	enddo
	Rnorm=(Rnorm**r1d2)
c	if (iter.gt.15) then
c          call prgerror(1,'sf_vonmises_isotropic_elastoplasticity',
c     *                  ' local iterations>15!')
c       stop
c	endif
c
      call decompLR_matrix(dR,P,7,0)
	if (Rnorm.le.tol) goto 200	!>>>>>>>>>>>>converged

c	solve the matrix system for Cp, g:

	call solve_matrix   (dR,P,R,dvar,7)
c	update Cp and g
          g = g+dvar(7)
c	if(g.lt.0.d0) then
c          call prgerror(2,'sf_vonmises_isotropic_elastoplasticity',
c     *                  ' negative plastic multiplier!')
c	 stop
c	endif
	   alpha=alphan+sqr2d3*g
          do i=4, 6
            dvar(i) = r1d2 * dvar(i)
          enddo
          call set_s_ap(Cp,dvar)

c	update stress and Rg
	call inverse_s(invCp, Cp)
c	d(Cp^-1)/dCp=-(Cp^-1)(Cp)(Cp^-1)
	call diff_ATA_m(diCpdCp,invCp)


      call piola2(S,dSdC,dSdCp,
     * C,invC,diCdC,Cp,invCp,diCpdCp,matData)


c	compute Rg, f and their derivatives wrt C, Cp and g
	call residg( Rg, dRgdCp, dRgdC, dRgdg, 
     * fun, dfdCp,  dfdC,
     * S, dSdC, dSdCp, C, Cp, alpha, matData, htype)

	enddo	!NEWTON!!!!!!!!!!!!!!!!!!!!!!!!!!!!	ITERATION end
	
	if (iter.eq.maxiter) then
	error=-1
	return
	endif

c	construct spatial ep tangent modulus

200	continue

	call inverse_matrix(dR,imtx77,P,7)      ! 

		do i=1,6
		do j=1,7
		imtx67(i,j)=imtx77(i,j)
		enddo
		enddo

	do i=1,6
		do j=1,6
		mtx76(i,j)=dRpdC(i,j)
		enddo
	enddo
	do i=1,6
		mtx76(7,i)=dRgdC(i)
	enddo

	call mult_matrix_p (tmp41,imtx67,mtx76,6,7,6)

	call mult_matrix_m (ccref,dSdCp,tmp41,6,6,6)
	call set_4_ap (ccref, dSdC)          ! 


	call set_4_f(ccref,ccref,r2)        ! ccref=2(dSdC-dSdCp:tmp43)OK

	endif              ! PLASTIC


      call pushfwrd_s(stre,F,S)
      call pushfwrd_4(cc,F,ccref)
      fact = 1.d0 / det_u(F)

      do i=1, 6
        stre(i) = stre(i) * fact
        do j=1, 6
          cc(i,j) = cc(i,j) * fact
        enddo
      enddo
      
c	store converged internal variables

100	continue
	do i=1,6
	iv2(i)=Cp(i)
	enddo
	iv2(7)=alpha
	

	return
	end
	
	
		  subroutine residp( Rp, dRpdCp, dRpdC, dRpdg, 
     * fun, dfdCp,  dfdC, diCpdCp, 
     * Cp, invCp, invCpn, g, iter, level )
      implicit none
      double precision 	Rp(6), dRpdCp(6,6), dRpdC(6,6),  
     * dRpdg(6), fun(6), dfdCp(6,6),  dfdC(6,6), 
     * invCp(6), Cp(6), invCpn(6), rct,
     * invCpf(3,3), diCpdCp(6,6), temp43(6,6),
     * stemp(6), utemp(3,3), temp41(6,6), temp42(6,6),
     * term(6), dtermdC(6,6), dtermdCp(6,6),dtermdg(6), 
     * nterm(6), dntermdC(6,6), dntermdCp(6,6),dntermdg(6), 
     * fact1, g, r1d2, fact2, r3, r2, level, r1

	  integer iter,i, jct

      data r3, r1d2, r1, r2 / 3.d0, 0.5d0, 1.d0, 2.d0 /

  	  call mult_u_ss_p ( invCpf, invCp, fun)

c	Rp(6)
	  fact1=g*g/r2
	  call set_s_p( Rp, Cp) !	Rp=	Cp
	  call mult_u_ss_p ( utemp, invCpn, Cp)	!utemp=Cpn^(-1) Cp
	  call mult_s_su_p ( stemp, Cp, utemp)	!stemp=Cp Cpn^(-1) Cp
	  call set_s_am( Rp, stemp) !	Rp=	Cp - CpCpn^(-1)Cp
	  call set_s_af( Rp, fun, g) !	Rp=	Cp - CpCpn^(-1)Cp+(g)f
	  call mult_s_su_f (term,fun,invCpf,fact1) !	term= (g^2/2)(finvCpf) 
	  call set_s_ap(Rp,term) !	Rp=	Cp - CpCpn^(-1)Cp+(g)f + (g^2/2)(finvCpf) 

c	dRp/dC(6,6)
	  call set_4_f( dRpdC, dfdC, g) 	!	dRpdC= g.dfdC
	  call diff_TAT_p( temp41, invCpf) !	temp41=	d(fCp^(-1)f)/df
	  call mult_4_44_f( dtermdC, temp41, dfdC, fact1) !	temp42=	d(fCp^(-1)f)/df.dfdC
	  call set_4_ap( dRpdC, dtermdC) !	dRpdC= g.dfdC + (g^2/2)d(fCp^(-1)f)/df.dfdC OKOK

c	dRp/dCp(6,6)
	  call unit_4_p( dRpdCp) !	dRpdCp	= I
   	  call diff_TAT_p( temp42, utemp) ! temp42=	d(Cp Cpn^(-1) Cp)/dCp
	  call set_4_am( dRpdCp, temp42) !	dRpdCp	= I-d(Cp Cpn^(-1) Cp)/dCp
	  call set_4_af( dRpdCp, dfdCp, g) !	dRpdCp	=           +(g)dfdCp
	  call mult_4_44_f(  dtermdCp,temp41,dfdCp, fact1) ! dtermdCp=	(g^2/2)[d(fCp^(-1)f)/df:dfdCp
	  call diff_ATA_p( temp42, fun) ! temp42=d(fCp^(-1)f)/dCp^(-1)
	  call mult_4_44_af( dtermdCp,temp42,diCpdCp, fact1)		
	  call set_4_ap( dRpdCp, dtermdCp) ! + d(fCp^(-1)f)/dCp^(-1):dCp^(-1)/dCp] OKOK

c	dRp/dg(6)
	  call set_s_p( dRpdg, fun)				! dRpdg=f
	  call mult_s_su_p( dtermdg, fun, invCpf)		! dtermdg=fCp^(-1)f
	  call set_s_af( dRpdg, dtermdg, g)			! dRpdg=f+(g)fCp^(-1)f OKOK

	  if (level.gt.r2) then
	  jct=nint(level)
	  do i=3,jct
c	fact2=fact1*g/real(i)
	  fact2=g/real(i)
	  call mult_s_su_f ( nterm, term, invCpf, fact2)	!nterm=fact2(term.invCpf) 

	  call mult_u_ss_p ( utemp, term, invCp)
	  call diff_UT_p( temp41, utemp) !dntermdf(/fact2)
	  call diff_TU_p( temp42, invCpf) !dntermdterm(/fact2)
	  call diff_ATB_p( temp43,term, fun) !dntermdCp^(-1)(/fact2)

	  call mult_4_44_f( dntermdC, temp41, dfdC, fact2) !dntermdC= fact2*[dntermdf.dfdC	
	  call mult_4_44_af( dntermdC, temp42, dtermdC, fact2) !+dntermdterm.dtermdC]

 	  call mult_4_44_f( dntermdCp,temp43,diCpdCp, fact2)	!dntermdCp= fact2*[dntermdCp^(-1).dCp^(-1)dCp
	  call mult_4_44_af( dntermdCp,temp41,dfdCp, fact2)! +dntermdf.dfdCp
	  call mult_4_44_af( dntermdCp,temp42,dtermdCp, fact2)! +dntermdterm.dtermdCp]	

	  call mult_s_su_p ( dntermdg, dtermdg, invCpf) !dntermdg= dtermdg.invCpf

	  call set_s_p( term, nterm) 
	  call set_4_p( dtermdCp, dntermdCp)
	  call set_4_p( dtermdC, dntermdC)	
	  call set_s_p( dtermdg, dntermdg)

	  call set_s_ap( Rp, term)
	  call set_4_ap( dRpdCp, dtermdCp)
	  call set_4_ap( dRpdC, dtermdC)	
	  call set_s_af( dRpdg, dtermdg, fact1)

	  fact1=fact1*fact2

	  enddo
	  endif

	  return
	  end
	  
	  
	  	  subroutine residg(	Rg, dRgdCp, dRgdC, dRgdg, 
     *				fun, dfdCp,  dfdC,
     *				S, dSdC, dSdCp, C, Cp, alpha, d, htype)
      implicit none
      double precision 	Rg, dRgdCp(6),  dRgdC(6),  dRgdg, norm,
     * fun(6), dfdCp(6,6),  dfdC(6,6), 
     * S(6), dSdC(6,6), dSdCp(6,6), C(6), Cp(6), d(*),
     * stemp(6),stemp2(6), CS(3,3), dSCSdS(6,6), SdCSdC(6,6), 
     * temp4(6,6), CdSCpdS(6,6), SCp(3,3), 
     * trcs, trcscs, fact1, fact2, fact3, fact4, fact6,
     * traceAB_ss, traceAB_uu, sigma0, H, alpha,
     * r1, r2, r3, r1d2, r2d3, r1d3, sqr2d3, utemp(3,3),
     * sigmay, dsigmayda

	  integer i,j, htype

      data r3, r1d2, r1, r2 / 3.d0, 0.5d0, 1.d0, 2.d0 /

c	set parameters
c	r3=r1+r2
	  r1d3=r1/r3
   	  r2d3=r2/r3
	  sqr2d3=r2d3**r1d2
	  sigma0=d(3)
	  H=d(4)

	  do i=1,3
	  do j=1,3
	  utemp(i,j)=0.d0
	  enddo
	  enddo
	  call mult_u_ss_p (CS,C,S)			! CS
	  trcs=traceAB_ss(C,S)				! tr[CS]
	  trcscs=traceAB_uu(CS,CS)			! tr[CSCS]

	  call hard(sigmay, dsigmayda, sigma0, H, alpha, htype)

c	Rg
	  norm=(trcscs-r1d3*trcs*trcs)**r1d2		! norm=(tr[CSCS]-(1/3)tr[CS]^2)^(1/2)
	  Rg=norm-sqr2d3*sigmay
	
c	dRg/dC (=dnorm/dC)
	  fact1=-r2d3*trcs
  	  fact2=r1d2/norm
	  call mult_s_s4_p( stemp, C, dSdC)		! stemp=C:dSdC
	  call set_s_ap( stemp, S)			! stemp=S+C:dSdC
	  call set_s_f( dRgdC, stemp, fact1)		! dRgdC=-(2/3)*tr[CS](S+C:dSdC)
	  call mult_s_su_ap( dRgdC, S, CS)		! dRgdC=SCS-(2/3)*tr[CS](S+C:dSdC)
	  call diff_TAT_p( dSCSdS, CS)			! dSCSdS
	  call mult_s_s4_p( stemp, C, dSCSdS)		! stemp=C:dSCSdS
	  call mult_s_s4_p( stemp2, stemp, dSdC)	! stemp2=(C:dSCSdS)dSdC
	  call set_s_ap( dRgdC, stemp2)			! dRgdC=SCS-(2/3)*tr[CS](S+C:dSdC)+(C:dSCSdS)dSdC
	  call diff_ATA_p( SdCSdC, S)			! SdCSdC
	  call mult_s_s4_ap( dRgdC, C, SdCSdC)	! dRgdC=SCS-(2/3)*tr[CS](S+C:dSdC)+(C:dSCSdS)dSdC+C:SdCSdC
	  call set_s_f( dRgdC,dRgdC, fact2)	! dRgdC=(1/2norm){SCS-(2/3)*tr[CS](S+C:dSdC)+C:(dSCSdS:dSdC+SdCSdC} OKOK*

c	dRg/dCp (=dnorm/dCp)
	  call mult_s_s4_p( stemp, C, dSdCp)		! stemp=C:dSdCp
	  call set_s_f( dRgdCp, stemp, fact1)		! dRgdCp=-(2/3)*tr[CS](C:dSdCp)
	  call mult_s_s4_p( stemp, C, dSCSdS)		! stemp=C:dSCSdS
	  call mult_s_s4_p( stemp2, stemp, dSdCp)	! temp4=C:dSCSdS:dSdCp
	  call set_s_ap( dRgdCp, stemp2)		! dRgdCp=C:dSCSdS:dSdCp-(2/3)*tr[CS](C:dSdCp)
	  call set_s_f( dRgdCp, dRgdCp, fact2)	! dRgdCp=(1/2norm)[C:dSCSdS:dSdCp-(2/3)*tr[CS](C:dSdCp)] OKOK*
	
c	dRg/dg
	 dRgdg=-r2d3*dsigmayda					! dRgdg=-(2/3)H 

c	f 
	  fact3=-r1d3*trcs
	  fact4=r2/norm
	  fact6=-r1/norm
	  call mult_s_us_p (fun,CS,Cp)				!f= CS Cp
	  call set_s_af( fun, Cp, fact3)			!f= CS Cp-(1/3)tr[CS]Cp 
	  call set_s_f( fun,fun, fact4)				!f= (2/norm)[CS Cp-(1/3)tr[CS]Cp] 

c	df/dCp
	  call unit_4_f (dfdCp,fact3)				! dfdCp=-(1/3)tr[CS]I
	  call diff_ATB_p( CdSCpdS, C, Cp)			! CdSCpdS
	  call mult_4_44_ap( dfdCp, CdSCpdS, dSdCp)		! dfdCp=CdSCpdS.dSdCp-(1/3)tr[CS]I
	  call diff_UT_p( temp4, CS)				! temp4=dCSCpdCp
	  call set_4_ap (dfdCp,temp4)				! dfdCp=CdSCpdS.dSdCp-(1/3)tr[CS]I+dCSCpdCp
	  call mult_s_s4_p( stemp, C, dSdCp)			! stemp=C:dSdCp
	  call mult_4_ss_p( temp4, Cp, stemp)			! temp4=Cp:C:dSdCp
	  call set_4_af (dfdCp,temp4, -r1d3)			! dfdCp=(CdSCpdS-(1/3)Cp:C:dSdCp+dCSCpdCp-(1/3)tr[CS]I4
	  call set_4_f (dfdCp,dfdCp, fact4)			! dfdCp=(2/norm){CdSCpdS.dSdCp+dCSCpdCp-(1/3)tr[CS]{I4+Cp:C:dSdCp)}}
	  call mult_4_ss_p( temp4, fun, dRgdCp)		! temp4=f dRgdCp
	  call set_4_af (dfdCp,temp4, fact6)			! dfdCp=(2/norm){CdSCpdS.dSdCp+dCSCpdCp-(1/3)[tr[CS]I4+Cp:C:dSdCp]}
									!		-(1/norm)(f dRgdCp) 
c	df/dC
	  call mult_s_s4_p( stemp, C, dSdC)			! stemp=C:dSdC
	  call set_s_ap( stemp, S)				! stemp=C:dSdC+S
	  call mult_4_ss_p( dfdC, Cp, stemp )			! dfdC=Cp(C:dSdC+S)
	  call set_4_f (dfdC,dfdC, -r1d3)			! dfdC=(-1/3)Cp(C:dSdC+S)
	  call mult_4_44_ap( dfdC, CdSCpdS, dSdC)		! dfdC=(-1/3)Cp(C:dSdC+S)+CdSCpdS.dSdC
	  call mult_u_ss_p (SCp,S,Cp)				! SCp
	  call diff_TU_p( temp4, SCp)				! temp4=d(CSCp/dC)
	  call set_4_ap (dfdC,temp4)				! dfdC=d(CSCp/dC)-(1/3)Cp(C:dSdC+S)+CdSCpdS.dSdC
	  call mult_4_ss_p( temp4, fun, dRgdC)		! temp4=f.dRgdC
	  call set_4_af (dfdC,temp4, -r1d2)			! dfdC=d(CSCp/dC)-(1/3)Cp(C:dSdC+S)+CdSCpdS.dSdC-(1/2)f.dRgdC
	  call set_4_f (dfdC,dfdC, fact4)			! dfdC=(2/norm)dfdC 

	  return
	  end
	  
	  
	        subroutine sstrain(d,F,stre,cc, alphan, epn, alpha, ep)

      implicit none

      double precision d(*), F(3,*), stre(*), cc(6,*), 
     * mu, K , lamb, eps(6), ep(6), epn(6), alphan,
     * tdevstr(6),edevstr(6),s(6),sigmay,H,
     * sdir(6),trtot,norms,alpha,
     * lamda,beta,gamma,ftrial,dot_ss,
     * r1d2, r0, r1, r2, r3, r1d3, r2d3, 
     * sq2d3, r2mu, r3mu, r2d3mu, r4d3mu, sdirP(6,6)


      integer          i,j,ii

c
      data r0, r1d2, r1, r2/ 0.d0, 0.5d0, 1.d0, 2.d0 /
      r3 = 3.d0 
      r2d3 = 2.d0 / 3.d0
      r1d3 = 1.d0 / 3.d0
	  sq2d3=r2d3**r1d2

		do i=1,6
			do j=1,6
				cc(i,j)=0.d0
			enddo
		enddo
		

	
	
	
	

      mu     = d(1)
      K      = d(2)
	 sigmay=d(3)
	 H=d(4)
      lamb   = K - r2d3 * mu
	r2mu=r2*mu
	r3mu=r3*mu
	r2d3mu=r2d3*mu
	r4d3mu=r2d3mu+r2d3mu
	
c     strain tensor
      eps(1) = F(1,1) - r1
      eps(2) = F(2,2) - r1
      eps(3) = F(3,3) - r1
      eps(4) = r1d2 * (F(1,2) + F(2,1))
      eps(5) = r1d2 * (F(2,3) + F(3,2))
      eps(6) = r1d2 * (F(3,1) + F(1,3))
      
c.....calculate deviatoric elastic trial strain
	trtot=eps(1)+eps(2)+eps(3)
		do i=1,3
			tdevstr(i)=eps(i)-r1d3*trtot
			tdevstr(i+3)=eps(i+3)
		enddo	
c	
	do j=1,6
	edevstr(j)=tdevstr(j)-epn(j)
	enddo

c.....calculate deviatoric trial stress 
 
	do j=1,6
	s(j)=r2mu*edevstr(j)
	enddo

c.....norm of deviatoric trial stress
	norms=dot_ss(s,s)**r1d2

c.....Von Mises yield criterion
	ftrial=norms-sq2d3*(sigmay+H*alpha)
c
	if (ftrial.le.r0) then
c..........stress update
		do i=1,3
			stre(i)=s(i)+K*trtot
			stre(i+3)=s(i+3)
		enddo	
c..........tangent modulus update
	do i=1,3
		do j=1,3
			cc(i,j)=K-r2d3mu
		enddo
		cc(i,i)=cc(i,i)+r2mu
		cc(i+3,i+3)=cc(i+3,i+3)+mu
	enddo
	else
c...........elastoplastic stress update
		do j=1,6
		sdir(j)=s(j)/norms
		enddo
c
		lamda=ftrial/(r2mu*(r1+H/r3mu))

		do i=1,6
			s(i)=s(i)-r2mu*lamda*sdir(i)
		enddo
		do i=1,3
			stre(i)=s(i)+K*trtot
			stre(i+3)=s(i+3)
		enddo	

c............consistent tangent modulus update
		beta=r1-r2mu*lamda/norms
		gamma=r1/(r1+H/(r3mu))-(1-beta)
	call tenspro(sdir,sdirP)
	cc(1,1)=r4d3mu*beta+K
	cc(1,2)=K-r2d3mu*beta
	cc(1,3)=K-r2d3mu*beta
	cc(2,1)=K-r2d3mu*beta
	cc(2,2)=r4d3mu*beta+K
	cc(2,3)=K-r2d3mu*beta
	cc(3,1)=K-r2d3mu*beta
	cc(3,2)=K-r2d3mu*beta
	cc(3,3)=r4d3mu*beta+K
	cc(4,4)=mu*beta
	cc(5,5)=mu*beta
	cc(6,6)=mu*beta
	do i=1,6
		do j=1,6
			cc(i,j)=cc(i,j)-r2mu*gamma*sdirP(i,j)
		enddo
	enddo
c...........internal variables update
			alpha=alphan+lamda*sq2d3
		      
			      
			do j=1,6
			ep(j)=epn(j)+sdir(j)*lamda
			enddo
	endif
	


      return
      end
      
      
            subroutine piola2(S,dSdC,dSdCp,
     * C,invC,diCdC,Cp,invCp,diCpdCp,d)

      implicit none
      double precision 	S(6), dSdC(6,6), dSdCp(6,6),
     * invC(6), diCdC(6,6),C(6), det_s,
     * invCp(6), diCpdCp(6,6),Cp(6),
     * d(*), temp4(6,6),det2C,
     * mu, K, r2, r3, r1d2, r2d3, 
     * lamda, fact, fact1, lnj

      data r1d2, r2, r3  / 0.5d0, 2.d0 , 3.d0 /

c
c	set parameters
	  r2d3=r2/r3
	  det2C=det_s(C)**r1d2
	  lnj=log(det2C)
      mu    = d(1)
      K     = d(2)
	  lamda	= K - r2d3 * mu
	  fact1=lamda*lnj-mu

c	Compressible Neo Hookean
c	S=mu*Cp^-1+(lamda*lnJ-mu)*C^-1
        call set_s_f     (S, invCp, mu)
        call set_s_af    (S, invC, fact1)

c	dSdC=(lamda/2)C^-2+(lamda*lnJ-mu)*(-d(C^-1)/dC)
      fact = r1d2*lamda
      call set_4_f(dSdC,diCdC,fact1)
      call mult_4_ss_p(temp4,invC,invC)
	  call set_4_af(dSdC,temp4,fact)			!OKOK

c	dSdCp=mu*(-d(Cp^-1)/dCp)
	  call set_4_f (dSdCp, diCpdCp, mu)

c
	  return
	  end
