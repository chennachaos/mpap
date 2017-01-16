	subroutine s_hill_anisotropic_elastoplasticity
     * (matData,F,stre,cc,iv1,iv2,finiteInt,isw,error)
c
c	Infinitesimal/Finite Strain Elastoplasticity with Hill initial Anisotropy
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
c			 matData(8...13) - (s11,s22,s33,s12,s23,s31)/sigmaref
c
c=============================================================================

      implicit none
      double precision matData(*), F(3,*),
     * stre(*), cc(6,*), iv1(*), iv2(*),
     * Cp(6), Cpn(6), alpha, alphan, 
     * r2, r3, r1d2, liter

	integer    i, isw, finiteInt, error
	logical          finite

      data r3, r1d2, r2 / 3.d0, 0.5d0, 2.d0 /
      
         if (isw.ne.3) then
c        initialize internal variables
         call pzero(iv1,7)
         call pzero(iv2,7)
         return
         endif

        finite = (finiteInt.eq.1)
		if (finite) then
        call prgerror(1,'s_hill_anisotropic_elastoplasticity: ',
     *  'No finite option! ')
		 return     
		 endif

c	initialize Cp_n, alpha_n
	do i=1,6
	Cpn(i)=iv1(i)
	enddo
	alphan = iv1(7)
c	Cpn=Cp_n, alpha=alpha_n
	do i=1,6
	Cp(i)=Cpn(i)
	enddo
	alpha=alphan
	

      call sstrainani
     * (matData,F,stre,cc, alphan, Cpn, alpha, Cp, liter, error)

c	store converged internal variables
	do i=1,6
	iv2(i)=Cp(i)
	enddo
	iv2(7)=alpha

	return
	end
	
	
	
      subroutine sstrainani
     * (d,F,stre,cc, alphan, epn, alpha, ep,liter, error)

      implicit none

      double precision d(*), F(3,*), stre(*), cc(6,*), 
     * mu, K , lamb, eps(6), ep(6), epn(6),Pdev(6,6),
     * s(6),tstre(6),sigmay,sigma0,dsigmayda,H,g,inorm, 
     * norm, dnorm(6), d2norm(6,6),alpha, alphan,dot_ss,
     * ce(6,6), Rs(6), dRsds(6,6), dRsdg(6), idRsds(6,6), 
     * Rg, dRgds(6), dRgdg,
     * Rnorm, tol, R(7), dR(7,7), dvar(7),
     * stmp1(6),stmp2(6),tnorm,
     * r1d2, r0, r1, r2, r3, r1d3, r2d3, 
     * sq2d3, r2mu, r3mu,a44,a55,a66,a12,a23,a31,
     * s11,s22,s33,s12,s23,s31, tmp41(6,6),liter,
     * Rnorm0, Rnorm1, fact1, fact, Vnorm0(7),RFact, Vnorm1(7)


      integer          j, i, iter, maxiter, P(7), P6(6),htype, error

      logical lsearch

	lsearch=.true.

c
      data r0, r1d2, r1, r2/ 0.d0, 0.5d0, 1.d0, 2.d0 /
      r3 = 3.d0 
      r2d3 = 2.d0 / 3.d0
      r1d3 = 1.d0 / 3.d0
    	sq2d3=r2d3**r1d2

	  do i=1,6
		do j=1,6
			ce(i,j)=0.d0
			cc(i,j)=0.d0
			dRsds(i,j)=0.d0
			idRsds(i,j)=0.d0
			d2norm(i,j)=0.d0
		enddo
	  stre(i)=0.d0
	  s(i)=0.d0
	  stmp1(i)=0.d0
	  stmp2(i)=0.d0
	  Rs(i)=0.d0
	  dRsdg(i)=0.d0
	  dnorm(i)=0.d0
	  enddo
c
      mu     = d(1)
      K      = d(2)
	  sigma0=d(3)
	  H=d(4)
      tol     = d(5)
      maxiter = nint(d(6))
      htype = nint(d(7))
	  s11= 1.d0 /d(8)**r2
	  s22= 1.d0 /d(9)**r2
	  s33= 1.d0 /d(10)**r2
	  s12= 1.d0 /d(11)**r2
	  s23= 1.d0 /d(12)**r2
	  s31= 1.d0 /d(13)**r2
	  
c
      lamb   = K - r2d3 * mu
	r2mu=r2*mu
	r3mu=r3*mu
c
	a12=+s11+s22-s33
	a23=-s11+s22+s33
	a31=+s11-s22+s33
	a44=s12!*r1d3
	a55=s23!*r1d3
	a66=s31!*r1d3

c....deviatoric projection
	call unit_4_p (Pdev)
	Pdev(1,1)=r1d3*(a12+a31)
	Pdev(2,2)=r1d3*(a23+a12)
	Pdev(3,3)=r1d3*(a31+a23)
	Pdev(1,2)=-r1d3*a12
	Pdev(1,3)=-r1d3*a31
	Pdev(2,3)=-r1d3*a23
	Pdev(2,1)=Pdev(1,2)
	Pdev(3,1)=Pdev(1,3)
	Pdev(3,2)=Pdev(2,3)
	Pdev(4,4)=r1d2*a44
	Pdev(5,5)=r1d2*a55
	Pdev(6,6)=r1d2*a66

c     elastic tangent
	do i=1,3
		do j=1,3
			ce(i,j)=lamb
		enddo
	ce(i,i)=ce(i,i)+r2mu
	ce(i+3,i+3)=ce(i+3,i+3)+mu
	enddo

c     strain tensor
      eps(1) = F(1,1) - r1
      eps(2) = F(2,2) - r1
      eps(3) = F(3,3) - r1
      eps(4) = r1d2 * (F(1,2) + F(2,1))
      eps(5) = r1d2 * (F(2,3) + F(3,2))
      eps(6) = r1d2 * (F(3,1) + F(1,3))
c	
c	call set_s_am (eps,epn)
	do i=1,6
	eps(i)=eps(i)-epn(i)
	enddo

c.....calculate trial stress 

      call multx_s_4s_p(tstre,ce,eps)

c....initialize

	g=0.d0
	call set_s_p (stre,tstre)
			
c.....calculate trial norm
	call multx_s_s4_p (s,stre,Pdev)
	tnorm=(dot_ss(s,stre))**r1d2

c.....Hill yield criterion
	call hard(sigmay, dsigmayda, sigma0, H, alpha, htype)
	Rg=tnorm-sq2d3*sigmay
c
	if (Rg.le.r0) then
	call set_4_p (cc,ce)
c
	else
c
	do iter=1,maxiter

c
c	norm
	call multx_s_s4_p (s,stre,Pdev)
	norm=(dot_ss(s,stre))**r1d2
	inorm=1.d0/norm
c
c	dnorm/dstre
	call multx_s_4s_f (dnorm,Pdev,stre,inorm)
c
c	d^2norm/dstre^2
	call set_4_f (d2norm,Pdev,inorm)
	call mult_4_ss_af (d2norm,dnorm,dnorm,-inorm)
c
c
c	dRsdg=Cdnorm
	call multx_s_4s_p (dRsdg,ce,dnorm)

c	Rs=stre-tstre+gCdnorm
	call set_s_f (Rs,dRsdg,g)
	call set_s_ap (Rs,stre)
	call set_s_am (Rs,tstre)
c
c	dRsds=I+gCd2norm
	call unit_4_p (dRsds)
	call mult_4_44_af (dRsds,ce,d2norm,g)
c
c	Rg
	call hard(sigmay, dsigmayda, sigma0, H, alpha, htype)
	Rg=norm-sq2d3*sigmay
c
c	dRgds=dnorm
c
c	dRgdg
	dRgdg=-r2d3*dsigmayda
c
c	Construct and solve the matrix system for stres, g:
c	dR.dvar=R
	do i=1,6
		do j=1,6
		dR(i,j)=dRsds(i,j)
		enddo
	enddo
	do i=1,6
		dR(i,7)=dRsdg(i)
	enddo
	do i=1,6
		dR(7,i)=dnorm(i)
	enddo
		dR(7,7)=dRgdg
	do i=1,6
		R(i)=-Rs(i)
	enddo	
		R(7)=-Rg

c	check convergence on normR
	Rnorm=0.d0
	do i=1,7
	Rnorm=Rnorm+R(i)**r2
	enddo
	Rnorm=(Rnorm**r1d2)

	if(lsearch) then
		do i=1,7
		Vnorm1(i)=R(i)
		enddo
	endif

      call decompLR_matrix(dR,P,7,0)
	if (Rnorm.le.tol) goto 200	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>converged
	
c	solve the matrix system for dstress, dg:
	call solve_matrix   (dR,P,R,dvar,7)

          do i=4, 6
            dvar(i) = r1d2 * dvar(i)
          enddo

	fact=1.d0

	if(lsearch) then

	if (iter.le.2) then 
		Rfact=0.d0
		do i=1,7
		Vnorm0(i)=Vnorm1(i)
		enddo
		goto 300
	endif

		Rnorm0=0.d0
		do i=1,7
		Rnorm0=Rnorm0+dvar(i)*Vnorm0(i)
		enddo

		Rnorm1=0.d0
		do i=1,7
		Rnorm1=Rnorm1+dvar(i)*Vnorm1(i)
		enddo

		fact1=Rnorm0/Rnorm1

		if (abs(1.d0/fact1).lt.0.501) then
		Rfact=0.d0
		fact=1.d0
		do i=1,7
		Vnorm0(i)=Vnorm1(i)
		enddo
		goto 300

		else
	

		 if (fact1.lt.0.d0) then
		 fact=r1d2*fact1+((r1d2*fact1)**r2-fact1)**r1d2

		 else
		 fact=r1d2*fact1
		 endif
		Rfact=(1-fact)*Rnorm0+Rnorm1*fact**r2
		endif
		 if (abs(fact).lt.0.05)fact=r1d2	
	endif
300	continue

c	update stress and g
          g = g+fact*dvar(7)
          call set_s_af(stre,dvar,fact)
	alpha=alphan+sq2d3*g

	if (iter.eq.maxiter) then
c          write(*,*)'sf_vonmises_isotropic_elastoplasticity',
c    *                  ' Local NR iteration failed!'
	error=1
	return
	endif
	
	enddo	!NEWTON!!!!!!!!!!!!!!!!!!!!!!!!!!!!	ITERATION end


200	continue
	if(g.lt.0.d0) 
     * call prgerror(2,'sstrainani: ',
     * 'Negative Plastic Multiplier! ')

c............internal variables update

	liter=float(iter)
	call set_s_p (ep,epn)
	call set_s_af (ep,dnorm,g)

c............consistent tangent modulus update

      call decompLR_matrix(dRsds,P6,6,0)
	call inv_tensor(dRsds,idRsds,P6,6)

c	nnorm

	call mult_4_44_p (tmp41,idRsds,ce)

	call mult_s_s4_p (s,dnorm,tmp41)
	norm=dot_ss(s,dnorm)+r2d3*H
	inorm=-1.d0/norm

	call mult_s_4s_p (stmp1,tmp41,dnorm)
	call set_4_p (cc,tmp41)
	call mult_4_ss_af (cc, stmp1, s, inorm)

	endif

      return
      end




