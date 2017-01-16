	   subroutine hard(sigmay, dsigmayda, sigma0, H, alpha, htype)
      implicit none
      double precision 	sigmay, dsigmayda, sigma0, H, alpha,
     * sigmainf, del,
     * r1, r2, r3, r1d2

	  integer htype

      data r3, r1d2, r1, r2 / 3.d0, 0.5d0, 1.d0, 2.d0 /

     	if(htype.gt.1) goto 555
     	
c	Linear Hardening
	  sigmay=sigma0+H*alpha
	  dsigmayda=H	 
	  return				 
555	  continue
c	Exponential Hardening
	  sigmainf=0.715
	  del=-16.93
	  sigmay=sigma0+(sigmainf-sigma0)*(r1-exp(del*alpha))+H*alpha
	  dsigmayda=(sigma0-sigmainf)*del*exp(del*alpha)+H					 


	  return
	  end
