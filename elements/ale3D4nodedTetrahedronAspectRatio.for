
      integer function ale3D4nodedTetrahedronAspectRatio
     *                    (ALEdata,xl,s,p,ndm,nst,ndat)
      implicit none

      integer          ndm, nst, ndat,nen,nedge,jedge,
     *                 i, j, ii, i0, j0, jj, iedge,
     *                  nface, pow, err

      double precision ALEdata(*), xl(ndm,*), s(nst,*), p(*), 
     *                x(3,4), x0(3,4), r(3), k(3,3), 
     *                dedge(6,3,4), ddedge(6,3,4,3,4), 
     *                dW(6),ddW(6,6),
     *                A, dA(6),ddA(6,6),fact,rodri, 
     *                V, dV(6),ddV(6,6),tmp(6),
     *                D, dD(6),ddD(6,6),
     *                ri,dri(6), ddri(6,6),
     *                ro,dro(6), ddro(6,6),
     *                r1, r2, r3, r6

      data r1, r2, r3, r6 / 1.d0 , 2.d0, 3.d0 , 6.d0 /

c
c     LINEAR TETRAHEDRA,
c
c     OPTIMISE MESH BY MINIMISING   W = Sum (ro/ri)^2
c 
      
      err=0
      
      pow = 2

      nedge=6
      nface=4
      nen=4
      

      do i=1, nen
        do j=1, ndm
          x(j,i) = xl(j,i) + xl(j,i+nen)
          x0(j,i)=           xl(j,i+nen+nen)
        enddo
      enddo

 
         call comp_ADV(x,ndm,nen,nedge, nface,dedge, ddedge,
     *   A,dA,ddA,D,dD,ddD,V,dV,ddV,err)
     
       
     
      ro=D/(r6*V)
      ri=r3*V/A
      rodri=ro/ri


      do iedge=1, nedge 
      dro(iedge)=r1/(r6*V)*(dD(iedge)-D/V*dV(iedge))
      dri(iedge)=r3/A*(dV(iedge)-V/A*dA(iedge))
      enddo
      
      call tenspro6f(ddro, dV, dro, -r1/V)
      do iedge=1, nedge 
      tmp(iedge)=(V*dD(iedge)-D*dV(iedge))/V**r2
      enddo
      call tenspro6mf(ddro, ddro, tmp, dV, r1/(r6*V))
      
      call tenspro6f(ddri, dA, dri, -r1/A)
      do iedge=1, nedge 
      tmp(iedge)=(A*dV(iedge)-V*dA(iedge))/A**r2
      enddo
      call tenspro6mf(ddri, ddri, tmp, dA, r3/A)
      
      do iedge=1, nedge 
      do jedge=1, nedge
      ddro(iedge,jedge)=ddro(iedge,jedge)+
     *   r1/(r6*V)*(ddD(iedge,jedge)-D/V*ddV(iedge,jedge))
      ddri(iedge,jedge)=ddri(iedge,jedge)+
     *   r3/A*(ddV(iedge,jedge)-V/A*ddA(iedge,jedge))
      enddo
      enddo
      

      call comp_dW3d(ri,ro,dri,dro,ddri,ddro, pow,dW,ddW)
      

       do i=1,nen     
        i0 = (i-1)*ndm
        
        call pzero(r,ndm)
        do iedge=1, nedge
          r(1) = r(1) + dW(iedge) * dedge(iedge,1,i)
          r(2) = r(2) + dW(iedge) * dedge(iedge,2,i)
          r(3) = r(3) + dW(iedge) * dedge(iedge,3,i)
        enddo
 !          write(*,'(3(1x,g12.5)/)') (r(j),j=1,3)
 !         if (abs(r(1)).gt.10000.0.or.
 !    *        abs(r(2)).gt.10000.0.or.
 !    *        abs(r(3)).gt.10000.0)    err=-1


        do ii=1, ndm
          p(i0+ii) = p(i0+ii) - r(ii)
        enddo

        do j=1, nen
          j0 = (j-1)*ndm
          
          call pzero(k,ndm*ndm)
          do iedge=1, nedge 
            do jedge=1, nedge
              fact = ddW(iedge,jedge) * dedge(iedge,1,i)
               k(1,1) = k(1,1) + fact * dedge(jedge,1,j)
               k(1,2) = k(1,2) + fact * dedge(jedge,2,j)
               k(1,3) = k(1,3) + fact * dedge(jedge,3,j)
              fact = ddW(iedge,jedge) * dedge(iedge,2,i) 
               k(2,1) = k(2,1) + fact * dedge(jedge,1,j) 
               k(2,2) = k(2,2) + fact * dedge(jedge,2,j) 
               k(2,3) = k(2,3) + fact * dedge(jedge,3,j)
              fact = ddW(iedge,jedge) * dedge(iedge,3,i) 
               k(3,1) = k(3,1) + fact * dedge(jedge,1,j) 
               k(3,2) = k(3,2) + fact * dedge(jedge,2,j) 
               k(3,3) = k(3,3) + fact * dedge(jedge,3,j)
            enddo
            k(1,1) = k(1,1)+ dW(iedge) * ddedge(iedge,1,i,1,j)
            k(1,2) = k(1,2)+ dW(iedge) * ddedge(iedge,1,i,2,j)
            k(1,3) = k(1,3)+ dW(iedge) * ddedge(iedge,1,i,3,j)

            k(2,1) = k(2,1)+ dW(iedge) * ddedge(iedge,2,i,1,j)
            k(2,2) = k(2,2)+ dW(iedge) * ddedge(iedge,2,i,2,j)
            k(2,3) = k(2,3)+ dW(iedge) * ddedge(iedge,2,i,3,j)

            k(3,1) = k(3,1)+ dW(iedge) * ddedge(iedge,3,i,1,j)
            k(3,2) = k(3,2)+ dW(iedge) * ddedge(iedge,3,i,2,j)
            k(3,3) = k(3,3)+ dW(iedge) * ddedge(iedge,3,i,3,j)
          enddo

          do ii=1, ndm
            do jj=1, ndm
              s(i0+ii,j0+jj) = s(i0+ii,j0+jj) + k(ii,jj)
            enddo
          enddo
        enddo

      enddo

      ale3D4nodedTetrahedronAspectRatio = err
      
      if(err.ne.0) then 
      write(*,*) "err.ne.0", err
      write(*,*) "W=", rodri**r2
      write(*,'(6(1x,g12.5)/)') (dW(j),j=1,6)
      write(*,*) "ri=", ri
      write(*,'(6(1x,g12.5)/)') (dri(j),j=1,6)
      write(*,*) "ro=", ro
      write(*,'(6(1x,g12.5)/)') (dro(j),j=1,6)
      write(*,*) "V=", V
      write(*,'(6(1x,g12.5)/)') (dV(j),j=1,6)
      write(*,*) "A=", A
      write(*,'(6(1x,g12.5)/)') (dA(j),j=1,6)
      write(*,*) "D=", D
      write(*,'(6(1x,g12.5)/)') (dD(j),j=1,6)
      endif


      return

      end




c===========================================================================
       subroutine comp_ADV(x,ndm,nen,nedge,nface, dedge, ddedge,
     *   A,dA,ddA,D,dD,ddD,V,dV,ddV,err)
c===========================================================================

      implicit none

      integer          nedge,ndm,nen,
     *                 i, j, iedge,
     *                 nface, DAT(4,3), iface,
     *                 err, errA, errD, errV

      double precision x(3,4), edge(6), 
     *                  dedge(6,3,4), 
     *                 ddedge(6,3,4,3,4), 
     *                A, dA(6),ddA(6,6),
     *                V, dV(6),ddV(6,6),
     *                D, dD(6),ddD(6,6),
     *                tmpA, tmpdA(3),tmpddA(3,3), 
     *                 r0, r1, r2

      data r0 , r1, r2 / 0.d0 , 1.d0 , 2.d0/
      
      
   
         DAT(1,1)=1
         DAT(1,2)=2
         DAT(1,3)=3
         DAT(2,1)=1
         DAT(2,2)=5
         DAT(2,3)=6
         DAT(3,1)=2
         DAT(3,2)=6
         DAT(3,3)=4
         DAT(4,1)=3
         DAT(4,2)=4
         DAT(4,3)=5

         
      

      call pzero(dedge,nedge*ndm*nen)
      call pzero(ddedge,nedge*ndm*nen*ndm*nen)
      
        do iedge=1, nedge  
        
        if (iedge.gt.(nen-1)) then
        i = nen 
        j = iedge-(nen-1)
        else if (iedge.eq.1) then
        i = 2 
        j = 3
        else if (iedge.eq.2) then
        i = 3 
        j = 1
        else
        i = 1
        j = 2
        endif


       

        edge(iedge) =  sqrt((x(1,i)-x(1,j))**r2
     *                    + (x(2,i)-x(2,j))**r2
     *                    + (x(3,i)-x(3,j))**r2)
     
        dedge (iedge,1,i) =   (x(1,i)-x(1,j))/edge(iedge)
        dedge (iedge,1,j) = - (x(1,i)-x(1,j))/edge(iedge)
        dedge (iedge,2,i) =   (x(2,i)-x(2,j))/edge(iedge)
        dedge (iedge,2,j) = - (x(2,i)-x(2,j))/edge(iedge)
        dedge (iedge,3,i) =   (x(3,i)-x(3,j))/edge(iedge)
        dedge (iedge,3,j) = - (x(3,i)-x(3,j))/edge(iedge)

        ddedge(iedge,1,i,1,i) =( r1-dedge(iedge,1,i)*dedge(iedge,1,i))/
     *   edge(iedge)
        ddedge(iedge,1,i,1,j) =(-r1-dedge(iedge,1,i)*dedge(iedge,1,j))/
     *   edge(iedge)
        ddedge(iedge,1,i,2,i) =( r0-dedge(iedge,1,i)*dedge(iedge,2,i))/
     *   edge(iedge)
        ddedge(iedge,1,i,2,j) =( r0-dedge(iedge,1,i)*dedge(iedge,2,j))/
     *   edge(iedge)
        ddedge(iedge,1,i,3,i) =( r0-dedge(iedge,1,i)*dedge(iedge,3,i))/
     *   edge(iedge)
        ddedge(iedge,1,i,3,j) =( r0-dedge(iedge,1,i)*dedge(iedge,3,j))/
     *   edge(iedge)
     
        ddedge(iedge,1,j,1,i) =(-r1-dedge(iedge,1,j)*dedge(iedge,1,i))/
     *   edge(iedge)
        ddedge(iedge,1,j,1,j) =( r1-dedge(iedge,1,j)*dedge(iedge,1,j))/
     *   edge(iedge)
        ddedge(iedge,1,j,2,i) =( r0-dedge(iedge,1,j)*dedge(iedge,2,i))/
     *   edge(iedge)
        ddedge(iedge,1,j,2,j) =( r0-dedge(iedge,1,j)*dedge(iedge,2,j))/
     *   edge(iedge)
        ddedge(iedge,1,j,3,i) =( r0-dedge(iedge,1,j)*dedge(iedge,3,i))/
     *   edge(iedge)
        ddedge(iedge,1,j,3,j) =( r0-dedge(iedge,1,j)*dedge(iedge,3,j))/
     *   edge(iedge)
     
        ddedge(iedge,2,i,1,i) =( r0-dedge(iedge,2,i)*dedge(iedge,1,i))/
     *   edge(iedge)
        ddedge(iedge,2,i,1,j) =( r0-dedge(iedge,2,i)*dedge(iedge,1,j))/
     *   edge(iedge)
        ddedge(iedge,2,i,2,i) =( r1-dedge(iedge,2,i)*dedge(iedge,2,i))/
     *   edge(iedge)
        ddedge(iedge,2,i,2,j) =(-r1-dedge(iedge,2,i)*dedge(iedge,2,j))/
     *   edge(iedge)
        ddedge(iedge,2,i,3,i) =( r0-dedge(iedge,2,i)*dedge(iedge,3,i))/
     *   edge(iedge)
        ddedge(iedge,2,i,3,j) =( r0-dedge(iedge,2,i)*dedge(iedge,3,j))/
     *   edge(iedge)
     
        ddedge(iedge,2,j,1,i) =( r0-dedge(iedge,2,j)*dedge(iedge,1,i))/
     *   edge(iedge)
        ddedge(iedge,2,j,1,j) =( r0-dedge(iedge,2,j)*dedge(iedge,1,j))/
     *   edge(iedge)
        ddedge(iedge,2,j,2,i) =(-r1-dedge(iedge,2,j)*dedge(iedge,2,i))/
     *   edge(iedge)
        ddedge(iedge,2,j,2,j) =( r1-dedge(iedge,2,j)*dedge(iedge,2,j))/
     *   edge(iedge)
        ddedge(iedge,2,j,3,i) =( r0-dedge(iedge,2,j)*dedge(iedge,3,i))/
     *   edge(iedge)
        ddedge(iedge,2,j,3,j) =( r0-dedge(iedge,2,j)*dedge(iedge,3,j))/
     *   edge(iedge)
     
        ddedge(iedge,3,i,1,i) =( r0-dedge(iedge,3,i)*dedge(iedge,1,i))/
     *   edge(iedge)
        ddedge(iedge,3,i,1,j) =( r0-dedge(iedge,3,i)*dedge(iedge,1,j))/
     *   edge(iedge)
        ddedge(iedge,3,i,2,i) =( r0-dedge(iedge,3,i)*dedge(iedge,2,i))/
     *   edge(iedge)
        ddedge(iedge,3,i,2,j) =( r0-dedge(iedge,3,i)*dedge(iedge,2,j))/
     *   edge(iedge)
        ddedge(iedge,3,i,3,i) =( r1-dedge(iedge,3,i)*dedge(iedge,3,i))/
     *   edge(iedge)
        ddedge(iedge,3,i,3,j) =(-r1-dedge(iedge,3,i)*dedge(iedge,3,j))/
     *   edge(iedge)
     
        ddedge(iedge,3,j,1,i) =( r0-dedge(iedge,3,j)*dedge(iedge,1,i))/
     *   edge(iedge)
        ddedge(iedge,3,j,1,j) =( r0-dedge(iedge,3,j)*dedge(iedge,1,j))/
     *   edge(iedge)
        ddedge(iedge,3,j,2,i) =( r0-dedge(iedge,3,j)*dedge(iedge,2,i))/
     *   edge(iedge)
        ddedge(iedge,3,j,2,j) =( r0-dedge(iedge,3,j)*dedge(iedge,2,j))/
     *   edge(iedge)
        ddedge(iedge,3,j,3,i) =(-r1-dedge(iedge,3,j)*dedge(iedge,3,i))/
     *   edge(iedge)
        ddedge(iedge,3,j,3,j) =( r1-dedge(iedge,3,j)*dedge(iedge,3,j))/
     *   edge(iedge)
        
      enddo


      

c       Compute Surface Area
      A=r0
      call pzero(dA,nedge)
      call pzero(ddA,nedge*nedge)

      do iface=1, nface
      
      call comp_tA (edge(DAT(iface,1)),edge(DAT(iface,2)),
     *                edge(DAT(iface,3)),
     *                  tmpA, tmpdA, tmpddA,errA)
     
      if(errA.ne.0) err=errA
     
      A=A+tmpA
       do i=1, 3
      dA(DAT(iface,i))=dA(DAT(iface,i))+tmpdA(i)
         do j=1, 3
          ddA(DAT(iface,i),DAT(iface,j))=
     *    ddA(DAT(iface,i),DAT(iface,j)) + tmpddA(i,j)
         enddo
       enddo
      enddo
      

      
!Compute Delta(area of triangle with edges ad, be, cf)
      call  comp_tD(edge,D,dD,ddD,errD)
      if(errD.ne.0) err=errD
     
!Compute Volume
!      call pzero(dV,6)
!      call pzero(ddV,36)
!      v=1
      call comp_tV2(edge,V,dV,ddV,errV)
      if(errV.ne.0) err=errV
      
      return

      end


c===========================================================================
      subroutine comp_dW3d(ri,ro,dri,dro,ddri,ddro,pow,dWds,ddWds)
c---------------------------------------------------------------------------
      implicit none
      integer          i, j, pow
      double precision ri, ro, dri(6), dro(6), ddri(6,6), ddro(6,6), 
     *                 dWds(6), ddWds(6,6), fact, fact1, fact2, fact3,
     *                 fact6, W,  r1

      r1 = 1.d0 

      fact6 = - r1 / ri
      fact  = - ro * fact6
      W     = fact
      do i=2, pow
        W = W * fact
      enddo
      

      fact1 = W / ro
      fact2 = W * fact6

      fact3 = (float(pow) - r1) / W


      do i=1, 6
        dWds(i) =       fact1 * dro(i) + fact2 * dri(i)
        do j=1, i
          ddWds(i,j) =  fact3 * dWds(i) * dWds(j)
     *                + fact1 * ddro(i,j) + fact2 * ddri(i,j) 
     *                + fact6 * (dWds(i) * dri(j) + dWds(j) * dri(i))
        enddo
      enddo
      ddWds(1,2) = ddWds(2,1)
      ddWds(1,3) = ddWds(3,1)
      ddWds(1,4) = ddWds(4,1)
      ddWds(1,5) = ddWds(5,1)
      ddWds(1,6) = ddWds(6,1)
      ddWds(2,3) = ddWds(3,2)
      ddWds(2,4) = ddWds(4,2)
      ddWds(2,5) = ddWds(5,2)
      ddWds(2,6) = ddWds(6,2)
      ddWds(3,4) = ddWds(4,3)
      ddWds(3,5) = ddWds(5,3)
      ddWds(3,6) = ddWds(6,3)
      ddWds(4,5) = ddWds(5,4)
      ddWds(4,6) = ddWds(6,4)
      ddWds(5,6) = ddWds(6,5)

      return
      end

c===========================================================================

      




    
     

      
      
c===========================================================================
      subroutine comp_tA(s1,s2,s3,A,dA,ddA,err)
c---------------------------------------------------------------------------
      implicit none
      double precision s1,s2,s3, A, dA(3),ddA(3,3),perid2,
     *                 A2, fact, sgn,
     *                 r2, r3, r4, r6, r8
     
      integer err


      data r2, r3, r4, r6, r8
     *   / 2.d0, 3.d0, 4.d0, 6.d0, 8.d0 /
      

      perid2=(s1+s2+s3)/r2
      A2=perid2*(perid2-s1)*(perid2-s2)*(perid2-s3)
      
      err=0
       if (abs(A2).lt.1.d-14) then
       write(*,*) "A<0,  s1:", s1," s2:",s2," s3:",s3
       err=-1
       return
       endif
     
      sgn=A2/abs(A2)
      
      A=sqrt(abs(A2))

      fact=sgn/(r8*A)
      

       dA(1) = fact*s1*(-s1**r2 + s2**r2 + s3**r2)
       dA(2) = fact*s2*( s1**r2 - s2**r2 + s3**r2)
       dA(3) = fact*s3*( s1**r2 + s2**r2 - s3**r2)
       
       fact=sgn/(128.*A**r3)
       
        ddA(1,1) = fact * (s1**r6 + r3*s1**r2*(s2**r2 - s3**r2)**r2 - 
     *   r3*s1**r4*(s2**r2 + s3**r2) - (s2**r2 - s3**r2)**r2*(s2**r2 + 
     *   s3**r2))
        ddA(2,1) = fact * (r4*s1*s2*s3**r2*(s1**r2 + s2**r2 - s3**r2))
        ddA(3,1) = fact * (r4*s1*s2**r2*s3*(s1**r2 - s2**r2 + s3**r2))

        ddA(1,2) = fact * (r4*s1*s2*s3**r2*(s1**r2 + s2**r2 - s3**r2))
        ddA(2,2) = fact * (-s1**r6 + (s2**r2 - s3**r2)**r3 + 
     *   s1**r4*(r3*s2**r2 + s3**r2) + s1**r2*(-r3*s2**r4 - 
     *   r6*s2**r2*s3**r2 + s3**r4))
        ddA(3,2) = fact * (r4*s1**r2*s2*s3*(-s1**r2 + s2**r2 + s3**r2))
             
        ddA(1,3) = fact * (r4*s1*s2**r2*s3*(s1**r2 - s2**r2 + s3**r2))
        ddA(2,3) = fact * (r4*s1**r2*s2*s3*(-s1**r2 + s2**r2 + s3**r2))
        ddA(3,3) = fact * (-s1**r6 - (s2**r2 - s3**r2)**r3 + 
     *    s1**r4*(s2**r2 + r3*s3**r2) + s1**r2*(s2**r4 - 
     *    r6*s2**r2*s3**r2 - r3*s3**r4))

      return
      end
      
      
      
      
      
      
c===========================================================================
      subroutine comp_tD(s,D,dD,ddD,err)
c---------------------------------------------------------------------------
      implicit none
      double precision s(6), D, dD(6), ddD(6,6),  fact,   
     *                 r2, r3, r4, r6, r8, d2,
     *                 r1d2, sgn
    
      integer            err
      data  r2,   r3,   r4,   r6,   r1d2, r8
     *   / 2.d0, 3.d0, 4.d0, 6.d0, 0.5d0, 8.d0 /
     
     
      fact=(s(1)*s(4)+s(2)*s(5)+s(3)*s(6))*r1d2
      D2=fact*(fact-s(1)*s(4))*(fact-s(2)*s(5))
     *                *(fact-s(3)*s(6))

      
      err=0
       if (abs(D2).lt.1.d-14) then
        write(*,*) "D<0,  s1:", s(1)*s(4)," s2:",
     *     s(2)*s(5)," s3:",s(3)*s(6)
       err=-1
       return
       endif
     
      D=sqrt(abs(D2))
 
      sgn=D2/abs(D2)
      
      fact=sgn/(r8*D)
      
      dD(1) = fact * (s(1)*s(4)**r2*(-(s(1)**r2*s(4)**r2) 
     * + s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
     
      dD(2) = fact * (s(2)*s(5)**r2*(s(1)**r2*s(4)**r2 
     * - s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
     
      dD(3) = fact * (s(3)*s(6)**r2*(s(1)**r2*s(4)**r2 
     * + s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2))
     
      dD(4) = fact * (s(1)**r2*s(4)*(-(s(1)**r2*s(4)**r2) 
     * + s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
     
      dD(5) = fact * (s(2)**r2*s(5)*(s(1)**r2*s(4)**r2 
     * - s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
     
      dD(6) = fact * (s(3)**r2*s(6)*(s(1)**r2*s(4)**r2 
     * + s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2))
     
     
      fact=sgn/(128.*D**r3)

       ddD(1,1) = fact * (s(4)**r2*((s(1)**r2*s(4)**r2 
     *  - s(2)**r2*s(5)**r2)**r3 - s(3)**r2*(3*s(1)**r4*s(4)**r4 
     *  + r6*s(1)**r2*s(2)**r2*s(4)**r2*s(5)**r2  
     *  - s(2)**r4*s(5)**r4)*s(6)**r2 + s(3)**r4*(3*s(1)**r2*s(4)**r2 
     *  + s(2)**r2*s(5)**r2)*s(6)**r4 - s(3)**r6*s(6)**r6))
            
       ddD(2,1) = fact * (r4*s(1)*s(2)*s(3)**r2*s(4)**r2*s(5)**r2
     *  *s(6)**r2*(s(1)**r2*s(4)**r2 
     *  + s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2))
            
       ddD(3,1) = fact * (r4*s(1)*s(2)**r2*s(3)*s(4)**r2*s(5)**r2
     *  *s(6)**r2*(s(1)**r2*s(4)**r2 
     *  - s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
            
       ddD(4,1) = fact * (r2*s(1)*s(4)*((s(1)**r2*s(4)**r2 
     *  - s(2)**r2*s(5)**r2)**r3 
     *       - s(3)**r2*(r3*s(1)**r4*s(4)**r4 
     * + r2*s(1)**r2*s(2)**r2*s(4)**r2*s(5)**r2 
     * - s(2)**r4*s(5)**r4)*s(6)**r2 + s(3)**r4*(r3*s(1)**r2*s(4)**r2 
     * + s(2)**r2*s(5)**r2)*s(6)**r4 
     *  - s(3)**r6*s(6)**r6))
            
       ddD(5,1) = fact * (r4*s(1)*s(2)**r2*s(3)**r2*s(4)**r2
     *  *s(5)*s(6)**r2*(s(1)**r2*s(4)**r2 
     *   + s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2))
     
       ddD(6,1) = fact * (r4*s(1)*s(2)**r2*s(3)**r2*s(4)**r2*s(5)**r2
     *  *s(6)*(s(1)**r2*s(4)**r2 
     *   - s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
     
      
             
       ddD(1,2) = fact * (r4*s(1)*s(2)*s(3)**r2*s(4)**r2*s(5)**r2
     *  *s(6)**r2 *(s(1)**r2*s(4)**r2 + s(2)**r2*s(5)**r2 
     *  - s(3)**r2*s(6)**r2))
            
       ddD(2,2) = fact * (s(5)**r2*(-(s(1)**r6*s(4)**r6) 
     *  + (s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2)**r3 
     *  + s(1)**r4*s(4)**r4*(r3*s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2) 
     *  + s(1)**r2*s(4)**r2*(-r3*s(2)**r4*s(5)**r4 
     *  - r6*s(2)**r2*s(3)**r2*s(5)**r2*s(6)**r2 + s(3)**r4*s(6)**r4)))
            
       ddD(3,2) = fact * (r4*s(1)**r2*s(2)*s(3)*s(4)**r2*s(5)**r2
     *  *s(6)**r2*(-(s(1)**r2*s(4)**r2) + s(2)**r2*s(5)**r2 
     *  + s(3)**r2*s(6)**r2))
     
       ddD(4,2) = fact * (r4*s(1)**r2*s(2)*s(3)**r2*s(4)*s(5)**r2
     *  *s(6)**r2*(s(1)**r2*s(4)**r2 + s(2)**r2*s(5)**r2 
     *  - s(3)**r2*s(6)**r2))
     
       ddD(5,2) = fact * (r2*s(2)*s(5)*(-(s(1)**r6*s(4)**r6) 
     *  + (s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2)**r3 
     *  + s(1)**r4*s(4)**r4*(r3*s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2) 
     *  + s(1)**r2*s(4)**r2*(-r3*s(2)**r4*s(5)**r4 
     *  - r2*s(2)**r2*s(3)**r2*s(5)**r2*s(6)**r2 + s(3)**r4*s(6)**r4)))
     
       ddD(6,2) = fact * (r4*s(1)**r2*s(2)*s(3)**r2*s(4)**r2*s(5)**r2
     *  *s(6)*(-(s(1)**r2*s(4)**r2) + s(2)**r2*s(5)**r2 
     *  + s(3)**r2*s(6)**r2))


  
        ddD(1,3) = fact * (r4*s(1)*s(2)**r2*s(3)*s(4)**r2*s(5)**r2
     *  *s(6)**r2*(s(1)**r2*s(4)**r2 
     *  - s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
        
        ddD(2,3) = fact * (r4*s(1)**r2*s(2)*s(3)*s(4)**r2*s(5)**r2
     *   *s(6)**r2*(-(s(1)**r2*s(4)**r2) + s(2)**r2*s(5)**r2 
     *  + s(3)**r2*s(6)**r2))
        
        ddD(3,3) = fact * (s(6)**r2*(-(s(1)**r6*s(4)**r6) 
     *   - (s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2)**r3 
     *  + s(1)**r4*s(4)**r4*(s(2)**r2*s(5)**r2 
     *  + r3*s(3)**r2*s(6)**r2) + s(1)**r2*s(4)**r2*(s(2)**r4*s(5)**r4 
     *  - r6*s(2)**r2*s(3)**r2*s(5)**r2*s(6)**r2 
     *  - r3*s(3)**r4*s(6)**r4)))
     
       ddD(4,3) = fact * (r4*s(1)**r2*s(2)**r2*s(3)*s(4)
     *  *s(5)**r2*s(6)**r2*(s(1)**r2*s(4)**r2 - s(2)**r2*s(5)**r2 
     *  + s(3)**r2*s(6)**r2))
     
       ddD(5,3) = fact * (r4*s(1)**r2*s(2)**r2*s(3)*s(4)**r2*s(5)
     *  *s(6)**r2*(-(s(1)**r2*s(4)**r2) + s(2)**r2*s(5)**r2 
     *  + s(3)**r2*s(6)**r2))
     
       ddD(6,3) = fact * (r2*s(3)*s(6)*(-(s(1)**r6*s(4)**r6) 
     *  - (s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2)**r3 
     *  + s(1)**r4*s(4)**r4*(s(2)**r2*s(5)**r2 + r3*s(3)**r2*s(6)**r2) 
     *  + s(1)**r2*s(4)**r2*(s(2)**r4*s(5)**r4 
     *  - r2*s(2)**r2*s(3)**r2*s(5)**r2*s(6)**r2 
     *  - r3*s(3)**r4*s(6)**r4)))
       
      
           
      ddD(1,4) = fact * (r2*s(1)*s(4)*((s(1)**r2*s(4)**r2 
     * - s(2)**r2*s(5)**r2)**r3 - s(3)**r2*(r3*s(1)**r4*s(4)**r4 
     * + r2*s(1)**r2*s(2)**r2*s(4)**r2*s(5)**r2 
     * - s(2)**r4*s(5)**r4)*s(6)**r2 
     *     + s(3)**r4*(r3*s(1)**r2*s(4)**r2 
     * + s(2)**r2*s(5)**r2)*s(6)**r4 - s(3)**r6*s(6)**r6))
           
      ddD(2,4) = fact * (r4*s(1)**r2*s(2)*s(3)**r2*s(4)
     *  *s(5)**r2*s(6)**r2*(s(1)**r2*s(4)**r2 + s(2)**r2*s(5)**r2 
     * - s(3)**r2*s(6)**r2))
     
      ddD(3,4) = fact *  (r4*s(1)**r2*s(2)**r2*s(3)*s(4)*s(5)**r2
     * *s(6)**r2*(s(1)**r2*s(4)**r2 - s(2)**r2*s(5)**r2 
     *  + s(3)**r2*s(6)**r2))
     
       ddD(4,4) = fact * (s(1)**r2*((s(1)**r2*s(4)**r2 
     *  - s(2)**r2*s(5)**r2)**r3 - s(3)**r2*(r3*s(1)**r4*s(4)**r4 
     *  + r6*s(1)**r2*s(2)**r2*s(4)**r2*s(5)**r2 
     * - s(2)**r4*s(5)**r4)*s(6)**r2 + s(3)**r4*(r3*s(1)**r2*s(4)**r2 
     *  + s(2)**r2*s(5)**r2)*s(6)**r4 - s(3)**r6*s(6)**r6))
     
       ddD(5,4) = fact * (r4*s(1)**r2*s(2)**r2*s(3)**r2*s(4)*s(5)
     *  *s(6)**r2*(s(1)**r2*s(4)**r2 + s(2)**r2*s(5)**r2 
     * - s(3)**r2*s(6)**r2))
     
       ddD(6,4) = fact * (r4*s(1)**r2*s(2)**r2*s(3)**r2*s(4)
     *  *s(5)**r2*s(6)*(s(1)**r2*s(4)**r2 
     * - s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
   
   
     
      ddD(1,5) = fact * (r4*s(1)*s(2)**r2*s(3)**r2*s(4)**r2*s(5)
     *  *s(6)**r2*(s(1)**r2*s(4)**r2 + s(2)**r2*s(5)**r2 
     * - s(3)**r2*s(6)**r2))
           
      ddD(2,5) = fact *  (r2*s(2)*s(5)*(-(s(1)**r6*s(4)**r6) 
     * + (s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2)**r3 
     * + s(1)**r4*s(4)**r4*(r3*s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2) 
     * + s(1)**r2*s(4)**r2*(-r3*s(2)**r4*s(5)**r4 
     *  - r2*s(2)**r2*s(3)**r2*s(5)**r2*s(6)**r2 + s(3)**r4*s(6)**r4)))
           
      ddD(3,5) = fact *  (r4*s(1)**r2*s(2)**r2*s(3)*s(4)**r2*s(5)
     *  *s(6)**r2*(-(s(1)**r2*s(4)**r2) + s(2)**r2*s(5)**r2 
     *  + s(3)**r2*s(6)**r2))
     
     
      ddD(4,5) = fact *  (r4*s(1)**r2*s(2)**r2*s(3)**r2*s(4)*s(5)*
     *  s(6)**r2*(s(1)**r2*s(4)**r2 + s(2)**r2*s(5)**r2 
     * - s(3)**r2*s(6)**r2))
     
      ddD(5,5) = fact *  (s(2)**r2*(-(s(1)**r6*s(4)**r6) 
     * + (s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2)**r3 
     * + s(1)**r4*s(4)**r4*(r3*s(2)**r2*s(5)**r2 
     * + s(3)**r2*s(6)**r2) + s(1)**r2*s(4)**r2*(-r3*s(2)**r4*s(5)**r4 
     * - r6*s(2)**r2*s(3)**r2*s(5)**r2*s(6)**r2 + s(3)**r4*s(6)**r4)))
     
      ddD(6,5) = fact *  (r4*s(1)**r2*s(2)**r2*s(3)**r2*s(4)**r2*s(5)
     * *s(6)*(-(s(1)**r2*s(4)**r2) + s(2)**r2*s(5)**r2 
     * + s(3)**r2*s(6)**r2))
   
     
       ddD(1,6) = fact * (r4*s(1)*s(2)**r2*s(3)**r2*s(4)**r2*s(5)**r2
     *  *s(6)*(s(1)**r2*s(4)**r2 - s(2)**r2*s(5)**r2 
     * + s(3)**r2*s(6)**r2))
       
      ddD(2,6) = fact *  (r4*s(1)**r2*s(2)*s(3)**r2*s(4)**r2*s(5)**r2
     * *s(6)*(-(s(1)**r2*s(4)**r2) + s(2)**r2*s(5)**r2 
     * + s(3)**r2*s(6)**r2))
       
      ddD(3,6) = fact *  (r2*s(3)*s(6)*(-(s(1)**r6*s(4)**r6) 
     * - (s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2)**r3 
     * + s(1)**r4*s(4)**r4*(s(2)**r2*s(5)**r2 + r3*s(3)**r2*s(6)**r2) 
     * + s(1)**r2*s(4)**r2*(s(2)**r4*s(5)**r4 
     * - r2*s(2)**r2*s(3)**r2*s(5)**r2*s(6)**r2 
     * - r3*s(3)**r4*s(6)**r4)))
     
       ddD(4,6) = fact * (r4*s(1)**r2*s(2)**r2*s(3)**r2*s(4)*s(5)**r2*
     * s(6)*(s(1)**r2*s(4)**r2 - s(2)**r2*s(5)**r2 + s(3)**r2*s(6)**r2))
      
       ddD(5,6) = fact * (r4*s(1)**r2*s(2)**r2*s(3)**r2*s(4)**r2*s(5)
     *  *s(6)*(-(s(1)**r2*s(4)**r2) + s(2)**r2*s(5)**r2 
     *  + s(3)**r2*s(6)**r2))
     
       ddD(6,6) = fact * (s(3)**r2*(-(s(1)**r6*s(4)**r6) 
     *  - (s(2)**r2*s(5)**r2 - s(3)**r2*s(6)**r2)**r3 
     * + s(1)**r4*s(4)**r4*(s(2)**r2*s(5)**r2 + r3*s(3)**r2*s(6)**r2) 
     * + s(1)**r2*s(4)**r2*(s(2)**r4*s(5)**r4 
     * - r6*s(2)**r2*s(3)**r2*s(5)**r2*s(6)**r2 
     * - r3*s(3)**r4*s(6)**r4)))
       
     
      return
      end
      
      
      
c===========================================================================
      subroutine comp_tV(x,V,dV,ddV,err)
c---------------------------------------------------------------------------
      implicit none
      double precision x(3,4), V, dV(3,4), ddV(3,4,3,4), r1d6
      integer    err, i
      
      r1d6=1.d0/6.d0
      
          V=r1d6*(-
     *    x(1,3)*x(2,2)*x(3,1) + x(1,4)*x(2,2)*x(3,1) + 
     *  x(1,2)*x(2,3)*x(3,1) - 
     *    x(1,4)*x(2,3)*x(3,1) - x(1,2)*x(2,4)*x(3,1) + 
     *  x(1,3)*x(2,4)*x(3,1) + 
     *    x(1,3)*x(2,1)*x(3,2) - x(1,4)*x(2,1)*x(3,2) - 
     *  x(1,1)*x(2,3)*x(3,2) + 
     *    x(1,4)*x(2,3)*x(3,2) + x(1,1)*x(2,4)*x(3,2) - 
     *  x(1,3)*x(2,4)*x(3,2) - 
     *    x(1,2)*x(2,1)*x(3,3) + x(1,4)*x(2,1)*x(3,3) + 
     *  x(1,1)*x(2,2)*x(3,3) - 
     *    x(1,4)*x(2,2)*x(3,3) - x(1,1)*x(2,4)*x(3,3) + 
     *  x(1,2)*x(2,4)*x(3,3) + 
     *    x(1,2)*x(2,1)*x(3,4) - x(1,3)*x(2,1)*x(3,4) - 
     *  x(1,1)*x(2,2)*x(3,4) + 
     *    x(1,3)*x(2,2)*x(3,4) + x(1,1)*x(2,3)*x(3,4) - 
     *  x(1,2)*x(2,3)*x(3,4))
     
c      err=0
c      if (abs(V).lt.1.d-14) then
c       write(*,*) "V<0"
c       do i=1, 4
c       write(*,*) x(1,i),",",x(2,i),",",x(3,i)
c       write(*,*) 
c       enddo
c       err=-1
c       return
c       endif
     
     
      if (V .lt. 0.d0) then
        V = -V
        r1d6=-r1d6
      endif
  
      
       dV(1,1) = r1d6*(-x(2,3)*x(3,2) + x(2,4)*x(3,2) + x(2,2)*x(3,3) - 
     *                  x(2,4)*x(3,3) - x(2,2)*x(3,4) + x(2,3)*x(3,4))
       dV(1,2) = r1d6*( x(2,3)*x(3,1) - x(2,4)*x(3,1) - x(2,1)*x(3,3) + 
     *                  x(2,4)*x(3,3) + x(2,1)*x(3,4) - x(2,3)*x(3,4))
       dV(1,3) = r1d6*(-x(2,2)*x(3,1) + x(2,4)*x(3,1) + x(2,1)*x(3,2) - 
     *                  x(2,4)*x(3,2) - x(2,1)*x(3,4) + x(2,2)*x(3,4))
       dV(1,4) = r1d6*( x(2,2)*x(3,1) - x(2,3)*x(3,1) - x(2,1)*x(3,2) + 
     *                  x(2,3)*x(3,2) + x(2,1)*x(3,3) - x(2,2)*x(3,3))
     
       dV(2,1) = r1d6*( x(1,3)*x(3,2) - x(1,4)*x(3,2) - x(1,2)*x(3,3) + 
     *                  x(1,4)*x(3,3) + x(1,2)*x(3,4) - x(1,3)*x(3,4))
       dV(2,2) = r1d6*(-x(1,3)*x(3,1) + x(1,4)*x(3,1) + x(1,1)*x(3,3) - 
     *                  x(1,4)*x(3,3) - x(1,1)*x(3,4) + x(1,3)*x(3,4))
       dV(2,3) = r1d6*( x(1,2)*x(3,1) - x(1,4)*x(3,1) - x(1,1)*x(3,2) + 
     *                  x(1,4)*x(3,2) + x(1,1)*x(3,4) - x(1,2)*x(3,4))
       dV(2,4) = r1d6*(-x(1,2)*x(3,1) + x(1,3)*x(3,1) + x(1,1)*x(3,2) - 
     *                  x(1,3)*x(3,2) - x(1,1)*x(3,3) + x(1,2)*x(3,3))

       dV(3,1) = r1d6*(-x(1,3)*x(2,2) + x(1,4)*x(2,2) + x(1,2)*x(2,3) -
     *                  x(1,4)*x(2,3) - x(1,2)*x(2,4) + x(1,3)*x(2,4))
       dV(3,2) = r1d6*( x(1,3)*x(2,1) - x(1,4)*x(2,1) - x(1,1)*x(2,3) +
     *                  x(1,4)*x(2,3) + x(1,1)*x(2,4) - x(1,3)*x(2,4))
       dV(3,3) = r1d6*(-x(1,2)*x(2,1) + x(1,4)*x(2,1) + x(1,1)*x(2,2) - 
     *                  x(1,4)*x(2,2) - x(1,1)*x(2,4) + x(1,2)*x(2,4))
       dV(3,4) = r1d6*( x(1,2)*x(2,1) - x(1,3)*x(2,1) - x(1,1)*x(2,2) + 
     *                  x(1,3)*x(2,2) + x(1,1)*x(2,3) - x(1,2)*x(2,3))
       
       
       call pzero(ddV,3*4*3*4)
       
       ddV(1,2,2,1) = r1d6*(-x(3,3) + x(3,4))
       ddV(1,3,2,1) = r1d6*( x(3,2) - x(3,4))
       ddV(1,4,2,1) = r1d6*(-x(3,2) + x(3,3))
       ddV(1,1,2,2) = r1d6*(x(3,3) - x(3,4))
       ddV(1,3,2,2) = r1d6*(-x(3,1) + x(3,4))
       ddV(1,4,2,2) = r1d6*(x(3,1) - x(3,3))
       ddV(1,1,2,3) = r1d6*(-x(3,2) + x(3,4))
       ddV(1,2,2,3) = r1d6*(x(3,1) - x(3,4))
       ddV(1,4,2,3) = r1d6*(-x(3,1) + x(3,2))
       ddV(1,1,2,4) = r1d6*(x(3,2) - x(3,3))
       ddV(1,2,2,4) = r1d6*(-x(3,1) + x(3,3))
       ddV(1,3,2,4) = r1d6*(x(3,1) - x(3,2))
       ddV(1,2,3,1) = r1d6*(x(2,3) - x(2,4))
       ddV(1,3,3,1) = r1d6*(-x(2,2) + x(2,4))
       ddV(1,4,3,1) = r1d6*(x(2,2) - x(2,3))
       ddV(1,1,3,2) = r1d6*(-x(2,3) + x(2,4))
       ddV(1,3,3,2) = r1d6*(x(2,1) - x(2,4))
       ddV(1,4,3,2) = r1d6*(-x(2,1) + x(2,3))
       ddV(1,1,3,3) = r1d6*(x(2,2) - x(2,4))
       ddV(1,2,3,3) = r1d6*(-x(2,1) + x(2,4))
       ddV(1,4,3,3) = r1d6*(x(2,1) - x(2,2))
       ddV(1,1,3,4) = r1d6*(x(2,2) - x(2,4))
       ddV(1,2,3,4) = r1d6*(-x(2,1) + x(2,4))
       ddV(1,4,3,4) = r1d6*(x(2,1) - x(2,2))
       ddV(2,2,1,1) = r1d6*(x(3,3) - x(3,4))
       ddV(2,3,1,1) = r1d6*(-x(3,2) + x(3,4))
       ddV(2,4,1,1) = r1d6*(x(3,2) - x(3,3))
       ddV(2,1,1,2) = r1d6*(-x(3,3) + x(3,4))
       ddV(2,3,1,2) = r1d6*(x(3,1) - x(3,4))
       ddV(2,4,1,2) = r1d6*(-x(3,1) + x(3,3))
       ddV(2,1,1,3) = r1d6*(x(3,2) - x(3,4))
       ddV(2,2,1,3) = r1d6*(-x(3,1) + x(3,4))
       ddV(2,4,1,3) = r1d6*(x(3,1) - x(3,2))
       ddV(2,1,1,4) = r1d6*(-x(3,2) + x(3,3))
       ddV(2,2,1,4) = r1d6*(x(3,1) - x(3,3))
       ddV(2,3,1,4) = r1d6*(-x(3,1) + x(3,2))
       ddV(2,2,3,1) = r1d6*(-x(1,3) + x(1,4))
       ddV(2,3,3,1) = r1d6*(x(1,2) - x(1,4))
       ddV(2,4,3,1) = r1d6*(-x(1,2) + x(1,3))
       ddV(2,1,3,2) = r1d6*(x(1,3) - x(1,4))
       ddV(2,3,3,2) = r1d6*(-x(1,1) + x(1,4))
       ddV(2,4,3,2) = r1d6*(x(1,1) - x(1,3))
       ddV(2,1,3,3) = r1d6*(-x(1,2) + x(1,4))
       ddV(2,2,3,3) = r1d6*(x(1,1) - x(1,4))
       ddV(2,4,3,3) = r1d6*(-x(1,1) + x(1,2))
       ddV(2,1,3,4) = r1d6*( x(1,2) - x(1,3))
       ddV(2,2,3,4) = r1d6*(-x(1,1) + x(1,3))
       ddV(2,3,3,4) = r1d6*(x(1,1) - x(1,2))
       ddV(3,2,1,1) = r1d6*(-x(2,3) + x(2,4))
       ddV(3,3,1,1) = r1d6*(x(2,2) - x(2,4))
       ddV(3,4,1,1) = r1d6*(-x(2,2) + x(2,3))
       ddV(3,1,1,2) = r1d6*(x(2,3) - x(2,4))
       ddV(3,3,1,2) = r1d6*(-x(2,1) + x(2,4))
       ddV(3,4,1,2) = r1d6*(x(2,1) - x(2,3))
       ddV(3,1,1,3) = r1d6*(-x(2,2) + x(2,4))
       ddV(3,2,1,3) = r1d6*(x(2,1) - x(2,4))
       ddV(3,4,1,3) = r1d6*(-x(2,1) + x(2,2))
       ddV(3,1,1,4) = r1d6*(x(2,2) - x(2,3))
       ddV(3,2,1,4) = r1d6*(-x(2,1) + x(2,3))
       ddV(3,3,1,4) = r1d6*(x(2,1) - x(2,2))
       ddV(3,2,2,1) = r1d6*(x(1,3) - x(1,4))
       ddV(3,3,2,1) = r1d6*(-x(1,2) + x(1,4))
       ddV(3,4,2,1) = r1d6*(x(1,2) - x(1,3))
       ddV(3,1,2,2) = r1d6*(-x(1,3) + x(1,4))
       ddV(3,3,2,2) = r1d6*(x(1,1) - x(1,4))
       ddV(3,4,2,2) = r1d6*(-x(1,1) + x(1,3))
       ddV(3,1,2,3) = r1d6*(x(1,2) - x(1,4))
       ddV(3,2,2,3) = r1d6*(-x(1,1) + x(1,4))
       ddV(3,4,2,3) = r1d6*(x(1,1) - x(1,2))
       ddV(3,1,2,4) = r1d6*(-x(1,2) + x(1,3))
       ddV(3,2,2,4) = r1d6*(x(1,1) - x(1,3))
       ddV(3,3,2,4) = r1d6*(-x(1,1) + x(1,2))
     
      return
      end
      
      
c===========================================================================
      subroutine comp_tV2(s,V,dV,ddV,err)
c---------------------------------------------------------------------------
      implicit none
      integer iedge,nedge
      double precision s(6),s2(6), V, dV(6), ddV(6,6),   
     *                 r2, r3, r4, r5, r6, 
     *                 fact, sgn, V2
     
      integer err

      data r2,    r3,   r4,   r5,   r6
     *   / 2.d0, 3.d0, 4.d0, 5.d0, 6.0d0 /


      
      nedge=6

      do iedge=1, nedge
        s2(iedge) = s(iedge)**r2
      enddo
      
       V2=-s2(1)**r2*s2(4) - 
     -  s2(2)**r2*s2(5) - 
     -  s2(6)*(
     -  s2(3)**r2 + s2(4)*s2(5) - s2(3)*(s2(4) + s2(5) - s2(6)))+
     -  s2(2)*(
     -  s2(5)*(s2(4) - s2(5) + s2(6)) + s2(3)*(-s2(4) + s2(5) + s2(6)))+
     -  s2(1)*(
     -  s2(2)*(s2(4) + s2(5) - s2(6)) + 
     -  s2(3)*(s2(4) - s2(5) + s2(6)) + 
     -  s2(4)*(-s2(4) + s2(5) + s2(6)))
     
      err=0
      if (abs(V2).lt.1.d-14) then
       write(*,*) "V2=",V2
        do iedge=1, nedge
        write(*,*)  s(iedge)
        enddo
       err=-1
       return
       endif
     
      sgn=V2/abs(V2)
     
      V=Sqrt(abs(V2))/12.
     
      fact=sgn/(144.*V)
     
     
     
       dV(1) = fact *(s(1)*(-r2*s2(1)*s2(4) + s2(3)*s2(4) - s2(4)**r2 - 
     *  s2(3)*s2(5) + s2(4)*s2(5) + s2(3)*s2(6) + s2(4)*s2(6) + 
     -      s2(2)*(s2(4) + s2(5) - s2(6))))
     
      dV(2) = fact *(s(2)*(s2(1)*(s2(4) + s2(5) - s2(6)) + 
     * s2(5)*(-r2*s2(2) + s2(4) - s2(5) + s2(6)) + 
     -      s2(3)*(-s2(4) + s2(5) + s2(6))))
     
      dV(3) = fact * (s(3)*(s2(6)*(-r2*s2(3) + s2(4) + s2(5) - 
     * s2(6)) + s2(1)*(s2(4) - s2(5) + s2(6)) + 
     -      s2(2)*(-s2(4) + s2(5) + s2(6))))
     
      dV(4) = fact * (s(4)*(-s2(1)**r2 - (s2(3) - s2(5))*(s2(2) - 
     * s2(6)) +  s2(1)*(s2(2) + s2(3) - r2*s2(4) + s2(5) + s2(6))))
     
      dV(5) = fact * (s(5)*(-s2(2)**r2 + s2(1)*(s2(2) - s2(3) + s2(4)) +
     *   (s2(3) - s2(4))*s2(6) + s2(2)*(s2(3) + s2(4) - 
     * r2*s2(5) + s2(6))))
     
      dV(6) = fact * (s(6)*(-s2(3)**r2 + s2(3)*s2(4) + s2(1)*(-s2(2) + 
     * s2(3) + s2(4)) + s2(3)*s2(5) - s2(4)*s2(5) + s2(2)*(s2(3) + 
     * s2(5)) - r2*s2(3)*s2(6)))
     
     
     
       fact=sgn/(20736.*V**r3)


       !dV/d s(1)
         ddV(1,1) = fact *((-16.*s2(1)*(-r2*s2(1)*s2(4) + s2(3)*s2(4) - 
     *   s2(4)**r2 - s2(3)*s2(5) + s2(4)*s2(5) + s2(3)*s2(6) + 
     *  s2(4)*s2(6) + 
     -          s2(2)*(s2(4) + s2(5) - s2(6)))**r2 + 
     -  16.*(-r6*s2(1)*s2(4) + s2(3)*s2(4) - s2(4)**r2 - s2(3)*s2(5) + 
     -  s2(4)*s2(5) + s2(3)*s2(6) + 
     -         s2(4)*s2(6) + s2(2)*(s2(4) + s2(5) - s2(6)))*
     -       (-(s2(1)**r2*s2(4)) - s2(2)**r2*s2(5) - s2(6)*(s2(3)**r2 + 
     -  s2(4)*s2(5) - s2(3)*(s2(4) + s2(5) - s2(6))) + 
     -         s2(2)*(s2(5)*(s2(4) - s2(5) + s2(6)) + s2(3)*(-s2(4) + 
     - s2(5) + s2(6))) + 
     -         s2(1)*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) - 
     - s2(5) + s2(6)) + s2(4)*(-s2(4) + s2(5) + s2(6)))))/16.)
     
        ddV(2,1) = fact *(s(1)*s(2)*(s2(3)**r2 - s2(3)*s2(4) + 
     -  s2(1)*(s2(2) - s2(3) - s2(4)) - s2(3)*s2(5) + s2(4)*s2(5) - 
     -  s2(2)*(s2(3) + s2(5)) + r2*s2(3)*s2(6))*
     -    (s2(4)**r2 + (s2(5) - s2(6))**r2 - r2*s2(4)*(s2(5) + s2(6))))
     
        ddV(3,1) = fact *(-(s(1)*s(3)*(-s2(2)**r2 + s2(1)*(s2(2) - 
     -  s2(3) + s2(4)) + (s2(3) - s2(4))*s2(6) + s2(2)*(s2(3) + s2(4) -
     -  r2*s2(5) + s2(6)))*
     -      (s2(4)**r2 + (s2(5) - s2(6))**r2 - 
     -  r2*s2(4)*(s2(5) + s2(6)))))
     
        ddV(4,1) = fact *(s(1)*s(4)*(r2*s(1)**r6*s2(4) - 
     -   r2*s(2)**r6*s2(5) + s2(2)**r2*(s2(5)*(r5*s2(4) - r5*s2(5) + 
     -  s2(6)) + s2(3)*(-s2(4) + s2(5) + s2(6))) - 
     -      r3*s2(1)**r2*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) -
     -   s2(5) + s2(6)) + s2(4)*(-s2(4) + s2(5) + s2(6))) + 
     -      s2(2)*(s2(3)**r2*(-s2(4) + s2(5) + s2(6)) + 
     -  s2(5)*(-r3*s2(4)**r2 - r2*s2(5)**r2 + s2(5)*s2(6) + s2(6)**r2 + 
     -  s2(4)*(r5*s2(5) - r4*s2(6))) + 
     -         s2(3)*(r3*s2(4)**r2 + s2(5)**r2 + r6*s2(5)*s2(6) + 
     -  s2(6)**r2 - r4*s2(4)*(s2(5) + s2(6)))) + 
     -      s2(1)*(s2(2)**r2*(s2(4) + r5*s2(5) - s2(6)) + 
     -  s2(3)**r2*(s2(4) - s2(5) + r5*s2(6)) + 
     -         s2(2)*(-r3*s2(4)**r2 + r5*s2(5)**r2 - r4*s2(5)*s2(6) - 
     - s2(6)**r2 - r2*s2(4)*(s2(5) - r2*s2(6)) + 
     -  r4*s2(3)*(s2(4) - s2(5) - s2(6))) + 
     -         s2(4)*(r2*s2(4)**r2 + s2(5)**r2 + r4*s2(5)*s2(6) + 
     -  s2(6)**r2 - r3*s2(4)*(s2(5) + s2(6))) - 
     -         s2(3)*(r3*s2(4)**r2 + s2(5)**r2 + r4*s2(5)*s2(6) - 
     -  r5*s2(6)**r2 + s2(4)*(-r4*s2(5) + r2*s2(6)))) + 
     -      s2(6)*(-r2*s(3)**r6 + s2(3)**r2*(r5*s2(4) + s2(5) - 
     -  r5*s2(6)) + s2(4)*s2(5)*(r3*s2(4) - s2(5) - s2(6)) + 
     -         s2(3)*(-r3*s2(4)**r2 + s2(5)**r2 + s2(5)*s2(6) - 
     -  r2*s2(6)**r2 + s2(4)*(-r4*s2(5) + r5*s2(6))))))
     
        ddV(5,1) = fact *(-(s(1)*(s2(2)**r2 + (s2(3) - s2(4))**r2 - 
     -   r2*s2(2)*(s2(3) + s2(4)))*s(5)*
     -      (s2(6)*(-r2*s2(3) + s2(4) + s2(5) - s2(6)) + s2(1)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(2)*(-s2(4) + s2(5) + s2(6)))))
     
        ddV(6,1) = fact *(-(s(1)*(s2(2)**r2 + (s2(3) - s2(4))**r2 - 
     -  r2*s2(2)*(s2(3) + s2(4)))*s(6)*
     -      (s2(1)*(s2(4) + s2(5) - s2(6)) + s2(5)*(-r2*s2(2) + s2(4) - 
     -  s2(5) + s2(6)) + s2(3)*(-s2(4) + s2(5) + s2(6)))))
       !dV/d s(2)
       ddV(1,2) = fact *  (s(1)*s(2)*(s2(3)**r2 - s2(3)*s2(4) + 
     -  s2(1)*(s2(2) - s2(3) - s2(4)) - s2(3)*s2(5) + s2(4)*s2(5) - 
     -  s2(2)*(s2(3) + s2(5)) + 
     -      r2*s2(3)*s2(6))*(s2(4)**r2 + (s2(5) - s2(6))**r2 - 
     -  r2*s2(4)*(s2(5) + s2(6))))
     
       ddV(2,2) = fact *((-16.*s2(2)*(s2(1)*(s2(4) + s2(5) - s2(6)) + 
     -  s2(5)*(-r2*s2(2) + s2(4) - s2(5) + s2(6)) + 
     -  s2(3)*(-s2(4) + s2(5) + s2(6)))**r2 - 
     -      16.*(s2(1)*(s2(4) + s2(5) - s2(6)) + s2(5)*(-r6*s2(2) + 
     -  s2(4) - s2(5) + s2(6)) + s2(3)*(-s2(4) + s2(5) + s2(6)))*
     -       (s2(1)**r2*s2(4) + s2(2)**r2*s2(5) + s2(2)*(s2(3)*(s2(4) - 
     -  s2(5) - s2(6)) + s2(5)*(-s2(4) + s2(5) - s2(6))) + 
     -         s2(6)*(s2(3)**r2 + s2(4)*s2(5) - 
     -  s2(3)*(s2(4) + s2(5) - s2(6))) - 
     -         s2(1)*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(4)*(-s2(4) + s2(5) + s2(6)))))/16.)
     
       ddV(3,2) = fact *(s(2)*s(3)*(s2(4)**r2 + (s2(5) - s2(6))**r2 - 
     -  r2*s2(4)*(s2(5) + s2(6)))*
     -    (s2(1)**r2 + (s2(3) - s2(5))*(s2(2) - s2(6)) - s2(1)*(s2(2) + 
     -  s2(3) - r2*s2(4) + s2(5) + s2(6))))
     
       ddV(4,2) = fact *(-(s(2)*s(4)*(s2(1)**r2 + (s2(3) - s2(5))**r2 - 
     -  r2*s2(1)*(s2(3) + s2(5)))*
     -      (s2(6)*(-r2*s2(3) + s2(4) + s2(5) - s2(6)) + s2(1)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(2)*(-s2(4) + s2(5) + s2(6)))))
     
       ddV(5,2) = fact * (s(2)*s(5)*(-r2*s(1)**r6*s2(4) + 
     -  r2*s(2)**r6*s2(5) + r3*s2(2)**r2*(s2(3)*(s2(4) - s2(5) - 
     -  s2(6)) + s2(5)*(-s2(4) + s2(5) - s2(6))) + 
     -      s2(1)**r2*(s2(2)*(r5*s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) -
     -   s2(5) + s2(6)) + s2(4)*(-r5*s2(4) + r5*s2(5) + s2(6))) + 
     -      s2(6)*(-r2*s(3)**r6 + s2(3)**r2*(s2(4) + r5*s2(5) - 
     -  r5*s2(6)) - s2(4)*s2(5)*(s2(4) - r3*s2(5) + s2(6)) + 
     -         s2(3)*(s2(4)**r2 - r3*s2(5)**r2 + r5*s2(5)*s2(6) - 
     -  r2*s2(6)**r2 + s2(4)*(-r4*s2(5) + s2(6)))) + 
     -      s2(2)*(s2(3)**r2*(-s2(4) + s2(5) + r5*s2(6)) - 
     -  s2(3)*(s2(4)**r2 + r3*s2(5)**r2 + r2*s2(5)*s2(6) - 
     -  r5*s2(6)**r2 - r4*s2(4)*(s2(5) - s2(6))) + 
     -         s2(5)*(s2(4)**r2 + r2*s2(5)**r2 - r3*s2(5)*s2(6) + 
     -  s2(6)**r2 + s2(4)*(-r3*s2(5) + r4*s2(6)))) + 
     -      s2(1)*(-r3*s2(2)**r2*(s2(4) + s2(5) - s2(6)) + 
     -  s2(3)**r2*(s2(4) - s2(5) + s2(6)) + 
     -         s2(4)*(-r2*s2(4)**r2 - r3*s2(5)**r2 - r4*s2(5)*s2(6) + 
     -  s2(6)**r2 + s2(4)*(r5*s2(5) + s2(6))) - 
     -         s2(2)*(-r5*s2(4)**r2 + r3*s2(5)**r2 - r4*s2(5)*s2(6) + 
     -  s2(6)**r2 + r4*s2(3)*(s2(4) - s2(5) + s2(6)) + r2*s2(4)*(s2(5) +
     -   r2*s2(6))) + s2(3)*(s2(4)**r2 + r3*s2(5)**r2 - r4*s2(5)*s2(6) +
     -  s2(6)**r2 + s2(4)*(-r4*s2(5) + r6*s2(6))))))
     
       ddV(6,2) = fact *(s(2)*(s2(1)**r2 + (s2(3) - s2(5))**r2 - 
     -  r2*s2(1)*(s2(3) + s2(5)))*s(6)*
     -    (r2*s2(1)*s2(4) - s2(3)*s2(4) + s2(4)**r2 + s2(3)*s2(5) - 
     -  s2(4)*s2(5) - s2(3)*s2(6) - s2(4)*s2(6) - 
     -  s2(2)*(s2(4) + s2(5) - s2(6))))
       !dV/d s(3)
       ddV(1,3) = fact *  (-(s(1)*s(3)*(-s2(2)**r2 + s2(1)*(s2(2) - 
     -  s2(3) + s2(4)) + (s2(3) - s2(4))*s2(6) + s2(2)*(s2(3) + 
     -  s2(4) - r2*s2(5) + s2(6)))*
     -  (s2(4)**r2 + (s2(5) - s2(6))**r2 - r2*s2(4)*(s2(5) + s2(6)))))
     
       ddV(2,3) = fact * (s(2)*s(3)*(s2(4)**r2 + (s2(5) - s2(6))**r2 - 
     -  r2*s2(4)*(s2(5) + s2(6)))*
     -    (s2(1)**r2 + (s2(3) - s2(5))*(s2(2) - s2(6)) - 
     -  s2(1)*(s2(2) + s2(3) - r2*s2(4) + s2(5) + s2(6))))
     
       ddV(3,3) = fact * ((-16.*s2(3)*(s2(6)*(-r2*s2(3) + s2(4) + 
     -  s2(5) - s2(6)) + s2(1)*(s2(4) - s2(5) + s2(6)) + s2(2)*(-s2(4) +
     -  s2(5) + s2(6)))**r2 - 
     -      16.*(s2(6)*(-r6*s2(3) + s2(4) + s2(5) - s2(6)) + 
     -  s2(1)*(s2(4) - s2(5) + s2(6)) + s2(2)*(-s2(4) + s2(5) + s2(6)))*
     -       (s2(1)**r2*s2(4) + s2(2)**r2*s2(5) + s2(2)*(s2(3)*(s2(4) - 
     -  s2(5) - s2(6)) + s2(5)*(-s2(4) + s2(5) - s2(6))) + 
     -         s2(6)*(s2(3)**r2 + s2(4)*s2(5) - 
     -  s2(3)*(s2(4) + s2(5) - s2(6))) - 
     -         s2(1)*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(4)*(-s2(4) + s2(5) + s2(6)))))/16.)
     
       ddV(4,3) = fact *  (-(s(3)*s(4)*(s2(1)**r2 + (s2(2) - 
     -  s2(6))**r2 - r2*s2(1)*(s2(2) + s2(6)))*
     -      (s2(1)*(s2(4) + s2(5) - s2(6)) + s2(5)*(-r2*s2(2) + s2(4) - 
     -  s2(5) + s2(6)) + s2(3)*(-s2(4) + s2(5) + s2(6)))))
     
       ddV(5,3) = fact *  (s(3)*s(5)*(r2*s2(1)*s2(4) - s2(3)*s2(4) + 
     -  s2(4)**r2 + s2(3)*s2(5) - s2(4)*s2(5) - 
     -  s2(3)*s2(6) - s2(4)*s2(6) - 
     -      s2(2)*(s2(4) + s2(5) - s2(6)))*(s2(1)**r2 + (s2(2) - 
     -  s2(6))**r2 - r2*s2(1)*(s2(2) + s2(6))))
     
       ddV(6,3) = fact * (s(3)*s(6)*(-r2*s(1)**r6*s2(4) - 
     -  r2*s(2)**r6*s2(5) + s2(2)**r2*(s2(3)*(-s2(4) + r5*s2(5) + 
     -  s2(6)) + s2(5)*(s2(4) - r5*s2(5) + r5*s2(6))) + 
     -      s2(1)**r2*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(r5*s2(4) -
     -   s2(5) + s2(6)) + s2(4)*(-r5*s2(4) + s2(5) + r5*s2(6))) + 
     -    s2(6)*(r2*s(3)**r6 - s2(4)*s2(5)*(s2(4) + s2(5) - r3*s2(6)) - 
     -  r3*s2(3)**r2*(s2(4) + s2(5) - s2(6)) + 
     -         s2(3)*(s2(4)**r2 + s2(5)**r2 - r3*s2(5)*s2(6) + 
     -  r2*s2(6)**r2 + s2(4)*(r4*s2(5) - r3*s2(6)))) + 
     -      s2(2)*(r3*s2(3)**r2*(s2(4) - s2(5) - s2(6)) + 
     -  s2(5)*(s2(4)**r2 - r2*s2(5)**r2 + r5*s2(5)*s2(6) - 
     -  r3*s2(6)**r2 + s2(4)*(s2(5) - r4*s2(6))) - 
     -         s2(3)*(s2(4)**r2 - r5*s2(5)**r2 + r2*s2(5)*s2(6) + 
     -  r3*s2(6)**r2 + r4*s2(4)*(s2(5) - s2(6)))) + 
     -      s2(1)*(s2(2)**r2*(s2(4) + s2(5) - s2(6)) - 
     -  r3*s2(3)**r2*(s2(4) - s2(5) + s2(6)) + 
     -         s2(2)*(s2(4)**r2 + s2(5)**r2 - r4*s2(5)*s2(6) + 
     -  r3*s2(6)**r2 + s2(4)*(r6*s2(5) - r4*s2(6)) - r4*s2(3)*(s2(4) + 
     -  s2(5) - s2(6))) + 
     -         s2(3)*(r5*s2(4)**r2 - s2(5)**r2 + r4*s2(5)*s2(6) - 
     -  r3*s2(6)**r2 - r2*s2(4)*(r2*s2(5) + s2(6))) + 
     -         s2(4)*(-r2*s2(4)**r2 + s2(5)**r2 - r4*s2(5)*s2(6) - 
     -  r3*s2(6)**r2 + s2(4)*(s2(5) + r5*s2(6))))))
       !dV/d s(4)
       
       ddV(1,4) = fact *     (s(1)*s(4)*(r2*s(1)**r6*s2(4) - 
     -  r2*s(2)**r6*s2(5) + s2(2)**r2*(s2(5)*(r5*s2(4) - r5*s2(5) + 
     -  s2(6)) + s2(3)*(-s2(4) + s2(5) + s2(6))) - 
     -      r3*s2(1)**r2*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) -
     -   s2(5) + s2(6)) + s2(4)*(-s2(4) + s2(5) + s2(6))) + 
     -      s2(2)*(s2(3)**r2*(-s2(4) + s2(5) + s2(6)) + 
     *  s2(5)*(-r3*s2(4)**r2 - r2*s2(5)**r2 + s2(5)*s2(6) + s2(6)**r2 + 
     -  s2(4)*(r5*s2(5) - r4*s2(6))) + 
     -         s2(3)*(r3*s2(4)**r2 + s2(5)**r2 + r6*s2(5)*s2(6) + 
     -  s2(6)**r2 - r4*s2(4)*(s2(5) + s2(6)))) + 
     -      s2(1)*(s2(2)**r2*(s2(4) + r5*s2(5) - s2(6)) + 
     -  s2(3)**r2*(s2(4) - s2(5) + r5*s2(6)) + 
     -         s2(2)*(-r3*s2(4)**r2 + r5*s2(5)**r2 - r4*s2(5)*s2(6) - 
     -  s2(6)**r2 - r2*s2(4)*(s2(5) - r2*s2(6)) + r4*s2(3)*(s2(4) - 
     -  s2(5) - s2(6))) + 
     -         s2(4)*(r2*s2(4)**r2 + s2(5)**r2 + r4*s2(5)*s2(6) + 
     -  s2(6)**r2 - r3*s2(4)*(s2(5) + s2(6))) - 
     -         s2(3)*(r3*s2(4)**r2 + s2(5)**r2 + r4*s2(5)*s2(6) - 
     -  r5*s2(6)**r2 + s2(4)*(-r4*s2(5) + r2*s2(6)))) + 
     -      s2(6)*(-r2*s(3)**r6 + s2(3)**r2*(r5*s2(4) + s2(5) - 
     -  r5*s2(6)) + s2(4)*s2(5)*(r3*s2(4) - s2(5) - s2(6)) + 
     -         s2(3)*(-r3*s2(4)**r2 + s2(5)**r2 + s2(5)*s2(6) - 
     -  r2*s2(6)**r2 + s2(4)*(-r4*s2(5) + r5*s2(6))))))
     
        ddV(2,4) = fact *    (-(s(2)*s(4)*(s2(1)**r2 + (s2(3) - 
     -   s2(5))**r2 - r2*s2(1)*(s2(3) + s2(5)))*
     -      (s2(6)*(-r2*s2(3) + s2(4) + s2(5) - s2(6)) + s2(1)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(2)*(-s2(4) + s2(5) + s2(6)))))
     
        ddV(3,4) = fact *    (-(s(3)*s(4)*(s2(1)**r2 + (s2(2) - 
     -   s2(6))**r2 - r2*s2(1)*(s2(2) + s2(6)))*
     -      (s2(1)*(s2(4) + s2(5) - s2(6)) + s2(5)*(-r2*s2(2) + 
     -  s2(4) - s2(5) + s2(6)) + s2(3)*(-s2(4) + s2(5) + s2(6)))))
     
        ddV(4,4) = fact *  ((-16.*s2(4)*(s2(1)**r2 + (s2(3) - 
     -   s2(5))*(s2(2) - s2(6)) - s2(1)*(s2(2) + s2(3) - r2*s2(4) + 
     -  s2(5) + s2(6)))**r2 + 
     -      16.*(-s2(1)**r2 - (s2(3) - s2(5))*(s2(2) - s2(6)) + 
     -  s2(1)*(s2(2) + s2(3) - r6*s2(4) + s2(5) + s2(6)))*
     -       (-(s2(1)**r2*s2(4)) - s2(2)**r2*s2(5) - s2(6)*(s2(3)**r2 + 
     -  s2(4)*s2(5) - s2(3)*(s2(4) + s2(5) - s2(6))) + 
     -         s2(2)*(s2(5)*(s2(4) - s2(5) + s2(6)) + s2(3)*(-s2(4) + 
     -  s2(5) + s2(6))) + 
     -         s2(1)*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(4)*(-s2(4) + s2(5) + s2(6)))))/16.)
     
        ddV(5,4) = fact *  (s(4)*s(5)*(s2(3)**r2 - s2(3)*s2(4) + 
     -   s2(1)*(s2(2) - s2(3) - s2(4)) - s2(3)*s2(5) + s2(4)*s2(5) - 
     -  s2(2)*(s2(3) + s2(5)) + r2*s2(3)*s2(6))*
     -    (s2(1)**r2 + (s2(2) - s2(6))**r2 - r2*s2(1)*(s2(2) + s2(6))))
     
        ddV(6,4) = fact *  (-(s(4)*(s2(1)**r2 + (s2(3) - s2(5))**r2 - 
     -  r2*s2(1)*(s2(3) + s2(5)))*s(6)*
     -      (-s2(2)**r2 + s2(1)*(s2(2) - s2(3) + s2(4)) + (s2(3) - 
     -  s2(4))*s2(6) + s2(2)*(s2(3) + s2(4) - r2*s2(5) + s2(6)))))
       !dV/d s(5)
        ddV(1,5) = fact *  (-(s(1)*(s2(2)**r2 + (s2(3) - s2(4))**r2 - 
     -   r2*s2(2)*(s2(3) + s2(4)))*s(5)*
     -      (s2(6)*(-r2*s2(3) + s2(4) + s2(5) - s2(6)) + s2(1)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(2)*(-s2(4) + s2(5) + s2(6)))))
     
        ddV(2,5) = fact *(s(2)*s(5)*(-r2*s(1)**r6*s2(4) + 
     -   r2*s(2)**r6*s2(5) + r3*s2(2)**r2*(s2(3)*(s2(4) - s2(5) - 
     -  s2(6)) + s2(5)*(-s2(4) + s2(5) - s2(6))) + 
     -      s2(1)**r2*(s2(2)*(r5*s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) -
     -   s2(5) + s2(6)) + s2(4)*(-r5*s2(4) + r5*s2(5) + s2(6))) + 
     -      s2(6)*(-r2*s(3)**r6 + s2(3)**r2*(s2(4) + r5*s2(5) - 
     -  r5*s2(6)) - s2(4)*s2(5)*(s2(4) - r3*s2(5) + s2(6)) + 
     -         s2(3)*(s2(4)**r2 - r3*s2(5)**r2 + r5*s2(5)*s2(6) - 
     -  r2*s2(6)**r2 + s2(4)*(-r4*s2(5) + s2(6)))) + 
     -      s2(2)*(s2(3)**r2*(-s2(4) + s2(5) + r5*s2(6)) - 
     -  s2(3)*(s2(4)**r2 + r3*s2(5)**r2 + r2*s2(5)*s2(6) - 
     -  r5*s2(6)**r2 - r4*s2(4)*(s2(5) - s2(6))) + 
     -         s2(5)*(s2(4)**r2 + r2*s2(5)**r2 - r3*s2(5)*s2(6) + 
     -  s2(6)**r2 + s2(4)*(-r3*s2(5) + r4*s2(6)))) + 
     -      s2(1)*(-r3*s2(2)**r2*(s2(4) + s2(5) - s2(6)) + 
     -  s2(3)**r2*(s2(4) - s2(5) + s2(6)) + 
     -         s2(4)*(-r2*s2(4)**r2 - r3*s2(5)**r2 - r4*s2(5)*s2(6) + 
     -  s2(6)**r2 + s2(4)*(r5*s2(5) + s2(6))) - 
     -         s2(2)*(-r5*s2(4)**r2 + r3*s2(5)**r2 - r4*s2(5)*s2(6) + 
     -  s2(6)**r2 + r4*s2(3)*(s2(4) - s2(5) + s2(6)) + 
     -  r2*s2(4)*(s2(5) + r2*s2(6))) + 
     -         s2(3)*(s2(4)**r2 + r3*s2(5)**r2 - r4*s2(5)*s2(6) + 
     -  s2(6)**r2 + s2(4)*(-r4*s2(5) + r6*s2(6))))))
     
        ddV(3,5) = fact * (s(3)*s(5)*(r2*s2(1)*s2(4) - s2(3)*s2(4) + 
     -   s2(4)**r2 + s2(3)*s2(5) - s2(4)*s2(5) - 
     -  s2(3)*s2(6) - s2(4)*s2(6) - 
     -      s2(2)*(s2(4) + s2(5) - s2(6)))*(s2(1)**r2 + (s2(2) - 
     -  s2(6))**r2 - r2*s2(1)*(s2(2) + s2(6))))
     
        ddV(4,5) = fact * (s(4)*s(5)*(s2(3)**r2 - s2(3)*s2(4) + 
     -   s2(1)*(s2(2) - s2(3) - s2(4)) - s2(3)*s2(5) + s2(4)*s2(5) - 
     -  s2(2)*(s2(3) + s2(5)) + r2*s2(3)*s2(6))*
     -    (s2(1)**r2 + (s2(2) - s2(6))**r2 - r2*s2(1)*(s2(2) + s2(6))))
     
        ddV(5,5) = fact * ((-16.*s2(5)*(-s2(2)**r2 + s2(1)*(s2(2) - 
     -   s2(3) + s2(4)) + (s2(3) - s2(4))*s2(6) + s2(2)*(s2(3) + s2(4) -
     -   r2*s2(5) + s2(6)))**r2 - 
     -      16.*(-s2(2)**r2 + s2(1)*(s2(2) - s2(3) + s2(4)) + (s2(3) - 
     -  s2(4))*s2(6) + s2(2)*(s2(3) + s2(4) - r6*s2(5) + s2(6)))*
     -       (s2(1)**r2*s2(4) + s2(2)**r2*s2(5) + s2(2)*(s2(3)*(s2(4) - 
     -  s2(5) - s2(6)) + s2(5)*(-s2(4) + s2(5) - s2(6))) + 
     -         s2(6)*(s2(3)**r2 + s2(4)*s2(5) - 
     -  s2(3)*(s2(4) + s2(5) - s2(6))) - 
     -         s2(1)*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(4)*(-s2(4) + s2(5) + s2(6)))))/16.)
     
        ddV(6,5) = fact * ((s2(2)**r2 + (s2(3) - s2(4))**r2 - 
     -   r2*s2(2)*(s2(3) + s2(4)))*s(5)*s(6)*
     -    (s2(1)**r2 + (s2(3) - s2(5))*(s2(2) - s2(6)) - s2(1)*(s2(2) + 
     -  s2(3) - r2*s2(4) + s2(5) + s2(6))))
       !dV/d s(6)
        ddV(1,6) = fact *(-(s(1)*(s2(2)**r2 + (s2(3) - s2(4))**r2 - 
     -   r2*s2(2)*(s2(3) + s2(4)))*s(6)*
     -      (s2(1)*(s2(4) + s2(5) - s2(6)) + s2(5)*(-r2*s2(2) + s2(4) - 
     -  s2(5) + s2(6)) + s2(3)*(-s2(4) + s2(5) + s2(6)))))
     
        ddV(2,6) = fact *(s(2)*(s2(1)**r2 + (s2(3) - s2(5))**r2 - 
     -   r2*s2(1)*(s2(3) + s2(5)))*s(6)*
     -    (r2*s2(1)*s2(4) - s2(3)*s2(4) + s2(4)**r2 + s2(3)*s2(5) - 
     -  s2(4)*s2(5) - s2(3)*s2(6) - s2(4)*s2(6) - 
     -  s2(2)*(s2(4) + s2(5) - s2(6))))
     
        ddV(3,6) = fact *(s(3)*s(6)*(-r2*s(1)**r6*s2(4) - 
     -   r2*s(2)**r6*s2(5) + s2(2)**r2*(s2(3)*(-s2(4) + r5*s2(5) + 
     -  s2(6)) + s2(5)*(s2(4) - r5*s2(5) + r5*s2(6))) + 
     -      s2(1)**r2*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(r5*s2(4) -
     -   s2(5) + s2(6)) + s2(4)*(-r5*s2(4) + s2(5) + r5*s2(6))) + 
     -      s2(6)*(r2*s(3)**r6 - s2(4)*s2(5)*(s2(4) + s2(5) - 
     -  r3*s2(6)) - r3*s2(3)**r2*(s2(4) + s2(5) - s2(6)) + 
     -         s2(3)*(s2(4)**r2 + s2(5)**r2 - r3*s2(5)*s2(6) + 
     -  r2*s2(6)**r2 + s2(4)*(r4*s2(5) - r3*s2(6)))) + 
     -      s2(2)*(r3*s2(3)**r2*(s2(4) - s2(5) - s2(6)) + 
     -  s2(5)*(s2(4)**r2 - r2*s2(5)**r2 + r5*s2(5)*s2(6) - 
     -  r3*s2(6)**r2 + s2(4)*(s2(5) - r4*s2(6))) - 
     -         s2(3)*(s2(4)**r2 - r5*s2(5)**r2 + r2*s2(5)*s2(6) + 
     -  r3*s2(6)**r2 + r4*s2(4)*(s2(5) - s2(6)))) + 
     -      s2(1)*(s2(2)**r2*(s2(4) + s2(5) - s2(6)) - 
     -  r3*s2(3)**r2*(s2(4) - s2(5) + s2(6)) + 
     -         s2(2)*(s2(4)**r2 + s2(5)**r2 - r4*s2(5)*s2(6) + 
     -  r3*s2(6)**r2 + s2(4)*(r6*s2(5) - r4*s2(6)) - 
     -  r4*s2(3)*(s2(4) + s2(5) - s2(6))) + 
     -         s2(3)*(r5*s2(4)**r2 - s2(5)**r2 + r4*s2(5)*s2(6) - 
     -  r3*s2(6)**r2 - r2*s2(4)*(r2*s2(5) + s2(6))) + 
     -         s2(4)*(-r2*s2(4)**r2 + s2(5)**r2 - r4*s2(5)*s2(6) - 
     -  r3*s2(6)**r2 + s2(4)*(s2(5) + r5*s2(6))))))
     
        ddV(4,6) = fact *(-(s(4)*(s2(1)**r2 + (s2(3) - s2(5))**r2 - 
     -   r2*s2(1)*(s2(3) + s2(5)))*s(6)*
     -      (-s2(2)**r2 + s2(1)*(s2(2) - s2(3) + s2(4)) + (s2(3) - 
     -  s2(4))*s2(6) + s2(2)*(s2(3) + s2(4) - r2*s2(5) + s2(6)))))
     
        ddV(5,6) = fact * ((s2(2)**r2 + (s2(3) - s2(4))**r2 - 
     -   r2*s2(2)*(s2(3) + s2(4)))*s(5)*s(6)*
     -    (s2(1)**r2 + (s2(3) - s2(5))*(s2(2) - s2(6)) - s2(1)*(s2(2) + 
     -  s2(3) - r2*s2(4) + s2(5) + s2(6))))
     
        ddV(6,6) = fact *((-16.*s2(6)*(s2(3)**r2 - s2(3)*s2(4) + 
     -   s2(1)*(s2(2) - s2(3) - s2(4)) - s2(3)*s2(5) + s2(4)*s2(5) - 
     *  s2(2)*(s2(3) + s2(5)) + 
     -          r2*s2(3)*s2(6))**r2 + 16.*(-s2(3)**r2 + s2(3)*s2(4) + 
     -  s2(1)*(-s2(2) + s2(3) + s2(4)) + s2(3)*s2(5) - s2(4)*s2(5) + 
     -         s2(2)*(s2(3) + s2(5)) - r6*s2(3)*s2(6))*(-
     -  (s2(1)**r2*s2(4)) - s2(2)**r2*s2(5) - 
     -         s2(6)*(s2(3)**r2 + s2(4)*s2(5) - 
     -  s2(3)*(s2(4) + s2(5) - s2(6))) + 
     -         s2(2)*(s2(5)*(s2(4) - s2(5) + s2(6)) + 
     -  s2(3)*(-s2(4) + s2(5) + s2(6))) + 
     -         s2(1)*(s2(2)*(s2(4) + s2(5) - s2(6)) + s2(3)*(s2(4) - 
     -  s2(5) + s2(6)) + s2(4)*(-s2(4) + s2(5) + s2(6)))))/16.)

      return
      end
      
      
      
      subroutine ale3D4nodedTetrahedronAspectRatiotest
     *                    (ALEdata,xl,s,p,ndm,nst,ndat)
      implicit none

      integer          ndm, nst, ndat,nen,nedge,jedge,
     *                 i, j, ii, i0, j0, jj, iedge,
     *                  nface, pow, err

      double precision ALEdata(*), xl(ndm,*), s(nst,*), p(*), 
     *                x(3,4), x0(3,4), r(3), k(3,3), 
     *                dedge(6,3,4), ddedge(6,3,4,3,4), 
     *                dW(6),ddW(6,6),
     *                A, dA(6),ddA(6,6),fact,rodri, 
     *                V, dV(6),ddV(6,6),tmp(6),
     *                D, dD(6),ddD(6,6),
     *                ri,dri(6), ddri(6,6),
     *                ro,dro(6), ddro(6,6),
     *                r1, r2, r3, r6, ratio

      data r1, r2, r3, r6 / 1.d0 , 2.d0, 3.d0 , 6.d0 /

      
      err=0
      
      pow = 2

      nedge=6
      nface=4
      nen=4
      
      ratio=ALEdata(1)
      
      
      do i=1, nen
        do j=1, ndm
          x(j,i) = xl(j,i) + xl(j,i+nen)
          x0(j,i)=           xl(j,i+nen+nen)
        enddo
      enddo

 
         call comp_ADV(x0,ndm,nen,nedge, nface,dedge, ddedge,
     *   A,dA,ddA,D,dD,ddD,V,dV,ddV,err)
     
       
     
      ro=D/(r6*V)
      ri=r3*V/A
      rodri=ro/ri
      
      
      if(rodri.gt.ratio) then
      write(*,*) rodri, " is too big, use better mesh!"
c     call prgerror(1,'ale3D4nodedTetrahedronAspectRatiotest',
c    *   'ro/ri too big, use better mesh!')
      endif



      return

      end
