      subroutine get_b(shape,u,grad_u,r,s,b)
c     ******************************************************************
c     *                                                                *
c     * B[i,j] = r[i]r[j] - (r[k]s[k])^2 * s[i]s[j]                    *
c     *                                                                *
c     * b(1) = B[1,1]                                                  *
c     * b(2) = B[2,2]                                                  *
c     * b(3) = B[1,2] = B[2,1]                                         *
c     *                                                                *
c     ******************************************************************

      implicit none

      integer          i
      double precision c, norm_r, norm_u
      double precision b(3), r(*), s(*), u(*)
      double precision grad_u(2,*), shape(3,*)

      norm_u = sqrt(u(1) * u(1) + u(2) * u(2))
      if (norm_u.lt.1.0d-16) then
        s(1) = 0.0d0
        s(2) = 0.0d0
      else
        s(1) = u(1) / norm_u
        s(2) = u(2) / norm_u
      endif

      r(1) = grad_u(1,1) * s(1) + grad_u(2,1) * s(2)
      r(2) = grad_u(1,2) * s(1) + grad_u(2,2) * s(2)

      norm_r = sqrt(r(1) * r(1) + r(2) * r(2))
      if (norm_r.lt.1.0d-16) then
        r(1) = 0.0d0
        r(2) = 0.0d0
      else
        r(1) = r(1) / norm_r
        r(2) = r(2) / norm_r
      endif

      c = (r(1) * s(1) + r(2) * s(2)) * (r(1) * s(1) + r(2) * s(2))

      b(1) = r(1) * r(1) - c * s(1) * s(1)
      b(2) = r(2) * r(2) - c * s(2) * s(2)
      b(3) = r(1) * r(2) - c * s(1) * s(2)

      return

      end



      subroutine get_dbdu(r,s,dbdu)
c     ******************************************************************
c     *                                                                *
c     *  dbdu(1,1) = dB[1,1]/dU[1]          dbdu(1,2) = dB[1,1]/dU[2]  *
c     *  dbdu(2,1) = dB[2,2]/dU[1]          dbdu(2,2) = dB[2,2]/dU[2]  *
c     *  dbdu(3,1) = dB[1,2]/dU[1]          dbdu(3,2) = dB[1,2]/dU[2]  *
c     *                                                                *
c     ******************************************************************

      implicit none

      double precision c, d, fact1, fact2, r2
      double precision r(2), s(2), u(3)
      double precision dbdu(3,2), drdu(2,2), dsdu(2,2), shape(3,3)

      data r2 / 2.0d0 /

      call get_drdu(shape,u,s,drdu)
      call get_dsdu(shape,u,s,dsdu)

      c = (r(1)*s(1)+r(2)*s(2))
      d = c * c

      fact1 = r2 * c * (drdu(1,1) * s(1) + drdu(2,1) * s(2) +
     &                  r(1) * dsdu(1,1) + r(2) * dsdu(2,1))
      fact2 = r2 * c * (drdu(1,2) * s(1) + drdu(2,2) * s(2) +
     &                  r(1) * dsdu(1,2) + r(2) * dsdu(2,2))

      dbdu(1,1) = (drdu(1,1)*r(1)+r(1)*drdu(1,1))-fact1*s(1)*s(1)-
     &            d*(dsdu(1,1)*s(1)+s(1)*dsdu(1,1))
      dbdu(1,2) = (drdu(1,2)*r(1)+r(1)*drdu(1,2))-fact2*s(1)*s(1)-
     &            d*(dsdu(1,2)*s(1)+s(1)*dsdu(1,2))
      dbdu(2,1) = (drdu(2,1)*r(2)+r(2)*drdu(2,1))-fact1*s(2)*s(2)-
     &            d*(dsdu(2,1)*s(2)+s(2)*dsdu(2,1))
      dbdu(2,2) = (drdu(2,2)*r(2)+r(2)*drdu(2,2))-fact2*s(2)*s(2)-
     &            d*(dsdu(2,2)*s(2)+s(2)*dsdu(2,2))
      dbdu(3,1) = (drdu(1,1)*r(2)+r(1)*drdu(2,1))-fact1*s(1)*s(2)-
     &            d*(dsdu(1,1)*s(2)+s(1)*dsdu(2,1))
      dbdu(3,2) = (drdu(1,2)*r(2)+r(1)*drdu(2,2))-fact2*s(1)*s(2)-
     &            d*(dsdu(1,2)*s(2)+s(1)*dsdu(2,2))

*     write(*,*)dbdu(1,1),dbdu(1,2)
*     write(*,*)dbdu(2,1),dbdu(2,2)
*     write(*,*)dbdu(3,1),dbdu(3,2)

      return

      end



      subroutine get_a(b,grad_u,a)
c     ******************************************************************
c     *                                                                *
c     *  a(1) => A[1,1] = B[1,k]*grad_u[k,1]                           *
c     *  a(2) => A[2,1] = B[2,k]*grad_u[k,1]                           *
c     *  a(3) => A[1,2] = B[1,k]*grad_u[k,2]                           *
c     *  a(4) => A[2,2] = B[2,k]*grad_u[k,2]                           *
c     *                                                                *
c     ******************************************************************

      double precision a(*), b(*)
      double precision grad_u(2,*)

      a(1) = b(1)*grad_u(1,1) + b(3)*grad_u(2,1)
      a(2) = b(3)*grad_u(1,1) + b(2)*grad_u(2,1)
      a(3) = b(1)*grad_u(1,2) + b(3)*grad_u(2,2)
      a(4) = b(3)*grad_u(1,2) + b(2)*grad_u(2,2)

      return

      end



      subroutine get_dadu(shape,grad_u,b,dbdu,dadu)
c     *******************************************************************
c     *                                                                 *
c     *                                                                 *
c     *******************************************************************

      integer          i
      double precision b(3)
      double precision dadu(4,2), dbdu(3,2), grad_u(2,2), shape(3,3)

      dadu(1,1) = dbdu(1,1)*grad_u(1,1)+dbdu(3,1)*grad_u(2,1)
      dadu(1,2) = dbdu(1,2)*grad_u(1,1)+dbdu(3,2)*grad_u(2,1)
      dadu(2,1) = dbdu(3,1)*grad_u(1,1)+dbdu(2,1)*grad_u(2,1)
      dadu(2,2) = dbdu(3,2)*grad_u(1,1)+dbdu(2,2)*grad_u(2,1)
      dadu(3,1) = dbdu(1,1)*grad_u(1,2)+dbdu(3,1)*grad_u(2,2)
      dadu(3,2) = dbdu(1,2)*grad_u(1,2)+dbdu(3,2)*grad_u(2,2)
      dadu(4,1) = dbdu(3,1)*grad_u(1,2)+dbdu(2,1)*grad_u(2,2)
      dadu(4,2) = dbdu(3,2)*grad_u(1,2)+dbdu(2,2)*grad_u(2,2)

      do i = 1, 3
        dadu(1,1) = dadu(1,1)+b(1)*shape(1,i)
        dadu(1,2) = dadu(1,2)+b(3)*shape(1,i)
        dadu(2,1) = dadu(2,1)+b(3)*shape(1,i)
        dadu(2,2) = dadu(2,2)+b(2)*shape(1,i)
        dadu(3,1) = dadu(3,1)+b(1)*shape(2,i)
        dadu(3,2) = dadu(3,2)+b(3)*shape(2,i)
        dadu(4,1) = dadu(4,1)+b(3)*shape(2,i)
        dadu(4,2) = dadu(4,2)+b(2)*shape(2,i)
      enddo

      return

      end




      subroutine get_drdu(shape,u,s,drdu)
c     ******************************************************************
c     *                                                                *
c     *           A[i]                  d || u ||                      *
c     *  r[i] = --------- where A[i] = -----------                     *
c     *          || A ||                 d x[i]                        *
c     *                                                                *
c     *                                                                *
c     *  drdu(i,j) = dr[i]/du[j]                                       *
c     ******************************************************************

      implicit none

      integer          i
      double precision fact1, fact2, invna, invnu, norm_a, norm_u
      double precision a(2), dnadu(2), s(*), u(*)
      double precision dadu(2,2), drdu(2,2), grad_u(2,2), shape(3,*)
c
c     a(i) = d||u||/dx_i
c
      a(1) = grad_u(1,1) * s(1) + grad_u(2,1) * s(2)
      a(2) = grad_u(1,2) * s(1) + grad_u(2,2) * s(2)

      norm_a = sqrt(a(1)*a(1)+a(2)*a(2))
c
c     dadu(i,j) = da[i]/du[j]
c
      call pzero(dadu,4)

      norm_u = sqrt(u(1) * u(1) + u(2) * u(2))

      if (norm_u.lt.1.0d-16) then
        invnu = 0.0d0
      else
        invnu = 1.0d0 / norm_u
      endif

      fact1 = grad_u(1,1) * s(1) + grad_u(2,1) * s(2)
      fact2 = grad_u(1,2) * s(1) + grad_u(2,2) * s(2)
      do i = 1, 3
        dadu(1,1) = dadu(1,1) + shape(1,i) * s(1) +
     &              invnu * shape(3,i) * (grad_u(1,1) - fact1 * s(1))
        dadu(1,2) = dadu(1,2) + shape(1,i) * s(2) +
     &              invnu * shape(3,i) * (grad_u(2,1) - fact1 * s(2))
        dadu(2,1) = dadu(2,1) + shape(2,i) * s(1) +
     &              invnu * shape(3,i) * (grad_u(1,2) - fact2 * s(1))
        dadu(2,2) = dadu(2,2) + shape(2,i) * s(2) +
     &              invnu * shape(3,i) * (grad_u(2,2) - fact2 * s(2))
      enddo
c
c     dnadu(i) = d||a||/du[i]
c
      if (norm_a.lt.1.0d-16) then
        invna = 0.0d0
      else
        invna = 1.0d0 / norm_a
      endif

      dnadu(1) = invna * (dadu(1,1) * a(1) + dadu(2,1) * a(2))
      dnadu(2) = invna * (dadu(1,2) * a(1) + dadu(2,2) * a(2))
c
c     dr[i]/du[j]
c
      drdu(1,1) = invna * (dadu(1,1) - invna * dnadu(1))
      drdu(1,2) = invna * (dadu(1,2) - invna * dnadu(2))
      drdu(2,1) = invna * (dadu(2,1) - invna * dnadu(1))
      drdu(2,2) = invna * (dadu(2,2) - invna * dnadu(2))

      return

      end



      subroutine get_dsdu(shape,u,s,dsdu)
c     ******************************************************************
c     *                                                                *
c     *  dsdu(1,1) = dS[1]/dU[1]                                       *
c     *  dsdu(1,2) = dS[1]/dU[2]                                       *
c     *  dsdu(2,1) = dS[2]/dU[1]                                       *
c     *  dsdu(2,2) = dS[2]/dU[2]                                       *
c     *                                                                *
c     ******************************************************************

      implicit none

      integer          i
      double precision fact, invnu, norm_u
      double precision s(2), u(3)
      double precision dsdu(2,2), shape(3,3)

      norm_u  = sqrt(u(1) * u(1) + u(2) * u(2))
      if (norm_u.lt.1.0d-16) then
        invnu = 0.0d0
      else
        invnu = 1.0d0 / norm_u
      endif

      call pzero(dsdu,4)

      do i = 1, 3
        fact = invnu * shape(3,i)
        dsdu(1,1) = dsdu(1,1) + fact * (1.0d0 - s(1) * s(1))
        dsdu(1,2) = dsdu(1,2) - fact * s(1) * s(2)
        dsdu(2,1) = dsdu(2,1) - fact * s(2) * s(1)
        dsdu(2,2) = dsdu(2,2) + fact * (1.0d0 - s(2) * s(2))
      enddo

      return

      end
