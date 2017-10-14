
c
c     integer function length(strg)
c
c     logical function fileOK(name,l)
c
c     integer function next_key(list,nlist,iof)
c
c     subroutine strip(strg) 
c     logical stripOK(strg)
c
c     subroutine replace(strg,dummy,txt,l,n)
c
c     subroutine presskey(key)
c
c     subroutine write_s(A)
c     subroutine write_u(A)
c     subroutine write_4(X)
c


c==============================================================================

      integer function length(strg)
c------------------------------------------------------
c     return length of strg
c
      implicit none
      integer j
      character strg*(*)
      j = len(strg)
c      do while ((j.gt.0).and.(.not.(
c     *  ((strg(j:j).ge.'a').and.(strg(j:j).le.'z')).or.
c     *  ((strg(j:j).ge.'A').and.(strg(j:j).le.'Z')).or.
c     *  ((strg(j:j).ge.'1').and.(strg(j:j).le.'9')).or.
c     *   (strg(j:j).eq.'0').or.
c     *   (strg(j:j).eq.'.').or.
c     *   (strg(j:j).eq.'_'))))
      do while (j.gt.0.and.(strg(j:j).eq.' '.or.strg(j:j).eq.''))
        j = j - 1
      enddo
      length = j
      return
      end
c==============================================================================

      logical function fileOK(name,mxl)
c---------------------------------------------------------
c     check whether name is an admissible file name
c
      implicit none
      integer mxl, length, i, ll
      character name*(*)
      
      ll = length(name)
      if (ll.gt.mxl.or.ll.lt.1) then
        fileOK = .false.
      else
        fileOK = .true.
        do i=1, ll
          if (.not.(
     *      ((name(i:i).ge.'a').and.(name(i:i).le.'z')).or.
     *      ((name(i:i).ge.'A').and.(name(i:i).le.'Z')).or.
     *      ((name(i:i).ge.'1').and.(name(i:i).le.'9')).or.
     *       (name(i:i).eq.'0').or.
     *       (name(i:i).eq.'.').or.
     *       (name(i:i).eq.'_'))) 
     *      fileOK = .false.
        enddo
      endif

      return
      end
c==============================================================================

      integer function next_key(list,nlist,iof)
c----------------------------------------------------
c     read next key from file and return id number
c
      implicit none
      integer   nlist, iof, ios, length, i, l
      character strg*4, list(nlist)*4
      logical   pcomp

      i = nlist + 1
      read(iof,'(A)',iostat=ios) strg
      do while ((ios.eq.0).and.(i.gt.nlist))
        i = 1
        l = length(list(i))
        do while (.not.(pcomp(strg,list(i),l)).and.(i.le.nlist))
          i = i + 1
          l = length(list(i))
        enddo
        if (i.gt.nlist) read(iof,'(A)',iostat=ios) strg
      enddo

      if (ios.eq.0) then
        next_key = i
      else
        next_key = 0
      endif

      return
      end
c==============================================================================

      subroutine strip(strg)
c----------------------------------------------------
c     remove comments and empty parts at beginning and end from strg
c
      implicit none
      integer   l, ll, length
      character strg*(*)

c...  remove comments
      ll = length(strg)
      l  = 1
      do while ((strg(l:l).ne.'!').and.(strg(l:l).ne.'#').and.(l.le.ll))
        l = l + 1
      enddo  
      strg = strg(1:l-1)
      ll = l - 1

c...  remove empty space at the beginning
      l = 1
      do while (strg(l:l).eq.' '.and.l.le.ll)
        l = l + 1
      enddo  
      if (l.gt.ll) then
        strg = ''
      else
        strg = strg(l:ll)
      endif

c...  remove empty space at the end
      strg = strg(1:length(strg))
      
      return

      end
c==============================================================================

      logical function stripOK(strg)
c----------------------------------------------------
c     remove comments from strg
c     return .false. if clean strg is empty 
c
      implicit none
      integer   length
      character strg*(*)

      stripOK = .true.

      call strip(strg)

      if (length(strg).eq.0) stripOK = .false.

      return

      end
c==============================================================================

      logical function mypcomp(a,b)
c----------------------------------------------------
c     compare two character strings a and b
c
      implicit  none
      integer   la, lb, length
      logical   pcomp
      character a*(*), b*(*)

      mypcomp = .false.

      la = length(a)
      lb = length(b)

c      write(*,*) la, lb
c      write(*,*) a(1:la)
c      write(*,*) b(1:lb)

      if (la.ne.lb) return

      mypcomp = pcomp(a,b,la)

      end
c==============================================================================

      logical function pcomp(a,b,n)
c----------------------------------------------------
c     compare two character strings a and b
c
      implicit  none
      integer   n, inc, i, ia, ib, iau, izu

      character a*(*),b*(*)

      pcomp = .false.

      iau = ichar('A')
      izu = ichar('Z')

      inc = iau - ichar('a')

      do i=1, n

        ia = ichar(a(i:i))
        ib = ichar(b(i:i))

        if (ia.ge.iau.and.ia.le.izu) ia = ia - inc
        if (ib.ge.iau.and.ib.le.izu) ib = ib - inc
        
        if (ia.ne.ib) return
      end do

      pcomp = .true.

      end
c==============================================================================

      subroutine replace(strg,dummy,txt,l,n)
c------------------------------------------------------
c     replace dummy(1..n) in strg by txt(1..n)
c
      implicit none
      integer   l, n, length,
     *          ls, lt, ld, i, il, j, k
      character strg*(*), dummy*(*), txt*(*)

      ls = length(strg)

      do i=1, n

        il = (i - 1) * l

        ld = length(dummy(il+1:il+l))
        lt = length(txt(il+1:il+l))

        j = 1

        do while (j.le.ls-ld)

          do while ((strg(j:j-1+ld).ne.dummy(il+1:il+l))
     *      .and.(j.le.ls-ld+1))
            j = j + 1
          enddo

          if (j.le.ls-ld+1) then

            ls = ls - ld + lt

            if (lt.gt.ld) then
              do k=ls, j+lt, -1
                strg(k:k) = strg(k-lt+ld:k-lt+ld)
              enddo
            elseif (ld.gt.lt) then
              do k=j+lt, ls
                strg(k:k) = strg(k+ld-lt:k+ld-lt)
              enddo
            endif

            do k=j, j+lt-1
              strg(k:k) = txt(il+k-j+1:il+k-j+1)
            enddo

            if (lt-ld.lt.0) then
              do k=ls+1, ls+ld-lt
                strg(k:k) = ''
              enddo
            endif

          endif

        enddo

      enddo

      return
      end
c==============================================================================

      subroutine presskey(key)
c------------------------------------------------------
c     wait until 'key' has been pressed
c 
c     key = 1  -> ENTER
c           2  -> Esc
c           3  -> x - ENTER
c           4  -> q - ENTER
c
      implicit none

      integer   key

      character ch

      save

      read(*,'(a)') ch





      return
      end
c==============================================================================


      subroutine write_s(A)
      implicit none
      double precision A(6)
      write(*,'(/3(1x,f6.3))') A(1), A(4), A(6)
      write(*,'( 3(1x,f6.3))') A(4), A(2), A(5)
      write(*,'( 3(1x,f6.3))') A(6), A(5), A(3)
      return
      end

      subroutine write_u(A)
      implicit none
      double precision A(3,3)
      write(*,'(/3(1x,f8.3))') A(1,1), A(1,2), A(1,3)
      write(*,'( 3(1x,f8.3))') A(2,1), A(2,2), A(2,3)
      write(*,'( 3(1x,f8.3))') A(3,1), A(3,2), A(3,3)
      return
      end

      subroutine write_4(X)
      implicit none
      integer          i, j
      double precision X(6,6)
      write(*,*)
      do i=1, 6
        write(*,'(6(1x,f8.3))') (X(i,j),j=1,6)
      enddo
      return
      end










