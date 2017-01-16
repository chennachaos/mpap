C *******************************************************************
C COPYRIGHT (c) 1995 Timothy A. Davis and Council for the Central
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
C Licence, see http://hsl.rl.ac.uk/acuk/cou.html
C
C Please note that for a UK ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between AEA
C    Technology plc and the Licensee on suitable terms and conditions,
C    which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither CCLRC nor AEA Technology plc shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C                    Laboratory of the Research Councils
C######DATE 30 November 1995
C## Initialization of XOUT and IOUT added to subroutine MA38DD.  20/1/98
C## In cases where MA38ED is not called this avoids a possibly wrong
C## decision on storage requirements after label 9000 in MA38DD.
C## Initialization of XRMAX moved to top of subroutine MA38DD.   14/1/99
C## Small modifications to commenst on iout/xout in MA38DD and MA38ED.
C#######################################################################
C## MA38: double precision version
C## user-callable routines are listed first
C#######################################################################

        SUBROUTINE MA38ID (KEEP, CNTL, ICNTL)
        INTEGER            KEEP(20)
        DOUBLE PRECISION   CNTL(10)
        INTEGER            ICNTL(20)
C=== MA38ID ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C## End of user documentation for MA38ID ###############################
C=======================================================================
C=======================================================================
        DOUBLE PRECISION FD05AD
        EXTERNAL         FD05AD
C=======================================================================
C=======================================================================
        INTEGER I
        DOUBLE PRECISION ZERO, TENTH, TWO
        PARAMETER (TENTH = 0.1D0, TWO = 2.0D0, ZERO = 0.0D0)
C=======================================================================
C=======================================================================
        ICNTL (1) = 6
        ICNTL (2) = 6
        ICNTL (3) = 2
        ICNTL (4) = 1
        ICNTL (5) = 4
        ICNTL (6) = 0
        ICNTL (7) = 16
        ICNTL (8) = 0
        DO 10 I = 9, 20
           ICNTL (I) = 0
10      CONTINUE
        CNTL (1) = TENTH
        CNTL (2) = TWO
        CNTL (3) = FD05AD (1)
        DO 30 I = 4, 10
           CNTL (I) = ZERO
30      CONTINUE
        KEEP (6) = 0
        KEEP (7) = 64
        KEEP (8) = 1
        DO 20 I = 9, 20
           KEEP (I) = 0
20      CONTINUE
        RETURN
        END
        SUBROUTINE MA38AD (N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     *                     INDEX, KEEP, CNTL, ICNTL, INFO, RINFO)
        INTEGER          N, NE, JOB
        LOGICAL          TRANSA
        INTEGER          LVALUE, LINDEX
        DOUBLE PRECISION VALUE(LVALUE)
        INTEGER          INDEX(LINDEX), KEEP(20)
        DOUBLE PRECISION CNTL(10)
        INTEGER          ICNTL(20), INFO(40)
        DOUBLE PRECISION RINFO(20)
C=== MA38AD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C## End of user documentation for MA38AD ###############################
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC MAX, MIN
        EXTERNAL  MA38DD, MA38KD, MA38ND, MA38YD
C=======================================================================
C=======================================================================
        INTEGER I, NZ, LUX1, LUI1, IUSE, XUSE, LUIR1, NZOFF, NBLKS
        LOGICAL PRESRV
        DOUBLE PRECISION
     $          ZERO
        PARAMETER (ZERO = 0.0D0)
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        DO 10 I = 1, 40
           INFO (I) = 0
10      CONTINUE
        DO 20 I = 1, 20
           RINFO (I) = ZERO
20      CONTINUE
        KEEP (1) = 0
        KEEP (2) = 0
        KEEP (3) = 0
        KEEP (4) = 0
        KEEP (5) = 0
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        CALL MA38YD (1, 1,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          RINFO(1), RINFO(1), 1, RINFO(1), 1)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IUSE = 0
        XUSE = 0
        INFO (5) = NE
        INFO (6) = NE
        IF (N .LT. 1) THEN
           CALL MA38ND (1, ICNTL, INFO, -1, -1)
           GO TO 9000
        ENDIF
        IF (NE .LT. 1) THEN
           CALL MA38ND (1, ICNTL, INFO, -2, -1)
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        NZ = NE
        IUSE = 2*N+1 + MAX (2*NZ, N+1) + NZ
        XUSE = 2*NZ
        INFO (18) = IUSE
        INFO (20) = XUSE
        INFO (19) = IUSE
        INFO (21) = XUSE
        IF (LINDEX .LT. IUSE) THEN
           CALL MA38ND (1, ICNTL, INFO, -3, IUSE)
        ENDIF
        IF (LVALUE .LT. XUSE) THEN
           CALL MA38ND (1, ICNTL, INFO, -4, XUSE)
        ENDIF
        IF (INFO (1) .LT. 0) THEN
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        CALL MA38KD (N, NZ, TRANSA, VALUE, LVALUE, INFO, ICNTL,
     $     INDEX, LINDEX-(2*N+1), INDEX(LINDEX-2*N), INDEX(LINDEX-N), 1)
        IF (INFO (1) .LT. 0) THEN
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IUSE = NZ + (N+1)
        XUSE = NZ
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        PRESRV = JOB .EQ. 1
        IF (PRESRV) THEN
           IUSE = IUSE + (N+1)
CFPP$ NODEPCHK L
           DO 30 I = 1, N+1
              INDEX (NZ+N+1+I) = INDEX (I)
30         CONTINUE
           CALL MA38DD (N, NZ, INDEX (NZ+N+2),
     $          VALUE (NZ+1), LVALUE-NZ,
     $          INDEX (NZ+2*N+3), LINDEX-(NZ+2*N+2),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, CNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX (N+2), VALUE, N, NZ, KEEP, NE)
           IF (INFO (1) .LT. 0) THEN
              GO TO 9000
           ENDIF
           LUX1 = LUX1 + NZ
           LUI1 = LUI1 + (NZ+2*N+2)
           LUX1 = LUX1 - NZ
           LUI1 = LUI1 - (NZ+N+1)
           DO 40 I = NZ+N+1, 1, -1
              INDEX (LUI1+I-1) = INDEX (I)
40         CONTINUE
           DO 50 I = NZ, 1, -1
              VALUE (LUX1+I-1) = VALUE (I)
50         CONTINUE
        ELSE
           CALL MA38DD (N, NZ, INDEX,
     $          VALUE, LVALUE,
     $          INDEX (N+2), LINDEX-(N+1),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, CNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX, VALUE, 0, 1, KEEP, NE)
           IF (INFO (1) .LT. 0) THEN
              GO TO 9000
           ENDIF
           LUI1 = LUI1 + (N+1)
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (TRANSA) THEN
           INDEX (LINDEX-6) = 1
        ELSE
           INDEX (LINDEX-6) = 0
        ENDIF
        INDEX (LINDEX-5) = NZOFF
        INDEX (LINDEX-4) = NBLKS
        IF (PRESRV) THEN
           INDEX (LINDEX-3) = 1
        ELSE
           INDEX (LINDEX-3) = 0
        ENDIF
        INDEX (LINDEX-2) = NZ
        INDEX (LINDEX-1) = N
        INDEX (LINDEX) = NE
        LUIR1 = LUI1
        IF (PRESRV) THEN
           LUIR1 = LUIR1 + N+1 + NZ
        ENDIF
        IF (NBLKS .GT. 1) THEN
           LUIR1 = LUIR1 + NZOFF
        ENDIF
        KEEP (1) = LUX1
        KEEP (2) = LVALUE
        KEEP (3) = LUI1
        KEEP (4) = LUIR1
        KEEP (5) = LINDEX
        IUSE = LINDEX - LUI1 + 1
        XUSE = LVALUE - LUX1 + 1
        INFO (22) = INFO (22) + (LINDEX - LUIR1 + 1)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
9000    CONTINUE
        IF (INFO (1) .LT. 0) THEN
           KEEP (1) = 0
           KEEP (2) = 0
           KEEP (3) = 0
           KEEP (4) = 0
           KEEP (5) = 0
        ENDIF
        INFO (18) = MIN (LINDEX, MAX (INFO (18), IUSE))
        INFO (20) = MIN (LVALUE, MAX (INFO (20), XUSE))
        CALL MA38YD (1, 2,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          RINFO(1), RINFO(1), 1, RINFO(1), 1)
        RETURN
        END
        SUBROUTINE MA38BD (N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     *                     INDEX, KEEP, CNTL, ICNTL, INFO, RINFO)
        INTEGER          N, NE, JOB
        LOGICAL          TRANSA
        INTEGER          LVALUE, LINDEX
        DOUBLE PRECISION VALUE(LVALUE)
        INTEGER          INDEX(LINDEX), KEEP(20)
        DOUBLE PRECISION CNTL(10)
        INTEGER          ICNTL(20), INFO(40)
        DOUBLE PRECISION RINFO(20)
C=== MA38BD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C## End of user documentation for MA38BD ###############################
C=======================================================================
C=======================================================================
        INTRINSIC MAX, MIN
        EXTERNAL MA38KD, MA38ND, MA38PD, MA38YD
C=======================================================================
C=======================================================================
        INTEGER I, NZ, LUX1, LUI1, IUSE, XUSE, N1, NZ1, NBLKS,
     $          LIND2, LUIR1, LUSIZ, LUI2, RPERMP, CPERMP,
     $          OFFPP, LUBLPP, BLKPP, ON, NZOFF, IP2, IO, PRL
        LOGICAL PRESRV, BADLU
        DOUBLE PRECISION
     $          ZERO
        PARAMETER (ZERO = 0.0D0)
C=======================================================================
C=======================================================================
        IO = ICNTL (2)
        PRL = ICNTL (3)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        DO 10 I = 1, 40
           INFO (I) = 0
10      CONTINUE
        DO 20 I = 1, 20
           RINFO (I) = ZERO
20      CONTINUE
        KEEP (1) = 0
        KEEP (2) = 0
        KEEP (3) = 0
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        CALL MA38YD (2, 1,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          RINFO(1), RINFO(1), 1, RINFO(1), 1)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IUSE = 0
        XUSE = 0
        INFO (5) = NE
        INFO (6) = NE
        IF (N .LT. 1) THEN
           CALL MA38ND (2, ICNTL, INFO, -1, -1)
           GO TO 9000
        ENDIF
        IF (NE .LT. 1) THEN
           CALL MA38ND (2, ICNTL, INFO, -2, -1)
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        LUIR1 = KEEP (4)
        LUI2 = KEEP (5)
        LUSIZ = LUI2 - LUIR1 + 1
        BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1 .OR. LUI2.GT.LINDEX
        IF (BADLU) THEN
           CALL MA38ND (2, ICNTL, INFO, -7, 0)
           GO TO 9000
        ENDIF
        IF (2*NE .GT. LUIR1) THEN
           CALL MA38ND (2, ICNTL, INFO, -2, LUIR1/2)
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (LUI2 .LT. LINDEX) THEN
           DO 30 I = LINDEX, LINDEX - LUSIZ + 1, -1
              INDEX (I) = INDEX (I - LINDEX + LUI2)
30         CONTINUE
           LUIR1 = LINDEX - LUSIZ + 1
           KEEP (5) = LINDEX
           KEEP (4) = LUIR1
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        N1 = INDEX (LINDEX-1)
        NZ1 = INDEX (LINDEX-2)
        NBLKS = INDEX (LINDEX-4)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        RPERMP = (LINDEX-6) - N
        CPERMP = RPERMP - N
        IP2 = CPERMP - 1
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (NBLKS .GT. 1) THEN
           OFFPP = CPERMP - (N+1)
           BLKPP = OFFPP - (NBLKS+1)
           LUBLPP = BLKPP - (NBLKS)
           IP2 = LUBLPP - 1
           ON = N
        ELSE
           OFFPP = 1
           BLKPP = 1
           LUBLPP = 1
           ON = 0
        ENDIF
        BADLU = N .NE. N1 .OR. NZ1 .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $          NBLKS .LE. 0 .OR. NBLKS .GT. N
        IF (BADLU) THEN
           CALL MA38ND (2, ICNTL, INFO, -7, 0)
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        NZ = NE
        IUSE = 2*N+1 + MAX (2*NZ,N+1) + NZ + LUSIZ
        XUSE = 2*NZ
        INFO (18) = IUSE
        INFO (20) = XUSE
        INFO (19) = IUSE
        INFO (21) = XUSE
        INFO (23) = XUSE
        LIND2 = LUIR1 - 1
        IF (LINDEX .LT. IUSE) THEN
           CALL MA38ND (2, ICNTL, INFO, -3, IUSE)
        ENDIF
        IF (LVALUE .LT. XUSE) THEN
           CALL MA38ND (2, ICNTL, INFO, -4, XUSE)
        ENDIF
        IF (INFO (1) .LT. 0) THEN
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        CALL MA38KD (N, NZ, TRANSA, VALUE, LVALUE, INFO, ICNTL,
     $     INDEX, LIND2-(2*N+1), INDEX (LIND2-2*N), INDEX (LIND2-N), 2)
        IF (INFO (1) .LT. 0) THEN
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IUSE = NZ + (N+1) + LUSIZ
        XUSE = NZ
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        PRESRV = JOB .EQ. 1
        IF (PRESRV) THEN
           IUSE = IUSE + (N+1)
CFPP$ NODEPCHK L
           DO 40 I = 1, N+1
              INDEX (NZ+N+1+I) = INDEX (I)
40         CONTINUE
           CALL MA38PD (N, NZ, INDEX (NZ+N+2),
     $          VALUE (NZ+1), LVALUE-NZ,
     $          INDEX (NZ+2*N+3), LIND2-(NZ+2*N+2),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX (N+2), VALUE, N, NZ,
     $          INDEX (LUIR1), IP2 - LUIR1 + 1,
     $          INDEX (LUBLPP), INDEX (BLKPP), INDEX (OFFPP), ON,
     $          INDEX (CPERMP), INDEX (RPERMP), NE)
           IF (INFO (1) .LT. 0) THEN
              GO TO 9000
           ENDIF
           LUX1 = LUX1 + NZ
           LUI1 = LUI1 + (NZ+2*N+2)
           LUX1 = LUX1 - (NZ)
           LUI1 = LUI1 - (NZ+N+1)
           DO 50 I = NZ+N+1, 1, -1
              INDEX (LUI1+I-1) = INDEX (I)
50         CONTINUE
           DO 60 I = NZ, 1, -1
              VALUE (LUX1+I-1) = VALUE (I)
60         CONTINUE
        ELSE
           CALL MA38PD (N, NZ, INDEX,
     $          VALUE, LVALUE,
     $          INDEX (N+2), LIND2-(N+1),
     $          LUX1, LUI1, IUSE, XUSE, NZOFF, NBLKS,
     $          ICNTL, INFO, RINFO,
     $          PRESRV, INDEX, INDEX, VALUE, 0, 1,
     $          INDEX (LUIR1), IP2 - LUIR1 + 1,
     $          INDEX (LUBLPP), INDEX (BLKPP), INDEX (OFFPP), ON,
     $          INDEX (CPERMP), INDEX (RPERMP), NE)
           IF (INFO (1) .LT. 0) THEN
              GO TO 9000
           ENDIF
           LUI1 = LUI1 + (N+1)
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (TRANSA) THEN
           INDEX (LINDEX-6) = 1
        ELSE
           INDEX (LINDEX-6) = 0
        ENDIF
        INDEX (LINDEX-5) = NZOFF
        INDEX (LINDEX-4) = NBLKS
        IF (PRESRV) THEN
           INDEX (LINDEX-3) = 1
        ELSE
           INDEX (LINDEX-3) = 0
        ENDIF
        INDEX (LINDEX-2) = NZ
        INDEX (LINDEX-1) = N
        INDEX (LINDEX) = NE
        KEEP (1) = LUX1
        KEEP (2) = LVALUE
        KEEP (3) = LUI1
        KEEP (4) = LUIR1
        KEEP (5) = LINDEX
        IUSE = LINDEX - LUI1 + 1
        XUSE = LVALUE - LUX1 + 1
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
9000    CONTINUE
        IF (INFO (1) .LT. 0) THEN
           KEEP (1) = 0
           KEEP (2) = 0
           KEEP (3) = 0
           KEEP (4) = 0
           KEEP (5) = 0
        ENDIF
        INFO (18) = MIN (LINDEX, MAX (INFO (18), IUSE))
        INFO (19) = INFO (18)
        INFO (22) = INFO (19)
        INFO (20) = MIN (LVALUE, MAX (INFO (20), XUSE))
        CALL MA38YD (2, 2,
     $          N, NE, JOB, TRANSA, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          RINFO(1), RINFO(1), 1, RINFO(1), 1)
        RETURN
        END
        SUBROUTINE MA38CD (N, JOB, TRANSC, LVALUE, LINDEX, VALUE, INDEX,
     *                     KEEP, B, X, W, CNTL, ICNTL, INFO, RINFO)
        INTEGER          N, JOB
        LOGICAL          TRANSC
        INTEGER          LVALUE, LINDEX
        DOUBLE PRECISION VALUE(LVALUE)
        INTEGER          INDEX(LINDEX), KEEP(20)
        DOUBLE PRECISION B(N), X(N), W(*), CNTL(10)
        INTEGER          ICNTL(20), INFO(40)
        DOUBLE PRECISION RINFO(20)
C=== MA38CD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C## End of user documentation for MA38CD ###############################
C=======================================================================
C=======================================================================
        INTRINSIC MAX
        EXTERNAL   MA38JD, MA38ND, MA38YD
C=======================================================================
C=======================================================================
        INTEGER NBLKS, OFFIP, OFFXP, N1, NZ, NE, OFFPP, BLKPP, LUBLPP,
     $          APP, AN, ANZ, ON, LUI1, LUI2, LUX1, LUX2, AIP, AXP,
     $          CPERMP, RPERMP, NZOFF, IRSTEP, YP, LY, LW, SP, IP1, IP2,
     $          XP1, LUIR1, IO, PRL
        LOGICAL PRESRV, BADLU
        DOUBLE PRECISION
     $          ZERO
        PARAMETER (ZERO = 0.0D0)
C=======================================================================
C=======================================================================
        IO = ICNTL (2)
        PRL = ICNTL (3)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        INFO (1) = 0
        INFO (24) = 0
        RINFO (7) = ZERO
        RINFO (8) = ZERO
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IRSTEP = MAX (0, ICNTL (8))
        IF (IRSTEP .EQ. 0) THEN
           LW = 2*N
        ELSE
           LW = 4*N
        ENDIF
        CALL MA38YD (3, 1,
     $          N, NE, JOB, TRANSC, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          B, X, N, W, LW)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        LUX1 = KEEP (1)
        LUX2 = KEEP (2)
        LUI1 = KEEP (3)
        LUIR1 = KEEP (4)
        LUI2 = KEEP (5)
        BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1
     $     .OR. LUI2 .GT. LINDEX
     $     .OR. LUX1 .LE. 0 .OR. LUX1 .GT. LUX2 .OR. LUX2 .GT. LVALUE
     $     .OR. LUI1 .LE. 0 .OR. LUIR1 .LT. LUI1 .OR. LUIR1 .GT. LUI2
        IF (BADLU) THEN
           CALL MA38ND (3, ICNTL, INFO, -7, 0)
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        NE = INDEX (LUI2)
        N1 = INDEX (LUI2-1)
        NZ = INDEX (LUI2-2)
        PRESRV = INDEX (LUI2-3) .NE. 0
        NBLKS = INDEX (LUI2-4)
        NZOFF = INDEX (LUI2-5)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        RPERMP = (LUI2-6) - N
        CPERMP = RPERMP - N
        IP2 = CPERMP - 1
        XP1 = LUX1
        IP1 = LUI1
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (PRESRV) THEN
           APP = IP1
           AIP = APP + N+1
           IP1 = AIP + NZ
           AXP = XP1
           XP1 = AXP + NZ
           AN = N
           ANZ = NZ
        ELSE
           APP = 1
           AIP = 1
           AXP = 1
           AN = 1
           ANZ = 1
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (NBLKS .GT. 1) THEN
           OFFIP = IP1
           IP1 = IP1 + NZOFF
           OFFXP = XP1
           XP1 = XP1 + NZOFF
           OFFPP = CPERMP - (N+1)
           BLKPP = OFFPP - (NBLKS+1)
           LUBLPP = BLKPP - (NBLKS)
           IP2 = LUBLPP - 1
           ON = N
        ELSE
           OFFIP = 1
           OFFXP = 1
           OFFPP = 1
           BLKPP = 1
           LUBLPP = 1
           ON = 1
        ENDIF
        BADLU = N .NE. N1 .OR. NZ .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $     NBLKS .LE. 0 .OR. NBLKS .GT. N .OR.
     $     XP1 .GT. LUX2 .OR. NZOFF .LT. 0 .OR. IP1 .NE. LUIR1
        IF (BADLU) THEN
           CALL MA38ND (3, ICNTL, INFO, -7, 0)
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (IRSTEP .GT. 0 .AND. .NOT. PRESRV) THEN
           CALL MA38ND (3, ICNTL, INFO, 8, 0)
           IRSTEP = 0
        ENDIF
        IF (IRSTEP .GT. 0 .AND. (JOB .EQ. 1 .OR. JOB .EQ. 2)) THEN
           CALL MA38ND (3, ICNTL, INFO, 8, 1)
           IRSTEP = 0
        ENDIF
        IF (IRSTEP .EQ. 0) THEN
           YP = 1
           LY = 1
           SP = 1
           LW = 2*N
        ELSE
           YP = 2*N+1
           LY = N
           SP = 3*N+1
           LW = 4*N
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        CALL MA38JD (N, JOB, TRANSC, LUX2-XP1+1, VALUE (XP1),
     $     IP2-LUIR1+1, INDEX (LUIR1), B, X,
     $     W, W (N+1), LY, W (YP), W (SP),
     $     CNTL, INFO, RINFO, INDEX (CPERMP), INDEX (RPERMP),
     $     AN, ANZ, INDEX (APP), INDEX (AIP), VALUE (AXP),
     $     ON, MAX (1, NZOFF), INDEX (OFFPP), INDEX (OFFIP),
     $     VALUE (OFFXP), NBLKS, INDEX (LUBLPP), INDEX (BLKPP), IRSTEP)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
9000    CONTINUE
        CALL MA38YD (3, 2,
     $          N, NE, JOB, TRANSC, LVALUE, LINDEX, VALUE,
     $          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     $          B, X, N, W, LW)
        RETURN
        END
C#######################################################################
C##  non-user-callable routines:
C#######################################################################
        SUBROUTINE MA38DD (N, NZ, CP, XX, XSIZE, II, ISIZE, XTAIL,
     *          ITAIL, IUSE, XUSE, NZOFF, NBLKS, ICNTL, CNTL, INFO,
     *          RINFO, PRESRV, AP, AI, AX, AN, ANZ, KEEP, NE)
        INTEGER N, NZ, ISIZE, II(ISIZE), ICNTL(20), INFO(40),
     *          CP(N+1), XSIZE, XTAIL, ITAIL, IUSE, XUSE, AN, ANZ,
     *          AP(AN+1), AI(ANZ), KEEP(20), NZOFF, NBLKS, NE
        LOGICAL PRESRV
        DOUBLE PRECISION XX(XSIZE), CNTL(10), RINFO(20), AX(ANZ)
C=== MA38DD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC MAX
        EXTERNAL MA38ED, MA38HD, MA38MD, MA38ND
C=======================================================================
C=======================================================================
        INTEGER KN, NZDIA, BLKPP, LUBLPP, P, OFFIP, XHEAD, ROW,
     $          OFFXP, OFFPP, IHEAD, K1, K2, BLK, PRP, P2, CPERMP,
     $          RPERMP, NSGLTN, NPIV, MNZ, NSYM, K, COL, RMAX, CMAX,
     $          TOTNLU, XRMAX, XRUSE
        LOGICAL TRYBTF, IOUT, XOUT
        DOUBLE PRECISION
     $          ZERO, ONE, A
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        NBLKS = 1
        NZOFF = 0
        NZDIA = NZ
        NSGLTN = 0
        NPIV = 0
        RMAX = 1
        CMAX = 1
        TOTNLU = 0
        XOUT = .FALSE.
        IOUT = .FALSE.
        XRMAX = 2*NE
        IF (PRESRV) THEN
           IHEAD = 1
           XHEAD = 1
        ELSE
           IHEAD = NZ + 1
           XHEAD = NZ + 1
        ENDIF
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        ITAIL = ITAIL - (2*N+7)
        IUSE = IUSE + (2*N+7)
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = INFO (18)
        CPERMP = ITAIL
        RPERMP = CPERMP + N
        IF (IHEAD .GT. ITAIL) THEN
           GO TO 9000
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        TRYBTF = ICNTL (4) .EQ. 1
        IF (TRYBTF) THEN
           ITAIL = ITAIL - (N+1)
           OFFPP = ITAIL
           ITAIL = ITAIL - (5*N+1)
           P = ITAIL
           IUSE = IUSE + (6*N+2)
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = INFO (18)
           IF (PRESRV) THEN
              IF (IHEAD .GT. ITAIL) THEN
                 GO TO 9000
              ENDIF
              CALL MA38HD (AX, ANZ, AI, ANZ, N, NZ, NZDIA, NZOFF,
     $           NBLKS, CP, II (CPERMP), II (RPERMP), II(P), II(P+N),
     $           II (P+2*N), II (P+3*N), II (P+4*N), II (OFFPP),
     $           PRESRV)
           ELSE
              IHEAD = IHEAD + NZ
              XHEAD = XHEAD + NZ
              IUSE = IUSE + NZ
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (19) = INFO (18)
              INFO (21) = INFO (20)
              IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
                 GO TO 9000
              ENDIF
              CALL MA38HD (XX, 2*NZ, II, 2*NZ, N, NZ, NZDIA, NZOFF,
     $              NBLKS, CP, II (CPERMP), II (RPERMP), II(P), II(P+N),
     $              II (P+2*N), II (P+3*N), II (P+4*N), II (OFFPP),
     $              PRESRV)
              IHEAD = IHEAD - NZ
              XHEAD = XHEAD - NZ
              IUSE = IUSE - NZ
              XUSE = XUSE - NZ
           ENDIF
           IF (NBLKS .GT. 1) THEN
              BLKPP = OFFPP - (NBLKS+1)
              LUBLPP = BLKPP - (NBLKS)
              ITAIL = LUBLPP
              IUSE = IUSE - (6*N+2) + (2*NBLKS+N+2)
           ELSE
              ITAIL = (ISIZE + 1) - (2*N+7)
              IUSE = IUSE - (6*N+2)
           ENDIF
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        XRUSE = NZ
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (NBLKS .EQ. 1) THEN
           CALL MA38ED (CP, N, II (CPERMP), II (RPERMP), NZOFF,
     $          ITAIL, XTAIL, XX, XSIZE, XUSE, II, ITAIL-1, IUSE,
     $          ICNTL, CNTL, INFO, RINFO, NBLKS,
     $          AP, AI, AX, PRESRV, 1, AN, ANZ, II, KEEP,
     $          RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
           IF (IOUT .OR. XOUT) THEN
              GO TO 9000
           ENDIF
           IF (INFO (1) .LT. 0) THEN
              GO TO 9010
           ENDIF
           IHEAD = 1
           XHEAD = 1
           II (ITAIL) = 1
        ELSE
           PRP = OFFPP
           IF (PRESRV) THEN
              NZOFF = 0
CFPP$ NODEPCHK L
              DO 10 K = 1, N
                 II (PRP + II (RPERMP+K-1) - 1) = K
10            CONTINUE
           ENDIF
           DO 30 BLK = NBLKS, 1, -1
              K1 = II (BLKPP+BLK-1)
              K2 = II (BLKPP+BLK) - 1
              KN = K2-K1+1
              IF (.NOT. PRESRV) THEN
                 P = CP (K1)
                 CP (K2+1) = IHEAD
              ENDIF
              IF (KN .GT. 1) THEN
                 CALL MA38ED (CP (K1), KN,
     $              II (CPERMP+K1-1), II (RPERMP+K1-1), NZOFF,
     $              ITAIL, XTAIL, XX, XTAIL-1, XUSE, II, ITAIL-1,
     $              IUSE, ICNTL, CNTL, INFO, RINFO, NBLKS,
     $              AP, AI, AX, PRESRV, K1, AN, ANZ, II (PRP), KEEP,
     $              RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
                 IF (IOUT .OR. XOUT) THEN
                    GO TO 9000
                 ENDIF
                 IF (INFO (1) .LT. 0) THEN
                    GO TO 9010
                 ENDIF
                 IF (PRESRV) THEN
                    IHEAD = 1
                    XHEAD = 1
                 ELSE
                    IHEAD = P
                    XHEAD = P
                 ENDIF
                 II (LUBLPP+BLK-1) = ITAIL
              ELSE
                 A = ZERO
                 IF (PRESRV) THEN
                    COL = II (CPERMP + K1 - 1)
                    DO 20 P2 = AP (COL), AP (COL + 1) - 1
                       ROW = II (PRP + AI (P2) - 1)
                       IF (ROW .LT. K1) THEN
                          NZOFF = NZOFF + 1
                       ELSE
                          A = AX (P2)
                       ENDIF
20                  CONTINUE
                    IHEAD = 1
                    XHEAD = 1
                 ELSE IF (P .NE. IHEAD) THEN
                    A = XX (P)
                    IHEAD = P
                    XHEAD = P
                    IUSE = IUSE - 1
                    XUSE = XUSE - 1
                    XRUSE = XRUSE - 1
                 ENDIF
                 NSGLTN = NSGLTN + 1
                 IF (A .EQ. ZERO) THEN
                    A = ONE
                 ELSE
                    NPIV = NPIV + 1
                 ENDIF
                 XTAIL = XTAIL - 1
                 XUSE = XUSE + 1
                 XRUSE = XRUSE + 1
                 XRMAX = MAX (XRMAX, XRUSE)
                 INFO (20) = MAX (INFO (20), XUSE)
                 INFO (21) = MAX (INFO (21), XUSE)
                 IF (XHEAD .GT. XTAIL) THEN
                    GO TO 9000
                 ENDIF
                 II (LUBLPP+BLK-1) = -XTAIL
                 XX (XTAIL) = A
              ENDIF
30         CONTINUE
CFPP$ NODEPCHK L
           DO 40 P = LUBLPP, LUBLPP + NBLKS - 1
              IF (II (P) .GT. 0) THEN
                 II (II (P)) = II (II (P)) - XTAIL + 1
                 II (P) = II (P) - ITAIL + 1
              ELSE
                 II (P) = (-II (P)) - XTAIL + 1
              ENDIF
40         CONTINUE
           PRP = IHEAD
           IHEAD = IHEAD + N
           IUSE = IUSE + N
           IF (NBLKS .EQ. N) THEN
              ITAIL = ITAIL - 1
              IUSE = IUSE + 1
              P2 = ITAIL
           ENDIF
           ITAIL = ITAIL - NZOFF
           OFFIP = ITAIL
           XTAIL = XTAIL - NZOFF
           OFFXP = XTAIL
           IUSE = IUSE + NZOFF
           XUSE = XUSE + NZOFF
           XRUSE = XRUSE + NZOFF
           XRMAX = MAX (XRMAX, XRUSE)
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), IUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XUSE)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
              GO TO 9000
           ENDIF
           MNZ = NZOFF
           IF (PRESRV) THEN
              CALL MA38MD (CP, N, II (RPERMP), II (CPERMP), NZOFF,
     $          II (OFFPP), II (OFFIP), XX (OFFXP), II (PRP),
     $          ICNTL, AP, AI, AX, AN, ANZ, PRESRV, NBLKS, II (BLKPP),
     $          MNZ, 1, P)
           ELSE
              CALL MA38MD (CP, N, II (RPERMP), II (CPERMP), NZOFF,
     $          II (OFFPP), II (OFFIP), XX (OFFXP), II (PRP),
     $          ICNTL, AP, II, XX, 0, MNZ, PRESRV, 0, II(BLKPP),
     $          MNZ, 1, P)
           ENDIF
           IF (NBLKS .EQ. N) THEN
              II (P2) = 0
           ENDIF
           IHEAD = 1
           XHEAD = 1
           IUSE = IUSE - N
           IF (.NOT. PRESRV) THEN
              IUSE = IUSE - NZOFF
              XUSE = XUSE - NZOFF
           ENDIF
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
9000    CONTINUE
        IF (IOUT .OR. IHEAD .GT. ITAIL) THEN
           CALL MA38ND (1, ICNTL, INFO, -3, INFO (19))
        ENDIF
        IF (XOUT .OR. XHEAD .GT. XTAIL) THEN
           CALL MA38ND (1, ICNTL, INFO, -4, INFO (21))
        ENDIF
9010    CONTINUE
        INFO (4) = 0
        NZDIA = NZ - NZOFF
        INFO (5) = NZ
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (8) = NSGLTN
        INFO (9) = NBLKS
        INFO (12) = INFO (10) + INFO (11) + N + INFO (7)
        NSYM = 0
        IF (INFO (1) .GE. 0) THEN
           DO 50 K = 1, N
              IF (II (CPERMP+K-1) .EQ. II (RPERMP+K-1)) THEN
                 NSYM = NSYM + 1
              ENDIF
50         CONTINUE
        ENDIF
        INFO (16) = NSYM
        INFO (17) = INFO (17) + NPIV
        RINFO (1) = RINFO (4) + RINFO (5) + RINFO (6)
        IF (INFO (1) .GE. 0 .AND. INFO (17) .LT. N) THEN
           CALL MA38ND (1, ICNTL, INFO, 4, INFO (17))
        ENDIF
        IF (PRESRV) THEN
           INFO (22) = MAX (3*NE+2*N+1, NE+3*N+2,
     $                           2*NZ+4*N+10+RMAX+3*CMAX+4*TOTNLU)
        ELSE
           INFO (22) = MAX (3*NE+2*N+1, NE+3*N+2, 2*NZ+3*N+2,
     $                             NZ+3*N+ 9+RMAX+3*CMAX+4*TOTNLU)
        ENDIF
        INFO (23) = XRMAX
        RETURN
        END
        SUBROUTINE MA38ED (CP, N, CPERM, RPERM, NZOFF,
     *          ITAIL, XTAIL, XX, XSIZE, XUSE, II, ISIZE, IUSE,
     *          ICNTL, CNTL, INFO, RINFO, NBLKS,
     *          AP, AI, AX, PRESRV, K1, AN, ANZ, PR, KEEP,
     *          RMAX, CMAX, TOTNLU, XRMAX, XRUSE, IOUT, XOUT)
        INTEGER XSIZE, ISIZE, N, ICNTL(20), INFO(40), XUSE, IUSE,
     *          ITAIL, XTAIL, II(ISIZE), CP(N+1), CPERM(N), NZOFF,
     *          AN, ANZ, RPERM(N), AI(ANZ), AP(AN+1), K1, PR(AN),
     *          NBLKS, KEEP(20), RMAX, CMAX, TOTNLU, XRMAX, XRUSE
        LOGICAL PRESRV, IOUT, XOUT
        DOUBLE PRECISION XX(XSIZE), CNTL(10), RINFO(20), AX(ANZ)
C=== MA38ED ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C			and if the corresponding error status has not
C			yet been set (the error status is set in the caller,
C			MA38DD)
C			and if the corresponding error status has not
C			yet been set (the error status is set in the caller,
C			MA38DD)
C=======================================================================
C=======================================================================
        INTRINSIC MAX, SQRT
        EXTERNAL MA38FD
C=======================================================================
C=======================================================================
        INTEGER CP1, PC, PEND, PCOL, CDEG, COL, CSIZ, NZ, XP, IP, IS, P,
     $          DN, DSIZ, WRKSIZ, I, CLEN, D1, D2, N2, ROW, CSCAL
        PARAMETER (CSCAL = 9)
        DOUBLE PRECISION XN
C=======================================================================
C=======================================================================
        IOUT = .FALSE.
        XOUT = .FALSE.
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        D1 = KEEP (7)
        D2 = KEEP (8)
        XN = N
        XN = SQRT (XN)
        N2 = XN
        DSIZ = MAX (0, D1, D2 * N2)
        DN = 0
        IF (PRESRV) THEN
           IF (NBLKS .EQ. 1) THEN
              DO 10 COL = 1, N
                 IF (AP (COL+1) - AP (COL) .GT. DSIZ) THEN
                    DN = DN + 1
                 ENDIF
10            CONTINUE
           ELSE
              DO 40 COL = 1, N
                 CDEG = AP (CPERM (COL) + 1)- AP (CPERM (COL))
                 IF (CDEG .GT. DSIZ) THEN
                    CDEG = 0
                    DO 20 P = AP (CPERM (COL)), AP (CPERM (COL) + 1) -1
                       ROW = PR (AI (P))
                       IF (ROW .GE. K1) THEN
                          CDEG = CDEG + 1
                          IF (CDEG .GT. DSIZ) THEN
                             DN = DN + 1
                             GO TO 30
                          ENDIF
                       ENDIF
20                  CONTINUE
30                  CONTINUE
                 ENDIF
40            CONTINUE
           ENDIF
        ELSE
           DO 50 COL = 1, N
              IF (CP (COL+1) - CP (COL) .GT. DSIZ) THEN
                 DN = DN + 1
              ENDIF
50         CONTINUE
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (NBLKS .EQ. 1) THEN
           WRKSIZ = 8*N + 3*DN
        ELSE
           WRKSIZ = 10*N + 3*DN
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (PRESRV) THEN
           CP1 = 1
           XP = 1
           IP = 1 + WRKSIZ
           IF (NBLKS .EQ. 1) THEN
              NZ = ANZ
              IS = NZ + WRKSIZ + CSCAL*N
              IUSE = IUSE + IS
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (19) = MAX (INFO (19), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)
              IOUT = IS .GT. ISIZE
              XOUT = NZ .GT. XSIZE
              IF (IOUT .OR. XOUT) THEN
                 GO TO 9000
              ENDIF
              PC = IP
              DO 70 COL = 1, N
                 CP (COL) = PC - WRKSIZ
                 CDEG = AP (COL+1) - AP (COL)
                 CLEN = CDEG
                 CSIZ = CDEG + CSCAL
                 II (PC) = CSIZ
                 II (PC+1) = CDEG
                 II (PC+5) = 0
                 II (PC+6) = CLEN
                 II (PC+7) = 0
                 II (PC+8) = 0
                 II (PC+2) = XP
                 XP = XP + CDEG
                 PC = PC + CSCAL
                 P = AP (COL)
                 DO 60 I = 0, CDEG - 1
                    II (PC + I) = AI (P + I)
60               CONTINUE
                 PC = PC + CDEG
70            CONTINUE
              DO 80 P = 1, NZ
                 XX (P) = AX (P)
80            CONTINUE
           ELSE
              DO 100 COL = 1, N
                 PC = IP
                 CP (COL) = PC - WRKSIZ
                 IP = IP + CSCAL
                 IOUT = IP .GT. ISIZE
                 IF (IOUT) THEN
                    GO TO 9000
                 ENDIF
                 II (PC+2) = XP
                 CDEG = IP
                 DO 90 P = AP (CPERM (COL)), AP (CPERM (COL)+1)-1
                    ROW = PR (AI (P))
                    IF (ROW .GE. K1) THEN
                       IOUT = IP .GT. ISIZE
                       XOUT = XP .GT. XSIZE
                       IF (IOUT .OR. XOUT) THEN
                          GO TO 9000
                       ENDIF
                       II (IP) = ROW - K1 + 1
                       XX (XP) = AX (P)
                       IP = IP + 1
                       XP = XP + 1
                    ELSE
                       NZOFF = NZOFF + 1
                    ENDIF
90               CONTINUE
                 CDEG = IP - CDEG
                 CLEN = CDEG
                 CSIZ = CDEG + CSCAL
                 II (PC) = CSIZ
                 II (PC+1) = CDEG
                 II (PC+5) = 0
                 II (PC+6) = CLEN
                 II (PC+7) = 0
                 II (PC+8) = 0
100           CONTINUE
              NZ = XP - 1
              IS = NZ + WRKSIZ + CSCAL*N
              IUSE = IUSE + IS
              XUSE = XUSE + NZ
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (19) = MAX (INFO (19), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)
           ENDIF
           XRUSE = XRUSE + NZ
           XRMAX = MAX (XRMAX, XRUSE)
        ELSE
           CP1 = CP (1)
           NZ = CP (N+1) - CP1
           PC = CP1 + WRKSIZ + (NZ+CSCAL*N)
           IUSE = IUSE + WRKSIZ + CSCAL*N
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), IUSE)
           IOUT = PC .GT. ISIZE+1
           IF (IOUT) THEN
              GO TO 9000
           ENDIF
           XP = NZ + 1
           IP = NZ + CSCAL*N + 1
           PEND = CP (N+1)
           DO 120 COL = N, 1, -1
              PCOL = CP (COL)
              DO 110 P = PEND-1, PCOL, -1
                 PC = PC - 1
                 II (PC) = II (P)
110           CONTINUE
              PC = PC - CSCAL
              CDEG = PEND - PCOL
              CLEN = CDEG
              PEND = PCOL
              CSIZ = CDEG + CSCAL
              IP = IP - CSIZ
              CP (COL) = IP
              II (PC) = CSIZ
              II (PC+1) = CDEG
              II (PC+5) = 0
              II (PC+6) = CLEN
              II (PC+7) = 0
              II (PC+8) = 0
              XP = XP - CDEG
              II (PC+2) = XP
120        CONTINUE
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        XP = CP1
        IP = CP1 + WRKSIZ
        IF (NBLKS .EQ. 1) THEN
           CALL MA38FD (CP, NZ, N, 1, CPERM, RPERM, ITAIL, XTAIL,
     $          XX (XP), XSIZE-XP+1, II (IP), ISIZE-IP+1, ICNTL, CNTL,
     $          INFO, RINFO, .FALSE., IUSE, XUSE,
     $          RPERM, CPERM, II (CP1), II (CP1+N),
     $          II (CP1+2*N), II (CP1+3*N), II (CP1+4*N), II (CP1+5*N),
     $          II (CP1+6*N+DN), II (CP1+7*N+2*DN),
     $          DN, DSIZ, RMAX, CMAX, TOTNLU, XRMAX, XRUSE)
        ELSE
           CALL MA38FD (CP, NZ, N, N, CPERM, RPERM, ITAIL, XTAIL,
     $          XX (XP), XSIZE-XP+1, II (IP), ISIZE-IP+1, ICNTL, CNTL,
     $          INFO, RINFO, .TRUE., IUSE, XUSE,
     $          II (CP1), II (CP1+N), II (CP1+2*N), II (CP1+3*N),
     $          II (CP1+4*N), II (CP1+5*N), II (CP1+6*N), II (CP1+7*N),
     $          II (CP1+8*N+DN), II (CP1+9*N+2*DN),
     $          DN, DSIZ, RMAX, CMAX, TOTNLU, XRMAX, XRUSE)
        ENDIF
        IF (INFO (1) .LT. 0) THEN
           RETURN
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IUSE = IUSE - WRKSIZ
        ITAIL = ITAIL + IP - 1
        XTAIL = XTAIL + XP - 1
        II (ITAIL) = XTAIL
        RETURN
C=======================================================================
C=======================================================================
9000    CONTINUE
        RETURN
        END
        SUBROUTINE MA38FD (CP, NZ, N, PN, CPERM, RPERM, ITAIL, XTAIL,
     *          XX, XSIZE, II, ISIZE, ICNTL, CNTL, INFO, RINFO, PGIVEN,
     *          IUSE, XUSE, WIR, WIC, WPR, WPC, WM, HEAD,
     *          WJ, RP, WC, WR, DN, DSIZ,
     *          RMAX, CMAX, TOTNLU, XRMAX, XRUSE)
        INTEGER XSIZE, ISIZE, ICNTL(20), INFO(40), PN,
     *          ITAIL, XTAIL, NZ, N, II(ISIZE), CP(N+1), DN, DSIZ,
     *          RPERM(PN), CPERM(PN), WIR(N), WIC(N), WPR(N),
     *          WPC(N), WM(N), HEAD(N), RP(N+DN), WC(N+DN),
     *          WR(N+DN), IUSE, XUSE, WJ(N),
     *          RMAX, CMAX, TOTNLU, XRMAX, XRUSE
        LOGICAL PGIVEN
        DOUBLE PRECISION XX(XSIZE), CNTL(10), RINFO(20)
C=== MA38FD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTEGER IDAMAX
        INTRINSIC ABS, MAX, MIN
        EXTERNAL MA38GD, MA38ND, DGEMM, DGEMV, IDAMAX
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C=======================================================================
C=======================================================================
        INTEGER SWPCOL, SWPROW, FDIMC, K0, COLPOS, ROWPOS, ROW2, RDEG2,
     $          P, I, J, FFROW, PIVROW, PIVCOL, LUDEGR, LUDEGC, E1,
     $          FXP, LURP, LUCP, IP, NEXT, FFLEFR, PC, MNEXT, MPREV,
     $          FFLEFC, FEDEGR, FEDEGC, K, XUDP, XDP, XSP, XLP, S, COL2,
     $          BESTCO, COL, E, ROW, COST, SRCHED, PR, F1, RSCAN, REP,
     $          KLEFT1, FFSIZE, FFXP, W0, FFDIMR, FFDIMC, KLEFT, XLDP,
     $          EP, SCAN1, SCAN2, SCAN3, SCAN4, NZL, NZU, DEGC, CEP
        INTEGER MINDEG, NSRCH, NPIV, ESON, LUIP1, DNZ, IWORST, WXP,
     $          NB, LUPP, NLU, NSONS, INEED, XNEED, LDIMC, LXP, RLEN2,
     $          RSIZ, LSONS, SONLST, XHEAD, IHEAD, DELN, DLEN,
     $          SLIST, XP, LUIP, RDEG, CDEG1, PFREE, XFREE, CDEG2,
     $          F, CDEG, MTAIL, MHEAD, RSIZ2, CSIZ2, IP2, MAXDR, MAXDC,
     $          XS, IS, LUXP, FSP, FLP, FDP, JJ, USONS, NDN, P2,
     $          CSIZ, CELN, CLEN, RELN, RLEN, UXP, PC2, PR2
        INTEGER CNXT, CPRV, CXP, FLUIP, LUSONP, FLEFTR, FLEFTC,
     $          FMAXR, FMAXC, SLOTS, LIMIT, RSCAL, CSCAL, FSCAL, EXTRA,
     $          FMAX, MINMEM, DUMMY1, DUMMY2, DUMMY3, DUMMY4, W1
        LOGICAL SYMSRC, PFOUND, MOVELU, OKCOL, OKROW, BETTER
        DOUBLE PRECISION
     $          TOLER, MAXVAL, RELPT, GRO, ONE, ZERO, X
        PARAMETER (ONE = 1.0D0, ZERO = 0.0D0,
     $          RSCAL = 2, CSCAL = 9, FSCAL = 7,
     $          MINMEM = 24)
        INTEGER INTZER
C=======================================================================
C=======================================================================
        NSRCH = MAX (1, ICNTL (5))
        SYMSRC = ICNTL (6) .NE. 0
        NB = MAX (1, ICNTL (7))
        RELPT = MAX (ZERO, MIN (CNTL (1), ONE))
        GRO = MAX (ONE, CNTL (2))
        NDN = N + DN
        W0 = NDN + 2
        KLEFT = N
        NPIV = 0
        NLU = 0
        MINDEG = 1
        FMAX = 1
        IHEAD = NZ + CSCAL*N + 1
        XHEAD = NZ + 1
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1
        XFREE = -1
        PFREE = 0
        INFO (19) = MAX (INFO (19), IUSE+MINMEM)
        IF (IHEAD.GT.ITAIL.OR.ISIZE.LT.MINMEM.OR.XHEAD.GT.XTAIL) THEN
           GO TO 9000
        ENDIF
        BESTCO = 0
        LIMIT = N + 2*NDN
        LSONS = NDN + 1
        USONS = NDN + 1
        DO 10 I = 1, N
           WIR (I) = -1
           WIC (I) = -2
           HEAD (I) = 0
           WC (I) = 0
           WR (I) = 0
10      CONTINUE
        MHEAD = 0
        MTAIL = 0
        DO 20 COL = 1, N
           PC = CP (COL)
           CLEN = II (PC+6)
           IF (CLEN .GT. 0) THEN
              IF (MHEAD .EQ. 0) MHEAD = PC
              II (PC+4) = MTAIL
              II (PC+3) = 0
              IF (MTAIL .NE. 0) II (MTAIL+3) = PC
              MTAIL = PC
           ELSE
              II (PC+2) = 0
              II (PC+4) = 0
              II (PC+3) = 0
           ENDIF
20      CONTINUE
        E = N
        DNZ = 0
        DO 50 COL = 1, N
           PC = CP (COL)
           CLEN = II (PC+6)
           CEP = (PC+9)
           IF (CLEN .GT. DSIZ) THEN
              DNZ = DNZ + CLEN
              DO 30 IP = CEP, CEP + CLEN - 1
                 ROW = II (IP)
                 WR (ROW) = WR (ROW) + 1
30            CONTINUE
              E = E + 1
              EP = PC
              RP (E) = EP
              FDIMC = CLEN
              FLEFTC = CLEN
              FLEFTR = 1
              II (EP+1) = FDIMC
              II (EP+5) = FLEFTR
              II (EP+6) = FLEFTC
              WR (E) = W0-1
              WC (E) = W0-1
              LURP = (EP+8)
              II (LURP) = COL
              FMAX = MAX (FMAX, FLEFTC)
           ELSE
              DO 40 IP = CEP, CEP + CLEN - 1
                 ROW = II (IP)
                 WC (ROW) = WC (ROW) + 1
40            CONTINUE
           ENDIF
50      CONTINUE
        PR = IHEAD
        CSIZ = CSCAL + 2
        IS = (NZ + RSCAL*N + DNZ) + (DN * CSIZ)
        IHEAD = IHEAD + IS
        IUSE = IUSE + IS
        INEED = IUSE
        XNEED = XUSE
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .GT. ITAIL) THEN
           GO TO 9000
        ENDIF
        IF (DN .GT. 0) THEN
           EXTRA = MIN ((ITAIL - IHEAD) / DN, DSIZ + 6)
           CSIZ = CSIZ + EXTRA
           IS = DN * EXTRA
           IHEAD = IHEAD + IS
           IUSE = IUSE + IS
           INFO (18) = MAX (INFO (18), IUSE)
        ENDIF
        DO 60 ROW = 1, N
           RP (ROW) = PR
           REP  = (PR+2)
           RELN = WR (ROW)
           RLEN = WC (ROW)
           RSIZ = 2*RELN + RLEN + RSCAL
           II (PR) = RSIZ
           RDEG = RELN + RLEN
           II (PR+1) = RDEG
           WM (ROW) = REP
           PR = PR + RSIZ
60      CONTINUE
        PC = PR
        DO 80 E = N+1, N+DN
           EP = RP (E)
           LUCP = (EP+9)
           FDIMC = II (EP+1)
CFPP$ NODEPCHK L
           DO 70 F = 0, FDIMC - 1
              ROW = II (LUCP+F)
              II (WM (ROW)    ) = E
              II (WM (ROW) + 1) = F
              WM (ROW) = WM (ROW) + 2
70         CONTINUE
           LURP = (EP+8)
           COL = II (LURP)
           CP (COL) = PC
           II (PC) = CSIZ
           CDEG = FDIMC
           II (PC+1) = CDEG
           II (PC+2) = 0
           II (PC+4) = 0
           II (PC+3) = 0
           II (PC+5) = 1
           II (PC+6) = 0
           II (PC+7) = 0
           II (PC+8) = 0
           CEP = (PC+9)
           II (CEP  ) = E
           II (CEP+1) = 0
           PC = PC + CSIZ
80      CONTINUE
        DO 100 COL = 1, N
           PC = CP (COL)
           CEP = (PC+9)
           CLEN = II (PC+6)
CFPP$ NODEPCHK L
           DO 90 P = CEP, CEP + CLEN - 1
              ROW = II (P)
              II (WM (ROW)) = COL
              WM (ROW) = WM (ROW) + 1
90         CONTINUE
100     CONTINUE
        RINFO (2) = RINFO (2) + NZ
        DO 110 COL = N, 1, -1
           PC = CP (COL)
           CDEG = II (PC+1)
           IF (CDEG .LE. 0) THEN
              CDEG = -(N+2)
              II (PC+1) = CDEG
           ELSE
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) II (CP (CNXT)+8) = COL
              HEAD (CDEG) = COL
           ENDIF
110     CONTINUE
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        DO 1540 DUMMY1 = 1, N
        IF (NPIV .GE. N) THEN
           GO TO 2000
        ENDIF
C=======================================================================
C=======================================================================
        IF (MTAIL .NE. 0 .AND. II (MTAIL+6) .EQ. 0) THEN
           XP = II (MTAIL+2)
           XUSE = XUSE - (XHEAD - XP)
           XHEAD = XP
           IF (MTAIL .EQ. PFREE) THEN
              PFREE = 0
              XFREE = -1
           ENDIF
           MTAIL = II (MTAIL+4)
           IF (MTAIL .NE. 0) THEN
              II (MTAIL+3) = 0
           ELSE
              MHEAD = 0
           ENDIF
        ENDIF
C=======================================================================
C=======================================================================
        NSONS = 0
        SONLST = 0
        SRCHED = 0
        PIVCOL = 0
        SLIST = 0
        DO 255 DUMMY2 = 1, N
           COL = 0
           DO 140 CDEG = MINDEG, N
              COL = HEAD (CDEG)
              IF (COL .NE. 0) THEN
                 GO TO 150
              ENDIF
140        CONTINUE
           IF (COL .EQ. 0) THEN
              GO TO 260
           ENDIF
150        CONTINUE
           PC = CP (COL)
           CNXT = II (PC+7)
           IF (CNXT .NE. 0) THEN
              II (CP (CNXT)+8) = 0
           ENDIF
           HEAD (CDEG) = CNXT
           MINDEG = CDEG
           XS = CDEG
           IF (XS .GT. XTAIL-XHEAD) THEN
              INFO (15) = INFO (15) + 1
              INTZER = 0
              CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                     II, ISIZE, IHEAD, IUSE,
     $                     CP, RP, DN, N, WIR, WIC, WR, WC,
     $                     INTZER, INTZER, INTZER, INTZER, .FALSE.,
     $                     PFREE, XFREE, MHEAD, MTAIL, SLOTS)
              PC = CP (COL)
           ENDIF
           WXP = XHEAD
           XHEAD = XHEAD + XS
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (XHEAD .GT. XTAIL) THEN
              GO TO 9000
           ENDIF
           CDEG = 0
           CEP = (PC+9)
           CELN = II (PC+5)
           DO 190 IP = CEP, CEP + 2*CELN - 2, 2
              E = II (IP)
              F = II (IP+1)
              EP = RP (E)
              FDIMC = II (EP+1)
              FXP = II (EP+2)
              IF (E .LE. N) THEN
                 FLUIP = II (EP)
                 LUDEGC = II (FLUIP+3)
                 LUCP = (FLUIP + 7)
              ELSE
                 LUDEGC = FDIMC
                 LUCP = (EP+9)
              ENDIF
              XP = FXP + F * FDIMC
              CDEG1 = CDEG
              DO 160 P = LUCP, LUCP + LUDEGC - 1
                 ROW = II (P)
                 IF (ROW .GT. 0) THEN
                    IF (WIR (ROW) .LE. 0) THEN
                       CDEG = CDEG + 1
                       WM (CDEG) = ROW
                    ENDIF
                 ENDIF
160           CONTINUE
              DO 170 I = CDEG1+1, CDEG
                 ROW = WM (I)
                 WIR (ROW) = I
                 XX (WXP+I-1) = ZERO
170           CONTINUE
CFPP$ NODEPCHK L
              DO 180 J = 0, LUDEGC - 1
                 ROW = II (LUCP+J)
                 IF (ROW .GT. 0) THEN
                    XX (WXP + WIR (ROW) - 1) =
     $              XX (WXP + WIR (ROW) - 1) + XX (XP+J)
                 ENDIF
180           CONTINUE
190        CONTINUE
           CDEG1 = CDEG
           CLEN = II (PC+6)
           CSIZ = II (PC)
           IP = PC + CSIZ - CLEN
           CXP = II (PC+2)
CFPP$ NODEPCHK L
           DO 200 I = 0, CLEN - 1
              ROW = II (IP+I)
              WM (CDEG+1+I) = ROW
              XX (WXP+CDEG+I) = XX (CXP+I)
200        CONTINUE
           CDEG = CDEG + CLEN
           II (PC+1) = CDEG
           MAXVAL = ABS (XX (WXP - 1 + IDAMAX (CDEG, XX (WXP), 1)))
           RINFO (3) = RINFO (3) + CDEG
           TOLER = RELPT * MAXVAL
           RDEG = N+1
           IF (CDEG .NE. 0 .AND. MAXVAL .GT. ZERO) THEN
              IF (SYMSRC) THEN
                 ROW = COL
                 ROWPOS = WIR (ROW)
                 IF (ROWPOS .LE. 0) THEN
                    DO 210 I = CDEG1 + 1, CDEG1 + CLEN
                       IF (WM (I) .EQ. ROW) THEN
                          ROWPOS = I
                          GO TO 220
                       ENDIF
210                 CONTINUE
220                 CONTINUE
                 ENDIF
                 IF (ROWPOS .GT. 0) THEN
                    X = ABS (XX (WXP-1+ROWPOS))
                    IF (X .GE. TOLER .AND. X .GT. ZERO) THEN
                       PR = RP (ROW)
                       RDEG = II (PR+1)
                    ENDIF
                 ENDIF
              ENDIF
              IF (RDEG .EQ. N+1) THEN
                 ROW = N+1
                 DO 230 I = 1, CDEG
                    ROW2 = WM (I)
                    PR = RP (ROW2)
                    RDEG2 = II (PR+1)
                    BETTER = RDEG2 .LT. RDEG .OR.
     $                      (RDEG2 .EQ. RDEG .AND. ROW2 .LT. ROW)
                    IF (BETTER) THEN
                       X = ABS (XX (WXP-1+I))
                       IF (X.GE.TOLER .AND. X.GT.ZERO) THEN
                          ROW = ROW2
                          RDEG = RDEG2
                          ROWPOS = I
                       ENDIF
                    ENDIF
230              CONTINUE
              ENDIF
           ENDIF
           XHEAD = XHEAD - XS
           XUSE = XUSE - XS
           XNEED = XNEED - XS
           DO 240 I = 1, CDEG1
              WIR (WM (I)) = -1
240        CONTINUE
           IF (RDEG .EQ. N+1) THEN
              CDEG = -(N+2)
              II (PC+1) = CDEG
           ELSE
              SRCHED = SRCHED + 1
              II (PC+7) = SLIST
              SLIST = COL
              COST = (CDEG - 1) * (RDEG - 1)
              IF (PIVCOL .EQ. 0 .OR. COST .LT. BESTCO) THEN
                 FFLEFC = CDEG
                 DO 250 I = 1, FFLEFC-1
                    WPC (I) = WM (I)
250              CONTINUE
                 WPC (ROWPOS) = WM (FFLEFC)
                 PIVCOL = COL
                 PIVROW = ROW
                 BESTCO = COST
              ENDIF
           ENDIF
           IF (SRCHED .GE. NSRCH) THEN
              GO TO 260
           ENDIF
255     CONTINUE
260     CONTINUE
C=======================================================================
C=======================================================================
        IF (PIVCOL .EQ. 0) THEN
           K = N - NPIV + 1
           DO 270 COL = 1, N
              IF (CP (COL) .NE. 0) THEN
                 K = K - 1
                 WPC (K) = COL
                 CP (COL) = 0
              ENDIF
270        CONTINUE
           K = N - NPIV + 1
           DO 280 ROW = 1, NDN
              IF (ROW .GT. N) THEN
                 E = ROW
                 RP (E) = 0
              ELSE IF (RP (ROW) .NE. 0) THEN
                 RLEN = WC (ROW)
                 IF (RLEN .GE. 0 .AND. RLEN .LE. N) THEN
                    K = K - 1
                    WPR (K) = ROW
                    RP (ROW) = 0
                 ELSE IF (RLEN .NE. -(NDN+2)) THEN
                    E = ROW
                    EP = RP (ROW)
                    WR (E) = -(NDN+2)
                    WC (E) = -(NDN+2)
                    FLUIP = II (EP)
                    RP (E) = FLUIP
                 ENDIF
              ENDIF
280        CONTINUE
           GO TO 2000
        ENDIF
C=======================================================================
C=======================================================================
        DO 300 I = 1, SRCHED
           COL = SLIST
           PC = CP (COL)
           SLIST = II (PC+7)
           IF (COL .NE. PIVCOL) THEN
              CDEG = II (PC+1)
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) THEN
                 II (CP (CNXT)+8) = COL
              ENDIF
              HEAD (CDEG) = COL
              MINDEG = MIN (MINDEG, CDEG)
           ENDIF
300     CONTINUE
C=======================================================================
C=======================================================================
        PR = RP (PIVROW)
        FFLEFR = 0
        REP = (PR+2)
        RELN = WR (PIVROW)
        DO 330 IP = REP, REP + 2*RELN - 2, 2
           E = II (IP)
           EP = RP (E)
           IF (E .LE. N) THEN
              FLUIP = II (EP)
              LUCP = (FLUIP + 7)
              LUDEGR = II (FLUIP+2)
              LUDEGC = II (FLUIP+3)
              LURP = LUCP + LUDEGC
              F1 = FFLEFR
              DO 310 P = LURP, LURP + LUDEGR - 1
                 COL = II (P)
                 IF (COL .GT. 0) THEN
                    IF (WIC (COL) .EQ. -2) THEN
                       FFLEFR = FFLEFR + 1
                       WPR (FFLEFR) = COL
                    ENDIF
                 ENDIF
310           CONTINUE
              DO 320 I = F1+1, FFLEFR
                 WIC (WPR (I)) = 0
320           CONTINUE
           ELSE
              LURP = (EP+8)
              COL = II (LURP)
              IF (WIC (COL) .EQ. -2) THEN
                 FFLEFR = FFLEFR + 1
                 WPR (FFLEFR) = COL
                 WIC (COL) = 0
              ENDIF
           ENDIF
330     CONTINUE
        RSIZ = II (PR)
        RLEN = WC (PIVROW)
        DO 340 P = PR + RSIZ - RLEN, PR + RSIZ - 1
           COL = II (P)
           IF (WIC (COL) .EQ. -2) THEN
              FFLEFR = FFLEFR + 1
              WPR (FFLEFR) = COL
           ENDIF
340     CONTINUE
C=======================================================================
C=======================================================================
        FFROW = PIVROW
        E1 = PIVROW
        K = 1
        K0 = 0
        FFDIMR = MIN (KLEFT, INT (GRO * FFLEFR))
        FFDIMC = MIN (KLEFT, INT (GRO * FFLEFC))
        FMAXR = FFLEFR
        FMAXC = FFLEFC
        FFSIZE = FFDIMC * FFDIMR
        RSCAN = MAX (DSIZ, FFDIMR)
        DO 350 I = 1, FFLEFC - 1
           WIR (WPC (I)) = I - 1
350     CONTINUE
        DO 360 I = 1, FFLEFR
           WIC (WPR (I)) = (I - 1) * FFDIMC
360     CONTINUE
        COL = WPR (FFLEFR)
        COLPOS = (WIC (PIVCOL)/FFDIMC)+1
        WPR (COLPOS) = COL
        WIC (COL) = WIC (PIVCOL)
        WIC (PIVCOL) = (FFDIMR - 1) * FFDIMC
        WIR (PIVROW) = FFDIMC - 1
        FFLEFR = FFLEFR - 1
        FFLEFC = FFLEFC - 1
        IF (FFSIZE + FFDIMC .GT. XTAIL-XHEAD) THEN
           INFO (15) = INFO (15) + 1
           INTZER = 0
           CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                  II, ISIZE, IHEAD, IUSE,
     $                  CP, RP, DN, N, WIR, WIC, WR, WC,
     $                  INTZER, INTZER, INTZER, INTZER, .FALSE.,
     $                  PFREE, XFREE, MHEAD, MTAIL, SLOTS)
        ENDIF
        FFXP = XHEAD
        XHEAD = XHEAD + FFSIZE
        WXP = XHEAD
        XHEAD = XHEAD + FFDIMC
        XUSE = XUSE + FFSIZE + FFDIMC
        XNEED = XNEED + FFSIZE + FFDIMC
        INFO (20) = MAX (INFO (20), XUSE)
        INFO (21) = MAX (INFO (21), XNEED)
        IF (XHEAD .GT. XTAIL) THEN
           GO TO 9000
        ENDIF
        XRUSE = XRUSE + FFSIZE
        XRMAX = MAX (XRMAX, XRUSE)
        DO 380 J = 0, FFLEFR - 1
           DO 370 I = 0, FFLEFC - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
370        CONTINUE
380     CONTINUE
        DO 390 J = 0, FFLEFR - 1
           XX (FFXP + J*FFDIMC + FFDIMC-1) = ZERO
390     CONTINUE
        DO 400 I = 0, FFLEFC - 1
           XX (FFXP + (FFDIMR-1)*FFDIMC + I) = ZERO
400     CONTINUE
        XX (FFXP + (FFDIMR-1)*FFDIMC + FFDIMC-1) = ZERO
        DO 410 J = 1, FFLEFR
           PC = CP (WPR (J))
           CDEG = II (PC+1)
           IF (CDEG .GT. 0) THEN
              CNXT = II (PC+7)
              CPRV = II (PC+8)
              IF (CNXT .NE. 0) THEN
                 II (CP (CNXT)+8) = CPRV
              ENDIF
              IF (CPRV .NE. 0) THEN
                 II (CP (CPRV)+7) = CNXT
              ELSE
                 HEAD (CDEG) = CNXT
              ENDIF
           ENDIF
410     CONTINUE
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        SCAN1 = 0
        SCAN2 = 0
        SCAN3 = 0
        SCAN4 = 0
        DO 1395 DUMMY3 = 1, N
C=======================================================================
C=======================================================================
        KLEFT1 = KLEFT - 1
        ROW = PIVROW
        DO 440 J = SCAN1, FFLEFC
           IF (J .NE. 0) THEN
              ROW = WPC (J)
           ENDIF
           PR = RP (ROW)
           REP = (PR+2)
           RELN = WR (ROW)
CFPP$ NODEPCHK L
           DO 430 P = REP, REP + 2*RELN - 2, 2
              E = II (P)
              IF (WC (E) .LT. W0) THEN
                 EP = RP (E)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 WR (E) = FLEFTR + W0
                 WC (E) = FLEFTC + W0
              ENDIF
              WC (E) = WC (E) - 1
430        CONTINUE
440     CONTINUE
        COL = PIVCOL
        DO 460 J = SCAN2, FFLEFR
           IF (J .NE. 0) THEN
              COL = WPR (J)
           ENDIF
           PC = CP (COL)
           CELN = II (PC+5)
           CEP = (PC+9)
CFPP$ NODEPCHK L
           DO 450 P = CEP, CEP + 2*CELN - 2, 2
              E = II (P)
              IF (WR (E) .LT. W0) THEN
                 EP = RP (E)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 WR (E) = FLEFTR + W0
                 WC (E) = FLEFTC + W0
              ENDIF
              WR (E) = WR (E) - 1
450        CONTINUE
460     CONTINUE
        COL = PIVCOL
        DO 700 JJ = SCAN3, FFLEFR
           IF (JJ .NE. 0) THEN
              COL = WPR (JJ)
           ENDIF
           CDEG = 0
           DELN = 0
           PC = CP (COL)
           CEP = (PC+9)
           CELN = II (PC+5)
           IP2 = CEP + 2*CELN - 2
           XUDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
           DO 470 IP = CEP, IP2, 2
              E = II (IP)
              IF (WC (E) .GT. W0) THEN
                    CDEG = CDEG + (WC (E) - W0)
              ELSE
                 DELN = DELN + 1
                 WM (DELN) = IP
              ENDIF
470       CONTINUE
          IF (DELN .NE. 0) THEN
              P2 = IP2
              DO 480 I = DELN, 1, -1
                 E = II (WM (I)  )
                 F = II (WM (I)+1)
                 II (WM (I)  ) = II (P2  )
                 II (WM (I)+1) = II (P2+1)
                 II (P2  ) = E
                 II (P2+1) = F
                 P2 = P2 - 2
480           CONTINUE
              DO 670 IP = P2 + 2, IP2, 2
                 E = II (IP)
                 IF (WC (E) .LT. W0) THEN
                    GOTO 670
                 ENDIF
                 EP = RP (E)
                 FDIMC = II (EP+1)
                 FXP = II (EP+2)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 IF (E .LE. N) THEN
                    FLUIP = II (EP)
                    LUDEGR = II (FLUIP+2)
                    LUDEGC = II (FLUIP+3)
                    LUCP = (FLUIP + 7)
                    LURP = LUCP + LUDEGC
                    IF (WIR (E) .EQ. -1) THEN
                       WIR (E) = SONLST - N - 2
                       SONLST = E
                       NSONS = NSONS + 1
                    ENDIF
                 ELSE
                    LUDEGR = 1
                    LUDEGC = FDIMC
                    LUCP = (EP+9)
                    LURP = (EP+8)
                 ENDIF
                 IF (WR (E) .EQ. W0) THEN
                    IF (LUDEGC .EQ. FLEFTC) THEN
                       DO 490 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          WM (I+1) = WIR (ROW2)
490                    CONTINUE
                       IF (LUDEGR .EQ. FLEFTR) THEN
                          DO 510 J = 0, LUDEGR-1
                             COL2 = II (LURP+J)
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 500 I = 0, LUDEGC-1
                                XX (XDP + WM (I+1)) =
     $                          XX (XDP + WM (I+1)) +
     $                          XX (FXP + J*FDIMC + I)
500                          CONTINUE
510                       CONTINUE
                       ELSE
                          DO 530 J = 0, LUDEGR-1
                             COL2 = II (LURP+J)
                             IF (COL2 .GT. 0) THEN
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 520 I = 0, LUDEGC-1
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
520                             CONTINUE
                             ENDIF
530                       CONTINUE
                       ENDIF
                    ELSE
                       DEGC = 0
                       DO 540 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
                          ENDIF
540                    CONTINUE
                       IF (LUDEGR .EQ. FLEFTR) THEN
                          DO 560 J = 0, LUDEGR-1
                             COL2 = II (LURP+J)
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 550 I = 1, DEGC
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
550                          CONTINUE
560                       CONTINUE
                       ELSE
                          DO 580 J = 0, LUDEGR-1
                             COL2 = II (LURP+J)
                             IF (COL2 .GT. 0) THEN
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 570 I = 1, DEGC
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
570                             CONTINUE
                             ENDIF
580                       CONTINUE
                       ENDIF
                    ENDIF
                    WR (E) = -(NDN+2)
                    WC (E) = -(NDN+2)
                    IF (E .LE. N) THEN
                       RP (E) = FLUIP
                       II (EP) = FSCAL
                       INEED = INEED - FSCAL
                    ELSE
                       RP (E) = 0
                       II (EP) = FDIMC + CSCAL
                       INEED = INEED - (FDIMC + CSCAL)
                    ENDIF
                    II (EP+1) = -1
                    II (EP+6) = 0
                    MPREV = II (EP+4)
                    MNEXT = II (EP+3)
                    XNEED = XNEED - LUDEGR*LUDEGC
                    IF (MNEXT .NE. 0 .AND. II (MNEXT+6) .EQ. 0) THEN
                       MNEXT = II (MNEXT+3)
                       II (EP+3) = MNEXT
                       IF (MNEXT .NE. 0) THEN
                          II (MNEXT+4) = EP
                       ELSE
                          MTAIL = EP
                       ENDIF
                    ENDIF
                    IF (MPREV .NE. 0 .AND. II (MPREV+6) .EQ. 0) THEN
                       II (EP+2) = II (MPREV+2)
                       MPREV = II (MPREV+4)
                       II (EP+4) = MPREV
                       IF (MPREV .NE. 0) THEN
                          II (MPREV+3) = EP
                       ELSE
                          MHEAD = EP
                       ENDIF
                    ENDIF
                    IF (MNEXT .NE. 0) THEN
                       XS = II (MNEXT+2) - II (EP+2)
                    ELSE
                       XS = FFXP - II (EP+2)
                    ENDIF
                    IF (XS .GT. XFREE) THEN
                       XFREE = XS
                       PFREE = EP
                    ENDIF
                    XRUSE = XRUSE - LUDEGR*LUDEGC
                 ELSE IF (WR (E) - W0 .LE. FLEFTR/2) THEN
                    WC (E) = -USONS
                    USONS = E
                    IF (LUDEGC .EQ. FLEFTC) THEN
                       DO 590 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          WM (I+1) = WIR (ROW2)
590                    CONTINUE
                       DO 610 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN
                             IF (WIC (COL2) .GE. 0) THEN
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 600 I = 0, LUDEGC-1
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
600                             CONTINUE
                                II (LURP+J) = -COL2
                             ENDIF
                          ENDIF
610                    CONTINUE
                    ELSE
                       DEGC = 0
                       DO 620 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
                          ENDIF
620                    CONTINUE
                       DO 640 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN
                             IF (WIC (COL2) .GE. 0) THEN
                                XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                                DO 630 I = 1, DEGC
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
630                             CONTINUE
                                II (LURP+J) = -COL2
                             ENDIF
                          ENDIF
640                    CONTINUE
                    ENDIF
                    FLEFTR = WR (E) - W0
                    II (EP+5) = FLEFTR
                 ELSE
                    F = II (IP+1)
                    IF (LUDEGC .EQ. FLEFTC) THEN
CFPP$ NODEPCHK L
                       DO 650 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          XX (XUDP + WIR (ROW2)) =
     $                    XX (XUDP + WIR (ROW2)) +
     $                    XX (FXP + F*FDIMC + I)
650                    CONTINUE
                    ELSE
CFPP$ NODEPCHK L
                       DO 660 I = 0, LUDEGC-1
                          ROW2 = II (LUCP+I)
                          IF (ROW2 .GT. 0) THEN
                             XX (XUDP + WIR (ROW2)) =
     $                       XX (XUDP + WIR (ROW2)) +
     $                       XX (FXP + F*FDIMC + I)
                          ENDIF
660                    CONTINUE
                    ENDIF
                    II (EP+5) = FLEFTR - 1
                    II (LURP+F) = -COL
                 ENDIF
670           CONTINUE
              II (PC+5) = II (PC+5) - DELN
              INEED = INEED - 2*DELN
           ENDIF
           CLEN = II (PC+6)
           IF (CLEN .GT. 0) THEN
              CSIZ = II (PC)
              IP = PC + CSIZ - CLEN
              DLEN = 0
CFPP$ NODEPCHK L
              DO 680 I = 0, CLEN - 1
                 ROW = II (IP+I)
                 IF (WIR (ROW) .GE. 0) THEN
                    DLEN = DLEN + 1
                    WM (DLEN) = I
                 ENDIF
680           CONTINUE
              IF (DLEN .NE. 0) THEN
                 CXP = II (PC+2)
                 DO 690 J = 1, DLEN
                    I = WM (J)
                    ROW = II (IP+I)
                    XX (XUDP + WIR (ROW)) =
     $              XX (XUDP + WIR (ROW)) + XX (CXP+I)
                    II (IP +I) = II (IP +J-1)
                    XX (CXP+I) = XX (CXP+J-1)
690              CONTINUE
                 CLEN = CLEN - DLEN
                 CXP = CXP + DLEN
                 INEED = INEED - DLEN
                 XNEED = XNEED - DLEN
                 II (PC+6) = CLEN
                 IF (CLEN .NE. 0) THEN
                    II (PC+2) = CXP
                 ELSE
                    MPREV = II (PC+4)
                    MNEXT = II (PC+3)
                    IF (MNEXT .NE. 0 .AND. II (MNEXT+6) .EQ. 0) THEN
                       MNEXT = II (MNEXT+3)
                       II (PC+3) = MNEXT
                       IF (MNEXT .NE. 0) THEN
                          II (MNEXT+4) = PC
                       ELSE
                          MTAIL = PC
                       ENDIF
                    ENDIF
                    IF (MPREV .NE. 0 .AND. II (MPREV+6) .EQ. 0) THEN
                       II (PC+2) = II (MPREV+2)
                       MPREV = II (MPREV+4)
                       II (PC+4) = MPREV
                       IF (MPREV .NE. 0) THEN
                          II (MPREV+3) = PC
                       ELSE
                          MHEAD = PC
                       ENDIF
                    ENDIF
                    IF (PC .EQ. MHEAD) THEN
                       II (PC+2) = 1
                    ENDIF
                    IF (MNEXT .NE. 0) THEN
                       XS = II (MNEXT+2) - II (PC+2)
                    ELSE
                       XS = FFXP - II (PC+2)
                    ENDIF
                    IF (XS .GT. XFREE) THEN
                       XFREE = XS
                       PFREE = PC
                    ENDIF
                 ENDIF
              ENDIF
              CDEG = CDEG + CLEN
           ENDIF
           CDEG2 = II (PC+1)
           CDEG = MIN (KLEFT1 - FFLEFC, CDEG2, CDEG)
           II (PC+1) = CDEG
700     CONTINUE
710     CONTINUE
        IF (USONS .NE. NDN+1) THEN
           NEXT = -WC (USONS)
           WC (USONS) = W0
           USONS = NEXT
        GOTO 710
        ENDIF
        ROW = PIVROW
        DO 840 JJ = SCAN4, FFLEFC
           IF (JJ .NE. 0) THEN
              ROW = WPC (JJ)
           ENDIF
           RDEG = 0
           DELN = 0
           PR = RP (ROW)
           REP = (PR+2)
           RELN = WR (ROW)
           IP2 = REP + 2*RELN - 2
CFPP$ NODEPCHK L
           DO 720 IP = REP, IP2, 2
              E = II (IP)
              IF (WR (E) .GT. W0) THEN
                 RDEG = RDEG + (WR (E) - W0)
              ELSE
                 DELN = DELN + 1
                 WM (DELN) = IP
              ENDIF
720        CONTINUE
           IF (DELN .NE. 0) THEN
              P2 = IP2
              DO 730 I = DELN, 1, -1
                 E = II (WM (I)  )
                 F = II (WM (I)+1)
                 II (WM (I)  ) = II (P2  )
                 II (WM (I)+1) = II (P2+1)
                 II (P2  ) = E
                 II (P2+1) = F
                 P2 = P2 - 2
730           CONTINUE
              DO 810 IP = P2 + 2, IP2, 2
                 E = II (IP)
                 IF (WR (E) .LT. W0) THEN
                    GOTO 810
                 ENDIF
                 EP = RP (E)
                 FDIMC = II (EP+1)
                 FXP = II (EP+2)
                 FLEFTR = II (EP+5)
                 FLEFTC = II (EP+6)
                 IF (E .LE. N) THEN
                    FLUIP = II (EP)
                    LUDEGR = II (FLUIP+2)
                    LUDEGC = II (FLUIP+3)
                    LUCP = (FLUIP + 7)
                    LURP = LUCP + LUDEGC
                    IF (WIR (E) .EQ. -1) THEN
                       WIR (E) = SONLST - N - 2
                       SONLST = E
                       NSONS = NSONS + 1
                    ENDIF
                 ELSE
                    LUDEGR = 1
                    LUDEGC = FDIMC
                    LUCP = (EP+9)
                    LURP = (EP+8)
                 ENDIF
                 IF (WC (E) - W0 .LE. FLEFTC/2) THEN
                    WR (E) = -LSONS
                    LSONS = E
                    DEGC = 0
                    DO 740 I = 0, LUDEGC-1
                       ROW2 = II (LUCP+I)
                       IF (ROW2 .GT. 0) THEN
                          IF (WIR (ROW2) .GE. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW2)
                             II (LUCP+I) = -ROW2
                          ENDIF
                       ENDIF
740                 CONTINUE
                    IF (LUDEGR .EQ. FLEFTR) THEN
                       DO 760 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                          DO 750 I = 1, DEGC
                             XX (XDP + WM (I)) =
     $                       XX (XDP + WM (I)) +
     $                       XX (FXP + J*FDIMC + WJ (I))
750                       CONTINUE
760                    CONTINUE
                    ELSE
                       DO 780 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN
                             XDP = FFXP + WIC (COL2)
CFPP$ NODEPCHK L
                             DO 770 I = 1, DEGC
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
770                          CONTINUE
                          ENDIF
780                    CONTINUE
                    ENDIF
                    FLEFTC = WC (E) - W0
                    II (EP+6) = FLEFTC
                 ELSE
                    XLDP = FFXP + WIR (ROW)
                    F = II (IP+1)
                    IF (LUDEGR .EQ. FLEFTR) THEN
CFPP$ NODEPCHK L
                       DO 790 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          XX (XLDP + WIC (COL2)) =
     $                    XX (XLDP + WIC (COL2)) +
     $                    XX (FXP + J*FDIMC + F)
790                    CONTINUE
                    ELSE
CFPP$ NODEPCHK L
                       DO 800 J = 0, LUDEGR-1
                          COL2 = II (LURP+J)
                          IF (COL2 .GT. 0) THEN
                             XX (XLDP + WIC (COL2)) =
     $                       XX (XLDP + WIC (COL2)) +
     $                       XX (FXP + J*FDIMC + F)
                          ENDIF
800                    CONTINUE
                    ENDIF
                    II (EP+6) = FLEFTC - 1
                    II (LUCP+F) = -ROW
                 ENDIF
810           CONTINUE
              WR (ROW) = WR (ROW) - DELN
              INEED = INEED - 2*DELN
           ENDIF
           RLEN = WC (ROW)
           IF (RLEN .GT. 0) THEN
              IF (RLEN .LE. RSCAN) THEN
                 RSIZ = II (PR)
                 IP = PR + RSIZ - RLEN
                 DLEN = 0
CFPP$ NODEPCHK L
                 DO 820 P = IP, IP + RLEN - 1
                    COL = II (P)
                    IF (WIC (COL) .NE. -2) THEN
                       DLEN = DLEN + 1
                       WM (DLEN) = P
                    ENDIF
820              CONTINUE
                 IF (DLEN .NE. 0) THEN
                    DO 830 J = 1, DLEN
                       II (WM (J)) = II (IP+J-1)
830                 CONTINUE
                    RLEN = RLEN - DLEN
                    INEED = INEED - DLEN
                    WC (ROW) = RLEN
                 ENDIF
              ENDIF
              RDEG = RDEG + RLEN
           ENDIF
           RDEG2 = II (PR+1)
           RDEG = MIN (KLEFT1 - FFLEFR, RDEG2, RDEG)
           II (PR+1) = RDEG
840     CONTINUE
850     CONTINUE
        IF (LSONS .NE. NDN+1) THEN
           NEXT = -WR (LSONS)
           WR (LSONS) = W0
           LSONS = NEXT
        GOTO 850
        ENDIF
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        DO 1324 DUMMY4 = 1, N
C=======================================================================
C=======================================================================
        XDP = FFXP + (FFDIMR - K) * FFDIMC
        X = XX (XDP + FFDIMC - K)
        X = ONE / X
        DO 870 P = XDP, XDP + FFLEFC-1
           XX (P) = XX (P) * X
870     CONTINUE
        RINFO (4) = RINFO (4) + FFLEFC
C=======================================================================
C=======================================================================
        KLEFT = KLEFT - 1
        NPIV = NPIV + 1
        INFO (17) = INFO (17) + 1
        WPR (N-NPIV+1) = PIVROW
        WPC (N-NPIV+1) = PIVCOL
        WIR (PIVROW) = -1
        WIC (PIVCOL) = -1
        RLEN = WC (PIVROW)
        INEED = INEED - CSCAL - RSCAL - RLEN
        PR = RP (PIVROW)
        PC = CP (PIVCOL)
        II (PR+1) = -1
        II (PC+1) = -1
        RP (PIVROW) = 0
        CP (PIVCOL) = 0
C=======================================================================
C=======================================================================
        FEDEGC = FFLEFC
        FEDEGR = FFLEFR
        PFOUND = .FALSE.
        OKCOL = FFLEFC .GT. 0
        OKROW = .FALSE.
        IF (OKCOL) THEN
           COLPOS = 0
           PIVCOL = N+1
           CDEG = N+1
           DO 880 J = 1, FFLEFR
              COL = WPR (J)
              PC = CP (COL)
              CDEG2 = II (PC+1)
              BETTER = CDEG2 .GE. 0 .AND.
     $                (CDEG2 .LT. CDEG .OR.
     $                (CDEG2 .EQ. CDEG .AND. COL .LT. PIVCOL))
              IF (BETTER) THEN
                 CDEG = CDEG2
                 COLPOS = J
                 PIVCOL = COL
              ENDIF
880        CONTINUE
           OKCOL = COLPOS .NE. 0
        ENDIF
C=======================================================================
C=======================================================================
        IF (OKCOL) THEN
           PC = CP (PIVCOL)
           CLEN = II (PC+6)
           OKCOL = FEDEGC + CLEN .LE. FFDIMC
        ENDIF
        IF (OKCOL) THEN
           P = FFXP + (COLPOS - 1) * FFDIMC - 1
CFPP$ NODEPCHK L
           DO 890 I = 1, FFLEFC
              XX (WXP-1+I) = XX (P+I)
890        CONTINUE
           IF (K-K0 .GT. 0 .AND. FFLEFC .NE. 0) THEN
              CALL DGEMV ('N', FFLEFC, K-K0,
     $          -ONE, XX (FFXP + (FFDIMR - K) * FFDIMC)        ,FFDIMC,
     $                XX (FFXP + (COLPOS - 1) * FFDIMC + FFDIMC - K), 1,
     $           ONE, XX (WXP)                                      , 1)
              RINFO (3) = RINFO (3) + 2*FFLEFC*(K-K0)
           ENDIF
           CEP = (PC+9)
           CELN = II (PC+5)
           DO 930 IP = CEP, CEP + 2*CELN - 2, 2
              E = II (IP)
              F = II (IP+1)
              EP = RP (E)
              FLEFTC = II (EP+6)
              FDIMC = II (EP+1)
              FXP = II (EP+2)
              IF (E .LE. N) THEN
                 FLUIP = II (EP)
                 LUCP = (FLUIP + 7)
                 LUDEGC = II (FLUIP+3)
              ELSE
                 LUCP = (EP+9)
                 LUDEGC = FDIMC
              ENDIF
              XP = FXP + F * FDIMC
              F1 = FEDEGC
              DO 900 P = LUCP, LUCP + LUDEGC - 1
                 ROW = II (P)
                 IF (ROW .GT. 0) THEN
                    IF (WIR (ROW) .LT. 0) THEN
                       F1 = F1 + 1
                       WPC (F1) = ROW
                    ENDIF
                 ENDIF
900           CONTINUE
              OKCOL = F1 + CLEN .LE. FFDIMC
              IF (.NOT. OKCOL) THEN
                 GO TO 940
              ENDIF
              DO 910 I = FEDEGC+1, F1
                 ROW = WPC (I)
                 WIR (ROW) = I - 1
                 XX (WXP-1+I) = ZERO
910           CONTINUE
              FEDEGC = F1
CFPP$ NODEPCHK L
              DO 920 J = 0, LUDEGC - 1
                 ROW = II (LUCP+J)
                 IF (ROW .GT. 0) THEN
                    XX (WXP + WIR (ROW)) =
     $              XX (WXP + WIR (ROW)) + XX (XP+J)
                 ENDIF
920           CONTINUE
930        CONTINUE
940        CONTINUE
        ENDIF
C=======================================================================
C=======================================================================
        IF (OKCOL) THEN
           CSIZ = II (PC)
           IP = PC + CSIZ - CLEN
           CXP = II (PC+2)
CFPP$ NODEPCHK L
           DO 950 I = 0, CLEN - 1
              ROW = II (IP+I)
              WIR (ROW) = FEDEGC + I
              WPC (FEDEGC+1+I) = ROW
              XX  (WXP+FEDEGC+I) = XX (CXP+I)
950        CONTINUE
           FEDEGC = FEDEGC + CLEN
           CDEG = FEDEGC - FFLEFC
           II (PC+1) = CDEG
           MAXVAL = ABS (XX (WXP-1 + IDAMAX (FEDEGC, XX (WXP), 1)))
           RINFO (3) = RINFO (3) + FEDEGC
           TOLER = RELPT * MAXVAL
           RDEG = N+1
           IF (MAXVAL .GT. ZERO) THEN
              IF (SYMSRC) THEN
                 PIVROW = PIVCOL
                 ROWPOS = WIR (PIVROW) + 1
                 IF (ROWPOS .GT. 0 .AND. ROWPOS .LE. FFLEFC) THEN
                    X = ABS (XX (WXP-1+ROWPOS))
                    IF (X.GE.TOLER .AND. X.GT.ZERO) THEN
                       PR = RP (PIVROW)
                       RDEG = II (PR+1)
                    ENDIF
                 ENDIF
              ENDIF
              IF (RDEG .EQ. N+1) THEN
                 PIVROW = N+1
                 DO 960 I = 1, FFLEFC
                    ROW2 = WPC (I)
                    PR = RP (ROW2)
                    RDEG2 = II (PR+1)
                    BETTER = RDEG2 .LT. RDEG .OR.
     $                      (RDEG2 .EQ. RDEG .AND. ROW2 .LT. PIVROW)
                    IF (BETTER) THEN
                       X = ABS (XX (WXP-1+I))
                       IF (X.GE.TOLER .AND. X.GT.ZERO) THEN
                          PIVROW = ROW2
                          RDEG = RDEG2
                          ROWPOS = I
                       ENDIF
                    ENDIF
960              CONTINUE
              ENDIF
           ELSE
              CDEG = -(N+2)
              II (PC+1) = CDEG
           ENDIF
           OKROW = RDEG .NE. N+1
        ENDIF
C=======================================================================
C=======================================================================
        IF (OKROW) THEN
           PR = RP (PIVROW)
           REP = (PR+2)
           RELN = WR (PIVROW)
           DO 990 IP = REP, REP + 2*RELN - 2, 2
              E = II (IP)
              EP = RP (E)
              IF (E .LE. N) THEN
                 FLUIP = II (EP)
                 LUCP = (FLUIP + 7)
                 LUDEGR = II (FLUIP+2)
                 LUDEGC = II (FLUIP+3)
                 LURP = LUCP + LUDEGC
                 FLEFTR = II (EP+5)
                 OKROW = FLEFTR .LE. FFDIMR
                 IF (.NOT. OKROW) THEN
                    GO TO 1000
                 ENDIF
                 F1 = FEDEGR
                 DO 970 P = LURP, LURP + LUDEGR - 1
                    COL = II (P)
                    IF (COL .GT. 0) THEN
                       IF (WIC (COL) .EQ. -2) THEN
                          F1 = F1 + 1
                          WPR (F1) = COL
                       ENDIF
                    ENDIF
970              CONTINUE
                 OKROW = F1 .LE. FFDIMR
                 IF (.NOT. OKROW) THEN
                    GO TO 1000
                 ENDIF
                 DO 980 I = FEDEGR+1, F1
                    WIC (WPR (I)) = (I - 1) * FFDIMC
980              CONTINUE
                 FEDEGR = F1
              ELSE
                 LURP = (EP+8)
                 COL = II (LURP)
                 IF (WIC (COL) .EQ. -2) THEN
                    WIC (COL) = FEDEGR * FFDIMC
                    FEDEGR = FEDEGR + 1
                    WPR (FEDEGR) = COL
                    OKROW = FEDEGR .LE. FFDIMR
                    IF (.NOT. OKROW) THEN
                       GO TO 1000
                    ENDIF
                 ENDIF
              ENDIF
990        CONTINUE
1000       CONTINUE
        ENDIF
        IF (OKROW) THEN
           RLEN = WC (PIVROW)
           IF (RLEN .GT. 0) THEN
              F1 = FEDEGR
              RSIZ = II (PR)
              P2 = PR + RSIZ
              DO 1010 P = P2 - RLEN, P2 - 1
                 COL = II (P)
                 IF (WIC (COL) .EQ. -2) THEN
                    F1 = F1 + 1
                    WPR (F1) = COL
                 ENDIF
1010          CONTINUE
              RLEN2 = F1 - FEDEGR
              IF (RLEN2 .LT. RLEN) THEN
                 DO 1020 I = FEDEGR+1, F1
                    II (P2 - F1 + I - 1) = WPR (I)
1020             CONTINUE
                 INEED = INEED - (RLEN - RLEN2)
                 WC (PIVROW) = RLEN2
              ENDIF
              RDEG = F1 - FFLEFR
              II (PR+1) = RDEG
              OKROW = F1 .LE. FFDIMR
              IF (OKROW) THEN
                 DO 1030 I = FEDEGR+1, F1
                    WIC (WPR (I)) = (I - 1) * FFDIMC
1030             CONTINUE
                 FEDEGR = F1
              ENDIF
           ELSE
              RDEG = FEDEGR - FFLEFR
              II (PR+1) = RDEG
           ENDIF
        ENDIF
        PFOUND = OKROW .AND. OKCOL
        IF (.NOT. PFOUND) THEN
           MOVELU = K .GT. 0
           DO 1040 I = FFLEFR+1, FEDEGR
              WIC (WPR (I)) = -2
1040       CONTINUE
           FEDEGR = FFLEFR
           DO 1050 I = FFLEFC+1, FEDEGC
              WIR (WPC (I)) = -1
1050       CONTINUE
           FEDEGC = FFLEFC
        ELSE
           MOVELU = FEDEGC .GT. FFDIMC - K .OR. FEDEGR .GT. FFDIMR - K
        ENDIF
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        IF (K-K0 .GE. NB .OR. MOVELU) THEN
           CALL DGEMM ('N', 'N', FFLEFC, FFLEFR, K-K0,
     $          -ONE, XX (FFXP + (FFDIMR - K) * FFDIMC), FFDIMC,
     $                XX (FFXP +  FFDIMC - K)          , FFDIMC,
     $           ONE, XX (FFXP)                        , FFDIMC)
           RINFO (6) = RINFO (6) + 2*FFLEFC*FFLEFR*(K-K0)
           K0 = K
        ENDIF
C=======================================================================
C=======================================================================
        IF (MOVELU) THEN
           LUDEGR = FFLEFR
           LUDEGC = FFLEFC
           XS = K*LUDEGC + K*LUDEGR + K*K
           IS = 7 + LUDEGC + LUDEGR + NSONS
           IF (IS .GT. ITAIL-IHEAD .OR. XS .GT. XTAIL-XHEAD) THEN
              IF (IS .GT. ITAIL-IHEAD) THEN
                 INFO (14) = INFO (14) + 1
              ENDIF
              IF (XS .GT. XTAIL-XHEAD) THEN
                 INFO (15) = INFO (15) + 1
              ENDIF
              CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                     II, ISIZE, IHEAD, IUSE,
     $                     CP, RP, DN, N, WIR, WIC, WR, WC,
     $                     FFXP, FFSIZE, WXP, FFDIMC, .FALSE.,
     $                     PFREE, XFREE, MHEAD, MTAIL, SLOTS)
           ENDIF
           ITAIL = ITAIL - IS
           LUIP = ITAIL
           IUSE = IUSE + IS
           INEED = INEED + IS
           XTAIL = XTAIL - XS
           LUXP = XTAIL
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (19) = MAX (INFO (19), INEED)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
              GO TO 9000
           ENDIF
           XRUSE = XRUSE + XS
           XRMAX = MAX (XRMAX, XRUSE)
           II (LUIP) = LUXP
           II (LUIP+1) = K
           II (LUIP+2) = LUDEGR
           II (LUIP+3) = LUDEGC
           II (LUIP+4) = NSONS
           II (LUIP+5) = 0
           II (LUIP+6) = 0
           E = FFROW
           IF (E .EQ. E1) THEN
              LUIP1 = LUIP
           ENDIF
           WR (E) = -(NDN+2)
           WC (E) = -(NDN+2)
           LUCP = (LUIP + 7)
           DO 1060 I = 0, LUDEGC-1
              II (LUCP+I) = WPC (I+1)
1060       CONTINUE
           LURP = LUCP + LUDEGC
           DO 1070 I = 0, LUDEGR-1
              II (LURP+I) = WPR (I+1)
1070       CONTINUE
           LUSONP = LURP + LUDEGR
           IP = LUSONP
           E = SONLST
1080       CONTINUE
           IF (E .GT. 0) THEN
              EP = RP (E)
              IF (WC (E) .EQ. -(NDN+2)) THEN
                 II (IP) = E
              ELSE IF (WC (E) .EQ. W0) THEN
                 II (IP) = E + N
              ELSE IF (WR (E) .EQ. W0) THEN
                 II (IP) = E + 2*N
              ENDIF
              NEXT = WIR (E) + N + 2
              WIR (E) = -1
              E = NEXT
              IP = IP + 1
           GOTO 1080
           ENDIF
           NSONS = 0
           SONLST = 0
           LDIMC = K + LUDEGC
           XP = FFXP + (FFDIMR-1)*FFDIMC + FFDIMC-1
           DO 1100 J = 0, K-1
CFPP$ NODEPCHK L
              DO 1090 I = 0, K-1
                 XX (LUXP + J*LDIMC + I) = XX (XP - J*FFDIMC - I)
1090          CONTINUE
1100       CONTINUE
           IF (LUDEGC .NE. 0) THEN
              LXP = LUXP + K
              XP = FFXP + (FFDIMR-1)*FFDIMC
              DO 1120 J = 0, K-1
CFPP$ NODEPCHK L
                 DO 1110 I = 0, LUDEGC-1
                    XX (LXP + J*LDIMC + I) = XX (XP - J*FFDIMC + I)
1110             CONTINUE
1120          CONTINUE
           ENDIF
           IF (LUDEGR .NE. 0) THEN
              UXP = LUXP + K * LDIMC
              XP = FFXP + FFDIMC-1
              DO 1140 J = 0, LUDEGR-1
CFPP$ NODEPCHK L
                 DO 1130 I = 0, K-1
                    XX (UXP + J*K + I) = XX (XP + J*FFDIMC - I)
1130             CONTINUE
1140          CONTINUE
           ENDIF
           NLU = NLU + 1
           NZU = (K*(K-1)/2) + K*LUDEGC
           NZL = (K*(K-1)/2) + K*LUDEGR
           INFO (10) = INFO (10) + NZL
           INFO (11) = INFO (11) + NZU
           K = 0
           K0 = 0
           IF (PFOUND) THEN
              NSONS = 1
              E = FFROW
              WIR (E) = - N - 2
              SONLST = E
              RP (E) = LUIP
              FFROW = PIVROW
           ENDIF
        ENDIF
C=======================================================================
C=======================================================================
        IF (.NOT. PFOUND) THEN
           GO TO 1400
        ENDIF
C=======================================================================
C=======================================================================
        XSP = (COLPOS - 1) * FFDIMC
        XDP = (FFDIMR - K - 1) * FFDIMC
        FSP = FFXP + XSP
        FDP = FFXP + XDP
        IF (K-K0 .GT. 0 .AND. FFLEFC .NE. 0) THEN
           CALL DGEMV ('N', FFLEFC, K-K0,
     $          -ONE, XX (FDP + FFDIMC    ), FFDIMC,
     $                XX (FSP + FFDIMC - K), 1,
     $           ONE, XX (FSP             ), 1)
           RINFO (5) = RINFO (5) + 2*FFLEFC*(K-K0)
        ENDIF
        IF (FFLEFR .LT. FFDIMR - K) THEN
           XLP = (FFLEFR - 1) * FFDIMC
           IF (FFLEFR .EQ. COLPOS) THEN
CFPP$ NODEPCHK L
              DO 1160 I = 0, FFLEFC - 1
                 XX (FDP+I) = XX (FSP+I)
1160          CONTINUE
CFPP$ NODEPCHK L
              DO 1170 I = FFDIMC - K, FFDIMC - 1
                 XX (FDP+I) = XX (FSP+I)
1170          CONTINUE
           ELSE
              FLP = FFXP + XLP
CFPP$ NODEPCHK L
              DO 1190 I = 0, FFLEFC - 1
                 XX (FDP+I) = XX (FSP+I)
                 XX (FSP+I) = XX (FLP+I)
1190          CONTINUE
CFPP$ NODEPCHK L
              DO 1200 I = FFDIMC - K, FFDIMC - 1
                 XX (FDP+I) = XX (FSP+I)
                 XX (FSP+I) = XX (FLP+I)
1200          CONTINUE
              SWPCOL = WPR (FFLEFR)
              WPR (COLPOS) = SWPCOL
              WIC (SWPCOL) = XSP
           ENDIF
           IF (FEDEGR .NE. FFLEFR) THEN
              SWPCOL = WPR (FEDEGR)
              WPR (FFLEFR) = SWPCOL
              WIC (SWPCOL) = XLP
           ENDIF
        ELSE IF (COLPOS .NE. FFDIMR - K) THEN
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
           DO 1220 I = 0, FFLEFC - 1
              X = XX (FDP+I)
              XX (FDP+I) = XX (FSP+I)
              XX (FSP+I) = X
1220       CONTINUE
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
           DO 1230 I = FFDIMC - K, FFDIMC - 1
              X = XX (FDP+I)
              XX (FDP+I) = XX (FSP+I)
              XX (FSP+I) = X
1230       CONTINUE
           SWPCOL = WPR (FFDIMR - K)
           WPR (COLPOS) = SWPCOL
           WIC (SWPCOL) = XSP
        ENDIF
        WIC (PIVCOL) = XDP
        FEDEGR = FEDEGR - 1
        SCAN2 = FFLEFR
        FFLEFR = FFLEFR - 1
C=======================================================================
C=======================================================================
        XSP = ROWPOS - 1
        XDP = FFDIMC - K - 1
        FSP = FFXP + XSP
        FDP = FFXP + XDP
        IF (FFLEFC .LT. FFDIMC - K) THEN
           XLP = FFLEFC - 1
           IF (FFLEFC .EQ. ROWPOS) THEN
CFPP$ NODEPCHK L
              DO 1250 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC
                 XX (FDP+J) = XX (FSP+J)
1250          CONTINUE
CFPP$ NODEPCHK L
              DO 1260 J = (FFDIMR - K - 1) * FFDIMC,
     $                    (FFDIMR - 1) * FFDIMC, FFDIMC
                 XX (FDP+J) = XX (FSP+J)
1260          CONTINUE
           ELSE
              FLP = FFXP + XLP
CFPP$ NODEPCHK L
              DO 1280 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC
                 XX (FDP+J) = XX (FSP+J)
                 XX (FSP+J) = XX (FLP+J)
1280          CONTINUE
CFPP$ NODEPCHK L
              DO 1290 J = (FFDIMR - K - 1) * FFDIMC,
     $                    (FFDIMR - 1) * FFDIMC, FFDIMC
                 XX (FDP+J) = XX (FSP+J)
                 XX (FSP+J) = XX (FLP+J)
1290          CONTINUE
              SWPROW = WPC (FFLEFC)
              WPC (ROWPOS) = SWPROW
              WIR (SWPROW) = XSP
           ENDIF
           IF (FEDEGC .NE. FFLEFC) THEN
              SWPROW = WPC (FEDEGC)
              WPC (FFLEFC) = SWPROW
              WIR (SWPROW) = XLP
           ENDIF
        ELSE IF (ROWPOS .NE. FFDIMC - K) THEN
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
           DO 1310 J = 0, (FFLEFR - 1) * FFDIMC, FFDIMC
              X = XX (FDP+J)
              XX (FDP+J) = XX (FSP+J)
              XX (FSP+J) = X
1310       CONTINUE
CFPP$ NODEPCHK L
CFPP$ NOLSTVAL L
           DO 1320 J = (FFDIMR - K - 1) * FFDIMC,
     $                 (FFDIMR - 1) * FFDIMC, FFDIMC
              X = XX (FDP+J)
           XX (FDP+J) = XX (FSP+J)
                 XX (FSP+J) = X
1320       CONTINUE
           SWPROW = WPC (FFDIMC - K)
           WPC (ROWPOS) = SWPROW
           WIR (SWPROW) = XSP
        ENDIF
        WIR (PIVROW) = XDP
        FEDEGC = FEDEGC - 1
        SCAN1 = FFLEFC
        FFLEFC = FFLEFC - 1
        IF (K-K0 .GT. 0 .AND. FFLEFR .GT. 0) THEN
           CALL DGEMV ('T', K-K0, FFLEFR,
     $       -ONE, XX (FDP + 1)                    , FFDIMC,
     $             XX (FDP + (FFDIMR - K) * FFDIMC), FFDIMC,
     $        ONE, XX (FDP)                        , FFDIMC)
           RINFO (5) = RINFO (5) + 2*(K-K0)*FFLEFR
        ENDIF
C=======================================================================
C=======================================================================
        IF (FEDEGC .EQ. FFLEFC) THEN
           SCAN3 = FFLEFR + 1
        ELSE
           SCAN3 = 0
        ENDIF
        IF (FEDEGR .EQ. FFLEFR) THEN
           SCAN4 = FFLEFC + 1
        ELSE
           SCAN4 = 0
        ENDIF
C=======================================================================
C=======================================================================
        K = K + 1
        IF (FEDEGR .NE. FFLEFR .OR. FEDEGC .NE. FFLEFC) THEN
           GO TO 1325
        ENDIF
1324    CONTINUE
1325    CONTINUE
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        DO 1350 J = FFLEFR, FEDEGR - 1
           DO 1330 I = 0, FEDEGC - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
1330       CONTINUE
           DO 1340 I = FFDIMC - K, FFDIMC - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
1340       CONTINUE
1350    CONTINUE
CFPP$ NODEPCHK L
        DO 1380 I = FFLEFC, FEDEGC - 1
CFPP$ NODEPCHK L
           DO 1360 J = 0, FFLEFR - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
1360       CONTINUE
CFPP$ NODEPCHK L
           DO 1370 J = FFDIMR - K, FFDIMR - 1
              XX (FFXP + J*FFDIMC + I) = ZERO
1370       CONTINUE
1380    CONTINUE
        DO 1390 J = FFLEFR+1, FEDEGR
           PC = CP (WPR (J))
           CDEG = II (PC+1)
           IF (CDEG .GT. 0) THEN
              CNXT = II (PC+7)
              CPRV = II (PC+8)
              IF (CNXT .NE. 0) THEN
                 II (CP (CNXT)+8) = CPRV
              ENDIF
              IF (CPRV .NE. 0) THEN
                 II (CP (CPRV)+7) = CNXT
              ELSE
                 HEAD (CDEG) = CNXT
              ENDIF
           ENDIF
1390    CONTINUE
        FFLEFC = FEDEGC
        FFLEFR = FEDEGR
        FMAXR = MAX (FMAXR, FFLEFR + K)
        FMAXC = MAX (FMAXC, FFLEFC + K)
C=======================================================================
C=======================================================================
1395    CONTINUE
1400    CONTINUE
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        II (LUIP1+5) = FMAXR
        II (LUIP1+6) = FMAXC
        INFO (13) = INFO (13) + 1
        DO 1410 J = FFLEFR, 1, -1
           COL = WPR (J)
           PC = CP (COL)
           CDEG = II (PC+1)
           CDEG = MIN (KLEFT, CDEG + FFLEFC)
           IF (CDEG .GT. 0) THEN
              II (PC+1) = CDEG
              CNXT = HEAD (CDEG)
              II (PC+7) = CNXT
              II (PC+8) = 0
              IF (CNXT .NE. 0) THEN
                 II (CP (CNXT)+8) = COL
              ENDIF
              HEAD (CDEG) = COL
              MINDEG = MIN (MINDEG, CDEG)
           ENDIF
1410    CONTINUE
CFPP$ NODEPCHK L
        DO 1420 I = 1, FFLEFC
           ROW = WPC (I)
           PR = RP (ROW)
           RDEG = II (PR+1)
           RDEG = MIN (KLEFT, RDEG + FFLEFR)
           II (PR+1) = RDEG
1420    CONTINUE
        W1 = W0 + FMAX + 1
        IF (W1 .LE. W0) THEN
           W0 = NDN+2
           DO 1430 E = 1, N+DN
              IF (WR (E) .GT. NDN) THEN
                 WR (E) = W0-1
                 WC (E) = W0-1
              ENDIF
1430       CONTINUE
        ELSE
           W0 = W1
        ENDIF
        XUSE = XUSE - FFDIMC
        XNEED = XNEED - FFDIMC
        XHEAD = XHEAD - FFDIMC
        E = FFROW
        XS = FFLEFR * FFLEFC
        FMAX = MAX (FMAX, FFLEFR, FFLEFC)
        XRUSE = XRUSE - FFSIZE + XS
        IF (FFLEFR .LE. 0 .OR. FFLEFC .LE. 0) THEN
           RP (E) = LUIP
           XUSE = XUSE - FFSIZE
           XNEED = XNEED - FFSIZE
           XHEAD = FFXP
           DO 1440 I = 1, FFLEFR
              WIC (WPR (I)) = -2
1440       CONTINUE
           DO 1450 I = 1, FFLEFC
              WIR (WPC (I)) = -1
1450       CONTINUE
           GOTO 1540
        ENDIF
        IF (FSCAL .GT. ITAIL-IHEAD) THEN
           INFO (14) = INFO (14) + 1
           INTZER = 0
           CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                  II, ISIZE, IHEAD, IUSE,
     $                  CP, RP, DN, N, WIR, WIC, WR, WC,
     $                  FFXP, FFSIZE, INTZER, INTZER, .FALSE.,
     $                  PFREE, XFREE, MHEAD, MTAIL, SLOTS)
        ENDIF
        EP = IHEAD
        IHEAD = IHEAD + FSCAL
        IUSE = IUSE + FSCAL
        INEED = INEED + FSCAL
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .GT. ITAIL) THEN
           GO TO 9000
        ENDIF
        RP (E) = EP
        II (EP) = LUIP
        II (EP+5) = FFLEFR
        II (EP+6) = FFLEFC
        WR (E) = W0-1
        WC (E) = W0-1
        RINFO (2) = RINFO (2) + XS
        IF (XS .LE. XFREE) THEN
           XDP = II (PFREE+2)
           II (PFREE+2) = II (PFREE+2) + XS
           XFREE = XFREE - XS
           MPREV = II (PFREE+4)
           IF (XFREE .EQ. 0) THEN
              MNEXT = II (PFREE+3)
              PFREE = 0
              XFREE = -1
           ELSE
              MNEXT = PFREE
           ENDIF
           IF (MNEXT .NE. 0) THEN
              II (MNEXT+4) = EP
           ELSE
              MTAIL = EP
           ENDIF
           IF (MPREV .NE. 0) THEN
              II (MPREV+3) = EP
           ELSE
              MHEAD = EP
           ENDIF
           DO 1470 J = 0, FFLEFR - 1
CFPP$ NODEPCHK L
              DO 1460 I = 0, FFLEFC - 1
                 XX (XDP + J*FFLEFC + I) = XX (FFXP + J*FFDIMC + I)
1460          CONTINUE
1470       CONTINUE
           XHEAD = FFXP
           XUSE = XUSE - FFSIZE
           XNEED = XNEED - FFSIZE + XS
           FFDIMC = FFLEFC
           II (EP+1) = FFDIMC
           II (EP+2) = XDP
           II (EP+3) = MNEXT
           II (EP+4) = MPREV
        ELSE
           XNEED = XNEED - FFSIZE + XS
           XS = FFSIZE - (FFLEFC + (FFLEFR-1)*FFDIMC)
           XHEAD = XHEAD - XS
           XUSE = XUSE - XS
           II (EP+1) = FFDIMC
           II (EP+2) = FFXP
           II (EP+3) = 0
           II (EP+4) = MTAIL
           IF (MTAIL .EQ. 0) THEN
              MHEAD = EP
           ELSE
              II (MTAIL+3) = EP
           ENDIF
           MTAIL = EP
        ENDIF
        INEED = INEED + 2*(FFLEFR+FFLEFC)
        IWORST = INEED + LIMIT + CSCAL
        INFO (19) = MAX (INFO (19), IWORST)
        INFO (18) = MAX (INFO (18), IWORST)
        DO 1500 I = 1, FFLEFR
           COL = WPR (I)
           PC = CP (COL)
           CELN = II (PC+5)
           CSIZ = II (PC)
           CLEN = II (PC+6)
           WIC (COL) = -2
           IF (2*(CELN+1) + CLEN + CSCAL .GT. CSIZ) THEN
              IS = 2 * (CELN + 1) + CLEN
              IS = MIN (IS + MAX (16, IS), LIMIT)
              CSIZ2 = IS + CSCAL
              IF (CSIZ2 .GT. ITAIL-IHEAD) THEN
                 INFO (14) = INFO (14) + 1
                 INTZER = 0
                 CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                        II, ISIZE, IHEAD, IUSE,
     $                        CP, RP, DN, N, WIR, WIC, WR, WC,
     $                        INTZER, INTZER, INTZER, INTZER, .TRUE.,
     $                        PFREE, XFREE, MHEAD, MTAIL, SLOTS)
                 PC = CP (COL)
                 CSIZ = II (PC)
              ENDIF
              PC2 = IHEAD
              IHEAD = IHEAD + CSIZ2
              IUSE = IUSE + CSIZ2
              INFO (18) = MAX (INFO (18), IUSE)
              IF (IHEAD .GT. ITAIL) THEN
                 GO TO 9000
              ENDIF
CFPP$ NODEPCHK L
              DO 1480 J = 0, CSCAL + 2*CELN - 1
                 II (PC2+J) = II (PC+J)
1480          CONTINUE
CFPP$ NODEPCHK L
              DO 1490 J = 0, CLEN - 1
                 II (PC2+CSIZ2-CLEN+J) = II (PC+CSIZ-CLEN+J)
1490          CONTINUE
              IF (CLEN .GT. 0) THEN
                 MNEXT = II (PC2+3)
                 MPREV = II (PC2+4)
                 IF (MNEXT .NE. 0) THEN
                    II (MNEXT+4) = PC2
                 ELSE
                    MTAIL = PC2
                 ENDIF
                 IF (MPREV .NE. 0) THEN
                    II (MPREV+3) = PC2
                 ELSE
                    MHEAD = PC2
                 ENDIF
              ENDIF
              CP (COL) = PC2
              II (PC2) = CSIZ2
              II (PC+1) = -1
              II (PC+6) = 0
              PC = PC2
           ENDIF
           CEP = (PC+9)
           II (CEP + 2*CELN  ) = E
           II (CEP + 2*CELN+1) = I - 1
           II (PC+5) = CELN + 1
1500    CONTINUE
        DO 1530 I = 1, FFLEFC
           ROW = WPC (I)
           PR = RP (ROW)
           RSIZ = II (PR)
           RELN = WR (ROW)
           RLEN = WC (ROW)
           WIR (ROW) = -1
           IF (2*(RELN+1) + RLEN + RSCAL .GT. RSIZ) THEN
              IS = 2 * (RELN + 1) + RLEN
              IS = MIN (IS + MAX (16, IS), LIMIT)
              RSIZ2 = IS + RSCAL
              IF (RSIZ2 .GT. ITAIL-IHEAD) THEN
                 INFO (14) = INFO (14) + 1
                 INTZER = 0
                 CALL MA38GD (XX, XSIZE, XHEAD, XUSE,
     $                        II, ISIZE, IHEAD, IUSE,
     $                        CP, RP, DN, N, WIR, WIC, WR, WC,
     $                        INTZER, INTZER, INTZER, INTZER, .TRUE.,
     $                        PFREE, XFREE, MHEAD, MTAIL, SLOTS)
                 PR = RP (ROW)
                 RSIZ = II (PR)
              ENDIF
              PR2 = IHEAD
              IHEAD = IHEAD + RSIZ2
              IUSE = IUSE + RSIZ2
              INFO (18) = MAX (INFO (18), IUSE)
              IF (IHEAD .GT. ITAIL) THEN
                 GO TO 9000
              ENDIF
CFPP$ NODEPCHK L
              DO 1510 J = 0, RSCAL + 2*RELN - 1
                 II (PR2+J) = II (PR+J)
1510          CONTINUE
CFPP$ NODEPCHK L
              DO 1520 J = 0, RLEN - 1
                 II (PR2+RSIZ2-RLEN+J) = II (PR+RSIZ-RLEN+J)
1520          CONTINUE
              RP (ROW) = PR2
              II (PR2) = RSIZ2
              II (PR+1) = -1
              PR = PR2
           ENDIF
           REP = (PR+2)
           II (REP + 2*RELN  ) = E
           II (REP + 2*RELN+1) = I - 1
           WR (ROW) = RELN + 1
1530    CONTINUE
C=======================================================================
C=======================================================================
1540    CONTINUE
2000    CONTINUE
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        IUSE = IUSE - (IHEAD - 1)
        XUSE = XUSE - (XHEAD - 1)
        INEED = IUSE
        XNEED = XUSE
        IHEAD = 1
        XHEAD = 1
        IF (NLU .EQ. 0) THEN
           ITAIL = ISIZE
           XTAIL = XSIZE
           IUSE = IUSE + 1
           XUSE = XUSE + 1
           INEED = IUSE
           XNEED = XUSE
           IP = ITAIL
           XP = XTAIL
        ENDIF
        DO 2010 K = 1, N
           ROW = WPR (N-K+1)
           COL = WPC (N-K+1)
           WIR (K) = ROW
           WIC (K) = COL
2010    CONTINUE
        DO 2020 K = 1, N
           ROW = WIR (K)
           COL = WIC (K)
           WPR (ROW) = K
           WPC (COL) = K
2020    CONTINUE
        IF (PGIVEN) THEN
           DO 2030 ROW = 1, N
              WM (WPR (ROW)) = RPERM (ROW)
2030       CONTINUE
           DO 2040 ROW = 1, N
              RPERM (ROW) = WM (ROW)
2040       CONTINUE
           DO 2050 COL = 1, N
              WM (WPC (COL)) = CPERM (COL)
2050       CONTINUE
           DO 2060 COL = 1, N
              CPERM (COL) = WM (COL)
2060       CONTINUE
        ENDIF
        IS = NLU + 5
        LUIP1 = ITAIL
        ITAIL = ITAIL - IS
        IUSE = IUSE + IS
        INEED = IUSE
        INFO (18) = MAX (INFO (18), IUSE)
        INFO (19) = MAX (INFO (19), INEED)
        IF (IHEAD .LE. ITAIL) THEN
           II (ITAIL+1) = NLU
           II (ITAIL+2) = NPIV
           LUPP = ITAIL+5
           IF (NLU .EQ. 0) THEN
              II (IP) = 0
              XX (XP) = ZERO
           ENDIF
           S = 0
           MAXDR = 1
           MAXDC = 1
           DO 2100 K = 1, N
              E = WIR (K)
              LUIP = RP (E)
              IF (LUIP .GT. 0) THEN
                 S = S + 1
                 II (LUPP+S-1) = LUIP - LUIP1 + 1
                 LUXP = II (LUIP)
                 II (LUIP) = LUXP - XTAIL + 1
                 P = (LUIP + 7)
                 LUDEGC = II (LUIP+3)
                 MAXDC = MAX (MAXDC, LUDEGC)
                 DO 2070 J = 1, LUDEGC
                    II (P) = WPR (ABS (II (P)))
                    P = P + 1
2070             CONTINUE
                 LUDEGR = II (LUIP+2)
                 MAXDR = MAX (MAXDR, LUDEGR)
                 DO 2080 J = 1, LUDEGR
                    II (P) = WPC (ABS (II (P)))
                    P = P + 1
2080             CONTINUE
                 NSONS = II (LUIP+4)
                 DO 2090 J = 1, NSONS
                    ESON = II (P)
                    IF (ESON .LE. N) THEN
                       II (P) = WM (ESON)
                    ELSE IF (ESON .LE. 2*N) THEN
                       II (P) = WM (ESON-N) + N
                    ELSE
                       II (P) = WM (ESON-2*N) + 2*N
                    ENDIF
                    P = P + 1
2090             CONTINUE
                 WM (E) = S
              ENDIF
2100       CONTINUE
           CMAX = MAX (CMAX, MAXDC)
           RMAX = MAX (RMAX, MAXDR)
           TOTNLU = TOTNLU + NLU
           II (ITAIL+3) = MAXDC
           II (ITAIL+4) = MAXDR
           XRUSE = XRUSE - NZ
           RETURN
        ENDIF
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
9000    CONTINUE
        IF (IHEAD .GT. ITAIL .OR. ISIZE .LT. MINMEM) THEN
           CALL MA38ND (1, ICNTL, INFO, -3, INFO (19))
        ENDIF
        IF (XHEAD .GT. XTAIL) THEN
           CALL MA38ND (1, ICNTL, INFO, -4, INFO (21))
        ENDIF
        RETURN
        END
        SUBROUTINE MA38GD (XX, XSIZE, XHEAD, XUSE,
     *          II, ISIZE, IHEAD, IUSE,
     *          CP, RP, DN, N, WIR, WIC, WR, WC,
     *          FFXP, FFSIZE, WXP, FFDIMC, DOSLOT,
     *          PFREE, XFREE, MHEAD, MTAIL, SLOTS)
        INTEGER N, DN, ISIZE, II(ISIZE), IHEAD, RP(N+DN),
     *          CP(N+1), WIR(N), WIC(N), XSIZE, XUSE,
     *          IUSE, XHEAD, FFXP, FFSIZE, WXP,
     *          FFDIMC, WR(N+DN), WC(N+DN), PFREE, XFREE, MHEAD,
     *          MTAIL, SLOTS
        LOGICAL DOSLOT
        DOUBLE PRECISION XX (XSIZE)
C=== MA38GD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTEGER WHAT, FSIZ, ROW, COL, P, IDP, XDP, I, E, EP, FDIMC,
     $          LUDEGR, LUDEGC, J, PC, CELN, CLEN, RELN, RLEN,
     $          CSIZ1, CSIZ2, RSIZ1, RSIZ2, FLUIP, CXP, FXP, RDEG,
     $          CDEG, CSCAL, RSCAL, FSCAL
        PARAMETER (CSCAL = 9, RSCAL = 2, FSCAL = 7)
        LOGICAL SLOT
C=======================================================================
C=======================================================================
        SLOTS = 0
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
CFPP$ NODEPCHK L
        DO 10 COL = 1, N
           PC = CP (COL)
           IF (PC .NE. 0) THEN
              CDEG = II (PC+1)
              CP (COL) = CDEG
              II (PC+1) = COL+N
           ENDIF
10      CONTINUE
CFPP$ NODEPCHK L
        DO 20 ROW = 1, N
           P = RP (ROW)
           RLEN = WC (ROW)
           IF (P .EQ. 0) THEN
              CONTINUE
           ELSE IF (RLEN .GE. 0 .AND. RLEN .LE. N) THEN
              RDEG = II (P+1)
              RP (ROW) = RDEG
              II (P+1) = ROW+2*N
           ELSE IF (WR (ROW) .EQ. -(N+DN+2)) THEN
              CONTINUE
           ELSE
              FDIMC = II (P+1)
              RP (ROW) = FDIMC
              II (P+1) = ROW
           ENDIF
20      CONTINUE
CFPP$ NODEPCHK L
        DO 30 E = N+1, N+DN
           EP = RP (E)
           IF (EP .NE. 0) THEN
              FDIMC = II (EP+1)
              RP (E) = FDIMC
              II (EP+1) = E+2*N
           ENDIF
30      CONTINUE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        XDP = 1
        P = MHEAD
40      CONTINUE
        IF (P .NE. 0) THEN
           WHAT = II (P+1)
           IF (WHAT .GT. 3*N) THEN
              E = WHAT - 2*N
              FXP = II (P+2)
              II (P+2) = XDP
CFPP$ NODEPCHK L
              DO 50 J = 0, RP (E) - 1
                 XX (XDP+J) = XX (FXP+J)
50            CONTINUE
              XDP = XDP + RP (E)
           ELSE IF (WHAT .EQ. -1 .OR. II (P+6) .EQ. 0) THEN
              IF (II (P+4) .NE. 0) THEN
                 II (II (P+4)+3) = II (P+3)
              ELSE
                 MHEAD = II (P+3)
              ENDIF
              IF (II (P+3) .NE. 0) THEN
                 II (II (P+3)+4) = II (P+4)
              ELSE
                 MTAIL = II (P+4)
              ENDIF
           ELSE IF (WHAT .LE. N) THEN
              E = WHAT
              FXP = II (P+2)
              II (P+2) = XDP
              FLUIP = II (P)
              LUDEGR = II (FLUIP+2)
              LUDEGC = II (FLUIP+3)
              FDIMC = RP (E)
              IF (FDIMC .EQ. LUDEGC) THEN
CFPP$ NODEPCHK L
                 DO 60 I = 0, (LUDEGR * LUDEGC) - 1
                    XX (XDP+I) = XX (FXP+I)
60               CONTINUE
              ELSE
                 DO 80 J = 0, LUDEGR - 1
CFPP$ NODEPCHK L
                    DO 70 I = 0, LUDEGC - 1
                       XX (XDP + J*LUDEGC + I) = XX (FXP + J*FDIMC + I)
70                  CONTINUE
80               CONTINUE
                 RP (E) = LUDEGC
              ENDIF
              XDP = XDP + LUDEGR*LUDEGC
           ELSE IF (WHAT .LE. 2*N) THEN
              CXP = II (P+2)
              II (P+2) = XDP
              CLEN = II (P+6)
CFPP$ NODEPCHK L
              DO 90 J = 0, CLEN - 1
                 XX (XDP+J) = XX (CXP+J)
90            CONTINUE
              XDP = XDP + CLEN
           ENDIF
           P = II (P+3)
        GOTO 40
        ENDIF
        PFREE = 0
        XFREE = -1
        IF (FFXP .NE. 0) THEN
CFPP$ NODEPCHK L
           DO 100 I = 0, FFSIZE - 1
              XX (XDP+I) = XX (FFXP+I)
100        CONTINUE
           FFXP = XDP
           XDP = XDP + FFSIZE
        ENDIF
        IF (WXP .NE. 0) THEN
           WXP = XDP
           XDP = XDP + FFDIMC
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        P = 1
        IDP = P
110     CONTINUE
        IF (P .LT. IHEAD) THEN
           WHAT = II (P+1)
           IF (WHAT .GT. 3*N) THEN
              E = WHAT - 2*N
              FSIZ = RP (E) + CSCAL
              II (P+1) = RP (E)
              RP (E) = IDP
CFPP$ NODEPCHK L
              DO 120 I = 0, FSIZ - 1
                 II (IDP+I) = II (P+I)
120           CONTINUE
              IF (II (IDP+4) .NE. 0) THEN
                 II (II (IDP+4)+3) = IDP
              ELSE
                 MHEAD = IDP
              ENDIF
              IF (II (IDP+3) .NE. 0) THEN
                 II (II (IDP+3)+4) = IDP
              ELSE
                 MTAIL = IDP
              ENDIF
              P = P + FSIZ
              IDP = IDP + FSIZ
           ELSE IF (WHAT .EQ. -1) THEN
              P = P + II (P)
           ELSE IF (WHAT .GE. 1 .AND. WHAT .LE. N) THEN
              E = WHAT
              FDIMC = RP (E)
              II (P+1) = FDIMC
              RP (E) = IDP
CFPP$ NODEPCHK L
              DO 130 I = 0, FSCAL - 1
                 II (IDP+I) = II (P+I)
130           CONTINUE
              IF (II (IDP+4) .NE. 0) THEN
                 II (II (IDP+4)+3) = IDP
              ELSE
                 MHEAD = IDP
              ENDIF
              IF (II (IDP+3) .NE. 0) THEN
                 II (II (IDP+3)+4) = IDP
              ELSE
                 MTAIL = IDP
              ENDIF
              P = P + FSCAL
              IDP = IDP + FSCAL
           ELSE IF (WHAT .LE. 2*N) THEN
              CSIZ1 = II (P)
              COL = WHAT - N
              CELN = II (P+5)
              CLEN = II (P+6)
              CSIZ2 = 2*CELN + CLEN + CSCAL
              SLOT = DOSLOT .AND. WIC (COL) .GE. 0 .AND. P .GE. IDP+2
              IF (SLOT) THEN
                 CSIZ2 = CSIZ2 + 2
                 SLOTS = SLOTS + 2
              ENDIF
              CDEG = CP (COL)
              II (P+1) = CDEG
              CP (COL) = IDP
              II (P) = CSIZ2
CFPP$ NODEPCHK L
              DO 140 I = 0, CSCAL + 2*CELN - 1
                 II (IDP+I) = II (P+I)
140           CONTINUE
              IF (CLEN .GT. 0) THEN
                 IF (II (IDP+4) .NE. 0) THEN
                    II (II (IDP+4)+3) = IDP
                 ELSE
                    MHEAD = IDP
                 ENDIF
                 IF (II (IDP+3) .NE. 0) THEN
                    II (II (IDP+3)+4) = IDP
                 ELSE
                    MTAIL = IDP
                 ENDIF
              ENDIF
              P = P + CSIZ1 - CLEN
              IDP = IDP + CSCAL + 2*CELN
              IF (SLOT) THEN
                 IDP = IDP + 2
              ENDIF
CFPP$ NODEPCHK L
              DO 150 I = 0, CLEN - 1
                 II (IDP+I) = II (P+I)
150           CONTINUE
              P = P + CLEN
              IDP = IDP + CLEN
           ELSE
              RSIZ1 = II (P)
              ROW = WHAT - 2*N
              RELN = WR (ROW)
              RLEN = WC (ROW)
              RSIZ2 = 2*RELN + RLEN + RSCAL
              SLOT = DOSLOT .AND. WIR (ROW) .GE. 0 .AND. P .GE. IDP+2
              IF (SLOT) THEN
                 RSIZ2 = RSIZ2 + 2
                 SLOTS = SLOTS + 2
              ENDIF
              RDEG = RP (ROW)
              II (P+1) = RDEG
              RP (ROW) = IDP
              II (P) = RSIZ2
CFPP$ NODEPCHK L
              DO 160 I = 0, RSCAL + 2*RELN - 1
                 II (IDP+I) = II (P+I)
160           CONTINUE
              P = P + RSIZ1 - RLEN
              IDP = IDP + RSCAL + 2*RELN
              IF (SLOT) THEN
                 IDP = IDP + 2
              ENDIF
CFPP$ NODEPCHK L
              DO 170 I = 0, RLEN - 1
                 II (IDP+I) = II (P+I)
170           CONTINUE
              P = P + RLEN
              IDP = IDP + RLEN
           ENDIF
        GOTO 110
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IUSE = IUSE - (IHEAD - IDP)
        IHEAD = IDP
        XUSE = XUSE - (XHEAD - XDP)
        XHEAD = XDP
        RETURN
        END
        SUBROUTINE MA38HD (XX, XSIZE, II, ISIZE, N, NZ, NZDIA, NZOFF,
     *          NBLKS, CP, CPERM, RPERM, PR, PC,
     *          W, ZPERM, BP, OFFP,
     *          PRESRV)
        INTEGER N, NZ, ISIZE, II(ISIZE), NZDIA, NZOFF, NBLKS, CP(N+1),
     *          CPERM(N), RPERM(N), PR(N), PC(N), W(N), ZPERM(N),
     *          BP(N+1), OFFP(N+1), XSIZE
        LOGICAL PRESRV
        DOUBLE PRECISION XX (XSIZE)
C=== MA38HD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        EXTERNAL MC21BD, MC13ED
C=======================================================================
C=======================================================================
        INTEGER COL, NDIAG, I, PO, PB, BLK, P, ROW, K1, K
C=======================================================================
C=======================================================================
        NZDIA = NZ
        NZOFF = 0
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        DO 10 COL = 1, N
           W (COL) = CP (COL+1) - CP (COL)
10      CONTINUE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        CALL MC21BD(N, II, NZ, CP, W, ZPERM, NDIAG, OFFP, CPERM, PR, PC)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        DO 20 COL = 1, N
           OFFP (COL) = CP (ZPERM (COL))
           W (COL) = CP (ZPERM (COL)+1) - CP (ZPERM (COL))
20      CONTINUE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        CALL MC13ED(N, II, NZ, OFFP, W, RPERM, BP, NBLKS, CPERM, PR, PC)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (NBLKS .NE. 1) THEN
           DO 30 COL = 1, N
              CPERM (COL) = ZPERM (RPERM (COL))
30         CONTINUE
           IF (.NOT. PRESRV) THEN
              DO 40 K = 1, N
                 PC (CPERM (K)) = K
                 PR (RPERM (K)) = K
40            CONTINUE
              BP (NBLKS+1) = N+1
              DO 60 BLK = 1, NBLKS
                 DO 50 I = BP (BLK), BP (BLK+1)-1
                    W (I) = BP (BLK)
50               CONTINUE
60            CONTINUE
              PB = NZ + 1
              DO 80 COL = 1, N
                 ZPERM (COL) = PB
                 K1 = W (COL)
CFPP$ NODEPCHK L
                 DO 70 P = CP (CPERM (COL)), CP (CPERM (COL)+1)-1
                    ROW = PR (II (P))
                    IF (W (ROW) .EQ. K1) THEN
                       II (PB) = ROW - K1 + 1
                       XX (PB) = XX (P)
                       PB = PB + 1
                    ENDIF
70               CONTINUE
80            CONTINUE
              NZDIA = PB - (NZ + 1)
              NZOFF = NZ - NZDIA
              PO = 1
              DO 100 COL = 1, N
                 OFFP (COL) = PO
                 K1 = W (PC (COL))
CFPP$ NODEPCHK L
                 DO 90 P = CP (COL), CP (COL+1)-1
                    ROW = PR (II (P))
                    IF (W (ROW) .NE. K1) THEN
                       II (PO) = II (P)
                       XX (PO) = XX (P)
                       PO = PO + 1
                    ENDIF
90               CONTINUE
100           CONTINUE
              OFFP (N+1) = PO
              PB = NZ + 1
CFPP$ NODEPCHK L
              DO 110 I = 0, NZDIA - 1
                 II (PO+I) = II (PB+I)
                 XX (PO+I) = XX (PB+I)
110           CONTINUE
              DO 120 COL = 1, N
                 CP (COL) = ZPERM (COL) - NZDIA
120           CONTINUE
           ENDIF
           BP (NBLKS+1) = N+1
CFPP$ NODEPCHK L
           DO 130 BLK = NBLKS + 1, 1, -1
              BP (BLK + (N-NBLKS)) = BP (BLK)
130        CONTINUE
        ENDIF
        RETURN
        END
        SUBROUTINE MA38JD (N, JOB, TRANSC, LUXSIZ, LUX,
     *          LUISIZ, LUI, B, X, R, Z, LY, Y, S, CNTL, INFO,
     *          RINFO, CPERM, RPERM, AN, ANZ, AP, AI, AX, ON,
     *          NZOFF, OFFP, OFFI, OFFX, NBLKS, LUBLKP, BLKP, IRSTEP)
        INTEGER N, JOB, LUXSIZ, LUISIZ, LUI(LUISIZ), LY, IRSTEP,
     *          INFO(40), CPERM(N), RPERM(N), AN,
     *          ANZ, AP(AN+1), AI(ANZ), ON, NZOFF, OFFP(ON+1),
     *          OFFI(NZOFF), NBLKS, LUBLKP(NBLKS), BLKP(NBLKS+1)
        LOGICAL TRANSC
        DOUBLE PRECISION LUX(LUXSIZ), B(N), X(N), R(N), Z(N), Y(LY),
     *          S(LY), CNTL(10), RINFO(20), AX(ANZ),
     *          OFFX(NZOFF)
C=== MA38JD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC ABS, MAX
        INTEGER IDAMAX
        EXTERNAL MA38LD, MA38ND, MA38SD, MA38TD, MA38UD
C=======================================================================
C=======================================================================
        INTEGER NLU, I, BLK, K1, K2, KN, P, STEP, NPIV, J
        DOUBLE PRECISION
     $          ZERO, ONE, XNORM, TAU, NCTAU, OMEGA1, OMEGA2, D1,
     $          D2, OMEGA, OMLAST, OM1LST, OM2LST, TWO, EPS, MAXEPS,
     $          THOSND, A, AXX, R2, X2, Y2, Z2
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $          MAXEPS = TWO ** (-15), THOSND = 1000.0D0)
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        OMEGA = ZERO
        OMEGA1 = ZERO
        OMEGA2 = ZERO
        EPS = CNTL (3)
        IF (EPS .LE. ZERO .OR. EPS .GT. MAXEPS) THEN
           EPS = MAXEPS
        ENDIF
        NCTAU = THOSND * N * EPS
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (NBLKS .EQ. 1) THEN
           NLU = LUI (2)
           NPIV = LUI (3)
        ENDIF
C-----------------------------------------------------------------------
        IF (JOB .EQ. 1) THEN
C-----------------------------------------------------------------------
           IF (.NOT. TRANSC) THEN
              DO 10 I = 1, N
                 X (I) = B (RPERM (I))
10            CONTINUE
              IF (NBLKS .EQ. 1) THEN
                 CALL MA38LD (NLU, N, LUI(6), LUI(6+NLU), LUX,X,Z)
              ELSE
                 DO 20 BLK = 1, NBLKS
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .GT. 1) THEN
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL MA38LD (NLU, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                    ENDIF
20               CONTINUE
              ENDIF
           ELSE
              DO 30 I = 1, N
                 R (I) = B (I)
30            CONTINUE
              IF (NBLKS .EQ. 1) THEN
                 CALL MA38TD (NLU, NPIV, N, LUI(6), LUI(6+NLU), LUX,R,Z)
              ELSE
                 DO 40 BLK = 1, NBLKS
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .GT. 1) THEN
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL MA38TD (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF
40               CONTINUE
              ENDIF
              DO 50 I = 1, N
                 X (RPERM (I)) = R (I)
50            CONTINUE
           ENDIF
C-----------------------------------------------------------------------
        ELSE IF (JOB .EQ. 2) THEN
C-----------------------------------------------------------------------
           IF (TRANSC) THEN
              DO 60 I = 1, N
                 X (I) = B (CPERM (I))
60            CONTINUE
              IF (NBLKS .EQ. 1) THEN
                 CALL MA38SD (NLU, N, LUI(6), LUI(6+NLU), LUX,X,Z)
              ELSE
                 DO 100 BLK = 1, NBLKS
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    IF (KN .EQ. 1) THEN
                       X (K1) = X (K1) / LUX (LUBLKP (BLK))
                       R (K1) = X (K1)
                    ELSE
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL MA38SD (NLU, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                       DO 70 I = K1, K2
                          R (I) = X (I)
70                     CONTINUE
                       CALL MA38TD (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF
                    DO 90 I = K1, K2
                       R2 = R (I)
                       DO 80 P = OFFP (I), OFFP (I+1)-1
                          X (OFFI (P)) = X (OFFI (P)) - OFFX (P) * R2
80                     CONTINUE
90                  CONTINUE
100              CONTINUE
              ENDIF
           ELSE
              IF (NBLKS .EQ. 1) THEN
                 DO 110 I = 1, N
                    R (I) = B (I)
110              CONTINUE
                 CALL MA38UD (NLU, NPIV, N, LUI(6), LUI(6+NLU), LUX,R,Z)
              ELSE
                 DO 150 BLK = NBLKS, 1, -1
                    K1 = BLKP (BLK)
                    K2 = BLKP (BLK+1) - 1
                    KN = K2-K1+1
                    DO 130 I = K1, K2
                       X2 = ZERO
                       DO 120 P = OFFP (I), OFFP (I+1)-1
                          X2 = X2 + OFFX (P) * R (OFFI (P))
120                    CONTINUE
                       X (I) = X2
130                 CONTINUE
                    IF (KN .EQ. 1) THEN
                       R (K1) = (B (K1) - X (K1)) / LUX (LUBLKP (BLK))
                    ELSE
                       P = LUBLKP (BLK)
                       NLU = LUI (P+1)
                       NPIV = LUI (P+2)
                       CALL MA38LD (NLU, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), X (K1), Z)
                       DO 140 I = K1, K2
                          R (I) = B (I) - X (I)
140                    CONTINUE
                       CALL MA38UD (NLU, NPIV, KN, LUI (P+5),
     $                    LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                    ENDIF
150              CONTINUE
              ENDIF
              DO 160 I = 1, N
                 X (CPERM (I)) = R (I)
160           CONTINUE
           ENDIF
C-----------------------------------------------------------------------
        ELSE
C-----------------------------------------------------------------------
           DO 450 STEP = 0, IRSTEP
              IF (.NOT. TRANSC) THEN
                 IF (STEP .EQ. 0) THEN
                    DO 170 I = 1, N
                       R (I) = B (RPERM (I))
170                 CONTINUE
                 ELSE
                    DO 180 I = 1, N
                       Z (I) = B (I)
180                 CONTINUE
                    DO 200 I = 1, N
                       X2 = X (I)
                       DO 190 P = AP (I), AP (I+1) - 1
                          Z (AI (P)) = Z (AI (P)) - AX (P) * X2
190                    CONTINUE
200                 CONTINUE
                    DO 210 I = 1, N
                       R (I) = Z (RPERM (I))
210                 CONTINUE
                 ENDIF
                 IF (NBLKS .EQ. 1) THEN
                    CALL MA38LD (NLU, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                    CALL MA38UD (NLU, NPIV, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                 ELSE
                    DO 240 BLK = NBLKS, 1, -1
                       K1 = BLKP (BLK)
                       K2 = BLKP (BLK+1) - 1
                       KN = K2-K1+1
                       DO 230 I = K1, K2
                          R2 = R (I)
                          DO 220 P = OFFP (I), OFFP (I+1)-1
                             R2 = R2 - OFFX (P) * R (OFFI (P))
220                       CONTINUE
                          R (I) = R2
230                    CONTINUE
                       IF (KN .EQ. 1) THEN
                          R (K1) = R (K1) / LUX (LUBLKP (BLK))
                       ELSE
                          P = LUBLKP (BLK)
                          NLU = LUI (P+1)
                          NPIV = LUI (P+2)
                          CALL MA38LD (NLU, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                          CALL MA38UD (NLU, NPIV, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                       ENDIF
240                 CONTINUE
                 ENDIF
                 IF (STEP .EQ. 0) THEN
                    DO 250 I = 1, N
                       X (CPERM (I)) = R (I)
250                 CONTINUE
                 ELSE
                    DO 260 I = 1, N
                       X (CPERM (I)) = X (CPERM (I)) + R (I)
260                 CONTINUE
                 ENDIF
              ELSE
                 IF (STEP .EQ. 0) THEN
                    DO 270 I = 1, N
                       R (I) = B (CPERM (I))
270                 CONTINUE
                 ELSE
                    DO 280 I = 1, N
                       Z (I) = B (I)
280                 CONTINUE
                    DO 300 I = 1, N
                       Z2 = Z (I)
                       DO 290 P = AP (I), AP (I+1) - 1
                          Z2 = Z2 - AX (P) * X (AI (P))
290                    CONTINUE
                       Z (I) = Z2
300                 CONTINUE
                    DO 310 I = 1, N
                       R (I) = Z (CPERM (I))
310                 CONTINUE
                 ENDIF
                 IF (NBLKS .EQ. 1) THEN
                    CALL MA38SD (NLU, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                    CALL MA38TD (NLU, NPIV, N,LUI(6),LUI(6+NLU),LUX,R,Z)
                 ELSE
                    DO 340 BLK = 1, NBLKS
                       K1 = BLKP (BLK)
                       K2 = BLKP (BLK+1) - 1
                       KN = K2-K1+1
                       IF (KN .EQ. 1) THEN
                          R (K1) = R (K1) / LUX (LUBLKP (BLK))
                       ELSE
                          P = LUBLKP (BLK)
                          NLU = LUI (P+1)
                          NPIV = LUI (P+2)
                          CALL MA38SD (NLU, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                          CALL MA38TD (NLU, NPIV, KN, LUI (P+5),
     $                       LUI (P+5+NLU), LUX (LUI (P)), R (K1), Z)
                       ENDIF
                       DO 330 I = K1, K2
                          R2 = R (I)
                          DO 320 P = OFFP (I), OFFP (I+1)-1
                             R (OFFI (P)) = R (OFFI (P)) - OFFX (P) * R2
320                       CONTINUE
330                    CONTINUE
340                 CONTINUE
                 ENDIF
                 IF (STEP .EQ. 0) THEN
                    DO 350 I = 1, N
                       X (RPERM (I)) = R (I)
350                 CONTINUE
                 ELSE
                    DO 360 I = 1, N
                       X (RPERM (I)) = X (RPERM (I)) + R (I)
360                 CONTINUE
                 ENDIF
              ENDIF
              IF (IRSTEP .GT. 0) THEN
                 XNORM = ABS (X (IDAMAX (N, X, 1)))
                 DO 370 I = 1, N
                    R (I) = B (I)
                    Z (I) = ZERO
                    Y (I) = ZERO
370              CONTINUE
                 IF (.NOT. TRANSC) THEN
                    DO 390 J = 1, N
                       X2 = X (J)
CFPP$ NODEPCHK L
                       DO 380 P = AP (J), AP (J+1) - 1
                          I = AI (P)
                          A = AX (P)
                          AXX = A * X2
                          R (I) = R (I) -     (AXX)
                          Z (I) = Z (I) + ABS (AXX)
                          Y (I) = Y (I) + ABS (A)
380                    CONTINUE
390                 CONTINUE
                 ELSE
                    DO 410 I = 1, N
                       R2 = R (I)
                       Z2 = Z (I)
                       Y2 = Y (I)
CFPP$ NODEPCHK L
                       DO 400 P = AP (I), AP (I+1) - 1
                          J = AI (P)
                          A = AX (P)
                          AXX = A * X (J)
                          R2 = R2 -     (AXX)
                          Z2 = Z2 + ABS (AXX)
                          Y2 = Y2 + ABS (A)
400                    CONTINUE
                       R (I) = R2
                       Z (I) = Z2
                       Y (I) = Y2
410                 CONTINUE
                 ENDIF
                 OMLAST = OMEGA
                 OM1LST = OMEGA1
                 OM2LST = OMEGA2
                 OMEGA1 = ZERO
                 OMEGA2 = ZERO
                 DO 420 I = 1, N
                    TAU = (Y (I) * XNORM + ABS (B (I))) * NCTAU
                    D1 = Z (I) + ABS (B (I))
                    IF (D1 .GT. TAU) THEN
                       OMEGA1 = MAX (OMEGA1, ABS (R (I)) / D1)
                    ELSE IF (TAU .GT. ZERO) THEN
                       D2 = Z (I) + Y (I) * XNORM
                       OMEGA2 = MAX (OMEGA2, ABS (R (I)) / D2)
                    ENDIF
420              CONTINUE
                 OMEGA = OMEGA1 + OMEGA2
                 RINFO (7) = OMEGA1
                 RINFO (8) = OMEGA2
                 INFO (24) = STEP
                 IF (ONE + OMEGA .LE. ONE) THEN
                    RETURN
                 ENDIF
                 IF (STEP .GT. 0 .AND. OMEGA .GT. OMLAST / TWO) THEN
                    IF (OMEGA .GT. OMLAST) THEN
                       DO 430 I = 1, N
                          X (I) = S (I)
                          RINFO (7) = OM1LST
                          RINFO (8) = OM2LST
430                    CONTINUE
                    ENDIF
                    INFO (24) = STEP - 1
                    RETURN
                 ENDIF
                 DO 440 I = 1, N
                    S (I) = X (I)
440              CONTINUE
              ENDIF
450        CONTINUE
C-----------------------------------------------------------------------
        ENDIF
C-----------------------------------------------------------------------
        RETURN
        END
        SUBROUTINE MA38KD (N, NZ, TRANSA, XX, XSIZE, INFO, ICNTL,
     *                     II, ISIZE, W, WP, WHO)
        INTEGER ISIZE, II(ISIZE), N, NZ, W(N), WP(N+1), INFO(40),
     *          ICNTL(20), XSIZE,WHO
        DOUBLE PRECISION XX (XSIZE)
        LOGICAL TRANSA
C=== MA38KD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC MAX
        EXTERNAL MA38ND, MA38ZD
C=======================================================================
C=======================================================================
        INTEGER ROW, COL, PDEST, P, NZ1, PCOL, IP, XP, IO, PRL, NINVLD,
     $          NDUPL, I
        LOGICAL PR3
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IO = ICNTL (2)
        PRL = ICNTL (3)
        PR3 = PRL .GE. 3 .AND. IO .GE. 0
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        NINVLD = 0
        NDUPL = 0
        DO 10 COL = 1, N
           W (COL) = 0
10      CONTINUE
        NZ1 = NZ
        DO 20 P = NZ, 1, -1
           ROW = II (P)
           COL = II (NZ+P)
           IF (ROW.LT.1.OR.ROW.GT.N.OR.COL.LT.1.OR.COL.GT.N) THEN
              IF (PR3) THEN
                 CALL MA38ZD (WHO, 99, ROW, COL, XX(P), IO)
              ENDIF
              II (P)    = II (NZ1)
              II (NZ+P) = II (NZ+NZ1)
              XX (P)    = XX (NZ1)
              NZ1 = NZ1 - 1
           ELSE
              IF (TRANSA) THEN
                 W (ROW) = W (ROW) + 1
              ELSE
                 W (COL) = W (COL) + 1
              ENDIF
           ENDIF
20      CONTINUE
        NINVLD = NZ - NZ1
        IF (NINVLD .NE. 0) THEN
           CALL MA38ND (WHO, ICNTL, INFO, 1, NINVLD)
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        WP (1) = 1
        DO 30 I = 1, N
           WP (I+1) = WP (I) + W (I)
30      CONTINUE
        DO 40 I = 1, N
           W (I) = WP (I)
40      CONTINUE
        IP = MAX (2*NZ, N+1)
        XP = NZ
        IF (TRANSA) THEN
           DO 50 P = 1, NZ1
              ROW = II (P)
              COL = II (NZ+P)
              II (IP + W (ROW)) = COL
              XX (XP + W (ROW)) = XX (P)
              W (ROW) = W (ROW) + 1
50         CONTINUE
        ELSE
           DO 60 P = 1, NZ1
              ROW = II (P)
              COL = II (NZ+P)
              II (IP + W (COL)) = ROW
              XX (XP + W (COL)) = XX (P)
              W (COL) = W (COL) + 1
60         CONTINUE
        ENDIF
        NZ = NZ1
CFPP$ NODEPCHK L
        DO 70 P = 1, NZ
           II (N+1+P) = II (IP+P)
           XX (P) = XX (XP+P)
70      CONTINUE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        DO 80 ROW = 1, N
           W (ROW) = 0
80      CONTINUE
        PDEST = 1
        DO 100 COL = 1, N
           PCOL = PDEST
           DO 90 P = WP (COL), WP (COL+1)-1
              ROW = II (N+1+P)
              IF (W (ROW) .GE. PCOL) THEN
                 XX (W (ROW)) = XX (W (ROW)) + XX (P)
                 IF (PR3) THEN
                    IF (TRANSA) THEN
                       CALL MA38ZD (WHO, 98, COL, ROW, XX (P), IO)
                    ELSE
                       CALL MA38ZD (WHO, 98, ROW, COL, XX (P), IO)
                    ENDIF
                 ENDIF
              ELSE
                 W (ROW) = PDEST
                 IF (PDEST .NE. P) THEN
                    II (N+1+PDEST) = ROW
                    XX (PDEST) = XX (P)
                 ENDIF
                 PDEST = PDEST + 1
              ENDIF
90         CONTINUE
           WP (COL) = PCOL
100     CONTINUE
        WP (N+1) = PDEST
        NZ1 = PDEST - 1
        NDUPL = NZ - NZ1
        IF (NDUPL .NE. 0) THEN
           CALL MA38ND (WHO, ICNTL, INFO, 2, NDUPL)
        ENDIF
        NZ = NZ1
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        DO 110 COL = 1, N+1
           II (COL) = WP (COL)
110     CONTINUE
        INFO (2) = NDUPL
        INFO (3) = NINVLD
        INFO (5) = NZ
        INFO (6) = NZ
        INFO (7) = 0
        IF (NZ .EQ. 0) THEN
           CALL MA38ND (WHO, ICNTL, INFO, -2, -1)
        ENDIF
        RETURN
        END
        SUBROUTINE MA38LD (NLU, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, N, LUP(NLU), LUI(*)
        DOUBLE PRECISION LUX(*), X(N), W(N)
C=== MA38LD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        EXTERNAL DTRSV, DGEMV
C=======================================================================
C=======================================================================
        INTEGER I, K, S, LUIP, LUXP, LUK, LUDEGC, LUCP, LXP, ROW
        DOUBLE PRECISION
     $          ONE
        PARAMETER (ONE = 1.0D0)
C=======================================================================
C=======================================================================
        K = 0
        DO 40 S = 1, NLU
           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LXP    = LUXP + LUK
           IF (LUK .EQ. 1) THEN
              K = K + 1
CFPP$ NODEPCHK L
              DO 10 I = 1, LUDEGC
                 ROW = LUI (LUCP+I-1)
                 X (ROW) = X (ROW) - LUX (LXP+I-1) * X (K)
10            CONTINUE
           ELSE
              CALL DTRSV ('L', 'N', 'U', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)
              DO 20 I = 1, LUDEGC
                 ROW = LUI (LUCP+I-1)
                 W (I) = X (ROW)
20            CONTINUE
              CALL DGEMV ('N', LUDEGC, LUK, -ONE,
     $           LUX (LXP), LUDEGC + LUK, X (K+1), 1, ONE, W, 1)
              DO 30 I = 1, LUDEGC
                 ROW = LUI (LUCP+I-1)
                 X (ROW) = W (I)
30            CONTINUE
              K = K + LUK
           ENDIF
40      CONTINUE
        RETURN
        END
        SUBROUTINE MA38MD (W, N, RPERM, CPERM, NZOFF,
     *          OFFP, OFFI, OFFX, PR,
     *          ICNTL, MP, MI, MX, MN, MNZ, PRESRV, NBLKS, BLKP,
     *          ONZ, WHO, NBELOW)
        INTEGER N, NZOFF, W(N+1), RPERM(N), CPERM(N), ONZ,
     *          OFFP(N+1), OFFI(ONZ), PR(N), ICNTL(20), MN, MNZ,
     *          MP(MN+1), MI(MNZ), NBLKS, BLKP(NBLKS+1), WHO, NBELOW
        LOGICAL PRESRV
        DOUBLE PRECISION OFFX(ONZ), MX(MNZ)
C=== MA38MD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        EXTERNAL MA38ZD
C=======================================================================
C=======================================================================
        INTEGER ROW, COL, P, BLK, K, K1, K2, IO, PRL
        LOGICAL PR3
C=======================================================================
C=======================================================================
        IO = ICNTL (2)
        PRL = ICNTL (3)
        PR3 = PRL .GE. 3 .AND. IO .GE. 0
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
CFPP$ NODEPCHK L
        DO 10 K = 1, N
           PR (RPERM (K)) = K
10      CONTINUE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        W (1) = 1
        DO 20 ROW = 2, N
           W (ROW) = 0
20      CONTINUE
        NBELOW = 0
        IF (PRESRV) THEN
           DO 50 BLK = 1, NBLKS
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              DO 40 COL = K1, K2
CFPP$ NODEPCHK L
                 DO 30 P = MP (CPERM (COL)), MP (CPERM (COL)+1)-1
                    ROW = PR (MI (P))
                    IF (ROW .LT. K1) THEN
                       W (ROW) = W (ROW) + 1
                    ELSE IF (ROW .GT. K2 .AND. WHO .EQ. 2) THEN
                       IF (PR3) THEN
                          CALL MA38ZD (2, 96, MI (P), COL, MX (P), IO)
                       ENDIF
                       NBELOW = NBELOW + 1
                    ENDIF
30               CONTINUE
40            CONTINUE
50         CONTINUE
        ELSE
           DO 70 COL = 1, N
CFPP$ NODEPCHK L
              DO 60 P = OFFP (COL), OFFP (COL+1) - 1
                 ROW = PR (MI (P))
                 W (ROW) = W (ROW) + 1
60            CONTINUE
70         CONTINUE
        ENDIF
        DO 80 ROW = 2, N
           W (ROW) = W (ROW) + W (ROW-1)
80      CONTINUE
        W (N+1) = W (N)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (PRESRV) THEN
           DO 110 BLK = NBLKS, 1, -1
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              DO 100 COL = K2, K1, - 1
CFPP$ NODEPCHK L
                 DO 90 P = MP (CPERM (COL)), MP (CPERM (COL)+1)-1
                    ROW = PR (MI (P))
                    IF (ROW .LT. K1) THEN
                       W (ROW) = W (ROW) - 1
                       OFFI (W (ROW)) = COL
                       OFFX (W (ROW)) = MX (P)
                    ENDIF
90               CONTINUE
100           CONTINUE
110        CONTINUE
        ELSE
           DO 130 COL = N, 1, -1
CFPP$ NODEPCHK L
              DO 120 P = OFFP (CPERM (COL)), OFFP (CPERM (COL) + 1) - 1
                 ROW = PR (MI (P))
                 W (ROW) = W (ROW) - 1
                 OFFI (W (ROW)) = COL
                 OFFX (W (ROW)) = MX (P)
120           CONTINUE
130        CONTINUE
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        DO 140 ROW = 1, N+1
           OFFP (ROW) = W (ROW)
140     CONTINUE
        NZOFF = OFFP (N+1) - 1
        RETURN
        END
        SUBROUTINE MA38ND (WHO, ICNTL, INFO, ERROR, S)
        INTEGER WHO, ICNTL(20), INFO(40), ERROR, S
C=== MA38ND ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC MOD
        EXTERNAL MA38ZD
C=======================================================================
C=======================================================================
        LOGICAL BOTH
        DOUBLE PRECISION
     $          ZERO
        PARAMETER (ZERO = 0.0D0)
        INTEGER IOERR, PRL
C=======================================================================
C=======================================================================
        IOERR = ICNTL (1)
        PRL = ICNTL (3)
        IF (ERROR .LT. 0) THEN
           BOTH = (INFO (1) .EQ. -3 .AND. ERROR .EQ. -4) .OR.
     $            (INFO (1) .EQ. -4 .AND. ERROR .EQ. -3)
           IF (BOTH) THEN
              INFO (1) = -5
           ELSE
              INFO (1) = ERROR
           ENDIF
           IF (PRL .GE. 1) THEN
              CALL MA38ZD (WHO, ERROR, S, 0, ZERO, IOERR)
           ENDIF
        ELSE IF (ERROR .GT. 0) THEN
           IF (INFO (1) .GE. 0) THEN
              IF (MOD (INFO (1) / ERROR, 2) .EQ. 0) THEN
                 INFO (1) = INFO (1) + ERROR
              ENDIF
           ENDIF
           IF (PRL .GE. 2) THEN
              CALL MA38ZD (WHO, ERROR, S, 0, ZERO, IOERR)
           ENDIF
        ENDIF
        RETURN
        END
        SUBROUTINE MA38OD (XX, XSIZE, XHEAD, XUSE,
     *          LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     *          FFXP, FFSIZE, PFREE, XFREE)
        INTEGER LUI(*), NLU, FRDIMC(NLU+2), FRXP(NLU+2),
     *          FRNEXT(NLU+2), FRPREV(NLU+2), LUP(NLU),
     *          XSIZE, XUSE, XHEAD, FFXP, FFSIZE,
     *          PFREE, XFREE
        DOUBLE PRECISION XX(XSIZE)
C=== MA38OD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC ABS
C=======================================================================
C=======================================================================
        INTEGER XDP, I, E, FDIMC, LUDEGR, LUDEGC, J, FLUIP, FXP,
     $          MHEAD, MTAIL
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        MHEAD = NLU+1
        MTAIL = NLU+2
        XDP = FRXP (MHEAD)
        E = FRNEXT (MHEAD)
10      CONTINUE
        IF (E .NE. MTAIL) THEN
           FDIMC = FRDIMC (E)
           IF (FDIMC .EQ. 0) THEN
              FRNEXT (FRPREV (E)) = FRNEXT (E)
              FRPREV (FRNEXT (E)) = FRPREV (E)
           ELSE
              FXP = FRXP (E)
              FRXP (E) = XDP
              FLUIP = LUP (E)
              LUDEGR = ABS (LUI (FLUIP+2))
              LUDEGC = ABS (LUI (FLUIP+3))
              IF (FDIMC .EQ. LUDEGC) THEN
CFPP$ NODEPCHK L
                 DO 20 I = 0, (LUDEGR * LUDEGC) - 1
                    XX (XDP+I) = XX (FXP+I)
20               CONTINUE
              ELSE
                 DO 40 J = 0, LUDEGR - 1
CFPP$ NODEPCHK L
                    DO 30 I = 0, LUDEGC - 1
                       XX (XDP + J*LUDEGC + I) = XX (FXP + J*FDIMC + I)
30                  CONTINUE
40               CONTINUE
                 FRDIMC (E) = LUDEGC
              ENDIF
              XDP = XDP + LUDEGR*LUDEGC
           ENDIF
           E = FRNEXT (E)
        GOTO 10
        ENDIF
        FRXP (MTAIL) = XDP
        PFREE = 0
        XFREE = -1
        IF (FFXP .NE. 0) THEN
CFPP$ NODEPCHK L
           DO 50 I = 0, FFSIZE - 1
              XX (XDP+I) = XX (FFXP+I)
50         CONTINUE
           FFXP = XDP
           XDP = XDP + FFSIZE
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        XUSE = XUSE - (XHEAD - XDP)
        XHEAD = XDP
        RETURN
        END
        SUBROUTINE MA38PD (N, NZ, CP, XX, XSIZE, II, ISIZE, XTAIL,
     *          ITAIL, IUSE, XUSE, NZOFF, NBLKS, ICNTL, INFO,
     *          RINFO, PRESRV, AP, AI, AX, AN, ANZ, LUI, LUISIZ,
     *          LUBLKP, BLKP, OFFP, ON, CPERM, RPERM, NE)
        INTEGER N, NZ, ISIZE, II(ISIZE), ICNTL(20), INFO(40),
     *          CP(N+1), XSIZE, XTAIL, ITAIL, IUSE, XUSE, AN, ANZ,
     *          AP(AN+1), AI(ANZ), LUISIZ, LUI(LUISIZ), NBLKS,
     *          LUBLKP(NBLKS), BLKP(NBLKS+1), ON, OFFP(ON+1),
     *          CPERM(N), RPERM(N), NZOFF, NE
        LOGICAL PRESRV
        DOUBLE PRECISION XX(XSIZE), RINFO(20), AX(ANZ)
C=== MA38PD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC MAX
        EXTERNAL  MA38MD, MA38ND, MA38QD, MA38RD
C=======================================================================
C=======================================================================
        INTEGER I, NZDIA, P, IHEAD, NSGLTN, NSYM, WP, ARIP, ARXP, NPIV,
     $          WRKSIZ, NLU, PRP, MC, MR, DUMMY1, DUMMY2, NZ2, K, BLK,
     $          K1, K2, KN, NZBLK, COL, ROW, PRL, IO, LUIP, MNZ, ARNZ,
     $          XHEAD, OFFIP, OFFXP, NOUTSD, NBELOW, NZORIG, XRMAX
        DOUBLE PRECISION
     $          ZERO, ONE, A
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C=======================================================================
C=======================================================================
        IO = ICNTL (2)
        PRL = ICNTL (3)
        NZORIG = NZ
        IF (PRESRV) THEN
           IHEAD = 1
           XHEAD = 1
        ELSE
           IHEAD = NZ + 1
           XHEAD = NZ + 1
        ENDIF
        NZOFF = 0
        NZDIA = NZ
        NSGLTN = 0
        NPIV = 0
        NOUTSD = 0
        NBELOW = 0
        ITAIL = ISIZE + 1
        XTAIL = XSIZE + 1
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        XRMAX = 2*NE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (NBLKS .GT. 1 .AND. PRESRV) THEN
           DO 10 K = 1, N
              OFFP (RPERM (K)) = K
10         CONTINUE
        ELSE
           ITAIL = ITAIL - (2*N+1)
           IUSE = IUSE + 2*N+1
           PRP = ITAIL
           WP = PRP + N
           IUSE = IUSE + NZ
           XUSE = XUSE + NZ
           ARXP = XHEAD
           ARIP = IHEAD
           IHEAD = IHEAD + NZ
           XHEAD = XHEAD + NZ
           INFO (18) = MAX (INFO (18), IUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XUSE)
           IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
              GO TO 9000
           ENDIF
           IF (NBLKS .EQ. 1) THEN
              IF (PRESRV) THEN
                 CALL MA38QD (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $              II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $              ICNTL, AP, BLKP, AI, AX, OFFP, ON, NZ,
     $              0, N, NZ2, I)
              ELSE
                 CALL MA38QD (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $              II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $              ICNTL, CP, BLKP, II, XX, OFFP, ON, NZ,
     $              0, N, NZ2, I)
              ENDIF
           ELSE
              CALL MA38QD (PRESRV, N, NZ, CPERM, RPERM, II (PRP),
     $           II (WP), NBLKS, XX (ARXP), II (ARIP), NZOFF, NZDIA,
     $           ICNTL, CP, BLKP, II, XX, OFFP, ON, NZ,
     $           0, N, NZ2, NBELOW)
           ENDIF
           DO 20 I = 1, N+1
              CP (I) = II (WP+I-1)
20         CONTINUE
           IUSE = IUSE - (2*N+1)
           IF (.NOT. PRESRV) THEN
              XUSE = XUSE - NZ
              IUSE = IUSE - NZ
           ENDIF
           ITAIL = ISIZE + 1
           XTAIL = XSIZE + 1
           NZ = NZ2
           IHEAD = NZ + 1
           XHEAD = NZ + 1
        ENDIF
        INFO (5) = NZ
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (4) = NBELOW
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (NBLKS .EQ. 1) THEN
           NLU = LUI (2)
           MC = LUI (4)
           MR = LUI (5)
           WRKSIZ = 2*N + MR + 3*MC + 4*(NLU+2)
           ITAIL = ITAIL - WRKSIZ
           IUSE = IUSE + WRKSIZ
           P = ITAIL
           INFO (18) = MAX (INFO (18), IUSE)
           IF (IHEAD .GT. ITAIL) THEN
              GO TO 9000
           ENDIF
           CALL MA38RD (CP, NZ, N, XTAIL,
     $          XX, XSIZE, XUSE, II, CPERM, RPERM,
     $          ICNTL, INFO, RINFO, MC, MR,
     $          II (P), II (P+N), II (P+2*N), II (P+2*N+MR),
     $          II (P+2*N+MR+MC), II (P+2*N+MR+2*MC),
     $          II (P+2*N+MR+3*MC), II (P+2*N+MR+3*MC+(NLU+2)),
     $          II (P+2*N+MR+3*MC+2*(NLU+2)),
     $          II (P+2*N+MR+3*MC+3*(NLU+2)),
     $          NLU, LUI (6), LUI (NLU+6), NOUTSD,
     $          XRMAX)
           IF (INFO (1) .LT. 0) THEN
              GO TO 9010
           ENDIF
           IUSE = IUSE - WRKSIZ - NZ
           ITAIL = ITAIL + WRKSIZ
           LUI (1) = 1
           IHEAD = 1
           XHEAD = 1
        ELSE
           IF (PRESRV) THEN
              NZOFF = 0
           ENDIF
           DO 70 BLK = NBLKS, 1, -1
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
              KN = K2-K1+1
              A = ZERO
              IF (PRESRV) THEN
                 IF (KN .GT. 1) THEN
                    NZBLK = 0
                    DO 40 K = K1, K2
                       COL = CPERM (K)
CFPP$ NODEPCHK L
                       DO 30 P = AP (COL), AP (COL+1) - 1
                          ROW = OFFP (AI (P))
                          IF (ROW .LT. K1) THEN
                             NZOFF = NZOFF + 1
                          ELSE IF (ROW .LE. K2) THEN
                             NZBLK = NZBLK + 1
                          ENDIF
30                     CONTINUE
40                  CONTINUE
                    ITAIL = ITAIL - (KN+1)
                    WP = ITAIL
                    IHEAD = NZBLK + 1
                    XHEAD = NZBLK + 1
                    IUSE = IUSE + NZBLK + KN+1
                    XUSE = XUSE + NZBLK
                    XRMAX = MAX (XRMAX, XUSE)
                    INFO (18) = MAX (INFO (18), IUSE)
                    INFO (20) = MAX (INFO (20), XUSE)
                    INFO (21) = MAX (INFO (21), XUSE)
                    IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
                       GO TO 9000
                    ENDIF
                    CALL MA38QD (PRESRV, N, NZ, CPERM, RPERM, OFFP,
     $                 II (WP), NBLKS, XX, II, DUMMY1, DUMMY2,
     $                 ICNTL, AP, BLKP, AI, AX, OFFP, 0, NZBLK,
     $                 BLK, KN, NZ2, I)
                    DO 50 I = 0, KN
                       CP (K1+I) = II (WP+I)
50                  CONTINUE
                    IUSE = IUSE - (KN+1)
                    ITAIL = ITAIL + (KN+1)
                 ELSE
                    COL = CPERM (K1)
                    DO 60 P = AP (COL), AP (COL + 1) - 1
                       ROW = OFFP (AI (P))
                       IF (ROW .LT. K1) THEN
                          NZOFF = NZOFF + 1
                       ELSE IF (ROW .EQ. K1) THEN
                          A = AX (P)
                       ENDIF
60                  CONTINUE
                    IHEAD = 1
                    XHEAD = 1
                 ENDIF
              ELSE
                 IF (BLK .EQ. 1) THEN
                    CP (K2+1) = NZOFF + 1
                 ELSE
                    CP (K2+1) = CP (BLKP (BLK-1))
                 ENDIF
                 IHEAD = CP (K1)
                 XHEAD = IHEAD
                 IF (KN .EQ. 1) THEN
                    IF (CP (K1) .GT. CP (K1+1)) THEN
                       A = XX (CP (K1) - 1)
                       IHEAD = IHEAD - 1
                       XHEAD = XHEAD - 1
                       IUSE = IUSE - 1
                       XUSE = XUSE - 1
                    ENDIF
                 ENDIF
              ENDIF
              IF (KN .GT. 1) THEN
                 ARNZ = CP (K1) - 1
                 LUIP = LUBLKP (BLK)
                 NLU = LUI (LUIP+1)
                 MC = LUI (LUIP+3)
                 MR = LUI (LUIP+4)
                 WRKSIZ = 2*KN + MR + 3*MC + 4*(NLU+2)
                 ITAIL = ITAIL - WRKSIZ
                 IUSE = IUSE + WRKSIZ
                 P = ITAIL
                 INFO (18) = MAX (INFO (18), IUSE)
                 IF (IHEAD .GT. ITAIL) THEN
                    GO TO 9000
                 ENDIF
                 CALL MA38RD (CP (K1), ARNZ, KN, XTAIL,
     $                XX, XTAIL-1, XUSE, II, CPERM (K1), RPERM (K1),
     $                ICNTL, INFO, RINFO, MC, MR,
     $                II (P), II (P+KN), II (P+2*KN), II (P+2*KN+MR),
     $                II (P+2*KN+MR+MC), II (P+2*KN+MR+2*MC),
     $                II (P+2*KN+MR+3*MC), II (P+2*KN+MR+3*MC+(NLU+2)),
     $                II (P+2*KN+MR+3*MC+2*(NLU+2)),
     $                II (P+2*KN+MR+3*MC+3*(NLU+2)),
     $                NLU, LUI (LUIP+5), LUI (LUIP+NLU+5), NOUTSD,
     $                XRMAX)
                 IF (INFO (1) .LT. 0) THEN
                    GO TO 9010
                 ENDIF
                 IUSE = IUSE - WRKSIZ
                 ITAIL = ITAIL + WRKSIZ
                 LUI (LUIP) = XTAIL
                 IUSE = IUSE - (IHEAD - CP (K2+1))
                 IHEAD = CP (K2+1)
                 XHEAD = IHEAD
              ELSE
                 NSGLTN = NSGLTN + 1
                 IF (A .EQ. ZERO) THEN
                    A = ONE
                 ELSE
                    NPIV = NPIV + 1
                 ENDIF
                 XTAIL = XTAIL - 1
                 XUSE = XUSE + 1
                 XRMAX = MAX (XRMAX, XUSE)
                 INFO (20) = MAX (INFO (20), XUSE)
                 INFO (21) = MAX (INFO (21), XUSE)
                 IF (XHEAD .GT. XTAIL) THEN
                    GO TO 9000
                 ENDIF
                 XX (XTAIL) = A
                 LUBLKP (BLK) = -XTAIL
              ENDIF
70         CONTINUE
CFPP$ NODEPCHK L
           DO 80 BLK = 1, NBLKS
              IF (LUBLKP (BLK) .GT. 0) THEN
                 LUI (LUBLKP (BLK)) = LUI (LUBLKP (BLK)) - XTAIL + 1
              ELSE
                 LUBLKP (BLK) = (-LUBLKP (BLK)) - XTAIL + 1
              ENDIF
80         CONTINUE
           IF (PRESRV) THEN
              PRP = IHEAD
              IHEAD = IHEAD + N
              IUSE = IUSE + N
              ITAIL = ITAIL - NZOFF
              OFFIP = ITAIL
              XTAIL = XTAIL - NZOFF
              OFFXP = XTAIL
              IUSE = IUSE + NZOFF
              XUSE = XUSE + NZOFF
              XRMAX = MAX (XRMAX, XUSE)
              INFO (18) = MAX (INFO (18), IUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XUSE)
              IF (IHEAD .GT. ITAIL .OR. XHEAD .GT. XTAIL) THEN
                 GO TO 9000
              ENDIF
              MNZ = NZOFF
              IF (NZOFF .EQ. 0) THEN
                 OFFIP = 1
                 OFFXP = 1
              ENDIF
              CALL MA38MD (CP, N, RPERM, CPERM, NZOFF,
     $             OFFP, II (OFFIP), XX (OFFXP), II (PRP),
     $             ICNTL, AP, AI, AX, AN, ANZ, PRESRV, NBLKS, BLKP,
     $             MNZ, 2, NBELOW)
              IHEAD = 1
              XHEAD = 1
              IUSE = IUSE - N
           ELSE
              DO 90 I = NZOFF, 1, -1
                 II (ITAIL+I-NZOFF-1) = II (I)
                 XX (XTAIL+I-NZOFF-1) = XX (I)
90            CONTINUE
              IHEAD = 1
              XHEAD = 1
              ITAIL = ITAIL - NZOFF
              XTAIL = XTAIL - NZOFF
           ENDIF
        ENDIF
        DO 100 I = 1, LUISIZ
           LUI (I) = ABS (LUI (I))
100     CONTINUE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
9000    CONTINUE
        IF (IHEAD .GT. ITAIL) THEN
           CALL MA38ND (2, ICNTL, INFO, -3, INFO (18))
        ENDIF
        IF (XHEAD .GT. XTAIL) THEN
           CALL MA38ND (2, ICNTL, INFO, -4, INFO (21))
        ENDIF
9010    CONTINUE
        IUSE = IUSE - (N+1)
        INFO (4) = NOUTSD + NBELOW
        NZDIA = NZORIG - NZOFF - NOUTSD - NBELOW
        INFO (5) = NZOFF + NZDIA
        INFO (6) = NZDIA
        INFO (7) = NZOFF
        INFO (8) = NSGLTN
        INFO (9) = NBLKS
        INFO (12) = INFO (10) + INFO (11) + N + INFO (7)
        NSYM = 0
        DO 110 K = 1, N
           IF (CPERM (K) .EQ. RPERM (K)) THEN
              NSYM = NSYM + 1
           ENDIF
110     CONTINUE
        INFO (16) = NSYM
        INFO (17) = INFO (17) + NPIV
        RINFO (1) = RINFO (4) + RINFO (5) + RINFO (6)
        IF (INFO (4) .GT. 0) THEN
           CALL MA38ND (2, ICNTL, INFO, 1, -INFO (4))
        ENDIF
        IF (INFO (1) .GE. 0 .AND. INFO (17) .LT. N) THEN
           CALL MA38ND (2, ICNTL, INFO, 4, INFO (17))
        ENDIF
        INFO (23) = XRMAX
        RETURN
        END
        SUBROUTINE MA38QD (PRESRV, N, NZ, CPERM, RPERM, PR,
     *          W, NBLKS, ARX, ARI, NZOFF, NZDIA,
     *          ICNTL, MP, BLKP, MI, MX, OFFP, ON, NZBLK,
     *          CBLK, KN, NZ2, NBELOW)
        INTEGER N, NZ, CPERM(N), RPERM(N), PR(N), KN, W(KN+1),
     *          NBLKS, NZBLK, ARI(NZBLK), NZOFF, NZDIA, MP(N+1),
     *          MI(NZ), ON, ICNTL(20), BLKP(NBLKS+1), NZ2,
     *          OFFP(ON+1), CBLK, NBELOW
        LOGICAL PRESRV
        DOUBLE PRECISION ARX(NZBLK), MX(NZ)
C=== MA38QD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC MIN
        EXTERNAL MA38MD
C=======================================================================
C=======================================================================
        INTEGER I, P, ROW, COL, BLK, BASE, K1, K2, K, B1, B2, K0
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        NZOFF = 0
        NBELOW = 0
        IF (NBLKS .EQ. 1) THEN
           DO 10 K = 1, N
              PR (RPERM (K)) = K
10         CONTINUE
        ELSE IF (NBLKS .GT. 1 .AND. .NOT. PRESRV) THEN
           CALL MA38MD (W, N, RPERM, CPERM, NZOFF,
     $        OFFP, ARI, ARX, PR,
     $        ICNTL, MP, MI, MX, N, NZ, .TRUE., NBLKS, BLKP,
     $        NZ, 2, NBELOW)
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        DO 20 I = 1, KN+1
           W (I) = 0
20      CONTINUE
        BASE = NZOFF + 1
        IF (CBLK .NE. 0) THEN
           K0 = BLKP (CBLK) - 1
           B1 = CBLK
           B2 = CBLK
        ELSE
           K0 = 0
           B1 = 1
           B2 = NBLKS
        ENDIF
        DO 80 BLK = B1, B2
           IF (NBLKS .GT. 1) THEN
              K1 = BLKP (BLK)
              K2 = BLKP (BLK+1) - 1
           ELSE
              K1 = 1
              K2 = N
           ENDIF
           DO 40 COL = K1, K2
              DO 30 P = MP (CPERM (COL)), MP (CPERM (COL) + 1) - 1
                 ROW = PR (MI (P))
                 IF (ROW .GE. K1 .AND. ROW .LE. K2) THEN
                    I = MIN (ROW, COL) - K0
                    W (I) = W (I) + 1
                 ENDIF
30            CONTINUE
40         CONTINUE
           W (K2-K0+1) = W (K2-K0) + BASE
           DO 50 I = K2-K0, K1-K0+1, -1
              W (I) = W (I+1) + W (I-1)
50         CONTINUE
           W (K1-K0) = W (K1-K0+1)
           DO 70 COL = K1, K2
              DO 60 P = MP (CPERM (COL)), MP (CPERM (COL) + 1) - 1
                 ROW = PR (MI (P))
                 IF (ROW .GE. K1 .AND. ROW .LE. K2) THEN
                    IF (ROW .GE. COL) THEN
                       I = COL - K0 + 1
                       W (I) = W (I) - 1
                       ARI (W (I)) = ROW - K1 + 1
                       ARX (W (I)) = MX (P)
                    ELSE
                       I = ROW - K0 + 1
                       W (I) = W (I) - 1
                       ARI (W (I)) = -(COL - K1 + 1)
                       ARX (W (I)) = MX (P)
                    ENDIF
                 ENDIF
60            CONTINUE
70         CONTINUE
           BASE = W (K1-K0)
           W (K2-K0+1) = 0
80      CONTINUE
        W (KN+1) = NZOFF + 1
        NZDIA = BASE - NZOFF - 1
        NZ2 = NZOFF + NZDIA
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (.NOT. PRESRV) THEN
           DO 90 I = 1, NZ
              MI (I) = ARI (I)
              MX (I) = ARX (I)
90         CONTINUE
        ENDIF
        RETURN
        END
        SUBROUTINE MA38RD (CP, NZ, N, XTAIL, XX, XSIZE, XUSE, ARI,
     *          CPERM, RPERM, ICNTL, INFO, RINFO, MC, MR,
     *          WIR, WIC, WPR, WPC, WM, WJ, FRDIMC, FRXP, FRNEXT,
     *          FRPREV, NLU, LUP, LUI, NOUTSD, XRMAX)
        INTEGER XSIZE, ICNTL(20), INFO(40), N, CPERM(N), RPERM(N),
     *          XTAIL, NZ, ARI(NZ), CP(N+1), MR, MC, NOUTSD,
     *          WIR(N), WIC(N), WPR(MR), XRMAX, WPC(MC), WM(MC),
     *          NLU, FRDIMC(NLU+2), FRXP(NLU+2), XUSE, WJ(MC),
     *          FRNEXT(NLU+2), FRPREV(NLU+2), LUP(NLU), LUI(*)
        DOUBLE PRECISION XX(XSIZE), RINFO(20)
C=== MA38RD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC ABS, MAX
        EXTERNAL MA38ND, MA38OD, MA38ZD, DGEMV, DGEMM, DTRSV, DTRSM
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C=======================================================================
C=======================================================================
        INTEGER SWPCOL, SWPROW, FDIMC, K0, COLPOS, ROWPOS, PIVOT, FFPP,
     $          P, I, J, LUDEGR, LUDEGC, KPOS, SP, FFRP, FFCP, TYPE,
     $          FXP, LURP, LUCP, NEXT, FFLEFR, PREV, XHEAD, FDEGR,
     $          FFLEFC, K, XCDP, XDP, XSP, S, FDEGC, FLURP, FLUCP,
     $          COL, E, ROW, MHEAD, MTAIL, UXP, LUK, IO, FLUIP, LUSONP,
     $          FFSIZE, FFXP, FFDIMR, FFDIMC, XRDP, NPIV, NB, LUNSON,
     $          XNEED, LDIMR, LDIMC, LXP, PRL, XP, LUIP, PFREE, XFREE,
     $          XS, LUXP, FSP, FLP, FDP, DEGC, NZU, NZL, XRUSE
        LOGICAL PR3, ALLCOL, ALLROW
        DOUBLE PRECISION
     $          ONE, ZERO, X
        PARAMETER (ONE = 1.0D0, ZERO = 0.0D0)
C=======================================================================
C=======================================================================
        IO = ICNTL (2)
        PRL = ICNTL (3)
        NB = MAX (1, ICNTL (7))
        NPIV = 0
        XHEAD = CP (1)
        XTAIL = XSIZE + 1
        XNEED = XUSE
        XRUSE = XUSE
        XRMAX = MAX (XRMAX, XRUSE)
        MHEAD = NLU+1
        MTAIL = NLU+2
        XFREE = -1
        PFREE = 0
        PR3 = PRL .GE. 3 .AND. IO .GE. 0
        DO 10 I = 1, N
           WIR (I) = -1
           WIC (I) = -1
10      CONTINUE
        DO 20 E = 1, NLU+2
           FRDIMC (E) = 0
           FRXP (E) = 0
           FRNEXT (E) = 0
           FRPREV (E) = 0
20      CONTINUE
        FRNEXT (MHEAD) = MTAIL
        FRPREV (MTAIL) = MHEAD
        FRXP (MHEAD) = XHEAD
        FRXP (MTAIL) = XHEAD
        RINFO (2) = RINFO (2) + NZ
        FFLEFR = 0
        FFLEFC = 0
        FFSIZE = 0
        FFXP = XHEAD
C=======================================================================
C=======================================================================
        DO 600 S = 1, NLU
C=======================================================================
C=======================================================================
           LUIP = LUP (S)
           LUK = LUI (LUIP+1)
           LUDEGC = LUI (LUIP+3)
           LUDEGR = LUI (LUIP+2)
           LUNSON = LUI (LUIP+4)
           LUCP = (LUIP + 7)
           LURP = LUCP + LUDEGC
           LUSONP = LURP + LUDEGR
           LDIMC = LUK + LUDEGC
           LDIMR = LUK + LUDEGR
C=======================================================================
C=======================================================================
           IF (LUI (LUIP+6) .NE. 0) THEN
              DO 30 I = 1, FFLEFR
                 WIC (WPR (I)) = -1
30            CONTINUE
              DO 40 I = 1, FFLEFC
                 WIR (WPC (I)) = -1
40            CONTINUE
              XS = FFLEFR * FFLEFC
              IF (FFSIZE .NE. 0) THEN
                 XNEED = XNEED - (FFSIZE - XS)
                 XRUSE = XRUSE - (FFSIZE - XS)
                 INFO (13) = INFO (13) + 1
              ENDIF
              IF (FFLEFR .LE. 0 .OR. FFLEFC .LE. 0) THEN
                 XUSE = XUSE - (XHEAD - FRXP (MTAIL))
                 XHEAD = FRXP (MTAIL)
              ELSE
                 E = S - 1
                 RINFO (2) = RINFO (2) + XS
                 IF (XS .LE. XFREE) THEN
                    XFREE = XFREE - XS
                    IF (PFREE .EQ. MTAIL) THEN
                       PREV = FRPREV (MTAIL)
                       NEXT = MTAIL
                       XDP = FRXP (MTAIL)
                       FRXP (MTAIL) = XDP + XS
                    ELSE
                       PREV = PFREE
                       NEXT = FRNEXT (PFREE)
                       XDP = FRXP (NEXT) - XS
                       IF (XFREE .EQ. 0 .AND. PFREE .NE. MHEAD) THEN
                          PREV = FRPREV (PREV)
                          PFREE = 0
                          XFREE = -1
                       ENDIF
                    ENDIF
                    DO 60 J = 0, FFLEFR - 1
CFPP$ NODEPCHK L
                       DO 50 I = 0, FFLEFC - 1
                          XX (XDP+J*FFLEFC+I) = XX (FFXP+J*FFDIMC+I)
50                     CONTINUE
60                  CONTINUE
                    XUSE = XUSE - (XHEAD - FRXP (MTAIL))
                    XHEAD = FRXP (MTAIL)
                    FRXP (E) = XDP
                    FRDIMC (E) = FFLEFC
                 ELSE
                    XS = FFSIZE - (FFLEFC + (FFLEFR-1)*FFDIMC)
                    XHEAD = XHEAD - XS
                    XUSE = XUSE - XS
                    PREV = FRPREV (MTAIL)
                    NEXT = MTAIL
                    FRXP (MTAIL) = XHEAD
                    FRXP (E) = FFXP
                    FRDIMC (E) = FFDIMC
                 ENDIF
                 FRNEXT (PREV) = E
                 FRPREV (NEXT) = E
                 FRNEXT (E) = NEXT
                 FRPREV (E) = PREV
              ENDIF
              IF (PFREE .EQ. MTAIL) THEN
                 PFREE = 0
                 XFREE = -1
              ENDIF
              FFDIMC = LUI (LUIP+6)
              FFDIMR = LUI (LUIP+5)
              FFSIZE = FFDIMR * FFDIMC
              FFXP = 0
              IF (FFSIZE .GT. XTAIL-XHEAD) THEN
                 INFO (15) = INFO (15) + 1
                 CALL MA38OD (XX, XSIZE, XHEAD, XUSE,
     $              LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     $              FFXP, FFSIZE, PFREE, XFREE)
              ENDIF
              FFXP = XHEAD
              XHEAD = XHEAD + FFSIZE
              XUSE = XUSE + FFSIZE
              XNEED = XNEED + FFSIZE
              XRUSE = XRUSE + FFSIZE
              XRMAX = MAX (XRMAX, XRUSE)
              INFO (20) = MAX (INFO (20), XUSE)
              INFO (21) = MAX (INFO (21), XNEED)
              IF (XHEAD .GT. XTAIL) THEN
                 GO TO 9000
              ENDIF
              DO 70 P = FFXP, FFXP + FFSIZE - 1
                 XX (P) = ZERO
70            CONTINUE
              DO 80 K = 1, LUK
                 WIC (NPIV + K) = (LDIMR - K) * FFDIMC
                 WIR (NPIV + K) =  LDIMC - K
80            CONTINUE
              DO 90 I = 0, LUDEGR - 1
                 COL = LUI (LURP+I)
                 WIC (COL) = I * FFDIMC
                 WPR (I+1) = COL
90            CONTINUE
              DO 100 I = 0, LUDEGC - 1
                 ROW = LUI (LUCP+I)
                 WIR (ROW) = I
                 WPC (I+1) = ROW
100           CONTINUE
           ELSE
              DO 120 J = FFLEFR, LDIMR - 1
                 DO 110 I = 0, LDIMC - 1
                    XX (FFXP + J*FFDIMC + I) = ZERO
110              CONTINUE
120           CONTINUE
              DO 140 I = FFLEFC, LDIMC - 1
CFPP$ NODEPCHK L
                 DO 130 J = 0, FFLEFR - 1
                    XX (FFXP + J*FFDIMC + I) = ZERO
130              CONTINUE
140           CONTINUE
              DO 220 K = 1, LUK
                 PIVOT = NPIV + K
                 XSP = WIC (PIVOT)
                 KPOS = LDIMR - K + 1
                 XDP = (KPOS - 1) * FFDIMC
                 WIC (PIVOT) = XDP
                 IF (XSP .GE. 0) THEN
                    COLPOS = (XSP / FFDIMC) + 1
                    FSP = FFXP + XSP
                    FDP = FFXP + XDP
                    IF (FFLEFR .LT. KPOS) THEN
                       IF (FFLEFR .EQ. COLPOS) THEN
CFPP$ NODEPCHK L
                          DO 150 I = 0, LDIMC - 1
                             XX (FDP+I) = XX (FSP+I)
                             XX (FSP+I) = ZERO
150                       CONTINUE
                       ELSE
                          FLP = FFXP + (FFLEFR - 1) * FFDIMC
CFPP$ NODEPCHK L
                          DO 160 I = 0, LDIMC - 1
                             XX (FDP+I) = XX (FSP+I)
                             XX (FSP+I) = XX (FLP+I)
                             XX (FLP+I) = ZERO
160                       CONTINUE
                          SWPCOL = WPR (FFLEFR)
                          WPR (COLPOS) = SWPCOL
                          WIC (SWPCOL) = XSP
                       ENDIF
                    ELSE IF (COLPOS .NE. KPOS) THEN
CFPP$ NODEPCHK L
                       DO 180 I = 0, LDIMC - 1
                          X = XX (FDP+I)
                          XX (FDP+I) = XX (FSP+I)
                          XX (FSP+I) = X
180                    CONTINUE
                       SWPCOL = WPR (KPOS)
                       WPR (COLPOS) = SWPCOL
                       WIC (SWPCOL) = XSP
                    ENDIF
                    FFLEFR = FFLEFR - 1
                 ENDIF
                 XSP = WIR (PIVOT)
                 KPOS = LDIMC - K + 1
                 XDP = (KPOS - 1)
                 WIR (PIVOT) = XDP
                 IF (XSP .GE. 0) THEN
                    ROWPOS = XSP + 1
                    FSP = FFXP + XSP
                    FDP = FFXP + XDP
                    IF (FFLEFC .LT. KPOS) THEN
                       IF (FFLEFC .EQ. ROWPOS) THEN
CFPP$ NODEPCHK L
                          DO 190 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC
                             XX (FDP+J) = XX (FSP+J)
                             XX (FSP+J) = ZERO
190                       CONTINUE
                       ELSE
                          FLP = FFXP + (FFLEFC - 1)
CFPP$ NODEPCHK L
                          DO 200 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC
                             XX (FDP+J) = XX (FSP+J)
                             XX (FSP+J) = XX (FLP+J)
                             XX (FLP+J) = ZERO
200                       CONTINUE
                          SWPROW = WPC (FFLEFC)
                          WPC (ROWPOS) = SWPROW
                          WIR (SWPROW) = XSP
                       ENDIF
                    ELSE IF (ROWPOS .NE. KPOS) THEN
CFPP$ NODEPCHK L
                       DO 210 J = 0, (LDIMR - 1) * FFDIMC, FFDIMC
                          X = XX (FDP+J)
                          XX (FDP+J) = XX (FSP+J)
                          XX (FSP+J) = X
210                    CONTINUE
                       SWPROW = WPC (KPOS)
                       WPC (ROWPOS) = SWPROW
                       WIR (SWPROW) = XSP
                    ENDIF
                    FFLEFC = FFLEFC - 1
                 ENDIF
220           CONTINUE
              I = FFLEFR
              DO 230 P = LURP, LURP + LUDEGR - 1
                 COL = LUI (P)
                 IF (WIC (COL) .LT. 0) THEN
                    WIC (COL) = I * FFDIMC
                    I = I + 1
                    WPR (I) = COL
                 ENDIF
230           CONTINUE
              I = FFLEFC
              DO 240 P = LUCP, LUCP + LUDEGC - 1
                 ROW = LUI (P)
                 IF (WIR (ROW) .LT. 0) THEN
                    WIR (ROW) = I
                    I = I + 1
                    WPC (I) = ROW
                 ENDIF
240           CONTINUE
           ENDIF
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
           DO 260 K = 1, LUK
              I = NPIV + K
              XCDP = FFXP + WIC (I)
              XRDP = FFXP + WIR (I)
              DO 250 P = CP (I+1), CP (I) - 1
                 J = ARI (P)
                 IF (J .GT. 0) THEN
                    XP = XCDP + WIR (J)
                    IF (XP .LT. XCDP) THEN
                       NOUTSD = NOUTSD + 1
                       IF (PR3) THEN
                          ROW = RPERM (J)
                          COL = CPERM (I)
                          CALL MA38ZD (2, 97, ROW, COL, XX (P), IO)
                       ENDIF
                    ELSE
                       XX (XP) = XX (XP) + XX (P)
                    ENDIF
                 ELSE
                    XP = XRDP + WIC (-J)
                    IF (XP .LT. XRDP) THEN
                       NOUTSD = NOUTSD + 1
                       IF (PR3) THEN
                          ROW = RPERM (I)
                          COL = CPERM (-J)
                          CALL MA38ZD (2, 97, ROW, COL, XX (P), IO)
                       ENDIF
                    ELSE
                       XX (XP) = XX (XP) + XX (P)
                    ENDIF
                 ENDIF
250           CONTINUE
260        CONTINUE
           P = CP (NPIV + LUK + 1)
           XS = CP (NPIV + 1) - P
           FRXP (MHEAD) = P
           XNEED = XNEED - XS
           IF (XS .GT. XFREE) THEN
              XFREE = XS
              PFREE = MHEAD
           ENDIF
C=======================================================================
C=======================================================================
           DO 480 SP = LUSONP, LUSONP + LUNSON - 1
              E = LUI (SP)
              IF (E .LE. N) THEN
                 TYPE = 1
              ELSE IF (E .LE. 2*N) THEN
                 E = E - N
                 TYPE = 2
              ELSE
                 E = E - 2*N
                 TYPE = 3
              ENDIF
              FDIMC = FRDIMC (E)
              IF (FDIMC .NE. 0) THEN
                 FXP = FRXP (E)
                 FLUIP = LUP (E)
                 FDEGR = LUI (FLUIP+2)
                 FDEGC = LUI (FLUIP+3)
                 ALLCOL = FDEGR .GT. 0
                 ALLROW = FDEGC .GT. 0
                 FDEGR = ABS (FDEGR)
                 FDEGC = ABS (FDEGC)
                 FLUCP = (FLUIP + 7)
                 FLURP = FLUCP + FDEGC
                 IF (TYPE .EQ. 1) THEN
                    IF (ALLROW) THEN
                       DO 270 I = 0, FDEGC-1
                          ROW = LUI (FLUCP+I)
                          WM (I+1) = WIR (ROW)
270                    CONTINUE
                       IF (ALLCOL) THEN
                          DO 290 J = 0, FDEGR-1
                             COL = LUI (FLURP+J)
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 280 I = 0, FDEGC-1
                                XX (XDP + WM (I+1)) =
     $                          XX (XDP + WM (I+1)) +
     $                          XX (FXP + J*FDIMC + I)
280                          CONTINUE
290                       CONTINUE
                       ELSE
                          DO 310 J = 0, FDEGR-1
                             COL = LUI (FLURP+J)
                             IF (COL .GT. 0) THEN
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 300 I = 0, FDEGC-1
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
300                             CONTINUE
                             ENDIF
310                       CONTINUE
                       ENDIF
                    ELSE
                       DEGC = 0
                       DO 320 I = 0, FDEGC-1
                          ROW = LUI (FLUCP+I)
                          IF (ROW .GT. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
                          ENDIF
320                    CONTINUE
                       IF (ALLCOL) THEN
                          DO 340 J = 0, FDEGR-1
                             COL = LUI (FLURP+J)
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 330 I = 1, DEGC
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
330                          CONTINUE
340                       CONTINUE
                       ELSE
                          DO 360 J = 0, FDEGR-1
                             COL = LUI (FLURP+J)
                             IF (COL .GT. 0) THEN
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 350 I = 1, DEGC
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
350                             CONTINUE
                             ENDIF
360                       CONTINUE
                       ENDIF
                    ENDIF
                    FRDIMC (E) = 0
                    PREV = FRPREV (E)
                    NEXT = FRNEXT (E)
                    XNEED = XNEED - FDEGR*FDEGC
                    XRUSE = XRUSE - FDEGR*FDEGC
                    IF (FRDIMC (PREV) .LE. 0) THEN
                       FRNEXT (PREV) = NEXT
                       FRPREV (NEXT) = PREV
                       E = PREV
                       PREV = FRPREV (E)
                    ENDIF
                    IF (FRDIMC (NEXT) .LE. 0) THEN
                       FRXP (NEXT) = FRXP (E)
                       IF (E .LE. NLU) THEN
                          FRNEXT (PREV) = NEXT
                          FRPREV (NEXT) = PREV
                       ENDIF
                       E = NEXT
                       NEXT = FRNEXT (E)
                       IF (FRNEXT (MHEAD) .EQ. MTAIL) THEN
                          FRXP (MTAIL) = FRXP (MHEAD)
                       ENDIF
                    ENDIF
                    IF (NEXT .EQ. 0) THEN
                       XS = FFXP - FRXP (E)
                    ELSE
                       XS = FRXP (NEXT) - FRXP (E)
                    ENDIF
                    IF (XS .GT. XFREE) THEN
                       XFREE = XS
                       PFREE = E
                    ENDIF
                 ELSE IF (TYPE .EQ. 2) THEN
                    IF (ALLROW) THEN
                       DO 370 I = 0, FDEGC-1
                          ROW = LUI (FLUCP+I)
                          WM (I+1) = WIR (ROW)
370                    CONTINUE
                       DO 390 J = 0, FDEGR-1
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN
                             IF (WIC (COL) .GE. 0) THEN
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 380 I = 0, FDEGC-1
                                   XX (XDP + WM (I+1)) =
     $                             XX (XDP + WM (I+1)) +
     $                             XX (FXP + J*FDIMC + I)
380                             CONTINUE
                                LUI (FLURP+J) = -COL
                             ENDIF
                          ENDIF
390                    CONTINUE
                    ELSE
                       DEGC = 0
                       DO 400 I = 0, FDEGC-1
                          ROW = LUI (FLUCP+I)
                          IF (ROW .GT. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
                          ENDIF
400                    CONTINUE
                       DO 420 J = 0, FDEGR-1
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN
                             IF (WIC (COL) .GE. 0) THEN
                                XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                                DO 410 I = 1, DEGC
                                   XX (XDP + WM (I)) =
     $                             XX (XDP + WM (I)) +
     $                             XX (FXP + J*FDIMC + WJ (I))
410                             CONTINUE
                                LUI (FLURP+J) = -COL
                             ENDIF
                          ENDIF
420                    CONTINUE
                    ENDIF
                    LUI (FLUIP+2) = -FDEGR
                 ELSE
                    DEGC = 0
                    DO 430 I = 0, FDEGC-1
                       ROW = LUI (FLUCP+I)
                       IF (ROW .GT. 0) THEN
                          IF (WIR (ROW) .GE. 0) THEN
                             DEGC = DEGC + 1
                             WJ (DEGC) = I
                             WM (DEGC) = WIR (ROW)
                             LUI (FLUCP+I) = -ROW
                          ENDIF
                       ENDIF
430                 CONTINUE
                    IF (ALLCOL) THEN
                       DO 450 J = 0, FDEGR-1
                          COL = LUI (FLURP+J)
                          XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                          DO 440 I = 1, DEGC
                             XX (XDP + WM (I)) =
     $                       XX (XDP + WM (I)) +
     $                       XX (FXP + J*FDIMC + WJ (I))
440                       CONTINUE
450                    CONTINUE
                    ELSE
                       DO 470 J = 0, FDEGR-1
                          COL = LUI (FLURP+J)
                          IF (COL .GT. 0) THEN
                             XDP = FFXP + WIC (COL)
CFPP$ NODEPCHK L
                             DO 460 I = 1, DEGC
                                XX (XDP + WM (I)) =
     $                          XX (XDP + WM (I)) +
     $                          XX (FXP + J*FDIMC + WJ (I))
460                          CONTINUE
                          ENDIF
470                    CONTINUE
                    ENDIF
                    LUI (FLUIP+3) = -FDEGC
                 ENDIF
              ENDIF
480        CONTINUE
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
           K0 = 0
           FFLEFR = LDIMR
           FFLEFC = LDIMC
           FFCP = FFXP + FFLEFR * FFDIMC
           FFRP = FFXP + FFLEFC
           FFPP = FFXP + FFLEFC + FFLEFR * FFDIMC
           DO 500 K = 1, LUK
              IF (K-K0-2 .GT. 0) THEN
                 CALL DTRSV ('U', 'N', 'U', K-K0-1,
     $                         XX (FFPP         ), FFDIMC,
     $                         XX (FFPP - FFDIMC), 1)
                 RINFO (5) = RINFO (5) + (K-K0-2)*(K-K0-1)
              ENDIF
              IF (K-K0-1 .GT. 0) THEN
                 CALL DGEMV ('N', FFLEFC, K-K0-1,
     $                   -ONE, XX (FFCP         ), FFDIMC,
     $                         XX (FFPP - FFDIMC), 1,
     $                    ONE, XX (FFCP - FFDIMC), 1)
                 RINFO (5) = RINFO (5) + 2*FFLEFC*(K-K0-1)
              ENDIF
              FFCP = FFCP - FFDIMC
              FFRP = FFRP - 1
              FFPP = FFPP - FFDIMC - 1
              FFLEFR = FFLEFR - 1
              FFLEFC = FFLEFC - 1
              X = XX (FFPP)
              IF (ABS (X) .EQ. ZERO) THEN
                 GO TO 9010
              ENDIF
              X = ONE / X
              DO 490 P = FFCP, FFCP + FFLEFC - 1
                 XX (P) = XX (P) * X
490           CONTINUE
              RINFO (4) = RINFO (4) + FFLEFC
              INFO (17) = INFO (17) + 1
              IF (K-K0 .GE. NB .OR. K .EQ. LUK) THEN
                 CALL DTRSM ('L', 'U', 'N', 'U', K-K0, FFLEFR, ONE,
     $                      XX (FFPP), FFDIMC,
     $                      XX (FFRP), FFDIMC)
                 CALL DGEMM ('N', 'N', FFLEFC, FFLEFR, K-K0,
     $                -ONE, XX (FFCP ), FFDIMC,
     $                      XX (FFRP ), FFDIMC,
     $                 ONE, XX (FFXP), FFDIMC)
                 RINFO (6) = RINFO (6) + FFLEFR*(K-K0-1)*(K-K0)
     $                                         + 2*FFLEFC*FFLEFR*(K-K0)
                 K0 = K
              ENDIF
500        CONTINUE
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
           XS = LUK*LUDEGC + LUK*LUDEGR + LUK*LUK
           IF (XS .GT. XTAIL-XHEAD) THEN
              INFO (15) = INFO (15) + 1
              CALL MA38OD (XX, XSIZE, XHEAD, XUSE,
     $              LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     $              FFXP, FFSIZE, PFREE, XFREE)
           ENDIF
           XTAIL = XTAIL - XS
           LUXP = XTAIL
           XUSE = XUSE + XS
           XNEED = XNEED + XS
           XRUSE = XRUSE + XS
           XRMAX = MAX (XRMAX, XRUSE)
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
           IF (XHEAD .GT. XTAIL) THEN
              GO TO 9000
           ENDIF
           LUI (LUIP) = LUXP
           DO 510 I = 0, LUDEGC-1
              LUI (LUCP+I) = WPC (I+1)
510        CONTINUE
           DO 520 I = 0, LUDEGR-1
              LUI (LURP+I) = WPR (I+1)
520        CONTINUE
           XP = FFXP + (LDIMR-1)*FFDIMC + LDIMC-1
           DO 540 J = 0, LUK-1
CFPP$ NODEPCHK L
              DO 530 I = 0, LUK-1
                 XX (LUXP + J*LDIMC + I) = XX (XP - J*FFDIMC - I)
530           CONTINUE
540        CONTINUE
           IF (LUDEGC .NE. 0) THEN
              LXP = LUXP + LUK
              XP = FFXP + (LDIMR-1)*FFDIMC
              DO 560 J = 0, LUK-1
CFPP$ NODEPCHK L
                 DO 550 I = 0, LUDEGC-1
                    XX (LXP + J*LDIMC + I) = XX (XP - J*FFDIMC + I)
550              CONTINUE
560           CONTINUE
           ENDIF
           IF (LUDEGR .NE. 0) THEN
              UXP = LUXP + LUK * LDIMC
              XP = FFXP + LDIMC-1
              DO 580 J = 0, LUDEGR-1
CFPP$ NODEPCHK L
                 DO 570 I = 0, LUK-1
                    XX (UXP + J*LUK + I) = XX (XP + J*FFDIMC - I)
570              CONTINUE
580           CONTINUE
           ENDIF
           NZU = (LUK*(LUK-1)/2) + LUK*LUDEGC
           NZL = (LUK*(LUK-1)/2) + LUK*LUDEGR
           INFO (10) = INFO (10) + NZL
           INFO (11) = INFO (11) + NZU
           DO 590 PIVOT = NPIV + 1, NPIV + LUK
              WIR (PIVOT) = -1
              WIC (PIVOT) = -1
590        CONTINUE
           NPIV = NPIV + LUK
C=======================================================================
C=======================================================================
600     CONTINUE
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        IF (NPIV .LT. N) THEN
           IF (PR3) THEN
              DO 620 I = NPIV+1, N
                 DO 610 P = CP (I+1), CP (I) - 1
                    J = ARI (P)
                    IF (J .GT. 0) THEN
                       ROW = RPERM (J)
                       COL = CPERM (I)
                    ELSE
                       ROW = RPERM (I)
                       COL = CPERM (-J)
                    ENDIF
                    CALL MA38ZD (2, 95, ROW, COL, XX(P), IO)
610              CONTINUE
620           CONTINUE
           ENDIF
           NOUTSD = NOUTSD + (CP (NPIV+1) - CP (N+1))
        ENDIF
        IF (FFSIZE .NE. 0) THEN
           INFO (13) = INFO (13) + 1
        ENDIF
        XUSE = XUSE - (XHEAD - CP (N+1))
        XNEED = XUSE
        XHEAD = CP (N+1)
        IF (NLU .EQ. 0) THEN
           XTAIL = XSIZE
           XUSE = XUSE + 1
           XRUSE = XUSE
           XNEED = XUSE
           INFO (20) = MAX (INFO (20), XUSE)
           INFO (21) = MAX (INFO (21), XNEED)
        ENDIF
        IF (XHEAD .LE. XTAIL) THEN
           IF (NLU .EQ. 0) THEN
              XX (XTAIL) = ZERO
           ENDIF
           DO 630 S = 1, NLU
              LUIP = LUP (S)
              LUXP = LUI (LUIP)
              LUI (LUIP) = LUXP - XTAIL + 1
630        CONTINUE
           XRUSE = XUSE
           XRMAX = MAX (XRMAX, XRUSE)
           RETURN
        ENDIF
C=======================================================================
C=======================================================================
9000    CONTINUE
        CALL MA38ND (2, ICNTL, INFO, -4, INFO (21))
        RETURN
9010    CONTINUE
        CALL MA38ND (2, ICNTL, INFO, -6, 0)
        RETURN
        END
        SUBROUTINE MA38SD (NLU, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, N, LUP(NLU), LUI(*)
        DOUBLE PRECISION LUX(*), X(N), W(N)
C=== MA38SD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        EXTERNAL  DTRSV, DGEMV
C=======================================================================
C=======================================================================
        INTEGER I, K, S, LUIP, LUXP, LUK, LUDEGR, LUDEGC, LURP, UXP,
     $          LUCP, ROW
        DOUBLE PRECISION
     $          ONE
        PARAMETER (ONE = 1.0D0)
C=======================================================================
C=======================================================================
        K = 0
        DO 40 S = 1, NLU
           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGR = LUI (LUIP+2)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LURP   = LUCP + LUDEGC
           IF (LUK .EQ. 1) THEN
              K = K + 1
              X (K) = X (K) / LUX (LUXP)
              UXP = LUXP + LUDEGC + 1
CFPP$ NODEPCHK L
              DO 10 I = 1, LUDEGR
                 ROW = LUI (LURP+I-1)
                 X (ROW) = X (ROW) - LUX (UXP+I-1) * X (K)
10            CONTINUE
           ELSE
              UXP = LUXP + LUK * (LUDEGC + LUK)
              CALL DTRSV ('U', 'T', 'N', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)
              DO 20 I = 1, LUDEGR
                 ROW = LUI (LURP+I-1)
                 W (I) = X (ROW)
20            CONTINUE
              CALL DGEMV ('T', LUK, LUDEGR, -ONE,
     $           LUX (UXP), LUK, X (K+1), 1, ONE, W, 1)
              DO 30 I = 1, LUDEGR
                 ROW = LUI (LURP+I-1)
                 X (ROW) = W (I)
30            CONTINUE
              K = K + LUK
           ENDIF
40      CONTINUE
        RETURN
        END
        SUBROUTINE MA38TD (NLU, NPIV, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, NPIV, N, LUP(NLU), LUI(*)
        DOUBLE PRECISION LUX(*), X(N), W(N)
C=== MA38TD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        EXTERNAL DTRSV, DGEMV
C=======================================================================
C=======================================================================
        INTEGER J, K, S, LUIP, LUXP, LUK, LUDEGC, LUCP, LXP, COL
        DOUBLE PRECISION
     $          ONE
        PARAMETER (ONE = 1.0D0)
C=======================================================================
C=======================================================================
        K = NPIV
        DO 30 S = NLU, 1, -1
           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LXP    = LUXP + LUK
           IF (LUK .EQ. 1) THEN
CFPP$ NODEPCHK L
              DO 10 J = 1, LUDEGC
                 COL = LUI (LUCP+J-1)
                 X (K) = X (K) - LUX (LXP+J-1) * X (COL)
10            CONTINUE
              K = K - 1
           ELSE
              K = K - LUK
              DO 20 J = 1, LUDEGC
                 COL = LUI (LUCP+J-1)
                 W (J) = X (COL)
20            CONTINUE
              CALL DGEMV ('T', LUDEGC, LUK, -ONE,
     $           LUX (LXP), LUDEGC + LUK, W, 1, ONE, X (K+1), 1)
              CALL DTRSV ('L', 'T', 'U', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)
           ENDIF
30      CONTINUE
        RETURN
        END
        SUBROUTINE MA38UD (NLU, NPIV, N, LUP, LUI, LUX, X, W)
        INTEGER NLU, NPIV, N, LUP(NLU), LUI(*)
        DOUBLE PRECISION LUX(*), X(N), W(N)
C=== MA38UD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        EXTERNAL DTRSV, DGEMV
C=======================================================================
C=======================================================================
        INTEGER J, K, S, LUIP, LUXP, LUK, LUDEGR, LUDEGC, LURP, UXP,
     $          LUCP, COL
        DOUBLE PRECISION
     $          ONE
        PARAMETER (ONE = 1.0D0)
C=======================================================================
C=======================================================================
        K = NPIV
        DO 30 S = NLU, 1, -1
           LUIP   = LUP (S)
           LUXP   = LUI (LUIP)
           LUK    = LUI (LUIP+1)
           LUDEGR = LUI (LUIP+2)
           LUDEGC = LUI (LUIP+3)
           LUCP   = (LUIP + 7)
           LURP   = LUCP + LUDEGC
           UXP    = LUXP + LUK * (LUDEGC + LUK)
           IF (LUK .EQ. 1) THEN
CFPP$ NODEPCHK L
              DO 10 J = 1, LUDEGR
                 COL = LUI (LURP+J-1)
                 X (K) = X (K) - LUX (UXP+J-1) * X (COL)
10            CONTINUE
              X (K) = X (K) / LUX (LUXP)
              K = K - 1
           ELSE
              K = K - LUK
              DO 20 J = 1, LUDEGR
                 COL = LUI (LURP+J-1)
                 W (J) = X (COL)
20            CONTINUE
              CALL DGEMV ('N', LUK, LUDEGR, -ONE,
     $           LUX (UXP), LUK, W, 1, ONE, X (K+1), 1)
              CALL DTRSV ('U', 'N', 'N', LUK,
     $           LUX (LUXP), LUDEGC + LUK, X (K+1), 1)
           ENDIF
30      CONTINUE
        RETURN
        END
        SUBROUTINE MA38YD (WHO, WHERE,
     *          N, NE, JOB, TRANS, LVALUE, LINDEX, VALUE,
     *          INDEX, KEEP, CNTL, ICNTL, INFO, RINFO,
     *          B, X, LX, W, LW)
        INTEGER WHO, WHERE, N, NE, JOB, LVALUE, LINDEX, INDEX(LINDEX),
     *          KEEP(20), ICNTL(20), INFO(40), LX, LW
        DOUBLE PRECISION VALUE(LVALUE), CNTL(10), RINFO(20), B(LX),
     *          X(LX), W(LW)
        LOGICAL TRANS
C=== MA38YD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        INTRINSIC MIN
C=======================================================================
C=======================================================================
        LOGICAL TRANSA, TRANSC, PRLU, BADLU, SGLTON, PRESRV, SYMBOL
        INTEGER IO, PRL, PRN, K, LUI1, LUI2, LUX1, LUX2, ROW, COL,
     $          FACNE, FACN, NZ, FACJOB, NBLKS, NZOFF, FACTRA, CPERMP,
     $          RPERMP, APP, AXP, AIP, OFFIP, OFFXP, LUBLPP, OFFPP,
     $          BLKPP, P1, P2, P, BLK, K1, K2, KN, LUIIP, LUXXP, NPIV,
     $          NLU, E, LUK, LUPP, LUIP, LUXP, LUDEGR, LUDEGC, LUNSON,
     $          LUSONP, LUCP, LURP, I, J, NZCOL, NZROW, UXP, SON,
     $          PRMAX, LUDIMR, LUDIMC, MAXDR, MAXDC, LUIR1, IP1, IP2,
     $          XP1
        DOUBLE PRECISION
     $          ONE
        PARAMETER (PRMAX = 10, ONE = 1.0D0)
C=======================================================================
C=======================================================================
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IO  = ICNTL(2)
        PRL = ICNTL(3)
        IF (PRL .LT. 3 .OR. IO .LT. 0) THEN
           RETURN
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (WHO .EQ. 1) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'MA38AD input:'
              PRLU = .FALSE.
           ELSE
              WRITE (IO, 6) 'MA38AD output:'
              PRLU = .TRUE.
           ENDIF
        ELSE IF (WHO .EQ. 2) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'MA38BD input:'
              PRLU = .TRUE.
           ELSE
              WRITE (IO, 6) 'MA38BD output:'
              PRLU = .TRUE.
           ENDIF
        ELSE IF (WHO .EQ. 3) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'MA38CD input:'
              PRLU = .TRUE.
           ELSE
              WRITE (IO, 6) 'MA38CD output:'
              PRLU = .FALSE.
           ENDIF
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (WHERE .EQ. 1) THEN
           WRITE (IO, 1)  'Scalar arguments:'
           WRITE (IO, 1)  '   N:         ', N, ' : order of matrix A'
           IF (WHO .EQ. 3) THEN
              LUI2 = KEEP (5)
              TRANSA = .FALSE.
              IF (LUI2-6 .GE. 1 .AND. LUI2-6 .LE. LINDEX) THEN
                 TRANSA = INDEX (LUI2-6) .NE. 0
              ENDIF
              TRANSC = TRANS
              IF (.NOT. TRANSC) THEN
                 IF (JOB .EQ. 1) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve P''Lx=b'
                 ELSE IF (JOB .EQ. 2) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve UQ''x=b'
                 ELSE IF (.NOT. TRANSA) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve Ax=b (PAQ=LU was factorized)'
                 ELSE
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve A''x=b (PA''Q=LU was factorized)'
                 ENDIF
              ELSE
                 IF (JOB .EQ. 1) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve L''Px=b'
                 ELSE IF (JOB .EQ. 2) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve QU''x=b'
                 ELSE IF (.NOT. TRANSA) THEN
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve A''x=b (PAQ=LU was factorized)'
                 ELSE
                    WRITE (IO, 1) '   JOB:       ', JOB,
     $              ' : solve Ax=b (PA''Q=LU was factorized)'
                 ENDIF
              ENDIF
              IF (TRANSC) THEN
                 WRITE (IO, 1)
     $           '   TRANSC:          .true. : see JOB above '
              ELSE
                 WRITE (IO, 1)
     $           '   TRANSC:         .false. : see JOB above '
              ENDIF
           ELSE
              WRITE (IO, 1) '   NE:        ', NE,
     $        ' : entries in matrix A'
              IF (JOB .EQ. 1) THEN
                 WRITE (IO, 1) '   JOB:       ', JOB,
     $           ' : matrix A preserved'
              ELSE
                 WRITE (IO, 1) '   JOB:       ', JOB,
     $           ' : matrix A not preserved'
              ENDIF
              TRANSA = TRANS
              IF (TRANSA) THEN
                 WRITE (IO, 1)
     $           '   TRANSA:          .true. : factorize A transpose'
              ELSE
                 WRITE (IO, 1)
     $           '   TRANSA:         .false. : factorize A'
              ENDIF
           ENDIF
           WRITE (IO, 1) '   LVALUE:    ',LVALUE,
     $     ' : size of VALUE array'
           WRITE (IO, 1) '   LINDEX:    ',LINDEX,
     $     ' : size of INDEX array'
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (WHERE .EQ. 1) THEN
           WRITE (IO, 1)
     $     'Control parameters, normally initialized by MA38ID:'
           WRITE (IO, 1) '   ICNTL (1): ', ICNTL (1),
     $     ' : I/O unit for error and warning messages'
           WRITE (IO, 1) '   ICNTL (2): ', IO,
     $     ' : I/O unit for diagnostics'
           WRITE (IO, 1) '   ICNTL (3): ', PRL,
     $     ' : printing control'
           IF (WHO .EQ. 1) THEN
              IF (ICNTL (4) .EQ. 1) THEN
                 WRITE (IO, 1) '   ICNTL (4): ', ICNTL (4),
     $           ' : use block triangular form (BTF)'
              ELSE
                 WRITE (IO, 1) '   ICNTL (4): ', ICNTL (4),
     $           ' : do not permute to block triangular form (BTF)'
              ENDIF
              WRITE (IO, 1) '   ICNTL (5): ', ICNTL (5),
     $        ' : columns examined during pivot search'
              IF (ICNTL (6) .NE. 0) THEN
                 WRITE (IO, 1) '   ICNTL (6): ', ICNTL (6),
     $           ' : preserve symmetry'
              ELSE
                 WRITE (IO, 1) '   ICNTL (6): ', ICNTL (6),
     $           ' : do not preserve symmetry'
              ENDIF
           ENDIF
           IF (WHO .NE. 3) THEN
              WRITE (IO, 1) '   ICNTL (7): ', ICNTL (7),
     $        ' : block size for dense matrix multiply'
           ELSE
              WRITE (IO, 1) '   ICNTL (8): ', ICNTL (8),
     $        ' : maximum number of iterative refinement steps'
           ENDIF
           IF (WHO .EQ. 1) THEN
              WRITE (IO, 3) '   CNTL (1): ',CNTL (1),
     $        ' : relative pivot tolerance'
              WRITE (IO, 3) '   CNTL (2): ',CNTL (2),
     $        ' : frontal matrix growth factor'
              WRITE (IO, 1) '   KEEP (7):  ',KEEP(7),
     $        ' : dense row/col control, d1'
              WRITE (IO, 1) '   KEEP (8):  ',KEEP(8),
     $        ' : dense row/col control, d2'
           ELSE IF (WHO .EQ. 3) THEN
              WRITE (IO, 3) '   CNTL (3): ',CNTL(3),
     $        ' : machine epsilon'
           ENDIF
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (WHERE .NE. 1) THEN
           WRITE (IO, 1) 'Output information:'
           IF (INFO (1) .LT. 0) THEN
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : error occurred.'
           ELSE IF (INFO (1) .GT. 0) THEN
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : warning occurred'
           ELSE
              WRITE (IO, 1) '   INFO (1):  ', INFO (1),
     $        ' : no error or warning occurred'
           ENDIF
           IF (WHO .NE. 3) THEN
              WRITE (IO, 1) '   INFO (2):  ', INFO (2),
     $        ' : duplicate entries in A'
              WRITE (IO, 1) '   INFO (3):  ', INFO (3),
     $        ' : invalid entries in A (indices not in 1..N)'
              WRITE (IO, 1) '   INFO (4):  ', INFO (4),
     $        ' : invalid entries in A (not in prior pattern)'
              WRITE (IO, 1) '   INFO (5):  ', INFO (5),
     $        ' : entries in A after summing duplicates'
              WRITE (IO, 1)
     $  '                             and removing invalid entries'
              WRITE (IO, 1) '   INFO (6):  ', INFO (6),
     $        ' : entries in diagonal blocks of A'
              WRITE (IO, 1) '   INFO (7):  ', INFO (7),
     $        ' : entries in off-diagonal blocks of A'
              WRITE (IO, 1) '   INFO (8):  ', INFO (8),
     $        ' : 1-by-1 diagonal blocks in A'
              WRITE (IO, 1) '   INFO (9):  ', INFO (9),
     $        ' : diagonal blocks in A (>1 only if BTF used)'
              WRITE (IO, 1) '   INFO (10): ', INFO (10),
     $        ' : entries below diagonal in L'
              WRITE (IO, 1) '   INFO (11): ', INFO (11),
     $        ' : entries above diagonal in U'
              WRITE (IO, 1) '   INFO (12): ', INFO (12),
     $        ' : entries in L + U + offdiagonal blocks of A'
              WRITE (IO, 1) '   INFO (13): ', INFO (13),
     $        ' : frontal matrices'
              WRITE (IO, 1) '   INFO (14): ', INFO (14),
     $        ' : integer garbage collections'
              WRITE (IO, 1) '   INFO (15): ', INFO (15),
     $        ' : real garbage collections'
              WRITE (IO, 1) '   INFO (16): ', INFO (16),
     $        ' : diagonal pivots chosen'
              WRITE (IO, 1) '   INFO (17): ', INFO (17),
     $        ' : numerically valid pivots found in A'
              WRITE (IO, 1) '   INFO (18): ', INFO (18),
     $        ' : memory used in INDEX'
              WRITE (IO, 1) '   INFO (19): ', INFO (19),
     $        ' : minimum memory needed in INDEX'
              WRITE (IO, 1) '   INFO (20): ', INFO (20),
     $        ' : memory used in VALUE'
              WRITE (IO, 1) '   INFO (21): ', INFO (21),
     $        ' : minimum memory needed in VALUE'
              WRITE (IO, 1) '   INFO (22): ', INFO (22),
     $        ' : memory needed in INDEX for next call to MA38BD'
              WRITE (IO, 1) '   INFO (23): ', INFO (23),
     $        ' : memory needed in VALUE for next call to MA38BD'
           ELSE
              WRITE (IO, 1) '   INFO (24): ', INFO (24),
     $        ' : steps of iterative refinement taken'
           ENDIF
           IF (WHO .NE. 3) THEN
              WRITE (IO, 3) '   RINFO (1):', RINFO (1),
     $        ' : total BLAS flop count'
              WRITE (IO, 3) '   RINFO (2):', RINFO (2),
     $        ' : assembly flop count'
              WRITE (IO, 3) '   RINFO (3):', RINFO (3),
     $        ' : pivot search flop count'
              WRITE (IO, 3) '   RINFO (4):', RINFO (4),
     $        ' : Level-1 BLAS flop count'
              WRITE (IO, 3) '   RINFO (5):', RINFO (5),
     $        ' : Level-2 BLAS flop count'
              WRITE (IO, 3) '   RINFO (6):', RINFO (6),
     $        ' : Level-3 BLAS flop count'
           ELSE IF (LW .EQ. 4*N) THEN
              WRITE (IO, 3) '   RINFO (7):', RINFO (7),
     $        ' : sparse error estimate omega1'
              WRITE (IO, 3) '   RINFO (8):', RINFO (8),
     $        ' : sparse error estimate omega2'
           ENDIF
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (WHERE .EQ. 1 .AND. WHO .NE. 3) THEN
           IF (TRANSA) THEN
              IF (PRL .GE. 5) THEN
                 WRITE (IO, 1) 'The input matrix A transpose:'
                 WRITE (IO, 1) '   VALUE (1 ... ',NE,
     $           ' ): numerical values'
                 WRITE (IO, 1) '   INDEX (1 ... ',NE,
     $           ' ): column indices'
                 WRITE (IO, 1) '   INDEX (',NE+1,' ... ',2*NE,
     $           ' ): row indices'
              ENDIF
              WRITE (IO, 1)
     $        'Input matrix A transpose (entry: row, column, value):'
           ELSE
              IF (PRL .GE. 5) THEN
                 WRITE (IO, 1) 'The input matrix A:'
                 WRITE (IO, 1) '   VALUE (1 ... ',NE,
     $           ' ): numerical values'
                 WRITE (IO, 1) '   INDEX (1 ... ',NE,
     $           ' ): row indices'
                 WRITE (IO, 1) '   INDEX (',NE+1,' ... ',2*NE,
     $           ' ): column indices'
              ENDIF
              WRITE (IO, 1)
     $        'Input matrix A (entry: row, column, value):'
           ENDIF
           PRN = MIN (PRMAX, NE)
           IF (PRL .GE. 4) THEN
              PRN = NE
           ENDIF
           DO 20 K = 1, PRN
              IF (TRANSA) THEN
                 ROW = INDEX (K+NE)
                 COL = INDEX (K)
              ELSE
                 ROW = INDEX (K)
                 COL = INDEX (K+NE)
              ENDIF
              WRITE (IO, 2) K, ROW, COL, VALUE (K)
20         CONTINUE
           IF (PRN .LT. NE) THEN
              WRITE (IO, 7)
           ENDIF
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (PRLU .AND. INFO (1) .LT. 0) THEN
           WRITE (IO, 1)
     $     'LU factors not printed because of error flag, INFO (1) ='
     $     , INFO (1)
           PRLU = .FALSE.
        ENDIF
        IF (PRLU) THEN
           LUX1 = KEEP (1)
           LUX2 = KEEP (2)
           LUI1 = KEEP (3)
           LUIR1 = KEEP (4)
           LUI2 = KEEP (5)
           XP1 = LUX1
           IP1 = LUI1
           IP2 = LUI2
           SYMBOL = WHO .EQ. 2 .AND. WHERE .EQ. 1
           IF (PRL .GE. 5) THEN
              IF (SYMBOL) THEN
                 WRITE (IO, 1)
     $           'KEEP (4...5) gives the location of LU factors'
                 WRITE (IO, 1)
     $           '   which must be preserved for calls to MA38BD: '
              ELSE
                 WRITE (IO, 1)
     $           'KEEP (1...5) gives the location of LU factors'
                 WRITE (IO, 1)
     $           '   which must be preserved for calls to MA38CD: '
                 WRITE (IO, 1) '      VALUE ( KEEP (1): ', LUX1,
     $           ' ... KEEP (2): ', LUX2,' )'
                 WRITE (IO, 1) '      INDEX ( KEEP (3): ', LUI1,
     $           ' ... KEEP (5): ', LUI2,' )'
                 WRITE (IO, 1) '   and for calls to MA38BD: '
              ENDIF
              WRITE (IO, 1) '      INDEX ( KEEP (4): ',LUIR1,
     $        ' ... KEEP (5): ', LUI2,' )'
           ENDIF
           BADLU = LUIR1 .LE. 0 .OR. LUI2-6 .LT. LUIR1 .OR.
     $        LUI2 .GT. LINDEX
           IF (.NOT. SYMBOL) THEN
              BADLU = BADLU .OR. LUX1 .LE. 0 .OR.
     $        LUX1 .GT. LUX2 .OR. LUX2 .GT. LVALUE .OR. LUI1 .LE. 0 .OR.
     $        LUIR1 .LT. LUI1 .OR. LUIR1 .GT. LUI2
           ENDIF
           IF (BADLU) THEN
              FACNE  = 0
              FACN   = 0
              NZ     = 0
              FACJOB = 0
              NBLKS  = 0
              NZOFF  = 0
              FACTRA = 0
           ELSE
              FACNE  = INDEX (LUI2)
              FACN   = INDEX (LUI2-1)
              NZ     = INDEX (LUI2-2)
              FACJOB = INDEX (LUI2-3)
              NBLKS  = INDEX (LUI2-4)
              NZOFF  = INDEX (LUI2-5)
              FACTRA = INDEX (LUI2-6)
           ENDIF
           PRESRV = FACJOB .NE. 0
           TRANSA = FACTRA .NE. 0
           RPERMP = (LUI2-6) - (FACN)
           CPERMP = RPERMP - (FACN)
           IP2 = CPERMP - 1
           IF (PRL .GE. 5) THEN
              IF (SYMBOL) THEN
                 WRITE (IO, 1) 'Layout of LU factors in INDEX:'
              ELSE
                 WRITE (IO, 1)
     $           'Layout of LU factors in VALUE and INDEX:'
              ENDIF
           ENDIF
           IF (PRESRV) THEN
              APP = IP1
              AIP = APP + (FACN+1)
              IP1 = AIP + (NZ)
              AXP = XP1
              XP1 = XP1 + (NZ)
              IF (PRL .GE. 5 .AND. .NOT. SYMBOL) THEN
                 WRITE (IO, 1)'   preserved copy of original matrix:'
                 WRITE (IO, 1)'      INDEX ( ',APP,' ... ', AIP-1,
     $           ' ): column pointers'
                 WRITE (IO, 1)'      INDEX ( ',AIP,' ... ', IP1-1,
     $           ' ): row indices'
                 WRITE (IO, 1)'      VALUE ( ',AXP,' ... ', XP1-1,
     $           ' ): numerical values'
              ENDIF
           ELSE
              IF (PRL .GE. 5 .AND. .NOT. SYMBOL) THEN
                 WRITE (IO, 1) '   original matrix not preserved.'
              ENDIF
           ENDIF
           BADLU = BADLU .OR.
     $          N .NE. FACN .OR. NZ .LE. 0 .OR. LUIR1 .GT. IP2 .OR.
     $          NBLKS .LE. 0 .OR. NBLKS .GT. N
           IF (.NOT. SYMBOL) THEN
              BADLU = BADLU .OR. XP1 .GT. LUX2 .OR. NZOFF .LT. 0
           ENDIF
           IF (BADLU) THEN
              NBLKS = 0
           ENDIF
           IF (NBLKS .LE. 1) THEN
              IF (PRL .GE. 5) THEN
                 WRITE (IO, 1)
     $           '   collection of elements in LU factors:'
                 WRITE (IO, 1) '      INDEX ( ',LUIR1,' ... ', IP2,
     $           ' ): integer data'
                 IF (.NOT. SYMBOL) THEN
                    WRITE (IO, 1) '      VALUE ( ',XP1,' ... ', LUX2,
     $              ' ): numerical values'
                 ENDIF
              ENDIF
           ELSE
              OFFIP = IP1
              IP1 = IP1 + (NZOFF)
              OFFXP = XP1
              XP1 = XP1 + (NZOFF)
              OFFPP = CPERMP - (N+1)
              BLKPP = OFFPP - (NBLKS+1)
              LUBLPP = BLKPP - (NBLKS)
              IP2 = LUBLPP - 1
              BADLU = BADLU .OR. LUIR1 .GT. IP2
              IF (.NOT. SYMBOL) THEN
                 BADLU = BADLU .OR. IP1 .GT. IP2 .OR.
     $           XP1 .GT. LUX2 .OR. LUIR1 .NE. IP1
              ENDIF
              IF (PRL .GE. 5) THEN
                 WRITE (IO, 1)
     $           '   matrix permuted to upper block triangular form.'
                 IF (NZOFF .NE. 0 .AND. .NOT. SYMBOL) THEN
                    WRITE (IO, 1)'   entries not in diagonal blocks:'
                    WRITE (IO, 1)'      INDEX ( ',OFFIP,' ... ',
     $              LUIR1-1, ' ): row indices'
                    WRITE (IO, 1)'      VALUE ( ',OFFXP,' ... ',
     $              XP1-1, ' ): numerical values'
                 ENDIF
                 WRITE (IO, 1)
     $  '   collection of elements in LU factors of diagonal blocks:'
                 IF (LUIR1 .LE. LUBLPP-1) THEN
                    WRITE (IO, 1) '      INDEX ( ',LUIR1,' ... ',
     $              IP2, ' ): integer data'
                 ENDIF
                 IF (XP1 .LE. LUX2 .AND. .NOT. SYMBOL) THEN
                    WRITE (IO, 1) '      VALUE ( ',XP1,' ... ', LUX2,
     $              ' ): numerical values'
                 ENDIF
                 WRITE (IO, 1) '   other block triangular data:'
                 WRITE (IO, 1) '      INDEX ( ',LUBLPP,' ... ',
     $           BLKPP-1, ' ): pointers to block factors'
                 WRITE (IO, 1) '      INDEX ( ', BLKPP,' ... ',
     $           OFFPP-1, ' ): index range of blocks'
                 IF (.NOT. SYMBOL) THEN
                    WRITE (IO, 1) '      INDEX ( ', OFFPP,' ... ',
     $              LUI2-7,' ): off-diagonal row pointers'
                 ENDIF
              ENDIF
           ENDIF
           IF (PRL .GE. 5) THEN
              WRITE (IO, 1)
     $        '   permutation vectors (start at KEEP(4)-2*N-6):'
              WRITE (IO, 1) '      INDEX ( ',CPERMP,' ... ',RPERMP-1,
     $        ' ): column permutations'
              WRITE (IO, 1) '      INDEX ( ',RPERMP,' ... ',LUI2-7,
     $        ' ): row permutations'
              WRITE (IO, 1) '   other data in INDEX: '
              WRITE (IO, 1) '      INDEX ( ',LUI2-6,' ): ', FACTRA,
     $        ' : TRANSA MA38AD/MA38BD argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2-5,' ): ', NZOFF,
     $        ' : entries in off-diagonal part'
              WRITE (IO, 1) '      INDEX ( ',LUI2-4,' ): ', NBLKS,
     $        ' : number of diagonal blocks'
              WRITE (IO, 1) '      INDEX ( ',LUI2-3,' ): ', FACJOB,
     $        ' : JOB MA38AD/MA38BD argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2-2,' ): ', NZ,
     $        ' : entries in original matrix'
              WRITE (IO, 1) '      INDEX ( ',LUI2-1,' ): ', FACN,
     $        ' : N MA38AD/MA38BD argument'
              WRITE (IO, 1) '      INDEX ( ',LUI2  ,' ): ', FACNE,
     $        ' : NE MA38AD/MA38BD argument'
           ENDIF
           IF (.NOT. SYMBOL) THEN
              BADLU = BADLU .OR. IP1 .NE. LUIR1
           ENDIF
           IP1 = LUIR1
           IF (BADLU) THEN
              WRITE (IO, 1) 'LU factors uncomputed or corrupted.'
              PRESRV = .FALSE.
              NBLKS = 0
           ENDIF
           IF (PRESRV .AND. .NOT. SYMBOL) THEN
              WRITE (IO, 8)
              WRITE (IO, 1)
     $        'Preserved copy of original matrix (stored by column),'
              WRITE (IO, 1) 'one entry per line (row index, value):'
              DO 40 COL = 1, N
                 P1 = INDEX (APP-1 + COL)
                 P2 = INDEX (APP-1 + COL+1) - 1
                 WRITE (IO, 1) '   col: ', COL
                 IF (PRL .EQ. 3) THEN
                    P2 = MIN (PRMAX, P2)
                 ENDIF
                 DO 30 P = P1, P2
                    WRITE (IO, 5) INDEX (AIP-1 + P), VALUE (AXP-1 + P)
30               CONTINUE
                 IF (PRL .EQ. 3 .AND. P2 .GE. PRMAX) THEN
                    WRITE (IO, 7)
                    GO TO 50
                 ENDIF
40            CONTINUE
50            CONTINUE
           ENDIF
           IF (NBLKS .GT. 1 .AND. .NOT. SYMBOL) THEN
              WRITE (IO, 8)
              WRITE (IO, 1)
     $        'Entries not in diagonal blocks (stored by row):'
              WRITE (IO, 1) 'one entry per line (column index, value):'
              IF (NZOFF .EQ. 0) THEN
                 WRITE (IO, 1) '   (none)'
              ENDIF
              DO 70 ROW = 1, N
                 P1 = INDEX (OFFPP-1 + ROW)
                 P2 = INDEX (OFFPP-1 + ROW+1) - 1
                 IF (P2 .GE. P1) THEN
                    WRITE (IO, 1) '   row: ', ROW
                    IF (PRL .EQ. 3) THEN
                       P2 = MIN (PRMAX, P2)
                    ENDIF
                    DO 60 P = P1, P2
                       WRITE (IO, 5)
     $                 INDEX (OFFIP-1 + P), VALUE (OFFXP-1 + P)
60                  CONTINUE
                 ENDIF
                 IF (PRL .EQ. 3 .AND. P2 .GE. PRMAX) THEN
                    WRITE (IO, 7)
                    GO TO 80
                 ENDIF
70            CONTINUE
80            CONTINUE
           ENDIF
           WRITE (IO, 8)
           IF (NBLKS .GT. 0) THEN
              IF (SYMBOL) THEN
                 WRITE (IO, 1) 'Nonzero pattern of prior LU factors:'
              ELSE
                 WRITE (IO, 1) 'LU factors:'
              ENDIF
           ENDIF
           PRN = 0
           DO 200 BLK = 1, NBLKS
              IF (NBLKS .GT. 1) THEN
                 K1 = INDEX (BLKPP-1 + BLK)
                 K2 = INDEX (BLKPP-1 + BLK+1) - 1
                 KN = K2-K1+1
                 SGLTON = KN .EQ. 1
                 IF (SGLTON) THEN
                    LUXXP = XP1-1 + INDEX (LUBLPP-1 + BLK)
                 ELSE
                    LUIIP = IP1-1 + INDEX (LUBLPP-1 + BLK)
                 ENDIF
                 IF (BLK .GT. 1) THEN
                    WRITE (IO, 9)
                 ENDIF
              ELSE
                 SGLTON = .FALSE.
                 K1 = 1
                 K2 = N
                 KN = N
                 LUIIP = IP1
              ENDIF
              IF (SGLTON) THEN
                 IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
                    WRITE (IO, 7)
                    GO TO 210
                 ENDIF
                 PRN = PRN + 1
                 IF (SYMBOL) THEN
                    WRITE (IO, 1) 'Block: ', BLK,
     $              ' (singleton) at index : ', K1
                 ELSE
                    WRITE (IO, 4) 'Block: ', BLK,
     $              ' (singleton) at index : ', K1,' value: ',
     $              VALUE (LUXXP)
                 ENDIF
                 IF (PRL .GE. 5) THEN
                    WRITE (IO, 1) 'located in VALUE ( ', LUXXP,' )'
                 ENDIF
              ELSE
                 LUXXP = XP1-1 + INDEX (LUIIP)
                 NLU = INDEX (LUIIP+1)
                 NPIV = INDEX (LUIIP+2)
                 MAXDC = INDEX (LUIIP+3)
                 MAXDR = INDEX (LUIIP+4)
                 LUPP = LUIIP+5
                 IF (NBLKS .GT. 1) THEN
                    WRITE (IO, 1) 'Block: ',BLK,' first index: ',K1,
     $              ' last index: ',K2
                 ENDIF
                 IF (PRL .GE. 5) THEN
                    WRITE (IO, 1) 'elements: ', NLU, ' pivots: ', NPIV
                    WRITE (IO, 1) 'largest contribution block: ',
     $                         MAXDC, ' by ', MAXDR
                    WRITE (IO, 1)'located in INDEX ( ',LUIIP,' ... )'
                    IF (.NOT. SYMBOL) THEN
                       WRITE (IO, 1) 'and in VALUE ( ',LUXXP,' ... )'
                    ENDIF
                 ENDIF
                 LUIIP = LUPP + NLU
                 K = 0
                 DO 190 E = 1, NLU
                    LUIP = LUIIP-1 + INDEX (LUPP-1 + E)
                    LUXP = LUXXP-1 + INDEX (LUIP)
                    LUK  = INDEX (LUIP+1)
                    LUDEGR = INDEX (LUIP+2)
                    LUDEGC = INDEX (LUIP+3)
                    LUNSON = INDEX (LUIP+4)
                    LUDIMR = INDEX (LUIP+5)
                    LUDIMC = INDEX (LUIP+6)
                    LUCP = LUIP + 7
                    LURP = LUCP + LUDEGC
                    LUSONP = LURP + LUDEGR
                    IF (PRL .GE. 5) THEN
                       WRITE (IO, 1) '   e: ', E, ' pivots: ', LUK
                       WRITE (IO, 1) '   children in dag: ', LUNSON,
     $                 ' frontal matrix: ', LUDIMR, ' by ', LUDIMC
                       ENDIF
                    P = LUXP
                    DO 140 J = 1, LUK
                       COL = K+J
                       NZCOL = LUK-J+1+LUDEGC
                       WRITE (IO, 1) '      L, col: ', COL
                       PRN = PRN + 1
                       ROW = COL
                       IF (SYMBOL) THEN
                          WRITE (IO, 5) ROW
                       ELSE
                          WRITE (IO, 5) ROW, ONE
                       ENDIF
                       P = P + 1
                       DO 120 I = J+1, LUK
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF
                          PRN = PRN + 1
                          ROW = K+I
                          IF (SYMBOL) THEN
                             WRITE (IO, 5) ROW
                          ELSE
                             WRITE (IO, 5) ROW, VALUE (P)
                          ENDIF
                          P = P + 1
120                    CONTINUE
                       DO 130 I = 1, LUDEGC
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF
                          PRN = PRN + 1
                          ROW = INDEX (LUCP-1+I)
                          IF (SYMBOL) THEN
                             WRITE (IO, 5) ROW
                          ELSE
                             WRITE (IO, 5) ROW, VALUE (P)
                          ENDIF
                          P = P + 1
130                    CONTINUE
                       P = P + J
140                 CONTINUE
                    UXP = LUXP + LUK*(LUDEGC+LUK)
                    DO 170 I = 1, LUK
                       ROW = K+I
                       NZROW = LUK-I+1+LUDEGR
                       WRITE (IO, 1) '      U, row: ', ROW
                       P = LUXP + (I-1) + (I-1) * (LUDEGC+LUK)
                       DO 150 J = I, LUK
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF
                          PRN = PRN + 1
                          COL = K+J
                          IF (SYMBOL) THEN
                             WRITE (IO, 5) COL
                          ELSE
                             WRITE (IO, 5) COL, VALUE (P)
                          ENDIF
                          P = P + (LUDEGC+LUK)
150                    CONTINUE
                       P = UXP
                       DO 160 J = 1, LUDEGR
                          IF (PRL .EQ. 3 .AND. PRN .GE. PRMAX) THEN
                             WRITE (IO, 7)
                             GO TO 210
                          ENDIF
                          PRN = PRN + 1
                          COL = INDEX (LURP-1+J)
                          IF (SYMBOL) THEN
                             WRITE (IO, 5) COL
                          ELSE
                             WRITE (IO, 5) COL, VALUE (P)
                          ENDIF
                          P = P + LUK
160                    CONTINUE
                       UXP = UXP + 1
170                 CONTINUE
                    IF (PRL .GE. 5) THEN
                       DO 180 I = 1, LUNSON
                          PRN = PRN + 1
                          SON = INDEX (LUSONP-1+I)
                          IF (SON .LE. KN) THEN
                             WRITE (IO, 1) '      LUson: ', SON
                          ELSE IF (SON .LE. 2*KN) THEN
                             WRITE (IO, 1) '      Uson:  ', SON-KN
                          ELSE
                             WRITE (IO, 1) '      Lson:  ', SON-2*KN
                          ENDIF
180                    CONTINUE
                    ENDIF
                    K = K + LUK
190              CONTINUE
              ENDIF
200        CONTINUE
210        CONTINUE
           IF (.NOT. BADLU) THEN
              PRN = MIN (PRMAX, N)
              IF (PRL .GE. 4) THEN
                 PRN = N
              ENDIF
              WRITE (IO, 8)
              WRITE (IO, 1) 'Column permutations'
              DO 220 I = 1, PRN
                 WRITE (IO, 5) INDEX (CPERMP+I-1)
220           CONTINUE
              IF (PRN .LT. N) THEN
                 WRITE (IO, 7)
              ENDIF
              WRITE (IO, 8)
              WRITE (IO, 1) 'Row permutations'
              DO 230 I = 1, PRN
                 WRITE (IO, 5) INDEX (RPERMP+I-1)
230           CONTINUE
              IF (PRN .LT. N) THEN
                 WRITE (IO, 7)
              ENDIF
           ENDIF
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (WHO .EQ. 3) THEN
           WRITE (IO, 8)
           PRN = MIN (PRMAX, N)
           IF (PRL .GE. 4) THEN
              PRN = N
           ENDIF
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 1) 'W (1 ... ',LW,
     $        ' ), work vector: not printed'
              WRITE (IO, 1) 'B (1 ... ',N,' ), right-hand side: '
              DO 240 I = 1, PRN
                 WRITE (IO, 5) I, B (I)
240           CONTINUE
              IF (PRN .LT. N) THEN
                 WRITE (IO, 7)
              ENDIF
           ELSE
              IF (INFO (1) .LT. 0) THEN
                 WRITE (IO, 1) 'W (1 ... ',LW,' ), work vector, and'
                 WRITE (IO, 1) 'X (1 ... ',N, ' ), solution,'
                 WRITE (IO, 1)
     $           '   not printed because of error flag, INFO (1) = ',
     $           INFO (1)
              ELSE
                 IF (LW .EQ. 4*N) THEN
                    WRITE (IO, 1) 'W (1 ... ',N,' ), residual: '
                    DO 250 I = 1, PRN
                       WRITE (IO, 5) I, W (I)
250                 CONTINUE
                    IF (PRN .LT. N) THEN
                       WRITE (IO, 7)
                    ENDIF
                    WRITE (IO, 1) 'W (',N+1,' ... ',LW,
     $              ' ), work vector: not printed'
                 ELSE
                    WRITE (IO, 1) 'W (1 ... ',LW,
     $              ' ), work vector: not printed'
                 ENDIF
                 WRITE (IO, 1) 'X (1 ... ',N,' ), solution: '
                 DO 260 I = 1, PRN
                    WRITE (IO, 5) I, X (I)
260              CONTINUE
                 IF (PRN .LT. N) THEN
                    WRITE (IO, 7)
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        IF (WHO .EQ. 1) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'end of MA38AD input '
           ELSE
              WRITE (IO, 6) 'end of MA38AD output'
           ENDIF
        ELSE IF (WHO .EQ. 2) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'end of MA38BD input '
           ELSE
              WRITE (IO, 6) 'end of MA38BD output'
           ENDIF
        ELSE IF (WHO .EQ. 3) THEN
           IF (WHERE .EQ. 1) THEN
              WRITE (IO, 6) 'end of MA38CD input '
           ELSE
              WRITE (IO, 6) 'end of MA38CD output'
           ENDIF
        ENDIF
        RETURN
C=======================================================================
C=======================================================================
1       FORMAT (' ', A, :, I12, :, A, :, I12, :,
     $               A, :, I12, :, A, :, I12, :, A, :, I12)
2       FORMAT (' ', I12, ': ', I12, ' ', I12, ' ', D11.4)
3       FORMAT (' ', A, D13.4, A)
4       FORMAT (' ', A, I12, A, I12, A, D11.4)
5       FORMAT (' ', I12, :, ': ', D11.4)
6       FORMAT (' ', 59('='), A)
7       FORMAT ('    ...')
8       FORMAT (' ', 79 ('-'))
9       FORMAT (' ', 79 ('.'))
        END
        SUBROUTINE MA38ZD (WHO, ERROR, I, J, X, IO)
        INTEGER WHO, ERROR, I, J, IO
        DOUBLE PRECISION X
C=== MA38ZD ============================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
        IF (IO .LT. 0) THEN
           RETURN
        ENDIF
        IF (WHO .EQ. 1) THEN
           IF (ERROR .EQ. -1) THEN
              WRITE (IO, 1) 'MA38AD: N less than one.'
           ELSE IF (ERROR .EQ. -2) THEN
              WRITE (IO, 1) 'MA38AD: NE less than one.'
           ELSE IF (ERROR .EQ. -3) THEN
              WRITE (IO, 1)
     $        'MA38AD: LINDEX too small.  Must be greater than ', I
           ELSE IF (ERROR .EQ. -4) THEN
              WRITE (IO, 1)
     $        'MA38AD: LVALUE too small.  Must be greater than ', I
           ELSE IF (ERROR .EQ. 1) THEN
              WRITE (IO, 1) 'MA38AD: ', I,
     $        ' invalid entries ignored (out of range 1..N).'
           ELSE IF (ERROR .EQ. 2) THEN
              WRITE (IO, 1) 'MA38AD: ', I,' duplicate entries summed.'
           ELSE IF (ERROR .EQ. 4) THEN
              WRITE (IO, 1)
     $        'MA38AD: matrix is singular.  Only ', I, ' pivots found.'
           ELSE IF (ERROR .EQ. 99) THEN
              WRITE (IO, 2)
     $        'MA38AD: invalid entry (out of range 1..N):', I, J, X
           ELSE IF (ERROR .EQ. 98) THEN
              WRITE (IO, 2)
     $        'MA38AD: duplicate entry summed:', I, J, X
           ENDIF
        ELSE IF (WHO .EQ. 2) THEN
           IF (ERROR .EQ. -1) THEN
              WRITE (IO, 1) 'MA38BD: N less than one.'
           ELSE IF (ERROR .EQ. -2) THEN
              IF (I .LT. 0) THEN
                 WRITE (IO, 1) 'MA38BD: NE less than one.'
              ELSE
                 WRITE (IO, 1)
     $           'MA38BD: NE too large.  Must be less than ', I
              ENDIF
           ELSE IF (ERROR .EQ. -3) THEN
              WRITE (IO, 1)
     $        'MA38BD: LINDEX too small.  Must be greater than ', I
           ELSE IF (ERROR .EQ. -4) THEN
              WRITE (IO, 1)
     $        'MA38BD: LVALUE too small.  Must be greater than ', I
           ELSE IF (ERROR .EQ. -6) THEN
              WRITE (IO, 1) 'MA38BD: pivot order from MA38AD failed.'
           ELSE IF (ERROR .EQ. -7) THEN
              WRITE (IO, 1)
     $        'MA38BD: LU factors uncomputed or corrupted.'
           ELSE IF (ERROR .EQ. 1) THEN
              IF (I .GT. 0) THEN
                 WRITE (IO, 1) 'MA38BD: ', I,
     $           ' invalid entries ignored (out of range 1..N).'
              ELSE
                 WRITE (IO, 1) 'MA38BD: ',-I,
     $           ' invalid entries ignored (not in prior pattern).'
              ENDIF
           ELSE IF (ERROR .EQ. 2) THEN
              WRITE (IO, 1) 'MA38BD: ', I,' duplicate entries summed.'
           ELSE IF (ERROR .EQ. 4) THEN
              WRITE (IO, 1) 'MA38BD: matrix is singular.  Only ', I,
     $        ' pivots found.'
           ELSE IF (ERROR .EQ. 99) THEN
              WRITE (IO, 2)
     $        'MA38BD: invalid entry (out of range 1..N):', I, J, X
           ELSE IF (ERROR .EQ. 98) THEN
              WRITE (IO, 2)
     $        'MA38BD: duplicate entry summed:', I, J, X
           ELSE IF (ERROR .EQ. 97) THEN
              WRITE (IO, 2)
     $        'MA38BD: invalid entry (not in pattern of prior factors)',
     $        I, J, X
           ELSE IF (ERROR .EQ. 96) THEN
              WRITE (IO, 2)
     $        'MA38BD: invalid entry (below diagonal blocks):', I, J, X
           ELSE IF (ERROR .EQ. 95) THEN
              WRITE (IO, 2)
     $        'MA38BD: invalid entry (prior matrix singular):', I, J, X
           ENDIF
        ELSE IF (WHO .EQ. 3) THEN
           IF (ERROR .EQ. -7) THEN
              WRITE (IO, 1)
     $        'MA38CD: LU factors uncomputed or corrupted.'
           ELSE IF (ERROR .EQ. 8) THEN
              IF (I .EQ. 0) THEN
                 WRITE (IO, 1)
     $  'MA38CD: no iterative refinement: original matrix not preserved'
              ELSE
                 WRITE (IO, 1)
     $  'MA38CD: no iterative refinement: only for Ax=b or A''x=b'
              ENDIF
           ENDIF
        ENDIF
        RETURN
C=======================================================================
C=======================================================================
1       FORMAT (' ', A, :, I12, :, A)
2       FORMAT (' ', A, /,
     $          '    row: ', I12, ' col: ', I12,' value: ', D11.4)
        END
