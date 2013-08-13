      subroutine mtaper(pimult,n,nwin,npta,el,ta)
C     MTAPER -- Generate Slepian tapers for leakage-resistant spectral
C               estimation.
C
C     Assumes:
C        pimult - taper order (the multiple of pi of the desired taper, e.g.
C                2.5 pi tapers would have pimult value of 2.5)
C        n - number of points in time series
C        nwin - number of tapers wanted (nwin << n, note!)
C        npta - physical dimension of ta array (ta(npta,nwin) is dimension)
C
C     Returns:
C        el - eigenvalues/bandwidth retention factors of the nwin tapers
C             produced (real*8)
C        ta - taper values for nwin tapers ((ta(1:n,i),i=1,nwin) are values)
C             (real*8)
C
C     Code originally by Jeff Park (Yale Univ.), late '80s I think.
C     Hacked by G. Helffrich (Tokyo Inst. of Technology) to modularize by
C        removing parameter passing through common areas.
C     April 7, 2001.

      logical debug
      parameter (ntmax=20,nscr=2**14,debug=.false.)
      real*8 el(nwin),ta(npta,nwin)
      real*8 pi,ww,cs,eps,rlu,rlb
      real*8 dfac,drat,gamma,bh,ell,aa,tapsq
      real*8 z(nscr,8)
      integer ip(ntmax)

      data pi/3.14159265358979d0/

C     Sanity checking
      if (nwin.ge.n) pause '**MTAPER: nwin not << n! - invalid use.'
      if (ntmax.lt.nwin) then
         write(0,*) '**MTAPER: nwin too large, ',ntmax,
     x       ' max - recompile.'
         pause
      endif
      if (n .gt. nscr) then
         write(0,*) '**MTAPER: n too large, ',nscr,' max. - recompile.'
	 pause
      endif

      if (debug) then
         do i=1,nwin
	    do j=1,n
	       ta(j,i)=1.0
	    enddo
	    el(i)=1.0
	 enddo
	 return
      endif

C     Error tolerance for computed eigenvalues.
      eps=1.d-13

      ww=dble(pimult)/dble(n)
      cs=dcos(2.d0*pi*ww)

C     Set up matrix for eispack subroutine solution
      do i=0,n-1
        z(i+1,1)=-cs*(dble(n-1)/2.d0-dble(i))**2
        z(i+1,2)=-dble(i)*dble(n-i)/2.d0
        z(i+1,3)=z(i+1,2)**2
      end do
C     Eigenvalues returned to caller in el(1:nwin)
      call tridib(n,eps,z(1,1),z(1,2),z(1,3),rlb,rlu,1,nwin,el,ip,
     x       ierr,z(1,4),z(1,5))
      if (ierr .ne. 0) pause '**MTAPER: multiple eigenvalues...'
C     Eigenvectors returned to caller in ta(1:n,1:nwin)
      call tinvit(npta,n,z(1,1),z(1,2),z(1,3),nwin,el,ip,
     x       ta,ierr,z(1,4),z(1,5),z(1,6),z(1,7),z(1,8))
      if (ierr .ne. 0) pause '**MTAPER: unconverged eigenvector...'

C     We calculate the eigenvalues of the dirichlet-kernel problem
C     (i.e. the bandwidth retention factors) from Slepian 1978 asymptotic
C     formula, gotten from Thomson 1982 eq 2.5, and supplemented by the
C     asymptotic formula for k near 2n from Slepian 1978 eq 61
      dfac=n*pi*ww
      drat=8.d0*dfac
      dfac=4.d0*dsqrt(pi*dfac)*dexp(-2.d0*dfac)
      do k=1,nwin
        el(k)=1.d0-dfac
        dfac=dfac*drat/k  ! is this correct formula? yes,but fails as k -> 2n
      end do
      gamma=dlog(8.d0*n*dsin(2.d0*pi*ww))+0.5772156649d0
      ii = 1
      do k=1,nwin
        bh=-2.d0*pi*(n*ww-(k-1)/2.d0-.25d0)/gamma
        ell=1.d0/(1.d0+dexp(pi*bh))
        el(k)=max(ell,el(k))

C       Normalize the eigentapers to preserve power for a white process
C       i.e. they have rms value unity.  Also ensure that alternate tapers
C       change sign in their first derivatives, and first is positive (ii).
        tapsq=0.d0
        do i=1,n
	  tapsq=tapsq+ta(i,k)**2
        end do
c       aa=dsqrt(tapsq/n)
        aa=ii*sign(dsqrt(tapsq/n),ta(2,k)-ta(1,k))
        do i=1,n
	  ta(i,k)=ta(i,k)/aa
        end do
c       Comment out following line if you want all tapers with positive first
c       derivative at t=0
        ii = -ii
      end do
      end

      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
C
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      REAL*8 D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      REAL*8 U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,MACHEP
      REAL*8 DSQRT,DABS,DFLOAT
      INTEGER IND(M)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0D0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0D0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE;
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES;
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER;
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT:
C
C        ALL INPUT ARRAYS ARE UNALTERED;
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS;
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
C  FOR F_FLOATING DEC FORTRAN
C      DATA MACHEP/1.1D-16/
C  FOR G_FLOATING DEC FORTRAN
       DATA MACHEP/1.25D-15/
C  FOR IEEE-754 FLOATING REPRESENTATION
C      DATA MACHEP/2.22D-16/
C
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = 1.0D0 - E2(1)
      Q = 0
C     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX ::::::::::
  100 P = Q + 1
C
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140
  120 CONTINUE
C     :::::::::: FIND VECTORS BY INVERSE ITERATION ::::::::::
  140 TAG = TAG + 1
      S = 0
C
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
C     :::::::::: CHECK FOR ISOLATED ROOT ::::::::::
         XU = 1.0D0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0D0
         GO TO 870
  490    NORM = DABS(D(P))
         IP = P + 1
C
         DO 500 I = IP, Q
  500    NORM = NORM + DABS(D(I)) + DABS(E(I))
C     :::::::::: EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ::::::::::
         EPS2 = 1.0D-3 * NORM
         EPS3 = MACHEP * NORM
         UK = DBLE(Q-P+1)
         EPS4 = UK * EPS3
         UK = EPS4 / DSQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
C     :::::::::: LOOK FOR CLOSE OR COINCIDENT ROOTS ::::::::::
  510    IF (DABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3
C     :::::::::: ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ::::::::::
  520    V = 0.0D0
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (DABS(E(I)) .LT. DABS(U)) GO TO 540
C     :::::::::: WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ::::::::::
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0D0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0D0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
C
         IF (U .EQ. 0.0D0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0D0
         RV3(Q) = 0.0D0
C     :::::::::: BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ::::::::::
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
C     :::::::::: ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ::::::::::
         IF (GROUP .EQ. 0) GO TO 700
         J = R
C
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = 0.0D0
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = 0.0D0
C
         DO 720 I = P, Q
  720    NORM = NORM + DABS(RV6(I))
C
         IF (NORM .GE. 1.0D0) GO TO 840
C     :::::::::: FORWARD SUBSTITUTION ::::::::::
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0D0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     :::::::::: ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ::::::::::
  780    DO 820 I = IP, Q
            U = RV6(I)
C     :::::::::: IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ::::::::::
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
C
         ITS = ITS + 1
         GO TO 600
C     :::::::::: SET ERROR -- NON-CONVERGED EIGENVECTOR ::::::::::
  830    IERR = -R
         XU = 0.0D0
         GO TO 870
C     :::::::::: NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ::::::::::
  840    U = 0.0D0
C
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
C
         XU = 1.0D0 / DSQRT(U)
C
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0D0
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
      IF (Q .LT. N) GO TO 100
 1001 RETURN
C     :::::::::: LAST CARD OF TINVIT ::::::::::
      END

      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      REAL*8 D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      REAL*8 U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP
      REAL*8 DABS,DMAX1,DMIN1,DFLOAT
      INTEGER IND(M)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,
C     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,
C     USING BISECTION.
C
C     ON INPUT:
C
C        N IS THE ORDER OF THE MATRIX;
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY;
C
C        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED
C          EIGENVALUES;
C
C        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER
C          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.
C
C     ON OUTPUT:
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE;
C
C        D AND E ARE UNALTERED;
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO;
C
C        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED
C          EIGENVALUES;
C
C        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES
C          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER;
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE;
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER
C     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
C  FOR F_FLOATING DEC FORTRAN
C      DATA MACHEP/1.1D-16/
C  FOR G_FLOATING DEC FORTRAN
       DATA MACHEP/1.25D-15/
C  FOR IEEE-754 FLOATING REPRESENTATION
C      DATA MACHEP/2.22D-16/
C
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0D0
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES ::::::::::
      DO 40 I = 1, N
         X1 = U
         U = 0.0D0
         IF (I .NE. N) U = DABS(E(I+1))
         XU = DMIN1(D(I)-(X1+U),XU)
         X0 = DMAX1(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         IF (DABS(E(I)) .GT. MACHEP * (DABS(D(I)) + DABS(D(I-1))))
     X      GO TO 40
   20    E2(I) = 0.0D0
   40 CONTINUE
C
      X1 = DMAX1(DABS(XU),DABS(X0)) * MACHEP * DBLE(N)
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
C     :::::::::: DETERMINE AN INTERVAL CONTAINING EXACTLY
C                THE DESIRED EIGENVALUES ::::::::::
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5D0
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
C     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ::::::::::
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0D0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0D0
         V = 0.0D0
         IF (Q .EQ. N) GO TO 110
         U = DABS(E(Q+1))
         V = E2(Q+1)
  110    XU = DMIN1(D(Q)-(X1+U),XU)
         X0 = DMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0D0) GO TO 140
  120 CONTINUE
C
  140 X1 = DMAX1(DABS(XU),DABS(X0)) * MACHEP
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     :::::::::: CHECK FOR ISOLATED ROOT WITHIN INTERVAL ::::::::::
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * DBLE(Q-P+1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     :::::::::: FIND ROOTS BY BISECTION ::::::::::
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     :::::::::: LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ::::::::::
      K = M2
  250    XU = LB
C     :::::::::: FOR I=K STEP -1 UNTIL M1 DO -- ::::::::::
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     :::::::::: NEXT BISECTION STEP ::::::::::
  300    X1 = (XU + X0) * 0.5D0
         IF ((X0 - XU) .LE. (2.0D0 * MACHEP *
     X      (DABS(XU) + DABS(X0)) + DABS(EPS1))) GO TO 420
C     :::::::::: IN-LINE PROCEDURE FOR STURM SEQUENCE ::::::::::
  320    S = P - 1
         U = 1.0D0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0D0) GO TO 325
            V = DABS(E(I)) / MACHEP
            IF (E2(I) .EQ. 0.0D0) V = 0.0D0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0D0) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     :::::::::: REFINE INTERVALS ::::::::::
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     :::::::::: K-TH EIGENVALUE FOUND ::::::::::
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     :::::::::: ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ::::::::::
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     :::::::::: SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
C                EXACTLY THE DESIRED EIGENVALUES ::::::::::
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
C     :::::::::: LAST CARD OF TRIDIB ::::::::::
      END
