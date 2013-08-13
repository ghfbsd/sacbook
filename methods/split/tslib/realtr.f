c
c
c TITLE - REALTR = REAL TRANSFORM
c     FOURIER TRANSFORM OF REAL SERIES FROM OUTPUT OF FFT
c
c              IF ISN=1, THIS SUBROUTINE COMPLETES THE FOURIER TRANS-
c              FORM OF 2*N REAL DATA VALUES, WHERE THE ORIGINAL DATA
c              VALUES ARE STORED ALTERNATELY IN ARRAYS A AND B, AND ARE
c              FIRST TRANSFORMED BY A COMPLEX FOURIER TRANSFORM OF
c              DIMENSION N.
c              THE COSINE COEFFICIENTS ARE IN A(1),A(2),...A(N+1) AND
c              THE SINE COEFFICIENTS ARE IN B(1),B(2),...B(N+1).
c              A TYPICAL CALLING SEQUENCE IS
c              CALL FFT (A,B,N,-1)
c              CALL REALTR (A,B,N,1)
c              THE RESULTS SHOULD BE MULTIPLIED BY 0.5/N TO GIVE THE
c              USUAL SCALING OF COEFFICIENTS.
c              IF ISN=1, THE INVERSE TRANSFORMATION IS DONE, THE
c              FIRST STEP IN EVALUATING A REAL FOURIER SERIES.
c              A TYPICAL CALLING SEQUENCE IS
c              CALL REALTR (A,B,N,-1)
c              CALL FFT (A,B,N,-1)
c              THE RESULTS SHOULD BE MULTIPLIED BY 0.5 TO GIVE THE
c              USUAL SCALING, AND THE TIME DOMAIN RESULTS ALTERNATE
c              IN ARRAYS A AND B, I.E. A(1),B(1), A(2),B(2),
c              ...A(N),B(N).
c              THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE
c              COMPLEX ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO
c              TWO TO GIVE THE CORRECT INDEXING INCREMENT AND A(2)
c              USED TO PASS THE INITIAL ADDRESS FOR THE SEQUENCE OF
c              IMAGINARY VALUES,E.G.
c              CALL FFT(A,A(2),N,2)
c              CALL REALTR(A,A(2),N,2)
c              IN THIS CASE, THE COSINE AND SINE COEFFICIENTS
c              ALTERNATE IN A.
c
c INPUTS
c
c     A(I)            I=1, IABS(ISN)*N, IABS(ISN) CONTAINS THE REAL
c              PART OF THE FOURIER TRANSFORM OF A REAL DATA SERIES
c              TREATED AS COMPLEX.
c     B(I)            CONTAINS THE IMAGINARY PART CORRESPONDING TO A
c     N               THE NUMBER OF COMPLEX NUMBERS REPRESENTED IN
c              A AND B
c     ISN             IS NEGATIVE FOR THE INVERSE TRANSFORM AND
c              POSITIVE FOR THE TRANSFORM IABS(ISN) IS THE ADVANCE
c              FOR INDEXING THROUGH A OR B.
c
c
c OUTPUTS
c
c     A(1)     THE REAL PART OF THE FOURIER TRANSFORM OF A REAL DATA
c              SERIES
c     B(1)     THE IMAGINARY PART OF THE FOURIER TRANSFORM OF A
c              REAL DATA SERIES
c     F        = 0 IF NO ERRORS
c              = 1 IF N .LE. 0
c
c
c EXAMPLES
c
c  SUPPOSE  X  IS A REAL ARRAY WITH  N  ELEMENTS (N=2**L), AND
c  WE WISH THE DIGITAL FOURIER COSINE AND SINE TRANSFORMS OF IT;
c  THEN THESE CALLS ARE REQUIRED (USING PDP-11 ROUTINES)
c      CALL FFTTWO(X,N/2)
c      CALL REALTR(X,X(2),N/2,2)
c  NOW  X  CONTAINS THE COS AND SIN TRANSFORMS (UNNORMALIZED) ,
c  STARTING AT ZERO FREQUENCY AND ALTERNATING COS,SIN,COS,SIN,....
c  NOTE THAT  N  HAS BEEN HALVED IN THE CALLS TO FFTTWO  AND  REALTR
c
c ***** WARNING - THE TRANSFORMS OF THE NYQUIST FREQUENCY 
c      APPEAR IN THE POSITIONS  X(N+1),X(N+2)  SO THAT
c      THE ARRAY  X  MUST BE DIMENSIONED AT LEAST  X(N+2) IN
c      THE CALLING  PROGRAM   ***************
c
c
      subroutine realtr(a, b, n, isn)
      dimension a(1), b(1)
      real im
      if (n .le. 1) return 
      inc = iabs(isn)
      nk = (n * inc) + 2
      nh = nk / 2
      sd = (2. * atan(1.)) / float(n)
      cd = 2. * (sin(sd) ** 2)
      sd = sin(sd + sd)
      sn = 0.
      if (isn .gt. 0) goto 10
      cn = -1.
      sd = - sd
      goto 20
   10 cn = 1.
      a(nk - 1) = a(1)
      b(nk - 1) = b(1)
   20 do 30 j = 1, nh, inc
      k = nk - j
      aa = a(j) + a(k)
      ab = a(j) - a(k)
      ba = b(j) + b(k)
      bb = b(j) - b(k)
      re = (cn * ba) + (sn * ab)
      im = (sn * ba) - (cn * ab)
      b(k) = im - bb
      b(j) = im + bb
      a(k) = aa - re
      a(j) = aa + re
      aa = cn - ((cd * cn) + (sd * sn))
      sn = ((sd * cn) - (cd * sn)) + sn
c
   30 cn = aa
      return 
c
      end
