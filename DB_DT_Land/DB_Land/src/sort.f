      SUBROUTINE SORT(X,N,Y)
C
C     PURPOSE--THIS SUBROUTINE SORTS (IN ASCENDING ORDER)
C              THE N ELEMENTS OF THE SINGLE PRECISION VECTOR X
C              AND PUTS THE RESULTING N SORTED VALUES INTO THE
C              SINGLE PRECISION VECTOR Y.
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                OBSERVATIONS TO BE SORTED. 
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C     OUTPUT ARGUMENTS--Y      = THE SINGLE PRECISION VECTOR
C                                INTO WHICH THE SORTED DATA VALUES
C                                FROM X WILL BE PLACED.
C     OUTPUT--THE SINGLE PRECISION VECTOR Y
C             CONTAINING THE SORTED
C             (IN ASCENDING ORDER) VALUES
C             OF THE SINGLE PRECISION VECTOR X.
C     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS. 
C     RESTRICTIONS--THE DIMENSIONS OF THE VECTORS IL AND IU 
C                   (DEFINED AND USED INTERNALLY WITHIN
C                   THIS SUBROUTINE) DICTATE THE MAXIMUM
C                   ALLOWABLE VALUE OF N FOR THIS SUBROUTINE.
C                   IF IL AND IU EACH HAVE DIMENSION K,
C                   THEN N MAY NOT EXCEED 2**(K+1) - 1.
C                   FOR THIS SUBROUTINE AS WRITTEN, THE DIMENSIONS
C                   OF IL AND IU HAVE BEEN SET TO 36,
C                   THUS THE MAXIMUM ALLOWABLE VALUE OF N IS
C                   APPROXIMATELY 137 BILLION.
C                   SINCE THIS EXCEEDS THE MAXIMUM ALLOWABLE
C                   VALUE FOR AN INTEGER VARIABLE IN MANY COMPUTERS,
C                   AND SINCE A SORT OF 137 BILLION ELEMENTS
C                   IS PRESENTLY IMPRACTICAL AND UNLIKELY,
C                   THEN THERE IS NO PRACTICAL RESTRICTION
C                   ON THE MAXIMUM VALUE OF N FOR THIS SUBROUTINE.
C                   (IN LIGHT OF THE ABOVE, NO CHECK OF THE 
C                   UPPER LIMIT OF N HAS BEEN INCORPORATED
C                   INTO THIS SUBROUTINE.)
C     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN. 
C     COMMENT--THE SMALLEST ELEMENT OF THE VECTOR X
C              WILL BE PLACED IN THE FIRST POSITION
C              OF THE VECTOR Y,
C              THE SECOND SMALLEST ELEMENT IN THE VECTOR X
C              WILL BE PLACED IN THE SECOND POSITION
C              OF THE VECTOR Y, ETC.
C     COMMENT--THE INPUT VECTOR X REMAINS UNALTERED.
C     COMMENT--IF THE ANALYST DESIRES A SORT 'IN PLACE',
C              THIS IS DONE BY HAVING THE SAME
C              OUTPUT VECTOR AS INPUT VECTOR IN THE CALLING SEQUENCE. 
C              THUS, FOR EXAMPLE, THE CALLING SEQUENCE
C              CALL SORT(X,N,X)
C              IS ALLOWABLE AND WILL RESULT IN
C              THE DESIRED 'IN-PLACE' SORT.
C     COMMENT--THE SORTING ALGORTHM USED HEREIN
C              IS THE BINARY SORT.
C              THIS ALGORTHIM IS EXTREMELY FAST AS THE
C              FOLLOWING TIME TRIALS INDICATE.
C              THESE TIME TRIALS WERE CARRIED OUT ON THE
C              UNIVAC 1108 EXEC 8 SYSTEM AT NBS
C              IN AUGUST OF 1974.
C              BY WAY OF COMPARISON, THE TIME TRIAL VALUES
C              FOR THE EASY-TO-PROGRAM BUT EXTREMELY
C              INEFFICIENT BUBBLE SORT ALGORITHM HAVE
C              ALSO BEEN INCLUDED--
C              NUMBER OF RANDOM        BINARY SORT       BUBBLE SORT
C               NUMBERS SORTED
C                N = 10                 .002 SEC          .002 SEC
C                N = 100                .011 SEC          .045 SEC
C                N = 1000               .141 SEC         4.332 SEC
C                N = 3000               .476 SEC        37.683 SEC
C                N = 10000             1.887 SEC      NOT COMPUTED
C     REFERENCES--CACM MARCH 1969, PAGE 186 (BINARY SORT ALGORITHM
C                 BY RICHARD C. SINGLETON).
C               --CACM JANUARY 1970, PAGE 54.
C               --CACM OCTOBER 1970, PAGE 624.
C               --JACM JANUARY 1961, PAGE 41.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE--301-921-2315
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --NOVEMBER  1975. 
C
C---------------------------------------------------------------------
C
      DIMENSION X(1),Y(1)
      DIMENSION IU(36),IL(36) 
C
      IPR=6
C
C     CHECK THE INPUT ARGUMENTS FOR ERRORS
C
      IF(N.LT.1)GOTO50
      IF(N.EQ.1)GOTO55
      HOLD=X(1)
      DO60I=2,N
      IF(X(I).NE.HOLD)GOTO90
   60 CONTINUE
      WRITE(IPR, 9)HOLD
      DO61I=1,N
      Y(I)=X(I)
   61 CONTINUE
      RETURN
   50 WRITE(IPR,15) 
      WRITE(IPR,47)N
      RETURN
   55 WRITE(IPR,18) 
      Y(1)=X(1)
      RETURN
   90 CONTINUE
    9 FORMAT(1H ,108H***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUME
     1NT (A VECTOR) TO THE SORT   SUBROUTINE HAS ALL ELEMENTS = ,E15.8,6
     1H *****)
   15 FORMAT(1H , 91H***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE
     1 SORT   SUBROUTINE IS NON-POSITIVE *****)
   18 FORMAT(1H ,100H***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUME
     1NT TO THE SORT   SUBROUTINE HAS THE VALUE 1 *****)
   47 FORMAT(1H , 35H***** THE VALUE OF THE ARGUMENT IS ,I8   ,6H *****)
C
C-----START POINT-----------------------------------------------------
C
C     COPY THE VECTOR X INTO THE VECTOR Y
      DO100I=1,N
      Y(I)=X(I)
  100 CONTINUE
C
C     CHECK TO SEE IF THE INPUT VECTOR IS ALREADY SORTED
C
      NM1=N-1
      DO200I=1,NM1
      IP1=I+1
      IF(Y(I).LE.Y(IP1))GOTO200
      GOTO250
  200 CONTINUE
      RETURN
  250 M=1 
      I=1 
      J=N 
  305 IF(I.GE.J)GOTO370
  310 K=I 
      MID=(I+J)/2
      AMED=Y(MID)
      IF(Y(I).LE.AMED)GOTO320 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
  320 L=J 
      IF(Y(J).GE.AMED)GOTO340 
      Y(MID)=Y(J)
      Y(J)=AMED
      AMED=Y(MID)
      IF(Y(I).LE.AMED)GOTO340 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
      GOTO340
  330 Y(L)=Y(K)
      Y(K)=TT
  340 L=L-1
      IF(Y(L).GT.AMED)GOTO340 
      TT=Y(L)
  350 K=K+1
      IF(Y(K).LT.AMED)GOTO350 
      IF(K.LE.L)GOTO330
      LMI=L-I
      JMK=J-K
      IF(LMI.LE.JMK)GOTO360
      IL(M)=I
      IU(M)=L
      I=K 
      M=M+1
      GOTO380
  360 IL(M)=K
      IU(M)=J
      J=L 
      M=M+1
      GOTO380
  370 M=M-1
      IF(M.EQ.0)RETURN
      I=IL(M)
      J=IU(M)
  380 JMI=J-I
      IF(JMI.GE.11)GOTO310
      IF(I.EQ.1)GOTO305
      I=I-1
  390 I=I+1
      IF(I.EQ.J)GOTO370
      AMED=Y(I+1)
      IF(Y(I).LE.AMED)GOTO390 
      K=I 
  395 Y(K+1)=Y(K)
      K=K-1
      IF(AMED.LT.Y(K))GOTO395 
      Y(K+1)=AMED
      GOTO390
      END 