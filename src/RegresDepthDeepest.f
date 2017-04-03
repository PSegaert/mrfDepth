      SUBROUTINE SWEEPMEDRES(X,N,NP,A,MAXIT,ITER,MDEPAPPR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Algorithm to approximate the deepest regression in higher dimensions.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C X        : DOUBLE PRECISION (N,P) #1
C n        : integer(1)             #2
C np       : integer(1)             #3
C A        : DOUBLE PRECISION (p)   #4
C maxit    : integer(1)             #5
C iter     : integer(1)             #6
C MDEPAPPR : integer(1)             #7

C P.J. Rousseeuw, S. Van Aelst, An algorithm for deepest multiple regression, 
C in: Proceedings in Computational Statistics, 2000, Physica-Verlag, Heidelberg, 2000, p. 139.

  
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NP
      DOUBLE PRECISION X(N,NP),A(NP),Y(N)
      DOUBLE PRECISION XVAR(N),XMED(N)
      DOUBLE PRECISION B(NP,NP),DPMEDX(NP),XTWEE(N,NP)
      DOUBLE PRECISION RES(N),EPS,DpMEDIAN_REGDepth,dpmed,C
      DOUBLE PRECISION YMED
      INTEGER MAXIT,ITER,MDEPAPPR,index,nconv,mdep
      
      INTEGER I,J,K
      
      EPS=1.D-9

        DO 10 I=1,N
           DO 20 J=1,NP 
           XTWEE(I,J)=X(I,J)
 20     CONTINUE
 10        CONTINUE 


      DO 103 I=1,NP
         DO 104 J=1,NP
             B(I,J)=0 D0
  104    CONTINUE
  103 CONTINUE    

CC   Preprocessing: First the median of the x-variables is determined 

      DO 1000 J=1,(NP-1)
        DO 1001 I=1,N
           XVAR(I)=X(I,J)
1001    CONTINUE
        DpMEDX(J)=DpMEDIAN_REGDepth(N,XVAR)
1000  CONTINUE 
   
CC   Construction of the sweeping variables: Out of each variable X_j (j>1) 
CC   we sweep the variables  X_1,...,X_{j-1} by using regression through the
CC   origin.
 
      DO 108 J=2,(NP-1)
         DO 109 I=1,N
            XVAR(I)=X(I,J)
  109    CONTINUE
         DO 110 K=1,(J-1)
               DpMED=DpMEDIAN_REGDepth(N,XVAR)
               INDEX=0               
            DO 111 I=1,N
               DENOM=X(I,K)-DpMEDX(K)               
               IF ((ABS(DENOM)).GT.EPS) THEN
                  INDEX=INDEX+1
                  XMED(INDEX)=(XVAR(I)-DpMED)/DENOM
               ENDIF               
  111       CONTINUE
            B(K,J)=DpMEDIAN_REGDepth(INDEX,XMED)
            DO 112 I=1,N
               XVAR(I)=XVAR(I)-(B(K,J)*X(I,K))
  112       CONTINUE
  110    CONTINUE
         DO 113 I=1,N
            X(I,J)=XVAR(I)
  113    CONTINUE   
  108 CONTINUE
 
CC  We substract the median of the sweeping variables.

      DO 1016 J=1,(NP-1)
         DO 1017 I=1,N
            XVAR(I)=X(I,J)
 1017    CONTINUE
         DPMEDX(J)=DpMEDIAN_REGDepth(N,XVAR)
         DO 1018 I=1,N  
            X(I,J)=XVAR(I)-DPMEDX(J)
 1018    CONTINUE
 1016 CONTINUE  

CC  The actual sweeping. We repeat the procedure until convergence.
CC  The maximal number of iterations equals 100. 

      DO 114 I=1,N
         Y(I)=X(I,NP)
  114 CONTINUE      
      DO 115 J=1,NP
         A(J)=0 D0
  115 CONTINUE
      NCONV=0
      DO 116 K=1,MAXIT
         DO 117 J=1,(NP-1)
            C=0 D0
            YMED=DpMEDIAN_REGDepth(N,Y)
            INDEX=0 
            DO 118 I=1,N
               IF (ABS(X(I,J)).GT.EPS) THEN
                  INDEX=INDEX+1
                  XMED(INDEX)=(Y(I)-YMED)/X(I,J)
               ENDIF         
  118       CONTINUE
            C=DpMEDIAN_REGDepth(INDEX,XMED)
            A(J)=A(J)+C     
            DO 119 I=1,N
               Y(I)=Y(I)-C*(X(I,J)+DPMEDX(J))
  119       CONTINUE
            IF ((NCONV .EQ. 0) .AND. (ABS(C) .GT. EPS)) THEN
               NCONV=1
            ENDIF  
  117    CONTINUE                    
         IF (NCONV .EQ. 0) THEN
            ITER=K 
            GOTO 500
         ELSE 
            NCONV=0
         ENDIF
  116 CONTINUE            
      ITER=MAXIT
  500 A(NP)=DpMEDIAN_REGDepth(N,Y)    

CC   We determine the values of the parameters.

      DO 120 J=(NP-2),1,-1
         DO 121 K=(J+1),(NP-1)
            A(J)=A(J)-(A(K)*B(J,K))
  121    CONTINUE
  120 CONTINUE

      DO 300 I=1,N
         RES(I)=XTWEE(I,NP)-A(NP)
  300 CONTINUE
         DO 301 J=1,(NP-1)
         DO 302 I=1,N
            RES(I)=RES(I)-(XTWEE(I,J)*A(J))
  302    CONTINUE
  301 CONTINUE       
CC      write(*,*)(A(j),j=1,5)

CC We make np residuals zero, this means that we let the fit pass through 
CC np observations.

      CALL SECTIONPOINT(XTWEE,A,N,NP)

CC We determine the residuals of the final fit.

      DO 200 I=1,N
         Y(I)=XTWEE(I,NP)-A(NP)
  200 CONTINUE
      DO 201 J=1,(NP-1)
         DO 202 I=1,N
            Y(I)=Y(I)-(XTWEE(I,J)*A(J))
  202    CONTINUE
  201 CONTINUE            

CC  We determine Min_{i=1,..(p-1)}(rdepth(fit,(x_i,y_i)). 

      MDEPAPPR=N  
      DO 123 J=1,(NP-1)
         DO 124 I=1,N
            XVAR(I)=X(I,J)
            XMED(I)=Y(I)
 124     CONTINUE           
         CALL DSORT(XVAR,XMED,N,2)
         CALL RDEPTH(XVAR,XMED,N,MDEP) 
         IF (MDEP .LT. MDEPAPPR) THEN
            MDEPAPPR=MDEP
         ENDIF 
 123  CONTINUE         
       DO 125 J=2,(NP-1)
         DO 126 I=1,N
            XVAR(I)=XTWEE(I,J)
            XMED(I)=Y(I)
 126     CONTINUE 
         CALL DSORT(XVAR,XMED,N,2)
         CALL RDEPTH(XVAR,XMED,N,MDEP) 
         IF (MDEP .LT. MDEPAPPR) THEN
            MDEPAPPR=MDEP
         ENDIF 
 125  CONTINUE   
      END               
    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION DpMEDIAN_REGDepth(N,X)
CC Computes the median of a data set. DATA SET MUST NOT BE ORDERED.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,I
      DOUBLE PRECISION X(N),XN(N),DpMED,DpMEDIAN_REGDepth
      DOUBLE PRECISION FINDQ
      

      DO 10 I=1,N
         XN(I)=X(I)
 10   CONTINUE         
      IF ((2*INT(N/2)).EQ.N) THEN
        DpMED= FINDQ(XN,N,N/2)
        DpMED=(FINDQ(XN,N,(N/2)+1)+DpMED)/2.D0
      ELSE
        DpMED=FINDQ(XN,N,INT(N/2)+1)
      ENDIF
      DpMEDIAN_REGDepth=DpMED 
      RETURN 
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RDEPTH(X,RES,LENGTH,DEPTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DEPTH,RPLUS,RMIN
      DOUBLE PRECISION EPS
      DIMENSION X(LENGTH),RES(LENGTH),LNEG(LENGTH),LPOS(LENGTH)
      INTEGER NPOS,NNEG,I,K
      EPS=1.D-9
      
 
      DEPTH=LENGTH
      NPOS=0
      NNEG=0
      DO 300 I=1,LENGTH
         R=RES(I)
      IF (R .LT. -EPS) THEN
      LNEG(I)=1
      LPOS(I)=0
      NNEG=NNEG+1
      ELSE IF (R .GT. EPS) THEN
      LPOS(I)=1
      LNEG(I)=0
      NPOS=NPOS+1
      ELSE
      LNEG(I)=1
      LPOS(I)=1
      NNEG=NNEG+1
      NPOS=NPOS+1
      END IF
  300 CONTINUE

      LPLUS=0
      LMIN=0
      RPLUS=NPOS
      RMIN=NNEG
      DO 400 K=1,LENGTH
      LPLUS=LPLUS+LPOS(K)
      LMIN=LMIN+LNEG(K)
      RPLUS=RPLUS-LPOS(K)
      RMIN=RMIN-LNEG(K)
      IF(K .NE. LENGTH) THEN 
        IF(X(K).NE.X(K+1)) THEN 
          MINXI=MIN(LPLUS+RMIN,LMIN+RPLUS)
          IF (MINXI .LT. DEPTH) DEPTH=MINXI
        ENDIF
      ELSE
        MINXI=MIN(LPLUS+RMIN,LMIN+RPLUS)
        IF (MINXI .LT. DEPTH) DEPTH=MINXI      
      ENDIF
  400 CONTINUE
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DSORT (DX, DY, N, KFLAG)
C***BEGIN PROLOGUE  DSORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   DSORT sorts array DX and optionally makes the same interchanges in
C   array DY.  The array DX may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      DX - array of values to be sorted   (usually abscissas)
C      DY - array to be (optionally) carried along
C      N  - number of values in array DX to be sorted
C      KFLAG - control parameter
C            =  2  means sort DX in increasing order and carry DY along.
C            =  1  means sort DX in increasing order (ignoring DY)
C            = -1  means sort DX in decreasing order (ignoring DY)
C            = -2  means sort DX in decreasing order and carry DY along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891024  Changed category.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to DX,DY; changed
C           code to parallel SSORT. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  DSORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION DX(*), DY(*)
C     .. Local Scalars ..
      DOUBLE PRECISION R, T, TT, TTY, TY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
C      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  DSORT
      NN = N
C      IF (NN .LT. 1) THEN
C         CALL XERMSG ('SLATEC', 'DSORT',
C     +      'The number of values to be sorted is not positive.', 1, 1)
C         RETURN
C      ENDIF
C
      KK = ABS(KFLAG)
C      IF (KK.NE.1 .AND. KK.NE.2) THEN
C         CALL XERMSG ('SLATEC', 'DSORT',
C     +      'The sort control parameter, K, is not 2, 1, -1, or -2.', 2,
C     +      1)
C         RETURN
C      ENDIF
C
C     Alter array DX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            DX(I) = -DX(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort DX only
C
      M = 1
      I = 1
      J = NN
      R = 0.375D0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (DX(L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (DX(K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         GO TO 40
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = DX(I+1)
      IF (DX(I) .LE. T) GO TO 80
      K = I
C
   90 DX(K+1) = DX(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 90
      DX(K+1) = T
      GO TO 80
C
C     Sort DX and carry DY along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375D0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
      TY = DY(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (DX(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (DX(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = DX(I+1)
      TY = DY(I+1)
      IF (DX(I) .LE. T) GO TO 170
      K = I
C
  180 DX(K+1) = DX(K)
      DY(K+1) = DY(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 180
      DX(K+1) = T
      DY(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            DX(I) = -DX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SECTIONPOINT(XTWEE,STARTP,N,NP) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N,NP,I
      INTEGER MAXPP1
      DIMENSION XTWEE(N,NP),YS(N,NP)
      DIMENSION STARTP(NP),UI(NP),INDEX(NP),HVEC2(NP*NP)
      DOUBLE PRECISION NOMIN,EPS,R
      
      EPS=1.D-9

      MAXPP1=NP*NP
      DO 81 J=1,NP
         INDEX(J)=0
   81 CONTINUE     
      R=1 D12
      DO 90 I = 1,N
         RI=0
         DO 91 J=1,(NP-1)
            RI=RI-(STARTP(J)*XTWEE(I,J))
   91    CONTINUE
         RI=RI+XTWEE(I,NP)-STARTP(NP)
         IF (ABS(RI) .LT. ABS(R)) THEN
         R=RI
         INDEX(1)=I
         ENDIF
   90 CONTINUE                       
      STARTP(NP)=STARTP(NP)+R
      DO 92 K=1,NP-1
         ALPHA=1 D9
         CALL NEWDIR(XTWEE,N,NP,K,K+1,MAXPP1,YS,INDEX,UI,HVEC2)
         DO 93 L=1,N   
            DENOM=UI(NP)
            NOMIN=XTWEE(L,NP)-STARTP(NP)
            DO 94 LL=1,NP-1
               NOMIN=NOMIN-(XTWEE(L,LL)*STARTP(LL))
               DENOM=DENOM+(XTWEE(L,LL)*UI(LL))
   94       CONTINUE    
            IF ((ABS(NOMIN) .LT. EPS) .AND. (ABS(DENOM) .GT. EPS)) THEN
               ALPHA=0 D0
               INDEX(K+1)=L
            ELSE IF((ABS(NOMIN).GT.EPS).AND.(ABS(DENOM).GT.EPS)) THEN
               ALPHAI=NOMIN/DENOM
               IF (ABS(ALPHAI) .LT. ABS(ALPHA)) THEN
                  ALPHA=ALPHAI
                  INDEX(K+1)=L
               ENDIF 
            ENDIF
   93    CONTINUE
         DO 95 KK=1,NP
            STARTP(KK)=STARTP(KK)+(ALPHA*UI(KK))
   95    CONTINUE
   92 CONTINUE  
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        
      SUBROUTINE NEWDIR(XTWEE,N,NP,K,KK,MAXPP1,YS,INDEX,UI,HVEC2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N,NP,K,KK,LL,I,J,MAXPP1
      DIMENSION XTWEE(N,NP),YS(K,KK)
      DIMENSION UI(NP),INDEX(NP),HVEC2(MAXPP1)
      DO 1000, I=1,K
         YS(I,1)=1
         DO 1001 J=1,K-1
            YS(I,J+1)=XTWEE(INDEX(I),NP+J-K)    
 1001    CONTINUE
         YS(I,K+1)=-XTWEE(INDEX(I),NP-K)      
 1000 CONTINUE
      CALL EQUAT(YS,K,K+1,HVEC2,MAXPP1,K,1,NERR)
      DO 1002 L=1,NP-K-1
         UI(L)=0 D0
 1002 CONTINUE
      UI(NP-K)=1 D0
      DO 1003 LL=1,K-1
         UI(NP+LL-K)=YS(1+LL,1)           
 1003 CONTINUE  
      UI(NP)=YS(1,1)   
      END  


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




