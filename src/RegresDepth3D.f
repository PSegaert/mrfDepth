      SUBROUTINE RDEPTH3(U,V,W,N1,X1,X2,Y,N2,RDEP,FLAG)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  This routine computes the regression depth of planes (U,V,W)
C  in a 3-dimensional data set (X1,X2,Y) of size N2.
C  Missing values are not allowed.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C N1     : INTEGER (1)
C        : The number of planes to compute depth off
C U      : DOUBLE PRECISION (N1)            
C        : X1 slope of the planes of which the depth is computed.
C V      : DOUBLE PRECISION (N1)            
C        : X2 slope of the planes of which the depth is computed.
C W      : DOUBLE PRECISION (N1)            
C        : Intercept of the plane of which the depth is computed.
C N2     : INTEGER (1)
C        : Size of the data set
C X1     : DOUBLE PRECISION (N2)            
C        : Measurements for the first variable. 
C X2     : DOUBLE PRECISION (N2)           
C        : Measurements for the second variable. 
C Y      : DOUBLE PRECISION (N2)           
C        : Measurements for the response variable. 
C EPS    : DOUBLE PRECISION (1)
C        : Numerical Tolerance
C     
C OUTPUT
C RDEP   : DOUBLE PRECISION (N1)
C        : The regression depth of the planes	
C FLAG   : INTEGER (N1)
C        : Indicates whether all points line lie on the samen line (FLAG = 1) or not (FLAG = 0)	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      INTEGER I,J,N1,N2
      DOUBLE PRECISION U(N1),V(N1),W(N1)
      DOUBLE PRECISION X1(N2),X2(N2),Y(N2)
      DOUBLE PRECISION XX1(N2),XX2(N2),YY(N2)     
      DOUBLE PRECISION RDEP(N1)
      INTEGER FLAG(N1)
      DOUBLE PRECISION EPS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      EPS=1.D-8
      
C Actual Routine
      DO 10 I=1,N1
       DO 11 J=1,N2
        XX1(J)=X1(J)
        XX2(J)=X2(J)
        YY(J)=Y(J)
11     CONTINUE       
       CALL RDEPTH31(XX1,XX2,YY,N2,U(I),V(I),W(I),RDEP(I),FLAG(I),EPS)
10    CONTINUE

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE RDEPTH31(X1,X2,Y,N,U,V,W,RDEP,FLAG,EPS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  This routine computes the regression depth of planes (U,V,W)
C  in a 3-dimensional data set (X1,X2,Y) of size N2.
C  Missing values are not allowed.
C  This subroutine was described in:
C       Rousseeuw, P.J. and Hubert, M. (1999). 
C       Regression Depth
C       Journal of the American Statistical Association, 94, 388-402.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C X1     : DOUBLE PRECISION (N)            
C        : Measurements for the first variable. 
C X2     : DOUBLE PRECISION (N)           
C        : Measurements for the second variable. 
C Y      : DOUBLE PRECISION (N)
C        : Measurements for the response variable. 
C N      : INTEGER (1)
C        : Size of the data set
C U      : DOUBLE PRECISION             
C        : X1 slope of the plane of which the depth is computed.
C V      : DOUBLE PRECISION             
C        : X2 slope of the plane of which the depth is computed.
C W      : DOUBLE PRECISION             
C        : Intercept of the plane of which the depth is computed.
C EPS    : DOUBLE PRECISION             
C        : Numerical Precision
C     
C OUTPUT
C RDEP   : DOUBLE PRECISION 
C        : The regression depth of the planes	
C FLAG   : INTEGER 
C        : Indicates whether all points line lie on the samen line (FLAG = 1) or not (FLAG = 0)	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,MAXN
      DOUBLE PRECISION X1(N),X2(N),Y(N)
      DOUBLE PRECISION RDEP,U,V,W,EPS,ALPHA(N),D
      INTEGER RESID(N),JRES(N)
      INTEGER NDEP,NNEGTOT,NPOSTOT,NDIM
      INTEGER FLAG
      
      INTEGER I
      
      MAXN=N
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
C  Compute all residuals, and find the total number of positive and
C  negative residuals.
C
      NNEGTOT=0
      NPOSTOT=0
      DO 10 I=1,N
         D=Y(I)-U*X1(I)-V*X2(I)-W
         IF (DABS(D).LE.EPS) THEN
            RESID(I)=0
         ELSEIF (D.GT.EPS) THEN
            RESID(I)=1
         ELSE 
            RESID(I)=-1
         ENDIF
         IF (RESID(I).LE.0) NNEGTOT=NNEGTOT+1
         IF (RESID(I).GE.0) NPOSTOT=NPOSTOT+1
 10   CONTINUE

      CALL STANDRDEPTH3(N,X1,X2,ALPHA,EPS)

      CALL RDEPTH31B(N,X1,X2,ALPHA,RESID,JRES,EPS,NDEP,
     +     NNEGTOT,NPOSTOT,NDIM)

      RDEP=(NDEP+0.D0)/(N+0.D0)
      FLAG=NDIM
      
      RETURN
 200  END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE STANDRDEPTH3(N,X1,X2,XN,EPS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine standardises and centers the data by use of the 
C coordinate wise median and mad.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N
      DOUBLE PRECISION X1(N),X2(N),XN(N),EPS

      CALL STANDRDEPTH31(N,X1,XN,EPS)
      CALL STANDRDEPTH31(N,X2,XN,EPS)

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE STANDRDEPTH31(N,X,XN,EPS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Subroutine for the standarisation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,I
      DOUBLE PRECISION X(N),XN(N),EPS
      DOUBLE PRECISION QLOC,QSCA,AVE,VAR,FINDQ

      DO 10 I=1,N
         XN(I)=X(I)
 10   CONTINUE         

      IF ((2*INT(N/2)).EQ.N) THEN
         QLOC=FINDQ(XN,N,N/2)
         QLOC=(FINDQ(XN,N,(N/2)+1)+QLOC)/2.D0
      ELSE
         QLOC=FINDQ(XN,N,INT(N/2)+1)
      ENDIF
      DO 30 I=1,N
         XN(I)=DABS(X(I)-QLOC)
 30   CONTINUE
      IF ((2*INT(N/2)).EQ.N) THEN
         QSCA=FINDQ(XN,N,N/2)
         QSCA=(FINDQ(XN,N,(N/2)+1)+QSCA)/2.D0
      ELSE
         QSCA=FINDQ(XN,N,INT(N/2)+1)
      ENDIF
      IF (DABS(QSCA).LT.EPS) THEN
         AVE=0.D0
         DO 40 I=1,N
            AVE=AVE+X(I)
 40      CONTINUE
         AVE=AVE/(N+0.D0)
         VAR=0.D0
         DO 50 I=1,N
            VAR=VAR+(X(I)-AVE)*(X(I)-AVE)
 50      CONTINUE  
         IF (N.NE.1) VAR=VAR/(N-1.D0)
         IF (DABS(VAR).LT.EPS) THEN
            QSCA=1.D0
         ELSE
            QSCA=DSQRT(VAR)
         ENDIF
      ENDIF
      DO 60 I=1,N
         X(I)=(X(I)-QLOC)/QSCA
 60   CONTINUE         

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE RDEPTH31B(N,X1,X2,ALPHA,RESID,JRES,EPS,NDEP,
     +     NNEGTOT,NPOSTOT,NDIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C To be called by RDEPTH31
C This subroutine performs the main part of the calculations.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N
      DOUBLE PRECISION X1(N),X2(N),ALPHA(N)
      DOUBLE PRECISION EPS,D,X1N,X2N,ANGLE
      DOUBLE PRECISION ALPHK,BETAK,P,P2
      INTEGER RESID(N),JRES(N)
      INTEGER I,NDIM,NDEP,NTPOS,NTNEG,NTNUL,NT
      INTEGER J,NN,NNEG,NNEGTOT,NNPOS,NPOSTOT
      INTEGER NU,ND,NNEGU,NPOSU,JA,JB,NN2,NPOS,NNNEG
      INTEGER LRPOS,LRNEG

      P=DACOS(-1.D0)
      P2=P*2.D0
      NDIM=2
C
C  Handle special cases where N is less or equal to 1.
C
      IF (N.LE.1) THEN
         NDEP=0
         IF ((N.EQ.1).AND.(RESID(1).EQ.0)) NDEP=1
         RETURN
      ENDIF
C
C  General case: initialize regression depth.
C
      NDEP=N
C
C  Main loop: use each object x(i) as rotation center.
C
      DO 100 I=1,N
C
C  Identify points which collapse with the rotation center:
C      NTNEG = ties with negative residual,
C      NTPOS = ties with positive residual,
C      NTNUL = ties with zero residual.
C  With the other points: construct the array ALPHA.
C      
         NTPOS=0
         NTNEG=0
         NTNUL=0
         NT=0
         DO 20 J=1,N
            D=DSQRT((X1(J)-X1(I))*(X1(J)-X1(I))+
     +           (X2(J)-X2(I))*(X2(J)-X2(I)))
            IF (D.LE.EPS) THEN
               IF (RESID(J).LT.0) NTNEG=NTNEG+1
               IF (RESID(J).GT.0) NTPOS=NTPOS+1
               IF (RESID(J).EQ.0) NTNUL=NTNUL+1
               NT=NT+1
            ELSE
               X1N=(X1(J)-X1(I))/D
               X2N=(X2(J)-X2(I))/D
               IF (ABS(X1N).GT.ABS(X2N)) THEN
                  IF (X1(J).GE.X1(I)) THEN
                     ALPHA(J-NT)=DASIN(X2N)
                     IF(ALPHA(J-NT).LT.0.D0) THEN
                        ALPHA(J-NT)=P2+ALPHA(J-NT)
                     ENDIF
                  ELSE
                     ALPHA(J-NT)=P-DASIN(X2N)
                  ENDIF
               ELSE
                  IF (X2(J).GE.X2(I)) THEN
                     ALPHA(J-NT)=DACOS(X1N)
                  ELSE
                     ALPHA(J-NT)=P2-DACOS(X1N)
                  ENDIF
               ENDIF
               IF (ALPHA(J-NT).GE.(P2-EPS)) ALPHA(J-NT)=0.D0
               JRES(J-NT)=RESID(J)
            ENDIF
 20      CONTINUE
C
C  Compute 
C      NN = the total number of points, different from the rotation center.
C      NNNEG = the total number of points with negative residual,
C              different from the rotation center,
C      NNPOS = the total number of points with positive residual,
C              different from the rotation center.
C
         NN=N-NT
         NNNEG=NNEGTOT-NTNEG-NTNUL
         NNPOS=NPOSTOT-NTPOS-NTNUL
C
C  If all ties have identical signs of residuals, they must not be 
C  taken into account (in that case the line can be translated a little
C  such that the center lies on the appropriate side), otherwise
C  the number of zero residuals is added to both the number of 
C  positive and negative residuals (a zero residual has always to be 
C  changed to make theta a nonfit).
C
         IF ((NTPOS.EQ.NT).OR.(NTNEG.EQ.NT)) THEN
            NTPOS=0
            NTNEG=0
         ELSE
            NTPOS=NTPOS+NTNUL
            NTNEG=NTNEG+NTNUL
         ENDIF
C
C  Special case: all points are identical.
C
         IF (NN.LT.1) THEN
            NDEP=MIN0(NTPOS,NTNEG)
            RETURN
         ENDIF
C
C  Sort the array ALPHA together with the residuals.
C
         CALL SORTRDEPTH3(ALPHA,JRES,NN)
C
C  Make smallest alpha equal to zero,
C  and compute NU = number of alpha < pi;
C              NNEGU = number of negative residuals up to alpha(nu);
C              NPOSU = number of positive residuals up to alpha(nu).
C
         ANGLE=ALPHA(1)
         NU=0
         ND=0
         NNEGU=0
         NPOSU=0
         DO 30 J=1,NN
            ALPHA(J)=ALPHA(J)-ANGLE
            IF (ALPHA(J).LT.(P-EPS)) THEN
               NU=NU+1
               IF (JRES(J).LE.0) NNEGU=NNEGU+1
               IF (JRES(J).GE.0) NPOSU=NPOSU+1
            ENDIF
            IF ((DABS(ALPHA(J)).LE.EPS).OR.
     +           (DABS(ALPHA(J)-P).LE.EPS)) ND=ND+1
 30      CONTINUE
         IF (ND.EQ.NN) NDIM=1

C
C  Mergesort the alpha with their antipodal angles beta,
C  and at the same time compute the regression depth.
C
         JA=1
         JB=1
         ALPHK=ALPHA(1)
         IF ((NU+1).LE.NN) THEN
            BETAK=ALPHA(NU+1)-P
         ELSE
            BETAK=ALPHA(NU+1-NN)+P
         ENDIF
         NN2=NN*2
         NPOS=NNPOS
         NNEG=NNNEG
         DO 40 J=1,NN2
            IF (ALPHK.LT.(BETAK+EPS)) THEN
               IF(JRES(JA).LE.0)NNEG=NNEG+1
               IF(JRES(JA).GE.0)NPOS=NPOS+1
               IF (JA.LT.NN) THEN
                  JA=JA+1
                  ALPHK=ALPHA(JA)
               ELSE
                  ALPHK=P2+1.D0
               ENDIF
            ELSE
               IF ((JB+NU).LE.NN) THEN
                  IF(JRES(JB+NU).GE.0)NPOSU=NPOSU+1
                  IF(JRES(JB+NU).LE.0)NNEGU=NNEGU+1
               ELSE
                  IF(JRES(JB+NU-NN).GE.0)NPOSU=NPOSU+1
                  IF(JRES(JB+NU-NN).LE.0)NNEGU=NNEGU+1
               ENDIF
               IF (NPOSU.EQ.(NNPOS+1)) THEN
                  NPOSU=1
                  NPOS=NPOS-NNPOS
               ENDIF
               IF (NNEGU.EQ.(NNNEG+1)) THEN
                  NNEGU=1
                  NNEG=NNEG-NNNEG
               ENDIF
               ANGLE=BETAK
               IF (JB.LT.NN) THEN
                  JB=JB+1
                  IF ((JB+NU).LE.NN) THEN
                     BETAK=ALPHA(JB+NU)-P
                  ELSE
                     BETAK=ALPHA(JB+NU-NN)+P
                  ENDIF
               ELSE
                  BETAK=P2+1.D0
               ENDIF
               IF (DABS(ANGLE-BETAK).GT.EPS) THEN
                  LRPOS=NPOS-NPOSU
                  LRNEG=NNNEG-(NNEG-NNEGU)
                  IF (NTPOS.LT.NTNEG) THEN
                     NDEP=MIN0(NDEP,LRPOS+LRNEG+NTPOS)
                  ELSE
                     NDEP=MIN0(NDEP,LRPOS+LRNEG+NTNEG)
                  ENDIF
                  LRPOS=NNPOS-LRPOS
                  LRNEG=NNNEG-LRNEG
                  IF (NTPOS.LT.NTNEG) THEN
                     NDEP=MIN0(NDEP,LRPOS+LRNEG+NTPOS)
                  ELSE
                     NDEP=MIN0(NDEP,LRPOS+LRNEG+NTNEG)
                  ENDIF
               ENDIF
            ENDIF
 40      CONTINUE         
 100  CONTINUE
      
      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE SORTRDEPTH3(B,RESID,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  Sorts an array B (of length N) in O(NlogN) time.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N
      DOUBLE PRECISION B(N),AMM,XX
      INTEGER RESID(N)
      INTEGER JLV(N),JRV(N)
      INTEGER JSS,JNDL,JR,JNC,J,JTWE,NRES
      
      JSS=1
      JLV(1)=1
      JRV(1)=N
  10  JNDL=JLV(JSS)
      JR=JRV(JSS)
      JSS=JSS-1
  20  JNC=JNDL
      J=JR
      JTWE=(JNDL+JR)/2
      XX=B(JTWE)
  30  IF (B(JNC).GE.XX) GOTO 40
      JNC=JNC+1
      GOTO 30
  40  IF (XX.GE.B(J)) GOTO 50
      J=J-1
      GOTO 40
  50  IF (JNC.GT.J) GOTO 60
      AMM=B(JNC)
      B(JNC)=B(J)
      B(J)=AMM
      NRES=RESID(JNC)
      RESID(JNC)=RESID(J)
      RESID(J)=NRES
      JNC=JNC+1
      J=J-1
  60  IF (JNC.LE.J) GOTO 30
      IF ((J-JNDL).LT.(JR-JNC)) GOTO 80
      IF (JNDL.GE.J) GOTO 70
      JSS=JSS+1
      JLV(JSS)=JNDL
      JRV(JSS)=J
  70  JNDL=JNC
      GOTO 100
  80  IF (JNC.GE.JR) GOTO 90
      JSS=JSS+1
      JLV(JSS)=JNC
      JRV(JSS)=JR
  90  JR=J
 100  IF (JNDL.LT.JR) GOTO 20
      IF (JSS.NE.0) GOTO 10
      RETURN
      END
