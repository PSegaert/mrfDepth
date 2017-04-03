      SUBROUTINE RDEPTH4(T,N1,X1,X2,X3,Y,N2,RDEP,FLAG)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  This routine computes the regression depth of
C  hyperplanes T in a 4-dimensional data set (X1,X2,X3,Y) of size N2.
C  Missing values are not allowed.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C N1     : INTEGER (1)
C        : The number of hyperplanes to compute the depth off.
C T      : DOUBLE PRECISION (N1x4)            
C        : Hyperplanes of which the depth is computed.
C        : The first 3 elements are the slopes in the directions of
C        : the 3 variables describing the X-space.
C N2     : INTEGER (1)
C        : Size of the data set
C X1     : DOUBLE PRECISION (N2)            
C        : Measurements for the first variable. 
C X2     : DOUBLE PRECISION (N2)           
C        : Measurements for the second variable. 
C X3     : DOUBLE PRECISION (N2)           
C        : Measurements for the third variable. 
C Y      : DOUBLE PRECISION (N2)           
C        : Measurements for the response variable.
C EPS    : DOUBLE PRECISION             
C        : Numerical Precision 
C     
C OUTPUT
C RDEP   : DOUBLE PRECISION (N1)
C        : The regression depth of the planes   
C FLAG   : INTEGER (N1)
C        : Indicates whether all points line lie on the same
C        : plane (FLAG = 2)
C        : line (FLAG = 1)
C        : or not (FLAG = 0)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      INTEGER N1,N2
      DOUBLE PRECISION T(N1,4)
      DOUBLE PRECISION X1(N2),X2(N2),X3(N2),Y(N2)
      DOUBLE PRECISION XX1(N2),XX2(N2),XX3(N2),YY(N2)
      DOUBLE PRECISION RDEP(N1)
      INTEGER FLAG(N1)
      DOUBLE PRECISION EPS

      INTEGER I,J
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      
      EPS=1.D-8

      DO 10 I=1,N1
        DO 11 J=1,N2
            XX1(J)=X1(J)
            XX2(J)=X2(J)
            XX3(J)=X3(J)
            YY(J)=Y(J)
11      CONTINUE 
        CALL RDEPTH41(T(I,:),XX1,XX2,XX3,YY,
     +                N2,RDEP(I),FLAG(I),EPS)
10    CONTINUE

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE RDEPTH41(T,X1,X2,X3,Y,N,RDEP,FLAG,EPS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  This routine computes the regression depth of
C  a hyperplane T in a 4-dimensional data set (X1,X2,X3,Y) of size N.
C  Missing values are not allowed.
C  This subroutine was described in:
C       Rousseeuw, P.J. and Struyf, A. (1998). 
C       Computing location depth and regression depth 
C       in higher dimensions
C       Statistics and Computing, 8, 193-203.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C T      : DOUBLE PRECISION (4)            
C        : Hyperplane of which the depth is computed.
C        : The first 3 elements are the slopes in the directions of
C        : the 3 variables describing the X-space.
C X1     : DOUBLE PRECISION (N)            
C        : Measurements for the first variable. 
C X2     : DOUBLE PRECISION (N)           
C        : Measurements for the second variable. 
C X3     : DOUBLE PRECISION (N)           
C        : Measurements for the third variable. 
C Y      : DOUBLE PRECISION (N)           
C        : Measurements for the response variable. 
C N      : INTEGER (1)
C        : Size of the data set
C EPS    : DOUBLE PRECISION             
C        : Numerical Precision
C     
C OUTPUT
C RDEP   : DOUBLE PRECISION (1)
C        : The regression depth of the planes   
C FLAG   : INTEGER (1)
C        : Indicates whether all points line lie on the same
C        : plane (FLAG = 2)
C        : line (FLAG = 1)
C        : or not (FLAG = 0)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X1(N),X2(N),X3(N),Y(N)
      DOUBLE PRECISION XN1(N),XN2(N),XN3(N)
      DOUBLE PRECISION RDEP,T(4),EPS,ALPHA(N),D
      INTEGER RESID(N),JRES(N)
      INTEGER FLAG

      MAXN=N
      
C
C  Compute all residuals, and find the total number of positive and
C  negative residuals.
C
      NNEGTOT=0
      NPOSTOT=0
      DO 10 I=1,N
         D=Y(I)-T(1)*X1(I)-T(2)*X2(I)-T(3)*X3(I)-T(4)
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

      CALL STANDRDEPTH4(N,X1,X2,X3,XN1,EPS)

      CALL RDEPTH41B(T,N,X1,X2,X3,XN1,XN2,XN3,Y,ALPHA,RESID,JRES,
     +     EPS,NDEP,NNEGTOT,NPOSTOT,NDIM)

      RDEP=(NDEP+0.D0)/(N+0.D0)      
      FLAG=NDIM
      
      RETURN
 200  END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE STANDRDEPTH4(N,X1,X2,X3,XN,EPS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine standardises and centers the data by use of the 
C coordinate wise median and mad.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X1(N),X2(N),X3(N),XN(N),EPS

      CALL STANDRDEPTH41(N,X1,XN,EPS,1)
      CALL STANDRDEPTH41(N,X2,XN,EPS,2)
      CALL STANDRDEPTH41(N,X3,XN,EPS,3)

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE STANDRDEPTH41(N,X,XN,EPS,J)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Subroutine for the standarisation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N),XN(N),EPS
      DOUBLE PRECISION QLOC,QSCA,AVE,VAR

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






      SUBROUTINE RDEPTH41B(T,N,X1,X2,X3,XN1,XN2,XN3,Y,ALPHA,RESID,JRES,
     +     EPS,NDEP,NNEGTOT,NPOSTOT,NDIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C To be called by RDEPTH41
C This subroutine performs the main part of the calculations.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X1(N),X2(N),X3(N),Y(N),ALPHA(N)
      DOUBLE PRECISION T(4),EPS,D,ANGLE
      DOUBLE PRECISION XN1(N),XN2(N),XN3(N),P,P2
      DOUBLE PRECISION ALPHK,BETAK,A(2,3),B(3,2)
      INTEGER RESID(N),JRES(N),NTNEG(3),NTPOS(3),NTNUL(3)

      P=DACOS(-1.D0)
      P2=P*2.D0
      NDIM=3
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
C  Main loops: use each line (x(ic),x(i)) as rotation center.
C
      DO 200 IC=1,N
         DO 15 I=1,N
            XN1(I)=X1(I)-X1(IC)
            XN2(I)=X2(I)-X2(IC)
            XN3(I)=X3(I)-X3(IC)
 15      CONTINUE
         DO 20 I=1,N
            IF ((DABS(XN1(I)).LE.EPS).AND.(DABS(XN2(I)).LE.EPS).AND.
     +           (DABS(XN3(I)).LE.EPS)) THEN
               GOTO 20
            ENDIF
C
C  Calculate the matrix of the orthogonal projection on the plane through 
C  x(ic), orthogonal to the line through x(ic) and x(i).
C  Let the third coordinate coincide with the line (x(ic),x(i)).
C
            IF (DABS(XN1(I)).GT.EPS) THEN
               B(2,1)=1.D0
               B(3,1)=1.D0
               B(1,1)=-(XN2(I)+XN3(I))/XN1(I)
            ELSEIF (DABS(XN2(I)).GT.EPS) THEN
               B(1,1)=1.D0
               B(3,1)=1.D0
               B(2,1)=-(XN1(I)+XN3(I))/XN2(I)
            ELSE
               B(1,1)=1.D0
               B(2,1)=1.D0
               B(3,1)=-(XN1(I)+XN2(I))/XN3(I)
            ENDIF
            B(1,2)=B(2,1)*XN3(I)-B(3,1)*XN2(I)
            B(2,2)=B(3,1)*XN1(I)-B(1,1)*XN3(I)
            B(3,2)=B(1,1)*XN2(I)-XN1(I)*B(2,1)
            
            A(1,1)=(B(2,2)*XN3(I)-XN2(I)*B(3,2))
            A(1,2)=-(B(1,2)*XN3(I)-XN1(I)*B(3,2))
            A(1,3)=(B(1,2)*XN2(I)-B(2,2)*XN1(I))
            A(2,1)=-(B(2,1)*XN3(I)-XN2(I)*B(3,1))
            A(2,2)=(B(1,1)*XN3(I)-XN1(I)*B(3,1))
            A(2,3)=-(B(1,1)*XN2(I)-XN1(I)*B(2,1))
C
C Compute the new planar coordinates for all points.
C Identify points which collapse with the rotation center:
C     NTNUL = real ties, 
C     NTPOS = the original point lies on the positive side of 
C             the projection plane,
C     NTNEG = the original point lies on the negative side of 
C             the projection plane.
C For each of these arrays:
C     first component = ties with negative residual,
C     second component = ties with zero residual,
C     third component = ties with positive residual.
C With the other points: construct the array ALPHA.
C      
            DO 25 J=1,3
               NTNUL(J)=0
               NTPOS(J)=0
               NTNEG(J)=0
 25         CONTINUE
            NT=0
            DO 30 J=1,N
               T(1)=XN1(J)*A(1,1)+XN2(J)*A(1,2)+XN3(J)*A(1,3)
               T(2)=XN1(J)*A(2,1)+XN2(J)*A(2,2)+XN3(J)*A(2,3)
               IF ((DABS(T(1)).LE.EPS).AND.
     +              (DABS(T(2)).LE.EPS)) THEN
                  NT=NT+1
                  D=XN1(J)*XN1(I)+XN2(J)*XN2(I)+XN3(J)*XN3(I)
                  IF (DABS(D).LE.EPS) THEN
                     IF (RESID(J).LT.0) NTNUL(1)=NTNUL(1)+1
                     IF (RESID(J).GT.0) NTNUL(3)=NTNUL(3)+1
                     IF (RESID(J).EQ.0) NTNUL(2)=NTNUL(2)+1
                  ELSEIF (D.GT.EPS) THEN 
                     IF (RESID(J).LT.0) NTPOS(1)=NTPOS(1)+1
                     IF (RESID(J).GT.0) NTPOS(3)=NTPOS(3)+1
                     IF (RESID(J).EQ.0) NTPOS(2)=NTPOS(2)+1
                  ELSE 
                     IF (RESID(J).LT.0) NTNEG(1)=NTNEG(1)+1
                     IF (RESID(J).GT.0) NTNEG(3)=NTNEG(3)+1
                     IF (RESID(J).EQ.0) NTNEG(2)=NTNEG(2)+1
                  ENDIF
               ELSE
                  D=DSQRT(T(1)*T(1)+T(2)*T(2))
                  ALPHK=T(1)/D
                  BETAK=T(2)/D
                  IF (ABS(ALPHK).GT.ABS(BETAK)) THEN
                     IF (T(1).GE.(0.D0-EPS)) THEN
                        ALPHA(J-NT)=DASIN(BETAK)
                        IF(ALPHA(J-NT).LT.0.D0) THEN
                           ALPHA(J-NT)=P2+ALPHA(J-NT)
                        ENDIF
                     ELSE
                        ALPHA(J-NT)=P-DASIN(BETAK)
                     ENDIF
                  ELSE
                     IF (T(2).GE.(0.D0-EPS)) THEN
                        ALPHA(J-NT)=DACOS(ALPHK)
                     ELSE
                        ALPHA(J-NT)=P2-DACOS(ALPHK)
                     ENDIF
                  ENDIF
                  IF (ALPHA(J-NT).GE.(P2-EPS)) ALPHA(J-NT)=0.D0
                  JRES(J-NT)=RESID(J)
               ENDIF
 30         CONTINUE
C
C  Compute 
C      NN = the total number of points, different from the rotation center.
C      NNNEG = the total number of points with negative residual,
C              different from the rotation center,
C      NNPOS = the total number of points with positive residual,
C              different from the rotation center.
C
            NN=N-NT
            NU=NTNUL(2)+NTNEG(2)+NTPOS(2)
            NNEGU=NTPOS(1)+NTNUL(1)+NTNEG(1)
            NPOSU=NTPOS(3)+NTNUL(3)+NTNEG(3)
            NNNEG=NNEGTOT-NNEGU-NU
            NNPOS=NPOSTOT-NPOSU-NU
C
C  If all ties have identical signs of residuals, they must not be 
C  taken into account (in that case the plane can be translated a little
C  such that the center lies on the appropriate side), otherwise
C  the number of zero residuals is added to both the number of 
C  positive and negative residuals (a zero residual has always to be 
C  changed to make theta a nonfit).
C
            IF ((NPOSU.EQ.NT).OR.(NNEGU.EQ.NT)) THEN
               NTPOS(1)=0
               NTNEG(1)=0
               NTNUL(1)=0
               NTPOS(3)=0
               NTNEG(3)=0
               NTNUL(3)=0
            ELSE
               NTNEG(1)=NTNEG(1)+NTNEG(2)
               NTPOS(1)=NTPOS(1)+NTPOS(2)
               NTNUL(1)=NTNUL(1)+NTNUL(2)
               NTNEG(3)=NTNEG(3)+NTNEG(2)
               NTPOS(3)=NTPOS(3)+NTPOS(2)
               NTNUL(3)=NTNUL(3)+NTNUL(2)
            ENDIF
C
C  Compute
C      NT = the minimal number of ties which residuals always have to 
C           change sign, even if small translations and rotations of
C           the plane (only changing the position of the ties) 
C           are performed.
C
            NT=MIN0(MIN0(NTPOS(3)+NTNUL(1)+NTNEG(1),
     +                   NTPOS(1)+NTNUL(3)+NTNEG(3)),
     +              MIN0(NTPOS(3)+NTNUL(3)+NTNEG(1),
     +                   NTPOS(1)+NTNUL(1)+NTNEG(3)))
            NT=MIN0(NT,MIN0(NPOSU+NU,NNEGU+NU))
C
C  Special case: all projections are identical.
C
            IF (NN.LT.1) THEN
               NDEP=MIN0(NDEP,NT)
               NDIM=1
               GOTO 20
            ENDIF
C
C  Sort the array ALPHA, together with the residuals.
C
            CALL SORTRDEPTH4(ALPHA,JRES,NN)
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
            DO 55 J=1,NN
               ALPHA(J)=ALPHA(J)-ANGLE
               IF (ALPHA(J).LT.(P-EPS)) THEN
                  NU=NU+1
                  IF (JRES(J).LE.0) NNEGU=NNEGU+1
                  IF (JRES(J).GE.0) NPOSU=NPOSU+1
               ENDIF
               IF ((DABS(ALPHA(J)).LE.EPS).OR.
     +              (DABS(ALPHA(J)-P).LE.EPS)) ND=ND+1
 55         CONTINUE
            IF (ND.EQ.NN) NDIM=2

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
            DO 60 J=1,NN2
               IF (ALPHK.LE.(BETAK+EPS)) THEN
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
                     NDEP=MIN0(NDEP,LRPOS+LRNEG+NT)
                     LRPOS=NNPOS-LRPOS
                     LRNEG=NNNEG-LRNEG
                     NDEP=MIN0(NDEP,LRPOS+LRNEG+NT)
                  ENDIF
               ENDIF
 60         CONTINUE
 20      CONTINUE
 200  CONTINUE

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE SORTRDEPTH4(B,RESID,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  Sorts an array B (of length N) in O(NlogN) time.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION B(N),AMM,XX
      INTEGER RESID(N)
      DIMENSION JLV(N),JRV(N)
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
