      SUBROUTINE HSDEP2(UU,VV,N1,XX,YY,N2,HDEP,SDEP)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is a wrapper to calculate the halfspace depths and
C simplicial depths of a set of bivariate points (U,V) given a 
C certain dataset (X,Y)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C U      : DOUBLE PRECISION (N1)            
C        : First coordinates of the points to calculate depth 
C V      : DOUBLE PRECISION (N1)           
C        : Second coordinates of the points to calculate depth 
C N1     : INTEGER (1)
C        : Number of points to calculate depth of
C X      : DOUBLE PRECISION (N2)            
C        : First coordinates of the points in dataset 
C Y      : DOUBLE PRECISION (N2)           
C        : Second coordinates of the points in dataset 
C N2     : INTEGER (1)
C        : Number of points in the dataset
C     
C OUTPUT
C HDEP   : DOUBLE PRECISION (N1)
C        : Halfspace depth of the points (U,V)  
C SDEP   : DOUBLE PRECISION (N1)
C        : Simplicial depth of the points (U,V) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
C INPUT
      INTEGER N1,N2
      DOUBLE PRECISION UU(N1),VV(N1),XX(N2),YY(N2)
C OUTPUT
      DOUBLE PRECISION HDEP(N1), SDEP(N1)
C WORKSPACE
      DOUBLE PRECISION U(N1),V(N1),X(N2),Y(N2)
      DOUBLE PRECISION X2(N2),Y2(N2)
      DOUBLE PRECISION xmean,ymean,xdev,ydev
      DOUBLE PRECISION BETA(N2)
      DOUBLE PRECISION EPS
      INTEGER I,J, F(N2)
      DOUBLE PRECISION DPF(N2)
      INTEGER JLV(N2),JRV(N2)
      INTEGER HDEPT
      DOUBLE PRECISION SDEPT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine

C Initialise constants
      EPS=1.D-8

C
C     Standardization of the data set.
C
      xmean=0.0
      ymean=0.0
      xdev=0.0
      ydev=0.0
      DO 2 I=1,N2
         xmean=xmean+xx(i)
         ymean=ymean+yy(i)
2     CONTINUE
      xmean=xmean/n2
      ymean=ymean/n2
      DO 3 I=1,N2
         xdev=xdev+((xx(i)-xmean)*(xx(i)-xmean))
         ydev=ydev+((yy(i)-ymean)*(yy(i)-ymean))
3     continue
      xdev=xdev/(N2-1)
      ydev=ydev/(N2-1)
      xdev=dsqrt(xdev)
      ydev=dsqrt(ydev)
      DO 5 I=1,N2
         if (xdev.gt.eps) THEN
            x(i)=(xx(i)-xmean)/xdev
         else
            x(i)=xx(i)
         endif 
         if (ydev.gt.eps) THEN
            y(i)=(yy(i)-ymean)/ydev
         else
            y(i)=yy(i)
         endif         
5     CONTINUE

      DO 6 I=1,N1
         if (xdev.gt.eps) THEN
            u(i)=(uu(i)-xmean)/xdev
         else
            u(i)=uu(i)
         endif 
         if (ydev.gt.eps) THEN
            v(i)=(vv(i)-ymean)/ydev
         else
            v(i)=vv(i)
         endif         
6     CONTINUE

      DO 123 i=1,N1
         DO 124 j=1,N2
            X2(J)=X(J)
            Y2(J)=Y(J)
  124    CONTINUE
        CALL HSDEP21(U(i),V(i),N2,X2,Y2,BETA,F,DPF,
     +              JLV,JRV,HDEPT,SDEPT)

        HDEP(i) = (HDEPT+0.d0)/(N2+0.d0)
        SDEP(i) = SDEPT
  123 CONTINUE

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE HSDEP21(U,V,N,X,Y,
     +                   BETA,F,DPF,JLV,JRV,HDEP,SDEP)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  This routine calculates the halfspace and the simplicial depth 
C  of a point theta = (U,V), given a certain dataset (X,Y).
C  This subroutine was described in:
C       Rousseeuw, P.J. and Ruts, I. (1996). AS 307:  Bivariate location 
C       depth, Applied Statistics (JRSS-C) 45, 516-526.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C U      : DOUBLE PRECISION (1) 
C        : First coordinate of the point to calculate depth 
C V      : DOUBLE PRECISION (1)           
C        : Second coordinate of the point to calculate depth 
C N      : INTEGER (1)
C        : Number of points in the dataset
C X      : DOUBLE PRECISION (N)            
C        : First coordinates of the points in dataset 
C Y      : DOUBLE PRECISION (N)           
C        : Second coordinates of the points in dataset 
C WORKSPACE
C BETA   : DOUBLE PRECISION (N)
C        : Vector being used to store angles
C F      : INTEGER (N)
C        : Vector being used in to count
C DPF    : DOUBLE PRECISION (N)
C        : Dummy-Vector being used in sorting
C JLV    : INTEGER (N)
C        : Dummy-Vector being used in sorting
C JRV    : INTEGER (N)
C        : Dummy-Vector being used in sorting  
C OUTPUT
C HDEP   : INTEGER (1)
C        : Halfspace depth of the point (U,V) 
C SDEP   : DOUBLE PRECISION (1)
C        : Simplicial depth of the point (U,V)  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
C     INPUT
      INTEGER N
      DOUBLE PRECISION U,V,X(N),Y(N),EPSI
            
C     WORKSPACE
      DOUBLE PRECISION BETA(N)
      INTEGER F(N),DummyVector(N)
      DOUBLE PRECISION DPF(N)
      INTEGER JLV(N),JRV(N)
      
C     OUTPUT
      INTEGER HDEP
      DOUBLE PRECISION SDEP

C     INTERNAL 
      DOUBLE PRECISION P,P2,D,XU,YU,ANGLE,ALPHK,BETAK
      INTEGER GI, NZ, NN, NU, JA, JB, NN2, NF, KI
      INTEGER I,J
      INTEGER(8) K
      INTEGER NUMH
      INTEGER(8) NUMS
      INTEGER(8) NBAD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Initialise constants
      EPSI=1.D-8
      P=DACOS(DBLE(-1.0))
      P2=P*2.0

C NZ keeps track of points coinciding with theta
      NZ=0
      NUMH=0
      NUMS=0
      HDEP=0.0
      SDEP=0.0
C
      IF (N.LT.1) RETURN
C Start the algorithm
C  Construct the array BETA.
      DO 10 I=1,N
          D=DSQRT((X(I)-U)*(X(I)-U)+(Y(I)-V)*(Y(I)-V))
C If distance between theta and point (x_i,y_i) is smaller
C than the tollerance level, consider them equal
          IF (D.LE.EPSI) THEN
              NZ=NZ+1
          ELSE
C Normalise the point
              XU=(X(I)-U)/D
              YU=(Y(I)-V)/D
C Calculate angle of normalised point Theta-(X_i,Y_i)
              IF (DABS(XU).GT.DABS(YU)) THEN
                  IF (X(I).GE.U) THEN
                      BETA(I-NZ)=DASIN(YU)
                      IF(BETA(I-NZ).LT.0.0) THEN
                          BETA(I-NZ)=P2+BETA(I-NZ)
                      ENDIF
                  ELSE
                      BETA(I-NZ)=P-DASIN(YU)
                  ENDIF
              ELSE
                  IF (Y(I).GE.V) THEN
                      BETA(I-NZ)=DACOS(XU)
                  ELSE
                      BETA(I-NZ)=P2-DACOS(XU)
                  ENDIF
              ENDIF
C If angle is close enought to 2*pi, consider it 2*pi.
              IF (BETA(I-NZ).GE.(P2-EPSI)) BETA(I-NZ)=0.0
          ENDIF
  10  CONTINUE
C The amount of points different from theta is stored in NN
      NN=N-NZ
C If all but one point coincides with theta,
C the depth of theta is known and equal to NZ.
      IF (NN.LE.1) GOTO 60
C  Sort the array BETA.
      DO 15 I=1,NN
      DPF(I)=DBLE(F(I))
15    CONTINUE
      CALL SORT(BETA,F,DummyVector,DPF,NN,JLV,JRV)
C  Check whether Theta=(U,V) lies outside the data cloud.
      ANGLE=BETA(1)-BETA(NN)+P2
      DO 20 I=2,NN
          ANGLE=DMAX1(ANGLE,(BETA(I)-BETA(I-1)))
  20  CONTINUE
      IF (ANGLE.GT.(P+EPSI)) GOTO 60
C  Make smallest BETA equal to zero,
C  and compute NU = number of BETA < PI.
      ANGLE=BETA(1)
      NU=0
      DO 30 I=1,NN
          BETA(I)=BETA(I)-ANGLE
          IF (BETA(I).LT.(P-EPSI)) NU=NU+1
  30  CONTINUE
      IF (NU.GE.NN) GOTO 60
C  Mergesort the BETA with their antipodal angles,
C  and at the same time update I, F(I), and NBAD.
      JA=1
      JB=1
      ALPHK=BETA(1)
      BETAK=BETA(NU+1)-P
      NN2=NN*2
      NBAD=0
      I=NU
      NF=NN
      DO 40 J=1,NN2
          IF ((ALPHK+EPSI).LT.BETAK) THEN
              NF=NF+1
              IF (JA.LT.NN) THEN
                  JA=JA+1
                  ALPHK=BETA(JA)
              ELSE
                  ALPHK=P2+1.0
              ENDIF
          ELSE
              I=I+1
              IF (I.EQ.(NN+1)) THEN
                  I=1
                  NF=NF-NN
              ENDIF
              F(I)=NF
              NBAD=NBAD+K((NF-I),2)
              IF (JB.LT.NN) THEN
                  JB=JB+1
                  IF ((JB+NU).LE.NN) THEN
                      BETAK=BETA(JB+NU)-P
                  ELSE
                      BETAK=BETA(JB+NU-NN)+P
                  ENDIF
              ELSE
                  BETAK=P2+1.0
              ENDIF
          ENDIF
  40  CONTINUE
      NUMS=K(NN,3)-NBAD
C  Computation of NUMH for halfspace depth.
      GI=0
      JA=1
      ANGLE=BETA(1)
      NUMH=MIN0(F(1),(NN-F(1)))
      DO 50 I=2,NN
          IF(BETA(I).LE.(ANGLE+EPSI)) THEN
              JA=JA+1
          ELSE
              GI=GI+JA
              JA=1
              ANGLE=BETA(I)
          ENDIF
          KI=F(I)-GI
          NUMH=MIN0(NUMH,MIN0(KI,(NN-KI)))
   50 CONTINUE
C  Adjust for the number NZ of data points equal to theta=(U,V):
   60 NUMS=NUMS + K(NZ,1)*K(NN,2)+K(NZ,2)*K(NN,1) + K(NZ,3)
      IF(N.GE.3) SDEP = (NUMS+0.0)/(K(N,3)+0.0)
      NUMH=NUMH+NZ
      HDEP=NUMH
C     End the subroutine

      RETURN
      END





CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      FUNCTION K(M,J)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Returns the value zero if M<=J; otherwise computes the
C number of combinations of J out of M.
C Only for M equal to 1,2 or 3.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER(8) K
      INTEGER M,J
      INTEGER(8) MM
C
      K=0
      IF (M.LT.J) THEN
          K=0
      ELSE
          MM=int8(M)
          IF (J.EQ.1) K=MM
          IF (J.EQ.2) K=(MM*(MM-1))/2
          IF (J.EQ.3) K=(MM*(MM-1)*(MM-2))/6
      ENDIF
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
