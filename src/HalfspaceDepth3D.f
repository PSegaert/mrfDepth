      SUBROUTINE HSDEP3(U,V,W,N1,X,Y,Z,N2,HDEP,ERR)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is a wrapper to calculate the halfspace depths 
C of a set of bivariate points (U,V,W) given a certain dataset (X,Y,Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C INPUT
C U      : DOUBLE PRECISION (N1)            
C        : First coordinates of the points to calculate depth 
C V      : DOUBLE PRECISION (N1)           
C        : Second coordinates of the points to calculate depth 
C W      : DOUBLE PRECISION (N1)           
C        : Third coordinates of the points to calculate depth 
C N1     : INTEGER (1)
C        : Number of points to calculate depth of
C X      : DOUBLE PRECISION (N2)            
C        : First coordinates of the points in dataset 
C Y      : DOUBLE PRECISION (N2)           
C        : Second coordinates of the points in dataset 
C Z      : DOUBLE PRECISION (N2)           
C        : Third coordinates of the points in dataset 
C N2     : INTEGER (1)
C        : Number of points in the dataset
C EPS    : DOUBLE PRECISION (1)
C        : Numerical precision
C OUTPUT
C HDEP   : REAL(8) (N1)
C        : Halfspace depth of the points (U,V)	
C ERR    : REAL(8) (N1)
C        : Diagnostic Information									
C          1 Indicates that all points and theta lie on the same line
C          2 Indicates that all points and theta lie on the same plane	
C          3 Indicates that all points and theta do not lie on a subspace	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
C     INPUT
      INTEGER N1,N2
      DOUBLE PRECISION U(N1),V(N1),W(N1),X(N2),Y(N2),Z(N2),EPS
      DOUBLE PRECISION X2(N2),Y2(N2),Z2(N2)
C     OUTPUT
      REAL(8) HDEP(N1), ERR(N1)

C     WORKSPACE
      DOUBLE PRECISION ALPHA(N2)
      DOUBLE PRECISION XN(N2),YN(N2)
      INTEGER I,J, F(N2),NDIM,NDEP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      EPS=1.D-8
        
      DO 10 I=1,N1
C For each variable calculate (X-med)/mad.
         DO 124 j=1,N2
            X2(J)=X(J)
            Y2(J)=Y(J)
            Z2(J)=Z(J)
  124    CONTINUE
          CALL STANDHSDEP3(N2,X2,Y2,Z2,U(I),V(I),W(I),XN,EPS)
C Calculate the three dimensional Tukey Depth of the points.
          CALL HSDEPTH31(N2,U(I),V(I),W(I),
     +    X2,Y2,Z2,ALPHA,F,XN,YN,EPS,NDIM,NDEP)
          HDEP(I)=(NDEP+0.0)/(N2+0.0)
C Set the corresponding flag.
          IF (NDIM.EQ.2) THEN 
            ERR(I) = 2
          ELSE IF (NDIM.EQ.1) THEN 
            ERR(I) = 1
          ELSE IF (NDIM.EQ.3) THEN 
            ERR(I) = 3
          ELSE 
            ERR(I) = -1
          END IF
10    CONTINUE
      
      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE HSDEPTH31(N,UU,VV,WW,X,Y,Z,ALPHA,F,XN,YN,EPS,NDIM,NDEP)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine calculates the halfspace depth of a point (U,V,W)  
C given a certain dataset (X,Y,Z).
C  This subroutine was described in:
C       Rousseeuw, P.J. and Struyf, A. (1998):  
C       Computing location depth and regression depth in higher dimensions,
C       Statistics and Computing Volume 8, pp 193-203.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C INPUT
C N      : INTEGER (1)
C        : Number of points in the dataset
C U      : DOUBLE PRECISION (1)            
C        : First coordinate of the point to calculate depth 
C V      : DOUBLE PRECISION (1)           
C        : Second coordinate of the point to calculate depth 
C W      : DOUBLE PRECISION (1)           
C        : Third coordinate of the point to calculate depth 
C X      : DOUBLE PRECISION (N)            
C        : First coordinates of the points in dataset 
C Y      : DOUBLE PRECISION (N)           
C        : Second coordinates of the points in dataset 
C Z      : DOUBLE PRECISION (N)           
C        : Third coordinates of the points in dataset 
C EPS    : DOUBLE PRECISION (1)
C        : Numerical precision
C OUTPUT
C NDIM   : INTEGER (1)
C        : Diagnostic Information							
C          0 Indicates no flag						
C          1 Indicates that all points and theta lie on the same line
C          2 Indicates that all points and theta lie on the same plane
C NDEP   : INTEGER (1)
C        : Halfspace depth of the point (U,V) multiplied by N.
C Workspace
C ALPHA  : DOUBLE PRECISION (N)
C        : Used to store angles in the algorithm
C F      : INTEGER (N) 
C        : Used to count in the algorithm
C XN     : DOUBLE PRECISION (N)
C        : Vector used in the algorithm for projected coordinates
C YN     : DOUBLE PRECISION (N)
C        : Vector used in the algorithm for projected coordinates
C	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C INPUT
      INTEGER N
      DOUBLE PRECISION U,V,W,X(N),Y(N),Z(N),EPS
      DOUBLE PRECISION UU,VV,WW

      DOUBLE PRECISION xmean,ymean,zmean,xdev,ydev,zdev

C OUTPUT
      INTEGER NDEP,NDIM
C WORKSPACE
      DOUBLE PRECISION ALPHA(N),XN(N),YN(N)
      DOUBLE PRECISION A(2,3),B(3,2),DP
      INTEGER F(N),NTNUL,NTPOS,NTNEG
      INTEGER I,J,NH
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine

C
C     Standardization of the data set.
C
      xmean=0.0
      ymean=0.0
      zmean=0.0
      xdev=0.0
      ydev=0.0
      zdev=0.0
      DO 2 I=1,N
         xmean=xmean+x(i)
         ymean=ymean+y(i)
         zmean=zmean+z(i)
2     CONTINUE
      xmean=xmean/n
      ymean=ymean/n
      zmean=zmean/n
      DO 3 I=1,n
         xdev=xdev+((x(i)-xmean)*(x(i)-xmean))
         ydev=ydev+((y(i)-ymean)*(y(i)-ymean))
         zdev=zdev+((z(i)-zmean)*(z(i)-zmean))
3     continue
      xdev=xdev/(n-1)
      ydev=ydev/(n-1)
      zdev=zdev/(n-1)
      xdev=dsqrt(xdev)
      ydev=dsqrt(ydev)
      zdev=dsqrt(zdev)
      DO 5 I=1,N
         if (xdev.gt.eps) THEN
            x(i)=(x(i)-xmean)/xdev
            U=((UU-xmean)/xdev)
         endif 
         if (ydev.gt.eps) THEN
            y(i)=(y(i)-ymean)/ydev
            V=((VV-ymean)/ydev)
         endif
         if (ydev.gt.eps) THEN
            z(i)=(z(i)-zmean)/zdev
            W=((WW-zmean)/zdev)
         endif         
5     CONTINUE

C Make theta the center of the dataset.
      DO 10 I=1,N
         X(I)=X(I)-U
         Y(I)=Y(I)-V
         Z(I)=Z(I)-W
 10   CONTINUE
      NDIM=3
C Handle special cases where N is less or equal to 1.
      IF (N.LE.1) THEN
         IF ((N.EQ.1) .AND. (DABS(X(1)).LE.EPS) .AND. 
     +        (DABS(Y(1)).LE.EPS) .AND. 
     +        (DABS(Z(1)).LE.EPS)) THEN
            NDEP=1
         ELSE
            NDEP=0
         ENDIF
         RETURN
      ENDIF
C General case: initialize halfspace depth.
      NDEP=N
C Loop over all lines (theta,x(i)).
      DO 20 I=1,N
         IF ((DABS(X(I)).LE.EPS).AND.(DABS(Y(I)).LE.EPS).AND.
     +        (DABS(Z(I)).LE.EPS)) THEN
            GOTO 20
         ENDIF
C Calculate the matrix of the orthogonal projection on the plane through 
C theta, orthogonal to the line through theta and x(i).
C Let the third coordinate coincide with the line (theta,x(i)).
         IF (DABS(X(I)).GT.EPS) THEN
            B(2,1)=1.D0
            B(3,1)=1.D0
            B(1,1)=-(Y(I)+Z(I))/X(I)
         ELSEIF (DABS(Y(I)).GT.EPS) THEN
            B(1,1)=1.D0
            B(3,1)=1.D0
            B(2,1)=-(X(I)+Z(I))/Y(I)
         ELSE
            B(1,1)=1.D0
            B(2,1)=1.D0
            B(3,1)=-(X(I)+Y(I))/Z(I)
         ENDIF
         B(1,2)=B(2,1)*Z(I)-B(3,1)*Y(I)
         B(2,2)=B(3,1)*X(I)-B(1,1)*Z(I)
         B(3,2)=B(1,1)*Y(I)-X(I)*B(2,1)   
         A(1,1)=(B(2,2)*Z(I)-Y(I)*B(3,2))
         A(1,2)=-(B(1,2)*Z(I)-X(I)*B(3,2))
         A(1,3)=(B(1,2)*Y(I)-B(2,2)*X(I))
         A(2,1)=-(B(2,1)*Z(I)-Y(I)*B(3,1))
         A(2,2)=(B(1,1)*Z(I)-X(I)*B(3,1))
         A(2,3)=-(B(1,1)*Y(I)-X(I)*B(2,1))
C Compute the new planar coordinates for all points.
C If a point collapses with theta, identify its position:
C     NTNUL = real ties, 
C     NTPOS = the original point lies on the positive side of 
C             the projection plane,
C     NTNEG = the original point lies on the negative side of 
C             the projection plane.
         NTNUL=0
         NTPOS=0
         NTNEG=0
         DO 30 J=1,N
            XN(J)=X(J)*A(1,1)+Y(J)*A(1,2)+Z(J)*A(1,3)
            YN(J)=X(J)*A(2,1)+Y(J)*A(2,2)+Z(J)*A(2,3)
            IF ((DABS(XN(J)).LE.EPS).AND.
     +           (DABS(YN(J)).LE.EPS)) THEN
               DP=X(J)*X(I)+Y(J)*Y(I)+Z(J)*Z(I)
               IF (DABS(DP).LE.EPS) THEN
                  NTNUL=NTNUL+1
               ELSEIF (DP.GT.EPS) THEN 
                  NTPOS=NTPOS+1
               ELSE 
                  NTNEG=NTNEG+1
               ENDIF
            ENDIF
 30      CONTINUE
         IF ((NTNUL+NTNEG+NTPOS).EQ.N) GOTO 50
C Compute the halfspace depth in two dimensions.
         CALL HSDEPTH31B(0.D0,0.D0,N,XN,YN,ALPHA,F,NH,NTPOS,NTNEG,
     +        NTNUL,EPS,NDIM)
C Update the three-dimensional halfspace depth.
         NDEP=MIN0(NDEP,NH)
 20   CONTINUE
      RETURN
C All points and theta lie on one line.
 50   NDEP=MIN0(NTNUL+NTPOS,NTNUL+NTNEG)
      NDIM=1

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE HSDEPTH31B(U,V,N,X,Y,
     +ALPHA,F,NH,NTPOS,NTNEG,NTNUL,EPS,NDIM)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine calculates is used to calculate the halfspace depth of
C a point (U,V,W)  given a certain dataset (X,Y,Z). Specifically it 
C calculates the depth of the points projected in two dimensions. 
C The two dimesional algorithm can not be used since we have to keep
C track of points that are equal after projection. 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C INPUT
C U      : DOUBLE PRECISION (1)            
C        : First coordinate of the projected point to calculate depth 
C V      : DOUBLE PRECISION (1)           
C        : Second coordinate of the projected point to calculate depth 
C N      : INTEGER (1)
C        : Number of points in the dataset
C X      : DOUBLE PRECISION (N)            
C        : First coordinates of the projected points in dataset 
C Y      : DOUBLE PRECISION (N)           
C        : Second coordinates of the projected points in dataset 
C NTNUL  : INTEGER (1)
C        : real ties 
C NTPOS  : INTEGER (1)
C        : the original point lies on the positive side of the projection plane
C NTNEG  : INTEGER (1)
C        : the original point lies on the negative side of the projection plane
C EPS    : DOUBLE PRECISION (1)
C        : Numerical precision
C OUTPUT
C NDIM   : INTEGER (1)
C        : Diagnostic Information							
C          0 Indicates no flag						
C          1 Indicates that all points and theta lie on the same line
C          2 Indicates that all points and theta lie on the same plane
C Workspace
C ALPHA  : DOUBLE PRECISION (N)
C        : Used in the algorithm
C F      : INTEGER (N) 
C        : Used to count in the algorithm
C NH     : INTEGER (1)
C        : Used in the algorithm 
C	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C INPUT
      INTEGER N
      DOUBLE PRECISION U,V,X(N),Y(N),EPS
      INTEGER NTNUL,NTPOS,NTNEG
C OUTPUT
      INTEGER NDIM      
C WORKSPACE      
      DOUBLE PRECISION ALPHA(N),P,P2,D,XU,YU,ANGLE,ALPHK,BETAK
      INTEGER F(N)
      INTEGER GI,NH,NUMH,NT,ND,I,NN,NU,JA,JB,NN2,NF,J,KI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine

      NUMH=0
      NH=0
      IF (N.LT.1) RETURN
      P=DACOS(-1.D0)
      P2=P*2.D0
      NT=0
      ND=0
C  Construct the array ALPHA.
      DO 10 I=1,N
          D=DSQRT((X(I)-U)*(X(I)-U)+(Y(I)-V)*(Y(I)-V))
          IF (D.LE.EPS) THEN
              NT=NT+1
          ELSE
              XU=(X(I)-U)/D
              YU=(Y(I)-V)/D
              IF (DABS(XU).GT.DABS(YU)) THEN
                  IF (X(I).GE.U) THEN
                      ALPHA(I-NT)=DASIN(YU)
                      IF(ALPHA(I-NT).LT.0.D0) THEN
                          ALPHA(I-NT)=P2+ALPHA(I-NT)
                      ENDIF
                  ELSE
                      ALPHA(I-NT)=P-DASIN(YU)
                  ENDIF
              ELSE
                  IF (Y(I).GE.V) THEN
                      ALPHA(I-NT)=DACOS(XU)
                  ELSE
                      ALPHA(I-NT)=P2-DACOS(XU)
                  ENDIF
              ENDIF
              IF (ALPHA(I-NT).GE.(P2-EPS)) ALPHA(I-NT)=0.D0
          ENDIF
  10  CONTINUE
      NN=N-NT
      IF (NN.LE.1) THEN
         NUMH=MIN0(NTNEG,NTPOS)
         GOTO 60
      ENDIF
C  Sort the array ALPHA.
      CALL SORT2(ALPHA,NN)
C  Check whether theta=(U,V) lies outside the data cloud.
      ANGLE=ALPHA(1)-ALPHA(NN)+P2
      DO 20 I=2,NN
          ANGLE=DMAX1(ANGLE,(ALPHA(I)-ALPHA(I-1)))
  20  CONTINUE
      IF (ANGLE.GT.(P+EPS)) THEN
         NUMH=MIN0(NTNEG,NTPOS)
         GOTO 60
      ENDIF
C  Make smallest alpha equal to zero,
C  and compute NU = number of alpha < pi.
      ANGLE=ALPHA(1)
      NU=0
      DO 30 I=1,NN
          ALPHA(I)=ALPHA(I)-ANGLE
          IF (ALPHA(I).LT.(P-EPS)) NU=NU+1
          IF ((DABS(ALPHA(I)).LE.EPS).OR.
     +         (DABS(ALPHA(I)-P).LE.EPS)) ND=ND+1
  30  CONTINUE
      IF (ND.EQ.NN) NDIM=2
      IF (NU.GE.NN) THEN
         NUMH=MIN0(NTNEG,NTPOS)
         GOTO 60
      ENDIF
C  Mergesort the alpha with their antipodal angles beta,
C  and at the same time update I, and F(I).
      JA=1
      JB=1
      ALPHK=ALPHA(1)
      BETAK=ALPHA(NU+1)-P
      NN2=NN*2
      I=NU
      NF=NN
      DO 40 J=1,NN2
          IF ((ALPHK+EPS).LT.BETAK) THEN
              NF=NF+1
              IF (JA.LT.NN) THEN
                  JA=JA+1
                  ALPHK=ALPHA(JA)
              ELSE
                  ALPHK=P2+1.D0
              ENDIF
          ELSE
              I=I+1
              IF (I.EQ.(NN+1)) THEN
                  I=1
                  NF=NF-NN
              ENDIF
              F(I)=NF
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
          ENDIF
  40  CONTINUE
C  Compute the halfspace depth.
C  Correct for ties (which where no ties in three dimensions)
C  by considering some small deviations of the planes (without changing
C  their intersection line with the plane which contains the 2D-data).
      GI=0
      JA=1
      ANGLE=ALPHA(1)
      NUMH=MIN0(MIN0(F(1)+NTNEG,F(1)+NTPOS),
     +     MIN0(NN-F(1)+NTNEG,NN-F(1)+NTPOS))
      DO 50 I=2,NN
         IF(ALPHA(I).LE.(ANGLE+EPS)) THEN
            JA=JA+1
         ELSE
            GI=GI+JA
            JA=1
            ANGLE=ALPHA(I)
         ENDIF
         KI=F(I)-GI
         NUMH=MIN0(NUMH,MIN0(MIN0(KI+NTNEG,KI+NTPOS),
     +        MIN0(NN-KI+NTNEG,NN-KI+NTPOS)))
 50   CONTINUE
C  Adjust for the number NTNUL of data points equal to theta.
 60   NH=NUMH+NTNUL
      RETURN
      END







CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







      SUBROUTINE SORT2(B,N)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Sorts an array B of length N in O(NlogN) time.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C N      : INTEGER (1)
C        : Number of elements in B
C B      : DOUBLE PRECISION (N)
C        : The vector to be sorted
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variable
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Input
      INTEGER N
      DOUBLE PRECISION B(N)
C Workspace
      DOUBLE PRECISION AMM,XX
      INTEGER JLV(N),JRV(N)
      INTEGER JSS,JNC,J,JTWE,JR,JNDL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
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







CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







      SUBROUTINE STANDHSDEP3(N,X,Y,Z,U,V,W,XN,EPS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Center the data by calculating (X-med)/mad for every variable
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C N      : INTEGER (1)
C        : Number of points in the dataset
C X      : DOUBLE PRECISION (N)            
C        : First coordinates of the points in dataset 
C Y      : DOUBLE PRECISION (N)           
C        : Second coordinates of the points in dataset 
C Z      : DOUBLE PRECISION (N)           
C        : Third coordinates of the points in dataset 
C U      : DOUBLE PRECISION (1)            
C        : First coordinate of the point to calculate depth 
C V      : DOUBLE PRECISION (1)           
C        : Second coordinate of the point to calculate depth 
C W      : DOUBLE PRECISION (1)           
C        : Third coordinate of the point to calculate depth
C XN     : DOUBLE PRECISION (N)
C        : Vector used in the algorithm 
C EPS    : DOUBLE PRECISION (1)
C        : Numerical precision
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variable
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C INPUT
      INTEGER N
      DOUBLE PRECISION X(N),Y(N),Z(N),U,V,W,XN(N),EPS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      CALL STANDHSDEP31(N,X,U,XN,EPS)
      CALL STANDHSDEP31(N,Y,V,XN,EPS)
      CALL STANDHSDEP31(N,Z,W,XN,EPS)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







      SUBROUTINE STANDHSDEP31(N,X,U,XN,EPS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Center the data by calculating (X-med)/mad for every variable
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C N      : INTEGER (1)
C        : Number of points in the dataset
C X      : DOUBLE PRECISION (N)            
C        : Coordinates to be centered
C U      : DOUBLE PRECISION (1)            
C        : First coordinate of the point to calculate depth 
C XN     : DOUBLE PRECISION (N)
C        : Vector used in the algorithm 
C EPS    : DOUBLE PRECISION (1)
C        : Numerical precision
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variable
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C INPUT
      INTEGER N
      DOUBLE PRECISION X(N),U,XN(N),EPS
C WORKSPACE
      DOUBLE PRECISION QLOC,QSCA,AVE,VAR,FINDQ
      INTEGER I
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
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
      U=(U-QLOC)/QSCA

      RETURN
      END
