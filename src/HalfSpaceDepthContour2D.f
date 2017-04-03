      subroutine iso2hdw(XX,YY,N,K,XCont,YCont,
     +                   empty,kount,dithered)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This subroutine calculates the isodepth contours
C for two dimensional dataset (X,Y) based on the halfspace depth.
C Missing values are not allowed.
C  This subroutine was described in:
C       Ruts, I. and Rousseeuw, P.J.:
C       Computing depth contours of bivariate point clouds
C       Computational Statistics and Data Analysis 23, 153-168 (1996).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C X        : DOUBLE PRECISION (N)
C          : First variable of  the data set.
C Y        : DOUBLE PRECISION (N)
C          : Second variable of  the data set.
C N        : INTEGER(1)
C          : The number of points in the data set
C K        : INTEGER (1)
C          : Depth of the contour to be computed
C          : Format is 1/n,2/n,...,n/n
C XCont    : DOUBLE PRECISION (N*(N-1)/2)
C          : First coordinates of the vertices in the contour.
C YCont    : DOUBLE PRECISION (N*(N-1)/2)
C          : Second coordinates of the vertices in the contour.
C Empty    : INTEGER (1)
C          : 1 = Contour is empty
C          : 0 = Contour is non-empty
C Kount    : INTEGER (1)
C          : The number of vertices in the contour.
C dithered : INTEGER (1)
C          : 1 = Dithering was applied
C          : 0 = No dithering was applied
C EPS      : DOUBLE PRECISION (1)
C          : Numerical Tolerance (10^-6)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
C INPUT
      INTEGER N, K
      DOUBLE PRECISION X(N),Y(N),XX(N),YY(N)
      DOUBLE PRECISION EPS

C     OUTPUT
      DOUBLE PRECISION XCont(N*(N-1)/2), YCont(N*(N-1)/2)
      INTEGER empty,KOUNT, dithered

C     WORKSPACE
      INTEGER I,M,MAXNUM,NUM,NCIRQ(N),MCIRQ(N)
      INTEGER JLV(N), JRV(N), IND1(N*(N-1)/2), IND2(N*(N-1)/2)
      DOUBLE PRECISION ANGLE(N*(N-1)/2)
      INTEGER KORNR(int(4*N*SQRT(REAL(N))+1),4)
      DOUBLE PRECISION fac
      INTEGER moredith,dithMemory,seed
      DOUBLE PRECISION xmean,ymean,xdev,ydev



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      M=((N*(N-1))/2)
      MAXNUM=int(4*N*SQRT(REAL(N))+1)
      fac=100000.0
      moredith=0
      dithMemory=0
      seed=256
      EPS=1.D-13

C
C     Standardization of the data set.
C
      xmean=0.0
      ymean=0.0
      xdev=0.0
      ydev=0.0
      DO 2 I=1,N
         xmean=xmean+xx(i)
         ymean=ymean+yy(i)
2     CONTINUE
      xmean=xmean/n
      ymean=ymean/n
      DO 3 I=1,N
         xdev=xdev+((xx(i)-xmean)*(xx(i)-xmean))
         ydev=ydev+((yy(i)-ymean)*(yy(i)-ymean))
3     continue
      xdev=xdev/(N-1)
      ydev=ydev/(N-1)
      xdev=dsqrt(xdev)
      ydev=dsqrt(ydev)
      DO 5 I=1,N
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

C Check whether dithering is needed
C Sort the data and compute the angles
10    CALL checkData2D(X,Y,N,fac,NCIRQ,MCIRQ,ANGLE,
     +JLV,JRV,IND1,IND2,dithered)

C The cases n<=3 are simple
      if ((n.eq.1).and.(k.eq.1)) goto 41
41    if (n.le.3) then
         do 42 i=1,n
         XCont(i)=x(i)
         YCont(i)=y(i)
42       continue
         EMPTY=0
         kount=n
         goto 200
      endif

      CALL isoFin98(X,Y,N,K,
     + NCIRQ,MCIRQ,JLV,JRV,IND1,IND2,ANGLE,KORNR,EMPTY,NUM,EPS)

C  Scan KORNR and store the coordinates of the vertices.
      KOUNT=0
      if (EMPTY.eq.1) then
          goto 200
      else
          CALL fillCont(X,Y,N,KORNR,MAXNUM,XCont,YCont,KOUNT,NUM,EPS)
      endif

      RETURN
200   end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      subroutine checkData2D(X,Y,N,fac,NCIRQ,MCIRQ,ANGLE,
     +JLV,JRV,IND1,IND2,dithered)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is to be called by the iso2hdw algorithm.
C It checks whether dithering is needed, sorts the points according
C to their x-values and computes the angles between the points.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER N
      DOUBLE PRECISION X(N), Y(N)
      DOUBLE PRECISION WX(N), WY(N), WX1(N), WY1(N)
      DOUBLE PRECISION ANGLE(N*(N-1)/2),D(N*(N-1)/2), ANG1
      INTEGER NCIRQ(N),MCIRQ(N), JLV(N), JRV(N)
      INTEGER IND1(N*(N-1)/2), IND2(N*(N-1)/2)
      DOUBLE PRECISION rand(2)
      INTEGER I,J,IV,L,M, LEFT
      INTEGER dithered

C     CONSTANTS
      INTEGER seed
      DOUBLE PRECISION PI,PI2,fac
      SEED=256
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      DO 5 I=1,N
         WX(I)=X(I)
         WY(I)=Y(I)
         wx1(i)=x(i)
         wy1(i)=y(i)
         NCIRQ(I)=I
         MCIRQ(I)=I
 5    CONTINUE
C  To test whether the data points are in general position,
C  we first check whether any two points coincide.
C  When we sort the points with respect to the x-coordinate, we
C  permute NCIRQ and MCIRQ in the same way, to initialize them.
      dithered = 0
1     CALL SORT(WX,NCIRQ,MCIRQ,WY,N,JLV,JRV)
      I=0
      J=0
10    I=I+1
      IF (I+1.GT.N) GOTO 15
      J=I+1
      IF (WX(I).NE.WX(J)) GOTO 10
      IF (WY(I).EQ.WY(J)) THEN
C         The points ',NCIRQ(I),' and ',NCIRQ(J),' coincide
         GOTO 210
      ELSE
         if (wy(i).lt.wy(j)) then
            iv=mcirq(j)
            mcirq(j)=mcirq(i)
            mcirq(i)=iv
         endif
         IF (J+1.LE.N) THEN
            L=J+1
            IF (WX(I).EQ.WX(L)) THEN
C              The points ',NCIRQ(I),', ',NCIRQ(J),' and 'NCIRQ(L),' lie on a vertical line.
               GOTO 210
            ELSE
               GOTO 10
            ENDIF
         ENDIF
      ENDIF
C  Compute all the angles formed by pairs of data points.
15    M=((N*(N-1))/2)
      L=1
      DO 20 I=1,N
      DO 25 J=I+1,N
      IF (X(I).EQ.X(J)) THEN
         ANGLE(L)=PI2
            ELSE
         ANGLE(L)=DATAN((Y(I)-Y(J))/(X(I)-X(J)))
         IF (ANGLE(L).le.0.0) ANGLE(L)=ANGLE(L)+PI
            ENDIF
      IND1(L)=I
      IND2(L)=J
      L=L+1
25    CONTINUE
20    CONTINUE
C  Sort all the angles and permute IND1 and IND2 in the same way.
C  To avoid using several SORT-routines, we will always permute
C  two integer arrays and one real array.
      CALL SORT(ANGLE,IND1,IND2,D,M,JLV,JRV)
C  Test whether any three points are collinear
      LEFT=1
30    ANG1=ANGLE(LEFT)
      DO 35 J=LEFT+1,M
      IF (ANGLE(J).GT.ANG1) THEN
      LEFT=J
      GOTO 30
         ELSE
            DO 36 I=LEFT,J-1
               IF ((IND1(I).EQ.IND1(J)).or.
     +            (IND1(I).EQ.IND2(J))) THEN
C The data are not in general position:' The points',IND1(J),', ',IND2(J),' and ',IND2(I),' are collinear.'
         GOTO 210
         ENDIF
         IF ((IND2(I).EQ.IND1(J)).or.
     +            (IND2(I).EQ.IND2(J))) THEN
C The data are not in general position:' The points',IND1(J),', ',IND2(J),' and ',IND1(I),' are collinear.'
         GOTO 210
         ENDIF
36          CONTINUE
         ENDIF
35    CONTINUE
      goto 37
C If the data are not in general position use dithering.
C If not end the routine
210   dithered=1
      do 211 i=1,n
         ncirq(i)=i
         mcirq(i)=i
         call norrandp(2,seed,rand)
         wx(i)=wx1(i)+rand(1)/fac
         x(i)=wx(i)
         wy(i)=wy1(i)+rand(2)/fac
         y(i)=wy(i)
211   continue
      goto 1

37    RETURN
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      subroutine isoFin98(X,Y,N,K,
     +NCIRQ,MCIRQ,JLV,JRV,IND1,IND2,ANGLE,KORNR,EMPTY,NUM,EPS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is to be called by the iso2hdw algorithm.
C It searches for the special k-dividers and find the candidate
C vertices for the contour.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER N,K
      DOUBLE PRECISION X(N),Y(N)
      INTEGER KOUNT
      INTEGER EMPTY
C     WORKSPACE
      INTEGER MAXNUM
      INTEGER NCIRQ(N),MCIRQ(N),NRANK(N),F(N)
      INTEGER JLV(N),JRV(N)
      INTEGER IND1(N*(N-1)/2),IND2(N*(N-1)/2)
      INTEGER KAND1(int(4*N*SQRT(REAL(N))+1))
      INTEGER KAND2(int(4*N*SQRT(REAL(N))+1))
      INTEGER KORNR(int(4*N*SQRT(REAL(N))+1),4)
      INTEGER KON,KONTROL,NDATA,NDK,HALT,halt2,jj,JFULL
      INTEGER IV,IW1,IW2,NEXT,JFLAG,NUM
      INTEGER HDEP1,HDEP2,HDEP3,HDEP4,HDEP5,I,J,L,M
      DOUBLE PRECISION BETA(N)
      DOUBLE PRECISION ANGLE(N*(N-1)/2),DPF(N)
      DOUBLE PRECISION ALPHA(int(4*N*SQRT(REAL(N))+1))
      DOUBLE PRECISION D(int(4*N*SQRT(REAL(N))+1))
      DOUBLE PRECISION xcord1,ycord1
      DOUBLE PRECISION XCORD,YCORD,ANG1,m1,m2
      DOUBLE PRECISION SDEP
C     CONSTANTS
      DOUBLE PRECISION PI,PI2,EPS
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
C     Define for use in the further algorithm
      M=N*(N-1)/2
      MAXNUM=int(4*N*SQRT(REAL(N))+1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
C
C  (Re)initialize NCIRQ and NRANK
C
      DO 45 I=1,N
      NCIRQ(I)=MCIRQ(I)
45    CONTINUE
      DO 50 I=1,N
      IV=NCIRQ(I)
      NRANK(IV)=I
50    CONTINUE
C
C  Let the line rotate from zero to ANGLE(1)
C
      KOUNT=1
      HALT=0
      if (angle(1).gt.pi2) then
         L=1
      CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +                    K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
         halt=1
      endif
      L=2
 60   KONTROL=0
C Check if line 2 is perpendicular to a projectionline with
C angle between 0 and angle(1)
      IF ((PI.LE.(ANGLE(L)+PI2)).AND.((ANGLE(L)-PI2).LT.ANGLE(1))) THEN
      CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +                    K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
      KONTROL=1
      ENDIF
      L=L+1
      IF (KONTROL.EQ.1) HALT=1
      IF ((L.EQ.M+1).AND.(KONTROL.EQ.1)) THEN
       JFLAG=1
       GOTO 79
      ENDIF
      IF (((HALT.EQ.1).AND.(KONTROL.EQ.0)).OR.(L.EQ.M+1)) THEN
       GOTO 70
      ELSE
       GOTO 60
      ENDIF
 70   if (l.gt.1) then
         JFLAG=L-1
      else
         jflag=m
      endif
      J=0
C
C  In case the first switch didn't occur between zero and ANGLE(1),
C  look for it between the following angles.
C
      IF ((L.EQ.M+1).AND.(KONTROL.EQ.0)) THEN
       HALT=0
         halt2=0
 73      J=J+1
         if (j.eq.m+1) j=1
         L=J+1
         if (l.eq.m+1) l=1
 75      KONTROL=0
      IF ((ANGLE(L)+PI2).LT.PI) THEN
      ANG1=ANGLE(L)+PI2
         ELSE
      ANG1=ANGLE(L)-PI2
         ENDIF
         if (j.eq.m) then
            jj=1
            if (halt2.eq.0) angle(1)=angle(1)+pi
         else
            jj=j+1
         endif
       IF ((ANGLE(J).LE.ANG1).AND.(ANG1.LT.ANGLE(jj))) THEN
         if (angle(1).gt.pi) angle(1)=angle(1)-pi
      CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +                       K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
      KONTROL=1
         ENDIF
         if (angle(1).gt.pi) angle(1)=angle(1)-pi
      IF (L.NE.M) THEN
      L=L+1
         ELSE
      L=1
         ENDIF
      IF (KONTROL.EQ.1) HALT=1
      IF ((HALT.EQ.1).AND.(KONTROL.EQ.0)) THEN
            if (halt2.eq.1) goto 101
            if (l.gt.1) then
            JFLAG=L-1
            else
            jflag=m
            endif
      GOTO 79
         ELSE
            IF (L.EQ.jj) THEN
               if (jj.eq.1) halt2=1
         GOTO 73
            ELSE
               GOTO 75
            ENDIF
         ENDIF
      ENDIF
C
C  The first switch has occurred. Now start looking for the next ones,
C  between the following angles.
C
79    DO 80 I=J+1,M-1
      L=JFLAG
90      KONTROL=0
      IF ((ANGLE(L)+PI2).LT.PI) THEN
      ANG1=ANGLE(L)+PI2
         ELSE
      ANG1=ANGLE(L)-PI2
         ENDIF
      IF ((ANGLE(I).LE.ANG1).AND.(ANG1.LT.ANGLE(I+1))) THEN
      CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,
     +                  ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
      KONTROL=1
         ENDIF
      IF (KONTROL.EQ.0) THEN
      JFLAG=L
         ELSE
      IF (L.NE.M) THEN
         L=L+1
            ELSE
         L=1
            ENDIF
      GOTO 90
         ENDIF
 80   CONTINUE
      L=JFLAG
C
C  Finally, look for necessary switches between the last angle and zero.
C
100   KONTROL=0
      IF ((ANGLE(L)+PI2).LT.PI) THEN
      ANG1=ANGLE(L)+PI2
      ELSE
      ANG1=ANGLE(L)-PI2
      ENDIF
      IF ((ANGLE(M).LE.ANG1).AND.(ANG1.LT.PI)) THEN
      CALL ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,
     +               ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
      KONTROL=1
      ENDIF
      IF (KONTROL.EQ.1) THEN
         IF (L.NE.M) THEN
       L=L+1
         ELSE
       L=1
         ENDIF
         GOTO 100
      ENDIF
101      NUM=KOUNT-1
C
C  Sort the NUM special k-dividers.
C  Permute KAND1, KAND2 and D in the same way.
C
      CALL SORT(ALPHA,KAND1,KAND2,D,NUM,JLV,JRV)

      IW1=1
      IW2=2
      JFULL=0
      NDK=0
120   NDATA=0
C
C  Compute the intersection point.
C
      IF (DABS(-DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1))
     +         +DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2))).LT.EPS) THEN
      IW2=IW2+1
      IF (IW2.EQ.NUM+1) IW2=1
      GOTO 120
      ENDIF
      XCORD=(DCOS(ALPHA(IW2))*D(IW1)-DCOS(ALPHA(IW1))*D(IW2))
     + /(-DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1))
     +                   +DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2)))
      YCORD=(-DSIN(ALPHA(IW2))*D(IW1)+DSIN(ALPHA(IW1))*D(IW2))
     + /(-DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2))
     +                   +DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1)))
C
C  Test whether the intersection point is a data point.
C  If so, adjust IW1 and IW2.
C

      IF ((KAND1(IW1).EQ.KAND1(IW2)).OR.(KAND1(IW1).EQ.KAND2(IW2)))
     +     NDATA=KAND1(IW1)
      IF ((KAND2(IW1).EQ.KAND1(IW2)).OR.(KAND2(IW1).EQ.KAND2(IW2)))
     +     NDATA=KAND2(IW1)
      IF (NDATA.NE.0) THEN
        iv=0
125       NEXT=IW2+1
            iv=iv+1
      IF (NEXT.EQ.(NUM+1)) NEXT=1
          if (next.ne.iw1) then
          IF ((NDATA.EQ.KAND1(NEXT)).OR.(NDATA.EQ.KAND2(NEXT))) THEN
      IW2=IW2+1
      IF (IW2.EQ.(NUM+1)) IW2=1
      GOTO 125
          ENDIF
          endif
             if (iv.eq.(num-1)) then
             num=1
                  KORNR(1,1)=KAND1(IW1)
                  KORNR(1,2)=KAND2(IW1)
                  KORNR(1,3)=KAND1(IW2)
                  KORNR(1,4)=KAND2(IW2)
                kount=1
                goto 170
             endif
      ENDIF
      IF (IW2.EQ.NUM) THEN
      KON=1
      ELSE
      KON=IW2+1
      ENDIF
       if (kon.eq.iw1) kon=kon+1
       if (kon.eq.num+1) kon=1

C
C  Test whether the intersection point lies to the left of the special
C  k-divider which corresponds to ALPHA(KON). If so, compute its depth.
C
      IF ((DSIN(ALPHA(KON))*XCORD-DCOS(ALPHA(KON))*YCORD
     +     -D(KON)).LE.eps) THEN
      CALL HSDEP21(XCORD,YCORD,N,X,Y,
     +             BETA,F,DPF,JLV,JRV,HDEP1,SDEP)

      IF (HDEP1.EQ.K) NDK=1
      IF (HDEP1.NE.K) THEN
      CALL HSDEP21(XCORD-EPS*10,YCORD-EPS*10,N,X,Y,BETA,F,DPF,
     +         JLV,JRV,HDEP2,SDEP)
      CALL HSDEP21(XCORD+EPS*10,YCORD+EPS*10,N,X,Y,BETA,F,DPF,
     +         JLV,JRV,HDEP3,SDEP)
      CALL HSDEP21(XCORD-EPS*10,YCORD+EPS*10,N,X,Y,BETA,F,DPF,
     +         JLV,JRV,HDEP4,SDEP)
      CALL HSDEP21(XCORD+EPS*10,YCORD-EPS*10,N,X,Y,BETA,F,DPF,
     +         JLV,JRV,HDEP5,SDEP)

      IF ((NDK.EQ.0).AND.
     +    ((HDEP1.ge.K).OR.(HDEP2.ge.K).OR.(HDEP3.ge.K)
     +      .OR.(HDEP4.ge.K).OR.(HDEP5.ge.K))) THEN
      NDK=1
      ENDIF
      IF ((HDEP1.LT.K).AND.(HDEP2.LT.K)
     +   .AND.(HDEP3.LT.K).AND.(HDEP4.LT.K)
     +   .AND.(HDEP5.LT.K).AND.(NDK.EQ.1)) THEN
C
C  The intersection point is not the correct one,
C  try the next special k-divider.
C
      IW2=IW2+1
      IF (IW2.EQ.(NUM+1)) IW2=1
      GOTO 120
      ENDIF
      ENDIF
C
C  Store IW1 and IW2 in KORNR. If KORNR has already been filled, check whether
C  we have encountered this intersection point before.
C
      IF ((IW2.GT.IW1).AND.(JFULL.EQ.0)) THEN
      DO 130 I=IW1,IW2-1
         KORNR(I,1)=KAND1(IW1)
         KORNR(I,2)=KAND2(IW1)
         KORNR(I,3)=KAND1(IW2)
         KORNR(I,4)=KAND2(IW2)
130         CONTINUE
      ELSE
      IF (IW2.GT.IW1) THEN
         DO 140 I=IW1,IW2-1
      IF ((KORNR(I,1).EQ.KAND1(IW1)).AND.(KORNR(I,2).EQ.KAND2(IW1))
     +   .AND.(KORNR(I,3).EQ.KAND1(IW2)).AND.(KORNR(I,4).EQ.KAND2(IW2)))
     +    THEN
      GOTO 170
             ELSE
      m1=(y(kornr(i,2))-y(kornr(i,1)))/(x(kornr(i,2))-x(kornr(i,1)))
      m2=(y(kornr(i,4))-y(kornr(i,3)))/(x(kornr(i,4))-x(kornr(i,3)))
      if (m1.ne.m2) then
      xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-
     +        m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2)
      ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-
     +        m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2)
      endif
               if ((dabs(xcord1-xcord).le.eps).and.
     +             (dabs(ycord1-ycord).le.eps)) then
               goto 170
               endif

      KORNR(I,1)=KAND1(IW1)
      KORNR(I,2)=KAND2(IW1)
      KORNR(I,3)=KAND1(IW2)
      KORNR(I,4)=KAND2(IW2)
             ENDIF
140            CONTINUE
      ELSE
         JFULL=1
         DO 150 I=IW1,NUM
      KORNR(I,1)=KAND1(IW1)
      KORNR(I,2)=KAND2(IW1)
      KORNR(I,3)=KAND1(IW2)
      KORNR(I,4)=KAND2(IW2)
150            CONTINUE
         DO 160 I=1,IW2-1
       IF ((KORNR(I,1).EQ.KAND1(IW1)).AND.(KORNR(I,2).EQ.KAND2(IW1))
     +  .AND.(KORNR(I,3).EQ.KAND1(IW2)).AND.(KORNR(I,4).EQ.KAND2(IW2)))
     +       THEN
      GOTO 170
             ELSE
      m1=(y(kornr(i,2))-y(kornr(i,1)))/(x(kornr(i,2))-x(kornr(i,1)))
      m2=(y(kornr(i,4))-y(kornr(i,3)))/(x(kornr(i,4))-x(kornr(i,3)))
      if (m1.ne.m2) then
      xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-
     +        m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2)
      ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-
     +        m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2)
      endif
               if ((dabs(xcord1-xcord).le.eps).and.
     +             (dabs(ycord1-ycord).le.eps)) then
               goto 170
               endif

      KORNR(I,1)=KAND1(IW1)
      KORNR(I,2)=KAND2(IW1)
      KORNR(I,3)=KAND1(IW2)
      KORNR(I,4)=KAND2(IW2)
             ENDIF
160            CONTINUE
      ENDIF
         ENDIF
      ELSE
C
C  The intersection point is not the correct one,
C  try the next special k-divider.
C
      IW2=IW2+1
      IF (IW2.EQ.(NUM+1)) IW2=1
      GOTO 120
      ENDIF
C
C  Look for the next vertex of the convex figure.
C
      IW1=IW2
      IW2=IW2+1
      IF (IW2.EQ.(NUM+1)) IW2=1
      GOTO 120
200   CONTINUE
170   CONTINUE
      if (ndk.eq.0) then
          empty = 1
      else
          empty = 0
      endif
      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







      SUBROUTINE fillCont(X,Y,N,KORNR,MAXNUM,
     +                    XCont,YCont,KOUNT,NUM,EPS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is to be called by the iso2hdw algorithm.
C It finds vertices for the contour based on supplied candidates in
C the vector KORNR.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER N, MAXNUM, KOUNT, NUM
      DOUBLE PRECISION X(N), Y(N)
      INTEGER KORNR(MAXNUM,4)
      DOUBLE PRECISION XCORD,YCORD,E1,E2,F1,F2,G1,G2
      DOUBLE PRECISION XCORD1,YCORD1,XCORDP,YCORDP,BETA(N)
      DOUBLE PRECISION XCont(N*(N-1)/2), YCont(N*(N-1)/2)
      INTEGER I
      INTEGER HDEP1, JLV(N), JRV(N), F(N)
      DOUBLE PRECISION SDEP,DPF(N)
      DOUBLE PRECISION EPS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      KOUNT=0
      I=1
      E1=Y(KORNR(I,2))-Y(KORNR(I,1))
      F1=X(KORNR(I,1))-X(KORNR(I,2))
      G1=X(KORNR(I,1))*(Y(KORNR(I,2))-Y(KORNR(I,1)))
     +   -Y(KORNR(I,1))*(X(KORNR(I,2))-X(KORNR(I,1)))
      E2=Y(KORNR(I,4))-Y(KORNR(I,3))
      F2=X(KORNR(I,3))-X(KORNR(I,4))
      G2=X(KORNR(I,3))*(Y(KORNR(I,4))-Y(KORNR(I,3)))
     +   -Y(KORNR(I,3))*(X(KORNR(I,4))-X(KORNR(I,3)))
      XCORD=(-F2*G1+F1*G2)/(E2*F1-E1*F2)
      YCORD=(-E2*G1+E1*G2)/(E1*F2-E2*F1)
      XCont(1)=XCORD
      YCont(1)=YCORD
        xcord1=xcord
        ycord1=ycord
        xcordp=xcord
        ycordp=ycord
      KOUNT=KOUNT+1
      I=I+1
190   IF ((KORNR(I,1).EQ.KORNR(I-1,1)).AND.(KORNR(I,2).EQ.KORNR(I-1,2))
     +.AND.(KORNR(I,3).EQ.KORNR(I-1,3).AND.KORNR(I,4).EQ.KORNR(I-1,4)))
     +THEN
      I=I+1
      ELSE
        IF ((KORNR(I,1).EQ.KORNR(1,1)).AND.(KORNR(I,2).EQ.KORNR(1,2))
     +    .AND.(KORNR(I,3).EQ.KORNR(1,3).AND.KORNR(I,4).EQ.KORNR(1,4)))
     +  THEN
      GOTO 200
        ELSE
          E1=Y(KORNR(I,2))-Y(KORNR(I,1))
          F1=X(KORNR(I,1))-X(KORNR(I,2))
          G1=X(KORNR(I,1))*(Y(KORNR(I,2))-Y(KORNR(I,1)))
     +       -Y(KORNR(I,1))*(X(KORNR(I,2))-X(KORNR(I,1)))
          E2=Y(KORNR(I,4))-Y(KORNR(I,3))
          F2=X(KORNR(I,3))-X(KORNR(I,4))
          G2=X(KORNR(I,3))*(Y(KORNR(I,4))-Y(KORNR(I,3)))
     +       -Y(KORNR(I,3))*(X(KORNR(I,4))-X(KORNR(I,3)))
          XCORD=(-F2*G1+F1*G2)/(E2*F1-E1*F2)
          YCORD=(-E2*G1+E1*G2)/(E1*F2-E2*F1)
          if (((dabs(xcord-xcordp).lt.eps).and.
     +        (dabs(ycord-ycordp).lt.eps)).or.
     +        ((dabs(xcord-xcord1).lt.eps).and.
     +        (dabs(ycord-ycord1).lt.eps))) then
             i=i+1
          else
             xcordp=xcord
             ycordp=ycord
             KOUNT=KOUNT+1
         XCont(KOUNT)=XCORD
             YCont(KOUNT)=YCORD
         I=I+1
          endif
        call HSDEP21(xcord,ycord,n,x,y,
     +               beta,f,dpf,jlv,jrv,hdep1,SDEP)
        ENDIF
      ENDIF
      IF (I.NE.(NUM+1)) GOTO 190

200   CONTINUE
      RETURN
      END