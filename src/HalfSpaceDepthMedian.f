      SUBROUTINE HALFMED2D(XX,YY,N,tukmed,hdepth,dithered)
C
C  Version: 18 December 1998
C
C
C  If you want to run the program for a sample size which is larger than
C  the number which is now set for MAXN in the following PARAMETER
C  instruction, you can reset MAXN to the desired sample size.
C  MAXNUM should then be reset to the integer part of int(4*N*sqrt(N)+1).
C
      INTEGER N
      INTEGER NCIRQ(N),MCIRQ(N)
      INTEGER JLV(N*(N-1)/2),JRV(N*(N-1)/2),F(N)
      INTEGER IND1(N*(N-1)/2),IND2(N*(N-1)/2)
      INTEGER KORNR(int(4*N*SQRT(REAL(N))+1),4)
      INTEGER KOUNT,NUM
      INTEGER I,empty,ib,ie,le,M,MAXNUM,nceil
      DOUBLE PRECISION X(N),Y(N),WX(N),WY(N)
      DOUBLE PRECISION XX(N),YY(N)
      DOUBLE PRECISION xmean,ymean,xdev,ydev
      DOUBLE PRECISION X2(N),Y2(N)
      DOUBLE PRECISION ANGLE(N*(N-1)/2)
      DOUBLE PRECISION PI,PI2,EPS,xsum,ysum
      DOUBLE PRECISION XCont(N*(N-1)/2), YCont(N*(N-1)/2)
      DOUBLE PRECISION sum,tukmed(2),fac,wx1(N),wy1(N),rand(2)
      INTEGER HDEP,KSTAR
      double precision dpf(N),SDEP
      INTEGER HDEP1
      DOUBLE PRECISION HDEPTH,SDEP1, BETA(N)
      INTEGER dithered
      INTEGER seed

      M=N*(N-1)/2 
      MAXNUM=int(4*N*SQRT(REAL(N))+1)
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
      EPS=1.D-8
      fac=100000.0
      seed=256

      do 11 i=1,N
         ncirq(i)=0
         mcirq(i)=0
         f(i)=0
         BETA(I)=0.d0
         DPF(I)=0.d0                 
11    continue         
      do 12 i=1,M
         jlv(i)=0
         jrv(i)=0
         ind1(i)=0
         ind2(i)=0 
         ANGLE(I)=0 
         XCont(I)=0       
         YCont(I)=0
12    continue
      do 13 i=1,maxnum
         kornr(i,1)=0
         kornr(i,2)=0
         kornr(i,3)=0
         kornr(i,4)=0
13    continue 
      do 14 i=1,2
        tukmed(i)=0
        rand(i)=0
14    continue 
      kount=0
      num=0
      i=0
      empty=0
      ib=0
      ie=0
      le=0
      xsum=0
      ysum=0
      sum=0
      HDEP1=0
      HDEPTH=0.d0
      SDEP1=0.d0
      dithered=0
      
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
      DO 4 I=1,N
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
4     CONTINUE

C Initialise some vectors
      DO 5 I=1,N
        NCIRQ(I)=I
        MCIRQ(I)=I
 5    CONTINUE
        
      CALL checkData2D(X,Y,N,fac,NCIRQ,MCIRQ,ANGLE,
     +JLV,JRV,IND1,IND2,dithered)

      DO 6 I=1,N
        WX(I)=X(I)
        WY(I)=Y(I)         
 6    CONTINUE

C Calculation of the Tukey median.
      xsum=0
      ysum=0
      tukmed(1)=0
      tukmed(2)=0
C The case with n<=3 is trivial.
      if (n.le.3) then
         do 41 i=1,n
         xsum=xsum+x(i)
         ysum=ysum+y(i)
41       continue
         tukmed(1)=xsum/n
         tukmed(2)=ysum/n
         goto 610
      endif
      
      kstar=0
      do 38 i=1,n
         call HSDEP21(x(i),y(i),n,x,y,beta,f,dpf,jlv,jrv,hdep,SDEP)
         if (hdep.gt.kstar) kstar=hdep
38    continue

      ib=kstar-1
      ie=int(dble(n)/2)
171     le=ie-ib
        if (le.eq.0) goto 185
         CALL isoFin98(X,Y,N,ib+nceil(le,2),
     +   NCIRQ,MCIRQ,JLV,JRV,IND1,IND2,ANGLE,KORNR,EMPTY,NUM,EPS)
        if (empty.eq.1) ie=ib+nceil(le,2)
        if (empty.eq.0) ib=ib+nceil(le,2)
        if (le.eq.1) goto 185
        goto 171
       
185        CALL isoFin98(X,Y,N,ib,
     +   NCIRQ,MCIRQ,JLV,JRV,IND1,IND2,ANGLE,KORNR,EMPTY,NUM,EPS)

C  Scan KORNR and store the coordinates of the vertices.
      KOUNT=0
      CALL fillCont(X,Y,N,KORNR,MAXNUM,XCont,YCont,KOUNT,NUM,EPS)
       xsum=0
       ysum=0
       Do 204 I=1,kount
          wx(i)=XCont(i)
          wy(i)=YCont(i)
          xsum=xsum+XCont(i)
          ysum=ysum+YCont(i)
204    continue

C     Calculation of the center of gravity
         if (kount.gt.1) then
         do 205 i=1,kount
         wx(i)=wx(i)-(xsum/kount)
         wy(i)=wy(i)-(ysum/kount)
205      continue
         sum=0
         tukmed(1)=0
         tukmed(2)=0
         do 206 i=1,kount-1
         sum=sum+dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i))
         tukmed(1)=tukmed(1)+
     +  ((wx(i)+wx(i+1))*dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i)))   
         tukmed(2)=tukmed(2)+
     +  ((wy(i)+wy(i+1))*dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i)))   
206      continue
         sum=sum+dabs(wx(kount)*wy(1)-wx(1)*wy(kount))
         tukmed(1)=tukmed(1)+
     +  ((wx(kount)+wx(1))*dabs(wx(kount)*wy(1)-wx(1)*wy(kount)))   
         tukmed(2)=tukmed(2)+
     +  ((wy(kount)+wy(1))*dabs(wx(kount)*wy(1)-wx(1)*wy(kount)))   
         tukmed(1)=(tukmed(1)/(3*sum))+(xsum/kount)
         tukmed(2)=(tukmed(2)/(3*sum))+(ysum/kount)
         else
         tukmed(1)=xsum
         tukmed(2)=ysum
         endif
         call HSDEP21(tukmed(1),tukmed(2),N,X,Y,
     +                BETA,F,DPF,JLV,JRV,HDEP1,SDEP1)
         HDEPTH=(HDEP1+0.d0)/(N+0.d0)
         
         if (xdev.gt.eps) THEN
            tukmed(1)=tukmed(1)*xdev+xmean
         endif 
         if (ydev.gt.eps) THEN
            tukmed(2)=tukmed(2)*ydev+ymean
         endif   

      RETURN
610   END
      
