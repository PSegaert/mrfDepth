      SUBROUTINE BAGPLOTF(n,x,y,whisk,nsub,tukm,NChull,Chull,interpol,
     +                   num,datatyp,IndexOutliers,num3,datatyp2,
     +                   pxpy,boxpl,nointer,PDepths)
C
C  This subroutine computes a bagplot for any bivariate data set.
C
C  There is a choice between computing all the whiskers (take whisk equal
C  to 1), or few whiskers (only one whisker on each edge of the bag - take
C  whisk equal to 2), or no whiskers (take whisk equal to 3), or star-shaped
C  whiskers (take whisk equal to 4).
C  If the data set is not in general position, dithering is used.
C  Some adaptations are made for small n, and for large n (take subsets).
C  If all the points are (nearly) collinear, the bagplot reduces to the
C  usual boxplot for univariate data.
C  In this case, boxplot is set equal to 1 or 2.
C
C  Version: 15 December 1999
C
C
      implicit integer(i-n), double precision(a-h,o-z)
      INTEGER N,NSUB
      INTEGER MAXN,MAXNUM,MAXM,MAXNUM2,MAXM2
      INTEGER NCIRQ(N),MCIRQ(N),NRANK(N),F(N)
      INTEGER JLV(N*(N-1)/2),JRV(N*(N-1)/2),kstar
      INTEGER IND1(NSUB*(NSUB-1)/2),IND2(NSUB*(NSUB-1)/2),ind(N),jnd(N)
      INTEGER KORNR(int(4*NSUB*SQRT(REAL(NSUB))+1),4)
      INTEGER LEFT,KOUNT,kount1,NUM,num1,num2,num3,hdep
      INTEGER I,J,K,L,M,empty,ib,ie,le,tm,nc,jk,tel,tel2
      integer numdep(N),ja,jb,index(N),typ(N)
      integer indhulp,start,ii,istarx,nointer,indoutl(n)
      integer IndexOutliers(n)
      integer notingp,seed,starttel,a(N),ntot
      integer boxpl,whisk,dattm(2*N),numdattm,angi
      integer numk,numk1,NBP_NCEIL
      double precision X(N),Y(N),WX(NSUB*(NSUB-1))
      double precision WY(NSUB*(NSUB-1)),dpf(N)
      double precision ANGLE(NSUB*(NSUB-1)/2)
      double precision beta(N),wx1(N),wy1(N)
      double precision angz(2*N),angy1(N),angy2(N)
      double precision px(NSUB*(NSUB-1),4),py(NSUB*(NSUB-1),4)
      double precision gamma(NSUB*(NSUB-1))
      double precision lambda(2*N),lambdanc
      double precision wxjb1,wyjb1,wxja1,wyja1,hulp
      double precision PI,PI2,EPS,xcord1,ycord1,xsum,ysum
      double precision xcordp,ycordp,fac,xdev,ydev
      double precision XCORD,YCORD,E1,E2,F1,F2,G1,G2,ANG1
      double precision sum,tukmed(2),tukm(2),dist,c
      double precision star(N*3,2),ang,xmean,ymean
      double precision zori(N,2),datatyp(n,3)
      double precision interpol(n*2,2),pxpy(n,3)
      double precision datatyp2(n,2),xrand(N),yrand(2)
      double precision Chull(NSUB*(NSUB-1)/2,2)
      integer NChull
      double precision sdep
      double precision PDepths(N,2)
      INTEGER dummyVector(N)

      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
      EPS=1.D-8
      fac=100000.0
      SEED=256

      MAXN=N
      MAXNUM=int(4*N*SQRT(REAL(N))+1)
      MAXNUM2=int(4*NSUB*SQRT(REAL(NSUB))+1)
      MAXM=N*(N-1)/2
      MAXM2=NSUB*(NSUB-1)/2

      do 11 i=1,maxn
         ncirq(i)=0
         mcirq(i)=0
         nrank(i)=0
         f(i)=0
         ind(i)=0
         jnd(i)=0
         numdep(i)=0
         index(i)=0
         typ(i)=0
         a(i)=0
         dattm(i)=0
11    continue
      do 12 i=1,maxm
         jlv(i)=0
         jrv(i)=0
12    continue
      do 13 i=1,maxnum2
         kornr(i,1)=0
         kornr(i,2)=0
         kornr(i,3)=0
         kornr(i,4)=0
13    continue
      left=0
      kount=0
      kount1=0
      num=0
      num1=0
      num2=0
      num3=0
      hdep=0
      j=0
      l=0
      empty=0
      ib=0
      ie=0
      le=0
      tm=0
      nc=0
      jk=0
      tel=0
      tel2=0
      ja=0
      jb=0
      indhulp=0
      start=0
      ii=0
      istarx=0
      notingp=0
      starttel=0
      ntot=0
      boxpl=0
      numdattm=0
      angi=0
      nointer=0
      do 14 i=1,maxn
         wx(i)=0.0
         wy(i)=0.0
         dpf(i)=0.0
         beta(i)=0.0
         angz(i)=0.0
         angy1(i)=0.0
         angy2(i)=0.0
         px(i,1)=0.0
         px(i,2)=0.0
         px(i,3)=0.0
         px(i,4)=0.0
         py(i,1)=0.0
         py(i,2)=0.0
         py(i,3)=0.0
         py(i,4)=0.0
         wx1(i)=0.0
         wy1(i)=0.0
         lambda(i)=0.0
         gamma(i)=0.0
         zori(i,1)=0.0
         zori(i,2)=0.0
         xrand(i)=0.0
         PDepths(i,1)=0.0
         PDepths(i,2)=0.0
         IndexOutliers(i)=0
14    continue
      NChull = 0
      do 16 i=1,maxm2
         ind1(i)=0
         ind2(i)=0
         angle(i)=0.0
         Chull(i,1)=0.0
         Chull(i,2)=0.0
16    continue
      do 18 i=1,n*3
         star(i,1)=0.0
         star(i,2)=0.0
18    continue
      kstar=0
      lambdanc=0.0
      wxjb1=0.0
      wyjb1=0.0
      wxja1=0.0
      wyja1=0.0
      hulp=0.0
      xcord=0.0
      ycord=0.0
      xcord1=0.0
      ycord1=0.0
      xcordp=0.0
      ycordp=0.0
      xsum=0.0
      ysum=0.0
      e1=0.0
      e2=0.0
      f1=0.0
      f2=0.0
      g1=0.0
      g2=0.0
      ang1=0.0
      ang=0.0
      sum=0.0
      tukmed(1)=0.0
      tukmed(2)=0.0
      tukm(1)=0.0
      tukm(2)=0.0
      dist=0.0
      c=0.0
      yrand(1)=0.0
      yrand(2)=0.0
      if (n.lt.20) whisk=1
C
C     Standardization of the data set.
C
      xmean=0.0
      ymean=0.0
      xdev=0.0
      ydev=0.0
      DO 2 I=1,N
         xmean=xmean+x(i)
         ymean=ymean+y(i)
2     CONTINUE
      xmean=xmean/n
      ymean=ymean/n
      DO 3 I=1,n
         xdev=xdev+((x(i)-xmean)*(x(i)-xmean))
         ydev=ydev+((y(i)-ymean)*(y(i)-ymean))
3     continue
      xdev=xdev/(n-1)
      ydev=ydev/(n-1)
      xdev=dsqrt(xdev)
      ydev=dsqrt(ydev)
      DO 5 I=1,N
         if (xdev.gt.eps) x(i)=(x(i)-xmean)/xdev
         if (ydev.gt.eps) y(i)=(y(i)-ymean)/ydev
         zori(i,1)=x(i)
         zori(i,2)=y(i)
         WX(I)=X(I)
         WY(I)=Y(I)
         wx1(i)=x(i)
         wy1(i)=y(i)
         NCIRQ(I)=I
         MCIRQ(I)=I
 5    CONTINUE
C
C     If n is large, take a subset.
C

      if (n.gt.nsub) then
         ntot=n
         n=nsub
         call rdraw(a,ntot,seed,n)
         do 41 i=1,n
            x(i)=wx1(a(i))
            y(i)=wy1(a(i))
            wx(i)=x(i)
            wy(i)=y(i)
41       continue
      else
          nsub=n
          MAXNUM2=int(4*NSUB*SQRT(REAL(NSUB))+1)
      endif

      M=((N*(N-1))/2)

C
C     Test whether more than half of the points lie on a vertical line.
C
      boxpl=0
      numdep(1)=1
      tel=1
      i=2
      j=1
42    if (dabs(x(i)-x(j)).lt.eps) then
         numdep(j)=numdep(j)+1
         if (numdep(j).gt.idnint(dble(n/2))) then
            boxpl=1
            if (xdev.gt.eps) then
            interpol(1,1)=x(j)*xdev+xmean
            interpol(2,1)=x(i)*xdev+xmean
            else
            interpol(1,1)=x(j)
            interpol(2,1)=x(i)
            endif
            if (ydev.gt.eps) then
            interpol(1,2)=y(j)*ydev+ymean
            interpol(2,2)=y(i)*ydev+ymean
            else
            interpol(1,2)=y(j)
            interpol(2,2)=y(i)
            endif
            goto 610
         endif
         i=i+1
         if(i.gt.n) then
            goto 46
         endif
         j=1
         goto 42
      else
         j=j+1
         if (j.eq.i) then
            numdep(i)=1
            i=i+1
            j=1
            tel=tel+1
            if (tel.gt.idnint(dble(n/2))+1) goto 46
         endif
         goto 42
      endif
46    continue

      CALL checkData2D(X,Y,N,fac,NCIRQ,MCIRQ,ANGLE,
     +                 JLV,JRV,IND1,IND2,notingp)

      DO 47 I=1,N
         WX(I)=X(I)
         WY(I)=Y(I)
47    CONTINUE

C
C     If all points are collinear, use a univariate boxplot.
C
      boxpl=0
      angi=1
      if (angle(m)-angle(1).le.eps) then
         boxpl=1
         goto 27
      endif
      if ((angle(1).le.eps).and.(angle(m).ge.pi*2-eps)) then
         do 26 i=2,m-1
            if ((angle(i).le.eps).or.(angle(i).ge.pi*2-eps)) then
               boxpl=1
            else
               boxpl=0
               goto 27
            endif
26       continue
      endif
C
C     Count whether 50% of the data points are collinear
C
      boxpl=0
      ang=angle(1)
      angi=1
      tel=1
      do 126 i=2,m
         if (dabs(angle(i)-ang).le.eps) then
            tel=tel+1
            if (tel.gt.(n*(n-2)/8)) then
               boxpl=1
               goto 27
            endif
         else
            ang=angle(i)
            angi=i
            tel=1
         endif
126   continue

27    if (boxpl.eq.1) then
C     All points collinear
         if (xdev.gt.eps) then
            interpol(1,1)=x(ind1(angi))*xdev+xmean
            interpol(2,1)=x(ind2(angi))*xdev+xmean
         else
            interpol(1,1)=x(ind1(angi))
            interpol(2,1)=x(ind2(angi))
         endif
         if (ydev.gt.eps) then
            interpol(1,2)=y(ind1(angi))*ydev+ymean
            interpol(2,2)=y(ind2(angi))*ydev+ymean
         else
            interpol(1,2)=y(ind1(angi))
            interpol(2,2)=y(ind2(angi))
         endif
         goto 610
      endif

C
C      The data are in general position.
C
37    if (notingp.eq.1) then
      do 212 i=1,n
          numdep(i)=0
212   continue
      endif
      kstar=0
      do 138 i=1,n
         numdep(i)=0
138   continue
      do 38 i=1,n
         call HSDEP21(x(i),y(i),n,x,y,beta,f,dpf,jlv,jrv,hdep,SDEP)
         PDepths(i,1)=(hdep+0.d0)/(N+0.d0)
         PDepths(i,2)=SDEP
         if (hdep.gt.kstar) kstar=hdep
         numdep(hdep)=numdep(hdep)+1
38    continue
C
C     Calculation of the Tukey median.
C
      tm=0
      xsum=0
      ysum=0
      tukmed(1)=0
      tukmed(2)=0
      if (n.le.3) then
         do 40 i=1,n
         xsum=xsum+x(i)
         ysum=ysum+y(i)
      if (xdev.gt.eps) datatyp(i,1)=x(i)*xdev+xmean
      if (ydev.gt.eps) datatyp(i,2)=y(i)*ydev+ymean
40       continue
         tukm(1)=xsum/n
      if (xdev.gt.eps) tukm(1)=tukm(1)*xdev+xmean
         tukm(2)=ysum/n
      if (ydev.gt.eps) tukm(2)=tukm(2)*ydev+ymean
         goto 800
      endif
      empty=0
      ib=kstar
      ie=int(dble(n)/2)
180     le=ie-ib
        if (le.lt.0) then
           nointer=1
           le=0
        endif
        if (le.eq.0) goto 185
       CALL isoFin98(X,Y,N,ib+nbp_nceil(le,2),NCIRQ,MCIRQ,JLV,JRV,
     +               IND1,IND2,ANGLE,KORNR(1:MAXNUM2,:),EMPTY,NUM,EPS)
        if (empty.eq.1) ie=ib+nbp_nceil(le,2)
        if (empty.eq.0) ib=ib+nbp_nceil(le,2)
        if (le.eq.1) goto 185
        goto 180
185     CALL isoFin98(X,Y,N,ib,NCIRQ,MCIRQ,JLV,JRV,
     +               IND1,IND2,ANGLE,KORNR(1:MAXNUM2,:),EMPTY,NUM,EPS)

C
C     Scan KORNR and compute coordinates of the vertices.
C
186   KOUNT=0
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
         wx(kount+1)=xcord
         wy(kount+1)=ycord
      if (tm.eq.0) then
         xsum=xcord
         ysum=ycord
      endif
      xcord1=xcord
      ycord1=ycord
      xcordp=xcord
      ycordp=ycord
      KOUNT=KOUNT+1
      I=I+1
      if (num.eq.1) goto 195
190   IF ((KORNR(I,1).EQ.KORNR(I-1,1)).AND.(KORNR(I,2).EQ.KORNR(I-1,2))
     +.AND.(KORNR(I,3).EQ.KORNR(I-1,3).AND.KORNR(I,4).EQ.KORNR(I-1,4)))
     +THEN
      I=I+1
      ELSE
        IF ((KORNR(I,1).EQ.KORNR(1,1)).AND.(KORNR(I,2).EQ.KORNR(1,2))
     +    .AND.(KORNR(I,3).EQ.KORNR(1,3).AND.KORNR(I,4).EQ.KORNR(1,4)))
     +  THEN
      GOTO 195
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
     +       ((dabs(xcord-xcord1).lt.eps).and.
     +        (dabs(ycord-ycord1).lt.eps))) then
             i=i+1
          else
             xcordp=xcord
             ycordp=ycord
             wx(kount+1)=xcord
             wy(kount+1)=ycord
      if (tm.eq.0) then
             xsum=xsum+xcord
             ysum=ysum+ycord
      endif
      KOUNT=KOUNT+1
      I=I+1
          endif
        ENDIF
      ENDIF
      IF (I.NE.(NUM+1)) GOTO 190
195   if (tm.eq.2) goto 500
      if (tm.eq.1) goto 300
c
c     Calculation of the center of gravity.
c
      if (tm.eq.0) then
         NChull = kount
          if (kount.gt.1) then
         do 205 i=1,kount
         if (xdev.gt.eps) then
            Chull(i,1)=wx(i)*xdev+xmean
         else
            Chull(i,1)=wx(i)
         endif
         if (ydev.gt.eps) then
            Chull(i,2)=wy(i)*ydev+ymean
         else
            Chull(i,2)=wy(i)
         endif
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
         Chull(i,2)=xsum
         tukmed(2)=ysum
         Chull(i,2)=ysum
         endif
207      if (xdev.gt.eps) then
            xcord=tukmed(1)*xdev+xmean
         else
            xcord=tukmed(1)
         endif
         if (ydev.gt.eps) then
            ycord=tukmed(2)*ydev+ymean
         else
            ycord=tukmed(2)
         endif
            tukm(1)=xcord
            tukm(2)=ycord
      endif
      tm=1
      if (nointer.eq.1) goto 800
C
C     Calculation of correct value of k.
C
      nc=int(n/2)
      j=kstar+1
200   j=j-1
      if (numdep(kstar).le.nc) then
         numdep(kstar)=numdep(kstar)+numdep(j-1)
         goto 200
      endif
      k=j+1
      numk1=numdep(kstar)
      numk=numk1-numdep(k-1)
      lambdanc=dble(nc-numk)/dble(numk1-numk)
C
C     Calculation of the vertices of Dk.
C
       CALL isoFin98(X,Y,N,k,NCIRQ,MCIRQ,JLV,JRV,
     +               IND1,IND2,ANGLE,KORNR(1:MAXNUM2,:),EMPTY,NUM,EPS)

      if (empty.eq.0) goto 186
300   continue
      kount1=kount
      do 201 i=1,kount1
      wx1(i)=wx(i)
      wy1(i)=wy(i)
201   continue
C
C     Calculation of the vertices of Dk-1.
C
      tm=2
       CALL isoFin98(X,Y,N,k-1,NCIRQ,MCIRQ,JLV,JRV,
     +               IND1,IND2,ANGLE,KORNR(1:MAXNUM2,:),EMPTY,NUM,EPS)
      if (empty.eq.0) goto 186
500   continue
C
C     If Dk-1 is a line segment, use a univariate boxplot.
C
      if (dabs(wx(2)-wx(1)).gt.eps) then
         hulp=(wy(2)-wy(1))/(wx(2)-wx(1))
         do 506 i=3,kount
         if (dabs(wy(i)-wy(1)-hulp*(wx(i)-wx(1))).gt.eps*10) goto 505
506      continue
         boxpl=2
         goto 509
      else
         do 504 i=3,kount
         if (dabs(wx(i)-wx(1)).gt.eps*10) goto 505
504      continue
         boxpl=2
         goto 509
      endif
509   continue
      do 507 i=1,n
      if (xdev.gt.eps) datatyp(i,1)=x(i)*xdev+xmean
      if (ydev.gt.eps) datatyp(i,2)=y(i)*ydev+ymean
      datatyp(i,3)=4.0
507   continue
      do 508 i=1,kount
      if (xdev.gt.eps) interpol(i,1)=wx(i)*xdev+xmean
      if (ydev.gt.eps) interpol(i,2)=wy(i)*ydev+ymean
508   continue
      goto 610

505      if (n.ge.10) then
      do 202 i=1,kount1
      wx1(i)=wx1(i)-tukmed(1)
      wy1(i)=wy1(i)-tukmed(2)
202   continue
      do 203 i=1,kount
      wx(i)=wx(i)-tukmed(1)
      wy(i)=wy(i)-tukmed(2)
203   continue
      do 204 i=1,n
      x(i)=x(i)-tukmed(1)
      y(i)=y(i)-tukmed(2)
204   continue
      do 1204 i=1,ntot
      zori(i,1)=zori(i,1)-tukmed(1)
      zori(i,2)=zori(i,2)-tukmed(2)
1204   continue
C
C     Compute angles of the data points.
C
      numdattm=0
      do 400 i=1,n
      index(i)=i
      if ((dabs(x(i)).lt.eps).and.(dabs(y(i)).lt.eps)) then
         numdattm=numdattm+1
         dattm(numdattm)=i
         angz(i)=1000.0
      else
      dist=dsqrt(x(i)*x(i)+y(i)*y(i))
      xcord=x(i)/dist
      ycord=y(i)/dist
      if (dabs(xcord).gt.dabs(ycord)) then
         if (xcord.ge.0.0) then
            angz(i)=dasin(ycord)
            if (angz(i).lt.0.0) angz(i)=angz(i)+pi*2
         else
            angz(i)=pi-dasin(ycord)
         endif
      else
         if (ycord.ge.0.0) then
            angz(i)=dacos(xcord)
         else
            angz(i)=pi*2-dacos(xcord)
         endif
      endif
         if (angz(i).ge.(pi*2-eps)) angz(i)=0.0
      endif
400   continue
C     numdattm is the number of data points equal to the Tukey median.
      call SORT(angz,index,dummyVector,dpf,n,jlv,jrv)
      do 401 i=1,n
         indoutl(i)=index(i)
401   continue
      n=n-numdattm
C
C     Compute angles of the vertices of the depth region Dk.
C
      if (kount1.eq.1) goto 412
      do 410 i=1,kount1
      ind(i)=i
      dist=dsqrt(wx1(i)*wx1(i)+wy1(i)*wy1(i))
      if (dist.gt.eps) then
      xcord=wx1(i)/dist
      ycord=wy1(i)/dist
      if (dabs(xcord).gt.dabs(ycord)) then
         if (xcord.ge.0.0) then
            angy1(i)=dasin(ycord)
            if (angy1(i).lt.0.0) angy1(i)=angy1(i)+pi*2
         else
            angy1(i)=pi-dasin(ycord)
         endif
      else
         if (ycord.ge.0.0) then
            angy1(i)=dacos(xcord)
         else
            angy1(i)=pi*2-dacos(xcord)
         endif
      endif
         if (angy1(i).ge.(pi*2-eps)) angy1(i)=0.0
      else
      nointer=1
      endif
410   continue
      call SORT(angy1,ind,dummyVector,dpf,kount1,jlv,jrv)
C
C     Compute angles of the vertices of the depth region Dk-1.
C
412      do 510 i=1,kount
      jnd(i)=i
      dist=dsqrt(wx(i)*wx(i)+wy(i)*wy(i))
      if (dist.gt.eps) then
      xcord=wx(i)/dist
      ycord=wy(i)/dist
      if (dabs(xcord).gt.dabs(ycord)) then
         if (xcord.ge.0.0) then
            angy2(i)=dasin(ycord)
            if (angy2(i).lt.0.0) angy2(i)=angy2(i)+pi*2
         else
            angy2(i)=pi-dasin(ycord)
         endif
      else
         if (ycord.ge.0.0) then
            angy2(i)=dacos(xcord)
         else
            angy2(i)=pi*2-dacos(xcord)
         endif
      endif
         if (angy2(i).ge.(pi*2-eps)) angy2(i)=0.0
      else
      nointer=1
      endif
510   continue
      call SORT(angy2,jnd,dummyVector,dpf,kount,jlv,jrv)
C
C     Calculation of arrays px and py for Dk.
C
      if (kount1.eq.1) goto 431
      jk=0
      wx1(kount1+1)=wx1(1)
      wy1(kount1+1)=wy1(1)
      angy1(kount1+1)=angy1(1)
      ind(kount1+1)=ind(1)
      if (angz(1).lt.angy1(1)) j=kount1
      if (angz(1).ge.angy1(1)-eps) j=1
      do 429 i=1,n
420   if ((angz(i).ge.angy1(j+1)-eps).and.(jk.eq.0)) then
         j=j+1
         if (j.eq.kount1) jk=1
         if (j.eq.kount1+1) j=1
         goto 420
      endif
      if ((dabs(wx1(ind(j+1))-wx1(ind(j))).gt.eps).and.
     +    (dabs(x(index(i))).gt.eps)) then
      dist=(y(index(i))/x(index(i)))-
     + ((wy1(ind(j+1))-wy1(ind(j)))/(wx1(ind(j+1))-wx1(ind(j))))
      px(i,1)=(wy1(ind(j))*wx1(ind(j+1))-
     + wx1(ind(j))*wy1(ind(j+1)))/(wx1(ind(j+1))-wx1(ind(j)))
      px(i,1)=px(i,1)/dist
      py(i,1)=px(i,1)*y(index(i))/x(index(i))
      else
         if (dabs(wx1(ind(j+1))-wx1(ind(j))).le.eps) then
      px(i,1)=wx1(ind(j))
      py(i,1)=px(i,1)*y(index(i))/x(index(i))
         else
      px(i,1)=0.0
      py(i,1)=((wy1(ind(j))-wy1(ind(j+1)))*wx1(ind(j))/
     +        (wx1(ind(j+1))-wx1(ind(j))))+wy1(ind(j))
         endif
      endif
429   continue
C
C     Calculation of arrays px and py for Dk-1.
C
431      jk=0
      wx(kount+1)=wx(1)
      wy(kount+1)=wy(1)
      angy2(kount+1)=angy2(1)
      jnd(kount+1)=jnd(1)
      if (angz(1).lt.angy2(1)) j=kount
      if (angz(1).ge.angy2(1)-eps) j=1
      do 529 i=1,n
520   if ((angz(i).ge.angy2(j+1)-eps).and.(jk.eq.0)) then
         j=j+1
         if (j.eq.kount) jk=1
         if (j.eq.kount+1) j=1
         goto 520
      endif
      if ((dabs(wx(jnd(j+1))-wx(jnd(j))).gt.eps).and.
     +    (dabs(x(index(i))).gt.eps)) then
      dist=(y(index(i))/x(index(i)))-
     + ((wy(jnd(j+1))-wy(jnd(j)))/(wx(jnd(j+1))-wx(jnd(j))))
      px(i,2)=(wy(jnd(j))*wx(jnd(j+1))-
     + wx(jnd(j))*wy(jnd(j+1)))/(wx(jnd(j+1))-wx(jnd(j)))
      px(i,2)=px(i,2)/dist
      py(i,2)=px(i,2)*y(index(i))/x(index(i))
      else
         if (dabs(wx(jnd(j+1))-wx(jnd(j))).le.eps) then
      px(i,2)=wx(jnd(j))
      py(i,2)=px(i,2)*y(index(i))/x(index(i))
         else
      px(i,2)=0.0
      py(i,2)=((wy(jnd(j))-wy(jnd(j+1)))*wx(jnd(j))/
     +  (wx(jnd(j+1))-wx(jnd(j))))+wy(jnd(j))
         endif
      endif
529   continue
      do 540 i=1,n
      if (kount1.eq.1) then
         px(i,1)=wx1(1)
         py(i,1)=wy1(1)
      endif
540   continue
C
C     Mergesort of angy1 and angy2 to obtain gamma
C     and calculation of arrays px and py.
C
      if (kount1.eq.1) then
         do 599 i=1,kount
            gamma(i)=angy2(i)
            px(i,3)=wx1(i)
            py(i,3)=wy1(i)
            px(i,4)=wx(i)
            py(i,4)=wy(i)
599      continue
         goto 602
      endif
      ja=1
      jb=1
      i=1
600   if (angy1(ja).le.angy2(jb)) then
         if (ja.le.kount1) then
            gamma(i)=angy1(ja)
            px(i,3)=wx1(ind(ja))
            py(i,3)=wy1(ind(ja))
            if (jb.eq.1) then
               wxjb1=wx(jnd(kount))
               wyjb1=wy(jnd(kount))
            else
               wxjb1=wx(jnd(jb-1))
               wyjb1=wy(jnd(jb-1))
            endif
       if ((dabs(wx(jnd(jb))-wxjb1).gt.eps).and.
     +     (dabs(wx1(ind(ja))).gt.eps)) then
       dist=(wy1(ind(ja))/wx1(ind(ja)))-
     + ((wy(jnd(jb))-wyjb1)/(wx(jnd(jb))-wxjb1))
            px(i,4)=(wyjb1*wx(jnd(jb))-
     + wxjb1*wy(jnd(jb)))/(wx(jnd(jb))-wxjb1)
            px(i,4)=px(i,4)/dist
            py(i,4)=px(i,4)*wy1(ind(ja))/wx1(ind(ja))
        else
           if (dabs(wx(jnd(jb))-wxjb1).le.eps) then
               px(i,4)=wxjb1
               py(i,4)=px(i,4)*wy1(ind(ja))/wx1(ind(ja))
           else
               px(i,4)=0.0
               py(i,4)=((wyjb1-wy(jnd(jb)))*wxjb1/
     +       (wx(jnd(jb))-wxjb1))+wyjb1
           endif
        endif
            ja=ja+1
            i=i+1
         else
            angy1(ja)=angy1(ja)+10
         endif
      else
         if (jb.le.kount) then
            gamma(i)=angy2(jb)
            px(i,4)=wx(jnd(jb))
            py(i,4)=wy(jnd(jb))
            if (ja.eq.1) then
               wxja1=wx1(ind(kount1))
               wyja1=wy1(ind(kount1))
            else
               wxja1=wx1(ind(ja-1))
               wyja1=wy1(ind(ja-1))
            endif
       if ((dabs(wx1(ind(ja))-wxja1).gt.eps).and.
     +     (dabs(wx(jnd(jb))).gt.eps)) then
       dist=(wy(jnd(jb))/wx(jnd(jb)))-
     + ((wy1(ind(ja))-wyja1)/(wx1(ind(ja))-wxja1))
            px(i,3)=(wyja1*wx1(ind(ja))-
     + wxja1*wy1(ind(ja)))/(wx1(ind(ja))-wxja1)
            px(i,3)=px(i,3)/dist
            py(i,3)=px(i,3)*wy(jnd(jb))/wx(jnd(jb))
       else
          if (dabs(wx1(ind(ja))-wxja1).le.eps) then
             px(i,3)=wxja1
             py(i,3)=px(i,3)*wy(jnd(jb))/wx(jnd(jb))
          else
             px(i,3)=0.0
             py(i,3)=((wyja1-wy1(ind(ja)))*wxja1/
     +    (wx1(ind(ja))-wxja1))+wyja1
          endif
      endif
            jb=jb+1
            i=i+1
         else
            angy2(jb)=angy2(jb)+10
         endif
      endif
      if (i.le.kount1+kount) goto 600
            do 601 i=1,kount1+kount
601   continue
C
C     Interpolation of two arrays px and py.
C
602      c=3.0
       if (nointer.eq.1) then
          kount1=0
       endif
       do 605 i=1,kount1+kount
         if (nointer.eq.0) then
       wx(i)=lambdanc*px(i,4)+(1-lambdanc)*px(i,3)
       wy(i)=lambdanc*py(i,4)+(1-lambdanc)*py(i,3)
         endif
       if (xdev.gt.eps) then
          xcord=(wx(i)+tukmed(1))*xdev+xmean
       else
          xcord=wx(i)+tukmed(1)
       endif
       if (ydev.gt.eps) then
          ycord=(wy(i)+tukmed(2))*ydev+ymean
       else
          ycord=wy(i)+tukmed(2)
       endif
       interpol(i,1)=xcord
       interpol(i,2)=ycord
605    continue
       if (nointer.eq.1) then
C      No interpolation, too many points coincide.
           do 1605 i=1,n
           if (xdev.gt.eps) then
              xcord=(x(index(i))+tukmed(1))*xdev+xmean
           else
              xcord=x(index(i))+tukmed(1)
           endif
           if (ydev.gt.eps) then
              ycord=(y(index(i))+tukmed(2))*ydev+ymean
           else
              ycord=y(index(i))+tukmed(2)
           endif
              datatyp(i,1)=xcord
              datatyp(i,2)=ycord
              datatyp(i,3)=1.0
1605       continue
           do 1606 i=1,numdattm
           if (xdev.gt.eps) then
              xcord=(x(dattm(i))+tukmed(1))*xdev+xmean
           else
              xcord=x(dattm(i))+tukmed(1)
           endif
           if (ydev.gt.eps) then
              ycord=(y(dattm(i))+tukmed(2))*ydev+ymean
           else
              ycord=y(dattm(i))+tukmed(2)
           endif
               datatyp(n+i,1)=xcord
               datatyp(n+i,2)=ycord
               datatyp(n+i,3)=0.0
1606      continue
          goto 610
       endif
C
C      Repeat some calculations for the whole data set.
C
1398      if (ntot.gt.nsub) then
         n=ntot
         do 1399 i=1,n
         x(i)=zori(i,1)
         y(i)=zori(i,2)
1399      continue
      else
         goto 1531
      endif
C
C     Compute angles of the data points.
C
      do 1300 i=1,n
         x(i)=x(i)-tukmed(1)
         y(i)=y(i)-tukmed(2)
1300  continue
      numdattm=0
      do 1400 i=1,n
      index(i)=i
      if ((dabs(x(i)).lt.eps).and.(dabs(y(i)).lt.eps)) then
         angz(i)=1000.0
         numdattm=numdattm+1
         dattm(numdattm)=i
      else
      dist=dsqrt(x(i)*x(i)+y(i)*y(i))
      xcord=x(i)/dist
      ycord=y(i)/dist
      if (dabs(xcord).gt.dabs(ycord)) then
         if (xcord.ge.0.0) then
            angz(i)=dasin(ycord)
            if (angz(i).lt.0.0) angz(i)=angz(i)+pi*2
         else
            angz(i)=pi-dasin(ycord)
         endif
      else
         if (ycord.ge.0.0) then
            angz(i)=dacos(xcord)
         else
            angz(i)=pi*2-dacos(xcord)
         endif
      endif
         if (angz(i).ge.(pi*2-eps)) angz(i)=0.0
      endif
1400   continue
      call SORT(angz,index,dummyVector,dpf,n,jlv,jrv)
      do 1401 i=1,n
         indoutl(i)=index(i)
1401   continue
      n=n-numdattm
C
C     Calculation of arrays px and py for B.
C
1431      jk=0
      wx(kount+kount1+1)=wx(1)
      wy(kount+kount1+1)=wy(1)
      gamma(kount+kount1+1)=gamma(1)
      if (angz(1).lt.gamma(1)) j=kount+kount1
      if (angz(1).ge.gamma(1)-eps) j=1
      do 1529 i=1,n
1520   if ((angz(i).ge.gamma(j+1)-eps).and.(jk.eq.0)) then
         j=j+1
         if (j.eq.kount+kount1) jk=1
         if (j.eq.kount+kount1+1) j=1
         goto 1520
      endif
      if ((dabs(wx(j+1)-wx(j)).gt.eps).and.
     +    (dabs(x(index(i))).gt.eps)) then
      dist=(y(index(i))/x(index(i)))-
     + ((wy(j+1)-wy(j))/(wx(j+1)-wx(j)))
      px(i,1)=(wy(j)*wx(j+1)-
     + wx(j)*wy(j+1))/(wx(j+1)-wx(j))
      px(i,1)=px(i,1)/dist
      py(i,1)=px(i,1)*y(index(i))/x(index(i))
      else
         if (dabs(wx(j+1)-wx(j)).le.eps) then
      px(i,1)=wx(j)
      py(i,1)=px(i,1)*y(index(i))/x(index(i))
         else
      px(i,1)=0.0
      py(i,1)=((wy(j)-wy(j+1))*wx(j)/
     +  (wx(j+1)-wx(j)))+wy(j)
         endif
      endif
1529   continue
C
C      Decide on the type of each data point.
C
1531   c=3.0
       num=0
       num1=0
       num2=0
       num3=0
       do 1639 i=1,n
       if (ntot.le.nsub) then
          px(i,1)=lambdanc*px(i,2)+(1-lambdanc)*px(i,1)
          py(i,1)=lambdanc*py(i,2)+(1-lambdanc)*py(i,1)
       endif
       if (px(i,1)*px(i,1)+py(i,1)*py(i,1).gt.eps) then
          lambda(i)=dsqrt(x(index(i))*x(index(i))+
     + y(index(i))*y(index(i)))/
     + dsqrt(px(i,1)*px(i,1)+py(i,1)*py(i,1))
       else
          if ((dabs(x(index(i))).lt.eps).and.
     +        (dabs(y(index(i))).lt.eps)) then
             lambda(i)=0.0
          else
             lambda(i)=c+1.0
          endif
       endif
       if (lambda(i).lt.eps) then
         num=num+1
         typ(i)=0
         goto 639
       endif
       if (lambda(i).le.1+eps) then
         num1=num1+1
         typ(i)=1
         goto 639
       endif
       if (lambda(i).le.c+eps) then
         num2=num2+1
         typ(i)=2
         goto 639
       endif
       if (lambda(i).gt.c) then
         num3=num3+1
         typ(i)=3
         indoutl(i)=index(i)
         IndexOutliers(num3) = index(i)
       endif
639    if (xdev.gt.eps) then
          xcord=(x(index(i))+tukmed(1))*xdev+xmean
       else
          xcord=x(index(i))+tukmed(1)
       endif
       if (ydev.gt.eps) then
          ycord=(y(index(i))+tukmed(2))*ydev+ymean
       else
          ycord=y(index(i))+tukmed(2)
       endif
           datatyp(i,1)=xcord
           datatyp(i,2)=ycord
           datatyp(i,3)=dble(typ(i))
1639       continue
       if (numdattm.ne.0) then
             num=num+numdattm
       do 1644 i=1,numdattm
             typ(n+i)=0
             lambda(n+i)=0.0
             px(n+i,1)=x(dattm(i))
             py(n+i,1)=y(dattm(i))
       if (xdev.gt.eps) then
          xcord=(x(dattm(i))+tukmed(1))*xdev+xmean
       else
          xcord=x(dattm(i))+tukmed(1)
       endif
       if (ydev.gt.eps) then
          ycord=(y(dattm(i))+tukmed(2))*ydev+ymean
       else
          ycord=y(dattm(i))+tukmed(2)
       endif
           datatyp(n+i,1)=xcord
           datatyp(n+i,2)=ycord
           datatyp(n+i,3)=0.0
1644   continue
       endif

      if (whisk.eq.1) goto 1640
      if (whisk.eq.2) goto 1641
1640  do 640 i=1,n
      if (xdev.gt.eps) then
         xcord=(px(i,1)+tukmed(1))*xdev+xmean
      else
         xcord=px(i,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(py(i,1)+tukmed(2))*ydev+ymean
      else
         ycord=py(i,1)+tukmed(2)
      endif
          pxpy(i,1)=xcord
          pxpy(i,2)=ycord
          pxpy(i,3)=0.0
640   continue
      if (whisk.eq.1) goto 610
      if (whisk.eq.3) goto 610
C
C     Retain only one whisker for each edge.
C
 1641 if (ntot.gt.1000) then
         n=500
         call rdraw(a,ntot,seed,n)
         do 2642 i=1,n
         if (xdev.gt.eps) then
            xcord=(x(index(a(i)))+tukmed(1))*xdev+xmean
         else
            xcord=x(index(a(i)))+tukmed(1)
         endif
         if (ydev.gt.eps) then
            ycord=(y(index(a(i)))+tukmed(2))*ydev+ymean
         else
            ycord=y(index(a(i)))+tukmed(2)
         endif
          datatyp2(i,1)=xcord
          datatyp2(i,2)=ycord
2642     continue
         do 1642 i=1,n
            beta(i)=angz(a(i))
            dattm(i)=typ(a(i))
            dpf(i)=lambda(a(i))
            px(i,2)=px(a(i),1)
            py(i,2)=py(a(i),1)
1642     continue
         do 1643 i=1,n
            angz(i)=beta(i)
            typ(i)=dattm(i)
            lambda(i)=dpf(i)
            px(i,1)=px(i,2)
            py(i,1)=py(i,2)
1643     continue
      endif
      if (num2.ne.0) then
      tel=1
      tel2=1
      gamma(kount1+kount+1)=gamma(1)+(pi*2)
      i=1
      ii=1
649   if (typ(i).ne.2) then
         i=i+1
         goto 649
      endif
      if (angz(i).lt.gamma(1)) then
         i=i+1
         goto 649
      endif
      start=i
648   if (angz(i).ge.gamma(tel+1)) then
         tel=tel+1
      goto 648
      endif
      hulp=lambda(i)
      indhulp=i
      do 650 i=1,start
         angz(i)=angz(i)+(pi*2)
650   continue
      i=start+1
651   if (typ(i).ne.2) then
         i=i+1
         if (i.eq.n+1) i=1
         goto 651
      endif
      if (i.eq.start) then
         if (tel2.gt.0) then

      if (xdev.gt.eps) then
         xcord=(px(indhulp,1)+tukmed(1))*xdev+xmean
      else
         xcord=px(indhulp,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(py(indhulp,1)+tukmed(2))*ydev+ymean
      else
         ycord=py(indhulp,1)+tukmed(2)
      endif
          pxpy(ii,1)=xcord
          pxpy(ii,2)=ycord
          pxpy(ii,3)=indhulp
          ii=ii+1
         endif
         goto 660
      endif
      if (angz(i).lt.gamma(tel+1)) then
         tel2=tel2+1
         if (lambda(i).gt.hulp) then
            hulp=lambda(i)
            indhulp=i
         endif
         i=i+1
         if (i.eq.n+1) i=1
         goto 651
      else
      if (tel2.gt.0) then

      if (xdev.gt.eps) then
         xcord=(px(indhulp,1)+tukmed(1))*xdev+xmean
      else
         xcord=px(indhulp,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(py(indhulp,1)+tukmed(2))*ydev+ymean
      else
         ycord=py(indhulp,1)+tukmed(2)
      endif
          pxpy(ii,1)=xcord
          pxpy(ii,2)=ycord
          pxpy(ii,3)=indhulp
          ii=ii+1
      endif
        tel=tel+1
        if (tel.eq.kount1+kount+1) tel=1
        if (angz(i).ge.gamma(tel+1)) then
           tel=tel+1
           if (tel.eq.kount1+kount+1) tel=1
        endif
        hulp=lambda(i)
        indhulp=i
        tel2=0
        goto 651
      endif
660   continue
      do 661 i=1,start
         angz(i)=angz(i)-(pi*2)
661   continue
C
C     Calculation of star-shaped whiskers.
C
      if (ntot.gt.nsub) goto 800
      if (whisk.ne.4) goto 800

      angz(n+1)=angz(1)
      do 669 i=1,n
         dattm(i)=0
669   continue
      do 670 i=1,n
         if ((dabs(angz(i+1)-angz(i)).lt.eps).and.
     +       (typ(i+1).eq.2).and.(typ(i).eq.2)) then
            if (lambda(i+1).gt.lambda(i)) dattm(i)=1
            if (lambda(i+1).le.lambda(i)) dattm(i+1)=1
         endif
670   continue
      dattm(1)=dattm(n+1)
      do 675 i=1,start
         angz(i)=angz(i)+(pi*2)
675   continue
      tel=1
      j=1
      if (dattm(start).eq.0) then
         i=start
      else
676      i=start+1
         if (i.eq.n+1) i=1
         if (typ(i).ne.2) goto 676
         if (dattm(i).eq.1) goto 676
      endif
      start=i
      ii=i
680   if (angz(i)-(pi*2).ge.gamma(tel+1)) then
         tel=tel+1
      goto 680
      endif
      do 690 l=1,tel
         gamma(l)=gamma(l)+(pi*2)
690   continue
      starttel=tel
      star(j,1)=x(index(i))
      star(j,2)=y(index(i))
      istarx=1
      if (xdev.gt.eps) then
         xcord=(star(j,1)+tukmed(1))*xdev+xmean
      else
         xcord=star(j,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(star(j,2)+tukmed(2))*ydev+ymean
      else
         ycord=star(j,2)+tukmed(2)
      endif
         star(j,1)=xcord
         star(j,2)=ycord
700   ii=ii+1
      if (ii.eq.n+1) ii=1
      if (typ(ii).ne.2) goto 700
      if (dattm(ii).eq.1) goto 700
      if (angz(ii).lt.gamma(tel+1)) then
      j=j+1
      star(j,1)=x(index(ii))
      star(j,2)=y(index(ii))
      istarx=1
      if (xdev.gt.eps) then
         xcord=(star(j,1)+tukmed(1))*xdev+xmean
      else
         xcord=star(j,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(star(j,2)+tukmed(2))*ydev+ymean
      else
         ycord=star(j,2)+tukmed(2)
      endif
         star(j,1)=xcord
         star(j,2)=ycord
      i=ii
      if (i.eq.start) goto 800
      goto 700
      else
730   tel=tel+1
      if (tel.eq.kount+kount1+1) tel=1
      IF (star(j,1).EQ.X(index(ii))) THEN
      ANG=PI2
      ELSE
      ANG=DATAN((star(j,2)-Y(index(ii)))/
     +(star(j,1)-X(index(ii))))
      IF (ANG.LE.0.0) ANG=ANG+PI
      ENDIF
      if (ang.lt.pi2) then
         if (x(index(ii)).lt.star(j,1)) ang=ang+pi
      else
         if (x(index(ii)).gt.star(j,1)) ang=ang+pi
      endif
720   IF (DSIN(ang)*wx(tel)-DCOS(ang)*wy(tel)
     +    .le.dsin(ang)*star(j,1)-dcos(ang)*star(j,2)) THEN
         if ((ii.eq.start).and.(tel.eq.starttel))
     +         gamma(tel+1)=gamma(tel+1)+(2*pi)
         if ((angz(ii).ge.gamma(tel+1))) then
            tel=tel+1
      if (tel.eq.kount+kount1+1) tel=1
            goto 720
         endif
      j=j+1
      star(j,1)=x(index(ii))
      star(j,2)=y(index(ii))
      istarx=1
      if (xdev.gt.eps) then
         xcord=(star(j,1)+tukmed(1))*xdev+xmean
      else
         xcord=star(j,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(star(j,2)+tukmed(2))*ydev+ymean
      else
         ycord=star(j,2)+tukmed(2)
      endif
         star(j,1)=xcord
         star(j,2)=ycord
      i=ii
      if (i.eq.start) goto 800
      goto 700
      ELSE
      j=j+1
740   star(j,1)=wx(tel)
      star(j,2)=wy(tel)
      if ((istarx.eq.1).and.(gamma(tel+1).lt.angz(ii))) then
         IF (x(index(i)).EQ.wx(tel)) THEN
         ANG=PI2
         ELSE
         ANG=DATAN((y(index(i))-wy(tel))/(x(index(i))-wx(tel)))
         IF (ANG.LE.0.0) ANG=ANG+PI
         ENDIF
         if (ang.lt.pi2) then
            if (wx(tel).lt.x(index(i))) ang=ang+pi
         else
            if (wx(tel).gt.x(index(i))) ang=ang+pi
         endif
         IF (DSIN(ang)*wx(tel+1)-DCOS(ang)*wy(tel+1)
     +    .gt.dsin(ang)*wx(tel)-dcos(ang)*wy(tel)) THEN
         tel=tel+1
         goto 740
         ELSE
         ENDIF
      endif
      istarx=0
      if (xdev.gt.eps) then
         xcord=(star(j,1)+tukmed(1))*xdev+xmean
      else
         xcord=star(j,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(star(j,2)+tukmed(2))*ydev+ymean
      else
         ycord=star(j,2)+tukmed(2)
      endif
         star(j,1)=xcord
         star(j,2)=ycord
      goto 730
      ENDIF
      endif
      else
      nointer=2
      endif
      else
      do 801 i=1,n
      datatyp(i,1)=(x(i)+tukmed(1))*xdev+xmean
      datatyp(i,2)=(y(i)+tukmed(2))*ydev+ymean
      datatyp(i,3)=0.0
801   continue
      endif
800   continue
610   num=kount1+kount
      RETURN
      END

      INTEGER FUNCTION NBP_NCEIL(M,J)

      INTEGER M,J

      IF (MOD(M,J).EQ.0) THEN
         NBP_NCEIL=INT(dble(M)/J)
      ELSE
         NBP_NCEIL=NINT(dble(M)/J+0.5)
      ENDIF
      RETURN
      END
