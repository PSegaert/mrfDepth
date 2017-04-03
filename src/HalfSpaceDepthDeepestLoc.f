      SUBROUTINE HSDEPTH_DEEPEST(X,N,np,
     +maxdir,nstp,ntry,nalt,
     +dpstM,HDEP,AlgStopFlag,ndir,nstep)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is a wrapper to calculate the approximate deepest 
C location in an NP-dimensional data set X.
C This routine was described in
C       Struyf, A. and Rousseeuw, P.J.:  
C       High Dimensional computation of the deepest location
C       Computational Statistics and Data Analysis 34, 415-426 (2000).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C X           : DOUBLE PRECISION (NxNP)
C             : The full data set.
C N           : INTEGER (1)
C             : The number of points in the dataset X.
C NP          : INTEGER (1)
C             : The number of variables in the data set X.
C maxdir      : INTEGER (1)
C             : Max number of directions used in the approximation.
C nstp        : INTEGER (1)
C             : max number of steps
C ntry        : INTEGER (1)
C             : max number of trials
C nalt        : INTEGER (1)
C             : max number of steps tried to increase the depth along a 
C		      : given direction.
C DPSTM       : DOUBLE PRECISION (NP)
C             : Coordinates of the approximate deepest point
C HDEP        : DOUBLE PRECISION (1)
C             : Approximate depth of dpstM
C AlgStopFlag : INTEGER (1)
C             : Indicates the stoping criterium of the algorithm.
C NDIR        : INTEGER (1)
C             : Number of directions used
C NSTEP       : INTEGER (1)
C             : Number of step used
C EPS         : DOUBLE PRECISION (1)
C             : Numerical Tolerance (10^-8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      implicit double precision (a-h,o-z)
C     INPUT
      INTEGER N,np
      double precision X(N,np)
      INTEGER maxdir,nstp,ntry,nalt
C     OUTPUT
      double precision dpstM(np),HDEP
      INTEGER nddpst,AlgStopFlag,ndir,nstep
C     WORKSPACE
      double precision xn(N),stepsM(2*(np+2),np)
      double precision cov(np,np),ave(np),u(maxdir,np)
      double precision evecs(np,np),evals(np)
      double precision utx(N,maxdir),utxsort(N,maxdir),d1(np)
      double precision locsca(np,2)
      integer jsamp(np), i1(N),i2(N)
      integer J
      double precision EPS
      ndir=maxdir
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine  
      EPS=1.D-8

      call stand(N,np,n,np,x,xn,eps,locsca)
      call deepest(n,np,ndir,x,eps,nddpst,dpstM,
     +     stepsM,xn,jsamp,cov,ave,evals,evecs,u,utx,utxsort,i1,i2,d1,
     +        nstp,ntry,nalt,AlgStopFlag,nstep)
      do 1 j=1,np
         dpstM(j)=dpstM(j)*locsca(j,2)+locsca(j,1)
 1    continue
      
      HDEP=(nddpst+0.d0)/(N+0.d0)

      RETURN
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      subroutine deepest(n,np,ndir,x,eps,nddpst,dpstM,
     +     stepsM,xn,jsamp,cov,ave,evals,evecs,u,utx,utxsort,i1,i2,d1,
     +     nstp,ntry,nalt,AlgStopFlag,nstep)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  To be called by HSDEPTH_DEEPEST
C  This subroutine does the actual computations
C  This subroutine was described in:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit double precision (a-h,o-z)
      INTEGER N,NP,ndir,nstp,ntry,nstep
      double precision x(N,NP),eps,stepsM(2*(np+2),np),xn(n)
      double precision cov(np,np),ave(np),evecs(np,np),evals(np),d1(np)
      double precision utx(n,ndir),utxsort(n,ndir),dpstM(np),u(ndir,np)
      double precision dpmedian
      integer jsamp(np),i1(n),i2(n)
      REAL(8) ran(1)
      INTEGER AlgStopFlag
      INTEGER ISEED
      INTEGER ntave,nalt,nback,lj,mold
      INTEGER I,J,L,NDDPST,numh,nt,mindep
      INTEGER nceil2,indmold,nrankl,nrankg
      INTEGER nsa,ngen,ind,jdir,nmin,nmax,nsin,p,ji,nsamp
      INTEGER ierr,indm,nid,nidalt,ndstep,ndold,inalt
      double precision utj,dmax,utdpst

cc  initialize the random seed.
      ISEED=256
cc  handle special case where n is equal to 1.
      if (n.eq.1)then
         do 1 j=1,np
            dpstM(j)=x(1,j)
 1       continue
         nddpst=1
         return
      endif
cc  handle special case where np is equal to 1.
      if (np.eq.1) then
         do 2 l=1,n
            xn(l)=x(l,1)
 2       continue
         dpstM(1)=dpmedian(xn,n)
         numh=0
         nt=0
         do 3 l=1,n
            if (x(l,1).gt.(dpstM(1)+eps)) then
               numh=numh+1
            elseif (x(l,1).ge.(dpstM(1)-eps)) then
               nt=nt+1
            endif
 3       continue
         nddpst=min0(numh+nt,n-numh)
         return
      endif
cc  general case.
C-------------
C Step 0
C-------------
cc  initialize the minimal depth of the deepest location.
      mindep=nceil2(dble(n/(np+1)),eps)
C-------------
C Step 1
C-------------      
cc compute the coordinate-wise median, used as first approximation 
cc for the deepest location.
      do 4 j=1,np
         do 5 i=1,n
            xn(i)=x(i,j)
 5       continue
         dpstM(j)=dpmedian(xn,n)
 4    continue
C-------------
C Step 2
C-------------
cc  construct ndir unit vectors: these are the directions used to compute
cc  depths and moving directions.
cc  a) add coordinate axes
      do 6 l=1,np
         do 7 j=1,np
            if (j.eq.l) then
               u(l,j)=1.d0
            else
               u(l,j)=0.d0
            endif
 7       continue
 6    continue

cc  b) add vectors connecting data points with the coordinate-wise median
cc  ( first check whether all directions can be used 
cc  or random selection is required. )
      if (n.le.int(ndir/4)) then
         nsa=n
         ngen=0
      else
         nsa=int(ndir/4)
         ngen=1
      endif
      ind=np

      do 8 jdir=1,nsa
         if (ngen.eq.1) then
            call UNIRAN(1,ISEED,ran)
            i=int(n*ran(1)+1.d0)
            if(i.gt.n)i=n
            jsamp(1)=i
         else
            jsamp(1)=jdir
         endif
         utj=0.d0
         do 9 j=1,np
            u(ind+1,j)=x(jsamp(1),j)-dpstM(j)
            utj=utj+u(ind+1,j)*u(ind+1,j)
 9       continue
         utj=dsqrt(utj)
         if (utj.gt.eps) then
            ind=ind+1
            do 10 j=1,np
               u(ind,j)=u(ind,j)/utj
 10         continue
         endif
 8    continue

cc  c) add vectors connecting two data points
cc  ( first check whether all directions can be used 
cc  or random selection is required. )
      if ((n*(n-1)/2).le.(int(ndir/2)-ind)) then
         nsa=(n*(n-1))/2
         ngen=0
      else
         nsa=int(ndir/2)-ind
         ngen=1
      endif
      nmin=1
      nmax=2

      do 11 jdir=1,nsa
         if (ngen.eq.1) then
            call UNIRAN(1,ISEED,ran)
            i=int(n*ran(1)+1.d0)
            if(i.gt.n)i=n
            jsamp(1)=i
 12         call UNIRAN(1,ISEED,ran)
            l=int(n*ran(1)+1.d0)
            if(l.gt.n)l=n
            if(l.eq.jsamp(1)) goto 12
            jsamp(2)=l
         else
            jsamp(1)=nmin
            jsamp(2)=nmax
            if (nmax.eq.n) then
               nmin=nmin+1
               nmax=nmin+1
            else
               nmax=nmax+1
            endif
         endif
         utj=0.d0
         do 13 j=1,np
            u(ind+1,j)=x(jsamp(1),j)-x(jsamp(2),j)
            utj=utj+u(ind+1,j)*u(ind+1,j)
 13      continue
         utj=dsqrt(utj)
         if (utj.gt.eps) then
            ind=ind+1
            do 14 j=1,np
               u(ind,j)=u(ind,j)/utj
 14        continue
         endif
 11   continue

cc  d) add vectors perpendicular to hyperplanes through np data points
cc  ( first check whether all directions can be used 
cc  or random selection is required. )
      nsin=0
      dmax=1.d0
      if (np.gt.int(n/2)) then 
         p=n-np
      else
         p=np
      endif
      j=p
 15   dmax=dmax*dble(n-p+j)/dble(p-j+1)
      if (dmax.gt.dble(ndir-ind)) then
         nsa=ndir-ind
         ngen=1
         goto 21
      endif
      j=j-1
      if (j.ge.1) goto 15
      nsa=int(dmax)
      ngen=0

 21   do 100 jdir=1,nsa
         if (ngen.eq.0) then
            if (jdir.eq.1) then
               do 16 j=1,np
                  jsamp(j)=j
 16            continue
            else
               j=np
 17            if (jsamp(j).lt.(n-np+j)) then
                  jsamp(j)=jsamp(j)+1
                  do 18 ji=(j+1),np
                     jsamp(ji)=jsamp(ji-1)+1
 18               continue
               else
                  j=j-1
                  goto 17
               endif
            endif
         else
            call UNIRAN(1,ISEED,ran)
            i=int(n*ran(1)+1.d0)
            if(i.gt.n)i=n
            jsamp(1)=i
            nsamp=1
 19         call UNIRAN(1,ISEED,ran)
            l=int(n*ran(1)+1.d0)
            if(l.gt.n)l=n
            do 20 j=1,nsamp
               if(l.eq.jsamp(j)) goto 19
 20         continue
            nsamp=nsamp+1
            jsamp(nsamp)=l
            if (nsamp.lt.np)goto 19
         endif
cc  compute the covariance matrix of the sample.
         do 30 j=1,np
            ave(j)=0.d0
            do 40 i=1,np
               ave(j)=ave(j)+x(jsamp(i),j)
 40         continue
            ave(j)=ave(j)/np
 30      continue
         do 50 j=1,np
            do 60 l=1,j
               cov(j,l)=0.d0
               do 70 i=1,np
                  cov(j,l)=cov(j,l)+(x(jsamp(i),j)-ave(j))
     +                 *(x(jsamp(i),l)-ave(l))
 70            continue
               cov(j,l)=cov(j,l)/(np-1)
               cov(l,j)=cov(j,l)
 60         continue
 50      continue
cc  compute the eigenvalues and corresponding eigenvectors 
cc  of the covariance matrix.
         call eigen(np,np,cov,evals,evecs,d1,ave,ierr)
         if (ierr.ne.0) then
C            write(*,700)ierr
            nsin=nsin+1
            goto 100
         endif
         if (evals(1).gt.eps) then
C            write(*,800)jdir
            nsin=nsin+1
            goto 100
         endif
cc  test for singularity of the sample.
         if (evals(2).le.eps) then
            nsin=nsin+1
         endif
cc  determine the direction orthogonal to the sample.
         nt=0
         do 80 j=1,np
            if (dabs(evecs(j,1)).le.eps) nt=nt+1
 80      continue
         if (nt.eq.np) then
C           write(*,900)jdir
            nsin=nsin+1
            goto 100
         endif
         ind=ind+1
         do 90 j=1,np
            u(ind,j)=evecs(j,1)
 90      continue
 100  continue   
cc search the deepest location.
cc initialize.
C-------------
C Step 4
C-------------
      ndir=ind
      do 101 jdir=1,ndir
         do 102 i=1,n
            utx(i,jdir)=0.d0
            do 103 j=1,np
               utx(i,jdir)=utx(i,jdir)+u(jdir,j)*x(i,j)
 103        continue
            xn(i)=utx(i,jdir)
 102        continue
         call sortLoc(xn,n,i1,i2)
         do 104 i=1,n
            utxsort(i,jdir)=xn(i)
 104     continue
 101  continue
      do 105 j=1,np
         stepsM(1,j)=dpstM(j)
 105  continue
C-------------
C Step 3
C-------------
      indM=1
      nid=0
      nidalt=0
      ndstep=-1
      ndold=-1
      nt=0
      nstep=0
      inalt=2
C      write(*,*)ndir,' directions will be evaluated.'

C-------------
C Step 5
C-------------
cc Start iterations.
 110  ndold=max0(ndstep,ndold)
      nstep=nstep+1
      ndstep=n+1
      do 120 l=1,ndir
cc Compute rank of new deepest in direction l.
         utj=0.d0
         do 130 j=1,np
            utj=utj+u(l,j)*stepsM(indM,j)
 130     continue
         do 140 i=1,n
            xn(i)=utxsort(i,l)
 140     continue
         call irank(utj,xn,n,eps,nrankl,nrankg)
         if (nrankg.lt.ndstep) then
            ndstep=nrankg
            do 150 j=1,np
               ave(j)=u(l,j)
 150        continue
            lj=l
            nt=1
         elseif (nrankg.eq.ndstep) then
            do 160 j=1,np
               ave(j)=ave(j)+u(l,j)
 160        continue            
            nt=nt+1
         endif
         if (nrankl.lt.ndstep) then
            ndstep=nrankl
            do 170 j=1,np
               ave(j)=-u(l,j)
 170        continue
            lj=-l
            nt=1
         elseif (nrankl.eq.ndstep) then
            do 180 j=1,np
               ave(j)=ave(j)-u(l,j)
 180        continue            
            nt=nt+1
         endif
 120  continue
      ntave=0
      do 190 j=1,np
         ave(j)=ave(j)/(nt+0.d0)
         if (dabs(ave(j)).lt.eps) ntave=ntave+1
 190  continue            
      if (ntave .eq. np) then
         do 195 j=1,np
            if (lj.gt.0) ave(j)=u(lj,j)
            if (lj.lt.0) ave(j)=-u(-lj,j)
 195     continue
      endif
C-------------
C Step 6
C-------------
      if (ndstep.ge.ndold) then
cc New deepest point found.
         nddpst=ndstep
C        c)
         do 200 j=1,np
            dpstM(j)=stepsM(indM,j)
 200     continue
C        d)
         if (inalt.ne.2) inalt=2
C        a)
         if (ndold.eq.ndstep) then
            nid=nid+1
            nidalt=nidalt+1
C        b)
         else 
            nid=0
            nidalt=0
         endif

      elseif (inalt.ne.2) then 
cc Alternative directions loop is active, but a bad direction was tried, 
cc hence try another direction.
         nid=nid+1
         goto 220
      else
         nidalt=nidalt+1
         nid=nid+1
      endif
C-------------
C Step 7
C-------------
      if (nidalt.ge.nalt) then
cc Save present configuration, and enter alternative directions loop.
         nback=ndstep
         do 210 j=1,np
            evals(j)=ave(j)
            evecs(j,1)=stepsM(indM,j)
 210     continue
         goto 220
      endif
      
      goto 300

cc Try some alternative directions, connecting the latest approximation
cc of the deepest location with one of the earlier approximations.
 220  nidalt=0
      if (inalt.eq.2*(np+2)) then
cc Restore old configuration.
         ndstep=nback
         if (indM.eq.2*(np+2)) then
            indM=1
         else
            indM=indM+1
         endif
         do 230 j=1,np
            ave(j)=evals(j)
            stepsM(indM,j)=evecs(j,1)
 230     continue            
         if (ndold.ge.ndstep) nid=nid+1
         inalt=2
         goto 300
      endif
      
      if (indM.eq.1) then
         indM=2*(np+2)
      else
         indM=indM-1
      endif
      if ((indM+inalt).gt.2*(np+2)) then
         mold=indM-2*(np+2)+inalt
      else
         mold=indM+inalt
      endif
      do 245 j=1,np
         ave(j)=stepsM(mold,j)-stepsM(indM,j)
 245  continue
      inalt=inalt+1
      goto 400
C-------------
C Step 8
C-------------      
 300  if (nid.ge.ntry) then
         AlgStopFlag=0
         return
      endif
C-------------
C Step 9
C-------------
cc Take a step in the computed direction.
 400  if (nstep.ge.nstp) then
         AlgStopFlag=1
         return
      endif
C-------------
C Step 10
C-------------
      utdpst=0.d0
      do 240 j=1,np
         utdpst=utdpst+ave(j)*ave(j)
 240  continue
      utdpst=dsqrt(utdpst)
      do 250 j=1,np
         ave(j)=ave(j)/utdpst
 250  continue
      utdpst=0.d0
      do 410 j=1,np
         utdpst=utdpst+ave(j)*stepsM(indM,j)
 410  continue
      do 420 i=1,n
         xn(i)=0.d0
         do 430 j=1,np
            xn(i)=xn(i)+ave(j)*x(i,j)
 430     continue
 420  continue
      call sortLoc(xn,n,i1,i2)
      call irank(utdpst,xn,n,eps,nrankl,nrankg)
      
      if (ndstep.ge.int(n/2)) then
         return
      elseif (nrankg.lt.mindep) then
         utdpst=-(utdpst-xn(n+1-mindep))
      else
         utdpst=-(utdpst-xn(n-nrankg))
      endif
      indMold=indM
      if (indM.eq.2*(np+2)) then
         indM=1
      else
         indM=indM+1
      endif   
      do 440 j=1,np
         stepsM(indM,j)=stepsM(indMold,j)+utdpst*ave(j)
 440  continue
      goto 110

      RETURN
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      subroutine sortLoc(b,n,jlv,jrv)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  Sorts an array b of length n in o(nlogn) time.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit double precision (a-h,o-z)
      integer n
      double precision b(n),amm,xx
      integer jlv(n),jrv(n)
      integer jss,jndl,jr,jnc,j,jtwe
      jss=1
      jlv(1)=1
      jrv(1)=n
  10  jndl=jlv(jss)
      jr=jrv(jss)
      jss=jss-1
  20  jnc=jndl
      j=jr
      jtwe=(jndl+jr)/2
      xx=b(jtwe)
  30  if (b(jnc).ge.xx) goto 40
      jnc=jnc+1
      goto 30
  40  if (xx.ge.b(j)) goto 50
      j=j-1
      goto 40
  50  if (jnc.gt.j) goto 60
      amm=b(jnc)
      b(jnc)=b(j)
      b(j)=amm
      jnc=jnc+1
      j=j-1
  60  if (jnc.le.j) goto 30
      if ((j-jndl).lt.(jr-jnc)) goto 80
      if (jndl.ge.j) goto 70
      jss=jss+1
      jlv(jss)=jndl
      jrv(jss)=j
  70  jndl=jnc
      goto 100
  80  if (jnc.ge.jr) goto 90
      jss=jss+1
      jlv(jss)=jnc
      jrv(jss)=jr
  90  jr=j
 100  if (jndl.lt.jr) goto 20
      if (jss.ne.0) goto 10
      return
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      subroutine irank(u,aw,n,eps,indle,indge)
      implicit double precision (a-h,o-z)
      integer N
      double precision aw(n),u,eps
      integer indge,indle,indl,indg,j,imax,imin

      if (u.lt.(aw(1)-eps)) then
         indge=n
         indle=0
         return
      elseif (u.le.(aw(1)+eps)) then
         indge=n
         indle=1
         indl=1
         goto 200
      endif
      if (u.gt.(aw(n)+eps)) then
         indge=0
         indle=n
         return
      elseif (u.ge.(aw(n)-eps)) then
         indge=1
         indle=n
         indg=n
         goto 100
      endif
      imin=1
      imax=n
 50   if ((imax-imin).eq.1) then
         indge=n-imin
         indle=imin
         return
      endif
      j=int((imax+imin)/2)
      if (u.lt.(aw(j)-eps)) then
         imax=j
      elseif (u.gt.(aw(j)+eps)) then
         imin=j
      else
         indge=n-j+1
         indle=j
         indl=j
         indg=j
         goto 100         
      endif
      goto 50
 100  if (dabs(aw(indg-1)-u).le.eps) then
         indge=indge+1
         indg=indg-1
         goto 100
      endif
      if (indle.eq.n) return
 200  if (dabs(aw(indl+1)-u).le.eps) then
         indle=indle+1
         indl=indl+1
         goto 200
      endif
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      function nceil2(u,eps)
      implicit double precision (a-h,o-z)
      integer nceil2
      double precision u,eps

      nceil2=int(u)
      if (dabs(dble(nceil2-u)).gt.eps) nceil2=nceil2+1
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      subroutine stand(maxn,maxp,n,np,x,xn,eps,locsca)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Standardises a dataset.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit double precision (a-h,o-z)
      integer n,maxn,maxp
      double precision x(maxn,maxp),xn(n),eps,findq
      double precision qloc,qsca,ave,var,locsca(maxp,2)
      integer i,np,j,jn
      
      jn=0
      do 10 j=1,np
         do 20 i=1,n
            xn(i)=x(i,j)
 20      continue
         if ((2*int(n/2)).eq.n) then
            qloc=findq(xn,n,n/2)
            qloc=(findq(xn,n,(n/2)+1)+qloc)/2.d0
         else
            qloc=findq(xn,n,int(n/2)+1)
         endif
         do 30 i=1,n
            xn(i)=dabs(x(i,j)-qloc)
 30      continue
         if ((2*int(n/2)).eq.n) then
            qsca=findq(xn,n,n/2)
            qsca=(findq(xn,n,(n/2)+1)+qsca)/2.d0
         else
            qsca=findq(xn,n,int(n/2)+1)
         endif
         if (dabs(qsca).lt.eps) then
            ave=0.d0
            do 40 i=1,n
               ave=ave+x(i,j)
 40         continue
            ave=ave/(n+0.d0)
            var=0.d0
            do 50 i=1,n
               var=var+(x(i,j)-ave)*(x(i,j)-ave)
 50         continue  
            if (n.ne.1) var=var/(n-1.d0)
            if (dabs(var).lt.eps) then
               if (np.ne.1) then
                  np=np-1
                  goto 10
               endif
            else
               qsca=dsqrt(var)
            endif
         endif
         jn=jn+1
         locsca(jn,1)=qloc
         locsca(jn,2)=qsca
         do 60 i=1,n
            x(i,jn)=(x(i,j)-qloc)/qsca
 60      continue         
 10   continue

      return
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      function dpmedian(aw,ncas)
cc  Finds the median of the array aw of length ncas.
      implicit double precision (a-h,o-z)
      integer ncas
      double precision dpmedian,findq
      double precision aw(ncas),qloc

      if ((2*int(ncas/2)).eq.ncas) then
         qloc=findq(aw,ncas,ncas/2)
         qloc=(findq(aw,ncas,(ncas/2)+1)+qloc)/2.d0
      else
         qloc=findq(aw,ncas,int(ncas/2)+1)
      endif
      dpmedian=qloc
      end






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






