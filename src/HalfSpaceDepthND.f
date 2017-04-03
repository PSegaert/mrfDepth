      SUBROUTINE HSDEPNP(T,N1,X,N2,NP,NDIR,HDEP,NSIN,ERR,
     +                     NReduVar,ReduVar,UsedNP)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is a wrapper to calculate the approximate the halfspace
C depth of NP-dimensional points in a matrix T of size N1xNP 
C in an NP-dimensional data set X of size N2xNP.
C Missing values are not allowed.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C T        : DOUBLE PRECISION (N1xNP)            
C          : Coordinates of the points to calculate depth 
C N1       : integer(1)
C          : Number of observations to calculate depth of
C X        : DOUBLE PRECISION (N2xNP)
C          : The full data set.
C N2       : integer(1)
C          : The number of points in the dataset
C NP       : integer(1)
C          : The number of variables in the data set.
C NDIR     : integer(1)
C          : The number of samples to draw.
C HDEP     : DOUBLE PRECISION (N1)
C          : Final approximation for the halfspace depth.
C NSIN     : DOUBLE PRECISION (N1)
C          : Integer indicating if how many samples were singular
C ERR      : DOUBLE PRECISION (N1)
C          : Flag indicating if all samples were singular
C NReduVar : integer (1)
C          : Flag indicating how many variables have variance zero
C ReduVar  : integer (NP)
C          : Vector of which the first NReduVar elements contain
C          : the variables with zero variance
C UsedNP   : integer (N1)
C          : The final amount of dimension used.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
C INPUT
      INTEGER N1,N2,NP
      DOUBLE PRECISION T(N1,NP),X(N2,NP),HDEP(N1)
      DOUBLE PRECISION XTEMP(N2,NP)
      INTEGER NDIR
      DOUBLE PRECISION EPS
C OUTPUT
      INTEGER NSIN(N1), ERR(N1) 
      INTEGER NReduVar,ReduVar(NP),UsedNP(N1)
C WORKSPACE
      INTEGER I,J           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine 
      EPS=1.D-8    
      DO 10 I=1,N1
      DO 5 J=1,N2
        XTEMP(J,:)=X(J,:)
 5    CONTINUE
        CALL HSDEPNP1(XTEMP,N2,NP,NDIR,T(I,:),EPS,HDEP(I),
     +                NSIN(I),ERR(I),NReduVar,ReduVar,UsedNP(I))
 10   CONTINUE
      
      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE HSDEPNP1(X,N,NP,NDIR,T,EPS,HDEP,NSIN,ERR,
     +                     NReduVar,ReduVar,UsedNP)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is calculates the approximate the halfspace
C depth of NP-dimensional points in a vector T of size NP 
C in an NP-dimensional data set X of size NxNP.
C Actual calculations are done in subroutine HSDEPNP2.
C Missing values are not allowed.
C  This subroutine was described in:
C       Rousseeuw, P.J and Struyf,A.:  
C       Computing location depth and regression depth in higher dimensions.
C       Statistics and Computing 8, 193-203 (1998).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C X      : DOUBLE PRECISION (NxNP)
C        : The full data set.
C N      : integer(1)
C        : The number of points in the dataset
C NP     : integer(1)
C        : The number of variables in the data set.
C NDIR   : integer(1)
C        : The number of samples to draw.
C T      : DOUBLE PRECISION (NP)            
C        : Coordinates of the points to calculate depth 
C EPS    : DOUBLE PRECISION
C        : Numerical Precision
C HDEP   : DOUBLE PRECISION 
C        : Final approximation for the halfspace depth.
C NSIN   : DOUBLE PRECISION
C        : Integer indicating how many samples were singular
C ERR    : DOUBLE PRECISION
C        : Flag indicating if all samples were singular
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NP,MAXN,MAXP,NNP
      DOUBLE PRECISION X(N,NP),T(NP),HDEP,XN(N)
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP),R(NP)
      DOUBLE PRECISION AVE(NP),COV(NP,NP)
      INTEGER NDEP,NDIR
      DOUBLE PRECISION EPS
      INTEGER JSAMP(NP)
      INTEGER NSIN, ERR
      INTEGER NReduVar,ReduVar(NP),UsedNP

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine   
      UsedNP=0
C     Initialize old variables in compliance
      MAXN=N
      MAXP=NP
      NNP=NP

C     Center the data around point T
      CALL STAND_HSDEPTHND(MAXN,MAXP,N,NP,X,T,XN,EPS,NDEP,
     +                     NReduVar,ReduVar)

      IF (NP.NE.NNP) THEN
         NNP=NP
         IF ((NDEP.EQ.0).OR.(NP.EQ.0)) then
            UsedNP=MAXP
            GOTO 10
         ENDIF
      ENDIF
C     Initialize NDEP
      NDEP=N
      CALL HSDEPNP2(N,NP,NNP,NDIR,MAXN,MAXP,X,JSAMP,T,R,
     +     EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN)
      UsedNP=NNP
 10   IF (NSIN.EQ.NDIR) THEN
         ERR=1
      ELSEIF (NSIN.GT.(0-EPS)) THEN
         ERR=0
      ELSE 
         ERR=-1
      ENDIF
      HDEP=(NDEP+0.D0)/(N+0.D0)
      RETURN
 300  END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC








      SUBROUTINE HSDEPNP2(N,NP,NNP,NDIR,MAXN,MAXP,X,JSAMP,T,R,
     +     EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is to be called by HSDEPNP1.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER MAXN,MAXP,NP,NDEP,NSIN
      DOUBLE PRECISION X(MAXN,MAXP),T(NP),R(NP),EPS
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP),COV(NP,NP),AVE(NP)
      INTEGER JSAMP(NP)
      INTEGER N,J,NNP,NUMH,NT,L
      INTEGER NDIR,NNP1,IERR

      DOUBLE PRECISION BETA(N),DPF(N),SDEP,WX(N),WY(N)
      INTEGER F(N),JLV(N),JRV(N),NDIM

C  Initialize the number of singular samples.
      NSIN=0
C  Handle special case where N is equal to 1.
      IF (N.EQ.1)THEN
         IF (N.EQ.1) THEN
            DO 10 J=1,NP
               IF (DABS(X(1,J)-T(J)).GT.EPS) GOTO 15
 10         CONTINUE
            NDEP=1
            RETURN
         ENDIF
 15      NDEP=0
         RETURN
      ENDIF
C  Handle special case where NNP is equal to 1.
 25   IF (NNP.EQ.1) THEN
         NUMH=0
         NT=0
         DO 20 L=1,N
            IF (X(L,1).GT.(T(1)+EPS)) THEN
               NUMH=NUMH+1
            ELSEIF (X(L,1).GE.(T(1)-EPS)) THEN
               NT=NT+1
            ENDIF
 20      CONTINUE
         NDEP=MIN0(NUMH+NT,N-NUMH)
         RETURN
      ENDIF
      IF (NNP.EQ.2) THEN
        CALL HSDEP21(T(1),T(2),N,X(:,1),X(:,2),
     +               BETA,F,DPF,JLV,JRV,NDEP,SDEP)
            GOTO 50
      ENDIF
      IF (NNP.EQ.3) THEN
        CALL HSDEPTH31(N,T(1),T(2),T(3),X(:,1),X(:,2),X(:,3),
     +               BETA,F,WX,WY,EPS,NDIM,NDEP)
            NNP=NDIM
            GOTO 50
      ENDIF

C  General case: call subroutine DEP.
      CALL DEP(N,NNP,NDIR,MAXN,MAXP,X,JSAMP,T,R,
     +     EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN)
C  If all points and theta are identified as lying on the same hyperplane,
C  reduce the dimension of the data set by projection on that hyperplane,
C  and compute the depth on the reduced data set.
      IF (NSIN.EQ.(-1)) THEN
         NSIN=0
         NNP1=NNP
         NNP=NNP-1
         CALL REDUCE(N,NNP,NNP1,MAXN,MAXP,X,T,R,EVECS,JSAMP,IERR)
         IF (IERR.LT.0) THEN
            GOTO 50
         ENDIF
         GOTO 25
      ENDIF

 50   RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC








      SUBROUTINE DEP(N,NP,NDIR,MAXN,MAXP,X,JSAMP,T,R,
     +     EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is to be called by HSDEPNP2.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER MAXN,MAXP,N,NP
      DOUBLE PRECISION X(MAXN,MAXP),T(NP),R(NP)
      DOUBLE PRECISION K,KT,EPS
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP),COV(NP,NP),AVE(NP)
      INTEGER JSAMP(NP)
      REAL(8) RAN(1)
      INTEGER ISEED
      INTEGER NSAMP,L,J,NSIN,NDEP,NRAN,NDIR
      INTEGER I,IERR,NT,NUMH
C  Initialize halfspace depth and random seed.
      NSIN=0
      NDEP=N
      ISEED = 256
      DO 100 NRAN=1,NDIR
C  Draw a random sample of size np.  
         CALL UNIRAN(1,ISEED,RAN)
         I=N*RAN(1)+1.
         IF(I.GT.N)I=N
         JSAMP(1)=I
         NSAMP=1
 20      CALL UNIRAN(1,ISEED,RAN)
         L=N*RAN(1)+1.
         IF(L.GT.N)L=N
         DO 30 J=1,NSAMP
            IF(L.EQ.JSAMP(J)) GOTO 20
 30      CONTINUE
         NSAMP=NSAMP+1
         JSAMP(NSAMP)=L
         IF (NSAMP.LT.NP)GOTO 20
C  Compute the covariance matrix of the sample.
         DO 40 J=1,NP
            AVE(J)=0.D0
            DO 50 I=1,NP
                  AVE(J)=AVE(J)+X(JSAMP(I),J)
 50         CONTINUE
            AVE(J)=AVE(J)/NP
 40      CONTINUE
         DO 60 J=1,NP
            DO 70 L=1,J
               COV(J,L)=0.D0
               DO 80 I=1,NP
                  COV(J,L)=COV(J,L)+(X(JSAMP(I),J)-AVE(J))
     +                 *(X(JSAMP(I),L)-AVE(L))
 80            CONTINUE
               COV(J,L)=COV(J,L)/(NP-1)
               COV(L,J)=COV(J,L)
 70         CONTINUE
 60      CONTINUE
C  Compute the eigenvalues and corresponding eigenvectors 
C  of the covariance matrix.
         CALL EIGEN(NP,NP,COV,EVALS,EVECS,R,AVE,IERR)
         IF (IERR.NE.0) THEN
C            WRITE(*,200)IERR
            NSIN=NSIN+1
            GOTO 100
         ENDIF
         IF (EVALS(1).GT.EPS) THEN
            NSIN=NSIN+1
            GOTO 100
         ENDIF
C  Test for singularity of the sample.
         IF (EVALS(2).LE.EPS) THEN
            NSIN=NSIN+1
         ENDIF
C  Project all points on the line through theta with direction given by 
C  the eigenvector of the smallest eigenvalue, i.e. the direction 
C  orthogonal on the hyperplane given by the np-subset.
C  Compute the one-dimensional halfspace depth of theta on this line.
         KT=0.D0
         NT=0
         DO 90 J=1,NP
            IF (DABS(EVECS(J,1)).LE.EPS) THEN
               NT=NT+1
            ELSE
               KT=KT+T(J)*EVECS(J,1)
            ENDIF
 90      CONTINUE
         IF (NT.EQ.NP) THEN
            NSIN=NSIN+1
            GOTO 100
         ENDIF

         NUMH=0
         NT=0
         DO 95 L=1,N
            K=0.D0
            DO 96 J=1,NP
               K=K+EVECS(J,1)*X(L,J)
 96         CONTINUE
            K=K-KT
            IF (K .GT. EPS) THEN
               NUMH=NUMH+1
            ELSEIF (K .GE. (0.D0-EPS)) THEN
               NT=NT+1
            ENDIF
 95      CONTINUE
C  If all projections collapse with theta, return to reduce the dimension.
         IF (NT.EQ.N) THEN
            NSIN=-1
            RETURN
         ENDIF
C  Update the halfspace depth.
         NDEP=MIN0(NDEP,MIN0(NUMH+NT,N-NUMH))
 100  CONTINUE
      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC








      SUBROUTINE STAND_HSDEPTHND(MAXN,MAXP,N,NP,X,T,XN,EPS,NDEP,
     +                     NReduVar,ReduVar)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine is to be called by HSDEPNP1 and centers the data X
C around P.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER MAXN,MAXP,N,NP
      DOUBLE PRECISION X(MAXN,MAXP),T(NP),XN(N),EPS
      DOUBLE PRECISION QLOC,QSCA,AVE,VAR
      DOUBLE PRECISION FINDQ
      INTEGER JN,J,I,NDEP
      INTEGER NReduVar,ReduVar(NP)

      NReduVar=0

      JN=0
      DO 10 J=1,NP
         ReduVar(J)=0
         DO 20 I=1,N
            XN(I)=X(I,J)
 20      CONTINUE
         IF ((2*INT(N/2)).EQ.N) THEN
            QLOC=FINDQ(XN,N,N/2)
            QLOC=(FINDQ(XN,N,(N/2)+1)+QLOC)/2.D0
         ELSE
            QLOC=FINDQ(XN,N,INT(N/2)+1)
         ENDIF
         DO 30 I=1,N
            XN(I)=DABS(X(I,J)-QLOC)
 30      CONTINUE
         IF ((2*INT(N/2)).EQ.N) THEN
            QSCA=FINDQ(XN,N,N/2)
            QSCA=(FINDQ(XN,N,(N/2)+1)+QSCA)/2.D0
         ELSE
            QSCA=FINDQ(XN,N,INT(N/2)+1)
         ENDIF
         IF (DABS(QSCA).LT.EPS) THEN
            AVE=0.D0
            DO 40 I=1,N
               AVE=AVE+X(I,J)
 40         CONTINUE
            AVE=AVE/(N+0.D0)
            VAR=0.D0
            DO 50 I=1,N
               VAR=VAR+(X(I,J)-AVE)*(X(I,J)-AVE)
 50         CONTINUE  
            IF (N.NE.1) VAR=VAR/(N-1.D0)
            IF (DABS(VAR).LT.EPS) THEN
               IF (DABS(T(J)-X(1,J)).GT.EPS) NDEP=0
               NP=NP-1
               NReduVar=NReduVar+1
               ReduVar(NReduVar)=J
               GOTO 10
            ELSE
               QSCA=DSQRT(VAR)
            ENDIF
         ENDIF
         JN=JN+1
         DO 60 I=1,N
            X(I,JN)=(X(I,J)-QLOC)/QSCA
 60      CONTINUE         
         T(JN)=(T(J)-QLOC)/QSCA
 10   CONTINUE

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
