      SUBROUTINE RDEPTHND(T,N1,X,N2,NP,NDIR,RDEP,NSIN,ERR)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  This subroutine computes an approximation for the regression depth
C  of fits T in an (NP)-dimensional data set X of size N2.
C  Missing values are not allowed.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Description of Arguments
C
C INPUT
C N1     : INTEGER (1)
C        : The number of hyperplanes to compute the depth off.
C T      : DOUBLE PRECISION (N1xNP)            
C        : Hyperplanes of which the depth is computed.
C        : The first NP-1 elements are the slopes in 
C        : the directions of the variables describing the X-space.
C X      : DOUBLE PRECISION (N2xNP)            
C        : The full data set. 
C N2     : INTEGER (1)
C        : Size of the data set
C NP     : INTEGER (1)
C        : The number of dimensions.
C NDIR   : INTEGER (1)
C        : The number of samples to draw.
C EPS    : DOUBLE PRECISION             
C        : Numerical Precision 
C     
C OUTPUT
C RDEP   : DOUBLE PRECISION (N1)
C        : The approximate regression depth of the planes	
C NSIN   : INTEGER (N1)
C        : Number of singular directions
C ERR    : INTEGER (N1)
C        : Indicating if all directions were singular (ERR=1) or not (ERR = 0)
	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      INTEGER N1,N2,NP,NDIR,I,J
      DOUBLE PRECISION X(N2,NP),XX(N2,NP),RDEP(N1),T(N1,NP)
      DOUBLE PRECISION EPS
      INTEGER NSIN(N1), ERR(N1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      
      EPS=1.D-8
      
      DO 10 I=1,N1
         DO 11 J=1,N2
         XX(J,:)=X(J,:)
11      CONTINUE  
       CALL RDEPTH_APPR1(T(I,:),XX,N2,NP-1,NDIR,
     +                   RDEP(I),EPS,NSIN(I),ERR(I))
 10   CONTINUE
      
      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE RDEPTH_APPR1(T,X,N,NP,NDIR,RDEP,EPS,NSIN,ERR)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  This subroutine computes an approximation for the regression depth
C  of fits T in an (NP)-dimensional data set X of size N2.
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
C T      : DOUBLE PRECISION (NP)            
C        : Hyperplane of which the depth is computed.
C        : The first NP-1 elements are the slopes in 
C        : the directions of the variables describing the X-space.
C X      : DOUBLE PRECISION (NxNP)            
C        : The full data set. 
C N      : INTEGER (1)
C        : Size of the data set
C NP     : INTEGER (1)
C        : The number of dimensions.
C NDIR   : INTEGER (1)
C        : The number of samples to draw.
C EPS    : DOUBLE PRECISION             
C        : Numerical Precision 
C     
C OUTPUT
C RDEP   : DOUBLE PRECISION (N1)
C        : The approximate regression depth of the planes	
C NSIN   : INTEGER (N1)
C        : Number of singular directions
C ERR    : INTEGER (N1)
C        : Indicating if all directions were singular (ERR=1) or not (ERR = 0)
	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Declaration of Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER MAXN,MAXP,MAXP1
      INTEGER N,NP,NNP,NDIR
      DOUBLE PRECISION X(N,NP+1),RDEP,T(NP+1),R(NP),XN(N)
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP)
      DOUBLE PRECISION AVE(NP),COV(NP,NP),EPS,D
      INTEGER RESID(N),JRES(N),JSAMP(N)
      INTEGER NDEP,NNEGTOT,NPOSTOT
      INTEGER NSIN, ERR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine
      MAXN=N
      MAXP=NP
      MAXP1=MAXP+1
C
C  Compute all residuals, and find the total number of positive and
C  negative residuals.
C
      NNEGTOT=0
      NPOSTOT=0
      DO 10 I=1,N
         D=X(I,NP+1)
         DO 15 J=1,NP
            D=D-T(J)*X(I,J)
 15      CONTINUE
         D=D-T(NP+1)
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

      CALL STAND_RDEPTH_APPR(MAXN,MAXP1,N,NP,X,XN,EPS)
      NNP=NP

      CALL RDEPTH_APPR_A(N,NP,NNP,NDIR,MAXN,MAXP1,X,R,RESID,JRES,XN,
     +          JSAMP,EPS,EVECS,EVALS,COV,AVE,NDEP,NSIN,NNEGTOT,NPOSTOT)

      IF (NSIN.EQ.NDIR) THEN
         ERR=1
      ELSEIF (NSIN.GT.(0.d0-EPS)) THEN
         ERR=0
      ELSE
         ERR=-1
      ENDIF
         
      RDEP=(NDEP+0.D0)/(N+0.D0)
      
      RETURN
 300  END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE STAND_RDEPTH_APPR(MAXN,MAXP1,N,NP,X,XN,EPS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C This routine standardises and centers the data by use of the 
C coordinate wise median and mad.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(MAXN,MAXP1),XN(N),EPS
      DOUBLE PRECISION QLOC,QSCA,AVE,VAR

      JN=0
      DO 10 J=1,NP
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
               NP=NP-1
               GOTO 10
            ELSE
               QSCA=DSQRT(VAR)
            ENDIF
         ENDIF
         JN=JN+1
         DO 60 I=1,N
            X(I,JN)=(X(I,J)-QLOC)/QSCA
 60      CONTINUE         
 10   CONTINUE

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE RDEPTH_APPR_A(N,NP,NNP,NDIR,MAXN,MAXP1,X,R,RESID,JRES,
     +     XN,JSAMP,EPS,EVECS,EVALS,COV,AVE,NDEP,NSIN,
     +     NNEGTOT,NPOSTOT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C To be called by RDEPTH_APPR1
C This subroutine performs the main part of the calculations.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(MAXN,MAXP1),R(NP),XN(N)
      DOUBLE PRECISION EPS,D
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP),COV(NP,NP),AVE(NP)
      INTEGER RESID(N),JRES(N),JSAMP(NP)
C
C  Initialize the number of singular samples.
C
      NSIN=0
C
C  Handle special case where N is equal to 1.
C
      IF (N.LE.1) THEN
         NDEP=0
         IF ((N.EQ.1).AND.(RESID(1).EQ.0)) NDEP=1
         RETURN
      ENDIF
C
C  Handle special case where NNP is equal to 1.
C
 25   IF (NNP.EQ.1) THEN
         CALL SORT_RDEPTH_APPR(X,RESID,N)
         NDEP=N
         NPOS=0
         NNEG=0
         DO 17 L=1,N
            IF (RESID(L).LE.0) NNEG=NNEG+1
            IF (RESID(L).GE.0) NPOS=NPOS+1
            IF (L.EQ.N) THEN
               D=1.D0
            ELSE 
               D=DABS(X(L,1)-X(L+1,1))
            ENDIF
            IF (D.GT.EPS) THEN
               NDEP=MIN0(NDEP,NPOS+NNEGTOT-NNEG)
               NDEP=MIN0(NDEP,NNEG+NPOSTOT-NPOS)
            ENDIF
 17      CONTINUE
         RETURN
      ENDIF
C
C  General case: call subroutine DEP.
C
      NP1=NNP+1
      CALL DEP_RDEPTH_APPR(N,NNP,NP1,NDIR,MAXN,MAXP1,
     +                     X,JSAMP,R,NNEGTOT,NPOSTOT,RESID,JRES,
     +                     XN,EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN)
C
C  If all points are identified as lying on the same hyperplane,
C  reduce the dimension of the data set by projection on that hyperplane,
C  and compute the depth on the reduced data set.
C
      IF (NSIN.EQ.(-1)) THEN
         NSIN=0
         NNP1=NNP
         NNP=NNP-1
         CALL REDUCE_RDEPTH_APPR(N,NNP,NNP1,
     +                            MAXN,MAXP1,X,R,EVECS,JSAMP,IERR)
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






      SUBROUTINE DEP_RDEPTH_APPR(N,NP,NP1,NDIR,MAXN,MAXP1,
     +                           X,JSAMP,R,NNEGTOT,NPOSTOT,RESID,JRES,
     +                           XN,EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C To be called by RDEPTH_APPR_A
C This subroutine performs the main part of the calculations.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(MAXN,MAXP1),XN(N)
      DOUBLE PRECISION EPS,R(NP),K,D,RAN(1)
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP),COV(NP,NP),AVE(NP)
      INTEGER RESID(N),JRES(N),JSAMP(NP)
C
C  Initialize regression depth and random seed.
C
      NDEP=N
      NRUN=0
      DO 100 NRAN=1,NDIR
C     
C  Draw a random sample of size np.
C     
         CALL UNIRAN(1,ISEED,RAN)
         I=N*RAN(1)+1.D0
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
C     
C  Compute the covariance matrix of the sample.
C     
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
C     
C  Compute the eigenvalues and corresponding eigenvectors 
C  of the covariance matrix.
C     
         CALL EIGEN(NP,NP,COV,EVALS,EVECS,R,AVE,IERR)
         IF (IERR.NE.0) THEN
            NSIN=NSIN+1
            GOTO 100
         ENDIF
         IF (EVALS(1).GT.EPS) THEN
            NSIN=NSIN+1
            GOTO 100
         ENDIF
C     
C  Test for singularity of the sample.
C     
         IF (EVALS(2).LE.EPS) THEN
            NSIN=NSIN+1
         ENDIF
C
C  Project all points on a line with direction given by 
C  the eigenvector of the smallest eigenvalue, i.e. the direction 
C  orthogonal on the hyperplane given by the np-subset.
C         
         KT=0.D0
         NT=0
         DO 90 J=1,NP
            IF (DABS(EVECS(J,1)).LE.EPS) NT=NT+1
 90      CONTINUE         
         IF (NT.EQ.NP) THEN
            NSIN=NSIN+1
            GOTO 100
         ENDIF

         NT=1
         DO 110 L=1,N
            K=0.D0
            DO 120 J=1,NP
               K=K+EVECS(J,1)*X(L,J)
 120        CONTINUE
            IF (L.EQ.1) THEN
               D=K
            ELSEIF (ABS(K-D).LE.EPS) THEN
               NT=NT+1
            ENDIF
            XN(L)=K
            JRES(L)=RESID(L)
 110     CONTINUE
C
C  If all projections collapse, return to reduce the dimension.
C
         IF (NT.EQ.N) THEN
            NSIN=-1
            RETURN
         ENDIF
C
C  Compute the one-dimensional regression depth of the fit,
C  and update the NP-dimensional depth.
C
         CALL SORT_RDEPTH_APPR(XN,JRES,N)
         NPOS=0
         NNEG=0
         DO 130 L=1,N
            IF (JRES(L).LE.0) NNEG=NNEG+1
            IF (JRES(L).GE.0) NPOS=NPOS+1
            IF (L.EQ.N) THEN
               D=1.D0
            ELSE 
               D=DABS(XN(L)-XN(L+1))
            ENDIF
            IF (D.GT.EPS) THEN
               NDEP=MIN0(NDEP,NPOS+NNEGTOT-NNEG)
               NDEP=MIN0(NDEP,NNEG+NPOSTOT-NPOS)
            ENDIF
 130     CONTINUE
 100  CONTINUE

      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE SORT_RDEPTH_APPR(B,RESID,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  Sorts an array B (of length N) in O(NlogN) time.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION B(N),AMM,XX
      INTEGER RESID(N)
      DIMENSION JLV(10000),JRV(10000)
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






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC






      SUBROUTINE REDUCE_RDEPTH_APPR(N,NNP,NNP1,MAXN,
     +  MAXP1,X,R,EVECS,W,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(MAXN,MAXP1),R(NNP1),S
      DOUBLE PRECISION EVECS(NNP1,NNP1)
      INTEGER W(NNP)
C
C  Invert matrix of base vectors EVECS.
C
      CALL VERT(EVECS,NNP+1,NNP+1,W,IERR)
      IF (IERR.LT.0) RETURN
C
C  Compute new NNP-dimensional coordinates for all points.
C      
      DO 30 IO=1,N
         DO 31 I=2,NNP+1
            R(I-1)=X(IO,1)*EVECS(I,1)
            DO 32 J=2,NNP+1
               R(I-1)=R(I-1)+X(IO,J)*EVECS(I,J)
 32         CONTINUE
 31      CONTINUE
         DO 33 I=1,NNP
            X(IO,I)=R(I)
 33      CONTINUE
 30   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
