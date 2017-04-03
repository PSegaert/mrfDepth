CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE VERT(V,LV,N,W,IERR)
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module VERT from package NAPACK.
* Retrieved from NETLIB on Wed Feb 19 03:31:44 1997.
* ======================================================================
C
C      ________________________________________________________
C     |                                                        |
C     |                INVERT A GENERAL MATRIX                 |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         V     --ARRAY CONTAINING MATRIX                |
C     |                                                        |
C     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
C     |                                                        |
C     |         W     --INTEGER WORK ARRAY WITH AT LEAST N-1   |
C     |                      ELEMENTS                          |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         V     --INVERSE                                |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |________________________________________________________|
C
      DOUBLE PRECISION V(LV,1),S,T
      INTEGER W(1),I,J,K,L,M,N,P
      IF ( N .EQ. 1 ) GOTO 110
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
C     ---------------------------------------
C     |*** FIND PIVOT AND START ROW SWAP ***|
C     ---------------------------------------
      P = L
      IF ( M .GT. N ) GOTO 30
      S = DABS(V(L,L))
      DO 20 I = M,N
           T = DABS(V(I,L))
           IF ( T .LE. S ) GOTO 20
           P = I
           S = T
20    CONTINUE
      W(L) = P
30    S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0. ) GOTO 120
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      V(L,L) = -1.
      S = 1./S
      DO 40 I = 1,N
40         V(I,L) = -S*V(I,L)
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0. ) GOTO 50
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      IF ( K .EQ. 0 ) GOTO 70
      DO 60 I = 1,K
60         V(I,J) = V(I,J) + T*V(I,L)
70    V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO 80 I = M,N
80         V(I,J) = V(I,J) + T*V(I,L)
      GOTO 50
C     -----------------------
C     |*** PIVOT COLUMNS ***|
C     -----------------------
90    L = W(K)
      DO 100 I = 1,N
           T = V(I,L)
           V(I,L) = V(I,K)
100        V(I,K) = T
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
      RETURN
110   IF ( V(1,1) .EQ. 0. ) GOTO 120
      V(1,1) = 1./V(1,1)
      RETURN
120   IERR=-1
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine eigen(nm,n,a,w,z,fv1,fv2,ierr)
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module RS from package EISPACK.
* Retrieved from NETLIB on Wed Nov 27 07:41:24 1996.
* ======================================================================
c
      implicit double precision (a-h,o-z)
      integer n,nm,ierr
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)

c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tql2.
c           the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end

      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE REDUCE(N,NNP,NNP1,MAXN,MAXP,X,T,R,EVECS,W,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(MAXN,MAXP),T(NNP1),R(NNP1)
      DOUBLE PRECISION EVECS(NNP1,NNP1)
      INTEGER W(NNP)

      IERR=0
C
C  Invert matrix of base vectors EVECS.
C
      CALL VERT(EVECS,NNP+1,NNP+1,W,IERR)
      IF (IERR.LT.0) RETURN
C
C  Compute new NNP-dimensional coordinates for all points and theta.
C      
      DO 30 I=2,NNP+1
         R(I-1)=T(1)*EVECS(I,1)
         DO 31 J=2,NNP+1
            R(I-1)=R(I-1)+T(J)*EVECS(I,J)
 31      CONTINUE
 30   CONTINUE
      DO 32 I=1,NNP
         T(I)=R(I)
 32   CONTINUE

      DO 40 IO=1,N
         DO 41 I=2,NNP+1
            R(I-1)=X(IO,1)*EVECS(I,1)
            DO 42 J=2,NNP+1
               R(I-1)=R(I-1)+X(IO,J)*EVECS(I,J)
 42         CONTINUE
 41      CONTINUE
         DO 43 I=1,NNP
            X(IO,I)=R(I)
 43      CONTINUE
 40   CONTINUE
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC    ------------------------------------------------------------------
CC     *EQUAT* : SOLVES A SYSTEM OF LINEAR EQUATIONS.
CC    ------------------------------------------------------------------
      SUBROUTINE EQUAT(AM,MAXP,MAXP1,HVEC,MAXPP1,NA,NB,NERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AM(MAXP,MAXP1)
      DOUBLE PRECISION HVEC(MAXPP1),TURN,SWAP,DETER
      JDM=MAXP
      DETER=1.0D0
      N=NA
      JMAT=N+NB
      JNK=0
      DO 10 J=1,JMAT
      JNK=(J-1)*MAXP
      DO 10 NC=1,MAXP
      JNK=JNK+1
      HVEC(JNK)=AM(NC,J)
 10   CONTINUE
      NZNDE=N-1
      LCLPL=-JDM
      DO 120 JHFD=1,N
      TURN=0.D0
      LCLPL=LCLPL+JDM+1
      JDEL=LCLPL+N-JHFD
      DO 40 JNCB=LCLPL,JDEL
      IF(DABS(HVEC(JNCB))-DABS(TURN)) 40,40,30
 30   TURN=HVEC(JNCB)
      LDEL=JNCB
 40   CONTINUE
      IF(DABS(TURN).LE.1D-8) GOTO 170
 50   IF(LDEL-LCLPL) 60,80,60
 60   DETER=-DETER
      LDEL=LDEL-JDM
      JNCB=LCLPL-JDM
      DO 70 JNCC=JHFD,JMAT
      LDEL=LDEL+JDM
      JNCB=JNCB+JDM
      SWAP=HVEC(JNCB)
      HVEC(JNCB)=HVEC(LDEL)
 70   HVEC(LDEL)=SWAP
 80   DETER=DETER*TURN
      IF (JHFD.EQ.N)  GOTO 120
      TURN=1./TURN
      JNCB=LCLPL+1
      DO 90 JNCC=JNCB,JDEL
 90   HVEC(JNCC)=HVEC(JNCC)*TURN
      JNCD=LCLPL
      JROW=JHFD+1
      DO 110 JNCB=JROW,N
      JNCD=JNCD+1
      JNCE=LCLPL
      JNCF=JNCD
      DO 100 JNCC=JROW,JMAT
      JNCE=JNCE+JDM
      JNCF=JNCF+JDM
 100  HVEC(JNCF)=HVEC(JNCF)-HVEC(JNCE)*HVEC(JNCD)
 110  CONTINUE
 120  CONTINUE
      NERR=0
      NEQA=N+1
      JBEGX=NZNDE*JDM+1
      DO 150 JNC=NEQA,JMAT
      JBEGX=JBEGX+JDM
      JENDX=JBEGX+N
      JBEGC=N*JDM+1
      JENDC=JBEGC+NZNDE
      DO 140 JNCB=1,NZNDE
      JENDX=JENDX-1
      JBEGC=JBEGC-JDM
      JENDC=JENDC-JDM-1
      HVEC(JENDX)=HVEC(JENDX)/HVEC(JENDC+1)
      SWAP=HVEC(JENDX)
      JNCD=JBEGX-1
      DO 130 JNCC=JBEGC,JENDC
      JNCD=JNCD+1
 130  HVEC(JNCD)=HVEC(JNCD)-HVEC(JNCC)*SWAP
 140  CONTINUE
      HVEC(JBEGX)=HVEC(JBEGX)/HVEC(1)
 150  CONTINUE
      JNC=-JDM
      JBEGX=NZNDE*JDM+1
      JENDX=JBEGX+NZNDE
      DO 160 JNCB=NEQA,JMAT
      JBEGX=JBEGX+JDM
      JENDX=JENDX+JDM
      JNC=JNC+JDM
      JNCD=JNC
      DO 160 JNCC=JBEGX,JENDX
      JNCD=JNCD+1
 160  HVEC(JNCD)=HVEC(JNCC)
      GOTO 180
 170  NERR=-1
 180  JNK=0
      DO 190 J=1,JMAT
      DO 190 NC=1,MAXP
      JNK=JNK+1
      AM(NC,J)=HVEC(JNK)
 190  CONTINUE
      RETURN
      END
