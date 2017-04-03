      SUBROUTINE SORT(B,I1,I2,R,N,JLV,JRV)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  Sorts a double precision array B of length N and permutes two 
C  integer arrays I1 and I2 and one double precision array R in
C  the same way.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER N,I1(N),I2(N),H1,H2, JSS
      DOUBLE PRECISION B(N),XX,AMM
      DOUBLE PRECISION R(N),H3
      INTEGER JLV(N),JRV(N)
      INTEGER JR, JNC, J, JNDL, JTWE
      JSS=1
      JLV(1)=1
      JRV(1)=N
 10   JNDL=JLV(JSS)
      JR=JRV(JSS)
      JSS=JSS-1
 20   JNC=JNDL
      J=JR
      JTWE=(JNDL+JR)/2
      XX=B(JTWE)
 30   IF (B(JNC).GE.XX) GOTO 40
      JNC=JNC+1
      GOTO 30
 40   IF (XX.GE.B(J)) GOTO 50
      J=J-1
      GOTO 40
 50   IF (JNC.GT.J) GOTO 60
      AMM=B(JNC)
      H1=I1(JNC)
      H2=I2(JNC)
      H3=R(JNC) 
      B(JNC)=B(J)
      I1(JNC)=I1(J)
      I2(JNC)=I2(J)
      R(JNC)=R(J)
      B(J)=AMM
      I1(J)=H1
      I2(J)=H2
      R(J)=H3
      JNC=JNC+1
      J=J-1
 60   IF (JNC.LE.J) GOTO 30
      IF ((J-JNDL).LT.(JR-JNC)) GOTO 80
      IF (JNDL.GE.J) GOTO 70
      JSS=JSS+1
      JLV(JSS)=JNDL
      JRV(JSS)=J
 70   JNDL=JNC
      GOTO 100
 80   IF (JNC.GE.JR) GOTO 90
      JSS=JSS+1
      JLV(JSS)=JNC
      JRV(JSS)=JR
 90   JR=J
100   IF (JNDL.LT.JR) GOTO 20
      IF (JSS.NE.0) GOTO 10
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      function findq(aw,ncas,k)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C Finds the k-th order statistic of the array aw of length ncas.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer ncas,l,lr,k,jnc,j
      double precision findq
      double precision aw(ncas)
      double precision ax,wa
      l=1
      lr=ncas
 20   if(l.ge.lr) goto 90
      ax=aw(k)
      jnc=l
      j=lr
 30   if(jnc.gt.j) goto 80
 40   if(aw(jnc).ge.ax) goto 50
      jnc=jnc+1
      goto 40
 50   if(aw(j).le.ax) goto 60
      j=j-1
      goto 50
 60   if(jnc.gt.j) goto 70
      wa=aw(jnc)
      aw(jnc)=aw(j)
      aw(j)=wa
      jnc=jnc+1
      j=j-1
 70   goto 30
 80   if(j.lt.k) l=jnc
      if(k.lt.jnc) lr=j
      goto 20
 90   findq=aw(k)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,
     +           ALPHA,ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine Description:
C  Updates NCIRQ and NRANK, detects the special k-dividers and stores 
C  their angles and the constant terms of their equations.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      INTEGER N,M, MAXNUM
      INTEGER NCIRQ(N),NRANK(N),IND1(M),IND2(M)
      INTEGER KAND1(MAXNUM),KAND2(MAXNUM)
      INTEGER KOUNT,K,L,IV,IV1,IV2,D1,D2
      DOUBLE PRECISION X(N),Y(N),ANGLE(M),D(M)
      DOUBLE PRECISION ALPHA(MAXNUM),DUM,PI,PI2
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Actual Routine    
      D1=IND1(L)
      IV1=NRANK(D1)
      D2=IND2(L)
      IV2=NRANK(D2)
      IV=NCIRQ(IV1)
      NCIRQ(IV1)=NCIRQ(IV2)
      NCIRQ(IV2)=IV
      IV=IV1
      NRANK(D1)=IV2
      NRANK(D2)=IV
	 IF (((IV1.EQ.K).AND.(IV2.EQ.(K+1)))
     +      .OR.((IV2.EQ.K).AND.(IV1.EQ.(K+1)))
     +      .OR.((IV1.EQ.(N-K)).AND.(IV2.EQ.(N-K+1))) 
     +      .OR.((IV2.EQ.(N-K)).AND.(IV1.EQ.(N-K+1)))) THEN
	    IF (ANGLE(L).LT.PI2) THEN
	       DUM=ANGLE(L)+PI2
            ELSE
	       DUM=ANGLE(L)-PI2
            ENDIF
            IF (((IV1.EQ.K).AND.(IV2.EQ.(K+1)))
     +         .OR.((IV2.EQ.K).AND.(IV1.EQ.(K+1)))) THEN
	       IF (DUM.LE.PI2) THEN
		  ALPHA(KOUNT)=ANGLE(L)+PI
               ELSE
		  ALPHA(KOUNT)=ANGLE(L)
               ENDIF
            ENDIF
	    IF (((IV1.EQ.(N-K)).AND.(IV2.EQ.(N-K+1)))
     +        .OR.((IV2.EQ.(N-K)).AND.(IV1.EQ.(N-K+1)))) THEN
	       IF (DUM.LE.PI2) THEN
		  ALPHA(KOUNT)=ANGLE(L)
               ELSE
		  ALPHA(KOUNT)=ANGLE(L)+PI
               ENDIF
            ENDIF
	    KAND1(KOUNT)=IND1(L)
	    KAND2(KOUNT)=IND2(L)
	    D(KOUNT)=DSIN(ALPHA(KOUNT))*X(IND1(L))
     +                -DCOS(ALPHA(KOUNT))*Y(IND1(L))
            KOUNT=KOUNT+1
         ENDIF
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION NCEIL(M,J)
      integer nceil,M,J
      IF (MOD(M,J).EQ.0) THEN
         NCEIL=INT(REAL(M)/J)
      ELSE
         NCEIL=NINT(REAL(M)/J+0.5)
      ENDIF
      RETURN
      END
