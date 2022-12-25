C
      SUBROUTINE TRANFR
C
C FORMAL SOLVES THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION 'SOURCE'.
C 'ERROR' IS THE RESULTING ERROR IN THE DEFINITION OF THE CONTINUUM
C SCATTERING SOURCE FUNCTION. TRANSFR CALCULATES THE MATRIX ELEMENTS
C OF THE PROBLEM. INTENSITIES AT TAU=0 ARE RETURNED IN /CSURF/.
C 80.08.05 *NORD*
C
C Updated 04-jun-89/aake:  All mu-quadratures are over even functions
C of mu, and use TRQUA2.
C
      INCLUDE 'spectrum.inc'
C
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP),
     &  FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  NIMPAC,KIMPAC(NRAYS),PIMPAC(NRAYS),MMU(NDP),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP)
      COMMON /CSPHER/NCORE,DIFLOG,RADIUS,RR(NDP)
      COMMON /CSURF/HSURF,YSURF(NRAYS)
      COMMON /ROSSC/ROSS(NDP),CDUMM(NDP) /RHOC/RHO(NDP)
      COMMON /TRDBUG/IDEBUG
C
      DIMENSION FUN(NRAYS),DER(NRAYS*2),DMU(NRAYS+1)
      DIMENSION TAULOG(NDP),XMU1(NRAYS)
C
C MU LOOP
      DO 131 I=1,NIMPAC
      NTAU=KIMPAC(I)
      NTAU1=NTAU-1
C
C CALCULATE DTAUT ALONG THE RAY
      ZOLD=0.0
      DO 100 K=1,NTAU
      Z=SQRT(RR(NTAU-K+1)**2-PIMPAC(I)**2)
      XMU(I,NTAU-K+1)=-Z/RR(NTAU-K+1)
      IF (K.EQ.1) GO TO 100
      DZ=Z-ZOLD
      DZDR=DZ/(RR(NTAU-K+1)-RR(NTAU-K+2))
      DTAUT(NTAU-K+2)=DZDR*0.5*(X(NTAU-K+1)+S(NTAU-K+1)
     & +X(NTAU-K+2)+S(NTAU-K+2))*(TAU(NTAU-K+2)-TAU(NTAU-K+1))
100   ZOLD=Z
      TAUT(1)=DZDR*(X(1)+S(1))*TAU(1)
      DO 101 K=2,NTAU
101   TAUT(K)=TAUT(K-1)+DTAUT(K)
C
C SAVE THE TAU SCALE FOR THE RADIAL RAY (I=1).  THIS IS USED IN
C THE EXTRA/INTERPOLATION OF THE PFEAU FOR MU=0., FURTHER DOWN.
C
      IF (I.EQ.1) THEN
        DO 102 K=1,NTAU
          TAULOG(K)=ALOG(TAUT(K))
cccc          print*,k,tau(k),taut(k),taulog(k)
102     continue
      ENDIF
C
C K=1
      A=1./DTAUT(2)
      B=A**2
      SP2(1,I)=1.+2.*A
      SP3(1,I)=-2.*B
      EX(I)=TAUT(1)*(1.-0.5*TAUT(1)*(1.-TAUT(1)/3.)*(1.-taut(1)/4.))
      IF (TAUT(1).GT.0.01) then
        EX(I)=1.-EXP(-TAUT(1))
      endif
      SP2(1,I)=SP2(1,I)/(1.+2.*A*EX(I))
      SP3(1,I)=SP3(1,I)/(1.+2.*A*EX(I))
C
C K=2,NTAU-1
      DO 110 K=2,NTAU1
      DTAUC=0.5*(DTAUT(K)+DTAUT(K+1))
C...      AD(K,I)=0.1666667*DTAUT(K)/DTAUC
C...      BD(K,I)=0.1666667*DTAUT(K+1)/DTAUC
      AD(K,I)=0.
      BD(K,I)=0.
      SP1(K,I)=-1./(DTAUT(K)*DTAUC)+AD(K,I)
      SP2(K,I)=1.
110   SP3(K,I)=-1./(DTAUT(K+1)*DTAUC)+BD(K,I)
C
C K=NTAU
      AD(NTAU,I)=0.0
      BD(NTAU,I)=0.0
      SP1(NTAU,I)=0.0
      SP2(NTAU,I)=1.0
      SP3(NTAU,I)=0.0
      IF (I.LE.NCORE) GO TO 120
      AD(NTAU,I)=0.3333333
      SP1(NTAU,I)=0.3333333-2./DTAUT(NTAU)**2
120   CONTINUE
C
C ELIMINATE SUBDIAGONAL
      DO 130 K=1,NTAU1
      SP1(K,I)=-SP1(K+1,I)/(SP2(K,I)-SP3(K,I))
      SP2(K+1,I)=SP2(K+1,I)+SP1(K,I)*SP2(K,I)
130   SP2(K,I)=SP2(K,I)-SP3(K,I)
131   CONTINUE
C
C FIND A GOOD RAY FOR TRANSC.  THE XMU VALUES ARE INCREASING FROM
C -1.0 TOWARDS 0.0, AND WE THUS ALWAYS FIND A VALUE FOR ISCAT.
C
C      IMAX=MIN0(MMU(JTAU0),NCORE)
C
C NCORE MUST BE SMALLER THAN MMU(K) FOR ALL K, SO EFFECTIVELY THE
C MIN0(MMU(),NCORE) IS EQUAL TO NCORE.  SHOULD ONE REQUIRE ISCAT TO
C BE A CORE RAY?  PROBABLY NOT, SINCE WITH A LARGE TAUM (ALLOWED),
C THE ISCAT RAY IS FORCED TOWARDS SMALLER IMPACT PARAMETERS, WHICH
C DEGRADES PERFORMANCE.
C
      IMAX=MMU(JTAU0)
      TMP1=1.
      DO 132 I=1,IMAX
C
C      IF (XMU(I,JTAU0).LT.-0.577) ISCAT=I
C
C  Better to look for the ray which is closest to the Eddington angle
C
	TMP2=ABS(XMU(I,JTAU0)+0.577)
	IF (TMP2.LT.TMP1) THEN
	  TMP1=TMP2
	  ISCAT=I
	ENDIF
132   CONTINUE
      RETURN
C---------------------------------------------------------------------
C
      ENTRY FORMAL
C
C MU LOOP
      DO 170 I=1,NIMPAC
      NTAU=KIMPAC(I)
      NTAU1=NTAU-1
C
C INITIATE
      P(1)=SOURCE(1)
ccc      print*,'p(','1',') ',p(1)
      DO 140 K=2,NTAU
      P(K)=(1.-AD(K,I)-BD(K,I))*SOURCE(K)+AD(K,I)*SOURCE(K-1)+
     & BD(K,I)*SOURCE(min(K+1,ntau))
ccc      print*,'p(',k,') ',p(k)
140   continue
      IF(I.LE.NCORE) P(JTAU1)=SOURCE(JTAU1)+XMU(I,JTAU1)**2*XH(JTAU1)
ccc      print*,'p(jtau1) ',jtau1,p(jtau1)
C
C ACCUMULATE RIGHT HAND SIDE
      DO 150 K=1,NTAU1
150   P(K+1)=P(K+1)+SP1(K,I)*P(K)
ccc        print*,'accumulated p'
ccc      do k=1,ntau
ccc        print*,k,p(k)
ccc      enddo
C
C BACKSUBSTITUTE
      PFEAU(I,NTAU)=P(NTAU)/SP2(NTAU,I)
ccc       print*,i,ntau,pfeAU(i,ntau),p(ntau),sp2(ntau,i)
      DO 160 K=1,NTAU1
      PFEAU(I,NTAU-K)=(P(NTAU-K)-
     & SP3(NTAU-K,I)*PFEAU(I,NTAU-K+1))/SP2(NTAU-K,I)
       if (k.eq.1) then
ccc       print*,i,ntau-k,pfeAU(i,ntau-k)
       endif
      IF (PFEAU(I,NTAU-K).LE.0.0) GO TO 230
160   CONTINUE
C
C END MU LOOP
      YSURF(I)=2.*(1.-EX(I))*PFEAU(I,1)+EX(I)**2*SOURCE(1)
      FUN(I)=-XMU(I,1)*(PFEAU(I,1)-SOURCE(1)*EX(I))
      IF (YSURF(I).LE.0.0) GO TO 231
170   CONTINUE
C
C  INTERPOLATE TO PFEAU AT MU=0 FOR THOSE K THAT HAVE NO RAY
*  THIS IS THE ORIGINAL CODE, WHICH EXTRAPOLATES IN MU, FOR CONSTANT K.
*  THIS WAS FOUND TO PRODUCE BAD RESULTS FOR THE FLUX, WHICH IS A
*  DERIVATIVE OF XK, WHICH IS THE SECOND MOMENT OF THE FUNCTION
*  BEING EXTRAPOLATED TO ZERO.
*      DO 181 K=1,JTAU1
*      II=MMU(K)
*      IF (KIMPAC(II).EQ.K) GO TO 181
*      PX=-XMU(II-2,K)/(XMU(II-1,K)-XMU(II-2,K))
*      QX=1.-PX
*      PFEAU(II,K)=EXP(ALOG(PFEAU(II-2,K))*QX+ALOG(PFEAU(II-1,K))*PX)
*181   CONTINUE
*
*  INTERPOLATE/EXTRAPOLATE IN K FOR PFEAU AT MU=0
*  THIS SHOWS WHAT IS INTENDED TO BE DONE, IF USING AN EXTERNAL
*  INTERPOLATION ROUTINE.
*       DO 181 I=1,NIMPAC
*         K=KIMPAC(I)
*         F(I)=PFEAU(I,K)
*         X(I)=ALOG(TAU(K))
*181   CONTINUE
*      DO 182 K=1,JTAU1
*        XX(K)=ALOG(TAU(K))
*182   CONTINUE
*      CALL INTERP(NIMPAC,X,F,JTAU1,XX,FF)
*      DO 183 K=1,JTAU1
*        I=MMU(K)
*        PFEAU(I,K)=FF(K)
*183   CONITNUE
C
C  DO IT INLINE:  EXTRAPOLATE/INTERPOLATE THE PFEAU IN ALOG(TAU), FOR MU=0.
C   START WITH 3 OUTERMOST POINTS (YOU PROBABLY NOTICED THAT KIMPAC(NIMPAC)=4)
C   JUMP OVER THE RAYS (WHERE PFEAU(MU=0.) IS ALREADY CALCULATED), AND
C   INTERPOLATE THE REST OF THE TIME.
C
      IF (NCORE.EQ.NIMPAC) THEN
C
C  TEMPORARY SECURITY TO AVOID DIVIDE BY ZERO ERRORS IN COMPUTATION OF PX,
C  WHEN THE OPACITY IS SO LARGE THAT THERE ARE ONLY NCORE RAYS.
C  IN THAT CASE, ONE USES THE OLD INTERPOLATION ROUTINE
C
      print*,'NCORE = NIMPAC. INTERPOLATION IN MU USED'
      DO 1810 K=1,JTAU1
        II=MMU(K)
        IF (KIMPAC(II).EQ.K) GO TO 1810
        PX=-XMU(II-2,K)/(XMU(II-1,K)-XMU(II-2,K))
        QX=1.-PX
        PFEAU(II,K)=EXP(ALOG(PFEAU(II-2,K))*QX+ALOG(PFEAU(II-1,K))*PX)
      if(pfeau(ii,k).lt.0) 
     &       print*, 'k',k,ii,pfeau(ii,k),' <0 pfeau'
1810   CONTINUE

      ELSE
C
C HERE IS THE NORMAL ONE.  EXTRAPOLATE FOR K<4, THEN INTERPOLATE
C
      I=NIMPAC
      DO 181 K=1,JTAU1
        IF (K.NE.KIMPAC(I)) THEN
          PX=(TAULOG(K)-TAULOG(KIMPAC(I)))/
     &       (TAULOG(KIMPAC(I-1))-TAULOG(KIMPAC(I)))
          QX=1.-PX
*
ccc      print*,ALOG(PFEAU(MMU(KIMPAC(I)),KIMPAC(I)))*QX+
ccc     &     ALOG(PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)))*PX 
ccc          if ( ALOG(PFEAU(MMU(KIMPAC(I)),KIMPAC(I)))*QX+
ccc     &     ALOG(PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)))*PX .gt.-20.) then
cc       print*,mmu(kimpac(i)),kimpac(i),PFEAU(MMU(KIMPAC(I)),KIMPAC(I)),
cc     &        mmu(kimpac(i-1)),kimpac(i-1),
cc     &        PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1))
*
           PFEAU(MMU(K),K)=sngl(dEXP(
     &       dLOG(dble(PFEAU(MMU(KIMPAC(I)),KIMPAC(I))))*QX+
     &       dLOG(dble(PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1))))*PX))
ccc          else
ccc           pfeau(mmu(k),k)=1.e-15
ccc          endif
*
        ENDIF
        IF (K.EQ.(KIMPAC(I-1)-1)) THEN
          I=I-1
        ENDIF
181   CONTINUE

      ENDIF

C
C CALCULATE MEAN INTENSITY USING CUBIC SPLINE WITH CENTRAL DIFFERENCE
      DO 190 K=1,JTAU1
        XJ(K)=TRQUA2(MMU(K),XMU(1,K),PFEAU(1,K),DMU,DER)
        ERROR(K)=(XJ(K)*S(K)+BPLAN(K)*X(K))/(X(K)+S(K))-SOURCE(K)
190   CONTINUE
      RETURN
C--------------------------------------------------------------------
C
      ENTRY TRMOM
C
C FLUX AT TAU(1)
      XH(1)=TRQUA2(MMU(1),XMU,FUN,DMU,DER)
C
C  CALCULATE SECOND MOMENT XK
C
      DO 201 K=1,JTAU1
      NMU=MMU(K)
      DO 200 I=1,NMU
200     FUN(I)=PFEAU(I,K)*XMU(I,K)**2
201   XK(K)=TRQUA2(NMU,XMU(1,K),FUN,DMU,DER)
C
C CALCULATE FIRST MOMENT, XH, FROM MOMENT RELATION
      DO 211 K=JTAU1+1,JTAU
211   XH(K)=( XK(K)-XK(K-1)+
     &       (XJ(K)+XJ(K-1)-3.*(XK(K)+XK(K-1)))*
     &       (RR(K-1)-RR(K))/(RR(K)+RR(K-1))    )*2.0/
     &          ( (TAU(K)-TAU(K-1)) * (X(K)+S(K)+X(K-1)+S(K-1)) )
C
C CALCULATE FIRST MOMENT BY USING R=DP/DTAU.  THIS IS MORE ACCURATE
C IN THE OPTICALLY THIN PARTS.
      ZOLD=0.0
      DO 212 K=2,JTAU1
        NMU=MMU(K)
        DO 213 I=1,NMU-1
          DZDR=(XMU(I,K)*RR(K)-XMU(I,K-1)*RR(K-1))/(RR(K-1)-RR(K))
          DTAU=DZDR*0.5*(X(K-1)+S(K-1)+X(K)+S(K))*(TAU(K)-TAU(K-1))
          XMU1(I)=
     &      -(XMU(I,K)*RR(K)+XMU(I,K-1)*RR(K-1))/(RR(K)+RR(K-1))
          FUN(I)=XMU1(I)*(PFEAU(I,K)-PFEAU(I,K-1))/DTAU
213     CONTINUE
        FUN(NMU)=0.
        XMU1(NMU)=0.
        XH(K)=-TRQUA2(NMU,XMU1,FUN,DMU,DER)
212   CONTINUE
C
C SURFACE FLUX
      NMU=MMU(1)
      PX=-XMU(NMU-2,1)/(XMU(NMU-1,1)-XMU(NMU-2,1))
      QX=1.-PX
      YSURF(NMU)=EXP(ALOG(YSURF(NMU-2))*QX+ALOG(YSURF(NMU-1))*PX)
      DO 220 I=1,NMU
220   FUN(I)=-XMU(I,1)*YSURF(I)
      HSURF=0.5*TRQUA2(NMU,XMU,FUN,DMU,DER)

cc      Print*,'check surf flux:',hsurf/xh(1)
      RETURN
C------------------------------------------------------------------
C
C EMERGENCY EXIT
230   KK=NTAU-K
      IDEBUG=2
      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
232   FORMAT('0NON-POSITIVE RESULT AT I,K,J0,J1,N=',5I3)
      GO TO 233
231   KK=0
      IDEBUG=3
      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
233   PRINT 237,NCORE,ISCAT,DIFLOG,RADIUS,EX(I)
237   FORMAT(' NCORE,ISCAT,DIFLOG,RADIUS,EX(I)=',2I3,3G12.3)
CCCCC      WRITE (13,*) JTAU,TAU,X,S,BPLAN,RADIUS,RR,RHO,ROSS
!!      PRINT 236,TAU,X,S,BPLAN,RR,RHO,ROSS,SOURCE
!!     & ,(YSURF(I),(PFEAU(I,K),K=1,39),I=1,NMU)
236   FORMAT('0TAU=',10(/10E12.4)/'0X=',10(/10E12.4)/'0S=',10(/10E12.4)
     & /'0BPLAN=',10(/10E12.4)/'0RR=',10(/10E12.4)/'0RHO=',10(/10E12.4)
     & /'0ROSS=',10(/10E12.4)/'0SOURCE=',10(/10E12.4)
     & /'0YSURF,PFEAU='/(10E12.4))
CCCCC      PRINT 238
238   FORMAT(' DEBUG INFORMATION WRITTEN ON UNIT 13')
      RETURN
      END
