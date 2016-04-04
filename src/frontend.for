      SUBROUTINE POTENTIAL(Z, ZS, A)
C     Get Z, ZS and A and create the grid and spline interpolation
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*12 FILEN
      PARAMETER (NDIM=10000)
      PARAMETER (SL=137.03599976D0,PI=3.1415926535897932D0)
      DIMENSION R0(NDIM),RV0(NDIM)
C  ****  Coulomb wave function parameters.
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C
C  ****  Read potential parameters.
C
   10 CONTINUE
      ALPHA=ABS(ALPHA)
C
C  ****  Potential radial grid.
C
      RNN=35.0D0/MAX(ALPHA,1.0D0)
      NV=600
      CALL SGRID(R0,RV0,RNN,1.0D-6,1.0D0,NV,NDIM)
C
      DO I=1,NV
        RV0(I)=Z+ZS*EXP(-ALPHA*R0(I))
      ENDDO

      CALL ERRSPL(ERR,R0,RV0,NV)
      WRITE(6,1001) ERR*100.0D0
 1001 FORMAT(1X,'*** Estimated error of the interpolating spline =',
     1  1P,E9.2,' %.')
      IF(ERR.GT.1.0D-1) WRITE(6,1002)
 1002 FORMAT(/2X,'**** The interpolation error seems to be too',
     1  ' large.'/7X,'The function R*V(R) should be tabulated',
     2  ' in a denser grid.')

      CALL VINT(R0,RV0,NV)
      END
      
      SUBROUTINE SCHROED_BOUND(N, L, EPS, E)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INTEGER, INTENT(IN) :: N, L
      DOUBLE PRECISION , INTENT(IN) :: EPS
      DOUBLE PRECISION, INTENT(OUT) :: E
      DOUBLE PRECISION TRUE_EPS
      PARAMETER (NDIM=10000) 
      PARAMETER (SL=137.03599976D0,PI=3.1415926535897932D0)
      DIMENSION R0(NDIM),RV0(NDIM)
C  ****  Coulomb wave function parameters.
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER

      TRUE_EPS=MAX(EPS,1.0D-15)
      E=-Z**2/(2.0D0*N*N)
      NGP=1000
      RN=1000.00
      CALL SGRID(RAD,RV0,RN,1.0D-6,5.0D0,NGP,NDIM)
      CALL SBOUND(E,TRUE_EPS,TRUE_EPS,N,L)
      IF(IER.NE.0) WRITE(6, *) 'There has been an error!'
      END

      SUBROUTINE SCHROED_FREE(E, L, EPS, IPS, CPS, META)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INTEGER, INTENT(IN) :: L
      DOUBLE PRECISION , INTENT(IN) :: E, EPS
      DOUBLE PRECISION, INTENT(OUT) :: IPS, CPS, META
      DOUBLE PRECISION TRUE_EPS
      PARAMETER (NDIM=10000) 
      PARAMETER (SL=137.03599976D0,PI=3.1415926535897932D0)
      DIMENSION R0(NDIM),RV0(NDIM)
C  ****  Coulomb wave function parameters.
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
    
      TRUE_EPS=MAX(EPS,1.0D-15)
      NGP=1000
      WAVEL=2.0D0*PI/SQRT(E+E)
      RN=0.025D0*WAVEL*DBLE(NGP)
      DRN=WAVEL/20.0D0
      CALL SGRID(RAD,RV0,RN,1.0D-6,DRN,NGP,NDIM)
      CALL SFREE(E,EPS,PHASE,L)
      IF(IER.NE.0) WRITE(6, *) 'There has been an error'
      IPS = PHASE
      CPS = DELTA
      META = ETA
      END

C

      SUBROUTINE DIRAC_BOUND(N, K, EPS, E)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INTEGER, INTENT(IN) :: N, K
      DOUBLE PRECISION , INTENT(IN) :: EPS
      DOUBLE PRECISION, INTENT(OUT) :: E
      DOUBLE PRECISION TRUE_EPS
      PARAMETER (NDIM=10000) 
      PARAMETER (SL=137.03599976D0,PI=3.1415926535897932D0)
      DIMENSION R0(NDIM),RV0(NDIM)
C  ****  Coulomb wave function parameters.
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      TRUE_EPS=MAX(EPS,1.0D-15)
      E=-Z**2/(2.0D0*N*N)
      NGP=1000
      RN=1000.00
      CALL SGRID(RAD,RV0,RN,1.0D-6,5.0D0,NGP,NDIM)
      CALL DBOUND(E,TRUE_EPS,TRUE_EPS,N,K)
      IF(IER.NE.0) WRITE(6, *) 'There has been an error!'
      END


      SUBROUTINE DIRAC_FREE(E, L, EPS, IPS, CPS, META)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INTEGER, INTENT(IN) :: L
      DOUBLE PRECISION , INTENT(IN) :: E, EPS
      DOUBLE PRECISION, INTENT(OUT) :: IPS, CPS, META
      DOUBLE PRECISION TRUE_EPS
      PARAMETER (NDIM=10000) 
      PARAMETER (SL=137.03599976D0,PI=3.1415926535897932D0)
      DIMENSION R0(NDIM),RV0(NDIM)
C  ****  Coulomb wave function parameters.
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
    
      TRUE_EPS=MAX(EPS,1.0D-15)
      NGP=1000
      WAVEL=2.0D0*PI/SQRT(E+E)
      RN=0.025D0*WAVEL*DBLE(NGP)
      DRN=WAVEL/20.0D0
      CALL SGRID(RAD,RV0,RN,1.0D-6,DRN,NGP,NDIM)
      CALL DFREE(E,EPS,PHASE,L)
      IF(IER.NE.0) WRITE(6, *) 'There has been an error'
      IPS = PHASE
      CPS = DELTA
      META = ETA
      END

      
      
C   PROGRAM DEMORAD (demo for subroutine package RADIALF)
C
C     This program solves the radial wave equation for partially
C  screened Coulomb fields of the form R*V(R)=Z+ZS*EXP(-A*R).
C
C                                           Francesc Salvat.  May, 2000.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*12 FILEN
      PARAMETER (NDIM=10000)
      PARAMETER (SL=137.03599976D0,PI=3.1415926535897932D0)
      DIMENSION R0(NDIM),RV0(NDIM)
C  ****  Coulomb wave function parameters.
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C
C  ****  Read potential parameters.
C
   10 CONTINUE
      WRITE(6,*) ' Potential function: R*V(R)=Z+ZS*EXP(-A*R)'
      WRITE(6,*) ' Enter Z, ZS and A ...'
      READ(5,*) Z,ZS,ALPHA
      ALPHA=ABS(ALPHA)
C
C  ****  Potential radial grid.
C
      RNN=35.0D0/MAX(ALPHA,1.0D0)
      NV=600
      CALL SGRID(R0,RV0,RNN,1.0D-6,1.0D0,NV,NDIM)
C
      DO I=1,NV
        RV0(I)=Z+ZS*EXP(-ALPHA*R0(I))
      ENDDO

C-----------------------------------------------------------------------
C  ----  Spline interpolation test.
      CALL ERRSPL(ERR,R0,RV0,NV)
      WRITE(6,1001) ERR*100.0D0
 1001 FORMAT(1X,'*** Estimated error of the interpolating spline =',
     1  1P,E9.2,' %.')
      IF(ERR.GT.1.0D-1) WRITE(6,1002)
 1002 FORMAT(/2X,'**** The interpolation error seems to be too',
     1  ' large.'/7X,'The function R*V(R) should be tabulated',
     2  ' in a denser grid.')
C-----------------------------------------------------------------------

      CALL VINT(R0,RV0,NV)
C
c      WRITE(6,1003)
c 1003 FORMAT(/2X,'The user grid can be read from a file (single',
c     1  ' column, increasing radii)'/2X,'or determined automati',
c     2  'cally. If you wish to use a prepared grid, enter'/2X,
c     1  'the file name, otherwise type ''n'' ...')
c      READ(5,'(A12)') FILEN


c      IF(FILEN.NE.'N'.AND.FILEN.NE.'n') THEN
C        IGRID=1
C        OPEN(7,FILE=FILEN)
C        DO I=1,NDIM
C          READ(7,*,END=20) RAD(I)
C          NGP=I
C        ENDDO
C        CLOSE(UNIT=7)
c     ELSE
c       IGRID=0
c     ENDIF
C
      IGRID=0

      OPEN(8,FILE='res.res')
   20 CONTINUE
      WRITE(6,*) '  '
      WRITE(6,*) ' Select one option ...'
      WRITE(6,*) '   1: Schrodinger bound state,   2: Sch',
     1  'rodinger free state,'
      WRITE(6,*) '   3: Dirac bound state,         4: Dir',
     2  'ac free state.'
      READ(5,*) IOPT
      WRITE(6,*) '  '
C
C  ************  Schrodinger equation. Bound state.
C
      IF(IOPT.EQ.1) THEN
        WRITE(6,*) ' Enter N, L and EPS ...'
        READ(5,*) N,L,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(N.LT.1.OR.L.GE.N) GO TO 20
        IF(IGRID.EQ.0) THEN
          NGP=1000
          RN=1000.0D0
          CALL SGRID(RAD,RV0,RN,1.0D-6,5.0D0,NGP,NDIM)
        ENDIF
        E=-Z**2/(2.0D0*N*N)
        CALL SBOUND(E,EPS,EPS,N,L)
        IF(IER.NE.0) GO TO 20
        WRITE(6,1100) Z,ZS,ALPHA,N,L,EPS,E
        WRITE(8,1100) Z,ZS,ALPHA,N,L,EPS,E
 1100 FORMAT(/1X,1P,'****  Schrodinger eq. ',
     1  'Potential function: R*V(R)=Z+ZS*EXP(-A*R)'
     2  /7X,'Z=',E13.6,', ZS=',E13.6,', A=',E13.6
     3  /7X,'Bound state: N=',I4,', L=',I4,'   (EPS=',E8.1,')'
     4  /7X,'Binding energy=',E22.15)
C
C  ************  Schrodinger equation. Free state.
C
      ELSE IF(IOPT.EQ.2) THEN
        WRITE(6,*) ' Enter E, L and EPS ...'
        READ(5,*) E,L,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(E.LT.0.0D0.OR.L.LT.0) GO TO 20
        IF(IGRID.EQ.0) THEN
          NGP=1000
          WAVEL=2.0D0*PI/SQRT(E+E)
          RN=0.025D0*WAVEL*DBLE(NGP)
          DRN=WAVEL/20.0D0
          CALL SGRID(RAD,RV0,RN,1.0D-6,DRN,NGP,NDIM)
        ENDIF
        CALL SFREE(E,EPS,PHASE,L)
        IF(IER.NE.0) GO TO 20
        WRITE(6,1200) Z,ZS,ALPHA,E,L,EPS,PHASE,DELTA,ETA
        WRITE(8,1200) Z,ZS,ALPHA,E,L,EPS,PHASE,DELTA,ETA
 1200 FORMAT(/1X,1P,'****  Schrodinger eq. ',
     1  'Potential function: R*V(R)=Z+ZS*EXP(-A*R)'
     2  /7X,'Z=',E13.6,', ZS=',E13.6,', A=',E13.6
     3  /7X,'Free state: E=',E13.6,', L=',I4,'   (EPS=',E8.1,')'
     4  /7X,'  Inner phase shift=',E22.15,
     5  /7X,'Coulomb phase shift=',E22.15,'   (ETA=',E13.6,')')
        WRITE(6,1201) ILAST,RAD(ILAST)
 1201 FORMAT(7X,'Matching radius: RAD(',I4,')=',1P,E22.15)
C
C  ************  Dirac equation. Bound state.
C
      ELSE IF(IOPT.EQ.3) THEN
        WRITE(6,*) ' Enter N, K and EPS ...'
        READ(5,*) N,K,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(N.LT.1.OR.K.EQ.0.OR.K.GE.N.OR.K.LT.-N) GO TO 20
        IF(IGRID.EQ.0) THEN
          NGP=1000
          RN=1000.0D0
          CALL SGRID(RAD,RV0,RN,1.0D-6,5.0D0,NGP,NDIM)
        ENDIF
        E=-Z**2/(2.0D0*N*N)
        CALL DBOUND(E,EPS,EPS,N,K)
        IF(IER.NE.0) GO TO 20
        WRITE(6,1300) Z,ZS,ALPHA,N,K,EPS,E
        WRITE(8,1300) Z,ZS,ALPHA,N,K,EPS,E
 1300 FORMAT(/1X,1P,'****  Dirac equation. ',
     1  'Potential function: R*V(R)=Z+ZS*EXP(-A*R)'
     2  /7X,'Z=',E13.6,', ZS=',E13.6,', A=',E13.6
     3  /7X,'Bound state: N=',I4,', K=',I4,'   (EPS=',E8.1,')'
     4  /7X,'Binding energy=',E22.15)
C
C  ************  Dirac equation. Free state.
C
      ELSE IF(IOPT.EQ.4) THEN
        WRITE(6,*) ' Enter E, K and EPS ...'
        READ(5,*) E,K,EPS
        EPS=MAX(EPS,1.0D-15)
        IF(E.LT.0.0D0.OR.K.EQ.0) GO TO 20
        IF(K.LT.0) THEN
          L=-K-1
        ELSE
          L=K
        ENDIF
        IF(IGRID.EQ.0) THEN
          NGP=1000
          WAVEL=2.0D0*PI/SQRT(E*(2.0D0+E/SL**2))
          RN=0.025D0*WAVEL*DBLE(NGP)
          DRN=WAVEL/20.0D0
          CALL SGRID(RAD,RV0,RN,1.0D-6,DRN,NGP,NDIM)
          IF(IER.NE.0) STOP 'Error in the grid definition (DF).'
        ENDIF
        CALL DFREE(E,EPS,PHASE,K)
        IF(IER.NE.0) GO TO 20
        WRITE(6,1400) Z,ZS,ALPHA,E,K,EPS,PHASE,DELTA,ETA
        WRITE(8,1400) Z,ZS,ALPHA,E,K,EPS,PHASE,DELTA,ETA
 1400 FORMAT(/1X,1P,'****  Dirac equation. ',
     1  'Potential function: R*V(R)=Z+ZS*EXP(-A*R)'
     2  /7X,'Z=',E13.6,', ZS=',E13.6,', A=',E13.6
     3  /7X,'Free state: E=',E13.6,', K=',I4,'   (EPS=',E8.1,')'
     4  /7X,'  Inner phase shift=',E22.15,
     5  /7X,'Coulomb phase shift=',E22.15,'   (ETA=',E13.6,')')
        WRITE(6,1401) ILAST,RAD(ILAST)
 1401 FORMAT(7X,'Matching radius: RAD(',I4,')=',1P,E22.15)
      ELSE
        GO TO 10
      ENDIF
C
C  ****  Radial wave functions are written in file 'waves.dat'.
C
      OPEN(7,FILE='waves.dat')
      WRITE(7,1500)
 1500 FORMAT(1X,'#  Radial wave functions calculated by RADIAL.',
     1  /1X,'#',7X,'R',14X,'P(R)',12X,'Q(R)')
      DO I=1,NGP
        WRITE(7,'(1X,1P,3E16.8)') RAD(I),P(I),Q(I)
      ENDDO
      CLOSE(UNIT=7)
      GO TO 20
      END
