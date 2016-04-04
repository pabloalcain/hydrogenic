CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C         RRRRR     AA    DDDDD   IIII    AA    L       FFFFFF         C
C         R    R   A  A   D    D   II    A  A   L       F              C
C         R    R  A    A  D    D   II   A    A  L       F              C
C         RRRRR   AAAAAA  D    D   II   AAAAAA  L       FFFF           C
C         R  R    A    A  D    D   II   A    A  L       F              C
C         R   R   A    A  DDDDB   IIII  A    A  LLLLLL  F              C
C                                                                      C
C                                                    (version 2005).   C
C                                                                      C
C  Numerical solution of the Schrodinger (S) and Dirac (D) radial      C
C  wave equations. Cubic spline interpolation + power series method.   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  The present package differs from the RADIAL package [F. Salvat, J.M.
C  Fernandez-Varea and W. Williamson, Jr, Comp. Phys. Commun. 90 (1995)
C  151-168], which is distributed by the CPC Program Library. The
C  following changes were introduced on the indicated dates:
C
C  --- January 6, 1997 (F. Salvat): The sign of the Dirac small radial
C  Q-function was reversed to conform with Rose's phase convention. The
C  manual was modified accordingly.
C
C  --- January 21, 1997 (F. Salvat): A correction was introduced to make
C  sure that the matching radius for Dirac bound states is beyond the
C  outer inflexion point of the P-function. Appendix C, on the nodes of
C  the radial functions, was added to the manual.
C
C  --- August 31, 1997 (F. Salvat): Errors were found in eqs. (143),
C  (163) and (164) of the manual. The program was corrected accordingly.
C
C  --- December 9, 1997 (F. Salvat): In the case of free states and pure
C  Coulomb potentials, the matching radius in the original routines was
C  found to be too small. Numerical results were accurate, but the inner
C  phase shift was equal to PI instead of zero and, therefore, radial
C  functions were negative for small R. A safer (larger) value has been
C  adopted.
C
C  --- December 10, 1997 (F. Salvat): An error was found, and corrected,
C  in subroutine DCOUL. This subroutine gave incorrect Dirac Coulomb
C  phase shifts when the argument R was less than the turning point.
C  Programming details of the Coulomb functions were also improved.
C
C  --- July 27, 1999 (J.M. Fernandez-Varea): Function BESJN was declared
C  as an external function to avoid conflict with some compilers that
C  use the same name for an intrinsic function.
C
C  --- April 30, 2000 (F. Salvat): The numerical value of the inverse
C  fine structure constant was updated to SL=137.03599976D0 (CODATA
C  Recommended Values of the Fundamental Constants: 1998).
C
C  --- September 26, 2001 (F. Salvat): The R_infty radial value, which
C  determines the starting point for the inward integration of bound
C  states (see eq. (170) of the manual), was increased by modifying the
C  value of the parameter TRINF in subroutines SINW and DINW.
C
C  --- August 28, 2002 (F. Salvat): The matching radius R_m for free
C  states was allowed to take smaller values. This reduces the calcula-
C  tion time when the outer grid points are widely spaced and allows the
C  calculation of free states with larger angular momenta.
C
C  --- April 24, 2005 (F. Salvat): The source file was reformatted. The
C  labels in final statements of DO loops were replaced by ENDDO
C  statements. Statements in DO loops and in IF()THEN ... ENDIF blocks
C  were indented for better readability. The global sign of the radial
C  functions for free states was redefined to make sure that the
C  magnitude of numerical inner phase shifts tends to zero when the
C  angular momentum l or j tends to infinity. This change was also
C  indicated in the manual.
C
C  --- May 4, 2005 (F. Salvat): The subroutines SGRID (which defines
C  radial grids with prescribed spacing) and ERRSPL (which estimates
C  spline-interpolation errors) have been added to this package.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C  References:
C  [1] F. Salvat and R. Mayol.
C      'Accurate numerical solution of Schrodinger and Dirac wave
C      equations for central fields'.
C      Comput. Phys. Commun. 62 (1991) 65-79.
C  [2] F. Salvat, J. M. Fernandez-Varea and W. Williamson, Jr.
C      'Accurate numerical solution of the radial Schrodinger and Dirac
C      wave equations'.
C      Comput. Phys. Commun. 90 (1995) 151-158.
C  [3] F. Salvat, J. M. Fernandez-Varea and W. Williamson, Jr.
C      'RADIALF: a FORTRAN subroutine package for the solution of the
C      radial Schrodinger and Dirac wave equations'.
C      Internal report, University of Barcelona, 2005.
C      This document is the manual of the RADIALF code system. Its PDF
C      file is included in the distribution package.
C
C
C  It is assumed that the (central) potential energy V(R) is such that
C  the function R*V(R) is finite for all R and tends to constant values
C  when R tends to zero and to infinity.
C
C****   All quantities are in atomic Hartree units.
C  For electrons and positrons
C     unit of length = A0 = 5.291772083D-11 m (= Bohr radius),
C     unit of energy = E0 = 27.2113834 eV (= Hartree energy).
C  For particles of mass 'M' (in units of the electron mass) the atomic
C  units of length and energy are
C     unit of length = A0/M,
C     unit of energy = M*E0.
C
C
C     The calling sequence from the main program is:
C
C****   CALL VINT(R,RV,NV)
C
C  This is an initialisation routine. It determines the natural cubic
C  spline that interpolates the table of values of the function R*V(R)
C  provided by the user.
C   Input arguments:
C     R(I) ..... input potential grid points (repeated values are
C                interpreted as discontinuities).
C     RV(I) .... R(I) times the potential energy at R=R(I).
C     NV ....... number of points in the table (.LE.NDIM).
C  The R(I) grid _must_ include the origin (R=0), and extend up to
C  radial distances for which the function R*V(R) reaches its (constant)
C  asymptotic value. The input grid points must be in non-decreasing
C  order, i.e. R(I+1).GE.R(I).
C
C****   CALL SBOUND(E,EPS,DELL,N,L) or DBOUND(E,EPS,DELL,N,K)
C
C  These subroutines solve the radial wave equations for bound states.
C   Input arguments:
C     E ........ estimated binding energy (a good initial estimate
C                speeds up the calculation).
C     EPS ...... global tolerance, i.e. allowed relative error in the
C                summation of the radial function series.
C     DELL ..... eigenvalue tolerance, i.e. the relative error of the
C                computed eigenvalue is required to be less than DELL
C                (a convenient value is DELL=EPS).
C     N ........ principal quantum number.
C     L ........ orbital angular momentum quantum number.
C     K ........ relativistic angular momentum quantum number.
C                (note: 0.LE.L.LE.N-1, -N.LE.K.LE.N-1)
C   Output argument:
C     E ........ binding energy.
C
C****   CALL SFREE(E,EPS,PHASE,L) or DFREE(E,EPS,PHASE,K)
C
C  These subroutines solve the radial wave equations for free states.
C   Input arguments:
C     E ........ kinetic energy.
C     EPS ...... global tolerance, i.e. allowed relative error in the
C                summation of the radial function series.
C     L ........ orbital angular momentum quantum number.
C     K ........ relativistic angular momentum quantum number.
C                (note: L.GE.0, K.NE.0)
C   Output arguments:
C     PHASE .... inner phase shift (in radians), caused by the short
C                range component of the potential. For modified Coulomb
C                potential, we have
C                       total phase shift = PHASE + DELTA
C                where DELTA is the Coulomb phase shift (that is
C                delivered through the common block
C                       COMMON/OCOUL/RK,ETA,DELTA  ).
C
C
C  The values of the radial functions are delivered through the common
C  block
C     COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C  with NDIM=10000 (if a larger number of grid points is needed, edit
C  this source file and change the value of the parameter NDIM in _all_
C  subroutines ). The grid of radii RAD(I), where the radial wave
C  functions are tabulated, can be arbitrarily selected by the user. The
C  quantities
C     RAD(I) ... user's radial grid,
C     NGP ...... number of grid points (.LE.NDIM),
C  must be defined before calling the solution subroutines. The radial
C  points are sorted in increasing order by the solution subroutines (to
C  avoid confusion, it is advisable to order the input grid in this
C  way). The output quantities are
C     P(I) ..... value of the radial function P(R) at the I-th grid
C                point, R=RAD(I).
C     Q(I) ..... value of the radial function Q(R) at the I-th grid
C                point (= P'(R) for Schrodinger particles).
C     ILAST .... *** Bound states: for R.GT.RAD(ILAST), P(R) and Q(R)
C                are set equal to 0.0D0.
C                *** Free states: for R.GT.RAD(ILAST), P(R) and Q(R)
C                are obtained in terms of the regular and irregular
C                asymptotic Coulomb functions as
C                  P(R)=COS(PHASE)*FU(R)+SIN(PHASE)*GU(R)
C                  Q(R)=COS(PHASE)*FL(R)+SIN(PHASE)*GL(R)
C                where FU, GU and FL, GL are given by subroutines SCOUL
C                and DCOUL with Z=RV(NV). When the absolute value of
C                RV(NV) is less than EPS, Z is set equal to zero, so
C                that the functions FU, GU, FL and GL then reduce to
C                spherical Bessel functions of integer order (given by
C                function BESJN).
C     IER ...... error code. A value larger than zero indicates that
C                some fatal error has been found during the calculation.
C
C
C  Bound state wave functions are normalised to unity. The adopted
C  normalisation for free states is such that P(R) oscillates with unit
C  amplitude in the asymptotic region.
C
C
C****   Error codes (and tentative solutions...):
C    0 .... everything is OK.
C    1 .... EMIN.GE.0 in 'BOUND' (Use a denser grid. If the error
C           persists then probably such a bound state does not exist).
C    2 .... E=0 in 'BOUND' (Probably this bound state does not exist).
C    3 .... RAD(NGP) too small in 'BOUND' (Extend the grid to larger
C           radii).
C    4 .... several zeros of P(R) in a single interval in 'BOUND' (Use
C           a denser grid).
C    5 .... E out of range in 'BOUND' (Accumulated round-off errors?).
C    6 .... RV(NGP)<<0 OR E>0 in 'BOUND' (Check the input potential
C           values).
C    7 .... E.LE.0 in 'FREE'.
C    8 .... RAD(NGP) too small in 'FREE' (Extend the grid to larger
C           radii).
C
C  The program stops when the input quantum numbers are out of range.
C
C  NOTE: The present source file implements the theory described in the
C  accompanying manual, ref. [3] (see above). Numbers in parenthesis in
C  comment lines indicate the relevant equations in that manual.
C
C  *********************************************************************
C                        SUBROUTINE VINT
C  *********************************************************************
      SUBROUTINE VINT(R,RV,NV)
C
C     Natural cubic spline interpolation for R*V(R) from the input radii
C  and potential values (128).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/RGRID/X(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      DIMENSION R(NDIM),RV(NDIM)
C
      IF(R(1).LT.0.0D0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in VINT: R(1).LT.0.')
        STOP
      ENDIF
      IF(NV.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in VINT: input potential grid with ',
     1    'more than ',I5,' data points.')
        STOP
      ENDIF
      R(1)=0.0D0
C
      IO=0
      I=0
      K=0
    1 I=I+1
      K=K+1
      X(K)=R(I)
      Y(K)=RV(I)
      IF(I.EQ.NV) GO TO 2
      IF(R(I).LT.R(I+1)-1.0D-12) GO TO 1
    2 CONTINUE
C
      CALL SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,K)
C
      K=K-1
      DO J=1,K
        IO=IO+1
        RG(IO)=X(J)
        RVG(IO)=Y(J)
        VA(IO)=A(J)
        VB(IO)=B(J)
        VC(IO)=C(J)
        VD(IO)=D(J)
      ENDDO
      IF(I.LT.NV) THEN
        K=0
        GO TO 1
      ENDIF
C  ****  An extra point is added to the grid, and R*V(R) is set equal
C        to RV(NV) for R.GE.R(NV)
      IO=IO+1
      RG(IO)=X(K+1)
      RVG(IO)=Y(K+1)
      VA(IO)=RVG(IO)
      VB(IO)=0.0D0
      VC(IO)=0.0D0
      VD(IO)=0.0D0
      NVT=IO+1
      RG(NVT)=2.0D0*RG(IO)
      RVG(NVT)=RVG(IO)
      VA(NVT)=RVG(IO)
      VB(NVT)=0.0D0
      VC(NVT)=0.0D0
      VD(NVT)=0.0D0
      RETURN
      END
C  *********************************************************************
C                         SUBROUTINE SBOUND
C  *********************************************************************
      SUBROUTINE SBOUND(E,EPS,DELL,N,L)
C
C     This subroutine solves the Schrodinger radial equation for bound
C  states.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
C  ****  Set IWR=1 to print partial results.
      PARAMETER (IWR=0)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
C
      IER=0
      IF(N.LT.1) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in SBOUND: N.LT.1.')
        STOP
      ENDIF
      IF(L.LT.0) THEN
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in SBOUND: L.LT.0.')
        STOP
      ENDIF
      IF(L.GE.N) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in SBOUND: L.GE.N.')
        STOP
      ENDIF
C  ****  Radial quantum number.
      NR=N-L-1
C
      DELL=ABS(DELL)
      IF(E.GT.-1.0D-1) E=-1.0D-1
      FL1=0.5D0*L*(L+1)
C
C  ****  Merge the 'RG' and 'RAD' grids,
C
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      DO 1 I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RG,RLOC,NVT,J)
        TST=MIN(ABS(RLOC-RG(J)),ABS(RLOC-RG(J+1)))
        IF(TST.LT.T) GO TO 1
        NRT=NRT+1
        R(NRT)=RLOC
        IND(NRT)=J
    1 CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Minimum of the effective radial potential. (165)
C
      EMIN=1.0D0
      DO I=2,NRT
        RN=R(I)
        J=IND(I)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EMIN=MIN(EMIN,(RVN+FL1/RN)/RN)
      ENDDO
      IF(EMIN.GT.-1.0D-35) THEN
        IER=1
        WRITE(6,1001)
 1001   FORMAT(1X,'*** Error 1 (in SBOUND): EMIN.GE.0.'/5X,
     1    '(Use a denser grid. If the error persists then probably',
     2    /6X,'such a bound state does not exist).')
        RETURN
      ENDIF
      RN=R(NRT)
      J=IND(NRT)
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
      ESUP=(RVN+FL1/RN)/RN
      IF(ESUP.GT.0.0D0) ESUP=0.0D0
      IF(L.EQ.0) THEN
        TEHYDR=-(R(1)/N)**2
        EMIN=MIN(E,TEHYDR,-10.0D0)
      ENDIF
      IF(E.GT.ESUP.OR.E.LT.EMIN) E=0.5D0*(ESUP+EMIN)
      EMAX=ESUP
      ICMIN=0
      ICMAX=0
C
C  ************  New shot.
C
    2 CONTINUE
      IF(E.GT.-1.0D-16) THEN
        IER=2
        WRITE(6,1002)
 1002   FORMAT(1X,'*** Error 2 (in SBOUND): E=0.',
     1    /5X,'(Probably this bound state does not exist).')
        RETURN
      ENDIF
C  ****  Outer turning point. (165)
      DO K=2,NRT
        IOTP=NRT+2-K
        RN=R(IOTP)
        J=IND(IOTP)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EKIN=E-(RVN+FL1/RN)/RN
        IF(EKIN.GT.0.0D0) GO TO 3
      ENDDO
    3 CONTINUE
      IOTP=IOTP+1
      IF(IOTP.GT.NRT-1) THEN
        IER=3
        WRITE(6,1003)
 1003   FORMAT(1X,'*** Error 3 (in SBOUND): RAD(NGP) too small.'
     1    /5X,'(Extend the grid to larger radii).')
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL SOUTW(E,EPS,L,NR,NZERO,IOTP)
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN
        IER=4
        WRITE(6,1004)
 1004   FORMAT(1X,'*** Error 4 (in SBOUND): Several zeros of P(R)',
     1    /5X,'in a single interval (Use a denser grid).')
        RETURN
      ENDIF
C
C  ****  Too many nodes.
      IF(NZERO.GT.NR) THEN
        IF(ICMIN.EQ.0) EMIN=EMIN-2.0D0*(EMAX-EMIN)
        EMAX=E
        ICMAX=1
        E=0.5D0*(EMIN+E)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,L
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  Too few nodes.
      IF(NZERO.LT.NR) THEN
        ICMIN=1
        EMIN=E
        E=0.5D0*(E+EMAX)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,L
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  The correct number of nodes has been found.
      PO=P(IOTP)
      QO=Q(IOTP)
C
C  ****  Inward solution.
C
      CALL SINW(E,EPS,L,IOTP)
      IF(IER.GT.0) RETURN
C  ****  Matching of the outward and inward solutions.
      FACT=PO/P(IOTP)
      DO I=IOTP,ILAST
        P(I)=P(I)*FACT
        Q(I)=Q(I)*FACT
      ENDDO
      QI=Q(IOTP)
      RLAST=R(ILAST)
C  ****  Normalisation.
      CALL SPLINE(R,P,A,B,C,D,0.0D0,0.0D0,ILAST)
      CALL INTEG2(R,A,B,C,D,0.0D0,RLAST,SUM,ILAST)
C  ****  Eigenvalue correction. (175)
      IF(SUM.LT.1.0D-15) SUM=1.0D0
      DE=PO*(QO-QI)/(SUM+SUM)
      EP=E+DE
C
      IF(DE.LT.0.0D0) THEN
        ICMAX=1
        EMAX=E
      ENDIF
C
      IF(DE.GT.0.0D0) THEN
        ICMIN=1
        EMIN=E
      ENDIF
C
      IF(ICMIN.EQ.0.AND.EP.LT.EMIN) THEN
        EMIN=1.1D0*EMIN
        IF(EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      ENDIF
      IF(ICMIN.EQ.1.AND.EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      IF(ICMAX.EQ.1.AND.EP.GT.EMAX) EP=0.5D0*(E+EMAX)
      IF(EP.GT.ESUP) EP=0.5D0*(ESUP+E)
C
      IF(IWR.NE.0) THEN
        WRITE(6,2000) N,L
 2000   FORMAT(/2X,'Subroutine SBOUND.   N =',I3,'   L =',I3)
        WRITE(6,2001) NR,NZERO,IOTP,NRT
 2001   FORMAT(2X,'NR =',I3,'   NZERO =',I3,'   IOTP = ',I5,
     1    '   NGP =',I5)
        WRITE(6,2002) EP
 2002   FORMAT(2X,'E new = ',1P,D22.15)
        WRITE(6,2003) E,DE
 2003   FORMAT(2X,'E old = ',1P,D22.15,'   DE = ',D11.4)
        WRITE(6,2004) EMIN,EMAX
 2004   FORMAT(2X,'EMIN = ',1P,D12.5,'   EMAX = ',D12.5)
      ENDIF
C
      IF(EP.GE.ESUP.AND.ABS(E-ESUP).LT.DELL*ABS(ESUP)) THEN
        IER=5
        WRITE(6,1005)
 1005   FORMAT(1X,'*** Error 5 (in SBOUND): E out of range.'/5X,
     1    '(Accumulated round-off errors?).')
        RETURN
      ENDIF
      EO=E
      E=EP
      IF(MIN(ABS(DE),ABS(E-EO)).GT.ABS(E*DELL)) GO TO 2
C  ****  Normalisation
      FACT=1.0D0/SQRT(SUM)
      DO I=1,ILAST
        P(I)=P(I)*FACT
        Q(I)=Q(I)*FACT
        IF(ABS(P(I)).LT.1.0D-35) P(I)=0.0D0
        IF(ABS(Q(I)).LT.1.0D-35) Q(I)=0.0D0
      ENDDO
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(R,RLOC,NRT,J)
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          PIO(I)=P(J)
          QIO(I)=Q(J)
        ELSE
          PIO(I)=P(J+1)
          QIO(I)=Q(J+1)
        ENDIF
      ENDDO
C
      IF(ABS(P(NRT)).GT.1.0D-5*ABS(P(IOTP))) THEN
        IER=3
        WRITE(6,1003)
      ENDIF
C
      RLOC=R(ILAST)
      CALL FINDI(RAD,RLOC,NGP,ILAST)
      ILAST=ILAST+1
      RETURN
      END
C  *********************************************************************
C                         SUBROUTINE DBOUND
C  *********************************************************************
      SUBROUTINE DBOUND(E,EPS,DELL,N,K)
C
C     This subroutine solves the Dirac radial equation for bound states.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG,
     1  SL=137.03599976D0)
C  ****  Set IWR=1 to print partial results.
      PARAMETER (IWR=0)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
C
      IER=0
      IF(N.LT.1) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DBOUND: N.LT.1.')
        STOP
      ENDIF
      IF(K.EQ.0) THEN
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in DBOUND: K.EQ.0.')
        STOP
      ENDIF
      IF(K.LT.-N) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in DBOUND: K.LT.-N.')
        STOP
      ENDIF
      IF(K.GE.N) THEN
        WRITE(6,2104)
 2104   FORMAT(1X,'*** Error in DBOUND: K.GE.N.')
        STOP
      ENDIF
C  ****  Orbital angular momentum quantum number. (15)
      IF(K.LT.0) THEN
        L=-K-1
        ELSE
        L=K
      ENDIF
C  ****  Radial quantum number.
      NR=N-L-1
C
      DELL=ABS(DELL)
      IF(E.GT.-1.0D-1) E=-1.0D-1
      FL1=0.5D0*L*(L+1)
C
C  ****  Merge the 'RG' and 'RAD' grids,
C
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
      R(I)=RG(I)
      IND(I)=I
      ENDDO
      NRT=NVT
      DO 1 I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RG,RLOC,NVT,J)
        TST=MIN(ABS(RLOC-RG(J)),ABS(RLOC-RG(J+1)))
        IF(TST.LT.T) GO TO 1
        NRT=NRT+1
        R(NRT)=RLOC
        IND(NRT)=J
    1 CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Minimum of the effective radial potential. (165)
C
      EMIN=1.0D0
      DO I=2,NRT
        RN=R(I)
        J=IND(I)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EMIN=MIN(EMIN,(RVN+FL1/RN)/RN)
      ENDDO
      IF(EMIN.GT.-1.0D-35) THEN
        IER=1
        WRITE(6,1001)
 1001   FORMAT(1X,'*** Error 1 (in DBOUND): EMIN.GE.0.'/5X,
     1    '(Use a denser grid. If the error persists then probably',
     2    /6X,'such a bound state does not exist).')
        RETURN
      ENDIF
      RN=R(NRT)
      J=IND(NRT)
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
      ESUP=(RVN+FL1/RN)/RN
      IF(ESUP.GT.0.0D0) ESUP=0.0D0
      IF(L.EQ.0) THEN
        TEHYDR=-(R(1)/N)**2
        EMIN=MIN(E,TEHYDR,-10.0D0)
      ENDIF
      IF(EMIN.LT.-SL*SL) EMIN=-SL*SL
      IF(E.GT.ESUP.OR.E.LT.EMIN) E=0.5D0*(ESUP+EMIN)
      EMAX=ESUP
      ICMIN=0
      ICMAX=0
C
C  ************  New shot.
C
    2 CONTINUE
      IF(E.GT.-1.0D-16) THEN
        IER=2
        WRITE(6,1002)
 1002   FORMAT(1X,'*** Error 2 (in DBOUND): E=0.',
     1    /5X,'(probably this bound state does not exist).')
        RETURN
      ENDIF
C  ****  Outer turning point. (165)
      DO J=2,NRT
        IOTP=NRT+2-J
        RN=R(IOTP)
        JJ=IND(IOTP)
        RVN=VA(JJ)+RN*(VB(JJ)+RN*(VC(JJ)+RN*VD(JJ)))
        EKIN=E-(RVN+FL1/RN)/RN
        IF(EKIN.GT.0.0D0) GO TO 3
      ENDDO
    3 CONTINUE
      IOTP=IOTP+1
      IF(IOTP.GT.NRT-1) THEN
        IER=3
        WRITE(6,1003)
 1003   FORMAT(1X,'*** Error 3 (in DBOUND): RAD(NGP) too small.'
     1    /5X,'(Extend the grid to larger radii).')
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL DOUTW(E,EPS,K,NR,NZERO,IOTP)
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN
        IER=4
        WRITE(6,1004)
 1004   FORMAT(1X,'*** Error 4 (in DBOUND): Several zeros of P(R)',
     1    /5X,'in a single interval (Use a denser grid).')
        RETURN
      ENDIF
C
C  ****  Too many nodes.
      IF(NZERO.GT.NR) THEN
      IF(ICMIN.EQ.0) EMIN=EMIN-2.0D0*(EMAX-EMIN)
      EMAX=E
      ICMAX=1
      E=0.5D0*(EMIN+E)
      IF(IWR.NE.0) THEN
        WRITE(6,2000) N,K
        WRITE(6,2001) NR,NZERO,IOTP,NRT
        WRITE(6,2002) E
        WRITE(6,2004) EMIN,EMAX
      ENDIF
C
      IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
        IER=5
        WRITE(6,1005)
        RETURN
      ENDIF
C
      GO TO 2
      ENDIF
C  ****  Too few nodes.
      IF(NZERO.LT.NR) THEN
        ICMIN=1
        EMIN=E
        E=0.5D0*(E+EMAX)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,K
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  The correct number of nodes has been found.
      PO=P(IOTP)
      QO=Q(IOTP)
C
C  ****  Inward solution.
C
      CALL DINW(E,EPS,K,IOTP)
      IF(IER.GT.0) RETURN
C  ****  Matching of the outward and inward solutions.
      FACT=PO/P(IOTP)
      DO I=IOTP,NRT
        P(I)=P(I)*FACT
        Q(I)=Q(I)*FACT
      ENDDO
      QI=Q(IOTP)
      RLAST=R(ILAST)
C  ****  Normalisation.
      CALL SPLINE(R,P,A,B,C,D,0.0D0,0.0D0,ILAST)
      CALL INTEG2(R,A,B,C,D,0.0D0,RLAST,SUMP,ILAST)
      CALL SPLINE(R,Q,A,B,C,D,0.0D0,0.0D0,ILAST)
      CALL INTEG2(R,A,B,C,D,0.0D0,RLAST,SUMQ,ILAST)
      SUM=SUMP+SUMQ
C  ****  Eigenvalue correction. (179)
      IF(SUM.LT.1.0D-15) SUM=1.0D0
      DE=SL*PO*(QO-QI)/SUM
      EP=E+DE
C
      IF(DE.LT.0.0D0) THEN
        ICMAX=1
        EMAX=E
      ENDIF
C
      IF(DE.GT.0.0D0) THEN
        ICMIN=1
        EMIN=E
      ENDIF
C
      IF(ICMIN.EQ.0.AND.EP.LT.EMIN) THEN
        EMIN=1.1D0*EMIN
        IF(EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      ENDIF
      IF(ICMIN.EQ.1.AND.EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      IF(ICMAX.EQ.1.AND.EP.GT.EMAX) EP=0.5D0*(E+EMAX)
      IF(EP.GT.ESUP) EP=0.5D0*(ESUP+E)
C
      IF(IWR.NE.0) THEN
        WRITE(6,2000) N,K
 2000   FORMAT(/2X,'Subroutine DBOUND.   N =',I3,'   K =',I3)
        WRITE(6,2001) NR,NZERO,IOTP,NRT
 2001   FORMAT(2X,'NR =',I3,'   NZERO =',I3,'   IOTP = ',I5,
     1    '   NGP =',I5)
        WRITE(6,2002) EP
 2002   FORMAT(2X,'E new = ',1P,D22.15)
        WRITE(6,2003) E,DE
 2003   FORMAT(2X,'E old = ',1P,D22.15,'   DE = ',D11.4)
        WRITE(6,2004) EMIN,EMAX
 2004   FORMAT(2X,'EMIN = ',1P,D12.5,'   EMAX = ',D12.5)
      ENDIF
C
      IF(EP.GT.ESUP.AND.ABS(E-ESUP).LT.DELL*ABS(ESUP)) THEN
        IER=5
        WRITE(6,1005)
 1005   FORMAT(1X,'*** Error 5 (in DBOUND): E out of range.'/5X,
     1    '(Accumulated round-off errors?).')
        RETURN
      ENDIF
      EO=E
      E=EP
      IF(MIN(ABS(DE),ABS(E-EO)).GT.ABS(E*DELL)) GO TO 2
C  ****  Normalisation.
      FACT=1.0D0/SQRT(SUM)
      DO I=1,ILAST
        P(I)=P(I)*FACT
        Q(I)=Q(I)*FACT
        IF(ABS(P(I)).LT.1.0D-35) P(I)=0.0D0
        IF(ABS(Q(I)).LT.1.0D-35) Q(I)=0.0D0
      ENDDO
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(R,RLOC,NRT,J)
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          PIO(I)=P(J)
          QIO(I)=Q(J)
        ELSE
          PIO(I)=P(J+1)
          QIO(I)=Q(J+1)
        ENDIF
      ENDDO
C
      IF(ABS(P(NRT)).GT.1.0D-5*ABS(P(IOTP))) THEN
        IER=3
        WRITE(6,1003)
      ENDIF
C
      RLOC=R(ILAST)
      CALL FINDI(RAD,RLOC,NGP,ILAST)
      ILAST=ILAST+1
      RETURN
      END
C  *********************************************************************
C                         SUBROUTINE SFREE
C  *********************************************************************
      SUBROUTINE SFREE(E,EPS,PHASE,L)
C
C     This subroutine solves the Schrodinger radial equation for free
C  states.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG,
     1  PI=3.1415926535897932D0,TPI=PI+PI,PIH=PI/2.0D0)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
      COMMON/OCOUL/RK,ETA,DELTA
      EXTERNAL BESJN
      ETA=0.0D0
      DELTA=0.0D0
C
      IER=0
      IF(L.LT.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in SFREE: L.LT.0.')
        STOP
      ENDIF
      FL1=0.5D0*L*(L+1)
C
      IF(E.LE.0.0D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7  (in SFREE): E.LE.0.')
        RETURN
      ENDIF
      RK=SQRT(E+E)
C
C  ****  Merge the 'RG' and 'RAD' grids,
C
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      DO 1 I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RG,RLOC,NVT,J)
        TST=MIN(ABS(RLOC-RG(J)),ABS(RLOC-RG(J+1)))
        IF(TST.LT.T) GO TO 1
        NRT=NRT+1
        R(NRT)=RLOC
        IND(NRT)=J
    1 CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      ZINF=RV(NVT)
      IF(ABS(ZINF).LT.EPS) THEN
C  ****  Finite range potentials.
        ILAST=NRT+1
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          T=EPS*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(RVN).GT.T) GO TO 2
          BNL1=BESJN(2,L+1,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 2  ! Test cutoff.
          BNL=BESJN(2,L,X)
          BJL=BESJN(1,L,X)
          BJL1=BESJN(1,L+1,X)
          ILAST=IL
          PA(ILAST)=X*BJL
          PB(ILAST)=-X*BNL
          QA(ILAST)=RK*((L+1.0D0)*BJL-X*BJL1)
          QB(ILAST)=-RK*((L+1.0D0)*BNL-X*BNL1)
        ENDDO
    2   CONTINUE
        IF(ILAST.EQ.NRT+1) THEN
          IER=8
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Error 8  (in SFREE): RAD(NGP) too small.'
     1      /5X,'(Extend the grid to larger radii).')
          RETURN
        ENDIF
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(ZINF)
        ILAST=NRT+1
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          IF(ABS(RVN-ZINF).GT.TAS) GO TO 3
          CALL SCOUL(ZINF,E,L,RN,P0,Q0,P1,Q1,ERR)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          PA(ILAST)=P0
          PB(ILAST)=P1
          QA(ILAST)=Q0
          QB(ILAST)=Q1
        ENDDO
    3   CONTINUE
        IF(ILAST.EQ.NRT+1) THEN
          IER=8
          WRITE(6,1008)
          RETURN
        ENDIF
      ENDIF
C
C  ****  Outward solution.
C
      CALL SOUTW(E,EPS,L,1,NZERO,ILAST)
C
C  ****  Phase shift. (187)
C
      IF(ABS(P(ILAST)).GT.EPS) THEN
        RATIO=Q(ILAST)/P(ILAST)
        PHASE=DATAN2(RATIO*PA(ILAST)-QA(ILAST),
     1              QB(ILAST)-RATIO*PB(ILAST))
        CD=COS(PHASE)
        SD=SIN(PHASE)
        RNORM=(CD*PA(ILAST)+SD*PB(ILAST))/P(ILAST)
      ELSE
        PHASE=DATAN2(-PA(ILAST),PB(ILAST))
        CD=COS(PHASE)
        SD=SIN(PHASE)
        RNORM=(CD*QA(ILAST)+SD*QB(ILAST))/Q(ILAST)
      ENDIF
C  ****  The calculated phase shift is reduced to the interval
C        (-PI/2,+PI/2) by subtracting or adding PI.
      IF(PHASE.GT.PIH) THEN
        PHASE=PHASE-PI
        CD=-CD
        SD=-SD
        RNORM=-RNORM
      ELSE IF(PHASE.LT.-PIH) THEN
        PHASE=PHASE+PI
        CD=-CD
        SD=-SD
        RNORM=-RNORM
      ENDIF
C  ****  Normalised wave function. (188)
      DO I=1,ILAST
        P(I)=RNORM*P(I)
        Q(I)=RNORM*Q(I)
        IF(ABS(P(I)).LT.1.0D-35) P(I)=0.0D0
        IF(ABS(Q(I)).LT.1.0D-35) Q(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          P(I)=CD*PA(I)+SD*PB(I)
          Q(I)=CD*QA(I)+SD*QB(I)
          IF(ABS(P(I)).LT.1.0D-35) P(I)=0.0D0
          IF(ABS(Q(I)).LT.1.0D-35) Q(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(R,RLOC,NRT,J)
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          PIO(I)=P(J)
          QIO(I)=Q(J)
        ELSE
          PIO(I)=P(J+1)
          QIO(I)=Q(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RAD,RLOC,NGP,ILAST)
      ILAST=ILAST+1
      IF(ILAST.GT.NGP) ILAST=NGP
      RETURN
      END
C  *********************************************************************
C                         SUBROUTINE DFREE
C  *********************************************************************
      SUBROUTINE DFREE(E,EPS,PHASE,K)
C
C     This subroutine solves the Dirac radial equation for free states.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG,
     1  SL=137.03599976D0,PI=3.1415926535897932D0,
     2  TPI=PI+PI,PIH=PI/2.0D0)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
      COMMON/OCOUL/RK,ETA,DELTA
      EXTERNAL BESJN
      ETA=0.0D0
      DELTA=0.0D0
C
      IER=0
      IF(K.EQ.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DFREE: K.EQ.0.')
        STOP
      ENDIF
C
      IF(E.LE.0.0D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** ERROR 7  (IN DFREE): E.LE.0.')
        RETURN
      ENDIF
C  ****  Orbital angular momentum quantum number. (15)
      IF(K.LT.0) THEN
        L=-K-1
        KSIGN=1
      ELSE
        L=K
        KSIGN=-1
      ENDIF
      FL1=0.5D0*L*(L+1)
      RK=SQRT(E*(E+2.0D0*SL*SL))/SL
C
C  ****  Merge the 'RG' and 'RAD' grids,
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in DFREE: Input potential grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      DO 1 I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RG,RLOC,NVT,J)
        TST=MIN(ABS(RLOC-RG(J)),ABS(RLOC-RG(J+1)))
        IF(TST.LT.T) GO TO 1
        NRT=NRT+1
        R(NRT)=RLOC
        IND(NRT)=J
    1 CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      ZINF=RV(NVT)
      IF(ABS(ZINF).LT.EPS) THEN
C  ****  Finite range potentials.
        FACTOR=SQRT(E/(E+2.0D0*SL*SL))
        ILAST=NRT+1
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          T=EPS*RN*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(RVN).GT.T) GO TO 2
          BNL=BESJN(2,L,X)
          IF(ABS(BNL).GT.1.0D6) GO TO 2  ! Test cutoff.
          BNL1=BESJN(2,L+KSIGN,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 2  ! Test cutoff.
          BJL=BESJN(1,L,X)
          BJL1=BESJN(1,L+KSIGN,X)
          ILAST=IL
          PA(ILAST)=X*BJL
          PB(ILAST)=-X*BNL
          QA(ILAST)=-FACTOR*KSIGN*X*BJL1
          QB(ILAST)=FACTOR*KSIGN*X*BNL1
        ENDDO
    2   CONTINUE
        IF(ILAST.EQ.NRT+1) THEN
          IER=8
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Error 8  (in DFREE): RAD(NGP) too small.'
     1      /5X,'(Extend the grid to larger radii).')
          RETURN
        ENDIF
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(ZINF)
        ILAST=NRT+1
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          IF(ABS(RVN-ZINF).GT.TAS) GO TO 3
          CALL DCOUL(ZINF,E,K,RN,P0,Q0,P1,Q1,ERR)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          PA(ILAST)=P0
          PB(ILAST)=P1
          QA(ILAST)=Q0
          QB(ILAST)=Q1
        ENDDO
    3   CONTINUE
        IF(ILAST.EQ.NRT+1) THEN
          IER=8
          WRITE(6,1008)
          RETURN
        ENDIF
      ENDIF
C
C  ****  Outward solution.
C
      CALL DOUTW(E,EPS,K,1,NZERO,ILAST)
C
C  ****  Phase shift. (187)
C
      RM=R(ILAST)
      IL=IND(ILAST-1)
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
      FG=(E-VF+2.0D0*SL*SL)/SL
      PO=P(ILAST)
      POP=-K*PO/RM+FG*Q(ILAST)
      IL=IND(ILAST)
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
      FG=(E-VF+2.0D0*SL*SL)/SL
      PIA=PA(ILAST)
      PIAP=-K*PIA/RM+FG*QA(ILAST)
      PIB=PB(ILAST)
      PIBP=-K*PIB/RM+FG*QB(ILAST)
C
      IF(ABS(PO).GT.EPS) THEN
        RATIO=POP/PO
        PHASE=DATAN2(RATIO*PIA-PIAP,PIBP-RATIO*PIB)
        CD=COS(PHASE)
        SD=SIN(PHASE)
        RNORM=(CD*PIA+SD*PIB)/PO
      ELSE
        PHASE=DATAN2(-PIA,PIB)
        CD=COS(PHASE)
        SD=SIN(PHASE)
        RNORM=(CD*PIAP+SD*PIBP)/POP
      ENDIF
C  ****  The calculated phase shift is reduced to the interval
C        (-PI/2,+PI/2) by subtracting or adding PI.
      IF(PHASE.GT.PIH) THEN
        PHASE=PHASE-PI
        CD=-CD
        SD=-SD
        RNORM=-RNORM
      ELSE IF(PHASE.LT.-PIH) THEN
        PHASE=PHASE+PI
        CD=-CD
        SD=-SD
        RNORM=-RNORM
      ENDIF
C  ****  Normalised wave function. (188)
      DO I=1,ILAST
        P(I)=RNORM*P(I)
        Q(I)=RNORM*Q(I)
        IF(ABS(P(I)).LT.1.0D-35) P(I)=0.0D0
        IF(ABS(Q(I)).LT.1.0D-35) Q(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          P(I)=CD*PA(I)+SD*PB(I)
          Q(I)=CD*QA(I)+SD*QB(I)
          IF(ABS(P(I)).LT.1.0D-35) P(I)=0.0D0
          IF(ABS(Q(I)).LT.1.0D-35) Q(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(R,RLOC,NRT,J)
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          PIO(I)=P(J)
          QIO(I)=Q(J)
        ELSE
          PIO(I)=P(J+1)
          QIO(I)=Q(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RAD,RLOC,NGP,ILAST)
      ILAST=ILAST+1
      IF(ILAST.GT.NGP) ILAST=NGP
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SOUTW
C  *********************************************************************
      SUBROUTINE SOUTW(E,EPS,L,NR,NZERO,IOTP)
C
C     Outward solution of the Schrodinger radial equation for a  piece-
C  wise cubic potential. Power series method.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS
      COMMON/NZT/NZMAX
      NZERO=0
      NZMAX=0
      AL=L
      N1=IOTP-1
C
      P(1)=0.0D0
      Q(1)=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PI=P(I)
        QI=Q(I)
        CALL SCH(E,AL,EPS)
        NZERO=NZERO+NCHS
        IF(NCHS.GT.NZMAX) NZMAX=NCHS
        IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN
        P(I+1)=PF
        Q(I+1)=QF
        IF(I.EQ.1) GO TO 1
C  ****  Renormalisation.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO K=1,I
          P(K)=P(K)*FACT
          Q(K)=Q(K)*FACT
          ENDDO
        ENDIF
    1 CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SINW
C  *********************************************************************
      SUBROUTINE SINW(E,EPS,L,IOTP)
C
C     Inward solution of the Schrodinger radial equation for a piece-
C  wise cubic potential. Power series method.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG,
     1  TRINF=22500.0D0)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS
      AL=L
C  ****  WKB solution at the outer grid point. (168)
      N=NRT
    1 N1=IND(N-1)
      RN=R(N)
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))
      CMU=2.0D0*RN*(RVN-E*RN)+AL*(AL+1)
      IF(CMU.LE.0.0D0) THEN
        IER=6
        WRITE(6,1006)
 1006   FORMAT(1X,'*** Error 6 (in SBOUND): RV(NGP)<<0 OR E>0.',
     1    /5X,'(Check the input potential values).')
        RETURN
      ENDIF
C  ****  Practical infinity. (170)
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN
        CRAT=1.0D0-RN*(RVN+RN*RVNP-2*E*RN)/(CMU*CMU)
        P(N)=1.0D0
        Q(N)=(0.5D0/RN)*CRAT-SQRT(CMU)/RN
        ILAST=N
      ELSE
        P(N)=0.0D0
        Q(N)=0.0D0
        N=N-1
        GO TO 1
      ENDIF
C
      N1=N-IOTP
      DO J=1,N1
        I=N-J
        I1=I+1
        RA=R(I1)
        RB=R(I)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PI=P(I1)
        QI=Q(I1)
        CALL SCH(E,AL,EPS)
        P(I)=PF
        Q(I)=QF
C  ****  Renormalisation.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO K=I1,N
            P(K)=P(K)*FACT
            Q(K)=Q(K)*FACT
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DOUTW
C  *********************************************************************
      SUBROUTINE DOUTW(E,EPS,K,NR,NZERO,IOTP)
C
C     Outward solution of the Dirac radial equation for a piecewise
C  cubic potential. Power series method.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,NSUM
      COMMON/NZT/NZMAX
      NZERO=0
      NZMAX=0
      AK=K
      IF(E.LT.0.0D0) THEN
        N1=NRT
      ELSE
        N1=IOTP-1
      ENDIF
C
      P(1)=0.0D0
      Q(1)=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PI=P(I)
        QI=Q(I)
        CALL DIR(E,AK,EPS)
        NZERO=NZERO+NCHS
        IF(NCHS.GT.NZMAX) NZMAX=NCHS
        IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN
        P(I+1)=PF
        Q(I+1)=QF
        IF(E.LT.0.0D0) THEN
C  ****  TCONV is the product of P and its second derivative at the
C        I-th grid point (positive if P is convex).
          TCONV=2.0D0*CA(3)*PI
          IF(I.GE.IOTP.AND.TCONV.GT.1.0D-15) THEN
            IOTP=I+1
            RETURN
          ENDIF
        ENDIF
        IF(I.EQ.1) GO TO 1
C  ****  Renormalisation.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO J=1,I
            P(J)=P(J)*FACT
            Q(J)=Q(J)*FACT
          ENDDO
        ENDIF
    1 CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DINW
C  *********************************************************************
      SUBROUTINE DINW(E,EPS,K,IOTP)
C
C     Inward solution of the Dirac radial equation for a piecewise cubic
C  potential. Power series method.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG,
     1  TRINF=22500.0D0,SL=137.03599976D0)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS
C  ****  Orbital angular momentum quantum number. (15)
      IF(K.LT.0) THEN
        L=-K-1
      ELSE
        L=K
      ENDIF
      AK=K
      AL=L
C  ****  WKB solution at the outer grid point. (177)
      N=NRT
      FACT=(E+2.0D0*SL*SL)/(SL*SL)
    1 N1=IND(N-1)
      RN=R(N)
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))
      CMU=FACT*RN*(RVN-E*RN)+AL*(AL+1)
      IF(CMU.LE.0.0D0) THEN
        IER=6
        WRITE(6,1006)
 1006   FORMAT(1X,'*** Error 6 (in DBOUND): RV(NGP)<<0 OR E>0.',
     1    /5X,'(Check the input potential values).')
        RETURN
      ENDIF
C  ****  Practical infinity. (170)
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN
        CRAT=(0.5D0-SQRT(CMU))/RN-0.25D0*FACT*(RVN+RN*RVNP
     1      -2.0D0*RN*E)/CMU
        P(N)=1.0D0
        Q(N)=SL*(CRAT+AK/RN)/(E+2.0D0*SL*SL)
        ILAST=N
      ELSE
        P(N)=0.0D0
        Q(N)=0.0D0
        N=N-1
        GO TO 1
      ENDIF
C
      N1=N-IOTP
      DO J=1,N1
        I=N-J
        I1=I+1
        RA=R(I1)
        RB=R(I)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PI=P(I1)
        QI=Q(I1)
        CALL DIR(E,AK,EPS)
        P(I)=PF
        Q(I)=QF
C  ****  Renormalisation.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO M=I1,N
            P(M)=P(M)*FACT
            Q(M)=Q(M)*FACT
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCH
C  *********************************************************************
      SUBROUTINE SCH(E,AL,EPS)
C
C     This subroutine solves the Schrodinger radial equation for a
C  a central potential V(R) such that
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C  Given the boundary conditions (i.e. the value of the radial function
C  and its derivative) at RA, the solution in the interval between RA
C  and RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AL .................... orbital angular momentum quantum number.
C
C  Input (common POTEN):
C     RV0, RV1, RV2, RV3 .... potential parameters.
C
C  Input-output (common SINOUT):
C     RA, RB ................ interval end points (input),
C     PI, QI ................ values of the radial function and its
C                             derivative at RA (input),
C     PF, QF ................ values of the radial function and its
C                             derivative at RB (output),
C     RLN ................... LOG of the re-normalising factor,
C     NSTEP ................. number of steps,
C     NCHS .................. number of zeros of P(R) in (RA,RB).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS
      COMMON/SSAVE/P0,Q0,P1,Q1,CA(60),R0,R1,NSUM
      NCHS=0
      RLN=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
        DIRECT=-1.0D0
      ELSE
        DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      P1=PI
      Q1=QI
    1 CONTINUE
      R0=R1
      P0=P1
      Q0=Q1
    2 CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL SCH0(E,AL,EPS)
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
    3 CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
    4 CONTINUE
      NSTEP=NSTEP+1
      TST=ABS(P1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalisation.
        RLN=RLN+LOG(TST)
        P1=P1/TST
        Q1=Q1/TST
      ENDIF
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      PF=P1
      QF=Q1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCH0
C  *********************************************************************
      SUBROUTINE SCH0(E,AL,EPS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Overflow level.
      PARAMETER (OVER=1.0D15)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SSAVE/P0,Q0,P1,Q1,CA(60),R0,R1,NSUM
C
      RVE=RV1-E
      IF(R0.GT.1.0D-10) GO TO 2
C
C  ****  First interval. (134-143)
C
      S=AL+1
      U0=AL*S
      U1=2*RV0*R1
      U2=2*RVE*R1**2
      U3=2*RV2*R1**3
      U4=2*RV3*R1**4
      UT=U0+U1+U2+U3+U4
C
      CA(1)=1.0D0
      CA(2)=U1*CA(1)/((S+1)*S-U0)
      CA(3)=(U1*CA(2)+U2*CA(1))/((S+2)*(S+1)-U0)
      CA(4)=(U1*CA(3)+U2*CA(2)+U3*CA(1))
     1     /((S+3)*(S+2)-U0)
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      Q1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      P2P1=S*(S-1)*CA(1)+(S+1)*S*CA(2)+(S+2)*(S+1)*CA(3)
     1     +(S+3)*(S+2)*CA(4)
C
      DO I=5,60
        K=I-1
        CA(I)=(U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)+U4*CA(I-4))
     1       /((S+K)*(S+K-1)-U0)
        P1=P1+CA(I)
        DQ1=(S+K)*CA(I)
        Q1=Q1+DQ1
        P2P1=P2P1+(S+K-1)*DQ1
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(P2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(CA(I))
        T2=ABS(R1*R1*(P2P1-UT*P1))
        TST1=EPS*MAX(ABS(P1),ABS(Q1)/I)
        TST2=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 1
      ENDDO
C  ****  Renormalisation. (143)
    1 CONTINUE
      NSUM=K+1
      Q1=Q1/(ABS(P1)*R1)
      P1=P1/ABS(P1)
      RETURN
C
C  ****  Middle region. (134-139)
C
    2 CONTINUE
      H=R1-R0
      H2=H*H
C
      RHO=H/R0
      U0=AL*(AL+1)+2*R0*(RV0+R0*(RVE+R0*(RV2+R0*RV3)))
      U1=2*(RV0+R0*(2*RVE+R0*(3*RV2+R0*4*RV3)))*H
      U2=2*(RVE+R0*(3*RV2+R0*6*RV3))*H2
      U3=2*(RV2+R0*4*RV3)*H2*H
      U4=2*RV3*H2*H2
      UT=U0+U1+U2+U3+U4
C
      CA(1)=P0
      CA(2)=Q0*H
      CA(3)=RHO*RHO*U0*CA(1)/2
      CA(4)=RHO*(RHO*(U0*CA(2)+U1*CA(1))-4*CA(3))/6
      CAK=(U0-2)*CA(3)+U1*CA(2)+U2*CA(1)
      CA(5)=RHO*(RHO*CAK-12*CA(4))/12
      CAK=(U0-6)*CA(4)+U1*CA(3)+U2*CA(2)+U3*CA(1)
      CA(6)=RHO*(RHO*CAK-24*CA(5))/20
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)+CA(5)+CA(6)
      Q1=CA(2)+2*CA(3)+3*CA(4)+4*CA(5)+5*CA(6)
      P2P1=2*CA(3)+6*CA(4)+12*CA(5)+20*CA(6)
C
      DO I=7,60
        K=I-1
        CAK=(U0-(K-2)*(K-3))*CA(I-2)+U1*CA(I-3)+U2*CA(I-4)
     1     +U3*CA(I-5)+U4*CA(I-6)
        CA(I)=RHO*(RHO*CAK-2*(K-1)*(K-2)*CA(K))/(K*(K-1))
        P1=P1+CA(I)
        DQ1=K*CA(I)
        Q1=Q1+DQ1
        P2P1=P2P1+K*(K-1)*CA(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(P2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(CA(I))
        T2=ABS(R1*R1*P2P1-H2*UT*P1)
        TST1=EPS*MAX(ABS(P1),ABS(Q1)/I)
        TST2=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 3
      ENDDO
C
    3 NSUM=K+1
      Q1=Q1/H
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DIR
C  *********************************************************************
      SUBROUTINE DIR(E,AK,EPS)
C
C     This subroutine solves the Dirac radial equation for a central
C  potential V(R) such that
C             R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C  Given the boundary conditions (i.e. the value of the large and small
C  radial functions) at RA, the solution in the interval between RA and
C  RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AK .................... relativistic angular momentum quantum
C                             number.
C
C  Input (common POTEN):
C     RV0, RV1, RV2, RV3 .... potential parameters.
C
C  Input-output (common DINOUT):
C     RA, RB ................ interval end points (input),
C     PI, QI ................ values of the large and small radial
C                             functions at RA (input),
C     PF, QF ................ values of the large and small radial
C                             functions at RB (output),
C     RLN ................... LOG of the re-normalising factor,
C     EPS ................... estimate of the global error in PF and QF,
C     NSTEP ................. number of steps,
C     NCHS .................. number of zeros of P(R) in (RA,RB).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,NSUM
      NCHS=0
      RLN=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
      DIRECT=-1.0D0
      ELSE
      DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      P1=PI
      Q1=QI
    1 CONTINUE
      R0=R1
      P0=P1
      Q0=Q1
    2 CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL DIR0(E,AK,EPS)
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
    3 CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
    4 CONTINUE
      NSTEP=NSTEP+1
      TST=ABS(P1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalisation.
        RLN=RLN+LOG(TST)
        P1=P1/TST
        Q1=Q1/TST
      ENDIF
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      PF=P1
      QF=Q1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DIR0
C  *********************************************************************
      SUBROUTINE DIR0(E,AK,EPS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Speed of light and overflow level.
      PARAMETER (SL=137.03599976D0,OVER=1.0D15)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,NSUM
C
      ISIG=1
      IF(AK.GT.0.0D0) ISIG=-1
      H=R1-R0
      H2=H*H
      RVE=RV1-E
C
      IF(R0.GT.1.0D-10) GO TO 4
C
C  ****  First interval. (147-164)
C
      U0=RV0/SL
      U1=RVE*R1/SL
      U2=RV2*R1**2/SL
      U3=RV3*R1**3/SL
      UT=U0+U1+U2+U3
      UQ=UT-2*SL*R1
      UH=U1-2*SL*R1
      IF(ABS(U0).LT.1.0D-10) GO TO 1
C
C  ****  U0.NE.0. (155-159)
      S=SQRT(AK*AK-U0*U0)
      DS=S+S
      CA(1)=1.0D0
      CB(1)=-(S+AK)/U0
      CAI=U1*CA(1)
      CBI=UH*CB(1)
      CA(2)=(-U0*CAI-(S+1-AK)*CBI)/(DS+1)
      CB(2)=((S+1+AK)*CAI-U0*CBI)/(DS+1)
      CAI=U1*CA(2)+U2*CA(1)
      CBI=UH*CB(2)+U2*CB(1)
      CA(3)=(-U0*CAI-(S+2-AK)*CBI)/(2*(DS+2))
      CB(3)=((S+2+AK)*CAI-U0*CBI)/(2*(DS+2))
      P1=CA(1)+CA(2)+CA(3)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)
      Q1=CB(1)+CB(2)+CB(3)
      QP1=S*CB(1)+(S+1)*CB(2)+(S+2)*CB(3)
C
      DO I=4,60
        K=I-1
        CAI=U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)
        CBI=UH*CB(K)+U2*CB(I-2)+U3*CB(I-3)
        CA(I)=(-U0*CAI-(S+K-AK)*CBI)/(K*(DS+K))
        CB(I)=((S+K+AK)*CAI-U0*CBI)/(K*(DS+K))
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+K)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=1. (160,161)
    1 CONTINUE
      IF(ISIG.LT.0) GO TO 2
      S=ABS(AK)
      DS1=S+S+1
      CA(1)=1.0D0
      CB(1)=U1*CA(1)/DS1
      CA(2)=0.0D0
      CB(2)=U2*CA(1)/(DS1+1)
      CA(3)=-UH*CB(1)/2
      CB(3)=(U1*CA(3)+U3*CA(1))/(DS1+2)
      CA(4)=-(UH*CB(2)+U2*CB(1))/3
      CB(4)=(U1*CA(4)+U2*CA(3))/(DS1+3)
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=(S+1)*CB(1)+(S+2)*CB(2)+(S+3)*CB(3)
C
      DO I=5,60
        K=I-1
        CA(I)=-(UH*CB(I-2)+U2*CB(I-3)+U3*CB(I-4))/K
        CB(I)=(U1*CA(I)+U2*CA(K)+U3*CA(I-2))/(DS1+K)
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+I)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=-1. (162,163)
    2 CONTINUE
      S=ABS(AK)+1
      DS1=S+ABS(AK)
      IF(UH.GT.0.0D0) THEN
      CB(1)=-1.0D0
      ELSE
      CB(1)=1.0D0
      ENDIF
      CA(1)=-UH*CB(1)/DS1
      CB(2)=0.0D0
      CA(2)=-U2*CB(1)/(DS1+1)
      CB(3)=U1*CA(1)/2
      CA(3)=-(UH*CB(3)+U3*CB(1))/(DS1+2)
      CB(4)=(U1*CA(2)+U2*CA(1))/3
      CA(4)=-(UH*CB(4)+U2*CB(3))/(DS1+3)
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=(S-1)*CB(1)+S*CB(2)+(S+1)*CB(3)
C
      DO I=5,60
        K=I-1
        CB(I)=(U1*CA(I-2)+U2*CA(I-3)+U3*CA(I-4))/K
        CA(I)=-(UH*CB(I)+U2*CB(K)+U3*CB(I-2))/(DS1+K)
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+K-1)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
C  ****  Renormalisation. (164)
    3 CONTINUE
      NSUM=K+1
      Q1=Q1/ABS(P1)
      P1=P1/ABS(P1)
      RETURN
C
C  ****  Middle region. (148-152)
C
    4 CONTINUE
      RHO=H/R0
      U0=(RV0+R0*(RVE+R0*(RV2+R0*RV3)))/SL
      U1=(RVE+R0*(2*RV2+R0*3*RV3))*H/SL
      U2=(RV2+R0*3*RV3)*H2/SL
      U3=RV3*H*H2/SL
      UB=U0-2*SL*R0
      UH=U1-2*SL*H
      UT=U0+U1+U2+U3
      UQ=UT-2*SL*R1
C
      CA(1)=P0
      CB(1)=Q0
      CA(2)=-RHO*(AK*CA(1)+UB*CB(1))
      CB(2)=RHO*(AK*CB(1)+U0*CA(1))
      CA(3)=-RHO*((AK+1)*CA(2)+UB*CB(2)+UH*CB(1))/2
      CB(3)=RHO*((AK-1)*CB(2)+U0*CA(2)+U1*CA(1))/2
      CA(4)=-RHO*((AK+2)*CA(3)+UB*CB(3)+UH*CB(2)+U2*CB(1))/3
      CB(4)=RHO*((AK-2)*CB(3)+U0*CA(3)+U1*CA(2)+U2*CA(1))/3
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=CA(2)+2*CA(3)+3*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=CB(2)+2*CB(3)+3*CB(4)
C
      DO I=5,60
        K=I-1
        CA(I)=-RHO*((AK+K-1)*CA(K)+UB*CB(K)+UH*CB(I-2)+U2*CB(I-3)
     1       +U3*CB(I-4))/K
        CB(I)=RHO*((AK-K+1)*CB(K)+U0*CA(K)+U1*CA(I-2)+U2*CA(I-3)
     1       +U3*CA(I-4))/K
        P1=P1+CA(I)
        PP1=PP1+K*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+K*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 5
      ENDDO
C
    5 CONTINUE
      NSUM=K+1
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC         Coulomb and Bessel functions         CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE SCOUL
C  *********************************************************************
      SUBROUTINE SCOUL(Z,E,L,R,F,FP,G,GP,ERR)
C
C     This subroutine computes radial Schrodinger-Coulomb wave functions
C  for free states.
C
C  **** All quantities in atomic units.
C
C  Input arguments:
C     Z ........ potential strength, i.e. value of R*V(R) (assumed
C                constant),
C     E ........ particle kinetic energy (positive),
C     L ........ orbital angular momentum quantum number (.GE.0),
C     R ........ radial distance (positive).
C
C  Output arguments:
C     F, FP .... regular Schrodinger-Coulomb function and its
C                derivative,
C     G, GP .... irregular Schrodinger-Coulomb function and its
C                derivative,
C     ERR ...... accuracy of the computed functions (relative
C                uncertainty).
C
C  Output through common /OCOUL/:
C     WAVNUM ... wave number,
C     ETA ...... Sommerfeld's parameter,
C     DELTA .... Coulomb phase shift (modulus 2*PI).
C
C     Radial functions are normalised so that, for large R, they
C  oscillate with unit amplitude.
C
C     Other subprograms required: subroutines FCOUL and SUM2F0,
C                                 and function CLGAM.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      COMMON/OCOUL/WAVNUM,ETA,DELTA
      EXTERNAL BESJN
C
C  ****  Parameters.
C
      WAVNUM=SQRT(E+E)
      IF(ABS(Z).GT.0.00001D0) THEN
        ETA=Z/WAVNUM
        ICAL=0
      ELSE
        ETA=0.0D0
        ICAL=1
      ENDIF
      RLAMB=L
      X=WAVNUM*R
C
      IF(E.LT.0.0001D0.OR.L.LT.0) THEN
        F=0.0D0
        FP=0.0D0
        G=1.0D35
        GP=-1.0D35
        DELTA=0.0D0
        ERR=1.0D0
        IF(E.LT.0.0001D0) WRITE(6,2101)
 2101   FORMAT(1X,'*** ERROR IN SCOUL: E IS TOO SMALL.')
        IF(L.LT.0) WRITE(6,2102)
 2102   FORMAT(1X,'*** ERROR IN SCOUL: L.LT.0.')
        RETURN
      ENDIF
      IF(ICAL.EQ.1) GO TO 1
C
C  ************  Coulomb functions.
C
      CALL FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)
      FP=FP*WAVNUM
      GP=GP*WAVNUM
      DELTA=DELTAC(ETA,RLAMB)
      IF(ERR.GE.1.0D-6) THEN
        F=0.0D0
        FP=0.0D0
        G=1.0D35
        GP=-1.0D35
        ERR=1.0D0
      ENDIF
      RETURN
C
C  ************  Z=0. Spherical Bessel functions.
C
    1 CONTINUE
      F=X*BESJN(1,L,X)
      FP=((L+1)*BESJN(1,L,X)-X*BESJN(1,L+1,X))*WAVNUM
      G=-X*BESJN(2,L,X)
      GP=-((L+1)*BESJN(2,L,X)-X*BESJN(2,L+1,X))*WAVNUM
      DELTA=0.0D0
      ERR=0.0D0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DCOUL
C  *********************************************************************
      SUBROUTINE DCOUL(Z,E,K,R,FU,FL,GU,GL,ERR)
C
C     This subroutine computes radial Dirac-Coulomb wave functions for
C  free states.
C
C  **** All quantities in atomic units.
C
C  Input arguments:
C     Z ........ potential strength, i.e. value of R*V(R) (assumed
C                constant),
C     E ........ particle kinetic energy (positive),
C     K ........ angular momentum quantum number KAPPA (.NE.0),
C     R ........ radial distance (positive).
C
C  Output arguments:
C     FU, FL ... upper and lower components of the regular Dirac-
C                Coulomb function,
C     GU, GL ... upper and lower components of the irregular Dirac-
C                Coulomb function,
C     ERR ...... accuracy of the computed functions (relative
C                uncertainty).
C
C  Output through common /OCOUL/:
C     WAVNUM ... wave number,
C     ETA ...... Sommerfeld's parameter,
C     DELTA .... Coulomb phase shift (modulus 2*PI).
C
C     Radial functions are normalised so that, for large r, the upper
C  component function oscillates with unit amplitude.
C
C     Other subprograms required: subroutines FCOUL and SUM2F0,
C                                 and function CLGAM.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI,
     1  SL=137.03599976D0,SL2=SL*SL,TSL2=SL2+SL2,ALPHA=1.0D0/SL)
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C
      IF(ABS(Z).GT.0.00001D0) THEN
        ZETA=Z*ALPHA
        ICAL=0
      ELSE
        ZETA=0.0D0
        ICAL=1
      ENDIF
      RLAMBS=K*K-ZETA*ZETA
      RLAMB=SQRT(RLAMBS)
      PC=SQRT(E*(E+TSL2))
      WAVNUM=PC/SL
      X=WAVNUM*R
C
      IF(E.LT.0.0001D0.OR.K.EQ.0) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERR=1.0D0
        DELTA=0.0D0
        IF(E.LT.0.0001D0) WRITE(6,2101)
 2101   FORMAT(1X,'*** ERROR IN DCOUL: E IS TOO SMALL.')
        IF(K.EQ.0) WRITE(6,2102)
 2102   FORMAT(1X,'*** ERROR IN DCOUL: K.EQ.0.')
        RETURN
      ENDIF
      IF(ICAL.EQ.1) GO TO 1
C
C  ****  Parameters.
C
      RLAMB1=RLAMB-1.0D0
      W=E+SL2
      ETA=ZETA*W/PC
      RLA=SQRT(RLAMBS+ETA*ETA)
      P1=K+RLAMB
      P2=RLAMB*SL2-K*W
      RNUR=ZETA*(W+SL2)
      RNUI=-P1*PC
      RNU=DATAN2(RNUI,RNUR)
      RNORM=1.0D0/(SQRT(RNUR*RNUR+RNUI*RNUI)*RLAMB)
C
C  ****  Coulomb phase shift.
C
      IF(K.GT.0) THEN
        L=K
      ELSE
        L=-K-1
      ENDIF
      DELTA0=DELTAC(ETA,RLAMB1)
      DELTA=RNU-(RLAMB-L-1)*PIH+DELTA0
      IF(Z.LT.0.0D0.AND.K.LT.0) THEN
        RNORM=-RNORM
        DELTA=DELTA-PI
      ENDIF
      IF(DELTA.GE.0.0D0) THEN
        DELTA=DMOD(DELTA,TPI)
      ELSE
        DELTA=-DMOD(-DELTA,TPI)
      ENDIF
C
C  ****  Coulomb functions.
C
      CALL FCOUL(ETA,RLAMB1,X,FM1,FPM1,GM1,GPM1,ERR)
      IF(ERR.GT.1.0D-6) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERR=1.0D0
        RETURN
      ENDIF
      SLA=(RLAMB/X)+(ETA/RLAMB)
      F=RLAMB*(SLA*FM1-FPM1)/RLA
      G=RLAMB*(SLA*GM1-GPM1)/RLA
C
C  ****  Dirac-Coulomb radial wave functions.
C
      Q2=P1*P2*RNORM
      Q1=RLA*PC*RNORM
      P1=P1*Q1
      Q1=ZETA*Q1
      P2=ZETA*P2*RNORM
C
      FU=P1*F+P2*FM1
      FL=-Q1*F-Q2*FM1
      GU=P1*G+P2*GM1
      GL=-Q1*G-Q2*GM1
      RETURN
C
C  ****  Z=0. Spherical Bessel functions.
C
    1 CONTINUE
      RLAMB=ABS(K)
      CALL FCOUL(0.0D0,RLAMB,X,F,FP,G,GP,ERR)
      DELTA=0.0D0
      IF(ERR.GE.1.0D-6) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERR=1.0D0
        RETURN
      ENDIF
      FM1=(RLAMB*F/X)+FP
      GM1=(RLAMB*G/X)+GP
      FACT=SQRT(E/(E+TSL2))
      IF(K.LT.0) THEN
        FU=FM1
        FL=-FACT*F
        GU=GM1
        GL=-FACT*G
      ELSE
        FU=F
        FL=FACT*FM1
        GU=G
        GL=FACT*GM1
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FCOUL
C  *********************************************************************
      SUBROUTINE FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)
C
C     Calculation of Coulomb functions for real ETA, RLAMB.GT.-1 and X
C  larger than, or of the order of XTP0 (the turning point for RLAMB=0).
C  Steed's continued fraction method is combined with several recursion
C  relations and an asymptotic expansion. The output value ERR=1.0D0
C  indicates that the adopted evaluation algorithm is not applicable
C  (X is too small).
C
C  Input arguments:
C     ETA ...... Sommerfeld's parameter,
C     RLAMB .... angular momentum,
C     X ........ variable (=wave number times radial distance).
C
C  Output arguments:
C     F, FP .... regular function and its derivative,
C     G, GP .... irregular function and its derivative,
C     ERR ...... relative numerical uncertainty. a value of the
C                order of 10**(-N) means that the calculated
C                functions are accurate to N decimal figures.
C                The maximum accuracy attainable with double
C                precision arithmetic is about 1.0D-15.
C
C     Other subprograms required: subroutine SUM2F0 and
C                                 functions DELTAC and CLGAM.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI,
     1  EPS=1.0D-16,TOP=1.0D5,NTERM=1000)
C
      IF(RLAMB.LT.-0.999D0) THEN
        WRITE(6,'(1X,''*** Error in FCOUL: RLAMB.LT.-0.999'')')
        STOP
      ENDIF
      IF(X.LT.EPS) GO TO 5
C
C  ****  Numerical constants.
C
      CI=DCMPLX(0.0D0,1.0D0)
      CI2=2.0D0*CI
      CIETA=CI*ETA
      X2=X*X
      ETA2=ETA*ETA
C
C  ****  Turning point (XTP). (44)
C
      IF(RLAMB.GE.0.0D0) THEN
        XTP=ETA+SQRT(ETA2+RLAMB*(RLAMB+1.0D0))
      ELSE
        XTP=EPS
      ENDIF
      ERRS=10.0D0
      IF(X.LT.XTP) GO TO 1
C
C  ************  Asymptotic expansion. (71-75)
C
C  ****  Coulomb phase-shift.
      DELTA=DELTAC(ETA,RLAMB)
C
      CPA=CIETA-RLAMB
      CPB=CIETA+RLAMB+1.0D0
      CPZ=CI2*X
      CALL SUM2F0(CPA,CPB,CPZ,C2F0,ERR1)
      CQA=CPA+1.0D0
      CQB=CPB+1.0D0
      CALL SUM2F0(CQA,CQB,CPZ,C2F0P,ERR2)
      C2F0P=CI*C2F0P*CPA*CPB/(2.0D0*X2)
C  ****  Functions.
      THETA=X-ETA*LOG(2.0D0*X)-RLAMB*PIH+DELTA
      IF(THETA.GT.1.0D4) THETA=DMOD(THETA,TPI)
      CEITH=CDEXP(CI*THETA)
      CGIF=C2F0*CEITH
      G=CGIF
      F=-CI*CGIF
C  ****  Derivatives.
      CGIFP=(C2F0P+CI*(1.0D0-ETA/X)*C2F0)*CEITH
      GP=CGIFP
      FP=-CI*CGIFP
C  ****  Global uncertainty. The Wronskian may differ from 1 due
C        to truncation and round off errors.
      ERR=MAX(ERR1,ERR2,ABS(G*FP-F*GP-1.0D0))
      IF(ERR.LE.EPS) RETURN
      ERRS=ERR
C
C  ************  Steed's continued fraction method.
C
    1 CONTINUE
      CIETA2=CIETA+CIETA
      ETAX=ETA*X
C
C  ****  Continued fraction for F. (60-70)
C
      INULL=0
      RLAMBN=RLAMB+1.0D0
      A1=-(RLAMBN+1.0D0)*(RLAMBN**2+ETA2)*X/RLAMBN
      B0=(RLAMBN/X)+(ETA/RLAMBN)
      B1=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)
      FA3=B0
      FA2=B0*B1+A1
      FB3=1.0D0
      FB2=B1
      RF=FA3
C
      DO N=2,NTERM
        RFO=RF
        DAF=ABS(RF)
        RLAMBN=RLAMB+N
        AN=-(RLAMBN**2-1.0D0)*(RLAMBN**2+ETA2)*X2
        BN=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)
        FA1=FA2*BN+FA3*AN
        FB1=FB2*BN+FB3*AN
C
        TST=ABS(FB1)
        IF(TST.LT.1.0D-25) THEN
          IF(INULL.GT.0) STOP
          INULL=1
          FA3=FA2
          FA2=FA1
          FB3=FB2
          FB2=FB1
          RF=RFO
        ELSE
          FA3=FA2/TST
          FA2=FA1/TST
          FB3=FB2/TST
          FB2=FB1/TST
          RF=FA2/FB2
          IF(ABS(RF-RFO).LT.EPS*DAF) GO TO 2
        ENDIF
      ENDDO
    2 CONTINUE
      IF(DAF.GT.1.0D-25) THEN
        ERRF=ABS(RF-RFO)/DAF
      ELSE
        ERRF=EPS
      ENDIF
      IF(ERRF.GT.ERRS) THEN
        ERR=ERRS
        IF(ERR.GT.1.0D-6) GO TO 5
        RETURN
      ENDIF
C
C  ****  Downward recursion for F and FP. Only if RLAMB.GT.1 and
C        X.LT.XTP. (48,49)
C
      RLAMB0=RLAMB
      IF(X.GE.XTP.OR.RLAMB0.LT.1.0D0) THEN
        ISHIFT=0
        XTPC=XTP
        RFM=0.0D0
      ELSE
        FT=1.0D0
        FTP=RF
        IS0=RLAMB0+1.0D-6
        TST=X*(X-2.0D0*ETA)
        RL1T=0.0D0
        DO I=1,IS0
          ETARL0=ETA/RLAMB0
          RL=SQRT(1.0D0+ETARL0**2)
          SSL=(RLAMB0/X)+ETARL0
          RLAMB0=RLAMB0-1.0D0
          FTO=FT
          FT=(SSL*FT+FTP)/RL
          FTP=SSL*FT-RL*FTO
          IF(FT.GT.1.0D10) THEN
            FTP=FTP/FT
            FT=1.0D0
          ENDIF
          RL1T=RLAMB0*(RLAMB0+1.0D0)
          IF(TST.GT.RL1T) THEN
            ISHIFT=I
            GO TO 3
          ENDIF
        ENDDO
        ISHIFT=IS0
    3   CONTINUE
        XTPC=ETA+SQRT(ETA2+RL1T)
        RFM=FTP/FT
      ENDIF
C
C  ****  Continued fraction for P+CI*Q with RLAMB0. (76-79)
C
      INULL=0
      CAN=CIETA-ETA2-RLAMB0*(RLAMB0+1.0D0)
      CB0=X-ETA
      CBN=2.0D0*(X-ETA+CI)
      CFA3=CB0
      CFA2=CB0*CBN+CAN
      CFB3=1.0D0
      CFB2=CBN
      CPIQ=CFA3
C
      DO N=2,NTERM
        CPIQO=CPIQ
        DAPIQ=CDABS(CPIQ)
        CAN=CAN+CIETA2+(N+N-2)
        CBN=CBN+CI2
        CFA1=CFA2*CBN+CFA3*CAN
        CFB1=CFB2*CBN+CFB3*CAN
        TST=CDABS(CFB1)
C
        IF(TST.LT.1.0D-25) THEN
          IF(INULL.GT.0) STOP
          INULL=1
          CFA3=CFA2
          CFA2=CFA1
          CFB3=CFB2
          CFB2=CFB1
          CPIQ=CPIQO
        ELSE
          CFA3=CFA2/TST
          CFA2=CFA1/TST
          CFB3=CFB2/TST
          CFB2=CFB1/TST
          CPIQ=CFA2/CFB2
          IF(CDABS(CPIQ-CPIQO).LT.EPS*DAPIQ) GO TO 4
        ENDIF
      ENDDO
    4 CONTINUE
      IF(DAPIQ.GT.1.0D-25) THEN
        ERRPIQ=CDABS(CPIQ-CPIQO)/DAPIQ
      ELSE
        ERRPIQ=EPS
      ENDIF
      IF(ERRPIQ.GT.ERRS) THEN
        ERR=ERRS
        IF(ERR.GT.1.0D-6) GO TO 5
        RETURN
      ENDIF
      CPIQ=CI*CPIQ/X
C
      RP=CPIQ
      RQ=-CI*CPIQ
      IF(RQ.LE.1.0D-25) GO TO 5
      ERR=MAX(ERRF,ERRPIQ)
C
C  ****  Inverting Steed's transformation. (57,58)
C
      IF(ISHIFT.LT.1) THEN
        RFP=RF-RP
        F=SQRT(RQ/(RFP**2+RQ**2))
        IF(FB2.LT.0.0D0) F=-F
        FP=RF*F
        G=RFP*F/RQ
        GP=(RP*RFP-RQ**2)*F/RQ
        IF(X.LT.XTP.AND.G.GT.TOP*F) GO TO 5
      ELSE
        RFP=RFM-RP
        FM=SQRT(RQ/(RFP**2+RQ**2))
        G=RFP*FM/RQ
        GP=(RP*RFP-RQ**2)*FM/RQ
        IF(X.LT.XTPC.AND.G.GT.TOP*FM) GO TO 5
C  ****  Upward recursion for G and GP (if ISHIFT.GT.0). (50,51)
        DO I=1,ISHIFT
          RLAMB0=RLAMB0+1.0D0
          ETARL0=ETA/RLAMB0
          RL=SQRT(1.0D0+ETARL0**2)
          SSL=(RLAMB0/X)+ETARL0
          GO=G
          G=(SSL*GO-GP)/RL
          GP=RL*GO-SSL*G
          IF(G.GT.1.0D35) GO TO 5
        ENDDO
        W=RF*G-GP
        F=1.0D0/W
        FP=RF/W
      ENDIF
C  ****  The Wronskian may differ from 1 due to round off errors.
      ERR=MAX(ERR,ABS(FP*G-F*GP-1.0D0))
      IF(ERR.LT.1.0D-6) RETURN
C
    5 CONTINUE
      F=0.0D0
      FP=0.0D0
      G=1.0D35
      GP=-1.0D35
      ERR=1.0D0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SUM2F0
C  *********************************************************************
      SUBROUTINE SUM2F0(CA,CB,CZ,CF,ERR)
C
C     Summation of the 2F0(CA,CB;CS) hypergeometric asymptotic series.
C  The positive and negative contributions to the real and imaginary
C  parts are added separately to obtain an estimate of rounding errors.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-16,ACCUR=0.5D-15,NTERM=75)
      RRP=1.0D0
      RRN=0.0D0
      RIP=0.0D0
      RIN=0.0D0
      CDF=1.0D0
      ERR2=0.0D0
      ERR3=1.0D0
      AR=0.0D0
      AF=0.0D0
      DO I=1,NTERM
        J=I-1
        CDF=CDF*(CA+J)*(CB+J)/(I*CZ)
        ERR1=ERR2
        ERR2=ERR3
        ERR3=CDABS(CDF)
        IF(ERR1.GT.ERR2.AND.ERR2.LT.ERR3) GO TO 1
        AR=CDF
        IF(AR.GT.0.0D0) THEN
          RRP=RRP+AR
        ELSE
          RRN=RRN+AR
        ENDIF
        AI=DCMPLX(0.0D0,-1.0D0)*CDF
        IF(AI.GT.0.0D0) THEN
          RIP=RIP+AI
        ELSE
          RIN=RIN+AI
        ENDIF
        CF=DCMPLX(RRP+RRN,RIP+RIN)
        AF=CDABS(CF)
        IF(AF.GT.1.0D25) THEN
          CF=0.0D0
          ERR=1.0D0
          RETURN
        ENDIF
        IF(ERR3.LT.1.0D-25*AF.OR.ERR3.LT.EPS) THEN
           ERR=EPS
           RETURN
        ENDIF
      ENDDO
C  ****  Round off error.
    1 CONTINUE
      TR=ABS(RRP+RRN)
      IF(TR.GT.1.0D-25) THEN
        ERRR=(RRP-RRN)*ACCUR/TR
      ELSE
        ERRR=1.0D0
      ENDIF
      TI=ABS(RIP+RIN)
      IF(TI.GT.1.0D-25) THEN
        ERRI=(RIP-RIN)*ACCUR/TI
      ELSE
        ERRI=1.0D0
      ENDIF
C  ****  ... and truncation error.
      IF(AR.GT.1.0D-25) THEN
        ERR=MAX(ERRR,ERRI)+ERR2/AF
      ELSE
        ERR=MAX(ERRR,ERRI)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION DELTAC
C  *********************************************************************
      FUNCTION DELTAC(ETA,RLAMB)
C
C     Calculation of Coulomb phase shift (modulus 2*PI). (47)
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)
      CI=DCMPLX(0.0D0,1.0D0)
C  ****  Coulomb phase-shift.
      DELTAC=-CI*CLGAM(RLAMB+1.0D0+CI*ETA)
      IF(DELTAC.GE.0.0D0) THEN
        DELTAC=DMOD(DELTAC,TPI)
      ELSE
        DELTAC=-DMOD(-DELTAC,TPI)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION CLGAM
C  *********************************************************************
      FUNCTION CLGAM(CZ)
C
C     This function gives LOG(GAMMA(CZ)) for complex arguments.
C
C   Ref.: M. Abramowitz and I.A. Stegun, 'Handbook of Mathematical
C         Functions'. Dover, New York (1974). PP 255-257.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      CZA=CZ
      ICONJ=0
      AR=CZA
      CLGAM=36.84136149D0
      IF(CDABS(CZA).LT.1.0D-16) RETURN
C
      AI=CZA*DCMPLX(0.0D0,-1.0D0)
      IF(AI.GT.0.0D0) THEN
        ICONJ=0
      ELSE
        ICONJ=1
        CZA=DCONJG(CZA)
      ENDIF
C
      CZFAC=1.0D0
      CZFL=0.0D0
    1 CONTINUE
      CZFAC=CZFAC/CZA
      IF(CDABS(CZFAC).GT.1.0D8) THEN
        CZFL=CZFL+CDLOG(CZFAC)
        CZFAC=1.0D0
      ENDIF
      CZA=CZA+1.0D0
      AR=CZA
      IF(CDABS(CZA).LT.1.0D-16) RETURN
      IF(CDABS(CZA).GT.15.0D0.AND.AR.GT.0.0D0) GO TO 2
      GO TO 1
C  ****  Stirling's expansion of CDLOG(GAMMA(CZA)).
    2 CONTINUE
      CZI2=1.0D0/(CZA*CZA)
      CZS=(43867.0D0/244188.0D0)*CZI2
      CZS=(CZS-3617.0D0/122400.0D0)*CZI2
      CZS=(CZS+1.0D0/156.0D0)*CZI2
      CZS=(CZS-691.0D0/360360.0D0)*CZI2
      CZS=(CZS+1.0D0/1188.0D0)*CZI2
      CZS=(CZS-1.0D0/1680.0D0)*CZI2
      CZS=(CZS+1.0D0/1260.0D0)*CZI2
      CZS=(CZS-1.0D0/360.0D0)*CZI2
      CZS=(CZS+1.0D0/12.0D0)/CZA
      CLGAM=(CZA-0.5D0)*CDLOG(CZA)-CZA+9.1893853320467274D-1+CZS
     1     +CZFL+CDLOG(CZFAC)
      IF(ICONJ.EQ.1) CLGAM=DCONJG(CLGAM)
      RETURN
      END
C  *********************************************************************
C                         FUNCION BESJN
C  *********************************************************************
      FUNCTION BESJN(JY,N,X)
C
C     This function computes the spherical Bessel functions of the first
C  kind and spherical Bessel functions of the second kind (also known as
C  spherical Neumann functions) for real positive arguments.
C
C  Input arguments:
C        JY ...... kind: 1(Bessel) or 2(Neumann).
C        N ....... order (integer).
C        X ....... argument (real and positive).
C
C  Ref.: M. Abramowitz and I.A. Stegun, 'Handbook of Mathematical
C        Functions'. Dover, New York (1974). pp 435-478.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      IF(X.LT.0) THEN
        WRITE(6,1000)
 1000   FORMAT(1X,'*** Negative argument in function BESJN.')
        STOP
      ENDIF
C  ****  Order and phase correction for Neumann functions.
C        Abramowitz and Stegun, eq. 10.1.15.
      IF(JY.EQ.2) THEN
        NL=-N-1
        IPH=2*MOD(ABS(N),2)-1
      ELSE
        NL=N
        IPH=1
      ENDIF
C  ****  Selection of calculation mode.
      IF(NL.LT.0) GO TO 5
      IF(X.GT.1.0D0*NL) GO TO 3
      XI=X*X
      IF(XI.GT.NL+NL+3.0D0) GO TO 2
C  ****  Power series for small arguments and positive orders.
C        Abramowitz and Stegun, eq. 10.1.2.
      F1=1.0D0
      IP=1
      IF(NL.NE.0) THEN
        DO I=1,NL
          IP=IP+2
          F1=F1*X/IP
        ENDDO
      ENDIF
      XI=0.5D0*XI
      BESJN=1.0D0
      PS=1.0D0
      DO I=1,1000
        IP=IP+2
        PS=-PS*XI/(I*IP)
        BESJN=BESJN+PS
        IF(ABS(PS).LT.1.0D-18*ABS(BESJN)) GO TO 1
      ENDDO
    1 BESJN=IPH*F1*BESJN
      RETURN
C  ****  Miller's method for positive orders and intermediate arguments.
C        Abramowitz and Stegun, eq. 10.1.19.
    2 XI=1.0D0/X
      F2=0.0D0
      F3=1.0D-35
      IP=2*(NL+31)+3
      DO I=1,31
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D30) THEN
          F2=F2/F3
          F3=1.0D0
        ENDIF
      ENDDO
      BESJN=1.0D0
      F2=F2/F3
      F3=1.0D0
      DO I=1,NL
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D30) THEN
          BESJN=BESJN/F3
          F2=F2/F3
          F3=1.0D0
        ENDIF
      ENDDO
      BESJN=IPH*XI*SIN(X)*BESJN/F3
      RETURN
C  ****  Recurrence relation for arguments greater than order.
C        Abramowitz and Stegun, eq. 10.1.19.
    3 XI=1.0D0/X
      F3=XI*SIN(X)
      IF(NL.EQ.0) GO TO 4
      F2=F3
      F3=XI*(F2-COS(X))
      IF(NL.EQ.1) GO TO 4
      IP=1
      DO I=2,NL
        F1=F2
        F2=F3
        IP=IP+2
        F3=IP*XI*F2-F1
      ENDDO
    4 BESJN=IPH*F3
      RETURN
C  ****  Recurrence relation for negative orders.
C        Abramowitz and Stegun, eq. 10.1.19.
    5 NL=ABS(NL)
      IF(X.LT.7.36D-1*(NL+1)*1.0D-35**(1.0D0/(NL+1))) THEN
        BESJN=-1.0D35
        RETURN
      ENDIF
      XI=1.0D0/X
      F3=XI*SIN(X)
      F2=XI*(F3-COS(X))
      IP=3
      DO I=1,NL
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D35) THEN
          BESJN=-1.0D35
          RETURN
        ENDIF
      ENDDO
      BESJN=IPH*F3
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC          Cubic spline interpolation          CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE SPLINE
C  *********************************************************************
      SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)
C
C  Cubic spline interpolation of tabulated data.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     Y(I) (I=1:N) ... corresponding function values.
C     S1,SN .......... second derivatives at X(1) and X(N). The natural
C                      spline corresponds to taking S1=SN=0.
C     N .............. number of grid points.
C  Output:
C     A(1:N),B(1:N),C(1:N),D(1:N) ... spline coefficients.
C
C  The interpolating cubic polynomial in the I-th interval, from X(I) to
C  X(I+1), is
C               P(x) = A(I)+x*(B(I)+x*(C(I)+x*D(I)))
C
C  Reference: M.J. Maron, 'Numerical Analysis: a Practical Approach',
C             MacMillan Publ. Co., New York, 1982.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      IF(N.LT.4) THEN
        WRITE(6,10) N
   10   FORMAT(5X,'Spline interpolation cannot be performed with',
     1    I4,' points. Stop.')
        STOP 'SPLINE. N is less than 4.'
      ENDIF
      N1=N-1
      N2=N-2
C  ****  Auxiliary arrays H(=A) and DELTA(=D).
      DO I=1,N1
        IF(X(I+1)-X(I).LT.1.0D-13) THEN
          WRITE(6,11)
   11     FORMAT(5X,'Spline X values not in increasing order. Stop.')
          STOP 'SPLINE. X values not in increasing order.'
        ENDIF
        A(I)=X(I+1)-X(I)
        D(I)=(Y(I+1)-Y(I))/A(I)
      ENDDO
C  ****  Symmetric coefficient matrix (augmented).
      DO I=1,N2
        B(I)=2.0D0*(A(I)+A(I+1))
        K=N1-I+1
        D(K)=6.0D0*(D(K)-D(K-1))
      ENDDO
      D(2)=D(2)-A(1)*S1
      D(N1)=D(N1)-A(N1)*SN
C  ****  Gauss solution of the tridiagonal system.
      DO I=2,N2
        R=A(I)/B(I-1)
        B(I)=B(I)-R*A(I)
        D(I+1)=D(I+1)-R*D(I)
      ENDDO
C  ****  The SIGMA coefficients are stored in array D.
      D(N1)=D(N1)/B(N2)
      DO I=2,N2
        K=N1-I+1
        D(K)=(D(K)-A(K)*D(K+1))/B(K-1)
      ENDDO
      D(N)=SN
C  ****  Spline coefficients.
      SI1=S1
      DO I=1,N1
        SI=SI1
        SI1=D(I+1)
        H=A(I)
        HI=1.0D0/H
        A(I)=(HI/6.0D0)*(SI*X(I+1)**3-SI1*X(I)**3)
     1      +HI*(Y(I)*X(I+1)-Y(I+1)*X(I))
     2      +(H/6.0D0)*(SI1*X(I)-SI*X(I+1))
        B(I)=(HI/2.0D0)*(SI1*X(I)**2-SI*X(I+1)**2)
     1      +HI*(Y(I+1)-Y(I))+(H/6.0D0)*(SI-SI1)
        C(I)=(HI/2.0D0)*(SI*X(I+1)-SI1*X(I))
        D(I)=(HI/6.0D0)*(SI1-SI)
      ENDDO
C  ****  Natural cubic spline for X.GT.X(N).
      FN=Y(N)
      FNP=B(N1)+X(N)*(2.0D0*C(N1)+X(N)*3.0D0*D(N1))
      A(N)=FN-X(N)*FNP
      B(N)=FNP
      C(N)=0.0D0
      D(N)=0.0D0
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FINDI
C  *********************************************************************
      SUBROUTINE FINDI(X,XC,N,I)
C
C  Finds the interval (X(I),X(I+1)) that contains the value XC.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     XC ............. point to be located.
C     N  ............. number of grid points.
C  Output:
C     I .............. interval index.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N)
C
      IF(XC.GT.X(N)) THEN
        I=N
        RETURN
      ENDIF
      IF(XC.LT.X(1)) THEN
        I=1
        RETURN
      ENDIF
      I=1
      I1=N
    1 IT=(I+I1)/2
      IF(XC.GT.X(IT)) THEN
        I=IT
      ELSE
        I1=IT
      ENDIF
      IF(I1-I.GT.1) GO TO 1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE INTEG
C  *********************************************************************
      SUBROUTINE INTEG(X,A,B,C,D,XL,XU,SUM,N)
C
C  Computes the integral of a cubic spline function.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     A(1:N),B(1:N),C(1:N),D(1:N) ... spline coefficients.
C     N  ............. number of grid points.
C     XL ............. lower limit of the integral.
C     XU ............. upper limit of the integral.
C  Output:
C     SUM ............ value of the integral.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),A(N),B(N),C(N),D(N)
C  ****  Set integration limits in increasing order.
      IF(XU.GT.XL) THEN
        XLL=XL
        XUU=XU
        SIGN=1.0D0
      ELSE
        XLL=XU
        XUU=XL
        SIGN=-1.0D0
      ENDIF
C  ****  Check integral limits.
      IF(XLL.LT.X(1).OR.XUU.GT.X(N)) THEN
        WRITE(6,10)
   10   FORMAT(5X,'Integral limits out of range. Stop.')
C       STOP 'INTEG. Integral limits out of range.'
      ENDIF
C  ****  Find involved intervals.
      SUM=0.0D0
      CALL FINDI(X,XLL,N,IL)
      CALL FINDI(X,XUU,N,IU)
C
      IF(IL.EQ.IU) THEN
C  ****  Only a single interval involved.
        X1=XLL
        X2=XUU
        SUM=X2*(A(IL)+X2*((B(IL)/2)+X2*((C(IL)/3)+X2*D(IL)/4)))
     1     -X1*(A(IL)+X1*((B(IL)/2)+X1*((C(IL)/3)+X1*D(IL)/4)))
      ELSE
C  ****  Contributions from several intervals.
        X1=XLL
        X2=X(IL+1)
        SUM=X2*(A(IL)+X2*((B(IL)/2)+X2*((C(IL)/3)+X2*D(IL)/4)))
     1     -X1*(A(IL)+X1*((B(IL)/2)+X1*((C(IL)/3)+X1*D(IL)/4)))
        IL=IL+1
        DO I=IL,IU
          X1=X(I)
          X2=X(I+1)
          IF(I.EQ.IU) X2=XUU
          SUMP=X2*(A(I)+X2*((B(I)/2)+X2*((C(I)/3)+X2*D(I)/4)))
     1        -X1*(A(I)+X1*((B(I)/2)+X1*((C(I)/3)+X1*D(I)/4)))
          SUM=SUM+SUMP
        ENDDO
      ENDIF
      SUM=SIGN*SUM
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE INTEG2
C  *********************************************************************
      SUBROUTINE INTEG2(X,A,B,C,D,XL,XU,SUM,N)
C
C  Computes the integral of a squared cubic spline function.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     A(1:N),B(1:N),C(1:N),D(1:N) ... spline coefficients.
C     N  ............. number of grid points.
C     XL ............. lower limit of the integral.
C     XU ............. upper limit of the integral.
C  Output:
C     SUM ............ value of the integral.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),A(N),B(N),C(N),D(N)
C  ****   Set integration limits in increasing order.
      IF(XU.GT.XL) THEN
        XLL=XL
        XUU=XU
        SIGN=1.0D0
      ELSE
        XLL=XU
        XUU=XL
        SIGN=-1.0D0
      ENDIF
C  ****  Check integral limits.
      IF(XLL.LT.X(1).OR.XUU.GT.X(N)) THEN
        WRITE(6,10)
   10   FORMAT(5X,'Integral limits out of range. Stop.')
C       STOP 'INTEG. Integral limits out of range.'
      ENDIF
C  ****  Find involved intervals.
      SUM=0.0D0
      CALL FINDI(X,XLL,N,IL)
      CALL FINDI(X,XUU,N,IU)
C
      IF(IL.EQ.IU) THEN
C  ****  Only a single interval involved.
        X1=XLL
        X2=XUU
        SUM=X2*(A(IL)*(A(IL)+X2*B(IL))
     1     +X2*X2*(((2*A(IL)*C(IL)+B(IL)**2)/3)
     2     +X2*(((B(IL)*C(IL)+A(IL)*D(IL))/2)
     3     +X2*(((2*B(IL)*D(IL)+C(IL)**2)/5)
     4     +X2*D(IL)*((C(IL)/3)+X2*D(IL)/7)))))
     5     -X1*(A(IL)*(A(IL)+X1*B(IL))
     6     +X1*X1*(((2*A(IL)*C(IL)+B(IL)**2)/3)
     7     +X1*(((B(IL)*C(IL)+A(IL)*D(IL))/2)
     8     +X1*(((2*B(IL)*D(IL)+C(IL)**2)/5)
     9     +X1*D(IL)*((C(IL)/3)+X1*D(IL)/7)))))
        IF(SUM.LT.0.0D0) SUM=0.0D0
      ELSE
C  ****  Contributions from several intervals.
        X1=XLL
        X2=X(IL+1)
        SUM=X2*(A(IL)*(A(IL)+X2*B(IL))
     1     +X2*X2*(((2*A(IL)*C(IL)+B(IL)**2)/3)
     2     +X2*(((B(IL)*C(IL)+A(IL)*D(IL))/2)
     3     +X2*(((2*B(IL)*D(IL)+C(IL)**2)/5)
     4     +X2*D(IL)*((C(IL)/3)+X2*D(IL)/7)))))
     5     -X1*(A(IL)*(A(IL)+X1*B(IL))
     6     +X1*X1*(((2*A(IL)*C(IL)+B(IL)**2)/3)
     7     +X1*(((B(IL)*C(IL)+A(IL)*D(IL))/2)
     8     +X1*(((2*B(IL)*D(IL)+C(IL)**2)/5)
     9     +X1*D(IL)*((C(IL)/3)+X1*D(IL)/7)))))
        IF(SUM.LT.0.0D0) SUM=0.0D0
        IL=IL+1
        DO I=IL,IU
          X1=X(I)
          X2=X(I+1)
          IF(I.EQ.IU) X2=XUU
          SUMP=X2*(A(I)*(A(I)+X2*B(I))
     1        +X2*X2*(((2*A(I)*C(I)+B(I)**2)/3)
     2        +X2*(((B(I)*C(I)+A(I)*D(I))/2)
     3        +X2*(((2*B(I)*D(I)+C(I)**2)/5)
     4        +X2*D(I)*((C(I)/3)+X2*D(I)/7)))))
     5        -X1*(A(I)*(A(I)+X1*B(I))
     6        +X1*X1*(((2*A(I)*C(I)+B(I)**2)/3)
     7        +X1*(((B(I)*C(I)+A(I)*D(I))/2)
     8        +X1*(((2*B(I)*D(I)+C(I)**2)/5)
     9        +X1*D(I)*((C(I)/3)+X1*D(I)/7)))))
          IF(SUMP.LT.0.0D0) SUMP=0.0D0
          SUM=SUM+SUMP
        ENDDO
      ENDIF
      SUM=SIGN*SUM
      RETURN
      END
C  *********************************************************************
C                        SUBROUTINE ERRSPL
C  *********************************************************************
      SUBROUTINE ERRSPL(ERR,X,Y,N)
C
C     This subroutine estimates the error introduced by the natural
C  cubic spline interpolation in a table X(I),Y(I) (I=1,...,N). The
C  interpolation error in the vicinity of X(K) is approximated by the
C  difference between Y(K) and the value obtained from the spline that
C  interpolates the table with the K-th point removed. ERR is the
C  largest relative error along the table.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NDIM=10000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/STORE/F(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      DIMENSION X(NDIM),Y(NDIM),R(NDIM)
C
      OPEN(13,FILE='errspl.out')
      WRITE(13,1001)
 1001 FORMAT(1X,'# test of the spline interpolation (output fro',
     1  'm subroutine ERRSPL)',/1X,'#',8X,'x',16X,'y(x)',12X,
     1  'y_int(x)',7X,'rel.dif.')
      ERR=0.0D0
      N1=N-1
      DO I=2,N1
        DO J=1,N1
          IF(J.LT.I) THEN
            R(J)=X(J)
            F(J)=Y(J)
          ELSE
            R(J)=X(J+1)
            F(J)=Y(J+1)
          ENDIF
        ENDDO
        CALL SPLINE(R,F,A,B,C,D,0.0D0,0.0D0,N1)
        RC=X(I)
        YI=A(I-1)+RC*(B(I-1)+RC*(C(I-1)+RC*D(I-1)))
        IF(ABS(Y(I)).GT.1.0D-3) THEN
          ERRP=1.0D0-YI/Y(I)
        ELSE
          ERRP=YI-Y(I)
        ENDIF
        ERR=MAX(ERR,ABS(ERRP))
        WRITE(13,1002) X(I),Y(I),YI,ERRP
 1002   FORMAT(1X,1P,3E18.10,E12.4)
      ENDDO
      CLOSE(13)
      RETURN
      END
C  *********************************************************************
C                        SUBROUTINE SGRID
C  *********************************************************************
      SUBROUTINE SGRID(R,DR,RN,R2,DRN,N,NMAX)
C
C  This subroutine sets up a radial grid R(I) (I=1:N) such that
C    1) A*R(I)+B*DLOG(R(I)+C)+D=I.
C    2) R(1)=0, R(N)=RN, R(2)=R2 and R(N-1)=RN-DRN (approx.).
C    3) The spacing between consecutive grid points, R(I+1)-R(I),
C       increases with I and is always less than DRN.
C
C  Input arguments:
C    RN ..... outer grid point (the grid extends from 0 up to RN).
C             RN must be greater than 1.0E-5
C    R2  .... R(2) (controls the grid spacing at small radii).
C             R2 must be less than 1.0E-2 and less than RN.
C    DRN .... R(N)-R(N-1) (controls the grid spacing at large radial
C             distances).
C    N ...... tentative number of grid points (it may be increased to
C             meet conditions 2 and 3).
C    NMAX ... physical dimensions of the arrays R(.) and DR(.); N cannot
C             exceed NMAX.
C
C  Output arguments:
C    N ...... number of grid points. It may be greater than the input
C             value.
C    R(1:N) ... radial grid points.
C    DR(1:N) ... values of the derivative of R with respect to I, which
C             are required to evaluate integrals using, e.g., Simpson's
C             method.
C
C     To describe the radial wave functions of bound electrons of an
C  atom or positive ion in its ground state configuration, the following
C  values of the input parameters should be adequate: RN of the order of
C  50, R2 about 1.0E-5 or smaller, DRN about 0.5, N=750 or larger.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION R(NMAX),DR(NMAX)
C  ****  Set IWR=1 to print the output radial grid.
      IWR=0
      IER=0
C
 1000 FORMAT(' RN=',1P,E13.6,', R2=',E13.6,', DRN=',E13.6,', N=',I5)
      RRN=RN
      IF(RRN.LT.1.0D-5) THEN
        WRITE(6,'('' *** Error (SGRID): RN is .LT. 1.E-5.'')')
        WRITE(6,1000) RRN,R2,DRN,N
        STOP 'RN is too small.'
      ENDIF
C
      RR2=R2
      IF(RR2.LT.1.0D-7) RR2=1.0D-7
      IF(RR2.GE.1.0D-2) THEN
        WRITE(6,1000) RRN,RR2,DRN,N
        WRITE(6,'('' *** Warning (SGRID): R2 is .GE. 1.E-2.'')')
        RR2=1.0D-2
        WRITE(6,1000) RRN,RR2,DRN,N
      ENDIF
C
      IF(RR2.GE.RRN) THEN
        WRITE(6,'('' *** Error (SGRID): R2 is .GE. RN.'')')
        WRITE(6,1000) RRN,RR2,DRN,N
        STOP 'R2 is greater than RN.'
      ENDIF
C
      RDRN=DRN
      IF(RDRN.LE.RR2) THEN
        WRITE(6,1000) RRN,RR2,RDRN,N
        WRITE(6,'('' *** Warning (SGRID): DRN is .LE. R2.'')')
        RDRN=MAX(2.0D0*(RRN-RR2)/MAX(N,10),RR2)
        WRITE(6,1000) RRN,RR2,RDRN,N
      ENDIF
C
      NR=MAX(N,10)
      TST=5.0D0*(RRN-RR2)/NR
      IF(RDRN.GT.TST) RDRN=TST
      NLOW=(RRN/RDRN)+10.0D0
      IF(NR.LT.NLOW) THEN
        WRITE(6,1000) RRN,RR2,RDRN,NR
        WRITE(6,'('' *** Warning (SGRID): NR is .LT. NLOW.'')')
        WRITE(6,'('' NLOW ='',I8)') NLOW
        NR=NLOW
        WRITE(6,1000) RRN,RR2,RDRN,NR
      ENDIF
C
      IF(NR.GT.NMAX) THEN
        WRITE(6,'('' *** Error (SGRID): NR is .GT. NMAX.'')')
        WRITE(6,1000) RRN,RR2,RDRN,NR
        WRITE(6,'('' NMAX ='',I8)') NMAX
        STOP 'NMAX is too small.'
      ENDIF
C
      NHIGH=0.5D0*((RRN/RR2)+(RRN/RDRN))-20.0D0
      IF(NR.GT.NHIGH) THEN
        A=(NR-1)/RRN
        AA=1.0D0/A
        IF(AA.LT.RDRN.AND.AA.LT.RR2) THEN
          B=0.0D0  ! Linear grid.
          C=1.0D0
          D=0.0D0
          GO TO 3
        ENDIF
        WRITE(6,'('' *** Error (SGRID): NR is .GT. NHIGH.'')')
        WRITE(6,'('' NHIGH ='',I8)') NHIGH
        WRITE(6,1000) RRN,RR2,RDRN,NR
        STOP 'Inconsistent grid definition.'
      ENDIF
C
C  ****  Grid parameters (189).
C
      AG=((NR-1)*(RR2/RRN)-1.0D0)*RDRN/(RDRN-RR2)
      IF(AG.LT.-1.0D0.OR.AG.GT.-0.500167D0) THEN
        A=(NR-1)/RRN
        AA=1.0D0/A
        IF(AA.LT.RDRN.AND.AA.LT.RR2) THEN
          B=0.0D0  ! Linear grid.
          C=1.0D0
          D=0.0D0
          GO TO 3
        ENDIF
        WRITE(6,1000) RRN,RR2,RDRN,NR
        WRITE(6,'('' *** Error (SGRID): AG is out of range.'')')
        WRITE(6,'('' AG ='',1P,E13.6)') AG
        STOP 'Inconsistent grid definition.'
      ENDIF
C
      XL=0.0D0
      XU=1000.0D0
    1 X=0.5D0*(XL+XU)
      F=(1.0D0+X)*(X*LOG((1.0D0+X)/X)-1.0D0)
      IF(F.GT.AG) THEN
        XU=X
      ELSE
        XL=X
      ENDIF
      IF(XU-XL.GT.1.0D-15) GO TO 1
C
      C=X*RRN
      B=(C/RRN)*(C+RRN)*(RDRN-RR2)/(RDRN*RR2)
      A=(C-RR2*B)/(C*RR2)
      D=1.0D0-B*LOG(C)
C
      R(1)=0.0D0
      DR(1)=(R(1)+C)/(A*(R(1)+C)+B)
      RR=1.0D-35
      DO I=2,NR
        RL=RR
        RU=RRN
    2   RR=0.5D0*(RU+RL)
        FR=A*RR+B*LOG(RR+C)+D-DBLE(I)
        IF(FR.GT.0.0D0) THEN
          RU=RR
        ELSE
          RL=RR
        ENDIF
        IF(RU-RL.GT.1.0D-15*RR) GO TO 2
        R(I)=RR
        DR(I)=(RR+C)/(A*(RR+C)+B)
        IF(DR(I).LT.DR(I-1)) THEN
          WRITE(6,'('' *** Error (SGRID): The grid spacing does'',
     1      '' not increase with I.'')')
          IER=1
          IWR=1
        ENDIF
      ENDDO
      GO TO 4
C
    3 CONTINUE
      R(1)=0.0D0
      DR(1)=AA
      DO I=2,NR
        R(I)=(I-1)*AA
        DR(I)=AA
      ENDDO
C
C  ****  Print the grid in a file.
C
    4 CONTINUE
      N=NR
      IF(IWR.EQ.1) THEN
        OPEN(17,FILE='sgrid.dat')
        WRITE(17,2001)
 2001   FORMAT('#  Radial grid:  A*R(I)+B*LOG(R(I)+C)+D=I',/'#')
        WRITE(17,2002) RRN,N,R(2),RDRN
 2002   FORMAT('#  ',1P,'RN =',E13.6,',    N =',I5,
     1    /'#  ','R2 =',E13.6,',  DRN =',E13.6,/'#')
        WRITE(17,2003) A,B,C,D
 2003   FORMAT('#  ',1P,' A =',E13.6,',    B =',E13.6,
     1    /'#  ',' C =',E13.6,',    D =',E13.6,/'#')
        WRITE(17,2004)
 2004   FORMAT('#   I',7X,'R(I)',10X,'DR(I)',/'# ',33('-'))
        DO I=1,N
          WRITE(17,'(1X,I5,1P,2E14.6)') I,R(I),DR(I)
        ENDDO
        IF(IER.EQ.1) WRITE(6,2005)
 2005   FORMAT(/'#   The grid spacing decreases with I. STOP.',/)
        CLOSE(17)
      ENDIF
      IF(IER.EQ.1) STOP 'The grid spacing decreases with I.'
C
      RETURN
      END
