SUBROUTINE hurr_in
  USE hurr_vars
  IMPLICIT NONE

  PDS    = 1013.
  FC     = 5.
  CD     = 1.1
  CECD   = 1.
  CD1    = 4.
  XVL    = 200.
  TAUR   = 12.
  RADMAX = 2.
  VTERM  = 7.
  TMIN   = 133.

  ! initial conditions
  RMAX   = 100.
  RO     = 600.
  VMAX   = 12.
  TMID   = 0.
  RSST   = 100.

  ! numerical parameters
  DT     = 20.
  NS     = 5
  EPS    = 0.1
  ALPHA  = 0.01
  RB     = 1500.
  ZB     = 25.
  NSPONGE= 5
  ETIME  = 10.
  TAVE   = 2.
  PLTIME = 2.5
  TIMMAX = 5.
  TIMEPL = 0.1
  ROG    = 200.
  ZOG    = 20.
  
  RETURN
END SUBROUTINE hurr_in

! ===================================================

SUBROUTINE double_timer
  USE hurr_vars
  IMPLICIT NONE

!  CALL param_setup
  ! THIS SHOULD BE DONE AT THE COMPLETION OF 
  ! THE FIRST INTEGRATION IN FWD MODEL
  DTL = 2. * DT
  NS  = 2  * NS
  DTL2 = .5 * DTL
  DO I = 1, M
     DTLR(I)  = 0.5 * DTL*RDR/RS(I)
  END DO
  DO J = 1 , N
     DTLZ(J)  = 0.5 * DTL*RDZ/RHOT(J)
  END DO

  RETURN
END SUBROUTINE double_timer

SUBROUTINE param_setup
  USE hurr_vars
  IMPLICIT NONE

  DZTEMP=ZB*1000./FLOAT(N)
  CALL INTERPOLATE(N,PDS,DZTEMP,TB,QVB,TBS,VMAXTH)
  
  TBS=(TBS+273)*(1000/PDS)**0.286
  TBT=2.0*TB(N)-TB(N-1)
  
  F=FC*1.0E-5
  CD=CD*0.001
  CD1=CD1*1.0E-5
  CE=CD*CECD
  RADT=1.0/(3600.0*TAUR)
  RADMAX=2.*DT*RADMAX/(3600.*24.0)
  RMAX=RMAX*1000.0
  RO=RO*1000.0
  RSST=RSST*1000.0
  RB=1000.0*RB
  ZB=1000.0*ZB
  ISTOP=ETIME*24.0*3600.0/DT
  NZSP=N-NSPONGE+1
  ISTART=1
  DR  = RB / FLOAT(M)
  DZ  = ZB / FLOAT(N)
  DTS = DT / FLOAT(NS)
  DTL = DT
  DZ2 = DZ * DZ
  DR2 = DR * DR
  RDR = 1. / DR
  RDZ = 1. / DZ
  RDR2= 1. / DR2
  RDZ2= 1. / DZ2
  PI  = 4.0 * ATAN(1.0)
  C2  = 90000.0
  CP  = 1005.0
  CV  = 718.0
  RD  = 287.0
  XKAPPA= RD / CP
  XLDCP = 2500.0
  G     = 9.81
  A1    = 7.5 * LOG(10.0)
  PNS   = (PDS/1000.) ** XKAPPA
  TAMB  = TBS*PNS
  ESSS   =  6.11 * EXP(A1* (TAMB-273.0)/(TAMB-36.0))
  QSP   = 0.622 * ESSS/(PDS-ESSS)
  QVBS  = QVB(1)
  TVBS  = TBS * (1.0 + 0.61 * QVBS)
  CC1   = 0.61 * G * DZ /(2.0* CP * TVBS)/(1.0 + 0.61 *QVBS)
  CC2   = 1.0 - G * DZ / (2.0* CP * TVBS * PNS)
  CC5   = -TBS / PNS
  CC3   =  G * DZ / CP
  XHL2  = 0.04 * DR2
  XVL2  = XVL * XVL
  TIME  = DT * (FLOAT(ISTART) - 1.0)
  
  !
  !        DEFINE GRID ARRAYS
  !
  DO I = 1, MP1
     R(I) = (FLOAT(I) - 1.0) * DR
     IF(I /= MP1)  RS(I) = (FLOAT(I) - 0.5) * DR
  END DO
  DO J = 1, NP1
     Z(J) = (FLOAT(J) - 1.0) * DZ
     IF( J /= NP1)  ZS(J) = (FLOAT(J) - 0.5) * DZ
  END DO
  
  !
  !        BASE STATE (HYDROSTATIC EQUATION FOR PN)
  !
  PN(1) = PNS - 0.5*CC3/(0.5*(TB(1)+TBS)*(1.0+0.61*.5*(QVB(1)+QVBS)) )
  PD(1) = 1000.0 * PN(1)**(1.0/XKAPPA)
  DO J = 2, N
     PN(J) = PN(J-1) - CC3/(0.5*(TB(J)+TB(J-1)) &
          * (1.0 + 0.61*0.5*(QVB(J)+QVB(J-1))) )
     TBAR=0.5*(PN(J)*TB(J)+PN(J-1)*TB(J-1))
     IF(TBAR <= TMIN) THEN
        PN(J)=PN(J-1)*EXP(-CC3/TMIN)
        TB(J)=TMIN/PN(J)
     END IF
     PD(J) = 1000.0 * PN(J)**(1.0/XKAPPA)
  END DO
  TBAR= 0.5*(PN(N-1)*TB(N-1)+PN(N)*TB(N))
  PNT = PN(N) -0.5*CC3/(0.5*(TB(N)+TBT)*(1.0+0.61*0.5*(QVB(N)+QVBT)))
  IF(TBAR <= TMIN)THEN
     PNT=PN(N)*EXP(-0.5*CC3/TMIN)
  END IF
  PDT = 1000.0 * PNT ** ( 1.0 / XKAPPA )
  
  DO J = N, 2, -1
     ES9 = 6.11 * EXP( A1* ( PN(J)*TB(J) - 273.0)/( PN(J)*TB(J) - 36.0))
     IF(( PD(J)-ES9 ) > ES9) THEN
        QS(J) = 0.622*ES9/(PD(J)-ES9)
     ELSE
        QS(J) = 0.622
     END IF
  END DO
  
  !
  !       ARRAYS FOR SMALL TIME STEP AND BUOYANCY CALCULATION
  !
  DO J = 1 , N
     RC2(J)    = RDR*DTS* C2/(CP * TB(J)* (1.0 + 0.61 * QVB(J)))
     CPTDR(J)  = DTS*RDR * CP * TB(J)* (1.0 + 0.61 * QVB(J))
     RTB(J)    = 0.5  * G / TB(J)
     RQVB(J)   = 0.61 * G * 0.5 / (1.0 + 0.61 * QVB(J))
     RHOTVB(J) = (1000.0/RD) * PN(J) ** (CV /RD)
     RHOT(J)   = 100.0 * RHOTVB(J) / (TB(J) * (1.0 + 0.61*QVB(J)))
     ZC2(J)    = RDZ*RC2(J)/RHOTVB(J)/RDR
     IF( J==1 )  CYCLE
     PNW = 0.5 * (PN(J) + PN(J-1))
     TVW = 0.5 * (TB(J) + TB(J-1))*(1.0 + 0.305*(QVB(J)+QVB(J-1)))
     RHOTVW(J)=  (1000.0/RD) * PNW ** (CV /RD)
     RHOW(J)  = 100.0 * RHOTVW(J) / TVW
     RTBW(J) = G / (TB(J) + TB(J-1))
     RQVBW(J)= 0.305*G/(1.0 + 0.305*(QVB(J)+QVB(J-1)))
     CPTDZ(J)= DTS*RDZ*CP*0.5*(TB(J)+TB(J-1))*(1.0+0.305*(QVB(J)+QVB(J-1)))
  END DO
  RHOTVW(1  ) = (1000.0/RD) * PNS ** (CV/RD)
  RHOTVW(NP1) = (1000.0/RD) * PNT ** (CV/RD)
  RHOW(  1  ) = 100.0 * RHOTVW(1  ) / TVBS
  RHOW(  NP1) = 100.0 * RHOTVW(NP1) / TBT
  !
  !         SPONGE LAYER DAMPING COEFFICIENT
  !
  DO J = 1, N
     ZSP = (ZS(J) - ZS(NZSP))/(ZS(N) - ZS(NZSP))
     IF(ZSP < 0) TAU(J) = 0.0
     IF(ZSP >=0.AND.ZSP<=0.5) TAU(J)=-0.5*ALPHA*(1-COS(ZSP*PI))
     IF(ZSP > 0.5) TAU(J)=-0.5*ALPHA*(1.0 + PI*(ZSP-0.5))
  END DO
  
  !
  !         ARRAYS TO INCREASE EFFICIENCY
  !
  DTSF = 0.5 * DTS * F
  DTSG = 0.5 * DTS * G
  DTL2 = 0.5 * DTL
  DO I = 1, M
     DTLR(I)  = 0.5  * DTL * RDR /RS(I)
     DTSV(I)  = 0.5  * DTS  / RS(I)
     IF( I /= 1 )  DTSR(I)  = 0.25 * DTS * RDR / R(I)
  END DO
  
  DO J = 1, N
     DTLZ(J)  = 0.5  * DTL * RDZ / RHOT(J)
     DTSZW(J) = 0.25 * DTS * RDZ / RHOW(J)
     DTSRW(J) = 0.25 * DTS * RDR / RHOW(J)
  END DO
  DO J = 1 , N
     DO I = 2 , M
        DTSZ(I,J)  = 0.25 * DTS * RDZ /(RHOT(J) * R(I))
     END DO
  END DO
  
  !
  !        ARRAYS FOR SEMI - IMPLICIT SMALL TIME STEP
  !
  CC4=0.25*(1.0+EP)**2
  DO J = 2, N
     A(J)=    CC4*CPTDZ(J)*RHOTVW(J+1)*    ZC2(J)
     B(J)=1.0+CC4*CPTDZ(J)*RHOTVW(J  )*(ZC2(J)+ZC2(J-1))
     C(J)=    CC4*CPTDZ(J)*RHOTVW(J-1)*   ZC2(J-1)
  END DO
  
  E(1)=0
  D( 1:M,1) = 0
  D9(1:M,1) = 0
  
  DO J = 2, N
     E(J) = A(J)/(B(J)-C(J)*E(J-1))
  END DO

  SST(:)=TAMB+TMID/(1.+(RS(:)/RSST)**2)
  ESS(:)=6.11*EXP(A1*(SST(:)-273.)/(SST(:)-36.))
  
  RETURN
END SUBROUTINE param_setup

SUBROUTINE hurr_initial
  USE hurr_vars
  IMPLICIT NONE
  
  ZD=ZS(nzsp)
  
  E1 = 3
  D2 = 2 * RMAX/(RO + RMAX)
  DO J = 1, N
     DO I = 1, M
        D1 = 2 * RMAX/(RS(I) + RMAX)
        IF( RS(I) < RO ) VR = SQRT(VMAX**2 *(RS(I)/RMAX)**2 * &
             (D1**E1 - D2**E1)+0.25*F**2*RS(I)**2) - 0.5*F*RS(I)
        IF( RS(I) >= RO ) VR = 0
        V9(I,J) = VR * ZS(J)/ZS(NZTR)
        IF(ZS(J) > ZS(NZTR)) V9(I,J) = VR*(ZS(J) - ZD)/(ZS(NZTR) - ZD)
        IF(ZS(J) > ZS(NZSP)) V9(I,J) = 0.0
        V19(I,J)  = V9(I,J)
        QL9(I,J)  = 0
        QL19(I,J) = 0
     END DO
  END DO
  
  DO J = 1, N
     P9(M,J) = 0
     P19(M,J)=P9(M,J)
     DO I = M, 2, -1
        P9(I-1,J) = P9(I,J) - ( DTS / CPTDR(J) ) * 0.5 * &
             (V9(I,J)**2/RS(I) + V9(I-1,J)**2/RS(I-1) + F*(V9(I,J) + V9(I-1,J)))
        P19(I-1,J)= P9(I-1,J)
     END DO
  END DO
  
  DO I = 1, M
     DO J = 2, NM1
        T9(I,J) = TB(J) + 0.25 &
             *(CPTDZ(J+1)*(P9(I,J+1)-P9(I,J  )) &
             + CPTDZ(J  )*(P9(I,J  )-P9(I,J-1)) )/(DTS*RTB(J))
        T19(I,J)=T9(I,J)
     END DO
     T9(I,1) = TB(1) + 0.5*CPTDZ(2)*(P9(I,2)-P9(I,1  ))/(DTS*RTB(1))
     T9(I,N) = TB(N) + 0.5*CPTDZ(N)*(P9(I,N)-P9(I,NM1))/(DTS*RTB(N))
     T19(I,1)=T9(I,1)
     T19(I,N)=T9(I,N)
  END DO
  
  !
  !    GIVE 0'S TO THOSE ARRAY ELEMENTS THAT ARE OUTSIDE PHYSICAL BOUNDARIES (IE
  !    TO ELEMENTS THAT WERE ADDED TO ARRAY DEFINITIONS TO AVOID INDEX OVERFLOW)
  !
  V9( 0:MP1,0  )=0.0
  P9( 0:MP1,0  )=0.0
  T9( 0:MP1,0  )=0.0
  QV9(0:MP1,0  )=0.0
  QL9(0:MP1,0  )=0.0
  U9( 0:MP1,0  )=0.0
  W9( 0:MP1,0  )=0.0
  V9( 0:MP1,NP1)=0.0
  P9( 0:MP1,NP1)=0.0
  T9( 0:MP1,NP1)=0.0
  QV9(0:MP1,NP1)=0.0
  QL9(0:MP1,NP1)=0.0
  U9( 0:MP1,NP1)=0.0
  
  W9(0,NP1)     =0.0
  W9(MP1,NP1)   =0.0
  
  V9( 0  ,0:NP1)=0.0
  P9( 0  ,0:NP1)=0.0
  T9( 0  ,0:NP1)=0.0
  QV9(0  ,0:NP1)=0.0
  QL9(0  ,0:NP1)=0.0
  U9( 0  ,0:NP1)=0.0
  W9( 0  ,0:NP1)=0.0
  V9( MP1,0:NP1)=0.0
  P9( MP1,0:NP1)=0.0
  T9( MP1,0:NP1)=0.0
  QV9(MP1,0:NP1)=0.0
  QL9(MP1,0:NP1)=0.0
  W9( MP1,0:NP1)=0.0
  
  U9(MP1,0)     =0.0
  U9(MP1,NP1)   =0.0
  
  !
  !    ADJUST Qv SO THAT THETAe= constant across vortex
  !
  
  DO J = 1 , N
     QV9( 1:M, J) = QVB(J)
     QV19(1:M, J) = QV9(1:M,J)
  END DO
  
  U9( 1:MP1,1:N) = 0
  U19(1:MP1,1:N) = 0
  
  W9( 1:M,1:NP1) = 0
  W19(1:M,1:NP1) = 0
  
  UA9(1, 1:N) = 0
  
END SUBROUTINE hurr_initial


SUBROUTINE INTERPOLATE(NLEVELS,PDS,DZ,TZ,QZ,SST,VMAX)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NLEVELS
  REAL*8, INTENT(IN) :: PDS, DZ
  REAL*8, INTENT(OUT) :: VMAX, SST
  REAL*8, DIMENSION(NLEVELS), INTENT(OUT) :: TZ, QZ

  INTEGER, PARAMETER ::NA=60
  REAL*8 :: PZ(NLEVELS), TAZ(NLEVELS),QAZ(NLEVELS)
  REAL*8 :: Z(NA),T(NA),Q(NA),P(NA)

  REAL*8 :: AZ, EPSI, G, PBAR, PMIN, RD, TM, TM1, TV1, TV2, TVBAR

  INTEGER :: I, J, IFL, K, NP

  OPEN(UNIT=101,FILE='s.in',STATUS='OLD')
  RD=287.
  EPSI=1./0.622
  G=9.8
  READ(101,20)NP,SST
20 FORMAT(T15,I3,T55,F5.2,////)
  DO I=1,NP
     READ(101,30)P(I),T(I),Q(I)
30   FORMAT(3X,F6.1,T12,F8.3,T25,F10.6)
  END DO
  CLOSE(101)
  DO I=1,NP
     T(I)=T(I)+273.15
     Q(I)=0.001*Q(I)
  END DO
  Z(1)=0.0
  DO I=2,NP
     TV1=T(I-1)*(1.+EPSI*Q(I-1)-Q(I-1))
     TV2=T(I)*(1.+EPSI*Q(I)-Q(I))
     TVBAR=0.5*(TV1+TV2)
     PBAR=0.5*(P(I-1)+P(I))
     Z(I)=Z(I-1)+(RD*TVBAR/(G*PBAR))*(P(I-1)-P(I))
  END DO

  DO I=1,NLEVELS !=20
     AZ=0.5*DZ+DZ*FLOAT(I-1)
     DO J=1,NP !=46
        K=J
        IF(Z(J).GE.AZ) EXIT
     END DO

     TM1=T(K-1)*(1000./P(K-1))**(0.286)
     TM=T(K)*(1000./P(K))**(0.286)
     TZ(I)=TM1+(TM-TM1)*(AZ-Z(K-1))/(Z(K)-Z(K-1))
     QZ(I)=Q(K-1)+(Q(K)-Q(K-1))*(AZ-Z(K-1))/(Z(K)-Z(K-1))

     QAZ(I)=1000.*QZ(I)
     TAZ(I)=T(K-1)+(T(K)-T(K-1))*(AZ-Z(K-1))/(Z(K)-Z(K-1))-273.15
     PZ(I)=P(K-1)+(P(K)-P(K-1))*(AZ-Z(K-1))/(Z(K)-Z(K-1))
  END DO

!
  CALL PCMIN(SST,PDS,PZ,TAZ,QAZ,NLEVELS,NLEVELS-1,PMIN,VMAX,IFL)
!
  RETURN
END SUBROUTINE INTERPOLATE

! ==================================================================

SUBROUTINE PCMIN(SST,PSL,P,T,R,NA,N,PMIN,VMAX,IFL)
!
!   ***   This subroutine calculates the maximum wind speed        ***
!   ***             and mimimum central pressure                   ***
!   ***    achievable in tropical cyclones, given a sounding       ***
!   ***             and a sea surface temperature.                 ***
!
!  INPUT:   SST: Sea surface temperature in C
!
!           PSL: Sea level pressure (mb)
!
!           P,T,R: One-dimensional arrays of dimension NA
!             containing pressure (mb), temperature (C),
!             and mixing ratio (g/kg). The arrays MUST be
!             arranged so that the lowest index corresponds
!             to the lowest model level, with increasing index
!             corresponding to decreasing pressure. The temperature
!             sounding should extend to at least the tropopause and 
!             preferably to the lower stratosphere, however the
!             mixing ratios are not important above the boundary
!             layer. Missing mixing ratios can be replaced by zeros.
!
!           NA: The dimension of P,T and R
!
!           N:  The actual number of points in the sounding
!                (N is less than or equal to NA)
!
!  OUTPUT:  PMIN is the minimum central pressure, in mb
!
!           VMAX is the maximum surface wind speed, in m/s
!                  (reduced to reflect surface drag)
!
!           IFL is a flag: A value of 1 means OK; a value of 0
!              indicates no convergence (hypercane); a value of 2
!              means that the CAPE routine failed.
!
!-----------------------------------------------------------------------------
  REAL*8 :: T(NA), P(NA), R(NA)
!
!   ***   Adjustable constant: Ratio of C_k to C_D    ***
!
  CKCD=1.0
!
!   ***   Adjustable constant for buoyancy of displaced parcels:  ***
!   ***    0=Reversible ascent;  1=Pseudo-adiabatic ascent        ***
!
  SIG=1.0
!
!   ***  Adjustable switch: if IDISS = 0, no dissipative heating is   ***
!   ***     allowed; otherwise, it is                                 ***
!
  IDISS=1
!
!   ***   Normalize certain quantities   ***
!
  SSTK=SST+273.15
  ES0=6.112*EXP(17.67*SST/(243.5+SST))
  DO I=1,N
     R(I)=R(I)*0.001
     T(I)=T(I)+273.15
  END DO

!
!   ***   Default values   ***
!
  VMAX=0.0
  PMIN=PSL
  IFL=1
!
  NP=0
  PM=950.0
!
!   ***   Find environmental CAPE *** 
!
  TP=T(1)
  RP=R(1)
  PP=P(1)
  CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEA,TOA,IFLAG)
  IF(IFLAG.NE.1)IFL=2
!
!   ***   Begin iteration to find mimimum pressure   ***
!
100 CONTINUE
!
!   ***  Find CAPE at radius of maximum winds   ***
!
  TP=T(1)
  PP=PM
  RP=0.622*R(1)*PSL/(PM*(0.622+R(1))-R(1)*PSL)
  CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEM,TOM,IFLAG) 
  IF(IFLAG.NE.1)IFL=2
  RAT=SSTK/TOM
  IF(IDISS.EQ.0)RAT=1.0
!
!  ***  Find saturation CAPE at radius of maximum winds   ***
!
  TP=SSTK
  PP=PM
  RP=0.622*ES0/(PM-ES0)
  CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEMS,TOMS,IFLAG)
  IF(IFLAG.NE.1)IFL=2
!
!  ***  Initial estimate of minimum pressure   ***
!
  RS0=RP
  TV1=T(1)*(1.+R(1)/0.622)/(1.+R(1))
  TVAV=0.5*(TV1+SSTK*(1.+RS0/0.622)/(1.+RS0))
  CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
  PNEW=PSL*EXP(-CAT/(287.04*TVAV))
!
!   ***  Test for convergence   ***
!
  IF(ABS(PNEW-PM).GT.0.2)THEN
     PM=PNEW
     NP=NP+1
     IF(NP.GT.1000.OR.PM.LT.400.0)THEN
        PMIN=400.0
        IFL=0
        GOTO 900
     END IF
     GOTO 100
  ELSE
     CAT=CAPEM-CAPEA+CKCD*RAT*(CAPEMS-CAPEM)
     PMIN=PSL*EXP(-CAT/(287.04*TVAV))
  END IF
900 CONTINUE
  FAC=MAX(0.0,(CAPEMS-CAPEM))
  VMAX=SQRT(CKCD*RAT*FAC)
!
!   ***  Renormalize sounding arrays   ***
  !
  DO I=1,N
     R(I)=R(I)*1000.0
     T(I)=T(I)-273.15
  END DO
!
  RETURN
END SUBROUTINE PCMIN

! ========================================================

SUBROUTINE CAPE(TP,RP,PP,T,R,P,ND,N,SIG,CAPED,TO,IFLAG)
!
!     This subroutine calculates the CAPE of a parcel with pressure PP (mb), 
!       temperature TP (K) and mixing ratio RP (gm/gm) given a sounding
!       of temperature (T in K) and mixing ratio (R in gm/gm) as a function
!       of pressure (P in mb). ND is the dimension of the arrays T,R and P,
!       while N is the actual number of points in the sounding. CAPED is
!       the calculated value of CAPE and TO is the temperature at the
!       level of neutral buoyancy.  IFLAG is a flag
!       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
!       not run owing to improper sounding (e.g.no water vapor at parcel level).
!       IFLAG=2 indicates that routine did not converge.                 
!
  REAL*8 :: T(ND),R(ND),P(ND),TVRDIF(100)   
  REAL*8 :: NA
!
!   ***   Default values   ***
!      
  CAPED=0.0
  TO=T(1)
  IFLAG=1
!
!   ***   Check that sounding is suitable    ***
!
  IF(RP.LT.1.0E-6.OR.TP.LT.200.0)THEN
     IFLAG=0
     RETURN
  END IF
!
!   ***   Assign values of thermodynamic constants     ***
!
  CPD=1005.7
  CPV=1870.0
  CL=4190.0
  CPVMCL=2320.0
  RV=461.5
  RD=287.04
  EPS=RD/RV
  ALV0=2.501E6
!
!   ***  Define various parcel quantities, including reversible   ***
!   ***                       entropy, S.                         ***
!                           
  TPC=TP-273.15
  ESP=6.112*EXP(17.67*TPC/(243.5+TPC))
  EVP=RP*PP/(EPS+RP)
  RH=EVP/ESP
  ALV=ALV0-CPVMCL*TPC
  S=(CPD+RP*CL)*LOG(TP)-RD*LOG(PP-EVP)+ALV*RP/TP-RP*RV*LOG(RH)            
!
!   ***  Find lifted condensation pressure, PLCL   ***
!     
  CHI=TP/(1669.0-122.0*RH-TP)
  PLCL=PP*(RH**CHI)
!
!   ***  Begin updraft loop   ***
!
  NCMAX=0
  DO J=2,N
!
!    ***   Don't bother lifting parcel above 60 mb    ***
!
     IF(P(J).LT.59.0) CYCLE
!
!    ***  Parcel quantities below lifted condensation level   ***
! 
     IF(P(J).GE.PLCL) THEN
        TG=TP*(P(J)/PP)**(RD/CPD)
        RG=RP
!
!   ***   Calculate buoyancy   ***
!  
        TLVR=TG*(1.+RG/EPS)/(1.+RG)
        TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
     ELSE
!
!   ***  Parcel quantities above lifted condensation level  ***
! 
        TG=T(J)          
        TJC=T(J)-273.15 
        ES=6.112*EXP(17.67*TJC/(243.5+TJC)) 
        RG=EPS*ES/(P(J)-ES)
!
!   ***  Iteratively calculate lifted parcel temperature and mixing   ***
!   ***                ratio for reversible ascent                    ***
!
        NC=0
120     CONTINUE
        NC=NC+1
!
!   ***  Calculate estimates of the rates of change of the entropy    ***
!   ***           with temperature at constant pressure               ***
!  
        ALV=ALV0-CPVMCL*(TG-273.15)
        SL=(CPD+RP*CL+ALV*ALV*RG/(RV*TG*TG))/TG
        EM=RG*P(J)/(EPS+RG)
        SG=(CPD+RP*CL)*LOG(TG)-RD*LOG(P(J)-EM)+ALV*RG/TG
        IF(NC.LT.3)THEN
           AP=0.3
        ELSE
           AP=1.0
        END IF
        TGNEW=TG+AP*(S-SG)/SL  
!
!   ***   Test for convergence   ***
!
        IF(ABS(TGNEW-TG).GT.0.01) THEN
           TG=TGNEW
           TC=TG-273.15
           ENEW=6.112*EXP(17.67*TC/(243.5+TC))
!
!   ***   Bail out if things get out of hand   ***
!
           IF(NC.GT.500.OR.ENEW.GT.(P(J)-1.0))THEN
              IFLAG=2
              RETURN
           END IF
           RG=EPS*ENEW/(P(J)-ENEW)           
           GOTO 120
        END IF
        NCMAX=MAX(NC,NCMAX)
!
!   *** Calculate buoyancy   ***
!
        RMEAN=SIG*RG+(1.-SIG)*RP
        TLVR=TG*(1.+RG/EPS)/(1.+RMEAN)
        TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
     END IF
  END DO  !  DO J=2,N
!
!  ***  Begin loop to find NA, PA, and CAPE from reversible ascent ***
!
  NA=0.0
  PA=0.0
!
!   ***  Find maximum level of positive buoyancy, INB    ***
!
  INB=1
  DO J=N,1,-1
     IF(TVRDIF(J).GT.0.0)INB=MAX(INB,J)
  END DO
  IF(INB.EQ.1)RETURN
!
!   ***  Find positive and negative areas and CAPE  ***
!
  IF(INB.GT.1)THEN
     DO J=2,INB
        TVM=0.5*(TVRDIF(J)+TVRDIF(J-1))
        PMA=0.5*(P(J)+P(J-1))
        IF(TVM.LE.0.0)THEN
           NA=NA-RD*TVM*(P(J-1)-P(J))/PMA
        ELSE
           PA=PA+RD*TVM*(P(J-1)-P(J))/PMA
        END IF
     END DO

!
!   ***   Find residual positive area above INB and TO  ***
!
     PAT=0.0
     TO=T(INB)
     IF(INB.LT.N)THEN
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/ &
             (TVRDIF(INB)-TVRDIF(INB+1))
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB)
        TO=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/ &
             (P(INB)-P(INB+1))
     END IF
!
!   ***   Find CAPE  ***
!            
     CAPED=PA+PAT-NA
     CAPED=MAX(CAPED,0.0)
  END IF
!
  RETURN
END SUBROUTINE CAPE
