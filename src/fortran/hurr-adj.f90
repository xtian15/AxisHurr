SUBROUTINE hurr_adj
  USE hurr_vars
  IMPLICIT NONE

  inU9 = U9
  inV9 = V9
  inW9 = W9
  inP9 = P9
  inT9 = T9
  inQV9 = QV9
  inQL9 = QL9
  inU19 = U19
  inV19 = V19
  inW19 = W19
  inT19 = T19
  inP19 = P19
  inQV19 = QV19
  inQL19 = QL19

  DO J=N,1,-1
     DO I=M,1,-1
        QLTEMP   = QL1(I,J)
        QL1(I,J) = 0
        QVTEMP   = QV1(I,J)
        QV1(I,J) = 0
        PTEMP    = P1(I,J)
        P1( I,J) = 0
        TTEMP    = T1(I,J)
        T1( I,J) = 0
        WTEMP    = W1(I,J+1)
        W1(I,J+1)= 0
        VTEMP    = V1(I,J)
        V1( I,J) = 0
        UTEMP    = U1(I+1,J)
        U1(I+1,J)= 0
        
        QL1(I,J) = QL(I,J)
        QL( I,J) = 0
        QV1(I,J) = QV(I,J)
        QV( I,J) = 0
        T1(I,J)  = T(I,J)
        T( I,J)  = 0
        P1(I,J)  = P(I,J)
        P( I,J)  = 0
        W1(I,J+1)= W(I,J+1)
        W( I,J+1)= 0
        V1(I,J)  = V(I,J)
        V( I,J)  = 0
        U1(I+1,J)= U(I+1,J)
        U( I+1,J)= 0

        
        QL1(I,J) = EPS * QLTEMP + QL1(I,J)
        QL(I,J)  = QLTEMP
        QLTEMP   = 0

        QV1(I,J) = EPS * QVTEMP + QV1(I,J)
        QV(I,J)  = QVTEMP
        QVTEMP   = 0

        T1(I,J)  = EPS * TTEMP + T1(I,J)
        T(I,J)   = TTEMP
        TTEMP    = 0

        P1(I,J)  = EPS * PTEMP + P1(I,J)
        P(I,J)   = PTEMP
        PTEMP    = 0

        W1(I,J+1)= EPS * WTEMP + W1(I,J+1)
        W(I,J+1) = WTEMP
        WTEMP    = 0

        V1(I,J)  = EPS * VTEMP + V1(I,J)
        V(I,J)   = VTEMP
        VTEMP    = 0

        U1(I+1,J)= EPS * UTEMP + U1(I+1,J)
        U(I+1,J) = UTEMP
        UTEMP    = 0
     END DO
  END DO

  ! variables assigned with values include:
  ! QL1, QV1, T1, P1, V1 (1:M,1:N)
  ! QL , QV , T , P , V  (1:M,1:N)
  ! W1 , W (1:M,2:N+1)
  ! U1 , U (2:N+1,1:M)

  CALL FWD_NO_VADV

  ! backing up variable values that will be repetitively updated
  
  bkV19  = V19
  bkT19  = T19
  bkQV19 = QV19
  bkQL19 = QL19
  ! ADVANCE V, T, QV, QL
  DO J=N, 1, -1
     DO I=M, 1, -1
        ! ===== fwd =====
        V19(I,J)  = V19(I,J)  + VA9(I,J)
        T19(I,J)  = T19(I,J)  + TA9(I,J)
        QV19(I,J) = QV19(I,J) + QVA9(I,J)
        QL19(I,J) = QL19(I,J) + QLA9(I,J)
        IF(QV19(I,J) < 0) QV19(I,J)= 0
        IF(QL19(I,J) < 0) QL19(I,J)=0
        !
        !     CONDENSATION / EVAPORATION
        !
        PNST9= PN(J) + P19(I,J) - P19(M,1)
        PDST9= 1000.0 * PNST9**(1/XKAPPA)
        TEMP9= PNST9 * T19(I,J)        
        ES9  = 6.11 * EXP(A1*(TEMP9-273)/(TEMP9-36))
        QSS9 = 0.622
        IF(PDST9-ES9 > ES9) QSS9 = 0.622*ES9/(PDST9-ES9)
        
        !     NO EVAPORATION NOR CONDENSATION
        IF(QV19(I,J)>=QSS9 .OR. QL19(I,J)>0) THEN
           R19 =  1/(1 + XLDCP * 237 * A1 * QSS9/(TEMP9-36)**2)
           QVD9 = R19*(QV19(I,J)-QSS9)
           !     NO EVAPORATION
           IF ((QL19(I,J)+QVD9) < 0) THEN              
              !     EVAPORATE
              T19( I,J) = T19(I,J)  - XLDCP * QL19(I,J)/PNST9
              QV19(I,J) = QV19(I,J) + QL19(I,J)
              QL19(I,J) = 0
           ELSE
              !     CONDENSE
              T19(I,J)  = T19(I,J)  + XLDCP * QVD9/PNST9
              QV19(I,J) = QV19(I,J)         - QVD9
              QL19(I,J) = QL19(I,J)         + QVD9
           END IF
        END IF
        ! ===============
        
        IF(QL19(I,J) < 0) QL1(I,J) = 0
        IF(QV19(I,J) < 0) QV1(I,J) = 0
        
        ! ===== fwd =====
        V19(I,J)  = bkV19(I,J)  + VA9(I,J)
        T19(I,J)  = bkT19(I,J)  + TA9(I,J)
        QV19(I,J) = bkQV19(I,J) + QVA9(I,J)
        QL19(I,J) = bkQL19(I,J) + QLA9(I,J)
        IF(QV19(I,J) < 0) QV19(I,J)= 0
        IF(QL19(I,J) < 0) QL19(I,J)=0
        !
        !     CONDENSATION / EVAPORATION
        !
        PNST9= PN(J) + P19(I,J) - P19(M,1)
        PDST9= 1000.0 * PNST9**(1/XKAPPA)
        TEMP9= PNST9 * T19(I,J)
        ES9  = 6.11 * EXP(A1*(TEMP9-273)/(TEMP9-36))
        QSS9 = 0.622
        IF(PDST9-ES9 > ES9) QSS9 = 0.622*ES9/(PDST9-ES9)
        ! ===============

        ! NO EVAPORATION NOR CONDENSATION
        IF(QV19(I,J)>=QSS9 .OR. QL19(I,J)>0) THEN
           ! ===== fwd =====
           R19 =  1/(1 + XLDCP * 237 * A1 * QSS9/(TEMP9-36)**2)
           QVD9 = R19*(QV19(I,J)-QSS9)
           ! ===============

           IF ((QL19(I,J)+QVD9) < 0) THEN
              QL1(I,J) = 0
              QL1(I,J) = QV1(I,J) + QL1(I,J)
              PNST = XLDCP*QL19(I,J)/PNST9**2 * T1(I,J)
              QL1(I,J) = -XLDCP/PNST9 * T1(I,J) + QL1(I,J)
              QVD = 0
           ELSE
              QVD = QL1(I,J)
              QVD = -QV1(I,J) + QVD
              PNST = -XLDCP * QVD9/PNST9**2 * T1(I,J)
              QVD = XLDCP/PNST9 * T1(I,J) + QVD
           END IF

           QSS = -R19 * QVD
           QV1(I,J) = R19 * QVD + QV1(I,J)
           R1  = (QV19(I,J)-QSS9) * QVD
           QVD = 0

           TEMP = -1/(1 + XLDCP * 237 * A1 * QSS9/(TEMP9-36)**2)**2 &
                *(-2 *    XLDCP * 237 * A1 * QSS9/(TEMP9-36)**3)*R1
           QSS  = -1/(1 + XLDCP * 237 * A1 * QSS9/(TEMP9-36)**2)**2 &
                *(        XLDCP * 237 * A1       /(TEMP9-36)**2)*R1 + QSS
           R1   = 0
        ELSE
           PNST = 0
           QVD  = 0
           QSS  = 0
           R1   = 0
           TEMP = 0
        END IF

        ! --- CONDENSATION / EVAPORATION ---
        IF(PDST9-ES9 > ES9) THEN
           ES   =  0.622*ES9/(PDST9-ES9)**2 * QSS
           PDST = -0.622*ES9/(PDST9-ES9)**2 * QSS
           ES   =  0.622    /(PDST9-ES9)    * QSS + ES
           QSS  =  0
        ELSE
           ES = 0
           PDST = 0
        END IF
        QSS = 0

        TEMP = 6.11 * EXP(A1*(TEMP9-273)/(TEMP9-36)) * A1*237/(TEMP9-36)**2 &
             *ES + TEMP
        ES   = 0

        T1(I,J) = PNST9 * TEMP + T1(I,J)
        PNST    = T19(I,J) * TEMP + PNST
        TEMP    = 0

        PNST    = 1000.0 * (1/XKAPPA)*PNST9**(1/XKAPPA-1)*PDST + PNST
        PDST    = 0

        P1(M,1) = -PNST + P1(M,1)
        P1(I,J) =  PNST + P1(I,J)
        PNST    = 0

        ! ----------------------------------

        ! ===== fwd =====
        V19(I,J)  = bkV19(I,J)  + VA9(I,J)
        T19(I,J)  = bkT19(I,J)  + TA9(I,J)
        QV19(I,J) = bkQV19(I,J) + QVA9(I,J)
        QL19(I,J) = bkQL19(I,J) + QLA9(I,J)
        ! ===============
        IF(QL19(I,J) < 0) QL1(I,J) = 0
        IF(QV19(I,J) < 0) QV1(I,J) = 0

        QLA(I,J) = QL1(I,J)
        QVA(I,J) = QV1(I,J)
        TA( I,J) = T1( I,J)
        VA( I,J) = V1( I,J)
     END DO
  END DO

  ! variables assigned with values include:
  ! PNST, QVD, QSS, R1, TEMP, ES, PDST
  ! QLA, QVA, TA, VA (1:M,1:N)

  V19  = bkV19
  T19  = bkT19
  QV19 = bkQV19
  QL19 = bkQL19
  !
  !     OUTER BOUNDARY --> SMALL TIME STEP
  !
  UA=0
  DO J=N, 1, -1
     UA(MP1,J) = U1(MP1,J)+UA(MP1,J)
  END DO

  PS = 0
  D  = 0
  WS = 0
  WA = 0

  DO NSMALL = NS, 1, -1
     DO J=1, N
        DO I=M, 1, -1
           ! P1(I,J) = 
           W1(I,J  ) =  0.5*(1.0+EP)*ZC2(J)*RHOTVW(J  )*P1(I,J) + W1(I,J  )
           W1(I,J+1) = -0.5*(1.0+EP)*ZC2(J)*RHOTVW(J+1)*P1(I,J) + W1(I,J+1)
           PS(I,J  ) =  P1(I,J) + PS(I,J)
           P1(I,J  ) =  0
           ! W1(I,J) = 
           D(I,J)    = W1(I,J) + D(I,J)
           W1(I,J+1) = E(J)*W1(I,J) + W1(I,J+1)
           W1(I,J  ) = 0
        END DO
     END DO

     DO J=N, 2, -1
        DO I=M, 1, -1
           CFACD9=DTS*RDZ*CP*0.5* &
                (T19(I,J)+T19(I,J-1))*(1.0+0.305*(QV9(I,J)+QV9(I,J-1)))
           
           ! D(I,J) = 
           D( I,J-1) =  C(J)*E(J)/A(J)*D(I,J) + D(I,J-1)
           PS(I,J-1) =  CFACD9*0.5*(1+EP)*E(J)/A(J)*D(I,J) + PS(I,J-1)
           PS(I,J  ) = -CFACD9*0.5*(1+EP)*E(J)/A(J)*D(I,J) + PS(I,J  )
           CFACD     = -0.5*(1+EP)*(PS9NS(I,J,NSMALL)-PS9NS(I,J-1,NSMALL)) &
                *E(J)/A(J)*D(I,J)
           WS( I,J)  = E(J)/A(J)*D(I,J) + WS(I,J)
           D(I,J)    = 0

           ! CFACD = 
           QV(I,J-1) = DTS*RDZ*CP*0.5*(T19(I,J)+T19(I,J-1))*0.305*CFACD + QV(I,J-1)
           QV(I,J  ) = DTS*RDZ*CP*0.5*(T19(I,J)+T19(I,J-1))*0.305*CFACD + QV(I,J  )
           T1(I,J-1) = DTS*RDZ*CP*0.5*(1.0+0.305*(QV9(I,J)+QV9(I,J-1))) &
                *CFACD + T1(I,J-1)
           T1(I,J  ) = DTS*RDZ*CP*0.5*(1.0+0.305*(QV9(I,J)+QV9(I,J-1))) &
                *CFACD + T1(I,J  )
           CFACD=0
        END DO
     END DO

     DO J=N,1,-1
        DO I=M,1,-1
           ! PS(I,J) =
           W1(I,J  ) =  0.5*(1.0-EP)*ZC2(J)*RHOTVW(J  )*PS(I,J) + W1(I,J  )
           W1(I,J+1) = -0.5*(1.0-EP)*ZC2(J)*RHOTVW(J+1)*PS(I,J) + W1(I,J+1)
           U1(I  ,J) =  RC2(J)/RS(I) * R(I  )*PS(I,J) + U1(I  ,J)
           U1(I+1,J) = -RC2(J)/RS(I) * R(I+1)*PS(I,J) + U1(I+1,J)
           P1(I,J)   =  PS(I,J) + P1(I,J)
           PS(I,J)   = 0
        END DO
     END DO

     DO J=N,2,-1
        DO I=M,1,-1
           CFACWS9=DTS*RDZ*CP*0.5* &
                (T19(I,J)+T19(I,J-1)) * (1.0+0.305*(QV9(I,J)+QV9(I,J-1)))
           
           ! WS(I,J) = 
           WA(I,J)   =  WS(I,J) + WA(I,J)
           P1(I,J-1) =  0.5*(1-EP)*CFACWS9*WS(I,J) + P1(I,J-1)
           P1(I,J  ) = -0.5*(1-EP)*CFACWS9*WS(I,J) + P1(I,J  )
           CFACWS    = -0.5*(1-EP)*(P19NS(I,J,NSMALL)-P19NS(I,J-1,NSMALL)) &
                *WS(I,J)
           W1(I,J)   =  WS(I,J) + W1(I,J)
           WS(I,J)   =  0

           ! CFACWS = 
           QV(I,J-1) =  DTS*RDZ*CP*0.5*(T19(I,J)+T19(I,J-1))*0.305*CFACWS &
                + QV(I,J-1)
           QV(I,J  ) =  DTS*RDZ*CP*0.5*(T19(I,J)+T19(I,J-1))*0.305*CFACWS &
                + QV(I,J  )
           T1(I,J-1) =  DTS*RDZ*CP*0.5*(1.0+0.305*(QV9(I,J)+QV9(I,J-1))) &
                *CFACWS + T1(I,J-1)
           T1(I,J  ) =  DTS*RDZ*CP*0.5*(1.0+0.305*(QV9(I,J)+QV9(I,J-1))) &
                *CFACWS + T1(I,J  )
           CFACWS    = 0
        END DO
     END DO

     DO J=N,1,-1
        DO I=M,2,-1
           CFACU19=DTS*RDR*CP*0.5&
                *(T19(I,J)*(1.0+0.61*QV9(I,J))+T19(I-1,J)*(1.0+0.61*QV9(I-1,J)))
           
           ! U1(I,J) = 
           UA(I,J) = U1(I,J) + UA(I,J)
           P1(I-1,J) =  CFACU19*U1(I,J) + P1(I-1,J)
           P1(I  ,J) = -CFACU19*U1(I,J) + P1(I  ,J)
           CFACU1    = -(P19NS(I,J,NSMALL)-P19NS(I-1,J,NSMALL))*U1(I,J)

           ! CFACU1 = 
           QV(I-1,J) = DTS*RDR*CP*0.5*T19(I-1,J)*0.61*CFACU1 + QV(I-1,J)
           QV(I  ,J) = DTS*RDR*CP*0.5*T19(I  ,J)*0.61*CFACU1 + QV(I  ,J)
           T1(I-1,J) = DTS*RDR*CP*0.5*(1.0+0.61*QV9(I-1,J))*CFACU1 + T1(I-1,J)
           T1(I  ,J) = DTS*RDR*CP*0.5*(1.0+0.61*QV9(I  ,J))*CFACU1 + T1(I  ,J)
        END DO
     END DO
  END DO  ! NSMALL = NS, 1, -1

  ! TIME SMOOTHER
  DO J=N,1,-1
     DO I=M,1,-1
        QL1(I,J) = EPS*QL(I,J) + QL1(I,J)
        QL( I,J) = (1 - 2*EPS)*QL(I,J)

        QV1(I,J) = EPS*QV(I,J) + QV1(I,J)
        QV( I,J) = (1 - 2*EPS)*QV(I,J)

        T1( I,J) = EPS*T( I,J) + T1( I,J)
        T(  I,J) = (1 - 2*EPS)*T( I,J)

        P1( I,J) = EPS*P( I,J) + P1( I,J)
        P(  I,J) = (1 - 2*EPS)*P( I,J)

        W1(I,J+1)= EPS*W(I,J+1) + W1(I,J+1)
        W( I,J+1)= (1 - 2*EPS)*W(I,J+1)

        V1( I,J) = EPS*V( I,J) + V1( I,J)
        V(  I,J) = (1 - 2*EPS)*V( I,J)

        U1(I+1,J)= EPS*U(I+1,J) + U1(I+1,J)
        U( I+1,J)= (1 - 2*EPS)*U(I+1,J)
     END DO
  END DO

  U9 = inU9
  V9 = inV9
  W9 = inW9
  P9 = inP9
  T9 = inT9
  QV9 = inQV9
  QL9 = inQL9
  U19 = inU19
  V19 = inV19
  W19 = inW19
  T19 = inT19
  QV19 = inQV19
  QL19 = inQL19
  P19 = inP19
  CALL FWD_NO_SMOOTHER
  
  !
  !      RAIN FALL
  !

  ! BOUNDARY
  DO I = M,1,-1
     DUM(I,2)=0.
     IF((QL19(I,2) - 0.001) > 0) DUM(I,2)=VTERM
     DUM(I,1)=0
     IF((QL19(I,1) - 0.001) > 0) DUM(I,1)=VTERM
     QL1(I,1) = -2*DTLZ(1)*DUM(I,1)*RHOT(1)*QLA(I,1) + QL1(I,1)
     QL1(I,2) =  2*DTLZ(1)*DUM(I,2)*RHOT(2)*QLA(I,1) + QL1(I,2)
  END DO

  DO J = 1, N
     DO I = 1, M
        DUM(I,J)=0
        IF((QL9(I,J) - 0.001) > 0) DUM(I,J)=VTERM
     END DO
  END DO

  DO J=NM1, 2, -1
     DO I=M, 1, -1
        QL(I,J-1) = -DTLZ(J)*DUM(I,J-1)*RHOT(J-1)*QLA(I,J) + QL(I,J-1)
        QL(I,J+1) =  DTLZ(J)*DUM(I,J+1)*RHOT(J+1)*QLA(I,J) + QL(I,J+1)
     END DO
  END DO
  !       FORCING FOR  V, T, QV, QL  EQUATIONS

  ! BOUNDARY

  I=M
  UB=0
  QRAD=0
  DO J=N, 1, -1
     ! QLA(I,J) = 
     IF(UB9(J)>0) THEN
        U(M  ,J)=-DTLR(I)*(QL19(I  ,J  )-QL19(I-1,J  ))*R(M  )*QLA(I,J) + U(M  ,J)
        U(MP1,J)=-DTLR(I)*(QL19(I  ,J  )-QL19(I-1,J  ))*R(MP1)*QLA(I,J) + U(MP1,J)
     END IF
     QL(I,J-1) =  DTLZ(J)*RHOW(J  )*W9(I,J  )*QLA(I,J) + QL(I,J-1)
     QL(I,J  ) = -DTLZ(J)*RHOW(J  )*W9(I,J  )*QLA(I,J) + QL(I,J  )
     W( I,J  ) = -DTLZ(J)*RHOW(J  )*(QL9(I  ,J  )-QL9(I  ,J-1))*QLA(I,J) + W(I,J  )
     QL(I,J  ) =  DTLZ(J)*RHOW(J+1)*W9(I,J+1)*QLA(I,J) + QL(I,J  )
     QL(I,J+1) = -DTLZ(J)*RHOW(J+1)*W9(I,J+1)*QLA(I,J) + QL(I,J+1)
     W( I,J+1) = -DTLZ(J)*RHOW(J+1)*(QL9(I  ,J+1)-QL9(I  ,J  ))*QLA(I,J) + W(I,J+1)
     QL1(I-1,J)=  DTLR(I)*UB9(J)*QLA(I,J) + QL1(I-1,J)
     QL1(I  ,J)= -DTLR(I)*UB9(J)*QLA(I,J) + QL1(I  ,J)
     !UB(J)     = -DTLR(I)*(QL19(I  ,J  )-QL19(I-1,J  ))*QLA(I,J) + UB(J)

     ! QVA(I,J) = 
     IF(UB9(J)>0) THEN
        U(M  ,J)=-DTLR(I)*(QV19(I  ,J  )-QV19(I-1,J  ))*R(M  )*QVA(I,J) + U(M  ,J)
        U(MP1,J)=-DTLR(I)*(QV19(I  ,J  )-QV19(I-1,J  ))*R(MP1)*QVA(I,J) + U(MP1,J)
     END IF
     QV(I,J-1) =  DTLZ(J)*RHOW(J  )*W9(I,J  )*QVA(I,J) + QV(I,J-1)
     QV(I,J  ) = -DTLZ(J)*RHOW(J  )*W9(I,J  )*QVA(I,J) + QV(I,J  )
     W( I,J  ) = -DTLZ(J)*RHOW(J  )*(QV9(I  ,J  )-QV9(I  ,J-1))*QVA(I,J) + W(I,J)
     QV(I,J  ) =  DTLZ(J)*RHOW(J+1)*W9(I,J+1)*QVA(I,J) + QV(I,J  )
     QV(I,J+1) = -DTLZ(J)*RHOW(J+1)*W9(I,J+1)*QVA(I,J) + QV(I,J+1)
     W( I,J+1) = -DTLZ(J)*RHOW(J+1)*(QV9(I  ,J+1)-QV9(I  ,J  ))*QVA(I,J) + W(I,J+1)
     QV1(I-1,J)=  DTLR(I)*UB9(J)*QVA(I,J) + QV1(I-1,J)
     QV1(I  ,J)= -DTLR(I)*UB9(J)*QVA(I,J) + QV1(I  ,J)
     !UB(J)     = -DTLR(I)*(QV19(I  ,J  )-QV19(I-1,J  ))*QVA(I,J) + UB(J)

     ! TA(I,J) = 
     IF(UB9(J)>0) THEN
       U(M  ,J)=-DTLR(I)*(T19(I  ,J  )-T19(I-1,J  ))*R(M  )*TA(I,J)+U(M  ,J)
       U(MP1,J)=-DTLR(I)*(T19(I  ,J  )-T19(I-1,J  ))*R(MP1)*TA(I,J)+U(MP1,J)
     END IF
     !QRAD(I,J) =  TA(I,J) + QRAD(I,J)
     T(I,J-1)  =  DTLZ(J)*RHOW(J  )*W9(I,J  )*TA(I,J) + T(I,J-1)
     T(I,J  )  = -DTLZ(J)*RHOW(J  )*W9(I,J  )*TA(I,J) + T(I,J  )
     W(I,J  )  = -DTLZ(J)*RHOW(J  )*(T9(I  ,J  )-T9(I  ,J-1))*TA(I,J) + W(I,J)
     T(I,J  )  =  DTLZ(J)*RHOW(J+1)*W9(I,J+1)*TA(I,J) + T(I,J  )
     T(I,J+1)  = -DTLZ(J)*RHOW(J+1)*W9(I,J+1)*TA(I,J) + T(I,J+1)
     W(I,J+1)  = -DTLZ(J)*RHOW(J+1)*(T9(I  ,J+1)-T9(I  ,J  ))*TA(I,J) + W(I,J+1)
     T1(I-1,J) =  DTLR(I)*UB9(J)*TA(I,J) + T1(I-1,J)
     T1(I  ,J) = -DTLR(I)*UB9(J)*TA(I,J) + T1(I  ,J)
     !UB(J)     = -DTLR(I)*(T19(I  ,J  )-T19(I-1,J  ))*TA(I,J) + UB(J)

     ! VA(I,J) = 
     IF(UB9(J)>0) THEN
        U(M  ,J)=-DTLR(I)*(V19(I  ,J  )-V19(I-1,J  ))*R(M  )*VA(I,J) + U(M  ,J)
        U(MP1,J)=-DTLR(I)*(V19(I  ,J  )-V19(I-1,J  ))*R(MP1)*VA(I,J) + U(MP1,J)
     END IF
     U(I  ,J) = -DTL2*(F+ V9(I,J)/RS(I))/RS(I) * R(I  )*VA(I,J) + U(I  ,J)
     U(I+1,J) = -DTL2*(F+ V9(I,J)/RS(I))/RS(I) * R(I+1)*VA(I,J) + U(I+1,J)
     V(I  ,J) = -DTL2*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)*(1/RS(I))*VA(I,J) &
          + V(I,J)
     V(I,J-1) =  DTLZ(J)*RHOW(J  )*W9(I,J  )*VA(I,J) + V(I,J-1)
     V(I,J  ) = -DTLZ(J)*RHOW(J  )*W9(I,J  )*VA(I,J) + V(I,J  )
     W(I,J  ) = -DTLZ(J)*RHOW(J  )*(V9(I  ,J  )-V9(I  ,J-1))*VA(I,J) + W(I,J)
     V(I,J  ) =  DTLZ(J)*RHOW(J+1)*W9(I,J+1)*VA(I,J) + V(I,J  )
     V(I,J+1) = -DTLZ(J)*RHOW(J+1)*W9(I,J+1)*VA(I,J) + V(I,J+1)
     W(I,J+1) = -DTLZ(J)*RHOW(J+1)*(V9(I  ,J+1)-V9(I  ,J  ))*VA(I,J) + W(I,J+1)
     V1(I-1,J)=  DTLR(I)*UB9(J)*VA(I,J) + V1(I-1,J)
     V1(I  ,J)= -DTLR(I)*UB9(J)*VA(I,J) + V1(I  ,J)
     !UB(J)    = -DTLR(I)*(V19(I  ,J  )-V19(I-1,J  ))*VA(I,J) + UB(J)     
  END DO

  DO J=N, 1, -1
     IF(UB9(J) < 0) UB(J) = 0
     U(M  ,J) = R(M  )*UB(J) + U(M  ,J)
     U(MP1,J) = R(MP1)*UB(J) + U(MP1,J)
     UB(J)    = 0
  END DO

  ! FORCING FOR  V, T, QV, QL  EQUATIONS
  DO J=N, 1, -1
     DO I=MM1, 1, -1
        ! QLA(I,J) = 
        QL(I,J-1) =  DTLZ(J)*RHOW(J  )*W9(I,J  )*QLA(I,J) + QL(I,J-1)
        QL(I,J  ) = -DTLZ(J)*RHOW(J  )*W9(I,J  )*QLA(I,J) + QL(I,J  )
        W( I,J  ) = -DTLZ(J)*RHOW(J  )*(QL9(I  ,J  )-QL9(I  ,J-1))*QLA(I,J) &
             + W(I,J  )
        QL(I,J  ) =  DTLZ(J)*RHOW(J+1)*W9(I,J+1)*QLA(I,J) + QL(I,J  )
        QL(I,J+1) = -DTLZ(J)*RHOW(J+1)*W9(I,J+1)*QLA(I,J) + QL(I,J+1)
        W( I,J+1) = -DTLZ(J)*RHOW(J+1)*(QL9(I  ,J+1)-QL9(I  ,J  ))*QLA(I,J) &
             + W(I,J+1)
        QL(I-1,J) =  DTLR(I)*R(I  )*U9(I  ,J)*QLA(I,J) + QL(I-1,J)
        QL(I  ,J) = -DTLR(I)*R(I  )*U9(I  ,J)*QLA(I,J) + QL(I  ,J)
        U( I  ,J) = -DTLR(I)*R(I  )*(QL9(I  ,J  )-QL9(I-1,J  ))*QLA(I,J) + U(I  ,J)
        QL(I  ,J) =  DTLR(I)*R(I+1)*U9(I+1,J)*QLA(I,J) + QL(I  ,J)
        QL(I+1,J) = -DTLR(I)*R(I+1)*U9(I+1,J)*QLA(I,J) + QL(I+1,J)
        U( I+1,J) = -DTLR(I)*R(I+1)*(QL9(I+1,J  )-QL9(I  ,J  ))*QLA(I,J) + U(I+1,J)
        
        
        ! QVA(I,J) = 
        QV(I,J-1) =  DTLZ(J)*RHOW(J  )*W9(I,J  )*QVA(I,J) + QV(I,J-1)
        QV(I,J  ) = -DTLZ(J)*RHOW(J  )*W9(I,J  )*QVA(I,J) + QV(I,J  )
        W( I,J  ) = -DTLZ(J)*RHOW(J  )*(QV9(I  ,J  )-QV9(I  ,J-1))*QVA(I,J) &
             + W(I,J  )
        QV(I,J  ) =  DTLZ(J)*RHOW(J+1)*W9(I,J+1)*QVA(I,J) + QV(I,J  )
        QV(I,J+1) = -DTLZ(J)*RHOW(J+1)*W9(I,J+1)*QVA(I,J) + QV(I,J+1)
        W( I,J+1) = -DTLZ(J)*RHOW(J+1)*(QV9(I  ,J+1)-QV9(I  ,J  ))*QVA(I,J) &
             + W(I,J+1)
        QV(I-1,J) =  DTLR(I)*R(I  )*U9(I  ,J)*QVA(I,J) + QV(I-1,J)
        QV(I  ,J) = -DTLR(I)*R(I  )*U9(I  ,J)*QVA(I,J) + QV(I  ,J)
        U( I  ,J) = -DTLR(I)*R(I  )*(QV9(I  ,J  )-QV9(I-1,J  ))*QVA(I,J) &
             + U(I  ,J)
        QV(I  ,J) =  DTLR(I)*R(I+1)*U9(I+1,J)*QVA(I,J) + QV(I  ,J)
        QV(I+1,J) = -DTLR(I)*R(I+1)*U9(I+1,J)*QVA(I,J) + QV(I+1,J)
        U( I+1,J) = -DTLR(I)*R(I+1)*(QV9(I+1,J  )-QV9(I  ,J  ))*QVA(I,J) &
             + U(I+1,J)

        ! TA(I,J) =
        !QRAD(I,J) =  TA(I,J) + QRAD(I,J)
        T(I,J-1)  =  DTLZ(J)*RHOW(J  )*W9(I,J  )*TA(I,J) + T(I,J-1)
        T(I,J  )  = -DTLZ(J)*RHOW(J  )*W9(I,J  )*TA(I,J) + T(I,J  )
        W(I,J  )  = -DTLZ(J)*RHOW(J  )*(T9(I  ,J  )-T9(I  ,J-1))*TA(I,J) &
             + W(I,J  )
        T(I,J  )  =  DTLZ(J)*RHOW(J+1)*W9(I,J+1)*TA(I,J) + T(I,J  )
        T(I,J+1)  = -DTLZ(J)*RHOW(J+1)*W9(I,J+1)*TA(I,J) + T(I,J+1)
        W(I,J+1)  = -DTLZ(J)*RHOW(J+1)*(T9(I  ,J+1)-T9(I  ,J  ))*TA(I,J) &
             + W(I,J+1)
        T(I-1,J)  =  DTLR(I)*R(I  )*U9(I  ,J)*TA(I,J) + T(I-1,J)
        T(I  ,J)  = -DTLR(I)*R(I  )*U9(I  ,J)*TA(I,J) + T(I  ,J)
        U(I  ,J)  = -DTLR(I)*R(I  )*(T9(I  ,J  )-T9(I-1,J  ))*TA(I,J) + U(I  ,J)
        T(I  ,J)  =  DTLR(I)*R(I+1)*U9(I+1,J)*TA(I,J) + T(I  ,J)
        T(I+1,J)  = -DTLR(I)*R(I+1)*U9(I+1,J)*TA(I,J) + T(I+1,J)
        U(I+1,J)  = -DTLR(I)*R(I+1)*(T9(I+1,J  )-T9(I  ,J  ))*TA(I,J) + U(I+1,J)

        ! VA(I,J) = 
        U(I  ,J) = -DTL2*(F+ V9(I,J)/RS(I))/RS(I)*R(I  )*VA(I,J) + U(I  ,J)
        U(I+1,J) = -DTL2*(F+ V9(I,J)/RS(I))/RS(I)*R(I+1)*VA(I,J) + U(I+1,J)
        V(I  ,J) = -DTL2*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)/RS(I)*VA(I,J) &
             + V(I,J)
        V(I,J-1) =  DTLZ(J)*RHOW(J  )*W9(I,J  )*VA(I,J) + V(I,J-1)
        V(I,J  ) = -DTLZ(J)*RHOW(J  )*W9(I,J  )*VA(I,J) + V(I,J  )
        W(I,J  ) = -DTLZ(J)*RHOW(J  )*(V9(I  ,J  )-V9(I  ,J-1))*VA(I,J) &
             + W(I,J  )
        V(I,J  ) =  DTLZ(J)*RHOW(J+1)*W9(I,J+1)*VA(I,J) + V(I,J  )
        V(I,J+1) = -DTLZ(J)*RHOW(J+1)*W9(I,J+1)*VA(I,J) + V(I,J+1)
        W(I,J+1) = -DTLZ(J)*RHOW(J+1)*(V9(I  ,J+1)-V9(I  ,J  ))*VA(I,J) &
             + W(I,J+1)
        V(I-1,J) =  DTLR(I)*R(I  )*U9(I  ,J)*VA(I,J) + V(I-1,J)
        V(I  ,J) = -DTLR(I)*R(I  )*U9(I  ,J)*VA(I,J) + V(I  ,J)
        U(I  ,J) = -DTLR(I)*R(I  )*(V9(I  ,J  )-V9(I-1,J  ))*VA(I,J) + U(I  ,J)
        V(I  ,J) =  DTLR(I)*R(I+1)*U9(I+1,J)*VA(I,J) + V(I  ,J)
        V(I+1,J) = -DTLR(I)*R(I+1)*U9(I+1,J)*VA(I,J) + V(I+1,J)
        U(I+1,J) = -DTLR(I)*R(I+1)*(V9(I+1,J  )-V9(I  ,J  ))*VA(I,J) + U(I+1,J)
     END DO
  END DO

  ! FORCING FOR W EQUATION --> OUTER BOUNDARY

  I=M
  UBW=0
  DO J=N,2,-1
     IF(UBW9(J)>0) THEN
        W1(I-1,J)= DTSRW(J)*(RHOT(J)*(U9(MP1,J)+U9(M,J))+RHOT(J-1)*2.*U9(M,J-1))&
             *WA(I,J)+W1(I-1,J)
        W1(I  ,J)=-DTSRW(J)*(RHOT(J)*(U9(MP1,J)+U9(M,J))+RHOT(J-1)*2.*U9(M,J-1))&
             *WA(I,J)+W1(I  ,J)
        U(M,J-1) =-DTSRW(J)*(W19(I  ,J) - W19(I-1,J))*RHOT(J-1)*WA(I,J)*2.+U(M,J-1)
        U(M,J  ) =-DTSRW(J)*(W19(I  ,J) - W19(I-1,J))*RHOT(J  )*WA(I,J)   +U(M,J )
        U(MP1,J) =-DTSRW(J)*(W19(I  ,J) - W19(I-1,J))*RHOT(J  )*WA(I,J)   +U(MP1,J)
     END IF
     QL(I,J-1) = -DTSG * WA(I,J) + QL(I,J-1)
     QL(I,J  ) = -DTSG * WA(I,J) + QL(I,J  )
     QV(I,J-1) =  DTS*RQVBW(J)*WA(I,J) + QV(I,J-1)
     QV(I,J  ) =  DTS*RQVBW(J)*WA(I,J) + QV(I,J  )
     T( I,J-1) =  DTS*RTBW( J)*WA(I,J) + T( I,J-1)
     T( I,J  ) =  DTS*RTBW( J)*WA(I,J) + T( I,J  )
     W( I,J-1) =  DTSZW(J)*(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  ))*WA(I,J) &
          + W(I,J-1)
     W( I,J  ) = -DTSZW(J)*(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  ))*WA(I,J) &
          + W(I,J  )
     W( I,J  ) = -DTSZW(J)*(W9(I,J  )-W9(I,J-1))*RHOW(J  )*WA(I,J) + W(I,J  )
     W( I,J-1) = -DTSZW(J)*(W9(I,J  )-W9(I,J-1))*RHOW(J-1)*WA(I,J) + W(I,J-1)
     W( I,J  ) =  DTSZW(J)*(RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  ))*WA(I,J) &
          + W(I,J  )
     W( I,J+1) = -DTSZW(J)*(RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  ))*WA(I,J) &
          + W(I,J+1)
     W( I,J  ) = -DTSZW(J)*(W9(I,J+1)-W9(I,J  ))*RHOW(J  )*WA(I,J) + W(I,J  )
     W( I,J+1) = -DTSZW(J)*(W9(I,J+1)-W9(I,J  ))*RHOW(J+1)*WA(I,J) + W(I,J+1)
  END DO

  ! FORCING FOR W EQUATION
  DO J=N, 2, -1
     DO I=MM1, 1, -1
        ! WA(I,J) = 
        QL(I,J-1) = -DTSG*WA(I,J) + QL(I,J-1)
        QL(I,J  ) = -DTSG*WA(I,J) + QL(I,J  )
        QV(I,J-1) =  DTS*RQVBW(J)*WA(I,J) + QV(I,J-1)
        QV(I,J  ) =  DTS*RQVBW(J)*WA(I,J) + QV(I,J  )
        T( I,J-1) =  DTS*RTBW( J)*WA(I,J) + T( I,J-1)
        T( I,J  ) =  DTS*RTBW( J)*WA(I,J) + T( I,J  )
        W( I,J-1) =  DTSZW(J)*(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  )) &
             *WA(I,J) + W(I,J-1)
        W( I,J  ) = -DTSZW(J)*(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  )) &
             *WA(I,J) + W(I,J  )
        W( I,J  ) = -DTSZW(J)*(W9(I,J  )-W9(I,J-1))*RHOW(J  )*WA(I,J) + W(I,J  )
        W( I,J-1) = -DTSZW(J)*(W9(I,J  )-W9(I,J-1))*RHOW(J-1)*WA(I,J) + W(I,J-1)
        W( I,J  ) =  DTSZW(J)*(RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  )) &
             *WA(I,J) + W(I,J  )
        W( I,J+1) = -DTSZW(J)*(RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  )) &
             *WA(I,J) + W(I,J+1)
        W( I,J  ) = -DTSZW(J)*(W9(I,J+1)-W9(I,J  ))*RHOW(J  )*WA(I,J) + W(I,J  )
        W( I,J+1) = -DTSZW(J)*(W9(I,J+1)-W9(I,J  ))*RHOW(J+1)*WA(I,J) + W(I,J+1)

        W( I-1,J)   =  DTSRW(J)*(RHOT(J  )*U9(I  ,J)+RHOT(J-1)*U9(I  ,J-1)) &
             *WA(I,J) + W(I-1,J)
        W( I  ,J)   = -DTSRW(J)*(RHOT(J  )*U9(I  ,J)+RHOT(J-1)*U9(I  ,J-1)) &
             *WA(I,J) + W(I  ,J)
        U( I  ,J-1) = -DTSRW(J)*(W9(I  ,J)-W9(I-1,J))*RHOT(J-1)*WA(I,J) &
             + U(I  ,J-1)
        U( I  ,J  ) = -DTSRW(J)*(W9(I  ,J)-W9(I-1,J))*RHOT(J  )*WA(I,J) &
             + U(I  ,J  )
        W( I  ,J)   =  DTSRW(J)*(RHOT(J  )*U9(I+1,J)+RHOT(J-1)*U9(I+1,J-1)) &
             *WA(I,J) + W(I  ,J)
        W( I+1,J)   = -DTSRW(J)*(RHOT(J  )*U9(I+1,J)+RHOT(J-1)*U9(I+1,J-1)) &
             *WA(I,J) + W(I+1,J)
        U( I+1,J-1) = -DTSRW(J)*(W9(I+1,J)-W9(I  ,J))*RHOT(J-1)*WA(I,J) &
             + U(I+1,J-1)
        U( I+1,J  ) = -DTSRW(J)*(W9(I+1,J)-W9(I  ,J))*RHOT(J  )*WA(I,J) &
             + U(I+1,J  )
     END DO
  END DO

  ! FORCING U EQUATIONS --> OUTER BOUNDARY

  ! DRAG AT J=1

  ! UA(MP1,1) =
  IF(U19(MP1,1)**2 + V19(M,1)**2>0) THEN
     V1(M  ,1) = -CD*DTL*RDZ*U19(MP1,1)*0.5/SQRT(U19(MP1,1)**2 + V19(M,1)**2) &
          *2*V19(M  ,1)*UA(MP1,1) + V1(M  ,1)
     U1(MP1,1) = -CD*DTL*RDZ*U19(MP1,1)*0.5/SQRT(U19(MP1,1)**2 + V19(M,1)**2) &
          *2*U19(MP1,1)*UA(MP1,1) + U1(MP1,1)
     U1(MP1,1) = -CD*DTL*RDZ*SQRT(U19(MP1,1)**2 + V19(M,1)**2)*UA(MP1,1) + U1(MP1,1)
  END IF

  DO J=N, 1, -1
     IF(U9(MP1,J) + CSTAR > 0) THEN
        U1(M  ,J) =  (U9(MP1,J)+CSTAR)*DTL*RDR*UA(MP1,J) + U1(M  ,J)
        U1(MP1,J) = -(U9(MP1,J)+CSTAR)*DTL*RDR*UA(MP1,J) + U1(MP1,J)
        U( MP1,J) = -DTL*RDR*(U19(MP1,J)-U19(M,J))*UA(MP1,J) + U(MP1,J)
     END IF
     ! UA(MP1,J) =
     V(M,J) = DTL*F*UA(MP1,J) + V(M,J)
     V(M,J) = DTL*2*V9(M,J)/RS(M) * UA(MP1,J) + V(M,J)
     UA(MP1,J) = 0
  END DO

  ! FORCING FOR U EQUATION
  
  DO J=N, 1, -1
     DO I=M, 2, -1
        ! UA(I,J) = 
        V(I-1,J)   =  DTSF * UA(I,J) + V(I-1,J)
        V(I  ,J)   =  DTSF * UA(I,J) + V(I  ,J)
        V(I-1,J)   =  DTSV(I-1)*2*V9(I-1,J)*UA(I,J) + V(I-1,J)
        V(I  ,J)   =  DTSV(I  )*2*V9(I  ,J)*UA(I,J) + V(I  ,J)
        U(I,J-1)   =  DTSZ(I,J)*RHOW(J  )*(RS(I)*W9(I,J  )+RS(I-1)*W9(I-1,J  )) &
             *UA(I,J) + U(I,J-1)
        U(I,J  )   = -DTSZ(I,J)*RHOW(J  )*(RS(I)*W9(I,J  )+RS(I-1)*W9(I-1,J  )) &
             *UA(I,J) + U(I,J  )
        W(I-1,J  ) = -DTSZ(I,J)*RHOW(J  )*(U9(I,J)-U9(I,J-1))*RS(I-1)*UA(I,J) &
             + W(I-1,J  )
        W(I  ,J  ) = -DTSZ(I,J)*RHOW(J  )*(U9(I,J)-U9(I,J-1))*RS(I  )*UA(I,J) &
             + W(I  ,J  )
        U(I,J  )   =  DTSZ(I,J)*RHOW(J+1)*(RS(I)*W9(I,J+1)+RS(I-1)*W9(I-1,J+1)) &
             *UA(I,J) + U(I,J  )
        U(I,J+1)   = -DTSZ(I,J)*RHOW(J+1)*(RS(I)*W9(I,J+1)+RS(I-1)*W9(I-1,J+1)) &
             *UA(I,J) + U(I,J+1)
        W(I-1,J+1) = -DTSZ(I,J)*RHOW(J+1)*(U9(I,J+1)-U9(I,J))*RS(I-1)*UA(I,J) &
             + W(I-1,J+1)
        W(I  ,J+1) = -DTSZ(I,J)*RHOW(J+1)*(U9(I,J+1)-U9(I,J))*RS(I  )*UA(I,J) &
             + W(I  ,J+1)
        U(I-1,J)   =  DTSR(I)*(R(I  )*U9(I  ,J)+R(I-1)*U9(I-1,J))*UA(I,J) &
             + U(I-1,J)
        U(I  ,J)   = -DTSR(I)*(R(I  )*U9(I  ,J)+R(I-1)*U9(I-1,J))*UA(I,J) &
             + U(I  ,J)
        U(I-1,J)   = -DTSR(I)*(U9(I  ,J)-U9(I-1,J))*R(I-1)*UA(I,J) + U(I-1,J)
        U(I  ,J)   = -DTSR(I)*(U9(I  ,J)-U9(I-1,J))*R(I  )*UA(I,J) + U(I  ,J)
        U(I  ,J)   =  DTSR(I)*(R(I+1)*U9(I+1,J)+R(I  )*U9(I  ,J))*UA(I,J) &
             + U(I  ,J)
        U(I+1,J)   = -DTSR(I)*(R(I+1)*U9(I+1,J)+R(I  )*U9(I  ,J))*UA(I,J) &
             + U(I+1,J)
        U(I  ,J)   = -DTSR(I)*(U9(I+1,J)-U9(I  ,J))*R(I  )*UA(I,J) + U(I  ,J)
        U(I+1,J)   = -DTSR(I)*(U9(I+1,J)-U9(I  ,J))*R(I+1)*UA(I,J) + U(I+1,J)
     END DO
  END DO

  ! RADIATION R < 1 DEG/DAY
  TDIF=0
  DO J=N, 1, -1
     DO I=M, 1, -1
        IF(QRAD9(I,J) <-RADMAX) QRAD(I,J) = 0
        IF(QRAD9(I,J) > RADMAX) QRAD(I,J) = 0
        TDIF    = -DTL*RADT*QRAD(I,J) + TDIF
        QRAD(I,J) = 0
        T1(I,J) =  TDIF + T1(I,J)
        TDIF    = 0
     END DO
  END DO

  CALL DIFFUSE_ADJ

  DO I=M, 1, -1
     CERS(I)   = CD/CE * CDRS(I) + CERS(I)
     CDRS(I)   = 0
     IF(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2>0) THEN
        V1(I,1)   = 0.5*CD1/SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             *2.0*V19(I,1)*CERS(I) + V1(I,1)
        U1(I+1,1) = 0.5*CD1/SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             *0.25*2.0*(U19(I,1)+U19(I+1,1))*CERS(I) + U1(I+1,1)
        U1(I  ,1) = 0.5*CD1/SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             *0.25*2.0*(U19(I,1)+U19(I+1,1))*CERS(I) + U1(I  ,1)
     END IF
     CERS(I)   = 0
  END DO

  ! SURFACE TEMP AND WATER-VAPOR MIXING RATIO
  DO I=M, 1, -1
     ! ----- forward computing -----
     PNSURF9 = (P19(I,1)-P19(M,1))/(1.0-0.5*G*DZ/(CP*SST(I)) )
     PDST9   = 1000.0 * (PNS+PNSURF9)**(1.0/XKAPPA)
     ! ----- forward computing ----- end
     
     PNSURF   = -SST(I)/(PNS+PNSURF9)**2 * TSURF(I)
     TSURF(I) = 0
     PDST     = -0.622*ESS(I) /(PDST9-ESS(I))**2 * QSURF(I)
     QSURF(I) = 0
     PNSURF   =  1000.0 * (1.0/XKAPPA)*(PNS+PNSURF9)**(1.0/XKAPPA-1.0)*PDST + PNSURF
     PDST     = 0
     P1(M,1)  = -1/(1.0-0.5*G*DZ/(CP*SST(I)) )*PNSURF + P1(M,1)
     P1(I,1)  =  1/(1.0-0.5*G*DZ/(CP*SST(I)) )*PNSURF + P1(I,1)
     PNSURF   = 0
  END DO

  ! ----- restoring state variables -----
  U9 = inU9
  V9 = inV9
  W9 = inW9
  P9 = inP9
  T9 = inT9
  QV9 = inQV9
  QL9 = inQL9
  U19 = inU19
  V19 = inV19
  W19 = inW19
  T19 = inT19
  P19 = inP19
  QV19 = inQV19
  QL19 = inQL19
  ! ----- restoring state variables ----- end

END SUBROUTINE hurr_adj

! ==================================================================

SUBROUTINE DIFFUSE_ADJ
  USE hurr_vars
  IMPLICIT NONE

  ! ----- LOCAL DEPENDENT VARIABLES -----
  REAL*8 :: &
       TRR9(M  ,N  ), TTT9(M  ,N  ), TZZ9(M  ,NP1), TRZ9(MP1,NP1), &
       TRT9(MP1,N  ), TZT9(M  ,NP1), TR9( MP1,N  ), TZ9( M  ,NP1), &
       TAD9(M  ,N  ), QVR9(MP1,N  ), QVZ9(M  ,NP1), QLR9(MP1,N  ), &
       QLZ9(M  ,NP1), XKM9(M  ,NP1), UDR9(M  ,N  ), WDZ9(M  ,N  ), &
       VDR9(MP1,N  ), VDZ9(M  ,N  ), TE9( M  ,N  ), &
       XKMH9(M ,NP1), UZWR9(MP1, N), DEFH29(M,N  ), DEF29(M,N), DTDZ9(M,N)
  REAL*8 :: &
       TRR(M  ,N  ), TTT(M  ,N  ), TZZ(M  ,NP1), TRZ(MP1,NP1), &
       TRT(MP1,N  ), TZT(M  ,NP1), TR( MP1,N  ), TZ( M  ,NP1), &
       TAD(M  ,N  ), QVR(MP1,N  ), QVZ(M  ,NP1), QLR(MP1,N  ), &
       QLZ(M  ,NP1), XKM(M  ,NP1), UDR(M  ,N  ), WDZ(M  ,N  ), &
       VDR(MP1,N  ), VDZ(M  ,N  ), TE( M  ,N  ), &
       XKMH(M ,NP1), UZWR(MP1, N), DEFH2(M,N  ), DEF2(M,N), DTDZ(M,N)
  
  REAL*8 :: TT9, AA9, VABS9
  REAL*8 :: TT, AA, VABS
  ! -------------------------------------

!  the common block prototype, which used to be shared with the MAIN
!  COMMON/DIFF/ &
!       UA(MP1,N) , WA(M,NP1), VA(M,N)  , TA(M,N), QVA(M,N), QLA(M,N), &
!       U1(MP1,N) , V1(M,N)  , W1(M,NP1), T1(M,N), QV1(M,N), QL1(M,N), &
!       XKM(M,NP1), P1(M,N)  , R(MP1) ,RS(M), Z(NP1),ZS(N),TAD(M,N)  , &
!       RTBW(N)   , RQVBW(N) , TAU(N  ) , TB(N), QVB(N), PN(N), PD(N), &
!       TSURF(M)  , QSURF(M) , RDR, RDZ, RDR2, RDZ2 , DTS, DTL, DZ   , &
!       XVL2, XHL2, CD, CE, XLDCP, XKAPPA, A1, G, CERS(M), CDRS(M)

  INCLUDE 'diffuse-fwd.inc'
  
  !  DISSIPATIVE HEATING
  TAD =0
  UDR =0
  XKMH=0
  VDR =0
  CDRS=0
  UZWR=0
  XKM =0
  VDZ =0
  TZT =0

  DO I=M-1, 1, -1
     ! TA(I,1) = 
     TAD(I,1) = TA(I,1) + TAD(I,1)

     ! TAD(I,1) = 
     UDR(I,2)   = DTL/(PN(1)*1005.0)*XKMH9(I,2)*0.25 * (UDR9(I,1)+UDR9(I,2))*2 &
          *TAD(I,1) + UDR(I,2)
     UDR(I,1)   = DTL/(PN(1)*1005.0)*XKMH9(I,2)*0.25 * (UDR9(I,1)+UDR9(I,2))*2 &
          *TAD(I,1) + UDR(I,1)
     XKMH( I,2) = DTL/(PN(1)*1005.0)*0.25 * (UDR9(I,1)+UDR9(I,2))**2 &
          *TAD(I,1) + XKMH(I,2)
     VDR(I+1,2) = DTL/(PN(1)*1005.0) *XKMH9(I,2)*(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))*2*TAD(I,1) + VDR(I+1,2)
     VDR(I  ,2) = DTL/(PN(1)*1005.0) *XKMH9(I,2)*(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))*2*TAD(I,1) + VDR(I  ,2)
     VDR(I+1,1) = DTL/(PN(1)*1005.0) *XKMH9(I,2)*(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))*2*TAD(I,1) + VDR(I+1,1)
     VDR(I  ,1) = DTL/(PN(1)*1005.0) *XKMH9(I,2)*(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))*2*TAD(I,1) + VDR(I  ,1)
     XKMH( I,2) = DTL/(PN(1)*1005.0) *(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))**2*TAD(I,1) + XKMH(I,2)
     V1(I,1)    = DTL/(PN(1)*1005.0) &
          *RDZ*CDRS9(I)*(0.25*(U19(I,1)+U19(I+1,1))**2 + V19(I,1)**2)**0.5*1.5 &
          *2*V19(I,1)*TAD(I,1) + V1(I,1)
     U1(I+1,1)  = DTL/(PN(1)*1005.0) &
          *RDZ*CDRS9(I)*(0.25*(U19(I,1)+U19(I+1,1))**2 + V19(I,1)**2)**0.5*1.5 &
          *0.25*2*(U19(I,1)+U19(I+1,1))*TAD(I,1) + U1(I+1,1)
     U1(I  ,1)  = DTL/(PN(1)*1005.0) &
          *RDZ*CDRS9(I)*(0.25*(U19(I,1)+U19(I+1,1))**2 + V19(I,1)**2)**0.5*1.5 &
          *0.25*2*(U19(I,1)+U19(I+1,1))*TAD(I,1) + U1(I  ,1)
     CDRS(I)    = DTL /(PN(1)*1005.0)*RDZ &
          *(0.25*(U19(I,1)+U19(I+1,1))**2 + V19(I,1)**2)**1.5*TAD(I,1) + CDRS(I)
     DO J=N-1, 2, -1
        ! TA(I,J) =
        TAD(I,J) = TA(I,J) + TAD(I,J)

        ! TAD(I,J) =
        UDR(I,J+1)  = DTL/(PN(J)*1005.0)* XKMH9(I,J+1) &
             *0.25*(UDR9(I,J)+UDR9(I,J+1))*2*TAD(I,J) + UDR(I,J+1)
        UDR(I,J  )  = DTL/(PN(J)*1005.0)* XKMH9(I,J+1) &
             *0.25*(UDR9(I,J)+UDR9(I,J+1))*2*TAD(I,J) + UDR(I,J  )
        XKMH(I,J+1) = DTL/(PN(J)*1005.0)* 0.25*(UDR9(I,J)+UDR9(I,J+1))**2 &
             *TAD(I,J) + XKMH(I,J+1)
        VDR(I+1,J+1)= DTL/(PN(J)*1005.0)* XKMH9(I,J+1)*(1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))*2 &
             *TAD(I,J) + VDR(I+1,J+1)
        VDR(I  ,J+1)= DTL/(PN(J)*1005.0)* XKMH9(I,J+1)*(1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))*2 &
             *TAD(I,J) + VDR(I  ,J+1)
        VDR(I+1,J  )= DTL/(PN(J)*1005.0)* XKMH9(I,J+1)*(1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))*2 &
             *TAD(I,J) + VDR(I+1,J  )
        VDR(I  ,J  )= DTL/(PN(J)*1005.0)* XKMH9(I,J+1)*(1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))*2 &
             *TAD(I,J) + VDR(I  ,J  )
        XKMH(I,J+1) = DTL/(PN(J)*1005.0)* (1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))**2 &
             *TAD(I,J) + XKMH(I,J+1)
        UZWR(I+1,J) = DTL/(PN(J)*1005.0)* XKM9(I,J)*0.25* &
             (UZWR9(I,J) + UZWR9(I+1,J))*2*TAD(I,J) + UZWR(I+1,J)
        UZWR(I  ,J) = DTL/(PN(J)*1005.0)* XKM9(I,J)*0.25* &
             (UZWR9(I,J) + UZWR9(I+1,J))*2*TAD(I,J) + UZWR(I  ,J)
        XKM(I,J)    = DTL/(PN(J)*1005.0)* 0.25*(UZWR9(I,J) + UZWR9(I+1,J))**2 &
             *TAD(I,J) + XKM(I,J)
        VDZ(I,J+1)  = DTL/(PN(J)*1005.0)* TZT9(I,J+1)*TAD(I,J) + VDZ(I,J+1)
        TZT(I,J+1)  = DTL/(PN(J)*1005.0)* VDZ9(I,J+1)*TAD(I,J) + TZT(I,J+1)
     END DO
  END DO

  QLZ=0
  QLR=0
  QVZ=0
  QVR=0
  TZ =0
  TR =0
  TRT=0
  DO I=M, 1, -1
     J=N

     ! QLA(I,J) = 
     QL1(I,J)   =  DTL*TAU(J)*QLA(I,J) + QL1(I,J)
     QLZ(I,J)   = -DTL*RDZ* QLA(I,J) + QLZ(I,J  )
     QLZ(I,J+1) =  DTL*RDZ* QLA(I,J) + QLZ(I,J+1)
     QLR(I  ,J) = -DTL*RDR/RS(I) * R(I  )*QLA(I,J) + QLR(I  ,J)
     QLR(I+1,J) =  DTL*RDR/RS(I) * R(I+1)*QLA(I,J) + QLR(I+1,J)

     ! QVA(I,J) =
     QV1(I,J  ) =  DTL*TAU(J)*QVA(I,J) + QV1(I,J)
     QVZ(I,J  ) = -DTL*RDZ/(PD(J)/PN(J))*QVA(I,J) + QVZ(I,J  )
     QVZ(I,J+1) =  DTL*RDZ/(PD(J)/PN(J))*QVA(I,J) + QVZ(I,J+1)
     QVR(I  ,J) = -DTL*RDR/RS(I)*R(I  )* QVA(I,J) + QVR(I  ,J)
     QVR(I+1,J) =  DTL*RDR/RS(I)*R(I+1)* QVA(I,J) + QVR(I+1,J)

     ! TA(I,J) =
     T1(I,J)   =  DTL*TAU(J) * TA(I,J) + T1(I,J)
     TZ(I,J  ) = -DTL*RDZ/PD(J) * TA(I,J) + TZ(I,J  )
     TZ(I,J+1) =  DTL*RDZ/PD(J) * TA(I,J) + TZ(I,J+1)
     TR(I  ,J) = -DTL*RDR/RS(I)*R(I  )*TA(I,J) + TR(I  ,J)
     TR(I+1,J) =  DTL*RDR/RS(I)*R(I+1)*TA(I,J) + TR(I+1,J)

     ! VA(I,J) =
     V1(I,J)    =  DTL*TAU(J) * VA(I,J) + V1(I,J)
     TZT(I,J  ) = -DTL*RDZ*VA(I,J) + TZT(I,J  )
     TZT(I,J+1) =  DTL*RDZ*VA(I,J) + TZT(I,J+1)
     TRT(I  ,J) = -DTL*RDR/(RS(I)**2)*R(I  )**2*VA(I,J) + TRT(I  ,J)
     TRT(I+1,J) =  DTL*RDR/(RS(I)**2)*R(I+1)**2*VA(I,J) + TRT(I+1,J)
     DO J=N-1, 1, -1
        ! QLA(I,J) =
        QL1(I,J  ) =  DTL*TAU(J) * QLA(I,J) + QL1(I,J)
        QLZ(I,J  ) = -DTL*RDZ * QLA(I,J) + QLZ(I,J)
        QLZ(I,J+1) =  DTL*RDZ * QLA(I,J) + QLZ(I,J+1)
        QLR(I  ,J) = -DTL*RDR/RS(I)*R(I  )*QLA(I,J) + QLR(I  ,J)
        QLR(I+1,J) =  DTL*RDR/RS(I)*R(I+1)*QLA(I,J) + QLR(I+1,J)
        
        ! QVA(I,J) =
        QV1(I,J  ) =  DTL*TAU(J)*QVA(I,J) + QV1(I,J)
        QVZ(I,J  ) = -DTL*RDZ/(0.5*(PD(J)/PN(J)+PD(J+1)/PN(J+1)))*QVA(I,J) &
             + QVZ(I,J  )
        QVZ(I,J+1) =  DTL*RDZ/(0.5*(PD(J)/PN(J)+PD(J+1)/PN(J+1)))*QVA(I,J) &
             + QVZ(I,J+1)
        QVR(I  ,J) = -DTL*RDR/RS(I)*R(I  )*QVA(I,J) + QVR(I  ,J)
        QVR(I+1,J) =  DTL*RDR/RS(I)*R(I+1)*QVA(I,J) + QVR(I+1,J)
        
        ! TA(I,J)  =
        T1(I,J  ) =  DTL*TAU(J)*TA(I,J) + T1(I,J)
        TZ(I,J  ) = -DTL*RDZ/(0.5*(PD(J)+PD(J+1)))*TA(I,J) + TZ(I,J  )
        TZ(I,J+1) =  DTL*RDZ/(0.5*(PD(J)+PD(J+1)))*TA(I,J) + TZ(I,J+1)
        TR(I  ,J) = -DTL*RDR/RS(I)*R(I  ) * TA(I,J) + TR(I  ,J)
        TR(I+1,J) =  DTL*RDR/RS(I)*R(I+1) * TA(I,J) + TR(I+1,J)
        
        ! VA(I,J)  =
        V1( I,J  ) =  DTL*TAU(J) * VA(I,J) + V1(I,J  )
        TZT(I,J  ) = -DTL*RDZ * VA(I,J) + TZT(I,J  )
        TZT(I,J+1) =  DTL*RDZ * VA(I,J) + TZT(I,J+1)
        TRT(I  ,J) = -DTL*RDR/(RS(I)**2)*R(I  )**2*VA(I,J) + TRT(I  ,J)
        TRT(I+1,J) =  DTL*RDR/(RS(I)**2)*R(I+1)**2*VA(I,J) + TRT(I+1,J)
     END DO
  END DO

  TZZ=0
  TRZ=0
  DO J=N, 2, -1
     DO I=M, 1, -1
        ! WA(I,J) =
        W1( I,J)   =  DTS*0.5*(TAU(J) + TAU(J-1))*WA(I,J) + W1(I,J)
        TZZ(I,J-1) = -DTS*RDZ*WA(I,J) + TZZ(I,J-1)
        TZZ(I,J  ) =  DTS*RDZ*WA(I,J) + TZZ(I,J  )
        TRZ(I  ,J) = -DTS*RDR/RS(I)*R(I  )*WA(I,J) + TRZ(I  ,J)
        TRZ(I+1,J) =  DTS*RDR/RS(I)*R(I+1)*WA(I,J) + TRZ(I+1,J)
     END DO
  END DO

  TTT=0
  TRR=0
  DO J=N, 1, -1
     DO I=M, 2, -1
        ! UA(I,J) = 
        U1( I,J)   =  DTS*TAU(J) * UA(I,J) + U1(I,J)
        TRZ(I,J  ) = -DTS*RDZ*UA(I,J) + TRZ(I,J  )
        TRZ(I,J+1) =  DTS*RDZ*UA(I,J) + TRZ(I,J+1)
        TTT(I,J)   = -DTS/R(I)*UA(I,J) + TTT(I,J)
        TRR(I-1,J) = -DTS*RDR/R(I)*RS(I-1)*UA(I,J) + TRR(I-1,J)
        TRR(I  ,J) =  DTS*RDR/R(I)*RS(I  )*UA(I,J) + TRR(I  ,J)
     END DO
  END DO

  ! QL FLUX
  !    QLZ(M,NP1)
  DO I=M, 1, -1
     QLZ(I,NP1) = 0
     QLZ(I,1  ) = 0
     DO J=N, 2, -1
        ! QLZ(I,J) =
        QL1(I,J-1) = -XKM9(I,J)*RDZ * QLZ(I,J) + QL1(I,J-1)
        QL1(I,J  ) =  XKM9(I,J)*RDZ * QLZ(I,J) + QL1(I,J  )
        XKM(I,J)   =  RDZ*(QL19(I,J)-QL19(I,J-1)) * QLZ(I,J) + XKM(I,J)
     END DO
  END DO

  !    QLR(IMP1,N)
  DO J=N, 1, -1
     QLR(M,J) = R(M)/R(MP1) *QLR(MP1,J) + QLR(M,J)
     QLR(1,J) = 0
     DO I=M, 2, -1
        ! QLR(I,J) =
        XKMH(I-1,J  ) = 0.25 *RDR*(QL19(I,J)-QL19(I-1,J))*QLR(I,J) &
             + XKMH(I-1,J  )
        XKMH(I-1,J+1) = 0.25 *RDR*(QL19(I,J)-QL19(I-1,J))*QLR(I,J) &
             + XKMH(I-1,J+1)
        XKMH(I  ,J  ) = 0.25 *RDR*(QL19(I,J)-QL19(I-1,J))*QLR(I,J) &
             + XKMH(I  ,J  )
        XKMH(I  ,J+1) = 0.25 *RDR*(QL19(I,J)-QL19(I-1,J))*QLR(I,J) &
             + XKMH(I  ,J+1)
        QL1(I-1,J)    = -0.25*RDR &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))*QLR(I,J) &
             + QL1(I-1,J)
        QL1(I  ,J)    =  0.25*RDR &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))*QLR(I,J) &
             + QL1(I  ,J)
     END DO
  END DO

  ! QV FLUX
  !     QVZ(M,NP1)
  CERS=0
  QSURF=0
  DO I=M, 1, -1
     ! forward computing
     VABS9=SQRT(0.25*(U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2)
     ! end forward computing
     
     QVZ(I,NP1) =  0
     QVZ(I,1  ) =  PD(1)/PN(1) * QVZ(I,1)

     VABS       =  (QV19(I,1)-QSURF9(I))*CERS9(I)*QVZ(I,1)
     CERS(I)    =  (QV19(I,1)-QSURF9(I))*VABS9*QVZ(I,1) + CERS(I)
     QSURF(I)   = -CERS9(I)*VABS9 * QVZ(I,1) + QSURF(I)
     QV1(I,1)   =  CERS9(I)*VABS9 * QVZ(I,1) + QV1(I,1)

     ! VABS =
     V1(I  ,1) = 0.5/SQRT(0.25*(U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2) &
          *2*V19(I,1)*VABS + V1(I,1)
     U1(I  ,1) = 0.5/SQRT(0.25*(U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2) &
          *0.25*2*(U19(I+1,1)+U19(I,1))*VABS + U1(I  ,1)
     U1(I+1,1) = 0.5/SQRT(0.25*(U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2) &
          *0.25*2*(U19(I+1,1)+U19(I,1))*VABS + U1(I+1,1)
     DO J=N, 2, -1
        ! QVZ(I,J) =
        QVZ(I,J  ) =  PD(J)/PN(J) * QVZ(I,J)
        QV1(I,J-1) = -XKM9(I,J)*RDZ* QVZ(I,J) + QV1(I,J-1)
        QV1(I,J  ) =  XKM9(I,J)*RDZ* QVZ(I,J) + QV1(I,J  )
        XKM(I,J  ) =  RDZ*(QV19(I,J)-QV19(I,J-1))* QVZ(I,J) + XKM(I,J)
     END DO
  END DO

  !     QVR(MP1,N)
  DO J=N, 1, -1
     QVR(M,J) = R(M)/R(MP1)*QVR(MP1,J) + QVR(M,J)
     QVR(1,J) = 0
     DO I=M, 1, -1
        ! QVR(I,J) =
        XKMH(I-1,J  ) = 0.25 *RDR*(QV19(I,J)-QV19(I-1,J))*QVR(I,J) &
             + XKMH(I-1,J  )
        XKMH(I-1,J+1) = 0.25 *RDR*(QV19(I,J)-QV19(I-1,J))*QVR(I,J) &
             + XKMH(I-1,J+1)
        XKMH(I  ,J  ) = 0.25 *RDR*(QV19(I,J)-QV19(I-1,J))*QVR(I,J) &
             + XKMH(I  ,J  )
        XKMH(I  ,J+1) = 0.25 *RDR*(QV19(I,J)-QV19(I-1,J))*QVR(I,J) &
             + XKMH(I  ,J+1)
        QV1(I-1,J)    = -0.25*RDR &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))*QVR(I,J) &
             + QV1(I-1,J)
        QV1(I  ,J)    =  0.25*RDR &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))*QVR(I,J) &
             + QV1(I  ,J)
     END DO
  END DO

  ! TEMPERATURE FLUX
  !     TZ(M,NP1)
  TSURF=0
  DO I=M, 1, -1
     ! forward calculation
     VABS9=SQRT(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2)
     ! end forward calculation
     
     TZ(I,NP1) = 0
     TZ(I,1  ) = PD(1) * TZ(I,1)
     ! TZ(I,1) = 
     VABS      = (T19(I,1)-TSURF9(I))*CERS9(I)*TZ(I,1)
     CERS(I)   = (T19(I,1)-TSURF9(I))*VABS9 *  TZ(I,1) + CERS(I)
     TSURF(I)  =-CERS9(I)*VABS9               *TZ(I,1) + TSURF(I)
     T1(I,1)   = CERS9(I)*VABS9               *TZ(I,1) + T1(I,1)

     ! VABS =
     IF(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2>0) THEN
        V1(I  ,1) = 0.5/SQRT(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2) &
             *2*V19(I,1)*VABS + V1(I  ,1)
        U1(I  ,1) = 0.5/SQRT(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2) &
             *0.25*2*(U19(I+1,1)+U19(I,1))*VABS + U1(I  ,1)
        U1(I+1,1) = 0.5/SQRT(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2) &
             *0.25*2*(U19(I+1,1)+U19(I,1))*VABS + U1(I+1,1)
     END IF
     DO J=N, 2, -1
        TZ(I,J) = PD(J) * TZ(I,J)
        ! TZ(I,J) =
        T1(I,J-1) = -XKM9(I,J)*RDZ* TZ(I,J) + T1(I,J-1)
        T1(I,J  ) =  XKM9(I,J)*RDZ* TZ(I,J) + T1(I,J  )
        XKM(I,J)  =  RDZ*(T19(I,J)-T19(I,J-1))*TZ(I,J) + XKM(I,J)
     END DO
  END DO

  !     TR(MP1,N)
  DO J=N, 1, -1
     TR(M,J) = R(M)/R(MP1) * TR(MP1,J) + TR(M,J)
     TR(1,J) = 0
     DO I=M, 2, -1
        ! TR(I,J) =
        XKMH(I-1,J  ) = 0.25 * RDR*(T19(I,J)-T19(I-1,J))* TR(I,J) &
             + XKMH(I-1,J  )
        XKMH(I-1,J+1) = 0.25 * RDR*(T19(I,J)-T19(I-1,J))* TR(I,J) &
             + XKMH(I-1,J+1)
        XKMH(I  ,J  ) = 0.25 * RDR*(T19(I,J)-T19(I-1,J))* TR(I,J) &
             + XKMH(I  ,J  )
        XKMH(I  ,J+1) = 0.25 * RDR*(T19(I,J)-T19(I-1,J))* TR(I,J) &
             + XKMH(I  ,J+1)
        T1(I-1,J) = -0.25*RDR*(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))&
             *TR(I,J) + T1(I-1,J)
        T1(I  ,J) =  0.25*RDR*(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))&
             *TR(I,J) + T1(I  ,J)
     END DO
  END DO

  !     TZZ(M,N)
  WDZ=0
  DO J=N, 1, -1
     DO I=M, 1, -1
        ! TZZ(I,J) =
        WDZ(I,J  ) = (XKM9(I,J+1)+XKM9(I,J))*TZZ(I,J) + WDZ(I,J)
        XKM(I,J  ) = WDZ9(I,J)*TZZ(I,J) + XKM(I,J  )
        XKM(I,J+1) = WDZ9(I,J)*TZZ(I,J) + XKM(I,J+1)
     END DO
  END DO

  !     TZT(M,NP1)
  DO I=M, 1, -1
     TZT(I,NP1) = 0
     ! TZT(I,1) = 
     IF(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2>0) THEN
        V1(I  ,1) = CDRS9(I)*V19(I,1)*0.5 &
             /SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             *2*V19(I,1)* TZT(I,1) + V1(I,1)
        U1(I+1,1) = CDRS9(I)*V19(I,1)*0.5 &
             /SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             *0.25*2*(U19(I,1)+U19(I+1,1))* TZT(I,1) + U1(I+1,1)
        U1(I  ,1) = CDRS9(I)*V19(I,1)*0.5 &
             /SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             *0.25*2*(U19(I,1)+U19(I+1,1))* TZT(I,1) + U1(I  ,1)
        V1(I  ,1) = CDRS9(I)*SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             *TZT(I,1) + V1(I,1)
        CDRS(I)   = V19(I,1)*SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             *TZT(I,1) + CDRS(I)
     END IF
     DO J=N, 2, -1
        ! TZT(I,J) =
        VDZ(I,J) = XKM9(I,J)* TZT(I,J) + VDZ(I,J)
        XKM(I,J) = VDZ9(I,J)* TZT(I,J) + XKM(I,J)
     END DO
  END DO

  !     TRT(M,N)
  DO J=N, 1, -1
     TRT(M,J) = R(M)**2/R(MP1)**2 * TRT(MP1,J) + TRT(M,J)
     TRT(1,J) = 0
     DO I=M, 2, -1
        ! TRT(I,J) = 
        XKMH(I-1,J  ) = 0.25*VDR9(I,J)*TRT(I,J) + XKMH(I-1,J  )
        XKMH(I-1,J+1) = 0.25*VDR9(I,J)*TRT(I,J) + XKMH(I-1,J+1)
        XKMH(I  ,J  ) = 0.25*VDR9(I,J)*TRT(I,J) + XKMH(I  ,J  )
        XKMH(I  ,J+1) = 0.25*VDR9(I,J)*TRT(I,J) + XKMH(I  ,J+1)
        VDR(I,J)      = 0.25*(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))&
             *TRT(I,J) + VDR(I,J)
     END DO
  END DO

  !     TRZ(M,NP1)
  DO J=NP1, 1, -1
     ! TRZ(MP1,J) = 
     TRZ(M,J) = R(M)/R(MP1) * TRZ(MP1,J) + TRZ(M,J)
     TRZ(1,J) = 0
  END DO

  DO I=M, 2, -1
     TRZ(I,NP1) = 0
     ! TRZ(I,1) = 
     IF(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2>0) THEN
        V1(I-1,1) = 0.5*(CDRS9(I)+CDRS9(I-1))*U19(I,1)*0.5 &
             /SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
             *0.25*2*(V19(I,1)+V19(I-1,1))*TRZ(I,1) + V1(I-1,1)
        V1(I  ,1) = 0.5*(CDRS9(I)+CDRS9(I-1))*U19(I,1)*0.5 &
             /SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
             *0.25*2*(V19(I,1)+V19(I-1,1))*TRZ(I,1) + V1(I  ,1)
        U1(I  ,1) = 0.5*(CDRS9(I)+CDRS9(I-1))*U19(I,1)*0.5 &
             /SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
             *2*U19(I,1)*TRZ(I,1) + U1(I  ,1)
        U1(I  ,1) = 0.5*(CDRS9(I)+CDRS9(I-1)) &
             *SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
             *TRZ(I,1) + U1(I,1)
        CDRS(I-1) = 0.5*U19(I,1)*SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
             *TRZ(I,1) + CDRS(I-1)
        CDRS(I  ) = 0.5*U19(I,1)*SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
             *TRZ(I,1) + CDRS(I  )
     END IF
  END DO

  DO J=N, 2, -1
     DO I=M, 2, -1
        ! TRZ(I,J) = 
        UZWR(I,J)  = 0.5*(XKM9(I-1,J)+XKM9(I,J))*TRZ(I,J) + UZWR(I,J)
        XKM(I  ,J) = 0.5*UZWR9(I,J)*TRZ(I,J) + XKM(I  ,J)
        XKM(I-1,J) = 0.5*UZWR9(I,J)*TRZ(I,J) + XKM(I-1,J)
     END DO
  END DO

  !     TTT(M,N)
  DO J=N, 1, -1
     DO I=M, 2, -1
        ! TTT(I,J) =
        XKMH(I-1,J  ) = U19(I,J)/R(I)*0.5 * TTT(I,J) + XKMH(I-1,J  )
        XKMH(I-1,J+1) = U19(I,J)/R(I)*0.5 * TTT(I,J) + XKMH(I-1,J+1)
        XKMH(I  ,J  ) = U19(I,J)/R(I)*0.5 * TTT(I,J) + XKMH(I  ,J  )
        XKMH(I  ,J+1) = U19(I,J)/R(I)*0.5 * TTT(I,J) + XKMH(I  ,J+1)
        U1(I,J)       = 0.5/R(I) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))*TTT(I,J) &
             + U1(I,J)
     END DO
     TTT(1,J) = 0
  END DO

  !     TRR(M,N)
  DO J=N, 1, -1
     DO I=M, 1, -1
        ! TRR(I,J) = 0
        UDR(I,J)    = (XKMH9(I,J+1)+XKMH9(I,J))*TRR(I,J) + UDR(I,J)
        XKMH(I,J  ) = UDR9(I,J)*TRR(I,J) + XKMH(I,J  )
        XKMH(I,J+1) = UDR9(I,J)*TRR(I,J) + XKMH(I,J+1)
        TRR(I,J) = 0
     END DO
  END DO

  DO I=M, 1, -1
     XKMH(I,N)   = XKMH(I,NP1)+XKMH(I,N)
     XKMH(I,NP1) = 0
     XKMH(I,2)   = XKMH(I,1  )+XKMH(I,2)
     XKMH(I,1  ) = 0
     XKM( I,N)   = XKM( I,NP1)+XKM( I,N)
     XKM( I,NP1) = 0
     XKM( I,2)   = XKM( I,1  )+XKM( I,2)
     XKM( I,1  ) = 0
  END DO

  DEFH2=0
  DTDZ=0
  DEF2=0
  DO J=N, 2, -1
     DO I=M, 1, -1
        ! ----- fwd -----
        XKM9(I,J) = 0
        IF(DEF29(I,J) > DTDZ9(I,J)) XKM9(I,J)= XVL2 * SQRT(DEF29(I,J) - DTDZ9(I,J))
        XKMH9(I,J) = XHL2 * SQRT(DEFH29(I,J))
        IF(XKM9(I,J) >= 0.4*DZ**2/DTL) XKM9(I,J) = 0.4 *DZ**2/DTL
        IF(XKMH9(I,J) < XKM9(I,J)) XKMH9(I,J) = XKM9(I,J)
        ! ---------------
        
        IF(XKMH9(I,J) < XKM9(I,J)) THEN
           XKM(I,J) = XKMH(I,J)
           XKMH(I,J)= 0
        END IF
        IF(XKM9(I,J) >= 0.4*DZ**2/DTL) XKM(I,J) = 0
        IF(DEFH29(I,J)>0) &
             DEFH2(I,J) = XHL2 * 0.5/SQRT(DEFH29(I,J))*XKMH(I,J) + DEFH2(I,J)
        XKMH(I,J)=0
        IF(DEF29(I,J) > DTDZ9(I,J)) THEN
           ! XKM(I,J) = 
           DTDZ(I,J) = -XVL2*0.5/SQRT(DEF29(I,J)-DTDZ9(I,J))*XKM(I,J) + DTDZ(I,J)
           DEF2(I,J) =  XVL2*0.5/SQRT(DEF29(I,J)-DTDZ9(I,J))*XKM(I,J) + DEF2(I,J)
           XKM(I,J)  =  0
        END IF
        XKM(I,J) = 0
     END DO
  END DO

  TE=0
  DO J=N, 2, -1
     DO I=M, 1, -1
        IF((QL19(I,J) * QL19(I,J-1)) >= 1.0E-8 ) THEN
           ! --- fwd ---
           TT9 = 0.5 &
                *(T19(I,J  ) * ( PN(J  )+P19(I,J  )-P19(M,1)) &
                + T19(I,J-1) * ( PN(J-1)+P19(I,J-1)-P19(M,1)) )
           AA9=(   1 + 4362    *(QV19(I,J)+QV19(I,J-1)) / TT9   ) &
                / (1 + 6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2)) &
                * (2 / (TB(J)+TB(J-1)) )
           ! -----------
           
           ! DTDZ(I,J) = 
           QV1(I,J-1) =  RDZ*G*DTDZ(I,J) + QV1(I,J-1)
           QL1(I,J-1) =  RDZ*G*DTDZ(I,J) + QL1(I,J-1)
           QV1(I,J  ) = -RDZ*G*DTDZ(I,J) + QV1(I,J  )
           QL1(I,J  ) = -RDZ*G*DTDZ(I,J) + QL1(I,J  )
           TE( I,J-1) = -RDZ*G*AA9 * DTDZ(I,J) + TE(I,J-1)
           TE( I,J  ) =  RDZ*G*AA9 * DTDZ(I,J) + TE(I,J  )
           AA         =  RDZ*G*(TE9(I,J)-TE9(I,J-1))*DTDZ(I,J)
           DTDZ(I,J)  =  0

           ! AA =
           TT = ( 1 + 4362  *(QV19(I,J)+QV19(I,J-1))/ TT9   ) &
                /(1+6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2))**2 &
                *(2/(TB(J)+TB(J-1)))*2*6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**3) &
                *AA
           QV1(I,J-1) = &
                -(1 + 4362  *(QV19(I,J)+QV19(I,J-1))/ TT9   ) &
                /(1+6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2))**2 &
                *(2/(TB(J)+TB(J-1)))*6738953/(TT9**2) * AA + QV1(I,J-1)
           QV1(I,J  ) = &
                -(1 + 4362  *(QV19(I,J)+QV19(I,J-1))/ TT9   ) &
                /(1+6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2))**2 &
                *(2/(TB(J)+TB(J-1)))*6738953/(TT9**2) * AA + QV1(I,J  )
           TT =   -(4362   *(QV19(I,J)+QV19(I,J-1))/TT9**2) &
                /(1+6738953*(QV19(I,J)+QV19(I,J-1))/(TT9**2)) &
                *(2/(TB(J)+TB(J-1))) * AA + TT
           QV1(I,J-1) = (4362/TT9) &
                /(1+6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2)) &
                *(2/(TB(J)+TB(J-1))) * AA + QV1(I,J-1)
           QV1(I,J  ) = (4362/TT9) &
                /(1+6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2)) &
                *(2/(TB(J)+TB(J-1))) * AA + QV1(I,J  )

           ! TT
           P1(M,1  ) = -0.5*T19(I,J-1)*TT + P1(M,1  )
           P1(I,J-1) =  0.5*T19(I,J-1)*TT + P1(I,J-1)
           T1(I,J-1) =  0.5*(PN(J-1)+P19(I,J-1)-P19(M,1)) * TT + T1(I,J-1)
           P1(M,1  ) = -0.5*T19(I,J  )*TT + P1(M,1  )
           P1(I,J  ) =  0.5*T19(I,J  )*TT + P1(I,J  )
           T1(I,J  ) =  0.5*(PN(J  )+P19(I,J  )-P19(M,1)) * TT + T1(I,J  )
        END IF
        ! DTDZ(I,J) =
        QV1(I,J-1) = -RDZ*2.0*RQVBW(J)*DTDZ(I,J) + QV1(I,J-1)
        QV1(I,J  ) =  RDZ*2.0*RQVBW(J)*DTDZ(I,J) + QV1(I,J  )
        T1( I,J-1) = -RDZ*2.0*RTBW( J)*DTDZ(I,J) + T1( I,J-1)
        T1( I,J  ) =  RDZ*2.0*RTBW( J)*DTDZ(I,J) + T1( I,J  )
     END DO
  END DO

  DO J=N, 2, -1
     DO I=M, 1, -1
        ! DEF2(I,J) = 
        VDZ(I,J)    = 2*VDZ9(I,J)*DEF2(I,J) + VDZ(I,J)
        UZWR(I  ,J) = 0.5*2*UZWR9(I  ,J)*DEF2(I,J) + UZWR(I  ,J)
        UZWR(I+1,J) = 0.5*2*UZWR9(I+1,J)*DEF2(I,J) + UZWR(I+1,J)
        WDZ( I,J-1) = 2*WDZ9(I,J-1)*DEF2(I,J) + WDZ(I,J-1)
        WDZ( I,J  ) = 2*WDZ9(I,J  )*DEF2(I,J) + WDZ(I,J  )
        DEFH2(I,J ) = DEF2(I,J) + DEFH2(I,J)
     END DO
  END DO

  DO J=N, 2, -1
     DO I=M, 1, -1
        ! DEFH2(I,J) = 
        VDR(I  ,J-1) = 0.25*2*VDR9(I  ,J-1)*DEFH2(I,J) + VDR(I  ,J-1)
        VDR(I+1,J-1) = 0.25*2*VDR9(I+1,J-1)*DEFH2(I,J) + VDR(I+1,J-1)
        VDR(I  ,J  ) = 0.25*2*VDR9(I  ,J  )*DEFH2(I,J) + VDR(I  ,J  )
        VDR(I+1,J  ) = 0.25*2*VDR9(I+1,J  )*DEFH2(I,J) + VDR(I+1,J  )
        U1( I  ,J-1) = (0.25/(RS(I)**2))*2*(U19(I+1,J-1)+U19(I,J-1)) &
             *DEFH2(I,J) + U1( I  ,J-1)
        U1( I+1,J-1) = (0.25/(RS(I)**2))*2*(U19(I+1,J-1)+U19(I,J-1)) &
             *DEFH2(I,J) + U1( I+1,J-1)
        U1( I  ,J  ) = (0.25/(RS(I)**2))*2*(U19(I+1,J  )+U19(I,J  )) &
             *DEFH2(I,J) + U1( I  ,J  )
        U1( I+1,J  ) = (0.25/(RS(I)**2))*2*(U19(I+1,J  )+U19(I,J  )) &
             *DEFH2(I,J) + U1( I+1,J  )
        UDR(I  ,J-1) = 2*UDR9(I,J-1)*DEFH2(I,J) + UDR(I  ,J-1)
        UDR(I  ,J  ) = 2*UDR9(I,J  )*DEFH2(I,J) + UDR(I  ,J  )
     END DO
  END DO

  DO J=N, 2, -1
     DO I=M, 1, -1
        ! VDZ(I,J) = 
        V1(I,J-1) = -RDZ * VDZ(I,J) + V1(I,J-1)
        V1(I,J  ) =  RDZ * VDZ(I,J) + V1(I,J  )
     END DO
  END DO

  DO J=N, 2, -1
     UZWR(M,J) = UZWR(MP1,J) + UZWR(M,J)
     UZWR(1,J) = 0
     DO I=M, 2, -1
        ! UZWR(I,J) = 
        W1(I-1,J) = -RDR*UZWR(I,J) + W1(I-1,J)
        W1(I  ,J) =  RDR*UZWR(I,J) + W1(I  ,J)
        U1(I,J-1) = -RDZ*UZWR(I,J) + U1(I,J-1)
        U1(I,J  ) =  RDZ*UZWR(I,J) + U1(I,J  )
     END DO
  END DO
  
  DO J=N, 1, -1
     VDR(M,J) = VDR(MP1,J) + VDR(M,J)
     VDR(1,J) = 0
     DO I=M, 2, -1
        ! VDR(I,J) =
        V1(I-1,J) = -R(I)*RDR/RS(I-1) * VDR(I,J) + V1(I-1,J)
        V1(I  ,J) =  R(I)*RDR/RS(I  ) * VDR(I,J) + V1(I  ,J)
     END DO
  END DO

  DO J=N, 1, -1
     DO I=M, 1, -1
        ! WDZ(I,J) = 
        W1(I,J  ) = -RDZ*WDZ(I,J) + W1(I,J  )
        W1(I,J+1) =  RDZ*WDZ(I,J) + W1(I,J+1)
        ! UDR(I,J) =
        U1(I  ,J) = -RDR*UDR(I,J) + U1(I  ,J)
        U1(I+1,J) =  RDR*UDR(I,J) + U1(I+1,J)
     END DO
  END DO

  DO J=N, 1, -1
     DO I=M, 1, -1
        ! TE(I,J) = 
        P1(M,1) =  T19(I,J) &
             *(  XLDCP*QV19(I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1)))**2) &
             *T19(I,J)*TE(I,J) + P1(M,1)
        P1(I,J) = -T19(I,J) &
             *(  XLDCP*QV19(I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1)))**2) &
             *T19(I,J)*TE(I,J) + P1(I,J)
        T1(I,J) = -T19(I,J) &
             *(  XLDCP*QV19(I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1)))**2) &
             *(PN(J)+P19(I,J)-P19(M,1))*TE(I,J) + T1(I,J)
        QV1(I,J)=  T19(I,J) &
             *(  XLDCP          /(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1)))) &
             *TE(I,J) + QV1(I,J)
        T1(I,J) = (1+XLDCP*QV19(I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1)))) &
             *TE(I,J) + T1(I,J)
     END DO
  END DO

  RETURN
END SUBROUTINE DIFFUSE_ADJ

! ==================================================================

SUBROUTINE FWD_NO_VADV
  USE hurr_vars
  IMPLICIT NONE

  ! CONSIDER ADD A SECTION WHERE THE VALUES OF ALL
  !     STATE VARIABLES ARE RESTORED.
  
  DO I = 1 , M  
     PNSURF9 = (P19(I,1)-P19(M,1))/(1.0-0.5*G*DZ/(CP*SST(I)) )
     PDST9   = 1000.0 * (PNS+PNSURF9)**(1.0/XKAPPA)
     QSURF9(I) = 0.622*ESS(I) /(PDST9-ESS(I))
     TSURF9(I) = SST(I)/(PNS+PNSURF9)
  END DO
  
  DO I = 1, M
     CERS9(I)=CE + CD1*SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2)
     CDRS9(I)=CERS9(I)*CD/CE
  END DO

  ! IN SUBROUTINE DIFFUSE, VARIABLES WHOSE VALUES WILL BE CHANGED INCLUDE:
  ! UA, WA, VA, TA, QVA, QLA
  CALL DIFFUSE_FWD
  !
  !     RADIATION (R < 1 DEG/DAY)
  !
  DO J = 1, N
     DO I = 1, M
        TDIF9=T19(I,J)-TB(J)
        QRAD9(I,J) = -DTL*TDIF9*RADT
        IF(QRAD9(I,J) > RADMAX) THEN
           QRAD9(I,J) = RADMAX
        END IF
        IF(QRAD9(I,J) <-RADMAX) THEN
           QRAD9(I,J) = -RADMAX
        END IF
     END DO
  END DO
  !
  !        FORCING FOR U EQUATION 
  !
  DO J = 1, N
     DO I = 2, M
        UA9(I,J) = UA9(I,J) - DTSR(I)*( &
             ( R(I+1)*U9(I+1,J)+R(I  )*U9(I  ,J))*(U9(I+1,J)-U9(I  ,J)) &
             +(R(I  )*U9(I  ,J)+R(I-1)*U9(I-1,J))*(U9(I  ,J)-U9(I-1,J)) ) &
             - DTSZ(I,J) *( &
              RHOW(J+1)*(RS(I)*W9(I,J+1)+RS(I-1)*W9(I-1,J+1))*(U9(I,J+1)-U9(I,J)) &
             +RHOW(J  )*(RS(I)*W9(I,J  )+RS(I-1)*W9(I-1,J  ))*(U9(I,J)-U9(I,J-1)))&
             +DTSV(I) * V9(I,J)**2 + DTSV(I-1) * V9(I-1,J)**2 &
             +DTSF    *(V9(I,J) + V9(I-1,J))
     END DO
  END DO
  !
  !       OUTER BOUNDARY  
  !
  DO J = 1, N
     UA9(MP1,J) = DTL * (V9(M,J)**2/RS(M) + F*V9(M,J))
     IF(U9(MP1,J) + CSTAR > 0) THEN
        UA9(MP1,J)= UA9(MP1,J)-(U9(MP1,J)+CSTAR)*DTL*RDR*(U19(MP1,J)-U19(M,J))
     END IF
  END DO
  !
  !       DRAG AT J = 1   
  !
  UA9(MP1,1) = UA9(MP1,1) &
       - CD*DTL*RDZ*U19(MP1,1)*SQRT(U19(MP1,1)**2 + V19(M,1)**2)
  !
  !       FORCING FOR W EQUATION  
  !
  DO J = 2, N
     DO I = 1, MM1
        WA9(I,J) = WA9(I,J) -DTSRW(J)*( &
              (RHOT(J  )*U9(I+1,J)+RHOT(J-1)*U9(I+1,J-1))*(W9(I+1,J)-W9(I  ,J)) &
             +(RHOT(J  )*U9(I  ,J)+RHOT(J-1)*U9(I  ,J-1))*(W9(I  ,J)-W9(I-1,J)) ) &
             -DTSZW(J)*( &
              (RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  ))*(W9(I,J+1)-W9(I,J  )) &
             +(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  ))*(W9(I,J  )-W9(I,J-1)) ) &
             + DTS*RTBW( J) * ( T9(I,J) - TB(J)  + T9(I,J-1)  - TB(J-1) ) &
             + DTS*RQVBW(J) * (QV9(I,J) - QVB(J) + QV9(I,J-1) - QVB(J-1)) &
             - DTSG * (QL9(I,J) +QL9(I,J-1))
     END DO
  END DO
  !       
  !       OUTER BOUNDARY  
  !       
  DO J = 2 , N
     UBW9(J) = RHOT(J)*(U9(MP1,J)+U9(M,J)) + RHOT(J-1)*(U9(M,J-1)+U9(M,J-1))
     IF(UBW9(J) < 0) THEN
        UBW9(J) = 0
     END IF
  END DO

  I = M
  DO J = 2, N
     WA9(I,J) = WA9(I,J) &
          -DTSRW(J) * UBW9(J) * (W19(I  ,J) - W19(I-1,J)) &
          -DTSZW(J)*( &
          ( RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  ))*(W9(I,J+1)-W9(I,J  )) &
          +(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  ))*(W9(I,J  )-W9(I,J-1)) ) &
          + DTS*RTBW( J) * ( T9(I,J) - TB(J)  + T9(I,J-1)  - TB(J-1) ) &
          + DTS*RQVBW(J) * (QV9(I,J) - QVB(J) + QV9(I,J-1) - QVB(J-1)) &
          - DTSG * (QL9(I,J) +QL9(I,J-1))
  END DO
  !
  !       FORCING FOR  V, T, QV, QL  EQUATIONS
  !
  DO J = 1, N
     DO I = 1, MM1
        VA9(I,J) = VA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (V9(I+1,J  ) - V9(I  ,J  ))  &
             +          R(I  ) * U9(I  ,J) * (V9(I  ,J  ) - V9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (V9(I  ,J+1) - V9(I  ,J  ))  &
             +          RHOW(J  )*W9(I,J  )* (V9(I  ,J  ) - V9(I  ,J-1)) ) &
             -DTL2*(F+ V9(I,J)/RS(I))*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)
        
        TA9(I,J) = TA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (T9(I+1,J  ) - T9(I  ,J  )) &
             +          R(I  ) * U9(I  ,J) * (T9(I  ,J  ) - T9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (T9(I  ,J+1) - T9(I  ,J  )) &
             +          RHOW(J  )*W9(I,J  )* (T9(I  ,J  ) - T9(I  ,J-1)) )! &
             !+QRAD9(I,J)
        
        QVA9(I,J) = QVA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (QV9(I+1,J  ) - QV9(I  ,J  )) &
             +          R(I  ) * U9(I  ,J) * (QV9(I  ,J  ) - QV9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QV9(I  ,J+1) - QV9(I  ,J  )) &
             +          RHOW(J  )*W9(I,J  )* (QV9(I  ,J  ) - QV9(I  ,J-1)) )
        
        QLA9(I,J) = QLA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (QL9(I+1,J  ) - QL9(I  ,J  ))  &
             +          R(I  ) * U9(I  ,J) * (QL9(I  ,J  ) - QL9(I-1,J  )) ) & 
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QL9(I  ,J+1) - QL9(I  ,J  ))  &
             +          RHOW(J  )*W9(I,J  )* (QL9(I  ,J  ) - QL9(I  ,J-1)) )
     END DO
  END DO
  !
  !       OUTER BOUNDARY
  !
  DO J = 1, N
     UB9(J) = R(MP1)*U9(MP1,J) + R(M)*U9(M,J)
     IF(UB9(J) < 0) THEN
        UB9(J) = 0
     END IF
  END DO
  I = M     
  DO J = 1, N
     VA9(I,J) = VA9(I,J) &
          -DTLR(I)* UB9(J) * (V19(I  ,J  ) - V19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* ( V9(I  ,J+1) - V9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )* ( V9(I  ,J  ) - V9(I  ,J-1)) ) &
          -DTL2*(F+ V9(I,J)/RS(I))*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)
     
     TA9(I,J) = TA9(I,J) &
          -DTLR(I)* UB9(J) * (T19(I  ,J  ) - T19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)*(T9(I  ,J+1) - T9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )*(T9(I  ,J  ) - T9(I  ,J-1)) )! &
          !+QRAD9(I,J)
     
     QVA9(I,J) = QVA9(I,J) &
          -DTLR(I)* UB9(J) * (QV19(I  ,J  ) - QV19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QV9(I  ,J+1) - QV9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )* (QV9(I  ,J  ) - QV9(I  ,J-1)) )
     
     QLA9(I,J) = QLA9(I,J) &
          -DTLR(I)* UB9(J) * (QL19(I  ,J  ) - QL19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QL9(I  ,J+1) - QL9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )* (QL9(I  ,J  ) - QL9(I  ,J-1)) )
  END DO
  !
  !         RAIN FALL
  !
  DO J = 1, N
     DO I = 1, M
        DUM(I,J)=0
        IF((QL9(I,J) - 0.001) > 0) DUM(I,J)=VTERM
     END DO
  END DO

  DO J = 2, NM1
     DO I = 1, M
        QLA9(I,J) = QLA9(I,J) + DTLZ(J) * &
             (DUM(I,J+1)*RHOT(J+1)*QL9(I,J+1)-DUM(I,J-1)*RHOT(J-1)*QL9(I,J-1))
     END DO
  END DO
  
  DO I = 1 , M
     DUM(I,2)=0.
     IF((QL19(I,2) - 0.001) > 0) DUM(I,2)=VTERM
     DUM(I,1)=0
     IF((QL19(I,1) - 0.001) > 0) DUM(I,1)=VTERM
     QLA9(I,1) = QLA9(I,1) + 2*DTLZ(1)* &
          (DUM(I,2)*RHOT(2)*QL19(I,2) - DUM(I,1)*RHOT(1)*QL19(I,1 ))
  END DO
  !
  !        TIME SMOOTHER
  !
  DO J = 1, N
     DO I = 1, M
        U9(I+1,J)= U9(I+1,J)+ EPS * (U19(I+1,J)-2.*U9(I+1,J))
        V9(I,J)  = V9(I,J)  + EPS * (V19(I,J)  -2.*V9(I,J)  )
        W9(I,J+1)= W9(I,J+1)+ EPS * (W19(I,J+1)-2.*W9(I,J+1))
        P9(I,J)  = P9(I,J)  + EPS * (P19(I,J)  -2.*P9(I,J)  )
        T9(I,J)  = T9(I,J)  + EPS * (T19(I,J)  -2.*T9(I,J)  )
        QV9(I,J) = QV9(I,J) + EPS * (QV19(I,J) -2.*QV9(I,J) )
        QL9(I,J) = QL9(I,J) + EPS * (QL19(I,J) -2.*QL9(I,J) )
     END DO
  END DO
  !
  !         SMALL TIME STEP
  !
  
  IF(.NOT.ALLOCATED(PS9NS)) ALLOCATE(PS9NS(M,N,NS))
  IF(.NOT.ALLOCATED(P19NS)) ALLOCATE(P19NS(M,N,NS))
  P19NS(:,:,1) = P19
  DO NSMALL = 1, NS
     DO J = 1, N
        DO I = 2, M
           CFACU19=DTS*RDR*CP*0.5&
                *(T19(I,J)*(1.0+0.61*QV9(I,J))+T19(I-1,J)*(1.0+0.61*QV9(I-1,J)))
           U19(I,J)=U19(I,J)-CFACU19*(P19(I,J)-P19(I-1,J))+UA9(I,J)
        END DO
     END DO
     
     DO J = 2, N
        DO I = 1, M
           CFACWS9=DTS*RDZ*CP*0.5* &
                (T19(I,J)+T19(I,J-1)) * (1.0+0.305*(QV9(I,J)+QV9(I,J-1)))
           WS9(I,J)=W19(I,J)-0.5*(1-EP)*CFACWS9*(P19(I,J)-P19(I,J-1))+WA9(I,J)
        END DO
     END DO
     
     DO J = 1, N
        DO I = 1, M
           PS9(I,J)=P19(I,J) &
                -RC2(J)*(R(I+1)*U19(I+1,J) - R(I)*U19(I,J))/RS(I) &
                -0.5*(1.0-EP)*ZC2(J) &
                *(RHOTVW(J+1)*W19(I,J+1)-RHOTVW(J)*W19(I,J))
           PS9NS(I,J,NSMALL) = PS9(I,J)
        END DO
     END DO
     
     DO J = 2, N
        DO I = 1, M
           CFACD9=DTS*RDZ*CP*0.5* &
                (T19(I,J)+T19(I,J-1))*(1.0+0.305*(QV9(I,J)+QV9(I,J-1)))
           D9(I,J) = (WS9(I,J) - CFACD9*0.5*(1+EP)*(PS9(I,J)-PS9(I,J-1)) &
                + C(J) * D9(I,J-1)) * E(J) /A(J)
        END DO
     END DO
     
     DO J = N, 1, -1
        DO I = 1, M
           W19(I,J) = E(J) * W19(I,J+1) + D9(I,J)

           P19(I,J) = PS9(I,J) -0.5*(1.0+EP)*ZC2(J)* &
                (RHOTVW(J+1)*W19(I,J+1) - RHOTVW(J)*W19(I,J))
           IF(NSMALL<=NS-1) P19NS(I,J,NSMALL+1) = P19(I,J)
        END DO
     END DO
  END DO ! DO NSMALL = 1 , NS
  !
  !        OUTER BOUNDARY
  !
  DO J = 1 , N
     U19(MP1,J) = U19(MP1,J) + UA9(MP1,J)
  END DO
  !
  !        ADVANCE V , T , QV , QL
  !

!  END DO TIMELOOP

  RETURN
END SUBROUTINE FWD_NO_VADV

!============================================================================

SUBROUTINE FWD_NO_SMOOTHER
  USE hurr_vars
  IMPLICIT NONE
  
  ! CONSIDER ADD A SECTION WHERE THE VALUES OF ALL
  !     STATE VARIABLES ARE RESTORED.
  
  DO I = 1 , M  
     PNSURF9 = (P19(I,1)-P19(M,1))/(1.0-0.5*G*DZ/(CP*SST(I)) )
     PDST9   = 1000.0 * (PNS+PNSURF9)**(1.0/XKAPPA)
     QSURF9(I) = 0.622*ESS(I) /(PDST9-ESS(I))
     TSURF9(I) = SST(I)/(PNS+PNSURF9)
  END DO
  
  DO I = 1, M
     CERS9(I)=CE + CD1*SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2)
     CDRS9(I)=CERS9(I)*CD/CE
  END DO
  
  ! IN SUBROUTINE DIFFUSE, VARIABLES WHOSE VALUES WILL BE CHANGED INCLUDE:
  ! UA, WA, VA, TA, QVA, QLA
  CALL DIFFUSE_FWD
  !
  !     RADIATION (R < 1 DEG/DAY)
  !
  DO J = 1, N
     DO I = 1, M
        TDIF9=T19(I,J)-TB(J)
        QRAD9(I,J) = -DTL*TDIF9*RADT
        IF(QRAD9(I,J) > RADMAX) THEN
           QRAD9(I,J) = RADMAX
        END IF
        IF(QRAD9(I,J) <-RADMAX) THEN
           QRAD9(I,J) = -RADMAX
        END IF
     END DO
  END DO
  !
  !        FORCING FOR U EQUATION 
  !
  DO J = 1, N
     DO I = 2, M
        UA9(I,J) = UA9(I,J) - DTSR(I)*( &
             ( R(I+1)*U9(I+1,J)+R(I  )*U9(I  ,J))*(U9(I+1,J)-U9(I  ,J)) &
             +(R(I  )*U9(I  ,J)+R(I-1)*U9(I-1,J))*(U9(I  ,J)-U9(I-1,J)) ) &
             - DTSZ(I,J) *( &
              RHOW(J+1)*(RS(I)*W9(I,J+1)+RS(I-1)*W9(I-1,J+1))*(U9(I,J+1)-U9(I,J)) &
             +RHOW(J  )*(RS(I)*W9(I,J  )+RS(I-1)*W9(I-1,J  ))*(U9(I,J)-U9(I,J-1)))&
             +DTSV(I) * V9(I,J)**2 + DTSV(I-1) * V9(I-1,J)**2 &
             +DTSF    *(V9(I,J) + V9(I-1,J))
     END DO
  END DO
  !
  !       OUTER BOUNDARY  
  !
  DO J = 1, N
     UA9(MP1,J) = DTL * (V9(M,J)**2/RS(M) + F*V9(M,J))
     IF(U9(MP1,J) + CSTAR > 0) THEN
        UA9(MP1,J)= UA9(MP1,J)-(U9(MP1,J)+CSTAR)*DTL*RDR*(U19(MP1,J)-U19(M,J))
     END IF
  END DO
  !
  !       DRAG AT J = 1   
  !
  UA9(MP1,1) = UA9(MP1,1) &
       - CD*DTL*RDZ*U19(MP1,1)*SQRT(U19(MP1,1)**2 + V19(M,1)**2)
  !
  !       FORCING FOR W EQUATION  
  !
  DO J = 2, N
     DO I = 1, MM1
        WA9(I,J) = WA9(I,J) -DTSRW(J)*( &
              (RHOT(J  )*U9(I+1,J)+RHOT(J-1)*U9(I+1,J-1))*(W9(I+1,J)-W9(I  ,J)) &
             +(RHOT(J  )*U9(I  ,J)+RHOT(J-1)*U9(I  ,J-1))*(W9(I  ,J)-W9(I-1,J)) ) &
             -DTSZW(J)*( &
              (RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  ))*(W9(I,J+1)-W9(I,J  )) &
             +(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  ))*(W9(I,J  )-W9(I,J-1)) ) &
             + DTS*RTBW( J) * ( T9(I,J) - TB(J)  + T9(I,J-1)  - TB(J-1) ) &
             + DTS*RQVBW(J) * (QV9(I,J) - QVB(J) + QV9(I,J-1) - QVB(J-1)) &
             - DTSG * (QL9(I,J) +QL9(I,J-1))
     END DO
  END DO
  !       
  !       OUTER BOUNDARY  
  !       
  DO J = 2 , N
     UBW9(J) = RHOT(J)*(U9(MP1,J)+U9(M,J)) + RHOT(J-1)*(U9(M,J-1)+U9(M,J-1))
     IF(UBW9(J) < 0) THEN
        UBW9(J) = 0
     END IF
  END DO

  I = M
  DO J = 2, N
     WA9(I,J) = WA9(I,J) &
          -DTSRW(J) * UBW9(J) * (W19(I  ,J) - W19(I-1,J)) &
          -DTSZW(J)*( &
          ( RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  ))*(W9(I,J+1)-W9(I,J  )) &
          +(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  ))*(W9(I,J  )-W9(I,J-1)) ) &
          + DTS*RTBW( J) * ( T9(I,J) - TB(J)  + T9(I,J-1)  - TB(J-1) ) &
          + DTS*RQVBW(J) * (QV9(I,J) - QVB(J) + QV9(I,J-1) - QVB(J-1)) &
          - DTSG * (QL9(I,J) +QL9(I,J-1))
  END DO
  !
  !       FORCING FOR  V, T, QV, QL  EQUATIONS
  !
  DO J = 1, N
     DO I = 1, MM1
        VA9(I,J) = VA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (V9(I+1,J  ) - V9(I  ,J  ))  &
             +          R(I  ) * U9(I  ,J) * (V9(I  ,J  ) - V9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (V9(I  ,J+1) - V9(I  ,J  ))  &
             +          RHOW(J  )*W9(I,J  )* (V9(I  ,J  ) - V9(I  ,J-1)) ) &
             -DTL2*(F+ V9(I,J)/RS(I))*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)
        
        TA9(I,J) = TA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (T9(I+1,J  ) - T9(I  ,J  )) &
             +          R(I  ) * U9(I  ,J) * (T9(I  ,J  ) - T9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (T9(I  ,J+1) - T9(I  ,J  )) &
             +          RHOW(J  )*W9(I,J  )* (T9(I  ,J  ) - T9(I  ,J-1)) )! &
             !+QRAD9(I,J)
        
        QVA9(I,J) = QVA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (QV9(I+1,J  ) - QV9(I  ,J  )) &
             +          R(I  ) * U9(I  ,J) * (QV9(I  ,J  ) - QV9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QV9(I  ,J+1) - QV9(I  ,J  )) &
             +          RHOW(J  )*W9(I,J  )* (QV9(I  ,J  ) - QV9(I  ,J-1)) )
        
        QLA9(I,J) = QLA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (QL9(I+1,J  ) - QL9(I  ,J  ))  &
             +          R(I  ) * U9(I  ,J) * (QL9(I  ,J  ) - QL9(I-1,J  )) ) & 
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QL9(I  ,J+1) - QL9(I  ,J  ))  &
             +          RHOW(J  )*W9(I,J  )* (QL9(I  ,J  ) - QL9(I  ,J-1)) )
     END DO
  END DO
  !
  !       OUTER BOUNDARY
  !
  DO J = 1, N
     UB9(J) = R(MP1)*U9(MP1,J) + R(M)*U9(M,J)
     IF(UB9(J) < 0) THEN
        UB9(J) = 0
     END IF
  END DO
  I = M     
  DO J = 1, N
     VA9(I,J) = VA9(I,J) &
          -DTLR(I)* UB9(J) * (V19(I  ,J  ) - V19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* ( V9(I  ,J+1) - V9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )* ( V9(I  ,J  ) - V9(I  ,J-1)) ) &
          -DTL2*(F+ V9(I,J)/RS(I))*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)
     
     TA9(I,J) = TA9(I,J) &
          -DTLR(I)* UB9(J) * (T19(I  ,J  ) - T19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)*(T9(I  ,J+1) - T9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )*(T9(I  ,J  ) - T9(I  ,J-1)) )! &
          !+QRAD9(I,J)
     
     QVA9(I,J) = QVA9(I,J) &
          -DTLR(I)* UB9(J) * (QV19(I  ,J  ) - QV19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QV9(I  ,J+1) - QV9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )* (QV9(I  ,J  ) - QV9(I  ,J-1)) )
     
     QLA9(I,J) = QLA9(I,J) &
          -DTLR(I)* UB9(J) * (QL19(I  ,J  ) - QL19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QL9(I  ,J+1) - QL9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )* (QL9(I  ,J  ) - QL9(I  ,J-1)) )
  END DO
  !
  !         RAIN FALL
  !
  DO J = 1, N
     DO I = 1, M
        DUM(I,J)=0
        IF((QL9(I,J) - 0.001) > 0) DUM(I,J)=VTERM
     END DO
  END DO
  
  DO J = 2, NM1
     DO I = 1, M
        QLA9(I,J) = QLA9(I,J) + DTLZ(J) * &
             (DUM(I,J+1)*RHOT(J+1)*QL9(I,J+1)-DUM(I,J-1)*RHOT(J-1)*QL9(I,J-1))
     END DO
  END DO
  
  DO I = 1 , M
     DUM(I,2)=0.
     IF((QL19(I,2) - 0.001) > 0) DUM(I,2)=VTERM
     DUM(I,1)=0
     IF((QL19(I,1) - 0.001) > 0) DUM(I,1)=VTERM
     QLA9(I,1) = QLA9(I,1) + 2*DTLZ(1)* &
          (DUM(I,2)*RHOT(2)*QL19(I,2) - DUM(I,1)*RHOT(1)*QL19(I,1 ))
  END DO
  !
  !        TIME SMOOTHER
  !

!  END DO TIMELOOP


  RETURN
END SUBROUTINE FWD_NO_SMOOTHER
