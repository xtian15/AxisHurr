SUBROUTINE hurr_fwd
  use hurr_vars
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
  DO J=1,N
     DO I=1,M
        V19(I,J)  = V19(I,J)  + VA9(I,J)
        T19(I,J)  = T19(I,J)  + TA9(I,J)
        QV19(I,J) = QV19(I,J) + QVA9(I,J)
        QL19(I,J) = QL19(I,J) + QLA9(I,J)
        IF(QV19(I,J) < 0) QV19(I,J)= 0
        IF(QL19(I,J) < 0) QL19(I,J)=0
        !
        !     CONDENSATION / EVAPORATION
        !
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
        IF(QV19(I,J) < 0) QV19(I,J)=0
        IF(QL19(I,J) < 0) QL19(I,J)=0
     END DO
  END DO
  !
  !        TIME FLIP
  !
  DO J=1,N
     DO I=1,M
        UTEMP9    = U9(I+1,J) + EPS * U19(I+1,J)
        VTEMP9    = V9(I,J)   + EPS * V19(I,J)
        WTEMP9    = W9(I,J+1) + EPS * W19(I,J+1)
        PTEMP9    = P9(I,J)   + EPS * P19(I,J)
        TTEMP9    = T9(I,J)   + EPS * T19(I,J)
        QVTEMP9   = QV9(I,J)  + EPS * QV19(I,J)
        QLTEMP9   = QL9(I,J)  + EPS * QL19(I,J)

        U9(I+1,J) = U19(I+1,J)
        V9(I,J)   = V19(I,J)
        W9(I,J+1) = W19(I,J+1)
        P9(I,J)   = P19(I,J)
        T9(I,J)   = T19(I,J)
        QV9(I,J)  = QV19(I,J)
        QL9(I,J)  = QL19(I,J)
        U19(I+1,J)= UTEMP9
        V19(I,J)  = VTEMP9
        W19(I,J+1)= WTEMP9
        T19(I,J)  = TTEMP9
        P19(I,J)  = PTEMP9
        QV19(I,J) = QVTEMP9
        QL19(I,J) = QLTEMP9
     END DO
  END DO
  IF(ITIME == 1) THEN
     DTL = 2. * DT
     NS  = 2  * NS
     DTL2 = .5 * DTL
     DO I = 1, M
        DTLR(I)  = 0.5 * DTL*RDR/RS(I)
     END DO
     DO J = 1 , N
        DTLZ(J)  = 0.5 * DTL*RDZ/RHOT(J)
     END DO
  END IF
  
  !END DO TIMELOOP
  
  RETURN
END SUBROUTINE HURR_FWD

! =================================================
SUBROUTINE DIFFUSE_FWD
  use hurr_vars
  IMPLICIT NONE

  ! ----- LOCAL DEPENDENT VARIABLES -----
  REAL*8 :: &
       TRR9(M  ,N  ), TTT9(M  ,N  ), TZZ9(M  ,NP1), TRZ9(MP1,NP1), &
       TRT9(MP1,N  ), TZT9(M  ,NP1), TR9( MP1,N  ), TZ9( M  ,NP1), &
       TAD9(M  ,N  ), QVR9(MP1,N  ), QVZ9(M  ,NP1), QLR9(MP1,N  ), &
       QLZ9(M  ,NP1), XKM9(M  ,NP1), UDR9(M  ,N  ), WDZ9(M  ,N  ), &
       VDR9(MP1,N  ), VDZ9(M  ,N  ), TE9( M  ,N  ), &
       XKMH9(M ,NP1), UZWR9(MP1, N), DEFH29(M,N  ), DEF29(M,N), DTDZ9(M,N)
  
  REAL*8 :: TT9, AA9, VABS9
  ! -------------------------------------
  
  DO J = 1, N
     DO I = 1, M
        TE9(I,J)= T19(I,J) &
             *(1+XLDCP*QV19(I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1))))
     END DO
  END DO
  
  DO J = 1, N
     DO I = 1, M
        UDR9(I,J) = RDR * (U19(I+1,J) - U19(I,J))
        WDZ9(I,J) = RDZ * (W19(I,J+1) - W19(I,J))
     END DO
  END DO
  
  DO J = 1, N
     DO I = 2, M
        VDR9(I,J) = R(I)*RDR * (V19(I,J)/RS(I) - V19(I-1,J)/RS(I-1))
     END DO
     VDR9(1  ,J) = 0
     VDR9(MP1,J) = VDR9(M,J)
  END DO
  
  DO J = 2, N
     DO I = 2, M
        UZWR9(I,J) = RDZ*(U19(I,J)-U19(I,J-1)) + RDR*(W19(I,J)-W19(I-1,J))
     END DO
     UZWR9(1  ,J) = 0
     UZWR9(MP1,J) = UZWR9(M,J)
  END DO
  
  DO J = 2, N
     DO I = 1, M
        VDZ9(I,J) = RDZ * (V19(I,J) - V19(I,J-1))
     END DO
  END DO
  DO J = 2, N
     DO I = 1, M
        DEFH29(I,J) = UDR9(I,J)**2 + UDR9(I,J-1)**2 &
             + (0.25/(RS(I)**2)) &
             *((U19(I+1,J  )+U19(I,J  ))**2 + (U19(I+1,J-1)+U19(I,J-1))**2) &
             +0.25* &
             (VDR9(I+1,J  )**2+VDR9(I,J  )**2+VDR9(I+1,J-1)**2+VDR9(I,J-1)**2)
     END DO
  END DO
  
  DO J = 2, N
     DO I = 1, M
        DEF29(I,J) = DEFH29(I,J) + WDZ9(I,J)**2 + WDZ9(I,J-1)**2 &
             + 0.5*(UZWR9(I+1,J)**2 + UZWR9(I,J)**2) + VDZ9(I,J)**2
     END DO
  END DO
  
  DO J = 2, N
     DO I = 1, M
        DTDZ9(I,J) = RDZ*( &
             2.0 * RTBW( J) * (T19( I,J) - T19( I,J-1)) &
            +2.0 * RQVBW(J) * (QV19(I,J) - QV19(I,J-1)) )
        IF( (QL19(I,J) * QL19(I,J-1)) >= 1.0E-8 ) THEN
           TT9 = 0.5 &
                *(T19(I,J  ) * ( PN(J  )+P19(I,J  )-P19(M,1)) &
                + T19(I,J-1) * ( PN(J-1)+P19(I,J-1)-P19(M,1)) )
           AA9=(   1 + 4362    *(QV19(I,J)+QV19(I,J-1)) / TT9   ) &
                / (1 + 6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2)) &
                * (2 / (TB(J)+TB(J-1)) )
           DTDZ9(I,J) = RDZ* G * ( AA9 * (TE9(I,J) - TE9(I,J-1)) &
                  - (QL19(I,J)+QV19(I,J)-QL19(I,J-1)-QV19(I,J-1)) )
        END IF
     END DO
  END DO
  
  DO J = 2, N
     DO I = 1, M
        XKM9(I,J) = 0
        IF(DEF29(I,J) > DTDZ9(I,J)) THEN
           XKM9(I,J)= XVL2 * SQRT(DEF29(I,J) - DTDZ9(I,J))
        END IF
        XKMH9(I,J) = XHL2 * SQRT(DEFH29(I,J))
        IF(XKM9(I,J) >= 0.4*DZ**2/DTL) THEN
           XKM9(I,J) = 0.4 *DZ**2/DTL
        END IF
        IF(XKMH9(I,J) < XKM9(I,J)) THEN
           XKMH9(I,J) = XKM9(I,J)
        END IF
     END DO
  END DO
  
  DO I = 1, M
     XKM9( I,1  ) = XKM9( I,2)
     XKM9( I,NP1) = XKM9( I,N)
     XKMH9(I,1  ) = XKMH9(I,2)
     XKMH9(I,NP1) = XKMH9(I,N)
  END DO
  !
  !        CALCULATE STRESS
  !
  !        TRR(M,N)
  !
  DO J = 1 , N
     DO I = 1 , M
        TRR9(I,J) = (XKMH9(I,J+1) + XKMH9(I,J)) * UDR9(I,J)
     END DO
  END DO
  !
  !        TTT(M,N)
  !
  DO J = 1, N
     TTT9(1,J) = 0
     DO I = 2, M
        TTT9(I,J) = U19(I,J)/R(I) * 0.5 &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
  END DO
  !
  !        TRZ(M,NP1)
  !
  DO J=2, N
     DO I=2, M
        TRZ9(I,J)=0.5* (XKM9(I-1,J) + XKM9(I,J)) * UZWR9(I,J)
     END DO
  END DO

  DO I =2, M
     TRZ9(I,1)= 0.5 * (CDRS9(I)+CDRS9(I-1)) * U19(I,1) &
          * SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2)
     TRZ9(I,NP1) = 0
  END DO
  
  DO J = 1, NP1
     TRZ9(1  ,J) = 0
     TRZ9(MP1,J) = TRZ9(M,J) * R(M) / R(MP1)
  END DO
  !
  !        TRT(M,N)
  !
  DO J=1, N
     DO I=2, M
        TRT9(I,J) = 0.25*VDR9(I,J) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
     TRT9(1  ,J) = 0
     TRT9(MP1,J) = TRT9(M,J) * R(M)**2 / R(MP1)**2
  END DO
  !
  !        TZT(M,NP1)
  !
  DO I = 1, M
     DO J = 2, N
        TZT9(I,J) = XKM9(I,J)*VDZ9(I,J)
     END DO
     TZT9(I,1  )= CDRS9(I) * V19(I,1) &
          *SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2)
     TZT9(I,NP1)= 0
  END DO
  !
  !          TZZ(M,N)
  !
  DO J = 1, N
     DO I = 1, M
        TZZ9(I,J) = (XKM9(I,J+1) + XKM9(I,J)) * WDZ9(I,J)
     END DO
  END DO
  !
  !       TEMPERATURE FLUX
  !
  !          TR(MP1,N)
  !
  DO J=1,N
     DO I=2,M
        TR9(I,J) = 0.25 * RDR*(T19(I,J)-T19(I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
     TR9(1  ,J)=0
     TR9(MP1,J)= TR9(M,J) * R(M) / R(MP1)
  END DO
  !
  !         TZ(M,NP1)
  !
  DO I=1,M
     DO J=2,N
        TZ9(I,J) = XKM9(I,J) * RDZ * (T19(I,J) - T19(I,J-1))
        TZ9(I,J) = TZ9(I,J) * PD(J)
     END DO
     VABS9=SQRT(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2)
     TZ9(I,1  ) = (T19(I,1)-TSURF9(I))*CERS9(I)*VABS9
     TZ9(I,1  ) =  TZ9(I,1)*PD(1)
     TZ9(I,NP1) = 0
  END DO
  !
  !       QV FLUX
  !
  !       QVR(MP1,N)
  !
  DO J=1,N
     DO I=2,M
        QVR9(I,J)= 0.25 *RDR*(QV19(I,J)-QV19(I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
     QVR9(1  ,J) = 0
     QVR9(MP1,J) = QVR9(M,J) * R(M) / R(MP1)
  END DO
  !
  !         QVZ(M,NP1)
  !
  DO I=1,M
     DO J=2,N
        QVZ9(I,J) = XKM9(I,J) * RDZ * (QV19(I,J) - QV19(I,J-1))
        QVZ9(I,J) = QVZ9(I,J) * PD(J)/PN(J)
     END DO
     VABS9=SQRT(0.25*(U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2)
     QVZ9(I,1  ) = (QV19(I,1)-QSURF9(I))*CERS9(I)*VABS9
     QVZ9(I,1  ) =  QVZ9(I,1)*PD(1)/PN(1)
     QVZ9(I,NP1) =  0
  END DO
  !
  !       QL FLUX
  !
  !       QLR(MP1,N)
  !
  DO J=1, N
     DO I=2, M
        QLR9(I,J)= 0.25 *RDR*(QL19(I,J)-QL19(I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
     QLR9(1  ,J) = 0
     QLR9(MP1,J) = QLR9(M,J) * R(M) / R(MP1)
  END DO
  !
  !         QLZ(M,NP1)
  !
  DO I=1, M
     DO J=2, N
        QLZ9(I,J) = XKM9(I,J) * RDZ * (QL19(I,J) - QL19(I,J-1))
     END DO
     QLZ9(I,1  ) = 0
     QLZ9(I,NP1) = 0
  END DO
  DO J = 1, N
     DO I = 2, M
        UA9(I,J)= -DTS*( &
             - RDR * (RS(I)*TRR9(I,J) - RS(I-1)*TRR9(I-1,J)) / R(I) &
             + TTT9(I,J) / R(I) &
             - RDZ * (TRZ9(I,J+1) - TRZ9(I,J)) ) &
               + DTS * TAU(J) * U19(I,J)
     END DO
  END DO
  DO J = 2, N
     DO I = 1, M
        WA9(I,J) = -DTS*( &
             - RDR * (R(I+1)*TRZ9(I+1,J) - R(I)*TRZ9(I,J)) / RS(I) &
             - RDZ * (TZZ9(I,J) - TZZ9(I,J-1)) ) &
             + DTS * 0.5 * (TAU(J) + TAU(J-1)) * W19(I,J)
     END DO
  END DO
  
  DO I = 1 , M
     DO J = 1 , N-1
        VA9(I,J) = -DTL * ( &
             - RDR * (R(I+1)**2*TRT9(I+1,J)-R(I)**2*TRT9(I,J)) / (RS(I)**2) &
             - RDZ * (TZT9(I,J+1) - TZT9(I,J)) ) &
             + DTL * TAU(J) * V19(I,J)
        
        TA9(I,J) = -DTL * ( &
             - RDR * (R(I+1)*TR9(I+1,J) - R(I)*TR9(I,J)) / RS(I) &
             - RDZ * (TZ9(I,J+1) - TZ9(I,J))/(0.5*(PD(J)+PD(J+1))) ) &
             + DTL * TAU(J) * (T19(I,J) - TB(J))
        
        QVA9(I,J) = -DTL * ( &
             - RDR * (R(I+1)*QVR9(I+1,J) - R(I)*QVR9(I,J)) / RS(I) &
             - RDZ * (QVZ9(I,J+1) - QVZ9(I,J) )/ &
             (0.5*(PD(J)/PN(J)+PD(J+1)/PN(J+1))) ) &
             + DTL * TAU(J) * QV19(I,J)
        
        QLA9(I,J) = -DTL * ( &
             - RDR * (R(I+1)*QLR9(I+1,J) - R(I)*QLR9(I,J)) / RS(I) &
             - RDZ * (QLZ9(I,J+1) - QLZ9(I,J)) ) &
             + DTL * TAU(J) * QL19(I,J)
     END DO
     VA9(I,J) = -DTL * ( &
          - RDR * (R(I+1)**2*TRT9(I+1,J)-R(I)**2*TRT9(I,J)) / (RS(I)**2) &
          - RDZ * (TZT9(I,J+1) - TZT9(I,J)) ) &
          + DTL * TAU(J) * V19(I,J)
     
     TA9(I,J) = -DTL * ( &
          - RDR * (R(I+1)*TR9(I+1,J) - R(I)*TR9(I,J)) / RS(I) &
          - RDZ * (TZ9(I,J+1) - TZ9(I,J) )/PD(J) ) &
          + DTL * TAU(J)*(T19(I,J) - TB(J))
     
     QVA9(I,J) = -DTL * ( &
          - RDR * (R(I+1)*QVR9(I+1,J) - R(I)*QVR9(I,J) ) / RS(I) &
          - RDZ * (QVZ9(I,J+1) - QVZ9(I,J)) / (PD(J)/PN(J)) ) &
          + DTL * TAU(J) * QV19(I,J)
       
     QLA9(I,J) = -DTL * ( &
          - RDR * (R(I+1)*QLR9(I+1,J) - R(I)*QLR9(I,J)) / RS(I) &
          - RDZ * (QLZ9(I,J+1) - QLZ9(I,J)) ) &
          + DTL * TAU(J) * QL19(I,J)
  END DO
  !
  !       Dissipative Heating
  !
  DO I=1,M-1
     DO J=2,N-1
        TAD9(I,J) = DTL / (PN(J)*1005.0) &
             *(TZT9(I,J+1)*VDZ9(I,J+1) &
             + XKM9(I,J)*0.25*(UZWR9(I,J) + UZWR9(I+1,J))**2 &
             + XKMH9(I,J+1)*(1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))**2 &
             + XKMH9(I,J+1)*0.25*(UDR9(I,J)+UDR9(I,J+1))**2 )
        TA9(I,J) = TA9(I,J)+TAD9(I,J)
     END DO
     TAD9(I,1)=DTL /(PN(1)*1005.0) &
          *( RDZ*CDRS9(I)*(0.25*(U19(I,1)+U19(I+1,1))**2 + V19(I,1)**2)**1.5 &
          + XKMH9(I,2)*(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))**2 &
          + XKMH9(I,2)*0.25 * (UDR9(I,1)+UDR9(I,2))**2 )
     TA9(I,1) = TA9(I,1)+TAD9(I,1)
  END DO
  
  RETURN
END SUBROUTINE DIFFUSE_FWD
