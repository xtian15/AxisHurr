
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