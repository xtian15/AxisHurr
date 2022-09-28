SUBROUTINE hurr_tlm
  USE hurr_vars
  IMPLICIT NONE
  
  DO I = 1 , M
     PNSURF  = (P1( I,1)-P1( M,1))/(1.0-0.5*G*DZ/(CP*SST(I)) )
     PNSURF9 = (P19(I,1)-P19(M,1))/(1.0-0.5*G*DZ/(CP*SST(I)) )
     PDST    = 1000.0 * (1.0/XKAPPA)*(PNS+PNSURF9)**(1.0/XKAPPA-1.0)*PNSURF
     PDST9   = 1000.0 * (PNS+PNSURF9)**(1.0/XKAPPA)
     QSURF(I)  = 0.622*ESS(I) /(PDST9-ESS(I))**2 * (-PDST)
     QSURF9(I) = 0.622*ESS(I) /(PDST9-ESS(I))
     TSURF(I)  = SST(I)/(PNS+PNSURF9)**2 * (-PNSURF)
     TSURF9(I) = SST(I)/(PNS+PNSURF9)
  END DO

  DO I = 1, M
     IF(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2>0) THEN
        CERS(I) = 0.5*CD1/SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             * (0.25*2.0*(U19(I,1)+U19(I+1,1))*(U1(I,1)+U1(I+1,1)) &
             +       2.0*V19(I,1)*V1(I,1))
     ELSE
        CERS(I) = 0
     END IF
     CERS9(I)=CE + CD1*SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2)
     
     CDRS(I) =CERS(I) *CD/CE
     CDRS9(I)=CERS9(I)*CD/CE
  END DO
  
  ! IN SUBROUTINE DIFFUSE, VARIABLES WHOSE VALUES WILL BE CHANGED INCLUDE:
  ! UA, WA, VA, TA, QVA, QLA
  CALL DIFFUSE_TLM
  !
  !     RADIATION (R < 1 DEG/DAY)
  !
  DO J = 1, N
     DO I = 1, M
        TDIF =T1( I,J)
        TDIF9=T19(I,J)-TB(J)
        QRAD(I,J)  = -DTL*TDIF *RADT
        QRAD9(I,J) = -DTL*TDIF9*RADT
        IF(QRAD9(I,J) > RADMAX) THEN
           QRAD( I,J) = 0
           QRAD9(I,J) = RADMAX
        END IF
        IF(QRAD9(I,J) <-RADMAX) THEN
           QRAD( I,J) =  0
           QRAD9(I,J) = -RADMAX
        END IF
     END DO
  END DO
  !
  !        FORCING FOR U EQUATION
  !
  DO J = 1, N
     DO I = 2, M
        UA( I,J) = UA( I,J) - DTSR(I)*( &
             ( R(I+1)*U( I+1,J)+R(I  )*U( I  ,J))*(U9(I+1,J)-U9(I  ,J)) &
             +(R(I+1)*U9(I+1,J)+R(I  )*U9(I  ,J))*(U( I+1,J)-U( I  ,J)) &
             +(R(I  )*U( I  ,J)+R(I-1)*U( I-1,J))*(U9(I  ,J)-U9(I-1,J)) &
             +(R(I  )*U9(I  ,J)+R(I-1)*U9(I-1,J))*(U( I  ,J)-U( I-1,J)) ) &
             - DTSZ(I,J) *( &
             RHOW(J+1)*(RS(I)*W( I,J+1)+RS(I-1)*W( I-1,J+1))*(U9(I,J+1)-U9(I,J)) &
            +RHOW(J+1)*(RS(I)*W9(I,J+1)+RS(I-1)*W9(I-1,J+1))*(U( I,J+1)-U( I,J)) &
            +RHOW(J  )*(RS(I)*W( I,J  )+RS(I-1)*W( I-1,J  ))*(U9(I,J)-U9(I,J-1)) &
            +RHOW(J  )*(RS(I)*W9(I,J  )+RS(I-1)*W9(I-1,J  ))*(U( I,J)-U( I,J-1)))&
            +DTSV(I) * 2*V9(I,J)*V(I,J) + DTSV(I-1) * 2*V9(I-1,J)*V(I-1,J) &
            +DTSF    * (V(I,J) + V(I-1,J))
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
     UA( MP1,J) = DTL * (2*V9(M,J)*V(M,J)/RS(M) + F*V(M,J))
     UA9(MP1,J) = DTL * (V9(M,J)**2/RS(M) + F*V9(M,J))
     IF(U9(MP1,J) + CSTAR > 0) THEN
        UA( MP1,J)= UA( MP1,J)-(U( MP1,J)      )*DTL*RDR*(U19(MP1,J)-U19(M,J))&
             -(U9(MP1,J)+CSTAR)*DTL*RDR*(U1( MP1,J)-U1( M,J))
        UA9(MP1,J)= UA9(MP1,J)-(U9(MP1,J)+CSTAR)*DTL*RDR*(U19(MP1,J)-U19(M,J))
     END IF
  END DO
  !
  !       DRAG AT J = 1
  !
  IF(U19(MP1,1)**2 + V19(M,1)**2>0) THEN
     UA( MP1,1) = UA( MP1,1) &
          - CD*DTL*RDZ*U1( MP1,1)*SQRT(U19(MP1,1)**2 + V19(M,1)**2) &
          - CD*DTL*RDZ*U19(MP1,1)*0.5/SQRT(U19(MP1,1)**2 + V19(M,1)**2) &
          *(2*U19(MP1,1)*U1(MP1,1) + 2*V19(M,1)*V1(M,1))
  END IF
  UA9(MP1,1) = UA9(MP1,1) &
       - CD*DTL*RDZ*U19(MP1,1)*SQRT(U19(MP1,1)**2 + V19(M,1)**2)
  !
  !       FORCING FOR W EQUATION
  !
  DO J = 2, N
     DO I = 1, MM1
        WA( I,J) = WA( I,J) -DTSRW(J)*( &
             (RHOT(J  )*U( I+1,J)+RHOT(J-1)*U( I+1,J-1))*(W9(I+1,J)-W9(I  ,J)) &
            +(RHOT(J  )*U9(I+1,J)+RHOT(J-1)*U9(I+1,J-1))*(W( I+1,J)-W( I  ,J)) &
            +(RHOT(J  )*U( I  ,J)+RHOT(J-1)*U( I  ,J-1))*(W9(I  ,J)-W9(I-1,J)) &
            +(RHOT(J  )*U9(I  ,J)+RHOT(J-1)*U9(I  ,J-1))*(W( I  ,J)-W( I-1,J)) ) &
             -DTSZW(J)*( &
             (RHOW(J+1)*W( I,J+1)+RHOW(J  )*W( I  ,J  ))*(W9(I,J+1)-W9(I,J  )) &
            +(RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  ))*(W( I,J+1)-W( I,J  )) &
            +(RHOW(J-1)*W( I,J-1)+RHOW(J  )*W( I  ,J  ))*(W9(I,J  )-W9(I,J-1)) &
            +(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  ))*(W( I,J  )-W( I,J-1)) ) &
             + DTS*RTBW( J) * ( T( I,J)          + T( I,J-1)            ) &
             + DTS*RQVBW(J) * (QV( I,J)          + QV( I,J-1)           ) &
             - DTSG * (QL( I,J) +QL( I,J-1))
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

  DO J = 2, N
     UB( J) = RHOT(J)*(U( MP1,J)+U( M,J)) + RHOT(J-1)*(U( M,J-1)+U( M,J-1))
     UB9(J) = RHOT(J)*(U9(MP1,J)+U9(M,J)) + RHOT(J-1)*(U9(M,J-1)+U9(M,J-1))
     IF(UB9(J) < 0) THEN
        UB( J) = 0
        UB9(J) = 0
     END IF
  END DO

  I = M
  DO J = 2, N
     WA( I,J) = WA( I,J) &
          -DTSRW(J) * UB( J) * (W19(I  ,J) - W19(I-1,J)) &
          -DTSRW(J) * UB9(J) * (W1( I  ,J) - W1( I-1,J)) &
          -DTSZW(J)*( &
          ( RHOW(J+1)*W( I,J+1)+RHOW(J  )*W( I  ,J  ))*(W9(I,J+1)-W9(I,J  )) &
          +(RHOW(J+1)*W9(I,J+1)+RHOW(J  )*W9(I  ,J  ))*(W( I,J+1)-W( I,J  )) &
          +(RHOW(J-1)*W( I,J-1)+RHOW(J  )*W( I  ,J  ))*(W9(I,J  )-W9(I,J-1)) &
          +(RHOW(J-1)*W9(I,J-1)+RHOW(J  )*W9(I  ,J  ))*(W( I,J  )-W( I,J-1)) ) &
          + DTS*RTBW( J) * ( T( I,J)          + T( I,J-1)            ) &
          + DTS*RQVBW(J) * (QV( I,J)          + QV( I,J-1)           ) &
          - DTSG * (QL( I,J) +QL( I,J-1))
     WA9(I,J) = WA9(I,J) &
          -DTSRW(J) * UB9(J) * (W19(I  ,J) - W19(I-1,J)) &
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
        VA( I,J) = VA( I,J) &
             -DTLR(I)* (R(I+1) * U( I+1,J) * (V9(I+1,J  ) - V9(I  ,J  ))  &
             +          R(I+1) * U9(I+1,J) * (V( I+1,J  ) - V( I  ,J  ))  &
             +          R(I  ) * U( I  ,J) * (V9(I  ,J  ) - V9(I-1,J  )) &
             +          R(I  ) * U9(I  ,J) * (V( I  ,J  ) - V( I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W( I,J+1)* (V9(I  ,J+1) - V9(I  ,J  ))  &
             +          RHOW(J+1)*W9(I,J+1)* (V( I  ,J+1) - V( I  ,J  ))  &
             +          RHOW(J  )*W( I,J  )* (V9(I  ,J  ) - V9(I  ,J-1))  &
             +          RHOW(J  )*W9(I,J  )* (V( I  ,J  ) - V( I  ,J-1)) ) &
             -DTL2*(   V( I,J)/RS(I))*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)&
             -DTL2*(F+ V9(I,J)/RS(I))*(R(I+1)*U( I+1,J)+R(I)*U( I,J))/RS(I)
        VA9(I,J) = VA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (V9(I+1,J  ) - V9(I  ,J  ))  &
             +          R(I  ) * U9(I  ,J) * (V9(I  ,J  ) - V9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (V9(I  ,J+1) - V9(I  ,J  ))  &
             +          RHOW(J  )*W9(I,J  )* (V9(I  ,J  ) - V9(I  ,J-1)) ) &
             -DTL2*(F+ V9(I,J)/RS(I))*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)
        
        TA( I,J) = TA( I,J) &
             -DTLR(I)* (R(I+1) * U( I+1,J) * (T9(I+1,J  ) - T9(I  ,J  )) &
             +          R(I+1) * U9(I+1,J) * (T( I+1,J  ) - T( I  ,J  )) &
             +          R(I  ) * U( I  ,J) * (T9(I  ,J  ) - T9(I-1,J  )) &
             +          R(I  ) * U9(I  ,J) * (T( I  ,J  ) - T( I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W( I,J+1)* (T9(I  ,J+1) - T9(I  ,J  )) &
             +          RHOW(J+1)*W9(I,J+1)* (T( I  ,J+1) - T( I  ,J  )) &
             +          RHOW(J  )*W( I,J  )* (T9(I  ,J  ) - T9(I  ,J-1)) &
             +          RHOW(J  )*W9(I,J  )* (T( I  ,J  ) - T( I  ,J-1)) )! &
             !+QRAD( I,J)
        TA9(I,J) = TA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (T9(I+1,J  ) - T9(I  ,J  )) &
             +          R(I  ) * U9(I  ,J) * (T9(I  ,J  ) - T9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (T9(I  ,J+1) - T9(I  ,J  )) &
             +          RHOW(J  )*W9(I,J  )* (T9(I  ,J  ) - T9(I  ,J-1)) )! &
             !+QRAD9(I,J)
        
        QVA( I,J) = QVA( I,J) &
             -DTLR(I)* (R(I+1) * U( I+1,J) * (QV9(I+1,J  ) - QV9(I  ,J  )) &
             +          R(I+1) * U9(I+1,J) * (QV( I+1,J  ) - QV( I  ,J  )) &
             +          R(I  ) * U( I  ,J) * (QV9(I  ,J  ) - QV9(I-1,J  )) &
             +          R(I  ) * U9(I  ,J) * (QV( I  ,J  ) - QV( I-1,J  )) )&
             -DTLZ(J)* (RHOW(J+1)*W( I,J+1)* (QV9(I  ,J+1) - QV9(I  ,J  )) &
             +          RHOW(J+1)*W9(I,J+1)* (QV( I  ,J+1) - QV( I  ,J  )) &
             +          RHOW(J  )*W( I,J  )* (QV9(I  ,J  ) - QV9(I  ,J-1)) &
             +          RHOW(J  )*W9(I,J  )* (QV( I  ,J  ) - QV( I  ,J-1)) )
        QVA9(I,J) = QVA9(I,J) &
             -DTLR(I)* (R(I+1) * U9(I+1,J) * (QV9(I+1,J  ) - QV9(I  ,J  )) &
             +          R(I  ) * U9(I  ,J) * (QV9(I  ,J  ) - QV9(I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QV9(I  ,J+1) - QV9(I  ,J  )) &
             +          RHOW(J  )*W9(I,J  )* (QV9(I  ,J  ) - QV9(I  ,J-1)) )
        
        QLA( I,J) = QLA( I,J) &
             -DTLR(I)* (R(I+1) * U( I+1,J) * (QL9(I+1,J  ) - QL9(I  ,J  )) &
             +          R(I+1) * U9(I+1,J) * (QL( I+1,J  ) - QL( I  ,J  )) &
             +          R(I  ) * U( I  ,J) * (QL9(I  ,J  ) - QL9(I-1,J  )) &
             +          R(I  ) * U9(I  ,J) * (QL( I  ,J  ) - QL( I-1,J  )) ) &
             -DTLZ(J)* (RHOW(J+1)*W( I,J+1)* (QL9(I  ,J+1) - QL9(I  ,J  )) &
             +          RHOW(J+1)*W9(I,J+1)* (QL( I  ,J+1) - QL( I  ,J  )) &
             +          RHOW(J  )*W( I,J  )* (QL9(I  ,J  ) - QL9(I  ,J-1)) &
             +          RHOW(J  )*W9(I,J  )* (QL( I  ,J  ) - QL( I  ,J-1)) )
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
     UB( J) = R(MP1)*U( MP1,J) + R(M)*U( M,J)
     UB9(J) = R(MP1)*U9(MP1,J) + R(M)*U9(M,J)
     IF(UB9(J) < 0) THEN
        UB( J) = 0
        UB9(J) = 0
     END IF
  END DO

  I = M
  DO J = 1, N
     VA( I,J) = VA( I,J) &
          -DTLR(I)* UB( J) * (V19(I  ,J  ) - V19(I-1,J  )) &
          -DTLR(I)* UB9(J) * (V1( I  ,J  ) - V1( I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W( I,J+1)* ( V9(I  ,J+1) - V9(I  ,J  )) &
          +          RHOW(J+1)*W9(I,J+1)* ( V( I  ,J+1) - V( I  ,J  )) &
          +          RHOW(J  )*W( I,J  )* ( V9(I  ,J  ) - V9(I  ,J-1)) &
          +          RHOW(J  )*W9(I,J  )* ( V( I  ,J  ) - V( I  ,J-1)) ) &
          -DTL2*(   V( I,J)/RS(I))*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I) &
          -DTL2*(F+ V9(I,J)/RS(I))*(R(I+1)*U( I+1,J)+R(I)*U( I,J))/RS(I)
     VA9(I,J) = VA9(I,J) &
          -DTLR(I)* UB9(J) * (V19(I  ,J  ) - V19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* ( V9(I  ,J+1) - V9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )* ( V9(I  ,J  ) - V9(I  ,J-1)) ) &
          -DTL2*(F+ V9(I,J)/RS(I))*(R(I+1)*U9(I+1,J)+R(I)*U9(I,J))/RS(I)
     
     TA( I,J) = TA( I,J) &
          -DTLR(I)* UB( J) * (T19(I  ,J  ) - T19(I-1,J  )) &
          -DTLR(I)* UB9(J) * (T1( I  ,J  ) - T1( I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W( I,J+1)*(T9(I  ,J+1) - T9(I  ,J  )) &
          +          RHOW(J+1)*W9(I,J+1)*(T( I  ,J+1) - T( I  ,J  )) &
          +          RHOW(J  )*W( I,J  )*(T9(I  ,J  ) - T9(I  ,J-1)) &
          +          RHOW(J  )*W9(I,J  )*(T( I  ,J  ) - T( I  ,J-1)) )! &
          !+QRAD( I,J)
     TA9(I,J) = TA9(I,J) &
          -DTLR(I)* UB9(J) * (T19(I  ,J  ) - T19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)*(T9(I  ,J+1) - T9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )*(T9(I  ,J  ) - T9(I  ,J-1)) )! &
          !+QRAD9(I,J)
     
     QVA( I,J) = QVA( I,J) &
          -DTLR(I)* UB( J) * (QV19(I  ,J  ) - QV19(I-1,J  )) &
          -DTLR(I)* UB9(J) * (QV1( I  ,J  ) - QV1( I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W( I,J+1)* (QV9(I  ,J+1) - QV9(I  ,J  )) &
          +          RHOW(J+1)*W9(I,J+1)* (QV( I  ,J+1) - QV( I  ,J  )) &
          +          RHOW(J  )*W( I,J  )* (QV9(I  ,J  ) - QV9(I  ,J-1)) &
          +          RHOW(J  )*W9(I,J  )* (QV( I  ,J  ) - QV( I  ,J-1)) )
     QVA9(I,J) = QVA9(I,J) &
          -DTLR(I)* UB9(J) * (QV19(I  ,J  ) - QV19(I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W9(I,J+1)* (QV9(I  ,J+1) - QV9(I  ,J  )) &
          +          RHOW(J  )*W9(I,J  )* (QV9(I  ,J  ) - QV9(I  ,J-1)) )
     
     QLA( I,J) = QLA( I,J) &
          -DTLR(I)* UB( J) * (QL19(I  ,J  ) - QL19(I-1,J  )) &
          -DTLR(I)* UB9(J) * (QL1( I  ,J  ) - QL1( I-1,J  )) &
          -DTLZ(J)* (RHOW(J+1)*W( I,J+1)* (QL9(I  ,J+1) - QL9(I  ,J  )) &
          +          RHOW(J+1)*W9(I,J+1)* (QL( I  ,J+1) - QL( I  ,J  )) &
          +          RHOW(J  )*W( I,J  )* (QL9(I  ,J  ) - QL9(I  ,J-1)) &
          +          RHOW(J  )*W9(I,J  )* (QL( I  ,J  ) - QL( I  ,J-1)) )
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
        QLA( I,J) = QLA( I,J) + DTLZ(J) * &
             (DUM(I,J+1)*RHOT(J+1)*QL( I,J+1)-DUM(I,J-1)*RHOT(J-1)*QL( I,J-1))
        QLA9(I,J) = QLA9(I,J) + DTLZ(J) * &
             (DUM(I,J+1)*RHOT(J+1)*QL9(I,J+1)-DUM(I,J-1)*RHOT(J-1)*QL9(I,J-1))
     END DO
  END DO
  
  DO I = 1 , M
     DUM(I,2)=0.
     IF((QL19(I,2) - 0.001) > 0) DUM(I,2)=VTERM
     DUM(I,1)=0
     IF((QL19(I,1) - 0.001) > 0) DUM(I,1)=VTERM
     QLA( I,1) = QLA( I,1) + 2*DTLZ(1)* &
          (DUM(I,2)*RHOT(2)*QL1( I,2) - DUM(I,1)*RHOT(1)*QL1( I,1 ))
     QLA9(I,1) = QLA9(I,1) + 2*DTLZ(1)* &
          (DUM(I,2)*RHOT(2)*QL19(I,2) - DUM(I,1)*RHOT(1)*QL19(I,1 ))
  END DO
  
  !
  !        TIME SMOOTHER
  !
  DO J = 1, N
     DO I = 1, M
        U( I+1,J)= U( I+1,J)+ EPS * (U1( I+1,J)-2.*U( I+1,J))
        U9(I+1,J)= U9(I+1,J)+ EPS * (U19(I+1,J)-2.*U9(I+1,J))
        
        V( I,J)  = V( I,J)  + EPS * (V1( I,J)  -2.*V( I,J)  )
        V9(I,J)  = V9(I,J)  + EPS * (V19(I,J)  -2.*V9(I,J)  )
        
        W( I,J+1)= W( I,J+1)+ EPS * (W1( I,J+1)-2.*W( I,J+1))
        W9(I,J+1)= W9(I,J+1)+ EPS * (W19(I,J+1)-2.*W9(I,J+1))
        
        P( I,J)  = P( I,J)  + EPS * (P1( I,J)  -2.*P( I,J)  )
        P9(I,J)  = P9(I,J)  + EPS * (P19(I,J)  -2.*P9(I,J)  )
        
        T( I,J)  = T( I,J)  + EPS * (T1( I,J)  -2.*T( I,J)  )
        T9(I,J)  = T9(I,J)  + EPS * (T19(I,J)  -2.*T9(I,J)  )
        
        QV( I,J) = QV( I,J) + EPS * (QV1( I,J) -2.*QV( I,J) )
        QV9(I,J) = QV9(I,J) + EPS * (QV19(I,J) -2.*QV9(I,J) )
        
        QL( I,J) = QL( I,J) + EPS * (QL1( I,J) -2.*QL( I,J) )
        QL9(I,J) = QL9(I,J) + EPS * (QL19(I,J) -2.*QL9(I,J) )
     END DO
  END DO

  !
  !         SMALL TIME STEP
  !
  DO NSMALL = 1, NS
     DO J = 1, N
        DO I = 2, M
           CFACU1 =DTS*RDR*CP*0.5 * (&
                (T1( I,J)*(1.0+0.61*QV9(I,J))+T1 (I-1,J)*(1.0+0.61*QV9(I-1,J)))&
                +(T19(I,J)*(    0.61*QV( I,J))+T19(I-1,J)*(    0.61*QV( I-1,J))))
           CFACU19=DTS*RDR*CP*0.5&
                *(T19(I,J)*(1.0+0.61*QV9(I,J))+T19(I-1,J)*(1.0+0.61*QV9(I-1,J)))
           
           U1( I,J)=U1( I,J)-CFACU1 *(P19(I,J)-P19(I-1,J)) &
                -CFACU19*(P1( I,J)-P1( I-1,J))+UA( I,J)
           U19(I,J)=U19(I,J)-CFACU19*(P19(I,J)-P19(I-1,J))+UA9(I,J)
        END DO
     END DO

     DO J = 2, N
        DO I = 1, M
           CFACWS =DTS*RDZ*CP*0.5* (&
                (T1( I,J)+T1( I,J-1)) * (1.0+0.305*(QV9(I,J)+QV9(I,J-1))) &
              +(T19(I,J)+T19(I,J-1)) * (    0.305*(QV( I,J)+QV( I,J-1))) )
           CFACWS9=DTS*RDZ*CP*0.5* &
                (T19(I,J)+T19(I,J-1)) * (1.0+0.305*(QV9(I,J)+QV9(I,J-1)))
           WS( I,J)=W1( I,J)-0.5*(1-EP)*CFACWS *(P19(I,J)-P19(I,J-1)) &
                -0.5*(1-EP)*CFACWS9*(P1( I,J)-P1( I,J-1))+WA( I,J)
           WS9(I,J)=W19(I,J)-0.5*(1-EP)*CFACWS9*(P19(I,J)-P19(I,J-1))+WA9(I,J)
        END DO
     END DO

     DO J = 1, N
        DO I = 1, M
           PS( I,J)=P1( I,J) &
                -RC2(J)*(R(I+1)*U1( I+1,J) - R(I)*U1( I,J))/RS(I) &
                -0.5*(1.0-EP)*ZC2(J) &
                *(RHOTVW(J+1)*W1( I,J+1)-RHOTVW(J)*W1( I,J))
           PS9(I,J)=P19(I,J) &
                -RC2(J)*(R(I+1)*U19(I+1,J) - R(I)*U19(I,J))/RS(I) &
                -0.5*(1.0-EP)*ZC2(J) &
                *(RHOTVW(J+1)*W19(I,J+1)-RHOTVW(J)*W19(I,J))
        END DO
     END DO

     DO J = 2, N
        DO I = 1, M
           CFACD =DTS*RDZ*CP*0.5*( &
                (T1( I,J)+T1( I,J-1))*(1.0+0.305*(QV9(I,J)+QV9(I,J-1))) &
                +(T19(I,J)+T19(I,J-1))*(    0.305*(QV( I,J)+QV( I,J-1))) )
           CFACD9=DTS*RDZ*CP*0.5* &
                (T19(I,J)+T19(I,J-1))*(1.0+0.305*(QV9(I,J)+QV9(I,J-1)))
           D( I,J) = (WS( I,J) - CFACD *0.5*(1+EP)*(PS9(I,J)-PS9(I,J-1)) &
                - CFACD9*0.5*(1+EP)*(PS( I,J)-PS( I,J-1)) &
                + C(J) * D( I,J-1)) * E(J) /A(J)
           D9(I,J) = (WS9(I,J) - CFACD9*0.5*(1+EP)*(PS9(I,J)-PS9(I,J-1)) &
                + C(J) * D9(I,J-1)) * E(J) /A(J)
        END DO
     END DO

     DO J = N, 1, -1
        DO I = 1, M
           W1( I,J) = E(J) * W1( I,J+1) + D( I,J)
           W19(I,J) = E(J) * W19(I,J+1) + D9(I,J)
           P1( I,J) = PS( I,J) -0.5*(1.0+EP)*ZC2(J)* &
                (RHOTVW(J+1)*W1( I,J+1) - RHOTVW(J)*W1( I,J))
           P19(I,J) = PS9(I,J) -0.5*(1.0+EP)*ZC2(J)* &
                (RHOTVW(J+1)*W19(I,J+1) - RHOTVW(J)*W19(I,J))
        END DO
     END DO
  END DO ! DO NSMALL = 1 , NS
  !
  !        OUTER BOUNDARY
  !
  DO J = 1 , N
     U1( MP1,J) = U1( MP1,J) + UA( MP1,J)
     U19(MP1,J) = U19(MP1,J) + UA9(MP1,J)
  END DO
  !
  !        ADVANCE V , T , QV , QL
  !
  DO J=1,N
     DO I=1,M
        V1( I,J)  = V1( I,J)  + VA( I,J)
        V19(I,J)  = V19(I,J)  + VA9(I,J)
        
        T1( I,J)  = T1( I,J)  + TA( I,J)
        T19(I,J)  = T19(I,J)  + TA9(I,J)
        
        QV1( I,J) = QV1( I,J) + QVA( I,J)
        QV19(I,J) = QV19(I,J) + QVA9(I,J)
        
        QL1( I,J) = QL1( I,J) + QLA( I,J)
        QL19(I,J) = QL19(I,J) + QLA9(I,J)
        IF(QV19(I,J) < 0) THEN
           QV1( I,J)= 0
           QV19(I,J)= 0
        END IF
        IF(QL19(I,J) < 0) THEN
           QL1( I,J)=0
           QL19(I,J)=0
        END IF
        !
        !     CONDENSATION / EVAPORATION
        !
        !
        PNST =         P1( I,J) - P1( M,1)
        PNST9= PN(J) + P19(I,J) - P19(M,1)
        
        PDST = 1000.0 * (1/XKAPPA)*PNST9**(1/XKAPPA-1)*PNST
        PDST9= 1000.0 * PNST9**(1/XKAPPA)
        
        TEMP = PNST  * T19(I,J) + PNST9 * T1(I,J)
        TEMP9= PNST9 * T19(I,J)
        
        ES   = 6.11 * EXP(A1*(TEMP9-273)/(TEMP9-36)) &
             *A1*237/(TEMP9-36)**2*TEMP
        ES9  = 6.11 * EXP(A1*(TEMP9-273)/(TEMP9-36))
        
        QSS  = 0
        QSS9 = 0.622
        IF(PDST9-ES9 > ES9) THEN
           QSS  = 0.622*ES /(PDST9-ES9) - 0.622*ES9/(PDST9-ES9)**2*(PDST-ES)
           QSS9 = 0.622*ES9/(PDST9-ES9)
        END IF
        
        !     NO EVAPORATION NOR CONDENSATION
        IF(QV19(I,J)>=QSS9 .OR. QL19(I,J)>0) THEN
           R1  = -1/(1 + XLDCP * 237 * A1 * QSS9/(TEMP9-36)**2)**2 &
                * (      XLDCP * 237 * A1 * QSS /(TEMP9-36)**2 &
                - 2 *    XLDCP * 237 * A1 * QSS9/(TEMP9-36)**3*TEMP)
           R19 =  1/(1 + XLDCP * 237 * A1 * QSS9/(TEMP9-36)**2)
           
           QVD  = R1 *(QV19(I,J)-QSS9) + R19*(QV1( I,J)-QSS)
           QVD9 = R19*(QV19(I,J)-QSS9)
           
           !     NO EVAPORATION
           IF ((QL19(I,J)+QVD9) < 0) THEN
              
              !     EVAPORATE
              T1(  I,J) = T1( I,J)  - XLDCP * QL1( I,J)/PNST9 &
                   - XLDCP * QL19(I,J)/PNST9**2*(-PNST)
              T19( I,J) = T19(I,J)  - XLDCP * QL19(I,J)/PNST9
              
              QV1( I,J) = QV1( I,J) + QL1( I,J)
              QV19(I,J) = QV19(I,J) + QL19(I,J)
              
              QL1( I,J) = 0
              QL19(I,J) = 0
           ELSE
              !     CONDENSE
              T1( I,J)  = T1( I,J)  + XLDCP * QVD /PNST9 &
                   + XLDCP * QVD9/PNST9**2*(-PNST)
              T19(I,J)  = T19(I,J)  + XLDCP * QVD9/PNST9
              
              QV1( I,J) = QV1( I,J)         - QVD
              QV19(I,J) = QV19(I,J)         - QVD9
              QL1( I,J) = QL1( I,J)         + QVD
              QL19(I,J) = QL19(I,J)         + QVD9
           END IF
        END IF
        IF(QV19(I,J) < 0) THEN
           QV1( I,J)=0
           QV19(I,J)=0
        END IF
        IF(QL19(I,J) < 0) THEN
           QL1( I,J)=0
           QL19(I,J)=0
        END IF
     END DO
  END DO
  
  !
  !        TIME FLIP
  !
  DO J=1,N
     DO I=1,M
        UTEMP     = U( I+1,J) + EPS * U1( I+1,J)
        UTEMP9    = U9(I+1,J) + EPS * U19(I+1,J)
        VTEMP     = V( I,J)   + EPS * V1( I,J)
        VTEMP9    = V9(I,J)   + EPS * V19(I,J)
        WTEMP     = W( I,J+1) + EPS * W1( I,J+1)
        WTEMP9    = W9(I,J+1) + EPS * W19(I,J+1)
        PTEMP     = P( I,J)   + EPS * P1( I,J)
        PTEMP9    = P9(I,J)   + EPS * P19(I,J)
        TTEMP     = T( I,J)   + EPS * T1( I,J)
        TTEMP9    = T9(I,J)   + EPS * T19(I,J)
        QVTEMP    = QV( I,J)  + EPS * QV1( I,J)
        QVTEMP9   = QV9(I,J)  + EPS * QV19(I,J)
        QLTEMP    = QL( I,J)  + EPS * QL1( I,J)
        QLTEMP9   = QL9(I,J)  + EPS * QL19(I,J)
        
        
        U( I+1,J) = U1( I+1,J)
        U9(I+1,J) = U19(I+1,J)
        V( I,J)   = V1( I,J)
        V9(I,J)   = V19(I,J)
        W( I,J+1) = W1( I,J+1)
        W9(I,J+1) = W19(I,J+1)
        P( I,J)   = P1( I,J)
        P9(I,J)   = P19(I,J)
        T( I,J)   = T1( I,J)
        T9(I,J)   = T19(I,J)
        QV( I,J)  = QV1( I,J)
        QV9(I,J)  = QV19(I,J)
        QL( I,J)  = QL1( I,J)
        QL9(I,J)  = QL19(I,J)
        
        
        U1( I+1,J)= UTEMP
        U19(I+1,J)= UTEMP9
        V1( I,J)  = VTEMP
        V19(I,J)  = VTEMP9
        W1( I,J+1)= WTEMP
        W19(I,J+1)= WTEMP9
        T1( I,J)  = TTEMP
        T19(I,J)  = TTEMP9
        P1( I,J)  = PTEMP
        P19(I,J)  = PTEMP9
        QV1( I,J) = QVTEMP
        QV19(I,J) = QVTEMP9
        QL1( I,J) = QLTEMP
        QL19(I,J) = QLTEMP9
     END DO
  END DO
  
  RETURN
END SUBROUTINE hurr_tlm

! ==========================================================

SUBROUTINE DIFFUSE_TLM
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

  DO J = 1, N
     DO I = 1, M
        TE( I,J)= T1( I,J) &
             *(1+XLDCP*QV19(I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1)))) &
             +    T19(I,J) * ( &
              (  XLDCP*QV1( I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1)))) &
             -   XLDCP*QV19(I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1)))**2 &
             *(  (T1( I,J)*(PN(J)+P19(I,J)-P19(M,1))) &
               + (T19(I,J)*(      P1( I,J)-P1( M,1))) ) )
        TE9(I,J)= T19(I,J) &
             *(1+XLDCP*QV19(I,J)/(T19(I,J)*(PN(J)+P19(I,J)-P19(M,1))))
     END DO
  END DO

  DO J = 1, N
     DO I = 1, M
        UDR( I,J) = RDR * (U1( I+1,J) - U1( I,J))
        UDR9(I,J) = RDR * (U19(I+1,J) - U19(I,J))
        WDZ( I,J) = RDZ * (W1( I,J+1) - W1( I,J))
        WDZ9(I,J) = RDZ * (W19(I,J+1) - W19(I,J))
     END DO
  END DO

  DO J = 1, N
     DO I = 2, M
        VDR( I,J) = R(I)*RDR * (V1( I,J)/RS(I) - V1( I-1,J)/RS(I-1))
        VDR9(I,J) = R(I)*RDR * (V19(I,J)/RS(I) - V19(I-1,J)/RS(I-1))
     END DO
     VDR( 1  ,J) = 0
     VDR9(1  ,J) = 0
     VDR( MP1,J) = VDR( M,J)
     VDR9(MP1,J) = VDR9(M,J)
  END DO

  DO J = 2, N
     DO I = 2, M
        UZWR( I,J) = RDZ*(U1( I,J)-U1( I,J-1)) + RDR*(W1( I,J)-W1( I-1,J))
        UZWR9(I,J) = RDZ*(U19(I,J)-U19(I,J-1)) + RDR*(W19(I,J)-W19(I-1,J))
     END DO
     UZWR( 1  ,J) = 0
     UZWR9(1  ,J) = 0
     UZWR( MP1,J) = UZWR( M,J)
     UZWR9(MP1,J) = UZWR9(M,J)
  END DO

  DO J = 2, N
     DO I = 1, M
        VDZ( I,J) = RDZ * (V1( I,J) - V1( I,J-1))
        VDZ9(I,J) = RDZ * (V19(I,J) - V19(I,J-1))
     END DO
  END DO
  DO J = 2, N
     DO I = 1, M
        DEFH2( I,J) = 2*UDR9(I,J)*UDR(I,J) + 2*UDR9(I,J-1)*UDR(I,J-1) &
             + (0.25/(RS(I)**2)) * (&
               2*(U19(I+1,J  )+U19(I,J  ))*(U1(I+1,J  )+U1(I,J  )) &
             + 2*(U19(I+1,J-1)+U19(I,J-1))*(U1(I+1,J-1)+U1(I,J-1)) ) &
             + 0.25* ( &
               2*VDR9(I+1,J  )*VDR(I+1,J  ) + 2*VDR9(I,J  )*VDR(I,J  ) &
             + 2*VDR9(I+1,J-1)*VDR(I+1,J-1) + 2*VDR9(I,J-1)*VDR(I,J-1) )
        DEFH29(I,J) = UDR9(I,J)**2 + UDR9(I,J-1)**2 &
             + (0.25/(RS(I)**2)) &
             *((U19(I+1,J  )+U19(I,J  ))**2 + (U19(I+1,J-1)+U19(I,J-1))**2) &
             +0.25* &
             (VDR9(I+1,J  )**2+VDR9(I,J  )**2+VDR9(I+1,J-1)**2+VDR9(I,J-1)**2)
     END DO
  END DO

  DO J = 2, N
     DO I = 1, M
        DEF2( I,J) = DEFH2( I,J) &
             + 2*WDZ9(I,J)*WDZ(I,J) + 2*WDZ9(I,J-1)*WDZ(I,J-1) + 0.5 &
             *(2*UZWR9(I+1,J)*UZWR(I+1,J) + 2*UZWR9(I,J)*UZWR(I,J))  &
             + 2*VDZ9(I,J)*VDZ(I,J)
        DEF29(I,J) = DEFH29(I,J) + WDZ9(I,J)**2 + WDZ9(I,J-1)**2 &
             + 0.5*(UZWR9(I+1,J)**2 + UZWR9(I,J)**2) + VDZ9(I,J)**2
     END DO
  END DO

  DO J = 2, N
     DO I = 1, M
        DTDZ( I,J) = RDZ*( &
             2.0 * RTBW( J) * (T1(  I,J) - T1(  I,J-1)) &
            +2.0 * RQVBW(J) * (QV1( I,J) - QV1( I,J-1)) )
        DTDZ9(I,J) = RDZ*( &
             2.0 * RTBW( J) * (T19( I,J) - T19( I,J-1)) &
            +2.0 * RQVBW(J) * (QV19(I,J) - QV19(I,J-1)) )
        IF( (QL19(I,J) * QL19(I,J-1)) >= 1.0E-8 ) THEN
           TT  = 0.5 &
                *(T1( I,J  ) * ( PN(J  )+P19(I,J  )-P19(M,1)) &
                + T19(I,J  ) * (         P1( I,J  )-P1( M,1)) &
                + T1( I,J-1) * ( PN(J-1)+P19(I,J-1)-P19(M,1)) &
                + T19(I,J-1) * (         P1( I,J-1)-P1( M,1)) )
           TT9 = 0.5 &
                *(T19(I,J  ) * ( PN(J  )+P19(I,J  )-P19(M,1)) &
                + T19(I,J-1) * ( PN(J-1)+P19(I,J-1)-P19(M,1)) ) 
           AA =(       4362    *(QV1( I,J)+QV1( I,J-1))/TT9 &
                      -4362    *(QV19(I,J)+QV19(I,J-1))/TT9**2*TT ) &
                / (1 + 6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2)) &
                * (2 / (TB(J)+TB(J-1)) ) &
              -(   1 + 4362    *(QV19(I,J)+QV19(I,J-1)) / TT9   ) &
                / (1 + 6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2))**2 &
                * (    6738953 *(QV1( I,J)+QV1( I,J-1))/(TT9**2) &
                    -2*6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**3)*TT) &
                * (2 / (TB(J)+TB(J-1)) )
           AA9=(   1 + 4362    *(QV19(I,J)+QV19(I,J-1)) / TT9   ) &
                / (1 + 6738953 *(QV19(I,J)+QV19(I,J-1))/(TT9**2)) &                
                * (2 / (TB(J)+TB(J-1)) )           
           DTDZ( I,J) = RDZ* G * ( AA  * (TE9(I,J) - TE9(I,J-1)) &
                                  +AA9 * (TE( I,J) - TE( I,J-1)) &
                - (QL1( I,J)+QV1( I,J)-QL1( I,J-1)-QV1( I,J-1)) )
           DTDZ9(I,J) = RDZ* G * ( AA9 * (TE9(I,J) - TE9(I,J-1)) &
                - (QL19(I,J)+QV19(I,J)-QL19(I,J-1)-QV19(I,J-1)) )
        END IF
     END DO
  END DO

  DO J = 2, N
     DO I = 1, M
        XKM( I,J) = 0
        XKM9(I,J) = 0

        IF(DEF29(I,J) > DTDZ9(I,J)) THEN
           XKM( I,J)= XVL2 * 0.5/SQRT(DEF29(I,J) - DTDZ9(I,J)) &
                *(DEF2(I,J) - DTDZ(I,J))
           XKM9(I,J)= XVL2 * SQRT(DEF29(I,J) - DTDZ9(I,J))
        END IF

        IF(DEFH29(I,J)>0) THEN
           XKMH(I,J) = XHL2 * 0.5/SQRT(DEFH29(I,J))*DEFH2(I,J)
        ELSE
           XKMH(I,J) = 0
        END IF
        XKMH9(I,J) = XHL2 * SQRT(DEFH29(I,J))

        IF(XKM9(I,J) >= 0.4*DZ**2/DTL) THEN
           XKM( I,J) = 0
           XKM9(I,J) = 0.4 *DZ**2/DTL 
        END IF
        IF(XKMH9(I,J) < XKM9(I,J)) THEN
           XKMH( I,J) = XKM( I,J)
           XKMH9(I,J) = XKM9(I,J)
        END IF
     END DO
  END DO

  DO I = 1, M
     XKM(  I,1  ) = XKM(  I,2)
     XKM9( I,1  ) = XKM9( I,2)
     XKM(  I,NP1) = XKM(  I,N)
     XKM9( I,NP1) = XKM9( I,N)
     XKMH( I,1  ) = XKMH( I,2)
     XKMH9(I,1  ) = XKMH9(I,2)
     XKMH( I,NP1) = XKMH( I,N)
     XKMH9(I,NP1) = XKMH9(I,N)
  END DO
  !
  !        CALCULATE STRESS       
  !
  !        TRR(M,N)       
  !
  DO J = 1 , N  
     DO I = 1 , M  
        TRR( I,J) = (XKMH( I,J+1) + XKMH( I,J)) * UDR9(I,J) &
                   +(XKMH9(I,J+1) + XKMH9(I,J)) * UDR( I,J)
        TRR9(I,J) = (XKMH9(I,J+1) + XKMH9(I,J)) * UDR9(I,J)
     END DO
  END DO
  !
  !        TTT(M,N)
  !
  DO J = 1, N
     TTT( 1,J) = 0
     TTT9(1,J) = 0
     DO I = 2, M
        TTT( I,J) = U1( I,J)/R(I) * 0.5 &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J)) &
             +      U19(I,J)/R(I) * 0.5 &
             *(XKMH( I,J+1)+XKMH( I,J)+XKMH( I-1,J+1)+XKMH( I-1,J))
        TTT9(I,J) = U19(I,J)/R(I) * 0.5 &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
  END DO
  !
  !        TRZ(M,NP1)
  !
  DO J=2, N
     DO I=2, M
        TRZ( I,J)=0.5* (XKM( I-1,J) + XKM( I,J)) * UZWR9(I,J) &
             +    0.5* (XKM9(I-1,J) + XKM9(I,J)) * UZWR( I,J)
        TRZ9(I,J)=0.5* (XKM9(I-1,J) + XKM9(I,J)) * UZWR9(I,J)
     END DO
  END DO

  DO I =2, M
     IF(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2>0) THEN
        TRZ( I,1)= 0.5 * (CDRS( I)+CDRS( I-1)) * U19(I,1) &
          * SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
          +     0.5 * (CDRS9(I)+CDRS9(I-1)) * U1( I,1) &
          * SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
          +   0.5 * (CDRS9(I)+CDRS9(I-1)) * U19(I,1) * 0.5 &
          / SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) &
          * (2*U19(I,1)*U1(I,1)+0.25*2*(V19(I,1)+V19(I-1,1))*(V1(I,1)+V1(I-1,1)) )
     ELSE
        TRZ(I,1)=0
     END IF
     TRZ9(I,1)= 0.5 * (CDRS9(I)+CDRS9(I-1)) * U19(I,1) &
          * SQRT(U19(I,1)**2+0.25*(V19(I,1)+V19(I-1,1))**2) 
     TRZ( I,NP1) = 0
     TRZ9(I,NP1) = 0
  END DO

  DO J = 1, NP1
     TRZ( 1  ,J) = 0
     TRZ9(1  ,J) = 0
     TRZ( MP1,J) = TRZ( M,J) * R(M) / R(MP1)
     TRZ9(MP1,J) = TRZ9(M,J) * R(M) / R(MP1)
  END DO
  !
  !        TRT(M,N)
  !
  DO J=1, N
     DO I=2, M
        TRT( I,J) = 0.25*VDR( I,J) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J)) &
             +      0.25*VDR9(I,J) &
             *(XKMH( I,J+1)+XKMH( I,J)+XKMH( I-1,J+1)+XKMH( I-1,J))
        TRT9(I,J) = 0.25*VDR9(I,J) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
     TRT( 1  ,J) = 0
     TRT9(1  ,J) = 0
     TRT( MP1,J) = TRT( M,J) * R(M)**2 / R(MP1)**2
     TRT9(MP1,J) = TRT9(M,J) * R(M)**2 / R(MP1)**2
  END DO
  !       
  !        TZT(M,NP1)     
  !       
  DO I = 1, M
     DO J = 2, N
        TZT( I,J) = XKM( I,J)*VDZ9(I,J) + XKM9(I,J)*VDZ(I,J)
        TZT9(I,J) = XKM9(I,J)*VDZ9(I,J)
     END DO
     IF(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2>0) THEN
        TZT( I,1  )= CDRS( I) * V19(I,1) &
             *SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             +       CDRS9(I) * V1( I,1) &
             *SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             +       CDRS9(I) * V19(I,1) * 0.5 &
             /SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2) &
             * (0.25*2*(U19(I,1)+U19(I+1,1))*(U1(I,1)+U1(I+1,1)) &
             + 2*V19(I,1)*V1(I,1) )
     ELSE
        TZT(I,1) = 0
     END IF
     TZT9(I,1  )= CDRS9(I) * V19(I,1) &
          *SQRT(0.25*(U19(I,1)+U19(I+1,1))**2+V19(I,1)**2)
     TZT( I,NP1)= 0
     TZT9(I,NP1)= 0
  END DO
  !
  !          TZZ(M,N)
  !
  DO J = 1, N
     DO I = 1, M
        TZZ( I,J) = (XKM( I,J+1) + XKM( I,J)) * WDZ9(I,J) &
             +      (XKM9(I,J+1) + XKM9(I,J)) * WDZ( I,J)
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
        TR( I,J) = 0.25 * RDR*(T1( I,J)-T1( I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J)) &
             +     0.25 * RDR*(T19(I,J)-T19(I-1,J)) &
             *(XKMH( I,J+1)+XKMH( I,J)+XKMH( I-1,J+1)+XKMH( I-1,J))
        TR9(I,J) = 0.25 * RDR*(T19(I,J)-T19(I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
     TR( 1  ,J)=0
     TR9(1  ,J)=0
     TR( MP1,J)= TR( M,J) * R(M) / R(MP1)
     TR9(MP1,J)= TR9(M,J) * R(M) / R(MP1)
  END DO
  !
  !         TZ(M,NP1)
  !
  DO I=1,M
     DO J=2,N
        TZ( I,J) = XKM( I,J) * RDZ * (T19(I,J) - T19(I,J-1)) &
             +     XKM9(I,J) * RDZ * (T1( I,J) - T1( I,J-1))
        TZ9(I,J) = XKM9(I,J) * RDZ * (T19(I,J) - T19(I,J-1))        
        TZ( I,J) = TZ( I,J) * PD(J)
        TZ9(I,J) = TZ9(I,J) * PD(J)
     END DO
     IF(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2>0) THEN
        VABS =0.5/SQRT(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2) &
             *(0.25*2*(U19(I+1,1)+U19(I,1))*(U1( I+1,1)+U1( I,1)) &
             + 2*V19(I,1)*V1(I,1) )
     ELSE
        VABS=0
     END IF
     VABS9=SQRT(0.25 * (U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2)
     TZ( I,1  ) = (T1( I,1)-TSURF( I))*CERS9(I)*VABS9 &
          +       (T19(I,1)-TSURF9(I))*CERS( I)*VABS9 &
          +       (T19(I,1)-TSURF9(I))*CERS9(I)*VABS
     TZ9(I,1  ) = (T19(I,1)-TSURF9(I))*CERS9(I)*VABS9
     TZ( I,1  ) =  TZ( I,1)*PD(1)
     TZ9(I,1  ) =  TZ9(I,1)*PD(1)
     TZ( I,NP1) = 0
     TZ9(I,NP1) = 0
  END DO
  !
  !       QV FLUX 
  !
  !       QVR(MP1,N)
  !
  DO J=1,N
     DO I=2,M
        QVR( I,J)= 0.25 *RDR*(QV1( I,J)-QV1( I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J)) &
             +     0.25 *RDR*(QV19(I,J)-QV19(I-1,J)) &
             *(XKMH( I,J+1)+XKMH( I,J)+XKMH( I-1,J+1)+XKMH( I-1,J))
        QVR9(I,J)= 0.25 *RDR*(QV19(I,J)-QV19(I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
     QVR( 1  ,J) = 0
     QVR9(1  ,J) = 0
     QVR( MP1,J) = QVR( M,J) * R(M) / R(MP1)
     QVR9(MP1,J) = QVR9(M,J) * R(M) / R(MP1)
  END DO
  !
  !         QVZ(M,NP1)
  !
  DO I=1,M
     DO J=2,N
        QVZ( I,J) = XKM( I,J) * RDZ * (QV19(I,J) - QV19(I,J-1)) &
             +      XKM9(I,J) * RDZ * (QV1( I,J) - QV1( I,J-1))
        QVZ9(I,J) = XKM9(I,J) * RDZ * (QV19(I,J) - QV19(I,J-1))
        QVZ( I,J) = QVZ( I,J) * PD(J)/PN(J)
        QVZ9(I,J) = QVZ9(I,J) * PD(J)/PN(J)
     END DO
     IF(0.25*(U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2>0) THEN
        VABS =0.5/SQRT(0.25*(U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2) &
             *(0.25*2*(U19(I+1,1)+U19(I,1))*(U1( I+1,1)+U1( I,1)) &
             + 2*V19(I,1)*V1(I,1) )
     ELSE
        VABS=0
     END IF
     VABS9=SQRT(0.25*(U19(I+1,1)+U19(I,1))**2 + V19(I,1)**2)
     QVZ( I,1  ) = (QV1( I,1)-QSURF( I))*CERS9(I)*VABS9 &
          +        (QV19(I,1)-QSURF9(I))*CERS( I)*VABS9 &
          +        (QV19(I,1)-QSURF9(I))*CERS9(I)*VABS
     QVZ9(I,1  ) = (QV19(I,1)-QSURF9(I))*CERS9(I)*VABS9
     QVZ( I,1  ) =  QVZ( I,1)*PD(1)/PN(1)
     QVZ9(I,1  ) =  QVZ9(I,1)*PD(1)/PN(1)
     QVZ( I,NP1) =  0
     QVZ9(I,NP1) =  0
  END DO
  !
  !       QL FLUX
  !
  !       QLR(MP1,N)
  !
  DO J=1, N
     DO I=2, M
        QLR( I,J)= 0.25 *RDR*(QL1( I,J)-QL1( I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J)) &
             +     0.25 *RDR*(QL19(I,J)-QL19(I-1,J)) &
             *(XKMH( I,J+1)+XKMH( I,J)+XKMH( I-1,J+1)+XKMH( I-1,J))
        QLR9(I,J)= 0.25 *RDR*(QL19(I,J)-QL19(I-1,J)) &
             *(XKMH9(I,J+1)+XKMH9(I,J)+XKMH9(I-1,J+1)+XKMH9(I-1,J))
     END DO
     QLR( 1  ,J) = 0
     QLR9(1  ,J) = 0
     QLR( MP1,J) = QLR( M,J) * R(M) / R(MP1)
     QLR9(MP1,J) = QLR9(M,J) * R(M) / R(MP1)
  END DO
  !
  !         QLZ(M,NP1)
  !
  DO I=1, M
     DO J=2, N
        QLZ( I,J) = XKM( I,J) * RDZ * (QL19(I,J) - QL19(I,J-1)) &
             +      XKM9(I,J) * RDZ * (QL1( I,J) - QL1( I,J-1))
        QLZ9(I,J) = XKM9(I,J) * RDZ * (QL19(I,J) - QL19(I,J-1))
     END DO
     QLZ( I,1  ) = 0
     QLZ9(I,1  ) = 0
     QLZ( I,NP1) = 0
     QLZ9(I,NP1) = 0
  END DO

  DO J = 1, N
     DO I = 2, M
        UA( I,J)= -DTS*( &
             - RDR * (RS(I)*TRR( I,J) - RS(I-1)*TRR( I-1,J)) / R(I) &
             + TTT( I,J) / R(I) &
             - RDZ * (TRZ( I,J+1) - TRZ( I,J)) ) &
             + DTS * TAU(J) * U1( I,J)

        UA9(I,J)= -DTS*( &
             - RDR * (RS(I)*TRR9(I,J) - RS(I-1)*TRR9(I-1,J)) / R(I) &
             + TTT9(I,J) / R(I) &
             - RDZ * (TRZ9(I,J+1) - TRZ9(I,J)) ) &
             + DTS * TAU(J) * U19(I,J)
     END DO
  END DO

  DO J = 2, N
     DO I = 1, M
        WA( I,J) = -DTS*( &
             - RDR * (R(I+1)*TRZ( I+1,J) - R(I)*TRZ( I,J)) / RS(I) &
             - RDZ * (TZZ( I,J) - TZZ( I,J-1)) ) &
             + DTS * 0.5 * (TAU(J) + TAU(J-1)) * W1( I,J)
        WA9(I,J) = -DTS*( &
             - RDR * (R(I+1)*TRZ9(I+1,J) - R(I)*TRZ9(I,J)) / RS(I) &
             - RDZ * (TZZ9(I,J) - TZZ9(I,J-1)) ) &
             + DTS * 0.5 * (TAU(J) + TAU(J-1)) * W19(I,J)
     END DO
  END DO

  DO I = 1 , M
     DO J = 1 , N-1
        VA( I,J) = -DTL * ( &
             - RDR * (R(I+1)**2*TRT( I+1,J)-R(I)**2*TRT( I,J)) / (RS(I)**2) &
             - RDZ * (TZT( I,J+1) - TZT( I,J)) ) &
             + DTL * TAU(J) * V1( I,J)
        VA9(I,J) = -DTL * ( &
             - RDR * (R(I+1)**2*TRT9(I+1,J)-R(I)**2*TRT9(I,J)) / (RS(I)**2) &
             - RDZ * (TZT9(I,J+1) - TZT9(I,J)) ) &
             + DTL * TAU(J) * V19(I,J)

        TA( I,J) = -DTL * ( &
             - RDR * (R(I+1)*TR( I+1,J) - R(I)*TR( I,J)) / RS(I) &
             - RDZ * (TZ( I,J+1) - TZ( I,J))/(0.5*(PD(J)+PD(J+1))) ) &
             + DTL * TAU(J) * T1( I,J)
        TA9(I,J) = -DTL * ( &
             - RDR * (R(I+1)*TR9(I+1,J) - R(I)*TR9(I,J)) / RS(I) &
             - RDZ * (TZ9(I,J+1) - TZ9(I,J))/(0.5*(PD(J)+PD(J+1))) ) &
             + DTL * TAU(J) * (T19(I,J) - TB(J))

        QVA( I,J) = -DTL * ( &
             - RDR * (R(I+1)*QVR( I+1,J) - R(I)*QVR( I,J)) / RS(I) &
             - RDZ * (QVZ( I,J+1) - QVZ( I,J) )/ &
             (0.5*(PD(J)/PN(J)+PD(J+1)/PN(J+1))) ) &
             + DTL * TAU(J) * QV1( I,J)
        QVA9(I,J) = -DTL * ( &
             - RDR * (R(I+1)*QVR9(I+1,J) - R(I)*QVR9(I,J)) / RS(I) &
             - RDZ * (QVZ9(I,J+1) - QVZ9(I,J) )/ &
             (0.5*(PD(J)/PN(J)+PD(J+1)/PN(J+1))) ) &
             + DTL * TAU(J) * QV19(I,J)

        QLA( I,J) = -DTL * ( &
             - RDR * (R(I+1)*QLR( I+1,J) - R(I)*QLR( I,J)) / RS(I) &
             - RDZ * (QLZ( I,J+1) - QLZ( I,J)) ) &
             + DTL * TAU(J) * QL1( I,J)
        QLA9(I,J) = -DTL * ( &
             - RDR * (R(I+1)*QLR9(I+1,J) - R(I)*QLR9(I,J)) / RS(I) &
             - RDZ * (QLZ9(I,J+1) - QLZ9(I,J)) ) &
             + DTL * TAU(J) * QL19(I,J)
     END DO
     VA( I,J) = -DTL * ( &
          - RDR * (R(I+1)**2*TRT( I+1,J)-R(I)**2*TRT( I,J)) / (RS(I)**2) &
          - RDZ * (TZT( I,J+1) - TZT( I,J)) ) &
          + DTL * TAU(J) * V1( I,J)
     VA9(I,J) = -DTL * ( &
          - RDR * (R(I+1)**2*TRT9(I+1,J)-R(I)**2*TRT9(I,J)) / (RS(I)**2) &
          - RDZ * (TZT9(I,J+1) - TZT9(I,J)) ) &
          + DTL * TAU(J) * V19(I,J)

     TA( I,J) = -DTL * ( &
          - RDR * (R(I+1)*TR( I+1,J) - R(I)*TR( I,J)) / RS(I) &
          - RDZ * (TZ( I,J+1) - TZ( I,J) )/PD(J) ) &
          + DTL * TAU(J)*T1( I,J)
     TA9(I,J) = -DTL * ( &
          - RDR * (R(I+1)*TR9(I+1,J) - R(I)*TR9(I,J)) / RS(I) &
          - RDZ * (TZ9(I,J+1) - TZ9(I,J) )/PD(J) ) &
          + DTL * TAU(J)*(T19(I,J) - TB(J))

     QVA( I,J) = -DTL * ( &
          - RDR * (R(I+1)*QVR( I+1,J) - R(I)*QVR( I,J) ) / RS(I) &
          - RDZ * (QVZ( I,J+1) - QVZ( I,J)) / (PD(J)/PN(J)) ) &
          + DTL * TAU(J) * QV1( I,J)
     QVA9(I,J) = -DTL * ( &
          - RDR * (R(I+1)*QVR9(I+1,J) - R(I)*QVR9(I,J) ) / RS(I) &
          - RDZ * (QVZ9(I,J+1) - QVZ9(I,J)) / (PD(J)/PN(J)) ) &
          + DTL * TAU(J) * QV19(I,J)

     QLA( I,J) = -DTL * ( &
          - RDR * (R(I+1)*QLR( I+1,J) - R(I)*QLR( I,J)) / RS(I) &
          - RDZ * (QLZ( I,J+1) - QLZ( I,J)) ) &
          + DTL * TAU(J) * QL1( I,J)
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
        TAD( I,J) = DTL / (PN(J)*1005.0) &
             *(TZT( I,J+1)*VDZ9(I,J+1) + TZT9(I,J+1)*VDZ( I,J+1) &
             + XKM( I,J)*0.25*(UZWR9(I,J) + UZWR9(I+1,J))**2 &
             + XKM9(I,J)*0.25*(UZWR9(I,J) + UZWR9(I+1,J))*2 &
             *                (UZWR( I,J) + UZWR( I+1,J)) &
             + XKMH( I,J+1)*(1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))**2 &
             + XKMH9(I,J+1)*(1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))*2 &
             *(VDR( I,J)+VDR( I+1,J)+VDR( I,J+1)+VDR( I+1,J+1)) &
             + XKMH( I,J+1)*0.25*(UDR9(I,J)+UDR9(I,J+1))**2 &
             + XKMH9(I,J+1)*0.25*(UDR9(I,J)+UDR9(I,J+1))*2 &
             *                   (UDR( I,J)+UDR( I,J+1)) )
        TAD9(I,J) = DTL / (PN(J)*1005.0) &
             *(TZT9(I,J+1)*VDZ9(I,J+1) &
             + XKM9(I,J)*0.25*(UZWR9(I,J) + UZWR9(I+1,J))**2 &
             + XKMH9(I,J+1)*(1.0/16.0) &
             *(VDR9(I,J)+VDR9(I+1,J)+VDR9(I,J+1)+VDR9(I+1,J+1))**2 &
             + XKMH9(I,J+1)*0.25*(UDR9(I,J)+UDR9(I,J+1))**2 )
        TA( I,J) = TA( I,J)+TAD( I,J)
        TA9(I,J) = TA9(I,J)+TAD9(I,J)        
     END DO
     TAD( I,1)=DTL /(PN(1)*1005.0) &
          *( RDZ*CDRS( I)*(0.25*(U19(I,1)+U19(I+1,1))**2 + V19(I,1)**2)**1.5 &
          +  RDZ*CDRS9(I)*(0.25*(U19(I,1)+U19(I+1,1))**2 + V19(I,1)**2)**0.5 &
          *1.5*(0.25*2*(U19(I,1)+U19(I+1,1))*(U1(I,1)+U1(I+1,1)) &
          + 2*V19(I,1)*V1(I,1)) &
          + XKMH( I,2)*(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))**2 &
          + XKMH9(I,2)*(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))*2 &
          *(VDR( I,1)+VDR( I+1,1)+VDR( I,2)+VDR( I+1,2)) &
          + XKMH( I,2)*0.25 * (UDR9(I,1)+UDR9(I,2))**2 &
          + XKMH9(I,2)*0.25 * (UDR9(I,1)+UDR9(I,2))*2*(UDR(I,1)+UDR(I,2)) )
     TAD9(I,1)=DTL /(PN(1)*1005.0) &
          *( RDZ*CDRS9(I)*(0.25*(U19(I,1)+U19(I+1,1))**2 + V19(I,1)**2)**1.5 &
          + XKMH9(I,2)*(1.0/16.0) &
          *(VDR9(I,1)+VDR9(I+1,1)+VDR9(I,2)+VDR9(I+1,2))**2 &
          + XKMH9(I,2)*0.25 * (UDR9(I,1)+UDR9(I,2))**2 )
     TA( I,1) = TA( I,1)+TAD( I,1)
     TA9(I,1) = TA9(I,1)+TAD9(I,1)
  END DO

  RETURN
END SUBROUTINE DIFFUSE_TLM
