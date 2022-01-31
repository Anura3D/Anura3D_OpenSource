*USER SUBROUTINES
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

       ! Arguments:
        !          I/O  Type
        !  PROPS    I   R()  : List with model parameters
        !  DSTRAN   I   R()  : Strain increment
        !  DDSDDE   O   R(,) : Material stiffness matrix
        !  STRESS  I/O  R()  : stresses
        !  STATEV  I/O  R()  : state variables
 
        !---  Local variables
      DIMENSION dSig(NTENS), Sig(NTENS)
      DIMENSION SigC(NTENS)
      DIMENSION SigE(NTENS),SigEQ(NTENS)
      Parameter (Ip=-1)
      DIMENSION xN1(3),xN2(3),xN3(3),Tmp1(3),Tmp2(3),Tmp3(3)
      IAPEX=0
      ITENS=0
      TauMax=1
        ! Contents of PROPS(2)
        !  1 : G       shear modulus
        !  2 : ENU     Poisson's ratio
        IntGlo= NPT
        G = PROPS(1)   
        ENU = PROPS(2)
        VNU = PROPS(2)
        one = 1.0d0
        two = 2.0d0 
        
        SPHI = PROPS(3)
        SPSI = PROPS(5)
        COHS = PROPS(4)
        TENS = PROPS(6)
        ! calculate elastic stress increment (DSig = elastic stiffness D * strain increment DEps)
        FAC = two * G / ( one - two * ENU )
        D1 = FAC * ( one - ENU )
        D2 = FAC * ENU
        DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
        dSig(1) = (D1 - D2) * DSTRAN(1) + D2 * DSTRANVOL
        dSig(2) = (D1 - D2) * DSTRAN(2) + D2 * DSTRANVOL
        dSig(3) = (D1 - D2) * DSTRAN(3) + D2 * DSTRANVOL
        dSig(4) = G * DSTRAN(4)
        if (NTENS == 6) then
        dSig(5) = G * DSTRAN(5)
        dSig(6) = G * DSTRAN(6)
        end if
      
        ! elastic stress
        SigE = STRESS + dSig
        Stress=SigE
        SigEQ = SigE
        SigC = SigE
        
        DDSDDE = 0.0
        DDSDDE(1:3,1:3) = D2
        DDSDDE(1,1) = D1
        DDSDDE(2,2) = D1
        DDSDDE(3,3) = D1
        DDSDDE(4,4) = G
        if (NTENS == 6) then
          DDSDDE(5,5) = G
          DDSDDE(6,6) = G
        end if
      

      
!------------ Calculate principle stresses and their direction

      call PrnSig(1,NTENS,SigE,xN1,xN2,xN3,Sig1,Sig2,Sig3,P,Q)

  901       format(a,6(1x,1p,e12.5),1x,i4)



!------------- evaluating the yield functions
      IF (SPHI >  0.999d0) SPHI=0.999d0
      F21 = 0.5d0*(SIG2-SIG1) + 0.5d0*(SIG2+SIG1)*SPHI - COHS
      F32 = 0.5d0*(SIG3-SIG2) + 0.5d0*(SIG3+SIG2)*SPHI - COHS
      F31 = 0.5d0*(SIG3-SIG1) + 0.5d0*(SIG3+SIG1)*SPHI - COHS
      FT1 = SIG1 - TENS
      FT2 = SIG2 - TENS
      FT3 = SIG3 - TENS
      IF (F31 > 1d-12 .OR. FT3 > 1d-12) GOTO 240
      IF (IPL  == 0) GOTO 360
!
!         ***POINTS CHANGING FROM PLASTIC TO ELASTIC STATE***
!
      IPL =0
      GOTO 360
!
!                   *** Plastic stress points ***
!
  240 IPL =1


     
      call PrnSig( 0,NTENS,SigEQ,Tmp1,Tmp2,Tmp3,SigV1,SigV2,SigV3,Dum1,Dum2)
   

      FF1 = 0.5d0*    (SIGV2-SIGV1) + 0.5d0*(SIGV2+SIGV1)*SPHI - COHS
      FF2 = 0.5d0*DABS(SIGV2-SIGV3) + 0.5d0*(SIGV2+SIGV3)*SPHI - COHS
      FF3 = 0.5d0*DABS(SIGV3-SIGV1) + 0.5d0*(SIGV3+SIGV1)*SPHI - COHS
      


      IF (FF1 > 1d-12 .OR. FF2 > 1d-12 .OR. FF3 > 1d-12) GOTO 250

!       ****** PLASTIC NEGATIVE POINTS ******
      IPL =-1
!
!         *** STRESSES ARE BROUGHT BACK TO THE YIELD SURFACE ***
!
  250 CONTINUE
      VNUI0 = VNU !AR(INTGLO)

      VNUI1 = 1 - 2*VNUI0
      VNUI2 = 1 + 2*VNUI0
      VNUIQ = VNUI2 / VNUI1
!
      DUM = SPSI / VNUI1
      PSIMIN = G *(-1+DUM)
      PSIMET = G *( 1+DUM)
      PSINU = 2*G *VNUI0*DUM
!
      HA= ( 1 - SPHI + SPSI - SPHI * SPSI) * SIG1 + 
     &      (-2 - 2/VNUI1*SPHI*SPSI) * SIG2 + 
     &     ( 1 - SPHI - SPSI + VNUIQ * SPHI * SPSI) * SIG3 + 
     &     2*(1 + SPSI) * COHS

      HB= ( 1 + SPHI + SPSI + VNUIQ * SPHI * SPSI) * SIG1 -
     &     ( 2 + 2/VNUI1*SPHI*SPSI) * SIG2 + 
     &     ( 1 + SPHI - SPSI - SPHI * SPSI) * SIG3 - 
     &     2*(1 - SPSI) * COHS

      HAB= (VNUI1+SPSI)*SIG1+ 
     &      (VNUI1-SPSI)*(SIG3-TENS)- 
     &      (VNUI1+SPSI)*(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
     
      HBA=(1-VNUI0)*SIG1-VNUI0*SIG3-(1-VNUI0)*(TENS*(1+SPHI)- 
     &      2*COHS)/(1-SPHI)+VNUI0*TENS

      HAO=(1-VNUI0)*SIG2-VNUI0*SIG3-VNUI1*TENS

      HOC=(1-VNUI0)*SIG1-VNUI0*SIG3-VNUI1*TENS

      HAA=(VNUI1+VNUI2*SPSI)*(SIG1-(TENS*(1+SPHI)-2*COHS)/ 
     &     (1-SPHI))+ 
     &     (VNUI1-SPSI)*(SIG2-TENS)+(VNUI1-SPSI)*(SIG3-TENS)

      HAAB=VNUI0*(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))- 
     &      (SIG2-TENS)+VNUI0*(SIG3-TENS)

      HAAO=(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))- 
     &      VNUI0*(SIG2-TENS)-VNUI0*(SIG3-TENS)

      HBB=(VNUI1+SPSI)*(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))+ 
     &     (VNUI1+SPSI)*(SIG2-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))+ 
     &     (VNUI1-VNUI2*SPSI)*(SIG3-TENS)

!        ***** DETERMINING RETURN AREA ****

      IAREA=0
      IASIGN=0

      IF (F31 > 0 .AND. HA >= 0 .AND. HB < 0 .AND. HAB < 0) THEN
        IAREA=2
        IASIGN=IASIGN+1
      ENDIF

      IF (F31 > 0 .AND. HB >= 0) THEN
        IASIGN=IASIGN+1
        IF (HBB < 0) THEN
          IAREA=1
          GOTO 260
        ELSE
          IAREA=8
        ENDIF
      ENDIF

      IF (F31 > 0 .AND. HA < 0) THEN
        IASIGN=IASIGN+1
        IF (HAA < 0) THEN

          IAREA=3
        ELSE
          IAREA=7
        ENDIF
      ENDIF

      IF (FT3 > 1d-12 .AND. 
     &     HBA >= 0 .AND. 
     &     HAO < 0 .AND. 
     &     HOC < 0) THEN
        IAREA=4
        IASIGN=IASIGN+1
      ENDIF

      IF (FT3 > 1d-12 .AND. HAB >= 0 .AND. HBA < 0) THEN
        IASIGN=IASIGN+1
        IF (HAAB > 0) THEN
          IAREA=5
        ELSE
          IF (IAREA == 0) IAREA=7
        ENDIF
      ENDIF

      IF (FT3 > 1d-12 .AND. HAO >= 0) THEN
        IASIGN=IASIGN+1
        IF (HAAO >= 0) THEN
          IAREA=6
        ELSE
          IF (IAREA == 0) IAREA=7
        ENDIF
      ENDIF

260   CONTINUE

      IF (IAREA == 0) THEN
        !call GiveWarning('NO ACTIVE RETURNING AREA! INTGLO=' // trim(String(INTGLO)))
      ENDIF


      DSP1=0
      DSP2=0
      DSP3=0

!
!         *** EXTENSION POINT ***
!
      IF (IAREA == 1) THEN
        A11 =       G *(1+SPHI*SPSI/VNUI1)
        A12 = 0.5d0*G *(1+SPHI+SPSI+VNUIQ*SPHI*SPSI)
        DETER  = A11*A11-A12*A12
        RLAM31 = (F31*A11-F32*A12) / DETER
        RLAM32 = (F32*A11-F31*A12) / DETER
!        IF (RLAM31 < 0) call GiveWarning('RLAM31 < 0  INTGLO' // trim(String(INTGLO)))
!        IF (RLAM32 < 0) call GiveWarning('RLAM32 < 0  INTGLO' // trim(String(INTGLO)))
        RLAM21 = 0
        DSP1 = RLAM31*PSIMIN + RLAM32*PSINU
        DSP2 = RLAM31*PSINU  + RLAM32*PSIMIN
        DSP3 = RLAM31*PSIMET + RLAM32*PSIMET
      ENDIF
!
!
!         *** REGULAR YIELD SURFACE ***
!
      IF (IAREA == 2) THEN
        A11 = G *(1+SPHI*SPSI/VNUI1)
        RLAM31 = F31 / A11
!        IF (RLAM31 < 0) call GiveWarning('RLAM31 < 0  INTGLO' // trim(String(INTGLO)))
        RLAM21 = 0
        RLAM32 = 0
        DSP1 = RLAM31*PSIMIN
        DSP2 = RLAM31*PSINU
        DSP3 = RLAM31*PSIMET
      END IF
!
!         *** COMPRESSION POINT ***
!
      IF (IAREA == 3) THEN
        A11 =       G *(1+SPHI*SPSI/VNUI1)
        A12 = 0.5d0*G *(1-SPHI-SPSI+VNUIQ*SPHI*SPSI)
        DETER  = A11*A11-A12*A12
        RLAM31 = (F31*A11-F21*A12) / DETER
        RLAM21 = (F21*A11-F31*A12) / DETER

!        IF (RLAM31 < 0) call GiveWarning('RLAM31 < 0  INTGLO' // trim(String(INTGLO)))
!        IF (RLAM21 < 0) call GiveWarning('RLAM21 < 0  INTGLO' // trim(String(INTGLO)))
        RLAM32 = 0
        DSP1 = RLAM31*PSIMIN + RLAM21*PSIMIN
        DSP2 = RLAM31*PSINU  + RLAM21*PSIMET
        DSP3 = RLAM31*PSIMET + RLAM21*PSINU
      END IF

      IF (IAREA == 4) THEN
        DUM=   2*G /VNUI1
        RLAMT3=FT3/(DUM*(1-VNUI0))
        IF (RLAMT3 < 0) Then
!          call WriteInLogFile('RLAMT3 < 0  INTGLO' // trim(String(INTGLO)))
!          call WriteInLogFile('FT3,RL3 = ' // trim(String(Ft3)) //' '// trim(String(RLamT3)))
!          call WriteInLogFile('G,vnui0 = ' // trim(String(G)) //' '// trim(String(vnui0)))
        End If
        DSP1 = DUM*RLAMT3*VNUI0
        DSP2 = DSP1
        DSP3 = DUM*RLAMT3*(1-VNUI0)
        ITENS=1
        IPL =2
      ENDIF

      IF (IAREA == 5) THEN
        A11 = VNUI1+SPHI*SPSI
        A12 = VNUI1+SPHI
        A21 = VNUI1+SPSI
        A22 = 2*(1-VNUI0)
        DUM = G /VNUI1
        DETER  = A11*A22-A12*A21

        RLAM31 = (F31*A22-FT3*A12) / DETER / DUM
        RLAMT3 = (FT3*A11-F31*A21) / DETER / DUM

!        IF (RLAM31 < 0) call WriteInLogFile('RLAM31 < 0  INTGLO' // trim(String(INTGLO)))
        IF (RLAMT3 < 0) Then
!          call WriteInLogFile('RLAMT3 < 0  (IA=5) INTGLO' // trim(String(INTGLO)))
!          call WriteInLogFile('F31 :' // trim(String(F31)))
!          call WriteInLogFile('FT3 :' // trim(String(FT3)))
!          call WriteInLogFile('A11 :' // trim(String(A11)))
!          call WriteInLogFile('A12 :' // trim(String(A12)))
!          call WriteInLogFile('A21 :' // trim(String(A21)))
!          call WriteInLogFile('A22 :' // trim(String(A22)))
!          call WriteInLogFile('HAB :' // trim(String(HAB)))
!          call WriteInLogFile('HBA :' // trim(String(HBA)))
        End If

        DSP1 = DUM*(RLAM31*(-VNUI1+SPSI)+RLAMT3*2*VNUI0)
        DSP2 = DUM*(RLAM31*2*VNUI0*SPSI+RLAMT3*2*VNUI0)
        DSP3 = DUM*(RLAM31*(VNUI1+SPSI)+RLAMT3*2*(1-VNUI0))
        ITENS=1
        IPL =2
      ENDIF

      IF (IAREA == 6) THEN
        RLAMT2 = (FT2*(1-VNUI0)-FT3*VNUI0)/(2*G)
        RLAMT3 = (FT3*(1-VNUI0)-FT2*VNUI0)/(2*G)

!        IF (RLAMT2 < 0) call WriteInLogFile('RLAMT2 < 0  INTGLO' // trim(String(INTGLO)))
        IF (RLAMT3 < 0) Then
!          call WriteInLogFile('RLAMT3 < 0  (IA=6) INTGLO' // trim(String(INTGLO)))
!         call WriteInLogFile('FT2 :' // trim(String(FT2)))
!          call WriteInLogFile('FT3 :' // trim(String(FT3)))
!          call WriteInLogFile('RL2 :' // trim(String(RL2)))
!          call WriteInLogFile('RL3 :' // trim(String(RL3)))
!          call WriteInLogFile('A21 :' // trim(String(A21)))
!          call WriteInLogFile('A22 :' // trim(String(A22)))
!          call WriteInLogFile('HAo :' // trim(String(HAo)))
!          call WriteInLogFile('HAAo :' // trim(String(HAAo)))
        End If

        DUM = 2*G/VNUI1
        DSP1 = DUM*VNUI0*(RLAMT2+RLAMT3)
        DSP2 = DUM*(RLAMT2*(1-VNUI0)+RLAMT3*   VNUI0 )
        DSP3 = DUM*(RLAMT2*   VNUI0 +RLAMT3*(1-VNUI0))

        ITENS=1
        IPL =2
      ENDIF

      IF (IAREA == 7) THEN
        DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
        DSP2 = SIG2-TENS
        DSP3 = SIG3-TENS
        ITENS=1
        IPL =2
      ENDIF

      IF (IAREA == 8) THEN
        DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
        DSP2 = SIG2-(SIG1-DSP1)
        DSP3 = SIG3-TENS
        ITENS=1
        IPL=2
      ENDIF

      IF (IAREA == 9) THEN
        DSP1 = SIG1-TENS
        DSP2 = SIG2-TENS
        DSP3 = SIG3-TENS
        IAPEX=1
        IPL=2
      ENDIF

!--------Check if the calculated principle stresses are on the yield surface
      SIG1R=SIG1-DSP1
      SIG2R=SIG2-DSP2
      SIG3R=SIG3-DSP3

      F21R = 0.5d0*(SIG2R-SIG1R) + 0.5d0*(SIG2R+SIG1R)*SPHI - COHS
      F32R = 0.5d0*(SIG3R-SIG2R) + 0.5d0*(SIG3R+SIG2R)*SPHI - COHS
      F31R = 0.5d0*(SIG3R-SIG1R) + 0.5d0*(SIG3R+SIG1R)*SPHI - COHS
      FT1R = SIG1R - TENS
      FT2R = SIG2R - TENS
      FT3R = SIG3R - TENS

      IF ( F21R > 1d-6.OR. 
     &      F32R > 1d-6.OR. 
     &      F31R > 1d-6.OR. 
     &      FT1R > 1d-6.OR. 
     &      FT2R > 1d-6.OR. 
     &      FT3R > 1d-6    ) THEN
        ICREC=0

        IF (IAREA == 5.AND.F32R > 1E-6) THEN
!------------ZONE 8
          ICREC=ICREC+1
          DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
          DSP2 = SIG2-(SIG1-DSP1)
          DSP3 = SIG3-TENS
        ENDIF

        IF (IAREA == 6.AND.FT1R > 1d-6) THEN
!------------ZONE 9
          ICREC=ICREC+1
          DSP1 = SIG1-TENS
          DSP2 = SIG2-TENS
          DSP3 = SIG3-TENS
          ITENS=0
          IAPEX=1
      ENDIF
!        IF (ICREC == 0) call GiveWarning('No correction is done!')
      ENDIF

!--------Check again if the calculated stresses are on the yield surface
      SIG1R=SIG1-DSP1
      SIG2R=SIG2-DSP2
      SIG3R=SIG3-DSP3

      F21 = 0.5d0*(SIG2R-SIG1R) + 0.5d0*(SIG2R+SIG1R)*SPHI - COHS
      F32 = 0.5d0*(SIG3R-SIG2R) + 0.5d0*(SIG3R+SIG2R)*SPHI - COHS
      F31 = 0.5d0*(SIG3R-SIG1R) + 0.5d0*(SIG3R+SIG1R)*SPHI - COHS
      FT1 = SIG1R - TENS
      FT2 = SIG2R - TENS
      FT3 = SIG3R - TENS

!      IF (F21 > 0.001) call GiveWarning('F21=' // trim(String(F21))//' STILL >0')
!      IF (F32 > 0.001) call GiveWarning('F32=' // trim(String(F32))//' STILL >0')
!      IF (F31 > 0.001) call GiveWarning('F31=' // trim(String(F31))//' STILL >0')
!      IF (FT1 > 0.001) call GiveWarning('FT1=' // trim(String(FT1))//' STILL >0')
!      IF (FT2 > 0.001) call GiveWarning('FT2=' // trim(String(FT2))//' STILL >0')
!      IF (FT3 > 0.001) call GiveWarning('FT3=' // trim(String(FT3))//' STILL >0')

      IF (SIG2R > SIG3R+0.01 .OR. SIG2R < SIG1R-0.01) THEN
!        call GiveWarning('ASSUMPTION NOT CORRECT!')
!        call GiveWarning('The original returning zone IAREA=' // trim(String(IAREA)))
      ENDIF

!
!         *** Computing Cartesian stress components ***
!
!      If (IntGlo == ip) then
!        call WriteInLogFile('SigR123 : ' // trim(String(Sig1R)) //' '//trim(String(Sig2R)) //' '//trim(String(Sig3R)))
!      endif

  903       format(a,6(1x,1p,e19.12),1x,i4)

      Call CarSig(Sig1R,Sig2R,Sig3R,xN1,xN2,xN3,NTENS,SigC)
      STRESS=SigC



!------- Calculate TAUMAX for checking accuracy on plastic point

!      If (IntGlo == ip) then
!        call WriteInLogFile('Sig123 : ' // trim(String(Sig1)) //' '//trim(String(Sig2)) //' '//trim(String(Sig3)))
!      endif

      TAUMAX=0.5d0*(SIG3-DSP3-SIG1+DSP1)
      taumax1=taumax
!     IF (TAUMAX <= -1E-6) call WriteInLogFile(' TAUMAX IS NEGATIVE!!!')
      TAUMAX = MAX(TAUMAX,COHS,0.5d0)
!      If (IntGlo == ip) then
!        call WriteInLogFile('DSig123 : ' // trim(String(DSP1)) //' '//trim(String(DSP2)) //' '//trim(String(DSP3)))
!        call WriteInLogFile('TauMax1 : ' // trim(String(TauMax1)))
!        call WriteInLogFile('Taumax  : ' // trim(String(Taumax)))
!      endif

360   CONTINUE

      RETURN

      End   
!-----------------------------------------------------------------------------     
      Subroutine MatTranspose (A,Ia,AT,IAt,N1,N2)
      implicit none
      integer :: Ia, N1, N2, IAt ! N1, N2 are the dimensions
      double precision :: A(Ia,*),AT(Iat,*)
      integer :: I, J

      Do I=1,N1
        Do J=1,N2
          AT(I,J)=A(J,I)
        End Do
      End Do

      End 