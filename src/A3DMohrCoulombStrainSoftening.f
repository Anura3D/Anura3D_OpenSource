
    module ModMohrCoulombStrainSoftening
    contains 
        
        SUBROUTINE ESM_MohrCoulombStrainSoftening(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
     DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)

      !implicit real(REAL_TYPE) (a-h, o-z) 
        implicit double precision (a-h, o-z) 
      integer :: NTENS, NSTATEV, NADDVAR, NPROPS, NPT, NOEL, IDSET, NUMBEROFPHASES
      double precision :: EUNLOADING, PLASTICMULTIPLIER
      CHARACTER*80 CMNAME         
     ! DIMENSION NPT(1),NOEL(1),IDSET(1),STRESS(NTENS),EUNLOADING(1),PLASTICMULTIPLIER(1),&
     !DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS),NUMBEROFPHASES(1)
       DIMENSION STRESS(NTENS),DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)

!---Local variables required in standard UMAT
        integer :: IStep, TimeStep
        double precision, dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
        double precision, dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
        double precision, dimension(:), allocatable :: stran
        double precision, dimension(:), allocatable :: time
        double precision, dimension(:), allocatable :: predef
        double precision, dimension(:), allocatable :: dpred    
        double precision, dimension(:), allocatable :: coords
        double precision, dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
        double precision, dimension(:,:), allocatable :: drot
        double precision, dimension(:,:), allocatable :: dfgrd0
        double precision, dimension(:,:), allocatable :: dfgrd1
        double precision :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
        double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
        double precision :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
        double precision :: pnewdt, dtime, temp, dtemp, celent
        double precision :: Value ! auxiliary variable holding any real valued number
        double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation  
        double precision :: VolumetricStress
        integer :: ITens
    
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

        
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
    
! Initialization
        Eunloading = 0.0
        PlasticMultiplier = 0.0
        !IDTask = 2
     
!Rename additional variables
        Porosity = AdditionalVar(1)
        WaterPressure = AdditionalVar(2)
        WaterPressure0 = AdditionalVar(3)
        GasPressure = AdditionalVar(4)
        GasPressure0 = AdditionalVar(5)
        DegreeSaturation = AdditionalVar(6)
        time(1) = AdditionalVar(7)   !TotalRealTime
        time(2) = AdditionalVar(8)   !OverallTotalTime
        dtime = AdditionalVar(9)     !TimeIncrement
        IStep = AdditionalVar(10)    
        TimeStep = AdditionalVar(11)   !Note: Very first time and load step: Istep=1 and TimeStep=1   
        
      IF((IStep==1).and.(TimeStep==1)) then !IDTask = 1
     
      !IF (IDTask == 1) then ! initialisation of state variables
       STATEV(1)=PROPS(3)
       STATEV(2)=PROPS(5)
       STATEV(3)=PROPS(7)
      END IF ! IDTask = 1
      
!---Call the UMAT      
          call umat_AM(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)
          
!---Evaluate mean effective stress           
        !VolumetricStress = stress(1) + stress(2) + stress(3)
        !if ((VolumetricStress>0.0) ) then!.or. (stress(1)>0.0) .or. (stress(2)>0.0) .or. (stress(3)>0.0)) then 
        !    do ITens = 1, NTENS 
        !        stress(ITens) = 0.0
        !    end do 
        !end if 
        
            
        

      
!---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
!---Always define this value to run the simulation

    ! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
    
          


        return

     end subroutine ESM_MohrCoulombStrainSoftening
        
        
        
        
        
        !*USER SUBROUTINES
      SUBROUTINE UMAT_AM(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
      RPL,DDSDDT,DRPLDE,DRPLDT,&
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      !INCLUDE 'ABA_PARAM.INC'
!
      implicit double precision (a-h, o-z) 
      CHARACTER*80 CMNAME
      !DIMENSION STRESS(NTENS),STATEV(NSTATEV),&
      !DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),&
      !STRAN(NTENS),DSTRAN(NTENS),TIME(2),&
      !PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      DIMENSION STRESS(NTENS),STATEV(NSTATEV),&
      DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),&
      STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
      PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      

! Arguments:
!          I/O  Type
!  PROPS    I   R()  : List with model parameters
!  DSTRAN   I   R()  : Strain increment
!  DDSDDE   O   R(,) : Material stiffness matrix
!  STRESS  I/O  R()  : stresses
!  STATEV  I/O  R()  : state variables
!
!
!---  Local variables
!
      Dimension :: DE(6,6), dSig(6), &
                  Sig(6), dEpsP(6), EpsP(6)

!
! Mohr-Coulomb Strain Softening model 
!
! Contents of PROPS(9) MCSS
!  1 : G       shear modulus
!  2 : ENU     Poisson's ratio
!  3 : cp      peak cohesion
!  4 : cr      residual cohesion
!  5 : phip    peak friction angle 
!  6 : phir    residual friction angle
!  7 : psip    peak dilation angle
!  8 : psir    residual dilation angle 
!  9 : factor  shape factor
!
        Rad  = 45d0 / datan(1d0)
!*
!* ... start correction routine
!*
        G      = PROPS(1)         ! shear modulus
        ENU    = PROPS(2)         ! Poisson's ratio
        cp     = PROPS(3)         ! peak cohesion
        cr     = PROPS(4)         ! residual cohesion
        phip   = PROPS(5)/Rad     ! peak friction angle (rad)
        phir   = PROPS(6)/Rad     ! residual friction angle (rad)
        psip   = PROPS(7)/Rad     ! peak dilation angle (rad)
        psir   = PROPS(8)/Rad     ! residual dilation angle (rad)
        factor = PROPS(9)         ! shape factor
        
        c    = STATEV(1)          ! cohesion 
        phi  = STATEV(2)          ! friction angle
        psi  = STATEV(3)          ! dilatancy angle
        Do i = 1,NTENS
        EpsP(i) = STATEV(3+i) 
        end do

        ipl     =   0
!*
        ! Fill elastic material matrix
        F1  = 2*G*(1-ENU)/(1-2*ENU)
        F2  = 2*G*( ENU )/(1-2*ENU)
        DE  = 0.0
        DE(1:3,1:3) = F2
        DE(1,1) = F1
        DE(2,2) = F1
        DE(3,3) = F1
        DE(4,4) = G
        DE(5,5) = G
        DE(6,6) = G
!*
        ! elastic stress increment
        Call MatVec( DE, 6, DSTRAN, 6, dSig)
        ! elastic stress
        Call AddVec( STRESS, dSig, 1d0, 1d0, 6, Sig )
        
        call MOHRStrainSoftening(IntGlo,F1,F2,G,cp,cr,phip,phir,psip,psir,factor,c,phi,psi,stress,dSig,EpsP,DSTRAN,dEpsP,Sig,IPL)

!*
!* ... stress state parameters update
!*
        Do i=1,NTENS
          STRESS(i) = Sig(i)
        End Do
        
        STATEV(1) = c   
        STATEV(2) = phi  
        STATEV(3) = psi   
        Do i = 1,NTENS
        STATEV(3+i) = EpsP(i)
        end do

!*
!* ... Tangent stiffness matrix to be returned (done by elastic stiffness)
!*
        G       =   PROPS(1)       ! G
        ENU     =   PROPS(2)       ! nu
        F1  = 2*G*(1-ENU)/(1-2*ENU)
        F2  = 2*G*( ENU )/(1-2*ENU)
        DDSDDE = 0.0
        DDSDDE(1:3,1:3) = F2
        DDSDDE(1,1) = F1
        DDSDDE(2,2) = F1
        DDSDDE(3,3) = F1
        DDSDDE(4,4) = G
        DDSDDE(5,5) = G
        DDSDDE(6,6) = G
!*
!* ... end UMAT routine
!*
      Return
      End  

!***********************************************************************
      Subroutine MOHRStrainSoftening(IntGlo,D1,D2,GG,cp,cr,phip,phir, &
       psip,psir,factor,c,phi,psi,Sig0,DSigE,EpsP,DEps,DEpsP,SigC,IPL)
      !**********************************************************************
      !
      ! Elastoplastic constitutive model with STRAIN SOFTENING, based on the 
      ! MOHR-COULOMB criterion (considering modifications of Abbo & Sloan (1995))
      ! Explicit MODIFIED EULER INTEGRATION SCHEME with automatic error control.
      ! Final correction of the yield surface drift (END OF STEP CORRECTION).
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i,n,m,it
      double precision :: F,F0,F2 !Evaluation of the Yield function
      double precision :: alpha !Elastic Strain proportion
      double precision :: SSTOL !Tolerance Relative Error
      double precision :: YTOL !Tolerance Relative Error of the yield function evaluation
      double precision :: SPTOL !Tolerance Softening parameters
      double precision :: Rn !Relative error
      double precision :: T,DT,T1,beta,DTmin !Substepping parameters
      double precision :: c1,phi1,psi1,c2,phi2,psi2
      double precision :: ctol,phitol,psitol !c,phi,psi tolerances
      double precision :: Dcr,Dphir,Dpsir !Diference between current and residial values
      double precision :: moduleEr,moduleSigDSig
      double precision :: EpsPEq,EpsPEq1,EpsPEq2 !Equivalent Plastic Deformation
      double precision :: DEpsPEq !Derivative Strain in function of Equivalent Plastic Deformation
      double precision :: p,J,Lode,S3TA
      double precision, dimension(3) :: Principal_stresses
      double precision, dimension(3) :: SigC_TC
      double precision, dimension(3,3) :: SigC_TC_mat
      double precision, dimension(3,3) :: SigC_cartesian
      integer :: ii, jj, kk, ll
      double precision, dimension(3,3) :: Principal_vectors
      double precision, dimension(6) :: SigYield, SigYield2
      double precision, dimension(6) :: DSigPP,DSigP1,DSigP2
      double precision, dimension(6) :: DEpsPP,DEpsPP1,DEpsPP2
      double precision, dimension(6) :: DEpsS,DEpsSS
      double precision, dimension(6) :: EpsP1,EpsP2
      double precision, dimension(6) :: DEpsPEqDPS,DEpsPEqDPS1
      double precision, dimension(6) :: sumSg,Er
      double precision, dimension(3) :: DSPDPEq,DSPDPEq1 !Variation of softening parameters (c,phi,psi) in function of plastic strain
      !In variables
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      double precision, intent(in) :: D1,D2,GG !Elastic Parameters
      double precision, intent(in) :: cp,cr,phip,phir,psip,psir,factor !Softening parameter
      double precision, intent(in), dimension(6) :: Sig0 !Initial Stress
      double precision, intent(in), dimension(6) :: DEps !Incremental total strain
      !Inout variables
      double precision, intent(inout):: c,phi,psi !cohesion,friction angle and dilatancy angle
      double precision, intent(inout), dimension(6) :: EpsP !Accumulated Plastic Strain
      double precision, intent(inout), dimension(6) :: SigC !Final Stress
      double precision, intent(inout), dimension(6) :: DSigE !Incremental Elastic Stress
      !Out variables
      integer, intent(out) :: IPL
      double precision, intent(out), dimension(6) :: DEpsP !Incremental plastic strain

      !Initialization
      DEpsPEq = 0.0d0
      EpsPEq = 0.0d0
      SigYield = 0.0d0
      DEpsP = 0.0d0
      F = 0.0d0
      it = 0

      if (c > cp.or.phi > phip.or.psi > psip) then
        c = min(c,cp)
        phi = min(phi,phip)
        psi = min(psi,psip)
      end if
      if (c < cr.or.phi < phir.or.psi < psir) then
        c = max(c,cr)
        phi = max(phi,phir)
        psi = max(psi,psir)
      end if

      !Tolerances
      SSTOL = 0.001d0 !Tolerance Relative Error (10-3 to 10-5)
      YTOL = 0.000001d0 !Tolerance Error on the Yield surface (10-6 to 10-9)
      SPTOL = 0.0001d0 !Tolerance Softening Parameters (0.0001d0)
      ctol = abs(cp-cr)*SPTOL
      phitol = abs(phip-phir)*SPTOL
      psitol = abs(psip-psir)*SPTOL
      DTmin = 0.000000001d0
      
      
      
      !call CalculateInvariants(IntGlo,SigC,p,J,Lode,S3TA)
      !if (p > 0) then !hardcoding T as zero
      !    call Get_EigenValues_EigenVectors(SigC, Principal_stresses, Principal_vectors)          
      !    !SigC = 0.0
      !    SigC_TC(3) = Principal_stresses(3) + ( - (Principal_stresses(1)+Principal_stresses(2)+Principal_stresses(3)) )/3
      !    SigC_TC(2) = Principal_stresses(2) + ( - (Principal_stresses(1)+Principal_stresses(2)+Principal_stresses(3)) )/3
      !    SigC_TC(1) = Principal_stresses(1) + ( - (Principal_stresses(1)+Principal_stresses(2)+Principal_stresses(3)) )/3
      !    SigC = [Sig_TC(1), Sig_TC(2), Sig_TC(3), 0, 0, 0]
      !end if
      ! START TENSION CUTOFF ---------------------------------------------------------------------------------------------------------------
      !call CalculateInvariants(IntGlo, SigC, p, J, Lode, S3TA)
      !
      !
      !if (p > 0) then
      !
      !    ! Compute principal stresses and directions
      !    call Get_EigenValues_EigenVectors(SigC, Principal_stresses, Principal_vectors)
      !
      !
      !    ! Adjust principal stresses for tension cutoff or other modifications
      !    SigC_TC(3) = Principal_stresses(3) - (Principal_stresses(1) + Principal_stresses(2) + Principal_stresses(3)) / 3.0
      !    SigC_TC(2) = Principal_stresses(2) - (Principal_stresses(1) + Principal_stresses(2) + Principal_stresses(3)) / 3.0
      !    SigC_TC(1) = Principal_stresses(1) - (Principal_stresses(1) + Principal_stresses(2) + Principal_stresses(3)) / 3.0
      !
      !
      !    ! Reassemble SigC_TC as a diagonal tensor
      !    do ii = 1, 3
      !        do jj = 1, 3
      !            SigC_TC_mat(ii, jj) = 0.0
      !        end do
      !        SigC_TC_mat(ii, ii) = SigC_TC(ii)
      !    end do
      !
      !
      !    ! Perform rotation back to Cartesian coordinates: SigC = R * SigC_TC * R^T
      !    do ii = 1, 3
      !        do jj = 1, 3
      !            SigC_cartesian(ii, jj) = 0.0
      !            do kk = 1, 3
      !                do ll = 1, 3
      !                    SigC_cartesian(ii, jj) = SigC_cartesian(ii, jj) + Principal_vectors(ii, kk) * SigC_TC_mat(kk, ll) * Principal_vectors(jj, ll)
      !                end do
      !            end do
      !        end do
      !    end do
      !    
      !    
      !    
      !    ! Assemble the rotated stress tensor in Voigt notation
      !    SigC(1) = SigC_cartesian(1, 1) ! ?_xx
      !    SigC(2) = SigC_cartesian(2, 2) ! ?_yy
      !    SigC(3) = SigC_cartesian(3, 3) ! ?_zz
      !    SigC(4) = SigC_cartesian(2, 3) ! ?_yz
      !    SigC(5) = SigC_cartesian(1, 3) ! ?_xz
      !    SigC(6) = SigC_cartesian(1, 2) ! ?_xy
      !    
      !    
      !end if
      
      ! END TENSION CUTOFF ---------------------------------------------------------------------------------------------------------------
      
      
      
      
      !Check the yield function value
      call DetermineYieldFunctionValue(IntGlo,SigC,c,phi,F)
      
      
      !If F<0 then the behaviour is elastic --> Return
      if (F <= YTOL) then
        IPL = 0
        ! Tension cutoff implementation 
      ! First -> compute mean effective stress 'p' and deviatoric stress 'J' for SigC
      !call CalculateInvariants(IntGlo,SigC,p,J,Lode,S3TA)
      
        return
      end if

      !If F>0, the behaviour is elastoplastic --> Continue
      Dcr = abs(c - cr)
      Dphir = abs(phi - phir)
      Dpsir = abs(psi - psir)
      !Check if we are in residual conditions or in softening conditions
      if (Dcr <= ctol.and.Dphir <= phitol.and.Dpsir <= psitol) then
        IPL = 1 !IPL=1 Residual conditions --> no changes of the strength parameters
        c = cr
        phi = phir
        psi = psir
      else
        IPL = 2 !IPL=2 Softening Conditions --> changes of the strength parameters
      end if

      !Determine the proportion (alpha) of the stress increment that lies within the yield function.
      !The PEGASUS ALGORITHM SCHEME FOR CONVENTIONAL ELASTOPLASTIC MODELS has been used
      call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)

      if (F0 < -YTOL) then !In this Time increment there is part of elastic behavior
        call DetermineElasticProportionPegasusMethod(IntGlo,Sig0,DSigE,DEps,c,phi,YTOL,alpha)
      else  
        alpha = 0.0d0 !In this time increment all is plastic behavior
      end if

      !Calculate the direction of the stress path--> missing
      !It is assumed that the direction is always outside the yield surface.

      !Determine the elastic portion of the stress increment
      DSigE = alpha * DSigE !Correct Incremental Elastic Stress

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Determine the plastic portion of the stress increment.
      !The method used is the MODIFIED EULER INTEGRATION SCHEME with error control
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Initialise parameters
      SigYield = Sig0 + DSigE !Sigma on the Yield surface
      DEpsS = (1.0d0-alpha) * DEps !Incremental Plastic Strain

      T = 0.0d0
      DT = 1.0d0

      !Start the plastification
      Do while (T <= 1.0d0)
        m = 0 !Counter
        Rn = 100

        call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)

        Do while (Rn > SSTOL.and.m < 1000)
            !1)Calculation of the portion of the plastic strain increment (DEpsPP)
            DEpsSS = DT * DEpsS !Portion of the plastic strain increment

            !Calculate a first estimate of the associated stress 
            !hardening/softening parameter changes
            call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
                EpsPEq,DSPDPEq)
            call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
            call DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c,phi,psi,SigYield,DEpsPEqDPS,DSPDPEq,DEpsSS,DSigP1,DEpsPP1)
            EpsP1 = EpsP + DEpsPP1

            call CalculateEpsPEq(EpsP1,EpsPEq1) !Determine Equivalent plastic Strain (EpsPEq)

            !if (IPL == 1) then !Residual conditions --> no changes of the strength parameters
            !    c1 = c
            !    phi1 = phi
            !    psi1 = psi
            !else !IPL=2 Softening Conditions --> changes of the strength parameters
                call CalculateSofteningParameters(EpsPEq1,factor,cp,cr,phip,phir,psip,psir,c1,phi1,psi1)
            !end if

            !2)Calculate a second estimate of the associated stress 
            !hardening/softening parameter changes
            SigYield2 = SigYield + DSigP1

            call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
                EpsPEq1,DSPDPEq1)
            call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP1,EpsPEq1,DEpsPEqDPS1)
            call DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c1,phi1,psi1,SigYield2,DEpsPEqDPS1,DSPDPEq1,DEpsSS,DSigP2,DEpsPP2)
            EpsP2 = EpsP + DEpsPP2

            call CalculateEpsPEq(EpsP2,EpsPEq2) !Determine Equivalent plastic Strain (EpsPEq)

            !if (IPL == 1) then !Residual conditions --> no changes of the strength parameters
            !    c2 = c
            !    phi2 = phi
            !    psi2 = psi
            !else  !IPL=2 Softening Conditions --> changes of the strength parameters
                call CalculateSofteningParameters(EpsPEq2,factor,cp,cr,phip,phir,psip,psir,c2,phi2,psi2)
            !end if

            !3)Obtain a more accurate modified Euler estimate of the changes in stress,
            !plastic strain and hardening/softening parameters
            DSigPP = 0.5d0 * (DSigP1 + DSigP2)

            !Calculation of the relative error
            Er = 0.5d0 * (DSigP1 - DSigP2)
            moduleEr = sqrt(Er(1)*Er(1)+Er(2)*Er(2)+ Er(3)*Er(3)+ Er(4)*Er(4)+Er(5)*Er(5)+Er(6)*Er(6))

            sumSg = SigYield + DSigPP
            moduleSigDSig = sqrt(sumSg(1)*sumSg(1) + sumSg(2)*sumSg(2) + sumSg(3)*sumSg(3)+ &
                                 sumSg(4)*sumSg(4) + sumSg(5)*sumSg(5) + sumSg(6)*sumSg(6))

            !Check the relative error (Rn) of the new stresses, with the defined tolerance (SSTOL)
            Rn = (moduleEr/moduleSigDSig)

            ! check whether decreasing of DT is possible, if not exit loop
            if (DT == DTmin) then
               exit
            end if
          
            !4)If Rn>SSTOL, the loop is not finished and the substep is recalculated smaller
            if (Rn > SSTOL) then
                beta = max (0.9d0*(sqrt(SSTOL/Rn)), 0.1d0)
                beta = min (beta, 1.1d0)
                DT = max (DT*beta, DTmin)
                m = m + 1 !Update counter
            end if

        end do

        !Update the accumulated stresses, plastic strain and softening parameters
        SigYield = SigYield + DSigPP
        DEpsPP = 0.5d0 * (DEpsPP1 + DEpsPP2)
        DEpsP = DEpsP + DEpsPP
        EpsP = EpsP + DEpsPP

        call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)

        call CalculateSofteningParameters(EpsPEq,factor,cp,cr,phip,phir,psip,psir,c,phi,psi)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!! END OF STEP CORRECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Check if we are on/in the yield surface, otherwise we are still outside (F>0)
        !and a correction is needed.
        call DetermineYieldFunctionValue(IntGlo,SigYield,c,phi,F)
        n=0 !Counter
        do while (abs(F) > YTOL.and.n < 20) !The correction is needed
            n = n + 1
            call CalculateEpsPEq(EpsP,EpsPEq)             !Determine Equivalent plastic Strain (EpsPEq)
            call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir, &
                EpsPEq,DSPDPEq)
            call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
            call EndOfStepCorrection(IntGlo,D1,D2,GG,IPL,F,SigYield,DSPDPEq,DEpsPEqDPS,EpsP,c,phi,psi)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !The substep is updated
        T1 = T + DT

        !If T1>1 the calculation is finished
        If (T1 >= 1d0) then
            SigC = SigYield   !Determine Final stresses
            
            ! Tension cutoff implementation 
            ! First -> compute mean effective stress 'p' and deviatoric stress 'J' for SigC
          !  call CalculateInvariants(IntGlo,SigC,p,J,Lode,S3TA)
          !  
            !if (p > 0) then !hardcoding T as zero
          !
          !      !call Get_EigenValues_EigenVectors(SigC, Principal_stresses, Principal_vectors)  
                !SigC = 0.0
          !!SigC_TC = 0.0
          !!SigC_TC(3) = SigC(3) + (0 - p)
          !!SigC_TC(2) = SigC(2) + (0 - p)
          !!SigC_TC(1) = SigC(1) + (0 - p)
          !      
           ! end if
            
            
            return
        end if

        !If T1<1, calculation of the next substep DT
        beta = min (0.9d0*(sqrt(SSTOL/Rn)), 1.1d0)
        if (m > 1) then ! the previous step failed
            beta = min (beta, 1.0d0)
            DT = beta * DT
            it = it+1
        else
            DT = beta * DT
            it = 0
        end if
        DT = max (DT, DTmin)
        DT = min (DT, 1.0d0-T1)
        T = T1
        
      end do  !If T=1 the loop is finished
      
      
      
      

      end subroutine MOHRStrainSoftening
 

      Subroutine DetermineElasticProportionPegasusMethod(IntGlo,Sig,DSig,DEps,c,phi,YTOL,alpha)
      !**********************************************************************
      !
      ! The PEGASUS METHOD method is used  
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: Its     
      double precision :: alpha0,alpha1,F0,F1,F
      double precision, dimension(6) :: Sig0,Sig1,SigNew
      !In variables
      double precision, intent(in), dimension(6) :: Sig, DSig
      double precision, intent(in), dimension(6) :: DEps
      double precision, intent(in) :: c,phi,YTOL
      integer, intent(in) :: IntGlo       !Global ID of Gauss point or particle
      !Out variables
      double precision, intent(out) :: alpha

      alpha0 = 0.0d0
      alpha1 = 1.0d0

      Sig0 = Sig + (alpha0*DSig) ! = Sig0
      Sig1 = Sig + (alpha1*DSig) ! = SigE
        
      call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)
      call DetermineYieldFunctionValue(IntGlo,Sig1,c,phi,F1)

      F=YTOL+1000
      Its = 0 !Counter

      do while (abs(F) > YTOL .and. Its < 1000)
        alpha = alpha1 - (F1*(alpha1-alpha0)/(F1-F0))
        SigNew = Sig + (alpha*DSig)
       
        call DetermineYieldFunctionValue(IntGlo,SigNew,c,phi,F)

        if ((F*F1) < 0.0d0) then
            alpha0 = alpha1
            F0 = F1
        else
            F0 = F1*F0/(F1+F)
        end if

        alpha1 = alpha
        F1 = F
        Its = Its + 1

      end do  
      if (Its >= 1000) then
            alpha = 0.0d0
      end if

      end subroutine DetermineElasticProportionPegasusMethod


      Subroutine CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      !**********************************************************************
      !
      ! Calcuation of the invariants (defined as Abbo & Sloan (1995))
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision :: Sx,Sy,Sz,SqTxy,SqTyz,SqTxz,suma,h1,h2,J2,J3
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00001 = 1.0D0
      double precision, parameter :: C00P16 = 0.166666666666666D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: CP3333 = 0.333333333333333D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0
      double precision, parameter :: TINY = 0.000000000000001D0
      !In variables
      double precision, intent(in), dimension(6) :: Sig
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      !Out variables
      double precision, intent(out) :: p,J,Lode,S3TA !Invariants

      p = C00000
      J = C00000
      Lode = C00000

      !Mean stress (p)
      p = CP3333 * (Sig(1) + Sig(2) + Sig(3))

      !Deviatoric stress (J)
      Sx = Sig(1) - p
      Sy = Sig(2) - p
      Sz = Sig(3) - p
      suma = (Sig(1)-Sig(2))*(Sig(1)-Sig(2))+(Sig(1)-Sig(3))*(Sig(1)-Sig(3))+(Sig(2)-Sig(3))*(Sig(2)-Sig(3))
      SqTxy =  Sig(4) * Sig(4)
      SqTyz =  Sig(5) * Sig(5)
      SqTxz =  Sig(6) * Sig(6)

      J2 = C00P16 * suma + SqTxy + SqTyz + SqTxz
      J3 = Sx*Sy*Sz + C00002 * Sig(4)*Sig(5)*Sig(6) - Sx*SqTyz - Sy*SqTxz - Sz*SqTxy
      J = SQRT(J2)

      !Lode's angle (Lode)
      if (J2 > 0.0) then !C00000

        h1 = -C00003/(C00002*C00IR3)
        h2 = J3/(J*J*J)
        S3TA = h1*h2
        if (S3TA < -C00001) then
            S3TA = -C00001
        else if (S3TA > C00001) then
            S3TA = C00001
      end if
        Lode = CP3333*asin(S3TA)
      else  !Special case of zero deviatoric stress
        Lode = C00000
        S3TA = C00000
      end if

      end subroutine CalculateInvariants

 
      Subroutine DetermineYieldFunctionValue(IntGlo,Sig,c,phi,F)
      !**********************************************************************
      !
      ! In this subroutine the yield function evaluated is a smooth hyperbolic approximation to the
      ! Mohr-Coulomb yield criterion (Abbo and Sloan, 1995).
      !
      ! The edges of the hexagonal pyramid and the tip have been smoothed.
      ! There are two parameters aSmooth (smoothes the tip) and ATTRAN(smoothes the edges)
      ! In this case aSmooth=0.0005*c*cot(phi) and LodeT=25º.
      ! If aSmooth=0 and LodeT=30º the original Mohr-Coulomb is obtained.
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision ::  p,J,Lode,S3TA !Invariants
      double precision ::  COH, SPHI, CPHI, COTPHI, STA, CTA, K, aSmooth, ASPHI2, SGN, A, B
      double precision, parameter :: C00001 = 1.0d0 !Parameters
      double precision, parameter :: C00003 = 3.0d0
      double precision, parameter :: C00P50 = 0.0005d0
      double precision, parameter :: C00000 = 0.0d0
      double precision, parameter :: C00IR3 = 0.577350269189626d0
      double precision, parameter :: C000P1 = 0.00000000001D0
      !Constants for rounded K function (for LodeT=25)
      !real(REAL_TYPE), parameter :: A1 = 1.432052062044227d0
      !real(REAL_TYPE), parameter :: A2 = 0.406941858374615d0
      !real(REAL_TYPE), parameter :: B1 = 0.544290524902313d0
      !real(REAL_TYPE), parameter :: B2 = 0.673903324498392d0
      !real(REAL_TYPE), parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !real(REAL_TYPE), parameter :: A1 = -138300705.446275
      !real(REAL_TYPE), parameter :: A2 = -138300706.472675
      !real(REAL_TYPE), parameter :: B1 = -138300706.3123
      !real(REAL_TYPE), parameter :: B2 = 0.192450089729875
      !real(REAL_TYPE), parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
      !In variables
      double precision, intent(in), dimension(6) :: Sig
      double precision, intent(in) :: c,phi
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle

      !Out variables
      double precision, intent(out) :: F

      F = C00000

      !Calculation of the invariants (p',J,Lode)
      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!Evaluation of the yield function with Smoothing!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Material parameters
      COH = c     !Cohesion
      SPHI = sin(phi) 
      CPHI = cos(phi)
      if (ABS(SPHI)>0.0) then 
      COTPHI = CPHI/SPHI
      else 
      COTPHI = 0.0    
      end if 
      
      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      if (abs(phi) == C00000) then
            ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      !Calculate K function
      if (abs(Lode) < ATTRAN) then
        STA = sin(Lode)
        CTA = cos(Lode)
        K = CTA - STA*SPHI*C00IR3
      else
        SGN = SIGN(C00001,Lode)
        A = A1 + A2*SGN*SPHI
        B = B1*SGN + B2*SPHI
        K = A - B*S3TA
      end if

      !Calculate value of Hyperbolic Yield function
      F = p*SPHI + sqrt(J*J*K*K+ASPHI2) - COH*CPHI
      
      end subroutine DetermineYieldFunctionValue


      Subroutine CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the yield function (F) and the plastic potencial punction (P).
      ! Based on Abbo & Sloan (1995)
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i
      double precision :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B,&
                           D, aSmooth, ASPHI2, SGN, T3TA, C3TA, J2, psi2
      double precision ::   K, dKdLode
      double precision :: SPSI, CPSI, TPSI, COTPSI, ASPSI2
      double precision :: i1, i2, Sx, Sy, Sz
      double precision :: DFDp,DFDJ,DFDLode !Derivatives F respect Invariants
      double precision :: DPDp,DPDJ,DPDLode !Derivatives P respect Invariants
      double precision :: C1, C2, C3
      double precision, dimension(6):: DpDSig,DJDSig,DJ3DSig !Derivatives Invariants

      double precision, parameter :: C00001 = 1.0D0 !Parameters
      double precision, parameter :: C000P5 = 0.5D0
      double precision, parameter :: C00P50 = 0.0005D0
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: C00004 = 4.0D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: CP3333 = 0.333333333333333D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0
      double precision, parameter :: C0R3I2 = 0.866025403784439D0
      double precision, parameter :: C000P1 = 0.000000000000001D0 
      double precision, parameter :: J0 = 0.001D0 
      !Constants for rounded K function (for LodeT=25)
      !real(REAL_TYPE), parameter :: A1 = 1.432052062044227d0
      !real(REAL_TYPE), parameter :: A2 = 0.406941858374615d0
      !real(REAL_TYPE), parameter :: B1 = 0.544290524902313d0
      !real(REAL_TYPE), parameter :: B2 = 0.673903324498392d0
      !real(REAL_TYPE), parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !real(REAL_TYPE), parameter :: A1 = -138300705.446275
      !real(REAL_TYPE), parameter :: A2 = -138300706.472675
      !real(REAL_TYPE), parameter :: B1 = -138300706.3123
      !real(REAL_TYPE), parameter :: B2 = 0.192450089729875
      !real(REAL_TYPE), parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
      !In variables
      double precision, intent(in) ::  c,phi,psi !Soft Parameters
      double precision, intent(in), dimension(6) :: Sig
      !Out variables
      double precision, intent(out), dimension(6) :: DFDSig, DPPDSig !Derivatives respect Sigma
      !Inout variables
      double precision, intent(inout) :: p,J,Lode,S3TA !Invariants

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! DFDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Material parameters
      COH = c !Cohesion
      SPHI = sin(phi) 
      CPHI = cos(phi)
      COTPHI = CPHI/SPHI
      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      if (abs(phi) == C00000) then
        ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      if (J == C00000) then
        J2 = C000P1
        J = sqrt(J2)
      else
        J2 = J*J
      end if

      CTA = cos(Lode)
      C3TA = CTA*(C00004*CTA*CTA-C00003)
      T3TA = S3TA/C3TA

      !Calculate K function and its derivative
      if (abs(Lode) < ATTRAN) then
        STA = S3TA/(C00004*CTA*CTA-C00001)
        K = CTA - STA*SPHI*C00IR3
        dKdLode =  - STA - C00IR3*SPHI*CTA
      else
        SGN = SIGN(C00001,Lode) ! It puts the Lode's sign to the number 1
        A = A1 + A2*SGN*SPHI
        B = B1*SGN + B2*SPHI
        K = A - B*S3TA
        dKdLode = - C00003*B*C3TA
      end if
      
      !Derivative Dp/DSig
      DpDSig(1) = CP3333
      DpDSig(2) = CP3333
      DpDSig(3) = CP3333
      DpDSig(4) = C00000
      DpDSig(5) = C00000
      DpDSig(6) = C00000
      
      !Derivative DJ/DSig
      i1 = C000P5/J
      if (J < 0.0001) then
        i1 = 0.0d0
      end if
      Sx = Sig(1)-p
      Sy = Sig(2)-p
      Sz = Sig(3)-p

      DJDSig(1) = i1 * Sx
      DJDSig(2) = i1 * Sy
      DJDSig(3) = i1 * Sz
      DJDSig(4) = i1 * C00002 * Sig(4)
      DJDSig(5) = i1 * C00002 * Sig(5)
      DJDSig(6) = i1 * C00002 * Sig(6)

      !Derivative DJ3/DSig
      i2 = CP3333*J*J
      DJ3DSig(1) = (Sy*Sz - Sig(5)*Sig(5) + i2)
      DJ3DSig(2) = (Sx*Sz - Sig(6)*Sig(6) + i2)
      DJ3DSig(3) = (Sx*Sy - Sig(4)*Sig(4) + i2)
      DJ3DSig(4) = C00002*(Sig(5)*Sig(6) - Sz*Sig(4))
      DJ3DSig(5) = C00002*(Sig(6)*Sig(4) - Sx*Sig(5))
      DJ3DSig(6) = C00002*(Sig(4)*Sig(5) - Sy*Sig(6))

      D = J*K/(sqrt(J2*K*K + ASPHI2))

      !C1F
      C1 = SPHI
      !C2F
      C2 = D*K - T3TA*D*dKdLode
      !C3F
      C3 = -C0R3I2*dKdLode*D/(J2*C3TA)
      
      !DFDSig!
      do i=1,6
        DFDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! DPPDSig = DFDSig (if associated Flow Rule)  !!!!!!!!!!!!!!!!!!!!!!
      !!!!! or
      !!!!! DPPDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (abs(J) < J0) then
        psi2 = phi - abs(J)*(phi - psi)/J0
      else
        psi2 = psi
      end if

      if (phi == psi2) then !If Associated Flow Rule, then DPPDSig = DFDSig
        DPPDSig = DFDSig

      else !If Non-Associated Flow Rule, then calculate...
        
          
          
          !Material parameters
          SPSI = sin(psi2) 
        CPSI = cos(psi2)
        
        if (abs(psi2) == C00000)then
            ASPSI2 = C00000
        else 
            COTPSI = CPSI/SPSI
        aSmooth = C00P50*COH*COTPSI !Smoothing parameter
        ASPSI2 = aSmooth*aSmooth*SPSI*SPSI
        end if
            
        
        
        
        
        

        !Calculate K function and its derivative
        if (abs(Lode) <= ATTRAN) then
            K = CTA - STA*SPSI*C00IR3
            dKdLode = - STA - C00IR3*SPSI*CTA
        else
            A = A1 + A2*SGN*SPSI
            B = B1*SGN + B2*SPSI
            K = A - B*S3TA
            dKdLode = - C00003*B*C3TA
        end if

        D = J*K/(sqrt(J*J*K*K + ASPSI2))

        !C1F
        C1 = SPSI
        !C2F
        C2 = D*K - T3TA*D*dKdLode
        !C3F
        C3 = -C0R3I2*dKdLode*D/(J2*C3TA)

        !DPPDSig
        do i=1,6
            DPPDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
        end do

      end if

      end subroutine CalculateDerivativesYieldFunctAndPlasticPotential


      Subroutine CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the yield function (F) with respect the strength parameters
      ! The strength parameters are: cohesion (COH) and friction angle (PHI)
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B,&
                         Denom, Num, aSmooth, ASPHI2, SGN
      double precision :: K, dKdPhi, dadc, dadPhi
      double precision, parameter :: C00001 = 1.0D0 !Parameters
      double precision, parameter :: C00P50 = 0.0005D0
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0
      double precision, parameter :: C000P1 = 0.00000000001D0
      !Constants for rounded K function (for LodeT=25)
      !real(REAL_TYPE), parameter :: A1 = 1.432052062044227d0
      !real(REAL_TYPE), parameter :: A2 = 0.406941858374615d0
      !real(REAL_TYPE), parameter :: B1 = 0.544290524902313d0
      !real(REAL_TYPE), parameter :: B2 = 0.673903324498392d0
      !real(REAL_TYPE), parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !real(REAL_TYPE), parameter :: A1 = -138300705.446275
      !real(REAL_TYPE), parameter :: A2 = -138300706.472675
      !real(REAL_TYPE), parameter :: B1 = -138300706.3123
      !real(REAL_TYPE), parameter :: B2 = 0.192450089729875
      !real(REAL_TYPE), parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians

      !In variables
      double precision, intent(in) :: p,J,Lode,S3TA !Invariants
      double precision, intent(in) :: c,phi !Soft Parameters
      !Out variables
      double precision, intent(out), dimension(2) :: DFDSP !Derivatives respect Soft Parameters


      !Material parameters
      COH = c !Cohesion
      SPHI = sin(phi) 
      CPHI = cos(phi)
      COTPHI = CPHI/SPHI

      !Calculate aSmooth and its derivatives
      if (abs(phi) == C00000) then
        COTPHI = C00000
        dadc = C00000
        dadPhi = C00000
      else
        dadc = C00P50*CPHI/SPHI
        dadPhi = - C00P50*COH/(SPHI*SPHI)
      end if
      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      if (abs(phi) == C00000) then
       ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      !Calculate K function and its derivatives
      if (abs(Lode) <= ATTRAN) then
        STA = sin(Lode)
        CTA = cos(Lode)
        K = CTA - STA*SPHI*C00IR3
        dKdPhi = - C00IR3*CPHI*STA
      else
        SGN = SIGN(C00001,Lode) !It puts the Lode's sign to the number 1
        A = A1 + A2*SGN*SPHI
        B = B1*SGN + B2*SPHI
        K = A - B*S3TA
        dKdPhi = A2*SGN*CPHI - B2*CPHI*S3TA
      end if

      !Operating..
      Denom = (sqrt(J*J*K*K + ASPHI2))
      Num =  J*J*K*dKdPhi + aSmooth*SPHI*SPHI*dadPhi + aSmooth*aSmooth*SPHI*CPHI

      !Derivative DF/Dc
      DFDSP(1) = aSmooth*SPHI*SPHI*dadc/Denom - CPHI

      !Derivative DF/Dphi
      DFDSP(2) = p*CPHI + Num/Denom + COH*SPHI

      if (J <= C00000) then
        DFDSP(1) = - CPHI
        DFDSP(2) = p*CPHI + COH*SPHI
      end if

      end subroutine CalculateDerivativesYieldFunctSofteningParameters


      subroutine CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,&
                              phip,phir,psip,psir,EpsPEq,DSPDPEq)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the strength parameters with respect
      ! the equivalent plastic shear strain
      !
      !**********************************************************************

      implicit none

      !In Variables
      double precision, intent(in) :: EpsPEq
      double precision, intent(in) :: factor,cp,cr,phip,phir,psip,psir
      !Out Variables
      double precision, intent(out), dimension(3):: DSPDPEq
     
      !Derivative Cohesion respect Equivalent Plastic Strain (Dc/DPEq)
      DSPDPEq(1) = -factor * (cp - cr) * (exp(-factor*EpsPEq))
      !Derivative Friction angle respect Equivalent Plastic Strain (Dphi/DPEq)
      DSPDPEq(2) = -factor * (phip - phir) * (exp(-factor*EpsPEq))
      !Derivative Dilatancy angle respect Equivalent Plastic Strain (Dpsi/DPEq)
      DSPDPEq(3) = -factor * (psip - psir) * (exp(-factor*EpsPEq))

      end subroutine CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain


      subroutine CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the equivalent plastic shear strain
      ! with respect the plastic strain
      !
      !**********************************************************************

      implicit none

      !Local Variables
      double precision :: k1, k2, k3
      double precision :: EpsPM
      double precision, dimension(3) :: EpsDev
      !In Variables
      double precision, intent(in), dimension(6) :: EpsP
      double precision, intent(in) :: EpsPEq
      !Out Variables
      double precision, intent(out), dimension(6):: DEpsPEqDPS

      !k1 = 2.0d0/(3.0d0*EpsPEq) 
      if (EpsPEq < 0.00000000001d0) then
        k1 = 0.0d0
      else 
        k1 = 2.0d0/(3.0d0*EpsPEq)  
      end if
      k2 = k1 * 1.0d0/3.0d0
      k3 = k1 * 2.0d0

      EpsPM = k2 * (EpsP(1) + EpsP(2) + EpsP(3))
      EpsDev(1) = EpsP(1)-EpsPM
      EpsDev(2) = EpsP(2)-EpsPM
      EpsDev(3) = EpsP(3)-EpsPM

      DEpsPEqDPS(1) = k2 * ( 2.0d0*EpsDev(1) - EpsDev(2) - EpsDev(3))
      DEpsPEqDPS(2) = k2 * (-EpsDev(1) + 2.0d0*EpsDev(2) - EpsDev(3))
      DEpsPEqDPS(3) = k2 * (-EpsDev(1) - EpsDev(2) + 2.0d0*EpsDev(3))
      DEpsPEqDPS(4) = k3 * EpsP(4)
      DEpsPEqDPS(5) = k3 * EpsP(5)
      DEpsPEqDPS(6) = k3 * EpsP(6)

      end subroutine CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain


      subroutine CalculateEpsPEq(EpsP,EpsPEq)
      !**********************************************************************
      !
      ! Calculation of the equivalent plastic shear strain 
      !
      !**********************************************************************
      
      implicit none

      !Local variables
      double precision:: EpsPM, C1, C2
      double precision, dimension(3) :: EpsDev
      !In variables
      double precision, intent(in), dimension(6) :: EpsP
      !Out variables
      double precision, intent(out) :: EpsPEq
      
      !EpsPEq = ((2/3)ep:ep)^(1/2), ep is the deviatoric plastic strain
      
      EpsPM = (1.0d0/3.0d0) * (EpsP(1) + EpsP(2) + EpsP(3))
      EpsDev(1) = EpsP(1)-EpsPM
      EpsDev(2) = EpsP(2)-EpsPM
      EpsDev(3) = EpsP(3)-EpsPM
      C1 = 2.0d0/3.0d0
      C2 = C1 * 2.0d0
      
      EpsPEq = sqrt((C1*EpsDev(1)*EpsDev(1)) + (C1*EpsDev(2)*EpsDev(2)) +  (C1*EpsDev(3)*EpsDev(3)) +&
                     (C2*EpsP(4)*EpsP(4)) + (C2*EpsP(5)*EpsP(5)) + (C2*EpsP(6)*EpsP(6)))

      end subroutine CalculateEpsPEq
      

      !Subroutine CalculateIncrementSofteningParameters(DSPDPEq,DEpsPEqDPS,DEpsP,Dh)
      !!**********************************************************************
      !!
      !! Calculation of the increment of the strenght parameters due to the softening
      !!
      !!**********************************************************************
      !
      !implicit none
      !
      !!Local variables
      !real(REAL_TYPE) :: k
      !!In variables
      !real(REAL_TYPE), intent(in), dimension(3) :: DSPDPEq
      !real(REAL_TYPE), intent(in), dimension(6) :: DEpsPEqDPS
      !real(REAL_TYPE), intent(in), dimension(6) :: DEpsP
      !!Out variables
      !real(REAL_TYPE), intent(out), dimension(3) :: Dh
      !
      !
      !k = DEpsPEqDPS(1)*DEpsP(1) + DEpsPEqDPS(2)*DEpsP(2) + DEpsPEqDPS(3)*DEpsP(3) + 
      !*       DEpsPEqDPS(4)*DEpsP(4) + DEpsPEqDPS(5)*DEpsP(5) + DEpsPEqDPS(6)*DEpsP(6)
      
      
      !Dh(1) = DSPDPEq(1)*k
      !Dh(2) = DSPDPEq(2)*k
      !Dh(3) = DSPDPEq(3)*k
      
      !Dh(1) = min (Dh(1) , 0.0d0)
      !Dh(2) = min (Dh(2) , 0.0d0)
      !Dh(3) = min (Dh(3) , 0.0d0)
      
      !end subroutine CalculateIncrementSofteningParameters


      Subroutine CalculateSofteningParameters(EpsPEq,factor,cp,cr,phip,phir,psip,psir,c,phi,psi)
      !**********************************************************************
      !
      ! Calculation of strenght parameters (c, phi, psi)
      !
      !**********************************************************************

      implicit none

      !In variables
      double precision, intent(in) :: EpsPEq,factor,cp,cr,phip,phir,psip,psir
      !Out variables
      double precision, intent(out) :: c,phi,psi  

      c = cr + (cp-cr)*exp(-factor*EpsPEq) 
      phi = phir + (phip-phir)*exp(-factor*EpsPEq) 
      psi = psir + (psip-psir)*exp(-factor*EpsPEq) 

      end subroutine CalculateSofteningParameters


      Subroutine DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c,phi,psi,Sig,DEpsPEqDPS,DSPDPEq,DEps,DSig,DEpsP)
      !**********************************************************************
      !
      ! Calculation of the stress increment and plastic strain increment
      !
      !         dSig = Dep * dEps
      !         dEpsP = Lambda * DPDSig
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i,k
      double precision :: A,Ai,Denom,Fact,LambdaNum,Lambda
      double precision :: p,J,Lode,S3TA !Invariants
      double precision, dimension(6,6) :: Num,Num1,Prod
      double precision, dimension(6) :: Denom1
      double precision, dimension(6) :: DPPDSig !Derivatives Plastic potential respect net stress
      double precision, dimension(6) :: DFDSig !Derivatives Yield function respect net stress
      double precision, dimension(2) :: DFDSP !Derivatives Yield function respect Soft Parameters
      double precision, dimension(6,6) :: Dep !Elastoplastic Constitutive Matrix 
      !In Variables
      double precision, intent(in) :: c,phi,psi !Softening parameters
      double precision, intent(in) :: D1,D2,GG !Elastic parameters
      double precision, intent(in), dimension(6):: DEpsPEqDPS
      double precision, intent(in), dimension(6) :: Sig
      double precision, intent(in), dimension(3) :: DSPDPEq !Derivatives respect Equivalent Plastic Strain
      double precision, intent(in), dimension(6) :: DEps
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      !Out Variables
      double precision, intent(out), dimension(6) :: DSig
      double precision, intent(out), dimension(6) :: DEpsP

      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      call CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
      call CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)

      !Parameter A (H = -A --> A>0 softening / A<0 hardening)
      A = 0.0d0
      Ai = (DFDSP(1)*DSPDPEq(1) + DFDSP(2)*DSPDPEq(2))
      do i=1,6
      A = A + Ai * DEpsPEqDPS(i) * DPPDSig(i)
      end do

      !Elastoplastic Constitutive Matrix (Dep)
      do i=1,6
        do k=1,6
            Prod(i,k) =  DPPDSig(i) * DFDSig(k)
        end do
      end do

      Num1(1,1) = D1*Prod(1,1) + D2*Prod(2,1) + D2*Prod(3,1)
      Num1(1,2) = D1*Prod(1,2) + D2*Prod(2,2) + D2*Prod(3,2)
      Num1(1,3) = D1*Prod(1,3) + D2*Prod(2,3) + D2*Prod(3,3)
      Num1(1,4) = D1*Prod(1,4) + D2*Prod(2,4) + D2*Prod(3,4)
      Num1(1,5) = D1*Prod(1,5) + D2*Prod(2,5) + D2*Prod(3,5)
      Num1(1,6) = D1*Prod(1,6) + D2*Prod(2,6) + D2*Prod(3,6)

      Num1(2,1) = D2*Prod(1,1) + D1*Prod(2,1) + D2*Prod(3,1)
      Num1(2,2) = D2*Prod(1,2) + D1*Prod(2,2) + D2*Prod(3,2)
      Num1(2,3) = D2*Prod(1,3) + D1*Prod(2,3) + D2*Prod(3,3)
      Num1(2,4) = D2*Prod(1,4) + D1*Prod(2,4) + D2*Prod(3,4)
      Num1(2,5) = D2*Prod(1,5) + D1*Prod(2,5) + D2*Prod(3,5)
      Num1(2,6) = D2*Prod(1,6) + D1*Prod(2,6) + D2*Prod(3,6)

      Num1(3,1) = D2*Prod(1,1) + D2*Prod(2,1) + D1*Prod(3,1)
      Num1(3,2) = D2*Prod(1,2) + D2*Prod(2,2) + D1*Prod(3,2)
      Num1(3,3) = D2*Prod(1,3) + D2*Prod(2,3) + D1*Prod(3,3)
      Num1(3,4) = D2*Prod(1,4) + D2*Prod(2,4) + D1*Prod(3,4)
      Num1(3,5) = D2*Prod(1,5) + D2*Prod(2,5) + D1*Prod(3,5)
      Num1(3,6) = D2*Prod(1,6) + D2*Prod(2,6) + D1*Prod(3,6)

      Num1(4,1) = GG*Prod(4,1)
      Num1(4,2) = GG*Prod(4,2)
      Num1(4,3) = GG*Prod(4,3)
      Num1(4,4) = GG*Prod(4,4)
      Num1(4,5) = GG*Prod(4,5)
      Num1(4,6) = GG*Prod(4,6)

      Num1(5,1) = GG*Prod(5,1)
      Num1(5,2) = GG*Prod(5,2)
      Num1(5,3) = GG*Prod(5,3)
      Num1(5,4) = GG*Prod(5,4)
      Num1(5,5) = GG*Prod(5,5)
      Num1(5,6) = GG*Prod(5,6)

      Num1(6,1) = GG*Prod(6,1)
      Num1(6,2) = GG*Prod(6,2)
      Num1(6,3) = GG*Prod(6,3)
      Num1(6,4) = GG*Prod(6,4)
      Num1(6,5) = GG*Prod(6,5)
      Num1(6,6) = GG*Prod(6,6)



      Num(1,1) = D1*Num1(1,1) + D2*Num1(1,2) + D2*Num1(1,3)
      Num(1,2) = D2*Num1(1,1) + D1*Num1(1,2) + D2*Num1(1,3)
      Num(1,3) = D2*Num1(1,1) + D2*Num1(1,2) + D1*Num1(1,3)
      Num(1,4) = GG*Num1(1,4)
      Num(1,5) = GG*Num1(1,5)
      Num(1,6) = GG*Num1(1,6)

      Num(2,1) = D1*Num1(2,1) + D2*Num1(2,2) + D2*Num1(2,3)
      Num(2,2) = D2*Num1(2,1) + D1*Num1(2,2) + D2*Num1(2,3)
      Num(2,3) = D2*Num1(2,1) + D2*Num1(2,2) + D1*Num1(2,3)
      Num(2,4) = GG*Num1(2,4)
      Num(2,5) = GG*Num1(2,5)
      Num(2,6) = GG*Num1(2,6)

      Num(3,1) = D1*Num1(3,1) + D2*Num1(3,2) + D2*Num1(3,3)
      Num(3,2) = D2*Num1(3,1) + D1*Num1(3,2) + D2*Num1(3,3)
      Num(3,3) = D2*Num1(3,1) + D2*Num1(3,2) + D1*Num1(3,3)
      Num(3,4) = GG*Num1(3,4)
      Num(3,5) = GG*Num1(3,5)
      Num(3,6) = GG*Num1(3,6)

      Num(4,1) = D1*Num1(4,1) + D2*Num1(4,2) + D2*Num1(4,3)
      Num(4,2) = D2*Num1(4,1) + D1*Num1(4,2) + D2*Num1(4,3)
      Num(4,3) = D2*Num1(4,1) + D2*Num1(4,2) + D1*Num1(4,3)
      Num(4,4) = GG*Num1(4,4)
      Num(4,5) = GG*Num1(4,5)
      Num(4,6) = GG*Num1(4,6)

      Num(5,1) = D1*Num1(5,1) + D2*Num1(5,2) + D2*Num1(5,3)
      Num(5,2) = D2*Num1(5,1) + D1*Num1(5,2) + D2*Num1(5,3)
      Num(5,3) = D2*Num1(5,1) + D2*Num1(5,2) + D1*Num1(5,3)
      Num(5,4) = GG*Num1(5,4)
      Num(5,5) = GG*Num1(5,5)
      Num(5,6) = GG*Num1(5,6)

      Num(6,1) = D1*Num1(6,1) + D2*Num1(6,2) + D2*Num1(6,3)
      Num(6,2) = D2*Num1(6,1) + D1*Num1(6,2) + D2*Num1(6,3)
      Num(6,3) = D2*Num1(6,1) + D2*Num1(6,2) + D1*Num1(6,3)
      Num(6,4) = GG*Num1(6,4)
      Num(6,5) = GG*Num1(6,5)
      Num(6,6) = GG*Num1(6,6)



      Denom1(1) = DFDSig(1)*D1 + DFDSig(2)*D2 + DFDSig(3)*D2
      Denom1(2) = DFDSig(1)*D2 + DFDSig(2)*D1 + DFDSig(3)*D2
      Denom1(3) = DFDSig(1)*D2 + DFDSig(2)*D2 + DFDSig(3)*D1
      Denom1(4) = DFDSig(4)*GG
      Denom1(5) = DFDSig(5)*GG
      Denom1(6) = DFDSig(6)*GG

      Denom =   Denom1(1)*DPPDSig(1) + Denom1(2)*DPPDSig(2) + &
                  Denom1(3)*DPPDSig(3) + Denom1(4)*DPPDSig(4) + &
             Denom1(5)*DPPDSig(5) + Denom1(6)*DPPDSig(6) - A

      Fact = 1d0/Denom

      !Dep
      Dep(1,1) = D1 - Fact*Num(1,1)
      Dep(1,2) = D2 - Fact*Num(1,2)
      Dep(1,3) = D2 - Fact*Num(1,3)
      Dep(1,4) = -Fact*Num(1,4)
      Dep(1,5) = -Fact*Num(1,5)
      Dep(1,6) = -Fact*Num(1,6)

      Dep(2,1) = D2 - Fact*Num(2,1)
      Dep(2,2) = D1 - Fact*Num(2,2)
      Dep(2,3) = D2 - Fact*Num(2,3)
      Dep(2,4) = -Fact*Num(2,4)
      Dep(2,5) = -Fact*Num(2,5)
      Dep(2,6) = -Fact*Num(2,6)

      Dep(3,1) = D2 - Fact*Num(3,1)
      Dep(3,2) = D2 - Fact*Num(3,2)
      Dep(3,3) = D1 - Fact*Num(3,3)
      Dep(3,4) = -Fact*Num(3,4)
      Dep(3,5) = -Fact*Num(3,5)
      Dep(3,6) = -Fact*Num(3,6)

      Dep(4,1) = -Fact*Num(4,1)
      Dep(4,2) = -Fact*Num(4,2)
      Dep(4,3) = -Fact*Num(4,3)
      Dep(4,4) = GG - Fact*Num(4,4)
      Dep(4,5) = -Fact*Num(4,5)
      Dep(4,6) = -Fact*Num(4,6)

      Dep(5,1) = -Fact*Num(5,1)
      Dep(5,2) = -Fact*Num(5,2)
      Dep(5,3) = -Fact*Num(5,3)
      Dep(5,4) = -Fact*Num(5,4)
      Dep(5,5) = GG - Fact*Num(5,5)
      Dep(5,6) = -Fact*Num(5,6)

      Dep(6,1) = -Fact*Num(6,1)
      Dep(6,2) = -Fact*Num(6,2)
      Dep(6,3) = -Fact*Num(6,3)
      Dep(6,4) = -Fact*Num(6,4)
      Dep(6,5) = -Fact*Num(6,5)
      Dep(6,6) = GG - Fact*Num(6,6)

      !!!!!!!!! Calculate Plastic multipliler(Lambda)!!!!!!!!!!!!!!!!!
      LambdaNum =   Denom1(1)*DEps(1) + Denom1(2)*DEps(2) + &
                   Denom1(3)*DEps(3) + Denom1(4)*DEps(4) + &
                   Denom1(5)*DEps(5) + Denom1(6)*DEps(6) 
      Lambda =  LambdaNum/Denom

      !!!!!!!!! Determine DSig --> (DSig = Dep*dEps) !!!!!!!!!!!
      do i=1,6
        DSig(i) = 0.0d0
        do k=1,6
            DSig(i) =  DSig(i) + Dep(i,k) * DEps(k)
        end do
      end do

      !!!!!!!!! Determine DEpsP --> (DEpsP = Lambda*DPDSig) !!!!!!!!!!!!
      do i=1,6
        DEpsP(i) = Lambda * DPPDSig(i)
      end do

      end subroutine DetermineDSigAndDEpsP


      subroutine EndOfStepCorrection(IntGlo,D1,D2,GG,IPL,F,Sig,DSPDPEq,DEpsPEqDPS,EpsP,c,phi,psi)
      !**********************************************************************
      !
      ! Final correction of the yield surface drift (END OF STEP CORRECTION).
      ! The stresses, the plastic strain and the strength parameters are corrected.
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i
      double precision :: p,J,Lode,S3TA !Invariants
      double precision :: Lambda,param,c2,phi2,psi2,F2
      double precision :: Denom,A,Ai
      double precision, dimension(2) :: DFDSP
      double precision, dimension(6) :: DPPDSig,DFDSig,Sig2,DEpsP,EpsP2
      double precision, dimension(6) :: Denom1
      double precision, dimension(3) :: Dh
      !In Variables
      integer, intent(in) :: IntGlo,IPL !Global ID of Gauss point or particle
      double precision, intent(in):: D1,D2,GG
      double precision, intent(in), dimension(3) :: DSPDPEq !Derivatives respect Equivalent Plastic Strain
      double precision, intent(in), dimension(6) :: DEpsPEqDPS !Derivatives respect Equivalent Plastic Strain
      !InOut Variables
      double precision, intent(inout):: c,phi,psi
      double precision, intent(inout), dimension(6) :: Sig
      double precision, intent(inout), dimension(6) :: EpsP
      double precision, intent(inout):: F

      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      call CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
      call CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)

      !Parameter A (hardening/softening parameter)
      A = 0.0d0
      Ai = (DFDSP(1)*DSPDPEq(1) + DFDSP(2)*DSPDPEq(2))
      do i=1,6
        A = A + Ai * DEpsPEqDPS(i) * DPPDSig(i)
      end do

      Denom1(1) = DPPDSig(1)*D1 + DPPDSig(2)*D2 + DPPDSig(3)*D2
      Denom1(2) = DPPDSig(1)*D2 + DPPDSig(2)*D1 + DPPDSig(3)*D2
      Denom1(3) = DPPDSig(1)*D2 + DPPDSig(2)*D2 + DPPDSig(3)*D1
      Denom1(4) = DPPDSig(4)*GG
      Denom1(5) = DPPDSig(5)*GG
      Denom1(6) = DPPDSig(6)*GG

      Denom = Denom1(1)*DFDSig(1) + Denom1(2)*DFDSig(2) + &
             Denom1(3)*DFDSig(3) + Denom1(4)*DFDSig(4) + &
              Denom1(5)*DFDSig(5) + Denom1(6)*DFDSig(6) - A

      Lambda = F/Denom !factor correction

      Sig2 = Sig - Lambda * Denom1 ! Sig2 = Sig + fact * Denom1 Stress corrected
      DEpsP = Lambda * DPPDSig
      EpsP2 = EpsP + DEpsP

      if (IPL == 1)then
        Dh = 0.0d0
      else
        param = DEpsPEqDPS(1) * DEpsP(1) + DEpsPEqDPS(2) * DEpsP(2) + DEpsPEqDPS(3) * DEpsP(3) + &
               DEpsPEqDPS(4) * DEpsP(4) + DEpsPEqDPS(5) * DEpsP(5) + DEpsPEqDPS(6) * DEpsP(6)
        Dh(1) = min (DSPDPEq(1)*param, 0.0d0)
        Dh(2) = min (DSPDPEq(2)*param, 0.0d0)
        Dh(3) = min (DSPDPEq(3)*param, 0.0d0)
      end if

      c2 = c + Dh(1)
      phi2 = phi + Dh(2)
      psi2 = psi + Dh(3)

      call DetermineYieldFunctionValue(IntGlo,Sig2,c2,phi2,F2)
      
      if ((abs(F2) > abs(F)).or.(Denom == 0.0d0)) then !NormalCorrectionScheme
        Denom = 0.0d0
        Denom = DFDSig(1)*DFDSig(1) + DFDSig(2)*DFDSig(2) + &
                 DFDSig(3)*DFDSig(3) + DFDSig(4)*DFDSig(4) + &
                 DFDSig(5)*DFDSig(5) + DFDSig(6)*DFDSig(6)

        Lambda = F/Denom
        Sig = Sig - Lambda * DFDSig
        DEpsP = Lambda * DPPDSig
        EpsP = EpsP + DEpsP
        call DetermineYieldFunctionValue(IntGlo,Sig,c,phi,F)
      else
        Sig = Sig2
        EpsP = EpsP2
        c = c2
        phi = phi2
        psi = psi2
        F = F2
      end if

      end subroutine EndOfStepCorrection


      subroutine CalculatePrincipalStresses(IntGlo,Sig,SigPrin)
      !**********************************************************************
      !
      ! Implemented in the frame of the MPM project.
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision, dimension(3) :: xN1,xN2,xN3
      double precision :: Sig1,Sig2,Sig3,p,q
      !In Variables
      integer, intent(in) :: IntGlo ! Global ID of Gauss point or particle
      double precision, intent(in), dimension(6) :: Sig
      !Out Variables
      double precision, intent(out), dimension(6) :: SigPrin

      call PrincipalSig(1,Sig,xN1,xN2,xN3,Sig1,Sig2,Sig3,P,Q)

      If (Sig1 >= Sig2.and.Sig2 >= Sig3) then
        SigPrin(1) = Sig1
        SigPrin(2) = Sig2
        SigPrin(3) = Sig3
      else if (Sig1 >= Sig3.and.Sig3 >= Sig2) then
        SigPrin(1) = Sig1
        SigPrin(2) = Sig3
        SigPrin(3) = Sig2
      else if (Sig3 >= Sig1.and.Sig1 >= Sig2) then
        SigPrin(1) = Sig3
        SigPrin(2) = Sig1
        SigPrin(3) = Sig2
      else if (Sig3 >= Sig2.and.Sig2 >= Sig1) then
        SigPrin(1) = Sig3
        SigPrin(2) = Sig2
        SigPrin(3) = Sig1
      else if (Sig2 >= Sig1.and.Sig1 >= Sig3) then
        SigPrin(1) = Sig2
        SigPrin(2) = Sig1
        SigPrin(3) = Sig3
      else if (Sig2 >= Sig3.and.Sig3 >= Sig1) then
        SigPrin(1) = Sig2
        SigPrin(2) = Sig3
        SigPrin(3) = Sig1
      end if

        SigPrin(4) = 0.0d0
        SigPrin(5) = 0.0d0
        SigPrin(6) = 0.0d0

      end subroutine CalculatePrincipalStresses


      Subroutine PrincipalSig(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit double precision (A-H,O-Z)
      Dimension S(*),xN1(*),xN2(*),xN3(*)

      If (iOpt.Eq.1) Then
        Call Eig_3(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
      Else
        Call Eig_3a(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
      End If
      Return
      End
      
      
      Subroutine Eig_3(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit double precision (A-H,O-Z)
      Dimension St(6),A(3,3),V(3,3),xN1(3),xN2(3),xN3(3)
!     *          xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
      !
      
      !A = 0.0
      
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! Set V to unity matrix
      V(1,1) = 1
      V(2,1) = 0
      V(3,1) = 0

      V(1,2) = 0
      V(2,2) = 1
      V(3,2) = 0

      V(1,3) = 0
      V(2,3) = 0
      V(3,3) = 1


      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      it = 0
      itmax = 50
      Do While ( it.Lt.itMax .And. abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
!     *           abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
          If (a(ip,iq) .Ne. 0.0) Then
              !if (abs(a(ip,iq))<1e-20) then 
              !    a(ip,iq)=1e-20
              !end if 

              tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
              
              
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            v1p=c*v(1,ip)-s*v(1,iq)
            v2p=c*v(2,ip)-s*v(2,iq)
            v3p=c*v(3,ip)-s*v(3,iq)
            v(1,iq)=s*v(1,ip)+c*v(1,iq)
            v(2,iq)=s*v(2,ip)+c*v(2,iq)
            v(3,iq)=s*v(3,ip)+c*v(3,iq)
            v(1,ip)=v1p
            v(2,ip)=v2p
            v(3,ip)=v3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      ! Sort eigenvalues S1 <= S2 <= S3
      is1 = 1
      is2 = 2
      is3 = 3
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
        it  = is3
        is3 = is2
        is2 = it
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      Do i=1,3
        xN1(i) = v(i,is1) ! first  column
        xN2(i) = v(i,is2) ! second column
        xN3(i) = v(i,is3) ! third  column
      End Do
      Return
      End ! Eig_3

      
      Subroutine Eig_3a(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      Implicit double precision (A-H,O-Z)
      Dimension St(6),A(3,3)   !  V(3,3),xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      If (iOpt.Eq.1) Tol = 1d-50*abs_max_s
      it=0
      itmax = 50
      Do While ( it.lt.itmax .And.&
                abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )

        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
          If (a(ip,iq) .Ne. 0.0) Then         ! ongelijk nul ?
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      Return
      End ! Eig_3a


      Subroutine MatVec(xMat,IM,Vec,N,VecR)
!C***********************************************************************
!C
!C     Calculate VecR = xMat*Vec
!C
!C I   xMat  : (Square) Matrix (IM,*)
!C I   Vec   : Vector
!C I   N     : Number of rows/colums
!C O   VecR  : Resulting vector
!C
!C***********************************************************************
      Implicit double precision (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
!C***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine MatVec

      
      Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
!C***********************************************************************
!C
!C     Calculate VecR() = R1*Vec1()+R2*Vec2()
!C
!C I   Vec1,
!C I   Vec2  : Vectors
!C I   R1,R2 : Multipliers
!C I   N     : Number of rows
!C O   VecR  : Resulting vector
!C
!C***********************************************************************
      Implicit double precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
!C***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine AddVec

      subroutine Get_EigenValues_EigenVectors(Sig0, Principal_stresses, Principal_vectors)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Obtains the maximum angle of principal stress rotation      !!!
      !!! A generalization is proposed for 3D cases in which the      !!!
	!!! the rotation is measure from hydrostatic axis in x,y,z      !!!
	!!! coordinates, to hydrostatic axis in x', y', z' coordinates  !!!
      !!!																!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none 
      ! Local
      double precision :: TENSOR(3,3)
 	double precision:: fv1(3),fv2(3), hyd_new(3)
	integer :: I, ierr, J
    
      ! Input
      double precision, dimension(6), intent(in) :: Sig0
    
      ! Output 
      double precision, dimension(3) :: Principal_stresses
      double precision, dimension(3,3) :: Principal_vectors
	
	
	!____ First assemble Tensor for principal stress calc________________________________________
	do I= 1,3
		TENSOR(I,I)=Sig0(I)		
	enddo
	TENSOR(1,2)=Sig0(4)
	TENSOR(2,1)=Sig0(4)
 	TENSOR(1,3)=Sig0(5)
	TENSOR(3,1)=Sig0(5) 
	TENSOR(2,3)=Sig0(6)
	TENSOR(3,2)=Sig0(6)
	!____________________________________________________________________________________________
	!____Call subroutine to find principal vectors_______________________________________________
	call rs(3,3,TENSOR,Principal_stresses,1,Principal_vectors,fv1,fv2,ierr)
	!____________________________________________________________________________________________
      !___Find new hydrostatic axis vector_________________________________________________________
	!hyd_new=0.0d0
	!do I=1,3
	!	do J=1,3
	!		hyd_new(I)=hyd_new(I)+Principal_vectors(I,J)
	!	enddo		
	!enddo
	! normalize the new hydrostatic axis vector__________________________________________________
	!call UnitaryVector(hyd_new, 3, hyd_new)
	!!____________________________________________________________________________________________
	!! Obtain angle of PSR _______________________________________________________________________
	!!Calculate hyd_new dot (1/sqrt3, 1/sqrt3, 1/sqrt3)
	!DirectionalCosine=0.0d0
	!do I=1,3
	!	DirectionalCosine=DirectionalCosine+(hyd_new(I)/(3.0d0**0.5))
	!enddo
	!Angle_PSR=acos(DirectionalCosine)
	!!____________________________________________________________________________________________
	!!Now Compute the change of PSR angle_________________________________________________________
	!AngleOfDifferentialPSR=Angle_PSR-Angle_PSR_0
	!____________________________________________________________________________________________
      end subroutine

	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!						EISPACK SUBROUTINES (This is a linear Algebra efficient package)
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)

      integer:: n,nm,ierr,matz
      double precision:: a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)

      !   this subroutine calls the recommended sequence of
      !   subroutines from the eigensystem subroutine package (eispack)
      !   to find the eigenvalues and eigenvectors (if desired)
      !   of a real symmetric matrix.
      !c
      !   on input
      !c
      !      nm  must be set to the row dimension of the two-dimensional
      !      array parameters as declared in the calling program
      !      dimension statement.
      !c
      !      n  is the order of the matrix  a.
      !c
      !      a  contains the real symmetric matrix.
      !c
      !      matz  is an integer variable set equal to zero if
      !      only eigenvalues are desired.  otherwise it is set to
      !      any non-zero integer for both eigenvalues and eigenvectors.
      !c
      !   on output
      !c
      !      w  contains the eigenvalues in ascending order.
      !c
      !      z  contains the eigenvectors if matz is not zero.
      !c
      !      ierr  is an integer output variable set equal to an error
      !         completion code described in the documentation for tqlrat
      !         and tql2.  the normal completion code is zero.
      !c
      !      fv1  and  fv2  are temporary storage arrays.
      !c
      !   questions and comments should be directed to burton s. garbow,
      !   mathematics and computer science div, argonne national laboratory
      !c
      !   this version dated august 1983.
      !c
      !   ------------------------------------------------------------------
      !c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50

   10 if (matz .ne. 0) go to 20
!c     .......... find eigenvalues only ..........
      !call  tred1(nm,n,a,w,fv1,fv2)
!*  tqlrat encounters catastrophic underflow on the Vax
!*     call  tqlrat(n,w,fv2,ierr)
      !call  tql1(n,w,fv1,ierr)
      go to 50
!c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
	end	
	

	
      subroutine tred2(nm,n,a,d,e,z)

      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale

!c     this subroutine is a translation of the algol procedure tred2,
!c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!c
!c     this subroutine reduces a real symmetric matrix to a
!c     symmetric tridiagonal matrix using and accumulating
!c     orthogonal similarity transformations.
!c
!c     on input
!c
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c
!c        n is the order of the matrix.
!c
!c        a contains the real symmetric input matrix.  only the
!c          lower triangle of the matrix need be supplied.
!c
!c     on output
!c
!c        d contains the diagonal elements of the tridiagonal matrix.
!c
!c        e contains the subdiagonal elements of the tridiagonal
!c          matrix in its last n-1 positions.  e(1) is set to zero.
!c
!c        z contains the orthogonal transformation matrix
!c          produced in the reduction.
!c
!c        a and z may coincide.  if distinct, a is unaltered.
!c
!c     questions and comments should be directed to burton s. garbow,
!c     mathematics and computer science div, argonne national laboratory
!c
!c     this version dated august 1983.
!c
!c     ------------------------------------------------------------------
!c
      do 100 i = 1, n

         do 80 j = i, n
   80    z(j,i) = a(j,i)

         d(i) = a(n,i)
  100 continue

      if (n .eq. 1) go to 510
!c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
!c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))

         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)

         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue

         go to 290

  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue

         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
!c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0

         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220

            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue

  220       e(j) = g
  240    continue
!c     .......... form p ..........
         f = 0.0d0

         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue

         hh = f / (h + h)
!c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
!c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)

            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)

            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue

  290    d(i) = h
  300 continue
!c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380

         do 330 k = 1, l
  330    d(k) = z(k,i) / h

         do 360 j = 1, l
            g = 0.0d0

            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)

            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue

  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0

  500 continue

  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue

      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
	end
	
	
	
	subroutine tql2(nm,n,d,e,z,ierr)

      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2!,pythag
!c
!c     this subroutine is a translation of the algol procedure tql2,
!c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!c     wilkinson.
!c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!c
!c     this subroutine finds the eigenvalues and eigenvectors
!c     of a symmetric tridiagonal matrix by the ql method.
!c     the eigenvectors of a full symmetric matrix can also
!c     be found if  tred2  has been used to reduce this
!c     full matrix to tridiagonal form.
!c
!c     on input
!c
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c
!c        n is the order of the matrix.
!c
!c        d contains the diagonal elements of the input matrix.
!c
!c        e contains the subdiagonal elements of the input matrix
!c          in its last n-1 positions.  e(1) is arbitrary.
!c
!c        z contains the transformation matrix produced in the
!c          reduction by  tred2, if performed.  if the eigenvectors
!c          of the tridiagonal matrix are desired, z must contain
!c          the identity matrix.
!c
!c      on output
!c
!c        d contains the eigenvalues in ascending order.  if an
!c          error exit is made, the eigenvalues are correct but
!c          unordered for indices 1,2,...,ierr-1.
!c
!c        e has been destroyed.
!c
!c        z contains orthonormal eigenvectors of the symmetric
!c          tridiagonal (or full) matrix.  if an error exit is made,
!c          z contains the eigenvectors associated with the stored
!c          eigenvalues.
!c
!c        ierr is set to
!c          zero       for normal return,
!c          j          if the j-th eigenvalue has not been
!c                     determined after 30 iterations.
!c
!c     calls pythag for  dsqrt(a*a + b*b) .
!c
!c     questions and comments should be directed to burton s. garbow,
!c     mathematics and computer science div, argonne national laboratory
!c
!c     this version dated august 1983.
!c
!c     ------------------------------------------------------------------
!c
      ierr = 0
      if (n .eq. 1) go to 1001

      do 100 i = 2, n
  100 e(i-1) = e(i)

      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0

      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
!c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
!c     .......... e(n) is always zero, so there is no exit
!c                through the bottom of the loop ..........
  110    continue

  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
!c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145

         do 140 i = l2, n
  140    d(i) = d(i) - h

  145    f = f + h
!c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue

  200    continue

         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)

         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue

         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p

         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue

  300 continue

      go to 1001
!c     .......... set error -- no convergence to an
!c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
	end
	
	
	double precision function pythag(a,b)
      double precision a,b
!c
!c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end

 

        
        
        
        
        
        
        
        
        
        
       ! *USER SUBROUTINES
!      SUBROUTINE UMAT_AM(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
!      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
!      !INCLUDE 'ABA_PARAM.INC'
!
!      CHARACTER*80 CMNAME
!      DIMENSION STRESS(NTENS),STATEV(NSTATEV), &
!      DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),&
!      STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
!      PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!
!       ! Arguments:
!        !          I/O  Type
!        !  PROPS    I   R()  : List with model parameters
!        !  DSTRAN   I   R()  : Strain increment
!        !  DDSDDE   O   R(,) : Material stiffness matrix
!        !  STRESS  I/O  R()  : stresses
!        !  STATEV  I/O  R()  : state variables
! 
!        !---  Local variables
!      DIMENSION dSig(NTENS), Sig(NTENS)
!      DIMENSION SigC(NTENS)
!      DIMENSION SigE(NTENS),SigEQ(NTENS)
!      Parameter (Ip=-1) ! I added an L here??
!      DIMENSION xN1(3),xN2(3),xN3(3),Tmp1(3),Tmp2(3),Tmp3(3)
!      IAPEX=0
!      ITENS=0
!      TauMax=1
!        ! Contents of PROPS(2)
!        !  1 : G       shear modulus
!        !  2 : ENU     Poisson's ratio
!        IntGlo= NPT
!        G = PROPS(1)   
!        ENU = PROPS(2)
!        VNU = PROPS(2)
!        one = 1.0d0
!        two = 2.0d0 
!        
!        SPHI = PROPS(3)
!        SPSI = PROPS(5)
!        COHS = PROPS(4)
!        TENS = PROPS(6)
!        ! calculate elastic stress increment (DSig = elastic stiffness D * strain increment DEps)
!        FAC = two * G / ( one - two * ENU )
!        D1 = FAC * ( one - ENU )
!        D2 = FAC * ENU
!        DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
!        dSig(1) = (D1 - D2) * DSTRAN(1) + D2 * DSTRANVOL
!        dSig(2) = (D1 - D2) * DSTRAN(2) + D2 * DSTRANVOL
!        dSig(3) = (D1 - D2) * DSTRAN(3) + D2 * DSTRANVOL
!        dSig(4) = G * DSTRAN(4)
!        if (NTENS == 6) then
!        dSig(5) = G * DSTRAN(5)
!        dSig(6) = G * DSTRAN(6)
!        end if
!      
!        ! elastic stress
!        SigE = STRESS + dSig
!        Stress=SigE
!        SigEQ = SigE
!        SigC = SigE
!        
!        DDSDDE = 0.0
!        DDSDDE(1:3,1:3) = D2
!        DDSDDE(1,1) = D1
!        DDSDDE(2,2) = D1
!        DDSDDE(3,3) = D1
!        DDSDDE(4,4) = G
!        if (NTENS == 6) then
!          DDSDDE(5,5) = G
!          DDSDDE(6,6) = G
!        end if
!      
!
!      
!!------------ Calculate principle stresses and their direction
!
!      call PrnSig(1,NTENS,SigE,xN1,xN2,xN3,Sig1,Sig2,Sig3,P,Q)
!
!  901       format(a,6(1x,1p,e12.5),1x,i4)
!
!
!
!!------------- evaluating the yield functions
!      IF (SPHI >  0.999d0) SPHI=0.999d0
!      F21 = 0.5d0*(SIG2-SIG1) + 0.5d0*(SIG2+SIG1)*SPHI - COHS
!      F32 = 0.5d0*(SIG3-SIG2) + 0.5d0*(SIG3+SIG2)*SPHI - COHS
!      F31 = 0.5d0*(SIG3-SIG1) + 0.5d0*(SIG3+SIG1)*SPHI - COHS
!      FT1 = SIG1 - TENS
!      FT2 = SIG2 - TENS
!      FT3 = SIG3 - TENS
!      IF (F31 > 1d-12 .OR. FT3 > 1d-12) GOTO 240
!      IF (IPL  == 0) GOTO 360
!!
!!         ***POINTS CHANGING FROM PLASTIC TO ELASTIC STATE***
!!
!      IPL =0
!      GOTO 360
!!
!!                   *** Plastic stress points ***
!!
!  240 IPL =1
!
!
!     
!      call PrnSig( 0,NTENS,SigEQ,Tmp1,Tmp2,Tmp3,SigV1,SigV2,SigV3,Dum1,Dum2)
!   
!
!      FF1 = 0.5d0*    (SIGV2-SIGV1) + 0.5d0*(SIGV2+SIGV1)*SPHI - COHS
!      FF2 = 0.5d0*DABS(SIGV2-SIGV3) + 0.5d0*(SIGV2+SIGV3)*SPHI - COHS
!      FF3 = 0.5d0*DABS(SIGV3-SIGV1) + 0.5d0*(SIGV3+SIGV1)*SPHI - COHS
!      
!
!
!      IF (FF1 > 1d-12 .OR. FF2 > 1d-12 .OR. FF3 > 1d-12) GOTO 250
!
!!       ****** PLASTIC NEGATIVE POINTS ******
!      IPL =-1
!!
!!         *** STRESSES ARE BROUGHT BACK TO THE YIELD SURFACE ***
!!
!  250 CONTINUE
!      VNUI0 = VNU !AR(INTGLO)
!
!      VNUI1 = 1 - 2*VNUI0
!      VNUI2 = 1 + 2*VNUI0
!      VNUIQ = VNUI2 / VNUI1
!!
!      DUM = SPSI / VNUI1
!      PSIMIN = G *(-1+DUM)
!      PSIMET = G *( 1+DUM)
!      PSINU = 2*G *VNUI0*DUM
!!
!      HA= ( 1 - SPHI + SPSI - SPHI * SPSI) * SIG1 + &
!           (-2 - 2/VNUI1*SPHI*SPSI) * SIG2 + &
!          ( 1 - SPHI - SPSI + VNUIQ * SPHI * SPSI) * SIG3 + &
!          2*(1 + SPSI) * COHS
!
!      HB= ( 1 + SPHI + SPSI + VNUIQ * SPHI * SPSI) * SIG1 -&
!          ( 2 + 2/VNUI1*SPHI*SPSI) * SIG2 + &
!          ( 1 + SPHI - SPSI - SPHI * SPSI) * SIG3 - &
!          2*(1 - SPSI) * COHS
!
!      HAB= (VNUI1+SPSI)*SIG1+ &
!           (VNUI1-SPSI)*(SIG3-TENS)- &
!           (VNUI1+SPSI)*(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
!     
!      HBA=(1-VNUI0)*SIG1-VNUI0*SIG3-(1-VNUI0)*(TENS*(1+SPHI)- &
!           2*COHS)/(1-SPHI)+VNUI0*TENS
!
!      HAO=(1-VNUI0)*SIG2-VNUI0*SIG3-VNUI1*TENS
!
!      HOC=(1-VNUI0)*SIG1-VNUI0*SIG3-VNUI1*TENS
!
!      HAA=(VNUI1+VNUI2*SPSI)*(SIG1-(TENS*(1+SPHI)-2*COHS)/ &
!          (1-SPHI))+ &
!          (VNUI1-SPSI)*(SIG2-TENS)+(VNUI1-SPSI)*(SIG3-TENS)
!
!      HAAB=VNUI0*(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))- &
!           (SIG2-TENS)+VNUI0*(SIG3-TENS)
!
!      HAAO=(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))- &
!           VNUI0*(SIG2-TENS)-VNUI0*(SIG3-TENS)
!
!      HBB=(VNUI1+SPSI)*(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))+ &
!          (VNUI1+SPSI)*(SIG2-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))+ &
!          (VNUI1-VNUI2*SPSI)*(SIG3-TENS)
!
!!        ***** DETERMINING RETURN AREA ****
!
!      IAREA=0
!      IASIGN=0
!
!      IF (F31 > 0 .AND. HA >= 0 .AND. HB < 0 .AND. HAB < 0) THEN
!        IAREA=2
!        IASIGN=IASIGN+1
!      ENDIF
!
!      IF (F31 > 0 .AND. HB >= 0) THEN
!        IASIGN=IASIGN+1
!        IF (HBB < 0) THEN
!          IAREA=1
!          GOTO 260
!        ELSE
!          IAREA=8
!        ENDIF
!      ENDIF
!
!      IF (F31 > 0 .AND. HA < 0) THEN
!        IASIGN=IASIGN+1
!        IF (HAA < 0) THEN
!
!          IAREA=3
!        ELSE
!          IAREA=7
!        ENDIF
!      ENDIF
!
!      IF (FT3 > 1d-12 .AND. &
!          HBA >= 0 .AND. &
!          HAO < 0 .AND. &
!          HOC < 0) THEN
!        IAREA=4
!        IASIGN=IASIGN+1
!      ENDIF
!
!      IF (FT3 > 1d-12 .AND. HAB >= 0 .AND. HBA < 0) THEN
!        IASIGN=IASIGN+1
!        IF (HAAB > 0) THEN
!          IAREA=5
!        ELSE
!          IF (IAREA == 0) IAREA=7
!        ENDIF
!      ENDIF
!
!      IF (FT3 > 1d-12 .AND. HAO >= 0) THEN
!        IASIGN=IASIGN+1
!        IF (HAAO >= 0) THEN
!          IAREA=6
!        ELSE
!          IF (IAREA == 0) IAREA=7
!        ENDIF
!      ENDIF
!
!260   CONTINUE
!
!      IF (IAREA == 0) THEN
!        !call GiveWarning('NO ACTIVE RETURNING AREA! INTGLO=' // trim(String(INTGLO)))
!      ENDIF
!
!
!      DSP1=0
!      DSP2=0
!      DSP3=0
!
!!
!!         *** EXTENSION POINT ***
!!
!      IF (IAREA == 1) THEN
!        A11 =       G *(1+SPHI*SPSI/VNUI1)
!        A12 = 0.5d0*G *(1+SPHI+SPSI+VNUIQ*SPHI*SPSI)
!        DETER  = A11*A11-A12*A12
!        RLAM31 = (F31*A11-F32*A12) / DETER
!        RLAM32 = (F32*A11-F31*A12) / DETER
!!        IF (RLAM31 < 0) call GiveWarning('RLAM31 < 0  INTGLO' // trim(String(INTGLO)))
!!        IF (RLAM32 < 0) call GiveWarning('RLAM32 < 0  INTGLO' // trim(String(INTGLO)))
!        RLAM21 = 0
!        DSP1 = RLAM31*PSIMIN + RLAM32*PSINU
!        DSP2 = RLAM31*PSINU  + RLAM32*PSIMIN
!        DSP3 = RLAM31*PSIMET + RLAM32*PSIMET
!      ENDIF
!!
!!
!!         *** REGULAR YIELD SURFACE ***
!!
!      IF (IAREA == 2) THEN
!        A11 = G *(1+SPHI*SPSI/VNUI1)
!        RLAM31 = F31 / A11
!!        IF (RLAM31 < 0) call GiveWarning('RLAM31 < 0  INTGLO' // trim(String(INTGLO)))
!        RLAM21 = 0
!        RLAM32 = 0
!        DSP1 = RLAM31*PSIMIN
!        DSP2 = RLAM31*PSINU
!        DSP3 = RLAM31*PSIMET
!      END IF
!!
!!         *** COMPRESSION POINT ***
!!
!      IF (IAREA == 3) THEN
!        A11 =       G *(1+SPHI*SPSI/VNUI1)
!        A12 = 0.5d0*G *(1-SPHI-SPSI+VNUIQ*SPHI*SPSI)
!        DETER  = A11*A11-A12*A12
!        RLAM31 = (F31*A11-F21*A12) / DETER
!        RLAM21 = (F21*A11-F31*A12) / DETER
!
!!        IF (RLAM31 < 0) call GiveWarning('RLAM31 < 0  INTGLO' // trim(String(INTGLO)))
!!        IF (RLAM21 < 0) call GiveWarning('RLAM21 < 0  INTGLO' // trim(String(INTGLO)))
!        RLAM32 = 0
!        DSP1 = RLAM31*PSIMIN + RLAM21*PSIMIN
!        DSP2 = RLAM31*PSINU  + RLAM21*PSIMET
!        DSP3 = RLAM31*PSIMET + RLAM21*PSINU
!      END IF
!
!      IF (IAREA == 4) THEN
!        DUM=   2*G /VNUI1
!        RLAMT3=FT3/(DUM*(1-VNUI0))
!        IF (RLAMT3 < 0) Then
!!          call WriteInLogFile('RLAMT3 < 0  INTGLO' // trim(String(INTGLO)))
!!          call WriteInLogFile('FT3,RL3 = ' // trim(String(Ft3)) //' '// trim(String(RLamT3)))
!!          call WriteInLogFile('G,vnui0 = ' // trim(String(G)) //' '// trim(String(vnui0)))
!        End If
!        DSP1 = DUM*RLAMT3*VNUI0
!        DSP2 = DSP1
!        DSP3 = DUM*RLAMT3*(1-VNUI0)
!        ITENS=1
!        IPL =2
!      ENDIF
!
!      IF (IAREA == 5) THEN
!        A11 = VNUI1+SPHI*SPSI
!        A12 = VNUI1+SPHI
!        A21 = VNUI1+SPSI
!        A22 = 2*(1-VNUI0)
!        DUM = G /VNUI1
!        DETER  = A11*A22-A12*A21
!
!        RLAM31 = (F31*A22-FT3*A12) / DETER / DUM
!        RLAMT3 = (FT3*A11-F31*A21) / DETER / DUM
!
!!        IF (RLAM31 < 0) call WriteInLogFile('RLAM31 < 0  INTGLO' // trim(String(INTGLO)))
!        IF (RLAMT3 < 0) Then
!!          call WriteInLogFile('RLAMT3 < 0  (IA=5) INTGLO' // trim(String(INTGLO)))
!!          call WriteInLogFile('F31 :' // trim(String(F31)))
!!          call WriteInLogFile('FT3 :' // trim(String(FT3)))
!!          call WriteInLogFile('A11 :' // trim(String(A11)))
!!          call WriteInLogFile('A12 :' // trim(String(A12)))
!!          call WriteInLogFile('A21 :' // trim(String(A21)))
!!          call WriteInLogFile('A22 :' // trim(String(A22)))
!!          call WriteInLogFile('HAB :' // trim(String(HAB)))
!!          call WriteInLogFile('HBA :' // trim(String(HBA)))
!        End If
!
!        DSP1 = DUM*(RLAM31*(-VNUI1+SPSI)+RLAMT3*2*VNUI0)
!        DSP2 = DUM*(RLAM31*2*VNUI0*SPSI+RLAMT3*2*VNUI0)
!        DSP3 = DUM*(RLAM31*(VNUI1+SPSI)+RLAMT3*2*(1-VNUI0))
!        ITENS=1
!        IPL =2
!      ENDIF
!
!      IF (IAREA == 6) THEN
!        RLAMT2 = (FT2*(1-VNUI0)-FT3*VNUI0)/(2*G)
!        RLAMT3 = (FT3*(1-VNUI0)-FT2*VNUI0)/(2*G)
!
!!        IF (RLAMT2 < 0) call WriteInLogFile('RLAMT2 < 0  INTGLO' // trim(String(INTGLO)))
!        IF (RLAMT3 < 0) Then
!!          call WriteInLogFile('RLAMT3 < 0  (IA=6) INTGLO' // trim(String(INTGLO)))
!!         call WriteInLogFile('FT2 :' // trim(String(FT2)))
!!          call WriteInLogFile('FT3 :' // trim(String(FT3)))
!!          call WriteInLogFile('RL2 :' // trim(String(RL2)))
!!          call WriteInLogFile('RL3 :' // trim(String(RL3)))
!!          call WriteInLogFile('A21 :' // trim(String(A21)))
!!          call WriteInLogFile('A22 :' // trim(String(A22)))
!!          call WriteInLogFile('HAo :' // trim(String(HAo)))
!!          call WriteInLogFile('HAAo :' // trim(String(HAAo)))
!        End If
!
!        DUM = 2*G/VNUI1
!        DSP1 = DUM*VNUI0*(RLAMT2+RLAMT3)
!        DSP2 = DUM*(RLAMT2*(1-VNUI0)+RLAMT3*   VNUI0 )
!        DSP3 = DUM*(RLAMT2*   VNUI0 +RLAMT3*(1-VNUI0))
!
!        ITENS=1
!        IPL =2
!      ENDIF
!
!      IF (IAREA == 7) THEN
!        DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
!        DSP2 = SIG2-TENS
!        DSP3 = SIG3-TENS
!        ITENS=1
!        IPL =2
!      ENDIF
!
!      IF (IAREA == 8) THEN
!        DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
!        DSP2 = SIG2-(SIG1-DSP1)
!        DSP3 = SIG3-TENS
!        ITENS=1
!        IPL=2
!      ENDIF
!
!      IF (IAREA == 9) THEN
!        DSP1 = SIG1-TENS
!        DSP2 = SIG2-TENS
!        DSP3 = SIG3-TENS
!        IAPEX=1
!        IPL=2
!      ENDIF
!
!!--------Check if the calculated principle stresses are on the yield surface
!      SIG1R=SIG1-DSP1
!      SIG2R=SIG2-DSP2
!      SIG3R=SIG3-DSP3
!
!      F21R = 0.5d0*(SIG2R-SIG1R) + 0.5d0*(SIG2R+SIG1R)*SPHI - COHS
!      F32R = 0.5d0*(SIG3R-SIG2R) + 0.5d0*(SIG3R+SIG2R)*SPHI - COHS
!      F31R = 0.5d0*(SIG3R-SIG1R) + 0.5d0*(SIG3R+SIG1R)*SPHI - COHS
!      FT1R = SIG1R - TENS
!      FT2R = SIG2R - TENS
!      FT3R = SIG3R - TENS
!
!      IF ( F21R > 1d-6.OR. &
!           F32R > 1d-6.OR. &
!           F31R > 1d-6.OR. &
!           FT1R > 1d-6.OR. &
!           FT2R > 1d-6.OR. &
!           FT3R > 1d-6    ) THEN
!        ICREC=0
!
!        IF (IAREA == 5.AND.F32R > 1E-6) THEN
!!------------ZONE 8
!          ICREC=ICREC+1
!          DSP1 = SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
!          DSP2 = SIG2-(SIG1-DSP1)
!          DSP3 = SIG3-TENS
!        ENDIF
!
!        IF (IAREA == 6.AND.FT1R > 1d-6) THEN
!!------------ZONE 9
!          ICREC=ICREC+1
!          DSP1 = SIG1-TENS
!          DSP2 = SIG2-TENS
!          DSP3 = SIG3-TENS
!          ITENS=0
!          IAPEX=1
!      ENDIF
!!        IF (ICREC == 0) call GiveWarning('No correction is done!')
!      ENDIF
!
!!--------Check again if the calculated stresses are on the yield surface
!      SIG1R=SIG1-DSP1
!      SIG2R=SIG2-DSP2
!      SIG3R=SIG3-DSP3
!
!      F21 = 0.5d0*(SIG2R-SIG1R) + 0.5d0*(SIG2R+SIG1R)*SPHI - COHS
!      F32 = 0.5d0*(SIG3R-SIG2R) + 0.5d0*(SIG3R+SIG2R)*SPHI - COHS
!      F31 = 0.5d0*(SIG3R-SIG1R) + 0.5d0*(SIG3R+SIG1R)*SPHI - COHS
!      FT1 = SIG1R - TENS
!      FT2 = SIG2R - TENS
!      FT3 = SIG3R - TENS
!
!!      IF (F21 > 0.001) call GiveWarning('F21=' // trim(String(F21))//' STILL >0')
!!      IF (F32 > 0.001) call GiveWarning('F32=' // trim(String(F32))//' STILL >0')
!!      IF (F31 > 0.001) call GiveWarning('F31=' // trim(String(F31))//' STILL >0')
!!      IF (FT1 > 0.001) call GiveWarning('FT1=' // trim(String(FT1))//' STILL >0')
!!      IF (FT2 > 0.001) call GiveWarning('FT2=' // trim(String(FT2))//' STILL >0')
!!      IF (FT3 > 0.001) call GiveWarning('FT3=' // trim(String(FT3))//' STILL >0')
!
!      IF (SIG2R > SIG3R+0.01 .OR. SIG2R < SIG1R-0.01) THEN
!!        call GiveWarning('ASSUMPTION NOT CORRECT!')
!!        call GiveWarning('The original returning zone IAREA=' // trim(String(IAREA)))
!      ENDIF
!
!!
!!         *** Computing Cartesian stress components ***
!!
!!      If (IntGlo == ip) then
!!        call WriteInLogFile('SigR123 : ' // trim(String(Sig1R)) //' '//trim(String(Sig2R)) //' '//trim(String(Sig3R)))
!!      endif
!
!  903       format(a,6(1x,1p,e19.12),1x,i4)
!
!      Call CarSig(Sig1R,Sig2R,Sig3R,xN1,xN2,xN3,NTENS,SigC)
!      STRESS=SigC
!
!
!
!!------- Calculate TAUMAX for checking accuracy on plastic point
!
!!      If (IntGlo == ip) then
!!        call WriteInLogFile('Sig123 : ' // trim(String(Sig1)) //' '//trim(String(Sig2)) //' '//trim(String(Sig3)))
!!      endif
!
!      TAUMAX=0.5d0*(SIG3-DSP3-SIG1+DSP1)
!      taumax1=taumax
!!     IF (TAUMAX <= -1E-6) call WriteInLogFile(' TAUMAX IS NEGATIVE!!!')
!      TAUMAX = MAX(TAUMAX,COHS,0.5d0)
!!      If (IntGlo == ip) then
!!        call WriteInLogFile('DSig123 : ' // trim(String(DSP1)) //' '//trim(String(DSP2)) //' '//trim(String(DSP3)))
!!        call WriteInLogFile('TauMax1 : ' // trim(String(TauMax1)))
!!        call WriteInLogFile('Taumax  : ' // trim(String(Taumax)))
!!      endif
!
!360   CONTINUE
!
!      RETURN
!
!      End subroutine  UMAT_AM
     
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

      End subroutine MatTranspose
        
        
        
        
      !  subroutine PrnSig(IOpt, ntens, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
      !!-------------------------------------------------------------------
      !!
      !!  Function: calculate principal stresses and directions
      !!            from cartesian stress vector
      !!
      !!  IOpt            I   I     flag to calculate principal direction (IOpt = 1)
      !!  IntGlo          I   I     global ID of Gauss point or particle 
      !!  S               I   R()   cartesian stress
      !!  xN1, xN2, xN3   O   R()   principal direction
      !!  S1, S2, S3      O   R     principal stress
      !!  P               O   R     isotropic stress (positive for tension)
      !!  Q               O   R     deviatoric stress
      !! 
      !!-------------------------------------------------------------------
      !
      !implicit none
      !
      !  ! arguments
      !  integer, intent(in) :: IOpt, ntens
      !  real(REAL_TYPE), intent(in), dimension(ntens) :: S ! size of stress tensor (6 for 3D, 4 for 2D)
      !  real(REAL_TYPE), intent(out), dimension(3) :: xN1, xN2, xN3
      !  real(REAL_TYPE), intent(out) :: S1, S2, S3, P, Q
      !
      !  if (IOpt == 1) then
      !    call Eig_3(0,ntens,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! Calculate principal direction
      !  else
      !    call Eig_3a(0,ntens,S,S1,S2,S3,P,Q) ! Do not calculate principal direction
      !  end if
      !
      !end subroutine PrnSig
        
        
        
        
!        subroutine Eig_3(iOpt, ntens, St, xN1, xN2, xN3, S1, S2, S3, P, Q)
!      !-------------------------------------------------------------------
!      !
!      !  Function: calculate principal stresses and directions 
!      !            from cartesian stress vector
!      !
!      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
!      ! 
!      !  IOpt            I   I     flag for output writing (IOpt = 1) 
!      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
!      !  xN1, xN2, xN3   O   R()   principal direction
!      !  S1, S2, S3      O   R     principal stress
!      !  P               O   R     isotropic stress (positive for tension)
!      !  Q               O   R     deviatoric stress
!      ! 
!      !-------------------------------------------------------------------
!
!      implicit none
!      
!        ! arguments
!        integer, intent(in) :: IOpt, ntens
!        real(REAL_TYPE), intent(in), dimension(ntens) :: St
!        real(REAL_TYPE), intent(out), dimension(3) :: xN1, xN2, xN3 ! is (3) for 2D and 3D
!        real(REAL_TYPE), intent(out) :: S1, S2, S3, P, Q
!        
!        ! local variables
!        real(REAL_TYPE), dimension(3, 3) :: A, V ! is (3,3) for 2D and 3D
!        real(REAL_TYPE) :: abs_max_s, tol, tau, sign_tau, t, c, s, temp1, temp2, temp3
!        integer :: i, k, it, itmax, ip, iq, iS1, iS2, iS3
!
!        ! Put cartesian stress vector into matrix A
!        select case(ntens)
!          case(4)
!              A = 0.0
!              A(1,1) = St(1) ! xx
!              A(1,2) = St(4) ! xy = yx
!              A(2,1) = St(4) ! xy = yx
!              A(2,2) = St(2) ! yy
!              A(3,3) = St(3) ! zz
!          
!              ! Set V to unity matrix
!              V = 0.0
!              V(1,1) = 1
!              V(2,2) = 1
!              V(3,3) = 1
!              
!          case(6)
!            A(1,1) = St(1) ! xx
!            A(1,2) = St(4) ! xy = yx
!            A(1,3) = St(6) ! zx = xz
!
!            A(2,1) = St(4) ! xy = yx
!            A(2,2) = St(2) ! yy
!            A(2,3) = St(5) ! zy = yz
!
!            A(3,1) = St(6) ! zx = xz
!            A(3,2) = St(5) ! zy = yz
!            A(3,3) = St(3) ! zz
!          
!           ! Set V to unity matrix
!            V(1,1) = 1
!            V(2,1) = 0
!            V(3,1) = 0
!
!            V(1,2) = 0
!            V(2,2) = 1
!            V(3,2) = 0
!
!            V(1,3) = 0
!            V(2,3) = 0
!            V(3,3) = 1
!        end select 
!        
!        ! get maximum value of cartesian stress vector
!        abs_max_s = 0.0
!        do i = 1, ntens
!          if (abs(St(i)) > abs_max_s) abs_max_s = abs(St(i))
!        end do
!      
!        ! set tolerance
!        tol = 1d-16 * abs_max_s
!        
!        ! get principal stresses and directions iteratively
!        it = 0
!        itmax = 50
!        do while ( (it < itmax) .and.  &
!                 (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) > tol) )
!        
!          it = it + 1
!          do k = 1,3
!                
!            if (k == 1) then
!              ip = 1
!              iq = 2
!            else if (k ==2) then
!              ip = 2
!              iq = 3
!            else
!              ip = 1
!              iq = 3
!            end if
!            
!            if (abs(A(ip,iq)) > 1d-50) then
!                
!              tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
!              if (tau >= 0.0) then
!                sign_tau = 1.0
!              else
!                sign_tau = -1.0
!              end if
!              
!              t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
!              c = 1.0 / sqrt(1.0 + t**2)
!              s = t * c
!              
!              temp1 = c * A(1, ip) - s * A(1, iq)
!              temp2 = c * A(2, ip) - s * A(2, iq)
!              temp3 = c * A(3, ip) - s * A(3, iq)
!              A(1, iq) = s * A(1, ip) + c * A(1, iq)
!              A(2, iq) = s * A(2, ip) + c * A(2, iq)
!              A(3, iq) = s * A(3, ip) + c * A(3, iq)
!              A(1, ip) = temp1
!              A(2, ip) = temp2
!              A(3, ip) = temp3
!
!              temp1 = c * V(1, ip) - s * V(1, iq)
!              temp2 = c * V(2, ip) - s * V(2, iq)
!              temp3 = c * V(3, ip) - s * V(3, iq)
!              V(1, iq) = s * V(1, ip) + c * V(1, iq)
!              V(2, iq) = s * V(2, ip) + c * V(2, iq)
!              V(3, iq) = s * V(3, ip) + c * V(3, iq)
!              V(1, ip) = temp1
!              V(2, ip) = temp2
!              V(3, ip) = temp3
!
!              temp1 = c * A(ip, 1) - s * A(iq, 1)
!              temp2 = c * A(ip, 2) - s * A(iq, 2)
!              temp3 = c * A(ip, 3) - s * A(iq, 3)
!              A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
!              A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
!              A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
!              A(ip, 1) = temp1
!              A(ip, 2) = temp2
!              A(ip, 3) = temp3
!            end if ! A(ip,iq)<>0
!            
!          end do ! k
!          
!          ! optional output writing
!          if (iOpt == 1) then
!!            call GiveMessage('error: ' // trim(String(abs(A(1, 2)) + abs(A(2, 3)) + abs(A(1, 3)))))
!          end if
!          
!        end do ! while
!     
!        ! get principal stresses from diagonal of A
!        S1 = A(1, 1)
!        S2 = A(2, 2)
!        S3 = A(3, 3)
!      
!        ! derived invariants
!        P = (S1 + S2 + S3) / 3.
!        Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )
!
!        ! Sort eigenvalues S1 <= S2 <= S3
!        iS1 = 1
!        iS2 = 2
!        iS3 = 3
!        
!        if (S1 > S2) then
!          t   = S2
!          S2  = S1
!          S1  = t
!          it  = iS2
!          iS2 = iS1
!          iS1 = it
!        end if
!        
!        if (S2 > S3) then
!          t   = S3
!          S3  = S2
!          S2  = t
!          it  = iS3
!          iS3 = iS2
!          iS2 = it
!        end if
!        
!        if (S1 > S2) then
!          t   = S2
!          S2  = S1
!          S1  = t
!          it  = iS2
!          iS2 = iS1
!          iS1 = it
!        end if
!        
!        ! get corresponding principal directions from V
!        do i = 1, 3 ! is 3 for 2D and 3D
!          xN1(i) = V(i, is1)
!          xN2(i) = V(i, is2)
!          xN3(i) = V(i, is3)
!        end do
!      
!      end subroutine Eig_3
!
!      
!      subroutine Eig_3a(iOpt, ntens, St, S1, S2, S3, P, Q)
!      !-------------------------------------------------------------------
!      !
!      !  Function: calculate principal stresses from cartesian stress vector
!      !
!      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
!      ! 
!      !  IOpt            I   I     flag for output writing (IOpt = 1) 
!      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
!      !  S1, S2, S3      O   R     principal stress
!      !  P               O   R     isotropic stress (positive for tension)
!      !  Q               O   R     deviatoric stress
!      ! 
!      !-------------------------------------------------------------------
!
!      
!      implicit none
!      
!        ! arguments
!        integer, intent(in) :: IOpt, ntens
!        real(REAL_TYPE), intent(in), dimension(ntens) :: St ! is (6) for 3D and (4) for 2D
!        real(REAL_TYPE), intent(out) :: S1, S2, S3, P, Q
!        
!        ! local variables
!        integer :: IDim
!        real(REAL_TYPE), dimension(:,:), allocatable :: A ! (3,3) for 2D and 3D
!        real(REAL_TYPE) :: abs_max_s, tol, tau, sign_tau, t, c, s, temp1, temp2, temp3
!        integer :: i, k, it, itmax, ip, iq
!        
!        IDim = 3
!        allocate(A(IDim,IDim)) 
!        
!        select case (ntens)       
!        case(4) ! for 2D
!          A = 0.0
!          A(1,1) = St(1) ! xx
!          A(1,2) = St(4) ! xy = yx
!          A(2,1) = St(4) ! xy = yx
!          A(2,2) = St(2) ! yy
!          A(3,3) = St(3) ! zz
!        case(6) ! for 3D
!          A(1,1) = St(1) ! xx
!          A(1,2) = St(4) ! xy = yx
!          A(1,3) = St(6) ! zx = xz
!
!          A(2,1) = St(4) ! xy = yx
!          A(2,2) = St(2) ! yy
!          A(2,3) = St(5) ! zy = yz
!
!          A(3,1) = St(6) ! zx = xz
!          A(3,2) = St(5) ! zy = yz
!          A(3,3) = St(3) ! zz
!        end select
!
!        ! get maximum value of cartesian stress vector
!        abs_max_s = 0.0
!        do i = 1,ntens
!          if (abs(St(i)) > abs_max_s) abs_max_s = abs(St(i))
!        end do
!      
!        ! set tolerance
!        tol = 1d-20 * abs_max_s
!        
!        ! get principal stresses and directions iteratively
!        it = 0
!        itmax = 50
!        do while ( (it < itmax) .and.  &
!                   (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) > tol) )
!        
!          it = it + 1
!          do k = 1,3
!            if (k == 1) then
!              ip = 1
!              iq = 2
!            else if (k ==2) then
!              ip = 2
!              iq = 3
!            else
!              ip = 1
!              iq = 3
!            end if
!            
!            if (abs(A(ip,iq)) > 1d-50) then
!                
!              tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
!              if (tau >= 0.0) then
!                sign_tau = 1.0
!              else
!                sign_tau = -1.0
!              end if
!              
!              t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
!              c = 1.0 / sqrt(1.0 + t**2)
!              s = t * c
!
!              temp1 = c * A(1, ip) - s * A(1, iq)
!              temp2 = c * A(2, ip) - s * A(2, iq)
!              temp3 = c * A(3, ip) - s * A(3, iq)
!              A(1, iq) = s * A(1, ip) + c * A(1, iq)
!              A(2, iq) = s * A(2, ip) + c * A(2, iq)
!              A(3, iq) = s * A(3, ip) + c * A(3, iq)
!              A(1, ip) = temp1
!              A(2, ip) = temp2
!              A(3, ip) = temp3
!
!              temp1 = c * A(ip, 1) - s * A(iq, 1)
!              temp2 = c * A(ip, 2) - s * A(iq, 2)
!              temp3 = c * A(ip, 3) - s * A(iq, 3)
!              A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
!              A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
!              A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
!              A(ip, 1) = temp1
!              A(ip, 2) = temp2
!              A(ip, 3) = temp3
!  
!            end if ! A(ip,iq)<>0
!            
!          end do ! k
!          
!          ! optional output writing
!          if (iOpt == 1) then
!!            call GiveMessage('error: ' // trim(String(abs(A(1, 2)) + abs(A(2, 3)) + abs(A(1, 3)))))
!          end if
!          
!        end do ! while
!     
!        ! get principal stresses from diagonal of A
!        S1 = A(1, 1)
!        S2 = A(2, 2)
!        S3 = A(3, 3)
!      
!        ! derived invariants
!        P = (S1 + S2 + S3) / 3.
!        Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )
!
!        ! Sort eigenvalues S1 <= S2 <= S3
!        if (S1 > S2) then
!          t   = S2
!          S2  = S1
!          S1  = t
!        end if
!        
!        if (S2 > S3) then
!          t   = S3
!          S3  = S2
!          S2  = t
!        end if
!        
!        if (S1 > S2) then
!          t   = S2
!          S2  = S1
!          S1  = t
!        end if
!        
!      end subroutine Eig_3a
        
        subroutine CarSig(S1, S2, S3, xN1, xN2, xN3, ntens, Stress)
      !----------------------------------------------------------
      !
      !  Function:  Returns the Cartesian stresses using 
      !             principal stress and principal direction
      !
      !  S1, S2, S3      I   R     principal stress
      !  xN1, xN2, xN3   I   R()   principal direction
      !  Stress          O   R()   cartesian stress
      !
      !----------------------------------------------------------
      
      implicit none
      
        ! arguments
        double precision, intent(in) :: S1, S2, S3
        double precision, intent(in), dimension(3) :: xN1, xN2, xN3
        integer, intent(in):: ntens
        double precision, intent(out), dimension(ntens) :: Stress
        
        ! local variables
        integer :: I
        integer :: IDim
        double precision, dimension(:, :), allocatable :: SM, T, TT, STT
        
        IDim = 3
        allocate(SM(IDim, IDim), T(IDim, IDim), TT(IDim, IDim), STT(IDim, IDim))

        ! put principal directions in matrix T
        do I = 1,3
            T(I,1) = xN1(I)
            T(I,2) = xN2(I)
            T(I,3) = xN3(I)
            TT(1,I) = T(I,1)
            TT(2,I) = T(I,2)
            TT(3,I) = T(I,3)            
        end do
      


        ! put principal stresses in matrix SM
        SM = 0.0
        SM(1,1) = S1
        SM(2,2) = S2
        SM(3,3) = S3

        ! calculate matrix multiplication T.S.TT
        call MatMat(SM, IDim, TT,  IDim, IDim, IDim, IDim , STT, IDim)
        call MatMat(T,  IDim, STT, IDim, IDim, IDim, IDim , SM,  IDim)
       
        ! Extract cartesian stress vector from stress matrix
        do I = 1, IDim
          Stress(I) = SM(I, I)
        end do        
        
        ! Extract cartesian stress vector from stress matrix
        Stress(4) = SM(2, 1)  
        if (ntens == 6) then
          Stress(5) = SM(3, 2)
          Stress(6) = SM(3, 1)
        end if
        
        end subroutine CarSig
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    

!!***********************************************************************
!      Subroutine MOHRStrainSoftening_PM4Sand(IntGlo,D1,D2,GG,cp,cr,phip,phir, &
!       psip,psir,factor,c,phi,psi,Sig0,DSigE,EpsP,DEps,DEpsP,SigC,IPL)
!      !**********************************************************************
!      !
!      ! Elastoplastic constitutive model with STRAIN SOFTENING, based on the 
!      ! MOHR-COULOMB criterion (considering modifications of Abbo & Sloan (1995))
!      ! Explicit MODIFIED EULER INTEGRATION SCHEME with automatic error control.
!      ! Final correction of the yield surface drift (END OF STEP CORRECTION).
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local variables
!      integer :: i,n,m,it
!      real(REAL_TYPE) :: F,F0,F2 !Evaluation of the Yield function
!      real(REAL_TYPE) :: alpha !Elastic Strain proportion
!      real(REAL_TYPE) :: SSTOL !Tolerance Relative Error
!      real(REAL_TYPE) :: YTOL !Tolerance Relative Error of the yield function evaluation
!      real(REAL_TYPE) :: SPTOL !Tolerance Softening parameters
!      real(REAL_TYPE) :: Rn !Relative error
!      real(REAL_TYPE) :: T,DT,T1,beta,DTmin !Substepping parameters
!      real(REAL_TYPE) :: c1,phi1,psi1,c2,phi2,psi2
!      real(REAL_TYPE) :: ctol,phitol,psitol !c,phi,psi tolerances
!      real(REAL_TYPE) :: Dcr,Dphir,Dpsir !Diference between current and residial values
!      real(REAL_TYPE) :: moduleEr,moduleSigDSig
!      real(REAL_TYPE) :: EpsPEq,EpsPEq1,EpsPEq2 !Equivalent Plastic Deformation
!      real(REAL_TYPE) :: DEpsPEq !Derivative Strain in function of Equivalent Plastic Deformation
!      real(REAL_TYPE) :: p,J,Lode,S3TA
!      real(REAL_TYPE), dimension(3) :: Principal_stresses
!      real(REAL_TYPE), dimension(3) :: SigC_TC
!      real(REAL_TYPE), dimension(3,3) :: Principal_vectors
!      real(REAL_TYPE), dimension(6) :: SigYield, SigYield2
!      real(REAL_TYPE), dimension(6) :: DSigPP,DSigP1,DSigP2
!      real(REAL_TYPE), dimension(6) :: DEpsPP,DEpsPP1,DEpsPP2
!      real(REAL_TYPE), dimension(6) :: DEpsS,DEpsSS
!      real(REAL_TYPE), dimension(6) :: EpsP1,EpsP2
!      real(REAL_TYPE), dimension(6) :: DEpsPEqDPS,DEpsPEqDPS1
!      real(REAL_TYPE), dimension(6) :: sumSg,Er
!      real(REAL_TYPE), dimension(3) :: DSPDPEq,DSPDPEq1 !Variation of softening parameters (c,phi,psi) in function of plastic strain
!      !In variables
!      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
!      real(REAL_TYPE), intent(in) :: D1,D2,GG !Elastic Parameters
!      real(REAL_TYPE), intent(in) :: cp,cr,phip,phir,psip,psir,factor !Softening parameter
!      real(REAL_TYPE), intent(in), dimension(6) :: Sig0 !Initial Stress
!      real(REAL_TYPE), intent(in), dimension(6) :: DEps !Incremental total strain
!      !Inout variables
!      real(REAL_TYPE), intent(inout):: c,phi,psi !cohesion,friction angle and dilatancy angle
!      real(REAL_TYPE), intent(inout), dimension(6) :: EpsP !Accumulated Plastic Strain
!      real(REAL_TYPE), intent(inout), dimension(6) :: SigC !Final Stress
!      real(REAL_TYPE), intent(inout), dimension(6) :: DSigE !Incremental Elastic Stress
!      !Out variables
!      integer, intent(out) :: IPL
!      real(REAL_TYPE), intent(out), dimension(6) :: DEpsP !Incremental plastic strain
!
!      !Initialization
!      DEpsPEq = 0.0d0
!      EpsPEq = 0.0d0
!      SigYield = 0.0d0
!      DEpsP = 0.0d0
!      F = 0.0d0
!      it = 0
!
!      if (c > cp.or.phi > phip.or.psi > psip) then
!        c = min(c,cp)
!        phi = min(phi,phip)
!        psi = min(psi,psip)
!      end if
!      if (c < cr.or.phi < phir.or.psi < psir) then
!        c = max(c,cr)
!        phi = max(phi,phir)
!        psi = max(psi,psir)
!      end if
!
!      !Tolerances
!      SSTOL = 0.01d0 !Tolerance Relative Error (10-3 to 10-5)
!      YTOL = 0.0001d0 !Tolerance Error on the Yield surface (10-6 to 10-9)
!      SPTOL = 0.01d0 !Tolerance Softening Parameters (0.0001d0)
!      ctol = abs(cp-cr)*SPTOL
!      phitol = abs(phip-phir)*SPTOL
!      psitol = abs(psip-psir)*SPTOL
!      DTmin = 0.000000001d0
!      
!      !Check the yield function value
!      call DetermineYieldFunctionValue(IntGlo,SigC,c,phi,F)
!      
!      
!      !If F<0 then the behaviour is elastic --> Return
!      if (F <= YTOL) then
!        IPL = 0
!      !  ! Tension cutoff implementation 
!      !! First -> compute mean effective stress 'p' and deviatoric stress 'J' for SigC
!      !call CalculateInvariants(IntGlo,SigC,p,J,Lode,S3TA)
!      !if (3*p > 0) then !hardcoding T as zero
!      !    call Get_EigenValues_EigenVectors(SigC, Principal_stresses, Principal_vectors)          
!      !    SigC = 0.0
!      !    !SigC_TC(3) = Principal_stresses(3) + (0 - (Principal_stresses(1)+Principal_stresses(2)+Principal_stresses(3)) )/3
!      !    !SigC_TC(2) = Principal_stresses(2) + (0 - (Principal_stresses(1)+Principal_stresses(2)+Principal_stresses(3)) )/3
!      !    !SigC_TC(1) = Principal_stresses(1) + (0 - (Principal_stresses(1)+Principal_stresses(2)+Principal_stresses(3)) )/3
!      !    !SigC = [Sig_TC(1), Sig_TC(2), Sig_TC(3), 0, 0, 0]
!      !end if
!        return
!      end if
!
!      !If F>0, the behaviour is elastoplastic --> Continue
!      Dcr = abs(c - cr)
!      Dphir = abs(phi - phir)
!      Dpsir = abs(psi - psir)
!      !Check if we are in residual conditions or in softening conditions
!      if (Dcr <= ctol.and.Dphir <= phitol.and.Dpsir <= psitol) then
!        IPL = 1 !IPL=1 Residual conditions --> no changes of the strength parameters
!        c = cr
!        phi = phir
!        psi = psir
!      else
!        IPL = 2 !IPL=2 Softening Conditions --> changes of the strength parameters
!      end if
!
!      !Determine the proportion (alpha) of the stress increment that lies within the yield function.
!      !The PEGASUS ALGORITHM SCHEME FOR CONVENTIONAL ELASTOPLASTIC MODELS has been used
!      call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)
!
!      if (F0 < -YTOL) then !In this Time increment there is part of elastic behavior
!        call DetermineElasticProportionPegasusMethod(IntGlo,Sig0,DSigE,DEps,c,phi,YTOL,alpha)
!      else  
!        alpha = 0.0d0 !In this time increment all is plastic behavior
!      end if
!
!      !Calculate the direction of the stress path--> missing
!      !It is assumed that the direction is always outside the yield surface.
!
!      !Determine the elastic portion of the stress increment
!      DSigE = alpha * DSigE !Correct Incremental Elastic Stress
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !Determine the plastic portion of the stress increment.
!      !The method used is the MODIFIED EULER INTEGRATION SCHEME with error control
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      !Initialise parameters
!      SigYield = Sig0 + DSigE !Sigma on the Yield surface
!      DEpsS = (1.0d0-alpha) * DEps !Incremental Plastic Strain
!
!      T = 0.0d0
!      DT = 1.0d0
!
!      !Start the plastification
!      Do while (T <= 1.0d0)
!        m = 0 !Counter
!        Rn = 100
!
!        call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)
!
!        Do while (Rn > SSTOL.and.m < 1000)
!            !1)Calculation of the portion of the plastic strain increment (DEpsPP)
!            DEpsSS = DT * DEpsS !Portion of the plastic strain increment
!
!            !Calculate a first estimate of the associated stress 
!            !hardening/softening parameter changes
!            call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
!                EpsPEq,DSPDPEq)
!            call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
!            call DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c,phi,psi,SigYield,DEpsPEqDPS,DSPDPEq,DEpsSS,DSigP1,DEpsPP1)
!            EpsP1 = EpsP + DEpsPP1
!
!            call CalculateEpsPEq(EpsP1,EpsPEq1) !Determine Equivalent plastic Strain (EpsPEq)
!
!            !if (IPL == 1) then !Residual conditions --> no changes of the strength parameters
!            !    c1 = c
!            !    phi1 = phi
!            !    psi1 = psi
!            !else !IPL=2 Softening Conditions --> changes of the strength parameters
!                call CalculateSofteningParameters(EpsPEq1,factor,cp,cr,phip,phir,psip,psir,c1,phi1,psi1)
!            !end if
!
!            !2)Calculate a second estimate of the associated stress 
!            !hardening/softening parameter changes
!            SigYield2 = SigYield + DSigP1
!
!            call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
!                EpsPEq1,DSPDPEq1)
!            call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP1,EpsPEq1,DEpsPEqDPS1)
!            call DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c1,phi1,psi1,SigYield2,DEpsPEqDPS1,DSPDPEq1,DEpsSS,DSigP2,DEpsPP2)
!            EpsP2 = EpsP + DEpsPP2
!
!            call CalculateEpsPEq(EpsP2,EpsPEq2) !Determine Equivalent plastic Strain (EpsPEq)
!
!            !if (IPL == 1) then !Residual conditions --> no changes of the strength parameters
!            !    c2 = c
!            !    phi2 = phi
!            !    psi2 = psi
!            !else  !IPL=2 Softening Conditions --> changes of the strength parameters
!                call CalculateSofteningParameters(EpsPEq2,factor,cp,cr,phip,phir,psip,psir,c2,phi2,psi2)
!            !end if
!
!            !3)Obtain a more accurate modified Euler estimate of the changes in stress,
!            !plastic strain and hardening/softening parameters
!            DSigPP = 0.5d0 * (DSigP1 + DSigP2)
!
!            !Calculation of the relative error
!            Er = 0.5d0 * (DSigP1 - DSigP2)
!            moduleEr = sqrt(Er(1)*Er(1)+Er(2)*Er(2)+ Er(3)*Er(3)+ Er(4)*Er(4)+Er(5)*Er(5)+Er(6)*Er(6))
!
!            sumSg = SigYield + DSigPP
!            moduleSigDSig = sqrt(sumSg(1)*sumSg(1) + sumSg(2)*sumSg(2) + sumSg(3)*sumSg(3)+ &
!                                 sumSg(4)*sumSg(4) + sumSg(5)*sumSg(5) + sumSg(6)*sumSg(6))
!
!            !Check the relative error (Rn) of the new stresses, with the defined tolerance (SSTOL)
!            Rn = (moduleEr/moduleSigDSig)
!
!            ! check whether decreasing of DT is possible, if not exit loop
!            if (DT == DTmin) then
!               exit
!            end if
!          
!            !4)If Rn>SSTOL, the loop is not finished and the substep is recalculated smaller
!            if (Rn > SSTOL) then
!                beta = max (0.9d0*(sqrt(SSTOL/Rn)), 0.1d0)
!                beta = min (beta, 1.1d0)
!                DT = max (DT*beta, DTmin)
!                m = m + 1 !Update counter
!            end if
!
!        end do
!
!        !Update the accumulated stresses, plastic strain and softening parameters
!        SigYield = SigYield + DSigPP
!        DEpsPP = 0.5d0 * (DEpsPP1 + DEpsPP2)
!        DEpsP = DEpsP + DEpsPP
!        EpsP = EpsP + DEpsPP
!
!        call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)
!
!        call CalculateSofteningParameters(EpsPEq,factor,cp,cr,phip,phir,psip,psir,c,phi,psi)
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !!!!!!!!!!!!!!!!!!!!!!! END OF STEP CORRECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !Check if we are on/in the yield surface, otherwise we are still outside (F>0)
!        !and a correction is needed.
!        call DetermineYieldFunctionValue(IntGlo,SigYield,c,phi,F)
!        n=0 !Counter
!        do while (abs(F) > YTOL.and.n < 10) !The correction is needed
!            n = n + 1
!            call CalculateEpsPEq(EpsP,EpsPEq)             !Determine Equivalent plastic Strain (EpsPEq)
!            call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir, &
!                EpsPEq,DSPDPEq)
!            call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
!            call EndOfStepCorrection(IntGlo,D1,D2,GG,IPL,F,SigYield,DSPDPEq,DEpsPEqDPS,EpsP,c,phi,psi)
!        end do
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        !The substep is updated
!        T1 = T + DT
!
!        !If T1>1 the calculation is finished
!        If (T1 >= 1d0) then
!            SigC = SigYield   !Determine Final stresses
!            return
!        end if
!
!        !If T1<1, calculation of the next substep DT
!        beta = min (0.9d0*(sqrt(SSTOL/Rn)), 1.1d0)
!        if (m > 1) then ! the previous step failed
!            beta = min (beta, 1.0d0)
!            DT = beta * DT
!            it = it+1
!        else
!            DT = beta * DT
!            it = 0
!        end if
!        DT = max (DT, DTmin)
!        DT = min (DT, 1.0d0-T1)
!        T = T1
!        
!      end do  !If T=1 the loop is finished
!      
!      ! Tension cutoff implementation 
!      ! First -> compute mean effective stress 'p' and deviatoric stress 'J' for SigC
!      call CalculateInvariants(IntGlo,SigC,p,J,Lode,S3TA)
!      if (3*p > 0) then !hardcoding T as zero
!          call Get_EigenValues_EigenVectors(SigC, Principal_stresses, Principal_vectors)  
!          SigC = 0.0
!          !SigC_TC = 0.0
!          !SigC_TC(3) = SigC(3) + (0 - p)
!          !SigC_TC(2) = SigC(2) + (0 - p)
!          !SigC_TC(1) = SigC(1) + (0 - p)
!      end if
!      
!      
!
!      end subroutine MOHRStrainSoftening_PM4Sand
! 
!
!      Subroutine DetermineElasticProportionPegasusMethod_PM4Sand(IntGlo,Sig,DSig,DEps,c,phi,YTOL,alpha)
!      !**********************************************************************
!      !
!      ! The PEGASUS METHOD method is used  
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local variables
!      integer :: Its     
!      real(REAL_TYPE) :: alpha0,alpha1,F0,F1,F
!      real(REAL_TYPE), dimension(6) :: Sig0,Sig1,SigNew
!      !In variables
!      real(REAL_TYPE), intent(in), dimension(6) :: Sig, DSig
!      real(REAL_TYPE), intent(in), dimension(6) :: DEps
!      real(REAL_TYPE), intent(in) :: c,phi,YTOL
!      integer, intent(in) :: IntGlo       !Global ID of Gauss point or particle
!      !Out variables
!      real(REAL_TYPE), intent(out) :: alpha
!
!      alpha0 = 0.0d0
!      alpha1 = 1.0d0
!
!      Sig0 = Sig + alpha0*DSig ! = Sig0
!      Sig1 = Sig + alpha1*DSig ! = SigE
!        
!      call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)
!      call DetermineYieldFunctionValue(IntGlo,Sig1,c,phi,F1)
!
!      F=YTOL+1000
!      Its = 0 !Counter
!
!      do while (abs(F) > YTOL.and.Its < 1000)
!        alpha = alpha1 - F1*(alpha1-alpha0)/(F1-F0)
!        SigNew = Sig + alpha*DSig
!       
!        call DetermineYieldFunctionValue(IntGlo,SigNew,c,phi,F)
!
!        if ((F*F1) < 0.0d0) then
!            alpha0 = alpha1
!            F0 = F1
!        else
!            F0 = F1*F0/(F1+F)
!        end if
!
!        alpha1 = alpha
!        F1 = F
!        Its = Its + 1
!
!      end do  
!      if (Its >= 1000) then
!            alpha = 0.0d0
!      end if
!
!      end subroutine DetermineElasticProportionPegasusMethod_PM4Sand
!
!
!      Subroutine CalculateInvariants_PM4Sand(IntGlo,Sig,p,J,Lode,S3TA)
!      !**********************************************************************
!      !
!      ! Calcuation of the invariants (defined as Abbo & Sloan (1995))
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local variables
!      real(REAL_TYPE) :: Sx,Sy,Sz,SqTxy,SqTyz,SqTxz,suma,h1,h2,J2,J3
!      real(REAL_TYPE), parameter :: C00000 = 0.0D0
!      real(REAL_TYPE), parameter :: C00001 = 1.0D0
!      real(REAL_TYPE), parameter :: C00P16 = 0.166666666666666D0
!      real(REAL_TYPE), parameter :: C00002 = 2.0D0
!      real(REAL_TYPE), parameter :: C00003 = 3.0D0
!      real(REAL_TYPE), parameter :: CP3333 = 0.333333333333333D0
!      real(REAL_TYPE), parameter :: C00IR3 = 0.577350269189626D0
!      real(REAL_TYPE), parameter :: TINY = 0.000000000000001D0
!      !In variables
!      real(REAL_TYPE), intent(in), dimension(6) :: Sig
!      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
!      !Out variables
!      real(REAL_TYPE), intent(out) :: p,J,Lode,S3TA !Invariants
!
!      p = C00000
!      J = C00000
!      Lode = C00000
!
!      !Mean stress (p)
!      p = CP3333 * (Sig(1) + Sig(2) + Sig(3))
!
!      !Deviatoric stress (J)
!      Sx = Sig(1) - p
!      Sy = Sig(2) - p
!      Sz = Sig(3) - p
!      suma = (Sig(1)-Sig(2))*(Sig(1)-Sig(2))+(Sig(1)-Sig(3))*(Sig(1)-Sig(3))+(Sig(2)-Sig(3))*(Sig(2)-Sig(3))
!      SqTxy =  Sig(4) * Sig(4)
!      SqTyz =  Sig(5) * Sig(5)
!      SqTxz =  Sig(6) * Sig(6)
!
!      J2 = C00P16 * suma + SqTxy + SqTyz + SqTxz
!      J3 = Sx*Sy*Sz + C00002 * Sig(4)*Sig(5)*Sig(6) - Sx*SqTyz - Sy*SqTxz - Sz*SqTxy
!      J = SQRT(J2)
!
!      !Lode's angle (Lode)
!      if (J2 > 0.0001) then !C00000
!
!        h1 = -C00003/(C00002*C00IR3)
!        h2 = J3/(J*J*J)
!        S3TA = h1*h2
!        if (S3TA < -C00001) then
!            S3TA = -C00001
!        else if (S3TA > C00001) then
!            S3TA = C00001
!      end if
!        Lode = CP3333*asin(S3TA)
!      else  !Special case of zero deviatoric stress
!        Lode = C00000
!        S3TA = C00000
!      end if
!
!      end subroutine CalculateInvariants_PM4Sand
!
! 
!      Subroutine DetermineYieldFunctionValue_PM4Sand(IntGlo,Sig,c,phi,F)
!      !**********************************************************************
!      !
!      ! In this subroutine the yield function evaluated is a smooth hyperbolic approximation to the
!      ! Mohr-Coulomb yield criterion (Abbo and Sloan, 1995).
!      !
!      ! The edges of the hexagonal pyramid and the tip have been smoothed.
!      ! There are two parameters aSmooth (smoothes the tip) and ATTRAN(smoothes the edges)
!      ! In this case aSmooth=0.0005*c*cot(phi) and LodeT=25º.
!      ! If aSmooth=0 and LodeT=30º the original Mohr-Coulomb is obtained.
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local variables
!      real(REAL_TYPE) ::  p,J,Lode,S3TA !Invariants
!      real(REAL_TYPE) ::  COH, SPHI, CPHI, COTPHI, STA, CTA, K, aSmooth, ASPHI2, SGN, A, B
!      real(REAL_TYPE), parameter :: C00001 = 1.0d0 !Parameters
!      real(REAL_TYPE), parameter :: C00003 = 3.0d0
!      real(REAL_TYPE), parameter :: C00P50 = 0.0005d0
!      real(REAL_TYPE), parameter :: C00000 = 0.0d0
!      real(REAL_TYPE), parameter :: C00IR3 = 0.577350269189626d0
!      real(REAL_TYPE), parameter :: C000P1 = 0.00000000001D0
!      !Constants for rounded K function (for LodeT=25)
!      !real(REAL_TYPE), parameter :: A1 = 1.432052062044227d0
!      !real(REAL_TYPE), parameter :: A2 = 0.406941858374615d0
!      !real(REAL_TYPE), parameter :: B1 = 0.544290524902313d0
!      !real(REAL_TYPE), parameter :: B2 = 0.673903324498392d0
!      !real(REAL_TYPE), parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
!      !Constants for rounded K function (for LodeT=29.5)
!      real(REAL_TYPE), parameter :: A1 = 7.138654723242414d0
!      real(REAL_TYPE), parameter :: A2 = 6.112267270920612d0
!      real(REAL_TYPE), parameter :: B1 = 6.270447753139589d0
!      real(REAL_TYPE), parameter :: B2 = 6.398760841429403d0
!      real(REAL_TYPE), parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
!      !Constants for rounded K function (for LodeT=30)
!      !real(REAL_TYPE), parameter :: A1 = -138300705.446275
!      !real(REAL_TYPE), parameter :: A2 = -138300706.472675
!      !real(REAL_TYPE), parameter :: B1 = -138300706.3123
!      !real(REAL_TYPE), parameter :: B2 = 0.192450089729875
!      !real(REAL_TYPE), parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
!      !In variables
!      real(REAL_TYPE), intent(in), dimension(6) :: Sig
!      real(REAL_TYPE), intent(in) :: c,phi
!      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
!
!      !Out variables
!      real(REAL_TYPE), intent(out) :: F
!
!      F = C00000
!
!      !Calculation of the invariants (p',J,Lode)
!      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !!Evaluation of the yield function with Smoothing!!!!!!!!!!!!!!!!!!!!!!!!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !Material parameters
!      COH = c     !Cohesion
!      SPHI = sin(phi) 
!      CPHI = cos(phi)
!      COTPHI = CPHI/SPHI
!      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
!      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
!      if (abs(phi) == C00000) then
!            ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
!      end if
!
!      !Calculate K function
!      if (abs(Lode) < ATTRAN) then
!        STA = sin(Lode)
!        CTA = cos(Lode)
!        K = CTA - STA*SPHI*C00IR3
!      else
!        SGN = SIGN(C00001,Lode)
!        A = A1 + A2*SGN*SPHI
!        B = B1*SGN + B2*SPHI
!        K = A - B*S3TA
!      end if
!
!      !Calculate value of Hyperbolic Yield function
!      F = p*SPHI + sqrt(J*J*K*K+ASPHI2) - COH*CPHI
!      
!      end subroutine DetermineYieldFunctionValue_PM4Sand
!
!
!      Subroutine CalculateDerivativesYieldFunctAndPlasticPotential_PM4Sand(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
!      !**********************************************************************
!      !
!      ! Calculation of the derivatives of the yield function (F) and the plastic potencial punction (P).
!      ! Based on Abbo & Sloan (1995)
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local variables
!      integer :: i
!      real(REAL_TYPE) :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B,&
!                           D, aSmooth, ASPHI2, SGN, T3TA, C3TA, J2, psi2
!      real(REAL_TYPE) ::   K, dKdLode
!      real(REAL_TYPE) :: SPSI, CPSI, TPSI, COTPSI, ASPSI2
!      real(REAL_TYPE) :: i1, i2, Sx, Sy, Sz
!      real(REAL_TYPE) :: DFDp,DFDJ,DFDLode !Derivatives F respect Invariants
!      real(REAL_TYPE) :: DPDp,DPDJ,DPDLode !Derivatives P respect Invariants
!      real(REAL_TYPE) :: C1, C2, C3
!      real(REAL_TYPE), dimension(6):: DpDSig,DJDSig,DJ3DSig !Derivatives Invariants
!
!      real(REAL_TYPE), parameter :: C00001 = 1.0D0 !Parameters
!      real(REAL_TYPE), parameter :: C000P5 = 0.5D0
!      real(REAL_TYPE), parameter :: C00P50 = 0.0005D0
!      real(REAL_TYPE), parameter :: C00000 = 0.0D0
!      real(REAL_TYPE), parameter :: C00003 = 3.0D0
!      real(REAL_TYPE), parameter :: C00004 = 4.0D0
!      real(REAL_TYPE), parameter :: C00002 = 2.0D0
!      real(REAL_TYPE), parameter :: CP3333 = 0.333333333333333D0
!      real(REAL_TYPE), parameter :: C00IR3 = 0.577350269189626D0
!      real(REAL_TYPE), parameter :: C0R3I2 = 0.866025403784439D0
!      real(REAL_TYPE), parameter :: C000P1 = 0.000000000000001D0 
!      real(REAL_TYPE), parameter :: J0 = 0.001D0 
!      !Constants for rounded K function (for LodeT=25)
!      !real(REAL_TYPE), parameter :: A1 = 1.432052062044227d0
!      !real(REAL_TYPE), parameter :: A2 = 0.406941858374615d0
!      !real(REAL_TYPE), parameter :: B1 = 0.544290524902313d0
!      !real(REAL_TYPE), parameter :: B2 = 0.673903324498392d0
!      !real(REAL_TYPE), parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
!      !Constants for rounded K function (for LodeT=29.5)
!      real(REAL_TYPE), parameter :: A1 = 7.138654723242414d0
!      real(REAL_TYPE), parameter :: A2 = 6.112267270920612d0
!      real(REAL_TYPE), parameter :: B1 = 6.270447753139589d0
!      real(REAL_TYPE), parameter :: B2 = 6.398760841429403d0
!      real(REAL_TYPE), parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
!      !Constants for rounded K function (for LodeT=30)
!      !real(REAL_TYPE), parameter :: A1 = -138300705.446275
!      !real(REAL_TYPE), parameter :: A2 = -138300706.472675
!      !real(REAL_TYPE), parameter :: B1 = -138300706.3123
!      !real(REAL_TYPE), parameter :: B2 = 0.192450089729875
!      !real(REAL_TYPE), parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
!      !In variables
!      real(REAL_TYPE), intent(in) ::  c,phi,psi !Soft Parameters
!      real(REAL_TYPE), intent(in), dimension(6) :: Sig
!      !Out variables
!      real(REAL_TYPE), intent(out), dimension(6) :: DFDSig, DPPDSig !Derivatives respect Sigma
!      !Inout variables
!      real(REAL_TYPE), intent(inout) :: p,J,Lode,S3TA !Invariants
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !!!!! DFDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      !Material parameters
!      COH = c !Cohesion
!      SPHI = sin(phi) 
!      CPHI = cos(phi)
!      COTPHI = CPHI/SPHI
!      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
!      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
!      if (abs(phi) == C00000) then
!        ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
!      end if
!
!      if (J == C00000) then
!        J2 = C000P1
!        J = sqrt(J2)
!      else
!        J2 = J*J
!      end if
!
!      CTA = cos(Lode)
!      C3TA = CTA*(C00004*CTA*CTA-C00003)
!      T3TA = S3TA/C3TA
!
!      !Calculate K function and its derivative
!      if (abs(Lode) < ATTRAN) then
!        STA = S3TA/(C00004*CTA*CTA-C00001)
!        K = CTA - STA*SPHI*C00IR3
!        dKdLode =  - STA - C00IR3*SPHI*CTA
!      else
!        SGN = SIGN(C00001,Lode) ! It puts the Lode's sign to the number 1
!        A = A1 + A2*SGN*SPHI
!        B = B1*SGN + B2*SPHI
!        K = A - B*S3TA
!        dKdLode = - C00003*B*C3TA
!      end if
!      
!      !Derivative Dp/DSig
!      DpDSig(1) = CP3333
!      DpDSig(2) = CP3333
!      DpDSig(3) = CP3333
!      DpDSig(4) = C00000
!      DpDSig(5) = C00000
!      DpDSig(6) = C00000
!      
!      !Derivative DJ/DSig
!      i1 = C000P5/J
!      if (J < 0.0001) then
!        i1 = 0.0d0
!      end if
!      Sx = Sig(1)-p
!      Sy = Sig(2)-p
!      Sz = Sig(3)-p
!
!      DJDSig(1) = i1 * Sx
!      DJDSig(2) = i1 * Sy
!      DJDSig(3) = i1 * Sz
!      DJDSig(4) = i1 * C00002 * Sig(4)
!      DJDSig(5) = i1 * C00002 * Sig(5)
!      DJDSig(6) = i1 * C00002 * Sig(6)
!
!      !Derivative DJ3/DSig
!      i2 = CP3333*J*J
!      DJ3DSig(1) = (Sy*Sz - Sig(5)*Sig(5) + i2)
!      DJ3DSig(2) = (Sx*Sz - Sig(6)*Sig(6) + i2)
!      DJ3DSig(3) = (Sx*Sy - Sig(4)*Sig(4) + i2)
!      DJ3DSig(4) = C00002*(Sig(5)*Sig(6) - Sz*Sig(4))
!      DJ3DSig(5) = C00002*(Sig(6)*Sig(4) - Sx*Sig(5))
!      DJ3DSig(6) = C00002*(Sig(4)*Sig(5) - Sy*Sig(6))
!
!      D = J*K/(sqrt(J2*K*K + ASPHI2))
!
!      !C1F
!      C1 = SPHI
!      !C2F
!      C2 = D*K - T3TA*D*dKdLode
!      !C3F
!      C3 = -C0R3I2*dKdLode*D/(J2*C3TA)
!      
!      !DFDSig!
!      do i=1,6
!        DFDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
!      end do
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !!!!! DPPDSig = DFDSig (if associated Flow Rule)  !!!!!!!!!!!!!!!!!!!!!!
!      !!!!! or
!      !!!!! DPPDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      if (abs(J) < J0) then
!        psi2 = phi - abs(J)*(phi - psi)/J0
!      else
!        psi2 = psi
!      end if
!
!      if (phi == psi2) then !If Associated Flow Rule, then DPPDSig = DFDSig
!        DPPDSig = DFDSig
!
!      else !If Non-Associated Flow Rule, then calculate...
!        !Material parameters
!          SPSI = sin(psi2) 
!        CPSI = cos(psi2)
!        if (abs(psi2) == C00000)then
!            ASPSI2 = C00000
!        
!        else
!            
!        COTPSI = CPSI/SPSI
!        aSmooth = C00P50*COH*COTPSI !Smoothing parameter
!        ASPSI2 = aSmooth*aSmooth*SPSI*SPSI
!        
!        end if
!        
!        
!
!        !Calculate K function and its derivative
!        if (abs(Lode) <= ATTRAN) then
!            K = CTA - STA*SPSI*C00IR3
!            dKdLode = - STA - C00IR3*SPSI*CTA
!        else
!            A = A1 + A2*SGN*SPSI
!            B = B1*SGN + B2*SPSI
!            K = A - B*S3TA
!            dKdLode = - C00003*B*C3TA
!        end if
!
!        D = J*K/(sqrt(J*J*K*K + ASPSI2))
!
!        !C1F
!        C1 = SPSI
!        !C2F
!        C2 = D*K - T3TA*D*dKdLode
!        !C3F
!        C3 = -C0R3I2*dKdLode*D/(J2*C3TA)
!
!        !DPPDSig
!        do i=1,6
!            DPPDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
!        end do
!
!      end if
!
!      end subroutine CalculateDerivativesYieldFunctAndPlasticPotential_PM4Sand
!
!
!      Subroutine CalculateDerivativesYieldFunctSofteningParameters_PM4Sand(p,J,Lode,S3TA,c,phi,DFDSP)
!      !**********************************************************************
!      !
!      ! Calculation of the derivatives of the yield function (F) with respect the strength parameters
!      ! The strength parameters are: cohesion (COH) and friction angle (PHI)
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local variables
!      real(REAL_TYPE) :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B,&
!                         Denom, Num, aSmooth, ASPHI2, SGN
!      real(REAL_TYPE) :: K, dKdPhi, dadc, dadPhi
!      real(REAL_TYPE), parameter :: C00001 = 1.0D0 !Parameters
!      real(REAL_TYPE), parameter :: C00P50 = 0.0005D0
!      real(REAL_TYPE), parameter :: C00000 = 0.0D0
!      real(REAL_TYPE), parameter :: C00003 = 3.0D0
!      real(REAL_TYPE), parameter :: C00002 = 2.0D0
!      real(REAL_TYPE), parameter :: C00IR3 = 0.577350269189626D0
!      real(REAL_TYPE), parameter :: C000P1 = 0.00000000001D0
!      !Constants for rounded K function (for LodeT=25)
!      !real(REAL_TYPE), parameter :: A1 = 1.432052062044227d0
!      !real(REAL_TYPE), parameter :: A2 = 0.406941858374615d0
!      !real(REAL_TYPE), parameter :: B1 = 0.544290524902313d0
!      !real(REAL_TYPE), parameter :: B2 = 0.673903324498392d0
!      !real(REAL_TYPE), parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
!      !Constants for rounded K function (for LodeT=29.5)
!      real(REAL_TYPE), parameter :: A1 = 7.138654723242414d0
!      real(REAL_TYPE), parameter :: A2 = 6.112267270920612d0
!      real(REAL_TYPE), parameter :: B1 = 6.270447753139589d0
!      real(REAL_TYPE), parameter :: B2 = 6.398760841429403d0
!      real(REAL_TYPE), parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
!      !Constants for rounded K function (for LodeT=30)
!      !real(REAL_TYPE), parameter :: A1 = -138300705.446275
!      !real(REAL_TYPE), parameter :: A2 = -138300706.472675
!      !real(REAL_TYPE), parameter :: B1 = -138300706.3123
!      !real(REAL_TYPE), parameter :: B2 = 0.192450089729875
!      !real(REAL_TYPE), parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
!
!      !In variables
!      real(REAL_TYPE), intent(in) :: p,J,Lode,S3TA !Invariants
!      real(REAL_TYPE), intent(in) :: c,phi !Soft Parameters
!      !Out variables
!      real(REAL_TYPE), intent(out), dimension(2) :: DFDSP !Derivatives respect Soft Parameters
!
!
!      !Material parameters
!      COH = c !Cohesion
!      SPHI = sin(phi) 
!      CPHI = cos(phi)
!      COTPHI = CPHI/SPHI
!
!      !Calculate aSmooth and its derivatives
!      if (abs(phi) == C00000) then
!        COTPHI = C00000
!        dadc = C00000
!        dadPhi = C00000
!      else
!        dadc = C00P50*CPHI/SPHI
!        dadPhi = - C00P50*COH/(SPHI*SPHI)
!      end if
!      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
!      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
!      if (abs(phi) == C00000) then
!       ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
!      end if
!
!      !Calculate K function and its derivatives
!      if (abs(Lode) <= ATTRAN) then
!        STA = sin(Lode)
!        CTA = cos(Lode)
!        K = CTA - STA*SPHI*C00IR3
!        dKdPhi = - C00IR3*CPHI*STA
!      else
!        SGN = SIGN(C00001,Lode) !It puts the Lode's sign to the number 1
!        A = A1 + A2*SGN*SPHI
!        B = B1*SGN + B2*SPHI
!        K = A - B*S3TA
!        dKdPhi = A2*SGN*CPHI - B2*CPHI*S3TA
!      end if
!
!      !Operating..
!      Denom = (sqrt(J*J*K*K + ASPHI2))
!      Num =  J*J*K*dKdPhi + aSmooth*SPHI*SPHI*dadPhi + aSmooth*aSmooth*SPHI*CPHI
!
!      !Derivative DF/Dc
!      DFDSP(1) = aSmooth*SPHI*SPHI*dadc/Denom - CPHI
!
!      !Derivative DF/Dphi
!      DFDSP(2) = p*CPHI + Num/Denom + COH*SPHI
!
!      if (J <= C00000) then
!        DFDSP(1) = - CPHI
!        DFDSP(2) = p*CPHI + COH*SPHI
!      end if
!
!      end subroutine CalculateDerivativesYieldFunctSofteningParameters_PM4Sand
!
!
!      subroutine CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain_(factor,cp,cr,&
!                              phip,phir,psip,psir,EpsPEq,DSPDPEq) !!PM4Sand
!      !**********************************************************************
!      !
!      ! Calculation of the derivatives of the strength parameters with respect
!      ! the equivalent plastic shear strain
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !In Variables
!      real(REAL_TYPE), intent(in) :: EpsPEq
!      real(REAL_TYPE), intent(in) :: factor,cp,cr,phip,phir,psip,psir
!      !Out Variables
!      real(REAL_TYPE), intent(out), dimension(3):: DSPDPEq
!     
!      !Derivative Cohesion respect Equivalent Plastic Strain (Dc/DPEq)
!      DSPDPEq(1) = -factor * (cp - cr) * (exp(-factor*EpsPEq))
!      !Derivative Friction angle respect Equivalent Plastic Strain (Dphi/DPEq)
!      DSPDPEq(2) = -factor * (phip - phir) * (exp(-factor*EpsPEq))
!      !Derivative Dilatancy angle respect Equivalent Plastic Strain (Dpsi/DPEq)
!      DSPDPEq(3) = -factor * (psip - psir) * (exp(-factor*EpsPEq))
!
!      end subroutine CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain_ !PM4Sand
!
!
!      subroutine CalculateDerivativesEqPlasticStrainRespectPlasticStrain_(EpsP,EpsPEq,DEpsPEqDPS) !PM4Sand
!      !**********************************************************************
!      !
!      ! Calculation of the derivatives of the equivalent plastic shear strain
!      ! with respect the plastic strain
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local Variables
!      real(REAL_TYPE) :: k1, k2, k3
!      real(REAL_TYPE) :: EpsPM
!      real(REAL_TYPE), dimension(3) :: EpsDev
!      !In Variables
!      real(REAL_TYPE), intent(in), dimension(6) :: EpsP
!      real(REAL_TYPE), intent(in) :: EpsPEq
!      !Out Variables
!      real(REAL_TYPE), intent(out), dimension(6):: DEpsPEqDPS
!
!      
!      if (EpsPEq < 0.00000000001d0) then
!        k1 = 0.0d0
!      else 
!          k1 = 2.0d0/(3.0d0*EpsPEq) ! I moved this here
!      end if
!      k2 = k1 * 1.0d0/3.0d0
!      k3 = k1 * 2.0d0
!
!      EpsPM = k2 * (EpsP(1) + EpsP(2) + EpsP(3))
!      EpsDev(1) = EpsP(1)-EpsPM
!      EpsDev(2) = EpsP(2)-EpsPM
!      EpsDev(3) = EpsP(3)-EpsPM
!
!      DEpsPEqDPS(1) = k2 * ( 2.0d0*EpsDev(1) - EpsDev(2) - EpsDev(3))
!      DEpsPEqDPS(2) = k2 * (-EpsDev(1) + 2.0d0*EpsDev(2) - EpsDev(3))
!      DEpsPEqDPS(3) = k2 * (-EpsDev(1) - EpsDev(2) + 2.0d0*EpsDev(3))
!      DEpsPEqDPS(4) = k3 * EpsP(4)
!      DEpsPEqDPS(5) = k3 * EpsP(5)
!      DEpsPEqDPS(6) = k3 * EpsP(6)
!
!      end subroutine CalculateDerivativesEqPlasticStrainRespectPlasticStrain_!PM4Sand
!
!
!      subroutine CalculateEpsPEq_PM4Sand(EpsP,EpsPEq)
!      !**********************************************************************
!      !
!      ! Calculation of the equivalent plastic shear strain 
!      !
!      !**********************************************************************
!      
!      implicit none
!
!      !Local variables
!      real(REAL_TYPE):: EpsPM, C1, C2
!      real(REAL_TYPE), dimension(3) :: EpsDev
!      !In variables
!      real(REAL_TYPE), intent(in), dimension(6) :: EpsP
!      !Out variables
!      real(REAL_TYPE), intent(out) :: EpsPEq
!      
!      !EpsPEq = ((2/3)ep:ep)^(1/2), ep is the deviatoric plastic strain
!      
!      EpsPM = (1.0d0/3.0d0) * (EpsP(1) + EpsP(2) + EpsP(3))
!      EpsDev(1) = EpsP(1)-EpsPM
!      EpsDev(2) = EpsP(2)-EpsPM
!      EpsDev(3) = EpsP(3)-EpsPM
!      C1 = 2.0d0/3.0d0
!      C2 = C1 * 2.0d0
!      
!      EpsPEq = sqrt(C1*EpsDev(1)*EpsDev(1) + C1*EpsDev(2)*EpsDev(2) +  C1*EpsDev(3)*EpsDev(3) +&
!                     C2*EpsP(4)*EpsP(4) + C2*EpsP(5)*EpsP(5) + C2*EpsP(6)*EpsP(6))
!
!      end subroutine CalculateEpsPEq_PM4Sand
!      
!
!      !Subroutine CalculateIncrementSofteningParameters(DSPDPEq,DEpsPEqDPS,DEpsP,Dh)
!      !!**********************************************************************
!      !!
!      !! Calculation of the increment of the strenght parameters due to the softening
!      !!
!      !!**********************************************************************
!      !
!      !implicit none
!      !
!      !!Local variables
!      !real(REAL_TYPE) :: k
!      !!In variables
!      !real(REAL_TYPE), intent(in), dimension(3) :: DSPDPEq
!      !real(REAL_TYPE), intent(in), dimension(6) :: DEpsPEqDPS
!      !real(REAL_TYPE), intent(in), dimension(6) :: DEpsP
!      !!Out variables
!      !real(REAL_TYPE), intent(out), dimension(3) :: Dh
!      !
!      !
!      !k = DEpsPEqDPS(1)*DEpsP(1) + DEpsPEqDPS(2)*DEpsP(2) + DEpsPEqDPS(3)*DEpsP(3) + 
!      !*       DEpsPEqDPS(4)*DEpsP(4) + DEpsPEqDPS(5)*DEpsP(5) + DEpsPEqDPS(6)*DEpsP(6)
!      
!      
!      !Dh(1) = DSPDPEq(1)*k
!      !Dh(2) = DSPDPEq(2)*k
!      !Dh(3) = DSPDPEq(3)*k
!      
!      !Dh(1) = min (Dh(1) , 0.0d0)
!      !Dh(2) = min (Dh(2) , 0.0d0)
!      !Dh(3) = min (Dh(3) , 0.0d0)
!      
!      !end subroutine CalculateIncrementSofteningParameters
!
!
!      Subroutine CalculateSofteningParameters_PM4Sand(EpsPEq,factor,cp,cr,phip,phir,psip,psir,c,phi,psi)
!      !**********************************************************************
!      !
!      ! Calculation of strenght parameters (c, phi, psi)
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !In variables
!      real(REAL_TYPE), intent(in) :: EpsPEq,factor,cp,cr,phip,phir,psip,psir
!      !Out variables
!      real(REAL_TYPE), intent(out) :: c,phi,psi  
!
!      c = cr + (cp-cr)*exp(-factor*EpsPEq) 
!      phi = phir + (phip-phir)*exp(-factor*EpsPEq) 
!      psi = psir + (psip-psir)*exp(-factor*EpsPEq) 
!
!      end subroutine CalculateSofteningParameters_PM4Sand
!
!
!      Subroutine DetermineDSigAndDEpsP_PM4Sand(IntGlo,D1,D2,GG,c,phi,psi,Sig,DEpsPEqDPS,DSPDPEq,DEps,DSig,DEpsP)
!      !**********************************************************************
!      !
!      ! Calculation of the stress increment and plastic strain increment
!      !
!      !         dSig = Dep * dEps
!      !         dEpsP = Lambda * DPDSig
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local variables
!      integer :: i,k
!      real(REAL_TYPE) :: A,Ai,Denom,Fact,LambdaNum,Lambda
!      real(REAL_TYPE) :: p,J,Lode,S3TA !Invariants
!      real(REAL_TYPE), dimension(6,6) :: Num,Num1,Prod
!      real(REAL_TYPE), dimension(6) :: Denom1
!      real(REAL_TYPE), dimension(6) :: DPPDSig !Derivatives Plastic potential respect net stress
!      real(REAL_TYPE), dimension(6) :: DFDSig !Derivatives Yield function respect net stress
!      real(REAL_TYPE), dimension(2) :: DFDSP !Derivatives Yield function respect Soft Parameters
!      real(REAL_TYPE), dimension(6,6) :: Dep !Elastoplastic Constitutive Matrix 
!      !In Variables
!      real(REAL_TYPE), intent(in) :: c,phi,psi !Softening parameters
!      real(REAL_TYPE), intent(in) :: D1,D2,GG !Elastic parameters
!      real(REAL_TYPE), intent(in), dimension(6):: DEpsPEqDPS
!      real(REAL_TYPE), intent(in), dimension(6) :: Sig
!      real(REAL_TYPE), intent(in), dimension(3) :: DSPDPEq !Derivatives respect Equivalent Plastic Strain
!      real(REAL_TYPE), intent(in), dimension(6) :: DEps
!      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
!      !Out Variables
!      real(REAL_TYPE), intent(out), dimension(6) :: DSig
!      real(REAL_TYPE), intent(out), dimension(6) :: DEpsP
!
!      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
!      call CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
!      call CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)
!
!      !Parameter A (H = -A --> A>0 softening / A<0 hardening)
!      A = 0.0d0
!      Ai = (DFDSP(1)*DSPDPEq(1) + DFDSP(2)*DSPDPEq(2))
!      do i=1,6
!      A = A + Ai * DEpsPEqDPS(i) * DPPDSig(i)
!      end do
!
!      !Elastoplastic Constitutive Matrix (Dep)
!      do i=1,6
!        do k=1,6
!            Prod(i,k) =  DPPDSig(i) * DFDSig(k)
!        end do
!      end do
!
!      Num1(1,1) = D1*Prod(1,1) + D2*Prod(2,1) + D2*Prod(3,1)
!      Num1(1,2) = D1*Prod(1,2) + D2*Prod(2,2) + D2*Prod(3,2)
!      Num1(1,3) = D1*Prod(1,3) + D2*Prod(2,3) + D2*Prod(3,3)
!      Num1(1,4) = D1*Prod(1,4) + D2*Prod(2,4) + D2*Prod(3,4)
!      Num1(1,5) = D1*Prod(1,5) + D2*Prod(2,5) + D2*Prod(3,5)
!      Num1(1,6) = D1*Prod(1,6) + D2*Prod(2,6) + D2*Prod(3,6)
!
!      Num1(2,1) = D2*Prod(1,1) + D1*Prod(2,1) + D2*Prod(3,1)
!      Num1(2,2) = D2*Prod(1,2) + D1*Prod(2,2) + D2*Prod(3,2)
!      Num1(2,3) = D2*Prod(1,3) + D1*Prod(2,3) + D2*Prod(3,3)
!      Num1(2,4) = D2*Prod(1,4) + D1*Prod(2,4) + D2*Prod(3,4)
!      Num1(2,5) = D2*Prod(1,5) + D1*Prod(2,5) + D2*Prod(3,5)
!      Num1(2,6) = D2*Prod(1,6) + D1*Prod(2,6) + D2*Prod(3,6)
!
!      Num1(3,1) = D2*Prod(1,1) + D2*Prod(2,1) + D1*Prod(3,1)
!      Num1(3,2) = D2*Prod(1,2) + D2*Prod(2,2) + D1*Prod(3,2)
!      Num1(3,3) = D2*Prod(1,3) + D2*Prod(2,3) + D1*Prod(3,3)
!      Num1(3,4) = D2*Prod(1,4) + D2*Prod(2,4) + D1*Prod(3,4)
!      Num1(3,5) = D2*Prod(1,5) + D2*Prod(2,5) + D1*Prod(3,5)
!      Num1(3,6) = D2*Prod(1,6) + D2*Prod(2,6) + D1*Prod(3,6)
!
!      Num1(4,1) = GG*Prod(4,1)
!      Num1(4,2) = GG*Prod(4,2)
!      Num1(4,3) = GG*Prod(4,3)
!      Num1(4,4) = GG*Prod(4,4)
!      Num1(4,5) = GG*Prod(4,5)
!      Num1(4,6) = GG*Prod(4,6)
!
!      Num1(5,1) = GG*Prod(5,1)
!      Num1(5,2) = GG*Prod(5,2)
!      Num1(5,3) = GG*Prod(5,3)
!      Num1(5,4) = GG*Prod(5,4)
!      Num1(5,5) = GG*Prod(5,5)
!      Num1(5,6) = GG*Prod(5,6)
!
!      Num1(6,1) = GG*Prod(6,1)
!      Num1(6,2) = GG*Prod(6,2)
!      Num1(6,3) = GG*Prod(6,3)
!      Num1(6,4) = GG*Prod(6,4)
!      Num1(6,5) = GG*Prod(6,5)
!      Num1(6,6) = GG*Prod(6,6)
!
!
!
!      Num(1,1) = D1*Num1(1,1) + D2*Num1(1,2) + D2*Num1(1,3)
!      Num(1,2) = D2*Num1(1,1) + D1*Num1(1,2) + D2*Num1(1,3)
!      Num(1,3) = D2*Num1(1,1) + D2*Num1(1,2) + D1*Num1(1,3)
!      Num(1,4) = GG*Num1(1,4)
!      Num(1,5) = GG*Num1(1,5)
!      Num(1,6) = GG*Num1(1,6)
!
!      Num(2,1) = D1*Num1(2,1) + D2*Num1(2,2) + D2*Num1(2,3)
!      Num(2,2) = D2*Num1(2,1) + D1*Num1(2,2) + D2*Num1(2,3)
!      Num(2,3) = D2*Num1(2,1) + D2*Num1(2,2) + D1*Num1(2,3)
!      Num(2,4) = GG*Num1(2,4)
!      Num(2,5) = GG*Num1(2,5)
!      Num(2,6) = GG*Num1(2,6)
!
!      Num(3,1) = D1*Num1(3,1) + D2*Num1(3,2) + D2*Num1(3,3)
!      Num(3,2) = D2*Num1(3,1) + D1*Num1(3,2) + D2*Num1(3,3)
!      Num(3,3) = D2*Num1(3,1) + D2*Num1(3,2) + D1*Num1(3,3)
!      Num(3,4) = GG*Num1(3,4)
!      Num(3,5) = GG*Num1(3,5)
!      Num(3,6) = GG*Num1(3,6)
!
!      Num(4,1) = D1*Num1(4,1) + D2*Num1(4,2) + D2*Num1(4,3)
!      Num(4,2) = D2*Num1(4,1) + D1*Num1(4,2) + D2*Num1(4,3)
!      Num(4,3) = D2*Num1(4,1) + D2*Num1(4,2) + D1*Num1(4,3)
!      Num(4,4) = GG*Num1(4,4)
!      Num(4,5) = GG*Num1(4,5)
!      Num(4,6) = GG*Num1(4,6)
!
!      Num(5,1) = D1*Num1(5,1) + D2*Num1(5,2) + D2*Num1(5,3)
!      Num(5,2) = D2*Num1(5,1) + D1*Num1(5,2) + D2*Num1(5,3)
!      Num(5,3) = D2*Num1(5,1) + D2*Num1(5,2) + D1*Num1(5,3)
!      Num(5,4) = GG*Num1(5,4)
!      Num(5,5) = GG*Num1(5,5)
!      Num(5,6) = GG*Num1(5,6)
!
!      Num(6,1) = D1*Num1(6,1) + D2*Num1(6,2) + D2*Num1(6,3)
!      Num(6,2) = D2*Num1(6,1) + D1*Num1(6,2) + D2*Num1(6,3)
!      Num(6,3) = D2*Num1(6,1) + D2*Num1(6,2) + D1*Num1(6,3)
!      Num(6,4) = GG*Num1(6,4)
!      Num(6,5) = GG*Num1(6,5)
!      Num(6,6) = GG*Num1(6,6)
!
!
!
!      Denom1(1) = DFDSig(1)*D1 + DFDSig(2)*D2 + DFDSig(3)*D2
!      Denom1(2) = DFDSig(1)*D2 + DFDSig(2)*D1 + DFDSig(3)*D2
!      Denom1(3) = DFDSig(1)*D2 + DFDSig(2)*D2 + DFDSig(3)*D1
!      Denom1(4) = DFDSig(4)*GG
!      Denom1(5) = DFDSig(5)*GG
!      Denom1(6) = DFDSig(6)*GG
!
!      Denom =   Denom1(1)*DPPDSig(1) + Denom1(2)*DPPDSig(2) + &
!                  Denom1(3)*DPPDSig(3) + Denom1(4)*DPPDSig(4) + &
!             Denom1(5)*DPPDSig(5) + Denom1(6)*DPPDSig(6) - A
!
!      Fact = 1d0/Denom
!
!      !Dep
!      Dep(1,1) = D1 - Fact*Num(1,1)
!      Dep(1,2) = D2 - Fact*Num(1,2)
!      Dep(1,3) = D2 - Fact*Num(1,3)
!      Dep(1,4) = -Fact*Num(1,4)
!      Dep(1,5) = -Fact*Num(1,5)
!      Dep(1,6) = -Fact*Num(1,6)
!
!      Dep(2,1) = D2 - Fact*Num(2,1)
!      Dep(2,2) = D1 - Fact*Num(2,2)
!      Dep(2,3) = D2 - Fact*Num(2,3)
!      Dep(2,4) = -Fact*Num(2,4)
!      Dep(2,5) = -Fact*Num(2,5)
!      Dep(2,6) = -Fact*Num(2,6)
!
!      Dep(3,1) = D2 - Fact*Num(3,1)
!      Dep(3,2) = D2 - Fact*Num(3,2)
!      Dep(3,3) = D1 - Fact*Num(3,3)
!      Dep(3,4) = -Fact*Num(3,4)
!      Dep(3,5) = -Fact*Num(3,5)
!      Dep(3,6) = -Fact*Num(3,6)
!
!      Dep(4,1) = -Fact*Num(4,1)
!      Dep(4,2) = -Fact*Num(4,2)
!      Dep(4,3) = -Fact*Num(4,3)
!      Dep(4,4) = GG - Fact*Num(4,4)
!      Dep(4,5) = -Fact*Num(4,5)
!      Dep(4,6) = -Fact*Num(4,6)
!
!      Dep(5,1) = -Fact*Num(5,1)
!      Dep(5,2) = -Fact*Num(5,2)
!      Dep(5,3) = -Fact*Num(5,3)
!      Dep(5,4) = -Fact*Num(5,4)
!      Dep(5,5) = GG - Fact*Num(5,5)
!      Dep(5,6) = -Fact*Num(5,6)
!
!      Dep(6,1) = -Fact*Num(6,1)
!      Dep(6,2) = -Fact*Num(6,2)
!      Dep(6,3) = -Fact*Num(6,3)
!      Dep(6,4) = -Fact*Num(6,4)
!      Dep(6,5) = -Fact*Num(6,5)
!      Dep(6,6) = GG - Fact*Num(6,6)
!
!      !!!!!!!!! Calculate Plastic multipliler(Lambda)!!!!!!!!!!!!!!!!!
!      LambdaNum =   Denom1(1)*DEps(1) + Denom1(2)*DEps(2) + &
!                   Denom1(3)*DEps(3) + Denom1(4)*DEps(4) + &
!                   Denom1(5)*DEps(5) + Denom1(6)*DEps(6) 
!      Lambda =  LambdaNum/Denom
!
!      !!!!!!!!! Determine DSig --> (DSig = Dep*dEps) !!!!!!!!!!!
!      do i=1,6
!        DSig(i) = 0.0d0
!        do k=1,6
!            DSig(i) =  DSig(i) + Dep(i,k) * DEps(k)
!        end do
!      end do
!
!      !!!!!!!!! Determine DEpsP --> (DEpsP = Lambda*DPDSig) !!!!!!!!!!!!
!      do i=1,6
!        DEpsP(i) = Lambda * DPPDSig(i)
!      end do
!
!      end subroutine DetermineDSigAndDEpsP_PM4Sand
!
!
!      subroutine EndOfStepCorrection_PM4Sand(IntGlo,D1,D2,GG,IPL,F,Sig,DSPDPEq,DEpsPEqDPS,EpsP,c,phi,psi)
!      !**********************************************************************
!      !
!      ! Final correction of the yield surface drift (END OF STEP CORRECTION).
!      ! The stresses, the plastic strain and the strength parameters are corrected.
!      !
!      !**********************************************************************
!
!      implicit none
!
!      !Local variables
!      integer :: i
!      real(REAL_TYPE) :: p,J,Lode,S3TA !Invariants
!      real(REAL_TYPE) :: Lambda,param,c2,phi2,psi2,F2
!      real(REAL_TYPE) :: Denom,A,Ai
!      real(REAL_TYPE), dimension(2) :: DFDSP
!      real(REAL_TYPE), dimension(6) :: DPPDSig,DFDSig,Sig2,DEpsP,EpsP2
!      real(REAL_TYPE), dimension(6) :: Denom1
!      real(REAL_TYPE), dimension(3) :: Dh
!      !In Variables
!      integer, intent(in) :: IntGlo,IPL !Global ID of Gauss point or particle
!      real(REAL_TYPE), intent(in):: D1,D2,GG
!      real(REAL_TYPE), intent(in), dimension(3) :: DSPDPEq !Derivatives respect Equivalent Plastic Strain
!      real(REAL_TYPE), intent(in), dimension(6) :: DEpsPEqDPS !Derivatives respect Equivalent Plastic Strain
!      !InOut Variables
!      real(REAL_TYPE), intent(inout):: c,phi,psi
!      real(REAL_TYPE), intent(inout), dimension(6) :: Sig
!      real(REAL_TYPE), intent(inout), dimension(6) :: EpsP
!      real(REAL_TYPE), intent(inout):: F
!
!      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
!      call CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
!      call CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)
!
!      !Parameter A (hardening/softening parameter)
!      A = 0.0d0
!      Ai = (DFDSP(1)*DSPDPEq(1) + DFDSP(2)*DSPDPEq(2))
!      do i=1,6
!        A = A + Ai * DEpsPEqDPS(i) * DPPDSig(i)
!      end do
!
!      Denom1(1) = DPPDSig(1)*D1 + DPPDSig(2)*D2 + DPPDSig(3)*D2
!      Denom1(2) = DPPDSig(1)*D2 + DPPDSig(2)*D1 + DPPDSig(3)*D2
!      Denom1(3) = DPPDSig(1)*D2 + DPPDSig(2)*D2 + DPPDSig(3)*D1
!      Denom1(4) = DPPDSig(4)*GG
!      Denom1(5) = DPPDSig(5)*GG
!      Denom1(6) = DPPDSig(6)*GG
!
!      Denom = Denom1(1)*DFDSig(1) + Denom1(2)*DFDSig(2) + &
!             Denom1(3)*DFDSig(3) + Denom1(4)*DFDSig(4) + &
!              Denom1(5)*DFDSig(5) + Denom1(6)*DFDSig(6) - A
!
!      Lambda = F/Denom !factor correction
!
!      Sig2 = Sig - Lambda * Denom1 ! Sig2 = Sig + fact * Denom1 Stress corrected
!      DEpsP = Lambda * DPPDSig
!      EpsP2 = EpsP + DEpsP
!
!      if (IPL == 1)then
!        Dh = 0.0d0
!      else
!        param = DEpsPEqDPS(1) * DEpsP(1) + DEpsPEqDPS(2) * DEpsP(2) + DEpsPEqDPS(3) * DEpsP(3) + &
!               DEpsPEqDPS(4) * DEpsP(4) + DEpsPEqDPS(5) * DEpsP(5) + DEpsPEqDPS(6) * DEpsP(6)
!        Dh(1) = min (DSPDPEq(1)*param, 0.0d0)
!        Dh(2) = min (DSPDPEq(2)*param, 0.0d0)
!        Dh(3) = min (DSPDPEq(3)*param, 0.0d0)
!      end if
!
!      c2 = c + Dh(1)
!      phi2 = phi + Dh(2)
!      psi2 = psi + Dh(3)
!
!      call DetermineYieldFunctionValue(IntGlo,Sig2,c2,phi2,F2)
!      
!      if ((abs(F2) > abs(F)).or.(Denom == 0.0d0)) then !NormalCorrectionScheme
!        Denom = 0.0d0
!        Denom = DFDSig(1)*DFDSig(1) + DFDSig(2)*DFDSig(2) + &
!                 DFDSig(3)*DFDSig(3) + DFDSig(4)*DFDSig(4) + &
!                 DFDSig(5)*DFDSig(5) + DFDSig(6)*DFDSig(6)
!
!        Lambda = F/Denom
!        Sig = Sig - Lambda * DFDSig
!        DEpsP = Lambda * DPPDSig
!        EpsP = EpsP + DEpsP
!        call DetermineYieldFunctionValue(IntGlo,Sig,c,phi,F)
!      else
!        Sig = Sig2
!        EpsP = EpsP2
!        c = c2
!        phi = phi2
!        psi = psi2
!        F = F2
!      end if
!
!      end subroutine EndOfStepCorrection_PM4Sand
        
        
        end module ModMohrCoulombStrainSoftening