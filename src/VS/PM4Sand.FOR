    module PM4Sand
    !----------------------------------------------
    ! PM4Sand stress-strain model to capture the cyclic
    ! plastic deformation of soils. 
    !
    !----------------------------------------------
   
    ! Modules to be used in the PM4Sand material --> Temporary
    use ModElementEvaluationTETRA
    use ModElementEvaluationTRI
    use ModElementEvaluationQUAD
    !use ModElementEvaluation
    use ModString
    use ModFeedback
    use ModGlobalConstants
    use ModCounters
    use ModGeometryMath
    use ModNURBS

    implicit none

    contains
    
    
    
    
    
    
    
    
    
    
    
    
    
    !---------------------------------------------------------------------------------------
        
        
!██████╗░███╗░░░███╗░░██╗██╗░██████╗░█████╗░███╗░░██╗██████╗░
!██╔══██╗████╗░████║░██╔╝██║██╔════╝██╔══██╗████╗░██║██╔══██╗
!██████╔╝██╔████╔██║██╔╝░██║╚█████╗░███████║██╔██╗██║██║░░██║
!██╔═══╝░██║╚██╔╝██║███████║░╚═══██╗██╔══██║██║╚████║██║░░██║
!██║░░░░░██║░╚═╝░██║╚════██║██████╔╝██║░░██║██║░╚███║██████╔╝
!╚═╝░░░░░╚═╝░░░░░╚═╝░░░░░╚═╝╚═════╝░╚═╝░░╚═╝╚═╝░░╚══╝╚═════╝░
        
        
        ! PM4SAND
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        SUBROUTINE ESM_AM_PM4SAND(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
            DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)
        
        implicit none 
       
       ! input parameters 
       character(len = 80) :: CMNAME ! constitutive model name
        
       integer(INTEGER_TYPE), intent(in) :: NPT ! point number
       integer(INTEGER_TYPE), intent(in) :: NOEL ! element number
       integer(INTEGER_TYPE), intent(in) :: IDSET ! ??
       integer(INTEGER_TYPE), intent(in) :: NTENS ! tensor order --> 6 in 3D, 4 in 2D
       integer(INTEGER_TYPE), intent(in) :: NUMBEROFPHASES ! number of phases
       
       integer(INTEGER_TYPE), intent(in) :: NSTATEV ! number of state variables
       integer(INTEGER_TYPE), intent(in) :: NADDVAR ! number of additional variables 
       integer(INTEGER_TYPE), intent(in) :: NPROPS ! number of properties 
       
       ! output parameters 
       
       ! inout parameters 
       real(REAL_TYPE), intent(inout) :: PLASTICMULTIPLIER ! scalar plastic multiplier 
       
       real(REAL_TYPE), dimension(NSTATEV), intent(inout) :: STATEV ! actual state variables
       real(REAL_TYPE), dimension(NADDVAR), intent(inout) :: ADDITIONALVAR ! actual additional variables 
       real(REAL_TYPE), dimension(NPROPS), intent(inout) :: PROPS ! actual properties 
       
       real(REAL_TYPE), dimension(NTENS), intent(inout) :: STRESS ! final stress tensor
       real(REAL_TYPE), dimension(NTENS), intent(inout) :: DSTRAN ! strain increment
       
       real(REAL_TYPE), dimension(NTENS,NTENS), intent(inout) :: EUNLOADING ! unloading elastic matrix
       
        !---Local variables required in standard UMAT
        integer(INTEGER_TYPE) :: IStep, TimeStep
        
        real(REAL_TYPE), dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
        real(REAL_TYPE), dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
        real(REAL_TYPE), dimension(:), allocatable :: stran ! strain 
        real(REAL_TYPE), dimension(:), allocatable :: time ! time
        real(REAL_TYPE), dimension(:), allocatable :: predef
        real(REAL_TYPE), dimension(:), allocatable :: dpred    
        real(REAL_TYPE), dimension(:), allocatable :: coords
        real(REAL_TYPE), dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
        real(REAL_TYPE), dimension(:,:), allocatable :: drot
        real(REAL_TYPE), dimension(:,:), allocatable :: dfgrd0
        real(REAL_TYPE), dimension(:,:), allocatable :: dfgrd1
        
        real(REAL_TYPE) :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
        real(REAL_TYPE) :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
        real(REAL_TYPE) :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
        real(REAL_TYPE) :: pnewdt, dtime, temp, dtemp, celent
        real(REAL_TYPE) :: Value ! auxiliary variable holding any real valued number
        real(REAL_TYPE) :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation  

        integer(INTEGER_TYPE) :: ndi, nshr, layer, kspt, kstep, kinc     
        
        integer(INTEGER_TYPE) :: IDTask
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
    
        ! Initialization
        Eunloading = 0.0
        PlasticMultiplier = 0.0
        IDTask = 2
     
        !Rename additional variables --> 11 additional state variables
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
        
      
        ! set first call as 1... do we even need this if we have first call as a state variable
        IF ((IStep==1).and.(TimeStep==1)) IDTask = 1 !--> why do I need this.... I guess for initialization 
     
        ! we need to always do this every time step
      
        !! 2D strain assignement 
        !STATEV(44) = DSTRAN(1) ! mEpsilon_n
        !STATEV(45) = DSTRAN(2)! mEpsilon_n
        !STATEV(46) = DSTRAN(4)! mEpsilon_n
        !
        !! 2D plane stress assignement
        !STATEV(50) = STRESS(1) ! SIGMA_XX ! mSigma_n
        !STATEV(51) = STRESS(2) ! SIGMA_YY ! mSigma_n
        !STATEV(52) = STRESS(4) ! SIGMA_XY ! mSigma_n
        !
        !! initially mSigma_n = mSigma
        !STATEV(53) = STRESS(1) ! SIGMA_XX ! mSigma
        !STATEV(54) = STRESS(2) ! SIGMA_YY ! mSigma
        !STATEV(55) = STRESS(4) ! SIGMA_XY ! mSigma
            
        !---Call the UMAT
        call UMAT_AM_PM4SAND(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)


      
        !---Definition of Eunloading -> required to define the max time step
        Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
        
        !---Always define this value to run the simulation

        ! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
        
        return
     
       end subroutine ESM_AM_PM4SAND
        
        
        
        
        
        !*USER SUBROUTINES
      SUBROUTINE UMAT_AM_PM4SAND(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
      RPL,DDSDDT,DRPLDE,DRPLDT,&
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      implicit none
  
       ! input parameters 
       integer(INTEGER_TYPE), INTENT(IN) :: NPT, NTENS ! point number
       integer(INTEGER_TYPE), intent(in) :: NOEL ! element number
     
       real(REAL_TYPE), intent(in) :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
       real(REAL_TYPE), intent(in) :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
       real(REAL_TYPE), intent(in) :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
       real(REAL_TYPE), intent(in) :: pnewdt, dtime, temp, dtemp, celent
      
       character(len = 80), intent(in) :: CMNAME ! constitutive model name
        
       real(REAL_TYPE), dimension(NSTATEV), intent(inout) :: STATEV ! actual state variables 
        
       integer(INTEGER_TYPE), intent(in) :: NSTATEV ! number of state variables 
       
       integer(INTEGER_TYPE), intent(in) :: NPROPS ! number of properties 
       real(REAL_TYPE), dimension(NPROPS), intent(in) :: PROPS ! actual properties 
      
      real(REAL_TYPE), intent(in), dimension(:) :: STRESS ! size of this has been noted before
      
      real(REAL_TYPE), intent(inout), dimension(NTENS,NTENS) :: DDSDDE
      real(REAL_TYPE), intent(in), dimension(NTENS) :: DDSDDT
      real(REAL_TYPE), intent(in), dimension(NTENS) :: DRPLDE
      
      real(REAL_TYPE), intent(in), dimension(NTENS) :: STRAN
      real(REAL_TYPE), intent(in), dimension(NTENS) :: DSTRAN
      
      real(REAL_TYPE), dimension(NTENS) :: StrainIncrement
      
      real(REAL_TYPE), intent(in), dimension(2) :: TIME
      real(REAL_TYPE), intent(in), dimension(1) :: PREDEF
      real(REAL_TYPE), intent(in), dimension(1) :: DPRED
      
      real(REAL_TYPE), intent(in), dimension(3) :: COORDS
      real(REAL_TYPE), intent(in), dimension(3,3) :: DROT
      real(REAL_TYPE), intent(in), dimension(3,3) :: DFGRD0
      real(REAL_TYPE), intent(in), dimension(3,3) :: DFGRD1

      integer(INTEGER_TYPE), intent(in) :: ndi, nshr, layer, kspt, kstep, kinc     
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
    
      real(REAL_TYPE), dimension(3,3) :: DE_2D  !, intent(out)
      real(REAL_TYPE), dimension(3) :: dSig_2D !, intent(out)
      real(REAL_TYPE), dimension(3) :: Sig_2D !, intent(out)
      real(REAL_TYPE), dimension(3) :: dEpsP_2D !, intent(out)
        
      INTEGER(INTEGER_TYPE) :: Iter
      
   
      Iter = 0
     
    
    !subroutine FullConstructor(Dr, G0, PostShake, me2p, &
    !                       P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, phi_cv, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
    !                       Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp, & ! z_max, Ado, 
    !                       KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
    !                       dEpsilonE, Sigma_n, Sigma, &
    !                       FirstCall, Ce)
    
      call PM4SandFullConstructor(PROPS(1), PROPS(2), PROPS(4), PROPS(5), & ! props(3) remians unchanged
                     STATEV(1), STATEV(2), STATEV(3), STATEV(4), STATEV(5), STATEV(6), STATEV(7), STATEV(8), STATEV(9), STATEV(10), &
                     STATEV(11), STATEV(12), STATEV(13), STATEV(14), STATEV(15), STATEV(16), STATEV(17), STATEV(18), &
                     STATEV(19), STATEV(20), STATEV(21), STATEV(22), STATEV(23), STATEV(24), STATEV(25), STATEV(26), STATEV(27), &
                     STATEV(28), STATEV(29), STATEV(30), STATEV(31), STATEV(32), STATEV(33), STATEV(34), &!STATEV(35), STATEV(36), &
                     STATEV(35:37), STATEV(38:40), STATEV(41:43), &
                     STATEV(44:46), STATEV(47:49), STATEV(50:52), &
                     STATEV(53:55), STATEV(56:58), STATEV(59:61), &
                     STATEV(62:64), STATEV(65:67), STATEV(68:70), &
                     STATEV(71:73), STATEV(74:76), STATEV(77:79), &
                     STATEV(80:82), STATEV(83:85), STATEV(86:88), &
                     STATEV(89), STATEV(90), STATEV(91), &
                     STATEV(92:94), STATEV(95:97), STATEV(98:100), &
                     STATEV(101), &
                     DE_2D)  ! This references the in/out variable Sigma_n
    
    StrainIncrement = DSTRAN * (-1.0) ! geotechnical engineering sign convention (compression positive)
    
    !-------------------
    ! we obtain the strain, mEpsilon, from the element and then input it into the integrate matrix
    ! we will find the strain increment inside the the integrate subroutine
    !call setTrialStrain(StrainIncrement, &!mEpsilon, mEpsilon_n, &
    !    Sigma, Sigma_n, &
    !    rrr, &
    !    Mb, Md, K_p, &
    !    Fabric, Fabric_in, Fabric_in_n, Fabric_n, &
    !    Alpha, Alpha_n, Alpha_in, Alpha_in_n, Alpha_in_p, Alpha_in_p_n, &
    !    Alpha_in_true, Alpha_in_true_n, Alpha_in_max, Alpha_in_max_n, Alpha_in_min, Alpha_in_min_n, &
    !    zpeak, zcum, pzp, zxp, z_max, Ado, KK, GG, Mcur, pzpFlag, P_atm, Dr, Mc, nd, nb, mm, &
    !    EpsilonE, mEpsilonE_n) !initStress, 
    
    !--> need to have the right inputs and outputs
    !call PM4SandIntegrate(StrainIncrement, & !strain_from_element, mEpsilon_n,
    !    G0, hpo, PostShake, me2p, &    
    !    P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, phi_cv, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
    !    Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp,  & !z_max, Ado,
    !    KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
    !    dEpsilonE, Sigma_n, Sigma, &
    !    Ce)
    
    
    !subroutine PM4SandIntegrate(StrainIncrement, & !strain_from_element, mEpsilon_n,
    !G0, hpo, PostShake, me2p, &    
    !P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
    !Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp,  & 
    !KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
    !dEpsilonE, Sigma_n, Sigma, &
    !Ce)
    
    !call PM4SandIntegrate(StrainIncrement, & !strain_from_element, mEpsilon_n,
    !    PROPS(2), PROPS(3), PROPS(4), PROPS(5), &    
    !    STATEV(1), STATEV(2), STATEV(3), STATEV(4), STATEV(5), STATEV(6), STATEV(7), STATEV(8), STATEV(9), STATEV(10), &
    !    STATEV(12), STATEV(13), STATEV(14), STATEV(15), STATEV(16), STATEV(17), STATEV(18), & ! phi_cv
    !    STATEV(19), STATEV(20), STATEV(21), STATEV(22), STATEV(23), STATEV(24), STATEV(25), STATEV(26), STATEV(27), &
    !    STATEV(28), STATEV(29), STATEV(30), STATEV(31), STATEV(32), STATEV(33), STATEV(34), &
    !    STATEV(35:37), STATEV(38:40), STATEV(41:43), &
    !    STATEV(44:46), STATEV(47:49), STATEV(50:52), &
    !    STATEV(53:55), STATEV(56:58), STATEV(59:61), &
    !    STATEV(62:64), STATEV(65:67), STATEV(68:70), &
    !    STATEV(71:73), STATEV(74:76), STATEV(77:79), &
    !    STATEV(80:82), STATEV(83:85), STATEV(86:88), & 
    !    STATEV(89), STATEV(90), STATEV(91), &
    !    STATEV(92:94), STATEV(95:97), STATEV(98:100), &
    !    DE_2D)
    
    call PM4SandIntegrate(StrainIncrement, & !strain_from_element, mEpsilon_n,
        PROPS(2), PROPS(3), PROPS(4), PROPS(5), &    
        STATEV(1), STATEV(2), STATEV(3), STATEV(4), STATEV(5), &
        STATEV(6), STATEV(7), STATEV(8), STATEV(9), STATEV(10), &
        STATEV(12), STATEV(13), STATEV(14), STATEV(15), STATEV(16), &
        STATEV(17), STATEV(18), &! phi_cv
        STATEV(19), STATEV(20), STATEV(21), STATEV(22), STATEV(23), &
        STATEV(24), STATEV(25), STATEV(26), STATEV(27), &
        STATEV(28), STATEV(29), STATEV(30), STATEV(31), &
        STATEV(32), STATEV(33), STATEV(34),&
        STATEV(35:37), STATEV(38:40), STATEV(41:43),&
        STATEV(44:46), STATEV(47:49), STATEV(50:52),&
        STATEV(53:55), STATEV(56:58), STATEV(59:61),&
        STATEV(62:64), STATEV(65:67), STATEV(68:70),&
        STATEV(71:73), STATEV(74:76), STATEV(77:79), &
        STATEV(80:82), STATEV(83:85), STATEV(86:88), &
        STATEV(89), STATEV(90), STATEV(91), &
        STATEV(92:94), STATEV(95:97), STATEV(98:100), &
        DE_2D)
    
    
    
    
    
    !subroutine CommitState(StrainIncrement, &
    !    z_max, nu, G0, mm, Md, Mb, PostShake, Pmin, P_atm, Cgd, p_sedo, Fsed_min, me2p, Fabric, Fabric_in, &! input VoidRatio, 
    !    Alpha_in, Alpha_in_p, Alpha_in_true, Alpha_in_max, Alpha_in_min, & ! input
    !    KK, GG, Mcur, rrr, Alpha_in_n, Alpha_n, Alpha_in_p_n, Alpha_in_true_n, Alpha_in_max_n, Alpha_in_min_n, Sigma_n, !dFabric, & !NextVoidRatio, & ! output
    !    zcum, zpeak, Sigma, Alpha, Fabric_n, Fabric_in_n, &
    !    Ce) ! inout
    
    
    !call PM4SandCommitState(StrainIncrement, &
    !    STATEV(8), STATEV(12), PROPS(2), STATEV(17), STATEV(21), STATEV(20), PROPS(4), STATEV(25), STATEV(1), STATEV(13), STATEV(24), STATEV(23), PROPS(5), STATEV(46:48), STATEV(52:54), &! input VoidRatio, 
    !    STATEV(64:66), STATEV(70:72), STATEV(76:78), STATEV(82:84), STATEV(88:90), & ! input
    !    STATEV(34), STATEV(35)  , STATEV(36), STATEV(37:39), STATEV(61:63), STATEV(55:57), STATEV(67:69), STATEV(73:75), STATEV(79:81), STATEV(85:87), STATEV(97:100), & !dFabric, & !NextVoidRatio, & ! output
    !    STATEV(29), STATEV(28), STATEV(101:103), STATEV(58:60), STATEV(43:45), STATEV(49:51), &
    !    DE_2D) ! inout
    
    !DDSDDE(1,1) = DE_2D(1,1)
    !  DDSDDE(2,1) = DE_2D(2,1)
    !  DDSDDE(3,1) = DE_2D(2,1)
    !  DDSDDE(4,1) = 0.0
    !  DDSDDE(5,1) = 0.0
    !  DDSDDE(6,1) = 0.0
    !  
    !  DDSDDE(1,2) = DE_2D(1,2)
    !  DDSDDE(2,2) = DE_2D(2,2)
    !  DDSDDE(3,2) = DE_2D(1,2)
    !  DDSDDE(4,2) = 0.0
    !  DDSDDE(5,2) = 0.0
    !  DDSDDE(6,2) = 0.0
    !  
    !  DDSDDE(1,3) = DE_2D(2,1)
    !  DDSDDE(2,3) = DE_2D(1,2)
    !  DDSDDE(3,3) = DE_2D(2,2)
    !  DDSDDE(4,3) = 0.0
    !  DDSDDE(5,3) = 0.0
    !  DDSDDE(6,3) = 0.0
    !  
    !  DDSDDE(1,4) = 0.0
    !  DDSDDE(2,4) = 0.0
    !  DDSDDE(3,4) = 0.0
    !  DDSDDE(4,4) = DE_2D(3,3)
    !  DDSDDE(5,4) = 0.0
    !  DDSDDE(6,4) = 0.0
    !  
    !  DDSDDE(1,5) = 0.0
    !  DDSDDE(2,5) = 0.0
    !  DDSDDE(3,5) = 0.0
    !  DDSDDE(4,5) = 0.0
    !  DDSDDE(5,5) = DE_2D(3,3)
    !  DDSDDE(6,5) = 0.0
    !  
    !  DDSDDE(1,6) = 0.0
    !  DDSDDE(2,6) = 0.0
    !  DDSDDE(3,6) = 0.0
    !  DDSDDE(4,6) = 0.0
    !  DDSDDE(5,6) = 0.0
    !  DDSDDE(6,6) = DE_2D(3,3)
    !  
    !  
    !  
    !  
    !  
    !  ! change the sign convention based on the tension +ve convention
    !  STRESS(1) = STATEV(98)
    !  STRESS(2) = STATEV(99)
    !  STRESS(3) = 0.0!STATEV(98)
    !  STRESS(4) = STATEV(100)
    !  STRESS(5) = 0.0
    !  STRESS(6) = 0.0
    
    call PM4SandCommitState(StrainIncrement, &
        STATEV(8), STATEV(12), PROPS(2), STATEV(17), STATEV(21), &
        STATEV(20), PROPS(4), STATEV(25), STATEV(1), STATEV(13), &
        STATEV(24), &
        STATEV(23), PROPS(5), STATEV(44:46), STATEV(50:52), &! input VoidRatio, 
       STATEV(62:64), STATEV(68:70), STATEV(74:76), STATEV(80:82), &
      STATEV(86:88),  &! input
        STATEV(32), STATEV(33)  , STATEV(34), STATEV(35:37), &
      STATEV(59:61), STATEV(53:55), STATEV(65:67), STATEV(71:73), &
      STATEV(77:79), STATEV(83:85), STATEV(95:97),  &!dFabric, & !NextVoidRatio, & ! output
        STATEV(29), STATEV(28), STATEV(98:100), STATEV(56:58), &!47:49), ! new stress is in 98:100
      STATEV(41:43), STATEV(47:49), &
        DE_2D, &
      STATEV(19)) ! inout
    
    
      End subroutine 
    
    
    
    !---------------------------------------------------------------------------------------
    
    
      
      
      
      
      
      
      
      
      
      
      
      
      
    
    
    subroutine PM4SandMaterial(oData) ! input the properties and the state variables  !pData, 
    !----------------------------------------------
    ! To have the variables that are needed as constants 
    ! in the PM4Sand implementation. 
    !----------------------------------------------
    implicit none 
    
    !real(REAL_TYPE), dimension(3) :: pData ! this contains the properties of PM4Sand 
    ! -pData(1) = D_r --> apparent relative density (variable controlling dilatancy and )
    ! -pData(2) = G_0 --> shear modulus coefficient (small-strain shear modulus which can be calculated from shear wave velocity)
    ! -pData(3) = h_p0 --> contraction rate parameter (variable adjusts contraction rates which can be adjusted to get desiered CSR)
    real(REAL_TYPE), dimension(107), intent(inout) :: oData ! this contains the initial state variables (no fabric parameters)
    !real(REAL_TYPE), dimension(3), intent(in) :: pData ! this contains the initial state variables (no fabric parameters)
    !logical :: FirstCall
    
    ! oData has 24 state variables that contains all the information about the state variables 
    
    ! I am not sure if we need this subroutine --> we do... BUT ONLY THE FIRST TIME WHEN WE NEED TO INITIALIZE THOSE PARAMETERS
    ! we need the primary inputs here stored in some sort of properties parameters
    !pData(1) = Dr
    !pData(2) = G0
    !pData(3) = hp0
    
    !m_e_init
    !m_Mc
    
    if (oData(25) == 1) then  ! if first call is yes 
        
    !-------------------------------------------if first time to initialize
    ! initialization of list of parameters
    ! these are standard parameters    
    oData(1) = 101.3    !// P_atm  in kPa ! Numbering according to the manual
	
    oData(2) = -1       !// h0          ! parameter 1: -Adjusts the ratio of plastic modulus to elastic modulus. 
                                        !              -Default value is h0 = (0.25+Dr)/2
                                        !              -Minimum value of 0.30
                                        !              -Provides reasonable G/Gmax and damping rationships for the 
                                        !               default value of G0.
                                        !              -May require adjustment in combination with any adjustments to G0.        
	
    oData(3) = 0.8      !// emax        ! parameter 2: -Maximum void ratio that affects the computation of density.
                                        !              -Default value is 0.8. 
	
    oData(4) = 0.5      !// emin        ! parameter 2: -Minimum void ratio that affects the computation of density.
                                        !              -Default value is 0.5. 
	
    oData(5) = 0.5      !// nb          ! parameter 3: -Controls dilatancy and thus also the peak friction angles. 
                                        !              -Default value is 0.5. --> Dense of crit
                                        !              -Default value is 0.5/4. --> Loose of crit      
	
    oData(6) = 0.1      !// nd          ! parameter 4: -Controls the stress-ratio at which contraction changes to dilation,
                                        !               which is often referred to as phase transformation. 
                                        !              -Produces a phase transformation angle slightly smaller than phi_cv.  
                                        !              -Default value is 0.1. --> Dense of crit
                                        !              -Default value is 4*0.1. --> Loose of crit
	
    oData(7) = -1       !// Ado         ! parameter 5: -Default value is based on Bolton's dilatancy relationship at times of 
                                        !               initialization. 
                                        !              -Typical values will be between 1.2 and 1.5.
    
	oData(8) = -1       !// z_max       ! parameter 6: -Default value is computed at the time of initialization:
                                        !               z_max = 0.70 exp(-6.1 ksi_R0) <=20
                                        !              -Maximum value is 20.
                                        !              -if ksi_R is initially 0.0, z_max = 0.7.
                                        !              -z_max increases when ksi_R increases (dense critical states) to reach a 
                                        !               maximum value of 20.
                                        !              -May require varying if relationship between Dr and cyclic strength is significantly
                                        !               different from that implied by liquefaction correlations of Idriss and Boulanger (2008).
    
	oData(9) = 250.0    !// cz          ! parameter 7: -Default value is 250. 
                                        !              -Controls strain levels at which fabric effects become important. 
    
	oData(10) = -1       !// ceps         ! parameter 8: -Default value varies with Dr. 
                                        !              -Value is 5.0 for Dr<35%. --> Loose
                                        !              -Value linearly decreases to its minimum value of 1.0 at Dr=75% --> Dense
	
    oData(11) = 2 * sin(33.0 * (3.14159265359/180.0))    
                                        !// m_Mc      ! parameter 9: -Default value is for phi_cv of 33 degrees. 
    
	oData(12) = 0.3     !// nu          ! parameter 10: -Default value is 0.30. 
                                        !               -k0 would correspond to k0=nu/(1-nu) --> k0 value = 0.43
    
	oData(13) = 2.0     !// Cgd         ! parameter 11: -Default value is 2.0.
                                        !               -G0 degrades with increasing plastic deviatoric strains (z_cum). 
                                        !               -Maximum degredation approaches a factor of 1/Cgd
    
	oData(14) = -1      !// C_DR        ! parameter 
	
    oData(15) = -1      !// Ckaf        ! parameter 12: -Default value varies with Dr, as:
                                        !               Ckaf = 5 + (220*(Dro - 0.26)^3)
                                        !               -Ckaf = 4.0 for Dr<10% increases to its maximum value of 
                                        !                35.0 at Dr=77%
                                        !               -Controls the effect that sustained static shear stresses have on 
                                        !                plastic modulus. 
	
    oData(16) = 10.0    !// Q           ! parameter 13: -Default vlaue is 10.0 for quartzitic sand per Bolton (1986)
	
    oData(17) = 1.5     !// R           ! parameter 14: -Default value is 1.5. 
                                        !               -This is a slight increase from Bolton (1986) to lower the CSL 
                                        !                to better match results for direct simple shear loading.
	
    oData(18) = 0.01    !// m           ! parameter 15: -Default value is 0.01.
                                        !               -Provides reasonable modeling and numerical stability. 
	
    oData(19) = -1       !// Fsed_min   !  post liq stuff
	oData(20) = -1       !//p_sdeo      ! post liq stuff
	oData(21) = 5		 !// IntScheme  ! this should be a choice 
	oData(22) = 0		 !// TanType    ! ??
	oData(23) = 1.0e-8	 !// TolF       ! error threshold (yield surface)
	oData(24) = 1.0e-8	 !// TolR       ! error threshold --> do we even need this??
    
    oData(25) = 0        !// First call ! --> set to zero after first call be been done to overwrite the state variables
    oData(26) = 0        !// Post shake switch
    
    oData(27) = oData(1) / 200.0      !m_Pmin !//
    oData(28) = oData(26) * 5.0       !m_Pmin2 !//
    
    oData(29) = -1           !// mpzpFlag
    
    oData(30) = -1           !// me2p
    
    oData(31) = -1           !// mDGamma !-->
    oData(32) = -1           !// mDGamma_n
    
    oData(33) = -1           !// mK
    oData(34) = -1           !// mG
    
    oData(35) = -1           !// mVoidRatio
    
    oData(36) = -1           !// mKp
    
    oData(37) = -1           !// mzcum
    oData(38) = -1           !// mzpeak
    oData(39) = -1           !// mpzp
    oData(40) = -1           !// mzxp
    !oData() = -1        --> zxpPk?    
    
    oData(41) = -1           !// mMb --> bounding surface stress ratio
    oData(42) = -1           !// mMd --> dilatancy surface stress ratio
    oData(43) = -1           !// mMcur --> current stress ratio
    
    
    ! replace these mEpsilon_n and mEpsilon with dStran
    !oData(44) = -1           !// mEpsilon_n(1) --> old 
    !oData(45) = -1           !// mEpsilon_n(2)
    !oData(46) = -1           !// mEpsilon_n(3)
    
    !oData(47) = -1           !// mEpsilon(1) --> new
    !oData(48) = -1           !// mEpsilon(2)
    !oData(49) = -1           !// mEpsilon(3)
    
    !oData(50) = -1           !// mSigma_n(1) --> old
    !oData(51) = -1           !// mSigma_n(2)
    !oData(52) = -1           !// mSigma_n(3)
    !
    !oData(53) = -1           !// mSigma(1) --> new
    !oData(54) = -1           !// mSigma(2)
    !oData(55) = -1           !// mSigma(3)
    
    oData(56) = -1           !// mSigma_b(1) ----> difference
    oData(57) = -1           !// mSigma_b(2)
    oData(58) = -1           !// mSigma_b(3)
    
    oData(59) = -1           !// mEpsilonE_n(1) --> old
    oData(60) = -1           !// mEpsilonE_n(2)
    oData(61) = -1           !// mEpsilonE_n(3)
    
    oData(62) = -1           !// mEpsilonE(1) --> new
    oData(63) = -1           !// mEpsilonE(2)
    oData(64) = -1           !// mEpsilonE(3)
    
    oData(65) = -1           !//mAlpha_n(1) --> old
    oData(66) = -1           !//mAlpha_n(2)
    oData(67) = -1           !//mAlpha_n(3)
    
    oData(68) = -1           !//mAlpha(1) --> new        
    oData(69) = -1           !//mAlpha(2)         
    oData(70) = -1           !//mAlpha(3)        
    
    oData(71) = -1           !//mAlpha_in_n(1) --> old 
    oData(72) = -1           !//mAlpha_in_n(2)
    oData(73) = -1           !//mAlpha_in_n(3)
    
    oData(74) = -1           !//mAlpha_in(1) --> new
    oData(75) = -1           !//mAlpha_in(2)
    oData(76) = -1           !//mAlpha_in(3)
    
    oData(77) = -1           !//mAlpha_in_p_n(1) --> old
    oData(78) = -1           !//mAlpha_in_p_n(2)
    oData(79) = -1           !//mAlpha_in_p_n(3)
    
    oData(80) = -1           !// mAlpha_in_p(1) --> new
    oData(81) = -1            !// mAlpha_in_p(2)
    oData(82) = -1            !// mAlpha_in_p(3)
    
    oData(83) = -1           !//mAlpha_in_true_n(1) --> old
    oData(84) = -1           !//mAlpha_in_true_n(2)
    oData(85) = -1           !//mAlpha_in_true_n(3)
    
    oData(86) = -1            !// mAlpha_in_true(1) --> new
    oData(87) = -1            !// mAlpha_in_true(2)
    oData(88) = -1            !// mAlpha_in_true(3)
    
    oData(89) = -1           !//mAlpha_in_max_n(1) --> old
    oData(90) = -1           !//mAlpha_in_max_n(2)
    oData(91) = -1           !//mAlpha_in_max_n(3)
    
    oData(92) = -1            !// mAlpha_in_max(1) --> new
    oData(93) = -1            !// mAlpha_in_max(2)
    oData(94) = -1            !// mAlpha_in_max(3)
    
    oData(95) = -1           !//mAlpha_in_min_n(1) --> old
    oData(96) = -1           !//mAlpha_in_min_n(2)
    oData(97) = -1           !//mAlpha_in_min_n(3)
    
    oData(98) = -1            !// mAlpha_in_min(1) --> new
    oData(99) = -1            !// mAlpha_in_min(2)
    oData(100) = -1            !// mAlpha_in_min(3)
    
    oData(101) = -1           !//mFabric_n(1)      --> old
    oData(102) = -1           !//mFabric_n(2)      
    oData(103) = -1           !//mFabric_n(3)
    
    oData(104) = -1           !//mFabric(1)        --> new
    oData(105) = -1           !//mFabric(2)        
    oData(106) = -1           !//mFabric(3)       
    
    oData(107) = -1           !//mFabric_in_n(1) --> old
    oData(108) = -1           !//mFabric_in_n(2)
    oData(109) = -1           !//mFabric_in_n(3)
    
    oData(110) = -1            !// mFabric_in(1) --> new
    oData(111) = -1            !// mFabric_in(2)
    oData(112) = -1            !// mFabric_in(3)
    
    else 
        
        !just keep the original oData from the MP
    
    end if
    
    
    end subroutine PM4SandMaterial
    
    
    
    
    
    
    
    
    
    !subroutine FullConstructor(Dr, G0, P_atm, emax, emin, nb, nd, Ado, z_max, cz, ceps, phi_cv, nu, Cgd, Cdr, Ckaf, QQ_Bolton, RR_Bolton, mm, me2p, & ! Input variables 
    !                        h0, Mc, nu, Cgd, Cdr, Ckaf, QQ_Bolton, RR_Bolton, mm, Fsed_min, p_sedo, & ! Output variables 
    !                        Pmin, Pmin2, Mb, Md, K_p, zpeak, zcum, pzp, zxp, z_max, Ado, & ! Output variables 
    !                        KK, GG, Mcur, pzpFlag, rrr, Sigma_b, & ! Output variables 
    !                        Fabric_n, Fabric, Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, Alpha_in_p, & ! Output variables 
    !                        Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, Ce, & ! Output variables 
    !                        PostShake, Fsed_min_out, p_sedo_out, & ! Output postshake 
    !                        TolF, TolR, VoidRatio, & ! output variables
    !                        Sigma_n) ! Inout variables
    subroutine PM4SandFullConstructor(Dr, G0, PostShake, me2p, &
                           P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, phi_cv, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
                           Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp, & ! z_max, Ado, 
                           KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
                           dEpsilonE, Sigma_n, Sigma, &
                           FirstCall, Ce)
    !----------------------------------------------
    ! To populate the stress-strain variables needed 
    ! 
    !---------------------------------------------- 
    implicit none 
    
    ! I would say we call this subroutine first to initialize the properties
    
    ! input variables
    real(REAL_TYPE), intent(in) :: Dr       ! primiary parameter 1 ! prop 1
    real(REAL_TYPE), intent(in) :: G0       ! primiary parameter 2 ! prop 2
    
    ! prop 3 --> hpo
    
    real(REAL_TYPE), intent(in) :: PostShake ! prop 4 logical, intent(in) :: 
    real(REAL_TYPE), intent(in) :: me2p ! prop 5 logical, intent(in) :: 
    
    ! output variables 
    real(REAL_TYPE), intent(out) :: P_atm    ! statev 1
    real(REAL_TYPE), intent(out) :: h0       ! statev 2              ! secondary parameter 1
    real(REAL_TYPE), intent(out) :: emax     ! statev 3                       ! secondary parameter 2
    real(REAL_TYPE), intent(out) :: emin     ! statev 4                        ! secondary parameter 2
    real(REAL_TYPE), intent(out) :: nb       ! statev 5                        ! secondary parameter 3
    real(REAL_TYPE), intent(out) :: nd       ! statev 6                        ! secondary parameter 4
    real(REAL_TYPE), intent(out) :: Ado      ! statev 7                        ! secondary parameter 5
    real(REAL_TYPE), intent(out) :: z_max    ! statev 8                        ! secondary parameter 6
    real(REAL_TYPE), intent(out) :: cz       ! statev 9                        ! secondary parameter 7
    real(REAL_TYPE), intent(out) :: ceps     ! statev 10                          ! secondary parameter 8
    real(REAL_TYPE), intent(out) :: phi_cv   ! statev 11                        ! secondary parameter 9
    real(REAL_TYPE), intent(out) :: nu       ! statev 12                        ! secondary parameter 10
    real(REAL_TYPE), intent(out) :: Cgd      ! statev 13                        ! secondary parameter 11

    real(REAL_TYPE), intent(out) :: Ckaf     ! statev 14                        ! secondary parameter 12

    real(REAL_TYPE), intent(out) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
    real(REAL_TYPE), intent(out) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
    
    real(REAL_TYPE), intent(out) :: mm         ! statev 17                             ! secondary parameter 15

    real(REAL_TYPE), intent(out) :: VoidRatio   ! statev 18

    real(REAL_TYPE), intent(out) :: Mc         ! statev 19                      ! secondary parameter 1
    real(REAL_TYPE), intent(out) :: Mb         ! statev 20 
    real(REAL_TYPE), intent(out) :: Md         ! statev 21
    
    real(REAL_TYPE), intent(out) :: Cdr        ! statev 22                       ! secondary parameter 
    
    real(REAL_TYPE), intent(out) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
    real(REAL_TYPE), intent(out) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
    
    real(REAL_TYPE), intent(out) :: Pmin       ! statev 25 
    real(REAL_TYPE), intent(out) :: Pmin2      ! statev 26

    real(REAL_TYPE), intent(out) :: K_p        ! statev 27
    
    real(REAL_TYPE), intent(out) :: zpeak      ! statev 28
    real(REAL_TYPE), intent(out) :: zcum       ! statev 29
    
    real(REAL_TYPE), intent(out) :: pzp        ! statev 30 
    real(REAL_TYPE), intent(out) :: zxp        ! statev 31
    
    real(REAL_TYPE), intent(out) :: KK         ! statev 32
    real(REAL_TYPE), intent(out) :: GG         ! statev 33  
    real(REAL_TYPE), intent(out) :: Mcur       ! statev 34
    
    real(REAL_TYPE), intent(out), dimension(3) :: rrr   ! statev 35:37
    real(REAL_TYPE), intent(out), dimension(3) :: Sigma_b ! statev 38:40
    
    real(REAL_TYPE), intent(out), dimension(3) :: Fabric_n ! statev 41:43
    real(REAL_TYPE), intent(out), dimension(3) :: Fabric ! statev 44:46
    
    real(REAL_TYPE), intent(out), dimension(3) :: Fabric_in_n ! statev 47:49
    real(REAL_TYPE), intent(out), dimension(3) :: Fabric_in ! statev 50:52
    
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_n ! statev 53:55
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha ! statev 56:58
    
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_n ! statev 59:61
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in ! statev 62:64
    
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_p_n ! statev 65:67
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_p ! statev 68:70
    
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_true_n ! statev 71:73
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_true ! statev 74:76
    
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_max_n ! statev 77:79
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_max ! statev 80:82
    
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_min_n ! statev 83:85
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha_in_min ! statev 86:88
    
    real(REAL_TYPE), intent(out) :: pzpFlag ! statev 89 ! statev has to be same type throughout
    
    real(REAL_TYPE), intent(out) :: TolF ! statev 90
	real(REAL_TYPE), intent(out) :: TolR ! statev 91 
    
    real(REAL_TYPE), intent(out), dimension(3,3) :: Ce
    
    real(REAL_TYPE), intent(out), dimension(3) :: dEpsilonE ! statev 92:94
    
    real(REAL_TYPE), intent(inout), dimension(3) :: Sigma_n ! statev 95:97                 ! current strain      
    real(REAL_TYPE), intent(out), dimension(3) :: Sigma ! statev 98:100
    
    ! inout variables
    real(REAL_TYPE), intent(inout) :: FirstCall ! statev 101                 ! current strain      
    
    ! local variables 
    real(REAL_TYPE) :: h0_1
    real(REAL_TYPE) :: h0_2
    
    real(REAL_TYPE) :: ksi
    real(REAL_TYPE) :: Mcut
    real(REAL_TYPE) :: Mfin
    
    real(REAL_TYPE) :: p0
    
    real(REAL_TYPE), dimension(3) :: initStress 
    real(REAL_TYPE), dimension(3) :: I1  
    
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    
    ! I need to double check if I am initializing from a linear elastic constitutive model
    Sigma = 0
    
    ! initialize identity matrix 
    I1(1) = 1 
    I1(2) = 1
    I1(3) = 0
    
    if (FirstCall == 1.0) then 
        
        ! initialize atmospheric pressure
        if (P_atm < 0) then 
            P_atm = 101.3 !kPa
        end if 
        
        ! initialize parameter 1
        if (h0 < 0) then
            h0_1 = 0.3
            h0_2 = (0.25+Dr)/2
            h0 = max(h0_1, h0_2)
        end if 
    
        if (emax < 0) then
            emax = 0.8
        end if 
    
        if (emin < 0) then
            emin = 0.5
        end if
    
        if (nb < 0) then
            nb = 0.5
        end if
    
        if (nd < 0) then
            nd = 0.1
        end if 

        if (cz < 0) then
            cz = 250
        end if 
    
        ! initialize parameter 8
        ! this is called ce in opensees but here it is renamed to prevent confusion
        if (Dr < 0) then
            if (Dr > 0.75) then 
                ceps = 0.2
            else if (Dr < 0.55) then 
                ceps = 0.5 
            else
                ceps = 0.5 - ((Dr - 0.55) * 1.5)
            end if
        end if 
    
    
        ! initialize parameter 9        
        if (phi_cv < 0) then
            phi_cv = 33.0
        end if 
        
        if (Mc < 0) then 
            Mc = 2.0 * sin(phi_cv * (3.14159265359/180.0) )
        end if 
    
        ! initialize parameter 10
        if (nu < 0) then
            nu = 0.3
        end if 
    
        ! initialize parameter 11
        if (Cgd < 0) then
            Cgd = 2.0
        end if 
    
        ! initialize parameter 
        if (Cdr < 0) then
            Cdr = (5 + ( 25 * (Dr - 0.35)) )
            Cdr = min(Cdr, 10.0)
        end if 
    
        ! initialize parameter 12
        if (Ckaf < 0) then
            
            Ckaf = (5.0 + 220.0 * (Dr - 0.26)**3)
            
            if (Ckaf > 35) then 
                Ckaf = 35.0
            end if
            
            if (Ckaf < 4) then 
                Ckaf = 4.0
            end if
        
        end if 
        
        ! initialize parameter 13
        if (QQ_Bolton < 0) then
            QQ_Bolton = 10.0
        end if 
        
        ! initialize parameter 14
        if (RR_Bolton < 0) then
            RR_Bolton = 1.5
        end if 
    
        ! initialize parameter 15
        if (mm < 0) then
            mm = 0.01
        end if 
        
        ! initialize parameter 
        if (Fsed_min < 0) then
            Fsed_min = (0.03 * exp(2.6 * Dr))
            Fsed_min = min(Fsed_min, 0.99)
        end if 
        
        ! initialize parameter 
        if (p_sedo < 0) then
            p_sedo = (P_atm / 5.0)
        end if 
    
        ! store initial stress
        initStress = Sigma_n
        
        ! find p0 
        call GetTrace(initStress, p0)
        p0 = 0.5 * p0 ! mean effective stress 
    
        ! minimum p' 
        if (Pmin < 0) then
            Pmin = max(p0/200, P_atm/200)
        end if
        
        ! p_min for stress
        if (Pmin2 < 0) then
            Pmin2 = Pmin * 10 ! not sure why this is 10 --> this is correct see table 3.2
        end if 
        
        ! check if p0 is less than m_Pmin --> tension cutoff aspects
        if (p0 < Pmin) then
        
            ! initial p is small, set p to p_min and store the difference(mSigmab), the difference
		    ! will be added to the stress returned to element
            Sigma_n = Pmin * I1 ! change stress
            Sigma_b = initStress - Sigma_n ! stress difference 
            
            p0 = Pmin ! corresponds to the minimum --> this should be m_Pmin2, added this negative here --> removed... follow geotech convention
            
            Alpha = 0 ! initialize to zero
            Alpha_n = 0 ! initialize to zero
        
        else
        
            Sigma_n = initStress ! same stress
            Sigma_b = 0.0 ! stress difference is zero
            
            call GetDevPart(initStress, Alpha_n) ! get the deviatoric part of the stress
            Alpha_n = Alpha_n/p0 ! calculate old Alpha which is the deviatoric part normalized by the mean effective stress
    
        end if
    
        ! calculate the relative state parameter
        call GetKsi(Dr, p0, RR_Bolton, QQ_Bolton, Pmin, P_atm, ksi)
        
        ! this is where we initialize m_z_max
        if (z_max < 0) then 
            z_max = min(0.7*exp(-6.1*ksi), 20.0)
        end if
    
        ! bounding and dilatancy surface variations depends on whether it is dense of loose
        if (ksi < 0) then  !DENSE ! --> I switched the signs because it was wrong in the Opensees implementation 
            
            ! dense of critical
            Mb = Mc * exp(-1.0 * nb * ksi)
            Md = Mc * exp(       nd * ksi)
        
            if (Mb > 2.0) then 
                ! Warning, Mb is larger than 2, using Ado = 1.5.
                Ado = 1.5
            else 
                Ado = 2.5 * (asin(Mb/2.0) - asin(Mc/2.0)) / (Mb - Md) !--> equation withou fabric effects
            end if 
        
        else !LOOSE
        
            ! loose of critical 
            Mb = Mc * exp(-1.0 * (nb/4.0) * ksi)
            Md = Mc * exp(       (nd*4.0) * ksi)
            
            if (Ado < 0) then 
                Ado = 1.24
            end if 
        
        end if  
        
        
        ! check if initial stresses are inside bounding and dilatancy surface--------
        Mcut = max(Mb, Md) ! max(Mb, Md)
        
        call GetDevPart(Sigma_n, GetDevPart_result)
        call GetNorm_Contr(GetDevPart_result, Mfin)
        Mfin = sqrt(2.0) * Mfin
        Mfin = Mfin / p0
        
        ! check that initial stresses are inside the bounding surface (or dilatancy surface if it is greater) and compute the 
        ! committed back-stress and stress ratio tensors from the stress tensor --> TABLE 3.2
        if (Mfin > Mcut) then !--> not sure why we do this... if outside bounding or dilatancy
        
            ! current stress ratio calculation
            rrr = (Sigma_n - (p0*I1)) * (1/p0) * (Mcut/Mfin) !--> current stress ratio
            ! initial stress outside bounding/dilatancy surface, scale shear stress and store the difference(mSigma_b),
		    ! the difference will be added to the stress returned to element to maintain global equilibrium
        
            Sigma_n = (p0*I1) + (rrr*p0)
            Sigma_b = initStress - Sigma_n
            Alpha_n = rrr * (Mcut - mm) * (1/Mcut)
        
        end if
        
        !---------------------------------------------------------------------------
        zcum = 0.0
        
        ! 6. Calculate the initial values of elastic shear modulus, elastic bulk modulus, plastic modulus, dilatancy... dilatancy not included
        ! calculate the elastic moduli (mK and mG)
        call GetElasticModuli_(Sigma_n, zcum, z_max, nu, G0, Md, Mb, &
            PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, Fsed_min, me2p) ! mK, mG, mMcur, mzcum are outputs
        
        ! calculate elastic stiffness (3x3)
        call GetStiffness(KK, GG, Ce) ! maybe we don't need this....
        
        ! why are we doing this.... I think we are just initializing the plastic modulus 
        K_p = 100 * GG !--> plastic modulus
        
        Alpha = Alpha_n
        
        ! initialize to zero
        Alpha_in = 0.0
        Alpha_in_n = 0.0
        Alpha_in_p = 0.0
        Alpha_in_p_n = 0.0
        
        ! initialize to mAlpha_n --> current stress ratio
        Alpha_in_true = Alpha_n ! alpha_in 
        Alpha_in_true_n = Alpha_n
	    Alpha_in_max = Alpha_n ! alpha_inMax
	    Alpha_in_max_n = Alpha_n 
	    Alpha_in_min = Alpha_n ! alpha_inMin
	    Alpha_in_min_n = Alpha_n
        
        ! initialize fabric to zero... we start modifying the fabric at large strain
        Fabric = 0.0
        Fabric_n = 0.0
        Fabric_in = 0.0
        Fabric_in_n = 0.0
        
        ! elastic strain increments
        dEpsilonE = 0 
        
        ! initialize fabric terms
        zpeak = z_max / 100000.0 ! z_peak
        pzp = max(p0, Pmin) / 100.0 ! p_zp
        zxp = 0.0 ! zxp set to zero 
        pzpFlag = 1.0!.true. 
        
        ! establish tolerances
        if (TolF < 0) then 
            TolF = 1.0e-8	 !// TolF       ! error threshold (yield surface)
        end if 
    
        if (TolR < 0) then 
            TolR = 1.0e-8	 !// TolR       ! error threshold --> do we even need this??
        end if
        
        if (VoidRatio < 0) then 
            VoidRatio = emax - ((emax - emin) * Dr) !e_init
        end if 
    
    end if ! first call
    
    end subroutine PM4SandFullConstructor
    
    
    subroutine PM4SandCommitState(StrainIncrement, &
        z_max, nu, G0, mm, Md, Mb, PostShake, Pmin, P_atm, Cgd, p_sedo, Fsed_min, me2p, Fabric, Fabric_in, &! input VoidRatio, 
        Alpha_in, Alpha_in_p, Alpha_in_true, Alpha_in_max, Alpha_in_min, & ! input
        KK, GG, Mcur, rrr, Alpha_in_n, Alpha_n, Alpha_in_p_n, Alpha_in_true_n, Alpha_in_max_n, Alpha_in_min_n, Sigma_n, & !dFabric, & !NextVoidRatio, & ! output
        zcum, zpeak, Sigma, Alpha, Fabric_n, Fabric_in_n, &
        Ce, Mc) ! inout
    !-------------------------------------------------
    ! To commit stress state and state parameters
    !
    !-------------------------------------------------
    implicit none 
    
    ! input 
    real(REAL_TYPE), dimension(3), intent(in) ::  StrainIncrement
    
    real(REAL_TYPE), intent(in) ::  z_max
    real(REAL_TYPE), intent(in) ::  nu
    real(REAL_TYPE), intent(in) ::  G0
    real(REAL_TYPE), intent(in) ::  mm
    
    real(REAL_TYPE), intent(in) ::  Mc
    real(REAL_TYPE), intent(in) ::  Md
    real(REAL_TYPE), intent(in) ::  Mb
    real(REAL_TYPE), intent(in) ::  PostShake
    real(REAL_TYPE), intent(in) ::  Pmin
    real(REAL_TYPE), intent(in) ::  P_atm
    real(REAL_TYPE), intent(in) ::  Cgd
    real(REAL_TYPE), intent(in) ::  p_sedo
    real(REAL_TYPE), intent(in) ::  Fsed_min
    real(REAL_TYPE), intent(in) ::  me2p
    !real(REAL_TYPE), intent(in) ::  VoidRatio
    
    real(REAL_TYPE), dimension(3), intent(in) ::  Fabric
    real(REAL_TYPE), dimension(3), intent(in) ::  Fabric_in
    
    real(REAL_TYPE), dimension(3), intent(in) ::  Alpha_in
    real(REAL_TYPE), dimension(3), intent(in) ::  Alpha_in_p
    real(REAL_TYPE), dimension(3), intent(in) ::  Alpha_in_true
    real(REAL_TYPE), dimension(3), intent(in) ::  Alpha_in_max
    real(REAL_TYPE), dimension(3), intent(in) ::  Alpha_in_min
    
    
    ! output 
    real(REAL_TYPE), intent(out) ::  KK
    real(REAL_TYPE), intent(out) ::  GG
    real(REAL_TYPE), intent(out) ::  Mcur
    
    real(REAL_TYPE), dimension(3), intent(out) ::  rrr
    
    real(REAL_TYPE), dimension(3), intent(out) ::  Alpha_in_n
    real(REAL_TYPE), dimension(3), intent(out) ::  Alpha_n
    real(REAL_TYPE), dimension(3), intent(out) ::  Alpha_in_p_n
    real(REAL_TYPE), dimension(3), intent(out) ::  Alpha_in_true_n
    real(REAL_TYPE), dimension(3), intent(out) ::  Alpha_in_max_n
    real(REAL_TYPE), dimension(3), intent(out) ::  Alpha_in_min_n
    
    real(REAL_TYPE), dimension(3), intent(out) ::  Sigma_n
    
    
    real(REAL_TYPE), dimension(3,3), intent(out) ::  Ce
    
    !real(REAL_TYPE), intent(out) ::  NextVoidRatio
    
    ! inout 
    real(REAL_TYPE), intent(inout) ::  zcum
    real(REAL_TYPE), intent(inout) ::  zpeak
    real(REAL_TYPE), dimension(3), intent(inout) ::  Sigma
    real(REAL_TYPE), dimension(3), intent(inout) ::  Alpha
    
    real(REAL_TYPE), dimension(3), intent(inout) ::  Fabric_n
    real(REAL_TYPE), dimension(3), intent(inout) ::  Fabric_in_n
    
    ! local 
    real(REAL_TYPE), dimension(3) ::  dFabric !, intent(out)

    real(REAL_TYPE) :: dVolStrain
    real(REAL_TYPE) :: GetTrace_Sigma
    real(REAL_TYPE) :: pp
    real(REAL_TYPE) :: DoubleDot2_2_Contr_dFabric_dFabric
    real(REAL_TYPE) :: DoubleDot2_2_Contr_Fabric_Fabric
    
    real(REAL_TYPE), dimension(3) :: I1
    
    I1(1) = 1
    I1(2) = 1
    I1(3) = 0
    
    call GetElasticModuli_(Sigma, zcum, z_max, nu, G0, Md, Mb, PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, Fsed_min, me2p)
    !(Sigma, KK, GG, Mcur, zcum) !--> temporary subroutine
    
    call GetTrace(Sigma, GetTrace_Sigma) !--> temporary subroutine
    
    if ( (Mcur>Mb) .and. (me2p==1.0) ) then
        pp = 0.5 * GetTrace_Sigma ! real type 
        rrr = (Sigma - (pp*I1)) * (1/pp) * (Mb/Mcur)
        Sigma = (pp*I1) + (rrr*pp)
        Alpha = rrr * (Mb - mm)/Mb
    end if
    
    Alpha_in_n = Alpha_in
	Alpha_n = Alpha
	Alpha_in_p_n = Alpha_in_p
	Alpha_in_true_n = Alpha_in_true
	Alpha_in_max_n = Alpha_in_max
	Alpha_in_min_n = Alpha_in_min
	Sigma_n = Sigma
	!Epsilon_n = Epsilon
	!EpsilonE_n = EpsilonE
	dFabric = Fabric - Fabric_n
    
    ! update cumulated fabric 
    call DoubleDot2_2_Contr(dFabric, dFabric, DoubleDot2_2_Contr_dFabric_dFabric) !--> temporary subroutine
    zcum = zcum + sqrt(DoubleDot2_2_Contr_dFabric_dFabric/2.0)
     
    call DoubleDot2_2_Contr(Fabric, Fabric, DoubleDot2_2_Contr_Fabric_Fabric) !--> temporary subroutine
    zpeak = max( sqrt(DoubleDot2_2_Contr_Fabric_Fabric/2.0), zpeak)
    
    Fabric_n = Fabric
    Fabric_in_n = Fabric_in
    !DGamma_n = DGamma 
    
    !call GetTrace(StrainIncrement, dVolStrain) ! dVolStrain !--> temporary subroutine
    !NextVoidRatio = VoidRatio - ( (1 + VoidRatio) * dVolStrain )
    !
    call GetStiffness(KK, GG, Ce) !--> temporary subroutine
    !call GetElastoPlasticTangent(Sigma_n, Ce, RR, nn, Kp, Ce) !--> temporary subroutine
    !mCep_Consistent = mCe
    
    
    end subroutine PM4SandCommitState
                                       
                                       
              
                                       
                           
                           
                           
                           
                           
                           
     
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
      ! I NEED TO UNCOMMENT THE BELOWWWWWWWWWWWWWWWWWWWWWWWWWWWWW                     
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
    !subroutine commitState()
    !implicit none
    !
    !! input
    !real(REAL_TYPE), dimension(3), intent(in) :: Sigma
    !integer(INTEGER_TYPE), intent(in) :: me2p
    !
    !real(REAL_TYPE), intent(in) :: Mb
    !
    !real(REAL_TYPE), dimension(3), intent(in) :: Alpha
    !real(REAL_TYPE), dimension(3), intent(in) :: Alpha_in
    !real(REAL_TYPE), dimension(3), intent(in) :: Alpha_in_p
    !
    !real(REAL_TYPE), dimension(3), intent(in) :: Alpha_in_true
    !real(REAL_TYPE), dimension(3), intent(in) :: Alpha_in_max
    !real(REAL_TYPE), dimension(3), intent(in) :: Alpha_in_min
    !
    !
    !real(REAL_TYPE), dimension(3), intent(in) :: Fabric
    !real(REAL_TYPE), dimension(3), intent(in) :: Fabric_in
    !real(REAL_TYPE), intent(in) :: VoidRatio
    !
    !real(REAL_TYPE), intent(in) :: mm
    !
    !real(REAL_TYPE), intent(in) :: Pmin
    !real(REAL_TYPE), intent(in) :: P_atm
    !real(REAL_TYPE), intent(in) :: Pmin2
    !
    !real(REAL_TYPE), intent(in) :: z_max
    !
    !real(REAL_TYPE), intent(in) :: K_p
    !
    !! output
    !real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_n
    !real(REAL_TYPE), dimension(3), intent(out) :: Alpha_n
    !real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_p_n
    !
    !real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_min_n
    !real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_true_n
    !real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_max_n
    !
    !real(REAL_TYPE), dimension(3), intent(out) :: dFabric
    !
    !real(REAL_TYPE), intent(out) :: GG
    !real(REAL_TYPE), intent(out) :: KK 
    !
    !real(REAL_TYPE), dimension(3), intent(out) :: Sigma_n
    !
    !real(REAL_TYPE), dimension(3), intent(out) :: dEpsilonE_n
    !
    !real(REAL_TYPE), dimension(3), intent(out) :: Fabric_n
    !real(REAL_TYPE), dimension(3), intent(out) :: Fabric_in_n
    !
    !real(REAL_TYPE), intent(out) :: NextVoidRatio
    !
    !! inout
    !real(REAL_TYPE), dimension(3), intent(inout) :: rrr ! inout becuase of the me2p option
    !real(REAL_TYPE), dimension(3,3), intent(inout) :: Ce
    !
    !real(REAL_TYPE), intent(inout) :: zcum
    !
    !real(REAL_TYPE), intent(inout) :: zpeak 
    !
    !! Local variables
    !real(REAL_TYPE) :: pp
    !real(REAL_TYPE), dimension(3) :: nn
    !real(REAL_TYPE), dimension(3) :: RR
    ! 
    !real(REAL_TYPE) :: Mcur !local
    !
    !real(REAL_TYPE) :: dVolStrain 
    !real(REAL_TYPE), dimension(3) :: I1
    !real(REAL_TYPE), dimension(3) :: dEpsilonE
    !
    !real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    !
    !!real(REAL_TYPE) :: mzpeak
    !I1(1) = 1
    !I1(2) = 1
    !I1(3) = 0
    !
    !! Get elastic moduli
    !call GetElasticModuli_(Sigma, zcum, z_max, nu, G0, Md, Mb, &
    !    PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, Fsed_min, me2p)
    !
    !! Check condition for plasticity
    !if ((Mcur > Mb) .and. (me2p==1)) then
    !    call GetTrace(Sigma, pp)
    !    pp = 0.5 * pp
    !    rrr = (Sigma - (pp * I1)) * (Mb / (Mcur*pp) )
    !    Sigma = (pp * I1) + (rrr * pp)
    !    Alpha = rrr * (Mb - mm) / Mb
    !endif
    !
    !! Update state variables
    !Alpha_in_n = Alpha_in
    !Alpha_n = Alpha
    !Alpha_in_p_n = Alpha_in_p
    !Alpha_in_true_n = Alpha_in_true
    !Alpha_in_max_n = Alpha_in_max
    !Alpha_in_min_n = Alpha_in_min
    !Sigma_n = Sigma
    !dEpsilonE_n = dEpsilonE ! elastic strain increment
    !
    !! Calculate change in fabric
    !dFabric = Fabric - Fabric_n
    !
    !! Update cumulated fabric
    !call DoubleDot2_2_Contr(dFabric, dFabric, DoubleDot2_2_Contr_result)
    !zcum = zcum + sqrt(DoubleDot2_2_Contr_result / 2.0)
    !
    !call DoubleDot2_2_Contr(Fabric, Fabric, DoubleDot2_2_Contr_result)
    !zpeak = max(sqrt(DoubleDot2_2_Contr_result / 2.0), zpeak)
    !
    !! Update fabric variables
    !Fabric_n = Fabric
    !Fabric_in_n = Fabric_in
    !
    !! Update void ratio
    !call GetTrace(StrainIncrement, dVolStrain) !mEpsilon
    !NextVoidRatio = VoidRatio - (1 + VoidRatio) * dVolStrain!GetTrace(mEpsilon)
    !
    !! Calculate stiffness matrices
    !call GetStiffness(KK, GG, Ce)
    !call GetElastoPlasticTangent(Sigma_n, Ce, RR, nn, K_p, Ce)
    !
    !
    !end subroutine commitState


    !---------------------------------------------------------------------
    


    subroutine revertToLastCommit()
    ! Need to be added
    end subroutine revertToLastCommit

!---------------------------------------------------------------------


    !subroutine revertToStart()
    !! Added: C.McGann, U.Washington for InitialStateAnalysis
    !if (ops_InitialStateAnalysis) then
    !    ! Do nothing, keep state variables from last step
    !else
    !    ! Normal call for revertToStart (not initialStateAnalysis)
    !    call initialize(mSigma)
    !endif
    !
    !end subroutine revertToStart

!---------------------------------------------------------------------


    !subroutine getCopy(clone)
    !! Outputs
    !type(PM4Sand), pointer :: clone
    !! Create a clone
    !clone => PM4Sand()
    !clone = this
    !
    !end subroutine getCopy

!---------------------------------------------------------------------

    subroutine getType(materialType)
    ! Outputs
    character(LEN=*) :: materialType
    ! Set material type
    materialType = "PlaneStrain"
    end subroutine getType

!---------------------------------------------------------------------

    subroutine getOrder(materialOrder)
    ! Outputs
    integer, intent(out) :: materialOrder
    ! Set material order
    materialOrder = 3
    end subroutine getOrder

!---------------------------------------------------------------------

    !subroutine setResponse(argv, argc, output, responseObj)
    !! Inputs
    !character(LEN=*), dimension(*), intent(in) :: argv
    !integer, intent(in) :: argc
    !type(OPS_Stream), intent(inout) :: output
    !! Outputs
    !type(Response), pointer :: responseObj
    !! Set response based on input arguments
    !if (strcmp(argv(1), "stress") == 0 .or. strcmp(argv(1), "stresses") == 0) then
    !    responseObj => MaterialResponse(this, 1, this%getStress())
    !elseif (strcmp(argv(1), "strain") == 0 .or. strcmp(argv(1), "strains") == 0) then
    !    responseObj => MaterialResponse(this, 2, this%getStrain())
    !elseif (strcmp(argv(1), "state") == 0) then
    !    responseObj => MaterialResponse(this, 3, this%getState())
    !elseif (strcmp(argv(1), "alpha") == 0 .or. strcmp(argv(1), "backstressratio") == 0) then
    !    responseObj => MaterialResponse(this, 4, this%getAlpha())
    !elseif (strcmp(argv(1), "fabric") == 0) then
    !    responseObj => MaterialResponse(this, 5, this%getFabric())
    !elseif (strcmp(argv(1), "alpha_in") == 0 .or. strcmp(argv(1), "alphain") == 0) then
    !    responseObj => MaterialResponse(this, 6, this%getAlpha_in())
    !elseif (strcmp(argv(1), "trackers") == 0 .or. strcmp(argv(1), "tracker") == 0) then
    !    responseObj => MaterialResponse(this, 7, this%getTracker())
    !else
    !    responseObj => NULL()
    !endif
    !end subroutine setResponse
                                       
                                       
                                       
    
    !------------------------------------------------------------------------
    
!    subroutine getResponse(responseID, matInfo)
!    use Information
!    implicit none
!    ! Inputs
!    integer, intent(in) :: responseID
!    type(Information), intent(inout) :: matInfo
!    ! Local variables
!    integer :: res
!
!    select case (responseID)
!    case (-1)
!        res = -1
!    case (1)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getStress()
!        endif
!        res = 0
!    case (2)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getStrain()
!        endif
!        res = 0
!    case (3)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getState()
!        endif
!        res = 0
!    case (4)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getAlpha()
!        endif
!        res = 0
!    case (5)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getFabric()
!        endif
!        res = 0
!    case (6)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getAlpha_in()
!        endif
!        res = 0
!    case (7)
!        if (associated(matInfo%theVector)) then
!            matInfo%theVector = getTracker()
!        endif
!        res = 0
!    case default
!        res = -1
!    end select
!
!end subroutine getResponse
    
        
    !------------------------------------------------------------------------
    
    !subroutine sendSelf(commitTag, theChannel)
    !use Channel
    !implicit none
    !! Inputs
    !integer, intent(in) :: commitTag
    !type(Channel), intent(in) :: theChannel
    !! Local variables
    !integer :: res
    !integer, parameter :: dataSize = 101
    !real(REAL_TYPE) :: data(dataSize)
    !
    !! Assign data to vector
    !data(1) = real(this%getTag())
    !
    !data(2) = this%m_Dr
    !data(3) = this%m_G0
    !data(4) = this%m_hpo
    !data(5) = this%massDen
    !data(6) = this%m_P_atm
    !data(7) = this%m_h0
    !data(8) = this%m_emax
    !data(9) = this%m_emin
    !data(10) = this%m_e_init
    !data(11) = this%m_nb
    !data(12) = this%m_nd
    !data(13) = this%m_Ado
    !data(14) = this%m_cz
    !data(15) = this%m_ce
    !data(16) = this%m_Mc
    !data(17) = this%m_nu
    !data(18) = this%m_Cgd
    !data(19) = this%m_Cdr
    !data(20) = this%m_Ckaf
    !data(21) = this%m_Q
    !data(22) = this%m_R
    !data(23) = this%m_m
    !data(24) = this%m_z_max
    !data(25) = this%m_Fsed_min
    !data(26) = this%m_p_sedo
    !data(27) = this%m_FirstCall
    !data(28) = this%m_PostShake
    !
    !data(29) = this%mTolF
    !data(30) = this%mTolR
    !data(31) = this%mScheme
    !data(32) = this%mTangType
    !data(33) = this%m_Pmin
    !data(34) = this%m_Pmin2
    !data(35) = this%m_pzpFlag
    !data(36) = this%me2p
    !
    !data(37) = this%mDGamma
    !data(38) = this%mDGamma_n
    !data(39) = this%mK
    !data(40) = this%mG
    !data(41) = this%mVoidRatio
    !data(42) = this%mKp
    !data(43) = this%mzcum
    !data(44) = this%mzpeak
    !data(45) = this%mpzp
    !data(46) = this%mzxp
    !data(47) = this%mMb
    !data(48) = this%mMd
    !data(49) = this%mMcur
    !
    !data(50) = this%mEpsilon(1);     data(53) = this%mEpsilon_n(1);      data(56) = this%mSigma(1);     data(59) = this%mSigma_n(1);    data(62) = this%mSigma_b(1)
    !data(51) = this%mEpsilon(2);     data(54) = this%mEpsilon_n(2);      data(57) = this%mSigma(2);     data(60) = this%mSigma_n(2);    data(63) = this%mSigma_b(2)
    !data(52) = this%mEpsilon(3);     data(55) = this%mEpsilon_n(3);      data(58) = this%mSigma(3);     data(61) = this%mSigma_n(3);    data(64) = this%mSigma_b(3)
    !
    !data(65) = this%mEpsilonE(1);    data(68) = this%mEpsilonE_n(1);     data(71) = this%mAlpha(1);     data(74) = this%mAlpha_n(1);    data(77) = this%mAlpha_in_n(1)
    !data(66) = this%mEpsilonE(2);    data(69) = this%mEpsilonE_n(2);     data(72) = this%mAlpha(2);     data(75) = this%mAlpha_n(2);    data(78) = this%mAlpha_in_n(2)
    !data(67) = this%mEpsilonE(3);    data(70) = this%mEpsilonE_n(3);     data(73) = this%mAlpha(3);     data(76) = this%mAlpha_n(3);    data(79) = this%mAlpha_in_n(3)
    !
    !data(80) = this%mAlpha_in_p_n(1); data(83) = this%mAlpha_in_true_n(1); data(86) = this%mAlpha_in_max_n(1);   data(89) = this%mAlpha_in_min_n(1)
    !data(81) = this%mAlpha_in_p_n(2); data(84) = this%mAlpha_in_true_n(2); data(87) = this%mAlpha_in_max_n(2);   data(90) = this%mAlpha_in_min_n(2)
    !data(82) = this%mAlpha_in_p_n(3); data(85) = this%mAlpha_in_true_n(3); data(88) = this%mAlpha_in_max_n(3);   data(91) = this%mAlpha_in_min_n(3)
    !
    !data(92) = this%mFabric(1);        data(95) = this%mFabric_n(1);      data(98) = this%mFabric_in_n(1)
    !data(93) = this%mFabric(2);        data(96) = this%mFabric_n(2);      data(99) = this%mFabric_in_n(2)
    !data(94) = this%mFabric(3);        data(97) = this%mFabric_n(3);      data(100) = this%mFabric_in_n(3)
    !
    !! Send data vector via the channel
    !!res
    !
    !
    !end subroutine 
    !------------------------------------------------------------------------
    
!    subroutine recvSelf(commitTag, theChannel, theBroker)
!    use Channel
!    use FEM_ObjectBroker
!    implicit none
!    ! Inputs
!    integer, intent(in) :: commitTag
!    type(Channel), intent(in) :: theChannel
!    type(FEM_ObjectBroker), intent(in) :: theBroker
!    ! Local variables
!    integer :: res
!    integer, parameter :: dataSize = 101
!    real(REAL_TYPE) :: data(dataSize)
!
!    ! Receive data from the channel
!    res = theChannel%recvVector(this%getDbTag(), commitTag, data)
!    if (res < 0) then
!        print *, "WARNING: PM4Sand::recvSelf - failed to receive vector from channel"
!        return
!    endif
!
!    ! Set received data to class attributes
!    this%setTag(int(data(1)))
!
!    ! Set other attributes using received data
!    this%m_Dr = data(2)
!    this%m_G0 = data(3)
!    this%m_hpo = data(4)
!    this%massDen = data(5)
!    this%m_P_atm = data(6)
!    this%m_h0 = data(7)
!    this%m_emax = data(8)
!    this%m_emin = data(9)
!    this%m_e_init = data(10)
!    this%m_nb = data(11)
!    this%m_nd = data(12)
!    this%m_Ado = data(13)
!    this%m_cz = data(14)
!    this%m_ce = data(15)
!    this%m_Mc = data(16)
!    this%m_nu = data(17)
!    this%m_Cgd = data(18)
!    this%m_Cdr = data(19)
!    this%m_Ckaf = data(20)
!    this%m_Q = data(21)
!    this%m_R = data(22)
!    this%m_m = data(23)
!    this%m_z_max = data(24)
!    this%m_Fsed_min = data(25)
!    this%m_p_sedo = data(26)
!    this%m_FirstCall = data(27)
!    this%m_PostShake = data(28)
!
!    this%mTolF = data(29)
!    this%mTolR = data(30)
!    this%mScheme = data(31)
!    this%mTangType = data(32)
!    this%m_Pmin = data(33)
!    this%m_Pmin2 = data(34)
!    this%m_pzpFlag = data(36)
!    this%me2p = data(37)
!
!    this%mDGamma = data(38)
!    this%mDGamma_n = data(39)
!    this%mK = data(40)
!    this%mG = data(41)
!    this%mVoidRatio = data(42)
!    this%mKp = data(43)
!    this%mzcum = data(44)
!    this%mzpeak = data(45)
!    this%mpzp = data(46)
!    this%mzxp = data(47)
!    this%mMb = data(48)
!    this%mMd = data(49)
!    this%mMcur = data(50)
!
!    this%mEpsilon(1) = data(51);     this%mEpsilon_n(1) = data(52);      this%mSigma(1) = data(53);     this%mSigma_n(1) = data(54);    this%mSigma_b(1) = data(55)
!    this%mEpsilon(2) = data(56);     this%mEpsilon_n(2) = data(57);      this%mSigma(2) = data(58);     this%mSigma_n(2) = data(59);    this%mSigma_b(2) = data(60)
!    this%mEpsilon(3) = data(61);     this%mEpsilon_n(3) = data(62);      this%mSigma(3) = data(63);     this%mSigma_n(3) = data(64);    this%mSigma_b(3) = data(65)
!
!    this%mEpsilonE(1) = data(66);    this%mEpsilonE_n(1) = data(67);     this%mAlpha(1) = data(68);     this%mAlpha_n(1) = data(69);    this%mAlpha_in_n(1) = data(70)
!    this%mEpsilonE(2) = data(71);    this%mEpsilonE_n(2) = data(72);     this%mAlpha(2) = data(73);     this%mAlpha_n(2) = data(74);    this%mAlpha_in_n(2) = data(75)
!    this%mEpsilonE(3) = data(76);    this%mEpsilonE_n(3) = data(77);     this%mAlpha(3) = data(78);     this%mAlpha_n(3) = data(79);    this%mAlpha_in_n(3) = data(80)
!
!    this%mAlpha_in_p_n(1) = data(81); this%mAlpha_in_true_n(1) = data(82); this%mAlpha_in_max_n(1) = data(83); this%mAlpha_in_min_n(1) = data(84)
!    this%mAlpha_in_p_n(2) = data(85); this%mAlpha_in_true_n(2) = data(86); this%mAlpha_in_max_n(2) = data(87); this%mAlpha_in_min_n(2) = data(88)
!    this%mAlpha_in_p_n(3) = data(89); this%mAlpha_in_true_n(3) = data(90); this%mAlpha_in_max_n(3) = data(91); this%mAlpha_in_min_n(3) = data(92)
!
!    this%mFabric(1) = data(93);       this%mFabric_n(1) = data(94);       this%mFabric_in_n(1) = data(95)
!    this%mFabric(2) = data(96);       this%mFabric_n(2) = data(97);       this%mFabric_in_n(2) = data(98)
!    this%mFabric(3) = data(99);       this%mFabric_n(3) = data(100);      this%mFabric_in_n(3) = data(101)
!end subroutine recvSelf
    
    
    !------------------------------------------------------------------------
    
!    subroutine Print(s, flag)
!    use OPS_Stream
!    implicit none
!    ! Inputs
!    type(OPS_Stream), intent(inout) :: s
!    integer, intent(in) :: flag
!
!    ! Print information to the stream
!    call s%write("PM4Sand Material, tag: ")
!    call s%write(this%getTag())
!    call s%write("Type: ")
!    call s%write(this%getType())
!    call s%write(endln)
!
!end subroutine Print
    
    !------------------------------------------------------------------------
    
!    subroutine SetParameter(argv, argc, param, result)
!    implicit none
!    ! Inputs
!    character(len=*), dimension(*), intent(in) :: argv
!    integer, intent(in) :: argc
!    type(Parameter), intent(inout) :: param
!    ! Output
!    integer, intent(out) :: result
!
!    ! Local variables
!    integer :: theMaterialTag
!
!    ! Check if argc is less than 2
!    if (argc < 2) then
!        result = -1
!        return
!    end if
!
!    ! Convert argv[1] to integer
!    read(argv(1), *) theMaterialTag
!
!    ! Check if theMaterialTag matches the current object's tag
!    if (theMaterialTag == this%getTag()) then
!        ! Check argv[0] for different cases and add objects accordingly
!        if (trim(argv(0)) == "updateMaterialStage") then
!            print *, this%getTag(), "update Material Stage"
!            result = param%addObject(1, this)
!        else if (trim(argv(0)) == "materialState") then
!            result = param%addObject(5, this)
!        else if (trim(argv(0)) == "IntegrationScheme") then
!            result = param%addObject(2, this)
!        else if (trim(argv(0)) == "refShearModulus" .or. trim(argv(0)) == "ShearModulus") then
!            result = param%addObject(6, this)
!        else if (trim(argv(0)) == "poissonRatio") then
!            result = param%addObject(7, this)
!        else if (trim(argv(0)) == "FirstCall") then
!            result = param%addObject(8, this)
!        else if (trim(argv(0)) == "voidRatio") then
!            result = param%addObject(9, this)
!        else if (trim(argv(0)) == "PostShake") then
!            result = param%addObject(13, this)
!        end if
!    else
!        result = -1
!    end if
!
!end subroutine SetParameter
    
    !------------------------------------------------------------------------
!    subroutine UpdateParameter(responseID, info)
!    implicit none
!    ! Inputs
!    integer, intent(in) :: responseID
!    type(Information), intent(in) :: info
!
!    ! Local variables
!    real(REAL_TYPE) :: eps_v
!
!    ! Update material parameters based on responseID
!    if (responseID == 1) then
!        me2p = info%theInt
!    else if (responseID == 5) then
!        me2p = int(info%theDouble)
!    else if (responseID == 2) then
!        mScheme = int(info%theDouble)
!    else if (responseID == 6) then
!        m_G0 = info%theDouble
!    else if (responseID == 7) then
!        m_nu = info%theDouble
!    else if (responseID == 8) then
!        m_FirstCall = info%theInt
!        call initialize(mSigma_n)
!        print *, this%getTag(), "initialize"
!    else if (responseID == 9) then
!        eps_v = GetTrace(mEpsilon)
!        m_e_init = (info%theDouble + eps_v) / (1.0 - eps_v)
!    else if (responseID == 13) then
!        m_PostShake = 1
!        ! mElastFlag = 1
!        call GetElasticModuli(mSigma, mK, mG, mMcur, mzcum)
!        print *, this%getTag(), "activate post shaking reconsolidation"
!    else
!        ! Invalid responseID
!        return -1
!    end if
!
!    ! Return success
!    return 0
!
!end subroutine UpdateParameter
    
    
    
    !------------------------------------------------------------------------
 !   subroutine initialize(P_atm, Dr, nb, nd, mm, Mc, nu, G0, PostShake, me2p, & ! inputs
 !                     Pmin, Pmin2, Mb, Md, K_p, zpeak, zcum, pzp, zxp, z_max, Ado, & ! outputs
 !                     KK, GG, Mcur, pzpFlag, rrr, Sigma_b, Fabric_n, Fabric, & ! outputs
 !                     Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, & ! outputs
 !                     Alpha_in_p_n, Alpha_in_p, Alpha_in_true_n, Alpha_in_true, & ! outputs
 !                     Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, Ce, & ! outputs
 !                     Sigma_n) ! inouts
 !   
 !   implicit none
 !   ! Initialize PM4Sand Material: involves the initialization of some internal parameters
 !   
 !   ! input variables 
 !   real(REAL_TYPE), intent(in) :: P_atm
 !   real(REAL_TYPE), intent(in) :: Dr
 !   
 !   real(REAL_TYPE), intent(in) :: nb
 !   real(REAL_TYPE), intent(in) :: nd
 !   real(REAL_TYPE), intent(in) :: mm
 !   
 !   real(REAL_TYPE), intent(in) :: Mc
 !
 !   real(REAL_TYPE), intent(in) :: nu
 !   real(REAL_TYPE), intent(in) :: G0
 !   
 !   logical, intent(in) :: PostShake
 !   logical, intent(in) :: me2p
 !   
 !   
 !   
 !   
 !   
 !   
 !   ! output variables 
 !   real(REAL_TYPE), intent(out) :: Pmin
 !   real(REAL_TYPE), intent(out) :: Pmin2
 !   
 !   real(REAL_TYPE), intent(out) :: Mb
 !   real(REAL_TYPE), intent(out) :: Md
 !   real(REAL_TYPE), intent(out) :: K_p
 !   
 !   real(REAL_TYPE), intent(out) :: zpeak
 !   real(REAL_TYPE), intent(out) :: zcum
 !   real(REAL_TYPE), intent(out) :: pzp
 !   real(REAL_TYPE), intent(out) :: zxp
 !   real(REAL_TYPE), intent(out) :: z_max   
 !   
 !   real(REAL_TYPE), intent(out) :: Ado
 !   
 !   real(REAL_TYPE), intent(out) :: KK
 !   real(REAL_TYPE), intent(out) :: GG
 !   real(REAL_TYPE), intent(out) :: Mcur
 !   
 !   logical, intent(out) :: pzpFlag
 !   
 !   real(REAL_TYPE), intent(out), dimension(3) :: rrr
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Sigma_b
 !   
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Fabric_n
 !   real(REAL_TYPE), dimension(3), intent(out) :: Fabric
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Fabric_in_n
 !   real(REAL_TYPE), dimension(3), intent(out) :: Fabric_in 
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_n
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_n
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_p_n
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_p
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_true_n
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_true
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_max_n
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_max
 !   
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_min_n
 !   real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_min
 !      
 !   real(REAL_TYPE), dimension(3,3), intent(out) :: Ce
 !
 !   
 !   ! inout variables
 !   real(REAL_TYPE), dimension(3), intent(inout) :: Sigma_n
 !   
 !   
 !   ! Local parameters
 !   real(REAL_TYPE), dimension(3) :: initStress  ! parameter
 !   
 !   real(REAL_TYPE) :: ksi
 !   real(REAL_TYPE) :: Mcut
 !   real(REAL_TYPE) :: Mfin
 !   
 !   real(REAL_TYPE) :: p0
 !   
 !   real(REAL_TYPE), dimension(3) :: I1 ! identity matrix 
 !   
 !   real(REAL_TYPE) :: GetNorm_Contr_result
 !   real(REAL_TYPE) :: GetTrace_result
 !   real(REAL_TYPE) :: GetKsi_result
 !   real(REAL_TYPE), dimension(3) :: GetDevPart_result
 !   real(REAL_TYPE), dimension(3,3) :: GetStiffness_result
 !   
 !   
 !   
 !
 !   
 !   
 !   
 !   
 !
 !   
 !   
 !   
 !   
 !   
 !
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   
 !   ! store initial stress
 !   initStress = Sigma_n
 !   
 !   
 !   ! initialize identity matrix 
 !   I1(1) = 1 
 !   I1(2) = 1
 !   I1(3) = 0
 !   
 !   ! find p0 
 !   call GetTrace(initStress, p0)
 !   p0 = 0.5 * p0 ! mean effective stress 
 !   
 !   ! minimum p' 
 !   Pmin = max(p0/200, P_atm/200)
 !   
 !   ! p_min for stress
 !   Pmin2 = Pmin * 10 ! not sure why this is 10 --> this is correct see table 3.2
 !   
 !   
 !   ! check if p0 is less than m_Pmin --> tension cutoff aspects
 !   if (p0 < Pmin) then
 !       ! initial p is small, set p to p_min and store the difference(mSigmab), the difference
	!	! will be added to the stress returned to element
 !       Sigma_n = Pmin * I1 ! change stress
 !       Sigma_b = initStress - mSigma_n ! stress difference 
 !       p0 = Pmin ! corresponds to the minimum --> this should be m_Pmin2, added this negative here --> removed... follow geotech convention
 !       Alpha = 0 ! initialize to zero
 !       Alpha_n = 0 ! initialize to zero
 !   else
 !       Sigma_n = initStress ! same stress
 !       Sigma_b = 0.0 ! stress difference is zero
 !       call GetDevPart(initStress, GetDevPart_result) ! get the deviatoric part of the stress
 !       Alpha_n = GetDevPart_result/p0 ! calculate old Alpha which is the deviatoric part normalized by the mean effective stress
 !   end if
 !   
 !   ! calculate the relative state parameter
 !   call GetKsi(Dr, p0, RR_Bolton, QQ_Bolton, Pmin, P_atm, ksi)
 !   !ksi = GetKsi_result
 !   
 !   ! this is where we initialize m_z_max
 !   !if (z_max < 0) then 
 !    z_max = min(0.7*exp(-6.1*ksi), 20.0)
 !   !end if
 !   
 !   
 !   ! bounding and dilatancy surface variations depends on whether it is dense of loose
 !   if (ksi < 0) then  !DENSE ! --> I switched the signs because it was wrong in the Opensees implementation 
 !       ! dense of critical
 !       Mb = Mc * exp(-1.0 * nb * ksi)
 !       Md = Mc * exp(       nd * ksi)
 !       
 !           if (Mb > 2.0) then 
 !               ! Warning, Mb is larger than 2, using Ado = 1.5.
 !               Ado = 1.5
 !           else 
 !               Ado = 2.5 * (asin(Mb/2.0) - asin(Mc/2.0)) / (Mb - Md) !--> equation withou fabric effects
 !           end if 
 !           
 !       
 !   else !LOOSE
 !       ! loose of critical 
 !       Mb = Mc * exp(-1.0 * (nb/4.0) * ksi)
 !       Md = Mc * exp(       (nd*4.0) * ksi)
 !       
 !       !if (Ado < 0) then 
 !           Ado = 1.24
 !       !end if 
 !   end if 
 !   
 !   
 !   ! check if initial stresses are inside bounding and dilatancy surface--------
 !   Mcut = max(Mb, Md) ! max(Mb, Md)
 !   
 !   call GetDevPart(Sigma_n, GetDevPart_result)
 !   call GetNorm_Contr(GetDevPart_result, GetNorm_Contr_result)
 !   Mfin = sqrt(2.0) * GetNorm_Contr_result !--> I added a negative sign here
 !   Mfin = Mfin / p0  !--> Mfin calculation
 !   
 !   
 !   ! check that initial stresses are inside the bounding surface (or dilatancy surface if it is greater) and compute the 
 !   ! committed back-stress and stress ratio tensors from the stress tensor --> TABLE 3.2
 !   if (Mfin > Mcut) then !--> not sure why we do this... if outside bounding or dilatancy
 !       
 !       ! current stress ratio calculation
 !       rrr = (Sigma_n - (p0*I1)) * (1/p0) * (Mcut/Mfin) !--> current stress ratio
 !       ! initial stress outside bounding/dilatancy surface, scale shear stress and store the difference(mSigma_b),
	!	! the difference will be added to the stress returned to element to maintain global equilibrium
 !       
 !       Sigma_n = (p0*I1) + (rrr*p0)
 !       Sigma_b = initStress - Sigma_n
 !       Alpha_n = rrr * (Mcut - mm) * (1/Mcut)
 !       
 !   end if
 !   
 !   !---------------------------------------------------------------------------
 !   zcum = 0.0
 !   ! 6. Calculate the initial values of elastic shear modulus, elastic bulk modulus, plastic modulus, dilatancy... dilatancy not included
 !   ! calculate the elastic moduli (mK and mG)
 !   call GetElasticModuli_(Sigma_n, zcum, z_max, nu, G0, Md, Mb, &
 !       PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, Fsed_min, me2p)!, me2p) ! mK, mG, mMcur, mzcum are outputs
 !   
 !   ! calculate elastic stiffness (3x3)
 !   call GetStiffness(KK, GG, Ce)
 !
 !   !Cep = GetStiffness_result
 !   !Cep_consistent = GetStiffness_result
 !   
 !   ! why are we doing this.... I think we are just initializing the plastic modulus 
 !   K_p = 100 * GG !--> plastic modulus
 !   Alpha = Alpha_n
 !   
 !   ! initialize to zero
 !   Alpha_in = 0.0
 !   Alpha_in_n = 0.0
 !   Alpha_in_p = 0.0
 !   Alpha_in_p_n = 0.0
 !   
 !   ! initialize to mAlpha_n --> current stress ratio
 !   Alpha_in_true = Alpha_n ! alpha_in 
 !   Alpha_in_true_n = Alpha_n
	!Alpha_in_max = Alpha_n ! alpha_inMax
	!Alpha_in_max_n = Alpha_n 
	!Alpha_in_min = Alpha_n ! alpha_inMin
	!Alpha_in_min_n = Alpha_n
 !   
 !   ! initialize fabric to zero... we start modifying the fabric at large strain
 !   Fabric = 0.0
 !   Fabric_n = 0.0
 !   Fabric_in = 0.0
 !   Fabric_in_n = 0.0
 !   
 !   !mTracker = 0.0 ! --> not sure what this is doing 
 !   
 !   ! initialize fabric terms
 !   zpeak = z_max / 100000.0 ! z_peak
 !   pzp = max(p0, Pmin) / 100.0 ! p_zp
 !   zxp = 0.0 ! zxp set to zero 
 !   pzpFlag = .true. ! 
 !   
 !   ! zxpPk???
 !   
 !   ! we need to call setTrialStrain
 !   
 !   end subroutine initialize
 !   
    
    
    !------------------------------------------------------------------------
    
 !   subroutine setTrialStrain(StrainIncrement, & !strain_from_element, mEpsilon_n,
 !   Dr, G0, hpo, PostShake, me2p, &    
 !   P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, phi_cv, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
 !   Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp,  & !z_max, Ado,
 !   KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
 !   dEpsilonE, Sigma_n, Sigma, &
 !   Ce)
 !   
 !   implicit none
 !   
 !   ! input variables
 !   real(REAL_TYPE), intent(in) :: Dr       ! primiary parameter 1 ! prop 1
 !   real(REAL_TYPE), intent(in) :: G0       ! primiary parameter 2 ! prop 2
 !   real(REAL_TYPE), intent(in) :: hpo
 !   ! prop 3 --> hpo
 !   
 !   real(REAL_TYPE), intent(in) :: PostShake ! prop 4
 !   real(REAL_TYPE), intent(in) :: me2p ! prop 5
 !   
 !   ! output variables 
 !   real(REAL_TYPE), intent(in) :: P_atm    ! statev 1
 !   real(REAL_TYPE), intent(in) :: h0       ! statev 2              ! secondary parameter 1
 !   real(REAL_TYPE), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
 !   real(REAL_TYPE), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
 !   real(REAL_TYPE), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
 !   real(REAL_TYPE), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
 !   real(REAL_TYPE), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
 !   real(REAL_TYPE), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
 !   real(REAL_TYPE), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
 !   real(REAL_TYPE), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
 !   real(REAL_TYPE), intent(in) :: phi_cv   ! statev 11                        ! secondary parameter 9
 !   real(REAL_TYPE), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
 !   real(REAL_TYPE), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11
 !
 !   real(REAL_TYPE), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12
 !
 !   real(REAL_TYPE), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
 !   real(REAL_TYPE), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
 !   
 !   real(REAL_TYPE), intent(in) :: mm         ! statev 17                             ! secondary parameter 15
 !
 !   real(REAL_TYPE), intent(inout) :: VoidRatio   ! statev 18
 !
 !   real(REAL_TYPE), intent(inout) :: Mc         ! statev 19                      ! secondary parameter 1
 !   real(REAL_TYPE), intent(inout) :: Mb         ! statev 20 
 !   real(REAL_TYPE), intent(inout) :: Md         ! statev 21
 !   
 !   real(REAL_TYPE), intent(inout) :: Cdr        ! statev 22                       ! secondary parameter 
 !   
 !   real(REAL_TYPE), intent(inout) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
 !   real(REAL_TYPE), intent(inout) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
 !   
 !   real(REAL_TYPE), intent(in) :: Pmin       ! statev 25 
 !   real(REAL_TYPE), intent(in) :: Pmin2      ! statev 26
 !
 !   real(REAL_TYPE), intent(inout) :: K_p        ! statev 27
 !   
 !   real(REAL_TYPE), intent(in) :: zpeak      ! statev 28
 !   real(REAL_TYPE), intent(in) :: zcum       ! statev 29
 !   
 !   real(REAL_TYPE), intent(inout) :: pzp        ! statev 30 
 !   real(REAL_TYPE), intent(inout) :: zxp        ! statev 31
 !   !real(REAL_TYPE), intent(in) :: z_max      ! statev 32
 !   !
 !   !real(REAL_TYPE), intent(in) :: Ado        ! statev 33
 !   
 !   real(REAL_TYPE), intent(inout) :: KK         ! statev 34
 !   real(REAL_TYPE), intent(inout) :: GG         ! statev 35  
 !   real(REAL_TYPE), intent(inout) :: Mcur       ! statev 36
 !   
 !   
 !   real(REAL_TYPE), intent(inout), dimension(3) :: rrr   ! statev 37:39
 !   real(REAL_TYPE), intent(inout), dimension(3) :: Sigma_b ! statev 40:42
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Fabric_n ! statev 43:45
 !   real(REAL_TYPE), intent(out), dimension(3) :: Fabric ! statev 46:48
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Fabric_in_n ! statev 49:51
 !   real(REAL_TYPE), intent(out), dimension(3) :: Fabric_in ! statev 52:54
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Alpha_n ! statev 55:57
 !   real(REAL_TYPE), intent(out), dimension(3) :: Alpha ! statev 58:60
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_n ! statev 61:63
 !   real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in ! statev 64:66
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_p_n ! statev 67:69
 !   real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in_p ! statev 70:72
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_true_n ! statev 73:75
 !   real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in_true ! statev 76:78
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_max_n ! statev 79:81
 !   real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in_max ! statev 82:84
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_min_n ! statev 85:87
 !   real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in_min ! statev 88:90
 !   
 !   real(REAL_TYPE), intent(inout), dimension(3,3) :: Ce
 !   
 !   integer(INTEGER_TYPE), intent(inout) :: pzpFlag ! statev ! statev 91
 !   
 !   real(REAL_TYPE), intent(in) :: TolF ! statev 92
 !   
 !   ! = 1.0e-8	 !// TolF       ! error threshold (yield surface)
	!real(REAL_TYPE), intent(in) :: TolR ! statev 93 
 !   
 !   
 !   ! = 1.0e-8	 !// TolR       ! error threshold --> do we even need this??
 !   
 !   real(REAL_TYPE), intent(out), dimension(3) :: dEpsilonE ! statev 94:96
 !   
 !   real(REAL_TYPE), intent(out), dimension(3) :: Sigma ! statev 101:103
 !   
 !   !real(REAL_TYPE), intent(out), dimension(3,3) :: Ce !
 !   
 !   ! inout variables
	!real(REAL_TYPE), intent(in), dimension(3) :: Sigma_n ! statev 97:100                 ! current strain      
 !   
 !   
 !
 !   ! Input variables
 !   real(REAL_TYPE), dimension(3), intent(in) :: StrainIncrement!strain_from_element 
 !   ! --> Epsilon_xx, Epsilon_yy, Epsilon_xy
 !   !
 !   !! variables that need to be an input?
 !   !real(REAL_TYPE), dimension(3) :: mEpsilon
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mEpsilon_n
 !   !
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mEpsilonE, mEpsilonE_n
 !   !
 !   !
 !   !! variables 
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mSigma
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mSigma_n
 !   !
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: rrr
 !   !real(REAL_TYPE), intent(inout) :: mMb
 !   !real(REAL_TYPE), intent(inout) :: mMd
 !   !real(REAL_TYPE), intent(inout) :: mKp
 !   !
 !   !
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mFabric
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_in
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_in_n
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mFabric_n
 !   !
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_n
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_n
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_p
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_p_n
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_true
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_true_n
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_max
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_max_n
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_min
 !   !real(REAL_TYPE), dimension(3), intent(inout) :: mAlpha_in_min_n
 !   !
 !   !real(REAL_TYPE), intent(inout) :: mzpeak
 !   !real(REAL_TYPE), intent(inout) :: mzcum
 !   !real(REAL_TYPE), intent(inout) :: mpzp
 !   !real(REAL_TYPE), intent(inout) :: mzxp
 !   !real(REAL_TYPE), intent(inout) :: m_z_max     ! parameter 
 !   !real(REAL_TYPE), intent(inout) :: m_Ado
 !   !
 !   !
 !   !real(REAL_TYPE), intent(inout) :: mK
 !   !real(REAL_TYPE), intent(inout) :: mG
 !   !real(REAL_TYPE), intent(inout) :: mMcur
 !   !
 !   !logical, intent(inout) :: m_pzpFlag
 !   !
 !   !real(REAL_TYPE), intent(in) :: P_atm
 !   !real(REAL_TYPE), intent(in) :: m_Dr
 !   !real(REAL_TYPE), intent(in) :: m_Mc
 !   !real(REAL_TYPE), intent(in) :: m_nd
 !   !real(REAL_TYPE), intent(in) :: m_nb
 !   !real(REAL_TYPE), intent(in) :: m_m
 !   !
 !   ! Set trial strain and call integrate subroutine
 !   !mEpsilon = strain_from_element
 !   ! note that now compression is positive, tension is negative
 !   ! but I am not sure if this is appropriate
 !   !mEpsilon = mEpsilon * (-1.0) ! -1.0 is for geotechnical sign convention --> we do not want this.... maybe we do...
 !   
 !   
 !   !StrainIncrement = StrainIncrement * (-1.0)
 !   
 !   
 !   ! TABLE 3.1
 !   ! now we want to integrate
 !   !call integrate(StrainIncrement, & !mEpsilon, mEpsilon_n, &
 !   !    Sigma, Sigma_n, rrr, Mb, Md, K_p, Fabric, Fabric_in, Fabric_in_n, Fabric_n, &
 !   !    Alpha, Alpha_n, Alpha_in, Alpha_in_n, Alpha_in_p, Alpha_in_p_n, &
 !   !    Alpha_in_true, Alpha_in_true_n, Alpha_in_max, Alpha_in_max_n, Alpha_in_min, Alpha_in_min_n, &
 !   !    zpeak, zcum, pzp, zxp, z_max, Ado, KK, GG, Mcur, pzpFlag, P_atm, Dr, Mc, nd, nb, mm, &
 !   !    EpsilonE, EpsilonE_n)
 !   
 !   
 !   call integrate(StrainIncrement, & !strain_from_element, mEpsilon_n,
 !       Dr, G0, hpo, PostShake, me2p, &    
 !       P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, phi_cv, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
 !       Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp,  & !z_max, Ado,
 !       KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
 !       dEpsilonE, Sigma_n, Sigma, &
 !       Ce)
 !   
 !   !integrate(mEpsilon, mEpsilon_n, mSigma, mSigma_n, rrr, mMb, mMd, mKp, mFabric, mFabric_in, mFabric_in_n, mFabric_n, &
 !   !    mAlpha, mAlpha_n, mAlpha_in, mAlpha_in_n, mAlpha_in_p, mAlpha_in_p_n, mAlpha_in_true, mAlpha_in_true_n, mAlpha_in_max, mAlpha_in_max_n, mAlpha_in_min, mAlpha_in_min_n, &
 !   !    mzpeak, mzcum, mpzp, mzxp, m_z_max, m_Ado, mK, mG, mMcur, m_pzpFlag, m_P_atm, m_Dr, m_Mc, m_nd, m_nb, m_m, &
 !   !    mEpsilonE, mEpsilonE_n)
 !   
 !   end subroutine setTrialStrain
    
    
    
    
    
    !------------------------------------------------------------------------

    subroutine GetState(State)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(16), intent(out) :: State

    ! variables
    real(REAL_TYPE), dimension(3) :: mEpsilonE
    real(REAL_TYPE), dimension(3) :: mAlpha_n
    real(REAL_TYPE), dimension(3) :: mFabric_n
    real(REAL_TYPE), dimension(3) :: mAlpha_in_n
    real(REAL_TYPE) :: mVoidRatio
    real(REAL_TYPE) :: mDGamma_n
    real(REAL_TYPE) :: mG
    real(REAL_TYPE) :: mKp
    
    ! Assemble state parameters into the result vector
    State(1:3) = mEpsilonE
    State(4:6) = mAlpha_n
    State(7:9) = mFabric_n
    State(10:12) = mAlpha_in_n
    State(13) = mVoidRatio
    State(14) = mDGamma_n
    State(15) = mG
    State(16) = mKp


    end subroutine GetState

    !------------------------------------------------------------------------

subroutine GetAlpha(Alpha)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Alpha
    
    ! variables 
    real(REAL_TYPE) :: mAlpha_n

    ! Return alpha tensor
    Alpha = mAlpha_n

end subroutine GetAlpha

    !------------------------------------------------------------------------

    subroutine GetFabric(Fabric)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Fabric
    
    real(REAL_TYPE) :: mFabric_n 

    ! Return fabric tensor
    Fabric = mFabric_n

    end subroutine GetFabric

    !------------------------------------------------------------------------

    subroutine GetAlpha_in(Alpha_in)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in

    real(REAL_TYPE) :: mAlpha_in_n
    
    ! Return alpha_in tensor
    Alpha_in = mAlpha_in_n

    end subroutine GetAlpha_in

    !------------------------------------------------------------------------

    subroutine GetTracker(Tracker)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Tracker

    ! variables 
    real(REAL_TYPE) :: mTracker
    
    ! Return tracker vector
    Tracker = mTracker

    end subroutine GetTracker

    !------------------------------------------------------------------------


    subroutine GetKp(Kp)
    implicit none
    ! Output
    real(REAL_TYPE), intent(out) :: Kp

    ! variables 
    real(REAL_TYPE) :: mKp
    
    ! Return Kp
    Kp = mKp

    end subroutine GetKp

    !------------------------------------------------------------------------


    subroutine GetG(G)
    implicit none
    ! Output
    real(REAL_TYPE), intent(out) :: G

    ! variables 
    real(REAL_TYPE):: mG
    
    ! Return shear modulus
    G = mG

    end subroutine GetG

    !------------------------------------------------------------------------


    subroutine GetAlpha_in_p(Alpha_in_p)
    implicit none
    ! Output
    real(REAL_TYPE), dimension(3), intent(out) :: Alpha_in_p

    real(REAL_TYPE), dimension(3) :: mAlpha_in_p_n
    
    ! Return previous alpha_in tensor
    Alpha_in_p = mAlpha_in_p_n


    end subroutine GetAlpha_in_p

    !------------------------------------------------------------------------


    subroutine GetDGamma(DGamma)
    implicit none
    ! Output
    real(REAL_TYPE), intent(out) :: DGamma
    
    real(REAL_TYPE) :: mDGamma_n

    ! Return previous L
    DGamma = mDGamma_n

    end subroutine GetDGamma
    
    !------------------------------------------------------------------------

    !subroutine GetTangent(TangType, Ce, Cep, Cep_Consistent, Tangent)
    !implicit none
    !! Inputs
    !integer, intent(in) :: TangType
    !real(REAL_TYPE), dimension(3,3), intent(in) :: Ce!, Cep, Cep_Consistent
    !! Output
    !real(REAL_TYPE), dimension(3,3), intent(out) :: Tangent
    !
    !! Determine which tangent to return based on TangType
    !if (TangType == 0) then
    !    Tangent = Ce
    !else if (TangType == 1) then
    !    Tangent = Cep
    !else
    !    Tangent = Cep_Consistent
    !endif
    !
    !
    !end subroutine GetTangent


    !------------------------------------------------------------------------

    subroutine getInitialTangent(initialTangent)
    ! Input/output variables
    real(REAL_TYPE), dimension(6,6), intent(out) :: initialTangent
    
    ! variables 
    real(REAL_TYPE), dimension(6,6) :: mCe
    
    ! Return initial tangent matrix
    initialTangent = mCe
    
    end subroutine getInitialTangent

    

    !------------------------------------------------------------------------

    subroutine getStress(stress)
    ! Input/output variables
    real(REAL_TYPE), dimension(3), intent(out) :: stress
    
    ! variables 
    real(REAL_TYPE), dimension(3) :: mSigma
    real(REAL_TYPE), dimension(3) :: mSigma_b
    
    ! Calculate and return stress
    stress = mSigma + mSigma_b
    stress = stress * (-1.0)  ! -1.0 is for geotechnical sign convention
    
    end subroutine getStress

    
        
    !------------------------------------------------------------------------

    subroutine getStrain(strain)
    ! Input/output variables
    real(REAL_TYPE), dimension(3), intent(out) :: strain
    
    
    ! variables 
    real(REAL_TYPE), dimension(3) :: mEpsilon
    
    ! Calculate and return strain
    strain = mEpsilon
    strain = strain * (-1.0)  ! -1.0 is for geotechnical sign convention
    
    end subroutine getStrain
    
    !------------------------------------------------------------------------

    
    !subroutine GetStress(Stress, Sigma_r)
    !implicit none
    !! Inputs
    !real(REAL_TYPE), dimension(6), intent(in) :: Stress
    !! Output
    !real(REAL_TYPE), dimension(6), intent(out) :: Sigma_r
    !
    !! Compute stress
    !Sigma_r = Stress + mSigma_b
    !Sigma_r = Sigma_r * (-1.0)
    !
    !
    !end subroutine GetStress
    
    
    
    !-----------------------------------------------------------------------
    subroutine GetElasticStrain(Epsilon, EpsilonE_r)
    implicit none
    ! Inputs
    real(REAL_TYPE), dimension(6), intent(in) :: Epsilon
    ! Output
    real(REAL_TYPE), dimension(6), intent(out) :: EpsilonE_r

    ! Compute elastic strain
    EpsilonE_r = Epsilon * (-1.0)

    end subroutine GetElasticStrain
    !subroutine getElasticStrain(elasticStrain)
    !! Input/output variables
    !real*8, dimension(3), intent(out) :: elasticStrain
    !
    !! Calculate elastic strain
    !elasticStrain = mEpsilon
    !elasticStrain = elasticStrain * (-1.0)
    !
    !end subroutine getElasticStrain
    
    
    
    
    !------------------------------------------------------------------------------------
    !// ---------------------------------------------------------------------------------
    !/*************************************************************/
    !// Plastic Integrator
    !/*************************************************************/
    
    subroutine PM4SandIntegrate(StrainIncrement, & !strain_from_element, mEpsilon_n,
    G0, hpo, PostShake, me2p, &    
    P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
    Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp,  & 
    KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
    dEpsilonE, Sigma_n, Sigma, &
    Ce)
    
    implicit none
    
    ! input variables
    !real(REAL_TYPE), intent(in) :: Dr       ! primiary parameter 1 ! prop 1
    real(REAL_TYPE), intent(in) :: G0       ! primiary parameter 2 ! prop 2
    real(REAL_TYPE), intent(in) :: hpo ! prop 3
    
    real(REAL_TYPE), intent(in) :: PostShake ! prop 4
    real(REAL_TYPE), intent(in) :: me2p ! prop 5
 
    real(REAL_TYPE), intent(in) :: P_atm    ! statev 1
    real(REAL_TYPE), intent(in) :: h0       ! statev 2              ! secondary parameter 1
    real(REAL_TYPE), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
    real(REAL_TYPE), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
    real(REAL_TYPE), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
    real(REAL_TYPE), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
    real(REAL_TYPE), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
    real(REAL_TYPE), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
    real(REAL_TYPE), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
    real(REAL_TYPE), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
    !real(REAL_TYPE), intent(in) :: phi_cv   ! statev 11                        ! secondary parameter 9
    real(REAL_TYPE), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
    real(REAL_TYPE), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11

    real(REAL_TYPE), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12
    
    real(REAL_TYPE), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
    real(REAL_TYPE), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
    
    real(REAL_TYPE), intent(in) :: mm         ! statev 17                             ! secondary parameter 15

    
    real(REAL_TYPE), intent(in) :: Pmin       ! statev 25 
    real(REAL_TYPE), intent(in) :: Pmin2      ! statev 26

    
    real(REAL_TYPE), intent(in) :: zpeak      ! statev 28
    real(REAL_TYPE), intent(in) :: zcum       ! statev 29
    
    !real(REAL_TYPE), intent(in) :: z_max      ! statev 32
    
    !real(REAL_TYPE), intent(in) :: Ado        ! statev 33
    
    
    real(REAL_TYPE), intent(in), dimension(3) :: Fabric_n ! statev 43:45
    
    real(REAL_TYPE), intent(in), dimension(3) :: Fabric_in_n ! statev 49:51
    
    real(REAL_TYPE), intent(in), dimension(3) :: Alpha_n ! statev 55:57
    
    real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_n ! statev 61:63
    
    real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_p_n ! statev 67:69
    
    real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_true_n ! statev 73:75
    
    real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_max_n ! statev 79:81
    
    real(REAL_TYPE), intent(in), dimension(3) :: Alpha_in_min_n ! statev 85:87
    
    real(REAL_TYPE), intent(in) :: TolF ! statev 92
    
	real(REAL_TYPE), intent(in) :: TolR ! statev 93 
    
    ! output variables 
    real(REAL_TYPE), intent(out), dimension(3) :: Fabric ! statev 46:48
    real(REAL_TYPE), intent(out), dimension(3) :: Fabric_in ! statev 52:54
    real(REAL_TYPE), intent(out), dimension(3) :: Alpha ! statev 58:60

    
    real(REAL_TYPE), intent(out), dimension(3) :: dEpsilonE ! statev 94:96
    
    real(REAL_TYPE), intent(out), dimension(3) :: Sigma ! statev 101:103
    
    real(REAL_TYPE), intent(inout), dimension(3,3) :: Ce !
    
    ! inout variables
	real(REAL_TYPE), intent(in), dimension(3) :: Sigma_n ! statev 97:100                 ! current strain      
    

    ! Inout variables
    real(REAL_TYPE), dimension(3), intent(in) :: StrainIncrement!strain_from_element 
    
    real(REAL_TYPE), intent(inout) :: VoidRatio   ! statev 18

    real(REAL_TYPE), intent(inout) :: Mc         ! statev 19                      ! secondary parameter 1
    real(REAL_TYPE), intent(inout) :: Mb         ! statev 20 
    real(REAL_TYPE), intent(inout) :: Md         ! statev 21
    
    real(REAL_TYPE), intent(inout) :: Cdr        ! statev 22                       ! secondary parameter 
    
    real(REAL_TYPE), intent(inout) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
    real(REAL_TYPE), intent(inout) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
    
    real(REAL_TYPE), intent(inout) :: K_p        ! statev 27

    real(REAL_TYPE), intent(inout) :: pzp        ! statev 30 
    real(REAL_TYPE), intent(inout) :: zxp        ! statev 31
    
    real(REAL_TYPE), intent(inout) :: KK         ! statev 34
    real(REAL_TYPE), intent(inout) :: GG         ! statev 35  
    real(REAL_TYPE), intent(inout) :: Mcur       ! statev 36
    
    real(REAL_TYPE), intent(inout), dimension(3) :: rrr   ! statev 37:39
    real(REAL_TYPE), intent(inout), dimension(3) :: Sigma_b ! statev 40:42
    
    
    real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in ! statev 64:66

    real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in_p ! statev 70:72
    real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in_true ! statev 76:78
    
    
    real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in_max ! statev 82:84
    
    real(REAL_TYPE), intent(inout), dimension(3) :: Alpha_in_min ! statev 88:90
    
    !real(REAL_TYPE), intent(inout), dimension(3,3) :: Ce
    
    !integer(INTEGER_TYPE), intent(inout) :: pzpFlag ! statev ! statev 91
    real(REAL_TYPE), intent(inout) :: pzpFlag ! statev ! statev 91
    
    ! Local variables
    real(REAL_TYPE), dimension(3) :: n_tr
    real(REAL_TYPE), dimension(3) :: tmp0
    real(REAL_TYPE), dimension(3) :: tmp1
    real(REAL_TYPE), dimension(3) :: Alpha_Alpha_in_true
    
    real(REAL_TYPE) :: NextVoidRatio ! this needs to either be intent(out) or overwrie VoidRatio
    real(REAL_TYPE) :: pp
    real(REAL_TYPE) :: zxpTemp
    real(REAL_TYPE) :: dGamma
    
    integer(INTEGER_TYPE) :: ii
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE) :: GetNorm_Contr_result

    
    
    
    
    ! Assignments, new = old
    Alpha = Alpha_n ! set current alpha as the previous alpha
    Alpha_in = Alpha_in_n ! set initial alpha as the previous initial alpha !--> EFFECTIVE INITIAL
    Alpha_in_true = Alpha_in_true_n ! true initial
    Alpha_in_p = Alpha_in_p_n ! previous initial
    Alpha_in_max = Alpha_in_max_n ! maximum --> used for negative loading 
    Alpha_in_min = Alpha_in_min_n ! minimum --> used for positive loading
    Fabric = Fabric_n ! set the current fabric as the previous fabric
    Fabric_in = Fabric_in_n ! set current initial fabric as the previous initial fabric

    ! Calculate trial stress
    Sigma = Sigma_n + matmul(Ce, StrainIncrement)!the third component has got to be the GAMMMA ENGINEERING SHEAR STRAIN because 
    ! using the current mAlpha and the trial stress tmp0, get the normal to the yield surface in the stress ratio space
    
    ! find the normal to the yield surface
    call GetNormalToYield(Sigma, Alpha, n_tr) 
    
    
    ! Check loading reversal condition
    ! if ((DoubleDot2_2_Contr(mAlpha - mAlpha_in_true, n_tr) < 0.0) && me2p) {
    ! find the difference between current and initial stress ratios
    Alpha_Alpha_in_true = Alpha - Alpha_in_true
    
    ! dot the back stress ratio difference with the normal 
    call DoubleDot2_2_Contr(Alpha_Alpha_in_true, n_tr, DoubleDot2_2_Contr_result) ! daxn -->dot it with the n_tr --> look at table 3.1
    if ( (DoubleDot2_2_Contr_result < 0.0) .and. (me2p /= .false.) ) then !--> me2p will actually be 1
        
        ! this is a load reversal
        Alpha_in_true = Alpha
        Alpha_in_p = Alpha_in ! --> previous initial becomes the initial 
        Fabric_in = Fabric ! --> initial fabric set to current due to this being a load reversal
        
        ! Update pzp: mean effective stress at the time that zp achieves its greates value
        call GetTrace(Sigma_n, pp)
        pp = 0.5 * pp
        
        if (pp <= Pmin) then
            pp = Pmin
        endif
        
        call GetNorm_Contr(Fabric_n, zxpTemp)
        zxpTemp = zxpTemp * pp
        
        ! record the maximum value of zp and its corresponding mean effective stress
        if ( ((zxpTemp > zxp) .and. (pp > pzp)) .or. pzpFlag==1.0 ) then
            zxp = zxpTemp
            pzp = pp
            pzpFlag = 0.0!.false.
        end if
        
        ! Track initial back-stress ratio history
        do ii = 1, 3
            if (Alpha_in(ii) > 0.0) then ! --> if positive 
                ! Minimum positive value
                Alpha_in_min(ii) = min(Alpha_in_min(ii), Alpha(ii)) ! minimum back-stress ratio
            else ! --> if negative
                ! Maximum negative value
                Alpha_in_max(ii) = max(Alpha_in_max(ii), Alpha(ii)) ! maximum back-stress ratio
            endif
        enddo
        
        ! Update mAlpha_in based on loading direction
        if (Alpha(3) * Alpha_in_p(3) > 0) then
            
            do ii = 1, 3
                if (n_tr(ii) > 0.0) then
                    ! Positive loading direction
                    Alpha_in(ii) = max(0.0, Alpha_in_min(ii))
                else
                    ! Negative loading direction
                    Alpha_in(ii) = min(0.0, Alpha_in_max(ii))
                endif
            enddo
            
        else
            
            ! when would this happen?
            ! alpha_xy spans the positive and negative quadrants then avoid the min/max business
            Alpha_in = Alpha
        endif
    
    end if

    ! Force elastic response
    if (me2p == .false.) then !--> this will never be the case --> probably for testing the model assuming elastic behavior
        !subroutine elastic_integrator(CurStress, NextStrainIncrement, G0, nu, Pmin, P_atm, & ! input
        !NextElasticStrainIncrement, NextStress, NextAlpha, NextVoidRatio, GG, KK, aC) ! output 
        call elastic_integrator(Sigma_n, StrainIncrement, G0, nu, Pmin, P_atm, &
                            dEpsilonE, Sigma, Alpha, GG, KK, Ce, VoidRatio, NextVoidRatio)
         !--> maybe to initialize stresses in an elastic manner without changing constitutive model
    else
        ! me2p -- .true.
        ! ElastoPlastic response --> always go here
        ! Explicit schemes
        call explicit_integrator(StrainIncrement, Sigma_n, Alpha_n, Fabric_n, Alpha_in, & !mEpsilon_n,
                Alpha_in_p, dEpsilonE, Sigma, Alpha, Fabric, DGamma, NextVoidRatio, &
                GG, KK, Ce, &
                TolF,  &
                G0, hpo, PostShake, me2p, &    
                P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, & !phi_cv, 
                Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, &
                Fabric_in, &
                MCur)!, Cep, Cep_Consistent)
    end if

end subroutine PM4SandIntegrate
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !------------------------------------------------------------------------
    !/*************************************************************/
    !// Elastic Integrator
    !/*************************************************************/
    
    subroutine elastic_integrator(CurStress, NextStrainIncrement, G0, nu, Pmin, P_atm, & ! input
        NextElasticStrainIncrement, NextStress, NextAlpha, GG, KK, aC, VoidRatio, NextVoidRatio) ! output 
    
    ! assuming linear elasticity for stress initialization
    implicit none
    
    ! input variables
    real(REAL_TYPE), dimension(3), intent(in) :: NextStrainIncrement
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress
    real(REAL_TYPE), intent(in) :: G0
    real(REAL_TYPE), intent(in) :: nu
    real(REAL_TYPE), intent(in) :: Pmin
    real(REAL_TYPE), intent(in) :: P_atm
    
    real(REAL_TYPE), intent(in) :: VoidRatio

    ! Output variables
    real(REAL_TYPE), intent(out) :: GG
    real(REAL_TYPE), intent(out) :: KK 
    real(REAL_TYPE), dimension(3), intent(out) :: NextElasticStrainIncrement
    real(REAL_TYPE), dimension(3, 3), intent(out) :: aC
    
    real(REAL_TYPE), dimension(3), intent(out) :: NextStress
    real(REAL_TYPE), dimension(3), intent(out) :: NextAlpha
    
    real(REAL_TYPE), intent(out) :: NextVoidRatio
    
    ! Local variables
    real(REAL_TYPE) :: pp
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
    real(REAL_TYPE), dimension(3) :: DetDevPart_result
    
    ! Calculate the next elastic strain increment (assume all elastic)
    NextElasticStrainIncrement = NextStrainIncrement
    
    ! Calculate void ratio 
    call GetTrace(NextElasticStrainIncrement, GetTrace_result)
    NextVoidRatio = VoidRatio - ( (1 + VoidRatio) * GetTrace_result )
    
    ! Calculate elastic moduli
    call GetElasticModuli(CurStress, G0, nu, P_atm, Pmin, GG, KK)
    
    ! Calculate elastic matrix (3x3)
    call GetStiffness(KK, GG, aC)
    
    ! Calculate the next stress
    call DoubleDot4_2(aC, NextStrainIncrement, DoubleDot4_2_result)
    NextStress = CurStress + DoubleDot4_2_result
    
    ! Calculate mean effective stress
    call GetTrace(NextStress, pp)
    pp = 0.5*pp !sum(NextStress)
    
    ! just to avoid dividing by zero in case pp is very small 
    ! (fyi won't really harden if in elastic integrator state)
    if (pp > Pmin) then 
        call GetDevPart(NextStress, DetDevPart_result)
        NextAlpha = DetDevPart_result / pp
    end if

    end subroutine elastic_integrator
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !-------------------------------------------------------------------------
    !// ----------------------------------------------------------------------
    !/*************************************************************/
    !// Explicit Integrator
    !/*************************************************************/
    
    subroutine explicit_integrator(StrainIncrement, CurStress, CurAlpha, CurFabric, alpha_in, alpha_in_p,  & !CurStrain, NextStrain,CurElasticStrain,
                            NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, Ce, TolF, &
                            G0, hpo, PostShake, me2p, &    
                            P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
                            Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, &
                            Fabric_in, &
                            MCur) !, aCep, aCep_Consistent)
    
    ! explicit elastoplastic integrator that should work with ForwardEuler and ModifiedEuler
    
    implicit none
    
    ! Input variables
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress
    real(REAL_TYPE), dimension(3), intent(in) :: CurAlpha
    real(REAL_TYPE), dimension(3), intent(in) :: CurFabric
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in_p
    real(REAL_TYPE), dimension(3), intent(in) :: StrainIncrement
    real(REAL_TYPE), intent(in) :: TolF
    real(REAL_TYPE), intent(in) :: VoidRatio
    
    !real(REAL_TYPE), intent(in) :: emax 
    !real(REAL_TYPE), intent(in) :: emin 
    !real(REAL_TYPE), intent(in) :: zpeak 
    !real(REAL_TYPE), intent(in) :: zcum 
    real(REAL_TYPE), intent(inout) :: MCur
    
    real(REAL_TYPE), dimension(3), intent(in) :: Fabric_in

    !real(REAL_TYPE), intent(in) :: Pmin
    !real(REAL_TYPE), intent(in) :: QQ_Bolton 
    !real(REAL_TYPE), intent(in) :: RR_bolton
    !
    !real(REAL_TYPE), intent(in) :: cz
    
    
    !!!!!!!!!!!!!!!!!
    
    !real(REAL_TYPE), intent(in) :: Dr       ! primiary parameter 1 ! prop 1
    real(REAL_TYPE), intent(in) :: G0       ! primiary parameter 2 ! prop 2
    real(REAL_TYPE), intent(in) :: hpo ! prop 3
    
    real(REAL_TYPE), intent(in) :: PostShake ! prop 4
    real(REAL_TYPE), intent(in) :: me2p ! prop 5
 
    real(REAL_TYPE), intent(in) :: P_atm    ! statev 1
    real(REAL_TYPE), intent(in) :: h0       ! statev 2              ! secondary parameter 1
    real(REAL_TYPE), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
    real(REAL_TYPE), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
    real(REAL_TYPE), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
    real(REAL_TYPE), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
    real(REAL_TYPE), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
    real(REAL_TYPE), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
    real(REAL_TYPE), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
    real(REAL_TYPE), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
    !real(REAL_TYPE), intent(in) :: phi_cv   ! statev 11                        ! secondary parameter 9
    real(REAL_TYPE), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
    real(REAL_TYPE), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11
    
    real(REAL_TYPE), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12
    
    real(REAL_TYPE), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
    real(REAL_TYPE), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
    
    real(REAL_TYPE), intent(in) :: mm         ! statev 17                             ! secondary parameter 15

    
    real(REAL_TYPE), intent(in) :: Pmin       ! statev 25 
    real(REAL_TYPE), intent(in) :: Pmin2      ! statev 26

    
    real(REAL_TYPE), intent(in) :: zpeak      ! statev 28
    real(REAL_TYPE), intent(in) :: zcum       ! statev 29
    
    !real(REAL_TYPE), intent(in) :: z_max      ! statev 32
    
    !real(REAL_TYPE), intent(in) :: Ado        ! statev 33
    
    real(REAL_TYPE), intent(inout) :: Mc         ! statev 19                      ! secondary parameter 1
    real(REAL_TYPE), intent(inout) :: Mb         ! statev 20 
    real(REAL_TYPE), intent(inout) :: Md         ! statev 21
    
    real(REAL_TYPE), intent(inout) :: Cdr        ! statev 22                       ! secondary parameter 
    
    real(REAL_TYPE), intent(inout) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
    real(REAL_TYPE), intent(inout) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
    
    real(REAL_TYPE), intent(inout) :: K_p        ! statev 27

    
    ! Output variables
    real(REAL_TYPE), dimension(3), intent(out) :: NextElasticStrainIncrement
    real(REAL_TYPE), dimension(3), intent(out) :: NextStress
    real(REAL_TYPE), dimension(3), intent(out) :: NextAlpha
    real(REAL_TYPE), dimension(3), intent(out) :: NextFabric
    real(REAL_TYPE), dimension(3, 3), intent(out) :: Ce !out as it is calculated here
    real(REAL_TYPE), intent(out) :: NextL
    real(REAL_TYPE), intent(out) :: NextVoidRatio
    
    ! inout variables
    real(REAL_TYPE), intent(inout) :: GG
    real(REAL_TYPE), intent(inout) :: KK
    
    ! Local variables
    real(REAL_TYPE) :: a0
    real(REAL_TYPE) :: a1
    real(REAL_TYPE) :: elasticRatio ! previously named aa
    
    real(REAL_TYPE) :: ff
    real(REAL_TYPE) :: fn
    real(REAL_TYPE) :: dVolStrain
    
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: GetNorm_Contr_result
    
    real(REAL_TYPE), dimension(3) :: dSigma
    real(REAL_TYPE), dimension(3) :: dDevStrain 
    real(REAL_TYPE), dimension(3) :: nn
    real(REAL_TYPE), dimension(3) :: I1
    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
    
    
    
    ! establish numerical values
    I1(1) = 1
    I1(2) = 1
    I1(3) = 0
    
    ! calculate the next void ratio using the strain
    call GetTrace(StrainIncrement, dVolStrain) ! trace to find volumetric strain increment
    NextVoidRatio = VoidRatio - ( (1 + VoidRatio) * dVolStrain )
    
    ! calculate the strain increment assuming fully elastic strain increment
    NextElasticStrainIncrement = StrainIncrement
    ! calculate volumetric strain and deviatoric strain increments
    !dVolStrain = dVolStrain
    dDevStrain = StrainIncrement - (dVolStrain / 3.0) ! deviatoric strain increment
    
    ! calculate 3x3 stiffness matrix (aC)
    call GetStiffness(KK, GG, Ce) 

    ! calculate the stress increment
    dSigma = (KK * dVolStrain * I1) + (2 * GG * dDevStrain)
    
    ! Update the stress assuming an entirely elastic stress increment
    NextStress = CurStress + dSigma
    
    ! Calculate some parameters for plasticity
    !-------------------------------------------------------------------------
    !ff = GetF(NextStress, CurAlpha)
    ! Correct the stress
    ! I added this so that we can have the correct considertion of the two shear components
    call ToCovariant(NextStress, ToCovariant_result)
    call GetFYieldFunction(ToCovariant_result, CurAlpha, ff) ! ff is the final yield function evaluation
    !-------------------------------------------------------------------------
    !fn = GetF(CurStress, CurAlpha)
    ! Correct the stress
    ! I added this so that we can have the correct considertion of the two shear components
    call ToCovariant(CurStress, ToCovariant_result)
    call GetFYieldFunction(ToCovariant_result, CurAlpha, fn) ! fn is the initial yield function evaluation
    !-------------------------------------------------------------------------
    !nn = GetNormalToYield(NextStress, CurAlpha)
    call GetNormalToYield(NextStress, CurAlpha, nn) ! nn is the normal of the yield function at the projected stress
    !-------------------------------------------------------------------------
    ! Perform plasticity calculations
    if (ff <= TolF) then ! Pure elastic step (no transition)
        
        ! Pure elastic loading/unloading
        NextAlpha = CurAlpha
        NextFabric = CurFabric
        NextL = 0
        
        return
        
    elseif (fn < -TolF) then ! transition from elastic to plastic
        
        ! establish the inputs to the intersection factor
        a0 = 0.0
        a1 = 1.0 
        
        ! evaluate the intersection factor aa        
        call IntersectionFactor(CurStress, StrainIncrement, CurAlpha, Ce, a0, a1, elasticRatio, TolF) 
        
        NextElasticStrainIncrement = elasticRatio * StrainIncrement ! elastic portion of the strain
        call DoubleDot4_2(Ce, NextElasticStrainIncrement, dSigma) ! elastic portion of the stress
        
        ! Temporary commented exp_int
        !call exp_int(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain, &
        !    CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
        !    NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
        
        !call ForwardEuler(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain, &
        !    CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
        !    NextL, NextVoidRatio, GG, KK, aC)!, aCep, aCep_Consistent)
        
        ! input stresses correspond to the onset of plasticity
        call ForwardEuler(CurStress + dSigma, StrainIncrement-NextElasticStrainIncrement,  &
                          CurAlpha, CurFabric, alpha_in, alpha_in_p,  &
                            G0, hpo, PostShake, me2p, &    
                            P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
                            Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, &
                            Fabric_in, &
                          NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric, &
                          NextL, NextVoidRatio, GG, KK, Ce, &
                            MCur) !, aCep, aCep_Consistent)
        
        return
        
        !ForwardEuler(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, &
        !                         NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
    else if (abs(fn) < TolF) then ! Pure plastic step or elastic unloading followed by plastic loading
        
        call GetNormalToYield(CurStress, CurAlpha, nn)
        call DoubleDot2_2_Contr(nn, dSigma, DoubleDot2_2_Contr_result)
        call GetNorm_Contr(dSigma, GetNorm_Contr_result)
        
        !if (DoubleDot2_2_Contr(GetNormalToYield(CurStress, CurAlpha), dSigma) / (max(1.0, GetNorm_Contr(dSigma))) > (-sqrt(mTolF))) then
        !cos(theta) = (nn dot dSigma) / (|nn| |dSigma|) < 0 --> negative means angle is greater than 90 so unloading is occur
        
        if (DoubleDot2_2_Contr_result / (max(1.0, GetNorm_Contr_result)) > (-sqrt(TolF))) then ! if true then no unloading Sloan et al. 2001
            ! Pure plastic step
            
            !subroutine ForwardEuler(CurStress, StrainIncrement, CurAlpha, CurFabric, alpha_in, alpha_in_p,  & !CurStrain, NextStrain, CurElasticStrain,
            !                emax, emin, zpeak, zcum, MCur, Fabric_in, VoidRatio, Pmin, &     
            !                QQ_Bolton, RR_bolton, &
            !                NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, aC)!, aCep, aCep_Consistent) ! do we need all three of aC, aCep, aCep_Consistent
            !
            
            ! Temporary commented exp_int
            call ForwardEuler(CurStress, StrainIncrement, &!CurStrain, CurElasticStrain, &
                CurAlpha, CurFabric, alpha_in, alpha_in_p, &
                G0, hpo, PostShake, me2p, &    
                P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, &
                Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, &
                Fabric_in, &
                NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric, &
                NextL, NextVoidRatio, GG, KK, Ce, &
                            MCur) !, aCep, aCep_Consistent)
            
            return 
            !call exp_int(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, &
            !    NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
            
        else
            
            ! Elastic unloading followed by plastic loading
            call IntersectionFactor_Unloading(CurStress, StrainIncrement, CurAlpha, Ce, TolF, elasticRatio) !CurStrain, NextStrain,
            
            NextElasticStrainIncrement = elasticRatio * StrainIncrement !r atio*(NextStrain - CurStrain)
            call DoubleDot4_2(Ce, NextElasticStrainIncrement, dSigma)
            
            ! Temporary commented exp_int
            !call exp_int(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain, &
            !    CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
            !    NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
            
            call ForwardEuler(CurStress + dSigma, StrainIncrement-NextElasticStrainIncrement, &!CurStrain + dElasStrain, CurElasticStrain + dElasStrain, &
                                CurAlpha, CurFabric, alpha_in, alpha_in_p, &
                                G0, hpo, PostShake, me2p, &    
                            P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, & !phi_cv, 
                            Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, &
                            Fabric_in, &
                                NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric, & !NextStrain, 
                                NextL, NextVoidRatio, GG, KK, Ce, &
                            MCur)
            
            return 
            
        end if
        
    else
        
        ! Correct the stress
        ! I added this so that we can have the correct considertion of the two shear components
        call ToCovariant(CurStress, ToCovariant_result)
        
        call GetFYieldFunction(ToCovariant_result, CurAlpha,ff)
        
        ! Illegal stress state
        print *, "PM4Sand: Encountered an illegal stress state! Tag: "!, !getTag()
        print *, "            f = ", ff
        ! Temporary commented exp_int
        !call exp_int(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, &
        !    NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
    
    end if

end subroutine explicit_integrator
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !-------------------------------------------------------------------------
    !// ----------------------------------------------------------------------
    !/*************************************************************/
    !// Forward-Euler Integrator
    !/*************************************************************/
    subroutine ForwardEuler(CurStress, StrainIncrement, CurAlpha, CurFabric, alpha_in, alpha_in_p,  & !CurStrain, NextStrain, CurElasticStrain,
                            G0, hpo, PostShake, me2p, &    
                            P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, VoidRatio, & !phi_cv, 
                            Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, &
                            Fabric_in, &
                            NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, GG, KK, Ce, &
        MCur)!, aCep, aCep_Consistent) ! do we need all three of aC, aCep, aCep_Consistent
    
    
    implicit none
    
    ! input variables
    real(REAL_TYPE), intent(in), dimension(3) :: CurStress 
    real(REAL_TYPE), intent(in), dimension(3) :: StrainIncrement
    
    !real(REAL_TYPE), intent(in), dimension(3) :: CurElasticStrain
    real(REAL_TYPE), intent(in), dimension(3) :: CurAlpha 
    real(REAL_TYPE), intent(in), dimension(3) :: CurFabric 
    real(REAL_TYPE), intent(in), dimension(3) :: alpha_in 
    real(REAL_TYPE), intent(in), dimension(3) :: alpha_in_p 
    
    
    
    ! input variables
    !real(REAL_TYPE), intent(in) :: Dr       ! primiary parameter 1 ! prop 1
    real(REAL_TYPE), intent(in) :: G0       ! primiary parameter 2 ! prop 2
    real(REAL_TYPE), intent(in) :: hpo ! prop 3
    
    real(REAL_TYPE), intent(in) :: PostShake ! prop 4
    real(REAL_TYPE), intent(in) :: me2p ! prop 5
 
    real(REAL_TYPE), intent(in) :: P_atm    ! statev 1
    real(REAL_TYPE), intent(in) :: h0       ! statev 2              ! secondary parameter 1
    real(REAL_TYPE), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
    real(REAL_TYPE), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
    real(REAL_TYPE), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
    real(REAL_TYPE), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
    real(REAL_TYPE), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
    real(REAL_TYPE), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
    real(REAL_TYPE), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
    real(REAL_TYPE), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
    !real(REAL_TYPE), intent(in) :: phi_cv   ! statev 11                        ! secondary parameter 9
    real(REAL_TYPE), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
    real(REAL_TYPE), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11

    real(REAL_TYPE), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12

    real(REAL_TYPE), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
    real(REAL_TYPE), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
    
    real(REAL_TYPE), intent(in) :: mm         ! statev 17                             ! secondary parameter 15

    
    real(REAL_TYPE), intent(in) :: Pmin       ! statev 25 
    real(REAL_TYPE), intent(in) :: Pmin2      ! statev 26

    
    real(REAL_TYPE), intent(in) :: zpeak      ! statev 28
    real(REAL_TYPE), intent(in) :: zcum       ! statev 29
    
    !real(REAL_TYPE), intent(in) :: z_max      ! statev 32
    
    
    
    
    
    
    
    !real(REAL_TYPE), intent(in) :: emax
    !real(REAL_TYPE), intent(in) :: emin
    !
    !real(REAL_TYPE), intent(in) :: zpeak
    !real(REAL_TYPE), intent(in) :: zcum
    !real(REAL_TYPE), intent(in) :: MCur
    
    real(REAL_TYPE), intent(in), dimension(3) :: Fabric_in
    
    real(REAL_TYPE), intent(in) :: VoidRatio
    
    
    
    ! this needs to be intent(in)
    real(REAL_TYPE), dimension(3) :: Alpha_in_true ! statev 73:75 !intent(in),
    
    
    
    ! output variables
    real(REAL_TYPE), intent(out), dimension(3) :: NextElasticStrainIncrement
    
    real(REAL_TYPE), intent(out), dimension(3) :: NextAlpha
    real(REAL_TYPE), intent(out), dimension(3) :: NextFabric
    real(REAL_TYPE), intent(out) :: NextL
    real(REAL_TYPE), intent(out) :: NextVoidRatio
    real(REAL_TYPE), intent(out) :: GG
    real(REAL_TYPE), intent(out) :: KK
    real(REAL_TYPE), intent(out), dimension(3,3) :: Ce
    !real(REAL_TYPE), intent(out), dimension(3,3) :: aCep
    !real(REAL_TYPE), intent(out), dimension(3,3) :: aCep_Consistent
    
    ! inout variables
    real(REAL_TYPE), intent(inout), dimension(3) :: NextStress
    
    
    
    
    real(REAL_TYPE), intent(inout) :: Mc         ! statev 19                      ! secondary parameter 1
    real(REAL_TYPE), intent(inout) :: Mb         ! statev 20 
    real(REAL_TYPE), intent(inout) :: Md         ! statev 21
    
    real(REAL_TYPE), intent(inout) :: Cdr        ! statev 22                       ! secondary parameter 
    
    real(REAL_TYPE), intent(inout) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
    real(REAL_TYPE), intent(inout) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
    
        
    real(REAL_TYPE), intent(inout) :: Mcur       ! statev 36

    ! Local variables
    !real(REAL_TYPE) :: CurVoidRatio
    real(REAL_TYPE) :: CurDr
    real(REAL_TYPE) :: Cka
    real(REAL_TYPE) :: hh
    real(REAL_TYPE) :: pp
    real(REAL_TYPE) :: dVolStrain
    real(REAL_TYPE) :: DD
    real(REAL_TYPE) :: AlphaAlphaBDotN
    
    real(REAL_TYPE), dimension(3) :: nn
    real(REAL_TYPE), dimension(3) :: RR
    real(REAL_TYPE), dimension(3) :: alphaD
    
    real(REAL_TYPE), dimension(3) :: alphaB
    
    real(REAL_TYPE), dimension(3) :: dPStrain
    real(REAL_TYPE), dimension(3) :: bb
    real(REAL_TYPE), dimension(3) :: dDevStrain
    real(REAL_TYPE), dimension(3) :: rrr
    real(REAL_TYPE), dimension(3) :: dSigma
    real(REAL_TYPE), dimension(3) :: dAlpha
    real(REAL_TYPE), dimension(3) :: dFabric    
    
    
    real(REAL_TYPE), dimension(3) :: dSigmaP
    
    !real(REAL_TYPE) :: pp
    
    real(REAL_TYPE) :: pzp
    real(REAL_TYPE) :: K_p
    real(REAL_TYPE) :: DGamma
    integer(INTEGER_TYPE) :: DebugFlag
    
    real(REAL_TYPE) :: ksi
    real(REAL_TYPE) :: TolF
    
    ! temporary local variables 
    real(REAL_TYPE) :: temp4 !--> this temporary variable is used to check the denominator of L
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: DoubleDot2_2_Mixed_result
    real(REAL_TYPE) :: Macauley_result
    real(REAL_TYPE) :: Macauley_result2
    real(REAL_TYPE) :: MacauleyIndex_result
    real(REAL_TYPE) :: GetNorm_Contr_result
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    real(REAL_TYPE), dimension(3) :: ToContraviant_result
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
	
    ! local variables with numeric values 
    integer(INTEGER_TYPE), dimension(3) :: I1 ! identity matrix
    real(REAL_TYPE) :: two3 ! two thirds
    
    
    ! establish numeric values
    I1(1) = 1
    I1(2) = 1
    I1(3) = 0
    
    two3 = 0.666666666666666667
    
    ! Get elastic moduli
    ! you need 'NextStress' because of fabric aspects
    call GetElasticModuli_(NextStress, zcum, z_max, nu, G0, Md, Mb, &
        PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, Fsed_min, me2p) ! --> csr uses mMcur
    ! Calculate current void ratio
    !call GetTrace(CurStrain, GetTrace_result)
    !call GetTrace(StrainIncrement, GetTrace_result)
    ! calculate the current void ratio where m_e_init remains the same for a material point
    !CurVoidRatio = VoidRatio - (1.0 + VoidRatio) * GetTrace_result !sum(CurStrain) 
    ! m_e_init needs to be global to the entire PM4Sand module 
    ! Calculate current Dr
    CurDr = (emax - VoidRatio) / (emax - emin) 
    ! m_emax and m_emin needs to be global the entire PM4Sand module 
    
    ! Calculate mean effective stress
    call GetTrace(CurStress, pp)
    pp = 0.5 * pp!sum(CurStress)
    if (pp < Pmin) pp = Pmin ! Apply tension cutoff
    
    ! Calculate void ratio from NextStrain
    call GetTrace(StrainIncrement, dVolStrain) ! volumetric strain
    NextVoidRatio = VoidRatio - ( (1.0 + VoidRatio) * dVolStrain ) !sum(NextStrain)
    ! Calculate strain increment
    !// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
    !dStrain = NextStrain - CurStrain
    
    ! Calculate NextElasticStrain
    NextElasticStrainIncrement = StrainIncrement!dStrain ! CurElasticStrain +
    !// using NextStress instead of CurStress to get correct n
	!// 1: calculate state parameters and initialize them
    
    ! Calculate state parameters
    call GetStateDependent(NextStress, CurAlpha, alpha_in, alpha_in_p, alpha_in_true, CurFabric, Fabric_in, &
        GG, zcum, zpeak, pzp, Mcur, CurDr, &
        Mc, nd, nb, Pmin, Pmin2, P_atm, & 
        mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
        QQ_Bolton, RR_Bolton, &
        nn, alphaB, alphaD, bb, RR, &
        Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, &
        ksi)
    
    
    ! Calculate volumetric strain increment
    call GetTrace(StrainIncrement, dVolStrain)
    !// dVolStrain = GetTrace(NextStrain - CurStrain);
    !dVolStrain = GetTrace_result!sum(dStrain)
    ! Calculate deviatoric strain increment
	!// dDevStrain = (NextStrain - CurStrain) - dVolStrain / 3.0 * mI1;
    dDevStrain = StrainIncrement - (dVolStrain / 3.0)
    ! Calculate stress ratio
    call GetDevPart(NextStress, GetDevPart_result)
    rrr = GetDevPart_result / pp
    ! Calculate denominator term
    !double precision :: temp4
    call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
    temp4 = K_p + (2.0 * GG) - (KK * DD * DoubleDot2_2_Contr_result) !DoubleDot2_2_Contr(nn, rrr)
    !// if (temp4 < 0.0) {
	!// 	mKp = -0.5 * (2 * G - K* D *DoubleDot2_2_Contr(n, r));
	!// 	temp4 = mKp + 2 * G - K* D *DoubleDot2_2_Contr(n, r);
	!// 	h = 1.5 * mKp / (p * AlphaAlphaBDotN);
	!// }
    
    ! Check denominator
    if (abs(temp4) < small) then !--> prevent division by zero
    
        ! Neutral loading
        dSigma = 0.0
        dAlpha = 0.0
        dFabric = 0.0
        ! dPStrain = dDevStrain + dVolStrain * mI1
        dPStrain = StrainIncrement !--> assume entire strain increment is elastic
        
    else
    
        ! Calculate L
        call DoubleDot2_2_Mixed(nn, dDevStrain, DoubleDot2_2_Mixed_result)
        call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
        !NextL = (2.0 * GG * DoubleDot2_2_Mixed(nn, dDevStrain) - DoubleDot2_2_Contr(nn, rrr) * KK * dVolStrain) / temp4 ! Equation 31
        NextL = ( (2.0 * GG * DoubleDot2_2_Mixed_result) -  (DoubleDot2_2_Contr_result * KK * dVolStrain) ) / temp4 ! Equation 31
        
        ! Set mDGamma to NextL
        DGamma = NextL ! mDGamma is the plastic mulitplier
        if (NextL < 0.0) then
            
            ! If NextL is negative
            if (debugFlag) then
                write(*,*) "NextL is smaller than 0"
                write(*,*) "NextL = ", NextL
            end if
            
            call ToContraviant(dDevStrain, ToContraviant_result)
            dSigma = (2.0 * GG * ToContraviant_result) + (KK * dVolStrain * I1) !--> focusing only on the elastic increment if L is negative
            dAlpha = 0.0
            dFabric = 0.0
            dPStrain = 0.0
    
        else
            
            ! If NextL is non-negative
            ! Calculate dSigma
            call ToContraviant(dDevStrain, ToContraviant_result)
            call Macauley(NextL, Macauley_result)
            !dSigma = 2.0 * GG * ToContraviant(dDevStrain) + KK * dVolStrain * mI1 - Macauley(NextL) * &
            !     (2.0 * GG * nn + KK * DD * mI1)
            dSigma = (2.0 * GG * ToContraviant_result) + (KK * dVolStrain * I1)  & ! elastic component
                - ( Macauley_result * ((2.0 * GG * nn) + (KK * DD * I1)) ) ! plastic component
            
            ! Update fabric
            call DoubleDot2_2_Contr(alphaD - CurAlpha, nn, DoubleDot2_2_Contr_result)
            
            if (DoubleDot2_2_Contr_result < 0.0) then
                ! dilation
                
                ! Update fabric according to Equation 57
                dFabric = (z_max * nn) + CurFabric
                
                call Macauley( (zcum / (2.0*z_max)) - 1.0, Macauley_result)
                call Macauley(NextL, Macauley_result2)
                call MacauleyIndex(-DD, MacauleyIndex_result) ! DD is negative if dilation
                
                !//  dz in Equation 57 
                dFabric = (-1.0 * (cz / (1 + Macauley_result)) * Macauley_result2 * &
                       MacauleyIndex_result * (1/DD) * dFabric) ! I included the dilatancy aspect in this equation because DD was not included???
            end if
            
            ! Update alpha
            !// dAlpha = two3 * NextL * h * b;
			!// dPStrain = NextL * mIIco * R;
            call ToCovariant(RR, ToCovariant_result)
            
            dPStrain = ToCovariant_result * NextL ! --> plastic strain (3x1)
            dAlpha = bb * (two3 * NextL * hh) ! --> change in alpha (Hardening) (3x1)
            
        end if
    
    end if


    ! Update NextFabric, NextElasticStrain, NextStress, and NextAlpha
    !// NextFabric = CurFabric + dFabric;
    NextFabric = CurFabric + dFabric
    !// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
    NextElasticStrainIncrement = StrainIncrement - dPStrain ! remove plastic strain increment !NextStrain - CurStrain !CurElasticStrain +
    !// NextStress = CurStress + dSigma;
    NextStress = CurStress + dSigma 
    !// NextAlpha = CurAlpha + dAlpha;
    NextAlpha = CurAlpha + dAlpha
    
    
    ! Perform stress correction
    !call Stress_Correction_(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio)
    call stress_correction_(alpha_in, alpha_in_p, alpha_in_true, Fabric_in, CurFabric, CurStress, Pmin, P_atm, &
                             TolF, Mc, emax, emin, zcum, zpeak, z_max, pzp, mm, h0, hpo, QQ_Bolton, RR_Bolton, &
                             DGamma, alphaD, alphaB, dSigmaP, DD, Mb, Md, K_p, NextVoidRatio, Mcur, GG, KK, NextStress, &
                             NextAlpha, dSigma, &
                             nb, nd, &
                             Pmin2, &
                             Ado)
    
    end subroutine ForwardEuler

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !--------------------------------------------------------------------------------------
    !// -------------------------------------------------------------------------------------------------------
    !/*************************************************************/
    !// Integrator Constraining Maximum Strain Increment
    !/*************************************************************/
    
    subroutine MaxStrainInc(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, &
        GG, KK, aC)!, aCep, aCep_Consistent)
    
    implicit none
    
    ! Inputs
    real(REAL_TYPE), intent(in) :: CurStress(3), CurStrain(3), CurElasticStrain(3), CurAlpha(3), CurFabric(3), alpha_in(3), alpha_in_p(3), NextStrain(3)
    ! Outputs
    real(REAL_TYPE), intent(out) :: NextElasticStrain(3), NextStress(3), NextAlpha(3), NextFabric(3), NextL, NextVoidRatio, GG, KK
    real(REAL_TYPE), intent(out) :: aC(3,3)!, aCep(3,3), aCep_Consistent(3,3)
    
    ! Function pointer to the integration scheme
    ! Temporary commented exp_int
    !procedure(PM4Sand_IntegrationScheme), pointer :: exp_int !--> more work
    
    ! Local variables
    real(REAL_TYPE) :: StrainInc(3), maxInc
    integer(INTEGER_TYPE) :: ii, numSteps
    real(REAL_TYPE) :: cStress(3), cStrain(3), cAlpha(3), cFabric(3), cAlpha_in(3), cAlpha_in_p(3), cEStrain(3)
    real(REAL_TYPE) :: nStrain(3)
    real(REAL_TYPE), dimension(3,3) :: nCe, nCep, nCepC
    
    real(REAL_TYPE) :: maxStrainIncValue
    
    maxStrainIncValue = 1e-5 !maxStrainInc

    ! Assign function pointer based on integration scheme
    ! Temporary commented exp_int
    !select case (mScheme)
    !    case (INT_MAXSTR_FE)
    !        exp_int => ForwardEuler !PM4Sand_
    !    case (INT_MAXSTR_ME)
    !        exp_int => ModifiedEuler !PM4Sand_
    !    case default
    !        exp_int => ModifiedEuler !PM4Sand_
    !end select
    
    ! Compute maximum strain increment
    StrainInc = NextStrain - CurStrain
    maxInc = StrainInc(1)
    do ii = 2, 3
        if (abs(StrainInc(ii)) > abs(maxInc)) maxInc = StrainInc(ii)
    end do
    
    ! Apply maximum strain increment constraint
    if (abs(maxInc) > maxStrainIncValue) then
        numSteps = floor(abs(maxInc) / maxStrainIncValue) + 1 
        ! floor(A) returns the greatest integer less than or equal to A
        StrainInc = (NextStrain - CurStrain) / numSteps !real(, kind=8)

        ! Initialize temporary variables
        cStress = CurStress
        cStrain = CurStrain
        cAlpha = CurAlpha
        cFabric = CurFabric
        cAlpha_in = alpha_in
        cAlpha_in_p = alpha_in_p
        cEStrain = CurElasticStrain

        do ii = 1, numSteps
            nStrain = cStrain + StrainInc

            ! Call integration scheme subroutine --> pointer
            ! Temporary commented exp_int
            !call exp_int(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, cAlpha_in_p, &
            !             nStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
            !             NextL, NextVoidRatio, GG, KK, nCe, nCep, nCepC)

            ! Update temporary variables
            cStress = NextStress
            cStrain = nStrain
            cEStrain = NextElasticStrain
            cAlpha = NextAlpha
            cFabric = NextFabric
        end do
    else
        ! Call integration scheme subroutine without applying maximum strain increment constraint
        ! Temporary commented exp_int
        !call exp_int(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, &
        !             NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, &
        !             NextL, NextVoidRatio, GG, KK, aC, aCep, aCep_Consistent)
    end if


    end subroutine MaxStrainInc !PM4Sand_
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !---------------------------------------------------------------------------------------
    !// -------------------------------------------------------------------------------------------------------
    !/*************************************************************/
    !// Modified-Euler Integrator
    !/*************************************************************/
 !   subroutine ModifiedEuler(CurStress, StrainIncrement, CurAlpha, CurFabric, &
 !       alpha_in, alpha_in_p, &
 !       emax, emin, zpeak, zcum, MCur, Fabric_in, VoidRatio, Pmin, &   
 !       QQ_Bolton, RR_bolton, cz, &
 !       NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, &
 !       GG, KK, aC)!, &
 !       
 !       
 !       !subroutine ModifiedEuler(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, &
 !       !alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, &
 !       !GG, KK, aC)!, &
 !       !TolE)!, aCep, aCep_Consistent)
 !   
 !   implicit none
 !   
 !   ! Arguments
 !   real(REAL_TYPE), intent(in), dimension(3) :: CurStress
 !   real(REAL_TYPE), intent(in), dimension(3) :: StrainIncrement
 !   
 !   !CurStrain, 
 !   !CurElasticStrain, 
 !   real(REAL_TYPE), intent(in), dimension(3) :: CurAlpha
 !   real(REAL_TYPE), intent(in), dimension(3) :: CurFabric
 !   real(REAL_TYPE), intent(in), dimension(3) :: alpha_in
 !   real(REAL_TYPE), intent(in), dimension(3) :: alpha_in_p
 !   !real(REAL_TYPE), intent(in), dimension(3) :: NextStrain
 !   
 !   real(REAL_TYPE), intent(in) :: emax
 !   real(REAL_TYPE), intent(in) :: emin
 !   real(REAL_TYPE), intent(in) :: z_max
 !   
 !   real(REAL_TYPE), intent(in) :: Mcur
 !   real(REAL_TYPE), intent(in) :: zcum
 !   real(REAL_TYPE), intent(in) :: zpeak
 !   
 !   real(REAL_TYPE), intent(in) :: Pmin
 !   real(REAL_TYPE), intent(in) :: cz
 !   
 !   real(REAL_TYPE), intent(in) :: QQ_Bolton
 !   real(REAL_TYPE), intent(in) :: RR_Bolton
 !   
 !   real(REAL_TYPE), intent(in) :: K_p
 !   
 !   real(REAL_TYPE), intent(in) :: pzp
 !   
 !   real(REAL_TYPE), intent(in) :: mm
 !   
 !   real(REAL_TYPE), intent(in), dimension(3) :: Fabric_in
 !   
 !   
 !   
 !   real(REAL_TYPE), intent(out), dimension(3) :: NextElasticStrainIncrement
 !   real(REAL_TYPE), intent(out), dimension(3) :: NextStress
 !   real(REAL_TYPE), intent(out), dimension(3) :: NextAlpha
 !   real(REAL_TYPE), intent(out), dimension(3) :: NextFabric
 !   
 !   real(REAL_TYPE), intent(out) :: NextL 
 !   real(REAL_TYPE), intent(out) :: NextVoidRatio 
 !   real(REAL_TYPE), intent(out) :: GG
 !   real(REAL_TYPE), intent(out) :: KK
 !   
 !   real(REAL_TYPE), dimension(3, 3), intent(in) :: aC!, aCep, aCep_Consistent
 !    
 !   !real(REAL_TYPE), intent(in) :: TolE
 !   
 !   ! Local variables
 !   real(REAL_TYPE) :: NextDr
 !   real(REAL_TYPE) :: dVolStrain
 !   real(REAL_TYPE) :: pp
 !   real(REAL_TYPE) :: Cka
 !   real(REAL_TYPE) :: temp4
 !   real(REAL_TYPE) :: curStepError
 !   real(REAL_TYPE) :: qq_ModEuler 
 !   real(REAL_TYPE) :: stressNorm 
 !   real(REAL_TYPE) :: hh 
 !   real(REAL_TYPE) :: DD 
 !   real(REAL_TYPE) :: AlphaAlphaBDotN
 !   
 !   real(REAL_TYPE) :: TT
 !   real(REAL_TYPE) :: dT
 !   real(REAL_TYPE) :: dT_min
 !   
 !   real(REAL_TYPE) :: small != 1.0e-10
 !   real(REAL_TYPE), dimension(3) :: tmp0 
 !   real(REAL_TYPE), dimension(3) :: tmp1 
 !   real(REAL_TYPE), dimension(3) :: tmp2
 !   real(REAL_TYPE), dimension(3) :: nn 
 !   real(REAL_TYPE), dimension(3) :: R1 
 !   real(REAL_TYPE), dimension(3) :: R2 
 !   real(REAL_TYPE), dimension(3) :: alphaD
 !   real(REAL_TYPE), dimension(3) :: dDevStrain 
 !   real(REAL_TYPE), dimension(3) :: rrr
 !   real(REAL_TYPE), dimension(3) :: bb
 !   real(REAL_TYPE), dimension(3) :: dSigma1
 !   real(REAL_TYPE), dimension(3) :: dSigma2
 !   real(REAL_TYPE), dimension(3) :: dAlpha1
 !   real(REAL_TYPE), dimension(3) :: dAlpha2
 !   real(REAL_TYPE), dimension(3) :: dFabric1 
 !   real(REAL_TYPE), dimension(3) :: dFabric2 
 !   real(REAL_TYPE), dimension(3) :: dPStrain1
 !   real(REAL_TYPE), dimension(3) :: dPStrain2
 !   real(REAL_TYPE) :: DGamma
 !   
 !   !real(REAL_TYPE), intent(in) :: e_init
 !   
 !   real(REAL_TYPE) :: two3, mIIco ! These should be initialized with actual values
 !   
 !   ! Local variables for use with subroutines below 
 !   real(REAL_TYPE) :: GetTrace_result
	!real(REAL_TYPE) :: GetNorm_Contr_result
 !   
 !   real(REAL_TYPE) :: DoubleDot2_2_Contr_result
 !   real(REAL_TYPE) :: DoubleDot2_2_Mixed_result
 !   
 !   real(REAL_TYPE) :: Macauley_result
 !   real(REAL_TYPE) :: Macauley_result2
 !   real(REAL_TYPE) :: MacauleyIndex_result
 !   
 !   real(REAL_TYPE), dimension(3) :: GetDevPart_result
 !   real(REAL_TYPE), dimension(3) :: ToContraviant_result
 !   real(REAL_TYPE), dimension(3) :: ToCovariant_result
 !   
 !   real(REAL_TYPE), dimension(3) :: nStress
 !   real(REAL_TYPE), dimension(3) :: nFabric
 !   real(REAL_TYPE), dimension(3) :: nAlpha
 !   
 !   real(REAL_TYPE), dimension(3) :: I1
 !   
 !   real(REAL_TYPE), dimension(3) :: alphaD_NextAlpha
 !   
 !   I1(1) = 1
 !   I1(2) = 1
 !   I1(3) = 0
 !   
 !   TT = 0.0 ! start from zero
 !   dT = 1.0 ! start at an increment of one
 !   dT_min = 1e-4 
 !   TolE = 1e-5
 !   
 !   small = 1.0e-10
 !   
 !   two3 = 0.66666666667 !2.0/3.0
 !   !mIIco = 1.0
 !   
 !   ! initialize NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric
 !   NextElasticStrainIncrement = StrainIncrement 
 !   NextStress = CurStress
 !   NextAlpha = CurAlpha
 !   NextFabric = CurFabric
 !   
 !   ! obtain  KK and GG
 !   call GetElasticModuli_(NextStress, zcum, z_max, nu, G0, Md, Mb, &
 !       PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, Fsed_min, me2p)
 !
 !   ! Calculate p
 !   call GetTrace(CurStress, pp)
 !   pp = 0.5 * pp
 !
 !   ! Check if p < m_Pmin / 5.0 and adjust NextStress if needed
 !   if (pp < (Pmin / 5.0) ) then ! --> some sort of a tension cutoff
 !       !if (debugFlag) then
 !       !    write(*,*) "Tag = ", this%getTag(), " : p < pmin / 5, should not happen"
 !       !end if
 !
 !       call GetDevPart(NextStress, GetDevPart_result)
 !       NextStress = GetDevPart_result + ((Pmin / 5.0) * I1)
 !       
 !   end if
 !   
 !   ! Main loop
 !   do while (TT < 1.0) ! while less than one (step 7 in page 140 in Sloan et al., 2001)
 !       
 !       ! Calculate NextVoidRatio and NextDr
 !       tmp0 = NextStrainIncrement * TT ! portion of the strain 
 !       !tmp0 = CurStrain + tmp0
 !       
 !       call GetTrace(tmp0, GetTrace_result) ! --> volumetric strain
 !       NextVoidRatio = VoidRatio - ( (1.0 + VoidRatio) * GetTrace_result )
 !       
 !       NextDr = (emax - NextVoidRatio) / (emax - emin)
 !       
 !       ! Calculate dVolStrain and dDevStrain
 !       call GetTrace(StrainIncrement, dVolStrain)
 !       dVolStrain = dT * dVolStrain
 !       
 !       dDevStrain = ((StrainIncrement) * dT) - (I1 * (-dVolStrain / 3.0))
 !       
 !       call GetTrace(NextStress, pp)
 !       pp = 0.5 * pp
 !       
 !       
 !       ! Calculate Delta 1
 !       ! Implement the calculation of dSigma1, dAlpha1, dFabric1, and dPStrain1
 !       call GetStateDependent(NextStress, NextAlpha, alpha_in, alpha_in_p, alpha_in_true, NextFabric, Fabric_in, &
 !               GG, zcum, zpeak, pzp, Mcur, NextDr, &
 !               Mc, nd, nb, Pmin, Pmin2, P_atm, & 
 !               mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
 !               QQ_Bolton, RR_Bolton, &
 !               nn, alphaB, alphaD, bb, R1, Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, &
 !               ksi)
 !       
 !   
 !       call GetDevPart(NextStress, rrr)
 !       rrr = rrr / pp
 !
 !       call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
 !       temp4 = K_p + (2 * GG) - (KK * DD * DoubleDot2_2_Contr_result) 
 !       
 !       if (abs(temp4) < small) then ! <-- to prevent dividing by zero
 !   
 !           ! neutral loading
 !           dSigma1 = 0.0
 !           dAlpha1 = 0.0
 !           dFabric1 = 0.0
 !           dPStrain1 = tmp0
 !
 !       else
 !           
 !           call DoubleDot2_2_Mixed(nn, dDevStrain, DoubleDot2_2_Mixed_result)
 !           call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
 !
 !           !NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4
 !           NextL = ( (2 * GG * DoubleDot2_2_Mixed_result) - (DoubleDot2_2_Contr_result * KK * dVolStrain) ) / temp4
 !           
 !           if (NextL < 0) then
 !               !if (debugFlag) then
 !               !    write(*,*) "1 NextL is smaller than 0"
 !               !    write(*,*) "NextL = ", NextL
 !               !end if
 !               call ToContraviant(dDevStrain, ToContraviant_result)
 !       
 !               dSigma1 = (2 * GG * ToContraviant_result) + (KK * dVolStrain * I1) ! if L<0 then purely elastic
 !               dAlpha1 = 0.0
 !               dFabric1 = 0.0
 !               dPStrain1 = 0.0
 !   
 !           else
 !               
 !               ! plastic stress portion
 !               tmp0 = nn * (2.0 * GG)
 !               tmp1 = I1 * (KK * DD)
 !               tmp1 = tmp1 + tmp0
 !               
 !               call Macauley(NextL, Macauley_result)
 !               tmp1 = (-Macauley_result) * tmp1 ! plastic stress portion
 !               tmp2 = I1 * (KK * dVolStrain) ! volumetric elastic stress
 !               
 !               call ToContraviant(dDevStrain, ToContraviant_result)
 !               tmp3 = ToContraviant_result * (2.0 * GG) ! shear elastic stress
 !               dSigma1 = tmp3 + tmp2 + tmp1
 !                        ! dev    vol   plastic
 !
 !               ! update fabric				
 !               alphaD_NextAlpha = alphaD - NextAlpha
 !               call DoubleDot2_2_Contr(alphaD_NextAlpha, nn, DoubleDot2_2_Contr_result)
 !               
 !               if (DoubleDot2_2_Contr_result < 0.0) then !if (DoubleDot2_2_Contr(alphaD - NextAlpha, n) < 0.0)
 !               
 !                   ! update fabric terms for dilation 
 !                   
 !                   ! Equation 57
 !                   dFabric1 = (nn*z_max) + NextFabric
 !                               
 !                   !dFabric1 = dFabric1 * (-m_cz / (1 + Macauley(mzcum / (2.0 * m_z_max) - 1.0)) * Macauley(NextL) * MacauleyIndex(-D))
 !                   call Macauley((zcum / (2.0 * z_max)) - 1.0, Macauley_result)
 !                   call Macauley(NextL, Macauley_result2)
 !                   call MacauleyIndex(-DD, MacauleyIndex_result)
 !                   
 !                   dFabric1 = dFabric1 * ( (-cz / (1 + Macauley_result)) * Macauley_result2 * MacauleyIndex_result * (1/DD) )
 !                   
 !               end if
 !               
 !               call ToCovariant(R1, dPStrain1)
 !               dPStrain1 = dPStrain1 * NextL
 !               dAlpha1 = bb * (two3 * NextL * hh)
 !               
 !           end if
 !           
 !       end if
 !
 !       ! Calculate Delta 2
 !       ! Implement the calculation of dSigma2, dAlpha2, dFabric2, and dPStrain2
 !       ! tmp0 = NextStress; tmp0 += dSigma1;   ! tmp0 is NextStress + dSigma1 until line 1511
 !       tmp0 = NextStress + dSigma1
 !       
 !       call GetTrace(tmp0, pp)
 !       pp = 0.5 * pp
 !       
 !       if (pp < 0) then
 !           if (dT == dT_min) then
 !       
 !               !if (debugFlag) then
 !               !
 !               !    write(*,*) "Delta 1: p < 0"
 !               !
 !               !endif
 !               
 !               NextElasticStrainIncrement = StrainIncrement !CurElasticStrain + (NextStrain - CurStrain)
 !               NextStress = CurStress
 !               NextAlpha = CurAlpha
 !               NextFabric = CurFabric
 !       
 !               return
 !   
 !           endif
 !   
 !           dT = max(0.1 * dT, dT_min)
 !   
 !           cycle
 !
 !       end if
 !
 !       ! tmp1.Zero();  tmp1 += NextAlpha; tmp1 += dAlpha1;  ! tmp1 is NextAlpha + dAlpha1
 !       tmp1 = NextAlpha + dAlpha1
 !       
 !       ! tmp2.Zero();  tmp2 += NextFabric; tmp2 += dFabric1;  ! tmp2 is NextFabric + dFabric1
 !       tmp2 = NextFabric + dFabric1
 !       
 !       call GetStateDependent(tmp0, tmp1, alpha_in, alpha_in_p, alpha_in_true, tmp2, fabric_in, &
 !               GG, zcum, zpeak, pzp, Mcur, NextDr, &
 !               Mc, nd, nb, Pmin, Pmin2, P_atm, & 
 !               mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
 !               QQ_Bolton, RR_Bolton, &
 !               nn, alphaB, alphaD, bb, R2, &
 !               Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, &
 !               ksi)
 !       
 !       ! r = GetDevPart(NextStress + dSigma1) / p;
 !       call GetDevPart(tmp0, rrr)
 !       rrr = rrr / pp
 !       call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
 !       temp4 = K_p + (2 * GG) - (KK * DD * DoubleDot2_2_Contr_result)!DoubleDot2_2_Contr(nn, rr)
 !       
 !       if (abs(temp4) < small) then !--> prevent division by zero!
 !       
 !           ! neutral loading
 !           dSigma2 = 0
 !           dAlpha2 = 0
 !           dFabric2 = 0
 !           
 !           ! dPStrain2 = dDevStrain + dVolStrain * mI1;
 !           dPStrain2 = dPStrain1
 !
 !       else
 !           
 !           call DoubleDot2_2_Mixed(nn, dDevStrain, DoubleDot2_2_Mixed_result)
 !           call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
 !
 !           !NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4
 !           NextL = ((2 * GG * DoubleDot2_2_Mixed_result) - (DoubleDot2_2_Contr_result * KK * dVolStrain)) / temp4
 !           DGamma = NextL
 !           
 !           if (NextL < 0) then
 !               !if (debugFlag) then
 !               !    write(*,*) "2 NextL is smaller than 0"
 !               !    write(*,*) "NextL = ", NextL
 !               !endif
 !       
 !               call ToContraviant(dDevStrain, ToContraviant_result)
 !               
 !               dSigma2 = (2 * GG * ToContraviant_result) + (KK * dVolStrain * I1)
 !               dAlpha2 = 0
 !               dFabric2 = 0
 !               dPStrain2 = 0
 !               ! dSigma2.Zero();
 !               ! dPStrain2 = dPStrain1;
 !   
 !           
 !           else
 !       
 !               ! dSigma2 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextL) *
 !               !  (2.0 * G * n + K * D * mI1);
 !               tmp0 = nn * (2.0 * GG) ! shear part 
 !               tmp1 = I1 * (KK * DD) ! volumetric part
 !               tmp1 = tmp1 + tmp0
 !               
 !               call Macauley(NextL, Macauley_result)
 !               
 !               tmp1 = tmp1 * (-Macauley_result) ! plastic stress
 !               tmp2 = I1 * (KK * dVolStrain) ! elastic volumetric stress
 !               
 !               call ToContraviant(dDevStrain, ToContraviant_result)
 !               
 !               tmp3 = ToContraviant_result * (2.0 * GG) ! elastic shear stress
 !               dSigma2 = tmp3 + tmp2 + tmp1
 !               
 !               ! update fabric
 !               alphaD_NextAlpha = alphaD - NextAlpha - dAlpha1
 !               call DoubleDot2_2_Contr(alphaD_NextAlpha, nn, DoubleDot2_2_Contr_result)
 !               
 !               if (DoubleDot2_2_Contr_result < 0.0) then ! if (DoubleDot2_2_Contr(alphaD - (NextAlpha + dAlpha1), n) < 0.0) then
 !                   
 !                   ! Equation 57
 !                   ! dFabric2 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + NextFabric + dFabric1);
 !                   dFabric2 = nn * z_max
 !                   dFabric2 = dFabric2 + NextFabric + dFabric1
 !                   
 !                   Macauley_result = 0.0
 !                   Macauley_result2 = 0.0
 !                   MacauleyIndex_result = 0.0
 !                   
 !                   call Macauley( (zcum / (2.0*z_max)) - 1.0, Macauley_result)
 !                   call Macauley(NextL, Macauley_result2)
 !                   call MacauleyIndex(-DD, MacauleyIndex_result)
 !
 !                   ! I included (1/DD) because that's what is shown in the literature 
 !                   dFabric2 = dFabric2 * -1.0 * ( (cz / (1 + Macauley_result)) * Macauley_result2 * MacauleyIndex_result * (1/DD))
 !
 !               end if
 !               
 !               ! dPStrain2 = NextL * mIIco * R2
 !               call ToCovariant(R2,dPStrain2)
 !               dPStrain2 = dPStrain2 * NextL
 !               
 !               ! dAlpha2 = two3 * NextL * h * b
 !               dAlpha2 = bb * (two3 * NextL * hh)
 !           
 !           
 !           endif
 !
 !       endif
 !
 !       ! Note:
 !       ! -dSigma1 is calculated based on an initial increment of strain and stress
 !       ! -dSigma2 is calcualted based on an updated stress state (NextStress+dSigma1), 
 !       !  including fabric updates 
 !       
 !       ! Update nStress, nFabric, and nAlpha
 !       nStress = (0.5 * (dSigma1 + dSigma2)) + NextStress
 !       nFabric = (0.5 * (dFabric1 + dFabric2)) + NextFabric
 !       nAlpha  = (0.5 * (dAlpha1 + dAlpha2)) + NextAlpha
 !       
 !       ! Check p
 !       call GetTrace(nStress, pp)
 !       pp = 0.5 * pp !GetTrace(nStress)
 !
 !       ! Perform stress correction if needed
 !       if (pp < 0) then
 !           
 !           if (dT == dT_min) then
 !               
 !               ! Implement the stress correction
 !               NextElasticStrainIncrement = StrainIncrement
 !               NextStress = CurStress
 !               NextAlpha = CurAlpha
 !               NextFabric = CurFabric
 !               
 !               return
 !               
 !           end if
 !               
 !           dT = max(0.1 * dT, dT_min)
 !   
 !           cycle
 !           
 !       end if
 !
 !       ! Calculate current step error
 !       call GetNorm_Contr(NextStress, stressNorm)
 !       
 !       tmp0 = dSigma2 - dSigma1
 !       
 !       if (stressNorm < 0.5) then
 !           
 !           call GetNorm_Contr(tmp0, curStepError)
 !           
 !       else
 !           
 !           call GetNorm_Contr(tmp0, curStepError)
 !           curStepError = curStepError / (2 * stressNorm)
 !
 !       end if
 !
 !       ! Check if the current step error is within tolerance
 !       if (curStepError > TolE) then
 !           
 !           ! Adjust step size
 !           qq_ModEuler = max(0.8 * sqrt(TolE / curStepError), 0.1)
 !           
 !           if (dT == dT_min) then
 !               ! Implement the step size adjustment and stress correction
 !               NextElasticStrainIncrement = NextElasticStrainIncrement - (0.5 * (dPStrain1 + dPStrain2))
 !               NextStress = nStress
 !               NextAlpha = nAlpha
 !               
 !               call Stress_Correction(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio)
 !               !//Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, 0.5 * (R1 + R2), n, r);
 !               
 !               TT = TT + dT
 !               
 !           end if
 !
 !           dT = max(qq_ModEuler * dT, dT_min)
 !       
 !       else
 !           
 !           ! Implement the step size adjustment and update NextElasticStrain, NextStress, NextAlpha, and NextFabric
 !           NextElasticStrainIncrement = NextElasticStrainIncrement - (0.5 * (dPStrain1 + dPStrain2))
 !           NextStress = nStress
 !           NextAlpha = nAlpha
 !           NextFabric = nFabric
 !           
 !           !call Stress_Correction_(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio)
 !           !//Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, 0.5 * (R1 + R2), n, r);
 !           
 !           call stress_correction_(alpha_in, alpha_in_p, alpha_in_true, Fabric_in, CurFabric, CurStress, Pmin, P_atm, &
 !                            TolF, Mc, emax, emin, zcum, zpeak, z_max, pzp, mm, h0, hpo, QQ_Bolton, RR_Bolton, &
 !                            DGamma, alphaD, alphaB, dSigmaP, DD, Mb, Md, K_p, NextVoidRatio, Mcur, GG, KK, NextStress, &
 !                            NextAlpha, dSigma)
 !           
 !           TT = TT + dT
 !           qq_ModEuler = max(0.8 * sqrt(TolE / curStepError), 0.5)
 !           dT = max(qq_ModEuler * dT, dT_min)
 !           dT = min(dT, 1.0 - TT)
 !       
 !       end if
 !   
 !   end do
 !   
 !
 !       
 !       
 !       end subroutine ModifiedEuler
    
    
    
    
    
    
    
    !--------------------------------------------------------------
    !// -------------------------------------------------------------------------------------------------------
    !/*************************************************************/
    !// Runge-Kutta Integrator
    !/*************************************************************/
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !--------------------------------------------------------------
    !/*************************************************************/
    !//            Pegasus Iterations                             //
    !/*************************************************************/
    !subroutine IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, a0, a1, aa)
    subroutine IntersectionFactor(CurStress, StrainIncrement, CurAlpha, Ce, a0, a1, ElasticRatio, TolF)
    
    implicit none
    
    ! inputs
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress
    real(REAL_TYPE), dimension(3), intent(in) :: CurAlpha
    real(REAL_TYPE), dimension(3), intent(in) :: StrainIncrement
    
    real(REAL_TYPE), dimension(3,3), intent(in) :: Ce
    
    real(REAL_TYPE), intent(in) :: TolF
    
    ! outputs
    real(REAL_TYPE), intent(out) :: ElasticRatio
    
    ! inout
    real(REAL_TYPE), intent(inout) :: a0 
    real(REAL_TYPE), intent(inout) :: a1
    
    !local variables
    real(REAL_TYPE) :: ff 
    real(REAL_TYPE) :: f0 
    real(REAL_TYPE) :: f1
    
    real(REAL_TYPE), dimension(3) :: NextElasticStress
    
    real(REAL_TYPE), dimension(3) :: dSigma 
    real(REAL_TYPE), dimension(3) :: dSigma0 
    real(REAL_TYPE), dimension(3) :: dSigma1
    
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
    
    integer(INTEGER_TYPE) :: ii
    real(REAL_TYPE) :: GetFYieldFunction_result
    
    ! Constants
    integer(INTEGER_TYPE) :: maxIter! = 10
    real(REAL_TYPE) :: small! = 1.0e-10
    
    ! establish numerical values
    maxIter = 10
    small = 1.0e-10
    
    ! Calculate dSigma0 and f0
    dSigma0 = matmul(Ce, (a0 * StrainIncrement)) ! a0 stress increment
    
    ! temporary new stress for a0
    NextElasticStress = CurStress + dSigma0
    
    ! Correct the stress
    ! I added this so that we can have the correct considertion of the two shear components
    call ToCovariant(NextElasticStress, ToCovariant_result)
    
    ! evaluate yield function for this a0
    call GetFYieldFunction(ToCovariant_result, CurAlpha, f0) ! a0 yield function evaluation
    
    ! calculate dSigma1 and f1
    dSigma1 =  matmul(Ce, (a1 * StrainIncrement)) ! a1 stress increment

    ! temporary new stress for a1
    NextElasticStress = CurStress + dSigma1 
    
    ! Correct the stress
    ! I added this so that we can have the correct considertion of the two shear components
    call ToCovariant(NextElasticStress, ToCovariant_result)
    
    ! evaluate yield function for this a1
    call GetFYieldFunction(ToCovariant_result, CurAlpha, f1) ! a1 yield function evaluation

    ! Iterate using the Pegasus method
    do ii = 1, maxIter
        
        ! Calculate new a
        ElasticRatio = a1 - (f1 * (a1 - a0) / (f1 - f0) ) !--> interpolation
        !dSigma = mCe * (aa * strainInc)
        dSigma = matmul(Ce , (ElasticRatio * StrainIncrement))
        
        NextElasticStress = CurStress + dSigma !--> new stress
        
        ! Correct the stress
        ! I added this so that we can have the correct considertion of the two shear components
        ToCovariant_result = 0.0
        call ToCovariant(NextElasticStress, ToCovariant_result)
        
        ! evaluate yield function
        call GetFYieldFunction(ToCovariant_result, CurAlpha, ff)

        ! Check convergence
        if (abs(ff) < TolF) then
            exit
        end if

        ! Update bounds for the next iteration
        if ( (ff * f0) < 0) then
            ! upper bound
            a1 = ElasticRatio
            f1 = ff
        else
            ! lower bound and gradient
            f1 = f1 * f0 / (f0 + ff)
            a0 = ElasticRatio
            f0 = ff
        end if

        ! Check for non-convergence
        if (ii == maxIter) then
            !if (debugFlag) then
            !    write(*, *) "Didn't find alpha!"
            !endif
            !aa = 0.0
            ! SOME THING WRONG HAPPENED IF WE ARE HERE
            exit
        end if
        
    end do

    ! Ensure a is within bounds
    if (ElasticRatio > (1.0 - small)) ElasticRatio = 1.0
    if (ElasticRatio < small) ElasticRatio = 0.0
    

    end subroutine IntersectionFactor
    
    !----------------------------------------------------------------
    !/*************************************************************/
    !//      Pegasus Iterations  (ElastoPlastic Unloading)        //
    !/*************************************************************/
    
    subroutine IntersectionFactor_Unloading(CurStress, StrainIncrement, CurAlpha, Ce, TolF, ElasticRatio) !CurStrain, NextStrain,
    
    implicit none
    
    ! inputs
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress
    real(REAL_TYPE), dimension(3), intent(in) :: StrainIncrement
    real(REAL_TYPE), dimension(3), intent(in) :: CurAlpha
    real(REAL_TYPE), dimension(3,3), intent(in) :: Ce
    real(REAL_TYPE), intent(in) :: TolF ! = 1.0e-8
    
    ! outputs
    real(REAL_TYPE), intent(out) :: ElasticRatio
    
    ! Local variables
    real(REAL_TYPE) :: a0
    real(REAL_TYPE) :: a1 
    real(REAL_TYPE) :: da
    
    real(REAL_TYPE) :: ff
    real(REAL_TYPE) :: f0
    real(REAL_TYPE) :: f1
    real(REAL_TYPE) :: fs
    
    real(REAL_TYPE), dimension(3) :: dSigma
    real(REAL_TYPE), dimension(3) :: dSigma0 
    real(REAL_TYPE), dimension(3) :: dSigma1 
    real(REAL_TYPE), dimension(3) :: NextElasticStressPoint !tmp
    logical :: flag = .false.

    integer(INTEGER_TYPE) :: ii
    integer(INTEGER_TYPE) :: kk
    integer(INTEGER_TYPE) :: nSub! = 20
    integer(INTEGER_TYPE) :: maxIter! = 10
    
    real(REAL_TYPE) :: GetFYieldFunction_result
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
    
    ! initialize variables
    nSub = 20
    maxIter = 10
    
    ElasticRatio = 0.0
    a0 = 0.0 
    a1 = 1.0
    
    ! Correct the stress
    ! I added this so that we can have the correct considertion of the two shear components
    ToCovariant_result = 0.0
    call ToCovariant(CurStress, ToCovariant_result)
    
    ! Compute initial value of f0
    call GetFYieldFunction(ToCovariant_result, CurAlpha, f0)
    fs = f0
    
    ! Compute initial value of dSigma
    call DoubleDot4_2(Ce, strainIncrement, dSigma)

    ! Pegasus iterations loop
    do ii = 1, maxIter
        
        ! Compute increment factor
        da = (a1 - a0) / real(nSub)
        
        
        do kk = 1, nSub ! Iterate over subintervals
            
            ! evaluate elastic ratio
            ElasticRatio = a0 + (da * real(kk))
            
            ! Compute stress increment
            NextElasticStressPoint = (ElasticRatio * dSigma) + CurStress
            
            ! Correct the stress
            ! I added this so that we can have the correct considertion of the two shear components
            ToCovariant_result = 0.0
            call ToCovariant(NextElasticStressPoint, ToCovariant_result)
            
            ! Compute value of f
            call GetFYieldFunction(ToCovariant_result, CurAlpha, ff)

            ! Check if f exceeds tolerance
            if (ff > TolF) then
                
                a1 = ElasticRatio
                
                if (f0 < -TolF) then
                    f1 = ff
                    flag = .true.
                    exit
                else
                    a0 = 0.0
                    f0 = fs
                    exit
                end if
            
            else
            
                a0 = ElasticRatio
                f0 = ff
            
            end if
            
        end do ! Iterate over subintervals

        ! Check termination conditions
        if (flag) exit
        
        if (ii == maxIter) then
            
            !if (debugFlag) then
            !    write(*, *) "Didn't find alpha! - Unloading, a0 = ", a0, ", a1 = ", a1
            !endif
            ElasticRatio = 0.0
            return
            
        end if
        
    end do

    ! Output the result
    !if (debugFlag) then
    !    write(*, *) "Found alpha - Unloading, a0 = ", a0, ", a1 = ", a1
    !endif
    ! we update a0 and a1, then plug it in IntersectionFactor to get aa
    call IntersectionFactor(CurStress, StrainIncrement, CurAlpha, Ce, a0, a1, ElasticRatio, TolF) !CurStrain, NextStrain,
    

    end subroutine IntersectionFactor_Unloading
    
    
    
    
    
    
    
    
    !--------------------------------------------------------
    !/*************************************************************/
    !//            Stress Correction                              //
    !/*************************************************************/
    subroutine stress_correction_(alpha_in, alpha_in_p, alpha_in_true, Fabric_in, CurFabric, CurStress, Pmin, P_atm, &
                             TolF, Mc, emax, emin, zcum, zpeak, z_max, pzp, mm, h0, hpo, QQ_Bolton, RR_Bolton, &
                             DGamma, alphaD, alphaB, dSigmaP, DD, Mb, Md, K_p, NextVoidRatio, Mcur, GG, KK, NextStress, &
                             NextAlpha, dSigma, &
                             nb, nd, &
                             Pmin2, Ado)
    
    implicit none
    
    ! in variables 
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in !y 
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in_p !y
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in_true !y
    real(REAL_TYPE), dimension(3), intent(in) :: Fabric_in !y 
    real(REAL_TYPE), dimension(3), intent(in) :: CurFabric !y
    real(REAL_TYPE), dimension(3), intent(in) :: CurStress !y
    real(REAL_TYPE), intent(in) :: Pmin !y 
    real(REAL_TYPE), intent(in) :: P_atm !y
    
    real(REAL_TYPE), intent(in) :: TolF !y
    
    real(REAL_TYPE), intent(in) :: Mc  !y
    real(REAL_TYPE), intent(in) :: emax !y
    real(REAL_TYPE), intent(in) :: emin !y
    
    real(REAL_TYPE), intent(in) :: zcum !y 
    real(REAL_TYPE), intent(in) :: zpeak !y
    real(REAL_TYPE), intent(in) :: z_max !y
    
    real(REAL_TYPE), intent(in) :: pzp !y
    
    real(REAL_TYPE), intent(in) :: mm !y
    
    real(REAL_TYPE), intent(in) :: h0 !y
    
    real(REAL_TYPE), intent(in) :: hpo !y
    
    real(REAL_TYPE), intent(in) :: QQ_Bolton !y
    real(REAL_TYPE), intent(in) :: RR_Bolton !y
    
    real(REAL_TYPE), intent(in) :: DGamma ! this is practically the L value to calculate the plastic stresses !
    
    
    real(REAL_TYPE), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
    real(REAL_TYPE), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
    
    real(REAL_TYPE), intent(in) :: Pmin2      ! statev 26

    real(REAL_TYPE), intent(in) :: Ado
    
    ! output variable
    real(REAL_TYPE), dimension(3), intent(out) :: alphaD !y
    real(REAL_TYPE), dimension(3), intent(out) :: alphaB !y

    real(REAL_TYPE), dimension(3), intent(out) :: dSigmaP !y
    
    real(REAL_TYPE), intent(out) :: DD !y
    real(REAL_TYPE), intent(out) :: Mb !y 
    real(REAL_TYPE), intent(out) :: Md !y
    real(REAL_TYPE), intent(out) :: K_p !y
    
    
    ! inout variables
    real(REAL_TYPE), intent(inout) :: NextVoidRatio
    
    real(REAL_TYPE), intent(inout) :: Mcur
    
    real(REAL_TYPE), intent(inout) :: GG
    real(REAL_TYPE), intent(inout) :: KK !--> probably inout
    
    real(REAL_TYPE), dimension(3), intent(inout) :: NextStress
    real(REAL_TYPE), dimension(3), intent(inout) :: NextAlpha
    
    real(REAL_TYPE), dimension(3), intent(inout) :: dSigma
    
    
    ! Local variables
    real(REAL_TYPE), dimension(3) :: dfrOverdSigma
    real(REAL_TYPE), dimension(3) :: dfrOverdAlpha
    real(REAL_TYPE), dimension(3) :: nn !y
    real(REAL_TYPE), dimension(3) :: RR !y
    
    real(REAL_TYPE), dimension(3) :: bb !y
    real(REAL_TYPE), dimension(3) :: aBar
    real(REAL_TYPE), dimension(3) :: rrr
    
    real(REAL_TYPE), dimension(3) :: tmp0
    real(REAL_TYPE), dimension(3) :: tmp1
    
    real(REAL_TYPE), dimension(3) :: nAlpha !y
    real(REAL_TYPE), dimension(3) :: nStress !y
    
    
    real(REAL_TYPE) :: ksi
    
    real(REAL_TYPE) :: CurDr !y
    
    real(REAL_TYPE) :: Ckaf !y
    real(REAL_TYPE) :: Cdr !y
    real(REAL_TYPE) :: Ceps !y
    
    real(REAL_TYPE) :: lambda
    real(REAL_TYPE) :: Cka 
    real(REAL_TYPE) :: hh !y
    real(REAL_TYPE) :: pp
    real(REAL_TYPE) :: fr
    real(REAL_TYPE) :: AlphaAlphaBDotN !y
    
    
    real(REAL_TYPE) :: maxIter
    
    real(REAL_TYPE) :: fr_old
    real(REAL_TYPE) :: alpha_up
    real(REAL_TYPE) :: alpha_mid
    real(REAL_TYPE) :: alpha_down
    
    integer(INTEGER_TYPE) :: jj
    integer(INTEGER_TYPE) :: ii 
    
    real(REAL_TYPE) :: two3
    
    real(REAL_TYPE) :: GetTrace_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: DoubleDot2_2_Mixed_result
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result2
    
    real(REAL_TYPE) :: Macauley_result
    real(REAL_TYPE) :: Macauley_result2
    real(REAL_TYPE) :: MacauleyIndex_result
    
	real(REAL_TYPE) :: GetNorm_Contr_result
    real(REAL_TYPE) :: GetFYieldFunction_result
    
    real(REAL_TYPE), dimension(3) :: I1
    
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    real(REAL_TYPE), dimension(3) :: ToContraviant_result
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
	
    !real(REAL_TYPE), dimension(3) :: ToCovariant_result
    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
    
    real(REAL_TYPE), dimension(3,3) :: GetStiffness_result, Ce
    
    ! Constants
    I1(1) = 1
    I1(2) = 1
    I1(3) = 0
    
    maxIter = 25
    two3 = 0.6666666666667

    
    ! Compute p
    call GetTrace(NextStress, pp)
    pp = 0.5 * pp!GetTrace(NextStress, GetTrace_result)

    ! Check tension cutoff
    if (pp < (Pmin / 5.0) ) then
        
        ! I added this so that we can have the correct considertion of the two shear components
        call ToCovariant(NextStress, ToCovariant_result)
        
        call GetFYieldFunction(ToCovariant_result, NextAlpha, fr)
        
        if (fr < TolF) then
            ! Stress state inside yield surface
            NextStress = NextStress + ( ((Pmin / 5.0) - pp) * I1)
        else
            ! Stress state outside yield surface
            NextStress = (Pmin / 5.0) * I1
            NextStress(3) = 0.8 * Mc * Pmin / 5.0
            NextAlpha = 0.0
            NextAlpha(3) = 0.8 * Mc
            return
        endif
        
    else
        ! Correct the stress
        ! I added this so that we can have the correct considertion of the two shear components
        call ToCovariant(NextStress, ToCovariant_result)
        
        ! check the initial state of stress: is it in the yield surface?
        call GetFYieldFunction(ToCovariant_result, NextAlpha, fr) !this works because the strain 
        
        ! if you are inside the yield surface then RETURN: You are done!
        if (fr < TolF) RETURN
        
        ! if you are outside the yield surface: then you need to calculate residuals:
        ! these represent the difference between the actual stress state and the stress state predicted by the plasticity model
        ! Compute CurDr
        ! use next stress and next alpha
        CurDr = (emax - NextVoidRatio) / (emax - emin)
        nStress = NextStress
        nAlpha = NextAlpha
        
        ! Stress correction loop
        do ii = 1, maxIter
            
            ! Compute r
            call GetDevPart(nStress, rrr)
            rrr = rrr / pp  ! calculate the stress ratio for the next stress state 
            
            ! get the state dependent variables which are a function of nStress and nAlpha
            ! Call GetStateDependent --> need to be checked
            call GetStateDependent(nStress, nAlpha, alpha_in, alpha_in_p, alpha_in_true, CurFabric, Fabric_in, &
                GG, zcum, zpeak, pzp, Mcur, CurDr, &
                Mc, nd, nb, Pmin, Pmin2, P_atm, & 
                mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
                QQ_Bolton, RR_Bolton, &
                nn, alphaB, alphaD, bb, RR, &
                Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, &
                ksi)           
            
            ! Compute the increment in stress caused by plastic deformation: It is calculated as the product of the 
            ! material stiffness matrix (aC) get the 3x3 2D stiffness matrix 
            call GetStiffness(KK, GG, Ce) 

            ! Compute dSigmaP
            ! ToCovariant multiplies the shear strain to get the engineering shear strain 
            ! to be used to find the stress using elasticity theory
            call ToCovariant(RR, ToCovariant_result) !--> RR is the direction of the plastic strain
                                                     !--> RR' is the deviatoric component of RR
            
            ! dSigmaP = aC * PLASTIC STRAIN = aC : [<L>R] !  WE ARE CONSIDERING ENTIRE R
            call DoubleDot4_2(Ce, DGamma * ToCovariant_result, dSigmaP) 
            !                   3x3 matrix          scalar * RR
            
            ! Compute aBar --> dAlpha
            aBar = two3 * hh * bb 
            
            call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
            
            ! Compute dfrOverdSigma
            dfrOverdSigma = (-0.5 * DoubleDot2_2_Contr_result * I1) + nn
            ! Compute dfrOverdAlpha
            dfrOverdAlpha = -pp * nn
            
            ! Compute lambda
            call DoubleDot2_2_Contr(dfrOverdSigma, dSigmaP, DoubleDot2_2_Contr_result)
            call DoubleDot2_2_Contr(dfrOverdAlpha, aBar, DoubleDot2_2_Contr_result2)
            lambda = fr / (DoubleDot2_2_Contr_result - DoubleDot2_2_Contr_result2)
            
            ! Update NextStress and NextAlpha
            ! This is consistency condition and flow rule
            tmp0 = - (lambda * dSigmaP) + nStress ! next stress
            tmp1 = (lambda * aBar) + nAlpha ! next alpha
            
            ! Correct the stress
            ! I added this so that we can have the correct considertion of the two shear components
            call ToCovariant(tmp0, ToCovariant_result)
            
            ! evaluate the yield function
            call GetFYieldFunction(ToCovariant_result, tmp1, GetFYieldFunction_result)
            
            
            if (abs(GetFYieldFunction_result) < abs(fr)) then ! first if statement in the stress correction loop
                nStress = tmp0
                nAlpha = tmp1
            else
                call DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma, DoubleDot2_2_Contr_result)
                lambda = fr / DoubleDot2_2_Contr_result !DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma) 
                nStress = nStress - (lambda * dfrOverdSigma)
            end if
            
            ! Correct the stress
            ! I added this so that we can have the correct considertion of the two shear components
            ToCovariant_result = 0.0
            call ToCovariant(nStress, ToCovariant_result)
            
            ! Update fr
            call GetFYieldFunction(nStress, nAlpha, fr)
            
             if (abs(fr) < TolF) then ! SECOND IF in the stress correction loop --> go to step 8
                 ! corrected state
                NextStress = nStress ! update NextStress
                NextAlpha = nAlpha ! update NextAlpha
                return
             end if
             
            ! Update p
            call GetTrace(nStress, pp)
            pp = max(0.5 * pp, Pmin)
            
        end do
        
        ! this is an additional layer of correction
        ! Search for the optimum alpha_mid
        dSigma = NextStress - CurStress
        alpha_up = 1.0
        alpha_mid = 0.5
        alpha_down = 0.0
        tmp0 = (alpha_mid * dSigma) + CurStress !mSigma
        
        ! Correct the stress
        ! I added this so that we can have the correct considertion of the two shear components
        ToCovariant_result = 0.0
        call ToCovariant(tmp0, ToCovariant_result)
        
        call GetFYieldFunction(ToCovariant_result, NextAlpha, fr_old)
        
        do jj = 1, maxIter
            
            if (fr_old < 0.0) then
                alpha_down = alpha_mid
                alpha_mid = 0.5 * (alpha_up + alpha_mid)
            else
                alpha_up = alpha_mid
                alpha_mid = 0.5 * (alpha_down + alpha_mid)
            endif
            
            tmp0 = (alpha_mid * dSigma) + CurStress!mSigma
            
            ! Correct the stress
            ! I added this so that we can have the correct considertion of the two shear components
            call ToCovariant(tmp0, ToCovariant_result)
            
            !GetFYieldFunction_result = 0.0
            call GetFYieldFunction(tmp0, NextAlpha, fr_old)
                        
            if (abs(fr_old) < TolF) then
                NextStress = tmp0
                exit
            endif
            
        end do
        
    endif
    

end subroutine Stress_Correction_
    
    
    
    
    
    
    
    
    
    !---------------------------------------------------------------
    
    subroutine Stress_Correction(NextStress, NextAlpha, dAlpha, rrr, &    
        RR, nn, mm, TolF)
    
    implicit none
    
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(inout) :: NextStress 
    real(REAL_TYPE), dimension(3), intent(inout) :: NextAlpha 
    real(REAL_TYPE), dimension(3), intent(inout) :: dAlpha 
    real(REAL_TYPE), dimension(3), intent(inout) :: rrr
    
    real(REAL_TYPE), dimension(3), intent(in) :: RR 
    real(REAL_TYPE), dimension(3), intent(in) :: nn
    
    real(REAL_TYPE), intent(in) :: mm
    
    real(REAL_TYPE), intent(in) :: TolF
    
    ! Local variables
    real(REAL_TYPE), dimension(3) :: dfrOverdSigma
    real(REAL_TYPE) :: lambda
    integer(INTEGER_TYPE) :: maxIter, ii
    real(REAL_TYPE) :: ff
    
    ! Local variables 
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    real(REAL_TYPE) :: GetFYieldFunction_result
    
    real(REAL_TYPE), dimension(3) :: ToCovariant_result
    
    real(REAL_TYPE), dimension(3) :: I1
    real(REAL_TYPE) :: root12

    ! Constants
    maxIter = 50
    !mTolF = 1.0e-5
    
    I1(1) = 1
    I1(2) = 1
    I1(3) = 0
    
    root12 = 0.70710678118
    ! Correct the stress
    ! I added this so that we can have the correct considertion of the two shear components
    ToCovariant_result = 0.0
    call ToCovariant(NextStress, ToCovariant_result)
    call GetFYieldFunction(ToCovariant_result, NextAlpha, ff)
    

    ! Check if correction is needed
    if (ff < TolF) return

    ! Stress correction loop
    do ii = 1, maxIter
        ! Compute dfrOverdSigma
        call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
        dfrOverdSigma = nn - (0.5 * DoubleDot2_2_Contr_result * I1)

        ! Compute lambda
        call DoubleDot2_2_Contr(dfrOverdSigma, RR, DoubleDot2_2_Contr_result)
        lambda = ff / DoubleDot2_2_Contr_result!DoubleDot2_2_Contr(dfrOverdSigma, R)

        ! Update NextStress and NextAlpha
        NextStress = NextStress - (RR * lambda)
        NextAlpha = NextAlpha - (dAlpha * lambda)
        
        ! Correct the stress
        ! I added this so that we can have the correct considertion of the two shear components
        ToCovariant_result = 0.0
        call ToCovariant(NextStress, ToCovariant_result)
        

        ! Recompute f
        call GetFYieldFunction(NextStress, NextAlpha, ff)

        ! Check convergence
        if (abs(ff) < TolF) exit

        ! Check maximum iterations
        !if (ii == maxIter) then
        !    if (debugFlag) print *, "Still outside with f = ", ff
        !endif
    end do

    end subroutine Stress_Correction
    
     
    !---------------------------------------------------------------
    !/*************************************************************/
    !/*************************************************************/
    !//            MATERIAL SPECIFIC METHODS                      //
    !/*************************************************************/
    !/*************************************************************/
    !// Macauley() -------------------------------------------------
    
    subroutine Macauley(xx, xx_Macauley) 
    !--------------------------------------------
    ! To do Macauley bracket. If it is negative, set to zero.
    !--------------------------------------------
    implicit none 
    real(REAL_TYPE), intent(in) :: xx
    real(REAL_TYPE), intent(out) :: xx_Macauley
    
    !// Macauley bracket
    if (xx>0) then 
        xx_Macauley = xx
    else 
        xx_Macauley = 0
    end if 
    
    end subroutine Macauley
    
    
    
    !---------------------------------------------------
    
    subroutine MacauleyIndex(xx, MacauleyIndex_result) 
    !--------------------------------------------
    ! To do Macauley bracket. If it is negative, set to zero.
    !--------------------------------------------
    implicit none 
    real(REAL_TYPE), intent(in) :: xx
    real(REAL_TYPE), intent(out) :: MacauleyIndex_result
    
    !// Macauley bracket
    if (xx>0) then 
        MacauleyIndex_result = 1
    else 
        MacauleyIndex_result = 0
    end if 
    
    end subroutine MacauleyIndex
    
    
    
    
   
    
    !-----------------------------------------------------
    
    subroutine GetFYieldFunction(nStress, nAlpha, ff) 
    !-----------------------------------------------------
    ! To evaluate the yield function. 
    !-----------------------------------------------------
    implicit none 
    
    ! input parameters 
    real(REAL_TYPE), intent(in), dimension(3) :: nAlpha
    real(REAL_TYPE), intent(in), dimension(3) :: nStress
    
    ! local parameters
    real(REAL_TYPE), dimension(3) :: ss ! deviator stress
    real(REAL_TYPE) :: pp
    
    ! output parameters
    real(REAL_TYPE) :: ff ! scalar
    
    ! Local parameters 
    real(REAL_TYPE) :: GetTrace_result
    real(REAL_TYPE) :: GetNorm_Contr_result
    real(REAL_TYPE), dimension(3) :: GetDevPart_result
    
    real(REAL_TYPE) :: root12
    real(REAL_TYPE) :: m_m
    
    real(REAL_TYPE) :: m_pmin
    
    
    root12 = 0.70710678118
    
    !// PM4Sand's yield function
	!// s = nStress - p*I
	!//   = [s_xx	s_xy ;	=	[Sigma_xx-p		Sigma_xy   ;
	!//		s_xy	s_yy];		 Sigma_xy		Sigma_yy-p];
	!// f = [(s - p*alpha):(s - p*alpha)]^0.5 - sqrt(1/2)*p*m
    
    ! get mean effective stress 
    call GetTrace(nStress, pp)
    pp = 0.5 * pp
    
    ! get deviatoric part of the stress by subtracting the mean eff stress from the xx and yy components
    call GetDevPart(nStress, GetDevPart_result)
    ss = GetDevPart_result !GetDevPart(nStress)
    
    ! subtract the stress ratio
    ss = ss - (pp*nAlpha)
    
    ! dot ss by ss
    call GetNorm_Contr(ss, GetNorm_Contr_result)
    
    ! evaluate the yield function ff
    ff = GetNorm_Contr_result - (root12 * m_m * pp)
    
    end subroutine GetFYieldFunction
    
    
    !----------------------------------------------------
    
    subroutine GetKsi(Dr, pp, RR_Bolton, QQ_Bolton, Pmin, P_atm, ksi)
    !-------------------------------------------
    ! To calculate the relative state parameter 
    ! based on Bolton (1986) as used by BZ14. .
    !-------------------------------------------
    implicit none
    
    ! input parameters
    real(REAL_TYPE), intent(in) :: Dr
    real(REAL_TYPE), intent(in) :: pp
    
    real(REAL_TYPE), intent(in)  :: RR_Bolton ! secondary parameter to find Dr,cs
    real(REAL_TYPE), intent(in)  :: QQ_Bolton ! secondary parameter to find Dr,cs
    
    real(REAL_TYPE), intent(in) :: Pmin
    real(REAL_TYPE), intent(in)  :: P_atm ! atmospheric pressure
    
    ! output parameters 
    real(REAL_TYPE), intent(out) :: ksi
    
    ! inout parameters
    
    ! local parameters
    real(REAL_TYPE) :: pn
    
    
    ! store initial 
    pn = pp
    
    ! cut off
    if (pn<=Pmin) then 
        pn = Pmin
    end if 
    
    !// find ksi_R
    !// calculating the relative density based on Bolton (1986)
    ksi = (RR_Bolton / (QQ_Bolton - log(100.0 * pn / P_atm)) ) - Dr
    
    end subroutine GetKsi
    
    
    
    
    !-------------------------------------------------
    subroutine GetElasticModuli_(sigma, zcum, z_max, nu, G0, Md, Mb, PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, Fsed_min, me2p)
    !--------------------------------------------
    ! To calculate the elastic modulus of the soil
    !--------------------------------------------
    implicit none 
    
    ! input variables 
    real(REAL_TYPE), intent(in), dimension(3) :: sigma(3)
    real(REAL_TYPE), intent(in) :: zcum
    real(REAL_TYPE), intent(in) :: z_max
    real(REAL_TYPE), intent(in) :: nu 
    real(REAL_TYPE), intent(in) :: G0
    real(REAL_TYPE), intent(in) :: Md
    real(REAL_TYPE), intent(in) :: Mb
    real(REAL_TYPE), intent(in) :: PostShake
    real(REAL_TYPE), intent(in) :: Pmin
    real(REAL_TYPE), intent(in) :: P_atm
    real(REAL_TYPE), intent(in) :: Cgd
    
    real(REAL_TYPE), intent(in) :: me2p
    
    real(REAL_TYPE), intent(in) :: p_sedo
    real(REAL_TYPE), intent(in) :: Fsed_min
    
    ! output variables 
    real(REAL_TYPE), intent(out) :: KK
    real(REAL_TYPE), intent(out) :: GG
    real(REAL_TYPE), intent(out) :: Mcur
    
    ! local variables
    real(REAL_TYPE) :: pp 
    real(REAL_TYPE) :: p_sed
    real(REAL_TYPE) :: F_sed
    
    real(REAL_TYPE) :: pn
    real(REAL_TYPE) :: qn
    real(REAL_TYPE) :: Csr0
    real(REAL_TYPE) :: Csr
    real(REAL_TYPE) :: msr
    
    real(REAL_TYPE) :: two3
    
    real(REAL_TYPE) :: temp
    
    real(REAL_TYPE) :: Macauley_result
    
    ! establish numeric values associated with this subroutine
    two3 = 0.666666666667
    msr = 4.0
    Csr0 = 0.5
    
    call GetTrace(sigma, pn)
    pn = 0.5 * pn
    
    if (pn <= Pmin) then 
        pn = Pmin 
    end if 
    
    qn = 2 * sqrt( (0.5*(sigma(1) - sigma(2)))**2.0 + sigma(3)**2.0)
    
    Mcur = qn/pn
    Csr = 1 - ( Csr0 * min(1.0, (Mcur/Mb)**msr) ) 
    
    temp = zcum / z_max
    
    if (me2p == 0) then 
    
        GG = G0 * P_atm !--> we never use this --> you can use this in stress initialization
    
    else ! me2p == .true.
        
        GG = G0 * P_atm * sqrt(pn / P_atm) * Csr * (1 + temp) / (1 + temp * Cgd)     
        
        !// this equation is in the table 
        !// "Elastic deviatoric strain increment"
        
        if (PostShake==1.0) then ! consider changing this to logical (love is never logical...)
			
            !// reduce elastic shear modulus for post shaking consolidation
            call GetTrace(sigma, pp)
			pp = 0.5 * pp
            
            call Macauley(1 - (Mcur/Md), Macauley_result)
			p_sed = p_sedo * (zcum / (zcum + z_max)) * Macauley_result**0.25 !// Equation 88
			F_sed = min(Fsed_min + (1 - Fsed_min) * (pp / 20.0 / (p_sed + small)), 1.0) !// Equation 87
			GG = GG * F_sed !// Equation 85
            
		end if 
        
    end if 
    
    
    KK = two3 * ( (1+nu)/(1 - (2*nu)) ) * GG
    !// Equation 86 (kind of...) 
	!// In the table "Elastic volumetric strain increment"
    
    !//
	!// K = (2/3) * ( 2*(1+nu)/3(1-2*nu) ) * G "Elastic volumetric strain increment"
	!//
    
    end subroutine GetElasticModuli_
    
    
    !---------------------------------------------------
    
    
    subroutine GetElasticModuli(sigma, G0, nu, P_atm, Pmin, GG, KK)
    !---------------------------------------------------
    ! To get the elastic parameters
    !---------------------------------------------------
    implicit none 
    
    ! input
    real(REAL_TYPE), intent(in), dimension(3) :: sigma
    real(REAL_TYPE), intent(in) :: G0    
    real(REAL_TYPE), intent(in) :: nu
    real(REAL_TYPE), intent(in) :: P_atm
    real(REAL_TYPE), intent(in) :: Pmin
    
    !logical, intent(in) :: me2p

    ! output
    real(REAL_TYPE), intent(out) :: GG
    real(REAL_TYPE), intent(out) :: KK
    
    
    ! Local variables 
    real(REAL_TYPE) :: pn
    real(REAL_TYPE) :: two3
    real(REAL_TYPE) :: GetTrace_sigma
    
    ! numeric values in this subroutine
    two3 = 2.0/3.0
    
    ! get mean effective stress
    call GetTrace(sigma, GetTrace_sigma)
    pn = 0.5*GetTrace_sigma
    
    ! tension cutoff
    if (pn<=Pmin) then 
        pn = Pmin
    end if 
 
    !if (me2p == .false.) then 
        
    ! shear modulus
    GG = G0 * P_atm ! .false --> linear elastic.
    !else 
    !    GG = G0 * P_atm * sqrt(pn / P_atm) ! .true. --> nonlinear elastic (do not use this here)
    !end if 

    ! bulk modulus 
    KK = two3 * ( (1+nu)/(1-(2*nu)) ) * GG !--> two3 needs to be a global number
    
    end subroutine GetElasticModuli                 
                     
                     
                     
                     
    !-------------------------------------------------                
    subroutine GetStiffness(KK, GG, CC)
    
    implicit none
    
    ! input variables 
    real(REAL_TYPE), intent(in) :: KK, GG
    
    ! output variables 
    real(REAL_TYPE), intent(out) :: CC(3, 3) ! 3x3 2D plane strain elasticity matrix
    
    ! local variables
    real(REAL_TYPE) :: aa, bb

    ! compute coefficients of the elasticity matrix
    aa = KK + ((4.0 / 3.0) * GG)
    bb = KK - ((2.0 / 3.0) * GG)

    ! fill stiffness matrix
    CC(1, 1) = aa
    CC(2, 2) = aa
    CC(3, 3) = GG
    CC(1, 2) = bb
    CC(2, 1) = bb


    end subroutine GetStiffness              
    
    
    
    
    !------------------------------------------------                 
    !subroutine GetCompliance()                 
    !               
    !implicit none
    !
    !
    !
    !
    !
    !end subroutine 
    
    
    !// returns the compliance matrix in its covariant-covariant form
    subroutine GetCompliance(KK, GG, DD)
    implicit none
    ! Inputs
    real(REAL_TYPE), intent(in) :: KK, GG
    ! Output
    real(REAL_TYPE), dimension(3, 3), intent(out) :: DD
    ! Local variables
    real(REAL_TYPE) :: aa, bb, cc

    ! Compute components of compliance matrix
    aa = (KK + 4.0 / 3.0 * GG) / (4.0 * GG * KK + 4.0 / 3.0 * GG**2.0)
    bb = (KK - 2.0 / 3.0 * GG) / (4.0 * GG * KK + 4.0 / 3.0 * GG**2.0)
    cc = 1.0 / GG

    ! Assign values to the compliance matrix
    DD(1, 1) = aa
    DD(2, 2) = aa
    DD(3, 3) = cc
    DD(1, 2) = bb
    DD(2, 1) = bb

end subroutine GetCompliance
    
    !------------------------------------------------
                     
    !subroutine GetElastoPlasticTangent()                 
    !          
    !implicit none
    !                 
    !              
    !
    !
    !end subroutine GetElastoPlasticTangent
    
    
    
    
    
!    subroutine GetElastoPlasticTangent(NextStress, aCe, RR, nn, K_p, aCep)
!    
!    ! --> I dont know what this is doing.....
!    implicit none
!    ! Inputs
!    real(REAL_TYPE), dimension(3), intent(in) :: NextStress
!    real(REAL_TYPE), dimension(3, 3), intent(in) :: aCe
!    real(REAL_TYPE), dimension(3), intent(in) :: RR
!    real(REAL_TYPE), dimension(3), intent(in) :: nn
!    real(REAL_TYPE), intent(in) :: K_p
!    ! Output
!    real(REAL_TYPE), dimension(3, 3), intent(out) :: aCep
!    ! Local variables
!    real(REAL_TYPE) :: pp, temp3
!    real(REAL_TYPE), dimension(3) :: rrr, temp1, temp2
!
!    ! Local variables 
!    real(REAL_TYPE) :: GetTrace_result
!    real(REAL_TYPE), dimension(3) :: GetDevPart_result
!    real(REAL_TYPE), dimension(3) :: DoubleDot4_2_result
!    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
!    real(REAL_TYPE) :: m_pmin
!    
!    real(REAL_TYPE), dimension(3) :: DoubleDot2_4_result
!    
!    integer(INTEGER_TYPE), dimension(3) :: mI1
!    
!    real(REAL_TYPE), dimension(3,3) :: Dyadic2_2_result
!    
!    real(REAL_TYPE), dimension(3,3) :: mIIco
!    
!    mI1(1) = 1
!    mI1(2) = 1
!    mI1(3) = 1
!    
!    mIIco(1,1) = 1
!    mIIco(2,1) = 1
!    mIIco(3,1) = 1
!    
!    mIIco(1,2) = 1
!    mIIco(2,2) = 1
!    mIIco(3,2) = 1
!    
!    mIIco(1,3) = 1
!    mIIco(2,3) = 1
!    mIIco(3,4) = 1
!    
!    
!    
!    
!    
!    
!    
!    
!    ! Calculate mean effective stress
!    call GetTrace(NextStress, pp)
!    pp = 0.5 * pp
!    if (pp < Pmin) pp = Pmin
!
!    ! Calculate deviatoric part of stress and normalize
!    call GetDevPart(NextStress, rrr)
!    rrr = rrr / pp
!
!    ! Compute temporary vectors and matrices
!    call DoubleDot4_2(aCe, RR, DoubleDot4_2_result)
!    temp1 = DoubleDot4_2_result!DoubleDot4_2(aCe, RR, DoubleDot4_2_result)
!    
!    call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
!    call DoubleDot2_4((nn - 1.0 / 2.0 * DoubleDot2_2_Contr_result) * mI1, aCe * mIIco, DoubleDot2_4_result)
!    temp2 = DoubleDot2_4_result !DoubleDot2_4(n - 1.0 / 2.0 * DoubleDot2_2_Contr(n, r) * mI1, aCe * mIIco)
!    
!    call DoubleDot2_2_Contr(temp2, RR, DoubleDot2_2_Contr_result)
!    temp3 = DoubleDot2_2_Contr_result + K_p
!
!    ! Update aCep based on plasticity condition
!    if (temp3 < small) then
!        aCep = aCe
!    else
!        call Dyadic2_2(temp1, temp2, Dyadic2_2_result)
!        aCep = aCe - 1.0 / temp3 * Dyadic2_2_result!Dyadic2_2(temp1, temp2)
!    end if
!end subroutine GetElastoPlasticTangent
    
    
    
    
    
    
    
    
    !-----------------------------------------------
    subroutine GetNormalToYield(stress, alpha, nn)
    !-----------------------------------------------
    ! To calculate the normal to the yield surface
    !-----------------------------------------------
    implicit none 
    
    ! output variables
    real(REAL_TYPE), intent(out), dimension(3) :: nn
    
    ! input variables 
    real(REAL_TYPE), intent(in), dimension(3) :: stress, alpha
    
    ! local variables
    real(REAL_TYPE) :: normN 
    real(REAL_TYPE) :: pp
    real(REAL_TYPE) :: small
    real(REAL_TYPE) :: root12
    real(REAL_TYPE), dimension(3) :: GetDevPart_stress
    
    ! define numeric values
    root12 = 0.7071067811865475
    
    call GetTrace(stress, pp)
    pp = 0.5 * pp
    
    call GetDevPart(stress, GetDevPart_stress)
    
    if (abs(pp)<small) then
        
        nn(2) = root12
        
    else 
        
        nn = (alpha*-pp) + GetDevPart_stress
        call GetNorm_Contr(nn, normN)
        
        if (normN<small) then 
            normN = 1
        end if
        nn = nn * (1/normN)
        
    end if
    
    
    end subroutine GetNormalToYield
    
    
    
    
    
    
    
    
    !-------------------------------------------------------------------
    
  
    
    
    
    !*************************************************************
    ! GetStateDependent() ----------------------------------------
    ! Alsardi: we are calculating the state parameters for PM4Sand
    ! mainly:
    ! -Ksi
    ! -Surface Ms
    ! -Plastic modulus
    ! -Dlatancy parameters
    ! -R (plastic potential derivative) --> m in Alba's notes
    subroutine GetStateDependent(stress, alpha, alpha_in, alpha_in_p, alpha_in_true, fabric, fabric_in, &
        GG, zcum, zpeak, pzp, Mcur, CurDr, &
        Mc, nd, nb, Pmin, Pmin2, P_atm, & 
        mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
        QQ_Bolton, RR_Bolton, &
        nn, alphaB, alphaD, bb, RR, &
        Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, ksi)
    ! get all the state parameters associated with a stress state
    implicit none
    ! Input parameters
    real(REAL_TYPE), dimension(3), intent(in) :: stress
    real(REAL_TYPE), dimension(3), intent(in) :: alpha
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in_p
    real(REAL_TYPE), dimension(3), intent(in) :: alpha_in_true
    real(REAL_TYPE), dimension(3), intent(in) :: fabric
    real(REAL_TYPE), dimension(3), intent(in) :: fabric_in
    real(REAL_TYPE), intent(in) :: GG
    real(REAL_TYPE), intent(in) :: zcum
    real(REAL_TYPE), intent(in) :: zpeak
    real(REAL_TYPE), intent(in) :: pzp
    real(REAL_TYPE), intent(in) :: Mcur
    real(REAL_TYPE), intent(in) :: CurDr
    
    real(REAL_TYPE), intent(in) :: Mc
    real(REAL_TYPE), intent(in) :: nd
    real(REAL_TYPE), intent(in) :: nb
    real(REAL_TYPE), intent(in) :: Pmin
    real(REAL_TYPE), intent(in) :: Pmin2
    real(REAL_TYPE), intent(in) :: P_atm
    real(REAL_TYPE), intent(in) :: mm
    real(REAL_TYPE), intent(in) :: z_max
    real(REAL_TYPE), intent(in) :: h0
    real(REAL_TYPE), intent(in) :: Ckaf
    real(REAL_TYPE), intent(in) :: Ado
    real(REAL_TYPE), intent(in) :: Ceps
    real(REAL_TYPE), intent(in) :: hpo
    real(REAL_TYPE), intent(in) :: Cdr
    
    real(REAL_TYPE), intent(in) :: QQ_Bolton
    real(REAL_TYPE), intent(in) :: RR_Bolton
    
    ! Output parameters
    real(REAL_TYPE), dimension(3), intent(out) :: nn
    real(REAL_TYPE), dimension(3), intent(out) :: alphaB
    real(REAL_TYPE), dimension(3), intent(out) :: alphaD
    real(REAL_TYPE), dimension(3), intent(out) :: bb
    real(REAL_TYPE), dimension(3), intent(out) :: RR
    
    real(REAL_TYPE), intent(out) :: ksi
    
    real(REAL_TYPE), intent(out) :: Mb
    real(REAL_TYPE), intent(out) :: Md
    real(REAL_TYPE), intent(out) :: DD
    real(REAL_TYPE), intent(out) :: K_p
    real(REAL_TYPE), intent(out) :: Cka
    real(REAL_TYPE), intent(out) :: hh
    real(REAL_TYPE), intent(out) :: AlphaAlphaBDotN
    
    ! Local variables
    real(REAL_TYPE), dimension(3) :: alphaD_alpha 
    real(REAL_TYPE), dimension(3) :: alphaDr_alpha
    real(REAL_TYPE), dimension(3) :: alpha_Alpha_in
    real(REAL_TYPE), dimension(3) :: alpha_Alpha_in_true
    real(REAL_TYPE), dimension(3) :: alpha_Alpha_p
    real(REAL_TYPE), dimension(3) :: minusFabric
    
    real(REAL_TYPE) :: Czpk1
    real(REAL_TYPE) :: Czpk2
    real(REAL_TYPE) :: Cpzp
    real(REAL_TYPE) :: Cpzp2
    real(REAL_TYPE) :: Cg1
    real(REAL_TYPE) :: Ckp
    
    real(REAL_TYPE) :: AlphaAlphaInDotN
    real(REAL_TYPE) :: AlphaAlphaInTrueDotN
    real(REAL_TYPE) :: Czin1
    real(REAL_TYPE) :: Czin2
    real(REAL_TYPE) :: Crot1
    real(REAL_TYPE) :: Mdr
    real(REAL_TYPE) :: pp
    
    
    real(REAL_TYPE) :: Cpmin
    
    real(REAL_TYPE) :: Crot2
    real(REAL_TYPE) :: C_pmin2 
    real(REAL_TYPE) :: Cin 
    real(REAL_TYPE) :: Cdz
    
    real(REAL_TYPE) :: temp
    real(REAL_TYPE) :: Ad
    
    real(REAL_TYPE) :: hp 
    real(REAL_TYPE) :: Adc
    
    real(REAL_TYPE) :: Drot
    
    ! Local variables
    real(REAL_TYPE) :: Macauley_result 
    real(REAL_TYPE) :: GetTrace_result 
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result 
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result2

    real(REAL_TYPE) :: root12
    real(REAL_TYPE) :: two3
    real(REAL_TYPE) :: one3
    real(REAL_TYPE), dimension(3) :: I1
   
    ! identity matrix
    I1(1) = 1
    I1(2) = 1
    I1(3) = 0
    
    ! two thirds 
    two3 = 0.66666666667
    ! one third
    one3 = 0.33333333333
    
    ! sqrt(0.5)
    root12 = 0.70710678118
    
    ! Mean effective stress
    call GetTrace(stress, pp)
    pp = 0.5d0 * pp
    
    if (pp <= Pmin) pp = Pmin ! kind of a tension cutoff
    
    ! evaluate relative state parameter index based on Bolton (1986) relationship
    call GetKsi(CurDr, pp, RR_Bolton, QQ_Bolton, Pmin, P_atm, ksi) ! 
    
    ! get the normal to the yield surface
    call GetNormalToYield(stress, alpha, nn) 
    
    ! update bounding and dilatancy using the relative state parameter (ksi)
    if (ksi <= 0.0d0) then ! negative relative state parameter index
        ! dense of critical
        Mb = Mc * exp(-1.0d0 * nb * ksi) ! bounding ratio based on ksi using m_nb
        Md = Mc * exp(         nd * ksi) ! dilatancy ratio update based on ksi using m_nd
    else ! positive relative state parameter index
        ! loose of critical
        Mb = Mc * exp(-1.0d0 * (nb / 4.0d0) * ksi) ! bounding ratio based on ksi using m_nb/4
        Md = Mc * exp(         (nd * 4.0d0) * ksi) ! dilatancy ratio update based on ksi using m_nb*4
    end if
    
    ! update image back-stress ratios for bounding surface
    alphaB = nn * root12 * (Mb - mm)
    
    ! update image back-stress ratios for dilatancy surface
    alphaD = nn * root12 * (Md - mm)
    
    ! calculate plastic modulus terms
    ! C_zpk1 = z_peak / (z_cum + z_max/5)
    Czpk1 = zpeak / (zcum + (z_max / 5.0d0)) ! Equation 62 --> Effect of fabric on plastic modulus
    
    ! C_zpk2 = z_peak / (z_cum + z_max/100)
    Czpk2 = zpeak / (zcum + (z_max / 100.0d0)) ! Equation 63 --> Effect of fabric on plastic modulus
    
    if (Czpk2 > 1.0d0 - small) Czpk2 = 1.0d0 - small ! there is a limit on this
    
    ! C_pzp2 = (-<-(p_zp - p)>)/(-<-(p_zp - p)> + p_min)
    call Macauley((pzp - pp), Macauley_result)
    !Cpzp2 = Macauley((pzp - pp)) / (Macauley((pzp - pp)) + m_Pmin) ! Equation 64 --> Effect of fabric on plastic modulus
    Cpzp2 = Macauley_result / (Macauley_result + Pmin) ! Equation 64 --> Effect of fabric on plastic modulus

    Cg1 = h0 / 200.0d0 ! look at the table in "Plastic deviatoric strain increment" C_Gamma1 = h_o / 200
    Ckp = 2.0d0 ! look at the table in "Plastic deviatoric strain increment" C_Kp = 2
    
    ! terms for k_p calculation
    bb = alphaB - alpha ! difference between the bounding surface image back stress ratio and the back stress ratio
    
    call DoubleDot2_2_Contr(bb, nn, AlphaAlphaBDotN)
    
    ! differnce between the back stress ratio and the initial back stress ratio
    alpha_Alpha_in = alpha - Alpha_in ! alpha - alpha_in
    
    ! dot alpha_mAlpha_in with normal nn
    call DoubleDot2_2_Contr(alpha_Alpha_in, nn, DoubleDot2_2_Contr_result)
    call Macauley(DoubleDot2_2_Contr_result, AlphaAlphaInDotN) ! Macauly bracket in the Kp calculation
    
    Alpha_Alpha_in_true = alpha - Alpha_in_true ! Alpha - Alpha^true_in
    
    call DoubleDot2_2_Contr(alpha_Alpha_in_true, nn, DoubleDot2_2_Contr_result)
    call Macauley(DoubleDot2_2_Contr_result, AlphaAlphaInTrueDotN)
    ! = Macauley_result!Macauley(DoubleDot2_2_Contr(alpha_mAlpha_in_true, n)) ! dot product in the denominator of Cka
    
    Cka = 1.0d0 + ( (Ckaf / (1.0d0 + (2.5d0 * AlphaAlphaInTrueDotN)**2)) * Cpzp2 * Czpk1 ) ! look at table "Plastic deviatoric strain increment"
    ! updataed K_p formulation following PM4Sand V3.1. mAlpha_in is the apparent back-stress ratio.
    alpha_Alpha_p = alpha - alpha_in_p ! if alpha - alpha_in_p is small
    
    call DoubleDot2_2_Contr(alpha_Alpha_p, nn, DoubleDot2_2_Contr_result)
    call Macauley(AlphaAlphaBDotN, Macauley_result)
      
    if (abs(AlphaAlphaBDotN) < small) then
        ! adding this condition to avoid division by zero error
        hh = 1.0d10 ! part of Kp calculation
    elseif (DoubleDot2_2_Contr_result <= 0) then !--> alpha_mAlpha_p dot nn <=0
        ! LOAD REVERSAL
        hh = (1.5d0 * GG * h0) *sqrt(abs(AlphaAlphaBDotN)) /(pp *(exp(AlphaAlphaInDotN) - 1 + Cg1) ) * &
            Cka / (1 + Ckp * (zpeak / z_max) * Macauley_result * sqrt(1 - Czpk2)) ! part of Kp calculation 
        ! Equation 35 --> Hardening and the update of the back-stress ratio
        ! This corresponds to equation 3.31 in Chen's PhD thesis to prevent shear stress overshooting
        hh = hh * (AlphaAlphaInDotN + Cg1) / (AlphaAlphaInTrueDotN + Cg1) ! part of Kp calculation --> Cg1 to prevent dividing by zero
    else 
        ! NO LOAD REVERSAL
        hh = (1.5d0 * GG * h0) * sqrt(abs(AlphaAlphaBDotN)) /(pp * (exp(AlphaAlphaInDotN) - 1 + Cg1) ) * &
            Cka / (1 + Ckp * (zpeak / z_max) * Macauley_result * sqrt(1 - Czpk2)) ! part of Kp calculation
    end if
        
    ! after we evaluate h we plug it in K_p to find the plastic bulk modulus
    call DoubleDot2_2_Contr(bb, nn, DoubleDot2_2_Contr_result)
    K_p = two3 * hh * pp * DoubleDot2_2_Contr_result !DoubleDot2_2_Contr(b, n) ! Equation 34 --> Hardening and the update of the back-stress ratio

    !--- calculate dilatancy
    ! For Czin1 calculation in Equation 67
    call DoubleDot2_2_Contr(fabric_in, nn, DoubleDot2_2_Contr_result)
    call DoubleDot2_2_Contr(fabric, nn, DoubleDot2_2_Contr_result2)

    call Macauley(1.0 - exp(-2.0*abs((DoubleDot2_2_Contr_result - DoubleDot2_2_Contr_result2) / z_max)), Czin1)
    !Czin1 = Macauley(1.0 - exp(-2.0*abs((DoubleDot2_2_Contr(fabric_in, nn) - DoubleDot2_2_Contr(fabric, n)) / m_z_max))) ! Plastic volumetric strain increment... see in table
    ! = Macauley_result! Plastic volumetric strain increment... see in table
    ! rotated dilatancy surface
    minusFabric = fabric * (-1.0)
    call DoubleDot2_2_Contr(minusFabric, nn, DoubleDot2_2_Contr_result)
    call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
    Crot1 = max(1.0 + (2.0 * Macauley_result / (sqrt(2.0)*z_max)) * (1 - Czin1), 1.0) ! Plastic volumetric strain increment... see in table
    
    ! Equation 66 --> Effect of fabric on plastic volumetric dilation 
    Mdr = Md / Crot1 ! Equation 65 --> Effect of fabric on plastic volumetric dilation 
    
    ! Vector alphaDr = root12 * (Mdr - m_m) * n;
    !alphaDr_alpha = nn
    alphaDr_alpha = nn * (root12 * (Mdr - mm)) ! alphaDr
    alphaDr_alpha = alphaDr_alpha - alpha ! Equation 68
    
    
    ! Alpha_DR = sqrt(1/2) * (M^dr - m)n
    !alphaD_alpha = alphaD
    alphaD_alpha = alphaD - alpha
    call DoubleDot2_2_Contr(alphaDr_alpha, nn, DoubleDot2_2_Contr_result)
    
    if (DoubleDot2_2_Contr_result <= 0) then ! AlphaD - Alpha : nn --> negative --> DILATION 
        
        ! dilation
        if (pzp == 0.0) then
            Cpzp = 1.0 ! to prevent dividing by zero
        else
            Cpzp = 1.0 / (1.0 + (2.5 * pp/pzp)**5.0) ! plastic vol strain increment ! Equation 3.55
        end if
        
        Cpmin = 1.0 / (1.0 + (Pmin2/pp)**2) ! plastic vol strain increment ! Equation 3.56
        Czin2 = (1.0 + (Czin1 * (zcum - zpeak)) / (3.0*z_max)) / (1.0 + (3.0 * Czin1 * (zcum - zpeak)) / (3.0*z_max)) ! plastic vol strain increment ! Equation 3.60
        
        ! double temp = (1.0 - Macauley(DoubleDot2_2_Contr(-1.0 * fabric, n)) * root12 / zpeak) ** 3
        call DoubleDot2_2_Contr(minusFabric, nn, DoubleDot2_2_Contr_result)
        call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
        temp = (1.0 - Macauley_result * root12 / zpeak) ** 3.0
        
        Ad = Ado * Czin2 / (((zcum**2 / z_max) * temp * (Ceps**2) * Cpzp * Cpmin * Czin1) + 1.0) ! plastic vol strain increment ! Equation 72
        call DoubleDot2_2_Contr(alphaD_alpha, nn, DoubleDot2_2_Contr_result)
        
        ! DD_non-rot in equation 70
        DD = Ad * DoubleDot2_2_Contr_result!DoubleDot2_2_Contr(alphaD_alpha, n) ! D_non-rot ! plastic vol strain increment
        ! double Drot = Ad * Macauley(DoubleDot2_2_Contr(-1.0*fabric, n)) / (sqrt(2.0)*m_z_max) * DoubleDot2_2_Contr(alphaDr - alpha, n) / m_Cdr
        
        call DoubleDot2_2_Contr(minusFabric, nn, DoubleDot2_2_Contr_result)
        call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
        
        call DoubleDot2_2_Contr(alphaDr_alpha, nn, DoubleDot2_2_Contr_result2)
        !call Macauley(DoubleDot2_2_Contr_result, Macauley_result2)
        !Drot = Ad * Macauley(DoubleDot2_2_Contr(minusFabric, n)) / (sqrt(2.0) * m_z_max) * DoubleDot2_2_Contr(alphaDr_alpha, n) / m_Cdr ! D_rot ! plastic vol strain increment
        !Drot --> Equation 69
        Drot = Ad * ( Macauley_result / (sqrt(2.0) * z_max) ) * ( DoubleDot2_2_Contr_result2 / Cdr )! D_rot ! plastic vol strain increment
        
        ! Equation 69 --> Drot
        if (DD > Drot) then ! plastic vol strain increment --> D_non-rot < D_rot
            call Macauley(Mb - Mcur, Macauley_result)
            DD = DD + (Drot - DD) * Macauley_result / (Macauley_result + 0.01) ! Equation 78 ! slight descrepancy in the denominator
        end if
        
        ! during dilation at very low effective stresses (i.e., p<=2p_min), Di is constrained to ensure soil at 
        ! dense than critical states continue to be dliative 
        if ((Pmin <= pp) .and. (pp <= 2*Pmin)) then
            call Macauley(Mb - Md, Macauley_result)
            ! equation 3.61 in Long Chen's PhD thesis
            DD = min(DD, -3.5 * Ado * Macauley_result * (2*Pmin - pp)/ Pmin) ! Equation 78 ! Equation 3.61
        end if
        
        
    else
        
        ! contraction
        
        ! bound K_p to non - negative, following flac practice
        K_p = max(0.0, K_p)
        ! he is using power 2 instead of 2.5 in BZ14
        call Macauley(0.5 - ksi, Macauley_result)
        hp = hpo * exp(-0.7 + (7.0*(Macauley_result)**2.0)) ! Equation 53 --> Plastic volumetric strains - Contraction ! Equation 3.41
        Crot2 = 1.0 - Czpk2 ! Equation 83
        ! Equation 3.664 in Long PhD thesis
        Cdz = max((1.0 - Crot2 * sqrt(2.0) * zpeak / z_max) * (z_max / (z_max + Crot2 * zcum)), 1.0 / (1.0 + z_max/2.0)) ! Equation 3.64
        
        call DoubleDot2_2_Contr(fabric, nn, DoubleDot2_2_Contr_result)
        call Macauley(DoubleDot2_2_Contr_result, Macauley_result)
        
        ! Equation 3.63 in Long PhD thesis
        !Adc = m_Ado * (1.0 + Macauley(DoubleDot2_2_Contr(fabric, n))) / hp / Cdz
        Adc = Ado * (1.0 + Macauley_result) / (hp * Cdz) ! Equation 3.63

        !Cin = 2.0 * Macauley(DoubleDot2_2_Contr(fabric, n)) / sqrt(2.0) / m_z_max
        Cin = 2.0 * Macauley_result / (sqrt(2.0) * z_max) ! Equation 3.66

        ! D = min(Adc * (DoubleDot2_2_Contr(alpha - mAlpha_in, n) + Cin) ** 2, 1.5 * m_Ado) * &
        !     DoubleDot2_2_Contr(alphaD - alpha, n) / (DoubleDot2_2_Contr(alphaD - alpha, n) + 0.16)
        call DoubleDot2_2_Contr(alpha_Alpha_in, nn, DoubleDot2_2_Contr_result)
        call DoubleDot2_2_Contr(alphaD_alpha, nn, DoubleDot2_2_Contr_result2)
        DD = min(Adc*(DoubleDot2_2_Contr_result + Cin)**2.0, 1.5 * Ado) * &
            DoubleDot2_2_Contr_result2 / (DoubleDot2_2_Contr_result2 + 0.16) ! note that the threshold is enforced by using min function
        ! Apply a factor to D so it doesn't go very big when p is small
        ! Equation 3.67 in Long PhD Thesis --> updated in PM4Sand 3.3
        if (pp < Pmin * 2.0) then
            C_pmin2 = 0.0
        elseif (pp >= Pmin * 18.0) then
            C_pmin2 = 1.0
        else
            C_pmin2 = (pp - (2.0 * Pmin)) / (16.0 * Pmin)
        end if
        
        DD = DD * C_pmin2
        
    end if
    
    ! R = n + one3 * D * mI1
    RR = nn + ((one3 * DD)*I1) ! --> used in plastic deviatoric strain increment
    
    ! I do not see where we calculate the deviatoric component of RR to get RR'?
    
    end subroutine GetStateDependent
    
    
    
    
    
    
    !/*************************************************************/
    !/*************************************************************/
    !//            SYMMETRIC TENSOR OPERATIONS                    //
    !/*************************************************************/
    !/*************************************************************/
    !// In all the functions below, by contravariant tensors, we mean stress-like tensors
    !// and by covariant tensors we mean strain-like tensors
    subroutine GetTrace(vv, GetTrace_result)
    !// computes the trace of the input argument
    implicit none
    ! Arguments
    real, dimension(3), intent(in) :: vv
    real, intent(out) :: GetTrace_result
    
    ! initialize
    GetTrace_result = 0.0

    ! Check vector size
    if (size(vv) /= 3) then
        print *, "ERROR! PM4Sand::GetTrace requires vector of size(3)!"
        GetTrace_result = 0.0
        return
    endif
    
    ! Compute the trace
    GetTrace_result = vv(1) + vv(2)

    end subroutine GetTrace
    
    
    
    subroutine GetDevPart(Vector, DevVector)
    !-------------------------------------------
    ! To calculate the deviatoric part of a matrix.
    !-------------------------------------------
    implicit none
    
    ! input variable 
    real(REAL_TYPE), dimension(3), intent(in) :: Vector
    real(REAL_TYPE), dimension(3), intent(out) :: DevVector
    
    ! local variables 
    real(REAL_TYPE) :: pp
    
    ! initialize 
    DevVector = 0.0
    pp = 0.0
    
    ! Check vector size
    if (size(Vector) /= 3) then
        print *, "ERROR! PM4Sand::GetDevPart requires vector of size(3)!"
        !result = 0.0
        return
    endif
    
    call GetTrace(Vector, pp)
    
    DevVector = Vector 
    DevVector(1) = DevVector(1) - 0.5*pp
    DevVector(2) = DevVector(2) - 0.5*pp
    
    
    
    end subroutine GetDevPart
    
    
    
    
    
    subroutine DoubleDot2_2_Contr(v1, v2, DoubleDot2_2_Contr_result)
    !--------------------------------------------------
    ! // computes doubledot product for vector-vector arguments,
    !    both "contravariant".
    implicit none 
    
    ! the size of input vectors needs to be both 3
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(in) :: v2
    real(REAL_TYPE), intent(out) :: DoubleDot2_2_Contr_result
    integer(INTEGER_TYPE) :: ii
    
    ! initialize 
    DoubleDot2_2_Contr_result = 0.0

    ! Check vector sizes
    if (size(v1) /= 3 .or. size(v2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_2_Contr requires vectors of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    DoubleDot2_2_Contr_result = 0 
    
    do ii = 1, 3 ! 3 is for the size of the vector 
        !result = result + v1(i) * v2(i) + (i > 1) * v1(i) * v2(i)
        DoubleDot2_2_Contr_result = DoubleDot2_2_Contr_result + ( v1(ii) * v2(ii) ) !+ (ii>2)*( v1(ii) * v2(ii) )
        if (ii == 3) then 
            ! double count shear component according to cauchy stress tensor
            DoubleDot2_2_Contr_result = DoubleDot2_2_Contr_result + ( v1(ii) * v2(ii) )
        end if 
    end do 
    
    end subroutine DoubleDot2_2_Contr
    
    
    
    
    
    subroutine DoubleDot2_2_Cov(v1, v2, DoubleDot2_2_Cov_result)
    !--------------------------------------------------------
    ! // computes doubledot product for vector-vector arguments, both "covariant"
    !
    implicit none 
    ! the size of input vectors needs to be both 3
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(in) :: v2
    real(REAL_TYPE), intent(out) :: DoubleDot2_2_Cov_result
    integer(INTEGER_TYPE) :: ii
    
      
    ! Check vector sizes
    if (size(v1) /= 3 .or. size(v2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_2_Cov requires vectors of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    DoubleDot2_2_Cov_result = 0
    
    
    ! Compute the doubledot product
    do ii = 1, 3
        !result = result + v1(i) * v2(i) - (i > 2) * 0.5 * v1(i) * v2(i)
        DoubleDot2_2_Cov_result = DoubleDot2_2_Cov_result + (v1(ii) * v2(ii))
        if (ii == 3) then 
            DoubleDot2_2_Cov_result = DoubleDot2_2_Cov_result - (0.5 * v1(ii) * v2(ii))
        end if 
        
    end do 
    
    end subroutine DoubleDot2_2_Cov
    
    
    
    
    subroutine DoubleDot2_2_Mixed(v1, v2, DoubleDot2_2_Mixed_result)
    ! // computes doubledot product for vector-vector arguments, 
    ! one "covariant" and the other "contravariant" 
    implicit none 
    ! the size of input vectors needs to be both 3
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(in) :: v2
    real(REAL_TYPE), intent(out) :: DoubleDot2_2_Mixed_result
    integer(INTEGER_TYPE) :: ii
    
    ! Check vector sizes
    if (size(v1) /= 3 .or. size(v2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_2_Mixed requires vectors of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    DoubleDot2_2_Mixed_result = 0
    
    do ii = 1, 3
        DoubleDot2_2_Mixed_result = DoubleDot2_2_Mixed_result + ( v1(ii) * v2(ii) )
    end do
    
    end subroutine DoubleDot2_2_Mixed
    
    
    
    
    subroutine GetNorm_Contr(vv, GetNorm_Contr_result)
    ! // computes contravariant (stress-like) norm of input 6x1 tensor
    implicit none
    
    real(REAL_TYPE), dimension(3), intent(in) :: vv
    real(REAL_TYPE), intent(out) :: GetNorm_Contr_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Contr_result
    
    
    ! Check vector size
    if (size(vv) /= 3) then
        print *, "ERROR! PM4Sand::GetNorm_Contr requires vector of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    GetNorm_Contr_result = 0
    DoubleDot2_2_Contr_result = 0 
    
    call DoubleDot2_2_Contr(vv, vv, DoubleDot2_2_Contr_result)
    
    GetNorm_Contr_result = sqrt(DoubleDot2_2_Contr_result)
    
    end subroutine GetNorm_Contr
    
    
    
    
    
    
    
    
    
    
    subroutine GetNorm_Cov(vv, GetNorm_Cov_result)
    ! // computes covariant (strain-like) norm of input 6x1 tensor
    implicit none 
    
    real(REAL_TYPE), dimension(3), intent(in) :: vv
    real(REAL_TYPE), intent(out) :: GetNorm_Cov_result
    
    real(REAL_TYPE) :: DoubleDot2_2_Cov_result
    
    ! Check vector size
    if (size(vv) /= 3) then
        print *, "ERROR! PM4Sand::GetNorm_Cov requires vector of size(3)!"
        !result = 0.0
        return
    endif
    
    ! initialize integer
    GetNorm_Cov_result = 0
    DoubleDot2_2_Cov_result = 0 
    
    ! evaluate GetNorm_Cov_result using DoubleDot2_2_Cov
    call DoubleDot2_2_Cov(vv, vv, DoubleDot2_2_Cov_result)
    
    ! square root to find GetNorm_Cov
    GetNorm_Cov_result = sqrt(DoubleDot2_2_Cov_result)
    
    
    end subroutine GetNorm_Cov
    
    
    
    
    
    
    subroutine Dyadic2_2(v1, v2, Dyadic2_2_result)
    !// computes dyadic product for two vector-storage arguments
    !// the coordinate form of the result depends on the coordinate form of inputs
    ! https://en.wikipedia.org/wiki/Dyadics
    implicit none 
    
    ! input vectors v1 and v2
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(in) :: v2
    ! output matrix result Dyadic2_2_result
    real(REAL_TYPE), dimension(3,3), intent(out) :: Dyadic2_2_result
    
    ! local variables 
    integer(INTEGER_TYPE) :: jj, ii
    
    ! Check vector sizes
    if (size(v1) /= 3 .or. size(v2) /= 3) then
        print *, "ERROR! PM4Sand::Dyadic2_2 requires vector of size(3)!"
        return
    endif
    
    ! initialize matrix
    Dyadic2_2_result = 0.0
    
    do ii = 1, 3
        do jj = 1, 3
            
            Dyadic2_2_result(ii,jj) = v1(ii) * v2(jj)
            
        end do 
    end do 
    
    
    
    
    end subroutine Dyadic2_2
    
    
    
    
    subroutine DoubleDot4_2(m1, v1, DoubleDot4_2_result)
    !// computes doubledot product for matrix-vector arguments
    !// caution: second coordinate of the matrix should be in opposite variant form of vector
    implicit none 
    ! input variables 
    real(REAL_TYPE), intent(in), dimension(3,3) :: m1
    real(REAL_TYPE), intent(in), dimension(3) :: v1
    ! output variables 
    real(REAL_TYPE), intent(out), dimension(3) :: DoubleDot4_2_result!m1Timesv1
    ! local variables 
    integer(INTEGER_TYPE) :: ii
    integer(INTEGER_TYPE) :: jj
    
    ! Check vector size
    if (size(v1) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot4_2 requires vector of size(3)!"
        return
    endif

    ! Check matrix size
    if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot4_2 requires 3-by-3 matrix "
        return
    endif
    
    
    ! Compute the result
    do ii = 1, 3
        DoubleDot4_2_result(ii) = 0.0
        do jj = 1, 3
            DoubleDot4_2_result(ii) = DoubleDot4_2_result(ii) + ( m1(ii, jj) * v1(jj) )
        enddo
    enddo
    
    end subroutine DoubleDot4_2
    
    
    
    
    
    
    
    
    
    
    subroutine DoubleDot2_4(v1, m1, DoubleDot2_4_result)
    !// computes doubledot product for matrix-vector arguments
    !// caution: first coordinate of the matrix should be in opposite 
    !// variant form of vector
    implicit none
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3, 3), intent(in) :: m1
    real(REAL_TYPE), dimension(3), intent(out) :: DoubleDot2_4_result
    ! Local variables
    integer(INTEGER_TYPE) :: ii, jj

    ! Check vector size
    if (size(v1) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_4 requires vector of size(3)!"
        return
    endif

    ! Check matrix size
    if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot2_4 requires 3-by-3 matrix "
        return
    endif

    ! Compute the result
    do ii = 1, 3
        DoubleDot2_4_result(ii) = 0.0
        do jj = 1, 3
            DoubleDot2_4_result(ii) = DoubleDot2_4_result(ii) + ( m1(jj, ii) * v1(jj) )
        enddo
    enddo

    end subroutine DoubleDot2_4
    
    
    
    
    
    
    
    subroutine DoubleDot4_4(m1, m2, DoubleDot4_4_result)
    !// computes doubledot product for matrix-matrix arguments
    !// caution: second coordinate of the first matrix should be in opposite 
    !// variant form of the first coordinate of second matrix
    implicit none
    ! Arguments
    real(REAL_TYPE), dimension(3, 3), intent(in) :: m1, m2
    real(REAL_TYPE), dimension(3, 3), intent(out) :: DoubleDot4_4_result
    ! Local variables
    integer(INTEGER_TYPE) :: ii, jj, kk

    ! Check matrix sizes
    if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3 .or. &
        size(m2, 1) /= 3 .or. size(m2, 2) /= 3) then
        print *, "ERROR! PM4Sand::DoubleDot4_4 requires 3-by-3 matrices"
        return
    endif

    ! Compute the result
    do ii = 1, 3
        do jj = 1, 3
            DoubleDot4_4_result(ii, jj) = 0.0
            do kk = 1, 3
                DoubleDot4_4_result(ii, jj) = DoubleDot4_4_result(ii, jj) + ( m1(ii, kk) * m2(kk, jj) )
            enddo
        enddo
    enddo
    
    
    
    end subroutine DoubleDot4_4
    
    
    
    
    
    subroutine ToContraviant(v1, ToContraviant_result)
    !/*************************************************************/
    !// ToContraviant() ---------------------------------------------
    implicit none
    ! Arguments
    real(REAL_TYPE), dimension(3), intent(in) :: v1
    real(REAL_TYPE), dimension(3), intent(out) :: ToContraviant_result
    ! Local variables
    integer :: ii
    
    ! initialize
    ToContraviant_result = 0.0

    ! Check vector size
    if (size(v1) /= 3) then
        print *, "ERROR! PM4Sand::ToContraviant requires vector of size(3)!"
        return
    endif

    ! Copy input vector to result
    ToContraviant_result = v1

    ! Adjust second component of the result
    ToContraviant_result(3) = ToContraviant_result(3) * 0.5

    end subroutine ToContraviant
    
    
    
    
    
    subroutine ToCovariant(v1, ToCovariant_result)
    !/*************************************************************/
    !// ToCovariant() ---------------------------------------------
    implicit none
    ! Arguments
    real, dimension(3), intent(in) :: v1
    real, dimension(3), intent(out) :: ToCovariant_result
    ! Local variables
    integer :: ii

    ! Check vector size
    if (size(v1) /= 3) then
        print *, "ERROR! PM4Sand::ToCovariant requires vector of size(3)!"
        return
    endif

    ! Copy input vector to result
    ToCovariant_result = v1

    ! Adjust second component of the result
    ToCovariant_result(3) = ToCovariant_result(3) * 2.0


    end subroutine ToCovariant
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    end module 