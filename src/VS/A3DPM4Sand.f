
    
    
    
!      module PM4SandMaterialModule
!      !implicit none
!
!      contains
!    
!    
!    
!    !---------------------------------------------------------------------------------------
!        
!        
!!██████╗░███╗░░░███╗░░██╗██╗░██████╗░█████╗░███╗░░██╗██████╗░
!!██╔══██╗████╗░████║░██╔╝██║██╔════╝██╔══██╗████╗░██║██╔══██╗
!!██████╔╝██╔████╔██║██╔╝░██║╚█████╗░███████║██╔██╗██║██║░░██║
!!██╔═══╝░██║╚██╔╝██║███████║░╚═══██╗██╔══██║██║╚████║██║░░██║
!!██║░░░░░██║░╚═╝░██║╚════██║██████╔╝██║░░██║██║░╚███║██████╔╝
!!╚═╝░░░░░╚═╝░░░░░╚═╝░░░░░╚═╝╚═════╝░╚═╝░░╚═╝╚═╝░░╚══╝╚═════╝░
!        
!        
!        ! PM4SAND
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        
!        SUBROUTINE ESM_AM_PM4SAND(NPT,NOEL,IDSET,STRESS, &
!      EUNLOADING,PLASTICMULTIPLIER, &
!      DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,&
!                NPROPS,PROPS,NUMBEROFPHASES,NTENS)
!        
!      !!  implicit none 
!      !! 
!      !! ! input parameters 
!      !! character(len = 80) :: CMNAME ! constitutive model name
!      !! 
!      !! !integer(4) :: 4
!      !! !integer(4) :: 8
!      !! 
!      !! integer(4), intent(in) :: NPT ! point number
!      !! integer(4), intent(in) :: NOEL ! element number
!      !! integer(4), intent(in) :: IDSET ! ??
!      !! integer(4), intent(in) :: NTENS ! tensor order --> 6 in 3D, 4 in 2D
!      !! integer(4), intent(in) :: NUMBEROFPHASES ! number of phases
!      !! 
!      !! integer(4), intent(in) :: NSTATEV ! number of state variables
!      !! integer(4), intent(in) :: NADDVAR ! number of additional variables 
!      !! integer(4), intent(in) :: NPROPS ! number of properties 
!      !! 
!      !! ! output parameters 
!      !! 
!      !! ! inout parameters 
!      !! real(8), intent(inout) :: PLASTICMULTIPLIER ! scalar plastic multiplier 
!      !! 
!      !! real(8), dimension(NSTATEV), intent(inout) :: STATEV ! actual state variables
!      !! real(8), dimension(NADDVAR), intent(inout) :: &
!      !!ADDITIONALVAR ! actual additional variables 
!      !! real(8), dimension(NPROPS), intent(inout) :: PROPS ! actual properties 
!      !! 
!      !! real(8), dimension(NTENS), intent(inout) :: STRESS ! final stress tensor
!      !! real(8), dimension(NTENS), intent(inout) :: DSTRAN ! strain increment
!      !! 
!      !! real(8), dimension(NTENS,NTENS), intent(inout) :: &
!      !!EUNLOADING ! unloading elastic matrix
!      !! 
!      !!  !---Local variables required in standard UMAT
!      !!  integer(4) :: IStep, TimeStep
!      !!  
!      !!  real(8), dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
!      !!  real(8), dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
!      !!  real(8), dimension(:), allocatable :: stran ! strain 
!      !!  real(8), dimension(:), allocatable :: time ! time
!      !!  real(8), dimension(:), allocatable :: predef
!      !!  real(8), dimension(:), allocatable :: dpred    
!      !!  real(8), dimension(:), allocatable :: coords
!      !!  real(8), dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
!      !!  real(8), dimension(:,:), allocatable :: drot
!      !!  real(8), dimension(:,:), allocatable :: dfgrd0
!      !!  real(8), dimension(:,:), allocatable :: dfgrd1
!      !!  
!      !!  real(8) :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
!      !!  real(8) :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
!      !!  real(8) :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
!      !!  real(8) :: pnewdt, dtime, temp, dtemp, celent
!      !!  real(8) :: Value ! auxiliary variable holding any real valued number
!      !!  real(8) :: Porosity, WaterPressure, WaterPressure0, &
!      !!GasPressure, GasPressure0, DegreeSaturation  
!      !!
!      !!  integer(4) :: ndi, nshr, layer, kspt, kstep, kinc
!      !!  
!      !!  integer(4) :: IDTask
!      !!  
!      !!  allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2),&
!      !!predef(1), dpred(1),  &
!      !!coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), &
!      !!dfgrd1(3,3) )
!    
!              !4 = 4
!              !8 = 8
!              
!        
!        !implicit real(REAL_TYPE) (a-h, o-z) 
!        implicit double precision (a-h, o-z) 
!      integer :: NTENS, NSTATEV, NADDVAR, NPROPS, NPT, NOEL, IDSET, NUMBEROFPHASES
!      double precision :: EUNLOADING, PLASTICMULTIPLIER
!      CHARACTER*80 CMNAME         
!     ! DIMENSION NPT(1),NOEL(1),IDSET(1),STRESS(NTENS),EUNLOADING(1),PLASTICMULTIPLIER(1),&
!     !DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS),NUMBEROFPHASES(1)
!       DIMENSION STRESS(NTENS),DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
!
!!---Local variables required in standard UMAT
!        integer :: IStep, TimeStep
!        double precision, dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
!        double precision, dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
!        double precision, dimension(:), allocatable :: stran
!        double precision, dimension(:), allocatable :: time
!        double precision, dimension(:), allocatable :: predef
!        double precision, dimension(:), allocatable :: dpred    
!        double precision, dimension(:), allocatable :: coords
!        double precision, dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
!        double precision, dimension(:,:), allocatable :: drot
!        double precision, dimension(:,:), allocatable :: dfgrd0
!        double precision, dimension(:,:), allocatable :: dfgrd1
!        double precision :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
!        double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
!        double precision :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
!        double precision :: pnewdt, dtime, temp, dtemp, celent
!        double precision :: Value ! auxiliary variable holding any real valued number
!        double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation  
!        double precision :: VolumetricStress
!        integer :: ITens
!    
!        integer :: ndi, nshr, layer, kspt, kstep, kinc     
!
!        
!        
!        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
!              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
!              
!        ! Initialization
!        Eunloading = 0.0
!        PlasticMultiplier = 0.0
!        IDTask = 2
!     
!        !Rename additional variables --> 11 additional state variables
!        Porosity = AdditionalVar(1)
!        WaterPressure = AdditionalVar(2)
!        WaterPressure0 = AdditionalVar(3)
!        GasPressure = AdditionalVar(4)
!        GasPressure0 = AdditionalVar(5)
!        DegreeSaturation = AdditionalVar(6)
!        time(1) = AdditionalVar(7)   !TotalRealTime
!        time(2) = AdditionalVar(8)   !OverallTotalTime
!        dtime = AdditionalVar(9)     !TimeIncrement
!        IStep = AdditionalVar(10)    
!        TimeStep = AdditionalVar(11)   !Note: Very first time and load step: Istep=1 and TimeStep=1   
!        
!      
!        ! UMAT checked that it follows close format of the Um
!        ! set first call as 1... do we even need this if we have first call as a state variable
!        IF ((IStep==1).and.(TimeStep==1)) IDTask = 1 !--> why do I need this.... I guess for initialization 
!     
!        ! we need to always do this every time step
!      
!        !! 2D strain assignement 
!        !STATEV(44) = DSTRAN(1) ! mEpsilon_n
!        !STATEV(45) = DSTRAN(2)! mEpsilon_n
!        !STATEV(46) = DSTRAN(4)! mEpsilon_n
!        !
!        !! 2D plane stress assignement
!        !STATEV(50) = STRESS(1) ! SIGMA_XX ! mSigma_n
!        !STATEV(51) = STRESS(2) ! SIGMA_YY ! mSigma_n
!        !STATEV(52) = STRESS(4) ! SIGMA_XY ! mSigma_n
!        !
!        !! initially mSigma_n = mSigma
!        !STATEV(53) = STRESS(1) ! SIGMA_XX ! mSigma
!        !STATEV(54) = STRESS(2) ! SIGMA_YY ! mSigma
!        !STATEV(55) = STRESS(4) ! SIGMA_XY ! mSigma
!            
!        !---Call the UMAT
!        call UMAT_PM4Sand(stress, statev, ddsdde, sse, &!UMAT_AM_PM4SAND
!        spd, scd, rpl,&
!       ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
!       dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, &
!        nprops, coords, drot, pnewdt, celent, dfgrd0, &
!           dfgrd1, noel, npt, layer, kspt, kstep, kinc)
!
!
!      
!        !---Definition of Eunloading -> required to define the max time step
!        Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
!        
!        !---Always define this value to run the simulation
!
!        ! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
!        
!        return
!     
!       end subroutine ESM_AM_PM4SAND
!        
!        
!        
!        
!        
!        !*USER SUBROUTINES 
!      SUBROUTINE UMAT_PM4Sand(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,& !UMAT_AM_PM4SAND
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!      
!      implicit none
!  
!       ! input parameters 
!       integer(4), INTENT(IN) :: NPT, NTENS ! point number
!       integer(4), intent(in) :: NOEL ! element number
!     
!       real(8), intent(in) :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
!       real(8), intent(in) :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
!       real(8), intent(in) :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
!       real(8), intent(in) :: pnewdt, dtime, temp, dtemp, celent
!      
!       character(len = 80), intent(in) :: CMNAME ! constitutive model name
!        
!       real(8), dimension(NSTATEV), intent(inout) :: STATEV ! actual state variables 
!        
!       integer(4), intent(in) :: NSTATEV ! number of state variables 
!       
!       integer(4), intent(in) :: NPROPS ! number of properties 
!       real(8), dimension(NPROPS), intent(in) :: PROPS ! actual properties 
!      
!      real(8), intent(out), dimension(NTENS) :: STRESS ! size of this has been noted before
!      
!      real(8), intent(inout), dimension(NTENS,NTENS) :: DDSDDE
!      real(8), intent(in), dimension(NTENS) :: DDSDDT
!      real(8), intent(in), dimension(NTENS) :: DRPLDE
!      
!      real(8), intent(in), dimension(NTENS) :: STRAN
!      real(8), intent(in), dimension(NTENS) :: DSTRAN
!      
!      real(8), dimension(3) :: StrainIncrement
!      
!      real(8), intent(in), dimension(2) :: TIME
!      real(8), intent(in), dimension(1) :: PREDEF
!      real(8), intent(in), dimension(1) :: DPRED
!      
!      real(8), intent(in), dimension(3) :: COORDS
!      real(8), intent(in), dimension(3,3) :: DROT
!      real(8), intent(in), dimension(3,3) :: DFGRD0
!      real(8), intent(in), dimension(3,3) :: DFGRD1
!
!      integer(4), intent(in) :: ndi, nshr, layer, kspt, &
!      kstep, kinc     
!! Arguments:
!!          I/O  Type
!!  PROPS    I   R()  : List with model parameters
!!  DSTRAN   I   R()  : Strain increment
!!  DDSDDE   O   R(,) : Material stiffness matrix
!!  STRESS  I/O  R()  : stresses
!!  STATEV  I/O  R()  : state variables
!!
!!
!!---  Local variables
!    
!      real(8), dimension(3,3) :: DE_2D  !, intent(out)
!      real(8), dimension(3) :: dSig_2D !, intent(out)
!      real(8), dimension(3) :: Sig_2D !, intent(out)
!      real(8), dimension(3) :: dEpsP_2D !, intent(out)
!        
!      INTEGER(4) :: Iter, ii
!      
!      real(8), dimension(3) :: nn
!      real(8), dimension(3) :: RR
!      
!      
!      
!      
!      
!      
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
!      real(8), dimension(3) :: rrr
!      
!        ! PM4SAND VARIABLES 
!        real(8) :: Dr, G0, hpo
!        
!        real(8) :: P_atm, h0, emax, emin, nb, nd, Ado, &
!      z_max, cz, ceps, &
!      Mc, nu, C_GD, C_DR, Ckaf, QQ_Bolton, RR_Bolton, &
!      mm, Fsed_min, p_sedo, &
!      TolF, TolR, Pmin, Pmin2, & !FirstCall, PostShake, TanType, IntSchem
!      DGamma, DGamma_n, KK, GG, InitialVoidRatio, Kp, zcum, & !pzpFlag, me2p
!      zpeak, pzp, zxp, Mb, Md, Mcur, phi_cv, Cgd, Cdr, &
!      K_p, NextVoidRatio   
!        
!        
!        real(8), dimension(3) :: dEpsilonE, &
!       Epsilon, Epsilon_n, Sigma, Sigma_n, &
!      Sigma_b, EpsilonE, EpsilonE_n, Alpha, Alpha_n, Alpha_in_n, &  
!      Alpha_in_p_n, Alpha_in_true_n, Alpha_in_max_n, Alpha_in_min_n, & 
!      Fabric, Fabric_n, Fabric_in_n   
!        
!        ! these are not state variables but we need them 
!        real(8), dimension(3) :: Alpha_in,  &
!      Alpha_in_p, Alpha_in_true, Alpha_in_max, Alpha_in_min, Fabric_in
!        
!        
!        real(8), dimension(3,3) :: Ce
!        
!   
!        logical :: me2p, pzpFlag, FirstCall, PostShake
!   
!      Iter = 0
!     
!      
!      ! I commented this for now 5/17/2024
!      ! Geotechnical Engineering convention for 2D planestrain problem
!      StrainIncrement = 0.0
!      StrainIncrement(1) = DSTRAN(1) * (-1.0) ! geotechnical engineering sign convention (compression positive)
!      StrainIncrement(2) = DSTRAN(2) * (-1.0)
!      StrainIncrement(3) = DSTRAN(4) ! Epsilon = 0.5 * Gamma  !0.5*
!
!      
!      do ii = 1, size(DSTRAN)
!      if (DSTRAN(ii) /= DSTRAN(ii)) then
!        print *, "DSTRAN(", ii, ") is NaN!"
!      end if
!      end do
!      
!      !Dr
!      !G0
!      !PostShake
!      !me2p   
!      !P_atm 
!      !h0 
!      !emax 
!      !emin
!      !nb 
!      !nd 
!      !Ad, 
!      !z_max 
!      !cz 
!      !ceps
!      !phi_cv 
!      !nu
!      !Cgd
!      !Ckaf
!      !QQ_Bolton
!      !RR_Bolton
!      !mm
!      !VoidRatio
!      !Mc
!      !Mb
!      !Md
!      !Cdr
!      !Fsed_min
!      !p_sedo
!      !Pmin
!      !Pmin2
!      !K_p
!      !zpeak
!      !zcum
!      !pzp 
!      !zxp  ! z_max, Ado, 
!      !KK,
!      !GG
!      !Mcur 
!      !rrr
!      !Sigma_b
!      !Fabric_n
!      !Fabric, 
!      !Fabric_in_n, 
!      !Fabric_in, 
!      !Alpha_n, 
!      !Alpha, 
!      !Alpha_in_n, 
!      !Alpha_in, 
!      !Alpha_in_p_n, 
!      !Alpha_in_p, 
!      !Alpha_in_true_n, 
!      !Alpha_in_true, 
!      !Alpha_in_max_n, 
!      !Alpha_in_max, 
!      !Alpha_in_min_n,
!      !Alpha_in_min, 
!      !pzpFlag, 
!      !TolF, 
!      !TolR, 
!      !dEpsilonE, 
!      !Sigma_n, 
!      !Sigma, 
!      !FirstCall
!      !Ce
!      
!      ! read properties (5 items)
!      Dr = PROPS(1) ! relative density 
!      G0 = PROPS(2) ! small-strain shear stiffness
!      hpo = PROPS(3) ! contraction rate parameter
!      PostShake = PROPS(4) ! postshake switch
!      me2p = PROPS(5) ! switch to invoke plasticity
!      
!      
!      
!      ! read state variables (114 items)
!      P_atm = STATEV(1)
!      h0 = STATEV(2)
!      emax = STATEV(3)
!      emin = STATEV(4)
!      nb = STATEV(5)
!      nd = STATEV(6)
!      Ado = STATEV(7)
!      
!      z_max = STATEV(8)
!      cz =STATEV(9)
!      ceps = STATEV(10)
!      phi_cv = STATEV(11)
!      nu = STATEV(12)
!      
!      Cgd = STATEV(13)
!      Ckaf = STATEV(14)
!      QQ_Bolton = STATEV(15)
!      RR_Bolton = STATEV(16)
!      
!      mm = STATEV(17)
!      InitialVoidRatio = STATEV(18)
!      
!      
!      Mc = STATEV(19)
!      Mb = STATEV(20)
!      Md = STATEV(21)
!      
!      Cdr = STATEV(22) ! Cdr<5.0
!      
!      Fsed_min = STATEV(23)
!      p_sedo = STATEV(24)
!      Pmin = STATEV(25)
!      Pmin2 = STATEV(26)
!      K_p = STATEV(27)
!      zpeak = STATEV(28)
!      zcum = STATEV(29)
!      
!      
!      pzp =STATEV(30)
!      zxp = STATEV(31)
!      
!      KK = STATEV(32) 
!      GG = STATEV(33)
!      
!      Mcur  = STATEV(34)
!      rrr = STATEV(35:37)
!      Sigma_b = STATEV(38:40)
!      
!      Fabric_n  = STATEV(41:43)
!      Fabric  = STATEV(44:46)
!      Fabric_in_n  = STATEV(47:49)
!      Fabric_in  = STATEV(50:52)
!      
!      Alpha_n = STATEV(53:55)
!      Alpha = STATEV(56:58)
!      Alpha_in_n = STATEV(59:61)
!      Alpha_in  = STATEV(62:64)
!      Alpha_in_p_n  = STATEV(65:67)
!      
!      
!      Alpha_in_p = STATEV(68:70)
!      Alpha_in_true_n  = STATEV(71:73)
!      Alpha_in_true = STATEV(74:76)
!      Alpha_in_max_n  = STATEV(77:79)
!      Alpha_in_max = STATEV(80:82)
!      
!      Alpha_in_min_n  = STATEV(83:85)
!      Alpha_in_min  = STATEV(86:88)
!      
!      pzpFlag  = STATEV(89)
!      TolF = STATEV(90)
!      TolR = STATEV(91)
!      
!      dEpsilonE =  STATEV(92:94) 
!      
!      Sigma_n  = STATEV(95:97)
!      Sigma  = STATEV(98:100)
!      
!      FirstCall = STATEV(101)
!      
!      Epsilon = STATEV(102:104)
!      Epsilon_n = STATEV(105:107)
!      
!      EpsilonE = STATEV(108:110)
!      EpsilonE_n = STATEV(111:113)
!      
!      ! I added this so that we can always update this
!      NextVoidRatio = STATEV(114)
!      
!      Ce = DE_2D
!      
!      
!      if (FirstCall) then 
!          Sigma_n(1) = stress(1) * (-1) 
!          Sigma_n(2) = stress(2) * (-1) 
!          Sigma_n(3) = stress(4) 
!          
!      end if 
!      
!      
!      
!      call PM4SandFullConstructor(Dr, G0, PostShake, me2p, &    
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, phi_cv, nu, &
!      Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, &
!      pzp, zxp, & ! z_max, Ado, 
!      KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, &
!      Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, &
!      Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, &
!      Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, & 
!      dEpsilonE, Sigma_n, Sigma, &
!      FirstCall, Ce, &
!      Epsilon, Epsilon_n, &
!      EpsilonE, EpsilonE_n, & 
!      NextVoidRatio )
!      
!      
!      
!      Epsilon = Epsilon + StrainIncrement     
!      
!      
!      
!      call PM4SandIntegrate(Epsilon, Epsilon_n, &
!      EpsilonE, EpsilonE_n, &
!      G0, hpo, PostShake, me2p, & 
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, &
!      ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, &
!      zcum, pzp, zxp, & 
!      KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, &
!      Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, & 
!      Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, &
!      Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
!      dEpsilonE, Sigma_n, Sigma, &
!      Ce, &
!      nn, RR, & 
!      NextVoidRatio)
!      
!      
!      
!      
!      
!      call PM4SandCommitState(StrainIncrement, &
!      z_max, nu, G0, mm, Md, Mb, PostShake, Pmin, P_atm, Cgd, p_sedo, &
!      Fsed_min, me2p, Fabric, Fabric_in, & ! input VoidRatio, 
!      Alpha_in, Alpha_in_p, Alpha_in_true, Alpha_in_max, &
!      Alpha_in_min, & ! input
!      KK, GG, Mcur, rrr, Alpha_in_n, Alpha_n, Alpha_in_p_n, &
!      Alpha_in_true_n, Alpha_in_max_n, Alpha_in_min_n, Sigma_n, & !dFabric, & !NextVoidRatio, & ! output
!      zcum, zpeak, Sigma, Alpha, Fabric_n, Fabric_in_n, &
!      Ce, &
!      Mc , &
!      nn, RR, &
!      K_p, &
!      Epsilon, Epsilon_n, &
!      EpsilonE, EpsilonE_n) ! inout
!    
!      
!   
!      
!      !------------------------------------------------------------------
!      
!      
!      
!      
!      
!      ! update global STATEV variable 
!      STATEV(1) = P_atm
!      STATEV(2) = h0
!      STATEV(3) = emax
!      STATEV(4) = emin
!      STATEV(5) = nb
!      STATEV(6) = nd 
!      STATEV(7) = Ado 
!      
!      STATEV(8) = z_max
!      STATEV(9) = cz
!      STATEV(10) = ceps 
!      STATEV(11) = phi_cv 
!      STATEV(12) = nu
!      
!      STATEV(13) = Cgd 
!      STATEV(14) = Ckaf 
!      STATEV(15) = QQ_Bolton 
!      STATEV(16) = RR_Bolton 
!      
!      STATEV(17) = mm
!      STATEV(18) = InitialVoidRatio
!      
!      
!      STATEV(19) = Mc
!      STATEV(20) = Mb 
!      STATEV(21) = Md 
!      
!      STATEV(22) = Cdr
!      STATEV(23) = Fsed_min
!      STATEV(24) = p_sedo
!      STATEV(25) = Pmin
!      STATEV(26) = Pmin2
!      STATEV(27) = K_p
!      STATEV(28) = zpeak
!      STATEV(29) = zcum
!      
!      
!      STATEV(30) = pzp
!      STATEV(31) = zxp 
!      
!      STATEV(32) = KK 
!      STATEV(33) = GG 
!      
!      STATEV(34) = Mcur 
!      STATEV(35:37) = rrr 
!      STATEV(38:40) = Sigma_b
!      
!      STATEV(41:43) = Fabric_n 
!      STATEV(44:46) = Fabric 
!      STATEV(47:49) = Fabric_in_n 
!      STATEV(50:52) = Fabric_in 
!      
!      STATEV(53:55) = Alpha_n 
!      STATEV(56:58) = Alpha 
!      STATEV(59:61) = Alpha_in_n 
!      STATEV(62:64) = Alpha_in 
!      STATEV(65:67) = Alpha_in_p_n 
!      
!      
!      
!      STATEV(68:70) = Alpha_in_p
!      STATEV(71:73) = Alpha_in_true_n 
!      STATEV(74:76) = Alpha_in_true 
!      STATEV(77:79) = Alpha_in_max_n 
!      STATEV(80:82) = Alpha_in_max 
!      
!      STATEV(83:85) = Alpha_in_min_n
!      STATEV(86:88) = Alpha_in_min
!      
!      STATEV(89) = pzpFlag
!      STATEV(90) = TolF
!      STATEV(91) = TolR
!      
!      STATEV(92:94) = dEpsilonE
!      
!      STATEV(95:97) = Sigma_n
!      STATEV(98:100) = Sigma
!      
!      STATEV(101) = FirstCall
!      
!      
!      
!      
!      ! these were added later on
!      STATEV(102:104) = Epsilon
!      STATEV(105:107) = Epsilon_n 
!      
!      STATEV(108:110) = EpsilonE 
!      STATEV(111:113) = EpsilonE_n
!      
!      STATEV(114) = NextVoidRatio
!      
!      DE_2D = Ce
!      
!      
!      !------------------------------------------------------------------
!      ! change the sign convention based on the tension +ve convention
!      STRESS(1) = STATEV(98) * (-1.0)
!      STRESS(2) = STATEV(99) * (-1.0)
!      STRESS(3) = 0.0!STATEV(98)
!      STRESS(4) = STATEV(100)
!      STRESS(5) = 0.0
!      STRESS(6) = 0.0
!      
!      
!      !STRESS(1) = 100.0 * (-1.0)
!      !STRESS(2) = 100.0 * (-1.0)
!      !STRESS(3) = 0.0!STATEV(98)
!      !STRESS(4) = STATEV(100)
!      !STRESS(5) = 0.0
!      !STRESS(6) = 0.0
!      
!      
!      
!      
!      DDSDDE(1,1) = DE_2D(1,1)
!      DDSDDE(2,1) = DE_2D(2,1)
!      DDSDDE(3,1) = DE_2D(2,1)
!      DDSDDE(4,1) = 0.0
!      DDSDDE(5,1) = 0.0
!      DDSDDE(6,1) = 0.0
!      
!      DDSDDE(1,2) = DE_2D(1,2)
!      DDSDDE(2,2) = DE_2D(2,2)
!      DDSDDE(3,2) = DE_2D(1,2)
!      DDSDDE(4,2) = 0.0
!      DDSDDE(5,2) = 0.0
!      DDSDDE(6,2) = 0.0
!      
!      DDSDDE(1,3) = DE_2D(2,1)
!      DDSDDE(2,3) = DE_2D(1,2)
!      DDSDDE(3,3) = DE_2D(2,2)
!      DDSDDE(4,3) = 0.0
!      DDSDDE(5,3) = 0.0
!      DDSDDE(6,3) = 0.0
!      
!      DDSDDE(1,4) = 0.0
!      DDSDDE(2,4) = 0.0
!      DDSDDE(3,4) = 0.0
!      DDSDDE(4,4) = DE_2D(3,3)
!      DDSDDE(5,4) = 0.0
!      DDSDDE(6,4) = 0.0
!      
!      DDSDDE(1,5) = 0.0
!      DDSDDE(2,5) = 0.0
!      DDSDDE(3,5) = 0.0
!      DDSDDE(4,5) = 0.0
!      DDSDDE(5,5) = DE_2D(3,3)
!      DDSDDE(6,5) = 0.0
!      
!      DDSDDE(1,6) = 0.0
!      DDSDDE(2,6) = 0.0
!      DDSDDE(3,6) = 0.0
!      DDSDDE(4,6) = 0.0
!      DDSDDE(5,6) = 0.0
!      DDSDDE(6,6) = DE_2D(3,3)
!      
!      
!      End subroutine UMAT_PM4Sand
!    
!    
!    
!    !---------------------------------------------------------------------------------------
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
!      subroutine PM4SandMaterial(oData) ! input the properties and the state variables  !pData, 
!    !----------------------------------------------
!    ! To have the variables that are needed as constants 
!    ! in the PM4Sand implementation. 
!    !----------------------------------------------
!    
!      implicit none 
!    
!    !real(8), dimension(3) :: pData ! this contains the properties of PM4Sand 
!    ! -pData(1) = D_r --> apparent relative density (variable controlling dilatancy and )
!    ! -pData(2) = G_0 --> shear modulus coefficient (small-strain shear modulus which can be calculated from shear wave velocity)
!    ! -pData(3) = h_p0 --> contraction rate parameter (variable adjusts contraction rates which can be adjusted to get desiered CSR)
!    
!      real(8), dimension(112), intent(inout) :: oData ! this contains the initial state variables (no fabric parameters) !107
!    !real(8), dimension(3), intent(in) :: pData ! this contains the initial state variables (no fabric parameters)
!    !logical :: FirstCall
!    
!    ! oData has 24 state variables that contains all the information about the state variables 
!    
!    ! I am not sure if we need this subroutine --> we do... BUT ONLY THE FIRST TIME WHEN WE NEED TO INITIALIZE THOSE PARAMETERS
!    ! we need the primary inputs here stored in some sort of properties parameters
!    !pData(1) = Dr
!    !pData(2) = G0
!    !pData(3) = hp0
!    
!    !m_e_init
!    !m_Mc
!    
!    
!      if (oData(25) == 1) then  ! if first call is yes 
!        
!    !-------------------------------------------if first time to initialize
!    ! initialization of list of parameters
!    ! these are standard parameters    
!    
!          oData(1) = 101.3    !// P_atm  in kPa ! Numbering according to the manual
!	
!    
!          oData(2) = -1       !// h0          ! parameter 1: -Adjusts the ratio of plastic modulus to elastic modulus. 
!                                        !              -Default value is h0 = (0.25+Dr)/2
!                                        !              -Minimum value of 0.30
!                                        !              -Provides reasonable G/Gmax and damping rationships for the 
!                                        !               default value of G0.
!                                        !              -May require adjustment in combination with any adjustments to G0.        
!	
!    
!          oData(3) = 0.8      !// emax        ! parameter 2: -Maximum void ratio that affects the computation of density.
!                                        !              -Default value is 0.8. 
!	
!    
!          oData(4) = 0.5      !// emin        ! parameter 2: -Minimum void ratio that affects the computation of density.
!                                        !              -Default value is 0.5. 
!	
!    
!          oData(5) = 0.5      !// nb          ! parameter 3: -Controls dilatancy and thus also the peak friction angles. 
!                                        !              -Default value is 0.5. --> Dense of crit
!                                        !              -Default value is 0.5/4. --> Loose of crit      
!	
!    
!          oData(6) = 0.1      !// nd          ! parameter 4: -Controls the stress-ratio at which contraction changes to dilation,
!                                        !               which is often referred to as phase transformation. 
!                                        !              -Produces a phase transformation angle slightly smaller than phi_cv.  
!                                        !              -Default value is 0.1. --> Dense of crit
!                                        !              -Default value is 4*0.1. --> Loose of crit
!	
!    
!          oData(7) = -1       !// Ado         ! parameter 5: -Default value is based on Bolton's dilatancy relationship at times of 
!                                        !               initialization. 
!                                        !              -Typical values will be between 1.2 and 1.5.
!    
!	
!          oData(8) = -1       !// z_max       ! parameter 6: -Default value is computed at the time of initialization:
!                                        !               z_max = 0.70 exp(-6.1 ksi_R0) <=20
!                                        !              -Maximum value is 20.
!                                        !              -if ksi_R is initially 0.0, z_max = 0.7.
!                                        !              -z_max increases when ksi_R increases (dense critical states) to reach a 
!                                        !               maximum value of 20.
!                                        !              -May require varying if relationship between Dr and cyclic strength is significantly
!                                        !               different from that implied by liquefaction correlations of Idriss and Boulanger (2008).
!    
!	
!          
!          
!          oData(9) = 250.0    !// cz          ! parameter 7: -Default value is 250. 
!                                        !              -Controls strain levels at which fabric effects become important. 
!    
!	
!          oData(10) = -1       !// ceps         ! parameter 8: -Default value varies with Dr. 
!                                        !              -Value is 5.0 for Dr<35%. --> Loose
!                                        !              -Value linearly decreases to its minimum value of 1.0 at Dr=75% --> Dense
!	
!    
!          oData(11) = 2 * sin(33.0 * (3.14159265359/180.0))    
!                                        !// m_Mc      ! parameter 9: -Default value is for phi_cv of 33 degrees. 
!    
!	
!          oData(12) = 0.3     !// nu          ! parameter 10: -Default value is 0.30. 
!                                        !               -k0 would correspond to k0=nu/(1-nu) --> k0 value = 0.43
!    
!	
!          oData(13) = 2.0     !// Cgd         ! parameter 11: -Default value is 2.0.
!                                        !               -G0 degrades with increasing plastic deviatoric strains (z_cum). 
!                                        !               -Maximum degredation approaches a factor of 1/Cgd
!    
!	
!          oData(14) = -1      !// C_DR        ! parameter 
!	
!    
!          oData(15) = -1      !// Ckaf        ! parameter 12: -Default value varies with Dr, as:
!                                        !               Ckaf = 5 + (220*(Dro - 0.26)^3)
!                                        !               -Ckaf = 4.0 for Dr<10% increases to its maximum value of 
!                                        !                35.0 at Dr=77%
!                                        !               -Controls the effect that sustained static shear stresses have on 
!                                        !                plastic modulus. 
!	
!    
!          oData(16) = 10.0    !// Q           ! parameter 13: -Default vlaue is 10.0 for quartzitic sand per Bolton (1986)
!	
!    
!          oData(17) = 1.5     !// R           ! parameter 14: -Default value is 1.5. 
!                                        !               -This is a slight increase from Bolton (1986) to lower the CSL 
!                                        !                to better match results for direct simple shear loading.
!	
!    
!          oData(18) = 0.01    !// m           ! parameter 15: -Default value is 0.01.
!                                        !               -Provides reasonable modeling and numerical stability. 
!	
!    
!          oData(19) = -1       !// Fsed_min   !  post liq stuff
!	
!          oData(20) = -1       !//p_sdeo      ! post liq stuff
!	
!          oData(21) = 5		 !// IntScheme  ! this should be a choice 
!	
!          oData(22) = 0		 !// TanType    ! ??
!	
!          
!          oData(23) = 1.0e-8	 !// TolF       ! error threshold (yield surface)
!	
!          oData(24) = 1.0e-8	 !// TolR       ! error threshold --> do we even need this??
!    
!          oData(25) = 0        !// First call ! --> set to zero after first call be been done to overwrite the state variables
!          oData(26) = 0        !// Post shake switch
!    
!          oData(27) = oData(1) / 200.0      !m_Pmin !//
!          oData(28) = oData(26) * 5.0       !m_Pmin2 !//
!    
!          oData(29) = -1           !// mpzpFlag
!    
!          oData(30) = -1           !// me2p
!    
!          oData(31) = -1           !// mDGamma !-->
!          oData(32) = -1           !// mDGamma_n
!    
!          oData(33) = -1           !// mK
!          oData(34) = -1           !// mG
!    
!          oData(35) = -1           !// mVoidRatio
!    
!          oData(36) = -1           !// mKp
!    
!          oData(37) = -1           !// mzcum
!          oData(38) = -1           !// mzpeak
!          oData(39) = -1           !// mpzp
!          oData(40) = -1           !// mzxp
!          !oData() = -1        --> zxpPk?    
!    
!          oData(41) = -1           !// mMb --> bounding surface stress ratio
!          oData(42) = -1           !// mMd --> dilatancy surface stress ratio
!          oData(43) = -1           !// mMcur --> current stress ratio
!    
!    
!    ! replace these mEpsilon_n and mEpsilon with dStran
!    !oData(44) = -1           !// mEpsilon_n(1) --> old 
!    !oData(45) = -1           !// mEpsilon_n(2)
!    !oData(46) = -1           !// mEpsilon_n(3)
!    
!    !oData(47) = -1           !// mEpsilon(1) --> new
!    !oData(48) = -1           !// mEpsilon(2)
!    !oData(49) = -1           !// mEpsilon(3)
!    
!    !oData(50) = -1           !// mSigma_n(1) --> old
!    !oData(51) = -1           !// mSigma_n(2)
!    !oData(52) = -1           !// mSigma_n(3)
!    !
!    !oData(53) = -1           !// mSigma(1) --> new
!    !oData(54) = -1           !// mSigma(2)
!    !oData(55) = -1           !// mSigma(3)
!    
!          oData(56) = -1           !// mSigma_b(1) ----> difference
!          oData(57) = -1           !// mSigma_b(2)
!          oData(58) = -1           !// mSigma_b(3)
!    
!          oData(59) = -1           !// mEpsilonE_n(1) --> old
!          oData(60) = -1           !// mEpsilonE_n(2)
!          oData(61) = -1           !// mEpsilonE_n(3)
!    
!          oData(62) = -1           !// mEpsilonE(1) --> new
!          oData(63) = -1           !// mEpsilonE(2)
!          oData(64) = -1           !// mEpsilonE(3)
!    
!          oData(65) = -1           !//mAlpha_n(1) --> old
!          oData(66) = -1           !//mAlpha_n(2)
!          oData(67) = -1           !//mAlpha_n(3)
!    
!          oData(68) = -1           !//mAlpha(1) --> new        
!          oData(69) = -1           !//mAlpha(2)         
!          oData(70) = -1           !//mAlpha(3)        
!    
!          oData(71) = -1           !//mAlpha_in_n(1) --> old 
!          oData(72) = -1           !//mAlpha_in_n(2)
!          oData(73) = -1           !//mAlpha_in_n(3)
!    
!          oData(74) = -1           !//mAlpha_in(1) --> new
!          oData(75) = -1           !//mAlpha_in(2)
!          oData(76) = -1           !//mAlpha_in(3)
!    
!          oData(77) = -1           !//mAlpha_in_p_n(1) --> old
!          oData(78) = -1           !//mAlpha_in_p_n(2)
!          oData(79) = -1           !//mAlpha_in_p_n(3)
!    
!          oData(80) = -1           !// mAlpha_in_p(1) --> new
!          oData(81) = -1            !// mAlpha_in_p(2)
!          oData(82) = -1            !// mAlpha_in_p(3)
!    
!          oData(83) = -1           !//mAlpha_in_true_n(1) --> old
!          oData(84) = -1           !//mAlpha_in_true_n(2)
!          oData(85) = -1           !//mAlpha_in_true_n(3)
!    
!          oData(86) = -1            !// mAlpha_in_true(1) --> new
!          oData(87) = -1            !// mAlpha_in_true(2)
!          oData(88) = -1            !// mAlpha_in_true(3)
!    
!          oData(89) = -1           !//mAlpha_in_max_n(1) --> old
!          oData(90) = -1           !//mAlpha_in_max_n(2)
!          oData(91) = -1           !//mAlpha_in_max_n(3)
!    
!          oData(92) = -1            !// mAlpha_in_max(1) --> new
!          oData(93) = -1            !// mAlpha_in_max(2)
!          oData(94) = -1            !// mAlpha_in_max(3)
!    
!          oData(95) = -1           !//mAlpha_in_min_n(1) --> old
!          oData(96) = -1           !//mAlpha_in_min_n(2)
!          oData(97) = -1           !//mAlpha_in_min_n(3)
!    
!          oData(98) = -1            !// mAlpha_in_min(1) --> new
!          oData(99) = -1            !// mAlpha_in_min(2)
!          oData(100) = -1            !// mAlpha_in_min(3)
!    
!          oData(101) = -1           !//mFabric_n(1)      --> old
!          oData(102) = -1           !//mFabric_n(2)      
!          oData(103) = -1           !//mFabric_n(3)
!    
!          oData(104) = -1           !//mFabric(1)        --> new
!          oData(105) = -1           !//mFabric(2)        
!          oData(106) = -1           !//mFabric(3)       
!    
!          oData(107) = -1           !//mFabric_in_n(1) --> old
!          oData(108) = -1           !//mFabric_in_n(2)
!          oData(109) = -1           !//mFabric_in_n(3)
!    
!          oData(110) = -1            !// mFabric_in(1) --> new
!          oData(111) = -1            !// mFabric_in(2)
!          oData(112) = -1            !// mFabric_in(3)
!    
!    
!      else 
!        
!        !just keep the original oData from the MP
!    
!    
!      end if
!    
!    
!    
!    
!      end subroutine PM4SandMaterial
!    
!    
!    
!    
!    
!    
!    
!    
!      subroutine PM4SandFullConstructor(Dr, G0, PostShake, me2p, &    
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, phi_cv, nu, &
!      Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum, &
!      pzp, zxp, & ! z_max, Ado, 
!      KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, &
!      Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, &
!      Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, &
!      Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
!      dEpsilonE, Sigma_n, Sigma, &
!      FirstCall, Ce, &
!      Epsilon, Epsilon_n, &
!      EpsilonE, EpsilonE_n, &
!      NextVoidRatio)
!    !----------------------------------------------
!    ! To populate the stress-strain variables needed 
!    ! 
!    !---------------------------------------------- 
!    
!      implicit none 
!    
!      ! I would say we call this subroutine first to initialize the properties
!      ! input variables
!    
!      real(8), intent(in) :: Dr       ! primiary parameter 1 ! prop 1
!    
!      real(8), intent(in) :: G0       ! primiary parameter 2 ! prop 2
!    
!    
!      ! prop 3 --> hpo
!    
!    
!      logical, intent(in) :: PostShake ! prop 4 logical, intent(in) ::  !real(8)
!    
!      logical, intent(in) :: me2p ! prop 5 logical, intent(in) :: 
!    
!    
!      ! output variables 
!    
!      real(8), intent(inout) :: P_atm    ! statev 1
!    
!      real(8), intent(inout) :: h0       ! statev 2              ! secondary parameter 1
!    
!      real(8), intent(inout) :: emax     ! statev 3                       ! secondary parameter 2
!      real(8), intent(inout) :: emin     ! statev 4                        ! secondary parameter 2
!      real(8), intent(inout) :: nb       ! statev 5                        ! secondary parameter 3
!      real(8), intent(inout) :: nd       ! statev 6                        ! secondary parameter 4
!      real(8), intent(inout) :: Ado      ! statev 7                        ! secondary parameter 5
!      real(8), intent(inout) :: z_max    ! statev 8                        ! secondary parameter 6
!      real(8), intent(inout) :: cz       ! statev 9                        ! secondary parameter 7
!      real(8), intent(inout) :: ceps     ! statev 10                          ! secondary parameter 8
!      real(8), intent(inout) :: phi_cv   ! statev 11                        ! secondary parameter 9
!      real(8), intent(inout) :: nu       ! statev 12                        ! secondary parameter 10
!      real(8), intent(inout) :: Cgd      ! statev 13                        ! secondary parameter 11
!
!      real(8), intent(inout) :: Ckaf     ! statev 14                        ! secondary parameter 12
!
!      real(8), intent(inout) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
!      real(8), intent(inout) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
!    
!      real(8), intent(inout) :: mm         ! statev 17                             ! secondary parameter 15
!
!      real(8), intent(inout) :: InitialVoidRatio   ! statev 18
!
!      real(8), intent(out) :: NextVoidRatio   ! statev 18
!      
!      
!      real(8), intent(inout) :: Mc         ! statev 19                      ! secondary parameter 1
!      real(8), intent(inout) :: Mb         ! statev 20 
!      real(8), intent(inout) :: Md         ! statev 21
!    
!      real(8), intent(inout) :: Cdr        ! statev 22                       ! secondary parameter 
!    
!      real(8), intent(inout) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
!      real(8), intent(inout) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
!    
!      real(8), intent(inout) :: Pmin       ! statev 25 
!      real(8), intent(inout) :: Pmin2      ! statev 26
!
!      real(8), intent(inout) :: K_p        ! statev 27
!    
!      real(8), intent(inout) :: zpeak      ! statev 28
!      real(8), intent(inout) :: zcum       ! statev 29
!    
!      real(8), intent(inout) :: pzp        ! statev 30 
!      real(8), intent(inout) :: zxp        ! statev 31
!    
!      real(8), intent(inout) :: KK         ! statev 32
!      real(8), intent(inout) :: GG         ! statev 33  
!      real(8), intent(inout) :: Mcur       ! statev 34
!    
!      real(8), intent(inout), dimension(3) :: rrr   ! statev 35:37
!      real(8), intent(inout), dimension(3) :: Sigma_b ! statev 38:40
!    
!      real(8), intent(inout), dimension(3) :: Fabric_n ! statev 41:43
!      real(8), intent(inout), dimension(3) :: Fabric ! statev 44:46
!    
!      real(8), intent(inout), dimension(3) :: Fabric_in_n ! statev 47:49
!      real(8), intent(inout), dimension(3) :: Fabric_in ! statev 50:52
!    
!      real(8), intent(inout), dimension(3) :: Alpha_n ! statev 53:55
!      real(8), intent(inout), dimension(3) :: Alpha ! statev 56:58
!    
!      real(8), intent(inout), dimension(3) :: Alpha_in_n ! statev 59:61
!      real(8), intent(inout), dimension(3) :: Alpha_in ! statev 62:64
!    
!      real(8), intent(inout), dimension(3) :: Alpha_in_p_n ! statev 65:67
!      real(8), intent(inout), dimension(3) :: Alpha_in_p ! statev 68:70
!    
!      real(8), intent(inout), dimension(3) :: Alpha_in_true_n ! statev 71:73
!      real(8), intent(inout), dimension(3) :: Alpha_in_true ! statev 74:76
!    
!      real(8), intent(inout), dimension(3) :: Alpha_in_max_n ! statev 77:79
!      real(8), intent(inout), dimension(3) :: Alpha_in_max ! statev 80:82
!    
!      real(8), intent(inout), dimension(3) :: Alpha_in_min_n ! statev 83:85
!      real(8), intent(inout), dimension(3) :: Alpha_in_min ! statev 86:88
!    
!      logical, intent(inout) :: pzpFlag ! statev 89 ! statev has to be same type throughout !real(8)
!    
!      real(8), intent(inout) :: TolF ! statev 90
!	real(8), intent(inout) :: TolR ! statev 91 
!    
!      real(8), intent(inout), dimension(3,3) :: Ce
!    
!      real(8), intent(inout), dimension(3) :: dEpsilonE ! statev 92:94
!    
!      real(8), intent(inout), dimension(3) :: Sigma_n ! statev 95:97                 ! current strain      
!      real(8), intent(inout), dimension(3) :: Sigma ! statev 98:100
!    
!      ! inout variables
!      logical, intent(inout) :: FirstCall ! statev 101   !real(8)              ! current strain      
!    
!      ! local variables 
!      real(8) :: h0_1
!      real(8) :: h0_2
!    
!      real(8) :: ksi
!      real(8) :: Mcut
!      real(8) :: Mfin
!    
!      real(8) :: p0
!    
!      real(8), dimension(3) :: initStress 
!      real(8), dimension(3) :: I1  
!    
!      !real(8), dimension(3) :: GetDevPart_result
!      
!      
!      real(8), dimension(3), intent(inout) :: Epsilon
!      real(8), dimension(3), intent(inout) :: Epsilon_n 
!      
!      real(8), dimension(3), intent(inout) :: EpsilonE
!      real(8), dimension(3), intent(inout) :: EpsilonE_n 
!    
!      
!    
!      ! initialize identity matrix 
!      I1(1) = 1 
!      I1(2) = 1
!      I1(3) = 0
!    
!      ! check if we need to reinitialize state variables
!      if (FirstCall == .TRUE.) then 
!          
!         ! I need to double check if I am initializing from a linear elastic constitutive model
!         !Sigma = 0
!      
!        
!        ! initialize atmospheric pressure
!        if (P_atm <= 0) then 
!            P_atm = 101.3 !1.01325!100.0!101.325 !kPa
!        end if 
!        
!        ! initialize parameter 1
!        if (h0 < 0) then
!            h0_1 = 0.3
!            h0_2 = (0.25+Dr)/2.0
!            h0 = max(h0_1, h0_2)
!        end if 
!    
!        if (emax < 0) then
!            emax = 0.8!0.739 !0.8!0.7389!0.8
!        end if 
!    
!        if (emin < 0) then
!            emin = 0.5!0.492!0.511 !0.4915!0.5
!        end if
!    
!        if (nb < 0) then
!            nb = 0.50!0.7!0.5!0.7!0.5!0.7!0.5
!        end if
!    
!        if (nd < 0) then
!            nd = 0.1
!        end if 
!
!        if (cz < 0) then
!            !cz = 200!250!200!250!200!250!200!250
!            
!            !if (Dr >= 0.75) then 
!            !    cz = 125!0.225!0.2!1.0!0.2
!            !else if (Dr <= 0.35) then 
!            !    cz = 250!250!0.35!5.0 !0.5 
!            !else
!            !    !ceps = 0.5 - ((Dr - 0.55) * 1.5)
!            !    cz = -(125*Dr) + 243.75
!            !end if
!            
!            cz = 250
!            
!        end if 
!    
!        ! initialize parameter 8
!        ! this is called ce in opensees but here it is renamed to prevent confusion
!        ! these values are reported by Long Chen (2020)
!        if (Dr > 0) then
!            if (Dr >= 0.75) then 
!                ceps = 0.2!0.2!1.0!0.2
!            else if (Dr <= 0.55) then 
!                ceps = 0.5!0.35!5.0 !0.5 
!            else
!                ceps = 0.5 - ((Dr - 0.55) * 1.5)
!            end if
!        end if 
!    
!    
!        ! initialize parameter 9        
!        if (phi_cv < 0) then
!            phi_cv = 33.0!6.0!33.0!36!33.0!36!33.0!35.6!33.0
!        end if 
!        
!        if (Mc < 0) then 
!            Mc = 2.0 * sin(phi_cv * (3.14159265359/180.0) )
!        end if 
!    
!        ! initialize parameter 10
!        if (nu < 0) then
!            nu = 0.3!0.3 !0.0
!        end if 
!    
!        ! initialize parameter 11
!        if (Cgd < 0) then
!            Cgd = 2.0
!        end if 
!    
!        ! initialize parameter --> this is from Version 3.3
!        if (Cdr < 0) then
!            Cdr = (5.0 + ( 25.0 * (Dr - 0.35)) ) ! this becomes negative if Dr is less than 0.35
!            Cdr = min(Cdr, 10.0)
!        end if 
!    
!        ! initialize parameter 12
!        if (Ckaf < 0) then
!            
!            Ckaf = 5.0 + (220.0 * ((Dr - 0.26)**3.0))
!            
!            if (Ckaf > 35.0) then 
!                Ckaf = 35.0
!            end if
!            
!            if (Ckaf < 4.0) then 
!                Ckaf = 4.0
!            end if
!        
!        end if 
!        
!        ! initialize parameter 13
!        if (QQ_Bolton < 0) then
!            QQ_Bolton = 10.0 !9.5!10.0            !9.5!10.0
!        end if 
!        
!        ! initialize parameter 14
!        if (RR_Bolton < 0) then
!            RR_Bolton = 1.5!0.7!1.5!0.7!1.5
!        end if 
!    
!        ! initialize parameter 15
!        if (mm < 0) then
!            mm = 0.01
!        end if 
!        
!        ! initialize parameter 
!        if (Fsed_min < 0) then
!            Fsed_min = (0.03 * exp(2.6 * Dr))
!            Fsed_min = min(Fsed_min, 0.99)
!        end if 
!        
!        ! initialize parameter 
!        if (p_sedo < 0) then
!            p_sedo = (P_atm / 5.0)
!        end if 
!    
!        ! store initial stress
!        initStress = Sigma_n ! local parameter
!        
!        ! find p0 
!        ! Step 2a: calculate mean effective stresses
!        p0 = 0.5 * GetTrace(initStress) ! mean effective stress 
!    
!        ! minimum p' 
!        if (Pmin < 0) then
!            Pmin = max(p0/200.0, P_atm/200.0)
!        end if
!        
!        ! p_min for stress
!        if (Pmin2 < 0) then
!            Pmin2 = Pmin * 10 ! not sure why this is 10 --> this is correct see table 3.2
!        end if 
!        
!        ! check if p0 is less than m_Pmin --> tension cutoff aspects
!        ! step 2b: if tensile stresses
!        if (p0 < Pmin) then
!        
!            ! initial p is small, set p to p_min and store the difference(mSigmab), the difference
!		    ! will be added to the stress returned to element
!            Sigma_n = Pmin * I1 ! change stress
!            Sigma_b = initStress - Sigma_n ! stress difference 
!            
!            p0 = Pmin ! corresponds to the minimum --> this should be m_Pmin2, added this negative here --> removed... follow geotech convention
!            
!            Alpha = 0 ! initialize to zero
!            Alpha_n = 0 ! initialize to zero
!        
!        else
!        
!            Sigma_n = initStress ! same stress
!            Sigma_b = 0.0 ! stress difference is zero
!            
!            ! alpha gets updated every time step at the beginning 
!            Alpha_n = GetDevPart(initStress) ! get the deviatoric part of the stress
!            Alpha_n = Alpha_n/p0 ! calculate old Alpha which is the deviatoric part normalized by the mean effective stress
!    
!        end if
!    
!        ! Step 3 caluclate relative state parameter and subsequently calculate the bounding and dilatancy stress ratios 
!        ! and Ado (from input property Dr and secondary parameters R, Q, nb, nd)
!        ! calculate the relative state parameter
!        ksi = GetKsi(Dr, p0, RR_Bolton, QQ_Bolton, Pmin, P_atm)
!        
!        ! this is where we initialize m_z_max
!        if (z_max < 0) then 
!            z_max = min(0.7*exp(-6.1*ksi), 20.0)
!        end if
!    
!        ! bounding and dilatancy surface variations depends on whether it is dense of loose
!        if (ksi < 0) then  !DENSE ! --> I switched the signs because it was wrong in the Opensees implementation 
!            
!            ! dense of critical
!            Mb = Mc * exp(-1.0 * nb * ksi)
!            Md = Mc * exp(       nd * ksi)
!        
!            if (Ado < 0) then 
!                
!
!                if (Mb > 2.0) then 
!                    ! Warning, Mb is larger than 2, using Ado = 1.5.
!                    Ado = 1.5
!                else 
!                    Ado = 2.5 * &
!                     (asin(Mb/2.0) - asin(Mc/2.0)) / (Mb - Md) !--> equation without fabric effects
!                end if 
!                
!            end if 
!                !Mb
!        
!        else !LOOSE
!        
!            ! loose of critical 
!            Mb = Mc * exp(-1.0 * (nb/4.0) * ksi)
!            Md = Mc * exp(       (nd*4.0) * ksi)
!            
!            if (Ado < 0) then 
!                Ado = 1.24
!            end if 
!        
!        end if  
!        
!        
!        ! check if initial stresses are inside bounding and dilatancy surface
!        Mcut = max(Mb, Md) 
!        
!        Mfin = sqrt(2.0) * GetNorm_Contr(GetDevPart(Sigma_n))
!        Mfin = Mfin / p0
!        
!        ! check that initial stresses are inside the bounding surface (or dilatancy surface if it is greater) and compute the 
!        ! committed back-stress and stress ratio tensors from the stress tensor --> TABLE 3.2
!        if (Mfin > Mcut) then !--> if outside bounding or dilatancy
!        
!            ! current stress ratio calculation
!            rrr = (Sigma_n - (p0*I1)) * (1/p0) * (Mcut/Mfin) !--> current stress ratio
!            ! initial stress outside bounding/dilatancy surface, scale shear stress and store the difference(mSigma_b),
!		  ! the difference will be added to the stress returned to element to maintain global equilibrium
!        
!            Sigma_n = (p0*I1) + (rrr*p0)
!            Sigma_b = initStress - Sigma_n
!            Alpha_n = rrr * (Mcut - mm) * (1/Mcut)
!            
!        !else
!        !    
!        !    rrr = Alpha_n
!        
!        end if
!        
!        !---------------------------------------------------------------------------
!        ! zero cumulative fabric
!        zcum = 0.0
!        
!        ! 6. Calculate the initial values of elastic shear modulus, elastic bulk modulus, plastic modulus, dilatancy... dilatancy not included
!        ! calculate the elastic moduli (mK and mG)
!        call GetElasticModuli_(Sigma_n, zcum, z_max, nu, G0, Md, Mb, &
!         PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, &
!        p_sedo, Fsed_min, me2p, &
!        Mc) ! mK, mG, mMcur, mzcum are outputs
!        
!        ! calculate elastic stiffness (3x3)
!        Ce = GetStiffness(KK, GG) ! maybe we don't need this....
!        
!        ! Initializing the plastic modulus, K_p 
!        K_p = 100 * GG !--> plastic modulus
!        
!        Alpha = Alpha_n
!        
!        ! initialize to zero --> staying consistent with the Washington code
!        Alpha_in = 0.0 
!        Alpha_in_n = 0.0
!        Alpha_in_p = 0.0
!        Alpha_in_p_n = 0.0
!        !Alpha_in = Alpha_n
!        !Alpha_in_n = Alpha_n
!        !Alpha_in_p = Alpha_n
!        !Alpha_in_p_n = Alpha_n
!        
!        ! initialize to mAlpha_n --> current stress ratio
!        Alpha_in_true = Alpha_n ! alpha_in 
!        Alpha_in_true_n = Alpha_n
!	    Alpha_in_max = Alpha_n ! alpha_inMax
!	    Alpha_in_max_n = Alpha_n 
!	    Alpha_in_min = Alpha_n ! alpha_inMin
!	    Alpha_in_min_n = Alpha_n
!        
!        ! initialize fabric to zero... we start modifying the fabric at large strain
!        Fabric = 0.0
!        Fabric_n = 0.0
!        Fabric_in = 0.0
!        Fabric_in_n = 0.0
!        
!        ! elastic strain increments
!        dEpsilonE = 0 
!        
!        ! initialize fabric terms
!        zpeak = z_max / 100000.0 ! z_peak
!        pzp = max(p0, Pmin) / 100.0 ! p_zp
!        zxp = 0.0 ! zxp set to zero 
!        pzpFlag = .TRUE.!.true. 
!        
!        ! establish tolerances
!        if (TolF < 0) then 
!            TolF = 1.0e-10	 !// TolF       ! error threshold (yield surface)
!        end if 
!    
!        if (TolR < 0) then 
!            TolR = 1.0e-8	 !// TolR       ! error threshold --> do we even need this??
!        end if
!        
!        if (InitialVoidRatio < 0) then 
!            InitialVoidRatio = emax - ((emax - emin) * Dr) !e_init
!        end if 
!        
!        
!        if (NextVoidRatio < 0) then 
!            NextVoidRatio = 0.0 !emax - ((emax - emin) * Dr) !e_init
!        end if 
!    
!         !turn off first call 
!        FirstCall = 0.0!1.0!0.0
!        
!        ! we reset the strains to zero 
!        Epsilon = 0.0
!        Epsilon_n = 0.0
!        
!        EpsilonE = 0.0
!        EpsilonE_n = 0.0
!        
!        
!      
!        
!      end if ! first call
!    
!      end subroutine PM4SandFullConstructor
!    
!    
!      subroutine PM4SandCommitState(StrainIncrement, &
!       z_max, nu, G0, mm, Md, Mb, PostShake, Pmin, P_atm, Cgd, p_sedo, &
!      Fsed_min, me2p, Fabric, Fabric_in, &
!        Alpha_in, Alpha_in_p, Alpha_in_true, Alpha_in_max, &
!      Alpha_in_min, &
!        KK, GG, Mcur, rrr, Alpha_in_n, Alpha_n, Alpha_in_p_n, &
!      Alpha_in_true_n, Alpha_in_max_n, Alpha_in_min_n, Sigma_n, &
!        zcum, zpeak, Sigma, Alpha, Fabric_n, Fabric_in_n, &
!        Cep, &
!      Mc , &
!      nn, RR, &
!      Kp, &
!      Epsilon, Epsilon_n,&
!      EpsilonE, EpsilonE_n) ! inout
!    !-------------------------------------------------
!    ! To commit stress state and state parameters
!    !
!    !-------------------------------------------------
!    
!      implicit none 
!    
!    
!      ! input 
!      real(8), dimension(3), intent(in) ::  StrainIncrement
!      real(8), dimension(3), intent(in) ::  nn
!      real(8), dimension(3), intent(in) ::  RR
!    
!      real(8), intent(in) ::  z_max
!      real(8), intent(in) ::  nu
!      real(8), intent(in) ::  G0
!      real(8), intent(in) ::  mm
!      
!      real(8), intent(in) ::  Mc
!      real(8), intent(in) ::  Md
!      real(8), intent(in) ::  Mb
!      
!      logical, intent(in) ::  PostShake
!      real(8), intent(in) ::  Pmin
!      real(8), intent(in) ::  P_atm
!      real(8), intent(in) ::  Cgd
!      real(8), intent(in) ::  p_sedo
!      real(8), intent(in) ::  Fsed_min
!      logical, intent(in) ::  me2p
!      
!      real(8), dimension(3), intent(in) ::  Fabric
!      real(8), dimension(3), intent(in) ::  Fabric_in
!    
!      real(8), dimension(3), intent(in) ::  Alpha_in
!      real(8), dimension(3), intent(in) ::  Alpha_in_p
!      real(8), dimension(3), intent(in) ::  Alpha_in_true
!      real(8), dimension(3), intent(in) ::  Alpha_in_max
!      real(8), dimension(3), intent(in) ::  Alpha_in_min
!    
!      
!    
!      ! output 
!      real(8), intent(out) ::  KK
!      real(8), intent(out) ::  GG
!      real(8), intent(out) ::  Mcur
!    
!      real(8), dimension(3), intent(out) ::  rrr
!    
!      real(8), dimension(3), intent(out) ::  Alpha_in_n
!      real(8), dimension(3), intent(out) ::  Alpha_n
!      real(8), dimension(3), intent(out) ::  Alpha_in_p_n
!      real(8), dimension(3), intent(out) ::  Alpha_in_true_n
!      real(8), dimension(3), intent(out) ::  Alpha_in_max_n
!      real(8), dimension(3), intent(out) ::  Alpha_in_min_n
!    
!      real(8), dimension(3), intent(out) ::  Sigma_n
!    
!    
!      real(8), dimension(3,3) ::  Ce
!      real(8), dimension(3,3), intent(out) ::  Cep
!    
!    !real(8), intent(out) ::  NextVoidRatio
!    
!      ! inout 
!      real(8), intent(inout) ::  zcum
!      real(8), intent(inout) ::  zpeak
!      real(8), dimension(3), intent(inout) ::  Sigma
!      real(8), dimension(3), intent(inout) ::  Alpha
!    
!      real(8), dimension(3), intent(inout) ::  Fabric_n
!      real(8), dimension(3), intent(inout) ::  Fabric_in_n
!    
!      
!      
!      real(8), dimension(3), intent(in) ::  Epsilon
!      real(8), dimension(3), intent(out) ::  Epsilon_n
!      real(8), dimension(3), intent(in) ::  EpsilonE
!      real(8), dimension(3), intent(out) ::  EpsilonE_n
!      
!      ! local 
!      real(8), dimension(3) ::  dFabric !, intent(out)
!
!      real(8) :: dVolStrain
!      real(8) :: GetTrace_Sigma
!      real(8) :: pp, aa
!      real(8), intent(in) :: Kp
!      real(8) :: DoubleDot2_2_Contr_dFabric_dFabric
!      real(8) :: DoubleDot2_2_Contr_Fabric_Fabric
!    
!      real(8), dimension(3) :: I1
!    
!      I1(1) = 1
!      I1(2) = 1
!      I1(3) = 0
!    
!      call GetElasticModuli_(Sigma, zcum, z_max, nu, G0, Md, Mb, &
!      PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, &
!      p_sedo, Fsed_min, me2p, &
!      Mc)
!      !(Sigma, KK, GG, Mcur, zcum) !--> temporary subroutine
!    
!      !call GetTrace(Sigma, GetTrace_Sigma) !--> temporary subroutine
!    
!      if ( (Mcur>Mb) .and. (me2p==1.0) ) then
!        pp = 0.5 * GetTrace(Sigma) ! real type 
!        rrr = (Sigma - (pp*I1)) * (1/pp) * (Mb/Mcur)
!        Sigma = (pp*I1) + (rrr*pp)
!        Alpha = rrr * (Mb - mm)/Mb
!      end if
!    
!      Alpha_in_n = Alpha_in
!	Alpha_n = Alpha
!	Alpha_in_p_n = Alpha_in_p
!	Alpha_in_true_n = Alpha_in_true
!	Alpha_in_max_n = Alpha_in_max
!	Alpha_in_min_n = Alpha_in_min
!	Sigma_n = Sigma
!	Epsilon_n = Epsilon
!	EpsilonE_n = EpsilonE
!	dFabric = Fabric - Fabric_n
!    
!      ! update cumulated fabric 
!       ! call DoubleDot2_2_Contr(dFabric, dFabric, 
!      !& DoubleDot2_2_Contr_dFabric_dFabric) !--> temporary subroutine
!      zcum = zcum + sqrt(DoubleDot2_2_Contr(dFabric, dFabric)/2.0) !/2.0
!     
!      ! call DoubleDot2_2_Contr(Fabric, Fabric, 
!      !& DoubleDot2_2_Contr_Fabric_Fabric) !--> temporary subroutine
!      zpeak = max( sqrt(DoubleDot2_2_Contr(Fabric, Fabric)/2.0), zpeak)
!    
!      Fabric_n = Fabric
!      Fabric_in_n = Fabric_in
!      !DGamma_n = DGamma 
!    
!    !call GetTrace(StrainIncrement, dVolStrain) ! dVolStrain !--> temporary subroutine
!    !NextVoidRatio = VoidRatio - ( (1 + VoidRatio) * dVolStrain )
!    !
!      Ce = GetStiffness(KK, GG) !--> temporary subroutine
!      
!      !if (nn(1)**2.0 + nn(2)**2.0 + nn(3)**2.0 +nn(3)**2.0 > 0) then
!      !aa = nn(1)**2.0 + nn(2)**2.0 + nn(3)**2.0 +nn(3)**2.0
!      !end if 
!      !
!      !if (RR(1)**2.0 + RR(2)**2.0 + RR(3)**2.0 +RR(3)**2.0 > 0) then
!      !aa = RR(1)**2.0 + RR(2)**2.0 + RR(3)**2.0 +RR(3)**2.0
!      !end if 
!      
!      ! get tangent elastoplastic matrix
!      call GetElastoPlasticTangent(Sigma_n, Ce, RR, nn, Kp, Pmin, Cep) !--> temporary subroutine
!    !mCep_Consistent = mCe
!    
!      Cep = Ce
!    
!      end subroutine PM4SandCommitState
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
!    !---------------------------------------------------------------------
!    
!
!
!      subroutine revertToLastCommit()
!      ! Need to be added
!      end subroutine revertToLastCommit
!
!!---------------------------------------------------------------------
!
!
!    !subroutine revertToStart()
!    !if (ops_InitialStateAnalysis) then
!    !    ! Do nothing, keep state variables from last step
!    !else
!    !    ! Normal call for revertToStart (not initialStateAnalysis)
!    !    call initialize(mSigma)
!    !endif
!    !
!    !end subroutine revertToStart
!
!!---------------------------------------------------------------------
!
!      subroutine getType(materialType)
!      ! Outputs
!      character(LEN=*) :: materialType
!      ! Set material type
!      materialType = "PlaneStrain"
!      end subroutine getType
!
!!---------------------------------------------------------------------
!
!      subroutine getOrder(materialOrder)
!      ! Outputs
!      integer, intent(out) :: materialOrder
!      ! Set material order
!      materialOrder = 3
!      end subroutine getOrder
!
!    
!    
!    
!    !------------------------------------------------------------------------
!
!      subroutine GetState(State)
!      implicit none
!      ! Output
!      real(8), dimension(16), intent(out) :: State
!
!      ! variables
!      real(8), dimension(3) :: mEpsilonE
!      real(8), dimension(3) :: mAlpha_n
!      real(8), dimension(3) :: mFabric_n
!      real(8), dimension(3) :: mAlpha_in_n
!      real(8) :: mVoidRatio
!      real(8) :: mDGamma_n
!      real(8) :: mG
!      real(8) :: mKp
!    
!      ! Assemble state parameters into the result vector
!      State(1:3) = mEpsilonE
!      State(4:6) = mAlpha_n
!      State(7:9) = mFabric_n
!      State(10:12) = mAlpha_in_n
!      State(13) = mVoidRatio
!      State(14) = mDGamma_n
!      State(15) = mG
!      State(16) = mKp
!
!
!      end subroutine GetState
!
!    !------------------------------------------------------------------------
!
!      subroutine GetAlpha(Alpha)
!      implicit none
!      ! Output
!      real(8), dimension(3), intent(out) :: Alpha
!    
!      ! variables 
!      real(8) :: mAlpha_n
!
!      ! Return alpha tensor
!      Alpha = mAlpha_n
!
!      end subroutine GetAlpha
!
!    !------------------------------------------------------------------------
!
!      subroutine GetFabric(Fabric)
!      implicit none
!      ! Output
!      real(8), dimension(3), intent(out) :: Fabric
!    
!      real(8) :: mFabric_n 
!
!      ! Return fabric tensor
!      Fabric = mFabric_n
!
!      end subroutine GetFabric
!
!    !------------------------------------------------------------------------
!
!      subroutine GetAlpha_in(Alpha_in)
!      implicit none
!      ! Output
!      real(8), dimension(3), intent(out) :: Alpha_in
!
!      real(8) :: mAlpha_in_n
!    
!      ! Return alpha_in tensor
!      Alpha_in = mAlpha_in_n
!
!      end subroutine GetAlpha_in
!
!    !------------------------------------------------------------------------
!
!      subroutine GetTracker(Tracker)
!      implicit none
!      ! Output
!      real(8), dimension(3), intent(out) :: Tracker
!
!      ! variables 
!      real(8) :: mTracker
!    
!      ! Return tracker vector
!      Tracker = mTracker
!
!      end subroutine GetTracker
!
!    !------------------------------------------------------------------------
!
!
!      subroutine GetKp(Kp)
!      implicit none
!      ! Output
!      real(8), intent(out) :: Kp
!
!      ! variables 
!      real(8) :: mKp
!    
!      ! Return Kp
!      Kp = mKp
!
!      end subroutine GetKp
!
!    !------------------------------------------------------------------------
!
!
!      subroutine GetG(G)
!      implicit none
!      ! Output
!      real(8), intent(out) :: G
!
!      ! variables 
!      real(8):: mG
!    
!      ! Return shear modulus
!      G = mG
!
!      end subroutine GetG
!
!    !------------------------------------------------------------------------
!
!
!      subroutine GetAlpha_in_p(Alpha_in_p)
!      implicit none
!      ! Output
!      real(8), dimension(3), intent(out) :: Alpha_in_p
!
!      real(8), dimension(3) :: mAlpha_in_p_n
!    
!      ! Return previous alpha_in tensor
!      Alpha_in_p = mAlpha_in_p_n
!
!
!      end subroutine GetAlpha_in_p
!
!    !------------------------------------------------------------------------
!
!
!      subroutine GetDGamma(DGamma)
!      implicit none
!      ! Output
!      real(8), intent(out) :: DGamma
!    
!      real(8) :: mDGamma_n
!
!      ! Return previous L
!      DGamma = mDGamma_n
!
!      end subroutine GetDGamma
!    
!    !------------------------------------------------------------------------
!
!    !subroutine GetTangent(TangType, Ce, Cep, Cep_Consistent, Tangent)
!    !implicit none
!    !! Inputs
!    !integer, intent(in) :: TangType
!    !real(8), dimension(3,3), intent(in) :: Ce!, Cep, Cep_Consistent
!    !! Output
!    !real(8), dimension(3,3), intent(out) :: Tangent
!    !
!    !! Determine which tangent to return based on TangType
!    !if (TangType == 0) then
!    !    Tangent = Ce
!    !else if (TangType == 1) then
!    !    Tangent = Cep
!    !else
!    !    Tangent = Cep_Consistent
!    !endif
!    !
!    !
!    !end subroutine GetTangent
!
!
!    !------------------------------------------------------------------------
!
!      subroutine getInitialTangent(initialTangent)
!      ! Input/output variables
!      real(8), dimension(6,6), intent(out) :: initialTangent
!    
!      ! variables 
!      real(8), dimension(6,6) :: mCe
!    
!      ! Return initial tangent matrix
!      initialTangent = mCe
!    
!      end subroutine getInitialTangent
!
!    
!
!    !------------------------------------------------------------------------
!
!      subroutine getStress(stress)
!      ! Input/output variables
!      real(8), dimension(3), intent(out) :: stress
!    
!      ! variables 
!      real(8), dimension(3) :: mSigma
!      real(8), dimension(3) :: mSigma_b
!    
!      ! Calculate and return stress
!      stress = mSigma + mSigma_b
!      stress = stress * (-1.0)  ! -1.0 is for geotechnical sign convention
!    
!      end subroutine getStress
!
!    
!        
!    !------------------------------------------------------------------------
!
!      subroutine getStrain(strain)
!      ! Input/output variables
!      real(8), dimension(3), intent(out) :: strain
!    
!    
!      ! variables 
!      real(8), dimension(3) :: mEpsilon
!    
!      ! Calculate and return strain
!      strain = mEpsilon
!      strain = strain * (-1.0)  ! -1.0 is for geotechnical sign convention
!    
!      end subroutine getStrain
!    
!    !------------------------------------------------------------------------
!
!    
!    !subroutine GetStress(Stress, Sigma_r)
!    !implicit none
!    !! Inputs
!    !real(8), dimension(6), intent(in) :: Stress
!    !! Output
!    !real(8), dimension(6), intent(out) :: Sigma_r
!    !
!    !! Compute stress
!    !Sigma_r = Stress + mSigma_b
!    !Sigma_r = Sigma_r * (-1.0)
!    !
!    !
!    !end subroutine GetStress
!    
!    
!    
!    !-----------------------------------------------------------------------
!      subroutine GetElasticStrain(Epsilon, EpsilonE_r)
!      implicit none
!      ! Inputs
!      real(8), dimension(6), intent(in) :: Epsilon
!      ! Output
!      real(8), dimension(6), intent(out) :: EpsilonE_r
!
!      ! Compute elastic strain
!      EpsilonE_r = Epsilon * (-1.0)
!
!      end subroutine GetElasticStrain
!    !subroutine getElasticStrain(elasticStrain)
!    !! Input/output variables
!    !real*8, dimension(3), intent(out) :: elasticStrain
!    !
!    !! Calculate elastic strain
!    !elasticStrain = mEpsilon
!    !elasticStrain = elasticStrain * (-1.0)
!    !
!    !end subroutine getElasticStrain
!    
!    
!    
!    
!    !------------------------------------------------------------------------------------
!    !// ---------------------------------------------------------------------------------
!    !/*************************************************************/
!    !// Plastic Integrator
!    !/*************************************************************/
!    
!      subroutine PM4SandIntegrate(Epsilon, Epsilon_n, &
!      EpsilonE, EpsilonE_n,&
!      G0, hpo, PostShake, me2p,  &
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, &
!      ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, &
!      zcum, pzp, zxp,  &
!      KK, GG, Mcur, rrr, Sigma_b, Fabric_n, Fabric, Fabric_in_n, &
!      Fabric_in, Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, &
!      Alpha_in_p, Alpha_in_true_n, Alpha_in_true, Alpha_in_max_n, &
!      Alpha_in_max, Alpha_in_min_n, Alpha_in_min, pzpFlag, TolF, TolR, &
!      dEpsilonE, Sigma_n, Sigma, &
!      Ce, &
!      nn, RR, &
!      NextVoidRatio)
!    
!      implicit none
!    
!      ! input variables
!      real(8), intent(in) :: G0       ! primiary parameter 2 ! prop 2
!      real(8), intent(in) :: hpo ! prop 3
!    
!      logical, intent(in) :: PostShake ! prop 4 !real(8)
!      logical, intent(in) :: me2p ! prop 5 !real(8)
! 
!      real(8), intent(in) :: P_atm    ! statev 1
!      real(8), intent(in) :: h0       ! statev 2              ! secondary parameter 1
!      real(8), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
!      real(8), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
!      real(8), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
!      real(8), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
!      real(8), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
!      real(8), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
!      real(8), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
!      real(8), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
!      real(8), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
!      real(8), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11
!
!      real(8), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12
!    
!      real(8), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
!      real(8), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
!    
!      real(8), intent(in) :: mm         ! statev 17                             ! secondary parameter 15
!
!    
!      real(8), intent(in) :: Pmin       ! statev 25 
!      real(8), intent(in) :: Pmin2      ! statev 26
!
!    
!      real(8), intent(in) :: zpeak      ! statev 28
!      real(8), intent(in) :: zcum       ! statev 29
!    
!    
!      real(8), intent(in), dimension(3) :: Fabric_n ! statev 43:45
!    
!      real(8), intent(in), dimension(3) :: Fabric_in_n ! statev 49:51
!    
!      real(8), intent(in), dimension(3) :: Alpha_n ! statev 55:57
!    
!      real(8), intent(in), dimension(3) :: Alpha_in_n ! statev 61:63
!    
!      real(8), intent(in), dimension(3) :: Alpha_in_p_n ! statev 67:69
!    
!      real(8), intent(in), dimension(3) :: Alpha_in_true_n ! statev 73:75
!    
!      real(8), intent(in), dimension(3) :: Alpha_in_max_n ! statev 79:81
!    
!      real(8), intent(in), dimension(3) :: Alpha_in_min_n ! statev 85:87
!    
!      real(8), intent(in) :: TolF ! statev 92
!    
!	real(8), intent(in) :: TolR ! statev 93 
!    
!      ! output variables 
!      real(8), intent(out), dimension(3) :: Fabric ! statev 46:48
!      real(8), intent(out), dimension(3) :: Fabric_in ! statev 52:54
!      real(8), intent(out), dimension(3) :: Alpha ! statev 58:60
!
!    
!      real(8), intent(out), dimension(3) :: dEpsilonE ! statev 94:96
!    
!      real(8), intent(out), dimension(3) :: Sigma ! statev 101:103
!    
!      real(8), intent(inout), dimension(3,3) :: Ce !
!    
!    ! inout variables
!	real(8), intent(in), dimension(3) :: Sigma_n ! statev 97:100                 ! current strain      
!    
!
!      ! Inout variables
!      real(8), dimension(3) :: StrainIncrement!strain_from_element !, intent(in)
!      real(8), dimension(3), intent(in) :: Epsilon
!      real(8), dimension(3), intent(in) :: Epsilon_n
!      real(8), dimension(3), intent(inout) :: EpsilonE
!      real(8), dimension(3), intent(inout) :: EpsilonE_n
!    
!      real(8), intent(in) :: InitialVoidRatio   ! statev 18
!
!      real(8), intent(inout) :: Mc         ! statev 19                      ! secondary parameter 1
!      real(8), intent(inout) :: Mb         ! statev 20 
!      real(8), intent(inout) :: Md         ! statev 21
!    
!      real(8), intent(in) :: Cdr        ! statev 22                       ! secondary parameter 
!    
!      real(8), intent(inout) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
!      real(8), intent(inout) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
!    
!      real(8), intent(inout) :: K_p        ! statev 27
!
!      real(8), intent(inout) :: pzp        ! statev 30 
!      real(8), intent(inout) :: zxp        ! statev 31
!    
!      real(8), intent(inout) :: KK         ! statev 34
!      real(8), intent(inout) :: GG         ! statev 35  
!      real(8), intent(inout) :: Mcur       ! statev 36
!    
!      real(8), intent(inout), dimension(3) :: rrr   ! statev 37:39
!      real(8), intent(inout), dimension(3) :: Sigma_b ! statev 40:42
!    
!    
!      real(8), intent(inout), dimension(3) :: Alpha_in ! statev 64:66
!
!      real(8), intent(inout), dimension(3) :: Alpha_in_p ! statev 70:72
!      real(8), intent(inout), dimension(3) :: Alpha_in_true ! statev 76:78
!    
!    
!      real(8), intent(inout), dimension(3) :: Alpha_in_max ! statev 82:84
!    
!      real(8), intent(inout), dimension(3) :: Alpha_in_min ! statev 88:90
!    
!    
!    !integer(4), intent(inout) :: pzpFlag ! statev ! statev 91
!      logical, intent(inout) :: pzpFlag ! statev ! statev 91
!    
!    ! Local variables
!      real(8), dimension(3) :: n_tr
!      real(8), dimension(3) :: tmp0
!      real(8), dimension(3) :: tmp1
!      real(8), dimension(3) :: Alpha_Alpha_in_true
!    
!      real(8), intent(out) :: NextVoidRatio
!      real(8) :: pp
!      real(8) :: zxpTemp
!      real(8) :: dGamma
!    
!      integer(4) :: ii
!    
!      real(8) :: DoubleDot2_2_Contr_result
!      real(8) :: GetTrace_result
!      real(8) :: GetNorm_Contr_result
!
!      real(8), intent(out), dimension(3) :: nn
!      real(8), intent(out), dimension(3) :: RR
!    
!    
!    
!      ! Assignments, new = old
!      Alpha = Alpha_n ! set current alpha as the previous alpha
!      Alpha_in = Alpha_in_n ! set initial alpha as the previous initial alpha !--> EFFECTIVE INITIAL
!      Alpha_in_true = Alpha_in_true_n ! true initial
!      Alpha_in_p = Alpha_in_p_n ! previous initial
!      Alpha_in_max = Alpha_in_max_n ! maximum --> used for negative loading 
!      Alpha_in_min = Alpha_in_min_n ! minimum --> used for positive loading
!      Fabric = Fabric_n ! set the current fabric as the previous fabric
!      Fabric_in = Fabric_in_n ! set current initial fabric as the previous initial fabric
!
!      ! Step 4: calculate the trail elastic stress increment and trial elastic stress:
!      ! Calculate trial stress
!      StrainIncrement = Epsilon - Epsilon_n
!      Sigma = Sigma_n + matmul(Ce, StrainIncrement)!the third component has got to be the GAMMMA ENGINEERING SHEAR STRAIN because 
!     
!      ! using the current mAlpha and the trial stress tmp0, get the normal to the yield surface in the stress ratio space
!      ! find the normal to the yield surface
!      n_tr = GetNormalToYield(Sigma, Alpha) 
!    
!      ! Check loading reversal condition
!      ! find the difference between current and initial stress ratios
!      Alpha_Alpha_in_true = Alpha - Alpha_in_true
!    
!      ! dot the back stress ratio difference with the normal 
!      if ( (DoubleDot2_2_Contr(Alpha_Alpha_in_true, n_tr) < 0.0) .and. &
!      (me2p == .TRUE.) ) then !--> me2p will actually be 1
!        
!        ! this is a load reversal
!        Alpha_in_p = Alpha_in ! --> previous initial becomes the initial 
!        Alpha_in_true = Alpha
!        Fabric_in = Fabric ! --> initial fabric set to current due to this being a load reversal
!        
!        ! Update pzp: mean effective stress at the time that zp achieves its greates value
!        pp = GetTrace(Sigma_n) !, pp
!        pp = 0.5 * pp
!        
!        if (pp <= Pmin) then
!            pp = Pmin
!        endif
!        
!        zxpTemp = GetNorm_Contr(Fabric_n)
!        zxpTemp = zxpTemp * pp
!        
!        ! record the maximum value of zp and its corresponding mean effective stress, pzp
!        if ( ((zxpTemp > zxp) .and. (pp > pzp)) .or. (pzpFlag==.TRUE.)) &
!      then
!            zxp = zxpTemp
!            pzp = pp
!            pzpFlag = .FALSE.
!        end if
!        
!        ! Track initial back-stress ratio history
!        do ii = 1, 3
!            if (Alpha_in(ii) > 0.0) then ! --> if positive 
!                ! Minimum positive value
!                Alpha_in_min(ii) = min(Alpha_in_min(ii), Alpha(ii)) ! minimum back-stress ratio
!            else ! --> if negative
!                ! Maximum negative value
!                Alpha_in_max(ii) = max(Alpha_in_max(ii), Alpha(ii)) ! maximum back-stress ratio
!            endif
!        enddo
!        
!        ! Update mAlpha_in based on loading direction
!        if (Alpha(3) * Alpha_in_p(3) > 0) then
!            
!            do ii = 1, 3
!                if (n_tr(ii) > 0.0) then
!                    ! Positive loading direction
!                    Alpha_in(ii) = max(0.0, Alpha_in_min(ii))
!                else
!                    ! Negative loading direction
!                    Alpha_in(ii) = min(0.0, Alpha_in_max(ii))
!                endif
!            enddo
!            
!        else
!            
!            ! when would this happen?
!            ! alpha_xy spans the positive and negative quadrants then avoid the min/max business
!            Alpha_in = Alpha
!        endif
!    
!    
!        end if
!
!    
!        ! Force elastic response    
!        if (me2p == .FALSE.) then !--> this will never be the case --> probably for testing the model assuming elastic behavior
!
!      !       call elastic_integrator(Sigma_n, StrainIncrement, 
!      !& G0, nu, Pmin, P_atm, 
!      !&  dEpsilonE, Sigma, Alpha, GG, KK, Ce, VoidRatio, NextVoidRatio)
!         !--> maybe to initialize stresses in an elastic manner without changing constitutive model
!    
!            
!            call elastic_integrator(Sigma_n, Epsilon_n, & 
!       EpsilonE_n, Epsilon,&
!       EpsilonE, Sigma, Alpha, NextVoidRatio, GG, &
!       KK, Ce, G0, nu, Pmin,& 
!       P_atm, InitialVoidRatio, me2p)
!            
!            
!            else
!        ! me2p -- .true.
!        ! ElastoPlastic response --> always go here
!        ! Explicit schemes
!        call explicit_integrator(&
!      Sigma_n, &
!      Epsilon_n,&
!      EpsilonE_n,  &
!      Alpha_n, &
!      Fabric_n, Alpha_in,&
!      Alpha_in_p, Epsilon, EpsilonE, &
!      Sigma, Alpha, Fabric, DGamma, NextVoidRatio, &
!      GG, KK, Ce, &
!      TolF,  &
!      G0, hpo, PostShake, me2p,&
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, &
!      nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &!phi_cv, 
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum,&
!      Fabric_in,&
!       MCur, &
!      pzp, &
!      alpha_in_true, &
!      nn, RR  )
!    
!        end if
!        
!
!        end subroutine PM4SandIntegrate
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
!    !------------------------------------------------------------------------
!    !/*************************************************************/
!    !// Elastic Integrator
!    !/*************************************************************/
!    
!    
!        subroutine elastic_integrator(CurStress, CurStrain, &
!       CurElasticStrain, NextStrain, NextElasticStrain, &
!       NextStress, NextAlpha, NextVoidRatio, GG, KK, aC,    &
!       G0, nu, Pmin, P_atm, &! input
!       InitialVoidRatio, me2p) ! output 
!    
!    ! assuming linear elasticity for stress initialization
!      implicit none
!    
!    ! input variables
!      real(8), dimension(3), intent(in) :: CurStress
!      real(8), dimension(3), intent(in) :: CurStrain
!      real(8), dimension(3), intent(in) :: CurElasticStrain
!      
!      real(8), dimension(3), intent(in) :: NextStrain
!      real(8), dimension(3), intent(out) :: NextElasticStrain
!      real(8), dimension(3), intent(out) :: NextStress
!      real(8), dimension(3), intent(out) :: NextAlpha
!      
!      real(8), intent(out) :: NextVoidRatio
!      
!      real(8), intent(out) :: GG
!      real(8), intent(out) :: KK 
!      
!      real(8), dimension(3, 3), intent(out) :: aC
!      
!      real(8), intent(in) :: G0
!      real(8), intent(in) :: nu
!      real(8), intent(in) :: Pmin
!      real(8), intent(in) :: P_atm
!    
!      real(8), intent(in) :: InitialVoidRatio
!      
!      !real(8), intent(in) :: nu
!    
!      ! Output variables
!      
!    ! Local variables
!      real(8) :: pp
!      real(8) :: GetTrace_result
!      real(8), dimension(3) ::dStrain ! DoubleDot4_2_result
!      !real(8), dimension(3) :: DetDevPart_result
!      
!      logical, intent(in) :: me2p
!    
!    !! Calculate the next elastic strain increment (assume all elastic)
!    !  NextElasticStrainIncrement = NextStrainIncrement
!      dStrain = NextStrain - CurStrain
!      
!    
!    ! Calculate void ratio 
!      !call GetTrace(NextElasticStrainIncrement)!, GetTrace_result)
!      NextVoidRatio = InitialVoidRatio - &
!       ( (1 + InitialVoidRatio) * GetTrace(NextStrain))
!    
!      
!      
!      NextElasticStrain = CurElasticStrain + dStrain 
!      
!    ! Calculate elastic moduli
!      call GetElasticModuli(CurStress, KK, GG, G0, nu, P_atm, Pmin, &
!        me2p)
!      
!      
!      
!      ! call GetElasticModuli_(NextStress, zcum, z_max, nu, G0, Md, Mb, 
!      !&   PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, 
!      !&  Fsed_min, me2p, 
!      !& Mc) ! --> csr uses mMcur
!    
!    ! Calculate elastic matrix (3x3)
!      aC = GetStiffness(KK, GG)
!    
!    ! Calculate the next stress
!      !call DoubleDot4_2(aC, NextStrainIncrement, DoubleDot4_2_result)
!      NextStress = CurStress + DoubleDot4_2(aC, dStrain)
!    
!    ! Calculate mean effective stress
!      !call GetTrace(NextStress, pp)
!      pp = 0.5*GetTrace(NextStress) !sum(NextStress)
!    
!    ! just to avoid dividing by zero in case pp is very small 
!    ! (fyi won't really harden if in elastic integrator state)
!      if (pp > Pmin) then 
!        !call , DetDevPart_result)
!        NextAlpha = GetDevPart(NextStress) / pp
!      end if
!    
!      end subroutine elastic_integrator
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
!    !-------------------------------------------------------------------------
!    !// ----------------------------------------------------------------------
!    !/*************************************************************/
!    !// Explicit Integrator
!    !/*************************************************************/
!    
!      subroutine explicit_integrator(CurStress, &
!     & CurStrain, CurElasticStrain, CurAlpha, CurFabric, &
!     & Alpha_in, Alpha_in_p, NextStrain, NextElasticStrain, &
!     & NextStress, NextAlpha, NextFabric, NextL, &
!     & NextVoidRatio, &
!     & GG, KK, aC, &! --> up until here it is the same 
!     & TolF, &!--> these are extra needed variables 
!     & G0, hpo, PostShake, me2p, &
!     & P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, Cgd, &
!     & Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!     & Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, zpeak, zcum,&
!     & Fabric_in, &
!     & MCur, &
!     & pzp, &
!     & alpha_in_true, &
!     & nn, RR) !, aCep, aCep_Consistent)
!    !---------------------------------------------------------------------------------------
!    ! Description:   
!    ! explicit elastoplastic integrator that should work with ForwardEuler and ModifiedEuler
!    !---------------------------------------------------------------------------------------
!
!      implicit none
!    
!      ! Input variables ---------------------------------------------------------------------
!      real(8), dimension(3), intent(in) :: CurStress
!      real(8), dimension(3), intent(in) :: CurStrain
!      real(8), dimension(3), intent(in) :: CurElasticStrain
!      
!      real(8), dimension(3), intent(in) :: CurAlpha
!      real(8), dimension(3), intent(in) :: CurFabric
!      
!      real(8), dimension(3), intent(in) :: Alpha_in
!      real(8), dimension(3), intent(in) :: Alpha_in_p
!      
!      real(8), dimension(3), intent(in) :: NextStrain
!      
!      real(8), intent(in) :: TolF
!      
!      real(8), intent(in) :: G0       ! primiary parameter 2 ! prop 2
!      real(8), intent(in) :: hpo ! prop 3
!    
!      logical, intent(in) :: PostShake ! prop 4 !real(8)
!      logical, intent(in) :: me2p ! prop 5 !real(8)
! 
!      real(8), intent(in) :: P_atm    ! statev 1
!      real(8), intent(in) :: h0       ! statev 2              ! secondary parameter 1
!      real(8), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
!      real(8), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
!      real(8), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
!      real(8), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
!      real(8), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
!      real(8), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
!      real(8), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
!      real(8), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
!
!      real(8), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
!      real(8), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11
!    
!      real(8), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12
!    
!      real(8), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
!      real(8), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
!    
!      real(8), intent(in) :: mm         ! statev 17                             ! secondary parameter 15
!
!      real(8), intent(in) :: InitialVoidRatio
!
!    
!      real(8), intent(in) :: Pmin       ! statev 25 
!      real(8), intent(in) :: Pmin2      ! statev 26
!
!    
!      real(8), intent(in) :: zpeak      ! statev 28
!      real(8), intent(in) :: zcum       ! statev 29
!
!      real(8), dimension(3), intent(in) :: Fabric_in
!      real(8), intent(in) :: Cdr        ! statev 22                       ! secondary parameter 
!
!      
!      ! Output variables ---------------------------------------------------------------------
!      real(8), dimension(3), intent(out) :: NextElasticStrain
!      real(8), dimension(3), intent(out) :: NextStress
!      real(8), dimension(3), intent(out) :: NextAlpha
!      real(8), dimension(3), intent(out) :: NextFabric
!      real(8), intent(out) :: NextL
!      real(8), intent(out) :: NextVoidRatio
!      
!      real(8), dimension(3, 3), intent(out) :: aC !Ce !out as it is calculated here
!
!      real(8), dimension(3), intent(out) :: nn
!      real(8), dimension(3), intent(out) :: RR
!      
!      
!      ! inout variables ---------------------------------------------------------------------
!      real(8), intent(inout) :: GG
!      real(8), intent(inout) :: KK
!
!      real(8), intent(inout) :: MCur
!    
!      real(8), intent(inout) :: Mc         ! statev 19                      ! secondary parameter 1
!      real(8), intent(inout) :: Mb         ! statev 20 
!      real(8), intent(inout) :: Md         ! statev 21
!    
!    
!      real(8), intent(inout) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
!      real(8), intent(inout) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
!    
!      real(8), intent(inout) :: K_p        ! statev 27
!
!      real(8), intent(inout) :: pzp
!      
!      real(8), dimension(3), intent(inout) :: alpha_in_true
!        
!    
!      ! Local variables ---------------------------------------------------------------------
!      real(8) :: a0
!      real(8) :: a1
!      real(8) :: elasticRatio ! previously named aa
!      
!      real(8) :: ff
!      real(8) :: fn
!      real(8) :: dVolStrain
!      
!      real(8), dimension(3) :: dElasStrain
!    
!      real(8), dimension(3) :: dSigma
!      real(8), dimension(3) :: dDevStrain 
!      
!      real(8), dimension(3) :: I1, dStrain
!    
!       real(8) :: debug
!       
!       real(8) :: denom, RatioValue
!    
!       ! Establish identity matrix. 
!       I1(1) = 1
!       I1(2) = 1
!       I1(3) = 0
!    
!       ! Calculate the next void ratio using the strain. 
!       ! InitialVoidRatio has a constant value. 
!       ! NextVoidRatio is updated each time step. 
!       NextVoidRatio = InitialVoidRatio &
!                    - ( (1 + InitialVoidRatio) * GetTrace(NextStrain) )
!    
!      ! Calculate total strain increment. 
!      dStrain = NextStrain - CurStrain
!      
!      ! Calculate NextElasticStrain assuming the entire 
!      ! strain is elastic. 
!      NextElasticStrain = CurElasticStrain + dStrain
!      
!      ! Calculate the volumetric strain increment. 
!      dVolStrain = GetTrace(dStrain)
!      
!      ! Calculate deviatoriic strain increment. 
!      dDevStrain = dStrain - (I1*(dVolStrain / 3.0))
!      
!      ! Calculate 3x3 stiffness matrix (aC) in 2D. 
!      aC = GetStiffness(KK, GG)
!      
!      ! calculate the stress increment assuming elastic behavior. 
!      dSigma = (KK * dVolStrain * I1) + &
!                 (2.0 * GG * ToContraviant(dDevStrain))
!    
!      ! Update the stress assuming an entirely elastic stress increment. 
!      NextStress = CurStress + dSigma
!    
!      ff = GetFYieldFunction(NextStress, CurAlpha, mm) ! ff is the final yield function evaluation
!
!      fn = GetFYieldFunction(CurStress, CurAlpha, mm) ! fn is the initial yield function evaluation
!    
!      nn = GetNormalToYield(NextStress, CurAlpha) ! nn is the normal of the yield function at the projected stress
!    
!      !-------------------------------------------------------------------------
!      ! Perform plasticity calculations
!      
!      if (ff <= TolF) then ! Pure elastic step (no transition)
!        
!        ! Pure elastic loading/unloading
!        NextAlpha = CurAlpha !--> no change
!        NextFabric = CurFabric
!        NextL = 0
!        RR = nn    
!        return
!     
!      elseif (fn < -TolF) then ! transition from elastic to plastic
!
!          ! establish the inputs to the intersection factor
!          a0 = 0.0
!          a1 = 1.0 
!
!          ! evaluate the intersection factor so that we can 
!          ! know the elastic to plastic ratio        
!          elasticRatio = &
!             IntersectionFactor(CurStress, CurStrain, &
!                                NextStrain, CurAlpha, &
!                                aC, a0, a1, TolF, mm) 
!        
!     
!        ! elastic portion of the strain
!        dElasStrain = dStrain * elasticRatio
!        
!        ! elastic portion of the stress
!        dSigma = DoubleDot4_2(aC, dElasStrain) 
!        
!        call MaxStrainInc(CurStress + dSigma, &! ModifiedEuler !ForwardEuler !MaxStrainInc
!     &                    CurStrain + dElasStrain, &
!     &                    CurElasticStrain + dElasStrain,&
!     &                    CurAlpha, CurFabric, &
!     &                    alpha_in, alpha_in_p,&
!     &                    NextStrain, NextElasticStrain,    &
!     &                    NextStress, NextAlpha, NextFabric, &
!     &                    NextL, NextVoidRatio,    &
!     &                    GG, KK,&
!     &                    aC, &! --> we stop here in the Opensees implementation
!     &                    G0, hpo, PostShake, me2p, &
!     &                    P_atm, h0, emax, emin, nb, nd, Ado, z_max, &
!     &                    cz, ceps, nu, Cgd, Ckaf, &
!     &                    QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!     &                    Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, &
!     &                    Pmin2, K_p, zpeak, zcum, &
!     &                    Fabric_in, &
!     &                    MCur,&
!     &                    TolF, &
!     &                    pzp, &
!     &                    alpha_in_true,&
!     &                    nn, RR  )
!        
!        return
!        
! 
!      else if (abs(fn) < TolF) then ! Pure plastic step or elastic unloading followed by plastic loading
!            
!        ! check denominator    
!        if (GetNorm_Contr(dSigma) == 0.0d0) then
!          denom = 1.0d0
!        else 
!          denom = GetNorm_Contr(dSigma)
!        end if    
!        
!        
!        RatioValue = &
!          DoubleDot2_2_Contr(GetNormalToYield(CurStress, CurAlpha),&
!                             dSigma) / denom
!
!            
!        
!        if ( RatioValue > (-sqrt(TolF)) ) then ! if true then no unloading Sloan et al. 2001
!        
!            ! Pure plastic step
!            call MaxStrainInc(CurStress, &! ModifiedEuler ! ForwardEuler !MaxStrainInc
!                             CurStrain, &
!                             CurElasticStrain,    &
!                             CurAlpha, CurFabric, &
!                             alpha_in, alpha_in_p, &
!                             NextStrain, NextElasticStrain, &
!                             NextStress, NextAlpha, NextFabric,&
!                             NextL, NextVoidRatio, &
!                             GG, KK, &
!                             aC,   &
!                             G0, hpo, PostShake, me2p, &
!                             P_atm, h0, emax, emin, nb, nd, Ado, z_max,&
!                             cz, ceps, nu, Cgd, Ckaf, &
!                            QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!                             Mc, Mb, Md, Cdr, Fsed_min, p_sedo, &
!                             Pmin, Pmin2, K_p, &
!                             zpeak, &
!                             zcum, &
!                             Fabric_in, &
!                             MCur, &
!                             TolF,     &
!                             pzp, &
!                             alpha_in_true,     &
!                             nn, RR )
!            
!            return 
!            
!        else
!            
!            ! Elastic unloading followed by plastic loading
!            elasticRatio = IntersectionFactor_Unloading(CurStress, &
!                         CurStrain, &
!                         NextStrain, CurAlpha, aC, TolF, mm) !CurStrain, NextStrain,
!            
!            dElasStrain = dStrain * elasticRatio
!            
!            dSigma = DoubleDot4_2(aC, dElasStrain)
!            
!            call MaxStrainInc(CurStress + dSigma, & ! ModifiedEuler ! ForwardEuler !MaxStrainInc
!      CurStrain + dElasStrain,&
!      CurElasticStrain + dElasStrain,&
!      CurAlpha, CurFabric, alpha_in, alpha_in_p,     &     
!      NextStrain, NextElasticStrain, NextStress, &
!      NextAlpha, NextFabric, NextL, NextVoidRatio, &
!      GG, KK, aC,   &
!      G0, hpo, PostShake, me2p, &
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps,nu,&
!      Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
!      zpeak, &
!      zcum, &
!      Fabric_in, &
!      MCur, &
!      TolF,     &
!      pzp, &
!      alpha_in_true,    & 
!      nn, RR   )                
!            
!            return 
!            
!        end if
!        
!    
!        else
!        
!        ! Correct the stress
!        ff = GetFYieldFunction(CurStress, CurAlpha, mm)
!        
!        ! Illegal stress state
!        !print *, "PM4Sand: Encountered an illegal stress state! Tag: "!, !getTag()
!        !print *, "            f = ", ff
!         
!        call MaxStrainInc(CurStress,   &!MaxStrainInc
!      CurStrain, &
!      CurElasticStrain,    &
!      CurAlpha, CurFabric, alpha_in, alpha_in_p, &
!      NextStrain, NextElasticStrain, &
!      NextStress, NextAlpha,&
!      NextFabric, NextL, NextVoidRatio, &
!      GG, KK, aC,   &
!      G0, hpo, PostShake, me2p, &
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps,nu,&
!      Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
!      zpeak, &
!      zcum, &
!      Fabric_in, &
!      MCur, &
!      TolF,     &
!      pzp, &
!      alpha_in_true,     &
!      nn, RR     )
!        
!        
!        return
!    
!    
!        end if
!
!
!        end subroutine explicit_integrator
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
!    !-------------------------------------------------------------------------
!    !// ----------------------------------------------------------------------
!    !/*************************************************************/
!    !// Forward-Euler Integrator
!    !/*************************************************************/
!      subroutine ForwardEuler(CurStress, &
!      CurStrain, &
!      CurElasticStrain,  &
!      CurAlpha, &
!      CurFabric, alpha_in, alpha_in_p,&
!      NextStrain, NextElasticStrain, &
!      NextStress, NextAlpha, NextFabric,&
!      NextL, NextVoidRatio, &
!      GG, KK, aC,&
!      G0, hpo, PostShake, me2p,     &
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, &
!      cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm,&
!      InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, &
!      Pmin2, K_p, zpeak, zcum, &
!      Fabric_in, &
!      MCur, &
!      TolF, &
!      pzp, &
!      alpha_in_true, &
!      nn, RR)
!
!      implicit none
!    
!      ! input variables
!      real(8), intent(in), dimension(3) :: CurStress 
!
!      real(8), intent(in), dimension(3) :: CurStrain 
!      real(8), intent(in), dimension(3) :: NextStrain
!
!      real(8), intent(in), dimension(3) :: CurElasticStrain 
!      real(8), intent(inout), dimension(3) :: NextElasticStrain 
!
!      real(8), intent(in), dimension(3) :: CurAlpha 
!      real(8), intent(in), dimension(3) :: CurFabric 
!      real(8), intent(in), dimension(3) :: alpha_in 
!      real(8), intent(in), dimension(3) :: alpha_in_p 
!    
!    
!    
!      ! input variables
!      real(8), intent(in) :: G0       ! primiary parameter 2 ! prop 2
!      real(8), intent(in) :: hpo ! prop 3
!    
!      logical, intent(in) :: PostShake ! prop 4 !real(8)
!      logical, intent(in) :: me2p ! prop 5 !real(8)
! 
!      real(8), intent(in) :: P_atm    ! statev 1
!      real(8), intent(in) :: h0       ! statev 2              ! secondary parameter 1
!      real(8), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
!      real(8), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
!      real(8), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
!      real(8), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
!      real(8), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
!      real(8), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
!      real(8), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
!      real(8), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
!      real(8), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
!      real(8), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11
!
!      real(8), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12
!
!      real(8), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
!      real(8), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
!    
!      real(8), intent(in) :: mm         ! statev 17                             ! secondary parameter 15
!
!    
!      real(8), intent(in) :: Pmin       ! statev 25 
!      real(8), intent(in) :: Pmin2      ! statev 26
!
!    
!      real(8), intent(in) :: zpeak      ! statev 28
!      real(8), intent(in) :: zcum       ! statev 29
!    
!      real(8), intent(in), dimension(3) :: Fabric_in
!      
!      real(8), intent(in) :: InitialVoidRatio
!    
!      real(8), dimension(3), intent(in) :: Alpha_in_true ! statev 73:75 !intent(in),
!    
!    ! output variables    
!      real(8), intent(out), dimension(3) :: NextAlpha
!      real(8), intent(out), dimension(3) :: NextFabric
!      real(8), intent(out) :: NextL
!      real(8), intent(out) :: NextVoidRatio
!      real(8), intent(out) :: GG
!      real(8), intent(out) :: KK
!      real(8), intent(inout), dimension(3,3) :: aC
!    
!    ! inout variables
!      real(8), intent(inout), dimension(3) :: NextStress
!
!      real(8), dimension(3) :: dStrain
!    
!      real(8), intent(in) :: Mc         ! statev 19                      ! secondary parameter 1
!      real(8), intent(inout) :: Mb         ! statev 20 
!      real(8), intent(inout) :: Md         ! statev 21
!    
!      real(8), intent(in) :: Cdr        ! statev 22                       ! secondary parameter 
!    
!      real(8), intent(in) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
!      real(8), intent(in) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
!    
!      real(8), intent(inout) :: Mcur       ! statev 36
!
!    ! Local variables
!      real(8) :: CurVoidRatio
!      real(8) :: CurDr
!      real(8) :: Cka
!      real(8) :: hh
!      real(8) :: pp
!      real(8) :: dVolStrain
!      real(8) :: DD
!      real(8) :: AlphaAlphaBDotN
!    
!      real(8), dimension(3), intent(out) :: nn
!      real(8), dimension(3), intent(out) :: RR
!      real(8), dimension(3) :: alphaD
!    
!      real(8), dimension(3) :: alphaB
!    
!      real(8), dimension(3) :: dPStrain
!      real(8), dimension(3) :: bb
!      real(8), dimension(3) :: dDevStrain
!      real(8), dimension(3) :: rrr
!      real(8), dimension(3) :: dSigma
!      real(8), dimension(3) :: dAlpha
!      real(8), dimension(3) :: dFabric    
!    
!      real(8), dimension(3) :: dSigmaP
!    
!      real(8), intent(inout) :: pzp
!      real(8), intent(out) :: K_p
!      real(8) :: DGamma
!      integer(4) :: DebugFlag
!    
!      real(8) :: ksi
!      real(8), intent(in) :: TolF
!    
!    ! temporary local variables 
!      real(8) :: temp4 !--> this temporary variable is used to check the denominator of L
!	
!    ! local variables with numeric values 
!      integer(4), dimension(3) :: I1 ! identity matrix
!      real(8) :: two3 ! two thirds
!      real(8) :: small 
!      
!      integer(4) :: ii
!      
!      small = 1e-5
!    
!    ! establish numeric values
!      I1(1) = 1
!      I1(2) = 1
!      I1(3) = 0
!    
!      two3 = 0.666666666666666667
!    
!    ! Get elastic moduli
!    ! you need 'NextStress' because of fabric aspects
!      call GetElasticModuli_(NextStress, zcum, z_max, nu, G0, Md, Mb, &
!        PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, &
!       Fsed_min, me2p, &
!      Mc) ! --> csr uses mMcur
!      
!      ! Calculate current void ratio
!      CurVoidRatio = InitialVoidRatio - &
!      ( (1 + InitialVoidRatio)*GetTrace(CurStrain) )
!      
!    ! Calculate current Dr
!      CurDr = (emax - CurVoidRatio) / (emax - emin) 
!
!    ! Calculate mean effective stress
!      pp = GetTrace(CurStress)
!      pp = 0.5 * pp !sum(CurStress)
!      if (pp < Pmin) pp = Pmin ! Apply tension cutoff
!    
!    ! Calculate void ratio from NextStrain
!      NextVoidRatio = InitialVoidRatio - &
!       ( (1.0 + InitialVoidRatio) * GetTrace(NextStrain) )
!    
!      ! Calculate strain increment
!      !// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
!      dStrain = NextStrain - CurStrain
!    
!      ! Calculate NextElasticStrain
!      NextElasticStrain = &
!       CurElasticStrain + dStrain !StrainIncrement!dStrain ! CurElasticStrain +
!      !// using NextStress instead of CurStress to get correct n
!	!// 1: calculate state parameters and initialize them
!    
!    ! Calculate state parameters
!      call GetStateDependent(NextStress, CurAlpha, alpha_in, alpha_in_p,&
!      alpha_in_true, CurFabric, Fabric_in, &
!        GG, zcum, zpeak, pzp, Mcur, CurDr, &
!        Mc, nd, nb, Pmin, Pmin2, P_atm, &
!        mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
!        QQ_Bolton, RR_Bolton, &
!        nn, alphaB, alphaD, bb, RR, &
!        Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, &
!        ksi)
!    
!    
!      ! Calculate volumetric strain increment
!      !// dVolStrain = GetTrace(NextStrain - CurStrain);
!      dVolStrain = GetTrace(dStrain)
!    
!      ! Calculate deviatoric strain increment
!	!// dDevStrain = (NextStrain - CurStrain) - dVolStrain / 3.0 * mI1;
!      dDevStrain = dStrain - (I1*(dVolStrain / 3.0) )
!      
!      ! Calculate stress ratio
!      rrr = GetDevPart(NextStress) / pp
!      
!      ! Calculate denominator term
!      temp4 = K_p + (2.0 * GG) - (KK * DD * DoubleDot2_2_Contr(nn, rrr))
!
!    
!    ! Check denominator if it is small
!      if (abs(temp4) < small) then !--> prevent division by zero
!    
!        ! Neutral loading
!        dSigma = 0.0
!        dAlpha = 0.0
!        dFabric = 0.0
!        ! dPStrain = dDevStrain + dVolStrain * mI1
!        dPStrain = dStrain 
!        
!      else
!    
!        ! Calculate L
!        NextL = ( (2.0 * GG * DoubleDot2_2_Mixed(nn, dDevStrain)) -  &
!      (DoubleDot2_2_Contr(nn, rrr) * KK * dVolStrain) ) / temp4 ! Equation 31
!        
!        ! Set mDGamma to NextL
!        DGamma = NextL ! mDGamma is the plastic mulitplier
!        if (NextL < 0.0) then
!            
!            ! If NextL is negative
!            if (debugFlag) then
!                !write(*,*) "NextL is smaller than 0"
!                !write(*,*) "NextL = ", NextL
!            end if
!
!            dSigma = (2.0 * GG * ToContraviant(dDevStrain)) + &
!                     (KK * dVolStrain * I1) !--> focusing only on the elastic increment if L is negative
!            dAlpha = 0.0
!            dFabric = 0.0
!            dPStrain = 0.0
!    
!        else
!            
!            ! If NextL is non-negative
!            ! Calculate dSigma
!            !dSigma = 2.0 * GG * ToContraviant(dDevStrain) + KK * dVolStrain * mI1 - Macauley(NextL) * &
!            !     (2.0 * GG * nn + KK * DD * mI1)
!            dSigma = (2.0 * GG * ToContraviant(dDevStrain)) + &
!                 (KK * dVolStrain * I1)  & ! elastic component
!          - ( Macauley(NextL) * ((2.0 * GG * nn) + (KK * DD * I1)) ) ! plastic component
!            
!            ! Update fabric
!            if (DoubleDot2_2_Contr(alphaD - CurAlpha, nn) < 0.0) then
!                ! dilation
!                
!                ! Update fabric according to Equation 57
!                dFabric = (z_max * nn) + CurFabric
!                
!                !//  dz in Equation 57 
!                dFabric = (-1.0 * &
!               (cz / (1 + Macauley( (zcum / (2.0*z_max)) - 1.0) )) * &
!               Macauley(NextL) * MacauleyIndex(-DD) * dFabric)
!                end if 
!            
!
!                ! Update alpha
!                !// dPStrain = NextL * mIIco * R;
!                !// dAlpha = two3 * NextL * h * b;
!            dPStrain = ToCovariant(RR) * NextL ! --> plastic strain (3x1)
!            dAlpha = bb * (two3 * NextL * hh) ! --> change in alpha (Hardening) (3x1)
!            
!            end if
!    
!    
!        end if
!
!        ! Update NextFabric, NextElasticStrain, NextStress, and NextAlpha
!        !// NextFabric = CurFabric + dFabric;
!        NextFabric = CurFabric + dFabric
!        
!        !// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
!        NextElasticStrain = NextElasticStrain - dPStrain 
!    
!        !// NextStress = CurStress + dSigma;
!        NextStress = CurStress + dSigma 
!        
!        !// NextAlpha = CurAlpha + dAlpha;
!        NextAlpha = CurAlpha + dAlpha
!
!        ! Perform stress correction
!        do ii = 1, size(NextStress)
!            if (NextStress(ii) /= NextStress(ii)) then
!                print *, "NextStress(", ii, ") is NaN!"
!            end if
!        end do
!            
!             call stress_correction_(NextStress, NextAlpha, alpha_in,&
!           alpha_in_p, CurFabric, NextVoidRatio, &
!            alpha_in_true, &
!      Fabric_in, CurStress, Pmin, P_atm, &
!                             TolF, Mc, emax, emin, zcum, zpeak, z_max, &
!      pzp, mm, h0, hpo, Cdr,&
!      Ceps, Ckaf,        &
!            QQ_Bolton, RR_Bolton, &
!      DGamma, alphaD, alphaB, &
!       dSigmaP, &
!      DD, Mb, Md, K_p, Mcur, GG, KK,  &
!                              dSigma, &
!                             nb, nd, &
!                             Pmin2, &
!                             Ado, &
!      nn, RR   )
!            
!             
!             do ii = 1, size(NextStress)
!                 if (NextStress(ii) /= NextStress(ii)) then
!                     print *, "NextStress(", ii, ") is NaN!"
!                  end if
!              end do
!            
!    
!      end subroutine ForwardEuler
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
!    !--------------------------------------------------------------------------------------
!    !// -------------------------------------------------------------------------------------------------------
!    !/*************************************************************/
!    !// Integrator Constraining Maximum Strain Increment
!    !/*************************************************************/
!      subroutine MaxStrainInc(CurStress, &
!      CurStrain,&
!      CurElasticStrain,&
!      CurAlpha, CurFabric, alpha_in, alpha_in_p,&
!      NextStrain, NextElasticStrain,&
!      NextStress, NextAlpha,&
!      NextFabric, NextL, &
!      NextVoidRatio, &
!      GG, KK, aC,&
!      G0, hpo, PostShake, me2p, &
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, &
!      cz, ceps, nu, Cgd, Ckaf, &
!      QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, &
!      Pmin2, K_p, zpeak, zcum, &
!      Fabric_in,  &
!      MCur,&
!      TolF, &
!      pzp, &
!      alpha_in_true, &
!      nn, RR  ) 
!    
!        implicit none
!        
!        ! input variables
!      real(8), intent(in), dimension(3) :: CurStress 
!
!      real(8), intent(in), dimension(3) :: CurStrain 
!      real(8), intent(in), dimension(3) :: NextStrain
!
!      real(8), intent(in), dimension(3) :: CurElasticStrain 
!      real(8), intent(inout), dimension(3) :: NextElasticStrain 
!        
!      real(8), intent(in), dimension(3) :: CurAlpha 
!      real(8), intent(in), dimension(3) :: CurFabric 
!      real(8), intent(in), dimension(3) :: alpha_in 
!      real(8), intent(in), dimension(3) :: alpha_in_p 
!    
!      ! inout variables
!      real(8), intent(inout), dimension(3) :: NextStress
!    
!      real(8), intent(out), dimension(3) :: NextAlpha
!      real(8), intent(out), dimension(3) :: NextFabric
!      real(8), intent(out) :: NextL
!
!      real(8), intent(out) :: NextVoidRatio
!
!      real(8), intent(out) :: GG
!      real(8), intent(out) :: KK
!	   
!      real(8), intent(out), dimension(3,3) :: aC
!    
!      real(8), intent(in) :: G0       ! primiary parameter 2 ! prop 2
!      real(8), intent(in) :: hpo ! prop 3
!    
!      logical, intent(in) :: PostShake ! prop 4 !real(8)
!      logical, intent(in) :: me2p ! prop 5 !real(8)
!    
!      real(8), intent(in) :: P_atm    ! statev 1
!      real(8), intent(in) :: h0       ! statev 2              ! secondary parameter 1
!      real(8), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
!      real(8), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
!      real(8), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
!      real(8), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
!      real(8), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
!      real(8), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
!      real(8), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
!      real(8), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
!      real(8), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
!      real(8), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11
!    
!      real(8), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12
!    
!      real(8), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
!      real(8), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
!    
!      real(8), intent(in) :: mm         ! statev 17                             ! secondary parameter 15
!    
!      real(8), intent(in) :: InitialVoidRatio
!    
!      real(8), intent(in) :: Pmin       ! statev 25 
!      real(8), intent(in) :: Pmin2      ! statev 26	
!	   
!      real(8), intent(out) :: K_p
!    
!      real(8), intent(in) :: zpeak      ! statev 28
!      real(8), intent(in) :: zcum       ! statev 29
!	
!      real(8), intent(in), dimension(3) :: Fabric_in
!      
!      real(8), intent(inout) :: Mcur       ! statev 36
!    
!      real(8), intent(in) :: TolF
!      real(8), intent(inout) :: pzp
!      
!      real(8), intent(in) :: Fsed_min
!      
!      real(8), intent(inout) :: Md
!      real(8), intent(inout) :: Mb
!      real(8), intent(in) :: Mc
!      
!      real(8), intent(in) :: Cdr
!      
!      real(8), intent(in) :: p_sedo
!    
!      real(8), dimension(3), intent(in) :: Alpha_in_true
!      real(8), dimension(3), intent(out) :: nn
!      real(8), dimension(3), intent(out) :: RR
!      
!      !--------------------------------------------------------------------------
!            
!      ! Local variables
!      real(8) :: maxInc
!      real(8), dimension(3) :: StrainIncStep 
!      integer(4) :: ii, numSteps
!      real(8), dimension(3) :: cStress 
!      real(8), dimension(3) :: cStrain 
!      real(8), dimension(3) :: cAlpha 
!      real(8), dimension(3) :: cFabric 
!      real(8), dimension(3) :: cAlpha_in 
!      real(8), dimension(3) :: cAlpha_in_p 
!      real(8), dimension(3) :: cEStrain
!      real(8), dimension(3) :: nStrain
!      
!      real(8), dimension(3) :: StrainInc
!      
!      real(8) :: maxStrainIncValue
!      
!      real(8), dimension(3) :: TotalStrainIncStep
!      real(8), dimension(3) :: PlasticStrainIncStep
!      
!      !------------------------------------------------------------------------------------
!      maxStrainIncValue = 1.0e-7
!    
!      ! Compute maximum strain increment comparing the components
!      StrainInc = NextStrain - CurStrain
!      maxInc = StrainInc(1)
!      
!      do ii = 2, 3
!          
!        if (abs(StrainInc(ii)) > &
!        abs(maxInc)) then 
!            maxInc = StrainInc(ii)
!        end if
!     
!      end do
!    
!      ! Apply maximum strain increment constraint
!      if (abs(maxInc) > maxStrainIncValue) then
!        numSteps = floor( abs(maxInc) / maxStrainIncValue ) + 1
!        StrainInc = (NextStrain - CurStrain) / real(numSteps)
!        
!        ! floor(A) returns the greatest integer less than or equal to A
!        ! Initialize temporary variables
!        cStress = CurStress
!        cStrain = CurStrain
!        cAlpha = CurAlpha
!        cFabric = CurFabric
!        cAlpha_in = alpha_in
!        cAlpha_in_p = alpha_in_p
!        cEStrain = CurElasticStrain
!    
!        do ii = 1, (numSteps-1)
!            
!            nStrain = cStrain + StrainInc !--> ultimately this will become nStrain = CurStrain + StrainIncrement
!            
!            call ModifiedEuler(cStress, &
!                             cStrain, &
!                             cEStrain,&
!                             cAlpha, &
!                             cFabric, cAlpha_in, cAlpha_in_p,  &
!                             nStrain, NextElasticStrain, NextStress, &
!                             NextAlpha, NextFabric, &
!                             NextL, NextVoidRatio, &
!                             GG, KK,       &
!                             aC,        &
!                             G0, hpo, PostShake, me2p,     &
!                             P_atm, h0, emax, emin, nb, nd, Ado, z_max, &
!          cz, ceps, nu, Cgd, Ckaf, &
!           QQ_Bolton, RR_Bolton, mm, InitialVoidRatio,  &
!          Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, &
!          Pmin2, K_p, zpeak, zcum, &
!          Fabric_in, &
!          MCur, &
!          TolF, &
!          pzp, &
!          alpha_in_true, &
!          nn, RR     )
!    
!            ! Update temporary variables
!            cStress = NextStress
!            cStrain = nStrain
!            cEStrain = NextElasticStrain
!            cAlpha = NextAlpha
!            cFabric = NextFabric
!    
!            end do
!    
!      else
!    
!      call ModifiedEuler(CurStress,&
!      CurStrain,&
!      CurElasticStrain,&
!      CurAlpha, CurFabric, alpha_in, alpha_in_p,  &
!      NextStrain, NextElasticStrain, NextStress, &
!      NextAlpha, NextFabric, NextL, NextVoidRatio, &
!      GG, KK, aC, &
!      G0, hpo, PostShake, me2p, &
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, &
!      cz, ceps, nu, Cgd, Ckaf, &
!      QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, &
!      Pmin2, K_p, zpeak, zcum, &
!      Fabric_in, &
!      MCur,&
!      TolF, &
!      pzp, &
!      alpha_in_true, &
!      nn, RR  )
!    
!      end if
!      
!      end subroutine MaxStrainInc
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
!    !---------------------------------------------------------------------------------------
!    !// ------------------------------------------------------------------------------------
!    !/*************************************************************/
!    !// Modified-Euler Integrator
!    !/*************************************************************/
!      subroutine ModifiedEuler(CurStress, &
!      CurStrain, &
!      CurElasticStrain, &
!      CurAlpha,&
!      CurFabric,&
!      alpha_in, alpha_in_p, &
!      NextStrain, NextElasticStrain, &
!      NextStress, NextAlpha, NextFabric, &
!      NextL, NextVoidRatio, &
!      GG, KK, aC,&
!      G0, hpo, PostShake, me2p,     &
!      P_atm, h0, emax, emin, nb, nd, Ado, z_max, &
!      cz, ceps, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm,&
!      InitialVoidRatio,&
!      Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, &
!      Pmin2, K_p, zpeak, zcum, &
!      Fabric_in, &
!      MCur, &
!      TolF, &
!      pzp, &
!      alpha_in_true, &
!      nn, RR) 
!
!      implicit none
!    
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! input variables
!      real(8), intent(in), dimension(3) :: CurStress 
!      
!      real(8), intent(in), dimension(3) :: CurStrain
!      
!      real(8), intent(in), dimension(3) :: CurElasticStrain
!      real(8), intent(in), dimension(3) :: CurAlpha 
!      real(8), intent(in), dimension(3) :: CurFabric 
!      real(8), intent(in), dimension(3) :: alpha_in 
!      real(8), intent(in), dimension(3) :: alpha_in_p 
!    
!      real(8), intent(in), dimension(3) :: NextStrain
!      real(8), intent(inout), dimension(3) :: NextElasticStrain
!      
!      real(8), intent(inout), dimension(3) :: NextStress
!      
!      real(8), intent(out), dimension(3) :: NextAlpha
!      real(8), intent(out), dimension(3) :: NextFabric
!
!      real(8), intent(out) :: NextL
!      real(8), intent(out) :: NextVoidRatio
!      real(8), intent(out) :: GG
!      real(8), intent(out) :: KK
!      
!      real(8), dimension(3, 3), intent(in) :: aC
!
!      ! input variables
!      real(8), intent(in) :: G0       ! primiary parameter 2 ! prop 2
!      real(8), intent(in) :: hpo ! prop 3
!    
!      logical, intent(in) :: PostShake ! prop 4 !real(8)
!      logical, intent(in) :: me2p ! prop 5 !real(8)
! 
!      real(8), intent(in) :: P_atm    ! statev 1
!      real(8), intent(in) :: h0       ! statev 2              ! secondary parameter 1
!      real(8), intent(in) :: emax     ! statev 3                       ! secondary parameter 2
!      real(8), intent(in) :: emin     ! statev 4                        ! secondary parameter 2
!      real(8), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
!      real(8), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
!      real(8), intent(in) :: Ado      ! statev 7                        ! secondary parameter 5
!      real(8), intent(in) :: z_max    ! statev 8                        ! secondary parameter 6
!      real(8), intent(in) :: cz       ! statev 9                        ! secondary parameter 7
!      real(8), intent(in) :: ceps     ! statev 10                          ! secondary parameter 8
!
!      real(8), intent(in) :: nu       ! statev 12                        ! secondary parameter 10
!      real(8), intent(in) :: Cgd      ! statev 13                        ! secondary parameter 11
!
!      real(8), intent(in) :: Ckaf     ! statev 14                        ! secondary parameter 12
!
!      real(8), intent(in) :: QQ_Bolton  ! statev 15                             ! secondary parameter 13
!      real(8), intent(in) :: RR_Bolton  ! statev 16                             ! secondary parameter 14
!    
!      real(8), intent(in) :: mm         ! statev 17                             ! secondary parameter 15
!
!    
!      real(8), intent(in) :: Pmin       ! statev 25 
!      real(8), intent(in) :: Pmin2      ! statev 26
!
!    
!      real(8), intent(in) :: zpeak      ! statev 28
!      real(8), intent(in) :: zcum       ! statev 29
!    
!    
!      real(8), intent(in), dimension(3) :: Fabric_in
!      
!      real(8), intent(in) :: InitialVoidRatio
!    
!    
!    
!      ! this needs to be intent(in)
!      real(8), dimension(3), intent(in) :: Alpha_in_true ! statev 73:75 !intent(in),
!    
!    ! inout variables
!    
!      real(8), intent(in) :: Mc         ! statev 19                      ! secondary parameter 1
!      real(8), intent(inout) :: Mb         ! statev 20 
!      real(8), intent(inout) :: Md         ! statev 21
!    
!      real(8), intent(in) :: Cdr        ! statev 22                       ! secondary parameter 
!    
!      real(8), intent(in) :: Fsed_min   ! statev 23             ! secondary parameter ! Postshake 
!      real(8), intent(in) :: p_sedo     ! statev 24     ! secondary parameter 1! Postshake                       
!    
!        
!      real(8), intent(inout) :: Mcur       ! statev 36
!
!    ! Local variables
!      real(8) :: CurDr
!      real(8) :: Cka
!      real(8) :: hh
!      real(8) :: pp
!      real(8) :: dVolStrain
!      real(8) :: DD
!      real(8) :: AlphaAlphaBDotN
!    
!      real(8), dimension(3), intent(out) :: nn
!      real(8), dimension(3), intent(out) :: RR
!    
!      real(8), dimension(3) :: alphaB
!    
!      real(8), dimension(3) :: dPStrain
!      
!      
!      real(8), dimension(3) :: dSigma
!      real(8), dimension(3) :: dAlpha
!      real(8), dimension(3) :: dFabric    
!    
!    
!      real(8), dimension(3) :: dSigmaP
!    
!      !real(8) :: pp
!    
!      real(8), intent(inout) :: pzp
!      real(8) :: K_p
!      real(8) :: DGamma
!      integer(4) :: DebugFlag
!    
!      real(8) :: ksi
!      real(8), intent(in) :: TolF
!    
!    ! temporary local variables 
!      real(8) :: temp4 !--> this temporary variable is used to check the denominator of L
!	
!    ! local variables with numeric values 
!      integer(4), dimension(3) :: I1 ! identity matrix
!      real(8) :: two3 ! two thirds
!      real(8) :: small 
!
!      real(8) :: TolE !, intent(in)
!    
!    ! Local variables
!      real(8) :: NextDr
!
!      real(8) :: curStepError
!      real(8) :: qq_ModEuler 
!      real(8) :: stressNorm 
!    
!      real(8) :: TT
!      real(8) :: dT
!      real(8) :: dT_min
!    
!      real(8), dimension(3) :: tmp0 
!      real(8), dimension(3) :: tmp1 
!      real(8), dimension(3) :: tmp2
!      real(8), dimension(3) :: tmp3
!            
!      real(8), dimension(3) :: R1 
!      real(8), dimension(3) :: R2 
!      real(8), dimension(3) :: alphaD
!      real(8), dimension(3) :: dDevStrain 
!      real(8), dimension(3) :: rrr
!      real(8), dimension(3) :: bb
!      real(8), dimension(3) :: dSigma1
!      real(8), dimension(3) :: dSigma2
!      real(8), dimension(3) :: dAlpha1
!      real(8), dimension(3) :: dAlpha2
!      real(8), dimension(3) :: dFabric1 
!      real(8), dimension(3) :: dFabric2 
!      real(8), dimension(3) :: dPStrain1
!      real(8), dimension(3) :: dPStrain2
!    
!      real(8) :: mIIco ! These should be initialized with actual values !two3, 
!
!    
!      real(8), dimension(3) :: nStress
!      real(8), dimension(3) :: nFabric
!      real(8), dimension(3) :: nAlpha
!    
!    
!      real(8), dimension(3) :: alphaD_NextAlpha
!    
!      I1(1) = 1
!      I1(2) = 1
!      I1(3) = 0
!    
!      TT = 0.0 ! start from zero
!      dT = 1.0 ! start at an increment of one
!      dT_min = 1.0e-5 
!      TolE = 1.0e-6
!      
!      small = 1.0e-10
!      
!      two3 = 0.66666666667 !2.0/3.0
!    
!    ! initialize NextElasticStrainIncrement, NextStress, NextAlpha, NextFabric      
!      NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain)
!      
!      NextStress = CurStress
!      NextAlpha = CurAlpha
!      NextFabric = CurFabric
!    
!    ! obtain  KK and GG
!      call GetElasticModuli_(NextStress, zcum, z_max, nu, G0, Md, Mb, &
!        PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, p_sedo, Fsed_min, &
!      me2p, &
!      Mc) ! --> we don't really use Mc here so we can remove 
! 
!    ! Calculate p
!      pp = GetTrace(CurStress)
!      pp = 0.5 * pp
! 
!    ! Check if p < m_Pmin / 5.0 and adjust NextStress if needed
!      if (pp < (Pmin / 5.0) ) then ! --> some sort of a tension cutoff
!        !if (debugFlag) then
!        !    write(*,*) "Tag = ", this%getTag(), " : p < pmin / 5, should not happen"
!        !end if
! 
!        NextStress = GetDevPart(NextStress) + ((Pmin / 5.0) * I1)
!        
!      end if
!    
!    ! Main loop
!      do while (TT < 1.0) ! while less than one (step 7 in page 140 in Sloan et al., 2001)
!        
!        ! Calculate NextVoidRatio and NextDr
!        tmp0 = ( (NextStrain - CurStrain) * TT ) + CurStrain
!        
!        NextVoidRatio = InitialVoidRatio - ( (1.0 + InitialVoidRatio) *&
!                                             GetTrace(tmp0) )
!        
!        NextDr = (emax - NextVoidRatio) / (emax - emin)
!        
!        ! Calculate dVolStrain and dDevStrain
!        tmp0 = NextStrain - CurStrain
!        dVolStrain = dT * GetTrace(tmp0)
!        
!        dDevStrain = (I1 * (-dVolStrain / 3.0))
!        tmp0 = tmp0 * dT
!        dDevStrain = dDevStrain + tmp0
!        
!        pp = 0.5 * GetTrace(NextStress)
!        
!        
!        ! Calculate Delta 1
!        ! Implement the calculation of dSigma1, dAlpha1, dFabric1, and dPStrain1
!        call GetStateDependent(NextStress, NextAlpha, alpha_in, &
!       alpha_in_p, alpha_in_true, NextFabric, Fabric_in, &
!                GG, zcum, zpeak, pzp, Mcur, NextDr, &
!                Mc, nd, nb, Pmin, Pmin2, P_atm, &
!                mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr,&
!               QQ_Bolton, RR_Bolton,&
!               nn, alphaB, alphaD, bb, R1, Mb, Md, DD, K_p, Cka, hh, &
!        AlphaAlphaBDotN, &
!               ksi)
!        
!    
!        
!        rrr = GetDevPart(NextStress) / pp
! 
!        !call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
!        temp4 = K_p + (2 * GG) - (KK * DD * DoubleDot2_2_Contr(nn, rrr)) 
!        
!        if (abs(temp4) < small) then ! <-- to prevent dividing by zero
!    
!            ! neutral loading
!            dSigma1 = 0.0
!            dAlpha1 = 0.0
!            dFabric1 = 0.0
!            dPStrain1 = tmp0
! 
!        else
!            
!            NextL = ( (2 * GG * DoubleDot2_2_Mixed(nn, dDevStrain)) - &
!     (DoubleDot2_2_Contr(nn, rrr) * KK * dVolStrain) ) / temp4
!            
!            if (NextL < 0) then
!                !if (debugFlag) then
!                    
!                !write(*,*) "1 NextL is smaller than 0"
!                !write(*,*) "NextL = ", NextL
!                !end if
!                !call ToContraviant(dDevStrain, ToContraviant_result)
!        
!                dSigma1 = (2 * GG * ToContraviant(dDevStrain)) + &
!            (KK * dVolStrain * I1) ! if L<0 then purely elastic
!                dAlpha1 = 0.0
!                dFabric1 = 0.0
!                dPStrain1 = 0.0
!    
!            else
!                
!                ! plastic stress portion
!                tmp0 = nn * (2.0 * GG)
!                tmp1 = I1 * (KK * DD)
!                tmp1 = tmp1 + tmp0
!                
!                !call Macauley(NextL, Macauley_result)
!                tmp1 = (-Macauley(NextL)) * tmp1 ! plastic stress portion
!                tmp2 = I1 * (KK * dVolStrain) ! volumetric elastic stress
!                
!                !call ToContraviant(dDevStrain, ToContraviant_result)
!                tmp3 = ToContraviant(dDevStrain) * (2.0 * GG) ! shear elastic stress
!                dSigma1 = tmp3 + tmp2 + tmp1
!                         ! dev    vol   plastic
! 
!                ! update fabric				
!                alphaD_NextAlpha = alphaD - NextAlpha
!                
!                if (DoubleDot2_2_Contr(alphaD_NextAlpha, nn) < 0.0) then
!                
!                    ! update fabric terms for dilation 
!                     
!                    ! Equation 57
!                    dFabric1 = (nn*z_max) + NextFabric
!                        
!                    dFabric1 = dFabric1 * &
!         ( (-1.0 * cz / (1.0 + Macauley((zcum / (2.0 * z_max))-1.0))) * &
!                    Macauley(NextL) * MacauleyIndex(-DD)  )
!                    
!                end if
!                
!                dPStrain1 = ToCovariant(R1) * NextL
!                dAlpha1 = bb * (two3 * NextL * hh)
!                
!            end if
!            
!        end if
! 
!        ! Calculate Delta 2
!        ! Implement the calculation of dSigma2, dAlpha2, dFabric2, and dPStrain2
!        tmp0 = NextStress + dSigma1
!        
!        pp = 0.5 * GetTrace(tmp0)
!        
!        if (pp < 0) then
!            
!            if (dT == dT_min) then
!        
!                !if (debugFlag) then
!                !
!                !    write(*,*) "Delta 1: p < 0"
!                !
!                !endif
!                
!                NextElasticStrain = CurElasticStrain + &
!                                     (NextStrain - CurStrain)
!                NextStress = CurStress
!                NextAlpha = CurAlpha
!                NextFabric = CurFabric
!        
!                return
!    
!            endif
!    
!            dT = max(0.1 * dT, dT_min)
!    
!            cycle
! 
!        end if
! 
!        ! tmp1.Zero();  tmp1 += NextAlpha; tmp1 += dAlpha1;  ! tmp1 is NextAlpha + dAlpha1
!        tmp1 = NextAlpha + dAlpha1
!        
!        ! tmp2.Zero();  tmp2 += NextFabric; tmp2 += dFabric1;  ! tmp2 is NextFabric + dFabric1
!        tmp2 = NextFabric + dFabric1
!        
!        call GetStateDependent(tmp0, tmp1, alpha_in, alpha_in_p, &
!       alpha_in_true, tmp2, fabric_in, &
!               GG, zcum, zpeak, pzp, Mcur, NextDr, &
!               Mc, nd, nb, Pmin, Pmin2, P_atm, &
!               mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr,&
!               QQ_Bolton, RR_Bolton,&
!               nn, alphaB, alphaD, bb, R2, &
!               Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, &
!               ksi)
!        
!        ! r = GetDevPart(NextStress + dSigma1) / p;
!        rrr = GetDevPart(tmp0) / pp
!        temp4 = K_p + (2 * GG) - (KK * DD * DoubleDot2_2_Contr(nn, rrr))
!        
!        if (abs(temp4) < small) then !--> prevent division by zero!
!        
!            ! neutral loading
!            dSigma2 = 0
!            dAlpha2 = 0
!            dFabric2 = 0
!            
!            ! dPStrain2 = dDevStrain + dVolStrain * mI1;
!            dPStrain2 = dPStrain1
! 
!        else
!            
!            !NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4
!            NextL = ((2 * GG * DoubleDot2_2_Mixed(nn, dDevStrain)) - &
!            (DoubleDot2_2_Contr(nn, rrr) * KK * dVolStrain)) / temp4
!            DGamma = NextL
!            
!            if (NextL < 0) then
!                !if (debugFlag) then
!                    
!                !write(*,*) "2 NextL is smaller than 0"
!                !write(*,*) "NextL = ", NextL
!                
!                !endif
!        
!                dSigma2 = (2 * GG * ToContraviant(dDevStrain)) + &
!                 (KK * dVolStrain * I1)
!                dAlpha2 = 0
!                dFabric2 = 0
!                dPStrain2 = 0
!            
!            else
!        
!                ! dSigma2 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextL) *
!                !  (2.0 * G * n + K * D * mI1);
!                tmp0 = nn * (2.0 * GG) ! shear part 
!                tmp1 = I1 * (KK * DD) ! volumetric part
!                
!                tmp1 = tmp1 + tmp0
!                tmp1 = tmp1 * (-Macauley(NextL)) ! plastic stress
!                
!                tmp2 = I1 * (KK * dVolStrain) ! elastic volumetric stress
!                
!                
!                tmp3 = ToContraviant(dDevStrain) * (2.0 * GG) ! elastic shear stress
!                dSigma2 = tmp3 + tmp2 + tmp1
!                
!                ! update fabric
!                alphaD_NextAlpha = alphaD - NextAlpha - dAlpha1
!
!                if ( DoubleDot2_2_Contr(alphaD_NextAlpha, nn) < 0.0 ) &
!                                                                 then
!                    
!                    ! Equation 57
!                    ! dFabric2 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + NextFabric + dFabric1);
!                    dFabric2 = nn * z_max
!                    dFabric2 = dFabric2 + NextFabric + dFabric1
!                    
! 
!                    ! I included (1/DD) because that's what is shown in the literature 
!                    dFabric2 = dFabric2 *   &
!        ( -1.0 * (cz / (1.0 + Macauley( (zcum / (2.0*z_max)) - 1.0 ))) &
!                   * Macauley(NextL) * &
!                   MacauleyIndex(-DD) )
! 
!                end if
!                
!                ! dPStrain2 = NextL * mIIco * R2
!                dPStrain2 = ToCovariant(R2) * NextL
!                
!                ! dAlpha2 = two3 * NextL * h * b
!                dAlpha2 = bb * (two3 * NextL * hh)
!            
!            
!            endif
! 
!        endif
! 
!        ! Note:
!        ! -dSigma1 is calculated based on an initial increment of strain and stress
!        ! -dSigma2 is calcualted based on an updated stress state (NextStress+dSigma1), 
!        !  including fabric updates 
!        
!        ! Update nStress, nFabric, and nAlpha
!        nStress = (0.5 * (dSigma1 + dSigma2)) + NextStress
!        nFabric = (0.5 * (dFabric1 + dFabric2)) + NextFabric
!        nAlpha  = (0.5 * (dAlpha1 + dAlpha2)) + NextAlpha
!        
!        
!        
!        ! Check p
!        pp = 0.5 * GetTrace(nStress)
! 
!        ! Perform stress correction if needed
!        if (pp < 0) then
!            
!            if (dT == dT_min) then
!
!                !opserr << "Delta 2: p < 0";
!                ! Implement the stress correction
!                NextElasticStrain = CurElasticStrain + &
!                                     (NextStrain - CurStrain)
!                NextStress = CurStress
!                NextAlpha = CurAlpha
!                NextFabric = CurFabric
!                return
!                
!            end if
!                
!            dT = max(0.1 * dT, dT_min)    
!            cycle
!            
!        end if
! 
!        ! Calculate current step error
!        stressNorm = GetNorm_Contr(NextStress)
!        
!        tmp0 = dSigma2 - dSigma1
!        
!        if (stressNorm < 0.5) then
!            
!            curStepError = GetNorm_Contr(tmp0)
!            
!        else
!            
!            curStepError = GetNorm_Contr(tmp0) / (2 * stressNorm)
! 
!        end if
! 
!        ! Check if the current step error is within tolerance
!        if (curStepError > TolE) then
!            
!            ! Adjust step size
!            qq_ModEuler = max(0.8 * sqrt(TolE / curStepError), 0.1)
!            
!            if (dT == dT_min) then 
!                
!                !// opserr << "reached dT_min\n";
!			  !// NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
!                
!                ! Implement the step size adjustment and stress correction
!                tmp0 = 0.5*(dPStrain1 + dPStrain2)
!                
!                NextElasticStrain = NextElasticStrain - tmp0
!                NextStress = nStress
!                NextAlpha = nAlpha
!                
!                call stress_correction_(NextStress, NextAlpha, alpha_in,&
!                                 alpha_in_p, CurFabric, NextVoidRatio, &
!                                 alpha_in_true, &
!                                 Fabric_in, CurStress, Pmin, P_atm, &
!                             TolF, Mc, emax, emin, zcum, zpeak, z_max, &
!                                 pzp, mm, h0, hpo, Cdr,&
!                                 Ceps, Ckaf,        &
!                                 QQ_Bolton, RR_Bolton, &
!                                 DGamma, alphaD, alphaB, &
!                                 dSigmaP, &
!                                 DD, Mb, Md, K_p, Mcur, GG, KK,  &
!                                 dSigma, &
!                                 nb, nd, &
!                                 Pmin2, &
!                                 Ado, &
!                                 nn, RR   )
!                
!                TT = TT + dT
!                
!            end if
! 
!            dT = max(qq_ModEuler * dT, dT_min)
!        
!        else
!            
!            ! Implement the step size adjustment and update NextElasticStrain, NextStress, NextAlpha, and NextFabric
!            NextElasticStrain = NextElasticStrain - &
!           (0.5 * (dPStrain1 + dPStrain2))
!            NextStress = nStress
!            NextAlpha = nAlpha
!            NextFabric = nFabric
!            
!            call stress_correction_(NextStress, NextAlpha, alpha_in,&
!                                 alpha_in_p, CurFabric, NextVoidRatio, &
!                                 alpha_in_true, &
!                                 Fabric_in, CurStress, Pmin, P_atm,  &
!                             TolF, Mc, emax, emin, zcum, zpeak, z_max, &
!                                 pzp, mm, h0, hpo, Cdr,&
!                                 Ceps, Ckaf,        &
!                                 QQ_Bolton, RR_Bolton, &
!                                 DGamma, alphaD, alphaB, &
!                                 dSigmaP, &
!                                 DD, Mb, Md, K_p, Mcur, GG, KK,  &
!                                 dSigma, &
!                                 nb, nd, &
!                                 Pmin2, &
!                                 Ado, &
!                                 nn, RR   )
!            
!            TT = TT + dT
!            
!            if (curStepError == 0) then 
!                curStepError = 1e-15
!            end if 
!            
!            
!            qq_ModEuler = max(0.8 * sqrt(TolE / curStepError), 0.5)
!            dT = max(qq_ModEuler * dT, dT_min)
!            dT = min(dT, 1.0 - TT)
!        
!        end if
!    
!        end do
!    
! 
!        return
!        
!        
!        end subroutine ModifiedEuler
!    
!    
!    
!    
!    
!    
!    
!    !--------------------------------------------------------------
!    !// -------------------------------------------------------------------------------------------------------
!    !/*************************************************************/
!    !// Runge-Kutta Integrator
!    !/*************************************************************/
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
!    !--------------------------------------------------------------
!    !/*************************************************************/
!    !//            Pegasus Iterations                             //
!    !/*************************************************************/
!      function IntersectionFactor(CurStress, CurStrain, NextStrain, &
!      CurAlpha, Ce, a0, a1, TolF, mm) result(ElasticRatio)
!    
!      implicit none
!    
!    ! inputs
!      real(8), dimension(3), intent(in) :: CurStress
!      real(8), dimension(3), intent(in) :: CurAlpha
!      real(8), dimension(3) :: StrainIncrement !, intent(in)
!      real(8), dimension(3), intent(in) :: CurStrain
!      real(8), dimension(3), intent(in) :: NextStrain
!    
!      real(8), dimension(3,3), intent(in) :: Ce
!    
!      real(8), intent(in) :: TolF
!      
!      real(8), intent(in) :: mm
!    
!    ! outputs
!      real(8) :: ElasticRatio
!    
!    ! inout
!      real(8), intent(inout) :: a0 
!      real(8), intent(inout) :: a1
!    
!    !local variables
!      real(8) :: ff 
!      real(8) :: f0 
!      real(8) :: f1
!    
!      real(8), dimension(3) :: NextElasticStress
!    
!      real(8), dimension(3) :: dSigma 
!      real(8), dimension(3) :: dSigma0 
!      real(8), dimension(3) :: dSigma1
!    
!      integer(4) :: ii
!    
!    ! Constants
!      integer(4) :: maxIter! = 10
!      real(8) :: small! = 1.0e-10
!    
!    ! establish numerical values
!      maxIter = 100
!      small = 1.0e-10
!      
!      StrainIncrement = NextStrain - CurStrain
!    
!    ! Calculate dSigma0 and f0
!      dSigma0 = matmul(Ce, (a0 * StrainIncrement)) ! a0 stress increment
!    
!    ! temporary new stress for a0
!      NextElasticStress = CurStress + dSigma0
!    
!    ! evaluate yield function for this a0
!      f0 = GetFYieldFunction(NextElasticStress, CurAlpha, mm) ! a0 yield function evaluation
!    
!    ! calculate dSigma1 and f1
!      dSigma1 =  matmul(Ce, (a1 * StrainIncrement)) ! a1 stress increment
!
!    ! temporary new stress for a1
!      NextElasticStress = CurStress + dSigma1 
!    
!    ! evaluate yield function for this a1
!      f1 = GetFYieldFunction(NextElasticStress, CurAlpha, mm) ! a1 yield function evaluation
!
!    ! Iterate using the Pegasus method
!      do ii = 1, maxIter
!        
!        ! Calculate new a
!        ElasticRatio = a1 - (f1 * (a1 - a0) / (f1 - f0) ) !--> interpolation
!        !dSigma = mCe * (aa * strainInc)
!        dSigma = matmul(Ce , (ElasticRatio * StrainIncrement))
!        
!        NextElasticStress = CurStress + dSigma !--> new stress
!        
!        
!        ! evaluate yield function
!        ff = GetFYieldFunction(NextElasticStress, CurAlpha, mm)
!
!        ! Check convergence
!        if (abs(ff) < TolF) then
!            exit
!        end if
!
!        ! Update bounds for the next iteration
!        if ( (ff * f0) < 0) then
!            ! upper bound
!            a1 = ElasticRatio
!            f1 = ff
!        else
!            ! lower bound and gradient
!            f1 = f1 * f0 / (f0 + ff)
!            a0 = ElasticRatio
!            f0 = ff
!        end if
!
!        ! Check for non-convergence
!        if (ii == maxIter) then
!            !if (debugFlag) then
!            !    write(*, *) "Didn't find alpha!"
!            !endif
!            ElasticRatio = 0.0
!            ! SOME THING WRONG HAPPENED IF WE ARE HERE
!            exit
!        end if
!        
!      end do
!
!    ! Ensure a is within bounds
!      if (ElasticRatio > (1.0 - small)) ElasticRatio = 1.0
!      if (ElasticRatio < small) ElasticRatio = 0.0
!    
!
!      end function IntersectionFactor
!    
!    !----------------------------------------------------------------
!    !/*************************************************************/
!    !//      Pegasus Iterations  (ElastoPlastic Unloading)        //
!    !/*************************************************************/
!    
!      function IntersectionFactor_Unloading(CurStress, &
!      CurStrain, NextStrain, CurAlpha, Ce, TolF, & !StrainIncrement
!      mm) result(ElasticRatio) !CurStrain, NextStrain,
!    
!      implicit none
!    
!    ! inputs
!      real(8), dimension(3), intent(in) :: CurStress
!      real(8), dimension(3) :: StrainIncrement
!      
!      real(8), dimension(3), intent(in) :: CurStrain
!      real(8), dimension(3), intent(in) :: NextStrain
!
!      real(8), dimension(3), intent(in) :: CurAlpha
!      real(8), dimension(3,3), intent(in) :: Ce
!      real(8), intent(in) :: TolF ! = 1.0e-8 to 1.0e-10
!    
!      real(8), intent(in) :: mm
!      
!    ! outputs
!      real(8) :: ElasticRatio
!    
!    ! Local variables
!      real(8) :: a0
!      real(8) :: a1 
!      real(8) :: da
!    
!      real(8) :: ff
!      real(8) :: f0
!      real(8) :: f1
!      real(8) :: fs
!    
!      real(8), dimension(3) :: dSigma
!      real(8), dimension(3) :: dSigma0 
!      real(8), dimension(3) :: dSigma1 
!      real(8), dimension(3) :: NextElasticStressPoint !tmp
!      logical :: flag = .false.
!
!      integer(4) :: ii
!      integer(4) :: kk
!      integer(4) :: nSub
!      integer(4) :: maxIter
!    
!    ! initialize variables
!      nSub = 100
!      maxIter = 100
!    
!      ElasticRatio = 0.0
!      a0 = 0.0 
!      a1 = 1.0
!    
!      
!      StrainIncrement = NextStrain - CurStrain 
!    
!    ! Compute initial value of f0
!      f0 = GetFYieldFunction(CurStress, CurAlpha, mm)
!      fs = f0
!    
!    ! Compute initial value of dSigma
!      dSigma = DoubleDot4_2(Ce, StrainIncrement)
!
!    ! Pegasus iterations loop
!      do ii = 1, maxIter
!        
!        ! Compute increment factor
!        da = (a1 - a0) / real(nSub)
!        
!        
!        do kk = 1, nSub-1 ! Iterate over subintervals
!            
!            ! evaluate elastic ratio
!            ElasticRatio = a0 + da!(da * real(kk))
!            
!            ! Compute stress increment
!            NextElasticStressPoint = (ElasticRatio * dSigma) + CurStress
!            
!            ! Compute value of f
!            ff = GetFYieldFunction(NextElasticStressPoint, CurAlpha, &
!                                                              mm)
!
!            ! Check if f exceeds tolerance
!            if (ff > TolF) then
!                
!                a1 = ElasticRatio
!                
!                if (f0 < -TolF) then
!                    f1 = ff
!                    flag = .true.
!                    exit
!                else
!                    a0 = 0.0
!                    f0 = fs
!                    exit
!                end if
!            
!            else
!            
!                a0 = ElasticRatio
!                f0 = ff
!            
!            end if
!            
!        end do ! Iterate over subintervals
!
!        ! Check termination conditions
!        if (flag) exit
!        
!        if (ii == maxIter) then
!            
!            !if (debugFlag) then
!            !    write(*, *) "Didn't find alpha! - Unloading, a0 = ", a0, ", a1 = ", a1
!            !endif
!            ElasticRatio = 0.0
!            return
!            
!        end if
!        
!    
!      end do
!
!    ! Output the result
!    !if (debugFlag) then
!    !    write(*, *) "Found alpha - Unloading, a0 = ", a0, ", a1 = ", a1
!    !endif
!      
!    ! we update a0 and a1, then plug it in IntersectionFactor to get aa
!      ElasticRatio = &
!      IntersectionFactor(CurStress, CurStrain, NextStrain, &
!      CurAlpha, Ce, &!StrainIncrement
!      a0, a1, TolF, mm) !CurStrain, NextStrain,
!    
!
!    
!      end function IntersectionFactor_Unloading
!    
!    
!    
!    
!    
!    
!    
!    
!    !--------------------------------------------------------
!    !/*************************************************************/
!    !//            Stress Correction                              //
!    !/*************************************************************/      
!      subroutine stress_correction_(NextStress, NextAlpha, alpha_in,&
!      alpha_in_p, CurFabric, NextVoidRatio,&
!      alpha_in_true,&
!      Fabric_in,  CurStress, Pmin, P_atm, &
!      TolF, Mc, emax, emin, zcum, zpeak, z_max, pzp, mm, h0, hpo, Cdr,&
!      Ceps, Ckaf, &
!      QQ_Bolton, RR_Bolton, &
!      DGamma, alphaD, alphaB, dSigmaP, DD, Mb, Md, K_p,  &
!      Mcur, GG, KK,  &
!      dSigma, &
!      nb, nd, &
!      Pmin2, Ado, &
!      nn, RR)
!    
!      implicit none
!    
!    ! in variables 
!      real(8), dimension(3), intent(in) :: alpha_in 
!      real(8), dimension(3), intent(in) :: alpha_in_p
!      real(8), dimension(3), intent(in) :: alpha_in_true
!      real(8), dimension(3), intent(in) :: Fabric_in 
!      real(8), dimension(3), intent(in) :: CurFabric
!      real(8), dimension(3), intent(in) :: CurStress
!      real(8), intent(in) :: Pmin 
!      real(8), intent(in) :: P_atm
!    
!      real(8), intent(in) :: TolF
!    
!      real(8), intent(in) :: Mc
!      real(8), intent(in) :: emax
!      real(8), intent(in) :: emin
!    
!      real(8), intent(in) :: zcum 
!      real(8), intent(in) :: zpeak
!      real(8), intent(in) :: z_max
!    
!      real(8), intent(in) :: pzp
!    
!      real(8), intent(in) :: mm
!    
!      real(8), intent(in) :: h0
!    
!      real(8), intent(in) :: hpo
!    
!      real(8), intent(in) :: QQ_Bolton
!      real(8), intent(in) :: RR_Bolton
!    
!      real(8), intent(in) :: DGamma ! this is practically the L value to calculate the plastic stresses !
!    
!    
!      real(8), intent(in) :: nb       ! statev 5                        ! secondary parameter 3
!      real(8), intent(in) :: nd       ! statev 6                        ! secondary parameter 4
!    
!      real(8), intent(in) :: Pmin2      ! statev 26
!
!      real(8), intent(in) :: Ado
!    
!    ! output variable
!      real(8), dimension(3), intent(out) :: alphaD
!      real(8), dimension(3), intent(out) :: alphaB
!      
!      real(8), dimension(3), intent(out) :: dSigmaP
!    
!      real(8), intent(out) :: DD
!      real(8), intent(out) :: Mb 
!      real(8), intent(out) :: Md 
!      real(8), intent(out) :: K_p 
!    
!    
!    ! inout variables
!      real(8), intent(inout) :: NextVoidRatio
!    
!      real(8), intent(inout) :: Mcur
!    
!      real(8), intent(inout) :: GG
!      real(8), intent(inout) :: KK !--> probably inout
!    
!      real(8), dimension(3), intent(inout) :: NextStress
!      real(8), dimension(3), intent(inout) :: NextAlpha
!    
!      real(8), dimension(3), intent(inout) :: dSigma
!    
!    
!    ! Local variables
!      real(8), dimension(3) :: dfrOverdSigma
!      real(8), dimension(3) :: dfrOverdAlpha
!      real(8), dimension(3), intent(out) :: nn
!      real(8), dimension(3), intent(out) :: RR
!    
!      real(8), dimension(3) :: bb
!      real(8), dimension(3) :: aBar
!      real(8), dimension(3) :: rrr
!    
!      real(8), dimension(3) :: tmp0
!      real(8), dimension(3) :: tmp1
!    
!      real(8), dimension(3) :: nAlpha
!      real(8), dimension(3) :: nStress
!    
!      real(8) :: ksi
!    
!      real(8) :: CurDr
!    
!      real(8), intent(in) :: Ckaf
!      real(8), intent(in) :: Cdr
!      real(8), intent(in) :: Ceps
!    
!      real(8) :: lambda
!      real(8) :: Cka
!      real(8) :: hh
!      real(8) :: pp
!      real(8) :: fr
!      real(8) :: AlphaAlphaBDotN
!    
!    
!      real(8) :: maxIter
!    
!      real(8) :: fr_old
!      real(8) :: alpha_up
!      real(8) :: alpha_mid
!      real(8) :: alpha_down
!    
!      integer(4) :: jj
!      integer(4) :: ii 
!    
!      real(8) :: two3
!    
!      real(8), dimension(3) :: I1
!    
!      real(8), dimension(3,3) :: Ce
!    
!    ! Constants
!      I1(1) = 1
!      I1(2) = 1
!      I1(3) = 0
!    
!      maxIter = 200
!      two3 = 0.6666666666667
!
!    
!    ! Compute p
!      pp = GetTrace(NextStress)
!      pp = 0.5 * pp !GetTrace(NextStress, GetTrace_result)
!
!    ! Check tension cutoff
!      if (pp < (Pmin / 5.0) ) then
!        
!        ! I added this so that we can have the correct considertion of the two shear components
!        !call ToCovariant(NextStress, ToCovariant_result)
!        
!        fr = GetFYieldFunction(NextStress, NextAlpha, mm)
!        
!        if (fr < TolF) then
!            ! Stress state inside yield surface
!            NextStress = NextStress + ( ((Pmin / 5.0) - pp) * I1)
!        else
!            ! Stress state outside yield surface
!            NextStress = (Pmin / 5.0) * I1
!            NextStress(3) = 0.8 * Mc * Pmin / 5.0
!            NextAlpha = 0.0
!            NextAlpha(3) = 0.8 * Mc
!            return
!        endif
!        
!      else
!        
!        ! check the initial state of stress: is it in the yield surface?
!        fr = GetFYieldFunction(NextStress, NextAlpha, mm) !this works because the strain 
!        
!        ! if you are inside the yield surface then RETURN: You are done!
!        if (fr < TolF) then 
!        
!            RETURN
!        else 
!            
!        ! if you are outside the yield surface: then you need to calculate residuals:
!        ! these represent the difference between the actual stress state and the stress state predicted by the plasticity model
!        ! Compute CurDr
!        ! use next stress and next alpha
!        CurDr = (emax - NextVoidRatio) / (emax - emin)
!        nStress = NextStress
!        nAlpha = NextAlpha
!        
!        ! Stress correction loop
!        do ii = 1, maxIter
!            
!            ! Compute r
!            rrr = GetDevPart(nStress)
!            rrr = rrr / pp  ! calculate the stress ratio for the next stress state 
!            
!            ! get the state dependent variables which are a function of nStress and nAlpha
!            call GetStateDependent(nStress, nAlpha, alpha_in, &
!      alpha_in_p, alpha_in_true, CurFabric, Fabric_in, &
!               GG, zcum, zpeak, pzp, Mcur, CurDr, &
!               Mc, nd, nb, Pmin, Pmin2, P_atm, &
!               mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
!               QQ_Bolton, RR_Bolton, &
!               nn, alphaB, alphaD, bb, RR, &
!               Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, &
!               ksi)           
!            
!            ! Compute the increment in stress caused by plastic deformation: It is calculated as the product of the 
!            ! material stiffness matrix (aC) get the 3x3 2D stiffness matrix 
!            Ce = GetStiffness(KK, GG) 
!
!            ! Compute dSigmaP
!            !--> RR is the direction of the plastic strain
!            !--> RR' is the deviatoric component of RR
!            
!            ! dSigmaP = aC * PLASTIC STRAIN = aC : [<L>R] !  WE ARE CONSIDERING ENTIRE R
!            dSigmaP = DoubleDot4_2(Ce, DGamma * ToCovariant(RR)) 
!            !                   3x3 matrix          scalar * RR
!            
!            ! Compute aBar --> dAlpha
!            aBar = two3 * hh * bb 
!            
!            !call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
!            
!            ! Compute dfrOverdSigma
!            dfrOverdSigma = (-0.5 * DoubleDot2_2_Contr(nn, rrr) * I1)  &
!                              + nn
!            ! Compute dfrOverdAlpha
!            dfrOverdAlpha = -pp * nn
!            
!            ! Compute lambda
!            lambda = fr / &
!      (DoubleDot2_2_Contr(dfrOverdSigma, dSigmaP) - &
!           DoubleDot2_2_Contr(dfrOverdAlpha, aBar))
!            
!            ! Update NextStress and NextAlpha
!            ! This is consistency condition and flow rule
!            tmp0 = - (lambda * dSigmaP) + nStress ! next stress
!            tmp1 = (lambda * aBar) + nAlpha ! next alpha
!            
!
!            
!            ! evaluate the yield function
!            if (abs(GetFYieldFunction(tmp0, tmp1, mm)) < abs(fr)) then ! first if statement in the stress correction loop
!                !nStress = tmp0
!                nStress = nStress - (lambda*dSigmaP)
!                !nAlpha = tmp1
!                nAlpha = nAlpha + (lambda*aBar)
!            else
!              
!                lambda = &
!               fr / DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma) 
!                
!                nStress = nStress - (lambda * dfrOverdSigma)
!            end if
!            
!            
!            ! Update fr
!            fr = GetFYieldFunction(nStress, nAlpha, mm)
!            
!             if (abs(fr) < TolF) then ! SECOND IF in the stress correction loop --> go to step 8
!                 ! corrected state
!                NextStress = nStress ! update NextStress
!                NextAlpha = nAlpha ! update NextAlpha
!                return
!             end if
!             
!            ! Update p
!            pp = GetTrace(nStress)
!            pp = max(0.5 * pp, Pmin)
!            
!        end do
!        
!        ! this is an additional layer of correction
!        ! Search for the optimum alpha_mid
!        dSigma = NextStress - CurStress
!        alpha_up = 1.0
!        alpha_mid = 0.5
!        alpha_down = 0.0
!        tmp0 = (alpha_mid * dSigma) + CurStress
!        
!        fr_old = GetFYieldFunction(tmp0, NextAlpha, mm)
!        
!        do jj = 1, maxIter
!            
!            if (fr_old < 0.0) then
!                alpha_down = alpha_mid
!                alpha_mid = 0.5 * (alpha_up + alpha_mid)
!            else
!                alpha_up = alpha_mid
!                alpha_mid = 0.5 * (alpha_down + alpha_mid)
!            endif
!            
!            tmp0 = (alpha_mid * dSigma) + CurStress!mSigma
!            
!            ! Correct the stress
!            ! I added this so that we can have the correct considertion of the two shear components
!            !call ToCovariant(tmp0, ToCovariant_result)
!            
!            !GetFYieldFunction_result = 0.0
!            fr_old = GetFYieldFunction(tmp0, NextAlpha, mm)
!                        
!            if (abs(fr_old) < TolF) then
!                NextStress = (alpha_mid * dSigma) + CurStress!tmp0
!                exit
!            end if
!            
!        end do
!        
!        end if 
!        
!      endif
!    
!
!      end subroutine Stress_Correction_
!    
!    
!    
!    
!    
!    
!    
!    
!    
!    !---------------------------------------------------------------
!    
!    !  subroutine Stress_Correction(NextStress, NextAlpha, dAlpha, mm, !--> not used
!    ! & RR, nn, rrr, TolF) --> Potts and Zdrakovic
!    !
!    !  implicit none
!    !
!    !! Arguments
!    !  real(8), dimension(3), intent(inout) :: NextStress 
!    !  real(8), dimension(3), intent(inout) :: NextAlpha 
!    !  real(8), dimension(3), intent(inout) :: dAlpha 
!    !  real(8), dimension(3), intent(inout) :: rrr
!    !
!    !  real(8), dimension(3), intent(in) :: RR 
!    !  real(8), dimension(3), intent(in) :: nn
!    !
!    !  real(8), intent(in) :: mm
!    !
!    !  real(8), intent(in) :: TolF
!    !
!    !! Local variables
!    !  real(8), dimension(3) :: dfrOverdSigma
!    !  real(8) :: lambda
!    !  integer(4) :: maxIter, ii
!    !  real(8) :: ff
!    !
!    !! Local variables 
!    !  real(8) :: DoubleDot2_2_Contr_result
!    !  real(8) :: GetFYieldFunction_result
!    !
!    !  real(8), dimension(3) :: ToCovariant_result
!    !
!    !  real(8), dimension(3) :: I1
!    !  real(8) :: root12
!    !
!    !! Constants
!    !  maxIter = 1000
!    !!mTolF = 1.0e-5
!    !
!    !  I1(1) = 1
!    !  I1(2) = 1
!    !  I1(3) = 0
!    !
!    !  root12 = 0.70710678118
!    !! Correct the stress
!    !! I added this so that we can have the correct considertion of the two shear components
!    !  !ToCovariant_result = 0.0
!    !  !!call ToCovariant(NextStress, ToCovariant_result)
!    !  ff = GetFYieldFunction(NextStress, NextAlpha, mm)
!    !
!    !
!    !! Check if correction is needed
!    !  if (ff < TolF) return
!    !
!    !! Stress correction loop
!    !  do ii = 1, maxIter
!    !    ! Compute dfrOverdSigma
!    !    !call DoubleDot2_2_Contr(nn, rrr, DoubleDot2_2_Contr_result)
!    !    dfrOverdSigma = nn - (0.5 * DoubleDot2_2_Contr(nn, rrr) * I1)
!    !
!    !    ! Compute lambda
!    !    !call DoubleDot2_2_Contr(dfrOverdSigma, RR, 
!    !  !& DoubleDot2_2_Contr_result)
!    !    lambda = ff / DoubleDot2_2_Contr(dfrOverdSigma, RR) !DoubleDot2_2_Contr(dfrOverdSigma, R)
!    !
!    !    ! Update NextStress and NextAlpha
!    !    NextStress = NextStress - (RR * lambda)
!    !    NextAlpha = NextAlpha - (dAlpha * lambda)
!    !    
!    !    ! Correct the stress
!    !    ! I added this so that we can have the correct considertion of the two shear components
!    !    !ToCovariant_result = 0.0
!    !    !call ToCovariant(NextStress, ToCovariant_result)
!    !    
!    !
!    !    ! Recompute f
!    !    ff = GetFYieldFunction(NextStress, NextAlpha, mm)
!    !
!    !    ! Check convergence
!    !    if (abs(ff) < TolF) exit
!    !
!    !    ! Check maximum iterations
!    !    !if (ii == maxIter) then
!    !    !    if (debugFlag) print *, "Still outside with f = ", ff
!    !    !endif
!    !  end do
!    !
!    !  end subroutine Stress_Correction
!    
!     
!    !---------------------------------------------------------------
!    !/*************************************************************/
!    !/*************************************************************/
!    !//            MATERIAL SPECIFIC METHODS                      //
!    !/*************************************************************/
!    !/*************************************************************/
!    !// Macauley() -------------------------------------------------
!    
!      function Macauley(xx) result(xx_Macauley) 
!    !--------------------------------------------
!    ! To do Macauley bracket. If it is negative, set to zero.
!    !--------------------------------------------
!      implicit none 
!      real(8), intent(in) :: xx
!      real(8) :: xx_Macauley !, intent(out)
!    
!    !// Macauley bracket
!      if (xx>0) then 
!        xx_Macauley = xx
!      else 
!        xx_Macauley = 0
!      end if 
!    
!      end function Macauley
!    
!    
!    !---------------------------------------------------
!    
!      function MacauleyIndex(xx) result(MacauleyIndex_result) 
!    !--------------------------------------------
!    ! To do Macauley bracket. If it is negative, set to zero.
!    !--------------------------------------------
!      implicit none 
!      real(8), intent(in) :: xx
!      real(8) :: MacauleyIndex_result !, intent(out)
!    
!    !// Macauley bracket
!      if (xx>0) then 
!        MacauleyIndex_result = 1
!      else 
!        MacauleyIndex_result = 0
!      end if 
!    
!      end function MacauleyIndex
!    
!    
!    
!    
!   
!    
!    !-----------------------------------------------------
!    
!      function GetFYieldFunction(nStress, nAlpha, mm) result(ff) 
!    !-----------------------------------------------------
!    ! To evaluate the yield function. 
!    !-----------------------------------------------------
!      implicit none 
!    
!    ! input parameters 
!      real(8), intent(in), dimension(3) :: nAlpha
!      real(8), intent(in), dimension(3) :: nStress
!    
!    ! local parameters
!      real(8), dimension(3) :: ss ! deviator stress
!      real(8) :: pp
!    
!    ! output parameters
!      real(8) :: ff ! scalar
!    
!    ! Local parameters 
!      real(8) :: GetTrace_result
!      real(8) :: GetNorm_Contr_result
!      real(8), dimension(3) :: GetDevPart_result
!    
!      real(8) :: root12
!      real(8), intent(in) :: mm
!    
!      real(8) :: m_pmin
!    
!    
!      root12 = 0.70710678118
!    
!    !// PM4Sand's yield function
!	!// s = nStress - p*I
!	!//   = [s_xx	s_xy ;	=	[Sigma_xx-p		Sigma_xy   ;
!	!//		s_xy	s_yy];		 Sigma_xy		Sigma_yy-p];
!	!// f = [(s - p*alpha):(s - p*alpha)]^0.5 - sqrt(1/2)*p*m
!    
!    ! get mean effective stress 
!      pp = GetTrace(nStress)
!      pp = 0.5 * pp
!    
!    ! get deviatoric part of the stress by subtracting the mean eff stress from the xx and yy components
!      ss = GetDevPart(nStress)
!      !ss = GetDevPart_result !GetDevPart(nStress)
!    
!    ! subtract the stress ratio
!      ss = ss - (pp*nAlpha)
!    
!    ! dot ss by ss
!      !call GetNorm_Contr(ss, GetNorm_Contr_result)
!    
!    ! evaluate the yield function ff
!      ff = GetNorm_Contr(ss) - (root12 * mm * pp)
!    
!      end function GetFYieldFunction
!    
!    
!    !----------------------------------------------------
!    
!      function GetKsi(Dr, pp, RR_Bolton, QQ_Bolton, Pmin, P_atm) &
!        result(ksi)
!    !-------------------------------------------
!    ! To calculate the relative state parameter 
!    ! based on Bolton (1986) as used by BZ14. .
!    !-------------------------------------------
!      implicit none
!    
!    ! input parameters
!      real(8), intent(in) :: Dr
!      real(8), intent(in) :: pp
!    
!      real(8), intent(in)  :: RR_Bolton ! secondary parameter to find Dr,cs
!      real(8), intent(in)  :: QQ_Bolton ! secondary parameter to find Dr,cs
!    
!      real(8), intent(in) :: Pmin
!      real(8), intent(in)  :: P_atm ! atmospheric pressure
!    
!    ! output parameters 
!      real(8) :: ksi !, intent(out)
!    
!    ! inout parameters
!    
!    ! local parameters
!      real(8) :: pn
!    
!    
!    ! store initial 
!      pn = pp
!    
!    ! cut off
!      if (pn<=Pmin) then 
!        pn = Pmin
!      end if 
!    
!    !// find ksi_R
!    !// calculating the relative density based on Bolton (1986)
!      ksi = (RR_Bolton / (QQ_Bolton - log(100.0 * pn / P_atm)) ) - Dr
!    
!      end function GetKsi
!    
!    
!    
!    
!    !-------------------------------------------------
!      subroutine GetElasticModuli_(sigma, zcum, z_max, nu, G0, &
!      Md, Mb, PostShake, Pmin, P_atm, KK, GG, Mcur, &
!      Cgd, p_sedo, Fsed_min, me2p, &
!      Mc)
!    !--------------------------------------------
!    ! To calculate the elastic modulus of the soil
!    !--------------------------------------------
!      implicit none 
!    
!    ! input variables 
!      real(8), intent(in), dimension(3) :: sigma(3)
!      real(8), intent(in) :: zcum
!      real(8), intent(in) :: z_max
!      real(8), intent(in) :: nu 
!      real(8), intent(in) :: G0
!      real(8), intent(in) :: Md
!      real(8), intent(in) :: Mb
!      logical, intent(in) :: PostShake !real(8)
!      real(8), intent(in) :: Pmin
!      real(8), intent(in) :: P_atm
!      real(8), intent(in) :: Cgd
!    
!      logical, intent(in) :: me2p !real(8)
!    
!      real(8), intent(in) :: p_sedo
!      real(8), intent(in) :: Fsed_min
!      real(8), intent(in) :: Mc
!    
!    ! output variables 
!      real(8), intent(out) :: KK
!      real(8), intent(out) :: GG
!      real(8), intent(out) :: Mcur
!    
!    ! local variables
!      real(8) :: pp 
!      real(8) :: p_sed
!      real(8) :: F_sed
!      
!      real(8) :: pn
!      real(8) :: qn
!      real(8) :: Csr0
!      real(8) :: Csr
!      real(8) :: msr
!    
!      real(8) :: two3
!    
!      real(8) :: temp
!    
!      
!      real(8) :: small, temp_min
!    
!    ! establish numeric values associated with this subroutine
!      two3 = 0.666666666667
!      msr = 4.0
!      Csr0 = 0.5
!      small = 1e-6
!    
!      pn = GetTrace(sigma)
!      pn = 0.5 * pn
!    
!      if (pn <= Pmin) then 
!        pn = Pmin 
!      end if 
!    
!      qn = 2.0 * sqrt( ( (0.5*(sigma(1) - sigma(2)))**2.0) &
!       + sigma(3)**2.0)
!    
!      Mcur = qn/pn
! 
!      !Csr = 1 - ( Csr0 * min(1.0, (Mc/Mb)**msr) ) 
!      temp_min = MIN(1.0,((Mcur/Mb)**msr))  
!      Csr = 1 - ( Csr0 *  temp_min) !((Mcur/Mb)**msr) ) 
!
!      !if (Csr>1) Csr=1
!      
!    
!       temp = zcum / z_max
!    
!    
!    
!      if (me2p == .FALSE.) then 
!    
!        GG = G0 * P_atm !--> we never use this --> you can use this in stress initialization
!    
!    
!      else ! me2p == .true.
!          
!        GG = G0 * P_atm * sqrt(pn / P_atm) * Csr * &
!      (1.0 + temp) / (1.0 + (temp * Cgd))     
!        
!        
!        if (PostShake==.TRUE.) then ! consider changing this to logical (love is never logical...)
!			
!            !// reduce elastic shear modulus for post shaking consolidation
!              pp = GetTrace(sigma)
!			pp = 0.5 * pp
!            
!            !call Macauley(1 - (Mcur/Md), Macauley_result)
!			p_sed = p_sedo * (zcum / (zcum + z_max)) * &
!             Macauley(1 - (Mcur/Md))**0.25 !// Equation 88
!			F_sed = min(Fsed_min + (1 - Fsed_min) * &
!      ((pp / 20.0 / (p_sed + small))**2.0) , 1.0) !// Equation 87 !--> shouldn't there be a square on here 
!			GG = GG * F_sed !// Equation 85
!            
!		end if 
!        
!      end if 
!    
!
!      ! Bulk modulus
!      KK = two3 * ( (1+nu)/(1 - (2*nu)) ) * GG
!    
!      end subroutine GetElasticModuli_
!    
!    
!    !---------------------------------------------------
!    
!    
!      subroutine GetElasticModuli(sigma, KK, GG, G0, nu, P_atm, Pmin, &
!      me2p) !--> Fabric independent
!    !&  result(GG, KK)
!    !---------------------------------------------------
!    ! To get the elastic parameters
!    !---------------------------------------------------
!      implicit none 
!    
!    ! input
!      real(8), intent(in), dimension(3) :: sigma
!      real(8), intent(in) :: G0    
!      real(8), intent(in) :: nu
!      real(8), intent(in) :: P_atm
!      real(8), intent(in) :: Pmin
!    
!    !logical, intent(in) :: me2p
!    
!    ! output
!      real(8), intent(out) :: GG
!      real(8), intent(out) :: KK
!    
!    
!    ! Local variables 
!      real(8) :: pn
!      real(8) :: two3
!      real(8) :: GetTrace_sigma
!      
!      logical :: me2p
!    
!    ! numeric values in this subroutine
!      two3 = 2.0/3.0
!    
!    ! get mean effective stress
!      !call GetTrace(sigma, GetTrace_sigma)
!      pn = 0.5*GetTrace(sigma)
!    
!    ! tension cutoff
!      if (pn<=Pmin) then 
!        pn = Pmin
!      end if 
!    
!      if (me2p == .false.) then 
!          ! shear modulus
!          GG = G0 * P_atm ! .false --> linear elastic.
!          
!          else 
!              GG = G0 * P_atm * sqrt(pn / P_atm) ! .true. --> nonlinear elastic (do not use this here)
!          end if 
!          
!          !if (nu == 0.5) then 
!          !    nu = 0.4999
!          !end if 
!          
!    
!    ! bulk modulus 
!      KK = two3 * ( (1.0+nu)/(1.0-(2.0*nu)) ) * GG !--> two3 needs to be a global number
!    
!      end subroutine GetElasticModuli                 
!                     
!                     
!                     
!                     
!    !-------------------------------------------------                
!      function GetStiffness(KK, GG) result(CC)
!    
!      implicit none
!    
!    ! input variables 
!      real(8), intent(in) :: KK, GG
!    
!    ! output variables 
!      real(8) :: CC(3, 3) ! 3x3 2D plane strain elasticity matrix !, intent(out)
!    
!    ! local variables
!      real(8) :: aa, bb
!      
!      CC = 0.0
!
!    ! compute coefficients of the elasticity matrix
!      aa = KK + ((4.0 / 3.0) * GG)
!      bb = KK - ((2.0 / 3.0) * GG)
!
!    ! fill stiffness matrix
!      CC(1, 1) = aa
!      CC(2, 2) = aa
!      CC(3, 3) = GG
!      CC(1, 2) = bb
!      CC(2, 1) = bb
!
!
!      end function GetStiffness              
!    
!    
!    
!    
!    !// returns the compliance matrix in its covariant-covariant form
!      subroutine GetCompliance(KK, GG, DD)
!      implicit none
!    ! Inputs
!      real(8), intent(in) :: KK, GG
!    ! Output
!      real(8), dimension(3, 3), intent(out) :: DD
!    ! Local variables
!      real(8) :: aa, bb, cc
!
!    ! Compute components of compliance matrix
!      aa = (KK + 4.0 / 3.0 * GG) / (4.0 * GG * KK + 4.0 / 3.0 * GG**2.0)
!      bb = (KK - 2.0 / 3.0 * GG) / (4.0 * GG * KK + 4.0 / 3.0 * GG**2.0)
!      cc = 1.0 / GG
!
!    ! Assign values to the compliance matrix
!      DD(1, 1) = aa
!      DD(2, 2) = aa
!      DD(3, 3) = cc
!      DD(1, 2) = bb
!      DD(2, 1) = bb
!
!      end subroutine GetCompliance
!    
!    !------------------------------------------------
!                     
!       subroutine GetElastoPlasticTangent(NextStress, aCe, RR, nn, K_p, &
!      Pmin, &
!      aCep)                 
!               
!       ! --> I dont know what this is doing.....
!       implicit none
!       ! Inputs
!       real(8), dimension(3), intent(in) :: NextStress
!       real(8), dimension(3, 3), intent(in) :: aCe
!       real(8), dimension(3), intent(in) :: RR
!       real(8), dimension(3), intent(in) :: nn
!       real(8), intent(in) :: K_p
!       ! Output
!       real(8), dimension(3, 3), intent(out) :: aCep
!       ! Local variables
!       real(8) :: pp, temp3, small
!       real(8), dimension(3) :: rrr, temp1, temp2
!      
!       ! Local variables 
!       real(8) :: GetTrace_result
!       real(8), dimension(3) :: GetDevPart_result
!       real(8), dimension(3) :: DoubleDot4_2_result
!       real(8) :: DoubleDot2_2_Contr_result
!       real(8), intent(in) :: Pmin
!      
!       real(8), dimension(3) :: DoubleDot2_4_result
!       
!       integer(4), dimension(3) :: mI1
!      
!       real(8), dimension(3,3) :: Dyadic2_2_result
!      
!       real(8), dimension(3,3) :: mIIco
!       
!       small = 1e-10
!      
!       mI1(1) = 1.0
!       mI1(2) = 1.0
!       mI1(3) = 0.0
!      
!       mIIco(1,1) = 1.0
!       mIIco(2,1) = 0.0
!       mIIco(3,1) = 0.0
!      
!       mIIco(1,2) = 0.0
!       mIIco(2,2) = 1.0
!       mIIco(3,2) = 0.0
!      
!       mIIco(1,3) = 0.0
!       mIIco(2,3) = 0.0
!       mIIco(3,3) = 1.0
!      
!       ! Calculate mean effective stress
!       pp = GetTrace(NextStress)
!       pp = 0.5 * pp
!       if (pp < Pmin) pp = Pmin
!       
!       ! Calculate deviatoric part of stress and normalize
!       rrr = GetDevPart(NextStress)
!       rrr = rrr / pp
!      
!       ! Compute temporary vectors and matrices
!       temp1 = DoubleDot4_2(aCe, RR) 
!      
!       temp2 = & 
!     DoubleDot2_4((nn - 1.0 / 2.0 * DoubleDot2_2_Contr(nn, rrr)) * mI1, &
!         aCe * mIIco)
!      
!       temp3 = DoubleDot2_2_Contr(temp2, RR) + K_p
!      
!       ! Update aCep based on plasticity condition
!       if (temp3 < small) then
!         aCep = aCe
!       else
!         !Dyadic2_2_result = Dyadic2_2(temp1, temp2)
!         aCep = aCe - ((1.0 / temp3) * Dyadic2_2(temp1, temp2))
!         end if
!         
!         aCep(3,1) = 0.0
!         aCep(3,2) = 0.0
!         aCep(1,3) = 0.0
!         aCep(2,3) = 0.0
!      
!      
!       end subroutine GetElastoPlasticTangent
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
!    !-----------------------------------------------
!      function GetNormalToYield(stress, alpha) result(nn)
!    !-----------------------------------------------
!    ! To calculate the normal to the yield surface
!    !-----------------------------------------------
!      implicit none 
!    
!    ! output variables
!      real(8), dimension(3) :: nn !, intent(out)
!    
!    ! input variables 
!      real(8), intent(in), dimension(3) :: stress, alpha
!    
!    ! local variables
!      real(8) :: normN 
!      real(8) :: pp
!      real(8) :: small
!      real(8) :: root12
!      real(8), dimension(3) :: GetDevPart_stress
!    
!      ! small 
!      small = 1e-6
!      
!    ! define numeric values
!      root12 = 0.7071067811865475
!    
!      pp = GetTrace(stress)
!      pp = 0.5 * pp
!    
!      !call GetDevPart(stress, GetDevPart_stress)
!    
!      if (abs(pp)<small) then
!        nn(3) = root12 ! --> this was 2 and needs to be 3
!      else 
!        
!        nn = (alpha*-pp) + GetDevPart(stress)
!        normN = GetNorm_Contr(nn)
!        
!        if (normN<small) then 
!            normN = 1
!        end if
!        nn = nn * (1/normN)
!        
!      end if
!    
!    
!      end function GetNormalToYield
!    
!    
!    
!    
!    
!    
!    
!    
!    !-------------------------------------------------------------------
!    
!  
!    
!    
!    
!    !*************************************************************
!    ! GetStateDependent() ----------------------------------------
!    ! Alsardi: we are calculating the state parameters for PM4Sand
!    ! mainly:
!    ! -Ksi
!    ! -Surface Ms
!    ! -Plastic modulus
!    ! -Dlatancy parameters
!    ! -R (plastic potential derivative) --> m in Alba's notes
!      subroutine GetStateDependent(stress, alpha, alpha_in, alpha_in_p, &
!      alpha_in_true, fabric, fabric_in, &
!        GG, zcum, zpeak, pzp, Mcur, CurDr, &
!       Mc, nd, nb, Pmin, Pmin2, P_atm, &
!       mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
!       QQ_Bolton, RR_Bolton, &
!       nn, alphaB, alphaD, bb, RR, &
!       Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, ksi)
!    ! get all the state parameters associated with a stress state
!      implicit none
!    ! Input parameters
!      real(8), dimension(3), intent(in) :: stress
!      real(8), dimension(3), intent(in) :: alpha
!      real(8), dimension(3), intent(in) :: alpha_in
!      real(8), dimension(3), intent(in) :: alpha_in_p
!      real(8), dimension(3), intent(in) :: alpha_in_true
!      real(8), dimension(3), intent(in) :: fabric
!      real(8), dimension(3), intent(in) :: fabric_in
!      real(8), intent(in) :: GG
!      real(8), intent(in) :: zcum
!      real(8), intent(in) :: zpeak
!      real(8), intent(in) :: pzp
!      real(8), intent(in) :: Mcur
!      real(8), intent(in) :: CurDr
!    
!      real(8), intent(in) :: Mc
!      real(8), intent(in) :: nd
!      real(8), intent(in) :: nb
!      real(8), intent(in) :: Pmin
!      real(8), intent(in) :: Pmin2
!      real(8), intent(in) :: P_atm
!      real(8), intent(in) :: mm
!      real(8), intent(in) :: z_max
!      real(8), intent(in) :: h0
!      real(8), intent(in) :: Ckaf
!      real(8), intent(in) :: Ado
!      real(8), intent(in) :: Ceps
!      real(8), intent(in) :: hpo
!      real(8), intent(in) :: Cdr
!    
!      real(8), intent(in) :: QQ_Bolton
!      real(8), intent(in) :: RR_Bolton
!    
!    ! Output parameters
!      real(8), dimension(3), intent(out) :: nn
!      real(8), dimension(3), intent(out) :: alphaB
!      real(8), dimension(3), intent(out) :: alphaD
!      real(8), dimension(3), intent(out) :: bb
!      real(8), dimension(3), intent(out) :: RR
!    
!      real(8), intent(out) :: ksi
!    
!      real(8), intent(out) :: Mb
!      real(8), intent(out) :: Md
!      real(8), intent(out) :: DD
!      real(8), intent(out) :: K_p
!      real(8), intent(out) :: Cka
!      real(8), intent(out) :: hh
!      real(8), intent(out) :: AlphaAlphaBDotN
!    
!    ! Local variables
!      real(8), dimension(3) :: alphaD_alpha 
!      real(8), dimension(3) :: alphaDr_alpha
!      real(8), dimension(3) :: alpha_Alpha_in
!      real(8), dimension(3) :: alpha_Alpha_in_true
!      real(8), dimension(3) :: alpha_Alpha_p
!      real(8), dimension(3) :: minusFabric
!    
!      real(8) :: Czpk1
!      real(8) :: Czpk2
!      real(8) :: Cpzp
!      real(8) :: Cpzp2
!      real(8) :: Cg1
!      real(8) :: Ckp
!    
!      real(8) :: AlphaAlphaInDotN
!      real(8) :: AlphaAlphaInTrueDotN
!      real(8) :: Czin1
!      real(8) :: Czin2
!      real(8) :: Crot1
!      real(8) :: Mdr
!      real(8) :: pp
!    
!    
!      real(8) :: Cpmin
!    
!      real(8) :: Crot2
!      real(8) :: C_pmin2 
!      real(8) :: Cin 
!      real(8) :: Cdz
!    
!      real(8) :: temp
!      real(8) :: Ad
!    
!      real(8) :: hp 
!      real(8) :: Adc
!    
!      real(8) :: Drot
!    
!    ! Local variables
!      real(8) :: root12
!      real(8) :: two3
!      real(8) :: one3
!      real(8), dimension(3) :: I1
!      
!      real(8) :: small 
!      
!      small = 1e-6
!   
!    ! identity matrix
!      I1(1) = 1
!      I1(2) = 1
!      I1(3) = 0
!    
!    ! two thirds 
!      two3 = 0.66666666667
!    ! one third
!      one3 = 0.33333333333
!    
!    ! sqrt(0.5)
!      root12 = 0.70710678118
!      
!      if (Ceps <0 ) then
!          print *, "Ceps is less than zero!"
!      end if
!          
!          
!          
!    
!    ! Mean effective stress
!      pp = GetTrace(stress)
!      pp = 0.5d0 * pp
!    
!      if (pp <= Pmin) then 
!          pp = Pmin ! kind of a tension cutoff
!      end if 
!      
!    
!    ! evaluate relative state parameter index based on Bolton (1986) relationship
!      ksi = GetKsi(CurDr, pp, RR_Bolton, QQ_Bolton, Pmin, P_atm) ! 
!    
!    ! get the normal to the yield surface
!      nn = GetNormalToYield(stress, alpha) 
!    
!    ! update bounding and dilatancy using the relative state parameter (ksi)
!      if (ksi <= 0.0d0) then ! negative relative state parameter index
!        ! dense of critical
!        Mb = Mc * exp(-1.0d0 * nb * ksi) ! bounding ratio based on ksi using m_nb
!        Md = Mc * exp(         nd * ksi) ! dilatancy ratio update based on ksi using m_nd
!      else ! positive relative state parameter index
!        ! loose of critical
!        Mb = Mc * exp(-1.0d0 * (nb / 4.0d0) * ksi) ! bounding ratio based on ksi using m_nb/4
!        Md = Mc * exp(         (nd * 4.0d0) * ksi) ! dilatancy ratio update based on ksi using m_nb*4
!      end if
!    
!    ! update image back-stress ratios for bounding surface
!      alphaB = nn * root12 * (Mb - mm)
!    
!    ! update image back-stress ratios for dilatancy surface
!      alphaD = nn * root12 * (Md - mm)
!    
!    ! calculate plastic modulus terms
!    ! C_zpk1 = z_peak / (z_cum + z_max/5)
!      Czpk1 = zpeak / (zcum + (z_max / 5.0d0)) ! Equation 62 --> Effect of fabric on plastic modulus
!    
!    ! C_zpk2 = z_peak / (z_cum + z_max/100)
!      Czpk2 = zpeak / (zcum + (z_max / 100.0d0)) ! Equation 63 --> Effect of fabric on plastic modulus
!    
!      if (Czpk2 > (1.0d0 - small)) Czpk2 = 1.0d0 - small ! there is a limit on this
!    
!    ! C_pzp2 = (-<-(p_zp - p)>)/(-<-(p_zp - p)> + p_min)
!      !call Macauley((pzp - pp), Macauley_result) --> Effect of fabric on plastic modulus
!      Cpzp2 = Macauley(pzp - pp) / (Macauley(pzp - pp) + Pmin) ! Equation 64 --> Effect of fabric on plastic modulus
!
!      Cg1 = h0 / 200.0d0 ! look at the table in "Plastic deviatoric strain increment" C_Gamma1 = h_o / 200
!      Ckp = 2.0d0 ! look at the table in "Plastic deviatoric strain increment" C_Kp = 2
!    
!    ! terms for k_p calculation
!      bb = alphaB - alpha ! difference between the bounding surface image back stress ratio and the back stress ratio
!    
!      AlphaAlphaBDotN = DoubleDot2_2_Contr(bb, nn)
!    
!    ! differnce between the back stress ratio and the initial back stress ratio
!      alpha_Alpha_in = alpha - Alpha_in ! alpha - alpha_in
!    
!    ! dot alpha_mAlpha_in with normal nn
!      !call DoubleDot2_2_Contr(alpha_Alpha_in, nn, 
!      !& DoubleDot2_2_Contr_result)
!      AlphaAlphaInDotN = &
!       Macauley(DoubleDot2_2_Contr(alpha_Alpha_in, nn))!, 
!      !& AlphaAlphaInDotN) ! Macauly bracket in the Kp calculation
!    
!      Alpha_Alpha_in_true = alpha - Alpha_in_true ! Alpha - Alpha^true_in
!    
!      !call DoubleDot2_2_Contr(alpha_Alpha_in_true, nn, 
!      !& DoubleDot2_2_Contr_result)
!    
!      AlphaAlphaInTrueDotN = &
!      Macauley(DoubleDot2_2_Contr(alpha_Alpha_in_true, nn)) !, 
!      !& AlphaAlphaInTrueDotN)
!    ! = Macauley_result!Macauley(DoubleDot2_2_Contr(alpha_mAlpha_in_true, n)) ! dot product in the denominator of Cka
!    
!      Cka = 1.0d0 + &
!      ( (Ckaf / (1.0d0 + ((2.5d0 * AlphaAlphaInTrueDotN)**2.0) )) &
!      * Cpzp2 * Czpk1 ) ! look at table "Plastic deviatoric strain increment"
!    ! updataed K_p formulation following PM4Sand V3.1. mAlpha_in is the apparent back-stress ratio.
!    
!      alpha_Alpha_p = alpha - alpha_in_p ! if alpha - alpha_in_p is small
!    
!
!        
!      if (abs(AlphaAlphaBDotN) < small) then
!        ! adding this condition to avoid division by zero error
!        hh = 1.0e10 ! part of Kp calculation
!      elseif (DoubleDot2_2_Contr(alpha_Alpha_p, nn) <= 0) then !--> alpha_mAlpha_p dot nn <=0
!        ! LOAD REVERSAL
!        hh = (1.5 * GG * h0) * ( sqrt(abs(AlphaAlphaBDotN)) &
!      / (pp * AlphaAlphaBDotN *(exp(AlphaAlphaInDotN) - 1.0 + Cg1) )) * &!* AlphaAlphaBDotN
!            (Cka / (1 + (Ckp * (zpeak / z_max)) * &
!       Macauley(AlphaAlphaBDotN) * &
!            sqrt(1 - Czpk2))) ! part of Kp calculation 
!        ! Equation 35 --> Hardening and the update of the back-stress ratio
!        ! This corresponds to equation 3.31 in Chen's PhD thesis to prevent shear stress overshooting
!      !   hh = hh * (AlphaAlphaInDotN + Cg1) / 
!      !&  (AlphaAlphaInTrueDotN + Cg1) ! part of Kp calculation --> Cg1 to prevent dividing by zero
!        hh = hh * (AlphaAlphaInDotN + Cg1)/(AlphaAlphaInTrueDotN + Cg1) ! part of Kp calculation --> Cg1 to prevent dividing by zero --> Crev
!      else 
!        ! NO LOAD REVERSAL
!
!      !   
!        hh = (1.5 * GG * h0) * ( sqrt(abs(AlphaAlphaBDotN)) &
!      / (pp * AlphaAlphaBDotN *(exp(min(AlphaAlphaInDotN, 300.0d0)) - 1.0 + Cg1) )) * &!* AlphaAlphaBDotN
!            (Cka / (1 + (Ckp * (zpeak / z_max)) * &
!       Macauley(AlphaAlphaBDotN) * &
!            sqrt(1 - Czpk2))) ! part of Kp calculation 
!        
!      end if
!        
!    ! after we evaluate h we plug it in K_p to find the plastic bulk modulus
!      !call DoubleDot2_2_Contr(bb, nn, DoubleDot2_2_Contr_result)
!      K_p = two3 * hh * pp * DoubleDot2_2_Contr(bb, nn) !DoubleDot2_2_Contr_result !DoubleDot2_2_Contr(b, n) ! Equation 34 --> Hardening and the update of the back-stress ratio
!
!    !--- calculate dilatancy
!      Czin1 = Macauley(1.0 - &
!      exp(-2.0*abs((DoubleDot2_2_Contr(fabric_in, nn) - &
!      DoubleDot2_2_Contr(fabric, nn)) / z_max)))
!      
!
!      ! Plastic volumetric strain increment... see in table
!      ! rotated dilatancy surface
!      minusFabric = fabric * (-1.0)
!
!      Crot1 = &
!     max(1.0 + ((2.0 * Macauley(DoubleDot2_2_Contr(minusFabric, nn)) / &
!      (sqrt(2.0)*z_max)) * (1 - Czin1)), 1.0) ! Plastic volumetric strain increment... see in table
!    
!    ! Equation 66 --> Effect of fabric on plastic volumetric dilation 
!      Mdr = Md / Crot1 ! Equation 65 --> Effect of fabric on plastic volumetric dilation 
!
!    
!      alphaDr_alpha = nn * (root12 * (Mdr - mm)) ! alphaDr
!      alphaDr_alpha = alphaDr_alpha - alpha ! Equation 68
!    
!      alphaD_alpha = alphaD - alpha
!
!    
!    
!      if ( DoubleDot2_2_Contr(alphaDr_alpha, nn) <= 0) then ! AlphaD - Alpha : nn --> negative --> DILATION 
!        
!        ! dilation
!        if (pzp == 0.0) then
!            Cpzp = 1.0 ! to prevent dividing by zero
!        else
!            Cpzp = 1.0 / (1.0 + ((2.5 * pp/pzp)**5.0)) ! plastic vol strain increment ! Equation 3.55
!        end if
!        
!        Cpmin = 1.0 / (1.0 + (Pmin2/pp)**2) ! plastic vol strain increment ! Equation 3.56
!        Czin2 = (1.0 + (Czin1 * (zcum - zpeak)) / (3.0*z_max)) &
!      / (1.0 + ((3.0 * Czin1 * (zcum - zpeak)) / (3.0*z_max))) ! plastic vol strain increment ! Equation 3.60
!        
!
!        temp = (1.0 - &
!       ( root12*Macauley(DoubleDot2_2_Contr(minusFabric, nn)) &
!      / ( zpeak )))** 3.0 !root12
!        
!        Ad = Ado * Czin2 / &
!      ((((zcum**2.0) / z_max) * temp * &
!        (Ceps**2.0) * Cpzp * Cpmin * Czin1) + 1.0) ! plastic vol strain increment ! Equation 72 !1.0 ! I changed this to 0.9
!
!        !1.0
!        
!        ! DD_non-rot in equation 70
!        DD = Ad * DoubleDot2_2_Contr(alphaD_alpha, nn) !DoubleDot2_2_Contr(alphaD_alpha, n) ! D_non-rot ! plastic vol strain increment
!        ! double Drot = Ad * Macauley(DoubleDot2_2_Contr(-1.0*fabric, n)) / (sqrt(2.0)*m_z_max) * DoubleDot2_2_Contr(alphaDr - alpha, n) / m_Cdr
!        
!          !2.0 * 5.0 *
!        
!        if (Cdr < 5.0) then
!            print *, "Cdr < 5.0!"    
!        end if
!        
!        Drot = Ad * ( Macauley(DoubleDot2_2_Contr(minusFabric, nn)) &
!       / (sqrt(2.0) * z_max) ) * &
!      ( DoubleDot2_2_Contr(alphaDr_alpha, nn) / Cdr )! D_rot ! plastic vol strain increment
!        
!        ! Equation 69 --> Drot
!        if (DD > Drot) then ! plastic vol strain increment --> D_non-rot < D_rot
!            !call Macauley(Mb - Mcur, Macauley_result)
!            DD = DD + &
!     ((Drot - DD) * Macauley(Mb - Mcur) / (Macauley(Mb - Mcur) + 0.01))! Equation 78 ! slight descrepancy in the denominator
!        end if
!        
!        ! during dilation at very low effective stresses (i.e., p<=2p_min), Di is constrained to ensure soil at 
!        ! dense than critical states continue to be dliative 
!        if ((Pmin <= pp) .and. (pp <= 2*Pmin)) then
!            !call Macauley(Mb - Md, Macauley_result)
!            ! equation 3.61 in Long Chen's PhD thesis
!            DD = min(DD, -3.5 * Ado * Macauley(Mb - Md) * & !-3.5
!       ((2.0*Pmin) - pp)/ Pmin) ! Equation 78 ! Equation 3.61
!        end if
!        
!        
!    
!        else
!        
!        ! contraction
!        
!        ! bound K_p to non - negative, following flac practice
!        K_p = max(0.0, K_p)
!        ! he is using power 2 instead of 2.5 in BZ14
!        !call Macauley(0.5 - ksi, Macauley_result)
!        hp = hpo * exp(-0.7 + (7.0*(Macauley(0.5 - ksi))**2.0)) ! Equation 53 --> Plastic volumetric strains - Contraction ! Equation 3.41
!        Crot2 = 1.0 - Czpk2 ! Equation 83
!        ! Equation 3.664 in Long PhD thesis
!        Cdz = max((1.0 - ((Crot2 * sqrt(2.0) * zpeak)/ z_max)) & 
!      * (z_max / (z_max + (Crot2 * zcum))), 1.0 / (1.0 + (z_max/2.0))) ! Equation 3.64
!        
!        
!        ! Equation 3.63 in Long PhD thesis
!        !Adc = m_Ado * (1.0 + Macauley(DoubleDot2_2_Contr(fabric, n))) / hp / Cdz
!        Adc = Ado * (1.0 + Macauley(DoubleDot2_2_Contr(fabric, nn))) &
!      / (hp * Cdz) ! Equation 3.63
!
!        Cin = 2.0 * Macauley(DoubleDot2_2_Contr(fabric, nn)) &
!       / (sqrt(2.0) * z_max) ! Equation 3.66
!
!
!        DD = &
!       min(Adc * (DoubleDot2_2_Contr(alpha_Alpha_in, nn) + Cin)**2.0, &
!      1.5 * Ado) * &
!      DoubleDot2_2_Contr(alphaD_alpha, nn) / &
!       (DoubleDot2_2_Contr(alphaD_alpha, nn) + 0.20) !0.10) !0.16) ! note that the threshold is enforced by using min function !0.16 !0.32 !0.235 !1.0
!        ! Apply a factor to D so it doesn't go very big when p is small
!        ! Equation 3.67 in Long PhD Thesis --> updated in PM4Sand 3.3
!        if (pp < Pmin * 2.0) then
!            C_pmin2 = 0.0
!        elseif (pp >= Pmin * 18.0) then
!            C_pmin2 = 1.0
!        else
!            C_pmin2 = (pp - (2.0 * Pmin)) / (16.0 * Pmin)
!        end if
!        
!        DD = DD * C_pmin2
!        
!        end if
!        
!        if (DD /= DD) then
!            print *, "DD is NaN!"    
!        end if
!    
!    ! R = n + one3 * D * mI1
!      RR = nn + ((one3 * DD)*I1) ! --> used in plastic deviatoric strain increment
!    
!    ! I do not see where we calculate the deviatoric component of RR to get RR'? RR' is equal to n because associative
!    
!      end subroutine GetStateDependent
!    
!    
!    
!    
!    !/*************************************************************/
!    !/*************************************************************/
!    !//            SYMMETRIC TENSOR OPERATIONS                    //
!    !/*************************************************************/
!    !/*************************************************************/
!    !// In all the functions below, by contravariant tensors, we mean stress-like tensors
!    !// and by covariant tensors we mean strain-like tensors
!      function GetTrace(vv)  result(GetTrace_result)
!    !// computes the trace of the input argument
!      implicit none
!    ! Arguments
!      real(8), dimension(3), intent(in) :: vv
!      real(8) :: GetTrace_result
!    
!    ! initialize
!      GetTrace_result = 0.0
!
!    ! Check vector size
!      if (size(vv) /= 3) then
!        print *, "ERROR! PM4Sand::GetTrace requires vector of size(3)!"
!        GetTrace_result = 0.0
!        return
!      endif
!    
!    ! Compute the trace
!      GetTrace_result = vv(1) + vv(2)
!
!      end function GetTrace
!
!
!    
!    
!      function GetDevPart(Vector) result(DevVector)
!    !-------------------------------------------
!    ! To calculate the deviatoric part of a matrix.
!    !-------------------------------------------
!      implicit none
!    
!    ! input variable 
!      real(8), dimension(3), intent(in) :: Vector
!      real(8), dimension(3) :: DevVector
!    
!    ! local variables 
!      real(8) :: pp
!    
!    ! initialize 
!      DevVector = 0.0
!      pp = 0.0
!    
!    ! Check vector size
!      if (size(Vector) /= 3) then
!        print *, &
!      "ERROR! PM4Sand::GetDevPart requires vector of size(3)!"
!        !result = 0.0
!        return
!      endif
!    
!      !call GetTrace(Vector, pp)
!    
!      DevVector = Vector 
!      DevVector(1) = DevVector(1) - 0.5*GetTrace(Vector)
!      DevVector(2) = DevVector(2) - 0.5*GetTrace(Vector)
!    
!    
!    
!      end function GetDevPart
!    
!    
!        
!        
!    
!      function DoubleDot2_2_Contr(v1, v2) result(res)
!    !--------------------------------------------------
!    ! // computes doubledot product for vector-vector arguments,
!    !    both "contravariant".
!      implicit none 
!    
!    ! the size of input vectors needs to be both 3
!      real(8), dimension(3), intent(in) :: v1
!      real(8), dimension(3), intent(in) :: v2
!      real(8) :: res !, intent(out)
!      integer(4) :: ii
!    
!      ! initialize 
!      res = 0.0
!
!    ! Check vector sizes
!      if (size(v1) /= 3 .or. size(v2) /= 3) then
!        print *, &
!      "ERROR! PM4Sand::DoubleDot2_2_Contr requires vectors of size(3)!"
!        !result = 0.0
!        return
!      endif
!    
!      ! initialize integer
!      res = 0
!      
!      do ii = 1, 3 ! 3 is for the size of the vector 
!        !result = result + v1(i) * v2(i) + (i > 1) * v1(i) * v2(i)
!        res = res + &
!      ( v1(ii) * v2(ii) )
!        if (ii == 3) then 
!            ! double count shear component according to cauchy stress tensor
!            res = res + &
!       ( v1(ii) * v2(ii) )
!        end if 
!      end do 
!    
!      end function !DoubleDot2_2_Contr
!    
!    
!    
!      function DoubleDot2_2_Cov(v1, v2) result(DoubleDot2_2_Cov_result)
!    !--------------------------------------------------------
!    ! // computes doubledot product for vector-vector arguments, both "covariant"
!    !
!      implicit none 
!    ! the size of input vectors needs to be both 3
!      real(8), dimension(3), intent(in) :: v1
!      real(8), dimension(3), intent(in) :: v2
!      real(8) :: DoubleDot2_2_Cov_result !, intent(out)
!      integer(4) :: ii
!    
!      
!    ! Check vector sizes
!      if (size(v1) /= 3 .or. size(v2) /= 3) then
!        print *, &
!      "ERROR! PM4Sand::DoubleDot2_2_Cov requires vectors of size(3)!"
!        !result = 0.0
!        return
!      endif
!    
!    ! initialize integer
!      DoubleDot2_2_Cov_result = 0
!    
!    
!    ! Compute the doubledot product
!      do ii = 1, 3
!        !result = result + v1(i) * v2(i) - (i > 2) * 0.5 * v1(i) * v2(i)
!        DoubleDot2_2_Cov_result = DoubleDot2_2_Cov_result + &
!      (v1(ii) * v2(ii))
!        if (ii == 3) then 
!            DoubleDot2_2_Cov_result = &
!      DoubleDot2_2_Cov_result - (0.5 * v1(ii) * v2(ii))
!        end if 
!        
!      end do 
!    
!      end function DoubleDot2_2_Cov
!    
!    
!      function DoubleDot2_2_Mixed(v1, v2) &
!      result(DoubleDot2_2_Mixed_result)
!    ! // computes doubledot product for vector-vector arguments, 
!    ! one "covariant" and the other "contravariant" 
!      implicit none 
!    ! the size of input vectors needs to be both 3
!      real(8), dimension(3), intent(in) :: v1
!      real(8), dimension(3), intent(in) :: v2
!      real(8) :: DoubleDot2_2_Mixed_result !intent(out)
!      integer(4) :: ii
!    
!    ! Check vector sizes
!      if (size(v1) /= 3 .or. size(v2) /= 3) then
!        print *, &
!      "ERROR! PM4Sand::DoubleDot2_2_Mixed requires vectors of size(3)!"
!        !result = 0.0
!        return
!      endif
!    
!    ! initialize integer
!      DoubleDot2_2_Mixed_result = 0
!    
!      do ii = 1, 3
!        DoubleDot2_2_Mixed_result = &
!      DoubleDot2_2_Mixed_result + ( v1(ii) * v2(ii) )
!      end do
!    
!      end function DoubleDot2_2_Mixed
!    
!    
!      function GetNorm_Contr(vv) result(GetNorm_Contr_result)
!    ! // computes contravariant (stress-like) norm of input 6x1 tensor
!      implicit none
!    
!      real(8), dimension(3), intent(in) :: vv
!      real(8) :: GetNorm_Contr_result
!    
!    
!    
!    ! Check vector size
!      if (size(vv) /= 3) then 
!        print *, &
!      "ERROR! PM4Sand::GetNorm_Contr requires vector of size(3)!"
!        !result = 0.0
!        return
!      endif
!    
!      ! initialize integer
!      GetNorm_Contr_result = 0
!        
!      GetNorm_Contr_result = sqrt(DoubleDot2_2_Contr(vv, vv))
!    
!      end function GetNorm_Contr
!    
!    
!    
!    
!    
!    
!    
!      function GetNorm_Cov(vv) result(GetNorm_Cov_result)
!    ! // computes covariant (strain-like) norm of input 6x1 tensor
!      implicit none 
!    
!      real(8), dimension(3), intent(in) :: vv
!      real(8) :: GetNorm_Cov_result !, intent(out)
!    
!      !real(8) :: DoubleDot2_2_Cov_result
!    
!    ! Check vector size
!      if (size(vv) /= 3) then
!        print *, &
!       "ERROR! PM4Sand::GetNorm_Cov requires vector of size(3)!"
!        !result = 0.0
!        return
!      endif
!    
!    ! initialize integer
!      GetNorm_Cov_result = 0
!     
!    
!    ! square root to find GetNorm_Cov
!      GetNorm_Cov_result = sqrt(DoubleDot2_2_Cov(vv, vv))
!    
!    
!      end function GetNorm_Cov
!    
!    
!    
!    
!      function Dyadic2_2(v1, v2) result(Dyadic2_2_result)
!    !// computes dyadic product for two vector-storage arguments
!    !// the coordinate form of the result depends on the coordinate form of inputs
!    ! https://en.wikipedia.org/wiki/Dyadics
!      implicit none 
!    
!    ! input vectors v1 and v2
!      real(8), dimension(3), intent(in) :: v1
!      real(8), dimension(3), intent(in) :: v2
!    ! output matrix result Dyadic2_2_result
!      real(8), dimension(3,3) :: Dyadic2_2_result !, intent(out)
!    
!    ! local variables 
!      integer(4) :: jj, ii
!    
!    ! Check vector sizes
!      if (size(v1) /= 3 .or. size(v2) /= 3) then
!        print *, "ERROR! PM4Sand::Dyadic2_2 requires vector of size(3)!"
!        return
!      endif
!    
!    ! initialize matrix
!      Dyadic2_2_result = 0.0
!    
!      do ii = 1, 3
!        do jj = 1, 3
!            
!            Dyadic2_2_result(ii,jj) = v1(ii) * v2(jj)
!            
!        end do 
!      end do 
!    
!    
!      end function Dyadic2_2
!    
!    
!    
!    
!      function DoubleDot4_2(m1, v1) result(DoubleDot4_2_result)
!    !// computes doubledot product for matrix-vector arguments
!    !// caution: second coordinate of the matrix should be in opposite variant form of vector
!      implicit none 
!    ! input variables 
!      real(8), intent(in), dimension(3,3) :: m1
!      real(8), intent(in), dimension(3) :: v1
!    ! output variables 
!      real(8), dimension(3) :: DoubleDot4_2_result!m1Timesv1 ! intent(out),
!    ! local variables 
!      integer(4) :: ii
!      integer(4) :: jj
!    
!    ! Check vector size
!      if (size(v1) /= 3) then
!        print *, &
!      "ERROR! PM4Sand::DoubleDot4_2 requires vector of size(3)!"
!        return
!      endif
!
!    ! Check matrix size
!      if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3) then
!        print *, "ERROR! PM4Sand::DoubleDot4_2 requires 3-by-3 matrix "
!        return
!      endif
!    
!    
!    ! Compute the result
!      do ii = 1, 3
!        DoubleDot4_2_result(ii) = 0.0
!        do jj = 1, 3
!            DoubleDot4_2_result(ii) = DoubleDot4_2_result(ii) + &
!       ( m1(ii, jj) * v1(jj) )
!        enddo
!      enddo
!    
!      end function DoubleDot4_2
!    
!    
!    
!    
!    
!    
!      function DoubleDot2_4(v1, m1) result(DoubleDot2_4_result)
!    !// computes doubledot product for matrix-vector arguments
!    !// caution: first coordinate of the matrix should be in opposite 
!    !// variant form of vector
!      implicit none
!    ! Arguments
!      real(8), dimension(3), intent(in) :: v1
!      real(8), dimension(3, 3), intent(in) :: m1
!      real(8), dimension(3) :: DoubleDot2_4_result !, intent(out)
!    ! Local variables
!      integer(4) :: ii, jj
!
!    ! Check vector size
!      if (size(v1) /= 3) then
!        print *, &
!      "ERROR! PM4Sand::DoubleDot2_4 requires vector of size(3)!"
!        return
!      endif
!
!    ! Check matrix size
!      if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3) then
!        print *, "ERROR! PM4Sand::DoubleDot2_4 requires 3-by-3 matrix "
!        return
!      endif
!
!    ! Compute the result
!      do ii = 1, 3
!        DoubleDot2_4_result(ii) = 0.0
!        do jj = 1, 3
!            DoubleDot2_4_result(ii) = DoubleDot2_4_result(ii) + &
!      ( m1(jj, ii) * v1(jj) )
!        enddo
!      enddo
!
!      end function DoubleDot2_4
!    
!    
!    
!    
!
!
!      function DoubleDot4_4(m1, m2) result(DoubleDot4_4_result)
!    !// computes doubledot product for matrix-matrix arguments
!    !// caution: second coordinate of the first matrix should be in opposite 
!    !// variant form of the first coordinate of second matrix
!      implicit none
!    ! Arguments
!      real(8), dimension(3, 3), intent(in) :: m1, m2
!      real(8), dimension(3, 3) :: DoubleDot4_4_result !, intent(out)
!    ! Local variables
!      integer(4) :: ii, jj, kk
!
!    ! Check matrix sizes
!      if (size(m1, 1) /= 3 .or. size(m1, 2) /= 3 .or. &
!         size(m2, 1) /= 3 .or. size(m2, 2) /= 3) then
!        print *, "ERROR! PM4Sand::DoubleDot4_4 requires 3-by-3 matrices"
!        return
!      endif
!
!    ! Compute the result
!      do ii = 1, 3
!        do jj = 1, 3
!            DoubleDot4_4_result(ii, jj) = 0.0
!            do kk = 1, 3
!                DoubleDot4_4_result(ii, jj) = &
!      DoubleDot4_4_result(ii, jj) + ( m1(ii, kk) * m2(kk, jj) )
!            enddo
!        enddo
!      enddo
!    
!    
!    
!      end function DoubleDot4_4
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
!      function ToContraviant(v1) result(ToContraviant_result)
!       !/*************************************************************/
!       !// ToContraviant() ---------------------------------------------
!      implicit none
!      ! Arguments
!      real(8), dimension(3), intent(in) :: v1
!      real(8), dimension(3) :: ToContraviant_result !, intent(out)
!      ! Local variables
!      integer :: ii
!    
!      ! initialize
!      ToContraviant_result = 0.0
!
!      ! Check vector size
!       if (size(v1) /= 3) then
!         print *, &
!      "ERROR! PM4Sand::ToContraviant requires vector of size(3)!"
!        return
!      endif
!
!      ! Copy input vector to result
!      ToContraviant_result = v1
!
!       ! Adjust second component of the result --> getting back tensorial strain
!      ToContraviant_result(3) = ToContraviant_result(3) * 0.5
!
!       end function ToContraviant
!    
!    
!    
!      
!      function ToCovariant(v1) result(ToCovariant_result)
!      !/*************************************************************/
!      !// ToCovariant() ---------------------------------------------
!      implicit none
!       ! Arguments
!      real(8), dimension(3), intent(in) :: v1
!      real(8), dimension(3) :: ToCovariant_result
!       ! Local variables
!      integer(4) :: ii
!
!    ! Check vector size
!      if (size(v1) /= 3) then
!        print *, &
!      "ERROR! PM4Sand::ToCovariant requires vector of size(3)!"
!        return
!      endif
!
!    ! Copy input vector to result
!      ToCovariant_result = v1
!
!    ! Adjust second component of the result
!      ToCovariant_result(3) = ToCovariant_result(3) * 2.0
!
!
!      end function ToCovariant
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
!    end module PM4SandMaterialModule
!
!    
    
    
    !======================================================================

    
    
    !======================================================================
! PM4Sand Material Module
! Fortran port of PM4Sand (Version 3.1) by Boulanger & Ziotopoulou
! Reference C++ implementation: OpenSees PM4Sand by L.Chen, P.Arduino
!
! Corrections vs original Fortran:
!   1. Contraction threshold 0.20 -> 0.16 (matches C++)
!   2. C_pmin2 applied in dilation branch (was missing)
!   3. PostShake modulus reduction restored (was commented out)
!   4. Void ratio committed every step (was commented out)
!   5. nu not modified in GetElasticModuli_ (intent(in) bug fixed)
!   6. DSTRAN/PROPS declared intent(in) (were wrongly inout)
!   7. Mb/Md are out from GetStateDependent (recomputed each call)
!   8. No double-write of NextVoidRatio
!======================================================================

      module PM4SandMaterialModule
      !implicit none

      contains

!----------------------------------------------------------------------
! ESM WRAPPER  (called by the finite element driver)
!----------------------------------------------------------------------
      SUBROUTINE ESM_AM_PM4SAND(NPT, NOEL, IDSET, STRESS, &
        EUNLOADING, PLASTICMULTIPLIER, &
        DSTRAN, NSTATEV, STATEV, NADDVAR, ADDITIONALVAR, CMNAME, &
        NPROPS, PROPS, NUMBEROFPHASES, NTENS)

      !implicit none

      !character(len=80), intent(in)    :: CMNAME
      !integer(4),        intent(in)    :: NPT, NOEL, IDSET
      !integer(4),        intent(in)    :: NTENS, NUMBEROFPHASES
      !integer(4),        intent(in)    :: NSTATEV, NADDVAR, NPROPS
      !
      !real(8), intent(inout)           :: PLASTICMULTIPLIER
      !real(8), intent(inout), dimension(NSTATEV)     :: STATEV
      !real(8), intent(inout), dimension(NADDVAR)     :: ADDITIONALVAR
      !real(8), intent(inout), dimension(NPROPS)      :: PROPS
      !real(8), intent(inout), dimension(NTENS)       :: STRESS
      !real(8), intent(inout), dimension(NTENS)       :: DSTRAN
      !real(8), intent(inout) :: EUNLOADING !, dimension(NTENS,NTENS)
      !
      !! local UMAT-style variables
      !integer(4) :: IStep, TimeStep, IDTask
      !integer(4) :: ndi, nshr, layer, kspt, kstep, kinc
      !
      !real(8), dimension(NTENS)       :: ddsddt, drplde, stran, dpred
      !real(8), dimension(1)           :: predef
      !real(8), dimension(2)           :: time
      !real(8), dimension(3)           :: coords
      !real(8), dimension(3,3)         :: drot, dfgrd0, dfgrd1
      !real(8), dimension(NTENS,NTENS) :: ddsdde
      !
      !real(8) :: sse, spd, scd, rpl, drpldt, pnewdt
      !real(8) :: dtime, temp, dtemp, celent
      !real(8) :: Porosity, WaterPressure, WaterPressure0
      !real(8) :: GasPressure, GasPressure0, DegreeSaturation
      
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

      ! initialise outputs
      EUNLOADING      = 0.0d0
      PLASTICMULTIPLIER = 0.0d0
      IDTask          = 2

      ! unpack additional variables
      Porosity         = ADDITIONALVAR(1)
      WaterPressure    = ADDITIONALVAR(2)
      WaterPressure0   = ADDITIONALVAR(3)
      GasPressure      = ADDITIONALVAR(4)
      GasPressure0     = ADDITIONALVAR(5)
      DegreeSaturation = ADDITIONALVAR(6)
      time(1)          = ADDITIONALVAR(7)   ! TotalRealTime
      time(2)          = ADDITIONALVAR(8)   ! OverallTotalTime
      dtime            = ADDITIONALVAR(9)   ! TimeIncrement
      IStep            = int(ADDITIONALVAR(10))
      TimeStep         = int(ADDITIONALVAR(11))

      IF ((IStep == 1) .and. (TimeStep == 1)) IDTask = 1

      call UMAT_PM4Sand(STRESS, STATEV, ddsdde, sse, &
        spd, scd, rpl, &
        ddsddt, drplde, drpldt, stran, DSTRAN, time, dtime, temp, &
        dtemp, predef, dpred, CMNAME, ndi, nshr, NTENS, NSTATEV, &
        PROPS, NPROPS, coords, drot, pnewdt, celent, dfgrd0, &
        dfgrd1, NOEL, NPT, layer, kspt, kstep, kinc)

      ! set max diagonal for time-step control
      EUNLOADING = max(ddsdde(1,1), ddsdde(2,2), ddsdde(3,3))

      end subroutine ESM_AM_PM4SAND

!======================================================================
! UMAT
! Standard Abaqus UMAT interface.
! Reads STATEV/PROPS, calls constructor+integrate+commit, writes back.
!======================================================================
      SUBROUTINE UMAT_PM4Sand(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
        RPL, DDSDDT, DRPLDE, DRPLDT, &
        STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, &
        CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
        COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, &
        NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

      implicit none

      !--- dummy arguments (Abaqus convention) -----------------------
      character(len=80), intent(in)    :: CMNAME
      integer(4),        intent(in)    :: NPT, NTENS, NOEL
      integer(4),        intent(in)    :: NSTATEV, NPROPS
      integer(4),        intent(in)    :: NDI, NSHR, LAYER, KSPT
      integer(4),        intent(in)    :: KSTEP, KINC

      real(8), intent(in)              :: SSE, SPD, SCD, RPL
      real(8), intent(in)              :: DRPLDT, PNEWDT
      real(8), intent(in)              :: DTIME, TEMP, DTEMP, CELENT

      ! intent(in) — we only read these
      real(8), intent(in), dimension(NTENS)      :: DDSDDT, DRPLDE
      real(8), intent(in), dimension(NTENS)      :: STRAN
      real(8), intent(in), dimension(NTENS)      :: DSTRAN
      real(8), intent(in), dimension(2)          :: TIME
      real(8), intent(in), dimension(1)          :: PREDEF, DPRED
      real(8), intent(in), dimension(3)          :: COORDS
      real(8), intent(in), dimension(3,3)        :: DROT
      real(8), intent(in), dimension(3,3)        :: DFGRD0, DFGRD1
      real(8), intent(in), dimension(NPROPS)     :: PROPS

      ! intent(inout) — read at entry, written at exit
      real(8), intent(inout), dimension(NTENS)       :: STRESS
      real(8), intent(inout), dimension(NSTATEV)     :: STATEV
      real(8), intent(inout), dimension(NTENS,NTENS) :: DDSDDE

      !--- primary model parameters (PROPS) --------------------------
      real(8) :: Dr, G0, hpo
      logical :: PostShake, me2p

      !--- secondary parameters (STATEV 1-31) ------------------------
      real(8) :: P_atm, h0, emax, emin, nb, nd, Ado
      real(8) :: z_max, cz, ceps, phi_cv, nu
      real(8) :: Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm
      real(8) :: InitialVoidRatio
      real(8) :: Mc, Mb, Md, Cdr
      real(8) :: Fsed_min, p_sedo, Pmin, Pmin2
      real(8) :: K_p, zpeak, zcum, pzp, zxp

      !--- evolving scalar state (STATEV 32-34) ----------------------
      real(8) :: KK, GG, Mcur
      !real(8) :: TolF, TolR

      !--- evolving tensor state (STATEV 35-113) ---------------------
      real(8), dimension(3) :: rrr, Sigma_b
      real(8), dimension(3) :: Fabric_n, Fabric, Fabric_in_n, Fabric_in
      real(8), dimension(3) :: Alpha_n, Alpha
      real(8), dimension(3) :: Alpha_in_n,     Alpha_in
      real(8), dimension(3) :: Alpha_in_p_n,   Alpha_in_p
      real(8), dimension(3) :: Alpha_in_true_n, Alpha_in_true
      real(8), dimension(3) :: Alpha_in_max_n,  Alpha_in_max
      real(8), dimension(3) :: Alpha_in_min_n,  Alpha_in_min
      real(8), dimension(3) :: dEpsilonE
      real(8), dimension(3) :: Sigma_n, Sigma
      real(8), dimension(3) :: Epsilon, Epsilon_n
      real(8), dimension(3) :: EpsilonE, EpsilonE_n

      !--- scalar state (STATEV 114) ---------------------------------
      real(8) :: NextVoidRatio
      
      !--- tolerances (STATEV 90-91) ---------------------------------
      real(8) :: TolF, TolR

      !--- flags (STATEV 89, 101) ------------------------------------
      logical :: pzpFlag, FirstCall

      !--- working variables -----------------------------------------
      real(8), dimension(3)   :: StrainIncrement, nn, RR
      real(8), dimension(3,3) :: Ce, Cep
      integer(4) :: ii

      !================================================================
      ! 1.  READ PRIMARY PROPERTIES
      !================================================================
      Dr        = PROPS(1)
      G0        = PROPS(2)
      hpo       = PROPS(3)
      PostShake = (PROPS(4) /= 0.0d0)
      me2p      = (PROPS(5) /= 0.0d0)

      !================================================================
      ! 2.  READ STATE VARIABLES
      !================================================================
      P_atm            = STATEV(1)
      h0               = STATEV(2)
      emax             = STATEV(3)
      emin             = STATEV(4)
      nb               = STATEV(5)
      nd               = STATEV(6)
      Ado              = STATEV(7)
      z_max            = STATEV(8)
      cz               = STATEV(9)
      ceps             = STATEV(10)
      phi_cv           = STATEV(11)
      nu               = STATEV(12)
      Cgd              = STATEV(13)
      Ckaf             = STATEV(14)
      QQ_Bolton        = STATEV(15)
      RR_Bolton        = STATEV(16)
      mm               = STATEV(17)
      InitialVoidRatio = STATEV(18)
      Mc               = STATEV(19)
      Mb               = STATEV(20)
      Md               = STATEV(21)
      Cdr              = STATEV(22)
      Fsed_min         = STATEV(23)
      p_sedo           = STATEV(24)
      Pmin             = STATEV(25)
      Pmin2            = STATEV(26)
      K_p              = STATEV(27)
      zpeak            = STATEV(28)
      zcum             = STATEV(29)
      pzp              = STATEV(30)
      zxp              = STATEV(31)
      KK               = STATEV(32)
      GG               = STATEV(33)
      Mcur             = STATEV(34)
      rrr              = STATEV(35:37)
      Sigma_b          = STATEV(38:40)
      Fabric_n         = STATEV(41:43)
      Fabric           = STATEV(44:46)
      Fabric_in_n      = STATEV(47:49)
      Fabric_in        = STATEV(50:52)
      Alpha_n          = STATEV(53:55)
      Alpha            = STATEV(56:58)
      Alpha_in_n       = STATEV(59:61)
      Alpha_in         = STATEV(62:64)
      Alpha_in_p_n     = STATEV(65:67)
      Alpha_in_p       = STATEV(68:70)
      Alpha_in_true_n  = STATEV(71:73)
      Alpha_in_true    = STATEV(74:76)
      Alpha_in_max_n   = STATEV(77:79)
      Alpha_in_max     = STATEV(80:82)
      Alpha_in_min_n   = STATEV(83:85)
      Alpha_in_min     = STATEV(86:88)
      pzpFlag          = (STATEV(89) /= 0.0d0)
      TolF             = STATEV(90)
      TolR             = STATEV(91)
      dEpsilonE        = STATEV(92:94)
      Sigma_n          = STATEV(95:97)
      Sigma            = STATEV(98:100)
      FirstCall        = (STATEV(101) /= 0.0d0)
      Epsilon          = STATEV(102:104)
      Epsilon_n        = STATEV(105:107)
      EpsilonE         = STATEV(108:110)
      EpsilonE_n       = STATEV(111:113)
      NextVoidRatio    = STATEV(114)

      !================================================================
      ! 3.  STRAIN INCREMENT  (2D plane-strain, geotechnical convention)
      !     C++: mEpsilon = strain_from_element * (-1.0)
      !     So internal strains are compression-positive.
      !     DSTRAN is tension+ve from Abaqus, so we flip the sign.
      !     DSTRAN(4) is the engineering shear strain increment.
      !================================================================
      StrainIncrement(1) = -DSTRAN(1)
      StrainIncrement(2) = -DSTRAN(2)
      StrainIncrement(3) = -DSTRAN(4)

      ! NaN guard
      do ii = 1, NTENS
        if (DSTRAN(ii) /= DSTRAN(ii)) then
          print *, "UMAT PM4Sand: NaN in DSTRAN(", ii, &
        ") element", NOEL, " point", NPT
        end if
      end do
      
      if (FirstCall) then 
          Sigma_n(1) = stress(1) * (-1) 
          Sigma_n(2) = stress(2) * (-1) 
          Sigma_n(3) = stress(4) * (-1) 
          
      end if 

      !================================================================
      ! 4.  INITIALISE Ce FROM CURRENT MODULI  (used by constructor)
      !================================================================
      !Ce = GetStiffness(KK, GG)

      !================================================================
      ! 5.  CONSTRUCTOR  — initialises state on the very first call
      !     (FirstCall == .TRUE.); a no-op afterwards
      !================================================================
      call PM4SandFullConstructor( &
        Dr, G0, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, &
        phi_cv, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, &
        InitialVoidRatio, Mc, Mb, Md, Cdr, Fsed_min, p_sedo, &
        Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp, &
        KK, GG, Mcur, rrr, Sigma_b, &
        Fabric_n, Fabric, Fabric_in_n, Fabric_in, &
        Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, &
        Alpha_in_p, Alpha_in_true_n, Alpha_in_true, &
        Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, &
        pzpFlag, TolF, TolR, dEpsilonE, Sigma_n, Sigma, &
        FirstCall, Ce, &
        Epsilon, Epsilon_n, EpsilonE, EpsilonE_n, &
        NextVoidRatio)
      
      Ce = GetStiffness(KK, GG)
      
      

      !================================================================
      ! 6.  ADVANCE TOTAL STRAIN
      !================================================================
      Epsilon = Epsilon + StrainIncrement

      !================================================================
      ! 7.  PLASTIC INTEGRATION
      !================================================================
      call PM4SandIntegrate( &
        Epsilon, Epsilon_n, EpsilonE, EpsilonE_n, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, pzp, zxp, &
        KK, GG, Mcur, rrr, Sigma_b, &
        Fabric_n, Fabric, Fabric_in_n, Fabric_in, &
        Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, &
        Alpha_in_p, Alpha_in_true_n, Alpha_in_true, &
        Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, &
        pzpFlag, TolF, TolR, dEpsilonE, Sigma_n, Sigma, &
        Ce, nn, RR, NextVoidRatio)

      !================================================================
      ! 8.  COMMIT STATE
      !     Updates _n variables, fabric accumulators, void ratio,
      !     elastic moduli, and returns the updated elastic tangent Ce.
      !================================================================
      call PM4SandCommitState( &
        z_max, nu, G0, mm, Mc, Md, Mb, PostShake, Pmin, P_atm, &
        Cgd, p_sedo, Fsed_min, me2p, &
        Fabric, Fabric_in, &
        Alpha_in, Alpha_in_p, Alpha_in_true, Alpha_in_max, Alpha_in_min, &
        KK, GG, Mcur, rrr, &
        Alpha_in_n, Alpha_n, Alpha_in_p_n, &
        Alpha_in_true_n, Alpha_in_max_n, Alpha_in_min_n, Sigma_n, &
        zcum, zpeak, Sigma, Alpha, Fabric_n, Fabric_in_n, &
        Ce, nn, RR, K_p, &
        Epsilon, Epsilon_n, EpsilonE, EpsilonE_n, &
        InitialVoidRatio, NextVoidRatio)

      !================================================================
      ! 9.  WRITE BACK ALL STATE VARIABLES
      !================================================================
      STATEV(1)       = P_atm
      STATEV(2)       = h0
      STATEV(3)       = emax
      STATEV(4)       = emin
      STATEV(5)       = nb
      STATEV(6)       = nd
      STATEV(7)       = Ado
      STATEV(8)       = z_max
      STATEV(9)       = cz
      STATEV(10)      = ceps
      STATEV(11)      = phi_cv
      STATEV(12)      = nu
      STATEV(13)      = Cgd
      STATEV(14)      = Ckaf
      STATEV(15)      = QQ_Bolton
      STATEV(16)      = RR_Bolton
      STATEV(17)      = mm
      STATEV(18)      = InitialVoidRatio
      STATEV(19)      = Mc
      STATEV(20)      = Mb
      STATEV(21)      = Md
      STATEV(22)      = Cdr
      STATEV(23)      = Fsed_min
      STATEV(24)      = p_sedo
      STATEV(25)      = Pmin
      STATEV(26)      = Pmin2
      STATEV(27)      = K_p
      STATEV(28)      = zpeak
      STATEV(29)      = zcum
      STATEV(30)      = pzp
      STATEV(31)      = zxp
      STATEV(32)      = KK
      STATEV(33)      = GG
      STATEV(34)      = Mcur
      STATEV(35:37)   = rrr
      STATEV(38:40)   = Sigma_b
      STATEV(41:43)   = Fabric_n
      STATEV(44:46)   = Fabric
      STATEV(47:49)   = Fabric_in_n
      STATEV(50:52)   = Fabric_in
      STATEV(53:55)   = Alpha_n
      STATEV(56:58)   = Alpha
      STATEV(59:61)   = Alpha_in_n
      STATEV(62:64)   = Alpha_in
      STATEV(65:67)   = Alpha_in_p_n
      STATEV(68:70)   = Alpha_in_p
      STATEV(71:73)   = Alpha_in_true_n
      STATEV(74:76)   = Alpha_in_true
      STATEV(77:79)   = Alpha_in_max_n
      STATEV(80:82)   = Alpha_in_max
      STATEV(83:85)   = Alpha_in_min_n
      STATEV(86:88)   = Alpha_in_min
      STATEV(89)      = merge(1.0d0, 0.0d0, pzpFlag)
      STATEV(90)      = TolF
      STATEV(91)      = TolR
      STATEV(92:94)   = dEpsilonE
      STATEV(95:97)   = Sigma_n
      STATEV(98:100)  = Sigma
      STATEV(101)     = merge(1.0d0, 0.0d0, FirstCall)
      STATEV(102:104) = Epsilon
      STATEV(105:107) = Epsilon_n
      STATEV(108:110) = EpsilonE
      STATEV(111:113) = EpsilonE_n
      STATEV(114)     = NextVoidRatio

      !================================================================
      ! 10. ASSEMBLE OUTPUT STRESS TENSOR  (6-component, tension +ve)
      !     C++: mSigma_r = mSigma + mSigma_b; mSigma_r *= -1.0
      !     Sigma_b holds the correction stored at initialisation when
      !     p0 < Pmin or stress was scaled to the bounding surface.
      !     Must be added before the sign flip to tension+ve convention.
      !================================================================
      STRESS(1) = -1.0d0 * (Sigma(1) + Sigma_b(1))
      STRESS(2) = -1.0d0 * (Sigma(2) + Sigma_b(2))
      STRESS(3) =  0.0d0
      STRESS(4) = -1.0d0 * (Sigma(3) + Sigma_b(3))
      STRESS(5) =  0.0d0
      STRESS(6) =  0.0d0

      !================================================================
      ! 11. ASSEMBLE TANGENT STIFFNESS  (6x6 from 3x3 Ce)
      !     C++ default TangType=0 returns elastic Ce.
      !     Plane-strain out-of-plane column/row copy from yy entry.
      !================================================================
      
      !================================================================
      ! 11. COMPUTE CEP AND ASSEMBLE TANGENT STIFFNESS
      !     Returns elasto-plastic tangent so stress-control driver
      !     converges without overshooting near the bounding surface.
      !================================================================
      Cep = Ce!GetElastoPlasticTangent(Sigma, Ce, RR, nn, K_p, Pmin, mm)
      
      DDSDDE = 0.0d0

      DDSDDE(1,1) = Cep(1,1)
      DDSDDE(2,1) = Cep(2,1)
      DDSDDE(3,1) = Cep(2,1)   ! zz couples same as yy

      DDSDDE(1,2) = Cep(1,2)
      DDSDDE(2,2) = Cep(2,2)
      DDSDDE(3,2) = Cep(1,2)

      DDSDDE(1,3) = Cep(2,1)
      DDSDDE(2,3) = Cep(1,2)
      DDSDDE(3,3) = Cep(2,2)

      DDSDDE(4,4) = Cep(3,3)
      DDSDDE(5,5) = Cep(3,3)
      DDSDDE(6,6) = Cep(3,3)

      return
      END SUBROUTINE UMAT_PM4Sand

     
      
      
      
!======================================================================
! PM4SandFullConstructor
! Initialises all state variables on the very first call.
! Matches C++ PM4Sand::initialize(Vector initStress).
! After the first call, FirstCall is set .FALSE. and the routine
! is a no-op.
!======================================================================
      subroutine PM4SandFullConstructor( &
        Dr, G0, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, &
        phi_cv, nu, Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, &
        InitialVoidRatio, Mc, Mb, Md, Cdr, Fsed_min, p_sedo, &
        Pmin, Pmin2, K_p, zpeak, zcum, pzp, zxp, &
        KK, GG, Mcur, rrr, Sigma_b, &
        Fabric_n, Fabric, Fabric_in_n, Fabric_in, &
        Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, &
        Alpha_in_p, Alpha_in_true_n, Alpha_in_true, &
        Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, &
        pzpFlag, TolF, TolR, dEpsilonE, Sigma_n, Sigma, &
        FirstCall, Ce, &
        Epsilon, Epsilon_n, EpsilonE, EpsilonE_n, &
        NextVoidRatio)

      implicit none

      ! Primary inputs (never change)
      real(8), intent(in)    :: Dr, G0
      logical, intent(in)    :: PostShake, me2p

      ! Secondary parameters — inout because defaults are written
      real(8), intent(inout) :: P_atm, h0, emax, emin, nb, nd, Ado
      real(8), intent(inout) :: z_max, cz, ceps, phi_cv, nu
      real(8), intent(inout) :: Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm
      real(8), intent(inout) :: InitialVoidRatio
      real(8), intent(inout) :: Mc, Mb, Md, Cdr
      real(8), intent(inout) :: Fsed_min, p_sedo, Pmin, Pmin2
      real(8), intent(inout) :: K_p, zpeak, zcum, pzp, zxp
      real(8), intent(inout) :: KK, GG, Mcur
      real(8), intent(inout) :: TolF, TolR
      real(8), intent(inout) :: NextVoidRatio

      ! Tensor state — inout
      real(8), intent(inout), dimension(3)   :: rrr, Sigma_b
      real(8), intent(inout), dimension(3)   :: Fabric_n, Fabric
      real(8), intent(inout), dimension(3)   :: Fabric_in_n, Fabric_in
      real(8), intent(inout), dimension(3)   :: Alpha_n, Alpha
      real(8), intent(inout), dimension(3)   :: Alpha_in_n,    Alpha_in
      real(8), intent(inout), dimension(3)   :: &
        Alpha_in_p_n,  Alpha_in_p
      real(8), intent(inout), dimension(3)   :: &
        Alpha_in_true_n, Alpha_in_true
      real(8), intent(inout), dimension(3)   :: &
        Alpha_in_max_n, Alpha_in_max
      real(8), intent(inout), dimension(3)   :: &
        Alpha_in_min_n, Alpha_in_min
      real(8), intent(inout), dimension(3)   :: dEpsilonE
      real(8), intent(inout), dimension(3)   :: Sigma_n, Sigma
      real(8), intent(inout), dimension(3)   :: Epsilon, Epsilon_n
      real(8), intent(inout), dimension(3)   :: EpsilonE, EpsilonE_n
      real(8), intent(inout), dimension(3,3) :: Ce

      ! Flags — inout
      logical, intent(inout) :: pzpFlag, FirstCall

      ! Local
      real(8), dimension(3) :: initStress, I1
      real(8) :: p0, ksi, Mcut, Mfin
      real(8) :: pi

      pi = 3.14159265358979d0

      I1(1) = 1.0d0;  I1(2) = 1.0d0;  I1(3) = 0.0d0

      ! Only initialise on first call
      if (.not. FirstCall) return

      !--- atmospheric pressure
      if (P_atm < 0.0d0) P_atm = 101.325d0

      !--- h0
      if (h0 < 0.0d0) h0 = max(0.3d0, (0.25d0 + Dr) / 2.0d0)

      !--- void ratio limits
      if (emax < 0.0d0) emax = 0.8d0 !0.7389!0.8d0 !0.7389!
      if (emin < 0.0d0) emin = 0.5d0 !0.49!0.5d0 !0.4915

      !--- dilatancy parameters
      if (nb < 0.0d0) nb = 0.5d0 !0.7d0 !0.5d0
      if (nd < 0.0d0) nd = 0.1d0

      !--- fabric parameters
      if (cz < 0.0d0) cz = 250.0d0 !200.0d0!

      !--- ce (ceps) — matches C++ exactly
      if (ceps < 0.0d0) then
        if (Dr > 0.75d0) then
          ceps = 0.2d0
        else if (Dr < 0.55d0) then
          ceps = 0.5d0
        else
          ceps = 0.5d0 - (Dr - 0.55d0) * 1.5d0
        end if
      end if

      !--- critical state friction angle / Mc
      if (phi_cv < 0.0d0) phi_cv = 33.0d0!30d0 !35.6d0!
      if (Mc < 0.0d0) Mc = 2.0d0 * sin(phi_cv * pi / 180.0d0)

      !--- Poisson ratio
      if (nu < 0.0d0) nu = 0.3d0

      !--- Cgd
      if (Cgd < 0.0d0) Cgd = 2.0d0

      !--- Cdr — matches C++: (5+25*(Dr-0.35)), capped at 10
      if (Cdr < 0.0d0) then
        Cdr = 5.0d0 + 25.0d0 * (Dr - 0.35d0)
        Cdr = min(Cdr, 10.0d0)
      end if

      !--- Ckaf — matches C++
      if (Ckaf < 0.0d0) then
        Ckaf = 5.0d0 + 220.0d0 * (Dr - 0.26d0)**3
        Ckaf = max(4.0d0, min(35.0d0, Ckaf))
      end if

      !--- Bolton parameters
      if (QQ_Bolton < 0.0d0) QQ_Bolton = 10.0d0 !9.5d0!10.0d0
      if (RR_Bolton < 0.0d0) RR_Bolton = 1.5d0 !0.7d0!1.5d0

      !--- yield surface size
      if (mm < 0.0d0) mm = 0.01d0

      !--- post-shaking parameters
      if (Fsed_min < 0.0d0) then
        Fsed_min = min(0.03d0 * exp(2.6d0 * Dr), 0.99d0)
      end if
      if (p_sedo < 0.0d0) p_sedo = P_atm / 5.0d0

      !--- store initial stress and find p0
      initStress = Sigma_n
      p0 = 0.5d0 * (initStress(1) + initStress(2))

      !--- minimum p'
      if (Pmin < 0.0d0) Pmin = max(p0 / 200.0d0, P_atm / 200.0d0)
      if (Pmin2 < 0.0d0) Pmin2 = Pmin * 10.0d0

      !--- tension cutoff: if p0 < Pmin adjust stress
      if (p0 < Pmin) then
        Sigma_n = Pmin * I1
        Sigma_b = initStress - Sigma_n
        p0      = Pmin
        Alpha   = 0.0d0
        Alpha_n = 0.0d0
      else
        Sigma_n = initStress
        Sigma_b = 0.0d0
        Alpha_n = GetDevPart(initStress) / p0
      end if

      !--- relative state parameter
      ksi = GetKsi(Dr, p0, RR_Bolton, QQ_Bolton, Pmin, P_atm)

      !--- z_max default
      if (z_max < 0.0d0) z_max = min(0.7d0 * exp(-6.1d0 * ksi), 20.0d0)

      !--- bounding / dilatancy surfaces — matches C++ initialize()
      if (ksi < 0.0d0) then
        ! dense of critical
        Mb = Mc * exp(-1.0d0 * nb * ksi)
        Md = Mc * exp(        nd * ksi)
        if (Ado < 0.0d0) then
          if (Mb > 2.0d0) then
            Ado = 1.5d0
          else
            Ado = 2.5d0 * (asin(Mb/2.0d0) - asin(Mc/2.0d0)) &
        / (Mb - Md)
          end if
        end if
      else
        ! loose of critical
        Mb = Mc * exp(-1.0d0 * (nb/4.0d0) * ksi)
        Md = Mc * exp(        (nd*4.0d0)  * ksi)
        if (Ado < 0.0d0) Ado = 1.24d0
      end if

      !--- check initial stress vs bounding/dilatancy surface
      Mcut = max(Mb, Md)
      Mfin = sqrt(2.0d0) * GetNorm_Contr(GetDevPart(Sigma_n)) / p0

      if (Mfin > Mcut) then
        rrr     = (Sigma_n - p0*I1) / p0 * (Mcut/Mfin)
        Sigma_n = p0*I1 + rrr*p0
        Sigma_b = initStress - Sigma_n
        Alpha_n = rrr * (Mcut - mm) / Mcut
      end if

      !--- zero cumulated fabric
      zcum = 0.0d0

      !--- elastic moduli
      call GetElasticModuli_(Sigma_n, zcum, z_max, nu, G0, Md, Mb, &
        PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, &
        p_sedo, Fsed_min, me2p, Mc)

      Ce  = GetStiffness(KK, GG)
      K_p = 100.0d0 * GG

      !--- back-stress initialisation — matches C++
      Alpha = Alpha_n

      Alpha_in        = 0.0d0
      Alpha_in_n      = 0.0d0
      Alpha_in_p      = 0.0d0
      Alpha_in_p_n    = 0.0d0

      Alpha_in_true   = Alpha_n
      Alpha_in_true_n = Alpha_n
      Alpha_in_max    = Alpha_n
      Alpha_in_max_n  = Alpha_n
      Alpha_in_min    = Alpha_n
      Alpha_in_min_n  = Alpha_n

      !--- fabric initialisation
      Fabric      = 0.0d0
      Fabric_n    = 0.0d0
      Fabric_in   = 0.0d0
      Fabric_in_n = 0.0d0

      !--- elastic strain
      dEpsilonE = 0.0d0

      !--- fabric scalar trackers
      zpeak   = z_max / 100000.0d0
      pzp     = max(p0, Pmin) / 100.0d0
      zxp     = 0.0d0
      pzpFlag = .TRUE.

      !--- tolerances
      if (TolF < 0.0d0) TolF = 1.0d-10
      if (TolR < 0.0d0) TolR = 1.0d-8

      !--- initial void ratio
      if (InitialVoidRatio < 0.0d0) then
        InitialVoidRatio = emax - (emax - emin) * Dr
      end if
      NextVoidRatio = 0.0d0

      !--- strain tensors reset
      Epsilon   = 0.0d0
      Epsilon_n = 0.0d0
      EpsilonE  = 0.0d0
      EpsilonE_n = 0.0d0

      !--- turn off first call
      FirstCall = .FALSE.

      end subroutine PM4SandFullConstructor

!======================================================================
! PM4SandCommitState
! Called after integration to lock in converged state.
! Matches C++ commitState() exactly.
!
! Intent rules enforced:
!   intent(in)    — quantities read but not written here
!   intent(inout) — Sigma/Alpha may be modified by me2p cap;
!                   zcum/zpeak accumulated; Fabric_n/Fabric_in_n updated
!   intent(out)   — committed _n copies, moduli, tangent
!======================================================================
      subroutine PM4SandCommitState( &
        z_max, nu, G0, mm, Mc, Md, Mb, &
        PostShake, Pmin, P_atm, Cgd, p_sedo, Fsed_min, me2p, &
        Fabric, Fabric_in, &
        Alpha_in, Alpha_in_p, Alpha_in_true, Alpha_in_max, Alpha_in_min, &
        KK, GG, Mcur, rrr, &
        Alpha_in_n, Alpha_n, Alpha_in_p_n, &
        Alpha_in_true_n, Alpha_in_max_n, Alpha_in_min_n, Sigma_n, &
        zcum, zpeak, Sigma, Alpha, Fabric_n, Fabric_in_n, &
        Ce, nn, RR, K_p, &
        Epsilon, Epsilon_n, EpsilonE, EpsilonE_n, &
        InitialVoidRatio, NextVoidRatio)

      implicit none

      !--- intent(in) ------------------------------------------------
      real(8), intent(in) :: z_max, nu, G0, mm, Mc, Md, Mb
      real(8), intent(in) :: Pmin, P_atm, Cgd, p_sedo, Fsed_min
      real(8), intent(in) :: K_p, InitialVoidRatio
      logical, intent(in) :: PostShake, me2p
      real(8), intent(in), dimension(3) :: nn, RR
      real(8), intent(in), dimension(3) :: Fabric, Fabric_in
      real(8), intent(in), dimension(3) :: Alpha_in, Alpha_in_p
      real(8), intent(in), dimension(3) :: Alpha_in_true
      real(8), intent(in), dimension(3) :: Alpha_in_max, Alpha_in_min
      real(8), intent(in), dimension(3) :: Epsilon, EpsilonE

      !--- intent(out) -----------------------------------------------
      real(8), intent(out) :: KK, GG, Mcur, NextVoidRatio
      real(8), intent(inout), dimension(3) :: rrr
      real(8), intent(out), dimension(3) :: Alpha_in_n, Alpha_n
      real(8), intent(out), dimension(3) :: Alpha_in_p_n
      real(8), intent(out), dimension(3) :: Alpha_in_true_n
      real(8), intent(out), dimension(3) :: Alpha_in_max_n
      real(8), intent(out), dimension(3) :: Alpha_in_min_n
      real(8), intent(out), dimension(3) :: Sigma_n
      real(8), intent(out), dimension(3) :: Epsilon_n, EpsilonE_n
      real(8), intent(out), dimension(3,3) :: Ce

      !--- intent(inout) ---------------------------------------------
      real(8), intent(inout) :: zcum, zpeak
      real(8), intent(inout), dimension(3) :: Sigma, Alpha
      real(8), intent(inout), dimension(3) :: Fabric_n, Fabric_in_n

      !--- local -----------------------------------------------------
      real(8), dimension(3)   :: dFabric, I1
      real(8) :: pp, Mb_loc
      I1   = (/1.0d0, 1.0d0, 0.0d0/)

      !--- recompute elastic moduli at committed stress state
      !    C++: this->GetElasticModuli(mSigma, mK, mG, mMcur, mzcum)
      !    Mb is intent(in) here, use local copy since GetElasticModuli_
      !    receives it as intent(in) too.
      Mb_loc = Mb
      call GetElasticModuli_(Sigma, zcum, z_max, nu, G0, Md, Mb_loc, &
        PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, &
        p_sedo, Fsed_min, me2p, Mc)

      !--- me2p bounding surface cap
      !    C++: if (mMcur > mMb && me2p) { ... }
      if ((Mcur > Mb_loc) .and. me2p) then
        pp    = 0.5d0 * (Sigma(1) + Sigma(2))
        rrr   = (Sigma - pp*I1) / pp * (Mb_loc/Mcur)
        Sigma = pp*I1 + rrr*pp
        Alpha = rrr * (Mb_loc - mm) / Mb_loc
      end if

      !--- commit all tensor state
      Alpha_in_n      = Alpha_in
      Alpha_n         = Alpha
      Alpha_in_p_n    = Alpha_in_p
      Alpha_in_true_n = Alpha_in_true
      Alpha_in_max_n  = Alpha_in_max
      Alpha_in_min_n  = Alpha_in_min
      Sigma_n         = Sigma
      Epsilon_n       = Epsilon
      EpsilonE_n      = EpsilonE

      !--- fabric accumulation
      !    C++: mzcum += sqrt(DoubleDot2_2_Contr(dFabric,dFabric)/2)
      !         mzpeak = max(sqrt(...), mzpeak)
      dFabric = Fabric - Fabric_n

      zcum  = zcum + sqrt(max(0.0d0, &
        (dFabric(1)**2 + dFabric(2)**2 + 2.0d0*dFabric(3)**2) &
        / 2.0d0))

      zpeak = max(zpeak, sqrt(max(0.0d0, &
        (Fabric(1)**2 + Fabric(2)**2 + 2.0d0*Fabric(3)**2) &
        / 2.0d0)))

      Fabric_n    = Fabric
      Fabric_in_n = Fabric_in

      !--- void ratio  (FIX: was commented out in original Fortran)
      !    C++: mVoidRatio = m_e_init - (1+m_e_init)*GetTrace(mEpsilon)
      !    GetTrace(Epsilon) = Epsilon(1)+Epsilon(2) in 2D
      NextVoidRatio = InitialVoidRatio &
        - (1.0d0 + InitialVoidRatio) * (Epsilon(1) + Epsilon(2))

      !--- elastic tangent
      !    C++ TangType=0 always returns mCe (elastic), not mCep
      Ce = GetStiffness(KK, GG) ! can revisit for Cep if we want to do implicit

      end subroutine PM4SandCommitState


!======================================================================
! PM4SandIntegrate
! Top-level integration driver.  Matches C++ integrate().
! Handles load-reversal detection, elastic/elasto-plastic branching.
!======================================================================
      subroutine PM4SandIntegrate( &
        Epsilon, Epsilon_n, EpsilonE, EpsilonE_n, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, pzp, zxp, &
        KK, GG, Mcur, rrr, Sigma_b, &
        Fabric_n, Fabric, Fabric_in_n, Fabric_in, &
        Alpha_n, Alpha, Alpha_in_n, Alpha_in, Alpha_in_p_n, &
        Alpha_in_p, Alpha_in_true_n, Alpha_in_true, &
        Alpha_in_max_n, Alpha_in_max, Alpha_in_min_n, Alpha_in_min, &
        pzpFlag, TolF, TolR, dEpsilonE, Sigma_n, Sigma, &
        Ce, nn, RR, NextVoidRatio)

      implicit none

      !--- intent(in) — committed state read at start
      real(8), intent(in) :: G0, hpo
      logical, intent(in) :: PostShake, me2p
      real(8), intent(in) :: P_atm, h0, emax, emin, nb, nd, Ado
      real(8), intent(in) :: z_max, cz, ceps, nu
      real(8), intent(in) :: Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm
      real(8), intent(in) :: InitialVoidRatio
      real(8), intent(in) :: Pmin, Pmin2, zpeak, zcum
      real(8), intent(in) :: TolF, TolR
      real(8), intent(in), dimension(3) :: Sigma_n
      real(8), intent(in), dimension(3) :: Fabric_n, Fabric_in_n
      real(8), intent(in), dimension(3) :: Alpha_n
      real(8), intent(in), dimension(3) :: Alpha_in_n, Alpha_in_p_n
      real(8), intent(in), dimension(3) :: Alpha_in_true_n
      real(8), intent(in), dimension(3) :: Alpha_in_max_n
      real(8), intent(in), dimension(3) :: Alpha_in_min_n
      real(8), intent(in), dimension(3) :: Epsilon_n, EpsilonE_n

      !--- intent(inout) — evolving during integration
      real(8), intent(inout) :: Mc, Mb, Md
      real(8), intent(inout) :: Cdr, Fsed_min, p_sedo, K_p
      real(8), intent(inout) :: pzp, zxp, KK, GG, Mcur
      real(8), intent(inout) :: NextVoidRatio
      real(8), intent(inout), dimension(3) :: rrr, Sigma_b
      real(8), intent(inout), dimension(3) :: Fabric, Fabric_in
      real(8), intent(inout), dimension(3) :: Alpha, Alpha_in
      real(8), intent(inout), dimension(3) :: Alpha_in_p
      real(8), intent(inout), dimension(3) :: Alpha_in_true
      real(8), intent(inout), dimension(3) :: Alpha_in_max
      real(8), intent(inout), dimension(3) :: Alpha_in_min
      real(8), intent(inout), dimension(3) :: dEpsilonE
      real(8), intent(inout), dimension(3) :: Sigma
      real(8), intent(inout), dimension(3) :: Epsilon, EpsilonE
      real(8), intent(inout), dimension(3,3) :: Ce
      logical, intent(inout) :: pzpFlag

      !--- intent(out)
      real(8), intent(out), dimension(3) :: nn, RR

      !--- local
      real(8), dimension(3) :: n_tr, StrainIncrement
      real(8), dimension(3) :: Alpha_mAlpha_in_true
      real(8) :: pp, zxpTemp
      real(8) :: elasticRatio, dVolStrain, DGamma
      integer(4) :: ii

      !--- reset to committed values (C++: mAlpha = mAlpha_n; etc.)
      Alpha         = Alpha_n
      Alpha_in      = Alpha_in_n
      Alpha_in_true = Alpha_in_true_n
      Alpha_in_p    = Alpha_in_p_n
      Alpha_in_max  = Alpha_in_max_n
      Alpha_in_min  = Alpha_in_min_n
      Fabric        = Fabric_n
      Fabric_in     = Fabric_in_n

      !--- trial elastic stress
      StrainIncrement = Epsilon - Epsilon_n
      Sigma = Sigma_n + matmul(Ce, StrainIncrement)

      !--- trial normal to yield surface
      n_tr = GetNormalToYield(Sigma, Alpha)

      !--- load-reversal detection
      !    C++: if (DoubleDot2_2_Contr(mAlpha-mAlpha_in_true, n_tr)<0 && me2p)
      Alpha_mAlpha_in_true = Alpha - Alpha_in_true
      if ((DoubleDot2_2_Contr(Alpha_mAlpha_in_true, n_tr) < 0.0d0) &
        .and. me2p) then

        Alpha_in_p    = Alpha_in
        Alpha_in_true = Alpha
        Fabric_in     = Fabric

        ! update pzp
        pp = 0.5d0*(Sigma_n(1)+Sigma_n(2))
        if (pp <= Pmin) pp = Pmin

        zxpTemp = GetNorm_Contr(Fabric_n) * pp

        if (((zxpTemp > zxp) .and. (pp > pzp)) .or. pzpFlag) then
          zxp     = zxpTemp
          pzp     = pp
          pzpFlag = .FALSE.
        end if

        ! track min/max back-stress history
        do ii = 1, 3
          if (Alpha_in(ii) > 0.0d0) then
            Alpha_in_min(ii) = min(Alpha_in_min(ii), Alpha(ii))
          else
            Alpha_in_max(ii) = max(Alpha_in_max(ii), Alpha(ii))
          end if
        end do

        ! update effective initial back-stress
        if (Alpha(3) * Alpha_in_p(3) > 0.0d0) then
          do ii = 1, 3
            if (n_tr(ii) > 0.0d0) then
              Alpha_in(ii) = max(0.0d0, Alpha_in_min(ii))
            else
              Alpha_in(ii) = min(0.0d0, Alpha_in_max(ii))
            end if
          end do
        else
          Alpha_in = Alpha
        end if

      end if

      !--- elastic-only stage (me2p == .false.)
      if (.not. me2p) then
        call elastic_integrator(Sigma_n, Epsilon_n, EpsilonE_n, &
        Epsilon, EpsilonE, Sigma, Alpha, NextVoidRatio, &
        GG, KK, Ce, G0, nu, Pmin, P_atm, InitialVoidRatio, me2p)
        nn = GetNormalToYield(Sigma, Alpha)
        RR = nn
        return
      end if

      !--- elasto-plastic stage
      call explicit_integrator(Sigma_n, Epsilon_n, EpsilonE_n, &
        Alpha_n, Fabric_n, Alpha_in, Alpha_in_p, &
        Epsilon, EpsilonE, Sigma, Alpha, Fabric, DGamma, &
        NextVoidRatio, GG, KK, Ce, &
        TolF, G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, pzp, Alpha_in_true, &
        nn, RR)

      end subroutine PM4SandIntegrate

!======================================================================
! elastic_integrator
! Pure linear-elastic stress update.  Matches C++ elastic_integrator().
!======================================================================
      subroutine elastic_integrator(CurStress, CurStrain, &
        CurElasticStrain, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, &
        NextVoidRatio, GG, KK, aC, &
        G0, nu, Pmin, P_atm, InitialVoidRatio, me2p)

      implicit none

      real(8), intent(in),  dimension(3)   :: CurStress, CurStrain
      real(8), intent(in),  dimension(3)   :: &
        CurElasticStrain, NextStrain
      real(8), intent(in)                  :: G0, nu, Pmin, P_atm
      real(8), intent(in)                  :: InitialVoidRatio
      logical, intent(in)                  :: me2p

      real(8), intent(out), dimension(3)   :: NextElasticStrain
      real(8), intent(out), dimension(3)   :: NextStress, NextAlpha
      real(8), intent(out)                 :: NextVoidRatio, GG, KK
      real(8), intent(out), dimension(3,3) :: aC

      real(8), dimension(3) :: dStrain
      real(8) :: pp

      dStrain = NextStrain - CurStrain

      ! void ratio
      NextVoidRatio = InitialVoidRatio &
        - (1.0d0 + InitialVoidRatio) * (NextStrain(1)+NextStrain(2))

      NextElasticStrain = CurElasticStrain + dStrain

      call GetElasticModuli_simple(CurStress, KK, GG, G0, nu, &
        P_atm, Pmin, me2p)

      aC = GetStiffness(KK, GG)

      NextStress = CurStress + DoubleDot4_2(aC, dStrain)

      pp = 0.5d0 * (NextStress(1) + NextStress(2))
      if (pp > Pmin) then
        NextAlpha = GetDevPart(NextStress) / pp
      end if

      end subroutine elastic_integrator


!======================================================================
! explicit_integrator
! Selects integration sub-scheme and handles elastic/plastic/transition.
! Matches C++ explicit_integrator().
!======================================================================
      subroutine explicit_integrator(CurStress, CurStrain, &
        CurElasticStrain, CurAlpha, CurFabric, Alpha_in, Alpha_in_p, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, &
        NextFabric, NextL, NextVoidRatio, GG, KK, aC, &
        TolF, G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, pzp, Alpha_in_true, &
        nn, RR)

      implicit none

      real(8), intent(in), dimension(3)   :: CurStress, CurStrain
      real(8), intent(in), dimension(3)   :: CurElasticStrain
      real(8), intent(in), dimension(3)   :: CurAlpha, CurFabric
      real(8), intent(in), dimension(3)   :: Alpha_in, Alpha_in_p
      real(8), intent(in), dimension(3)   :: NextStrain
      real(8), intent(in), dimension(3)   :: Fabric_in, Alpha_in_true
      real(8), intent(in)                 :: TolF, G0, hpo
      logical, intent(in)                 :: PostShake, me2p
      real(8), intent(in)                 :: P_atm, h0, emax, emin
      real(8), intent(in)                 :: nb, nd, Ado, z_max
      real(8), intent(in)                 :: cz, ceps, nu, Cgd, Ckaf
      real(8), intent(in)                 :: QQ_Bolton, RR_Bolton, mm
      real(8), intent(in)                 :: InitialVoidRatio
      real(8), intent(in)                 :: Pmin, Pmin2
      real(8), intent(in)                 :: zpeak, zcum

      real(8), intent(inout)              :: Mc, Mb, Md, Cdr
      real(8), intent(inout)              :: Fsed_min, p_sedo, K_p
      real(8), intent(inout)              :: GG, KK, Mcur, pzp

      real(8), intent(out), dimension(3)   :: NextElasticStrain
      real(8), intent(out), dimension(3)   :: NextStress, NextAlpha
      real(8), intent(out), dimension(3)   :: NextFabric
      real(8), intent(out), dimension(3)   :: nn, RR
      real(8), intent(out)                 :: NextL, NextVoidRatio
      real(8), intent(out), dimension(3,3) :: aC

      ! local
      real(8), dimension(3)   :: dStrain, dSigma, dDevStrain, &
        dElasStrain
      real(8), dimension(3)   :: I1
      real(8) :: f, fn, elasticRatio, dVolStrain, denom, RatioValue
      real(8) :: small

      small = 1.0d-10
      I1    = (/1.0d0, 1.0d0, 0.0d0/)

      ! void ratio
      NextVoidRatio = InitialVoidRatio &
        - (1.0d0 + InitialVoidRatio) * (NextStrain(1)+NextStrain(2))

      dStrain           = NextStrain - CurStrain
      NextElasticStrain = CurElasticStrain + dStrain
      dVolStrain        = dStrain(1) + dStrain(2)
      dDevStrain        = dStrain - I1*(dVolStrain/3.0d0)

      aC     = GetStiffness(KK, GG)
      dSigma = KK*dVolStrain*I1 + 2.0d0*GG*ToContraviant(dDevStrain)
      NextStress = CurStress + dSigma

      f  = GetFYieldFunction(NextStress, CurAlpha, mm)
      fn = GetFYieldFunction(CurStress,  CurAlpha, mm)
      nn = GetNormalToYield(NextStress, CurAlpha)

      !--- CASE 1: pure elastic
      if (f <= TolF) then
        NextAlpha  = CurAlpha
        NextFabric = CurFabric
        NextL      = 0.0d0
        RR         = nn
        return
      end if

      !--- CASE 2: elastic-to-plastic transition
      if (fn < -TolF) then
        elasticRatio = IntersectionFactor(CurStress, CurStrain, &
        NextStrain, CurAlpha, aC, 0.0d0, 1.0d0, TolF, mm)
        dElasStrain = dStrain * elasticRatio
        dSigma      = DoubleDot4_2(aC, dElasStrain)
        call MaxStrainInc(CurStress+dSigma, CurStrain+dElasStrain, &
        CurElasticStrain+dElasStrain, CurAlpha, CurFabric, &
        Alpha_in, Alpha_in_p, NextStrain, NextElasticStrain, &
        NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, &
        GG, KK, aC, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, TolF, pzp, Alpha_in_true, &
        nn, RR)
        return
      end if

      !--- CASE 3: on yield surface
      if (abs(fn) < TolF) then

        if (GetNorm_Contr(dSigma) < small) then
          denom = 1.0d0
        else
          denom = GetNorm_Contr(dSigma)
        end if

        RatioValue = DoubleDot2_2_Contr( &
        GetNormalToYield(CurStress, CurAlpha), dSigma) / denom

        if (RatioValue > -sqrt(TolF)) then
          ! pure plastic
          call MaxStrainInc(CurStress, CurStrain, CurElasticStrain, &
        CurAlpha, CurFabric, Alpha_in, Alpha_in_p, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, &
        NextFabric, NextL, NextVoidRatio, GG, KK, aC, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, TolF, pzp, Alpha_in_true, &
        nn, RR)
        else
          ! elastic unloading then plastic
          elasticRatio = IntersectionFactor_Unloading(CurStress, &
        CurStrain, NextStrain, CurAlpha, aC, TolF, mm)
          dElasStrain = dStrain * elasticRatio
          dSigma      = DoubleDot4_2(aC, dElasStrain)
          call MaxStrainInc(CurStress+dSigma, CurStrain+dElasStrain, &
        CurElasticStrain+dElasStrain, CurAlpha, CurFabric, &
        Alpha_in, Alpha_in_p, NextStrain, NextElasticStrain, &
        NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, &
        GG, KK, aC, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, TolF, pzp, Alpha_in_true, &
        nn, RR)
        end if
        return
      end if

      !--- CASE 4: illegal stress state — recover
      !print *, "PM4Sand: illegal stress state, f=", &
      !  GetFYieldFunction(CurStress, CurAlpha, mm)
      call MaxStrainInc(CurStress, CurStrain, CurElasticStrain, &
        CurAlpha, CurFabric, Alpha_in, Alpha_in_p, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, &
        NextFabric, NextL, NextVoidRatio, GG, KK, aC, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, TolF, pzp, Alpha_in_true, &
        nn, RR)

      end subroutine explicit_integrator


!======================================================================
! MaxStrainInc
! Sub-steps strain increment if too large, then calls ModifiedEuler.
! Matches C++ MaxStrainInc() with INT_MAXSTR_ME.
!======================================================================
      subroutine MaxStrainInc(CurStress, CurStrain, CurElasticStrain, &
        CurAlpha, CurFabric, Alpha_in, Alpha_in_p, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, &
        NextFabric, NextL, NextVoidRatio, GG, KK, aC, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, TolF, pzp, Alpha_in_true, &
        nn, RR)

      implicit none

      real(8), intent(in), dimension(3)    :: CurStress, CurStrain
      real(8), intent(in), dimension(3)    :: CurElasticStrain
      real(8), intent(in), dimension(3)    :: CurAlpha, CurFabric
      real(8), intent(in), dimension(3)    :: Alpha_in, Alpha_in_p
      real(8), intent(in), dimension(3)    :: NextStrain
      real(8), intent(in), dimension(3)    :: Fabric_in, Alpha_in_true
      real(8), intent(in)                  :: TolF, G0, hpo
      logical, intent(in)                  :: PostShake, me2p
      real(8), intent(in)                  :: P_atm, h0, emax, emin
      real(8), intent(in)                  :: nb, nd, Ado, z_max
      real(8), intent(in)                  :: cz, ceps, nu, Cgd, Ckaf
      real(8), intent(in)                  :: QQ_Bolton, RR_Bolton, mm
      real(8), intent(in)                  :: InitialVoidRatio
      real(8), intent(in)                  :: Pmin, Pmin2, zpeak, zcum

      real(8), intent(inout)               :: Mc, Mb, Md, Cdr
      real(8), intent(inout)               :: Fsed_min, p_sedo, K_p
      real(8), intent(inout)               :: GG, KK, Mcur, pzp

      real(8), intent(inout), dimension(3)    :: NextElasticStrain
      real(8), intent(inout), dimension(3)    :: NextStress
      real(8), intent(out), dimension(3)      :: NextAlpha, NextFabric
      real(8), intent(out), dimension(3)      :: nn, RR
      real(8), intent(out)                    :: NextL, NextVoidRatio
      real(8), intent(out), dimension(3,3)    :: aC

      ! local
      real(8), dimension(3)   :: StrainInc, cStress, cStrain, cAlpha
      real(8), dimension(3)   :: cFabric, cEStrain, nStrain
      real(8), dimension(3,3) :: nCe
      real(8) :: maxInc, maxStrainIncValue
      real(8) :: nL, nVoid, nG, nK
      integer(4) :: ii, numSteps

      maxStrainIncValue = 1.0d-6
      
      StrainInc = NextStrain - CurStrain
      maxInc    = StrainInc(1)
      do ii = 2, 3
        if (abs(StrainInc(ii)) > abs(maxInc)) maxInc = StrainInc(ii)
      end do

      if (abs(maxInc) > maxStrainIncValue) then

        numSteps  = floor(abs(maxInc) / maxStrainIncValue) + 1
        StrainInc = (NextStrain - CurStrain) / real(numSteps, 8)

        cStress  = CurStress;    cStrain  = CurStrain
        cAlpha   = CurAlpha;     cFabric  = CurFabric
        cEStrain = CurElasticStrain

        do ii = 1, numSteps
          nStrain = cStrain + StrainInc
          call ModifiedEuler(cStress, cStrain, cEStrain, &
        cAlpha, cFabric, Alpha_in, Alpha_in_p, &
        nStrain, NextElasticStrain, NextStress, NextAlpha, &
        NextFabric, NextL, nVoid, nG, nK, nCe, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, TolF, pzp, Alpha_in_true, &
        nn, RR)
          cStress  = NextStress;  cStrain  = nStrain
          cEStrain = NextElasticStrain
          cAlpha   = NextAlpha;   cFabric  = NextFabric
        end do

        NextVoidRatio = nVoid; GG = nG; KK = nK; aC = nCe

      else

        call ModifiedEuler(CurStress, CurStrain, CurElasticStrain, &
        CurAlpha, CurFabric, Alpha_in, Alpha_in_p, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, &
        NextFabric, NextL, NextVoidRatio, GG, KK, aC, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, TolF, pzp, Alpha_in_true, &
        nn, RR)

      end if

      end subroutine MaxStrainInc


!======================================================================
! ModifiedEuler
! Adaptive Modified Euler integrator with error control.
! Matches C++ ModifiedEuler().
!======================================================================
      subroutine ModifiedEuler(CurStress, CurStrain, CurElasticStrain, &
        CurAlpha, CurFabric, Alpha_in, Alpha_in_p, &
        NextStrain, NextElasticStrain, NextStress, NextAlpha, &
        NextFabric, NextL, NextVoidRatio, GG, KK, aC, &
        G0, hpo, PostShake, me2p, &
        P_atm, h0, emax, emin, nb, nd, Ado, z_max, cz, ceps, nu, &
        Cgd, Ckaf, QQ_Bolton, RR_Bolton, mm, InitialVoidRatio, &
        Mc, Mb, Md, Cdr, Fsed_min, p_sedo, Pmin, Pmin2, K_p, &
        zpeak, zcum, Fabric_in, Mcur, TolF, pzp, Alpha_in_true, &
        nn, RR)

      implicit none

      real(8), intent(in), dimension(3)    :: CurStress, CurStrain
      real(8), intent(in), dimension(3)    :: CurElasticStrain
      real(8), intent(in), dimension(3)    :: CurAlpha, CurFabric
      real(8), intent(in), dimension(3)    :: Alpha_in, Alpha_in_p
      real(8), intent(in), dimension(3)    :: NextStrain
      real(8), intent(in), dimension(3)    :: Fabric_in, Alpha_in_true
      real(8), intent(in)                  :: TolF, G0, hpo
      logical, intent(in)                  :: PostShake, me2p
      real(8), intent(in)                  :: P_atm, h0, emax, emin
      real(8), intent(in)                  :: nb, nd, Ado, z_max
      real(8), intent(in)                  :: cz, ceps, nu, Cgd, Ckaf
      real(8), intent(in)                  :: QQ_Bolton, RR_Bolton, mm
      real(8), intent(in)                  :: InitialVoidRatio
      real(8), intent(in)                  :: Pmin, Pmin2, zpeak, zcum

      real(8), intent(inout)               :: Mc, Mb, Md, Cdr
      real(8), intent(inout)               :: Fsed_min, p_sedo
      real(8), intent(inout)               :: K_p, GG, KK, Mcur, pzp

      real(8), intent(inout), dimension(3)    :: NextElasticStrain
      real(8), intent(inout), dimension(3)    :: NextStress
      real(8), intent(out), dimension(3)      :: NextAlpha, NextFabric
      real(8), intent(out), dimension(3)      :: nn, RR
      real(8), intent(out)                    :: NextL, NextVoidRatio
      real(8), intent(out), dimension(3,3)    :: aC

      ! local
      real(8), dimension(3)   :: I1, dDevStrain, tmp0, tmp1, tmp2
      real(8), dimension(3)   :: nStress, nAlpha, nFabric
      real(8), dimension(3)   :: dSigma1, dSigma2, dAlpha1, dAlpha2
      real(8), dimension(3)   :: &
        dFabric1, dFabric2, dPStrain1, dPStrain2
      real(8), dimension(3)   :: R1, R2, alphaB, alphaD, bb
      real(8), dimension(3)   :: alphaD_NextAlpha, rrr_loc
      real(8) :: TT, dT, dT_min, TolE
      real(8) :: pp, NextDr, dVolStrain, temp4, small
      real(8) :: Cka, hh, AlphaAlphaBDotN, DD, ksi
      real(8) :: stressNorm, curStepError, qq
      real(8) :: two3, DGamma
      real(8) :: Mb_loc, Md_loc   ! local copies — Mb/Md intent(inout)

      I1     = (/1.0d0, 1.0d0, 0.0d0/)
      TT     = 0.0d0;  dT = 1.0d0;  dT_min = 1.0d-4;  TolE = 1.0d-5
      small  = 1.0d-10;  two3 = 2.0d0/3.0d0

      NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain)
      NextStress  = CurStress
      NextAlpha   = CurAlpha
      NextFabric  = CurFabric

      ! Use local Mb/Md so GetElasticModuli_ output does not
      ! overwrite the intent(inout) Mb/Md from the caller
      Mb_loc = Mb
      Md_loc = Md
      call GetElasticModuli_(NextStress, zcum, z_max, nu, G0, Md_loc, &
        Mb_loc, PostShake, Pmin, P_atm, KK, GG, Mcur, Cgd, &
        p_sedo, Fsed_min, me2p, Mc)
      Mb = Mb_loc
      Md = Md_loc
      aC = GetStiffness(KK, GG)

      pp = 0.5d0*(CurStress(1)+CurStress(2))
      if (pp < Pmin/5.0d0) then
        NextStress = GetDevPart(NextStress) + (Pmin/5.0d0)*I1
      end if

      do while (TT < 1.0d0)

        tmp0 = CurStrain + TT*(NextStrain-CurStrain)
        NextVoidRatio = InitialVoidRatio &
        - (1.0d0+InitialVoidRatio)*(tmp0(1)+tmp0(2))
        NextDr = (emax - NextVoidRatio) / (emax - emin)

        tmp0       = (NextStrain - CurStrain)
        dVolStrain = dT * (tmp0(1)+tmp0(2))
        dDevStrain = I1*(-dVolStrain/3.0d0) + tmp0*dT
        pp         = 0.5d0*(NextStress(1)+NextStress(2))

        !--- Delta 1
        call GetStateDependent(NextStress, NextAlpha, Alpha_in, &
        Alpha_in_p, Alpha_in_true, NextFabric, Fabric_in, &
        GG, zcum, zpeak, pzp, Mcur, NextDr, &
        Mc, nd, nb, Pmin, Pmin2, P_atm, mm, z_max, h0, Ckaf, &
        Ado, ceps, hpo, Cdr, QQ_Bolton, RR_Bolton, &
        nn, alphaB, alphaD, bb, R1, &
        Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, ksi)

        rrr_loc = GetDevPart(NextStress) / pp
        temp4 = K_p + 2.0d0*GG - KK*DD*DoubleDot2_2_Contr(nn,rrr_loc)

        dFabric1  = 0.0d0
        dPStrain1 = 0.0d0

        if (abs(temp4) < small) then
          dSigma1 = 0.0d0;  dAlpha1 = 0.0d0
          dPStrain1 = tmp0*dT
        else
          NextL = (2.0d0*GG*DoubleDot2_2_Mixed(nn,dDevStrain) &
        - DoubleDot2_2_Contr(nn,rrr_loc)*KK*dVolStrain) / temp4
          DGamma = NextL
          if (NextL < 0.0d0) then
            dSigma1 = 2.0d0*GG*ToContraviant(dDevStrain) &
        + KK*dVolStrain*I1
            dAlpha1 = 0.0d0
          else
            dSigma1 = 2.0d0*GG*ToContraviant(dDevStrain) &
        + KK*dVolStrain*I1 &
        - Macauley(NextL)*(2.0d0*GG*nn + KK*DD*I1)
            alphaD_NextAlpha = alphaD - NextAlpha
            if (DoubleDot2_2_Contr(alphaD_NextAlpha,nn) < 0.0d0) then
              dFabric1 = (nn*z_max + NextFabric) * &
        (-1.0d0*cz/(1.0d0+Macauley(zcum/(2.0d0*z_max)-1.0d0)) &
        *Macauley(NextL)*MacauleyIndex(-DD))
            end if
            dPStrain1 = ToCovariant(R1)*NextL
            dAlpha1   = bb*(two3*NextL*hh)
          end if
        end if

        !--- check p after Delta 1
        tmp0 = NextStress + dSigma1
        pp   = 0.5d0*(tmp0(1)+tmp0(2))
        if (pp < 0.0d0) then
          if (dT <= dT_min) then
            NextElasticStrain = CurElasticStrain+(NextStrain-CurStrain)
            NextStress = CurStress;  NextAlpha = CurAlpha
            NextFabric = CurFabric;  return
          end if
          dT = max(0.1d0*dT, dT_min);  cycle
        end if

        !--- Delta 2
        tmp1 = NextAlpha  + dAlpha1
        tmp2 = NextFabric + dFabric1
        call GetStateDependent(tmp0, tmp1, Alpha_in, &
        Alpha_in_p, Alpha_in_true, tmp2, Fabric_in, &
        GG, zcum, zpeak, pzp, Mcur, NextDr, &
        Mc, nd, nb, Pmin, Pmin2, P_atm, mm, z_max, h0, Ckaf, &
        Ado, ceps, hpo, Cdr, QQ_Bolton, RR_Bolton, &
        nn, alphaB, alphaD, bb, R2, &
        Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, ksi)

        rrr_loc = GetDevPart(tmp0) / pp
        temp4 = K_p + 2.0d0*GG - KK*DD*DoubleDot2_2_Contr(nn,rrr_loc)

        dFabric2  = 0.0d0
        dPStrain2 = 0.0d0

        if (abs(temp4) < small) then
          dSigma2 = 0.0d0;  dAlpha2 = 0.0d0
          dPStrain2 = dPStrain1
        else
          NextL = (2.0d0*GG*DoubleDot2_2_Mixed(nn,dDevStrain) &
        - DoubleDot2_2_Contr(nn,rrr_loc)*KK*dVolStrain) / temp4
          DGamma = NextL
          if (NextL < 0.0d0) then
            dSigma2 = 2.0d0*GG*ToContraviant(dDevStrain) &
        + KK*dVolStrain*I1
            dAlpha2 = 0.0d0
          else
            dSigma2 = 2.0d0*GG*ToContraviant(dDevStrain) &
        + KK*dVolStrain*I1 &
        - Macauley(NextL)*(2.0d0*GG*nn + KK*DD*I1)
            alphaD_NextAlpha = alphaD - NextAlpha - dAlpha1
            if (DoubleDot2_2_Contr(alphaD_NextAlpha,nn) < 0.0d0) then
              dFabric2 = (nn*z_max + NextFabric + dFabric1) * &
        (-1.0d0*cz/(1.0d0+Macauley(zcum/(2.0d0*z_max)-1.0d0)) &
        *Macauley(NextL)*MacauleyIndex(-DD))
            end if
            dPStrain2 = ToCovariant(R2)*NextL
            dAlpha2   = bb*(two3*NextL*hh)
          end if
        end if

        !--- average update
        nStress = NextStress + 0.5d0*(dSigma1+dSigma2)
        nAlpha  = NextAlpha  + 0.5d0*(dAlpha1+dAlpha2)
        nFabric = NextFabric + 0.5d0*(dFabric1+dFabric2)

        pp = 0.5d0*(nStress(1)+nStress(2))
        if (pp < 0.0d0) then
          if (dT <= dT_min) then
            NextElasticStrain = CurElasticStrain+(NextStrain-CurStrain)
            NextStress = CurStress;  NextAlpha = CurAlpha
            NextFabric = CurFabric;  return
          end if
          dT = max(0.1d0*dT, dT_min);  cycle
        end if

        !--- error control
        stressNorm = GetNorm_Contr(NextStress)
        tmp0 = dSigma2 - dSigma1
        if (stressNorm < 0.5d0) then
          curStepError = GetNorm_Contr(tmp0)
        else
          curStepError = GetNorm_Contr(tmp0) / (2.0d0*stressNorm)
        end if

        if (curStepError > TolE) then
            
            if (curStepError == 0) then 
                curStepError = 1e-15
            end if 
            
          qq = max(0.8d0*sqrt(TolE/curStepError), 0.1d0)
          if (dT <= dT_min) then
            NextElasticStrain = NextElasticStrain &
        - 0.5d0*(dPStrain1+dPStrain2)
            NextStress = nStress;  NextAlpha = nAlpha
            call stress_correction_(NextStress, NextAlpha, Alpha_in, &
        Alpha_in_p, CurFabric, NextVoidRatio, Alpha_in_true, &
        Fabric_in, CurStress, Pmin, P_atm, TolF, Mc, emax, emin, &
        zcum, zpeak, z_max, pzp, mm, h0, hpo, Cdr, ceps, Ckaf, &
        QQ_Bolton, RR_Bolton, DGamma, alphaD, alphaB, &
        DD, Mb, Md, K_p, Mcur, GG, KK, dSigma1, nb, nd, Pmin2, &
        Ado, nn, RR)
            TT = TT + dT
          end if
          dT = max(qq*dT, dT_min)
        else
          NextElasticStrain = NextElasticStrain &
        - 0.5d0*(dPStrain1+dPStrain2)
          NextStress  = nStress
          NextAlpha   = nAlpha
          NextFabric  = nFabric
          call stress_correction_(NextStress, NextAlpha, Alpha_in, &
        Alpha_in_p, CurFabric, NextVoidRatio, Alpha_in_true, &
        Fabric_in, CurStress, Pmin, P_atm, TolF, Mc, emax, emin, &
        zcum, zpeak, z_max, pzp, mm, h0, hpo, Cdr, ceps, Ckaf, &
        QQ_Bolton, RR_Bolton, DGamma, alphaD, alphaB, &
        DD, Mb, Md, K_p, Mcur, GG, KK, dSigma1, nb, nd, Pmin2, &
        Ado, nn, RR)
          TT = TT + dT
          
          if (curStepError == 0) then 
                curStepError = 1e-15
            end if 
          
          qq = max(0.8d0*sqrt(TolE/curStepError), 0.5d0)
          dT = max(qq*dT, dT_min)
          dT = min(dT, 1.0d0-TT)
        end if

      end do

      end subroutine ModifiedEuler

!======================================================================
! stress_correction_
! Pulls stress back onto yield surface.
! Matches C++ Stress_Correction() (first overload).
!======================================================================
      subroutine stress_correction_(NextStress, NextAlpha, Alpha_in, &
        Alpha_in_p, CurFabric, NextVoidRatio, Alpha_in_true, &
        Fabric_in, CurStress, Pmin, P_atm, TolF, Mc, emax, emin, &
        zcum, zpeak, z_max, pzp, mm, h0, hpo, Cdr, Ceps, Ckaf, &
        QQ_Bolton, RR_Bolton, DGamma, alphaD_in, alphaB_in, &
        DD_in, Mb, Md, K_p, Mcur, GG, KK, dSigma, &
        nb, nd, Pmin2, Ado, nn, RR)

      implicit none

      real(8), intent(in), dimension(3) :: Alpha_in, Alpha_in_p
      real(8), intent(in), dimension(3) :: Alpha_in_true, Fabric_in
      real(8), intent(in), dimension(3) :: CurFabric, CurStress
      real(8), intent(in) :: Pmin, P_atm, TolF, Mc, emax, emin
      real(8), intent(in) :: zcum, zpeak, z_max, pzp, mm
      real(8), intent(in) :: h0, hpo, Cdr, Ceps, Ckaf
      real(8), intent(in) :: QQ_Bolton, RR_Bolton, DGamma
      real(8), intent(in) :: nb, nd, Pmin2, Ado
      real(8), intent(in), dimension(3) :: alphaD_in, alphaB_in
      real(8), intent(in) :: DD_in

      real(8), intent(inout), dimension(3)   :: NextStress, NextAlpha
      real(8), intent(inout), dimension(3)   :: dSigma
      real(8), intent(inout)                 :: NextVoidRatio
      real(8), intent(inout)                 :: Mcur, GG, KK
      real(8), intent(inout)                 :: Mb, Md, K_p

      real(8), intent(out), dimension(3) :: nn, RR

      ! local
      real(8), dimension(3)   :: nStress, nAlpha, dSigmaP, dfrOverdSigma
      real(8), dimension(3)   :: dfrOverdAlpha, aBar, rrr_loc, bb
      real(8), dimension(3)   :: alphaD, alphaB, tmp0, tmp1
      real(8), dimension(3)   :: I1
      real(8), dimension(3,3) :: Ce
      real(8) :: pp, fr, lambda, CurDr, Cka, hh, ksi, AlphaAlphaBDotN
      real(8) :: DD, alpha_up, alpha_mid, alpha_down, fr_old, two3
      integer(4) :: ii, jj, maxIter

      I1     = (/1.0d0, 1.0d0, 0.0d0/)
      maxIter = 200!25
      two3   = 2.0d0/3.0d0

      pp = 0.5d0*(NextStress(1)+NextStress(2))

      !--- tension cutoff
      if (pp < Pmin/5.0d0) then
        fr = GetFYieldFunction(NextStress, NextAlpha, mm)
        if (fr < TolF) then
          NextStress = NextStress + (Pmin/5.0d0 - pp)*I1
        else
          NextStress    = (Pmin/5.0d0)*I1
          NextStress(3) = 0.8d0*Mc*Pmin/5.0d0
          NextAlpha     = 0.0d0
          NextAlpha(3)  = 0.8d0*Mc
          nn = GetNormalToYield(NextStress, NextAlpha)
          RR = nn
        end if
        return
      end if

      fr = GetFYieldFunction(NextStress, NextAlpha, mm)
      if (fr < TolF) then
        nn = GetNormalToYield(NextStress, NextAlpha)
        RR = nn
        return
      end if

      ! correction loop  (C++ maxIter = 25)
      CurDr  = (emax - NextVoidRatio) / (emax - emin)
      nStress = NextStress
      nAlpha  = NextAlpha

      do ii = 1, maxIter
        pp      = max(0.5d0*(nStress(1)+nStress(2)), Pmin)
        rrr_loc = GetDevPart(nStress) / pp

        call GetStateDependent(nStress, nAlpha, Alpha_in, &
        Alpha_in_p, Alpha_in_true, CurFabric, Fabric_in, &
        GG, zcum, zpeak, pzp, Mcur, CurDr, &
        Mc, nd, nb, Pmin, Pmin2, P_atm, mm, z_max, h0, Ckaf, &
        Ado, Ceps, hpo, Cdr, QQ_Bolton, RR_Bolton, &
        nn, alphaB, alphaD, bb, RR, &
        Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, ksi)

        Ce      = GetStiffness(KK, GG)
        dSigmaP = DoubleDot4_2(Ce, DGamma*ToCovariant(RR))
        aBar    = two3*hh*bb

        dfrOverdSigma = -0.5d0*DoubleDot2_2_Contr(nn,rrr_loc)*I1 + nn
        dfrOverdAlpha = -pp*nn

        lambda = fr / (DoubleDot2_2_Contr(dfrOverdSigma, dSigmaP) &
        -  DoubleDot2_2_Contr(dfrOverdAlpha,  aBar))

        tmp0 = nStress - lambda*dSigmaP
        tmp1 = nAlpha  + lambda*aBar

        if (abs(GetFYieldFunction(tmp0,tmp1,mm)) < abs(fr)) then
          nStress = nStress - lambda*dSigmaP
          nAlpha  = nAlpha  + lambda*aBar
        else
          lambda  = fr / DoubleDot2_2_Contr(dfrOverdSigma,dfrOverdSigma)
          nStress = nStress - lambda*dfrOverdSigma
        end if

        fr = GetFYieldFunction(nStress, nAlpha, mm)
        if (abs(fr) < TolF) then
          NextStress = nStress
          NextAlpha  = nAlpha
          nn = GetNormalToYield(NextStress, NextAlpha)
          RR = nn
          return
        end if
      end do

      ! bisection fallback on stress path from CurStress to NextStress
      ! C++: dSigma = NextStress - mSigma  (mSigma = CurStress here)
      dSigma     = NextStress - CurStress
      alpha_up   = 1.0d0
      alpha_mid  = 0.5d0
      alpha_down = 0.0d0
      tmp0       = alpha_mid*dSigma + CurStress
      fr_old     = GetFYieldFunction(tmp0, NextAlpha, mm)

      do jj = 1, maxIter
        if (fr_old < 0.0d0) then
          alpha_down = alpha_mid
          alpha_mid  = 0.5d0*(alpha_up + alpha_mid)
        else
          alpha_up   = alpha_mid
          alpha_mid  = 0.5d0*(alpha_down + alpha_mid)
        end if
        tmp0   = alpha_mid*dSigma + CurStress
        fr_old = GetFYieldFunction(tmp0, NextAlpha, mm)
        if (abs(fr_old) < TolF) then
          NextStress = tmp0
          exit
        end if
      end do

      nn = GetNormalToYield(NextStress, NextAlpha)
      RR = nn

      end subroutine stress_correction_


!======================================================================
! GetElasticModuli_
! Full fabric-dependent moduli including PostShake reduction.
! FIX: PostShake block was commented out in original Fortran.
! FIX: nu is intent(in) — use local copy nu_loc.
!======================================================================
      subroutine GetElasticModuli_(sigma, zcum, z_max, nu, G0, &
        Md, Mb, PostShake, Pmin, P_atm, &
        KK, GG, Mcur, Cgd, p_sedo, Fsed_min, me2p, Mc)

      implicit none

      real(8), intent(in), dimension(3) :: sigma
      real(8), intent(in)  :: zcum, z_max, nu, G0
      real(8), intent(in)  :: Md, Mb, Pmin, P_atm, Cgd
      real(8), intent(in)  :: p_sedo, Fsed_min, Mc
      logical, intent(in)  :: PostShake, me2p

      real(8), intent(out) :: KK, GG, Mcur

      real(8) :: pn, qn, Csr, temp, p_sed, F_sed
      real(8) :: nu_loc, two3, Csr0, small
      integer :: msr

      two3  = 2.0d0/3.0d0
      Csr0  = 0.5d0;  small = 1.0d-10;  msr = 4

      pn = 0.5d0*(sigma(1)+sigma(2))
      if (pn <= Pmin) pn = Pmin

      qn   = 2.0d0*sqrt((0.5d0*(sigma(1)-sigma(2)))**2 + sigma(3)**2)
      Mcur = qn / pn

      Csr = 1.0d0 - Csr0*min(1.0d0, (Mcur/Mb)**msr)
      temp = zcum / z_max

      if (.not. me2p) then
        GG = G0 * P_atm
      else
        GG = G0*P_atm*sqrt(pn/P_atm)*Csr &
        *(1.0d0+temp)/(1.0d0+temp*Cgd)

        !--- PostShake modulus reduction (FIX: restored from C++)
        if (PostShake) then
          pn    = 0.5d0*(sigma(1)+sigma(2))   ! without Pmin cutoff
          p_sed = p_sedo*(zcum/(zcum+z_max)) &
        * max(0.0d0, 1.0d0-Mcur/Md)**0.25d0
          F_sed = min(Fsed_min &
        + (1.0d0-Fsed_min)*pn/(20.0d0*(p_sed+small)), &
        1.0d0)
          GG = GG * F_sed
        end if
      end if

      ! FIX: nu is intent(in) — use local copy for incompressible guard
      nu_loc = nu
      if (nu_loc == 0.5d0) nu_loc = 0.4999d0
      KK = two3*((1.0d0+nu_loc)/(1.0d0-2.0d0*nu_loc))*GG

      end subroutine GetElasticModuli_


!----------------------------------------------------------------------
! GetElasticModuli_simple
! Fabric-independent moduli used in elastic_integrator.
!----------------------------------------------------------------------
      subroutine GetElasticModuli_simple(sigma, KK, GG, G0, nu, &
        P_atm, Pmin, me2p)

      implicit none

      real(8), intent(in), dimension(3) :: sigma
      real(8), intent(in)  :: G0, nu, P_atm, Pmin
      logical, intent(in)  :: me2p
      real(8), intent(out) :: KK, GG

      real(8) :: pn, two3, nu_loc

      two3 = 2.0d0/3.0d0
      pn   = max(0.5d0*(sigma(1)+sigma(2)), Pmin)

      if (.not. me2p) then
        GG = G0 * P_atm
      else
        GG = G0 * P_atm * sqrt(pn/P_atm)
      end if

      nu_loc = nu
      if (nu_loc == 0.5d0) nu_loc = 0.4999d0
      KK = two3*((1.0d0+nu_loc)/(1.0d0-2.0d0*nu_loc))*GG

      end subroutine GetElasticModuli_simple


!======================================================================
! GetStateDependent
! Computes all state-dependent quantities needed for integration.
! Matches C++ GetStateDependent() exactly.
!
! FIXES vs original Fortran:
!   - Contraction threshold: 0.16 (was 0.20)
!   - C_pmin2 applied in contraction branch only (matches C++)
!   - Mb/Md are intent(out): recomputed every call from ksi
!======================================================================
      subroutine GetStateDependent(stress, alpha, alpha_in, &
        alpha_in_p, alpha_in_true, fabric, fabric_in, &
        GG, zcum, zpeak, pzp, Mcur, CurDr, &
        Mc, nd, nb, Pmin, Pmin2, P_atm, &
        mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr, &
        QQ_Bolton, RR_Bolton, &
        nn, alphaB, alphaD, bb, RR, &
        Mb, Md, DD, K_p, Cka, hh, AlphaAlphaBDotN, ksi)

      implicit none

      real(8), intent(in), dimension(3) :: stress, alpha
      real(8), intent(in), dimension(3) :: alpha_in, alpha_in_p
      real(8), intent(in), dimension(3) :: alpha_in_true
      real(8), intent(in), dimension(3) :: fabric, fabric_in
      real(8), intent(in) :: GG, zcum, zpeak, pzp, Mcur, CurDr
      real(8), intent(in) :: Mc, nd, nb, Pmin, Pmin2, P_atm
      real(8), intent(in) :: mm, z_max, h0, Ckaf, Ado, Ceps, hpo, Cdr
      real(8), intent(in) :: QQ_Bolton, RR_Bolton

      real(8), intent(out), dimension(3) :: nn, alphaB, alphaD, bb, RR
      real(8), intent(out) :: ksi, Mb, Md, DD, K_p, Cka, hh
      real(8), intent(out) :: AlphaAlphaBDotN

      ! local
      real(8), dimension(3) :: alphaD_alpha, alphaDr_alpha
      real(8), dimension(3) :: alpha_mAlpha_in, alpha_mAlpha_in_true
      real(8), dimension(3) :: alpha_mAlpha_p, minusFabric, I1
      real(8) :: pp, Czpk1, Czpk2, Cpzp2, Cg1, Ckp
      real(8) :: AlphaAlphaInDotN, AlphaAlphaInTrueDotN
      real(8) :: Czin1, Crot1, Mdr, C_pmin2
      real(8) :: Cpzp, Cpmin, Czin2, temp, Ad, Drot, Crot2
      real(8) :: hp, Cdz, Adc, Cin
      real(8) :: root12, two3, one3, small, arg

      root12 = 0.7071067811865475d0
      two3   = 2.0d0/3.0d0
      one3   = 1.0d0/3.0d0
      small  = 1.0d-10
      I1     = (/1.0d0, 1.0d0, 0.0d0/)

      pp = 0.5d0*(stress(1)+stress(2))
      if (pp <= Pmin) pp = Pmin

      ksi = GetKsi(CurDr, pp, RR_Bolton, QQ_Bolton, Pmin, P_atm)
      nn  = GetNormalToYield(stress, alpha)

      !--- Mb and Md (recomputed every call from ksi)
      if (ksi <= 0.0d0) then
        arg = min(-nb*ksi, 700.0d0)
        Mb  = Mc * exp(arg)
        arg = max(nd*ksi, -700.0d0)
        Md  = Mc * exp(arg)
      else
        arg = max(-(nb/4.0d0)*ksi, -700.0d0)
        Mb  = Mc * exp(arg)
        arg = min((nd*4.0d0)*ksi, 700.0d0)
        Md  = Mc * exp(arg)
      end if

      alphaB = nn * (root12*(Mb - mm))
      alphaD = nn * (root12*(Md - mm))

      !--- plastic modulus terms
      Czpk1 = zpeak / (zcum + z_max/5.0d0)
      Czpk2 = zpeak / (zcum + z_max/100.0d0)
      if (Czpk2 > 1.0d0-small) Czpk2 = 1.0d0-small

      Cpzp2 = Macauley(pzp-pp) / (Macauley(pzp-pp) + Pmin)

      Cg1 = h0 / 200.0d0
      Ckp = 2.0d0

      bb = alphaB - alpha
      AlphaAlphaBDotN = DoubleDot2_2_Contr(bb, nn)

      alpha_mAlpha_in = alpha - alpha_in
      AlphaAlphaInDotN = &
        Macauley(DoubleDot2_2_Contr(alpha_mAlpha_in,nn))

      alpha_mAlpha_in_true = alpha - alpha_in_true
      AlphaAlphaInTrueDotN = &
        Macauley(DoubleDot2_2_Contr(alpha_mAlpha_in_true, nn))

      Cka = 1.0d0 + &
        (Ckaf/(1.0d0+(2.5d0*AlphaAlphaInTrueDotN)**2)) &
        * Cpzp2*Czpk1

      alpha_mAlpha_p = alpha - alpha_in_p

      if (abs(AlphaAlphaBDotN) < small) then
        hh = 1.0d10
      else if (DoubleDot2_2_Contr(alpha_mAlpha_p, nn) <= 0.0d0) then
        ! load reversal
        hh = 1.5d0*GG*h0/pp &
        / (exp(AlphaAlphaInDotN)-1.0d0+Cg1) &
        / sqrt(abs(AlphaAlphaBDotN)) &
        * Cka / (1.0d0+Ckp*(zpeak/z_max) &
        *Macauley(AlphaAlphaBDotN)*sqrt(1.0d0-Czpk2))
        hh = hh*(AlphaAlphaInDotN+Cg1)/(AlphaAlphaInTrueDotN+Cg1)
      else
        ! no reversal
        hh = 1.5d0*GG*h0/pp &
        / (exp(AlphaAlphaInDotN)-1.0d0+Cg1) &
        / sqrt(abs(AlphaAlphaBDotN)) &
        * Cka / (1.0d0+Ckp*(zpeak/z_max) &
        *Macauley(AlphaAlphaBDotN)*sqrt(1.0d0-Czpk2))
      end if

      K_p = two3*hh*pp*DoubleDot2_2_Contr(bb, nn)

      !--- dilatancy
      Czin1 = Macauley(1.0d0 - exp(-2.0d0*abs( &
        (DoubleDot2_2_Contr(fabric_in,nn) &
        -DoubleDot2_2_Contr(fabric,nn)) / z_max)))

      minusFabric = -1.0d0*fabric
      Crot1 = max(1.0d0 + &
        2.0d0*Macauley(DoubleDot2_2_Contr(minusFabric,nn)) &
        /(sqrt(2.0d0)*z_max)*(1.0d0-Czin1), 1.0d0)

      Mdr           = Md / Crot1
      alphaDr_alpha = nn*(root12*(Mdr-mm)) - alpha
      alphaD_alpha  = alphaD - alpha

      

      if (DoubleDot2_2_Contr(alphaDr_alpha, nn) <= 0.0d0) then

        !--- DILATION
        if (pzp == 0.0d0) then
          Cpzp = 1.0d0
        else
          Cpzp = 1.0d0/(1.0d0+(2.5d0*pp/pzp)**5)
        end if
        Cpmin = 1.0d0/(1.0d0+(Pmin2/pp)**2)
        Czin2 = (1.0d0+Czin1*(zcum-zpeak)/(3.0d0*z_max)) &
        / (1.0d0+3.0d0*Czin1*(zcum-zpeak)/(3.0d0*z_max))

        temp = (1.0d0 - Macauley(DoubleDot2_2_Contr(minusFabric,nn)) &
        *root12/zpeak)**3

        Ad = Ado*Czin2 / &
        (zcum**2/z_max*temp*Ceps**2*Cpzp*Cpmin*Czin1 + 1.0d0)

        DD = Ad * DoubleDot2_2_Contr(alphaD_alpha, nn)

        Drot = Ad * Macauley(DoubleDot2_2_Contr(minusFabric,nn)) &
        / (sqrt(2.0d0)*z_max) &
        * DoubleDot2_2_Contr(alphaDr_alpha, nn) / Cdr

        if (DD > Drot) then
          DD = DD + (Drot-DD)*Macauley(Mb-Mcur) &
        / (Macauley(Mb-Mcur)+0.01d0)
        end if

        if ((Pmin <= pp) .and. (pp <= 2.0d0*Pmin)) then
          DD = min(DD, -3.5d0*Ado*Macauley(Mb-Md) &
        *(2.0d0*Pmin-pp)/Pmin)
        end if

        ! FIX: C_pmin2 applied in dilation branch (was missing)
        !DD = DD * C_pmin2

      else

        !--- CONTRACTION
        K_p = max(0.0d0, K_p)

        if (ksi <= 0.5d0) then
          hp = hpo*exp(-0.7d0+7.0d0*Macauley(0.5d0-ksi)**2)
        else
          hp = hpo*exp(-0.7d0)
        end if

        Crot2 = 1.0d0 - Czpk2
        Cdz   = max((1.0d0-Crot2*sqrt(2.0d0)*zpeak/z_max) &
        *(z_max/(z_max+Crot2*zcum)), &
        1.0d0/(1.0d0+z_max/2.0d0))

        Adc = Ado*(1.0d0+Macauley(DoubleDot2_2_Contr(fabric,nn))) &
        / (hp*Cdz)
        Cin = 2.0d0*Macauley(DoubleDot2_2_Contr(fabric,nn)) &
        / (sqrt(2.0d0)*z_max)

        ! FIX: threshold 0.16 (was 0.20 in original Fortran)
        DD = min(Adc*(DoubleDot2_2_Contr(alpha_mAlpha_in,nn)+Cin)**2, &
        1.5d0*Ado) &
        * DoubleDot2_2_Contr(alphaD_alpha,nn) &
        / (DoubleDot2_2_Contr(alphaD_alpha,nn) + 0.16d0)!0.16d0)
        
        
        !--- C_pmin2 (applied in both branches — matches C++)
      
        if (pp < Pmin*2.0d0) then
            C_pmin2 = 0.0d0
        else if (pp >= Pmin*18.0d0) then
            C_pmin2 = 1.0d0
        else
            C_pmin2 = (pp - 2.0d0*Pmin) / (16.0d0*Pmin)
        end if
        
        
        

        DD = DD * C_pmin2

      end if

      RR = nn + one3*DD*I1

      end subroutine GetStateDependent


!======================================================================
! IntersectionFactor  (Pegasus method)
!======================================================================
      function IntersectionFactor(CurStress, CurStrain, NextStrain, &
        CurAlpha, Ce, a0_in, a1_in, TolF, mm) result(ElasticRatio)

      implicit none

      real(8), intent(in), dimension(3)   :: CurStress, CurStrain
      real(8), intent(in), dimension(3)   :: NextStrain, CurAlpha
      real(8), intent(in), dimension(3,3) :: Ce
      real(8), intent(in) :: a0_in, a1_in, TolF, mm

      real(8) :: ElasticRatio

      real(8), dimension(3) :: StrainInc, dSigma, tmp
      real(8) :: a0, a1, f, f0, f1, small
      integer(4) :: ii, maxIter

      maxIter = 800;  small = 1.0d-10
      a0 = a0_in;    a1 = a1_in

      StrainInc = NextStrain - CurStrain

      dSigma = DoubleDot4_2(Ce, a0*StrainInc)
      tmp    = CurStress + dSigma
      f0     = GetFYieldFunction(tmp, CurAlpha, mm)

      dSigma = DoubleDot4_2(Ce, a1*StrainInc)
      tmp    = CurStress + dSigma
      f1     = GetFYieldFunction(tmp, CurAlpha, mm)

      do ii = 1, maxIter
        ElasticRatio = a1 - f1*(a1-a0)/(f1-f0)
        dSigma = DoubleDot4_2(Ce, ElasticRatio*StrainInc)
        tmp    = CurStress + dSigma
        f      = GetFYieldFunction(tmp, CurAlpha, mm)
        if (abs(f) < TolF) exit
        if (f*f0 < 0.0d0) then
          a1 = ElasticRatio;  f1 = f
        else
          f1 = f1*f0/(f0+f)
          a0 = ElasticRatio;  f0 = f
        end if
        if (ii == maxIter) then
          ElasticRatio = 0.0d0;  exit
        end if
      end do

      if (ElasticRatio > 1.0d0-small) ElasticRatio = 1.0d0
      if (ElasticRatio < small)        ElasticRatio = 0.0d0

      end function IntersectionFactor


!======================================================================
! IntersectionFactor_Unloading
!======================================================================
      function IntersectionFactor_Unloading(CurStress, CurStrain, &
        NextStrain, CurAlpha, Ce, TolF, mm) result(ElasticRatio)

      implicit none

      real(8), intent(in), dimension(3)   :: CurStress, CurStrain
      real(8), intent(in), dimension(3)   :: NextStrain, CurAlpha
      real(8), intent(in), dimension(3,3) :: Ce
      real(8), intent(in) :: TolF, mm

      real(8) :: ElasticRatio

      real(8), dimension(3) :: StrainInc, dSigma, tmp
      real(8) :: a0, a1, da, f, f0, f1, fs
      integer(4) :: ii, kk, nSub, maxIter
      logical :: flag

      nSub = 200;  maxIter = 500!500,800
      a0 = 0.0d0;  a1 = 1.0d0;  flag = .false.

      StrainInc = NextStrain - CurStrain
      f0        = GetFYieldFunction(CurStress, CurAlpha, mm)
      fs        = f0
      dSigma    = DoubleDot4_2(Ce, StrainInc)

      do ii = 1, maxIter
        da = (a1-a0) / real(nSub, 8)
        do kk = 1, nSub-1
          ElasticRatio = a0 + da*real(kk, 8)
          tmp = CurStress + ElasticRatio*dSigma
          f   = GetFYieldFunction(tmp, CurAlpha, mm)
          if (f > TolF) then
            a1 = ElasticRatio
            if (f0 < -TolF) then
              f1 = f;  flag = .true.;  exit
            else
              a0 = 0.0d0;  f0 = fs;  exit
            end if
          else
            a0 = ElasticRatio;  f0 = f
          end if
        end do
        if (flag) exit
        if (ii == maxIter) then
          ElasticRatio = 0.0d0;  return
        end if
      end do

      ElasticRatio = IntersectionFactor(CurStress, CurStrain, &
        NextStrain, CurAlpha, Ce, a0, a1, TolF, mm)

      end function IntersectionFactor_Unloading


!======================================================================
! UTILITY FUNCTIONS
!======================================================================

      function Macauley(xx) result(res)
        implicit none
        real(8), intent(in) :: xx;  real(8) :: res
        res = max(0.0d0, xx)
      end function Macauley

      function MacauleyIndex(xx) result(res)
        implicit none
        real(8), intent(in) :: xx;  real(8) :: res
        res = merge(1.0d0, 0.0d0, xx > 0.0d0)
      end function MacauleyIndex

      function GetTrace(vv) result(res)
        implicit none
        real(8), intent(in), dimension(3) :: vv
        real(8) :: res
        res = vv(1) + vv(2)
      end function GetTrace

      function GetDevPart(vv) result(dev)
        implicit none
        real(8), intent(in), dimension(3) :: vv
        real(8), dimension(3) :: dev
        real(8) :: p
        p      = 0.5d0*(vv(1)+vv(2))
        dev    = vv
        dev(1) = dev(1) - p
        dev(2) = dev(2) - p
      end function GetDevPart

      function DoubleDot2_2_Contr(v1, v2) result(res)
      ! stress:stress double-dot — shear counts twice
        implicit none
        real(8), intent(in), dimension(3) :: v1, v2
        real(8) :: res
        res = v1(1)*v2(1) + v1(2)*v2(2) + 2.0d0*v1(3)*v2(3)
      end function DoubleDot2_2_Contr

      function DoubleDot2_2_Mixed(v1, v2) result(res)
      ! mixed (strain:stress) double-dot
        implicit none
        real(8), intent(in), dimension(3) :: v1, v2
        real(8) :: res
        res = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
      end function DoubleDot2_2_Mixed

      function GetNorm_Contr(vv) result(res)
        implicit none
        real(8), intent(in), dimension(3) :: vv
        real(8) :: res
        res = sqrt(DoubleDot2_2_Contr(vv, vv))
      end function GetNorm_Contr

      function ToContraviant(v1) result(res)
      ! covariant -> contravariant (halve shear)
        implicit none
        real(8), intent(in), dimension(3) :: v1
        real(8), dimension(3) :: res
        res    = v1
        res(3) = res(3) * 0.5d0
      end function ToContraviant

      function ToCovariant(v1) result(res)
      ! contravariant -> covariant (double shear)
        implicit none
        real(8), intent(in), dimension(3) :: v1
        real(8), dimension(3) :: res
        res    = v1
        res(3) = res(3) * 2.0d0
      end function ToCovariant

      function DoubleDot4_2(m1, v1) result(res)
        implicit none
        real(8), intent(in), dimension(3,3) :: m1
        real(8), intent(in), dimension(3)   :: v1
        real(8), dimension(3) :: res
        integer(4) :: ii, jj
        do ii = 1, 3
          res(ii) = 0.0d0
          do jj = 1, 3
            res(ii) = res(ii) + m1(ii,jj)*v1(jj)
          end do
        end do
      end function DoubleDot4_2

      function GetStiffness(KK_in, GG_in) result(CC)
        implicit none
        real(8), intent(in) :: KK_in, GG_in
        real(8), dimension(3,3) :: CC
        real(8) :: aa, bb
        CC = 0.0d0
        aa = KK_in + (4.0d0/3.0d0)*GG_in
        bb = KK_in - (2.0d0/3.0d0)*GG_in
        CC(1,1) = aa;  CC(2,2) = aa;  CC(3,3) = GG_in
        CC(1,2) = bb;  CC(2,1) = bb
      end function GetStiffness

      function GetKsi(Dr, pp, RR_Bolton, QQ_Bolton, Pmin, P_atm) &
        result(ksi)
        implicit none
        real(8), intent(in) :: Dr, pp, RR_Bolton, QQ_Bolton, Pmin, P_atm
        real(8) :: ksi
        real(8) :: pn
        pn  = max(pp, Pmin)
        ksi = RR_Bolton/(QQ_Bolton-log(100.0d0*pn/P_atm)) - Dr
      end function GetKsi

      function GetNormalToYield(stress, alpha) result(nn)
        implicit none
        real(8), intent(in), dimension(3) :: stress, alpha
        real(8), dimension(3) :: nn
        real(8) :: pp, normN
        real(8), parameter :: root12 = 0.7071067811865475d0
        real(8), parameter :: small  = 1.0d-10
        pp = 0.5d0*(stress(1)+stress(2))
        if (abs(pp) < small) then
          nn    = 0.0d0;  nn(3) = root12
        else
          nn    = GetDevPart(stress) - pp*alpha
          normN = GetNorm_Contr(nn)
          if (normN < small) normN = 1.0d0
          nn = nn / normN
        end if
      end function GetNormalToYield

      function GetFYieldFunction(nStress, nAlpha, m_m) result(ff)
        implicit none
        real(8), intent(in), dimension(3) :: nStress, nAlpha
        real(8), intent(in) :: m_m
        real(8) :: ff, pp
        real(8), dimension(3) :: ss
        real(8), parameter :: root12 = 0.7071067811865475d0
        pp = 0.5d0*(nStress(1)+nStress(2))
        ss = GetDevPart(nStress) - pp*nAlpha
        ff = GetNorm_Contr(ss) - root12*m_m*pp
      end function GetFYieldFunction
      
      
      
      
      
      
      
      
! ================================================================
! CHANGE 3: Add this function anywhere inside the module
!           (e.g. just before the end module statement,
!            alongside GetStiffness and other utility functions)
! ================================================================

      function GetElastoPlasticTangent(NextStress, aCe, RR, nn, K_p, &
        Pmin, mm) result(Cep)

      implicit none

      real(8), intent(in), dimension(3)   :: NextStress, RR, nn
      real(8), intent(in), dimension(3,3) :: aCe
      real(8), intent(in)                 :: K_p, Pmin, mm

      real(8), dimension(3,3) :: Cep

      real(8), dimension(3) :: temp1, temp2, rrr_loc
      real(8) :: temp3, pp, small
      integer(4) :: ii, jj

      small = 1.0d-10

      pp      = max(0.5d0*(NextStress(1)+NextStress(2)), Pmin)
      rrr_loc = GetDevPart(NextStress) / pp

      ! temp1 = Ce * R
      do ii = 1, 3
        temp1(ii) = aCe(ii,1)*RR(1) &
        + aCe(ii,2)*RR(2) &
        + aCe(ii,3)*RR(3)
      end do

      ! dF/dSigma = n - 0.5*(n:r)*I1
      temp2(1) = nn(1) - 0.5d0*DoubleDot2_2_Contr(nn, rrr_loc)
      temp2(2) = nn(2) - 0.5d0*DoubleDot2_2_Contr(nn, rrr_loc)
      temp2(3) = nn(3)

      ! temp2 = (dF/dSigma) * Ce  (Ce is symmetric)
      rrr_loc(1) = aCe(1,1)*temp2(1) &
        + aCe(2,1)*temp2(2) &
        + aCe(3,1)*temp2(3)
      rrr_loc(2) = aCe(1,2)*temp2(1) &
        + aCe(2,2)*temp2(2) &
        + aCe(3,2)*temp2(3)
      rrr_loc(3) = aCe(1,3)*temp2(1) &
        + aCe(2,3)*temp2(2) &
        + aCe(3,3)*temp2(3)
      temp2 = rrr_loc

      ! denominator = temp2:R + K_p
      temp3 = DoubleDot2_2_Contr(temp2, RR) + K_p

      ! Cep = Ce - outer(temp1, temp2) / temp3
      if (temp3 > small) then
        do ii = 1, 3
          do jj = 1, 3
            Cep(ii,jj) = aCe(ii,jj) &
        - temp1(ii)*temp2(jj) / temp3
          end do
        end do
      else
        Cep = aCe
      end if

      end function GetElastoPlasticTangent
     



!======================================================================
      end module PM4SandMaterialModule

    
    
    
    
    