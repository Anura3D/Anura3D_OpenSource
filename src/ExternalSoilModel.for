	!*****************************************************************************
    !                                       ____  _____  
    !           /\                         |___ \|  __ \ 
    !          /  \   _ __  _   _ _ __ __ _  __) | |  | |
    !         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
    !        / ____ \| | | | |_| | | | (_| |___) | |__| |
    !       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/ 
    !
    !
	!	Anura3D - Numerical modelling and simulation of large deformations 
    !   and soil–water–structure interaction using the material point method (MPM)
    !
    !	Copyright (C) 2024  Members of the Anura3D MPM Research Community 
    !   (See Contributors file "Contributors.txt")
    !
    !	This program is free software: you can redistribute it and/or modify
    !	it under the terms of the GNU Lesser General Public License as published by
    !	the Free Software Foundation, either version 3 of the License, or
    !	(at your option) any later version.
    !
    !	This program is distributed in the hope that it will be useful,
    !	but WITHOUT ANY WARRANTY; without even the implied warranty of
    !	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !	GNU Lesser General Public License for more details.
    !
    !	You should have received a copy of the GNU Lesser General Public License
    !	along with this program.  If not, see <https://www.gnu.org/licenses/>.
	!
    !*****************************************************************************


module ModExternalSoilModel
!**********************************************************************
!
!  Function: Contains the routines related to calling user-defined soil models in external DLLs.
!     $Revision: 9064 $
!     $Date: 2023-08-07 13:55 -0400 (Moore, 08 Aug 2023) $
!
!**********************************************************************
      
use ModMPMData
use ModGlobalConstants
use ModReadCalculationData
use ModReadMaterialData
use ModMPMInit
use user32
use kernel32
use ModMeshInfo
use ModLinearElasticity
use ModMohrCoulomb
use ModBingham

contains


subroutine StressSolid(IDpt, IDel, BMatrix,IEntityID)
!**********************************************************************
!
!    Function:  calculate stresses at material point using external soil models
!
!*********************************************************************        
    
implicit none
        
    integer(INTEGER_TYPE), intent(in) :: IDpt ! global integration/material point number
    integer(INTEGER_TYPE), intent(in) :: IDel ! global element number
    ! B-matrix at the considered integration point (here only used if ApplyObjectiveStress=TRUE)
    real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: BMatrix    
    integer(INTEGER_TYPE), intent(in) :: IEntityID ! entity ID (here only used if ApplyObjectiveStress=TRUE)
    
    ! local variables
    character(len=80) :: cmname
    integer(INTEGER_TYPE) :: I ! counter
    integer(INTEGER_TYPE) :: IDset ! ID of material parameter set
    integer(INTEGER_TYPE) :: ntens ! Dimension of stress vector to pass to ESM 
    integer(INTEGER_TYPE), parameter :: nAddVar = 12
    real(REAL_TYPE), dimension(NPROPERTIES) :: props ! array of material properties
    real(REAL_TYPE), dimension(nAddVar) :: AdditionalVar
    real(REAL_TYPE), dimension(MatParams(MaterialIDArray(IDpt))%UMATDimension) :: Stress, StrainIncr ! stress and strain increment in integration/material point
    real(REAL_TYPE), dimension(NTENSOR) :: Sig0, StressIncr, StressPrinc, TempStrainIncr, TempStrainIncrPrevious
    real(REAL_TYPE), dimension(NSTATEVAR) :: StateVar ! state parameters in integration/material
    real(REAL_TYPE) :: Eunloading, PlasticMultiplier
    character(len=64) :: NameModel ! name of the constitutive model
    logical :: IsUndrEffectiveStress
    real(REAL_TYPE) :: DSigWP ! Change of water pressure at integration point 
    real(REAL_TYPE) :: DSigGP ! Change of gas pressure at integration point 
    real(REAL_TYPE) :: Bulk ! Bulk modulus
    real(REAL_TYPE) :: DEpsVol ! Incremental volumetric strain (water)
    procedure(DUMMYESM), pointer :: ESM
    
    ! Creating a pointer for the External Soil models
    pointer (p, ESM)             

    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
    NameModel = MatParams(IDset)%MaterialModel ! name of constitutive model as specified in GOM-file
    ntens = MatParams(IDset)%UMATDimension  ! 2D or 3D formulation of the External soil model   
          
    ! get strain increments in integration/material point
    TempStrainIncr = GetEpsStep(Particles(IDpt)) ! incremental strain vector assigned to point
    
        
    StrainIncr = 0.0

    do I=1, NTENSOR
    StrainIncr(I) = StrainIncr(I) + TempStrainIncr(I)
    enddo 
        
    DEpsVol = StrainIncr(1) + StrainIncr(2) + StrainIncr(3) ! volumetric strain, valid for 2D and 3D
          
    IsUndrEffectiveStress = &
    !code version 2016 and previous
    ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(IDSet)%MaterialType)=='2-phase')) .or. &
    !code version 2017.1 and following
    (trim(MatParams(IDSet)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
          
    ! initalise water pressure (only needed for undrained analyses)
    DSigWP = 0.0d0
    DSigGP = 0.0d0
    ! for effective stress analysis
    if (IsUndrEffectiveStress) then
        if (Particles(IDpt)%Porosity > 0.0) then
        Bulk = Particles(IDpt)%BulkWater / Particles(IDpt)%Porosity ! kN/m2
        DSigWP = Bulk * DEpsVol
        else
        DSigWP = 0.0
        end if
        call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP)
    end if ! effective stress analysis
          
    ! get stresses in integration/material point      
    do I = 1, NTENSOR
    Sig0(I) = SigmaEff0Array(IDpt, I) ! get initial stress of current step assigned to point 
    end do
    Stress=0.0
    do I=1, NTENSOR
        Stress(I) = Stress(I) + Sig0(I) !To use 3D UMAT also for 2D formulations
    enddo 
          
    ! initialise state variables (only for very first time and load step)

    !!! Commented out on 08/07/23 - reason: Update from the 2023 release

    ! if ((CalParams%IStep == 1).and.(CalParams%TimeStep == 1)) then
    ! StateVar = MatParams(IDset)%ESM_Statvar_in ! array of NSTATEVAR initial value of the State Variables for the external soil model in UMAT/VUMAT format
    ! else 
    StateVar = ESMstatevArray(IDpt,:)
    ! end if 
          
    
    
    !call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP) !Note that the subroutine checks Cavitation Threshold & Gas Pressure
       
    ! Undrained effective stress pore pressure
    if (IsUndrEffectiveStress) then
    Particles(IDPt)%WaterPressure = Particles(IDPt)%WaterPressure + DSigWP
    end if 
    
    !get values of variables of interest for UMAT model
    AdditionalVar(1) = Particles(IDPt)%Porosity
    AdditionalVar(2) = Particles(IDPt)%WaterPressure
    AdditionalVar(3) = Particles(IDPt)%WaterPressure0 
    AdditionalVar(4) = Particles(IDPt)%GasPressure
    AdditionalVar(5) = Particles(IDPt)%GasPressure0
    AdditionalVar(6) = Particles(IDPt)%DegreeSaturation
    AdditionalVar(7) = CalParams%TotalRealTime
    AdditionalVar(8) = CalParams%OverallRealTime
    AdditionalVar(9) = CalParams%TimeIncrement
    AdditionalVar(10) = CalParams%IStep
    AdditionalVar(11) = CalParams%TimeStep
    AdditionalVar(12) = Particles(IDpt)%BulkWater
          
    ! get name of DLL
    cmname = MatParams(IDSet)%SoilModelDLL
    ! get material properties  
    props = MatParams(IDSet)%ESM_Solid ! This sets props to be the material parameters for an external soil model
         
    if (trim(NameModel)//char(0) == trim('linear_elasticity')//char(0)) then
    props(1) = Particles(IDpt)%ShearModulus ! shear modulus, G
    cmname = UMAT_LINEAR_ELASTICITY
    elseif (trim(NameModel)//char(0) == trim(ESM_MOHR_COULOMB_STANDARD)//char(0)) then
    props(1) = Particles(IDpt)%ShearModulus ! shear modulus, G
    props(2) = MatParams(IDSet)%PoissonRatio 
    props(3) = SIN(MatParams(IDSet)%FrictionAngle*(Pi/180.0)) 
    props(4) = Particles(IDpt)%CohesionCosPhi 
    props(5) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0))
    props(6) = MatParams(IDSet)%TensileStrength
    cmname = UMAT_MOHR_COULOMB_STANDARD
    endif    
    
    ! initialise UMAT
    
    ! Create the external soil model pointer
    ! This might cause problems if parallelization is attempted for the constitutive models. 
    ! This pointer might be referencing the same subroutine. I'm not sure how fortran treats function pointers
    ESM => MatParams(IDSet)%ESM_POINTER
    call ESM(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)

    ! save unloading stiffness in Particles array  
    Particles(IDpt)%ESM_UnloadingStiffness = Eunloading
                 
    if (IsUndrEffectiveStress) then
        Particles(IDpt)%BulkWater = AdditionalVar(12)
    end if
    
    call SetIPL(IDpt, IDel, int(PlasticMultiplier))


    
    ! to use objective stress definition
    if (CalParams%ApplyObjectiveStress) then ! Consider large deformation terms
    call Hill(IdEl, ELEMENTNODES, IncrementalDisplacementSoil(1:Counters%N, IEntityID),  &
                     ReducedDof, ElementConnectivities, BMatrix, Sig0(1:NTENSOR), Stress(1:NTENSOR), DEpsVol)
    end if ! objective stress            
            
    ! write new stresses to global array
    do I=1, NTENSOR
        StressIncr(I) = Stress(I) - Sig0(I)
    end do             
                               
    ! save updated state variables and in Particles array
    ESMstatevArray(IDpt,:) = StateVar
          
    call CalculatePrincipalStresses(IDpt, Stress(1:NTENSOR), StressPrinc)
    call AssignStressStrainToGlobalArrayESM(IDpt, NTENSOR, StressIncr, StressPrinc, StrainIncr)

    ! write plasticity state to global array
    !  call SetIPL(IDpt, IDel, int(StateVar(50)))
    if (CalParams%ApplyBulkViscosityDamping) then
    RateVolStrain(IDEl) = DEpsVol / CalParams%TimeIncrement
    call CalculateViscousDamping_interface(IDpt, IDEl)
    end if  

end subroutine StressSolid

        subroutine CalculateViscousDamping_interface(ParticleID, IEl)
        !**********************************************************************
        !
        !> Computes a pressure term introducing bulk viscosity damping to the equation of motion.
        !>
        !! \param[in] ParticleID ID of considered material point.
        !! \param[in] IEl ID of element of the considered material point.
        !! \param[in] DilationalWaveSpeed Current wave speed computed for the considered material point.
        !
        !*********************************************************************
        
        implicit none
          integer(INTEGER_TYPE), intent(in) :: ParticleID, IEl

          real(REAL_TYPE) :: ViscousDampingPressure = 0.0
          real(REAL_TYPE) :: Density = 0.0
          real(REAL_TYPE) :: ElementLMinLocal = 0.0
          real(REAL_TYPE) :: RateVolStrainLocal = 0.0
          real(REAL_TYPE) :: MaterialIndex = 0.0
          real(REAL_TYPE) :: DilationalWaveSpeed = 0.0
          logical :: IsUndrEffectiveStress

          if (.not.CalParams%ApplyBulkViscosityDamping) return

          MaterialIndex = MaterialIDArray(ParticleID) !l;aksdjf;lkajsdf
          
           IsUndrEffectiveStress = &
              !code version 2016 and previous
              ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
              !code version 2017.1 and following
              (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
           
          !if (CalParams%ApplyEffectiveStressAnalysis
           if (IsUndrEffectiveStress &
          .or.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3))) then
            Density = MatParams(MaterialIndex)%DensityMixture / 1000.0
          else
            Density = (1 - MatParams(MaterialIndex)%InitialPorosity) * MatParams(MaterialIndex)%DensitySolid / 1000.0
          end if

          ElementLMinLocal = ElementLMin(IEl)
          RateVolStrainLocal = RateVolStrain(IEl)

          
          call GetWaveSpeed(ParticleID, DilationalWaveSpeed)

          ViscousDampingPressure = CalParams%BulkViscosityDamping1 *  &
            Density * DilationalWaveSpeed * ElementLMinLocal * RateVolStrainLocal

          if ((RateVolStrainLocal < 0.0).and.(CalParams%BulkViscosityDamping2 > 0.0)) then
            ViscousDampingPressure = ViscousDampingPressure + &
              Density * (CalParams%BulkViscosityDamping2 * ElementLMinLocal * RateVolStrainLocal)**2
          end if

          Particles(ParticleID)%DBulkViscousPressure = ViscousDampingPressure

        end subroutine CalculateViscousDamping_interface
        
        !*************************************************************************************************************************************
        !*******************************Everything below this line has been added when compared to the original (Jonathan Moore)**************
        !*************************************************************************************************************************************
    
!
    Subroutine ESM_VPSS_MC(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER, DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)
        
              !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"ESM" :: ESM
              !---Updating these variables to follow the updated fortran format

                    implicit none 
                    !!!!! These intents may need to be changed to intent(in) !!!!!
                    integer, intent(in) :: nstatev, ntens, nprops
                    integer :: naddvar
                    integer, intent(in) :: npt
                    integer, intent(in) :: noel
                    integer, intent(in) :: idset
                    double precision, dimension(ntens), intent(inout) :: stress
                    double precision, intent(inout) :: eunloading, plasticmultiplier
                    double precision, dimension(ntens), intent(inout) :: dstran
                    double precision, dimension(nstatev), intent(inout) :: statev
                    double precision, dimension(naddvar) :: additionalvar 
                    character(len=80) :: cmname !This is a parameter do not set intent
                    double precision, dimension(nprops), intent(inout) :: props
                    integer:: numberofphases    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
              !---Local variables required in standard UMAT
                    integer :: IStep, TimeStep
                    double precision :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
                    double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
                    double precision, dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
                    double precision, dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
                    double precision                            :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
                    double precision, dimension(:), allocatable :: stran
                    double precision, dimension(:), allocatable :: time
                    double precision, dimension(:), allocatable :: predef
                    double precision, dimension(:), allocatable :: dpred    
                    double precision, dimension(:), allocatable :: coords
                    double precision, dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
                    double precision, dimension(:,:), allocatable :: drot
                    double precision, dimension(:,:), allocatable :: dfgrd0
                    double precision, dimension(:,:), allocatable :: dfgrd1
                    double precision :: pnewdt, dtime, temp, dtemp, celent
                    double precision :: Value ! auxiliary variable holding any real valued number
                    double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation  
                    integer          :: ndi, nshr, layer, kspt, kstep, kinc     
              
              !---Local variables defned by the user
              ! e.g. integer :: var_local	  
              !---User can define here additional variables     
              
                    allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1), &
                          coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
              
              !Initialization
                    Eunloading = 0.0
                    PlasticMultiplier = 0.0
                
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
					
              !Call the UMAT
                    call UMAT_VPSS_MC(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
                            dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, pnewdt, celent, dfgrd0, &      
                            dfgrd1, noel, npt, layer, kspt, kstep, kinc)
              
              
              !---Definition of Eunloading -> required to define the max time step
                    Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
              !---Always define this value to run the simulation
              
              ! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
            return
        end subroutine ESM_VPSS_MC
  
  
        Subroutine UMAT_VPSS_MC(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, &
				                DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,&
				                CMNAME, NDI, NSHR, NTENS, &	NSTATEV, PROPS, NPROPS, COORDS, &
				                PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, &
				                LAYER, KSPT, KSTEP, KINC)

             !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
                     implicit none
		            !Defining inputs

            ! Arguments:
            !          I/O  Type
            !  PROPS    I   R()  : List with model parameters
            !  DSTRAN   I   R()  : Strain increment
            !  DTIME    I   R()  : Time increment
            !  DDSDDE   O   R(,) : Material stiffness matrix
            !  STRESS  I/O  R()  : stresses
            !  STATEV  I/O  R()  : state variables
            !
         
            integer, intent(in) :: NSTATEV, NPROPS, NPT
            integer :: NTENS
            integer, intent(in) :: NOEL 
			real (Real_Type), dimension(NTENS), intent(inout) :: STRESS
            real(Real_Type), dimension(NSTATEV), intent(inout) :: STATEV
            real(Real_Type), dimension(NTENS, NTENS), intent(inout) :: DDSDDE
			real(Real_Type), intent(in) :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
            real(Real_Type), intent(in) :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
			real(Real_Type), dimension(NTENS), intent(in) :: DDSDDT
            real(Real_Type), dimension(NTENS), intent(in) :: DRPLDE
			real(Real_Type)                               :: DRPLDT
            real(Real_Type), dimension(NTENS), intent(in) :: STRAN
            real(Real_Type), dimension(NTENS), intent(in) :: DSTRAN
            real(Real_Type), dimension(2), intent(in) :: TIME
            real(Real_Type), dimension(1), intent(in) :: PREDEF
            real(Real_Type), dimension(1), intent(in) :: DPRED
            real(Real_Type), dimension(NPROPS), intent(in) :: PROPS
            real(Real_Type), dimension(3), intent(in) :: COORDS
            real(Real_Type), dimension(3,3), intent(in) :: DFGRD0
            real(Real_Type), dimension(3,3), intent(in) :: DFGRD1
            REAL(Real_Type), intent(in) :: PNEWDT,  TEMP, DTEMP, CELENT
			double precision, intent(in) :: DTIME
            character(len = 80), intent(in):: CMNAME
			integer, intent(in) :: NDI, NSHR, LAYER, KSPT, KSTEP, KINC  
      
            ! Local variables:
            !
            !  DE        : Linear Elastic constitutive matrix
            !  dSig	     : Stress increment vector
            !  Sig	     : Stress vector
            !  dEpsE     : Elastic strain increment vector
            !  dEpsP     : Plastic strain increment vector
            !  dEps      : Total strain increment vector
            !  EpsE      : Elastic strain vector
            !  EpsP      : Plastic strain vector
            !  Eps       : Total strain vector
            !  EpsRate	 : Total strain rate tensor
            !
            integer     			:: N_S, N_i, N_m_s, N_m_i
            real(Real_Type), dimension(6,6) :: DE
            real(Real_Type), dimension(6)   :: Sig
            real(Real_Type), dimension(6)   :: dEpsP
            real(Real_Type), dimension(6)   :: EpsP, ERate, dummy_Erate, dummy_Erate_2
            real(Real_Type), dimension(6) :: dEpsE, EpsE, dEps, Eps
            logical :: switch_smooth, switch_original, switch_yield
            real(Real_Type) :: G_0, enu, eM_tc, eN, D_min, eh, alpha_G, alpha_K, alpha_D, D_part, G_s, RefERate !SSMC props local variables, (props)
            real(Real_Type) :: G, bk, eta_y, DP, eI_coeff, Sum_rate  !SSMC state variables (statv)
            double precision :: Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, Error_Yield_max
            double precision :: F1, F2, bK_0
            integer :: i
            
            
        !
		! Viscoplastic NAMC model with smoothed outer surface in pi plane
		!
		! Contents of PROPS(10) NAMC with HSR
		!  1 : G_0				shear modulus
		!  2 : enu				Poisson's ratio
		!  3 : eM_tc            Critical stress ratio for triaxial compression
        !  4 : eN               Nova's vol coupling coefficient
		!  5 : D_min			Minimum dilation
		!  6 : eh				hardening parameter
		!  7 : alpha_G			Shear modulus viscosity factor
		!  8 : alpha_K			Bulk modulus viscosity factor
		!  9 : alpha_D			dilation viscosity
		! 10 : D_part			Particle diameter
		! 11 : G_s				Specific gravity      
		! 12 : RefERate			Reference strain rate
		! 13 : Switch_smooth	Boolean switch for activating strain rate smoothing
		! 14 : N_S				Degree of smoothening
        ! 15 : switch_original	Changes from Wang's to Olzak&Perzyna consistency 
        ! 16 : N_m_s            # of cycles to maintain constant strain rate

		G_0			= PROPS(1)         ! shear modulus
        enu			= PROPS(2)         ! Poisson's ratio
        eM_tc		= PROPS(3)         ! Critical stress ratio
        eN          = PROPS(4)         ! Nova's vol coupling coefficient
        D_min       = PROPS(5)         ! Minimum dilation
        eh			= PROPS(6)         ! hardening parameter
		alpha_G	    = PROPS(7)	       ! Shear modulus viscosity factor
		alpha_K	    = PROPS(8)	       ! Bulk modulus viscosity factor
		if ((alpha_K==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
			alpha_K= 2.5*alpha_G
		endif
		alpha_D		= PROPS(9)    ! dilation angle viscosity
		if ((alpha_D==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
			alpha_D= alpha_G
        endif
        D_part		= PROPS(10)			! Associated particle diameter, value in mm
        G_s  		= PROPS(11)			! Specific gravity
		RefERate	= PROPS(12)		    ! reference strain rate
		if (RefERate==0.0d0) then
			RefERate=2.5e-5
		endif
		call dbltobool(PROPS(13), switch_smooth)  ! switch for activating strain rate smoothening
		N_S=PROPS(14)							  ! Degree of smoothening
        call dbltobool(PROPS(15), switch_original)! (1 for Wang, 0 for Olzak&Perzyna)
        N_m_s       = PROPS(16)                   ! Number of cycles to maintain strain rate
        
        G			    = STATEV(1)			         ! shear modulus
		bK			    = STATEV(2)			         ! bulk modulus
		eta_y		    = STATEV(3)                  ! friction ratio
        Dp		        = STATEV(4)                  ! Dilation
        eI_coeff		= STATEV(5)                  ! Inertial coefficient

        call dbltobool(STATEV(6),switch_yield)       ! Point is yielding
		do i=1,6
			EpsP(i)     =STATEV(6+i)				 !Plastic strain               
		end do
		N_i             =STATEV(13)			         ! Current number of strain rate sums
		SUM_rate        =STATEV(14)			         ! Current sum of strain rates
        N_m_i           =STATEV(15)                  ! Current number of cycles strain rate held constant
        do i=1,6
            dummy_Erate(i) = STATEV(15+i)
        end do
        do i=1,6
            dummy_Erate_2(i) = STATEV(21+i)
        end do
        
		!_____Error control state parameters__________________________________________________
		Error_yield_1=0.0d0                                                                   !
        Error_yield_2=0.0d0                                                                   !
        Error_Euler_max=0.0d0															      !
		Error_Yield_last=0.0d0															      !
		Error_Yield_max=0.0d0																  !
		!							                                                          !
		!_____________________________________________________________________________________!		


		if (DTIME==0.0d0) then
			ERate= 0.0d0    ! Current strain rate
		else
			ERate= (1/DTIME)*DSTRAN ! Current strain rate
        end if
        
		bK_0= 2*G_0*(1+ENU)/(3*(1-2*ENU))

		!***********************************************************************************

		!Call the refined modified Euler algorithm
		call NAMC_HSR(NOEL, G_0, enu, eM_tc, eN, D_min, eh, alpha_G, alpha_K, alpha_D, D_part, G_s,&
                     switch_smooth, RefERate, N_S, switch_original,&
					 G, bK, eta_y, Dp, EpsP, eI_coeff, switch_yield, N_i, SUM_rate,&
					 DSTRAN, STRESS, Sig, Erate, DTIME,&
					 Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, &
					 Error_Yield_max, N_m_s, N_m_i, dummy_Erate, dummy_Erate_2)
		!************************************************************************************

		!************************************************************************************
		!Stress and state variables updating
		Do i=1,NTENS
          STRESS(i) = Sig(i)
        End Do

        STATEV(1) = G
        STATEV(2) = bK
        STATEV(3) = eta_y
		STATEV(4) = Dp
        STATEV(5) =	eI_coeff	
		STATEV(6) = logic2dbl(switch_yield)
		do i=1,6
			STATEV(6+i) = EpsP(i)
		end do			
		STATEV(13)=N_i
		STATEV(14)=SUM_rate
        STATEV(15)=N_m_i
        do i=1,6
            STATEV(15+i)=dummy_Erate(i)
        end do
        do i=1,6
            STATEV(21+i)=dummy_Erate_2(i)
        end do
		!_____Error control state parameters__________________________________________________
		!     Comment if not wanted                                                           !
		!_____________________________________________________________________________________!	
		STATEV(28)=Error_yield_1
        STATEV(29)=Error_yield_2
        STATEV(30)=Error_Euler_max
		STATEV(31)=Error_Yield_last
		STATEV(32)=Error_Yield_max
        
		!************************************************************************************

		!************************************************************************************
		!Tangent stiffness matrix to be returned done by elastic stiffness
        F1  = bK+(4*G/3)
        F2  = bK-(2*G/3)
        DDSDDE = 0.0
        DDSDDE(1:3,1:3) = F2
        DDSDDE(1,1) = F1
        DDSDDE(2,2) = F1
        DDSDDE(3,3) = F1
        DDSDDE(4,4) = G
        DDSDDE(5,5) = G
        DDSDDE(6,6) = G
		!*************************************************************************************

		!End of UMAT
		Return
    end subroutine UMAT_VPSS_MC
    
!_____________________________________________________________________________________________    
!##    ##    ###    ##     ##  ######     ##     ##  ######  ########  
!###   ##   ## ##   ###   ### ##    ##    ##     ## ##    ## ##     ## 
!####  ##  ##   ##  #### #### ##          ##     ## ##       ##     ## 
!## ## ## ##     ## ## ### ## ##          #########  ######  ########  
!##  #### ######### ##     ## ##          ##     ##       ## ##   ##   
!##   ### ##     ## ##     ## ##    ##    ##     ## ##    ## ##    ##  
!##    ## ##     ## ##     ##  ######     ##     ##  ######  ##     ##  
    
    subroutine NAMC_HSR(NOEL, G_0, nu, M_tc, N, D_min, h, alpha_G, alpha_K, alpha_D, D_part, G_s,&
                     switch_smooth, RefERate, N_S, switch_original, &
					 G, K, eta_y, Dp, EpsP, I_coeff, switch_yield, N_i, SUM_rate,&
					 dEps, Sig_0, Sig, Erate, DTIME,&
					 Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, &
					 Error_Yield_max, N_m_s, N_m_i, dummy_Erate, dummy_Erate_2)
    
      !*********************************************************************************
      !
      ! Elasto-plastic constitutive model with strain softening and strain rate effects,
	  ! Non associated MOHR-COULOMB 
      ! Explicit MODIFIED EULER INTEGRATION SCHEME with automatic error control.
      ! Final correction of the yield surface drift (END OF STEP CORRECTION).
      !
      !*********************************************************************************
		implicit none
        
		!input variables
		integer, intent(in) :: NOEL !Global ID of Gauss point or particle
        integer, intent(in) :: N_S, N_m_s
        double precision, intent(in)::  G_0, nu, M_tc, N, D_min, alpha_G, alpha_K, &
										alpha_D, D_part, G_s, &
										RefERate, DTIME
        double precision, dimension(6), intent(in)::  dEps
        double precision, dimension(6), intent(inout):: Sig_0
        
        logical, intent(in):: switch_smooth, switch_original           
        
		!output variables
        integer, intent(inout):: N_i, N_m_i
        
        double precision, intent(inout):: h, G, K, eta_y, Dp, I_coeff, &
                                          Error_yield_1, Error_yield_2, Error_Euler_max,&
										  Error_Yield_last, Error_Yield_max, SUM_rate
        double precision, dimension(6), intent(inout):: Sig, Erate, EpsP, dummy_Erate, dummy_Erate_2
        
        logical, intent(inout):: switch_yield
    
		!local variables
        integer:: counter, MAXITER, SubStepping_MaxIter
        double precision:: p, q, theta, M, I_act, I_0, Gu, Ku, eta_yu, Du, dI, &
                           G1, K1, eta_y1, Dp1, G2, K2, eta_y2, Dp2, &
                           p_t, q_t, dI_t, dI_TT, I_TT, dD1, dD2
        double precision:: epsq_rate, epsq_p, eps_v
        double precision:: F0, FT, alpha
        double precision:: FTOL, STOL, DTmin , LTOL, R_TT, qR, RTOL
        double precision:: dummyvar(3), D1, D2
        double precision:: DE(6,6), dSig_el(6), Sig_t(6), dEps_t(6), dEps_TT(6), &
                           Sig1(6), Sig2(6), dEpsp1(6), dEpsp2(6), dEpsp(6), &
                           dSig1(6), dSig2(6), Epspt(6)
        double precision:: T, DT
        
        logical:: ApplyStrainRateUpdates=.false., &!Check if the strain rate path crosses the reference line
                  IsUnloading, & !If true material unloads (not elastic unloading)
                  Failed=.false. !If rel residual is too high failed is true
        !print *, NOEL
		!______________________________________________________________________________
        ! Error tolerances
        FTOL=1.0e-9		!Yield surface tolerance
        STOL=1.0e-3		!Relative error tolerance
        DTmin=1.0e-9	!Minimum pseudo-time increment, originally Dtmin = 1.0e-9
		LTOL=0.01d0		!Tolerance for elastic unloading	  
		MAXITER=20		!Max. number of iterations
        RTOL = STOL * 1.0e-1
        !______________________________________________________________________________
		!______________________________________________________________________________
        ! Initialization of error trackers
        Error_yield_1   =0.0	!Drift after first approximation
        Error_yield_2   =0.0	!Drift after second approximation
		Error_Euler_max =0.0    !Max relative residual
		Error_Yield_max =0.0    !max drift after averaging
		Error_Yield_last=0.0    !max abs drift after drift correction		
        !______________________________________________________________________________
        !Initialization of state variables
        call Get_invariants(Sig_0, p, q, theta)
		call Get_strain_invariants(EpsP, eps_v, epsq_p)!plastic deviatoric strain
        I_0=RefERate
        !call Get_I_coeff(D_part, G_s, -100.0, RefERate, I_0)!Reference inertial coefficient
        call Get_M(M_tc, theta, M)!Get M
        if (G==0.0d0) then
            G=G_0
        endif
        if (K==0.0d0) then
            K=2*G_0*(1+nu)/(3*(1-2*nu))
        endif        
        !if (I_coeff==0.0d0) then
        !    call Get_I_coeff(D_part, G_s, p, RefERate, I_coeff)
        !endif

        if (Dp==0.0d0) then ! Since plastic deviatoric strain epsq_p will be zero in elastic Dp will be be zero 
            call Get_Dp(h, D_min, I_coeff, I_0, epsq_p, alpha_D, ApplyStrainRateUpdates, Dp)
        endif
        
        !print *, eta_y
        if (eta_y==0.0d0) then
            call Get_invariants(dEps, dummyvar(1), dummyvar(2), theta)
            call Get_M(M_tc, theta, M)
            eta_y=M-Dp*(1.0-N)
        endif
        !print *, eta_y

        !_____________________________________________________________________________
        !Evaluate yield function at initial stress-state variable
        
        call YieldFunction(q, p, eta_y, F0)
        
        !_____________________________________________________________________________
        !Compute the current deviatoric strain rate
        
        call TwoNormTensor_strain(Erate,6,dummyvar(1))!Norm of the strain rate
        TrackStrainRate = dummyvar(1) !Track the strain rate
        
		 !Compute a smoothed strain rate if switch_smooth is true
		if (switch_smooth) then
			N_i=N_i+1
            !if (Noel ==1) then
            !    print *, N_i
            !endif
			if (N_i<N_s) then !not enough values
				Sum_rate=sum_rate+dummyvar(1) !accumulate the strain rate
				dummyvar(2)=Sum_rate/n_i !takes the average
			else !enough values
				Sum_rate=Sum_rate*(1.0-1.0/N_s) !approximate sum without one term
				Sum_rate=Sum_rate+dummyvar(1)
				dummyvar(2)=Sum_rate/N_s !averaged strain rate
            endif
            if (dummyvar(1)==0.0d0) then !in case in first step this value is zero
                Erate=0.0d0
            else
                Erate=(dummyvar(2)/dummyvar(1))*Erate !corrected strain rate tensor
            endif
            call TwoNormTensor_strain(Erate,6,dummyvar(1)) ! recalculate the two norm so the updated erate value can be tracked
            TrackSmoothStrainRate = dummyvar(1) !track the smoothed strain rate
        endif
        
        !The purpose of the below code is to fix the value of Erate for a series of steps
        !if (switch_smooth) then
        !    N_i = N_i + 1
        !    N_m_i = N_m_i + 1
        !    if (N_i < N_s) then
        !        SUM_rate = SUM_rate + dummyvar(1) ! Accumulate the strain rate
        !        dummyvar(2) = SUM_rate/N_i ! Takes the average
        !        if (dummyvar(1) ==0.0d0) then ! in case in the first step this value is zero
        !            Erate = 0.0d0
        !        else
        !            Erate = (dummyvar(2)/dummyvar(1)) * Erate ! corrected strain rate tensor
        !        endif
        !    elseif (N_m_i < N_m_s) then !Enough values
        !        SUM_rate = SUM_rate + dummyvar(1) ! Accumulate the strain rate
        !        dummyvar(2) = SUM_rate/N_m_i ! Takes the average
        !        dummy_Erate = (dummyvar(2)/dummyvar(1)) * Erate ! corrected strain rate tensor
        !        
        !        !SUM_rate = SUM_rate * (1-1/N_s)
        !        !SUM_rate = SUM_rate + dummyvar(1)
        !        !dummyvar(2)=SUM_rate/N_s
        !        !dummy_Erate = (dummyvar(2)/dummyvar(1)) * Erate
        !
        !        if (N_m_i == N_s .or. N_m_i == 1) then
        !            dummy_Erate_2 = dummy_Erate 
        !        endif
        !        Erate = dummy_Erate_2
        !    else
        !        SUM_rate = SUM_rate * (1-1/N_m_i)
        !        SUM_rate = SUM_rate + dummyvar(1)
        !        dummyvar(2)=SUM_rate/N_m_i
        !        dummy_Erate = (dummyvar(2)/dummyvar(1)) * Erate
        !        Erate = dummy_Erate
        !        N_m_i = 0
        !        SUM_rate = 0
        !    endif
        !    call TwoNormTensor_strain(Erate,6,dummyvar(1)) ! Recalculate the two norm so the updated Erate value can be tracked
        !    TrackSmoothStrainRate = dummyvar(1) !Track the smoothed strain rate
        !endif
                
        !________________________________________________________________________________
        
        !Update state variables due to HSR
        !Store state variables
        Gu=G
        Ku=K
        eta_yu=eta_y
        Du=Dp        
        
        !print *, "Before"
        !print *, "I_coeff: ", I_coeff
        !print *, "I_0: ", I_0
        !print *, "ApplyStrainRateUpdates: ", ApplyStrainRateUpdates
        
        call Get_strain_invariants(Erate, dummyvar(1), epsq_rate)! deviatoric strain rate
        TrackShearStrainRate = epsq_rate ! Track shear strain rate
        call Get_I_coeff(D_part, G_s, p, epsq_rate, I_act)!actual inertial coefficient
        TrackInertialCoefficient = I_act ! Track the intertial coefficient
        !I_act = 0 ! Test to see if this is the reason the models get different results for no updates and D_part = 0
        dI=I_act-I_coeff !change in inertial coefficient
        call check4crossing(I_coeff, I_act, dI, I_0, ApplyStrainRateUpdates) !Check if update needed    
        
        !print *, "After"
        !print *, "I_coeff: ", I_coeff
        !print *, "I_act: ", I_act
        !print *, "dI", dI
        !print *, "I_0: ", I_0
        !print *, "ApplyStrainRateUpdates: ", ApplyStrainRateUpdates
        
        if (applystrainrateupdates) then !update
            !h=h*(i_act/i_0)**alpha_g
            call update_gk(g_0, nu, i_act, i_0, alpha_g, alpha_k, gu, ku) !updates the elastic properties
            call get_dp(h, d_min, i_act, i_0, epsq_p, alpha_d, applystrainrateupdates, du) !updates the dilation
            eta_yu=m-du*(1.0-n) !updates eta
        endif
        !_________________________________________________________________________________
        
		! Fill elastic material matrix
		D1  = Ku+(4*Gu/3)
		D2  = Ku-(2*Gu/3)
		DE  = 0.0
		DE(1:3,1:3) = D2
		DE(1,1) = D1
		DE(2,2) = D1
		DE(3,3) = D1
		DE(4,4) = Gu
		DE(5,5) = Gu
		DE(6,6) = Gu
               
        !Get elastic predictor stress
        call MatVec(DE,6,dEps,6,dSig_el)
        call AddVec(Sig_0, dSig_el, 1.0, 1.0, 6, Sig_t)
        !Get new invariant stresses
        call Get_invariants(Sig_t, p_t, q_t, theta)
        !Evaluate yield function
        call YieldFunction(q_t, p_t, eta_yu, FT)
        
        !Track variables
        !Trackq_t = q_t
        !Trackp_t = p_t
        TrackEta_y = eta_y ! Track eta_y for plotting
        TrackDp = Dp

        !___________________________________________________________________________
        !Now check elastic loading, unloading
		if (FT<-FTOL) then !Elastic behavior
           !Update state parameters
            G=Gu
            K=Ku
            eta_y=eta_yu
            Dp=Du
            I_coeff=I_act

            !Update stress
            Sig=Sig_t
            switch_yield=.false.
            
            !Track Loading Variables
            DeformCateg = 100
            !print *, FT
            TrackFT = FT
        !___________________________________________________________________________
		!*************************  Plastic behavior  ****************************** 
        !***************************************************************************
        !***************************************************************************
            
        !___________________________________________________________________________
        !***********Checking surface bisection and viscoplastic unloading***********
        !===========================================================================
        else !plastic behavior            
            if (F0<-FTOL) then !Elasto-plastic transition                
                call Newton_Raphson(G, K, eta_y, Dp, G_0, nu, M, M_tc, N, D_min, h, D_part, &
									G_s, epsq_p, I_coeff, I_act, I_0, alpha_K, alpha_G, alpha_D,&
									Gu, Ku, eta_yu, Du, dEps, MAXITER, FTOL,&
									F0, Sig_0, alpha)
                !Track Deformation category
                DeformCateg = 200
                

            else!pure plastic deformations or unloading
                call Check_Unloading(M_tc, eta_y, eta_yu, dI, Sig_0, dSig_el, LTOL, IsUnloading)  !Determines if is unloading path
                if (IsUnloading) then !Find new alpha
                    call Newton_Raphson(G, K, eta_y, Dp, G_0, nu, M, M_tc,N, D_min, h, D_part, &
										G_s, epsq_p, I_coeff, I_act, I_0, alpha_K, alpha_G, alpha_D,&
										Gu, Ku, eta_yu, Du, dEps, MAXITER, FTOL,&
										F0, Sig_0, alpha)
                    !Track Deformation category
                    DeformCateg = 300
                
                else !Pure plasticity
                    alpha=0.0d0
                    !Track Deformation category
                    DeformCateg = 400
                endif                
            endif 
        !______________________________________________________________________________________
            
        !______________________________________________________________________________________
        !*************************Substepping algorithm****************************************
        !======================================================================================
            Sig=Sig_0 !Initialize stress
            dEps_t=(1.0-alpha)*dEps !remaining deformation
            dI_t=(1.0-alpha)*dI !remaining inertial energy
            counter=0
            T=0.0d0 !pseudo-time
            DT=1.0d0 !pseudo-time increment
            
            I_TT=I_coeff
            !print *, "Noel:", NOEL
            !print *, "N_i:", N_i
            SubStepping_MaxIter = 0
            do while (T<1.0d0) !main loop
                !SubStepping_MaxIter = SubStepping_MaxIter+1
                !print * , "T: ", T !Print out T for substepping algorithm
                Sig_t=Sig!Store initial stress tensor
                I_coeff=I_TT
                dEps_TT=DT*dEps_t!substepping
                dI_TT=DT*dI_t!rate substepping
                I_TT=I_coeff+dI_TT
                
                !_________________________________________________________________________________
                !First approximation
                
                !Store parameters
                G1=G
                K1=K
                eta_y1=eta_y
                Dp1=Dp
                Sig1=Sig_t
                EpsPt=EpsP
                !call the Euler's algorithm with input parameters
                call Euler_Algorithm(G_0, nu, M_tc, M, N, D_min, h, D_part, G_s,&
                                     G1, K1, eta_y1, Dp1, &	
                                     erate, I_0, I_coeff, dI_TT, alpha_G, alpha_K, alpha_D, DTIME, DT, &
                                     switch_original, &
                                     Sig1, EpsPt, dEps_TT, dD1, dEpsp1, dSig1)
                !_________________________________________________________________________________
                
                !=================================================================================
                !Store max F1 comment if not needed
                call Get_invariants(Sig1, p, q, theta)
                call YieldFunction(q, p, eta_y1, FT)
				if (abs(FT)>abs(Error_yield_1)) Error_yield_1=abs(FT)
                !=================================================================================
                
                !_________________________________________________________________________________
                !second approximation
                
                !Store parameters
                G2=G1
                K2=K1
                eta_y2=eta_y1
                Dp2=Dp1
                Sig2=Sig1
                !call the Euler's algorithm with input parameters
                call Euler_Algorithm(G_0, nu, M_tc, M, N, D_min, h, D_part, G_s,&
                                     G2, K2, eta_y2, Dp2, &	
                                    erate, I_0, I_TT, dI_TT, alpha_G, alpha_K, alpha_D, DTIME, DT, &
                                     switch_original, &
                                     Sig2, EpsPt, dEps_TT, dD2, dEpsp2, dSig2)
                !_________________________________________________________________________________
                
                !=================================================================================
                !Store max F2 comment if not needed
                call Get_invariants(Sig2, p, q, theta)
                call YieldFunction(q, p, eta_y2, FT)
				if (abs(FT)>abs(Error_yield_2)) Error_yield_2=abs(FT)
                !=================================================================================
                
                !_________________________________________________________________________________
                !Compute averages and error
                
                Sig_t=Sig_t+0.5*(dSig1+dSig2)!updated stress
				Dp1=Dp+0.5*(dD1+dD2) !Updated dilation
				eta_y1=M-Dp1
				call TwoNormTensor((dSig1-dSig2), 6, dummyvar(1)) !||Delta Sigma||
                call TwoNormTensor(Sig_t, 6, dummyvar(2)) !||Sig_T+DT||
                R_TT=0.5*max(dummyvar(1)/dummyvar(2), abs(dD1-dD2)/abs(eta_y1)) !Relative residual error
                !________________________________________________________________________________
                
                !=================================================================================
                !Store max rel residual comment if not needed
                 ! if (R_TT>Error_Euler_max) Error_Euler_max= R_TT  				            
                !=================================================================================
                
                
                if ((R_TT>STOL).and.(counter<=MAXITER)) then!Step failed
                    counter=counter+1
                    qR=max((0.9*sqrt(STOL/R_TT)), 0.1)
                    DT=max(qr*DT, DTmin)
                    Failed=.true.
                    I_TT=I_coeff
                    
                else !successful step
                    !___________________________________________________________________________
                    !Update plastic strain, stress, and state variables
                    Sig=Sig_t
                    EpsP=EpsP+ 0.5*(dEpsp1+dEpsp2)
                    G=0.5*(G1+G2)
                    K=0.5*(K1+K2)
					call Get_strain_invariants(EpsP, eps_v, epsq_p)
					call Get_Dp(h, D_min, I_TT, I_0,  epsq_p, alpha_D, ApplyStrainRateUpdates, Dp)
                    eta_y=M-Dp*(1.0-N)
                    !__________________________________________________________________________
                    
                    !________________________________________________________________________________
                    !***********************Stress drift correction**********************************
					!________________________________________________________________________________
                    call Get_invariants(Sig, p, q, theta) !stress invariants
                    call YieldFunction(q, p, eta_y, F0) !Initial drift
                    !print *, "F0:", F0
                    !=================================================================================
                    !Store last F0, comment if not needed
					if (abs(F0)>abs(Error_Yield_max)) Error_Yield_max=abs(F0)                          
                    !=================================================================================
                    
                    call Stress_Drift_Correction(G_0, nu, M_tc, M, N, D_min, h, D_part, G_s,&
                                                G, K, eta_y, Dp, &	
                                                Erate, I_0, I_TT, I_coeff, dI_TT, alpha_G, alpha_K, alpha_D, DTIME, DT, &
                                                switch_original, MAXITER, F0, FTOL, &
                                                Sig, EpsP, dEps_TT)
                    !_________________________________________________________________________________
                    
                    !=================================================================================
                    !Store last F0 comment if not needed
					if (abs(F0)>abs(Error_Yield_last)) Error_Yield_last=abs(F0)                           
                    !=================================================================================
                    
                    !_________________________________________________________________________________
                    !Final substepping if needed
                    
                    !Check to see if R_TT is close to zero
                    if (abs(R_TT)< RTOL) then
                        qr = 1.1
                    else
                        qR=min((0.9*sqrt(STOL/R_TT)), 1.1)
                    endif
                    
                    if (Failed) then
                        qR=min(qR, 1.0)
                        Failed=.false.
                    endif
                    DT=qR*DT
                    T=T+DT
                    DT=max(DT, DTmin)
                    DT=min(DT, 1.0-T)
                    switch_yield=.true. !plastic point                    
                    !_________________________________________________________________________________                    
                end if               
            end do
            I_coeff=I_act
        endif
        !Track variables
        TrackShearModulus = Gu
        call Get_invariants(Sig, p, q, theta) ! Recalculate invariants so that they can be stored
        Trackq_t = q
        Trackp_t = p
    end subroutine NAMC_HSR
    !*******************************************************************************************
    
    
    !___________________________________________________________________________________________
      !     ##    ## ######## ##      ## ########  #######  ##    ## 
      !     ###   ## ##       ##  ##  ##    ##    ##     ## ###   ## 
      !     ####  ## ##       ##  ##  ##    ##    ##     ## ####  ## 
      !     ## ## ## ######   ##  ##  ##    ##    ##     ## ## ## ## 
      !     ##  #### ##       ##  ##  ##    ##    ##     ## ##  #### 
      !     ##   ### ##       ##  ##  ##    ##    ##     ## ##   ### 
      !     ##    ## ########  ###  ###     ##     #######  ##    ## 
      !########     ###    ########  ##     ##  ######   #######  ##    ## 
      !##     ##   ## ##   ##     ## ##     ## ##    ## ##     ## ###   ## 
      !##     ##  ##   ##  ##     ## ##     ## ##       ##     ## ####  ## 
      !########  ##     ## ########  #########  ######  ##     ## ## ## ## 
      !##   ##   ######### ##        ##     ##       ## ##     ## ##  #### 
      !##    ##  ##     ## ##        ##     ## ##    ## ##     ## ##   ### 
      !##     ## ##     ## ##        ##     ##  ######   #######  ##    ##
    
    subroutine Newton_Raphson(G, K, eta_y, Dp, G_0, nu, M, M_tc, No, D_min, h, D_par, &
							  G_s, eps_q, I_coeff, I_act, I_0, k_K, k_G, k_D,&
                              Gu, Ku, eta_yu, Du, dEps, MAXITER, FTOL, &
							  F0, Sig_0, alpha)
    !**************************************************************************************
	!  This subroutine determine the elastic proportion with consistent state and elastic *
	!  parameters                                                                         *
	!**************************************************************************************
    implicit none
    !input
    integer, intent(in):: MAXITER
    double precision, intent(in):: G_0, nu, M, M_tc, No, D_min, h, D_par, G_s, &
                                   eps_q, k_G, k_K, k_D, I_act, I_0, FTOL
    double precision, dimension(6), intent(in):: dEps
    !output
    double precision, intent(inout):: G, K, eta_y, Dp, I_coeff, &
									  Gu, Ku, eta_yu, Du, F0, Sig_0(6)
    double precision, intent(out):: alpha
    !local variables
    logical:: ApplyStrainRateUpdates
    integer:: n, i
    double precision:: FT, dG, dK, deta, dD, dI, &
                       dI_alpha, I_alpha, &
					   p_alpha, q_alpha, dummyvar, F_prime
    double precision:: dEps_alpha(6), dSig_alpha(6), Sig_alpha(6), &
                       n_vec(6), L, dSigdAlpha(6) 
    double precision:: D1, D2, DE(6,6), dDdG(6,6), dDdK(6,6), aux(6,6)
    
    !Initialize parameters
    FT=1000
    alpha=0.5d0
    n=0
    !Determine Changes in state parameters
    dG=Gu-G
    dK=Ku-K
    deta=eta_yu-eta_y
    dD=Du-Dp
    dI=I_act-I_coeff    
    do while ((abs(FT)>=FTOL).and.(n<=MAXITER))
        !Store variables
        Gu=G
        Ku=K
        eta_yu=eta_y
        Du=Dp
        F0=FT
        
        n=n+1
        !___________________________________________________________________
        !Compute trial strains and stresses
        dEps_alpha=alpha*dEps
        dI_alpha=alpha*dI
        I_alpha=I_coeff+dI_alpha
        !___________________________________________________________________
        !Evaluate rate crossing
        call check4crossing(I_coeff, I_alpha, dI_alpha,I_0, ApplyStrainRateUpdates)
        if (ApplyStrainRateUpdates) then !Update parameters
            call Update_GK(G_0, nu, I_alpha, I_0, k_G, k_K, Gu, Ku)
            call Get_Dp(h, D_min, I_alpha, I_0, eps_q, k_D, ApplyStrainRateUpdates, Du)
            eta_yu=M-Du*(1.0-No)
        endif
        !___________________________________________________________________
        !Update trial stress 
        !Ensemble elastic matrix
		D1  = Ku+(4*Gu/3)
		D2  = Ku-(2*Gu/3)
		DE  = 0.0
		DE(1:3,1:3) = D2
		DE(1,1) = D1
		DE(2,2) = D1
		DE(3,3) = D1
		DE(4,4) = Gu
		DE(5,5) = Gu
		DE(6,6) = Gu  
        call MatVec(DE,6,dEps_alpha,6,dSig_alpha)
        call AddVec(Sig_0, dSig_alpha, 1.0, 1.0, 6, Sig_alpha)
        !___________________________________________________________
        !Evaluate yield function
        call Get_invariants(Sig_alpha, p_alpha, q_alpha, dummyvar)
        call YieldFunction(q_alpha, p_alpha, eta_yu, FT)
        
        !Evaluate n=dF/dSig and L=dF/dXs
        call Get_dF_to_dSigma(M_tc, eta_yu, Sig_alpha, n_vec)!dF/dSig
        L=-p_alpha*(1.0-No) !L=dF/dXs        
        !____________________________________________________________
        !Evaluate F'(alpha)
        !Elastic matrix derivatives
        dDdK=0.0
		dDdK(1:3,1:3) = 1.0
		dDdG=0.0
		dDdG(1:3,1:3)=-2./3.
		dDdG(1,1)=4./3.
		dDdG(2,2)=4./3.
		dDdG(3,3)=4./3.
		dDdG(4,4)=1.0
		dDdG(5,5)=1.0
		dDdG(5,5)=1.0
		aux=alpha* (dK*dDdK+dG*dDdG)
		aux=DE+aux
		call MatVec(aux, 6, dEps, 6, dSigdAlpha)
        F_prime=0.0d0
		do i=1,6 !dot product
			F_prime=F_prime+n_vec(i)*dSigdAlpha(i)
        enddo
        F_prime=F_prime+L*dD!F'(alpha)
        !______________________________________________________________
        !Update alpha
        alpha=alpha-FT/F_prime    
        if ((alpha>1.0d0).or.(alpha<0.0d0)) then
			alpha=0.0d0
		endif  
    end do
    !Update variables to return
    G=Gu
    K=Ku
    eta_y=eta_yu
    Dp=Du
    I_coeff=I_alpha
    Sig_0=Sig_alpha    
    end subroutine Newton_Raphson
    !*********************************************************************************
    
!______________________________________________________________________________________    
!######## ##     ## ##       ######## ########      ######  
!##       ##     ## ##       ##       ##     ##    ##    ## 
!##       ##     ## ##       ##       ##     ##    ##       
!######   ##     ## ##       ######   ########      ######  
!##       ##     ## ##       ##       ##   ##            ## 
!##       ##     ## ##       ##       ##    ##     ##    ## 
!########  #######  ######## ######## ##     ##     ######  
!   ###    ##        ######    #######  ########  #### ######## ##     ## ##     ## 
!  ## ##   ##       ##    ##  ##     ## ##     ##  ##     ##    ##     ## ###   ### 
! ##   ##  ##       ##        ##     ## ##     ##  ##     ##    ##     ## #### #### 
!##     ## ##       ##   #### ##     ## ########   ##     ##    ######### ## ### ## 
!######### ##       ##    ##  ##     ## ##   ##    ##     ##    ##     ## ##     ## 
!##     ## ##       ##    ##  ##     ## ##    ##   ##     ##    ##     ## ##     ## 
!##     ## ########  ######    #######  ##     ## ####    ##    ##     ## ##     ## 
    
    subroutine Euler_Algorithm(G_0, nu, M_tc, M, No,  D_min, h, Dpart, Gs,&
                               G, K, eta_y, Dp, &	
                               erate, I_0, I, dI, k_G, k_K, k_D, dtime, DT, &
                               switch_original, &
                               Sig, EpsP, dEps, dD, dEpsp, dSig)
    !************************************************************************
	! Euler's algorithm (finite difference) integration						*
	! Returns:																*
	! Increment of stresses and change of state parameters                  *                           
	!************************************************************************  
    implicit none
    !input
    logical, intent(in):: switch_original
    double precision, intent(in):: G_0, nu, M_tc, M, No, D_min, h, Dpart, Gs,&
                                   erate(6), I_0, k_G, k_K, k_D, dtime, DT,&
                                   dEps(6)
    double precision, intent(inout):: I, dI
    !output
    double precision, intent(inout):: G, K, eta_y, Dp, EpsP(6) , Sig(6)
    double precision, intent(out):: dD, dSig(6), dEpsp(6)
    !local variables
    logical:: ApplyStrainRateUpdate
    double precision:: I_f, DE(6,6), D1, D2, dSig_el(6), num, lambda, &
                       n_vec(6), p, q, dummyvar, L, b, m_vec(6), a(6), &
					   epsq_p, epsv_p, Hard, den, dummyvec(6), Hvp, erate_q

    !________________________________________________________________________
    !Evaluate for strain rate updating
    I_f=I+dI
    call check4crossing(I,  I_f, dI, I_0, ApplyStrainRateUpdate)
    
    if (ApplyStrainRateUpdate) then
        call Update_GK(G_0, nu, I, I_0, k_G, k_K, G, K)
	endif
    !________________________________________________________________________
	
	!________________________________________________________________________
	!Compute invariants and derivatives
	call Get_invariants(Sig, p, q, dummyvar)
	call Get_strain_invariants(EpsP, epsv_p, epsq_p)
	call Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec) !n=dFdSig
	call Get_dP_to_dSigma(Dp, Sig, m_vec) !m=dP/dSig
    L = -p*(1.0-No) !L=dF/dXs
	call Get_dD_to_dEpsq(D_min, h, I_0, k_D, epsq_p, epsv_p, &
                         EpsP, I, ApplyStrainRateUpdate, a) !a=dD/dEpsq^p
    call Get_dD_to_dI(D_min, h, I_0, k_D, epsq_p, I, b) !b=dXs/dI in place of dXs/dErate
	!________________________________________________________________________
    
    !________________________________________________________________________
    !Compute elastic predictor	
	!Ensemble elastic matrix    
    D1  = K+(4*G/3)
	D2  = K-(2*G/3)
	DE  = 0.0
	DE(1:3,1:3) = D2
	DE(1,1) = D1
	DE(2,2) = D1
	DE(3,3) = D1
	DE(4,4) = G
	DE(5,5) = G
	DE(6,6) = G 
    
    !Calculate the predictor
    call MatVec(DE, 6, dEps, 6, dSig_el)
    !________________________________________________________________________
    
    !________________________________________________________________________
    !compute numerator
	call DotProduct_2(dSig_el, n_vec, 6, num)! dSig_el . n_vec  
    
    if ((ApplyStrainRateUpdate).and.(.not.switch_original)) then !Dashpot visc.
        num=num+L*b*dI
    endif
    !_______________________________________________________________________
    
    !_______________________________________________________________________
    !Compute denominator
    
    !Compute hardening modulus H=-L.a.m
    call DotProduct_2(a, m_vec, 6, dummyvar)
    Hard=-L*dummyvar !hardening modulus
    
    !compute Den=n.D.m +H
    call MatVec(DE, 6, m_vec, 6, dummyvec)
    call DotProduct_2(n_vec, dummyvec, 6, den)
    den=den+Hard
    
    if ((ApplyStrainRateUpdate).and.(switch_original)) then !original
    !compute viscous hardening modulus Hvp=-L.bv.m/(dt*DT)
		call Get_strain_invariants(erate,dummyvar,erate_q)
		call Get_dEpsq_to_dEps(erate_q,erate,dummyvec)
		dummyvec=b*Dpart*sqrt(Gs/abs(p))*dummyvec
        call DotProduct_2(dummyvec, m_vec, 6, dummyvar)
        Hvp=-L*dummyvar/(dtime*DT)
        den=den+Hvp
    endif
    !_____________________________________________________________________
    
    !_____________________________________________________________________
    !compute viscoplastic mult. lambda, dEpsp, dSig, and dD
    
    lambda=num/den! viscoelastic multiplier    
    dEpsp=lambda*m_vec! Flow rule
    dummyvec=dEps-dEpsp
    call MatVec(DE, 6, dummyvec, 6, dSig) !stress increment
    dD=0.0
    call DotProduct_2(a, dEpsp, 6, dD)!plastic hard/softening
    dD=dD+b*dI !rate hard/softening
    !_____________________________________________________________________
    
    !_____________________________________________________________________
    !Update EpsP, Sig, D, and eta_y
    EpsP=EpsP+dEpsp
    Sig=Sig+dSig
    Dp=Dp+dD
    eta_y=M-Dp*(1.0-No)
    !____________________________________________________________________    
    end subroutine Euler_Algorithm
!*************************************************************************
    
    
!________________________________________________________________________   
    
! ######  ######## ########  ########  ######   ######  
!##    ##    ##    ##     ## ##       ##    ## ##    ## 
!##          ##    ##     ## ##       ##       ##       
! ######     ##    ########  ######    ######   ######  
!      ##    ##    ##   ##   ##             ##       ## 
!##    ##    ##    ##    ##  ##       ##    ## ##    ## 
! ######     ##    ##     ## ########  ######   ######  
!########  ########  #### ######## ######## 
!##     ## ##     ##  ##  ##          ##    
!##     ## ##     ##  ##  ##          ##    
!##     ## ########   ##  ######      ##    
!##     ## ##   ##    ##  ##          ##    
!##     ## ##    ##   ##  ##          ##    
!########  ##     ## #### ##          ## 
    
    subroutine Stress_Drift_Correction(G_0, nu, M_tc, M, No, D_min, h, Dpart, Gs, &
                                       G, K, eta_y, Dp, &	
                                       erate, I_0, I_f, I, dI, k_G, k_K, k_D, dtime, DT, &
                                       switch_original, MAXITER, F0, FTOL, &
                                       Sig, EpsP, dEps)
    !************************************************************************
	! Corrects any drift caused during the Mod. Euler's procedure       	*
    ! Returns corrected stresses and state parameters                       *                 
	!************************************************************************
    implicit none
    !input
    logical, intent(in):: switch_original
    integer, intent(in):: MAXITER
    double precision, intent(in):: G_0, nu, M_tc, M, No, D_min, h, Dpart, Gs,&
                                   erate(6), I_0, k_G, k_K, k_D, dtime, DT,&
                                   FTOL, dEps(6)
    double precision, intent(inout):: I, I_f, dI
    !output
    double precision, intent(inout):: G, K, eta_y, Dp, EpsP(6) , Sig(6), F0
    !local variables
    logical:: ApplyStrainRateUpdate
    integer:: n
    double precision::epsq_p, epsv_p, n_vec(6), m_vec(6), L, a(6), b,&
                      DE(6,6), D1, D2, Den,  Hard, Hvp, dlambda, &
                      dummyvec(6), dummyvar, FC, &
                      Du, eta_yu, Sigu(6), dD, dSig(6), dEpsp(6),Epspu(6), &
                     p, q, erate_q
    !________________________________________________________________________
    !Evaluate for strain rate updating    
    call check4crossing(I, I_f, dI, I_0, ApplyStrainRateUpdate)
    
    if (ApplyStrainRateUpdate) then
        call Update_GK(G_0, nu, I_f, I_0, k_G, k_K, G, K)
    endif
    !________________________________________________________________________
    
    !________________________________________________________________________
    !Initialize values
    n=0 !counter
    FC=F0    
    !________________________________________________________________________
    
    do while ((abs(FC)>FTOL).and.(n<MAXITER))
        !___________________________________________________________________
        !update values
    	F0=FC
        n=n+1
        !___________________________________________________________________
        
        !________________________________________________________________________
	    !Compute invariants and derivatives
	    call Get_invariants(Sig, p, q, dummyvar)
	    call Get_strain_invariants(EpsP, epsv_p, epsq_p)
	    call Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec) !n=dFdSig
	    call Get_dP_to_dSigma(Dp, Sig, m_vec) !m=dP/dSig
        L=-p*(1.0-No) !L=dF/dXs
	    call Get_dD_to_dEpsq(D_min, h, I_0, k_D, epsq_p, epsv_p, &
                             EpsP, I_f, ApplyStrainRateUpdate, a) !a=dD/dEpsq^p
        call Get_dD_to_dI(D_min, h, I_0, k_D, epsq_p, I_f, b) !b=dXs/dI in place of dXs/dErate
	    !________________________________________________________________________
        
        !___________________________________________________________________
        !compute denominator
        
        !Ensemble elastic matrix    
        D1  = K+(4*G/3)
	    D2  = K-(2*G/3)
	    DE  = 0.0
	    DE(1:3,1:3) = D2
	    DE(1,1) = D1
	    DE(2,2) = D1
	    DE(3,3) = D1
	    DE(4,4) = G
	    DE(5,5) = G
	    DE(6,6) = G 
        
        !Compute hardening modulus H=-L.a.m
        call DotProduct_2(a, m_vec, 6, dummyvar)
        Hard=-L*dummyvar !hardening modulus
    
        !compute Den=n.D.m +H
        call MatVec(DE, 6, m_vec, 6, dummyvec)
        call DotProduct_2(n_vec, dummyvec, 6, den)
        den=den+Hard
    
        if ((ApplyStrainRateUpdate).and.(switch_original)) then !original
        !compute viscous hardening modulus Hvp=-L.bdI/dErate.m/(dt*DT)
            call Get_strain_invariants(erate,dummyvar,erate_q)
		    call Get_dEpsq_to_dEps(erate_q,erate,dummyvec)
		    dummyvec=b*Dpart*sqrt(Gs/abs(p))*dummyvec
            call DotProduct_2(dummyvec, m_vec, 6, dummyvar)
            Hvp=-L*dummyvar/(dtime*DT)
            den=den+Hvp
        endif        
        !___________________________________________________________________
        
        !___________________________________________________________________
        !Compute DeltaLambda=F0/den and updated Sig, D, and eta_y
        dlambda=F0/den      !error of plastic multiplier
        dEpsp=dlambda*m_vec !error in plastic strain tensor
		Epspu=EpsP-dEpsp
        call MatVec(DE, 6, dEpsp, 6, dSig) !stress drift
        dD=0.0
        call DotProduct_2(a, dEpsp, 6, dD)!plastic hard/softening
		!if (switch_original) dD=dD+b*dI
        
        !update
        Sigu=Sig-dSig
        Du=Dp-dD
        eta_yu=M-Du*(1.0-No)
        !__________________________________________________________________
        
        !__________________________________________________________________
        !Evaluate yield function
        call Get_invariants(Sigu, p, q, dummyvar)
        call YieldFunction(q, p, eta_yu, FC)
        !__________________________________________________________________
        !Evaluate change direction
        if (abs(FC)>abs(F0)) then
            call DotProduct_2(n_vec, n_vec, 6, den)
        	dlambda=F0/den
            dEpsp=dlambda*m_vec !error in plastic strain tensor
			Epspu=EpsP-dEpsp
            call MatVec(DE, 6, dEpsp, 6, dSig) !stress drift
            Sigu=Sig-dSig
            Du=Dp
            eta_yu=eta_y
		end if
        !__________________________________________________________________
        
        !__________________________________________________________________
        !Update the state parameters
		EpsP=Epspu
        Sig=Sigu
        Dp=Du
        eta_y=eta_yu
        !__________________________________________________________________
    end do
    
    
    end subroutine Stress_Drift_Correction
!******************************************************************************
    
    
    
    
    
!______________________________________________________________________________  
!##     ##    ###    ########  ########  
!##     ##   ## ##   ##     ## ##     ## 
!##     ##  ##   ##  ##     ## ##     ## 
!######### ##     ## ########  ##     ## 
!##     ## ######### ##   ##   ##     ## 
!##     ## ##     ## ##    ##  ##     ## 
!##     ## ##     ## ##     ## ########    
! ######   #######  ######## ########       ###    ##    ## ########  
!##    ## ##     ## ##          ##         ## ##   ###   ## ##     ## 
!##       ##     ## ##          ##        ##   ##  ####  ## ##     ## 
! ######  ##     ## ######      ##       ##     ## ## ## ## ##     ## 
!      ## ##     ## ##          ##       ######### ##  #### ##     ## 
!##    ## ##     ## ##          ##       ##     ## ##   ### ##     ## 
! ######   #######  ##          ##       ##     ## ##    ## ########  
!########  ######## ########  #### ##     ##    ###    ######## #### ##     ## 
!##     ## ##       ##     ##  ##  ##     ##   ## ##      ##     ##  ##     ## 
!##     ## ##       ##     ##  ##  ##     ##  ##   ##     ##     ##  ##     ## 
!##     ## ######   ########   ##  ##     ## ##     ##    ##     ##  ##     ## 
!##     ## ##       ##   ##    ##   ##   ##  #########    ##     ##   ##   ##  
!##     ## ##       ##    ##   ##    ## ##   ##     ##    ##     ##    ## ##   
!########  ######## ##     ## ####    ###    ##     ##    ##    ####    ###    
!########  ######  
!##       ##    ## 
!##       ##       
!######    ######  
!##             ## 
!##       ##    ## 
!########  ###### 
    
	subroutine Get_Dp(h, D_min, I, I_0, eps_q, k, ApplyRateUpdating, D)
	!*********************************************************************
	! Returns the dilation for current inertial coefficient and dev.     *
	! strain															 *
	!*********************************************************************
    implicit none
    logical, intent(in):: ApplyRateUpdating
    double precision, intent(in):: h, D_min, I, I_0, eps_q, k
    !out
    double precision, intent(out):: D
    !local variables
    double precision:: D_mm
    if (ApplyRateUpdating) then
		D_mm=D_min*(I/I_0)**k !strain/rate hardening
    else
        D_mm=D_min
    endif
    
    D=h*D_mm*eps_q*exp(1.0-h*eps_q) !hardening rule
    end subroutine Get_Dp
    
    subroutine Update_GK(G_0, nu, I, I_0, k_G, k_K, G, K)
	!*********************************************************************
	! Returns updated elastic modulus                                    *
	!																	 *
	!*********************************************************************
    implicit none
    !input
    double precision, intent(in):: G_0, nu, I, I_0, k_G, k_K
    !output
    double precision, intent(out):: G, K
    !local variables
    double precision:: K_0
    G=G_0*(I/I_0)**k_G! updated modulus
    
    K_0=2*G_0*(1+nu)/(3*(1-2*nu))!bulk modulus
    K=K_0*(I/I_0)**k_K! updated modulus
    end subroutine Update_GK
    
    
    subroutine Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec)
	!************************************************************************
	! Returns the derivative of the yield function with respect to the		*
	! stress tensor 														*
    ! n=dF/dSigma =dF/dp*dp/dSigma+ dF/dq*dq/dSigma +dF/dtheta*dtheta/dSigma*
    ! n is a (1X6) vector													*                           
	!************************************************************************  
    implicit none
    !input
    double precision, intent(in):: M_tc, eta_y, Sig(6)
    !output
    double precision, dimension(6):: n_vec
    !local variables
    double precision:: p, q, theta, pi=2.0*acos(0.0d0), &
                       J2, J3, dJ3dsig(6), dfdtheta, &
					   dpdsig(6), dqdsig(6), dev(6), dev2(6), &
					   TrS2, II(6), dthetadSig(6), COS_3THETA
    !Get the invariants
    call Get_invariants(Sig, p, q, theta)
    !Get dF/dp=eta_y and dF/dq=1
    !Get dF/dtheta
    dfdtheta=0.45*p*M_tc*((cos(1.5d0*theta+0.25d0*pi))**0.2)*sin(1.5d0*theta+0.25d0*pi)
    !___________________________________________________________________________
    !1) Get dp/dSig=1/3 Imat
    dpdsig=0.0d0
    dpdsig(1)=1.0/3.0
    dpdsig(2)=1.0/3.0
    dpdsig(3)=1.0/3.0
    
    !2) Get dq/dsig= 2 *dev/3*q
    dev=Sig
    dev(1)=dev(1)-p
    dev(2)=dev(2)-p
    dev(3)=dev(3)-p
    
    dqdSig=(3.0/(2.0*q))*dev
    
    !3) Get dtheta/dSigma= (1/3cos3theta) d/dsigma((J3/2) * (3/J2)^1.5)
    J2=(q**2)/3.0
    J3=dev(1)*dev(2)*dev(3)-dev(1)*dev(6)**2-dev(2)*dev(4)**2-dev(3)*dev(5)**2+2.0*dev(4)*dev(5)*dev(6)
    !Fill S.S
    dev2(1)=dev(1)**2+dev(4)**2+dev(5)**2
    dev2(2)=dev(2)**2+dev(4)**2+dev(6)**2
    dev2(3)=dev(3)**2+dev(5)**2+dev(6)**2
    dev2(4)=dev(4)*(dev(1)+dev(2))+dev(5)*dev(6)
    dev2(5)=dev(5)*(dev(1)+dev(3))+dev(4)*dev(6)
    dev2(6)=dev(6)*(dev(2)+dev(3))+dev(4)*dev(5)
    !Compute dJ3dSig
    TrS2=dev2(1)+dev2(2)+dev2(3)  
    II=0.0d0!Identity tensor
    II(1)=1.0
    II(2)=1.0
    II(3)=1.0
    dJ3dsig=dev2-(TrS2*II/3.0d0)
    !Compute dtheta/dsig
    
    dthetadSig=dJ3dsig-(1.5*J3/J2)*dev
    COS_3THETA=cos(3.0*theta)
    dthetadSig=(sqrt(3.0)/(2.0*COS_3THETA*J2**1.5))*dthetadSig
    !__________________________________________________________________
    !Get n_vec=dF/dSig
	n_vec=(eta_y*dpdsig)+dqdSig+(dfdtheta*dthetadSig) !n_vec=dF/dSig
    end subroutine Get_dF_to_dSigma
    
    
    subroutine Get_dD_to_dI(D_min, h, I_0, kD, eps_q, I, b)
    !************************************************************************
	! Returns the derivative of the Dilation with respect to the inertial	*
	! coefficient 															*
    ! b=dD/dI																*
    ! b is a scalar															*                           
	!************************************************************************ 
    implicit none
    !input
    double precision, intent(in):: D_min, h, I_0, kD, eps_q, I
    !output
    double precision, intent(out)::b
    !local variables
    
    b=h*D_min*eps_q*exp(1-h*eps_q)*kD*((I/I_0)**(kD-1.0))/I_0
       
    end subroutine Get_dD_to_dI
    
    subroutine Get_dP_to_dSigma(D, Sig, m_vec)
	!************************************************************************
	! Returns the derivative of the plastic potential function with respect *
	! to the stress tensor													*
    ! m=dP/dSigma =dP/dp*dp/dSigma+ dP/dq*dq/dSigma							*
    ! m is a (1X6) vector													*                           
	!************************************************************************  
    implicit none
    !input
    double precision, intent(in):: D, Sig(6)
    !output
    double precision, dimension(6):: m_vec
    !local variables
    double precision:: p, q, theta, pi=2.0*acos(0.0d0), &
                       J2, J3, dJ3dsig(6), dfdtheta, &
					   dpdsig(6), dqdsig(6), dev(6), dev2(6), &
					   TrS2, II(6), dthetadSig(6), COS_3THETA
    !Get the invariants
    call Get_invariants(Sig, p, q, theta)
    !Get dP/dp=-D and dF/dq=1
    !___________________________________________________________________________
    !1) Get dp/dSig=1/3 Imat
    dpdsig=0.0d0
    dpdsig(1)=1.0/3.0
    dpdsig(2)=1.0/3.0
    dpdsig(3)=1.0/3.0
    
    !2) Get dq/dsig= 2 *dev/3*q
    dev=Sig
    dev(1)=dev(1)-p
    dev(2)=dev(2)-p
    dev(3)=dev(3)-p
    
    dqdSig=(3.0/(2.0*q))*dev  
    
    !__________________________________________________________________
    !Get m_vec=dP/dSig
	m_vec=(-D*dpdsig)+dqdSig !m_vec=dP/dSig
    end subroutine Get_dP_to_dSigma
    
    
    
    subroutine Get_dD_to_dEpsq(D_min, h, I_0, k_D, epsq_p, epsv_p, &
						        EpsP, I, ApplyStrainRateUpdate, a)
   	!************************************************************************
	! Returns the derivative of the Dilation with respect to the plastic    *
	! strain																*
    ! a=dXs/dEpsp= dD/dEpsq * dEPsq/dEpsp									*
    ! a is a (1X6) vector													*                           
	!************************************************************************
    implicit none
	!input
	logical, intent(in):: ApplyStrainRateUpdate
	double precision, intent(in):: D_min, h, I_0, k_D, epsq_p, epsv_p, &
                                    EpsP(6), I
	!output
	double precision, intent(out):: a(6)
	!local variables
	double precision:: D, dDdEpsq_p, dev(6),dEpsq_pdEpsp(6) 
	
	!________________________________________________________________________
	!Get dD/dEpsq_p
	if (ApplyStrainRateUpdate) then
		D=D_min*(I/I_0)**k_D
	else
		D=D_min
	end if
	
	dDdEpsq_p=h*D*exp(1.0-h*epsq_p)*(1.0-h*epsq_p)
    !_______________________________________________________________________
    
    !_______________________________________________________________________
    !Get dEpsp_Q/dEpsp= 
	dev=EpsP
    dev(1)=dev(1)-(epsv_p/3.0)
    dev(2)=dev(2)-(epsv_p/3.0)
    dev(3)=dev(3)-(epsv_p/3.0) !deviatoric stress tensor
    
    if (epsq_p>0.0d0) then !in case of zero plastic strain
        dEpsq_pdEpsp=(2.0/(3.0*epsq_p))*dev
    else
        dEpsq_pdEpsp=0.0d0
    endif    
    !_______________________________________________________________________
    
    !_______________________________________________________________________
    !Get a=dXs/dEpsp
    
    a=dDdEpsq_p*dEpsq_pdEpsp
    !______________________________________________________________________    
	end subroutine Get_dD_to_dEpsq
	
	subroutine Get_dEpsq_to_dEps(Epsq, Eps, dEqdEpsq)
	!************************************************************************
	! Returns the derivative of the deviatoric strain with respect to the   *
	! deviatoric strain	tensor						     					*
    ! dEqdEpsq is a (1X6) vector											*                           
	!************************************************************************
	implicit none
	!input
	double precision, intent(in):: Epsq, Eps(6)
	!output
	double precision, intent(out):: dEqdEpsq(6)
	!local variables
	double precision:: evol, dev(6)
	
	evol=Eps(1)+Eps(2)+Eps(3)!vol strain
    
	dev=Eps
	dev(1)=dev(1)-evol/3.0
	dev(2)=dev(2)-evol/3.0
	dev(3)=dev(3)-evol/3.0 !deviatoric strain tensor
	
	if (Epsq>0.0d0) then !in case of zero plastic strain
        dEqdEpsq=(2.0/(3.0*Epsq))*dev
    else
        dEqdEpsq=0.0d0
    endif  	
	end subroutine Get_dEpsq_to_dEps
!**********************************************************************************
	
    
!_________________________________________________________________________________   
!######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######  
!##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ## 
!##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##       
!######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######  
!##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ## 
!##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ## 
!##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######  
 

    Subroutine Get_strain_invariants(Eps, Eps_v, Eps_q)
	!*********************************************************************
	! Takes the strain tensor and returns deviatoric and vol. strain     *
	!																	 *
	!********************************************************************* 
    implicit none
    !input
    double precision, dimension(6), intent(in):: Eps
    !output
    double precision, intent(out):: Eps_v, Eps_q
    !local variables
	double precision:: dev(6)
    Eps_v=Eps(1)+Eps(2)+Eps(3)! vol strain
    
	dev=Eps 
	dev(1)=dev(1)-(Eps_v/3.0)
	dev(2)=dev(2)-(Eps_v/3.0)
	dev(3)=dev(3)-(Eps_v/3.0)!deviatoric strain tensor
	
    call TwoNormTensor_strain(dev, 6, Eps_q)
    Eps_q=Eps_q*sqrt(2.0/3.0) ! dev strain
    end subroutine  Get_strain_invariants
    
    
    subroutine Get_invariants(Sig, p, q, theta)
	!*********************************************************************
	! Takes the stress tensor Sig and return invariants p, q, and theta  *
	!																	 *
	!*********************************************************************
	implicit none
	!input variables
	double precision, dimension(6), intent(in):: Sig
	!output variables
	double precision, intent(out)::p, q, theta
	!local variables
	double precision:: dev(6), J2, J3, sin3theta    

	p=(Sig(1)+Sig(2)+Sig(3))/3.0 !mean stress
	dev=Sig
	dev(1)=dev(1)-p !computes deviatoric stress tensor
	dev(2)=dev(2)-p
	dev(3)=dev(3)-p
	
    call TwoNormTensor(dev, 6, J2)
    J2=(J2**2)/2.0 !J_2 invariant
    q=sqrt(3*J2) ! deviatoric stress
    
    !J3 stress invariant
    J3=dev(1)*dev(2)*dev(3)-dev(1)*dev(6)**2-dev(2)*dev(4)**2-dev(3)*dev(5)**2+2.0*dev(4)*dev(5)*dev(6)
    
    !sin3theta
    if (J2>0.0d0) then
       sin3theta=0.5*J3*(3.0/J2)**(1.5d0) 
    else !Assume triaxial compression
       sin3theta=-1.0d0 
    endif
	if (sin3theta<-1.0) sin3theta=-1.0d0
    if (sin3theta>1.0) sin3theta=1.0d0
        
    
    theta=-asin(sin3theta)/3.0d0 !Lode's angle

    end subroutine Get_invariants
    
    
    subroutine Check_Unloading(M_tc, eta_y, eta_yu, dI, Sig, dSig,&
								LTOL, IsUnloading)
	!*********************************************************************
	! Returns true if stress path is viscoplastic unloading              *
	!																	 *
	!*********************************************************************  
    !input
    implicit none
    double precision, intent(in):: M_tc, eta_y, eta_yu, dI, &
                                   Sig(6), dSig(6), LTOL
    !output
    logical, intent(out):: IsUnloading
    !local variables
    double precision:: deta, n_vec(6), n_norm, Sig_norm,&
						dSIg_inner_n, beta, phi
    IsUnloading=.false.
    
    deta=eta_yu-eta_y!change in state parameter    
    call Get_dF_to_dSigma(M_tc, eta_yu, Sig, n_vec)!Normal to surface
    call TwoNormTensor(n_vec, 6, n_norm) !norm of n_vec
    call TwoNormTensor(dSig, 6, Sig_norm) !norm of dSig
    call TensorInnerProduct(dSig, n_vec, 6,dSIg_inner_n) !inner product
    
    beta=acos(deta/(n_norm*Sig_norm))!conical aperture is a plane for inviscid mat.
    phi=acos(dSIg_inner_n/(n_norm*Sig_norm))!angle between stress rate and normal
    
    if (phi-beta>LTOL) IsUnloading=.true. !condition for unloading
    
    end subroutine Check_Unloading
    
    
    subroutine Get_I_coeff(D_part, G_s, p, eps_rate, I)
	!*********************************************************************
	! Returns the inertial coefficient                                   *
	!																	 *
	!********************************************************************* 
    implicit none
    !input
    double precision, intent(in):: D_part, G_s, p, eps_rate
    !output
    double precision, intent(out):: I
    !local variables
    I=D_part*eps_rate*sqrt(G_s/abs(p))
    end subroutine Get_I_coeff
    
    subroutine Get_M(M_tc, theta, M)
	!*********************************************************************
	! Returns M															 *
	!																	 *
	!*********************************************************************
    implicit none
    !in
    double precision, intent(in):: M_tc, theta
    !out
    double precision, intent(out):: M
    !local
    double precision:: COS_VAL, pi=2*acos(0.0d0)
    COS_VAL=cos(1.5*theta+0.25*pi)
    M=M_tc*(1+0.25*COS_VAL**1.2)
    end subroutine Get_M
    
    
    subroutine YieldFunction(q, p, eta_y, F)
	!*********************************************************************
	! Returns the value of the yield function evaluated at q, p , eta    *
	!																	 *
	!*********************************************************************
    implicit none
    !in
    double precision, intent(in):: q, p, eta_y
    !out
    double precision, intent(out):: F
    !local variables
    
    F=q+eta_y*p !sign is due to compression being negative in UMAT
	end subroutine YieldFunction
!***********************************************************************************************
	
!_______________________________________________________________________________________________
!##     ##    ###    ######## ##     ## 
!###   ###   ## ##      ##    ##     ## 
!#### ####  ##   ##     ##    ##     ## 
!## ### ## ##     ##    ##    ######### 
!##     ## #########    ##    ##     ## 
!##     ## ##     ##    ##    ##     ## 
!##     ## ##     ##    ##    ##     ## 
    
    Subroutine TensorInnerProduct(TensorA, TensorB, N, Re)
	!***********************************************************************
	!
	!     Calculate 2NormTensor = sqrt(A:A)
	!
	! I   Tensor  : (Square or vector of dimension N)
	! I   N     :   Number of elements
	! O   2Norm : Resulting norm
	!
	!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension TensorA(N), TensorB(N)
	!***********************************************************************
	  X=N/2
      Re=0.0d0
	  Do I=1,X
		  Re=Re+TensorA(I)*TensorB(I)
	  end Do
	  Do I=X+1,N
		  Re=Re+2*(TensorA(I)*TensorB(I))
	  end do
    end subroutine TensorInnerProduct
    
    
    Subroutine TwoNormTensor(Tensor, N, TwoNorm)
	!***********************************************************************
	!
	!     Calculate 2NormTensor = sqrt(A:A)
	!
	! I   Tensor  : (Square or vector of dimension N)
	! I   N     :   Number of elements
	! O   2Norm : Resulting norm
	!
	!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension Tensor(N)
	!***********************************************************************
	  X=N/2
      TwoNorm=0.0d0
	  Do I=1,X
		  TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
	  end Do
	  Do I=X+1,N
		  TwoNorm=TwoNorm+2*(Tensor(I)*Tensor(I))
	  end do
	  TwoNorm=sqrt(TwoNorm)

    end subroutine TwoNormTensor
        
    Subroutine TwoNormTensor_strain(Tensor, N, TwoNorm)
	!***********************************************************************
	!
	!     Calculate 2NormTensor = sqrt(A:A)
	!
	! I   Tensor  : (Square or vector of dimension N)
	! I   N     :   Number of elements
	! O   2Norm : Resulting norm
	!
	!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension Tensor(N)
	!***********************************************************************
	  X=N/2
      TwoNorm=0.0d0
	  Do I=1,X
		  TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
	  end Do
	  Do I=X+1,N
		  TwoNorm=TwoNorm+0.5*(Tensor(I)*Tensor(I))!The convention in UMAT is to use engineering shear strains
	  end do
	  TwoNorm=sqrt(TwoNorm)

	end subroutine TwoNormTensor_strain
    
    
	Subroutine MatVec(xMat,IM,Vec,N,VecR)
	!***********************************************************************
	!
	!     Calculate VecR = xMat*Vec
	!
	! I   xMat  : (Square) Matrix (IM,*)
	! I   Vec   : Vector
	! I   N     : Number of rows/colums
	! O   VecR  : Resulting vector
	!
	!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
	!***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
    End Subroutine MatVec
    

	 Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
	!***********************************************************************
	!
	!     Calculate VecR() = R1*Vec1()+R2*Vec2()
	!
	! I   Vec1,
	! I   Vec2  : Vectors
	! I   R1,R2 : Multipliers
	! I   N     : Number of rows
	! O   VecR  : Resulting vector
	!
	!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
	!***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
    End Subroutine AddVec
    
    Subroutine DotProduct_2(VecA, VecB,N, Dp)
	!***********************************************************************
	!
	!     Calculate the dot product of A(Nx1) and B(1xN)
	!
	! I   VecA VecB  : Vectors
	! I   N     :   Dimension
	! O   Dp : Dot product
	!
	!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension VecA(N), VecB(N)
	!***********************************************************************
	  Dp=0.0d0
	  Do I=1,N
			  Dp=Dp+VecA(I)*VecB(I)
	  end do

	end subroutine DotProduct_2

	subroutine dbltobool(A,B)
	!******************************************************************
	! Takes a double which values are either 1.0 or 0.0 and returns a *
	! Boolean
	!******************************************************************
	implicit none
	double precision, intent(in):: A
	logical, intent(out):: B
	if (A<1.0) then
		B=.false.
	else
		B=.true.
	endif
	end subroutine dbltobool

	function logic2dbl(a)
	  logical, intent(in) :: a

	  if (a) then
		logic2dbl = 1.d0
	  else
		logic2dbl = 0.d0
	  end if
	end function logic2dbl
	
	Subroutine check4crossing(IErate0I, IErateI, dErate_eff,RateRef, Apply)
    !******************************************************************
	! determines if strain rate updating must occur                   *
	! Boolean                                                         *
	!******************************************************************
		implicit none
		double precision, intent(inout):: IErate0I, IErateI, dErate_eff
        double precision, intent(in)   :: RateRef
		logical:: cond1, cond2, cond3
		logical, intent(out)::Apply
		  Apply=.false.
		  if(IErate0I<=RateRef) IErate0I=RateRef
		  if (IErateI<=RateRef) IErateI=RateRef
		  cond1=(IErate0I==RateRef).and.(IErateI>RateRef)
		  cond2=(IErate0I>RateRef).and.(IErateI==RateRef)
		  dErate_eff=IErateI-IErate0I
		  cond3=(IErate0I>RateRef).and.(IErateI>RateRef)
		  if (cond1.or.cond2.or.cond3) Apply=.true.
	end Subroutine check4crossing
end module ModExternalSoilModel