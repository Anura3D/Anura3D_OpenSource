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
    !	Copyright (C) 2022  Members of the Anura3D MPM Research Community 
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
!     $Date: 2021-02-20 10:27:53 +0100 (sab, 20 feb 2021) $
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


    pointer (p, ESM)             
          
    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
    NameModel = MatParams(IDset)%MaterialModel ! name of constitutive model as specified in GOM-file
    ntens = MatParams(IDset)%UMATDimension  ! 2D or 3D formulation of the External soil model   
          
    ! get strain increments in integration/material point
    TempStrainIncr = GetEpsStep(Particles(IDpt)) ! incremental strain vector assigned to point
    
    if (CalParams%ApplyImplicitQuasiStatic) then
        if (CalParams%ImplicitIntegration%Iteration > 1) then
            do I = 1, NTENSOR
                TempStrainIncrPrevious(I) = GetEpsStepPreviousI(Particles(IDpt), I)
            end do
            
            TempStrainIncr = TempStrainIncr - TempStrainIncrPrevious
            
        end if
    end if
        
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
    if ((CalParams%IStep == 1).and.(CalParams%TimeStep == 1)) then
    StateVar = MatParams(IDset)%ESM_Statvar_in
    !StateVar = 0.0
    !StateVar(1) = 0.0                        !back stress, orientation of yield surface cone(11)
    !StateVar(2) = 0.0                        !back stress, orientation of yield surface cone(22)
    !StateVar(3) = 0.0                        !back stress, orientation of yield surface cone(33)
    !StateVar(4) = 0.0                        !back stress, orientation of yield surface cone(12)
    !StateVar(5) = 0.0                        !back stress, orientation of yield surface cone(13)
    !StateVar(6) = 0.0                        !back stress, orientation of yield surface cone(23)
    !StateVar(7) = 0.952                      !void ratio
    !StateVar(8) = 0.0                        !fabric tensor(11) -> non-zero corresponding to a=1/3 for statistically isotropic particle's orientation
    !StateVar(9) = 0.0                        !fabric tensor(22) -> non-zero corresponding to a=1/3 for statistically isotropic particle's orientation
    !StateVar(10) = 0.0                        !fabric tensor(33) -> non-zero corresponding to a=1/3 for statistically isotropic particle's orientation
    !StateVar(11) = 0.0                        !fabric tensor(12) -> zero
    !StateVar(12) = 0.0                        !fabric tensor(13) -> zero
    !StateVar(13) = 0.0                        !fabric tensor(23) -> zero
    !StateVar(14) = 0.0                        !not used
    !StateVar(15) = 0.0                        !alpha value at stress reversal points (discrete update) (11)
    !StateVar(16) = 0.0                        !alpha value at stress reversal points (discrete update) (22)
    !StateVar(17) = 0.0                        !alpha value at stress reversal points (discrete update) (33)
    !StateVar(18) = 0.0                        !alpha value at stress reversal points (discrete update) (12)
    !StateVar(19) = 0.0                        !alpha value at stress reversal points (discrete update) (13)
    !StateVar(20) = 0.0                        !alpha value at stress reversal points (discrete update) (23)
    !StateVar(21) = 0.0                        !not used
    !StateVar(22) = 0.0                        !not used
    !StateVar(23) = 0.0                        !not used
    !StateVar(24) = 0.0                        !not used
    !StateVar(25) = 0.0                        !not used
    !StateVar(26) = 0.0                        !not used
    !StateVar(27) = 0.0                        !not used
    !StateVar(28) = 0.0                        !not used
    !StateVar(29) = 0.0                        !excess pore pressure (undrained case)
    !StateVar(30) = -100.0                     !mean effective stress
    !StateVar(31) = 0.0                        !deviator stress
    !StateVar(32) = 0.0                        !Lode parameter (cos(3theta))
    !StateVar(33) = 0.0                        !suggested size of first time substep
    !StateVar(34) = 0.0                        !number of function evaluation
    !StateVar(35) = 0.0                        !not used
    !StateVar(36) = 0.0                        !not used
    else 
    StateVar = ESMstatevArray(IDpt,:)
    end if 
          
    
    
    !call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP) !Note that the subroutine checks Cavitation Threshold & Gas Pressure
    ! commented to keep water pressure as hydrostatic 
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
    props = MatParams(IDSet)%ESM_Solid
         
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
    p = GetProcAddress(MatParams(IDSet)%SoilModelDLLHandle, "ESM"C) ! Pointing to the ESM .dll 
    
    if (NameModel == ESM_ARB_Model_SaniSand) then 
           
        call ESM_Sanisand(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
 
    elseif (NameModel == ESM_ARB_Model_MohrCoulomb) then 
    
        call ESM_MohrCoulomb(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
        
    elseif (NameModel == ESM_ARB_Model_MohrCoulombStrainSoftening) then 
   
        call ESM_MohrCoulombStrainSoftening(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
        
    elseif (NameModel == ESM_ARB_Model_Hypoplasticity) then 
        
        call ESM_Hypoplasticity(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
        
    else
        
        call ESM(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
    
    end if 
        
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
    enddo             
                               
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

          MaterialIndex = MaterialIDArray(ParticleID)
          
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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        !░██████╗░█████╗░███╗░░██╗██╗░██████╗░█████╗░███╗░░██╗██████╗░
        !██╔════╝██╔══██╗████╗░██║██║██╔════╝██╔══██╗████╗░██║██╔══██╗
        !╚█████╗░███████║██╔██╗██║██║╚█████╗░███████║██╔██╗██║██║░░██║
        !░╚═══██╗██╔══██║██║╚████║██║░╚═══██╗██╔══██║██║╚████║██║░░██║
        !██████╔╝██║░░██║██║░╚███║██║██████╔╝██║░░██║██║░╚███║██████╔╝
        !╚═════╝░╚═╝░░╚═╝╚═╝░░╚══╝╚═╝╚═════╝░╚═╝░░╚═╝╚═╝░░╚══╝╚═════╝░
        
              Subroutine ESM_SaniSand(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER, &
     DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UDSM" :: UDSM
      implicit double precision (a-h, o-z) 
      CHARACTER*80 CMNAME     
      DIMENSION STRESS(NTENS), &
     DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
      !NPT(1),NOEL(1),IDSET(1),EUNLOADING(1),,PLASTICMULTIPLIER(1),NUMBEROFPHASES(1)

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
              
  
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

!---Local variables defned by the user
! e.g. integer :: var_local	  
!---User can define here additional variables     
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
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
          call umat_SaniSand(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)

      
!---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
!---Always define this value to run the simulation
    
! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
    
          


        return

     end subroutine ESM_SaniSand

        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
!c------------------------------------------------------------------------------
!c complete file suite for Dafalias & Manzari (2004) SANISAND model for sand
!c------------------------------------------------------------------------------
      subroutine umat_SaniSand(stress,statev,ddsdde,sse,spd,scd,&
       rpl,ddsddt,drplde,drpldt,&
       stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
       ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
       celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!c------------------------------------------------------------------------------
!c user subroutine for Abaqus 6.3
!c------------------------------------------------------------------------------
!c
!c Implemented constitutive law:
!c -----------------------------
!c Dafalias & Manzari(2004) SANISAND model for sand
!c
!c Ref: Dafalias, Y. F. and Manzari, M. T.  
!c      Simple plasticity sand model accounting for fabric change effects
!c      J. Engng. Mechanics, ASCE (2004), 130(6):622-634.
!c
!c ----------------------------------------------------------------------------
!c The string for the material name may contain 9 characters.
!c ----------------------------------------------------------------------------
!c Material constants:
!c   	
!c     ---------------------------------------------------------------------
!c     props(j)      
!c     ---------------------------------------------------------------------
!c        1     p_a        Atmospheric pressure 
!c        2     e0         Void ratio on CSL at p = 0  
!c        3     lambda     CSL parameter (e:p plane)
!c        4     xi         CSL parameter (e:p plane)
!c        5     M_c        Slope of CSL in q:p plane, TX compression
!c        6     M_e        Slope of CSL in q:p plane, TX extension
!c        7     mm         opening of yield surface cone
!c        8     G0         Shear modulus constant
!c        9     nu         Poisson's ratio
!c        10    h0         Plastic modulus constant
!c        11    c_h        Plastic modulus constant
!c        12    n_b        Plastic modulus constant
!c        13    A0         Dilatancy constant
!c        14    n_d        Dilatancy constant
!c        15    z_max      Fabric index constant
!c        16    c_z        Fabric index constant
!c        17    bulk_w     Pore water bulk modulus (undrained conditions)
!c	      18    p_tmult    shift of mean stress pt=p_tmult*p_a	
!c        19    initial value of void ratio
!c     ----------------------------------------------------------------------
!c
!c     Solution dependent state variables (statev):
!c     definition via sdvini
!c
!c     group 1: internal variables (14 variables)
!c
!c        1 ... alpha_11	  back stress, orientation of yield surface cone
!c        2 ... alpha_22
!c        3 ... alpha_33
!c        4 ... alpha_12
!c        5 ... alpha_13
!c        6 ... alpha_23
!c
!c        7 ... void       void ratio
!c
!c        8 ... Fab_11     fabric tensor z
!c        9 ... Fab_22
!c       10 ... Fab_33
!c       11 ... Fab_12
!c       12 ... Fab_13
!c       13 ... Fab_23
!c
!c       14 ... not used        
!c
!c     group 2: memory variables for shear reversal (SR) and other purposes
!c
!c       15 ... alpha_sr_11	 alpha value at stress reversal points (discrete update)
!c       16 ... alpha_sr_22   
!c       17 ... alpha_sr_33   
!c       18 ... alpha_sr_12
!c       19 ... alpha_sr_13
!c       20 ... alpha_sr_23
!c
!c       21 ... not used
!c       22 ... not used
!c       23 ... not used
!c       24 ... not used
!c       25 ... not used
!c       26 ... not used
!c       27 ... not used
!c
!c       28 ... not used 
!c
!c     group 3: variables saved for post processing or other purposes
!c
!c       29 ... pore	    excess pore pressure (undrained case)
!c       30 ... p	    mean effective stress
!c       31 ... q	    deviator stress
!c       32 ... z	    Lode parameter (cos(3theta))
!c       33 ... dtsub	suggested size of first time substep
!c       34 ... nfev	    number of function evaluation
!c       35 ... not used
!c       36 ... not used
!c
!c Authors: 
!c     C. Tamagnini (tamag@unipg.it)
!c     Dipartimento di Ingegneria Civile e Ambientale 
!c     Università degli Studi di Perugia, Italy
!c
!c     M. Martinelli
!c     Dipartimento di Ingegneria Strutturale e Geotecnica
!c     Università di Roma "La Sapienza", Italy
!c
!c     C. Miriano
!c     Dipartimento di Ingegneria Strutturale e Geotecnica
!c     Università di Roma "La Sapienza", Italy
!c
!c Modifications: D. Masin, 2015
!c
!c NOTES: 
!c     - sign convention for stress and strain: tension and extension positive
!c     - stress and strain sign convention changed upon entering the SP algorithm
!c     - tangent stiffness operator evaluated according to two alternative options
!c       selected setting the logical flag "cons_lin":
!c       cons_lin.eq.0  -> numerical linearization via
!c                            direct perturbation of dstran
!c       cons_lin.eq.1 -> continuum tangent stiffness
!c                            (not optimal for full N-R iterative solver)
!c       cons_lin.eq.2 -> elastic tangent stiffness
!c                            (even less optimal for full N-R iterative solver, but sometimes more stable)
!c
!c
!c Last change: 4/2013  
!c
!c----------------------------------------------------------------------------
!c
      implicit none
!c
      character*80 cmname
!c
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,&
      layer, kspt, kstep, kinc
!c
      double precision stress(ntens), statev(nstatv),&
       ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),&
       stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),&
       props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      double precision sse, spd, scd, rpl, drpldt, dtime, temp,& 
       dtemp, pnewdt, celent
!c
!c ... 1. nasvdim    = maximum number of additional state variables
!c     2. tolintT    = prescribed error tolerance for the adaptive 
!c                     substepping scheme
!c     3. maxnint    = maximum number of time substeps allowed.
!c                     If the limit is exceeded abaqus is forced to reduce 
!c                     the overall time step size (cut-back) 
!c     4. DTmin      = minimum substeps size allowed.
!c                     If the limit is exceeded abaqus is forced to reduce 
!c                     the overall time step size (cut-back)
!c     5. perturb    = perturbation par. for computation of Jacobian matrices
!c     6. nfasv      = number of first additional state variable in statev field 
!c     7. prsw       = switch for printing information
!c
!c ... declaration of local variables
!c
      integer prsw,elprsw,cons_lin,abaqus,chiara,check_ff,drcor
!c
      integer i,error,maxnint,nfev,mario_DT_test,inittension
      integer nparms,nasvdim,nfasv,nydim,nzdim
      integer nasvy,nasvz,nyact,nzact,plastic,testing
!c
      !double precision dot_vect
!c	
      double precision parms(nprops),theta,tolintT,dtsub,DTmin,perturb
      double precision sig_n(6),sig_np1(6),DDtan(6,6),pore,PI
      double precision deps_np1(6),depsv_np1,norm_D2,norm_D,tolintTtest
      double precision eps_n(6),epsv_n,alphayield(6)
      double precision norm_deps2,norm_deps,pp,qq,cos3t,ddum
      double precision zero,tol_f,fact_thres,p_thres,stran_lim,eps_debug
      double precision p_atm,ptshift,phimob,tol_f_test,youngel,nuel
      double precision avoid,apsi,aec,fyield,Mb!yf_DM,psi_void_DM,
      double precision dummy,sdev(6),I1,alpha(6),cM,tau(6),gth,etanorm
      double precision sinphinorm
!c
      parameter (nasvdim  = 36)
      parameter (nydim    = 6+14)
      parameter (nzdim    = 14)
      parameter (tolintT  = 1.00d-3)
      parameter (tolintTtest = 1.0d-2) 
!c
      parameter (maxnint  = 50000)
      parameter (DTmin    = 1.0d-18)
      parameter (perturb  = 1.0d-4)
      parameter (nfasv    = 1)
      parameter (prsw     = 0)
      parameter (cons_lin = 1)
	parameter (abaqus = 0)
!c chiara
 	parameter (eps_debug = 0.9d-3)
!c
      parameter (zero = 0.0d0)
      parameter (PI = 3.14159265358979323846264338327950288)
      parameter (fact_thres=0.000000001d0)
!c
!c ... additional state variables
!c
      double precision  asv1(nydim-6),asv2(nzdim)
!c
!c ... solution vector (stresses, additional state variables)
!c
      double precision  y(nydim),y_n(nydim),z(nzdim),z_n(nzdim)
!c
!c      common /z_nct_errcode/error
!c	common /z_tolerance/tol_f
!c	common /z_check_yeld/check_ff
!c	common /z_drift_correction/drcor
!c	common /z_threshold_pressure/p_thres
!	
!c
        tol_f=1.0d-6
        tol_f_test=1.0d-6
	check_ff=0
	drcor=1
	plastic=0
	phimob=0.0d0
	ptshift=0.0d0
!c
!c
!c ... Error Management:
!c     ----------------
!c     error =  0 ... no problem in time integration
!c     error =  1 ... problems in evaluation of the time rate, (e.g. undefined 
!c                    stress state), reduce time integration substeps
!c     error =  3 ... problems in time integration, reduce abaqus load increment 
!c                    (cut-back)
!c     error=10 ... severe error, terminate calculation
!c
      !error=0
    ndi = 3
!c
!c ... check problem dimensions
!c
      if (ndi.ne.3) then
!c
        write(6,*) 'ERROR: this UMAT can be used only for elements'
        write(6,*) '       with 3 direct stress/strain components'
        write(6,*) 'noel = ',noel
        error=10
!c
      endif
!c      open(unit=6,position='Append',file=
!c     .	'C:/users/david/data/Zhejiang/plaxis-sani/UMATdebug.txt')
!c
!c ... check material parameters and move them to array parms
!c
      nparms=nprops
      call check_parms_DM(props,parms,nparms)
!c
!c ... print informations about time integration, useful when problems occur
!c
	p_atm=parms(1)
	p_thres=fact_thres*p_atm
!c
      elprsw = 0
      if (prsw .ne. 0) then
!c
!c ... print only in some defined elements
!c
!c	 if ((noel.eq.101).and.(npt.eq.1)) elprsw = 1
!c
      endif
!c
!c ... define number of additional state variables
!c
      call define(nasvy,nasvz)
      nyact = 6 + nasvy
      nzact = nasvz
      if (nyact.gt.nydim) then
        write(6,*) 'ERROR:nydim too small; program terminated'
        error=10
      elseif (nzact.gt.nzdim) then
        write(6,*) 'ERROR:nzdim too small; program terminated'
        error=10
      endif
!c
!c ... suggested time substep size, and initial excess pore pressure
!c
      pore = statev(29)
!c
!c ... changes sign conventions for stresses and strains, compute 
!c     current effective stress tensor and volumetric strain
!c
      ptshift=parms(18)*parms(1)
      do i=1,3
          stress(i) = stress(i)-ptshift
      enddo

      call move_sig(stress,ntens,-1*ptshift,sig_n)
	
      call move_sig(stress,ntens,pore,sig_n)
      call move_eps(dstran,ntens,deps_np1,depsv_np1)
      call move_eps(stran,ntens,eps_n,epsv_n)
      
      norm_D2=dot_vect(2,deps_np1,deps_np1,6)
      norm_D=sqrt(norm_D2)
!c chiara
	if (eps_n(1).gt.eps_debug) then
		chiara=1
	end if
!c
!c ... initialise void ratio and yield surface inclination
!c
      if (statev(7) .lt. 0.001) then
            do i=1,6        
               alphayield(i)=zero
            end do
            call deviator(sig_n,alphayield,ddum,pp)
      	    avoid=0
            if(parms(19) .le. 5.0) then 
                   avoid=parms(19)
            else if(parms(19) .gt. 5.0) then
        	   apsi=parms(19)-10.0d0
        	   aec=parms(2)-parms(3)*(pp/parms(1))**parms(4)
        	   avoid=aec+apsi
            endif
            statev(7)=avoid
            do i=1,6        
               statev(i)=alphayield(i)/pp
               statev(i+14)=alphayield(i)/pp
            end do
      end if
!c
!c ... move vector of additional state variables into asv1(nydim) and asv2(nzdim)
!c
      do i=1,nasvy
        asv1(i) = statev(i-1+nfasv)
      enddo
!c
      do i=1,nasvz
        asv2(i) = statev(i-1+nfasv+nasvy)
      enddo
!c
!c --------------------
!c ... Time integration
!c --------------------
!c
      call iniyz(y,nydim,z,nzdim,asv1,nasvy,asv2,nasvz,sig_n,ntens)
      call push(y,y_n,nydim)
      call push(z,z_n,nzdim)
!c
      if (elprsw.ne.0) then
        write(6,*) '================================================='
        write(6,*) '         Call of UMAT - DM SANISAND model:       '
        write(6,*) '================================================='
        call wrista(3,y,nydim,deps_np1,dtime,coords,statev,nstatv,&
                   parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
      endif
!c
!c ... local integration using adaptive RKF-23 method with error control
!c
      dtsub = 1e-5
      if((dtsub.le.zero).or.(dtsub.gt.dtime)) then
        dtsub = dtime
      end if
      
      testing=0
      kstep = 1
      kinc = 1
!c     For use in PLAXIS, activate the following line
      !if(kstep.eq.1 .AND. kinc.eq.1) testing=1
!c     For use in ABAQUS, the line above should be inactive
	
      if(norm_D.eq.0) testing=2
!c     FEM asking for ddsdde only

      nfev = 0 ! initialisation
      error = 0

      if(testing.eq.1) then
        call rkf23_upd_DM(y,z,nyact,nasvy,nasvz,tolintTtest,maxnint,&
              DTmin,deps_np1,parms,nparms,nfev,elprsw,&
     	       mario_DT_test,&
              error,tol_f_test,check_ff,drcor,p_thres,plastic)
!c ... give original state if the model fails without substepping
          if(error.ne.0) then
            do i=1,nyact        
               y(i)=y_n(i)
            end do
            error=0
          end if
      else if(testing.eq.2) then
            do i=1,nyact        
                  y(i)=y_n(i)
            end do
!c ... Normal RKF23 integration
      else   !testing.eq.0
        call rkf23_upd_DM(y,z,nyact,nasvy,nasvz,tolintT,maxnint,&
              DTmin,deps_np1,parms,nparms,nfev,elprsw,&
     	       mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
      end if
!c
!c ... error conditions (if any)
!c
      mario_DT_test = 0
    !commented by Abdelrahman
	if(mario_DT_test.eq.1) then
	call wrista(4,y,nydim,deps_np1,dtime,coords,statev,nstatv,&
                 parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
    endif
!c
    !error = 0
      if(error.eq.3) then
!c
!c ... reduce abaqus load increment
!c
        call wrista(2,y,nydim,deps_np1,dtime,coords,statev,nstatv,&
                parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
        write(6,*) 'subroutine UMAT: reduce step size in ABAQUS'
        write(6,*) 'error 3 activated'

        if(abaqus.ne.0) then
            pnewdt = 0.25d0
        else
!c ... write a message and return the original state
           do i=1,nyact        
                  y(i)=y_n(i)
           end do     
        endif
!c			
!c        if(abaqus.eq.0) then
!c                write(6,*) 'analysis ended because number of time ' 
!c		  write(6,*) 'substeps exceeded maximum number allowed'
!c		  write(6,*) ' (maxnint)'
!c	      call xit_DM
!c	  endif
  
	return
!c
      elseif(error.eq.10) then
      	write(6,*) 'error 10 activated'
!c
        call wrista(2,y,nydim,deps_np1,dtime,coords,statev,nstatv,&
                   parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
        call xit_DM
!c
      endif
!c
!c ... update dtsub and nfev
!c
      if(dtsub.le.0.0d0) then 
      	dtsub = 0
      else if(dtsub.ge.dtime) then 
      	dtsub = dtime
      end if
      statev(33)=dtsub
      statev(34)=dfloat(nfev)
!c
!c ... computation of Jacobian matrix
!c
      error=0
      if(cons_lin.eq.0) then
!c
!c ... parameter of the numerical differentiation
!c     double precision
!c
        norm_deps2=dot_vect(2,deps_np1,deps_np1,ntens)
        norm_deps=dsqrt(norm_deps2)
        theta=perturb*max(norm_deps,1.0d-6)
!c
!c ... compute consistent tangent via numerical perturbation
!c
        call pert_DM(y_n,y,z,nyact,nasvy,nasvz,tolintT,maxnint,DTmin,&
            deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DDtan,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
      else
!c        
        call tang_stiff(y,z,nyact,nasvy,nasvz,parms,nparms,&
              DDtan,cons_lin,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
      endif
!c
	if(error.ne.0) then
	 write(6,*) 'some error in pert or tang_stiff, value=',error
!c
!c        call wrista(2,y,nydim,deps_np1,dtime,coords,statev,nstatv,&
!c                   parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
!c        call xit_DM
!c
      endif
!c
!c ... convert solution (stress + cons. tangent) to abaqus format
!c     update pore pressure and compute total stresses

      inittension=0
!c     just checks for NAN
      call check_RKF_DM(inittension,y,nyact,nasvy,parms,nparms)
      if (inittension.ne.0) then
           do i=1,nyact        
                  y(i)=y_n(i)
           end do        
      end if
	
!c
      call solout(stress,ntens,asv1,nasvy,asv2,nasvz,ddsdde,&
                 y,nydim,z,pore,depsv_np1,parms,nparms,DDtan)
!c
!c ... updated vector of additional state variables to abaqus statev vector
!c
      do i=1,nasvy
        statev(i-1+nfasv) = asv1(i) 
      end do
!c
      do i=1,nasvz
        statev(i-1+nfasv+nasvy) = asv2(i)
      enddo
!c
!c ... transfer additional information to statev vector
!c
      do i=1,6
        sig_np1(i)=y(i)
      end do
      call inv_sig(sig_np1,pp,qq,cos3t)
!c
      statev(29) = pore 
      statev(30) = pp
      statev(31) = qq
      statev(32) = cos3t
     
      cM=parms(6)/parms(5)
      alpha(1)=y(7)
      alpha(2)=y(8) 
      alpha(3)=y(9) 
      alpha(4)=y(10) 
      alpha(5)=y(11) 
      alpha(6)=y(12)
      call deviator(sig_np1,sdev,I1,pp)
      do i=1,6
        tau(i)=sdev(i)-pp*alpha(i)
      end do
      call lode_DM(tau,cM,cos3t,gth,dummy)
      if (pp<0.001) then 
          etanorm = 0
      else 
          etanorm=gth*qq/pp
      end if
      sinphinorm=3*etanorm/(6+etanorm)
      !if (sinphinorm<-1) then 
      !    sinphinorm=-1
      !end if 
      
      statev(33) = asin(sinphinorm)*180/PI
      statev(34) = nfev
      
      
      
!c      check that bounding surtface is not violated 
!c      if (noel.eq.324 .and. npt.eq.3) then
!c      	    fyield=yf_DM(y,nyact,parms,nparms)  
!c      	    apsi=psi_void_DM(statev(7),pp,parms,nparms)
!c      	    Mb=parms(5)*dexp(-parms(12)*apsi)
!c      	    write(6,*) 'fyield=',fyield
!c      	    write(6,*) 'psi=',apsi
!c      	    write(6,*) 'Mb=',Mb
!c      	    write(6,*) 'qq/pp*gth=',qq/pp*gth
!c      	    write(6,*) '---------------------'
!c      end if
      
      do i=1,3
          stress(i) = stress(i)+ptshift
      enddo
!c      close(6)
!c
!c -----------------------
!c End of time integration
!c -----------------------
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine alpha_th_DM(flag,n,gth,psi,parms,nparms,alpha)
!c-----------------------------------------------------------------------------
!c	calculate tensors:
!c
!c     alpha_c => flag = 1
!c     alpha_b => flag = 2
!c     alpha_d => flag = 3
!c
!c   Dafalias & Manzari (2004) SANISAND model for sands
!c
!c   variables allocated in parms
!c
!c   1     p_a        Atmospheric pressure 
!c   2     e0         Void ratio on CSL at p = 0  
!c   3     lambda     CSL parameter (e:p plane)
!c   4     xi         CSL parameter (e:p plane)
!c   5     M_c        Slope of CSL in q:p plane, TX compression
!c   6     M_e        Slope of CSL in q:p plane, TX extension
!c   7     mm         opening of yield surface cone
!c   8     G0         Shear modulus constant
!c   9     nu         Poisson's ratio
!c   10    h0         Plastic modulus constant
!c   11    c_h        Plastic modulus constant
!c   12    n_b        Plastic modulus constant
!c   13    A0         Dilatancy constant
!c   14    n_d        Dilatancy constant
!c   15    z_max      Fabric index constant
!c   16    c_z        Fabric index constant
!c   17    bulk_w     Pore water bulk modulus (undrained conditions)
!c
!c	written 10/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer flag,nparms,i
!c
      double precision n(6),gth,psi,parms(nparms),alpha(6)
      double precision M_c,mm,n_b,n_d
!c
      double precision M,alpha_th
      double precision two,three,sqrt23
!c
      data two,three/2.0d0,3.0d0/
!c
      sqrt23=dsqrt(two/three)
!c
!c ... recover material parameters
!c
      M_c=parms(5)
	mm=parms(7)
	n_b=parms(12)
	n_d=parms(14)
!c
!c ... select which alpha tensor to evaluate
!c
      if(flag.eq.1) then
!c
!c ... critical state cone
!c
        M=M_c
!c
      elseif(flag.eq.2) then 
!c
!c ... bounding surface cone
!c
        M=M_c*dexp(-n_b*psi)
!c
      else
!c
!c ... dilatancy cone
!c
        M=M_c*dexp(n_d*psi)
!c
      endif
!c
!c ... tensor alpha_ij
!c
      alpha_th=M*gth-mm 
!c
      do i=1,6
        alpha(i)=sqrt23*alpha_th*n(i)      
      end do
!c	
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine check_crossing(y,y_tr,n,parms,nparms,prod)
!c-----------------------------------------------------------------------------
!c
!c  computes 
!c
!c  prod := dsig_tr*grad(f)_k 
!c
!c  useful for checking if crossing of yield locus occurs whenever 
!c  f_k = 0 and f_tr > 0 
!c
!c  variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer i,n,nparms
!c
      !double precision dot_vect
!c
      double precision y(n),y_tr(n),parms(nparms)
      double precision P(6),P1(6),dsig_tr(6)
      double precision prod
!c
!c ... gradient of yield surface at state y_k
!c
      call grad_f_DM(y,n,parms,nparms,P,P1)
!c
      do i=1,6
        dsig_tr(i)=y_tr(i)-y(i)
      end do ! i
!c		  
      prod=dot_vect(1,P,dsig_tr,6)  
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine check_parms_DM(props,parms,nprops)
!c-----------------------------------------------------------------------------
!c     checks input material parameters for Dafalias & Manzari (2004)  
!c     SANISAND model for sand
!c
!c     Material constants:
!c   	
!c     ---------------------------------------------------------------------
!c     props(j)      
!c     ---------------------------------------------------------------------
!c        1     p_a        Atmospheric pressure 
!c        2     e0         Void ratio on CSL at p = 0  
!c        3     lambda     CSL parameter (e:p plane)
!c        4     xi         CSL parameter (e:p plane)
!c        5     M_c        Slope of CSL in q:p plane, TX compression
!c	   6     M_e        Slope of CSL in q:p plane, TX extension
!c        7     mm         opening of yield surface cone
!c        8     G0         Shear modulus constant
!c        9     nu         Poisson's ratio
!c        10    h0         Plastic modulus constant
!c        11    c_h        Plastic modulus constant
!c        12    n_b        Plastic modulus constant
!c        13    A0         Dilatancy constant
!c        14    n_d        Dilatancy constant
!c	   15    z_max      Fabric index constant
!c        16    c_z        Fabric index constant
!c        17    bulk_w     Pore water bulk modulus (undrained conditions)
!c     ---------------------------------------------------------------------
!c
!c     Solution dependent state variables (statev):
!c     definition via sdvini
!c
!c     group 1: internal variables (14 variables)
!c
!c        1 ... alpha_11	  back stress, orientation of yield surface cone
!c        2 ... alpha_22
!c        3 ... alpha_33
!c        4 ... alpha_12
!c        5 ... alpha_13
!c        6 ... alpha_23
!c
!c        7 ... void       void ratio
!c
!c        8 ... Fab_11     fabric tensor z
!c        9 ... Fab_22
!c       10 ... Fab_33
!c       11 ... Fab_12
!c       12 ... Fab_13
!c       13 ... Fab_23
!c
!c       14 ... not used        
!c
!c     group 2: memory variables for shear reversal (SR) and other purposes
!c
!c       15 ... alpha_sr_11	 alpha value at stress reversal points (discrete update)
!c       16 ... alpha_sr_22   
!c       17 ... alpha_sr_33   
!c       18 ... alpha_sr_12
!c       19 ... alpha_sr_13
!c       20 ... alpha_sr_23
!c
!c       21 ... not used
!c       22 ... not used
!c       23 ... not used
!c       24 ... not used
!c       25 ... not used
!c       26 ... not used
!c       27 ... not used
!c
!c       28 ... not used 
!c
!c     group 3: variables saved for post processing or other purposes
!c
!c       29 ... pore	    excess pore pressure (undrained case)
!c       30 ... p	    mean effective stress
!c       31 ... q	    deviator stress
!c       32 ... z	    Lode parameter (cos(3theta))
!c       33 ... dtsub	suggested size of first time substep
!c       34 ... nfev	    number of function evaluation
!c       35 ... not used
!c       36 ... not used
!c
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer nprops
!c
      double precision props(nprops),parms(nprops)

      double precision p_a,e0,lambda,xi,M_c,M_e,mm
      double precision G0,nu,h0,c_h,n_b,A0,n_d,z_max
      double precision c_z,bulk_w,sinphi,PI,sinphiext
!c
	  double precision zero
!c
      parameter(zero=0.0d0)
      parameter(PI=3.14159265358979323846264338327950288)
!c
!c ... recover material parameters and initial state info
!c
      p_a=props(1)
      e0=props(2)
      lambda=props(3)
	xi=props(4)
      M_c=props(5)
      M_e=props(6)
      mm=props(7)
      G0=props(8)
      nu=props(9) 
      h0=props(10)
      c_h=props(11)
      n_b=props(12)
      A0=props(13)
      n_d=props(14)
      z_max=props(15)
      c_z=props(16)
      bulk_w=props(17)
!c
!c ... move vector props into local vector parms
!c
      call push(props,parms,nprops)

      if(parms(5) .gt. 5) then
      	      sinphi=sin(parms(5)/180*PI)
      	      parms(5)=6*sinphi/(3-sinphi)
      else
      	      sinphi=3*parms(5)/(6+parms(5))
      end if
      if(parms(6) .gt. 5) then
      	      sinphiext=sin(parms(6)/180*PI)
      	      parms(6)=6*sinphiext/(3+sinphiext)
      else if ((parms(6) .le. 5) .and. (parms(6) .gt. 0.01)) then
      	      sinphiext=3*parms(6)/(6-parms(6))
      else
      	      parms(6)=parms(5)*(3-sinphi)/(3+sinphi)
      end if

      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine define(nasvy,nasvz)
!c-----------------------------------------------------------------------------
      implicit none 
      integer nasvy,nasvz
!c
!c number of additional state variables stored in vectors y and z
!c must 14 each (otherwise change nasvdim in umat)
!c
!c        components of ASV(i) stored in y (14 variables)
!c
!c        1 ... alpha_11	  back stress, orientation of yield surface cone
!c        2 ... alpha_22
!c        3 ... alpha_33
!c        4 ... alpha_12
!c        5 ... alpha_13
!c        6 ... alpha_23
!c        7 ... void       void ratio
!c        8 ... Fab_11	  fabric tensor
!c        9 ... Fab_22
!c       10 ... Fab_33
!c       11 ... Fab_12
!c       12 ... Fab_13
!c       13 ... Fab_23
!c       14 ... not used        
!c
!c       components of ASV(i) stored in z (14 variables)
!c
!c       15 ... alpha_sr_11	alpha value at stress reversal points (discrete update)
!c       16 ... alpha_sr_22  
!c       17 ... alpha_sr_33
!c       18 ... alpha_sr_12
!c       19 ... alpha_sr_13
!c       20 ... alpha_sr_23
!c       21 ... not used     
!c       22 ... not used
!c       23 ... not used
!c       24 ... not used
!c       25 ... not used
!c       26 ... not used
!c       27 ... not used
!c       28 ... not used 
!c
      nasvy = 14
      nasvz = 14
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c------------------------------------------------------------------------------
      subroutine deviator(t,s,trace,mean)
!c------------------------------------------------------------------------------
!c calculate deviator and trace of 2nd order tensor t(6)
!c
!c NOTE: Voigt notation is used with the following index conversion
!c
!c       11 -> 1
!c       22 -> 2
!c       33 -> 3
!c       12 -> 4
!c       13 -> 5
!c       23 -> 6
!c
!c------------------------------------------------------------------------------
!c
      implicit none
!c
      double precision t(6),s(6),trace,mean
      double precision one,three,onethird
!c
      data one,three/1.0d0,3.0d0/
!c
!c ... some constants
!c
      onethird=one/three
!c
!c ... trace and mean value
!c
      trace=t(1)+t(2)+t(3)
      mean=onethird*trace
!c
!c ... deviator stress
!c
      s(1)=t(1)-mean
      s(2)=t(2)-mean
      s(3)=t(3)-mean
      s(4)=t(4)
      s(5)=t(5)
      s(6)=t(6)
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c------------------------------------------------------------------------------
      double precision function distance_sanisand(alpha_k,alpha,n)
!c------------------------------------------------------------------------------
!c computes distance function
!c
!c     d = (alpha^k_{ij}-alpha_{ij})n_{ij}   (k=sr,b,d)
!c
!c Dafalias & Manzari (2004) SANISAND model for sands
!c
!c------------------------------------------------------------------------------
      implicit none
!c
      integer i
!c
      !double precision dot_vect
!c
      double precision alpha_k(6),alpha(6),n(6),delta(6)
!c
      do i=1,6
        delta(i)=alpha_k(i)-alpha(i)
      end do
!c
      distance_sanisand=dot_vect(1,delta,n,6)
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c------------------------------------------------------------------------------
      double precision function dot_vect(flag,a,b,n)
!c------------------------------------------------------------------------------
!c dot product of a 2nd order tensor, stored in Voigt notation
!c
!c flag = 1 -> vectors are stresses in Voigt notation
!c flag = 2 -> vectors are strains in Voigt notation
!c flag = 3 -> ordinary dot product between R^n vectors
!c------------------------------------------------------------------------------
      implicit none
      integer i,n,flag
      double precision a(n),b(n)
      double precision zero,half,one,two,coeff
!c
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
!c
      if(flag.eq.1) then
!c
!c ... stress tensor (or the like)
!c
        coeff=two
!c
      elseif(flag.eq.2) then
!c
!c ... strain tensor (or the like)
!c
        coeff=half
!c
      else
!c
!c ... standard vectors
!c
        coeff=one
!c	
      end if
!c
      dot_vect=zero
!c
      do i=1,n
        if(i.le.3) then
          dot_vect = dot_vect+a(i)*b(i)
        else
          dot_vect = dot_vect+coeff*a(i)*b(i)
        end if
      end do
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine drift_corr_DM(y,n,z,nasvz,parms,nparms,tol,switch2,&
      mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c  performs consistent drift correction (see Sloan et al. 2001)
!c  Dafalias & Manzari(2004) SANISAND model for sand
!c
!c  written 8/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
!c
      !double precision dot_vect
!c
	integer switch2,mario_DT_test
!c
      !external matmul
      !double precision yf_DM
!c 
      integer n,nasvz,nparms,i,n_drift,max_ndrift,switch
	integer iter, itermax
	
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres

!c
      double precision y(n),y0(n),y1(n), z(nasvz),parms(nparms)
      double precision gradf(6),gradf1(6),gradg(6),gradg1(6)
      double precision DDe(6,6),UU(6),VV(6),h_alpha(6),Kpm1,p1,pp1
      double precision f0,tol,zero,one,denom,fnm1,p,three,onethird,f0_p
	double precision factor,f1,p_atm
!c
      parameter(zero=0.0d0,one=1.0d0,three=3.0d0)
      parameter(max_ndrift=10000, itermax=1000)
!c
!c      common /z_nct_errcode/error	
!c
!c ... initialize constants and vectors
!c
      call push(y,y0,n)
!c
!c ... check if current state is inside the elastic nucleus
!c
!c Chiara Miriano 15 maggio 2009
!c	onethird=one/three
!c	p=(y0(1)+y0(2)+y0(3))*onethird
!c	if(p.lt.zero) then
!c		do i=1,3
!c			y0(i)=y0(i)-p
!c		end do
!c	p=(y0(1)+y0(2)+y0(3))*onethird
!c	end if
!c Chiara Miriano  15 maggio 2009
!c
      f0=yf_DM(y0,n,parms,nparms)
	onethird=one/three
	p=(y0(1)+y0(2)+y0(3))*onethird
	
      n_drift=0
	switch=0
    if (p<0.001) then 
        f0_p = 0
    else 
        f0_p=f0/p
    end if
    
!c	p_atm=parms(1)
	
!c	if(p.lt.(p_atm/100)) f0_p=f0_p/1000000
!c	f0_p=f0
	switch2=0
!c
      do while(f0_p.gt.tol)
!cccc		do while(f0.gt.tol)
!c
        fnm1=f0
!c
!c ... current state outside yield surface, correct it until f0 < ftol
!c
        n_drift=n_drift+1
!c
!c ... elastic stiffness and gradients of f and g
!c
        call el_stiff_DM(y0,n,parms,nparms,DDe,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
        call grad_f_DM(y0,n,parms,nparms,gradf,gradf1)
        call grad_g_DM(y0,n,parms,nparms,gradg,gradg1)
!c
!c ... vectors UU=DDe*gradg and VV=DDe*gradf		
!c
	  call matmul(DDe,gradg1,UU,6,6,1)
	  call matmul(DDe,gradf1,VV,6,6,1)
!c
!c ... hardening function h_alpha and plastic modulus (1/Kp)
!c
       call plast_mod_DM(y0,n,z,nasvz,parms,nparms,h_alpha,Kpm1,&
      switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
	if (switch2.gt.zero) return
!c
        if(one/Kpm1.le.zero) then 
!c          write(6,*) 'ERROR: subroutine DRIFT_CORR:'
          write(6,*) 'subcritical softening condition'
!c          write(6,*) 'Kp = ',one/Kpm1
!c          error=10
	   error=3
!c		switch2=1
         return
!cc		switch=1
        end if
!c
!c ... correction for stress (y(1):y(6))
!c	
	  if(switch.eq.0) then
		do i=1,6
			y1(i)=y0(i)-Kpm1*f0*UU(i)
		end do
!c
!c ... correction for hardening variable alpha (y(7):y(12))
!c
		do i=1,6
			y1(6+i)=y0(6+i)+Kpm1*f0*h_alpha(i)
		end do
		do i=13,n
			y1(i)=y0(i)
		end do	
!c
!c ... recompute drift at the new state
!c
		f0=yf_DM(y1,n,parms,nparms)
		if(f0.gt.fnm1) then
			switch=1
!c	switch2=1
	p1=(y1(1)+y1(2)+y1(3))*onethird
!c	write(*,*)'switch2=',switch2
!c	write(*,*)'p=',p
!c	write(*,*)'p1=',p1
!c	return
		else
			call push(y1,y0,n)
		end if		
	  else
!c
!c ... normal correction in place of consistent correction
!c
!c		call push(y,y0,n)
		call push(y0,y1,n)
		f0=yf_DM(y0,n,parms,nparms)
!c		f1=yf_DM(y1,n,parms,nparms)
		denom=dot_vect(1,gradf,gradf,6)
	factor=one
	f1=f0
!ccc	iter=0
!ccc	do while(f1.ge.f0)
!ccc		iter=iter+1
!ccc		if (iter.gt.itermax) then
!ccc		error=10
!ccc		endif
		do i=1,6
			y1(i)=y0(i)-f0*gradf(i)/denom/factor
		end do
		do i=13,n
			y1(i)=y0(i)
		end do
		pp1=(y1(1)+y1(2)+y1(3))*onethird
!ccc		factor=factor*2
			if(pp1.lt.zero)then
!c			write(*,*)'pp1<0=',pp1
				switch2=1
				return
			endif

		f1=yf_DM(y1,n,parms,nparms)
!ccc	enddo


!c		write(*,*) 'drift_corr: normal correction'
		call push(y1,y0,n)

!c		pause
	  end if
!c
!c ... recompute drift at the new state
!c
	  f0=yf_DM(y0,n,parms,nparms)
	  p=(y0(1)+y0(2)+y0(3))*onethird
	  f0_p=f0/p
!c	  if(p.lt.(p_atm/100)) f0_p=f0_p/1000000
!c	  f0_p=f0	  
!c
        if(n_drift.gt.max_ndrift) then
          write(6,*) 'ERROR: subroutine DRIFT_CORR:'
          write(6,*) 'too many iterations, increase tolerance'
          write(6,*) 'n_drift = ',n_drift
          write(6,*) 'drift = ',f0_p
!c          error=10
!c	   error=3
!c          return
	  f0_p=0
        end if		

	
!c ... bottom of while loop
!c
!c	  write(*,*) 'drift_corr:      f0 = ',f0
!c	  write(*,*) 'drift_corr: n_drift = ',n_drift
	end do
!c
!c ... return corrected stress and q into vector y
!c
      call push(y0,y,n)      
!c	
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine el_stiff_DM(y,n,parms,nparms,DDe,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c------------------------------------------------------------------------------
!c subroutine to compute elastic stiffness 
!c Dafalias &Manzari SANISAND model (2004)
!c
!c material parameters
!c
!c  1     p_a        Atmospheric pressure 
!c  2     e0         Void ratio on CSL at p = 0  
!c  3     lambda     CSL parameter (e:p plane)
!c  4     xi         CSL parameter (e:p plane)
!c  5     M_c        Slope of CSL in q:p plane, TX compression
!c  6     M_e        Slope of CSL in q:p plane, TX extension
!c  7     mm         opening of yield surface cone
!c  8     G0         Shear modulus constant
!c  9     nu         Poisson's ratio
!c  10    h0         Plastic modulus constant
!c  11    c_h        Plastic modulus constant
!c  12    n_b        Plastic modulus constant
!c  13    A0         Dilatancy constant
!c  14    n_d        Dilatancy constant
!c  15    z_max      Fabric index constant
!c  16    c_z        Fabric index constant
!c  17    bulk_w     Pore water bulk modulus (undrained conditions)
!c
!c variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c  NOTE: soil mechanics convention (compression positive)
!c        all stress and strain vectors are 6-dimensional
!c------------------------------------------------------------------------------
      implicit none
!c
      integer i,j,n,nparms
!c
      double precision y(n),parms(nparms)
      double precision p_a,G0,nu,ratio
      double precision sig1,sig2,sig3,p,void
      double precision coeff1,coeff2
      double precision Kt,Gt,fe
      double precision Id(6,6),IxI(6,6),DDe(6,6)
      double precision zero,half,one,two,three 
      double precision pp,p_thres_E,tenm3
!c
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres

      parameter(zero=1.0d0,half=0.5d0)
      parameter(one=1.0d0,two=2.0d0,three=3.0d0)
	parameter(p_thres_E=0.001d0)
!cc      parameter(tenm3=0.001d0)
!c
!c	common /z_threshold_pressure/p_thres
!c
!c ... initialize matrices
!c
      call pzero(Id,36)
      call pzero(IxI,36)
      call pzero(DDe,36)
!c
      Id(1,1)=one
      Id(2,2)=one
      Id(3,3)=one
      Id(4,4)=half
      Id(5,5)=half
      Id(6,6)=half
!c
      IxI(1,1)=one
      IxI(2,1)=one
      IxI(3,1)=one
      IxI(1,2)=one
      IxI(2,2)=one
      IxI(3,2)=one
      IxI(1,3)=one
      IxI(2,3)=one
      IxI(3,3)=one
!c
!c ... recover material parameters
!c
      p_a=parms(1)
      G0=parms(8)
      nu=parms(9) 
!c
!c ... recover state variables
!c
      sig1=y(1)
      sig2=y(2)
      sig3=y(3)
!c
      void=y(13)
!c
!c ... mean stress
!c
      p=(sig1+sig2+sig3)/three
!c
	pp=p
!ccc	if(p.lt.p_thres_E)then
!ccc		pp=p_thres_E
!ccc	end if
	if(p.lt.p_thres)then
		pp=p_thres
	end if
!c
!c ... max. shear modulus, tangent shear modulus Gt, tangent bulk modulus Kt
!c
      ratio=three*(one-two*nu)/(two*(one+nu))
      fe=(2.97d0-void)*(2.97d0-void)/(one+void)
      Gt=G0*p_a*fe*dsqrt(pp/p_a)
      Kt=Gt/ratio
!c
!c ... elastic stiffness, stored in matrix DDe(6,6)
!c
      coeff1=Kt-two*Gt/three
      coeff2=two*Gt
!c
      do i=1,6
        do j=1,6
          DDe(i,j)=coeff1*IxI(i,j)+coeff2*Id(i,j)
        end do
      end do
!c
      return
      end 
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
	  subroutine f_hypoelas_DM(y,n,parms,nparms,deps,F,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c
!c ... computes the function F(y) for (hypo)elastic processes
!c     Dafalias & Manzari SANISAND Model (2004)
!c
!c  variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c  variables allocated in vector z(nasvz):
!c
!c	z(1)	= alpha_sr(1)
!c	z(2)	= alpha_sr(2)
!c	z(3)	= alpha_sr(3)
!c	z(4)	= alpha_sr(4)
!c	z(5)	= alpha_sr(5)
!c	z(6)	= alpha_sr(6)
!c	z(7)	= not used
!c	z(8)	= not used
!c	z(9)	= not used
!c	z(10)	= not used
!c	z(11)	= not used
!c	z(12)	= not used
!c	z(13)	= not used
!c	z(14)	= not used
!c
!c  written by Tamagnini 10/2008
!c
!c-----------------------------------------------------------------------------
!c
      implicit none
!c
      !external matmul
!c
      integer n,m,nparms
!c
      double precision y(n),parms(nparms),deps(6)
      double precision depsv,void
      double precision F(n),De(6,6),dsig_e(6)
      double precision one
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
!c
      data one/1.0d0/
!c
      call pzero(F,n)
!c
!c ... void ratio and volum. strain increment
!c
      void = y(13)
      depsv=deps(1)+deps(2)+deps(3)
!c
!c ... elastic stiffness matrix
!c
      call el_stiff_DM(y,n,parms,nparms,De,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
      call matmul(De,deps,dsig_e,6,6,1)
!c
      F(1)=dsig_e(1)
      F(2)=dsig_e(2)
      F(3)=dsig_e(3)
      F(4)=dsig_e(4)
      F(5)=dsig_e(5)
      F(6)=dsig_e(6)
!c
      F(13)=-(one+void)*depsv
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine f_plas_DM(y,n,nasvy,z,nz,parms,nparms,deps,kRK,nfev,&
      switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c calculate coefficient kRK from current state (stored in y and z) and
!c strain increment deps
!c 
!c Dafalias & Manzari (2004) SANISAND model for sands
!c
!c variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c variables allocated in vector z(nasvz):
!c
!c	z(1)	= alpha_sr(1)
!c	z(2)	= alpha_sr(2)
!c	z(3)	= alpha_sr(3)
!c	z(4)	= alpha_sr(4)
!c	z(5)	= alpha_sr(5)
!c	z(6)	= alpha_sr(6)
!c	z(7)	= not used
!c	z(8)	= not used
!c	z(9)	= not used
!c	z(10)	= not used
!c	z(11)	= not used
!c	z(12)	= not used
!c	z(13)	= not used
!c	z(14)	= not used
!c
!c written 10/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer n,nz,nasvy,nparms,i,nfev
!c
	integer switch2,mario_DT_test
!c
      double precision y(n),z(nz),kRK(n),parms(nparms),deps(6)
      double precision F_sig(6),F_q(nasvy)
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
!c
	double precision zero
!c
	parameter(zero=0.0d0)
!c
!c      common /z_nct_errcode/error	
!c
!c ... update counter for the number of function f(y) evaluations
!c
      nfev=nfev+1
!c
!c ... initialize kRK
!c
      call pzero(kRK,n)
!c	
!c ... build F_sig(6) and F_q(nasv) vectors and move them into kRK
!c
      call get_F_sig_q(y,n,nasvy,z,nz,parms,nparms,deps,F_sig,F_q,&
      switch2,mario_DT_test,error)
		
	if(switch2.gt.zero) return

      if(error.eq.10) return
!c
      do i=1,6
        kRK(i)=F_sig(i)
      end do			 
!c	
      do i=1,nasvy
        kRK(6+i)=F_q(i)
      end do			 
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine get_F_sig_q(y,n,nasvy,z,nz,parms,nparms,deps,F_sig,F_q,&
      switch2,mario_DT_test,error)
!c-----------------------------------------------------------------------------
!c
!c  computes vectors F_sigma and F_q in F(y)
!c  Dafalias & Manzari (2004) SANISAND model for sands
!c
!c  variables allocated in vector y(n):
!c
!c variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c variables allocated in vector z(nasvz):
!c
!c	z(1)	= alpha_sr(1)
!c	z(2)	= alpha_sr(2)
!c	z(3)	= alpha_sr(3)
!c	z(4)	= alpha_sr(4)
!c	z(5)	= alpha_sr(5)
!c	z(6)	= alpha_sr(6)
!c	z(7)	= not used
!c	z(8)	= not used
!c	z(9)	= not used
!c	z(10)	= not used
!c	z(11)	= not used
!c	z(12)	= not used
!c	z(13)	= not used
!c	z(14)	= not used
!c
!c  written 10/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
      !external matmul
!c		
	integer switch2,mario_DT_test
!c 
      integer nparms,n,nasvy,nz
!c
      double precision y(n),z(nz),parms(nparms),deps(6)
      double precision Dep(6,6),HH(nasvy,6),F_sig(6),F_q(nasvy)
      double precision zero
	
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres

	parameter(zero=0.0d0)
!c
!c ... compute tangent operators
!c
      call get_tan_DM(y,n,nasvy,z,nz,parms,nparms,Dep,HH,switch2,&
      mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
	if(switch2.gt.zero) then
!c		write(*,*) 'get_tan - switch2>0'
		return
	endif
!c
!c ... compute F_sig=Dep*deps
!c
      call matmul(Dep,deps,F_sig,6,6,1)
!c
!c ... compute F_q=HH*deps
!c
      call matmul(HH,deps,F_q,nasvy,6,1)
!c	
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine get_tan_DM(y,ny,nasvy,z,nz,parms,nparms,Dep,Hep,&
      switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c computes matrices Dep and Hep
!c Dafalias & Manzari (2004) SANISAND model for sands
!c
!c variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c variables allocated in vector z(nasvz):
!c
!c	z(1)	= alpha_sr(1)
!c	z(2)	= alpha_sr(2)
!c	z(3)	= alpha_sr(3)
!c	z(4)	= alpha_sr(4)
!c	z(5)	= alpha_sr(5)
!c	z(6)	= alpha_sr(6)
!c	z(7)	= not used
!c	z(8)	= not used
!c	z(9)	= not used
!c	z(10)	= not used
!c	z(11)	= not used
!c	z(12)	= not used
!c	z(13)	= not used
!c	z(14)	= not used
!c
!c material properties allocated in vector parms(nparms):
!c
!c   1     p_a        Atmospheric pressure 
!c   2     e0         Void ratio on CSL at p = 0  
!c   3     lambda     CSL parameter (e:p plane)
!c   4     xi         CSL parameter (e:p plane)
!c   5     M_c        Slope of CSL in q:p plane, TX compression
!c   6     M_e        Slope of CSL in q:p plane, TX extension
!c   7     mm         opening of yield surface cone
!c   8     G0         Shear modulus constant
!c   9     nu         Poisson's ratio
!c   10    h0         Plastic modulus constant
!c   11    c_h        Plastic modulus constant
!c   12    n_b        Plastic modulus constant
!c   13    A0         Dilatancy constant
!c   14    n_d        Dilatancy constant
!c   15    z_max      Fabric index constant
!c   16    c_z        Fabric index constant
!c   17    bulk_w     Pore water bulk modulus (undrained conditions)
!c
!c  NOTE: stress and strain convention: compression positive
!c
!c  written 10/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
      !external matmul
!c 
      integer nparms,ny,nz,nasvy,i,j,switch2,iter,iter_max,switch4
	integer mario_DT_test
!c
      double precision distance!,psi_void_DM!dot_vect,
!c
      double precision y(ny),z(nz),parms(nparms)
      double precision De(6,6),Dep(6,6),Hep(nasvy,6),m(6)
      double precision LL(6),LL1(6),RR(6),RR1(6),U(6),V(6)
!c
      double precision p_a,e0,lambda,xi,M_c,M_e,cM,mm
      double precision G0,nu,h0,c_h,n_b
      double precision A0,n_d,z_max,c_z,bulk_w
!c
      double precision sig(6),alpha(6),void,Fab(6)
      double precision alpha_sr(6),alpha_b(6)
      double precision s(6),tau(6),n(6)
	double precision norm2,norm,I1,p,psi,cos3t,gth,dgdth
	double precision b0,d_sr,hh,db
      double precision Hplas,LDeR,Kp,Kpm1
      double precision mtrR,brack_mtrR,tol_ff,tol_dil,Hvs
      double precision h_alpha(6),h_fab(6),HH_alpha(6,6),HH_fab(6,6)
	double precision ff0,chvoid!yf_DM,
!c
      double precision zero,tiny,half,one,two,three,large,kappa
      double precision onethird,twothird
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
!c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(tiny=1.0d-15,large=1.0e15)
!c Heaviside function parameter
	parameter(kappa=3.0d2)
!c
!c      common /z_nct_errcode/error	
!c
      data m/1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0/
!c
	switch2=zero
	switch4=zero
	iter=0
	iter_max=1e3
!c
!c ... initialize constants and vectors
!c
      onethird=one/three
      twothird=two/three
      half=one/two
!c
      call pzero(Dep,36)
      call pzero(Hep,6*nasvy)
!c
!c ... recover material parameters
!c
      p_a=parms(1)
      e0=parms(2)
      lambda=parms(3)
	xi=parms(4)
      M_c=parms(5)
      M_e=parms(6)
      mm=parms(7)
      G0=parms(8)
      nu=parms(9) 
      h0=parms(10)
      c_h=parms(11)
      n_b=parms(12)
      A0=parms(13)
      n_d=parms(14)
      z_max=parms(15)
      c_z=parms(16)
      bulk_w=parms(17)
!c	  
	cM=M_e/M_c
!c
!c ... recover state variables
!c
      do i=1,6
        sig(i)=y(i)	
      end do !i
!c
      do i=1,6
        alpha(i)=y(6+i)	
      end do !i
!c
      void=y(13)
!c
      do i=1,6
        Fab(i)=y(13+i)	
      end do !i
!c
      do i=1,6
        alpha_sr(i)=z(i)	
      end do !i	  
!c
!c ... deviator stress and mean pressure
!c
      call deviator(sig,s,I1,p)
!c
!c ... stress ratio tensor and unit vector n
!c
      do i=1,6
        tau(i)=s(i)-p*alpha(i)
      end do ! i
!c
      norm2=dot_vect(1,tau,tau,6)
      norm=dsqrt(norm2)
      if(norm.lt.tiny) then
        norm=tiny
      endif
      do i=1,6
        n(i)=tau(i)/norm
!c		if(n(i).lt.tiny)then
!c		n(i)=zero
!c		endif
      end do
!c
!c ... elastic stiffness
!c
      call el_stiff_DM(y,ny,parms,nparms,De,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
!c ... gradient of yield function L and tensor V=De*L
!c
      call grad_f_DM(y,ny,parms,nparms,LL,LL1)
      call matmul(De,LL1,V,6,6,1)
!c
!c ... gradient of plastic potential R and tensor U=De*R
!c
      call grad_g_DM(y,ny,parms,nparms,RR,RR1)
      call matmul(De,RR1,U,6,6,1)
!c
!c ... plastic modulus functions b0 and hh
!c	  
	  if (dabs(p).gt.zero) then
	    chvoid=c_h*void
	    if(chvoid.ge.1) then
!c	     error=3
             chvoid=0.99999
	    end if
	    b0=G0*h0*(one-chvoid)/dsqrt(p/p_a)
	  else
	    b0=large
	  end if
!c
      d_sr=distance_sanisand(alpha,alpha_sr,n)
	  if (d_sr.lt.zero) then
		call push(alpha,alpha_sr,6)
!c		write(*,*) 'alpha updated'
	  end if
!c
	  if (d_sr.lt.tiny) then
		d_sr=tiny
	  end if
!c
	hh=b0/d_sr

!c
!c ... hardening function h_alpha	  
!c
      psi=psi_void_DM(void,p,parms,nparms)
      call lode_DM(tau,cM,cos3t,gth,dgdth)
      call alpha_th_DM(2,n,gth,psi,parms,nparms,alpha_b)
      db=distance_sanisand(alpha_b,alpha,n)
!c
      do i=1,6
        h_alpha(i)=twothird*hh*(alpha_b(i)-alpha(i))
      end do
!c
!c ... hardening function h_fab
!c
 	mtrR=-RR(1)-RR(2)-RR(3)
      brack_mtrR=half*(mtrR+dabs(mtrR))
!c
!c chiara heaviside function
!c
!cc	tol_dil=tol_ff*A0
!c
!c	if((-tol_dil.lt.mtrR).and.(tol_dil.gt.mtrR)) then
!c		Hvs=one/(1+exp(two*mtrR*kappa))
!c		bracK_mtrR=Hvs*mtrR
!c	endif
!c
      do i=1,6
	  h_fab(i)=-c_z*brack_mtrR*(z_max*n(i)+Fab(i))
      end do
!c
!c ... plastic moduli Hplas and Kp
!c
      Hplas=twothird*hh*p*db

	if(Hplas.gt.1e+15) then
!c		write(*,*)'Hplas'

	endif
!c
      LDeR=dot_vect(1,LL1,U,6)
!c
      Kp=LDeR+Hplas

	ff0=yf_DM(y,ny,parms,nparms)
!c .......................................................................
	
	if(mario_DT_test.eq.zero) then

		if(LDeR.lt.zero) then
		switch2=1
!c		write(*,*)'LDeR < 0'
		return
		endif

	

		if(Kp.lt.zero) then
		switch2=1
!c		write(*,*)'function get_tan: Kp < zero'
		return
		endif
	
	else
		

		if(LDeR.le.zero) then
		switch2=1
!c		error=3
!c		write(*,*)'subroutine get_tan_DM: LDeR < 0'
		return
		endif
	endif



	if(Kp.lt.zero)then
!c		write(6,*)'function get_tan: Kp < 0'
	error=3
	return
	endif
	call push(alpha_sr,z,6)
      Kpm1=one/Kp
!c	if (Kpm1.lt.tiny) then
!c		Kpm1=zero
!c	  end if
!c
!c ... elastoplastic stiffness matrix
!c
      do i=1,6
        do j=1,6
          Dep(i,j)=De(i,j)-Kpm1*U(i)*V(j)
        end do !j
      end do !i
!c
!c ... hardening tensor H_alpha
      do i=1,6
        do j=1,6
          HH_alpha(i,j)=Kpm1*h_alpha(i)*V(j)
        end do !j
      end do !i
!c
!c ... hardening tensor H_fab
!c		
      do i=1,6
        do j=1,6
          HH_fab(i,j)=Kpm1*h_fab(i)*V(j)
        end do !j
      end do !i
!c
!c ... Build tangent evolution matrix Hep(nasv,6) row-wise
!c
      do j=1,6
!c
        Hep(1,j) =HH_alpha(1,j)			! alpha(1)
        Hep(2,j) =HH_alpha(2,j)			! alpha(2)
        Hep(3,j) =HH_alpha(3,j)			! alpha(3)
        Hep(4,j) =HH_alpha(4,j)			! alpha(4)
        Hep(5,j) =HH_alpha(5,j)			! alpha(5)
        Hep(6,j) =HH_alpha(6,j)			! alpha(6)
        Hep(7,j) =-(one+void)*m(j)		! void
        Hep(8,j) =HH_fab(1,j)			! Fab(1)
        Hep(9,j) =HH_fab(2,j)			! Fab(2)
        Hep(10,j)=HH_fab(3,j)			! Fab(3)
        Hep(11,j)=HH_fab(4,j)			! Fab(4)
        Hep(12,j)=HH_fab(5,j)			! Fab(5)
        Hep(13,j)=HH_fab(6,j)			! Fab(6)
!c
      end do !j
!c	
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine grad_f_DM(y,ny,parms,nparms,gradf,gradf1)
!c------------------------------------------------------------------------------
!c ... subroutine to compute the gradient of yield function at state y 
!c     Dafalias & Manzari (2004) SNAISAND model for sands
!c
!c     P(i)  = grad(f) stress-like vector in Voigt notation
!c     P1(i) = grad(f) strain-like vector in Voigt notation
!c
!c ... variables allocated in parms
!c
!c  1     p_a        Atmospheric pressure 
!c  2     e0         Void ratio on CSL at p = 0  
!c  3     lambda     CSL parameter (e:p plane)
!c  4     xi         CSL parameter (e:p plane)
!c  5     M_c        Slope of CSL in q:p plane, TX compression
!c  6     M_e        Slope of CSL in q:p plane, TX extension
!c  7     mm         opening of yield surface cone
!c  8     G0         Shear modulus constant
!c  9     nu         Poisson's ratio
!c  10    h0         Plastic modulus constant
!c  11    c_h        Plastic modulus constant
!c  12    n_b        Plastic modulus constant
!c  13    A0         Dilatancy constant
!c  14    n_d        Dilatancy constant
!c  15    z_max      Fabric index constant
!c  16    c_z        Fabric index constant
!c  17    bulk_w     Pore water bulk modulus (undrained conditions)
!c
!c  variables allocated in vector y(n):
!c
!c       y(1) = sig(1)
!c       y(2) = sig(2)
!c       y(3) = sig(3)
!c       y(4) = sig(4)
!c       y(5) = sig(5)
!c       y(6) = sig(6)
!c       y(7) = alpha(1)
!c       y(8) = alpha(2)
!c       y(9) = alpha(3)
!c      y(10) = alpha(4)
!c      y(11) = alpha(5)
!c      y(12) = alpha(6)
!c      y(13) = void
!c      y(14) = Fab(1)
!c      y(15) = Fab(2)
!c      y(16) = Fab(3)
!c      y(17) = Fab(4)
!c      y(18) = Fab(5)
!c      y(19) = Fab(6)
!c      y(20) = not used
!c
!c  NOTE: soil mechanics convention (compression positive)
!c        all stress and strain vectors are 6-dimensional
!c------------------------------------------------------------------------------
      implicit none
!c
      !double precision dot_vect
!c
      integer ny,nparms,i
!c
      double precision parms(nparms),y(ny),gradf(6),gradf1(6),del(6)
      double precision mm,sig(6),s(6),r(6),I1,p
      double precision alpha(6),tau(6),n(6)
      double precision norm,norm2,v,vv
      double precision one,two,three,sqrt23,onethird,small
	double precision n1,n2
!c
      parameter(one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(small=1.0d-10)
	parameter(n1=0.816496580927739,n2=-0.40824829046385)
!c
      data del/1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0/
!c
      sqrt23=dsqrt(two/three)
      onethird=one/three
!c
      call pzero(n,6)
!c
!c ... recover material parameters
!c
      mm=parms(7)
!c
!c ... recover state variables
!c
      sig(1)=y(1)
      sig(2)=y(2)
      sig(3)=y(3)
      sig(4)=y(4)
      sig(5)=y(5)
      sig(6)=y(6)
!c
      alpha(1)=y(7)
      alpha(2)=y(8) 
      alpha(3)=y(9) 
      alpha(4)=y(10) 
      alpha(5)=y(11) 
      alpha(6)=y(12)
!c
!c ... deviator stress and mean pressure
!c
      call deviator(sig,s,I1,p)
!c
!c ... reduced stress tensor and unit vector
!c
      do i=1,6
        tau(i)=s(i)-p*alpha(i)
      end do ! i
!c
      norm2=dot_vect(1,tau,tau,6)
      norm=dsqrt(norm2)
!c
	if(norm.lt.small) then
	norm=small
	endif
!c
      do i=1,6
        n(i)=tau(i)/norm
      enddo
!c
!c	norm_n=dot_vect(1,n,n,6)
!c
!c	if(norm2.lt.small) then
!c	    n(1)=n1
!c		n(2)=n2
!c		n(3)=n2
!c      endif
!c
!c
!c ... coefficient V
!c
      if(dabs(p).lt.small) then
        do i=1,6
          r(i)=s(i)/small
        enddo
      else
        do i=1,6
          r(i)=s(i)/p
        enddo
      endif
      v=dot_vect(1,r,n,6)
      vv=-onethird*v
!c
!c ... gradient of f
!c
      do i=1,6
        gradf(i)=n(i)+vv*del(i)
        if(i.le.3) then
          gradf1(i)=gradf(i)
        else
          gradf1(i)=two*gradf(i)
        endif
      enddo
!c
      return
      end
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c-----------------------------------------------------------------------------
      subroutine grad_g_DM(y,ny,parms,nparms,gradg,gradg1)
!c------------------------------------------------------------------------------
!c ... subroutine to compute the gradient of plastic potential at state y 
!c     Dafalias & Manzari (2004) SANISAND model for sands
!c
!c     gradg(i)  = grad(g) stress-like vector in Voigt notation
!c     gradg1(i) = grad(g) strain-like vector in Voigt notation
!c
!c ... variables allocated in parms
!c
!c  1     p_a        Atmospheric pressure 
!c  2     e0         Void ratio on CSL at p = 0  
!c  3     lambda     CSL parameter (e:p plane)
!c  4     xi         CSL parameter (e:p plane)
!c  5     M_c        Slope of CSL in q:p plane, TX compression
!c  6     M_e        Slope of CSL in q:p plane, TX extension
!c  7     mm         opening of yield surface cone
!c  8     G0         Shear modulus constant
!c  9     nu         Poisson's ratio
!c  10    h0         Plastic modulus constant
!c  11    c_h        Plastic modulus constant
!c  12    n_b        Plastic modulus constant
!c  13    A0         Dilatancy constant
!c  14    n_d        Dilatancy constant
!c  15    z_max      Fabric index constant
!c  16    c_z        Fabric index constant
!c  17    bulk_w     Pore water bulk modulus (undrained conditions)
!c
!c  variables allocated in vector y(n):
!c
!c       y(1) = sig(1)
!c       y(2) = sig(2)
!c       y(3) = sig(3)
!c       y(4) = sig(4)
!c       y(5) = sig(5)
!c       y(6) = sig(6)
!c       y(7) = alpha(1)
!c       y(8) = alpha(2)
!c       y(9) = alpha(3)
!c      y(10) = alpha(4)
!c      y(11) = alpha(5)
!c      y(12) = alpha(6)
!c      y(13) = void
!c      y(14) = Fab(1)  (stress--like)
!c      y(15) = Fab(2)
!c      y(16) = Fab(3)
!c      y(17) = Fab(4)
!c      y(18) = Fab(5)
!c      y(19) = Fab(6)
!c      y(20) = not used
!c
!c  NOTE: soil mechanics convention (compression positive)
!c        all stress and strain vectors are 6-dimensional
!c------------------------------------------------------------------------------
      implicit none
!c
      double precision distance,psi_void!,psi_void_DM!dot_vect,
!c
      integer ny,nparms,i
!c
      double precision M_c,M_e,cM,A0
!c
      double precision parms(nparms),y(ny),gradg(6),gradg1(6)
      double precision sig(6),s(6),alpha(6),Fab(6),I1,p
      double precision n(6),n2(6),tau(6),Rdev(6)
      double precision Ad,alpha_d(6),dd	  
	double precision cos3t,gth,dgdth
      double precision void,psi,dil,dil3
!c
      double precision temp1,temp2,temp3,temp4
      double precision norm,norm2
      double precision zero,one,two,three,six
	double precision half,sqrt6,onethird,small,del(6)
!c
	integer chiara
!c
      parameter(half=0.5d0,one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0)
      parameter(zero=0.0d0,small=1.0d-10)
!c
      data del/1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0/
!c
      sqrt6=dsqrt(six)
      onethird=one/three
!c
      call pzero(n,6)
!c
!c ... recover material parameters
!c
      M_c=parms(5)
      M_e=parms(6)
      A0=parms(13)
!c
      cM=M_e/M_c
!c
!c ... recover state variables
!c
      sig(1)=y(1)
      sig(2)=y(2)
      sig(3)=y(3)
      sig(4)=y(4)
      sig(5)=y(5)
      sig(6)=y(6)
!c
      alpha(1)=y(7)
      alpha(2)=y(8) 
      alpha(3)=y(9) 
      alpha(4)=y(10) 
      alpha(5)=y(11) 
      alpha(6)=y(12)
!c
      void=y(13)
!c
      Fab(1)=y(14)
      Fab(2)=y(15)
      Fab(3)=y(16)
      Fab(4)=y(17)
      Fab(5)=y(18)
      Fab(6)=y(19)
!c
!c ... deviator stress and mean pressure
!c
      call deviator(sig,s,I1,p)
!c
!c ... stress ratio tensor, unit tensor n and tensor n^2 
!c
      do i=1,6
        tau(i)=s(i)-p*alpha(i)
      end do ! i
!c
      norm2=dot_vect(1,tau,tau,6)
      norm=dsqrt(norm2)
!c
      if(norm.lt.small) then
	  norm=small
      endif
      do i=1,6
          n(i)=tau(i)/norm
      enddo
!c
      n2(1)=n(1)*n(1)+n(4)*n(4)+n(5)*n(5)
      n2(2)=n(4)*n(4)+n(2)*n(2)+n(6)*n(6)
      n2(3)=n(6)*n(6)+n(5)*n(5)+n(3)*n(3)
      n2(4)=n(1)*n(4)+n(4)*n(2)+n(6)*n(5)
      n2(5)=n(5)*n(1)+n(6)*n(4)+n(3)*n(5)
      n2(6)=n(4)*n(5)+n(2)*n(6)+n(6)*n(3)
!c
!c ... state parameter psi
!c
      psi=psi_void_DM(void,p,parms,nparms)
!c
!c ... Lode angle; functions g(theta) and (1/g)dg/dtheta
!c
      call lode_DM(tau,cM,cos3t,gth,dgdth)
!c
!c ... vector Rdev
!c	C&M
      temp1=one+three*cos3t*dgdth
      temp2=-three*sqrt6*dgdth
      do i=1,6
          Rdev(i)=temp1*n(i)+temp2*(n2(i)-onethird*del(i))
      enddo
!c	  	  
!c ... dilatancy function
!c
	temp3=dot_vect(1,Fab,n,6)
	temp4=half*(temp3+dabs(temp3))
      Ad=A0*(one+temp4)	  
!c	  
      call alpha_th_DM(3,n,gth,psi,parms,nparms,alpha_d)
      dd = distance_sanisand(alpha_d,alpha,n)
!c Chiara
	if((psi.gt.zero).and.(dd.lt.zero)) then
		dd=zero
	endif
!c end Chiara
!c
      dil=Ad*dd
      dil3=onethird*dil
!c
      do i=1,6
        gradg(i)=Rdev(i)+dil3*del(i)
        if(i.le.3) then
          gradg1(i)=gradg(i)
        else
          gradg1(i)=two*gradg(i)
        endif
      enddo
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine iniyz(y,nydim,z,nzdim,qq1,nasvy,qq2,nasvz,sig,ntens)
!c-----------------------------------------------------------------------------
!c  initializes the vectors of state variables
!c
!c  variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c  variables allocated in vector z(nasvz):
!c
!c	z(1)	= alpha_sr(1)
!c	z(2)	= alpha_sr(2)
!c	z(3)	= alpha_sr(3)
!c	z(4)	= alpha_sr(4)
!c	z(5)	= alpha_sr(5)
!c	z(6)	= alpha_sr(6)
!c	z(7)	= not used
!c	z(8)	= not used
!c	z(9)	= not used
!c	z(10)	= not used
!c	z(11)	= not used
!c	z(12)	= not used
!c	z(13)	= not used
!c	z(14)	= not used
!c
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer i,nydim,nzdim,nasvy,nasvz,ntens
!c
      double precision y(nydim),z(nzdim)
      double precision qq1(nasvy),qq2(nasvz),sig(ntens)
!c
      call pzero(y,nydim)
      call pzero(z,nzdim)
!c
      do i=1,ntens
        y(i) = sig(i)
      enddo
!c
      do i=1,nasvy
        y(6+i) = qq1(i)
      enddo
!c
      do i=1,nasvz
        z(i) = qq2(i)
      enddo
!c
      return
      end
!c-----------------------------------------------------------------------------
	subroutine intersect_DM(y0,y1,y_star,n,parms,nparms,tol_ff,&
         xi,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c
!c ... finds the intersection point between the stress path 
!c     and the yield surface (Papadimitriou & Bouckovalas model)
!c     using Newton method
!c
!c  variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c  variables allocated in vector z(nasvz):
!c
!c	z(1)	= alpha_sr(1)
!c	z(2)	= alpha_sr(2)
!c	z(3)	= alpha_sr(3)
!c	z(4)	= alpha_sr(4)
!c	z(5)	= alpha_sr(5)
!c	z(6)	= alpha_sr(6)
!c	z(7)	= not used
!c	z(8)	= not used
!c	z(9)	= not used
!c	z(10)	= not used
!c	z(11)	= not used
!c	z(12)	= not used
!c	z(13)	= not used
!c	z(14)	= not used
!c
!c ... variables allocated in parms
!c
!c   1     p_a        Atmospheric pressure 
!c   2     e0         Void ratio on CSL at p = 0  
!c   3     lambda     CSL parameter (e:p plane)
!c   4     xi         CSL parameter (e:p plane)
!c   5     M_c        Slope of CSL in q:p plane, TX compression
!c   6     M_e        Slope of CSL in q:p plane, TX extension
!c   7     mm         opening of yield surface cone
!c   8     G0         Shear modulus constant
!c   9     nu         Poisson's ratio
!c   10    h0         Plastic modulus constant
!c   11    c_h        Plastic modulus constant
!c   12    n_b        Plastic modulus constant
!c   13    A0         Dilatancy constant
!c   14    n_d        Dilatancy constant
!c	15    z_max      Fabric index constant
!c   16    c_z        Fabric index constant
!c   17    bulk_w     Pore water bulk modulus (undrained conditions)
!c
!c  Tamagnini 10/2008
!c
!c-----------------------------------------------------------------------------
!c
      implicit none
!c
      integer n,nparms,maxiter,kiter,i,kiter_bis,bisect
!c
      !double precision dot_vect!yf_DM,
!c
      double precision parms(nparms),y0(n),y1(n),y_star(n),y05(n)
      double precision tol_ff,fy_star,err,dfdxi,dfdxi_m1,xi,fy05
	double precision dxi, xip1
      double precision sig0(6),sig1(6),dsig(6),P_star(6),P1_star(6)
      double precision zero,one,half,three,onethird
	double precision pp_star,low, fy11, fy00, xi_max, xi_i, pp05
	double precision y00(n),y11(n)
	
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
!c
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0,three=3.0d0)
	parameter(low=1.0d-10)
!c
!c      common /z_nct_errcode/error	
!c
      xi=one
      maxiter=200000
      kiter=0
	bisect=0
	kiter_bis=0
!c
      do i=1,6
        sig0(i)=y0(i)
        sig1(i)=y1(i)
        dsig(i)=sig1(i)-sig0(i)
      end do !i
!c
      call push(y1,y_star,n)
!c
      fy_star=yf_DM(y_star,n,parms,nparms)
	onethird=one/three
	pp_star=(y_star(1)+y_star(2)+y_star(3))*onethird
      err=dabs(fy_star/pp_star)
	if(pp_star.gt.one) err=dabs(fy_star)	
!c	
!c
!c	
	if(bisect.eq.0) then
!c
!c
!c ... start Newton iteration
!c
      do while ((err.gt.tol_ff).and.(bisect.eq.0))
!c
        kiter=kiter+1
!c
        call grad_f_DM(y_star,n,parms,nparms,P_star,P1_star)
!c
        dfdxi=dot_vect(1,P_star,dsig,6)
		if (dfdxi.lt.low) then
!ccc			dfdxi=low
!ccc			write(*,*)'dfdxi.lt.zero'
			bisect=1
		endif
        dfdxi_m1=one/dfdxi
!c
!c ... search direction
!c
        dxi=-dfdxi_m1*fy_star
        xip1=xi+dxi
!c
!c	Line search
!c
	  do while ((xip1.lt.zero).or.(xip1.gt.one))
		dxi=half*dxi
	    xip1=xi+dxi
	  end do
!c
	  xi=xip1
!c
!c	End Line search
!c
!c ... find new intersection point and yield function value
!c
        do i=1,n
          y_star(i)=y0(i)+xi*(y1(i)-y0(i))
        end do !i
!c
        fy_star=yf_DM(y_star,n,parms,nparms)
		if (fy_star.lt.zero) then
!c			write(*,*)'fy_star.lt.zero'
			bisect=1
		else
		onethird=one/three
		pp_star=(y_star(1)+y_star(2)+y_star(3))*onethird
		err=dabs(fy_star/pp_star)
		if(pp_star.gt.one) err=dabs(fy_star)	
!cccc		err=dabs(fy_star)
		endif
!c
	if (kiter.gt.maxiter+1) then
          write(6,*) 'ERROR: max no. of iterations exceeded'
          write(6,*) 'Subroutine INTERSECT_DM'
          write(6,*) 'err = ',err
!c          error=10
!c	   error=3
!c          return 
	err=0
        end if

	
      end do ! bottom of Newton iteration
!c
!c ... check that 0 < xi < 1 (intersection point between initial and final states)
!c
      if((xi.lt.zero).and.(xi.gt.one)) then 
!c
        write(6,*) 'ERROR: the intersection point found lies'
        write(6,*) '       outside the line connecting initial'
        write(6,*) '       and trial stress states'
        write(6,*) 'Subroutine INTERSECT_DM'
        write(6,*) 'xi = ',xi
        xi = zero
!c	    error=10
!c	    error=3
        return 
!c	
	endif
	endif

!c
	if(bisect.eq.1) then
!c ... start bisection method
!ccc	write(*,*) 'bisection method'
!c
!c		find f((a+b)/2)
	
	do i=1,n
		y00(i)=y0(i)
		y11(i)=y1(i)
	enddo
	fy00 =yf_DM(y00,n,parms,nparms)
	fy11 =yf_DM(y11,n,parms,nparms)
	do i=1,n
		y05(i)=y0(i)
	enddo
	pp05=(y05(1)+y05(2)+y05(3))*onethird	
	fy05 =yf_DM(y05,n,parms,nparms)
!cccc	err = dabs(fy05)
	err=abs(fy05/pp05)
	if(pp05.gt.one) err=dabs(fy05)	

	do while(err.gt.tol_ff)
		kiter_bis=kiter_bis+1
!c
		do i=1,6
			y05(i)=half*(y00(i)+y11(i))
		enddo	
	
		fy05 =yf_DM(y05,n,parms,nparms)
		pp05=(y05(1)+y05(2)+y05(3))*onethird
!cccc		err = dabs(fy05)
		err=abs(fy05/pp05)
		if(pp05.gt.one) err=dabs(fy05)	
	
		if(fy05.lt.zero) then
			call push(y05,y00,n)
		else
			call push(y05,y11,n)
		endif
		
		if (kiter_bis.gt.maxiter+1) then
			write(6,*) 'ERROR: max no. of iterations exceeded'
			write(6,*) 'Subroutine INTERSECT_DM - bisection'
			write(6,*) 'err = ',err
			err=0
!c          		error=10
!c	                error=3
!c          		return 
	  	endif
	enddo

	do i=1,n
		y_star(i)=y05(i)
	enddo
		
!c	xi= (y05(1)-y0(1))/(y1(1)-y0(1))

	xi_max=zero
	do i=1,6
		if((y1(i)-y0(i)).ne.zero) then
		  xi_i= (y05(i)-y0(i))/(y1(i)-y0(i))
		  if(xi_i.gt.xi_max) then
			xi_max = xi_i
		  endif
		endif
	enddo
	xi = xi_max
	
!c ... end bisection method	
!c
	
	endif

!c 
      return
      end
!c
!c------------------------------------------------------------------------------
      subroutine inv_sig(sig,pp,qq,cos3t)
!c------------------------------------------------------------------------------
!c calculate invariants of stress tensor
!c Dafalias &Manzari SANISAND model (2004)
!c
!c NOTE: Voigt notation is used with the following index conversion
!c
!c       11 -> 1
!c       22 -> 2
!c       33 -> 3
!c       12 -> 4
!c       13 -> 5
!c       23 -> 6
!c
!c------------------------------------------------------------------------------
!c
      implicit none
!c
      double precision sig(6),sdev(6),s2(6)
      double precision I1,J2bar,J2bar_sq,J3bar,trs2,trs3
      double precision pp,qq,cos3t,numer,denom
!c
      double precision zero,one,two,three
      double precision onethird,half,onept5,sqrt3,tiny
!c
      !double precision dot_vect
!c
      data zero,one,two,three/0.0d0,1.0d0,2.0d0,3.0d0/
      data tiny/1.0d-15/
!c
!c ... some constants
!c
      onethird=one/three
      half=one/two
      onept5=three/two
      sqrt3=dsqrt(three)
!c
!c ... trace and mean stress
!c
      I1=sig(1)+sig(2)+sig(3)
      pp=onethird*I1
!c
!c ... deviator stress
!c
      sdev(1)=sig(1)-pp
      sdev(2)=sig(2)-pp
      sdev(3)=sig(3)-pp
      sdev(4)=sig(4)
      sdev(5)=sig(5)
      sdev(6)=sig(6)
!c
!c ... second invariants
!c
      trs2=dot_vect(1,sdev,sdev,6)
      J2bar=half*trs2
      qq=dsqrt(onept5*trs2)
!c
!c ... components of (sdev_ij)(sdev_jk) (stress-like Voigt vector)
!c
      s2(1)=sdev(1)*sdev(1)+sdev(4)*sdev(4)+sdev(5)*sdev(5)
      s2(2)=sdev(4)*sdev(4)+sdev(2)*sdev(2)+sdev(6)*sdev(6)
      s2(3)=sdev(6)*sdev(6)+sdev(5)*sdev(5)+sdev(3)*sdev(3)
      s2(4)=sdev(1)*sdev(4)+sdev(4)*sdev(2)+sdev(6)*sdev(5)
      s2(5)=sdev(5)*sdev(1)+sdev(6)*sdev(4)+sdev(3)*sdev(5)
      s2(6)=sdev(4)*sdev(5)+sdev(2)*sdev(6)+sdev(6)*sdev(3)
!c	     
!c ... Lode angle
!c
      if(trs2.lt.tiny) then 
!c
        cos3t=one
!c		
      else
!c
        trs3=dot_vect(1,sdev,s2,6)
!c
        J3bar=onethird*trs3
        J2bar_sq=dsqrt(J2bar)
        numer=three*sqrt3*J3bar
        denom=two*(J2bar_sq**3)
        cos3t=numer/denom
        if(dabs(cos3t).gt.one) then
          cos3t=cos3t/dabs(cos3t)
        end if
!c
      end if 
!c
      return
      end
!c
!c------------------------------------------------------------------------------
      subroutine lode_DM(r,cM,cos3t,gth,dgdth)
!c------------------------------------------------------------------------------
!c calculate cos(3*theta) from deviatoric stress ratio tensor 
!c tau(i) = s(i)-p*alpha(i) stored in vector r(6)
!c
!c computes functions g(theta) and (1/g)dg/dtheta from:
!c a) Argyris function (Argyris=1) (as in the original paper)
!c b) Van Eekelen function (Argyris=0) (more appropriate for high friction angles)
!c
!c       gth   = g(theta)
!c       dgdth = (1/g)dg/dtheta
!c
!c NOTE: Voigt notation is used with the following index conversion (ABAQUS)
!c
!c       11 -> 1
!c       22 -> 2
!c       33 -> 3
!c       12 -> 4
!c       13 -> 5
!c       23 -> 6
!c       
!c SM stress convention: compression positive
!c
!c------------------------------------------------------------------------------
!c
      implicit none
!c
      integer Argyris
!c
      double precision r(6),r2(6)
      double precision trr2,trr3,J2bar,J3bar,J2bar_sq
      double precision cM,n_VE,n_VEm1,numer,denom,cos3t
!c
      double precision tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
	double precision alpha,beta,gth,dgdth
      double precision one,two,three
      double precision onethird,half,sqrt3,tiny
!c
      !double precision dot_vect
!c
      data one,two,three/1.0d0,2.0d0,3.0d0/
      data tiny,n_VE/1.0d-15,-0.25d0/
	data Argyris/0/
!c
!c ... some constants
!c
      onethird=one/three
      half=one/two
      sqrt3=dsqrt(three)
!c
!c ... second invariant
!c
      trr2=dot_vect(1,r,r,6)
      J2bar=half*trr2
!c
!c ... components of (r_ij)(r_jk) (stress-like Voigt vector)
!c
      r2(1)=r(1)*r(1)+r(4)*r(4)+r(5)*r(5)
      r2(2)=r(4)*r(4)+r(2)*r(2)+r(6)*r(6)
      r2(3)=r(6)*r(6)+r(5)*r(5)+r(3)*r(3)
      r2(4)=r(1)*r(4)+r(4)*r(2)+r(6)*r(5)
      r2(5)=r(5)*r(1)+r(6)*r(4)+r(3)*r(5)
      r2(6)=r(4)*r(5)+r(2)*r(6)+r(6)*r(3)
!c	     
!c ... Lode angle
!c
      if(trr2.lt.tiny) then 
!c
        cos3t=one
!c		
      else
!c
        trr3=dot_vect(1,r,r2,6)
!c
        J3bar=onethird*trr3
        J2bar_sq=dsqrt(J2bar)
        numer=three*sqrt3*J3bar
        denom=two*(J2bar_sq**3)
        cos3t=numer/denom
        if(dabs(cos3t).gt.one) then
          cos3t=cos3t/dabs(cos3t)
        end if
!c
      end if 
!c	     
!c ... g function and its derivative
!c
	  if (Argyris.ne.0) then
!c
!c ... Argyris function
!c
	  gth=two*cM/((one+cM)-(one-cM)*cos3t)
	  dgdth=(1-cM)*gth/(two*cM) 
!c
        else
!c
!c ... Van Eekelen function
!c
        n_VEm1=one/n_VE
!c
        tmp1=one/(two**n_VE)
	  tmp2=cM**n_VEm1
        tmp3=one+tmp2
        tmp4=one-tmp2
!c
	  alpha=tmp1*(tmp3**n_VE)
        beta=tmp4/tmp3
!c
        tmp5=(one+beta*cos3t)**n_VE
        tmp6=one+beta*cos3t
!c
        gth=alpha*tmp5
        dgdth=n_VE*beta/tmp6
!c
        end if
!c
      return
      end
!c
!c------------------------------------------------------------------------------
      subroutine matmul(a,b,c,l,m,n)
!c------------------------------------------------------------------------------
!c matrix multiplication
!c------------------------------------------------------------------------------
      implicit none
!c
      integer i,j,k,l,m,n
!c
      double precision a(l,m),b(m,n),c(l,n)
!c
      do i=1,l
        do j=1,n
          c(i,j) = 0.0d0
          do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          enddo
        enddo
      enddo
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine move_eps(dstran,ntens,deps,depsv)
!c-----------------------------------------------------------------------------
!c Move strain/strain increment stran/dstran into eps/deps, 
!c computes volumetric strain/strain increment and switches 
!c sign convention from solid to soil mechanics
!c
!c NOTE: 
!c   stran  = strain tensor (extension positive)
!c	eps    = strain tensor (compression positive)
!c	epsv   = vol. strain (compression positive)
!c   dstran = strain increment tensor (extension positive)
!c	deps   = strain increment tensor (compression positive)
!c	depsv  = vol. strain increment (compression positive)
!c
!c   eps/deps has always 6 components
!c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision deps(6),dstran(ntens),depsv
!c
      call pzero(deps,6)
!c
      do i=1,ntens
        deps(i) = -dstran(i)
      enddo
!c
      depsv=deps(1)+deps(2)+deps(3)
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine move_sig(stress,ntens,pore,sig)
!c-----------------------------------------------------------------------------
!c Computes effective stress from total stress (stress) and pore pressure (pore)
!c and switches sign convention from solid to soil mechanics
!c
!c NOTE: stress = total stress tensor (tension positive)
!c       pore   = exc. pore pressure (undrained conds., compression positive)
!c       sig    = effective stress (compression positive)
!c
!c       sig has always 6 components
!c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision sig(6),stress(ntens),pore
!c
      call pzero(sig,6)
!c
      do i=1,ntens
        if(i.le.3) then
          sig(i) = -stress(i)-pore
        else
          sig(i) = -stress(i)
        end if
      enddo
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine norm_res_DM(y_til,y_hat,ny,norm_R)
!c-----------------------------------------------------------------------------
!c  evaluate norm of residual vector Res=||y_hat-y_til||
!c  Dafalias & Manzari(2004) SANISAND model for sand
!c
!c  variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c  written 10/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
	  implicit none
!c 
      integer ny,i
!c
      double precision y_til(ny),y_hat(ny)
      double precision err(ny),norm_R2,norm_R
      double precision sig_hat(6),sig_til(6),del_sig(6)
      double precision alpha_hat(6),alpha_til(6),del_alpha(6)
	double precision Fab_hat(6),Fab_til(6),del_Fab(6)
      double precision void_hat,void_til,del_void
      double precision norm_sig2,norm_alpha2,norm_Fab2
      double precision norm_sig,norm_alp,norm_Fab
      double precision zero!dot_vect,
!c
      parameter(zero=0.0d0)
!c
      call pzero(err,ny)
!c
!c ... recover stress tensor and internal variables
!c
      do i=1,6
        sig_hat(i)=y_hat(i)
        sig_til(i)=y_til(i)
        del_sig(i)=dabs(sig_hat(i)-sig_til(i))
      end do
!c
      do i=1,6
        alpha_hat(i)=y_hat(6+i)
        alpha_til(i)=y_til(6+i)
        del_alpha(i)=dabs(alpha_hat(i)-alpha_til(i))
      end do
!c
      void_hat=y_hat(13)
      void_til=y_til(13)
      del_void=dabs(void_hat-void_til)
!c
      do i=1,6
        Fab_hat(i)=y_hat(13+i)
        Fab_til(i)=y_til(13+i)
        del_Fab(i)=dabs(Fab_hat(i)-Fab_til(i))
      end do
!c
!c ... relative error norms
!c
      norm_sig2=dot_vect(1,sig_hat,sig_hat,6)
      norm_alpha2=dot_vect(1,alpha_hat,alpha_hat,6)
	norm_Fab2=dot_vect(1,Fab_hat,Fab_hat,6)
      norm_sig=dsqrt(norm_sig2)
      norm_alp=dsqrt(norm_alpha2)
	norm_Fab=dsqrt(norm_Fab2)
!c
      if(norm_sig.gt.zero) then
        do i=1,6
          err(i)=del_sig(i)/norm_sig
        end do
      end if
!c
      if(norm_alp.gt.zero) then
        do i=1,6
          err(6+i)=del_alpha(i)/norm_alp
        end do
      end if
!c
      err(13)=del_void/void_hat
!c
!c      if(norm_Fab.gt.zero) then
!c        do i=1,6
!c          err(13+i)=del_Fab(i)/norm_Fab
!c        end do
!c      end if
!c
!c chiara 4 maggio 2010
!c
      do i=1,6
		if((Fab_til(i).ne.zero).and.(norm_Fab.gt.zero)) then
			err(13+i)=del_Fab(i)/norm_Fab
		end if
      end do    
!c
!c ... global relative error norm
!c
	norm_R2=dot_vect(3,err,err,ny)
      norm_R=dsqrt(norm_R2)
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine pert_DM(y_n,y_np1,z,n,nasvy,nasvz,err_tol,&
                        maxnint,DTmin,deps_np1,parms,&
                        nparms,nfev,elprsw,theta,ntens,DD,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c  compute numerically consistent tangent stiffness 
!c  Dafalias & Manzari (2004) SANISAND model for sand
!c
!c  written 10/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
!c 
      integer elprsw
!c
      integer ntens,jj,kk
      integer n,nasvy,nasvz,nparms,nfev
      integer maxnint,mario_DT_test
!c
      double precision y_n(n),y_np1(n),y_star(n),z(nasvz),parms(nparms)
      double precision err_tol
      double precision theta,DTmin
      double precision deps_np1(6),deps_star(6)
      double precision dsig(6),DD(6,6)
      double precision zero,three
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
!c
      parameter(zero=0.0d0,three=3.0d0)
!c
!c      common /z_nct_errcode/error
!c      common /z_plastic_flag/plastic
!c
!c ... initialize DD and y_star
!c 
      call pzero(DD,36)
      call pzero(y_star,n)
!c
      if(plastic.eq.0) then
!c
!c ... elastic process, DD = De (explicit)
!c
        call el_stiff_DM(y_np1,n,parms,nparms,DD,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
      else
!c
!c ... plastic process, DD computed using numerical perturbation
!c     loop over strain basis vectors
!c
	do jj=1,ntens
        call push(y_n,y_star,n)
        call push(deps_np1,deps_star,ntens)
!c
!c        do jj=1,ntens
!c
!c ... perturbed strain increment
!c
	      deps_star(jj)=deps_star(jj)+theta
!c
!c ... perturbed final state, stored in y_star
!c
          if(error.ne.10) then
            call rkf23_upd_DM(y_star,z,n,nasvy,nasvz,err_tol,maxnint,&
                DTmin,deps_star,parms,nparms,nfev,elprsw,&
     		 mario_DT_test,&
                error,tol_f,check_ff,drcor,p_thres,plastic)
          end if
!c
!c ... stiffness components in column jj (kk = row index, jj = column index)
!c
          do kk=1,ntens
            dsig(kk)=y_star(kk)-y_np1(kk)
            DD(kk,jj)=dsig(kk)/theta
          end do !kk
!c
        end do !jj
!c
      end if	
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine plast_mod_DM(y,ny,z,nz,parms,nparms,h_alpha,Kpm1,&
      switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c  computes vector h_hard and plastic modulus Kp
!c  Dafalias & Manzari(2004) SANISAND model for sand
!c
!c  variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c  variables allocated in vector z(nasvz):
!c
!c	z(1)	= alpha_sr(1)
!c	z(2)	= alpha_sr(2)
!c	z(3)	= alpha_sr(3)
!c	z(4)	= alpha_sr(4)
!c	z(5)	= alpha_sr(5)
!c	z(6)	= alpha_sr(6)
!c	z(7)	= not used
!c	z(8)	= not used
!c	z(9)	= not used
!c	z(10)	= not used
!c	z(11)	= not used
!c	z(12)	= not used
!c	z(13)	= not used
!c	z(14)	= not used
!c
!c  material properties allocated in vector parms(nparms):
!c
!c   1     p_a        Atmospheric pressure 
!c   2     e0         Void ratio on CSL at p = 0  
!c   3     lambda     CSL parameter (e:p plane)
!c   4     xi         CSL parameter (e:p plane)
!c   5     M_c        Slope of CSL in q:p plane, TX compression
!c   6     M_e        Slope of CSL in q:p plane, TX extension
!c   7     mm         opening of yield surface cone
!c   8     G0         Shear modulus constant
!c   9     nu         Poisson's ratio
!c   10    h0         Plastic modulus constant
!c   11    c_h        Plastic modulus constant
!c   12    n_b        Plastic modulus constant
!c   13    A0         Dilatancy constant
!c   14    n_d        Dilatancy constant
!c   15    z_max      Fabric index constant
!c   16    c_z        Fabric index constant
!c   17    bulk_w     Pore water bulk modulus (undrained conditions)
!c
!c
!c  NOTE: stress and strain convention: compression positive
!c
!c  written 10/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
      !external matmul
!c 
      integer nparms,ny,nz,i,iter,switch2,mario_DT_test
!c
      double precision distance,ref_db,psi_void!,psi_void_DM!dot_vect,
!c
      double precision y(ny),z(nz),parms(nparms)
      double precision De(6,6),h_alpha(6)
      double precision LL(6),LL1(6),RR(6),RR1(6),U(6),V(6)
!c
      double precision p_a,e0,lambda,xi,M_c,M_e,cM,mm
      double precision G0,nu,h0,c_h,n_b
      double precision A0,n_d,z_max,c_z,bulk_w
!c
      double precision sig(6),alpha(6),void,Fab(6)
      double precision alpha_sr(6),alpha_b(6)
      double precision s(6),tau(6),n(6),I1,p,psi,cos3t
	double precision b0,d_sr,hh,db,gth,dgdth
      double precision HHp,LDeR,Kp,Kpm1,norm2,norm,chvoid
!c
      double precision zero,tiny,one,two,three,large
      double precision twothird
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
!c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(tiny=1.0d-15,large=1.0e15)
!c
!c      common /z_nct_errcode/error	
!c
!c ... initialize constants and vectors
!c
      twothird=two/three
!c
!c ... recover material parameters
!c
      p_a=parms(1)
      e0=parms(2)
      lambda=parms(3)
	xi=parms(4)
      M_c=parms(5)
      M_e=parms(6)
      mm=parms(7)
      G0=parms(8)
      nu=parms(9) 
      h0=parms(10)
      c_h=parms(11)
      n_b=parms(12)
      A0=parms(13)
      n_d=parms(14)
      z_max=parms(15)
      c_z=parms(16)
      bulk_w=parms(17)
!c	  
	  cM=M_e/M_c
!c
	switch2=zero
	iter=0
!c
!c ... recover state variables
!c
      do i=1,6
        sig(i)=y(i)	
      end do !i
!c
      do i=1,6
        alpha(i)=y(6+i)	
      end do !i
!c
      void=y(13)
!c
      do i=1,6
        Fab(i)=y(13+i)	
      end do !i
!c
      do i=1,6
        alpha_sr(i)=z(i)	
      end do !i	  
!c
!c ... deviator stress and mean pressure
!c
      call deviator(sig,s,I1,p)
!c
!c ... stress ratio tensor and unit vector n
!c
      do i=1,6
        tau(i)=s(i)-p*alpha(i)
      end do ! i
!c
      norm2=dot_vect(1,tau,tau,6)
      norm=dsqrt(norm2)
      if(norm.lt.tiny) then
        norm=tiny
      endif
      do i=1,6
        n(i)=tau(i)/norm
!c		if(n(i).lt.tiny)then
!c		n(i)=zero
!c		endif
      end do
!c
!c ... elastic stiffness
!c
      call el_stiff_DM(y,ny,parms,nparms,De,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
!c ... gradient of yield function L and tensor V=De*L
!c
      call grad_f_DM(y,ny,parms,nparms,LL,LL1)
      call matmul(De,LL1,V,6,6,1)
!c
!c ... gradient of plastic potential R and tensor U=De*R
!c
      call grad_g_DM(y,ny,parms,nparms,RR,RR1)
      call matmul(De,RR1,U,6,6,1)
!c
!c ... plastic modulus functions b0 and hh
!c	  
	  if (dabs(p).gt.zero) then
	    chvoid=c_h*void
	    if(chvoid.ge.1) then
!c	     error=3
             chvoid=0.99999
	    end if
	    b0=G0*h0*(one-chvoid)/dsqrt(p/p_a)
	  else
	    b0=large
	  end if
!c
      d_sr=distance_sanisand(alpha,alpha_sr,n)
	  if (d_sr.lt.zero) then
		call push(alpha,alpha_sr,6)
!c		write(*,*) 'alpha updated'
	  end if
!c
	  if (d_sr.lt.tiny) then
		d_sr=tiny
	  end if
!c
	hh=b0/d_sr

!c
!c ... hardening function h_alpha	  
!c
      psi=psi_void_DM(void,p,parms,nparms)
      call lode_DM(tau,cM,cos3t,gth,dgdth)
      call alpha_th_DM(2,n,gth,psi,parms,nparms,alpha_b)
      db=distance_sanisand(alpha_b,alpha,n)
!c
      do i=1,6
        h_alpha(i)=twothird*hh*(alpha_b(i)-alpha(i))
      end do
!c
!c ... plastic moduli HHp and Kp
!c
      HHp=twothird*hh*p*db

	if(HHp.gt.1e+15) then
!c		write(*,*)'Hplas'

	endif
!c
      LDeR=dot_vect(1,LL1,U,6)
!c
      Kp=LDeR+HHp

!c .......................................................................
	
	if(mario_DT_test.eq.zero) then

		if(LDeR.lt.zero) then
		switch2=1
!c		write(*,*)'LDeR < 0'
		return
		endif

	

		if(Kp.lt.zero) then
		switch2=1
!c		write(*,*)'function plast_mod: Kp < zero'
		return
		endif
	
	else
		

		if(LDeR.le.zero) then
			switch2=1
!c		error=3
!c		write(*,*)'subroutine plast_mod_DM: LDeR < 0'
		return
		endif
	endif



	if(Kp.lt.zero)then
!c		write(6,*)'function plast_mod: Kp < zero'
	error=3
	return
	endif






	call push(alpha_sr,z,6)



      Kpm1=one/Kp
!c	
      return
      end
!c
!c------------------------------------------------------------------------------
      double precision function psi_void_DM(void,p,parms,nparms)
!c------------------------------------------------------------------------------
!c computes state parameter psi (pyknotropy factor)
!c Dafalias & Manzari (2004) SANISAND model for sands
!c------------------------------------------------------------------------------
      implicit none
!c
      integer nparms
!c
      double precision void,p,parms(nparms)
      double precision p_a,e0,lambda,xi,ec
!c
!c ... recover material parameters
!c
      p_a=parms(1)
      e0=parms(2)
      lambda=parms(3)
	xi=parms(4)
!c
      ec=e0-lambda*(p/p_a)**xi
      psi_void_DM=void-ec
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine push(a,b,n)
!c-----------------------------------------------------------------------------
!c push vector a into vector b
!c-----------------------------------------------------------------------------
      implicit none
      integer i,n
      double precision a(n),b(n) 
!c
      do i=1,n
        b(i)=a(i)
      enddo
!c
      return
      end
!c
      subroutine pzero(v,nn)
!c
!c-----[--.----+----.----+----.-----------------------------------------]
!c      Purpose: Zero real array of data

!c      Inputs:
!c         nn     - Length of array

!c      Outputs:
!c         v(*)   - Array with zero values
!c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer n,nn
      double precision v(nn)

      do n = 1,nn
        v(n) = 0.0d0
      end do ! n

      end
!c
!c-----------------------------------------------------------------------------
      subroutine rkf23_upd_DM(y,z,n,nasvy,nasvz,err_tol,maxnint,DTmin,&
                deps_np1,parms,nparms,nfev,elprsw,&
     		 mario_DT_test,&
                error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c  Dafalias & Manzari(2004) SANISAND model for sand
!c
!c  numerical solution of y'=f(y)
!c  explicit, adapive RKF23 scheme with local time step extrapolation
!c
!c  variables allocated in vector y(n):
!c
!c	y(1)	= sig(1)
!c	y(2)	= sig(2)
!c	y(3)	= sig(3)
!c	y(4)	= sig(4)
!c	y(5)	= sig(5)
!c	y(6)	= sig(6)
!c	y(7)	= alpha(1)
!c	y(8)	= alpha(2)
!c	y(9)	= alpha(3)
!c	y(10)	= alpha(4)
!c	y(11)	= alpha(5)
!c	y(12)	= alpha(6)
!c   y(13)   = void
!c   y(14)   = Fab(1)
!c   y(15)   = Fab(2)
!c   y(16)   = Fab(3)
!c   y(17)   = Fab(4)
!c   y(18)   = Fab(5)
!c   y(19)   = Fab(6)
!c   y(20)   = not used
!c
!c  variables allocated in vector z(nasvz):
!c
!c	z(1)	= alpha_sr(1)
!c	z(2)	= alpha_sr(2)
!c	z(3)	= alpha_sr(3)
!c	z(4)	= alpha_sr(4)
!
!c	z(5)	= alpha_sr(5)
!c	z(6)	= alpha_sr(6)
!c	z(7)	= not used
!c	z(8)	= not used
!c	z(9)	= not used
!c	z(10)	= not used
!c	z(11)	= not used
!c	z(12)	= not used
!c	z(13)	= not used
!c	z(14)	= not used
!c
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer elprsw,mario,switch2,switch3,mario2
	integer mario_DT, mario_DT_test
!c
      integer n,nasvy,nasvz,nparms,i,ksubst,kreject,nfev
      integer maxnint,attempt,maxnint_1
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
!c
      double precision y(n),z(nasvz),deps_np1(6)
      double precision parms(nparms),DTmin,err_tol,err_tol_1, err_tol_n
      double precision zero,half,one,two,three,four,six
      double precision ptnine,one6,one3,two3,temp,prod,pt1
	double precision z1(nasvz),deps_np1_star(6), z_k(nasvz)
!c
      double precision y_k(n),y_tr(n),y_star(n),y_k1(n)!yf_DM,
      double precision y_2(n),y_3(n),y_til(n),y_hat(n)
      double precision p_atm,tol_ff,ff_tr,ff_k
      double precision T_k,DT_k,xi
      double precision kRK_1(n),kRK_2(n),kRK_3(n)
      double precision norm_R,S_hull
      double precision Fab(6),dev_fab(5),I1,f_p,absfp2
      double precision pp,onethird,ptone,p_thres2,tol_ff1,pp_k,pp_tr
	double precision ff_k_pp_k,ff_tr_pp_tr,pp_3,pp_2,pp_hat,ten,min_y_tr
	double precision iter, pp_kk
!c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(four=4.0d0,six=6.0d0,half=0.5d0,ptnine=0.9d0)
      parameter(pt1=1.0d-3,ptone=0.1d0,ten=1.0d1)
!c      parameter(drcor=1)
!c
!c      common /z_nct_errcode/error
!c      common /z_plastic_flag/plastic
!c	common /z_tolerance/tol_f
!c	common /z_check_yield/check_ff
!c	common /z_drift_correction/drcor
!c	common /z_threshold_pressure/p_thres
!c
!c
!c ... initialize y_k vector and other variables
!c
      call pzero(y_k,n)
!c
      one6=one/six
      one3=one/three
      two3=two/three
!c
      plastic=0
	mario = 0
	mario_DT=0
	mario_DT_test=0
    iter = 0

!c
	iter=iter+1
!ccccc	write(*,*)'iter =', iter
!ccccc	if(iter.gt.145) then
!ccccc		write(*,*)'stop'
!ccccc	endif
!c
!c ... start of update process
!c
      call push(y,y_k,n)
	call push(z,z_k,nasvz)
!c
!c ... set tolerance for yield function
!c
      p_atm=parms(1)
      tol_ff=tol_f*p_atm
!c
!c ... yield function at current state, T = 0.0
!c
!c	if (check_ff) then
!c
	ff_k=yf_DM(y_k,n,parms,nparms)
	onethird=one/three
	pp_k=(y_k(1)+y_k(2)+y_k(3))*onethird
	
    if (pp_k<0.0001) then 
        ff_k_pp_k=0
    else 
        ff_k_pp_k=ff_k/pp_k
    end if
        
	if(pp_k.gt.one) ff_k_pp_k=ff_k	

!c
!c ... abort execution if initial state is outside the YS
!c
!c	if (check_ff) then
!c
	if (ff_k_pp_k.gt.tol_ff) then 
!cccc	if (ff_k.gt.tol_ff) then

!c		write(6,*) 'ERROR: initial state is outside the YS'
!c		write(6,*) 'Subroutine RKF23_UPDATE'
!c		write(*,*) 'f = ',ff_k_pp_k
!c		write(*,*) 'f = ',ff_k
!c 		error=10
!c	        error=3
!c 		return 
!c
                call drift_corr_DM(y_k,n,z1,nasvz,parms,nparms,tol_ff,&
                 switch2,mario_DT_test,&
                 error,tol_f,check_ff,drcor,p_thres,plastic)

	end if
!c	endif
!c
!c ... compute trial solution (single step) and update z
!c
	deps_np1_star=deps_np1
	call trial_state(y_k,n,parms,nparms,deps_np1_star,y_tr,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
	pp_tr=(y_tr(1)+y_tr(2)+y_tr(3))*onethird
	if((pp_k.gt.(p_thres+p_thres))) then
		do while(pp_tr.le.p_thres)
			call trial_state(y_k,n,parms,nparms,deps_np1_star,y_tr,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
			pp_tr=(y_tr(1)+y_tr(2)+y_tr(3))*onethird 
			deps_np1_star=deps_np1_star*half
!c			write(*,*) 'pp_k>p_thres'
		end do
	elseif((pp_k.le.(p_thres+p_thres))&
     .and.(pp_tr.gt.(p_thres+p_thres))) then
		deps_np1_star=deps_np1
		call trial_state(y_k,n,parms,nparms,deps_np1_star,y_tr,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
		pp_tr=(y_tr(1)+y_tr(2)+y_tr(3))*onethird
!c		write(*,*) 'pp_k<=p_thres and p_tr>p_thres'
	elseif((pp_k.le.(p_thres+p_thres))&
     .and.(pp_tr.le.(p_thres+p_thres))) then
		 call push(y_k,y_tr,n)
		write(*,*) 'pp_k<=p_thres and pp_tr<=p_thres'
	endif

		
		
!c
!c ... yield function at trial state
!c        
      ff_tr=yf_DM(y_tr,n,parms,nparms)
	pp_tr=(y_tr(1)+y_tr(2)+y_tr(3))*onethird
	
    if (pp_tr <0.001) then 
        ff_tr_pp_tr = 0
    else 
        ff_tr_pp_tr=ff_tr/pp_tr
    end if

    if(pp_tr.gt.one) ff_tr_pp_tr=ff_tr
!c
!c ... compute scalar product of dsig_tr and grad(f) 
!c     to check if crossing of yl occurs
!c
      call check_crossing(y_k,y_tr,n,parms,nparms,prod)
!c        
!c ... check whether plastic loading, elastic unloading or mixed ep loading
!c
      if (ff_tr_pp_tr.lt.tol_ff) then
!cccc		if (ff_tr.lt.tol_ff) then
!c
!c ... Case 1: Elastic unloading inside the yield function: 
!c             trial state is the final state
!c
        call push(y_tr,y_k,n)
!c        
      else
!c
!c ... Case 2: Some plastic loading occurs during the step
!c
!c ................................................................................................................
!c
	if(pp_tr.lt.p_thres) then
!c	(situation not admissible)
!c
!ccc		p_thres2=ptone*p_thres
!c
		write(*,*) 'Low mean pressure, p=',pp_k
!c
!c
!c........................................................................................................
!c
!c
	else
!c
        	if ((ff_k_pp_k.lt.(-tol_ff)).or.(prod.lt.zero)) then
!cccc			if ((ff_k.lt.(-tol_ff)).or.(prod.lt.zero)) then
!c
!c ... Case 2a: the initial part of the stress path is inside the YL;
!c              find the intersection point and update current state at the YL
		call intersect_DM(y_k,y_tr,y_star,n,parms,nparms,tol_ff,xi,&
             error,tol_f,check_ff,drcor,p_thres,plastic)
          		call push(y_star,y_k,n)
!c
!c		Chiara 27.08.09 beginning
!c
			plastic=1
!c
!c		Chiara 27.08.09 end    
!c        
        	else
!c
!c Case 2b: all the stress path lies outside the YL
!c
        	xi=zero
!c
!c		Chiara 27.08.09 beginning
!c
		plastic=1
!c
!c		Chiara 27.08.09 end    
!c
        	end if
!c
!c
!c initialize normalized time and normalized time step
!c
!c
        T_k=xi     
        DT_k=(one-xi)
        ksubst=0
        kreject=0
        nfev=0 
	  attempt=1
	  maxnint_1=maxnint
	  err_tol_1=err_tol 
	  err_tol_n=err_tol  
	  switch3=0           
!c
!c ... start substepping 
!c
!ccc	mario = zero
!c
        do while((T_k.lt.one).and.(mario.eq.zero)&
     .and.(mario_DT.eq.zero)) !**********************************
!c
!ccc.and.(attempt.ne.3)
          ksubst=ksubst+1
!c
!c ... write substepping info
!c
!cccc		if(iter.gt.145)then
!c			if(ksubst.gt.82830) then
!c          write(*,*) ksubst,T_k,DT_k,pp_hat
!c	endif
!cccc	endif
!c
!c ... check for maximum number of substeps
!c
!c 
          if((ksubst.gt.maxnint_1).or.(switch3.eq.1)) then
!c
			if(attempt.eq.1) then
					maxnint_1=2.0*maxnint
					err_tol_1=1000.0*err_tol
		write(*,*) 'I had to increase tolintT,', 'T_k=', T_k
!c					call push(y_hat,y_k,n)
!c					call push(z1,z,nasvz)
!c					call push(y_k,y,n)
     	   				attempt=2
					DT_k=1-T_k
!c			return
!c				
			elseif(attempt.eq.2) then
!ccc					write(*,*) 'number of substeps ',ksubst
!ccc          		 	write(*,*) 'is too big, step rejected'
	      				write(*,*) 'attempt number ',attempt
!ccc					write(*,*) 'T_k =',T_k		
!ccccc					error=3
!ccccc					write(*,*) 'mario1 = one'
					mario=one
					call push(z1,z,nasvz)
					call push(y_k,y,n)
!ccccc					attempt=3
                    return
			endif
		endif
!c
  
				
!c
!c ... build RK functions
!c
		call push(z_k,z1,nasvz)
!cc		call push(y_k,y_k1,nasvy)
		pp_kk=(y_k(1)+y_k(2)+y_k(3))*onethird

!c
          call f_plas_DM(y_k,n,nasvy,z1,nasvz,parms,nparms,&
                        deps_np1,kRK_1,nfev,switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
          if(error.eq.10) return

		if (switch2.eq.zero) then
!c
!c ... find y_2
!c
			temp=half*DT_k
!c
			do i=1,n
			y_2(i)=y_k(i)+temp*kRK_1(i)
			end do
!c
			pp_2=(y_2(1)+y_2(2)+y_2(3))*onethird
!c
			if(pp_2.gt.zero)then
!cc	if((y_2(1).gt.zero).and.(y_2(2).gt.zero).and.(y_2(3).gt.zero))then
!c		
			call f_plas_DM(y_2,n,nasvy,z1,nasvz,parms,nparms,&
     		               deps_np1,kRK_2,nfev,switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
			if(error.eq.10) return

			if (switch2.eq.zero) then
!c					
!c ... find y_3
!c
				do i=1,n
				y_3(i)=y_k(i)-DT_k*kRK_1(i)+two*DT_k*kRK_2(i)
				end do
!c
				pp_3=(y_3(1)+y_3(2)+y_3(3))*onethird
!c
				if(pp_3.gt.zero)then 
!cc		if((y_3(1).gt.zero).and.(y_3(2).gt.zero).and.(y_3(3).gt.zero))then
!c
					call f_plas_DM(y_3,n,nasvy,z1,nasvz,parms,nparms,&
                        deps_np1,kRK_3,nfev,switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
					if(error.eq.10) return
					if (switch2.eq.zero) then
!c		 		
!c ... approx. solutions of 2nd (y_til) and 3rd (y_hat) order
!c
						do i=1,n	
						y_til(i)=y_k(i)+DT_k*kRK_2(i)
						y_hat(i)=y_k(i)+DT_k*&
                    (one6*kRK_1(i)+two3*kRK_2(i)+one6*kRK_3(i))
						end do
!c
!c ... local error estimate
!c
						call norm_res_DM(y_til,y_hat,n,norm_R)
!c
!c ... time step size estimator according to Hull
!c	
						if (norm_R .eq. 0) then 
                            S_hull = 1e-5
                        else 
                            S_hull=ptnine*DT_k*(err_tol/norm_R)**one3
                        end if
                        
!c
						if(norm_R.eq.zero) then
!ccc							write(*,*)'RKF23: norm_R=0'
!ccc							error=3
						endif
!c
!c
!c
!c
      if ((norm_R.lt.err_tol).and.(attempt.ne.2).and.&
           (attempt.ne.3)) then 				
!c
!c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
!c
            

			pp_hat=(y_hat(1)+y_hat(2)+y_hat(1))*onethird
			if (pp_hat.lt.p_thres) then
				write(*,*) 'mario_2 = one'
!cccccc				error =3
				mario=one
			endif
!c
!            
!c ... correct drift from yield surface using Sloan algorithm
!c
				if(drcor.ne.0) then
	       call drift_corr_DM(y_hat,n,z1,nasvz,parms,nparms,tol_ff,&
                   switch2,mario_DT_test,&
                   error,tol_f,check_ff,drcor,p_thres,plastic)
				end if
				pp_hat=(y_hat(1)+y_hat(2)+y_hat(1))*onethird
			if (switch2.eq.zero) then
				call push(y_hat,y_k,n)
				call push(z1,z_k,nasvz)

				T_k=T_k+DT_k;
				DT_k=min(four*DT_k,S_hull)
				DT_k=min((one-T_k),DT_k)
			endif	
!c
!c
			end if
!c
!c
		if ((norm_R.lt.err_tol_1).and.(attempt.eq.2)&
        .and.(switch2.eq.zero)) then 				
!c
!c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
!c
            

			pp_hat=(y_hat(1)+y_hat(2)+y_hat(1))*onethird
			if (pp_hat.lt.p_thres) then
				write(*,*) 'mario_3 = one'
				mario=one
!cccccc				error =3
			endif
!c
            
!c ... correct drift from yield surface using Sloan algorithm
!c
!cc			if (switch2.eq.zero) then
				if(drcor.ne.0) then
	       call drift_corr_DM(y_hat,n,z1,nasvz,parms,nparms,tol_ff,&
      switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
				end if

			if (switch2.eq.zero) then
				call push(y_hat,y_k,n)
				call push(z1,z_k,nasvz)

				T_k=T_k+DT_k;
				DT_k=min(four*DT_k,S_hull)
				DT_k=min((one-T_k),DT_k)
			endif !switch2.eq.zero
!cc			endif !switch2.eq.zero	
!c
!c
		end if !(norm_R.lt.err_tol_1).and.(attempt.eq.2).and.(switch2.eq.zero)
!c
!c
!c
			if((norm_R.gt.err_tol).and.(attempt.ne.3)&
      .and.(switch2.eq.zero)) then
!c
!c ... substep is not accepted, recompute with new (smaller) substep size DT
!c
!c		  write(*,*) 'substep size ',DT_k,'  norm_R =',norm_R
!c
            	DT_k=max(DT_k/four,S_hull)
!c
!c ... check for minimum step size
!c
            		if(DT_k.lt.DTmin) then
             			write(*,*) 'substep size ',DT_k,&
     		         		' is too small, step rejected'
						DT_k= one - T_k
						mario2=1
!cc						err_tol_n=100.0*err_tol
						write(*,*) 'err_tol_n=100*err_tol=',err_tol_n
						switch3=1
!c						write(*,*) 'mario_4 = one'
!c						error =3
            		end if
!c					
      end if
!c
      mario2 = 0
		if((norm_R.lt.err_tol_n).and.(attempt.ne.3)&
     	.and.(switch2.eq.zero).and.(mario2.eq.one)) then
!c
!c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
!c
            

			pp_hat=(y_hat(1)+y_hat(2)+y_hat(1))*onethird
			if (pp_hat.lt.p_thres) then
				write(*,*) 'mario_31 = one; pp_hat<p_thres'
				mario=one
!cccccc				error=3
			endif
!c
            
!c ... correct drift from yield surface using Sloan algorithm
!c
!cc			if (switch2.eq.zero) then
				if(drcor.ne.0) then
	       call drift_corr_DM(y_hat,n,z1,nasvz,parms,nparms,tol_ff,&
      switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
				end if

			if (switch2.eq.zero) then
				call push(y_hat,y_k,n)
				call push(z1,z_k,nasvz)

				T_k=T_k+DT_k;
				DT_k=min(four*DT_k,S_hull)
				DT_k=min((one-T_k),DT_k)
			endif !switch2.eq.zero
!cc			endif !switch2.eq.zero
		mario2=zero	
!cc					
         end if
!c
!c
			if (attempt.eq.3) then
!c
            		write(*,*) 'Third attempt: solution accepted'
		  		
!c
!c ... correct drift from yield surface using Sloan algorithm
!c
		  	if(drcor.ne.0) then
				call drift_corr_DM(y_k,n,z_k,nasvz,parms,nparms,tol_ff,&
      switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
	      		endif

			call push(y_hat,y_k,n)
!c
            		T_k=T_k+DT_k;
!c
          	end if


!c
!c
!c ... if switch2=1 (if switch2=1 occurs in accepted solution part)
			if (switch2.ne.zero) then
				DT_k=DT_k/four
	            		if(DT_k.lt.DTmin) then
						DT_k= one - T_k
!ccccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - switch2>0, mario5 = one'
!ccccccccc						error =3
!ccccccccc						return
!cccccc						attempt=3
            		end if	
			endif
			
	else
!c ... if switch2=1 (in kRk_3)
		DT_k=DT_k/four
           		if(DT_k.lt.DTmin) then
							DT_k= one - T_k
!ccccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - in kRk_3, mario6 = one, pp_3=',pp_3
!ccccccccc						error =3
!ccccccccc						return
!cccccc						attempt=3
            		end if	
	endif
		
		else								
!cc ... if pp_3 lower than zero
            	DT_k=DT_k/four
				if(DT_k.lt.DTmin) then
					DT_k= one - T_k
!ccccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - pp_3<0, mario7=one, pp_2=',pp_2
!ccccccccc						error =3
!ccccccccc						return
!cccccc						attempt=3
            			end if
		endif

	else
!c ... if switch2=1 (in kRk_2)
		DT_k=DT_k/four
				if(DT_k.lt.DTmin) then
					DT_k= one - T_k
!cccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - kRk_2, mario8=one, pp_2=',pp_2
!cccccccc						error =3
!cccccccc						return
!cccccc						attempt=3
            			end if
	endif


	else
!c ... if pp_2 lower than zero
		DT_k=DT_k/four
			if(DT_k.lt.DTmin) then
					DT_k= one - T_k
!ccccccc							mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - pp_2<0 - mario9=one, pp_kk=',pp_kk
!cccccccc						error =3
!cccccccc						return
!cccccc						attempt=3
            			end if
	endif
!c
	else
!c ... if switch2=1 (in kRk_1)
		DT_k=DT_k/four
			if(DT_k.lt.DTmin) then
					DT_k= one - T_k
!cccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - kRk_1, mario10=one, pp_kk=',pp_kk
!c
!c
!ccccccccc						error =3
!ccccccccc							return
!cccccc						attempt=3
            			end if
	endif

!c
!c
!c
!c ... bottom of while loop
!c
        end do !*****************************************************
!c
       end if 
      end if
!c
	if(mario.eq.1) then !stop solution, keep current configuration
!cccccc					DT_k=1-T_k
					write(*,*) 'attempt number ',attempt
      				write(*,*) 'accept solution'
					call push(z_k,z,nasvz)
					call push(y_k,y,n)
!c	error=3
!ccc	attempt=1
	else if(mario_DT.eq.1) then	 !abort solution, keep previous configuration
					write(*,*) 'mario_DT'
!c					call push(z1,z,nasvz)
!c					call push(y_k,y,n)
	else		
!c
      call push(y_k,y,n)
	call push(z_k,z,nasvz)
!c
	endif
!c
		if(drcor.ne.0) then
	call drift_corr_DM(y,n,z,nasvz,parms,nparms,tol_ff,&
      switch2,mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
	    endif
		if (switch2.ne.zero) then
			write(*,*) 'RKF23 - drift_corr - switch2>0'
!cccccc			error=3
			return
		endif 
!c
      return
!c
 2000 format('Substep no.',i4,'- T_k=',d12.4,'- DT_k=',d12.4,&
     ' -pp_hat=',d12.4)
!c
      end
!c
!-----------------------------------------------------------------------------
      subroutine solout(stress,ntens,asv1,nasvy,asv2,nasvz,ddsdde,&
                       y,nydim,z,pore,depsv_np1,parms,nparms,DD)
!-----------------------------------------------------------------------------
! copy the vector of state variables to umat output
! Dafalias &Manzari SANISAND model (2004)
!
! modified 4/2008 (Tamagnini)
!
! NOTE: solid mechanics convention for stress and strain components
!       pore is always positive in compression
!
! depsv_np1 = vol. strain increment, compression positive 
! y(1:6)    = effective stress, compression positive 
! pore      = excess pore pressure, compression positive 
! stress    = total stress, tension positive
!
!-----------------------------------------------------------------------------
      implicit none
!
      integer nydim,nasvy,nasvz,nparms,ntens,i,j
!
      double precision y(nydim),z(nasvz),asv1(nasvy),asv2(nasvz)
      double precision stress(ntens),ddsdde(ntens,ntens),DD(6,6)
      double precision parms(nparms),bulk_w,pore,depsv_np1 
!
      bulk_w=parms(17)
!
! ... update excess pore pressure (if undrained cond.), compression positive
!
      pore=pore+bulk_w*depsv_np1
!
! updated total stresses (effective stresses stored in y(1:6))
!
      do i=1,ntens
		if (i.le.3) then
          stress(i) = -y(i)-pore
		else
          stress(i) = -y(i)
        end if
      enddo
!
! additional state variables
!
      do i=1,nasvy
		asv1(i) = y(6+i)
      enddo
!
      do i=1,nasvz
		asv2(i) = z(i)
      enddo
!
! consistent tangent stiffness
!
      do j=1,ntens
        do i=1,ntens
		if((i.le.3).and.(j.le.3)) then
          ddsdde(i,j) = DD(i,j)+bulk_w
		else
          ddsdde(i,j) = DD(i,j)
		end if        
        end do
      enddo
!
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine tang_stiff(y,z,n,nasvy,nasvz,parms,nparms,DD,cons_lin,&
                error,tol_f,check_ff,drcor,p_thres,plastic)
!c-----------------------------------------------------------------------------
!c  compute continuum tangent stiffness at the end of the step
!c  Dafalias & Manzari (2004) SANISAND model for sand
!c
!c  written 10/2008 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
!c 
      integer switch2,mario_DT_test
!c
      integer n,nasvy,nasvz,nparms,cons_lin
!c
      double precision y(n),z(nasvz),parms(nparms)
      double precision DD(6,6),HH(nasvy,6)
      double precision zero,three
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres

      parameter(zero=0.0d0,three=3.0d0)
!c
!c      common /z_nct_errcode/error
!c      common /z_plastic_flag/plastic
!c
!c ... initialize DD
!c 
      call pzero(DD,36)
!c
      if(plastic.eq.1 .and. cons_lin.eq.1) then
!c
!c ... plastic process
!c
	  call get_tan_DM(y,n,nasvy,z,nasvz,parms,nparms,DD,HH,switch2,&
         mario_DT_test,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
     
       else       
             call el_stiff_DM(y,n,parms,nparms,DD,&
              error,tol_f,check_ff,drcor,p_thres,plastic)

!c	
!c ... if switch2=0 (LDeR<0) try Elastic tangent matrix 
!c
! 	if(switch2.ne.zero) then
!c        error=3
! 		write(*,*) 'tang_stiff - switch2=1 - elastoplastic_matrix=0'
! 	endif
!c
      end if	
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine trial_state(y_k,n,parms,nparms,deps,y_tr,&
              error,tol_f,check_ff,drcor,p_thres,plastic)               
!c-----------------------------------------------------------------------------
!c
!c ... computes the trial stress state (freezing plastic flow) 
!c
!c     eps  = strain at the beginning of the step
!c     deps = strain increment for this step
!c
!c	NOTE:
!c
!c	mode = 1	single step Runge-Kutta 3rd order
!c	mode = 2	single step forward Euler
!c
!c-----------------------------------------------------------------------------
!c
      implicit none
!c
      integer n,m,nparms,mode,i
!c
      double precision y_k(n),parms(nparms),deps(6)
      double precision y_2(n),y_3(n),y_tr(n)
      double precision one,two,three,six
      double precision kRK_1(n),kRK_2(n),kRK_3(n)
      double precision DT_k,DTk05,DTk2,DTk6,DTk23
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
!c
      parameter(mode=1)
!c
      data one,two,three,six/1.0d0,2.0d0,3.0d0,6.0d0/
!c
      DT_k=one
!c
!c ... compute F_el
!c
      call f_hypoelas_DM(y_k,n,parms,nparms,deps,kRK_1,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
      if (mode.eq.1) then
!c
!c ... 3rd order RK - build F function
!c
        DTk05=DT_k/two
        DTk2=two*DT_k
        DTk6=DT_k/six
        DTk23=two*DT_k/three
!c
        do i=1,n
          y_2(i)=y_k(i)+DTk05*kRK_1(i)
        end do ! i
!c
        call f_hypoelas_DM(y_2,n,parms,nparms,deps,kRK_2,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
        do i=1,n
          y_3(i)=y_k(i)-DT_k*kRK_1(i)+DTk2*kRK_2(i)
        end do ! i
!c
        call f_hypoelas_DM(y_3,n,parms,nparms,deps,kRK_3,&
              error,tol_f,check_ff,drcor,p_thres,plastic)
!c
        do i=1,n
          y_tr(i)=y_k(i)+DTk6*kRK_1(i)+DTk23*kRK_2(i)+DTk6*kRK_3(i)
        end do ! i
!c
      else
!c
!c ... forward Euler
!c
        do i=1,n
          y_tr(i)=y_k(i)+DT_k*kRK_1(i)
        end do ! i
!c
      end if
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine wrista(mode,y,nydim,deps_np1,dtime,coords,statev,&
                nstatv,parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
!c-----------------------------------------------------------------------------
!c ... subroutine for managing output messages
!c
!c     mode
!c
!c     all = writes: kstep, kinc, noel, npt
!c     2   = writes also: error message,coords(3),parms(nparms),ndi,
!c           nshr,stress(nstress),deps(nstress),dtime,statev(nstatv)
!c     3   = writes also: stress(nstress),deps(nstress),dtime,statev(nstatv)
!c
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer mode,nydim,nstatv,nparms,noel,npt,ndi,nshr,kstep,kinc,i    
!c
      double precision y(nydim),statev(nstatv),parms(nparms)
      double precision deps_np1(6),coords(3),dtime
!c
!c ... writes for mode = 2
!c
      if (mode.eq.2) then
        write(6,*) '==================================================='
        write(6,*) '   ERROR: abaqus job failed during call of UMAT'
        write(6,*) '==================================================='
        write(6,*) ' state dump: '
        write(6,*) 
      endif
!c
!c ... writes for all mode values
!c
      !commented by Abdelrahman
	!if(mode.ne.4) then
 !     write(6,111) 'Step: ',kstep, 'increment: ',kinc,&
 !     'element: ', noel, 'Integration point: ',npt
 !     write(6,*) 
	!endif
!c
!c ... writes for mode = 2
!c
      !commented by Abdelrahman
      if (mode.eq.2) then
        write(6,*) 'Co-ordinates of material point:'
        write(6,104) 'x1 = ',coords(1),' x2 = ',coords(2),' x3 = ',&
         coords(3)
        write(6,*) 
        write(6,*) 'Material parameters:'
        write(6,*) 
        do i=1,nparms
          write(6,105) 'prop(',i,') = ',parms(i)
        enddo 
        write(6,*)
        write(6,102) 'No. of mean components:  ',ndi
        write(6,102) 'No. of shear components: ',nshr
        write(6,*)
      endif
!c
!c ... writes for mode = 2 or 3
!c
      if ((mode.eq.2).or.(mode.eq.3)) then
        write(6,*) 'Stresses:'
        write(6,*) 
        write(6,101) 'sigma(1) = ',y(1)
        write(6,101) 'sigma(2) = ',y(2)
        write(6,101) 'sigma(3) = ',y(3)
        write(6,101) 'sigma(4) = ',y(4)
        write(6,101) 'sigma(5) = ',y(5)
        write(6,101) 'sigma(6) = ',y(6)
        write(6,*) 
        write(6,*) 'Strain increment:'
        write(6,*) 
        write(6,101) 'deps_np1(1) = ',deps_np1(1)
        write(6,101) 'deps_np1(2) = ',deps_np1(2)
        write(6,101) 'deps_np1(3) = ',deps_np1(3)
        write(6,101) 'deps_np1(4) = ',deps_np1(4)
        write(6,101) 'deps_np1(5) = ',deps_np1(5)
        write(6,101) 'deps_np1(6) = ',deps_np1(6)
        write(6,*) 
        write(6,*) 'Time increment:'
        write(6,*) 
        write(6,108) 'dtime = ',dtime
        write(6,*) 
        write(6,*) 'Internal state variables:'
        write(6,*) 
        write(6,109) 'alpha_11     = ',statev(1)
        write(6,109) 'alpha_22     = ',statev(2)
        write(6,109) 'alpha_33     = ',statev(3)
        write(6,109) 'alpha_12     = ',statev(4)
        write(6,109) 'alpha_13     = ',statev(5)
        write(6,109) 'alpha_23     = ',statev(6)
        write(6,109) 'void         = ',statev(7)
        write(6,109) 'Fab_11       = ',statev(8)
        write(6,109) 'Fab_22       = ',statev(9)
        write(6,109) 'Fab_33       = ',statev(10)
        write(6,109) 'Fab_12       = ',statev(11)
        write(6,109) 'Fab_13       = ',statev(12)
        write(6,109) 'Fab_23       = ',statev(13)
        write(6,109) 'dummy        = ',statev(14)
        write(6,109) 'alpha_sr_11  = ',statev(15)
        write(6,109) 'alpha_sr_22  = ',statev(16)
        write(6,109) 'alpha_sr_33  = ',statev(17)
        write(6,109) 'alpha_sr_12  = ',statev(18)
        write(6,109) 'alpha_sr_13  = ',statev(19)
        write(6,109) 'alpha_sr_23  = ',statev(20)
        write(6,109) 'dummy        = ',statev(21)
        write(6,109) 'dummy        = ',statev(22)
        write(6,109) 'dummy        = ',statev(23)
        write(6,109) 'dummy        = ',statev(24)
        write(6,109) 'dummy        = ',statev(25)
        write(6,109) 'dummy        = ',statev(26)
        write(6,109) 'dummy        = ',statev(27)
        write(6,109) 'dummy        = ',statev(28)
        write(6,*) 
        write(6,*) '==================================================='
!c
      endif

	if (mode.eq.4) then
	      write(6,111) 'Step: ',kstep, 'increment: ',kinc,&
      'element: ', noel, 'Integration point: ',npt
        write(6,104) 'DT<DTmin:x1=',coords(1),' x2 = ',&
         coords(2),' x3 = ',coords(3)
!c	      write(6,*) 
	endif
!c
101   format(1X,a10,f50.44)
102   format(1X,a25,i1)
103   format(1X,a7,i5)
104   format(1X,3(a12,f10.4,2X))
105   format(1X,a5,i2,a4,f20.3)
106   format(1X,3(a9,f12.4,2X))
107   format(1X,3(a10,f12.4,2X))
108   format(1X,a8,f12.4)
109   format(1X,a10,f50.44)
110   format(1X,a5,f10.4)
111   format(1X,a6,i4,2X,a11,i4,2X,a9,i10,2X,a19,i4)
!c       
      return
      end
!c
!c-----------------------------------------------------------------------------
      double precision function yf_DM(y,ny,parms,nparms)
!c-----------------------------------------------------------------------------
!c     compute yield function of Dafalias & Manzari (2004) 
!c     SANISAND model for sands
!c-----------------------------------------------------------------------------
!c
      implicit none
!c
      integer ny,nparms
!c
      !double precision dot_vect
      double precision y(ny),parms(nparms)
      double precision mm,zero,one,two,three,sqrt23,norm2
      double precision sig(6),s(6),trace,p,alpha(6),sbar(6)
!c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
!c
!c ... some constants and material parameters
!c
      sqrt23=dsqrt(two/three)
      mm=parms(7)
!c   	
!c ... recover state variables
!c
      sig(1)=y(1)
      sig(2)=y(2)
      sig(3)=y(3)
      sig(4)=y(4)
      sig(5)=y(5)
      sig(6)=y(6)
!c
      alpha(1)=y(7)
      alpha(2)=y(8)
      alpha(3)=y(9)
      alpha(4)=y(10)
      alpha(5)=y(11)
      alpha(6)=y(12)
!c
!c ... mean stress and deviator stress
!c
      call deviator(sig,s,trace,p)
!c
      sbar(1)=s(1)-p*alpha(1)
      sbar(2)=s(2)-p*alpha(2)
      sbar(3)=s(3)-p*alpha(3)
      sbar(4)=s(4)-p*alpha(4)
      sbar(5)=s(5)-p*alpha(5)
      sbar(6)=s(6)-p*alpha(6)
!c
!c ... compute yield function
!c
      norm2=dot_vect(1,sbar,sbar,6)
      yf_DM=dsqrt(norm2)-sqrt23*mm*p
!c
      return
      end
!c
!c.....MODIFICATO
!c
	subroutine xit_DM()
	stop
	return
	end
!c
      
!c------------------------------------------------------------------------------
      subroutine inv_sig_full(sig,pp,qq,cos3t,I1,I2,I3)
!c------------------------------------------------------------------------------
!c calculate invariants of stress tensor
!c
!c NOTE: Voigt notation is used with the following index conversion
!c
!c       11 -> 1
!c       22 -> 2
!c    33 -> 3
!c       12 -> 4
!c       13 -> 5
!c       23 -> 6
!c
!c------------------------------------------------------------------------------
!c
      implicit none
!c
      double precision sig(6),sdev(6)
      double precision eta(6),eta_d(6),eta_d2(6)
      double precision xmin1,xmin2,xmin3
      double precision tretadev3,pp,qq,cos3t,I1,I2,I3
      double precision norm2,norm2sig,norm2eta,numer,denom
!c
      double precision half,one,two,three,six
      double precision onethird,threehalves,sqrt6,tiny
!c
      !double precision dot_vect
!c
      data half,one/0.5d0,1.0d0/
      data two,three,six/2.0d0,3.0d0,6.0d0/
      data tiny/1.0d-18/
!c
!c ... some constants
!c
      onethird=one/three
      threehalves=three/two
      sqrt6=dsqrt(six)
!c
!c ... trace and mean stress
!c
      I1=sig(1)+sig(2)+sig(3)
      pp=onethird*I1
!c
!c ... deviator stress
!c
      sdev(1)=sig(1)-pp
      sdev(2)=sig(2)-pp
      sdev(3)=sig(3)-pp
      sdev(4)=sig(4)
      sdev(5)=sig(5)
      sdev(6)=sig(6)
!c
!c ... normalized stress and dev. normalized stress
!c

      if(I1.ne.0) then
         eta(1)=sig(1)/I1
         eta(2)=sig(2)/I1
         eta(3)=sig(3)/I1
         eta(4)=sig(4)/I1
         eta(5)=sig(5)/I1
        eta(6)=sig(6)/I1
      else
        eta(1)=sig(1)/tiny
        eta(2)=sig(2)/tiny
        eta(3)=sig(3)/tiny
        eta(4)=sig(4)/tiny
        eta(5)=sig(5)/tiny
        eta(6)=sig(6)/tiny        
      end if
!c
      eta_d(1)=eta(1)-onethird
      eta_d(2)=eta(2)-onethird
      eta_d(3)=eta(3)-onethird
      eta_d(4)=eta(4)
      eta_d(5)=eta(5)
      eta_d(6)=eta(6)
!c
!c ... second invariants
!c
      norm2=dot_vect(1,sdev,sdev,6)
      norm2sig=dot_vect(1,sig,sig,6)
      norm2eta=dot_vect(1,eta_d,eta_d,6)
!c
      qq=dsqrt(threehalves*norm2)
      I2=half*(norm2sig-I1*I1)
!c
!c ... components of (eta_d_ij)(eta_d_jk)
!c
      eta_d2(1)=eta_d(1)*eta_d(1)+eta_d(4)*eta_d(4)+eta_d(5)*eta_d(5)
      eta_d2(2)=eta_d(4)*eta_d(4)+eta_d(2)*eta_d(2)+eta_d(6)*eta_d(6)
      eta_d2(3)=eta_d(6)*eta_d(6)+eta_d(5)*eta_d(5)+eta_d(3)*eta_d(3)
      eta_d2(4)=eta_d(1)*eta_d(4)+eta_d(4)*eta_d(2)+eta_d(6)*eta_d(5)
      eta_d2(5)=eta_d(5)*eta_d(1)+eta_d(6)*eta_d(4)+eta_d(3)*eta_d(5)
      eta_d2(6)=eta_d(4)*eta_d(5)+eta_d(2)*eta_d(6)+eta_d(6)*eta_d(3)
!c           
!c ... Lode angle
!c
      if(norm2eta.lt.tiny) then 
!c
        cos3t=-one
!c               
      else
!c
        tretadev3=dot_vect(1,eta_d,eta_d2,6)
!c
        numer=-sqrt6*tretadev3
        denom=(dsqrt(norm2eta))**3
        cos3t=numer/denom
        if(dabs(cos3t).gt.one) then
             cos3t=cos3t/dabs(cos3t)
        end if
!c
      end if 
!c
!c ... determinant
!c
      xmin1=sig(2)*sig(3)-sig(6)*sig(6)
      xmin2=sig(4)*sig(3)-sig(6)*sig(5)
      xmin3=sig(4)*sig(6)-sig(5)*sig(2)
!c
      I3=sig(1)*xmin1-sig(4)*xmin2+sig(5)*xmin3

!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine check_RKF_DM(error_RKF,y,ny,nasv,parms,nparms)
!c-----------------------------------------------------------------------------
!c Checks is RKF23 solout vector y is OK for hypoplasticity
!c-----------------------------------------------------------------------------
      implicit none
!c
        integer error_RKF,ny,nasv,i,nparms,testnan,iopt
!c
        double precision y(ny),parms(nparms)
        double precision sig(6),pmean,sig_star(6)
        double precision I1,I2,I3,pp,qq,cos3t
        double precision ptshift,minstress,sin2phim,tolerance
        double precision OCR,omega,fSBS,sensit,cos2phic
        double precision coparam,sin2phicco
!c
!c	check for NAN
 	testnan=0
        do i=1,ny
       	  call umatisnan_DM(y(i),testnan)
        end do
        if(testnan.eq.1) error_RKF=1
!c
      return
      end
!c
!c-----------------------------------------------------------------------------
      subroutine umatisnan_DM(chcknum,testnan)
!c-----------------------------------------------------------------------------
!c
!c  checks whether number is NaN
!c
!c-----------------------------------------------------------------------------
        double precision chcknum
        integer testnan

	    if (.not.(chcknum .ge. 0. .OR. chcknum .lt. 0.)) testnan=1        
	    if (chcknum .gt. 1.d30) testnan=1        
	    if (chcknum .lt. -1.d30) testnan=1        
 	    if (chcknum .ne. chcknum) testnan=1        
       
        return
        end         


     
     
     
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        !        
        !███╗░░░███╗░█████╗░██╗░░██╗██████╗░░░░░░░░█████╗░░█████╗░██╗░░░██╗██╗░░░░░░█████╗░███╗░░░███╗██████╗░
        !████╗░████║██╔══██╗██║░░██║██╔══██╗░░░░░░██╔══██╗██╔══██╗██║░░░██║██║░░░░░██╔══██╗████╗░████║██╔══██╗
        !██╔████╔██║██║░░██║███████║██████╔╝█████╗██║░░╚═╝██║░░██║██║░░░██║██║░░░░░██║░░██║██╔████╔██║██████╦╝
        !██║╚██╔╝██║██║░░██║██╔══██║██╔══██╗╚════╝██║░░██╗██║░░██║██║░░░██║██║░░░░░██║░░██║██║╚██╔╝██║██╔══██╗
        !██║░╚═╝░██║╚█████╔╝██║░░██║██║░░██║░░░░░░╚█████╔╝╚█████╔╝╚██████╔╝███████╗╚█████╔╝██║░╚═╝░██║██████╦╝
        !╚═╝░░░░░╚═╝░╚════╝░╚═╝░░╚═╝╚═╝░░╚═╝░░░░░░░╚════╝░░╚════╝░░╚═════╝░╚══════╝░╚════╝░╚═╝░░░░░╚═╝╚═════╝░
        
        Subroutine ESM_MohrCoulomb(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
     DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"ESM" :: ESM
      implicit double precision (a-h, o-z) 
      CHARACTER*80 CMNAME     
      DIMENSION STRESS(NTENS),&
     DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
      !NPT(1),NOEL(1),IDSET(1),EUNLOADING(1),PLASTICMULTIPLIER(1),,NUMBEROFPHASES(1)

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
              
  
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

!---Local variables defned by the user
! e.g. integer :: var_local	  
!---User can define here additional variables     
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
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
          call umat_MohrCoulomb(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,& 
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0,&
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)

      
!---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
!---Always define this value to run the simulation
    
! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
    
          


        return

     end subroutine ESM_MohrCoulomb

        
        
        
        
        
        
     
     
     
     
     !*USER SUBROUTINES
      SUBROUTINE UMAT_MohrCoulomb(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
      RPL,DDSDDT,DRPLDE,DRPLDT,&
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      !INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
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

      call PrnSig_MohrCoulomb(1,NTENS,SigE,xN1,xN2,xN3,Sig1,Sig2,Sig3,P,Q)

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
      IPL = 0
      IF (IPL  == 0) GOTO 360 !commented by Abdelrahman Alsardi
!
!         ***POINTS CHANGING FROM PLASTIC TO ELASTIC STATE***
!
      IPL =0
      GOTO 360
!
!                   *** Plastic stress points ***
!
  240 IPL =1


     
      call PrnSig_MohrCoulomb( 0,NTENS,SigEQ,Tmp1,Tmp2,Tmp3,SigV1,SigV2,SigV3,Dum1,Dum2)
   

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
      HA= ( 1 - SPHI + SPSI - SPHI * SPSI) * SIG1 + &
           (-2 - 2/VNUI1*SPHI*SPSI) * SIG2 + &
          ( 1 - SPHI - SPSI + VNUIQ * SPHI * SPSI) * SIG3 + &
          2*(1 + SPSI) * COHS

      HB= ( 1 + SPHI + SPSI + VNUIQ * SPHI * SPSI) * SIG1 -&
          ( 2 + 2/VNUI1*SPHI*SPSI) * SIG2 + &
          ( 1 + SPHI - SPSI - SPHI * SPSI) * SIG3 - &
          2*(1 - SPSI) * COHS

      HAB= (VNUI1+SPSI)*SIG1+ &
           (VNUI1-SPSI)*(SIG3-TENS)- &
           (VNUI1+SPSI)*(TENS*(1+SPHI)-2*COHS)/(1-SPHI)
     
      HBA=(1-VNUI0)*SIG1-VNUI0*SIG3-(1-VNUI0)*(TENS*(1+SPHI)- &
           2*COHS)/(1-SPHI)+VNUI0*TENS

      HAO=(1-VNUI0)*SIG2-VNUI0*SIG3-VNUI1*TENS

      HOC=(1-VNUI0)*SIG1-VNUI0*SIG3-VNUI1*TENS

      HAA=(VNUI1+VNUI2*SPSI)*(SIG1-(TENS*(1+SPHI)-2*COHS)/ &
          (1-SPHI))+ &
          (VNUI1-SPSI)*(SIG2-TENS)+(VNUI1-SPSI)*(SIG3-TENS)

      HAAB=VNUI0*(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))- &
           (SIG2-TENS)+VNUI0*(SIG3-TENS)

      HAAO=(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))- &
           VNUI0*(SIG2-TENS)-VNUI0*(SIG3-TENS)

      HBB=(VNUI1+SPSI)*(SIG1-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))+& 
          (VNUI1+SPSI)*(SIG2-(TENS*(1+SPHI)-2*COHS)/(1-SPHI))+& 
          (VNUI1-VNUI2*SPSI)*(SIG3-TENS)

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

      IF (FT3 > 1d-12 .AND. &
          HBA >= 0 .AND. &
          HAO < 0 .AND. &
          HOC < 0) THEN
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

      IF ( F21R > 1d-6.OR. &
           F32R > 1d-6.OR. &
           F31R > 1d-6.OR. &
           FT1R > 1d-6.OR. &
           FT2R > 1d-6.OR. &
           FT3R > 1d-6    ) THEN
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
     
     
     
     
     
     
     
      Subroutine MZEROR(R,K)
!C
!C***********************************************************************
!C
!C     Function: To make a real array R with dimension K to zero
!C
!C***********************************************************************
!C
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J) = 0.0D0
      End Do

      Return
      End


      Subroutine MZEROI(I,K)
!C
!C***********************************************************************
!C
!C     Function: To make an integre array I with Dimension K to zero
!C
!C***********************************************************************
!C
      Dimension I(*)

      Do J=1,K
        I(J)=0
      End Do

      Return
      End

      Subroutine SETRVAL(R,K,V)
!C
!C***********************************************************************
!C
!C     Function: To fill a real array R with Dimension K with value V
!C
!C***********************************************************************
!C
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J)=V
      End Do

      Return
      End

      Subroutine SETIVAL(I,K,IV)
!C
!C***********************************************************************
!C
!C     Function: To fill an integer array I with Dimension K with value IV
!C
!C***********************************************************************
!C
      Implicit Double Precision (A-H,O-Z)
      Dimension I(*)

      Do J=1,K
        I(J)=IV
      End Do

      Return
      End

      Subroutine COPYIVEC(I1,I2,K)
!C
!C***********************************************************************
!C
!C     Function: To copy an integer array I1 with Dimension K to I2
!C
!C***********************************************************************
!C
      Implicit Double Precision (A-H,O-Z)
      Dimension I1(*),I2(*)

      Do  J=1,K
        I2(J)=I1(J)
      End Do

      Return
      End

      Subroutine COPYRVEC(R1,R2,K)
!C
!C***********************************************************************
!C
!C     Function: To copy a Double array R1 with Dimension K to R2
!C
!C***********************************************************************
!C
      Implicit Double Precision (A-H,O-Z)
      Dimension R1(*),R2(*)

      Do J=1,K
        R2(J)=R1(J)
      End Do

      Return
      End


      Logical Function IS0ARR_MohrCoulomb(A,N)
!C
!C***********************************************************************
!C    Function :  To check whether a real array contains only zero values.
!C                When an array contains only zero's is might not need to be
!C                written to the XXX file.
!C                exit Function when first non-zero value occured or when
!C                all elements are checked and are zero.
!C
!C    Input:  A : array to be checked
!C            N : number of elements in array that should be checked
!C
!C    Output : .TRUE.  when all elements are 0
!C             .FALSE. when at least one element is not zero
!C
!C    Called by :  Subroutine TOBXX
!C
!C***********************************************************************
!C
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*)
      Is0Arr_MohrCoulomb=.False.
      Do I=1,N
        If ( A(I) .Ne. 0 ) Return
      End Do
      Is0Arr_MohrCoulomb=.True.
      Return
      End

      Logical Function IS0IARR_MohrCoulomb(IARR,N)
!C
!C***********************************************************************
!C    Function :  To check whether a integer array contains only zero values.
!C                Similar to IS0ARR
!C
!C    Input:  IARR : array to be checked
!C            N    : number of elements in array that should be checked
!C
!C    Output : .TRUE.  when all elements are 0
!C             .FALSE. when at least one element is not zero
!C
!C    Called by :  Subroutine TOBXX
!C
!C***********************************************************************
!C
      Implicit Double Precision (A-H,O-Z)
      Dimension IARR(*)

      Is0IArr_MohrCoulomb=.False.
      Do I=1,N
        If ( IARR(I) .Ne. 0 ) Return
      End Do
      Is0IArr_MohrCoulomb=.True.
      Return
      End
!C***********************************************************************
      Subroutine MulVec(V,F,K)
!C***********************************************************************
!C
!C     Function: To multiply a real vector V with dimension K by F
!C
!C***********************************************************************
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*)

      Do J=1,K
        V(J)=F*V(J)
      End Do

      Return
      End     ! Subroutine Mulvec
!C***********************************************************************
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
      Implicit Double Precision (A-H,O-Z)
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

!C***********************************************************************
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
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
!C***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine AddVec
!C
!C***********************************************************************
      Double Precision Function DInProd(A,B,N)
!C***********************************************************************
!C
!C     Returns the Inproduct of two vectors
!C
!C I   A,B  : Two vectors
!C I   N    : Used length of vectors
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*),B(*)
!C***********************************************************************

      X = 0
      Do I=1,N
        X = X + A(I)*B(I)
      End Do
      DInProd = X
      Return
      End     ! Function DInProd
!C
!C***********************************************************************
      Subroutine MatMat(xMat1,Id1,xMat2,Id2,nR1,nC2,nC1,xMatR,IdR)
!C***********************************************************************
!C
!C     Calculate xMatR = xMat1*xMat2
!C
!C I   xMat1 : Matrix (Id1,*)
!C I   xMat2 : Matrix (Id2,*)
!C I   nR1   : Number of rows in resulting matrix    (= No rows in xMat1)
!C I   nC2   : Number of columns in resulting matrix (= No cols in xMat2)
!C I   nC1   : Number of columns in matrix xMat1
!C             = Number  rows    in matrix xMat2
!C O   xMatR : Resulting matrix (IdR,*)
!C
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(Id1,*),xMat2(Id2,*),xMatR(IdR,*)
!C**********************************************************************

      Do I=1,nR1
        Do J=1,nC2
          X=0
          Do K=1,nC1
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMat

!C***********************************************************************
      Subroutine MatMatSq(n, xMat1, xMat2, xMatR)
!C***********************************************************************
!C
!C     Calculate xMatR = xMat1*xMat2 for square matrices, size n
!C
!C I   n     : Dimension of matrices
!C I   xMat1 : Matrix (n,*)
!C I   xMat2 : Matrix (n,*)
!C O   xMatR : Resulting matrix (n,*)
!C
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(n,*),xMat2(n,*),xMatR(n,*)
!C**********************************************************************

      Do I=1,n
        Do J=1,n
          X=0
          Do K=1,n
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMatSq

!C***********************************************************************
      Subroutine WriVal ( io, C , V )
!C***********************************************************************
!C
!C Write (Double) value to file unit io (when io>0)
!C
!C***********************************************************************
!C
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,V
    1 Format( A,3x, 1x,1p,e12.5)
      Return
      End
!C***********************************************************************
      Subroutine WriIVl ( io, C , I )
!C***********************************************************************
!C
!C Write (integer) value to file unit io (when io>0)
!C
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,I
    1 Format( A,3x, 1x,I6)
      Return
      End
!C***********************************************************************
      Subroutine WriIVc ( io, C , iV , n )
!C***********************************************************************
!C
!C Write (integer) vector to file unit io (when io>0)
!C
!C***********************************************************************
      Character C*(*)
      Dimension iV(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) (iv(i),i=1,n)
    1 Format( ( 2(3x,5i4) ) )
      Return
      End
!C***********************************************************************
      Subroutine WriVec ( io, C , V , n )
!C***********************************************************************
!C
!C Write (Double) vector to file unit io (when io>0)
!C 6 values per line
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      If (Len_Trim(C).Le.6) Then
        Write(io,2) C,( V(i),i=1,n)
      Else
        Write(io,*) C
        Write(io,1) ( V(i),i=1,n)
      End If
    1 Format( ( 2(1x, 3(1x,1p,e10.3) ) ) )
    2 Format( A, ( T7, 2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
!C***********************************************************************
      Subroutine WriVec5( io, C , V , n )
!C***********************************************************************
!C
!C Write (Double) vector to file unit io (when io>0)
!C 5 values per line
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) ( V(i),i=1,n)
    1 Format( 5(1x,1p,e12.5) )
      Return
      End
!C***********************************************************************
      Subroutine WriMat ( io, C , V , nd, nr, nc )
!C***********************************************************************
!C
!C Write (Double) matrix to file unit io (when io>0)
!C 6 values per line
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(nd,*)

      If (io.Le.0) Return

      Write(io,*) C
      Do j=1,nr
        Write(io,1) j,( V(j,i),i=1,nc)
      End Do
    1 Format(i4, (  T7,2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
!C***********************************************************************

      Subroutine MatInvPiv(Aorig,B,N)
      Implicit Double Precision (A-H,O-Z)
      Dimension Aorig(n,*), B(n,*),A(:,:)
      Allocatable :: A
      Allocate ( A(n,n) ) ! No error checking !!
      Call CopyRVec(AOrig, A, n*n )
      Call MZeroR(B,n*n)
      Do i=1,n
        B(i,i) = 1d0
      End Do
      Do I=1,n
        T=A(I,I)
        iPiv=i
        Do j=i+1,n
          If ( Abs(A(j,i)) .Gt. Abs(A(iPiv,i))  ) iPiv=j
        End Do
        If (iPiv.Ne.i) Then
          Do j=1,n
            x         = A( i  ,j)
            A( i  ,j) = A(iPiv,j)
            A(iPiv,j) = x
            x         = B( i  ,j)
            B( i  ,j) = B(iPiv,j)
            B(iPiv,j) = x
          End Do
          T=A(I,I)
        End If
        Do J=1,n
          A(I,J)=A(I,J)/T
          B(I,J)=B(I,J)/T
        End Do
        Do K=1,n
          If (K.Ne.I) Then
            T=A(K,I)
            Do J=1,n
              A(K,J)=A(K,J)-T*A(I,J)
              B(K,J)=B(K,J)-T*B(I,J)
            End Do
          End If
        End Do
      End Do
      DeAllocate ( A  )
      Return
      End ! MatinvPiv

!C***********************************************************************  
            subroutine PrnSig_MohrCoulomb(IOpt, ntens, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions
      !            from cartesian stress vector
      !
      !  IOpt            I   I     flag to calculate principal direction (IOpt = 1)
      !  IntGlo          I   I     global ID of Gauss point or particle 
      !  S               I   R()   cartesian stress
      !  xN1, xN2, xN3   O   R()   principal direction
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      ! 
      !-------------------------------------------------------------------
      
      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt, ntens
        double precision, intent(in), dimension(ntens) :: S ! size of stress tensor (6 for 3D, 4 for 2D)
        double precision, intent(out), dimension(3) :: xN1, xN2, xN3
        double precision, intent(out) :: S1, S2, S3, P, Q

        if (IOpt == 1) then
          call Eig_3(0,ntens,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! Calculate principal direction
        else
          call Eig_3a(0,ntens,S,S1,S2,S3,P,Q) ! Do not calculate principal direction
        end if
     
      end subroutine PrnSig_MohrCoulomb
!C***********************************************************************
      
      subroutine Eig_3(iOpt, ntens, St, xN1, xN2, xN3, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions 
      !            from cartesian stress vector
      !
      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
      ! 
      !  IOpt            I   I     flag for output writing (IOpt = 1) 
      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
      !  xN1, xN2, xN3   O   R()   principal direction
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      ! 
      !-------------------------------------------------------------------

      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt, ntens
        double precision, intent(in), dimension(ntens) :: St
        double precision, intent(out), dimension(3) :: xN1, xN2, xN3 ! is (3) for 2D and 3D
        double precision, intent(out) :: S1, S2, S3, P, Q
        
        ! local variables
        double precision, dimension(3, 3) :: A, V ! is (3,3) for 2D and 3D
        double precision :: abs_max_s, tol, tau, sign_tau, t, c, s, temp1, temp2, temp3
        integer :: i, k, it, itmax, ip, iq, iS1, iS2, iS3

        ! Put cartesian stress vector into matrix A
        select case(ntens)
          case(4)
              A = 0.0
              A(1,1) = St(1) ! xx
              A(1,2) = St(4) ! xy = yx
              A(2,1) = St(4) ! xy = yx
              A(2,2) = St(2) ! yy
              A(3,3) = St(3) ! zz
          
              ! Set V to unity matrix
              V = 0.0
              V(1,1) = 1
              V(2,2) = 1
              V(3,3) = 1
              
          case(6)
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
        end select 
        
        ! get maximum value of cartesian stress vector
        abs_max_s = 0.0
        do i = 1, ntens
          if (abs(St(i)) > abs_max_s) abs_max_s = abs(St(i))
        end do
      
        ! set tolerance
        tol = 1d-16 * abs_max_s
        
        ! get principal stresses and directions iteratively
        it = 0
        itmax = 50
        do while ( (it < itmax) .and.  &
                 (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) > tol) )
        
          it = it + 1
          do k = 1,3
                
            if (k == 1) then
              ip = 1
              iq = 2
            else if (k ==2) then
              ip = 2
              iq = 3
            else
              ip = 1
              iq = 3
            end if
            
            if (abs(A(ip,iq)) > 1d-50) then
                
              tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
              if (tau >= 0.0) then
                sign_tau = 1.0
              else
                sign_tau = -1.0
              end if
              
              t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
              c = 1.0 / sqrt(1.0 + t**2)
              s = t * c
              
              temp1 = c * A(1, ip) - s * A(1, iq)
              temp2 = c * A(2, ip) - s * A(2, iq)
              temp3 = c * A(3, ip) - s * A(3, iq)
              A(1, iq) = s * A(1, ip) + c * A(1, iq)
              A(2, iq) = s * A(2, ip) + c * A(2, iq)
              A(3, iq) = s * A(3, ip) + c * A(3, iq)
              A(1, ip) = temp1
              A(2, ip) = temp2
              A(3, ip) = temp3

              temp1 = c * V(1, ip) - s * V(1, iq)
              temp2 = c * V(2, ip) - s * V(2, iq)
              temp3 = c * V(3, ip) - s * V(3, iq)
              V(1, iq) = s * V(1, ip) + c * V(1, iq)
              V(2, iq) = s * V(2, ip) + c * V(2, iq)
              V(3, iq) = s * V(3, ip) + c * V(3, iq)
              V(1, ip) = temp1
              V(2, ip) = temp2
              V(3, ip) = temp3

              temp1 = c * A(ip, 1) - s * A(iq, 1)
              temp2 = c * A(ip, 2) - s * A(iq, 2)
              temp3 = c * A(ip, 3) - s * A(iq, 3)
              A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
              A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
              A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
              A(ip, 1) = temp1
              A(ip, 2) = temp2
              A(ip, 3) = temp3
            end if ! A(ip,iq)<>0
            
          end do ! k
          
          ! optional output writing
          if (iOpt == 1) then
!            call GiveMessage('error: ' // trim(String(abs(A(1, 2)) + abs(A(2, 3)) + abs(A(1, 3)))))
          end if
          
        end do ! while
     
        ! get principal stresses from diagonal of A
        S1 = A(1, 1)
        S2 = A(2, 2)
        S3 = A(3, 3)
      
        ! derived invariants
        P = (S1 + S2 + S3) / 3.
        Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )

        ! Sort eigenvalues S1 <= S2 <= S3
        iS1 = 1
        iS2 = 2
        iS3 = 3
        
        if (S1 > S2) then
          t   = S2
          S2  = S1
          S1  = t
          it  = iS2
          iS2 = iS1
          iS1 = it
        end if
        
        if (S2 > S3) then
          t   = S3
          S3  = S2
          S2  = t
          it  = iS3
          iS3 = iS2
          iS2 = it
        end if
        
        if (S1 > S2) then
          t   = S2
          S2  = S1
          S1  = t
          it  = iS2
          iS2 = iS1
          iS1 = it
        end if
        
        ! get corresponding principal directions from V
        do i = 1, 3 ! is 3 for 2D and 3D
          xN1(i) = V(i, is1)
          xN2(i) = V(i, is2)
          xN3(i) = V(i, is3)
        end do
      
      end subroutine Eig_3

      
      subroutine Eig_3a(iOpt, ntens, St, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses from cartesian stress vector
      !
      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
      ! 
      !  IOpt            I   I     flag for output writing (IOpt = 1) 
      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      ! 
      !-------------------------------------------------------------------

      
      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt, ntens
        double precision, intent(in), dimension(ntens) :: St ! is (6) for 3D and (4) for 2D
        double precision, intent(out) :: S1, S2, S3, P, Q
        
        ! local variables
        integer :: IDim
        double precision, dimension(:,:), allocatable :: A ! (3,3) for 2D and 3D
        double precision :: abs_max_s, tol, tau, sign_tau, t, c, s, temp1, temp2, temp3
        integer :: i, k, it, itmax, ip, iq
        
        IDim = 3
        allocate(A(IDim,IDim)) 
        
        select case (ntens)       
        case(4) ! for 2D
          A = 0.0
          A(1,1) = St(1) ! xx
          A(1,2) = St(4) ! xy = yx
          A(2,1) = St(4) ! xy = yx
          A(2,2) = St(2) ! yy
          A(3,3) = St(3) ! zz
        case(6) ! for 3D
          A(1,1) = St(1) ! xx
          A(1,2) = St(4) ! xy = yx
          A(1,3) = St(6) ! zx = xz

          A(2,1) = St(4) ! xy = yx
          A(2,2) = St(2) ! yy
          A(2,3) = St(5) ! zy = yz

          A(3,1) = St(6) ! zx = xz
          A(3,2) = St(5) ! zy = yz
          A(3,3) = St(3) ! zz
        end select

        ! get maximum value of cartesian stress vector
        abs_max_s = 0.0
        do i = 1,ntens
          if (abs(St(i)) > abs_max_s) abs_max_s = abs(St(i))
        end do
      
        ! set tolerance
        tol = 1d-20 * abs_max_s
        
        ! get principal stresses and directions iteratively
        it = 0
        itmax = 50
        do while ( (it < itmax) .and.  &
                   (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) > tol) )
        
          it = it + 1
          do k = 1,3
            if (k == 1) then
              ip = 1
              iq = 2
            else if (k ==2) then
              ip = 2
              iq = 3
            else
              ip = 1
              iq = 3
            end if
            
            if (abs(A(ip,iq)) > 1d-50) then
                
              tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
              if (tau >= 0.0) then
                sign_tau = 1.0
              else
                sign_tau = -1.0
              end if
              
              t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
              c = 1.0 / sqrt(1.0 + t**2)
              s = t * c

              temp1 = c * A(1, ip) - s * A(1, iq)
              temp2 = c * A(2, ip) - s * A(2, iq)
              temp3 = c * A(3, ip) - s * A(3, iq)
              A(1, iq) = s * A(1, ip) + c * A(1, iq)
              A(2, iq) = s * A(2, ip) + c * A(2, iq)
              A(3, iq) = s * A(3, ip) + c * A(3, iq)
              A(1, ip) = temp1
              A(2, ip) = temp2
              A(3, ip) = temp3

              temp1 = c * A(ip, 1) - s * A(iq, 1)
              temp2 = c * A(ip, 2) - s * A(iq, 2)
              temp3 = c * A(ip, 3) - s * A(iq, 3)
              A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
              A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
              A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
              A(ip, 1) = temp1
              A(ip, 2) = temp2
              A(ip, 3) = temp3
  
            end if ! A(ip,iq)<>0
            
          end do ! k
          
          ! optional output writing
          if (iOpt == 1) then
!            call GiveMessage('error: ' // trim(String(abs(A(1, 2)) + abs(A(2, 3)) + abs(A(1, 3)))))
          end if
          
        end do ! while
     
        ! get principal stresses from diagonal of A
        S1 = A(1, 1)
        S2 = A(2, 2)
        S3 = A(3, 3)
      
        ! derived invariants
        P = (S1 + S2 + S3) / 3.
        Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )

        ! Sort eigenvalues S1 <= S2 <= S3
        if (S1 > S2) then
          t   = S2
          S2  = S1
          S1  = t
        end if
        
        if (S2 > S3) then
          t   = S3
          S3  = S2
          S2  = t
        end if
        
        if (S1 > S2) then
          t   = S2
          S2  = S1
          S1  = t
        end if
        
      end subroutine Eig_3a
      
!C
!C***********************************************************************
      Logical Function LEqual(A,B,Eps)
!C***********************************************************************
!C
!C     Returns .TRUE.  when two real values are (almost) equal,
!C             .FALSE. otherwise
!C
!C I   A,B  : Two real values to be compared
!C I   Eps  : Toleration (Magnitude ~= 1E-5)
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
!C***********************************************************************
      LEqual =.True.
      If (A .Eq. B) Return
      If (DAbs(A-B) .LT. 0.5D0*Eps*( DAbs(A) + DAbs(B) + Eps ) )Return
      LEqual =.False.
      Return
      End     ! function LEqual
!C
!C***********************************************************************
      Subroutine CrossProd(xN1,xN2,xN3)
!C***********************************************************************
!C
!C     Returns cross product of xN1 and xN2
!C
!C I   xN1,xN2 : Two basic vectors
!C O   xN3     : Resulting vector
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*)
!C***********************************************************************

      xN3(1) = xN1(2)*xN2(3) - xN1(3)*xN2(2)
      xN3(2) = xN1(3)*xN2(1) - xN1(1)*xN2(3)
      xN3(3) = xN1(1)*xN2(2) - xN1(2)*xN2(1)

      Return
      End     ! Subroutine CrossProd
!C
!C***********************************************************************
      Double Precision Function ArcSin(X,ie)
!C***********************************************************************
!C
!C     Returns the Arc Sine of X
!C
!C I   X : Input value
!C
!C     Note : In stead of using default routine DASIN we use this one
!C            because �X� can be slightly beyond 1 and this will give
!C            a RTE using DASIN(X)
!C
!C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
!C***********************************************************************
      Ie=0
      S = (1-X*X)
!      If (S .Lt. -1E-10) Ie=1
!      If (S .Lt. -1E-10) Write(*,1) X,S
!      If (S .Lt. -1E-10) Write(2,1) X,S
    1 Format(' ArcSin(',1x,1p,e13.5e3,') , S =',1x,1p,e13.5e3)
      If (S.LT.0) S = 0
      S = DSQRT(S)
      ArcSin = DATan2(X,S)
      Return
      End     ! function ArcSin
!C
!C***********************************************************************
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
      subroutine setveclen(xn,n,xl)
!C**********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN(*)
      x=0
      do i=1,n
        x=x+xn(i)**2
      end do
      if (x.Ne.0) Then
        f=xl/sqrt(x)
        do i=1,3
          xn(i)=xn(i)*f
        end do
      end if
      return
      end ! setveclen

!C**********************************************************************
!C End Of file
!C**********************************************************************

     
     
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
     
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
 !███╗░░░███╗░█████╗░░██████╗░██████╗
 !████╗░████║██╔══██╗██╔════╝██╔════╝
 !██╔████╔██║██║░░╚═╝╚█████╗░╚█████╗░
 !██║╚██╔╝██║██║░░██╗░╚═══██╗░╚═══██╗
 !██║░╚═╝░██║╚█████╔╝██████╔╝██████╔╝
 !╚═╝░░░░░╚═╝░╚════╝░╚═════╝░╚═════╝░
      
      
            SUBROUTINE ESM_MohrCoulombStrainSoftening(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
     DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"ESM" :: ESM
      implicit double precision (a-h, o-z) 
      CHARACTER*80 CMNAME         
      DIMENSION STRESS(NTENS),&
     DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
      !NPT(1),NOEL(1),IDSET(1),EUNLOADING(1),PLASTICMULTIPLIER(1),NUMBEROFPHASES(1)

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

    
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

        
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
    
! Initialization
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
        
        IDTask = 0
        
      IF((IStep==1).and.(TimeStep==1)) IDTask = 1
     
      IF (IDTask == 1) then ! initialisation of state variables
       STATEV(1)=PROPS(3)
       STATEV(2)=PROPS(5)
       STATEV(3)=PROPS(7)
      END IF ! IDTask = 1
      
!---Call the UMAT      
          call umat_MohrCoulombStrainSoftening(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)


      
!---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
!---Always define this value to run the simulation

    ! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
    
          


        return

      end subroutine ESM_MohrCoulombStrainSoftening
      
      

      !*USER SUBROUTINES
      SUBROUTINE UMAT_MohrCoulombStrainSoftening(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
      RPL,DDSDDT,DRPLDE,DRPLDT,&
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      !INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME
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
      SSTOL = 0.0001d0 !Tolerance Relative Error (10-3 to 10-5)
      YTOL = 0.000001d0 !Tolerance Error on the Yield surface (10-6 to 10-9)
      SPTOL = 0.0001d0 !Tolerance Softening Parameters (0.0001d0)
      ctol = abs(cp-cr)*SPTOL
      phitol = abs(phip-phir)*SPTOL
      psitol = abs(psip-psir)*SPTOL
      DTmin = 0.000000001d0
      
      !Check the yield function value
      call DetermineYieldFunctionValue(IntGlo,SigC,c,phi,F)
      
      !If F<0 then the behaviour is elastic --> Return
      if (F <= YTOL) then
        IPL = 0
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
        do while (abs(F) > YTOL.and.n < 10) !The correction is needed
            n = n + 1
            call CalculateEpsPEq(EpsP,EpsPEq)             !Determine Equivalent plastic Strain (EpsPEq)
            call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,& 
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

      Sig0 = Sig + alpha0*DSig ! = Sig0
      Sig1 = Sig + alpha1*DSig ! = SigE
        
      call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)
      call DetermineYieldFunctionValue(IntGlo,Sig1,c,phi,F1)

      F=YTOL+1000
      Its = 0 !Counter

      do while (abs(F) > YTOL.and.Its < 1000)
        alpha = alpha1 - F1*(alpha1-alpha0)/(F1-F0)
        SigNew = Sig + alpha*DSig
       
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
      if (J2 > C00000) then

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
      ! In this case aSmooth=0.0005*c*cot(phi) and LodeT=25�.
      ! If aSmooth=0 and LodeT=30� the original Mohr-Coulomb is obtained.
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
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
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
      COTPHI = CPHI/SPHI
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
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
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
        COTPSI = CPSI/SPSI
        aSmooth = C00P50*COH*COTPSI !Smoothing parameter
        ASPSI2 = aSmooth*aSmooth*SPSI*SPSI
        if (abs(psi2) == C00000)then
            ASPSI2 = C00000
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
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians

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
      
      EpsPEq = sqrt(C1*EpsDev(1)*EpsDev(1) + C1*EpsDev(2)*EpsDev(2) +  C1*EpsDev(3)*EpsDev(3) +&
                     C2*EpsP(4)*EpsP(4) + C2*EpsP(5)*EpsP(5) + C2*EpsP(6)*EpsP(6))

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
      !double precision :: k
      !!In variables
      !double precision, intent(in), dimension(3) :: DSPDPEq
      !double precision, intent(in), dimension(6) :: DEpsPEqDPS
      !double precision, intent(in), dimension(6) :: DEpsP
      !!Out variables
      !double precision, intent(out), dimension(3) :: Dh
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
      Implicit Double Precision (A-H,O-Z)
      Dimension S(*),xN1(*),xN2(*),xN3(*)

      If (iOpt.Eq.1) Then
        Call Eig_3_MohrCoulombStrainSoftening(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
      Else
        Call Eig_3a_MohrCoulombStrainSoftening(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
      End If
      Return
      End
      
      
      Subroutine Eig_3_MohrCoulombStrainSoftening(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
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

      
      Subroutine Eig_3a_MohrCoulombStrainSoftening(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      Implicit Double Precision (A-H,O-Z)
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


      Subroutine MatVec_MohrCoulombStrainSoftening(xMat,IM,Vec,N,VecR)
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
      Implicit Double Precision (A-H,O-Z)
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

      
      Subroutine AddVec_MohrCoulombStrainSoftening(Vec1,Vec2,R1,R2,N,VecR)
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
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
!C***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine AddVec


      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
!      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      !██╗░██████╗░█████╗░  ███╗░░░███╗░█████╗░██████╗░███████╗██╗░░░░░
!      !██║██╔════╝██╔══██╗  ████╗░████║██╔══██╗██╔══██╗██╔════╝██║░░░░░
!      !██║╚█████╗░███████║  ██╔████╔██║██║░░██║██║░░██║█████╗░░██║░░░░░
!      !██║░╚═══██╗██╔══██║  ██║╚██╔╝██║██║░░██║██║░░██║██╔══╝░░██║░░░░░
!      !██║██████╔╝██║░░██║  ██║░╚═╝░██║╚█████╔╝██████╔╝███████╗███████╗
!      !╚═╝╚═════╝░╚═╝░░╚═╝  ╚═╝░░░░░╚═╝░╚════╝░╚═════╝░╚══════╝╚══════╝
! 
!!c------------------------------------------------------------------------------
!!c-----------------------------------------------------------------------------
!!c     University Karlsruhe Institute of Technology KIT, Germany
!!c     University del Norte, Colombia
!!c     Jan, 2018
!!c     Dr.Ing- William Fuentes
!!c     
!!c     UMAT for:
!!c     IS- HP Wolffersdorff      PROPS(17)=1
!!c     ISA-HP Wolffersdorff      PROPS(17)=0
!!c     
!!C  
!!c    Umat written in Voigt notation    
!!c------------------------------------------------------------------------------
!!c-----------------------------------------------------------------------------
!      subroutine umat_ISA(stress,statev,ddsdde,sse,spd,scd,&
!      rpl,ddsddt,drplde,drpldt,&
!      stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
!      ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
!      celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!      implicit none
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------
!      character*80 cmname
!      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,&
!      layer, kspt, kstep, kinc
!      double precision stress(ntens), statev(nstatv),&
!      ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),&
!      stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),&
!      props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
!      double precision sse, spd, scd, rpl, drpldt, dtime, temp, &
!      dtemp, pnewdt, celent
!	 
!      real*8 STRESSEl(ntens), DDSDDEEl(ntens,ntens), stressu(ntens)&
!      ,normSTRESSU 
!      integer FirstInc
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------        
!!
!!c   
!      STRESSEl(1:ntens)=statev(57:56+ntens)
!      STRESSU(1:ntens)=STATEV(63:62+ntens)	     
!      call normS(STRESSU,normSTRESSU, ntens)		  
!!c     Flag for first increment	      
!      if ((normSTRESSU==0.0d0)) then
!      STRESSU=STRESS   
!      call normS(STRESSU,normSTRESSU, ntens)	
!      endif !  ((kinc.le.1).and.(kstep.le.1))          
!
!!c   
!       if (normSTRESSU.ne.0.0d0) then
!      call umatISA(stressU,statev,ddsdde,sse,spd,scd,&
!      rpl,ddsddt,drplde,drpldt,&
!      stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
!      ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
!      celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!        else
!        STRESSU=0.0d0
!        DDSDDE=0.0d0		
!      endif !  (normSTRESSU.ne.0.0d0)
!
!      call  PhantomElastic(STRESSEl,STATEV,DDSDDEEl,SSE,SPD,SCD,&
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)     
!
!	 
!       DDSDDE=DDSDDE +DDSDDEEl
!       stress=stressu +STRESSEL 
!       statev(57:56+ntens)=STRESSEl(1:ntens)
!       statev(63:62+ntens)=STRESSU(1:ntens)
!	   
!      end subroutine umat_ISA
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------	  
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------
!      subroutine umatISA(stress,statev,ddsdde,sse,spd,scd,&
!      rpl,ddsddt,drplde,drpldt,&
!      stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
!      ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
!      celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!      implicit none
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------
!      character*80 cmname
!      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,&
!      layer, kspt, kstep, kinc
!      double precision stress(ntens), statev(nstatv),&
!      ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),&
!      stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),&
!      props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
!      double precision sse, spd, scd, rpl, drpldt, dtime, temp, &
!      dtemp, pnewdt, celent
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------
!      real*8 lambda, np, ne, ei0, lambdac, npc, ec0, nu, Mc, nd,fb0,&
!      R, mR,betaR, chi, cz, zfab, c, rK, Ad, zmax, epsF, chimin,&
!      hb(ntens),cb(ntens),zb(ntens), deveps,hbtrial(ntens),fyieldtrial,&
!      rb(ntens), rbiso(ntens), Y0max, Kd, Kh, term1, Gh,	&
!      EEr(ntens,ntens), EE(ntens,ntens), delta(ntens), DstranU(ntens),&
!      delta2(ntens,1), delta2T(1,ntens), rbiso2(ntens,1),&
!      rbiso2T(1,ntens), devStress(ntens), p, q, void, vec1(NTENS)&
!     ,Nb(ntens), DstranUU(ntens), vec1N(1,NTENS), vecN1(NTENS,1) &
!     ,cbmax(ntens), cbbar(ntens), dotgamma ,dt, tolsubt,nDstran&
!     ,hb0(ntens), cb0(ntens), Jmom(ntens,ntens), rho, hbb(NTENS)&  
!     ,term2, term3 , rec, ei0s, ec0s , ei, ec,nnew, fb, fd &
!     ,gb,g0, r0, rd, rbb(ntens), r0b(ntens), rdb(ntens),dotEps(ntens)&
!     ,ndev(ntens),FPTL,  nndev, mflow1(ntens), Ydegree, nY &
!     ,zitot, Y0e, Nhp(ntens),vecNN(NTENS,NTENS),Kw , trace,pw&
!     ,stress0(NTENS),pcutmin, Ymax, gs, ECoeff, betaR0, yh&
!     ,udevstress(NTENS), dfdsig(NTENS), zbb(NTENS), zdot(NTENS)&
!     ,mb, zmax0, fdamage, chi0, fsos, recZi, recZ, rmin, eps(ntens)&
!     ,ndevE(ntens),normrbiso, normepsdev, ENb(ntens), Emb(ntens)&
!     , TolY, eaccPar, chimax , fe0, hbU(NTENS)
!      integer maxnint, isub  ,nsub,firstinc, issub, Nssub,maxnint2&
!      ,Load ,isElast, isTension, isisa, i
!      real*8 sq2, sq3,exprho, pcut,tolsubbt, dts, fdt,FirstINCR &
!     ,normndev, gammaW, EEinv(ntens,ntens),rhomin, fd0, zfb
!     
!      real*8 phic, hs, ed0, alpha, beta, mt,hatT(ntens), hTd(ntens)&
!     ,Fm,Fm2, sinphic, a, ed , fe,trStress,expBauer,NormhatT2, fs&
!     ,n, Isym(ntens,ntens), nstran,normhb,uhb(ntens),vechh(NTENS,NTENS)&
!     ,vecEE(NTENS,NTENS) ,IsymES(ntens,ntens),  betahe, Bhe, Hhe&
!     , che, mflow(ntens), EE2(ntens,ntens),eacc,ndstranU, rho2&
!     ,  ismr, devStressU(ntens),normz, NZ(NTENS),EEhat(ntens,ntens)&
!      ,beta_hor, termz,termz0, termExp
!
!      parameter (sq2=1.4142135623730950488016887242097d0)
!      parameter (sq3=1.7320508075688772935274463415059d0)  
!!c            
!      parameter (tolsubt=1.0d-5, maxnint=5000, exprho=1, pcut=0.001, &
!      pcutmin=0.001d0, tolsubbt=5.0d-5,maxnint2=5000)
!
!
!	     isElast=0
!      if ((isElast==1)) then
!      call  UMATElastic(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)     
!      return     
!      endif ! if (isElast==1) then
!
!
!
!
!
!
!    
!      term1=dble(kstep*kinc)+kinc
!      if (term1==0) then
!      firstinc=1
!      else
!      firstinc=0
!      endif 
!      if (firstinc==1) then
!!c     Initialization intergranular strain      
!      call Initialasv(statev, nstatv, props, nprops, ntens)
!      endif     
!!c       
!!c	 
!!c------------------------------------------------------------------------------
!!c     1) Read parameters
!!c------------------------------------------------------------------------------
!!c
!      ismr=1.0d0
!      isisa=0
!      phic=props(1)
!      Kw=props(2)
!      hs=props(3)       
!      n=props(4)
!      ed0=props(5) 
!      ec0=props(6)       
!      ei0=props(7) 
!      alpha=props(8) 
!      beta=props(9)
!      mt=props(10)
!      mr=props(11)      
!      R=props(12)
!      betaR=props(13)
!      chi=props(14)
!      chi0=chi
!      chimax=props(15)
!      eaccPar=props(16)
!      isisa=int(props(17))
!      cz=(props(18))
!      beta_hor=(props(19))
!      zmax=1.0d0
!      
!      Mc=6.0d0*sin(phic)/(3.0d0-sin(phic))
!      c=3.0d0/(3.d0+Mc)      
! 	
!!c
!!c     Delta Kronecker
!      delta=0.0d0
!      delta(1:3)=1.0d0	 
!      delta2(:,1)=delta(:)	
!      delta2T=Transpose(delta2)	
!      
!      call Isym1(Isym,ntens)  
!      IsymES=Isym
!      IsymES(4:ntens,:)=2.0d0*Isym(4:ntens,:)           
!!c	 
!!c------------------------------------------------------------------------------
!!c     2) Read state variables
!!c------------------------------------------------------------------------------
!!c   
!      void=statev(1)            ! void ratio
!      pw=statev(2)              ! pore water pressure
!      hb(1:ntens)=statev(3:2+ntens)    ! intergranular strain
!      cb(1:ntens)=statev(9:8+ntens)    ! back-intergranular strain      
!      eps(1:ntens)=stran(:)   !statev(51:50+ntens) ! total strain
!      zb(:)=statev(15:14+ntens) ! tensor Z
!      eacc=statev(30)
!!c     Variables at the beginning
!      hb0=hb     
!      cb0=cb    
!      
!!c      For evaluation of CSR under cyclic loading.     
!!c      if (props(2).ge.1.0d0) then
!!c          term1=abs(eps(1))*100.0d0
!!c      if (term1.ge.props(2)) then
!!c          stop
!!c      endif
!!c      endif
!!c
!!c 
!!c
!!c------------------------------------------------------------------------------
!!c     3) Start subincrements
!!c------------------------------------------------------------------------------
!!c  	
!!c     Strain rate
!      if (dtime.ne.0.0d0) then
!      dotEps=Dstran/dtime
!      else
!      dotEps=0.0d0
!      endif
!!c     Norm of the strain increment      
!      call normE(Dstran,nDstran, ntens) 
!!c     Number of subincrements      
!      nsub = max(int(nDstran/tolsubt),1)
!      if (nsub.ge.maxnint) nsub=maxnint
!!c     Subincrement of strain       
!      dstranU=dstran/(DBLE(nsub))
!!c     Unit strain rate DstranUU         
!      Call unitE(DstranU,DstranUU, ntens)
!!c     Subincrement time            
!      dt=dtime/(DBLE(nsub))
!      DDSDDE=0.0d0          ! initialization of jacobian
!      
!      
!
!      
!      do isub=1, Nsub
!!c	  
!!c     From total to effective stress	  
!      stress=stress+pw*delta	
!      p=-1.0d0/3.0d0*trace(stress,ntens)
!      isTension=0
!      if (p.le.pcutmin) then
!      stress=-pcutmin*delta
!      isTension=1
!      endif
!	  
!      stress0=stress  ! initial effective stress 
!      hb0=hb    ! initial intergranular strain
!      cb0=cb    ! initial back-intergranular strain
!!c	   
!!c------------------------------------------------------------------------------
!!c     4) Elastic stiffness
!!c------------------------------------------------------------------------------
!!c
!!       
!!C     Relative stress
!      trStress=(trace(stress, ntens))
!      hatT=stress/trStress  
!!C     Relative deviator stress
!      hTd=hatT-1.0d0/3.0d0*delta 
!      
!
!      CALL pq(Stress,p,q, ntens) 
!      call dev(Stress,devStress, ntens)  ! deviator stress
!      rb=devstress/p                     ! stress ratio tensor
!      rbiso=rb/(sq2/sq3*Mc)              ! normalized stress ratio tensor
!      call normS(rbiso,normrbiso, ntens)
!           
!!C ---------------------------------------------------------------  
!!C     Material constant a
!      sinphic=sin(phic) !sin of phi_c
!      a=sq3*(3.d0-sinphic)/(2.0d0*sq2*sinphic)
!!C	Material constant F 
!
!      call getThetaS(devStress,gS, c,ntens)  
!      Fm=1.0d0 + (normrbiso/(gS))*(gS-1.0d0)
!      if (gS.gt.1.0d0) gS=1.0d0
!      if (gS.lt.c) gS=c         
!      Fm2=Fm**2.d0 
!!C ---------------------------------------------------------------       
!!C     Minimum, critical and maximum void ratio (Bauers law)  
!      expBauer=exp(-(-trStress/hs)**n)		
!      ed = ed0*expBauer
!      ec = ec0*expBauer
!      ei = ei0*expBauer
!!C ---------------------------------------------------------------          
!!C     Barotropy factor fb 
!      fb=(hs/n*(ei0/ec0)**beta*(1.d0+ei)/ei*(-trStress/hs)**(1.d0-n))* &
!      (3.d0+a**2.0d0-a*sq3*((ei0-ed0)/(ec0-ed0))**alpha)**(-1.0d0)
!!C     Picnotropy factor fe (2.70)        
!      fe0 =(ec/void )**beta
!!C     Density factor fd (2.71)        
!      fd =((void-ed)/(ec-ed))**alpha 
!      
!!C --------------------------------------------------------------- 
!!c     Extension for cyclic mobility !Fuentes 2018
!!c             
!      vec1=hb-cb
!      Call unitE(vec1,Nb, ntens)         
!      call doubleEE(zb,Nb,termz0, ntens)
!      termz=mb(-termz0)   !*mb(1.0d0-exp(-0.2d0*p))
!
!	  fe=(fe0-mb(fe0-1.0d0)*termz)
!
!       
!!C --------------------------------------------------------------- 
!!C     Tensor L 
!
!      call normS(hatT,term1, ntens) 
!      NormhatT2=(term1)**2.0d0 
!      fs=fb*fe/NormhatT2   
!      rbiso2(:,1)=hatT                  ! transpose
!      rbiso2T(1,:)=hatT  
!      
!	  !EEhat=(Fm2*Isym+a**2.d0*(matmul(rbiso2,rbiso2T)))   
!      EE=EEhat*fs
!
!!C     Density factor fd       
!      fd0 =((void-ed)/(ec-ed))**alpha 
!      fd=(fd0+mb(1.0d0-fd0)*termz)
!	  
!!c     Non-linear tensor    
!      Nhp=fd*fs*a*Fm*(hatT+hTd)
!
!      
!       
!!c------------------------------------------------------------------------------
!!c      CONVENTIONAL INTERGRANULAR STRAIN (Niemunis and Herle 1996)
!!c------------------------------------------------------------------------------
!!c       
!      if (isisa==1) then
!!c	   
!!c------------------------------------------------------------------------------
!!c     5a) Intergranular strain
!!c------------------------------------------------------------------------------
!!c 
!      Call unitE(hb,uhb, ntens)
!      CALL doubleEE(uhb,DstranUU,term1, ntens)        
!      CALL normE(hb,normhb, ntens)  	     
!      Load=0  ! loading 
!      if (term1.lt.0.0d0) Load=1 !unloading
! 	 
!      rho=normhb/R						
!
!!c     Fourth order tensor      
!      vecN1(:,1)=uhb        
!      vec1N(1,1:3)=uhb(1:3) 
!      vec1N(1,4:ntens)=uhb(4:ntens)/2.0d0 ! to contravariants
!      !vechh=(matmul(vecN1, vec1N)) !covariant-contravariant
!      
!      if (Load==1) then				
!!c     Unloading
!      hb=hb0+dstranU
!      else
!!c     Loading
!
!      vecNN=IsymES -vechh*(rho**betaR)
!      !hb=hb0+matmul(vecNN,dstranu)
!      endif	
!      
!      
!
!!c     Check 1
!      Call normE(hb, term3, ntens)     
!      term3=term3/R
!!c    Correction      
!      TolY=-1.0d-6*R 
!      if (dt.ne.0.0d0) then
!      if (term3.ge.(1.0d0+TolY)) then        
!      Call unitE(hb,uhb, ntens)  	  
!      hb=R*uhb     
!      endif ! term1, term2, term3
!      endif ! (dt.ne.0.0d0)
!  
!      
!      
!      CALL normE(hb,normhb, ntens)        
!      rho=normhb/R	      
!      if (rho.ge.1.0d0) then
!      rho=1.0d0
!      Call unitE(hb,uhb, ntens)      
!      hb=R*uhb
!      endif      
!!c	   
!!c------------------------------------------------------------------------------
!!c     5b) Jacobian 
!!c------------------------------------------------------------------------------
!!c
!      stress=stress0 	
!
!     
!!c     Jacobian 
!      term1=1.0d0
!      term3=1.0d0
!      
!      vecN1(:,1)=Nhp(:)
!      vec1N(1,1:3)=uhb(1:3)
!      vec1N(1,4:ntens)=uhb(4:ntens)/2.0d0      
!      !vecNN=matmul(vecN1, vec1N) ! contravatiant-contravariant
!
!      if (ismr==1.0d0) then
!      vecEE=((rho**chi)*mT+(1.0d0-rho**chi)*mR)*EE
!      
!      else
!      vecEE=mr*EE
!      endif
!      if (Load==1) then
!!c     Unloading      
!      term3=rho**chi*(mR-mT)
!      !Jmom=vecEE+term3*matmul(EE, vechh)      
!      else
!!c     Loading        
!      term3=rho**chi*(1.0d0-mT)
!      !Jmom=vecEE+term3*matmul(EE, vechh)+&
!      !rho**chi*vecNN
!      endif
!      yh=rho**chi
!!c	   
!!c------------------------------------------------------------------------------
!!c     6) ISA INTERGRANULAR STRAIN (Fuentes, 2014)
!!c------------------------------------------------------------------------------
!!c       
!      ELSEif (isisa==0) then
!      hbtrial=hb+DstranU
!!c     Yield function
!      vec1=hbtrial-cb
!      call normE(vec1,term1, ntens) 
!      fyieldtrial=term1-R/2.0d0
!      Load=0     
!      if ((fyieldtrial.lt.0.0d0)) Load=1
!      if (Load==1) then
!!c     Elastic	
!      hb=hb+DstranU
!      termz=0.0d0
!	  termz0=0.0d0
!      else !(fyield.le.0.0d0) then
!!c     Plastic
!!c     Normal to the yield surface Nb
!       
!           
!      vec1=hbtrial-cb
!      Call unitE(vec1,Nb, ntens)     
!      hbb=R*Nb
!      vec1=(hbb-hb)
!      ! check
!      Call normE(vec1, term1, ntens)
!      if (term1==0.0d0) then
!      vec1=Nb
!      endif
!      
!      Call unitE(vec1,vec1, ntens)  
!      Call unitE(hb,hbU, ntens)  
!      call doubleEE(vec1,hbU,term2, ntens)   
!      term2=abs(term2)    
!         betaR=beta_hor+(betaR-beta_hor)*term2*(1.0d0-termz)   
!                    
!!c     Projection of back-intergranular strain cbmax      
!      cbmax=R/2.0d0*DstranUU     
!!c     Hardening function back-intergranular strain      
!      cbbar=(cbmax-cb)/R*betaR  
!!c     Consistency parameter      
!      call doubleEE(Nb,dotEps,dotgamma, ntens) 
!      call doubleEE(cbbar,Nb,term1, ntens)
!      dotgamma=(dotgamma)/(1.0d0+term1) 
!!c      if (dotgamma.le.0.0d0) dotgamma=0.0d0
!!c     Update 
!      cb=cb0+betaR/R*(cbmax-cb0)*&
!      dotgamma*dt/(1.0d0+betaR/R*dotgamma*dt)       
!      
!      term1=dotgamma*dt/(R/2.0d0)
!          
!      hb=hb0+(dstranu-term1*(hb0-cb))/(1.0d0+term1) 
!                
!      
!!c     Normal to the yield surface Nb    
!      vec1=hb-cb
!      Call unitE(vec1,Nb, ntens)   
!!c     Back intergranular strain 
!      if (dotgamma.ne.0.0d0) then     
!      cb=hb-R/2.0d0*Nb
!      endif
!        
!           
!!c      hb=cb+R/2.0d0*Nb
!!c     Check 1
!      Call normE(cb, term1, ntens)
!      term1=term1-R/2.0d0
!!c     Check 2
!      Call doubleEE(cb,DstranUU,term2, ntens)  
!      term2=term2-R/2.0d0
!!c     Check 3
!      Call normE(hb, term3, ntens)     
!      term3=term3/R
!!c    Correction      
!      TolY=-1.0d-6*R 
!      if ((dt.ne.0.0d0).and.(dotgamma.ne.0.0d0)) then
!      if ((term1.ge.TolY).or.(term2.ge.TolY).or.&
!      (term3.ge.(1.0d0+TolY))) then        
!       
!      Call unitE(cb,vec1, ntens)
!      cb=R/2.0d0*vec1
!      vec1=hb-cb
!      Call unitE(vec1,Nb, ntens)   	  
!      hb=cb+R/2.0d0*Nb     
!      endif ! term1, term2, term3
!      endif ! (dt.ne.0.0d0)
!      vec1=hb-cb
!      Call unitE(vec1,Nb, ntens)  
!      endif
!
!!c     Projected intergranular strain
!      hbb=R*Nb        
!      vec1=hbb-hb
!      call normE(vec1,term1, ntens)    
!      term1=term1/(2.0d0*R) 
!      if (load==1) then
!      rho=0.0d0
!      else
!      rho=1.0d0-term1 !**exprho
!      endif
!      if (rho.gt.(1.0d0)) rho=1.0d0
!      if (rho.lt.(0.0d0)) rho=0.0d0	 
!              
!!c	   
!!c------------------------------------------------------------------------------
!!c     7) JACOBIAN
!!c------------------------------------------------------------------------------
!!c 
!      if (Load==1) then
!!c     Elastic	
!      Jmom=mR*EE
!      yh=0.0d0
!      else 
!!c     Plastic 
!
!!c     Factor to reduce plastic accumulation
!      call doubleEE(Nb,DstranUU,term2, ntens)
!      if (term2.le.0.0d0) term2=0.0d0 
!      chi=chi0+eacc*(chimax-chi0)
!      term3=rho**chi
!      yh=term2*term3          
!
!      vecN1(:,1)=Nhp(:)        
!      vec1N(1,:)=DstranUU(:) !wil1
!      vec1N(1,4:ntens)=vec1N(1,4:ntens)/2.0d0 ! stress-type
!      !vecNN=matmul(vecN1, vec1N)        
!     
!
!!c     Jacobian of the subincrement     
!
!!c     Factor of stiffness increase 
!     
!      term1=mR+(1.0d0-mR)*yh*ismr
!      
!     
!      call normE(hb,term2, ntens)
!      rho2=(term2/R)         
!
!      term2=rho**chi
!      term3=term2       
!      Jmom=term1*(EE+(term3)*vecNN)                 
!      endif ! Load==1
!     
!      else
!      
!      write(*,*) "intergranular strain model not recognized"
!      write(*,*) "set PROPS(16)=0 for Niemunis intergranular strain"    
!      write(*,*) "set PROPS(16)=1 for ISA intergranular strain"     
!      endif    ! ISISA==1
!!c	   
!!c------------------------------------------------------------------------------
!!c     7) Next stress and state variables
!!c------------------------------------------------------------------------------
!!c       
!      
!      
!!c     Next stress 
!      !stress=stress0+matmul(Jmom,DstranU)  
!
!!c     Tension cut 
!
!      p=-1.0d0/3.0d0*trace(stress,ntens)
!      isTension=0
!      if (p.le.pcutmin) then
!!C	 write(*,*) "tension="
!      stress=-pcutmin*delta
!      isTension=1
!      endif
!  
!      if ((Kw.gt.0.0d0)) then
!      pw=pw-Kw*trace(dstranu,ntens)*(1.0d0+void)/void 	  
!      Call dyadSS(delta,delta,vecNN, ntens) 
!      vecNN=vecNN*(1.0d0+void)/void*Kw       
!      Jmom=Jmom+vecNN
!!c      vec1=MATMUL(vecNN, DSTRANU)
!!c      term1=Kw*trace(dstranu,ntens)*(1.0d0+void)/void 
!      endif ! (Kw.gt.0.0d0)
! 
! 
!!c     From effective to total stress
!      stress=stress-pw*delta   
!!c     Void ratio      
!      void=void+trace(dstranu,ntens)*(1.0d0+void)            	 
!!c     Jacobian	
!      DDSDDE=DDSDDE+Jmom/(dble(nsub)) 
!!c     Other state variables
!      eps=eps+DSTRANU 
!      
! 
!      
!      
!            
!      call normE(DstranU, ndstranU, ntens)
!
!
!      eacc =eacc+(eaccPar/R)*&
!      ((1.0d0*(1.0d0-yh))-eacc)*ndstranU
!      ! wil2
!
!
!      zb=zb+cz*mb(q/p/(fd0*Mc*Fm)-1.0d0) &
!      *(zmax*Nb-zb)*ndstranU
!      
!
!      
!      call normE(zb, normz, ntens)    
!      if (normz.ge.zmax) then
!      zb=zmax*Nb
!      endif
!      enddo !isub=1, Nsub  
!
!!c     Return state variables
!      statev(1)=void
!      statev(2)=pw
!      statev(3:2+ntens)=hb(:)           ! intergranular strain
!      statev(9:8+ntens)=cb(:)           ! back-intergranular strain   
!      statev(15:14+ntens)=zb(:)         ! tensor Z           
!     
!!c     Variables for visualization      
!      statev(51:50+ntens)=stran(:) !eps(1:ntens)  ! total strain
!      statev(30)=eacc  
!      statev(32)=normz 
!      vec1=STRESS+pw*delta
!      CALL pq(vec1,p,q, ntens) 
!      statev(69)=p
!      statev(70)=q   
!      statev(71)=q/p/(Mc*gS)
!      statev(74)=void-ec 
!      statev(75)=rho2  	  
!
!      return
! 210  stop 'I cannot open outputfile'      
!      end subroutine umatISA
!!c------------------------------------------------------------------------------
!!c------------------------------------------------------------------------------	 
!
!      SUBROUTINE Initialasv(statev, nstatev, props, nprops, ntens)
!!c
!      implicit none
!      integer nasvdim, nprops, ntens,nstatev
!      real*8 statev(nstatev), props(nprops)
!      real*8 OC, hb(ntens), R, term1, cb(ntens)
!      OC=0.9d0 !=1 to begin with plastic conditions
!               !0<OC<1 to begin with elastic conditions
!!c
!       
!      hb(:)=statev(3:2+ntens)   ! intergranular strain
!      cb(:)=statev(9:8+ntens)   ! back-intergranular strain
!      R=props(12)
!      call normE(hb, term1, ntens)
!
!      
!      if (term1==0.0) then
!	  
!      hb=0.0d0
!      cb=0.0d0
!	  
!      hb(1)=-0.57735d0*R
!      hb(2)=-0.57735d0*R
!      hb(3)=-0.57735d0*R
!      
!      cb(1)=hb(1)+0.57735d0*R/2.0d0
!      cb(2)=hb(1)+0.57735d0*R/2.0d0
!      cb(3)=hb(1)+0.57735d0*R/2.0d0     
!      endif  !  
!
!      statev(3:2+ntens)=hb(:)   ! intergranular strain
!      statev(9:8+ntens)=cb(:)   ! back-intergranular strain   
!      
!      
!
!            
!      END SUBROUTINE Initialasv 
!!C ---------------------------------------------------------------  
!!C ---------------------------------------------------------------  
!!*USER SUBROUTINES
!!C     Heading of UMAT
!      SUBROUTINE UMATElastic(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!C
!      implicit none
!!c      INCLUDE 'ABA_PARAM.INC'
!!C
!!C --------------------------------------------------------------- 
!!C     Declarating UMAT vaiables and constants 
!      character*80 cmname
!      integer ntens, ndi, nshr, nstatev, nprops, noel, npt,&
!      layer, kspt, kstep, kinc
!      double precision stress(ntens), statev(nstatev),&
!       ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),&
!       stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),&
!       props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
!      double precision sse, spd, scd, rpl, drpldt, dtime, temp, &
!       dtemp, pnewdt, celent
!!C --------------------------------------------------------------- 
!!C     Declarating other variables
!      real*8 E,anu, G, K, ALAMDA, amu , ELMOD(NTENS,NTENS) 
!      integer   i
!!C ---------------------------------------------------------------        
!!C     Material parameters and constants 
!      E=10000. !young modulus
!      ANU=0.3 !poisson modulus
!      G=E/(2.0d0*(1.0d0+ANU))!shear modulus
!      K=E/(3.0d0*(1.0d0-2.0d0*ANU))!bulk modulus
!      ALAMDA=K-2.0d0/3.0d0*G !Lame constant
!      AMU=G !Lame constant
!!C --------------------------------------------------------------- 
!!C     Elastic modulus 
!      ELMOD=0.0d0
!      ELMOD(1,1)=ALAMDA+2.0d0*AMU
!      ELMOD(2,2)=ELMOD(1,1)
!      ELMOD(3,3)=ELMOD(1,1)
!      do i=4, ntens
!      ELMOD(i,i)=AMU
!      enddo
!      ELMOD(1,2)=ALAMDA
!      ELMOD(1,3)=ALAMDA
!      ELMOD(2,3)=ALAMDA
!      ELMOD(2,1)=ELMOD(1,2)
!      ELMOD(3,1)=ELMOD(1,3)
!      ELMOD(3,2)=ELMOD(2,3)          
!!C ---------------------------------------------------------------  
!!C     Next stress
!      !STRESS=STRESS+MATMUL(ELMOD,DSTRAN)
!!C     Consistent tangent operator      
!      DDSDDE=ELMOD
!      return
!!C ---------------------------------------------------------------           
!      END SUBROUTINE UMATElastic  
!	  
!      SUBROUTINE PhantomElastic(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!C
!      implicit none
!!c      INCLUDE 'ABA_PARAM.INC'
!!C
!!C --------------------------------------------------------------- 
!!C     Declarating UMAT vaiables and constants 
!      character*80 cmname
!      integer ntens, ndi, nshr, nstatev, nprops, noel, npt,&
!      layer, kspt, kstep, kinc
!      double precision stress(ntens), statev(nstatev),&
!       ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),&
!       stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),&
!       props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
!      double precision sse, spd, scd, rpl, drpldt, dtime, temp, &
!       dtemp, pnewdt, celent
!!C --------------------------------------------------------------- 
!!C     Declarating other variables
!      real*8 E,anu, G, K, ALAMDA, amu , ELMOD(NTENS,NTENS) 
!      integer   i
!!C ---------------------------------------------------------------        
!!C     Material parameters and constants 
!      E=20. !young modulus
!      ANU=0.45 !poisson modulus
!      G=E/(2.0d0*(1.0d0+ANU))!shear modulus
!      K=E/(3.0d0*(1.0d0-2.0d0*ANU))!bulk modulus
!      ALAMDA=K-2.0d0/3.0d0*G !Lame constant
!      AMU=G !Lame constant
!!C --------------------------------------------------------------- 
!!C     Elastic modulus 
!      ELMOD=0.0d0
!      ELMOD(1,1)=ALAMDA+2.0d0*AMU
!      ELMOD(2,2)=ELMOD(1,1)
!      ELMOD(3,3)=ELMOD(1,1)
!      do i=4, ntens
!      ELMOD(i,i)=AMU
!      enddo
!      ELMOD(1,2)=ALAMDA
!      ELMOD(1,3)=ALAMDA
!      ELMOD(2,3)=ALAMDA
!      ELMOD(2,1)=ELMOD(1,2)
!      ELMOD(3,1)=ELMOD(1,3)
!      ELMOD(3,2)=ELMOD(2,3)          
!!C ---------------------------------------------------------------  
!!C     Next stress
!      !STRESS=STRESS+MATMUL(ELMOD,DSTRAN)
!!C     Consistent tangent operator      
!      DDSDDE=ELMOD
!      return
!!C ---------------------------------------------------------------           
!      END SUBROUTINE PhantomElastic  	    
!!c     
!!c     KIT University, IBF Institute for soil mechanics and rock mechanics
!!c     Jan, 2015
!!c     Dr.Ing- William Fuentes
!!c     
!!c     Tensorial operation library for Voigt notation
!!c     
!!C     
!!C --------------------------------------------------------------- 
!!C      
!!c     
!!c     KIT University, IBF Institute for soil mechanics and rock mechanics
!!c     Jan, 2015
!!c     Dr.Ing- William Fuentes
!!c     
!!c     Tensorial operation library for Voigt notation
!!c     
!!C     
!!C --------------------------------------------------------------- 
!!C      
!       SUBROUTINE dev(A,deviator, ntens)       
!!C      returns the deviatoric portion of tensor A
!       integer ntens
!       real*8 A(ntens), traceA, deviator(ntens),trace
!       deviator=A
!       traceA=trace(A, ntens)
!       deviator(1:3)=deviator(1:3)-traceA/3.0d0
!       END SUBROUTINE dev 
!!C        
!!C --------------------------------------------------------------- 
!!C      
!       SUBROUTINE p1(A,p, ntens)       
!!C      returns the mean stress of tensor A
!       integer ntens
!       real*8 A(ntens),p
!       p=-1.0d0/3.0d0*(A(1)+A(2)+A(3))
!       END SUBROUTINE p1       
!!C        
!!C --------------------------------------------------------------- 
!!C
!       REAL*8 FUNCTION trace(A, ntens)
!!C      Returns trace  of tensor A
!       integer ntens
!       real*8  A(ntens)
!       trace=A(1)+A(2)+A(3)
!       END FUNCTION trace 
!!C        
!!C --------------------------------------------------------------- 
!!C
!       REAL*8 FUNCTION mb(term1)
!!C      Returns < > Mc-Caulay Brackets of term1
!       integer ntens
!       real*8  term1
!       mb=0.0d0
!       if (term1.ge.0.0d0) mb=term1
!       END FUNCTION mb        
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE normS(A,res, ntens) 
!!C      returns the norm of tensor A  (stress-type)   	   
!       integer ntens
!       real*8 A2(ntens), A(ntens), g1, g2, res
!       integer i              
!       A2 = A*A
!       g1=A2(1)+A2(2)+A2(3)
!       g2=0.0d0
!       do i=4, ntens 
!       g2=g2+2.0D0*(A2(i))
!       enddo
!       res=sqrt(g1+g2)
!       END SUBROUTINE normS 
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE normE(A,res, ntens) 
!!C      returns the norm of tensor A  (strain-type)          	   
!       integer ntens
!       real*8 A2(ntens), A(ntens), g1, g2, res
!       integer i
!       A2 = A*A
!       g1=A2(1)+A2(2)+A2(3)
!       g2=0.0d0
!       do i=4, ntens 
!       g2=g2+(A2(i))/4.0d0
!       enddo
!       res=sqrt(g1+2.0d0*g2)
!       END SUBROUTINE normE      
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE unitE(A,Ares, ntens) 
!!C      returns the unit tensor A  (strain-type)  	   
!       integer ntens
!       real*8 A(ntens), Ares(ntens), g1, g2, res
!       integer i 
!       
!       call normE(A,res, ntens) 
!       if (res==0.0d0) then
!       Ares=0.0d0
!       return
!       else
!       Ares=A/res
!       endif ! (res==0.0d0
!       return
!       END SUBROUTINE unitE  
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE unitS(A,Ares, ntens) 
!!C      returns the unit tensor A  (stress-type)  	   
!       integer ntens
!       real*8 A(ntens), Ares(ntens), g1, g2, res
!       integer i
!       
!       call normS(A,res, ntens) 
!       if (res==0.0d0) then
!       Ares=0.0d0
!       return
!       else
!       Ares=A/res
!       endif ! (res==0.0d0
!       return
!       END SUBROUTINE unitS                  
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE EtoS(A,Ares, ntens)
!!C      returns A (strain-type) in contravariant (stress-type) 	   
!       integer ntens
!       real*8 A(ntens), Ares(ntens), g1, g2, res
!       integer i              
!       
!       do i=4, 3+ntens
!       A(i)=A(i)/2.0d0
!       enddo
!       return
!       END SUBROUTINE EtoS    
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE StoE(A,Ares, ntens)
!!C      returns A (stress-type) in covariant (strain-type) 
!       integer ntens
!       real*8 A(ntens), Ares(ntens), g1, g2, res
!       integer i                
!       
!       do i=4, 3+ntens
!       A(i)=A(i)*2.0d0
!       enddo
!       return
!       END SUBROUTINE StoE   
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE doubleEE(A,A2,res, ntens) 
!!C      returns double contraction  between strain and strain  	   
!       integer ntens
!       real*8 A(ntens), A2(ntens),  res
!       integer i           
!       
!       res=0.0d0
!       do i=1,3
!       res=res+A(i)*A2(i)
!       enddo
!       do i=4,ntens
!       res=res+A(i)*A2(i)/2.0d0
!       enddo
!       return
!       END SUBROUTINE doubleEE         
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE doubleSS(A,A2,res, ntens)
!!C      returns double contraction  between stress and stress 	   
!       integer ntens
!       real*8 A(ntens), A2(ntens),  res
!       integer i            
!       
!       res=0.0d0
!       do i=1,3
!       res=res+A(i)*A2(i)
!       enddo
!       do i=4,ntens
!       res=res+A(i)*A2(i)*2.0d0
!       enddo
!       return
!       END SUBROUTINE doubleSS    
!!C        
!!C --------------------------------------------------------------- 
!!C      
!       SUBROUTINE dyadSS(A,B,Ares, ntens)       
!!C      returns the dyadic product of two stress variables
!       integer ntens
!       real*8 A(ntens), B(ntens), Ares(ntens, ntens)
!       integer i,j
!       Ares=0.0d0
!       do i=1, ntens
!       do j=1, ntens       
!       Ares(i,j)=Ares(i,j)+A(i)*B(j)
!       enddo
!       enddo
!       return
!       END SUBROUTINE dyadSS        
!                         
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE Emod(Elmod,K,amu,ntens)
!!C      returns the elastic modulus
!!c      K is bulk modulus and amu is shear modulus	   
!       implicit none  
!       integer ntens
!       real*8 K, amu
!       real*8 Elmod(ntens,ntens)
!       real*8 ALAMDA
!       integer i
!!c      not with ntens!!!!!!
!       ALAMDA=K-2.0d0*amu/3.0d0
!       
!       ELMOD=0.0d0
!       ELMOD(1,1)=ALAMDA+2.0d0*AMU
!       ELMOD(2,2)=ELMOD(1,1)
!       ELMOD(3,3)=ELMOD(1,1)
!       do i=4, ntens
!       ELMOD(i,i)=AMU
!       enddo
!       ELMOD(1,2)=ALAMDA
!       ELMOD(1,3)=ALAMDA
!       ELMOD(2,3)=ALAMDA
!       ELMOD(2,1)=ELMOD(1,2)
!       ELMOD(3,1)=ELMOD(1,3)
!       ELMOD(3,2)=ELMOD(2,3)     
!       return       
!       END SUBROUTINE Emod   
!       
!       
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE Isym1(Isym,ntens)      
!!c      Returns the 4 order unit tensor for symmetric tensors
!!c      Contravatiant-Contravatiant      
!       implicit none
!       integer ntens, i
!       real*8 Isym(ntens,ntens), alamda, amu,  K
!
!      AMU=1.0d0/2.0d0 !constant       
!      K=1.0d0/3.0d0
!      ALAMDA=K-2.0d0/3.0d0*AMU !constant
!
!!C --------------------------------------------------------------- 
!!C     Elastic modulus 
!      Isym=0.0d0
!      Isym(1,1)=ALAMDA+2.0d0*AMU
!      Isym(2,2)=Isym(1,1)
!      Isym(3,3)=Isym(1,1)
!	do i=4, ntens
!	Isym(i,i)=AMU
!	enddo
!      Isym(1,2)=ALAMDA
!      Isym(1,3)=ALAMDA
!      Isym(2,3)=ALAMDA
!      Isym(2,1)=Isym(1,2)
!      Isym(3,1)=Isym(1,3)
!      Isym(3,2)=Isym(2,3) 
!      return 
!      End subroutine Isym1        
!       
!!C        
!!C --------------------------------------------------------------- 
!!C
!!C        
!!C --------------------------------------------------------------- 
!!C
!
!      
!      
!       SUBROUTINE getThetaS(A,g, c,ntens) 
!!C      returns the function g of tensor A (stress type) depending on the lodes angle 
!       integer ntens
!       real*8 A(ntens),g,c, B(6),cos3theta
!       real*8 Au(ntens), sq6
!       parameter (sq6=2.4494897427831780981972840747059d0)
!       integer i               
!       B=0.0d0
!       B(1:3)=A(1:3)
!       B(4:ntens)=A(4:ntens)
!       call dev(B, B,6) 
!       call unitE(B,B, 6) 
!
!       cos3theta=B(1)**3.0d0+B(2)**3.0d0+B(3)**3.0d0 &
!      +3.0d0*B(3)*B(5)**2.0d0+3.0d0*B(1)*(B(4)**2.0d0+B(5)**2.0d0) &
!      +6.0d0*B(4)*B(5)*B(6) +3.d0*B(3)*B(6)**2.0 &
!     + 3.0d0*B(2)*(B(4)**2.0d0+B(6)**2.0)
!
!       cos3theta=-sq6*cos3theta
!       g=2.0d0*c/&
!     ((1.0d0+c)-(1.0d0-c)*cos3theta) 
!       if (g.lt.c) g=c
!       if (g.gt.1.0d0) g=1.0d0
!       return
!       END SUBROUTINE getThetaS    
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE getThetaE(A,g, c,ntens) 
!!C      returns the function g of tensor A (strain type) depending on the lodes angle	   
!       integer ntens
!       real*8 A(ntens),g,c, B(6),cos3theta
!       real*8 Au(ntens), sq6
!       parameter (sq6=2.4494897427831780981972840747059d0)
!       integer i
!!c      
!       B=0.0d0
!       B(1:3)=A(1:3)
!       B(4:ntens)=A(4:ntens)/2.0d0
!       call dev(B, B,6) 
!       call unitE(B,B, 6) 
!
!       cos3theta=B(1)**3.0d0+B(2)**3.0d0+B(3)**3.0d0 &
!      +3.0d0*B(3)*B(5)**2.0d0+3.0d0*B(1)*(B(4)**2.0d0+B(5)**2.0d0) &
!      +6.0d0*B(4)*B(5)*B(6) +3.d0*B(3)*B(6)**2.0 &
!     + 3.0d0*B(2)*(B(4)**2.0d0+B(6)**2.0)
!
!       cos3theta=-sq6*cos3theta
!       g=2.0d0*c/&
!     ((1.0d0+c)-(1.0d0-c)*cos3theta) 
!       if (g.lt.c) g=c
!       if (g.gt.1.0d0) g=1.0d0
!       return
!       END SUBROUTINE getThetaE 
!!C        
!!C --------------------------------------------------------------- 
!!C
!       SUBROUTINE pq(T,p,q, ntens) 
!!C      returns p and q invariants of the stress T  	   
!       integer ntens
!       real*8 T(ntens), p,q, trace, DevT(NTENS)
!       integer i
!       real*8 sq32
!       parameter (sq32=1.2247448713915890490986420373529d0)                
!       p=-1.0d0/3.0d0*Trace(T,ntens)
!       call dev(T,devT, ntens) 
!       call normS(devT, q, ntens)
!       q=sq32*q
!       return
!       END SUBROUTINE pq          
!                         
      
      
      
      
      
      
      
!      
!██╗░░██╗██╗░░░██╗██████╗░░█████╗░██████╗░██╗░░░░░░█████╗░░██████╗████████╗██╗░█████╗░██╗████████╗██╗░░░██╗
!██║░░██║╚██╗░██╔╝██╔══██╗██╔══██╗██╔══██╗██║░░░░░██╔══██╗██╔════╝╚══██╔══╝██║██╔══██╗██║╚══██╔══╝╚██╗░██╔╝
!███████║░╚████╔╝░██████╔╝██║░░██║██████╔╝██║░░░░░███████║╚█████╗░░░░██║░░░██║██║░░╚═╝██║░░░██║░░░░╚████╔╝░
!██╔══██║░░╚██╔╝░░██╔═══╝░██║░░██║██╔═══╝░██║░░░░░██╔══██║░╚═══██╗░░░██║░░░██║██║░░██╗██║░░░██║░░░░░╚██╔╝░░
!██║░░██║░░░██║░░░██║░░░░░╚█████╔╝██║░░░░░███████╗██║░░██║██████╔╝░░░██║░░░██║╚█████╔╝██║░░░██║░░░░░░██║░░░
!╚═╝░░╚═╝░░░╚═╝░░░╚═╝░░░░░░╚════╝░╚═╝░░░░░╚══════╝╚═╝░░╚═╝╚═════╝░░░░╚═╝░░░╚═╝░╚════╝░╚═╝░░░╚═╝░░░░░░╚═╝░░░

      Subroutine ESM_hypoplasticity(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER, &
     DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UDSM" :: UDSM
      implicit double precision (a-h, o-z) 
      CHARACTER*80 CMNAME     
      DIMENSION STRESS(NTENS), &
     DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
      !NPT(1),NOEL(1),IDSET(1),EUNLOADING(1),,PLASTICMULTIPLIER(1),NUMBEROFPHASES(1)

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
              
  
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

!---Local variables defned by the user
! e.g. integer :: var_local	  
!---User can define here additional variables     
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
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
          call umat_hypoplasticity(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)

      
!---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
!---Always define this value to run the simulation
    
! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
    
          


        return

     end subroutine ESM_hypoplasticity
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      ! Copyright (C)  2009  C. Tamagnini, E. Sellari, D. Masin, P.A. von Wolffersdorff
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
!  USA.

!c------------------------------------------------------------------------------
      subroutine umat_hypoplasticity(stress,statev,ddsdde,sse,spd,scd,&
       rpl,ddsddt,drplde,drpldt,&
       stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
       ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
       celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!c------------------------------------------------------------------------------
!c user subroutine for Abaqus
!c------------------------------------------------------------------------------
!c
!c	Author: D. Masin, based on RKF23 implementation by C. Tamagnini
!c
!c----------------------------------------------------------------------------
!c
      implicit none
!c
      character*80 cmname
!c
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,&
      layer, kspt, kstep, kinc, inittension
!c
      double precision stress(ntens), statev(nstatv),&
       ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),&
       stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),&
       props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      double precision sse, spd, scd, rpl, drpldt, dtime, temp, &
       dtemp, pnewdt, celent
!c
!c ... 1. nasvdim    = maximum number of additional state variables
!c     2. tolintT    = prescribed error tolerance for the adaptive 
!c                     substepping scheme
!c     3. maxnint    = maximum number of time substeps allowed.
!c                     If the limit is exceeded abaqus is forced to reduce 
!c                     the overall time step size (cut-back) 
!c     4. DTmin      = minimum substeps size allowed.
!c                     If the limit is exceeded abaqus is forced to reduce 
!c                     the overall time step size (cut-back)
!c     5. perturb    = perturbation parameter for numerical computation of Jacobian matrices
!c     6. nfasv      = number of first additional state variable in statev field 
!c     7. prsw       = switch for printing information
!c
!c ... declaration of local variables
!c
        logical prsw,elprsw
!c
      integer i,error,maxnint,nfev,testnan,maxninttest
        integer nparms,nasvdim,nfasv,nydim,nasv,nyact,testing
!c
        !double precision dot_vect_h
!c       
      double precision parms(nprops),theta,tolintT,dtsub,DTmin,perturb
      double precision sig_n(6),sig_np1(6),DDtan(6,6),pore
        double precision deps_np1(6),depsv_np1,norm_deps,tolintTtest
        double precision norm_deps2,pp,qq,cos3t,I1,I2,I3,norm_D2,norm_D
        double precision ameanstress,avoid,youngel,tdepel0,tdepel1,nuel
        double precision Eyoung0,Eyoung1,nu0,nu1
!c
      parameter (nasvdim = 15)
      parameter (nydim = 6+nasvdim)
!c       parameter (tolintT = 1.0d-3) ...orig value...
        parameter (tolintT = 1.0d-3) 
        parameter (tolintTtest = 1.0d-1) 
!c
!c       parameter (maxnint = 1000) ...orig value...
        parameter (maxnint = 10000)
        parameter (maxninttest = 1000)
        parameter (DTmin = 1.0d-17)
        parameter (perturb = 1.0d-5)
        parameter (nfasv = 1)
        parameter (prsw = .true.)

!c
!c ... additional state variables
!c
      double precision  asv(nasvdim)
!c
!c ... solution vector (stresses, additional state variables)
!c
      double precision  y(nydim),y_n(nydim),dy(nydim)
!c
!c
!c ... Error Management:
!c     ----------------
!c     error =  0 ... no problem in time integration
!c     error =  1 ... problems in evaluation of the time rate, (e.g. undefined 
!c                    stress state), reduce time integration substeps
!c     error =  3 ... problems in time integration, reduce abaqus load increment 
!c                    (cut-back)
!c     error = 10 ... severe error, terminate calculation
!c
      error=0
!c
!c ... check problem dimensions
!c
                
      if (ndi.ne.3) then
!c
                write(1,*) 'ERROR: this UMAT can be used only for elm.'
                write(1,*) 'with 3 direct stress/strain components'
                write(1,*) 'noel = ',noel
                error=10
!c
      endif
!c
!c ... check material parameters and move them to array parms(nparms)
!c
      call check_parms_h(props,nprops,parms,nparms,error)
!c
!c ... print informations about time integration, useful when problems occur
!c
      elprsw = .false.
      if (prsw) then
!c
!c ... print only in some defined elements
!c
                if ((noel.eq.101).and.(npt.eq.1)) elprsw = .false.
      endif
!c
!c ... define number of additional state variables
!c
      call define_h(nasv)
      nyact = 6 + nasv
      if (nyact.gt.nydim) then
          write(1,*) 'ERROR: nasvdim too small in UMAT'
          error=10
      endif
!c
!c ... suggested time substep size, and initial excess pore pressure
!c
      dtsub = statev(13)
      pore = -statev(8)
!c
!c ... initialise void ratio
!c
      if (statev(7) .lt. 0.001) then
       	    ameanstress=-(stress(1)+stress(2)+stress(3))/3
       	    avoid=0
            if(Props(16) .le. 10.0) then 
            	 if(ameanstress .lt. 0.001) then
            	   avoid=props(16)
            	 else
        	       avoid=props(16)*dexp(-(3*ameanstress/&
             	       Props(3))**props(4))
     			 end if
            else if(props(16) .gt. 10.0) then
                  avoid=props(16)-10.0
            endif
            statev(7)=avoid
      end if

!c
!c ... vector of additional state variables
!c
      asv = 0
      do i=1,nasv
        asv(i) = statev(i-1+nfasv)
      enddo
!c
!c ... compute volume strain increment and current effective stress tensor
!c
      do i=1,6        
            sig_n(i)=0
            deps_np1(i)=0
      end do
      call move_sig_h(stress,ntens,pore,sig_n)
      call move_eps_h(dstran,ntens,deps_np1,depsv_np1)

      norm_D2=dot_vect_h(2,deps_np1,deps_np1,6)!dot_vect_h(2,deps_np1,deps_np1,6)
      norm_D=sqrt(norm_D2)

!c ... check whether the strain rate from the ABAQUS is not NAN	  

      testnan=0
      call umatisnan_h(norm_D,testnan)
      if (testnan .eq. 1) then 
	     call wrista_h(3,y,nydim,deps_np1,dtime,coords,statev,nstatv,&
                   parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
	     write(1,*) 'Error in integration, noel ',noel
	     write(1,*) 'Try to decrease the global step size'
	     call xit_h
      end if
!c
!c --------------------
!c ... Time integration
!c --------------------
!c

      call iniy_h(y,nydim,nasv,ntens,sig_n,asv)
      call push_h(y,y_n,nydim)

!c ... check whether the initial state is not tensile
      inittension=0
      call check_RKF_h(inittension,y,nyact,nasv,parms,nparms)
!c
      if (elprsw) then
        write(1,*) '==================================================='
        write(1,*) 'Call of umat:'
        write(1,*) '==================================================='
        call wrista_h(3,y,nydim,deps_np1,dtime,coords,statev,nstatv,&
                  parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
      endif

!c ... Switch for elasticity in the case tensile stress is reached
      youngel=0
!c
!c ... local integration using adaptive RKF-23 method, consistent Jacobian and error estimation
!c
      if((dtsub.le.0.0d0).or.(dtsub.gt.dtime)) then
        dtsub = dtime
      endif
!c
      testing=0
      kstep = 1
      kinc = 1
!c     For use in PLAXIS, activate the following line
      !if(kstep.eq.1 .AND. kinc.eq.1) testing=1
!c     For use in ABAQUS EXPLICIT, activate the following line
      !if(kstep.eq.1 .AND. kinc.eq.1) testing=3
!c     For use in ABAQUS, the two lines above should be inactive
	
      if(norm_D.eq.0) testing=2
!c     FEM asking for ddsdde only

      nfev = 0 ! initialisation

      if(inittension.eq.0) then

      if(testing.eq.1) then
          call rkf23_update_h(y,nyact,nasv,dtsub,tolintTtest,&
                           maxninttest,DTmin,&
                           deps_np1,parms,nparms,nfev,elprsw,&
                           dtime,error)
!c ... give original state if the model fails without substepping
          if(error.eq.3) then
            do i=1,nyact        
               y(i)=y_n(i)
            end do
            error=0
          end if
      else if(testing.eq.2) then
            do i=1,nyact        
                  y(i)=y_n(i)
            end do      
      else if(testing.eq.3) then
      	temp=parms(10)
      	parms(10)=0
        call perturbate_h(y_n,y,nyact,nasv,dtsub,&
           tolintT,maxnint,DTmin,&
           deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DDtan,&
           dtime,error)      	
        parms(10)=temp
        youngel=-100
        nuel=0.3
        call calc_elasti_h(y,nyact,nasv,dtsub,tolintT,&
             maxnint,DTmin,&
            deps_np1,parms,nparms,nfev,elprsw,&
     	    dtime,DDtan,&
     	    youngel,nuel,error)
!c ... Normal RKF23 integration
      else   !inittension.eq.0 .and. testing.eq.0
          call rkf23_update_h(y,nyact,nasv,dtsub,tolintT,&
                           maxnint,DTmin,&
                           deps_np1,parms,nparms,nfev,&
                           elprsw,dtime,error)
      end if
!c
!c ... error conditions (if any)
!c
      if (error.eq.3) then
!c
!c          pnewdt = 0.25d0
!c
           write(1,*) 'UMAT: step rejected in element '&
     			,noel,' point ',npt
           call wrista_h(1,y,nydim,deps_np1,dtime,&
                     coords,statev,nstatv,&
                     parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
!c          call xit_h
!c          return
!c ...      do not do anything, we are the most likely close to the tensile region
           do i=1,nyact        
                  y(i)=y_n(i)
           end do
!c
      elseif (error.eq.10) then
!c
           call wrista_h(2,y,nydim,deps_np1,dtime,&
                     coords,statev,nstatv,&
                     parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
           call xit_h
      endif ! end error.eq.3

!c ... compute ddsdde

      call perturbate_h(y_n,y,nyact,nasv,dtsub,tolintT,maxnint,DTmin,&
           deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DDtan,&
           dtime,error)

      else ! inittension.ne.0
!c          we were initilly in the tensile stress, calc elastic
	    youngel=100
	    nuel=0.48
	    call calc_elasti_h(y,nyact,nasv,dtsub,tolintT,&
                           maxnint,DTmin,&
                           deps_np1,parms,nparms,nfev,elprsw,&
     			    dtime,DDtan,&
     			    youngel,nuel,error)
      endif ! end inittension.eq.0
!c
!c ... update dtsub and nfev
!c
      if(dtsub.le.0.0d0) then 
      	dtsub = 0
      else if(dtsub.ge.dtime) then 
      	dtsub = dtime
      end if
      statev(13)=dtsub
      statev(10)=dfloat(nfev)
!c ... convert solution (stress + cons. tangent) to abaqus format
!c     update pore pressure and compute total stresses 
!c
      call solout_h(stress,ntens,asv,nasv,ddsdde,&
                 y,nydim,pore,depsv_np1,parms,nparms,DDtan)
     
!c
!c ... updated vector of additional state variables to abaqus statev vector
!c
      do i=1,nasv
           statev(i-1+nfasv) = asv(i) 
      end do
!c
!c ... transfer additional information to statev vector
!c
      do i=1,6
           sig_np1(i)=y(i)
      end do
      pp=-(sig_np1(1)+sig_np1(2)+sig_np1(3))/3
!c
      statev(8) = -pore 
      statev(9) = pp

      if(inittension.eq.0) then
      call calc_statev_h(sig_np1,statev,parms,nparms,nasv,&
      	nstatv,deps_np1)
      end if

!c
!c -----------------------
!c End of time integration
!c -----------------------
!c
      return
      end
!c------------------------------------------------------------------------------
!c-----------------------------------------------------------------------------
      subroutine check_parms_h(props,nprops,parms,nparms,error)
!c-----------------------------------------------------------------------------
!c checks input material parameters 
!c
!c written 10/2004 (Tamagnini & Sellari)
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer nprops,nparms,i,error
!c
      double precision props(nprops),parms(nprops)
        double precision zero,one,four,pi,pi_deg
        double precision phi_deg,phi,hs,en,ed0,ec0,ei0,alpha,beta
        double precision m_R,m_T,r_uc,beta_r,chi,bulk_w,p_t
!c
        parameter(zero=0.0d0,one=1.0d0,four=4.0d0,pi_deg=180.0d0)
!c
        nparms=nprops
!c
      do i=1,nprops
                parms(i)=props(i)
      enddo
!c
!c ... recover material parameters
!c
        phi_deg=parms(1) ! critical friction angle -> yes -> 1 --> 33.1
        hs    =parms(3) ! yes -> 5 --> 4e6
        en    =parms(4) ! yes -> 6 ... this is n in Wichtmann et al. (2019) -> 0.27
        ed0   =parms(5) ! yes -> 4 -> 0.677
        ec0   =parms(6) ! yes -> 3 -> 1.054
        ei0   =parms(7) ! yes -> 2 -> 1.212
        alpha =parms(8) ! yes  -> 7 -> 0.14
        beta  =parms(9) ! yes -> 8 -> 2.5
        m_R=parms(10)  ! yes -> 10 -> 2.2
        m_T=parms(11) ! yes -> 11 -> 1.1
        r_uc=parms(12) ! -> 9 -> 1e-4
        beta_r=parms(13) ! yes -> 12 -> 0.1
        chi=parms(14) ! yes -> 13 -> 5.5
        bulk_w=parms(15) ! yes = 2.7e6 kPa      
        p_t=parms(2) ! shift of the mean effective stress -> axis translation due to cohesion ... can be zero 
        
!c
        pi=four*datan(one)
        phi=phi_deg*pi/pi_deg
        parms(1)=phi
!c
        if(phi.le.zero) then
!c       
                write(1,*) 'ERROR: subroutine CHECK_PARMS:'
                write(1,*) 'phi = ',phi
                error = 10
                return 
!c
        end if
!c
        if(m_R.lt.zero) then
!c       
                write(1,*) 'ERROR: subroutine CHECK_PARMS:'
                write(1,*) 'm_R = ',m_R
                error = 10 
                return 
!c
        end if
!c
        if(m_T.lt.zero) then
!c       
                write(1,*) 'ERROR: subroutine CHECK_PARMS:'
                write(1,*) 'm_T = ',m_T
                error = 10 
                return 
!c
        end if
!c
        if(r_uc.lt.zero) then
!c       
                write(1,*) 'ERROR: subroutine CHECK_PARMS:'
                write(1,*) 'r_uc = ',r_uc
                error = 10 
                return 
!c
        end if
!c
        if(beta_r.lt.zero) then
!c       
                write(1,*) 'ERROR: subroutine CHECK_PARMS:'
                write(1,*) 'beta_r = ',beta_r
                error = 10 
                return 
!c
        end if
!c
        if(chi.lt.zero) then
!c       
                write(1,*) 'ERROR: subroutine CHECK_PARMS:'
                write(1,*) 'chi = ',chi
                error = 10 
                return 
!c
        end if
!c 
        if(bulk_w.lt.zero) then
!c       
                write(1,*) 'ERROR: subroutine CHECK_PARMS:'
                write(1,*) 'bulk_w = ',bulk_w
                error = 10 
                return 
!c
        end if
!c 
        if(p_t.lt.zero) then
!c       
                write(1,*) 'ERROR: subroutine CHECK_PARMS:'
                write(1,*) 'p_t = ',p_t
                error = 10 
                return 
!c
        end if
!c 
      return
      end
!c-----------------------------------------------------------------------------
      subroutine define_h(nasv)
!c-----------------------------------------------------------------------------
      implicit none 
      integer nasv
!c
!c number of additional state variables 
!c must be less than  18 (otherwise change nasvdim in umat)
!c
!c    nasv(1) ... del_11  intergranular strain component
!c    nasv(2) ... del_22  intergranular strain component
!c    nasv(3) ... del_33  intergranular strain component
!c    nasv(4) ... del_12  intergranular strain component
!c    nasv(5) ... del_13  intergranular strain component
!c    nasv(6) ... del_23  intergranular strain component
!c    nasv(7) ... void    void ratio
!c
!c modified 6/2005 (Tamagnini, Sellari & Miriano)
!c
      nasv = 7
      return
      end
!c------------------------------------------------------------------------------
      double precision function dot_vect_h(flag,a,b,n)
!c------------------------------------------------------------------------------
!c dot product of a 2nd order tensor, stored in Voigt notation
!c created 10/2004 (Tamagnini & Sellari)
!c
!c flag = 1 -> vectors are stresses in Voigt notation
!c flag = 2 -> vectors are strains in Voigt notation
!c flag = 3 -> ordinary dot product between R^n vectors
!c------------------------------------------------------------------------------
      implicit none
        integer i,n,flag
      double precision a(n),b(n)
        double precision zero,half,one,two,coeff
!c
        parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
!c
        if(flag.eq.1) then
!c
!c ... stress tensor (or the like)
!c
                coeff=two
!c
        elseif(flag.eq.2) then
!c
!c ... strain tensor (or the like)
!c
                coeff=half
!c
        else
!c
!c ... standard vectors
!c
                coeff=one
!c       
        end if
!c
        dot_vect_h=zero
!c
        do i=1,n
                if(i.le.3) then
                      dot_vect_h = dot_vect_h+a(i)*b(i)
                else
                      dot_vect_h = dot_vect_h+coeff*a(i)*b(i)
                end if
        end do
!c
      return
      end
!c-----------------------------------------------------------------------------
      subroutine get_F_sig_q_h(sig,q,nasv,parms,nparms,&
               deps,F_sig,F_q,error)
!c-----------------------------------------------------------------------------
!c
!c  finds vectors F_sigma and F_q in F(y)
!c
!c  written 6/2005 (Tamagnini, Sellari & Miriano)
!c-----------------------------------------------------------------------------
        implicit none
        !double precision dot_vect_h
        
!c 
      integer nparms,nasv,ii
!c
        double precision sig(6),q(nasv),parms(nparms),deps(6)
        double precision MM(6,6),HH(nasv,6),F_sig(6),F_q(nasv)
        double precision LL(6,6),NN(6),norm_D,norm_D2
        integer istrain,error
!c
!c ... compute tangent operators
!c
		if(parms(10) .le. 0.5) then
			istrain=0 
		else 
			istrain=1
		end if

        call get_tan_h(deps,sig,q,nasv,parms,nparms,MM,&
             HH,LL,NN,istrain,error)
!c
!c ... compute F_sig=MM*deps
!c
		if (istrain .eq. 1) then
        		call matmul_h(MM,deps,F_sig,6,6,1)
        else 
        		call matmul_h(LL,deps,F_sig,6,6,1)
		        norm_D2=dot_vect_h(2,deps,deps,6)!dot_vect_h(2,deps,deps,6)
		        norm_D=sqrt(norm_D2)
                do ii=1,6
                     F_sig(ii)=F_sig(ii)+NN(ii)*norm_D
                end do
        end if
!c
!c ... compute F_q=HH*deps
!c
        call matmul_h(HH,deps,F_q,nasv,6,1)
!c       
        return
        end
!c-----------------------------------------------------------------------------
      subroutine get_tan_h(deps,sig,q,nasv,parms,nparms,MM,HH,&
     		 LL,NN,istrain,error)
!c-----------------------------------------------------------------------------
!c  computes matrices M and H for Masin hypoplastic model for clays
!c  version with intergranular strains
!c
!c  NOTE: stress and strain convention: tension and extension positive
!c
!c  written 6/2005 (Tamagnini & Sellari)
!c-----------------------------------------------------------------------------
        implicit none
!c 
      integer nparms,nasv,i,j,error
!c
        !double precision dot_vect_h
!c
        double precision sig(6),q(nasv),parms(nparms),deps(6)
        double precision eta(6),eta_dev(6),del(6),void,sig_star(6)
        double precision eta_del(6),eta_delta(6),eta_eps(6)
        double precision norm_del,norm_del2,norm_deps,norm_deps2,eta_dn2
        double precision pp,qq,cos3t,I1,I2,I3,tanpsi
        double precision a,a2,FF,fd,fs
        double precision num,den,aF,Fa2,eta_n2,norm_m,norm_m2
        double precision II(6,6),IU(6,6)
        double precision MM(6,6),HH(nasv,6),LL(6,6),NN(6),AA(6,6),m(6)
        integer istrain
        double precision m_dir(6),m_dir1(6),Leta(6),H_del(6,6),H_e(6)
        double precision load,rho
        double precision zero,tiny,half,one,two,three,six,eight,nine
        double precision onethird,sqrt3,twosqrt2,sqrt2,oneeight,ln2m1
        double precision temp1,temp2,temp3,temp4
        double precision phi,hs,en,ed0,ec0,ei0,alpha,beta,r_uc
        double precision m_R,m_T,beta_r,chi,bulk_w,p_t,sinphi,sinphi2
        double precision ec,ed,ei,bauer,fb,fe,sq2,sq3,sq6,az
!c
        parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0)
        parameter(tiny=1.0d-17,half=0.5d0,eight=8.0d0,nine=9.0d0)
        parameter(sq2=1.4142135623730951455d0,&
               sq3=1.7320508075688771931d0,&
               sq6=2.4494897427831778813d0)

!c
!c ... initialize constants and vectors
!c
        onethird=one/three
        sqrt3=dsqrt(three)
        twosqrt2=two*dsqrt(two)
        sqrt2=dsqrt(two)
        oneeight=one/eight
        onethird=one/three
        ln2m1=one/dlog(two)
!c
        do i=1,6
                do j=1,6
                        MM(i,j)=zero
                        LL(i,j)=zero
                        II(i,j)=zero
                        IU(i,j)=zero
                        H_del(i,j)=zero
                end do
                eta_del(i)=zero
                eta_delta(i)=zero
                eta_eps(i)=zero
        end do
!c
        do i=1,nasv
                do j=1,6
                        HH(i,j)=zero
                end do
        end do
!c
!c ... fourth order identity tensors in Voigt notation
!c
        II(1,1)=one
        II(2,2)=one
        II(3,3)=one
        II(4,4)=half
        II(5,5)=half
        II(6,6)=half
!c
        IU(1,1)=one
        IU(2,2)=one
        IU(3,3)=one
        IU(4,4)=one
        IU(5,5)=one
        IU(6,6)=one
!c
!c ... recover material parameters
!c
        phi	  =parms(1)
        hs    =parms(3)
        en    =parms(4)
        ed0   =parms(5)
        ec0   =parms(6)
        ei0   =parms(7)
        alpha =parms(8)
        beta  =parms(9)
        m_R=parms(10) 
        m_T=parms(11)
        r_uc=parms(12)
        beta_r=parms(13)
        chi=parms(14)
        bulk_w=parms(15)
        p_t=parms(2)
!c
        sinphi=dsin(phi)
        sinphi2=sinphi*sinphi

!c
!c ... recover internal state variables
!c
        del(1)=q(1)
        del(2)=q(2)
        del(3)=q(3)
        del(4)=q(4)
        del(5)=q(5)
        del(6)=q(6)
        void=q(7)
!c
!c ... axis translation due to cohesion (p_t>0)
!c
        sig_star(1)=sig(1)-p_t
        sig_star(2)=sig(2)-p_t
        sig_star(3)=sig(3)-p_t
        sig_star(4)=sig(4)
        sig_star(5)=sig(5)
        sig_star(6)=sig(6)
!c
!c ... strain increment and intergranular strain directions
!c
        norm_deps2=dot_vect_h(2,deps,deps,6)!dot_vect_h(2,deps,deps,6)
        norm_del2=dot_vect_h(2,del,del,6)!dot_vect_h(2,del,del,6)
        norm_deps=dsqrt(norm_deps2)
        norm_del=dsqrt(norm_del2)
!c
        if(norm_del.ge.tiny) then
!c
                do i=1,6
                        eta_del(i)=del(i)/norm_del
                end do
!c
        end if
!c
        eta_delta(1)=eta_del(1)
        eta_delta(2)=eta_del(2)
        eta_delta(3)=eta_del(3)
        eta_delta(4)=half*eta_del(4)
        eta_delta(5)=half*eta_del(5)
        eta_delta(6)=half*eta_del(6)
!c
        if(norm_deps.ge.tiny) then
!c
                do i=1,6
                        eta_eps(i)=deps(i)/norm_deps
                end do
!c
        end if
!c
!c ... auxiliary stress tensors
!c
        call inv_sig_h(sig_star,pp,qq,cos3t,I1,I2,I3)
!c
!c        if (pp.gt.tiny) then
!c
!c ... if mean stress is negative, return with MM = 0, HH = 0 and error = 10 (severe)
!c
!c                write(1,*) 'ERROR: subroutine GET_TAN:'
!c                write(1,*) 'Mean stress is positive (tension): p = ',pp
!c                error = 10
!c                return 
!c
!c        end if
!c
        
        eta(1)=sig_star(1)/I1
        eta(2)=sig_star(2)/I1
        eta(3)=sig_star(3)/I1
        eta(4)=sig_star(4)/I1
        eta(5)=sig_star(5)/I1
        eta(6)=sig_star(6)/I1   
!c
        eta_dev(1)=eta(1)-onethird
        eta_dev(2)=eta(2)-onethird
        eta_dev(3)=eta(3)-onethird
        eta_dev(4)=eta(4)
        eta_dev(5)=eta(5)
        eta_dev(6)=eta(6)
!c
!c ... functions a and F
!c
        eta_dn2=dot_vect_h(1,eta_dev,eta_dev,6)!dot_vect_h(1,eta_dev,eta_dev,6)
        tanpsi=sqrt3*dsqrt(eta_dn2)
        temp1=oneeight*tanpsi*tanpsi+&
         (two-tanpsi*tanpsi)/(two+sqrt2*tanpsi*cos3t)
        temp2=tanpsi/twosqrt2
!c
        a=sqrt3*(three-sin(phi))/(twosqrt2*sin(phi))
        a2=a*a
        FF=dsqrt(temp1)-temp2
!c
!c ... barotropy and pyknotropy functions
!c
	    bauer=dexp(-(-I1/hs)**en)
      	ed = ed0*bauer
      	ec = ec0*bauer
      	ei = ei0*bauer

      	temp1=three+a*a-a*sq3*((ei0-ed0)/(ec0-ed0))**alpha
      	if(temp1.lt.zero) stop 'factor fb not defined'
	    fb=hs/en/temp1*(one+ei)/ei*(ei0/ec0)**beta*(-I1/hs)**(one-en)
    	fe=(ec/void)**beta
            	
        fs=fb*fe
!c
		if(void.ge.ed) then
        	fd=((void-ed)/(ec-ed))**alpha
        else
        	fd=0
        end if
!c
!c
!c ... tensor L
!c
        eta_n2=dot_vect_h(1,eta,eta,6)!dot_vect_h(1,eta,eta,6)
        do i = 1,6
                do j=1,6
                        LL(i,j)=(II(i,j)*FF*FF+&
     	                  a2*eta(i)*eta(j))/eta_n2
                end do
        end do

!c
!c ... tensor NN
!c

       do i=1,6
         NN(i) = FF*a*(eta(i)+eta_dev(i))/eta_n2
       enddo
        
!c
!c ... BEGIN INTERGR. STRAIN
!c

        if(istrain .eq. 1) then
!c
!c ... loading function
!c
        load=dot_vect_h(2,eta_del,eta_eps,6)!dot_vect_h(2,eta_del,eta_eps,6)
!c
!c ... intergranular strain--related tensors
!c
        rho=norm_del/r_uc
!c
        if (rho.gt.one) then
                rho=one
        end if
!c
        call matmul_h(LL,eta_del,Leta,6,6,1)
!c
!c ... tangent stiffness M(sig,q,eta_eps)
!c
        temp1=((rho**chi)*m_T+(one-rho**chi)*m_R)*fs
!c
        if (load.gt.zero) then
!c    
                temp2=(rho**chi)*(one-m_T)*fs
                temp3=(rho**chi)*fs*fd
!c
                do i=1,6
                  do j=1,6
                    AA(i,j)=temp2*Leta(i)*eta_delta(j)&
                           +temp3*NN(i)*eta_delta(j)
                    MM(i,j)=temp1*LL(i,j)+AA(i,j)
                  end do
                end do
!c
        else
!c
                temp4=(rho**chi)*(m_R-m_T)*fs
!c
                do i=1,6
                  do j=1,6
                        AA(i,j)=temp4*Leta(i)*eta_delta(j)
                        MM(i,j)=temp1*LL(i,j)+AA(i,j)
                  end do
                end do
!c
        end if
!c
!c ... intergranular strain evolution function
!c     NOTE: H_del transforms a strain-like vector into a strain-like vector
!c           eta_del(i) instead of eta_delta(i)
!c           I = 6x6 unit matrix
!c
        if (load.gt.zero) then
!c
                do i=1,6
                  do j=1,6
                H_del(i,j)=IU(i,j)-(rho**beta_r)*eta_del(i)*eta_delta(j)
                  end do
                end do
!c
        else
!c
                do i=1,6
              H_del(i,i)=one
                end do
!c
        end if
!c
!c ... void ratio evolution function (tension positive)
!c
        do i=1,6 
                if (i.le.3) then
                  H_e(i)=one+void
                else
              H_e(i)=zero
                end if
        end do
!c
!c ... assemble hardening matrix
!c
        do i=1,nasv
                if (i.le.6) then
                        do j=1,6
                                HH(i,j)=H_del(i,j)
                        end do
                else
                        do j=1,6
                                HH(i,j)=H_e(j)
                        end do
                end if
        end do
!c       
!c ... end istrain
        else if (istrain .eq. 0) then
        
        do i=1,6 
                if (i.le.3) then
                  H_e(i)=one+void
                else
              H_e(i)=zero
                end if
        end do        
        do i=1,nasv
                if (i.le.6) then
                        do j=1,6
                                HH(i,j)=0
                        end do
                else
                        do j=1,6
                                HH(i,j)=H_e(j)
                        end do
                end if
        end do        
!c ... end istrain/noistrain switch        
        end if

        do i=1,6
           do j=1,6
                LL(i,j)=LL(i,j)*fs
           end do
           NN(i)=NN(i)*fs*fd
        end do        

        return
        end
!c-----------------------------------------------------------------------------
      subroutine iniy_h(y,nydim,nasv,ntens,sig,qq)
!c-----------------------------------------------------------------------------
!c initializes the vector of state variables
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer i,nydim,nasv,ntens
!c
      double precision y(nydim),qq(nasv),sig(ntens)
!c
      do i=1,nydim
        y(i) = 0
      enddo
!c
      do i=1,ntens
        y(i) = sig(i)
      enddo
!c
!c additional state variables
!c
      do i=1,nasv
        y(6+i) = qq(i)
      enddo
!c
      return
      end
!c------------------------------------------------------------------------------
      subroutine inv_eps_h(eps,eps_v,eps_s,sin3t)
!c------------------------------------------------------------------------------
!c calculate invariants of strain tensor
!c------------------------------------------------------------------------------
!c
      implicit none
!c
      integer i
!c
      double precision eps(6),edev(6),edev2(6),ev3
        double precision tredev3,eps_v,eps_s,sin3t
        double precision norm2,numer,denom
!c
      double precision zero,one,two,three,six
      double precision onethird,twothirds,sqrt6
!c
      data zero,one,two,three,six/0.0d0,1.0d0,2.0d0,3.0d0,6.0d0/
!c
!c ... some constants
!c
        onethird=one/three
        twothirds=two/three
        sqrt6=dsqrt(six)
!c
!c ... volumetric strain
!c
      eps_v=eps(1)+eps(2)+eps(3)
!c
      ev3=onethird*eps_v
!c
!c ... deviator strain
!c
        edev(1)=eps(1)-ev3
        edev(2)=eps(2)-ev3
        edev(3)=eps(3)-ev3
        edev(4)=eps(4)/two
        edev(5)=eps(5)/two
        edev(6)=eps(6)/two
!c
!c ... second invariant
!c
        norm2=edev(1)*edev(1)+edev(2)*edev(2)+edev(3)*edev(3)+&
           two*(edev(4)*edev(4)+edev(5)*edev(5)+edev(6)*edev(6))
!c
        eps_s=dsqrt(twothirds*norm2)
!c
!c ... components of (edev_ij)(edev_jk)
!c
        edev2(1)=edev(1)*edev(1)+edev(4)*edev(4)+edev(5)*edev(5)
        edev2(2)=edev(4)*edev(4)+edev(2)*edev(2)+edev(6)*edev(6)
        edev2(3)=edev(6)*edev(6)+edev(5)*edev(5)+edev(3)*edev(3)
        edev2(4)=two*(edev(1)*edev(4)+edev(4)*edev(2)+edev(6)*edev(5))
        edev2(5)=two*(edev(5)*edev(1)+edev(6)*edev(4)+edev(3)*edev(5))
        edev2(6)=two*(edev(4)*edev(5)+edev(2)*edev(6)+edev(6)*edev(3))
!c            
!c ... Lode angle
!c
        if(eps_s.eq.zero) then 
!c
                sin3t=-one
!c               
        else
!c
                tredev3=zero
                do i=1,6
                        tredev3=tredev3+edev(i)*edev2(i)
                end do
!c
                numer=sqrt6*tredev3
                denom=(dsqrt(norm2))**3
                sin3t=numer/denom
                if(dabs(sin3t).gt.one) then
                        sin3t=sin3t/dabs(sin3t)
                end if
!c
        end if 
!c
      return
      end
!c------------------------------------------------------------------------------
      subroutine inv_sig_h(sig,pp,qq,cos3t,I1,I2,I3)
!c------------------------------------------------------------------------------
!c calculate invariants of stress tensor
!c
!c NOTE: Voigt notation is used with the following index conversion
!c
!c       11 -> 1
!c       22 -> 2
!c    33 -> 3
!c       12 -> 4
!c       13 -> 5
!c       23 -> 6
!c
!c------------------------------------------------------------------------------
!c
      implicit none
!c
      double precision sig(6),sdev(6)
      double precision eta(6),eta_d(6),eta_d2(6)
      double precision xmin1,xmin2,xmin3
      double precision tretadev3,pp,qq,cos3t,I1,I2,I3
      double precision norm2,norm2sig,norm2eta,numer,denom
!c
      double precision half,one,two,three,six
      double precision onethird,threehalves,sqrt6,tiny
!c
      !double precision dot_vect_h
!c
      data half,one/0.5d0,1.0d0/
      data two,three,six/2.0d0,3.0d0,6.0d0/
      data tiny/1.0d-18/
!c
!c ... some constants
!c
      onethird=one/three
      threehalves=three/two
      sqrt6=dsqrt(six)
!c
!c ... trace and mean stress
!c
      I1=sig(1)+sig(2)+sig(3)
      pp=onethird*I1
!c
!c ... deviator stress
!c
      sdev(1)=sig(1)-pp
      sdev(2)=sig(2)-pp
      sdev(3)=sig(3)-pp
      sdev(4)=sig(4)
      sdev(5)=sig(5)
      sdev(6)=sig(6)
!c
!c ... normalized stress and dev. normalized stress
!c

      if(I1.ne.0) then
         eta(1)=sig(1)/I1
         eta(2)=sig(2)/I1
         eta(3)=sig(3)/I1
         eta(4)=sig(4)/I1
         eta(5)=sig(5)/I1
        eta(6)=sig(6)/I1
      else
        eta(1)=sig(1)/tiny
        eta(2)=sig(2)/tiny
        eta(3)=sig(3)/tiny
        eta(4)=sig(4)/tiny
        eta(5)=sig(5)/tiny
        eta(6)=sig(6)/tiny        
      end if
!c
      eta_d(1)=eta(1)-onethird
      eta_d(2)=eta(2)-onethird
      eta_d(3)=eta(3)-onethird
      eta_d(4)=eta(4)
      eta_d(5)=eta(5)
      eta_d(6)=eta(6)
!c
!c ... second invariants
!c
      norm2=dot_vect_h(1,sdev,sdev,6)!dot_vect_h(1,sdev,sdev,6)
      norm2sig=dot_vect_h(1,sig,sig,6)!dot_vect_h(1,sig,sig,6)
      norm2eta=dot_vect_h(1,eta_d,eta_d,6)!dot_vect_h(1,eta_d,eta_d,6)
!c
      qq=dsqrt(threehalves*norm2)
      I2=half*(norm2sig-I1*I1)
!c
!c ... components of (eta_d_ij)(eta_d_jk)
!c
      eta_d2(1)=eta_d(1)*eta_d(1)+eta_d(4)*eta_d(4)+eta_d(5)*eta_d(5)
      eta_d2(2)=eta_d(4)*eta_d(4)+eta_d(2)*eta_d(2)+eta_d(6)*eta_d(6)
      eta_d2(3)=eta_d(6)*eta_d(6)+eta_d(5)*eta_d(5)+eta_d(3)*eta_d(3)
      eta_d2(4)=eta_d(1)*eta_d(4)+eta_d(4)*eta_d(2)+eta_d(6)*eta_d(5)
      eta_d2(5)=eta_d(5)*eta_d(1)+eta_d(6)*eta_d(4)+eta_d(3)*eta_d(5)
      eta_d2(6)=eta_d(4)*eta_d(5)+eta_d(2)*eta_d(6)+eta_d(6)*eta_d(3)
!c           
!c ... Lode angle
!c
      if(norm2eta.lt.tiny) then 
!c
        cos3t=-one
!c               
      else
!c
        tretadev3=dot_vect_h(1,eta_d,eta_d2,6)!dot_vect_h(1,eta_d,eta_d2,6)
!c
        numer=-sqrt6*tretadev3
        denom=(dsqrt(norm2eta))**3
        cos3t=numer/denom
        if(dabs(cos3t).gt.one) then
             cos3t=cos3t/dabs(cos3t)
        end if
!c
      end if 
!c
!c ... determinant
!c
      xmin1=sig(2)*sig(3)-sig(6)*sig(6)
      xmin2=sig(4)*sig(3)-sig(6)*sig(5)
      xmin3=sig(4)*sig(6)-sig(5)*sig(2)
!c
      I3=sig(1)*xmin1-sig(4)*xmin2+sig(5)*xmin3

!c
      return
      end
!c------------------------------------------------------------------------------
      subroutine matmul_h(a,b,c,l,m,n)
!c------------------------------------------------------------------------------
!c matrix multiplication
!c------------------------------------------------------------------------------
      implicit none
!c
      integer i,j,k,l,m,n
!c
      double precision a(l,m),b(m,n),c(l,n)
!c
      do i=1,l
        do j=1,n
          c(i,j) = 0.0d0
          do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          enddo
        enddo
      enddo
!c
      return
      end
!c-----------------------------------------------------------------------------
      subroutine move_asv_h(asv,nasv,qq_n)
!c-----------------------------------------------------------------------------
!c move internal variables in vector qq_n and changes intergranular strain 
!c from continuum to soil mechanics convention
!c
!c NOTE: del has always 6 components
!c
!c written 6/2005 (Tamagnini, Sellari & Miriano)
!c-----------------------------------------------------------------------------
      implicit none
      integer nasv,i
      double precision asv(nasv),qq_n(nasv),zero 
!c
        parameter(zero=0.0d0)
!c
      do i=1,nasv
                qq_n(i)=zero
      enddo
!c
!c ... intergranular strain tensor stored in qq_n(1:6)
!c
      do i=1,6
                qq_n(i) = -asv(i)
      enddo
!c
!c ... void ratio stored in qq_n(7)
!c
        qq_n(7) = asv(7) 
!c
      return
      end
!c-----------------------------------------------------------------------------
      subroutine move_eps_h(dstran,ntens,deps,depsv)
!c-----------------------------------------------------------------------------
!c Move strain increment dstran into deps and computes 
!c volumetric strain increment
!c
!c NOTE: all strains negative in compression; deps has always 6 components
!c
!c written 7/2005 (Tamagnini, Sellari & Miriano)
!c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision deps(6),dstran(ntens),depsv
!c
      do i=1,ntens
                deps(i) = dstran(i)
      enddo
!c
        depsv=deps(1)+deps(2)+deps(3)
!c
      return
      end
!c-----------------------------------------------------------------------------
      subroutine move_sig_h(stress,ntens,pore,sig)
!c-----------------------------------------------------------------------------
!c computes effective stress from total stress (stress) and pore pressure (pore)
!c
!c NOTE: stress = total stress tensor (tension positive)
!c         pore   = exc. pore pressure (undrained conds., compression positive)
!c         sig    = effective stress (tension positive)
!c
!c       sig has always 6 components
!c
!c written 7/2005 (Tamagnini, Sellari & Miriano)
!c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision sig(6),stress(ntens),pore,zero 
!c
        parameter(zero=0.0d0)
!c
      do i=1,6
                sig(i)=zero
      enddo
!c
      do i=1,ntens
                if(i.le.3) then
                        sig(i) = stress(i)+pore
                else
                        sig(i) = stress(i)
                end if
      enddo
!c
      return
      end
!c-----------------------------------------------------------------------------
      subroutine norm_res_h(y_til,y_hat,ny,nasv,norm_R)
!c-----------------------------------------------------------------------------
!c  evaluate norm of residual vector Res=||y_hat-y_til||
!c
!c  written 6/2005 (Tamagnini, Sellari & Miriano)
!c-----------------------------------------------------------------------------
        implicit none
!c 
      integer ny,nasv,ng,k,i,testnan
!c
      double precision y_til(ny),y_hat(ny),void_til,void_hat,del_void
      double precision err(ny),norm_R2,norm_R
      double precision norm_sig2,norm_q2,norm_sig,norm_q
      double precision sig_hat(6),sig_til(6),del_sig(6)
      double precision q_hat(nasv),q_til(nasv),del_q(nasv)
      double precision zero!dot_vect_h,
!c
      parameter(zero=0.0d0)
!c
      ng=6*nasv
      k=42+nasv
!c
      do i=1,ny
              err(i)=zero
      end do
!c
!c ... recover stress tensor and internal variables
!c
      do i=1,6
                sig_hat(i)=y_hat(i)
                sig_til(i)=y_til(i)
                del_sig(i)=dabs(sig_hat(i)-sig_til(i))
      end do
!c
      do i=1,nasv-1
                q_hat(i)=y_hat(6+i)
                q_til(i)=y_til(6+i)
                del_q(i)=dabs(q_hat(i)-q_til(i))
      end do
!c
      void_hat=y_hat(6+nasv)
      void_til=y_til(6+nasv)
      del_void=dabs(void_hat-void_til)
!c
!c ... relative error norms
!c
      norm_sig2=dot_vect_h(1,sig_hat,sig_hat,6)!dot_vect_h(1,sig_hat,sig_hat,6)
      norm_q2=dot_vect_h(2,q_hat,q_hat,6)!dot_vect_h(2,q_hat,q_hat,6)
      norm_sig=dsqrt(norm_sig2)
      norm_q=dsqrt(norm_q2)
!c
      if(norm_sig.gt.zero) then
                do i=1,6
                        err(i)=del_sig(i)/norm_sig
                end do
      end if
!c
      if(norm_q.gt.zero) then
                do i=1,nasv-1
                err(6+i)=del_q(i)/norm_q
                end do
      end if
!c
      err(6+nasv)=del_void/void_hat
!c
!c ... global relative error norm
!c
      norm_R2=dot_vect_h(3,err,err,ny)!dot_vect_h(3,err,err,ny)
      norm_R=dsqrt(norm_R2)
!c
      testnan=0
      call umatisnan_h(norm_sig,testnan)
      call umatisnan_h(norm_q,testnan)
      call umatisnan_h(void_hat,testnan)
      if(testnan.eq.1) then
	norm_R=1.d20
      end if

      return
      end

!c-----------------------------------------------------------------------------
      subroutine perturbate_h(y_n,y_np1,n,nasv,dtsub,err_tol,maxnint,&
         DTmin,deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DD, dtime,&
         error)
!c-----------------------------------------------------------------------------
!c
!c  compute numerically consistent tangent stiffness
!c
!c  written 12/2005 (Tamagnini)
!c-----------------------------------------------------------------------------
      implicit none
!c 
      logical elprsw
!c
      integer ntens,jj,kk,i
      integer n,nasv,nparms,nfev
      integer maxnint,error
!c
      double precision y_n(n),y_np1(n),y_star(n),parms(nparms)
      double precision dtsub,err_tol,DTmin, dtime
      double precision theta,sig(6),q(nasv)
      double precision deps_np1(6),deps_star(6)
      double precision dsig(6),DD(6,6),HHtmp(nasv,6)
      double precision LL(6,6),NN(6)
      integer istrain
      double precision zero
!c
      parameter(zero=0.0d0)
!c
!c ... initialize DD and y_star
!c 
      if(parms(10) .le. 0.5) then
          istrain=0 
      else 
          istrain=1
      end if

      do kk=1,6
          do jj=1,6
              DD(kk,jj)=zero
          end do
      end do
      do i=1,6
          sig(i)=y_n(i)
      end do
      do i=1,nasv
          q(i)=y_n(6+i)
      end do
        
      call push_h(y_n,y_star,n)

      if(error.ne.10) then
          call get_tan_h(deps_np1,sig,q,nasv,parms,nparms,&
               	DD,HHtmp,LL,NN,istrain,error)                
        end if
        if(istrain .eq. 0) then
          do kk=1,6
                do jj=1,6
                       DD(kk,jj)=LL(kk,jj)
                end do
          end do
        else
          do kk=1,6
                do jj=1,6
                        DD(kk,jj)=parms(10)*LL(kk,jj)
                end do
          end do
        end if

        return
        end        
        
!c-----------------------------------------------------------------------------
      subroutine push_h(a,b,n)
!c-----------------------------------------------------------------------------
!c push vector a into vector b
!c-----------------------------------------------------------------------------
      implicit none
      integer i,n
      double precision a(n),b(n) 
!c
      do i=1,n
                b(i)=a(i)
      enddo
!c
      return
      end
!c-----------------------------------------------------------------------------
      subroutine rhs_h(y,ny,nasv,parms,nparms,deps,kRK,nfev,error)
!c-----------------------------------------------------------------------------
!c calculate coefficient kRK from current state y and strain increment deps
!c Masin hypoplastic model for clays with intergranular strains
!c
!c written 12/2005 (Tamagnini & Sellari)
!c-----------------------------------------------------------------------------
      implicit none
!c
        integer error,ny,nparms,nasv,i,nfev
!c
      double precision zero,one,two,four 
        double precision y(ny),kRK(ny),parms(nparms),deps(6)
        double precision sig(6),q(nasv)
        double precision F_sig(6),F_q(nasv)
!c
        parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)
!c
!c ... update counter for the number of function f(y) evaluations
!c
        nfev=nfev+1
!c
!c ... initialize kRK
!c
        do i=1,ny
                kRK(i)=zero
        end do
!c
!c ... recover current state variables (sig,q)                   
!c
        do i=1,6
                sig(i)=y(i)
        end do
!c
      do i=1,nasv
                q(i)=y(6+i)
        end do
!c       
!c ... build F_sig(6) and F_q(nasv) vectors and move them into kRK
!c
        call get_F_sig_q_h(sig,q,nasv,parms,nparms,deps,F_sig,F_q,error)
        if(error.eq.10) return
!c
        do i=1,6
!c
                kRK(i)=F_sig(i)
!c
        end do                   
!c       
        do i=1,nasv
!c
                kRK(6+i)=F_q(i)
!c
        end do                   
!c
      return
      end
!c-----------------------------------------------------------------------------
      subroutine rkf23_update_h(y,n,nasv,dtsub,err_tol,maxnint,DTmin,&
                             deps_np1,parms,nparms,nfev,elprsw,dtime,&
                             error)
!c-----------------------------------------------------------------------------
!c
!c  numerical solution of y'=f(y)
!c  explicit, adapive RKF23 scheme with local time step extrapolation
!c
!c  Tamagnini, Sellari & Miriano 6/2005
!c
!c-----------------------------------------------------------------------------
        implicit none
!c
        logical elprsw
!c
      integer n,nasv,nparms,i,ksubst,kreject,nfev
        integer maxnint,error,error_RKF
!c
      double precision y(n),parms(nparms),dtsub,err_tol,DTmin
        double precision zero,half,one,two,three,four,six
        double precision ptnine,onesixth,onethird,twothirds,temp
!c
        double precision deps_np1(6),y_k(n),y_2(n),y_3(n),y_til(n)
        double precision y_hat(n)
        double precision T_k,DT_k,dtime
        double precision kRK_1(n),kRK_2(n),kRK_3(n)
        double precision norm_R,S_hull
!c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(four=4.0d0,six=6.0d0,half=0.5d0,ptnine=0.9d0)
!c
!c ... initialize y_k vector and other variables
!c
        do i=1,n
                y_k(i)=zero
        end do
!c
        onesixth=one/six
        onethird=one/three
        twothirds=two/three
!c
!c ... start of update process
!c
                
        error_RKF=0
        T_k=zero      
        DT_k=dtsub/dtime
        ksubst=0
        kreject=0
        nfev=0
!c
        do i=1,n
                y_k(i)=y(i)
        end do
!c
!c ... start substepping 
!c
        do while(T_k.lt.one) 
!c
                ksubst=ksubst+1
!c
!c ... write substepping info
!c
!c               write(*,1234) ksubst,T_k,DT_k
!c1234           format('Substep no.',i4,' -- T_k = ',d12.4,' -- DT_k = ',d12.4)
!c
!c ... check for maximum number of substeps
!c
                if(ksubst.gt.maxnint) then
                       	write(1,*) 'number of substeps ',ksubst,&
                                  ' is too big, step rejected'
                        error=3
                        return
                end if          
!c
!c ... build RK functions
!c
                call check_RKF_h(error_RKF,y_k,n,nasv,parms,nparms)
                if(error_RKF.eq.1) then 
		  error=3
		  return
		else
		  call rhs_h(y_k,n,nasv,parms,nparms,deps_np1,kRK_1,nfev,error)
		end if
                if(error.eq.10) return
!c
!c ... find y_2
!c
                temp=half*DT_k
!c
                do i=1,n
                        y_2(i)=y_k(i)+temp*kRK_1(i)
                end do

!c               
                call check_RKF_h(error_RKF,y_2,n,nasv,parms,nparms)
                if(error_RKF.eq.1) then 
		  error=3
		  return
		else
		  call rhs_h(y_2,n,nasv,parms,nparms,deps_np1,kRK_2,nfev,error)
		end if
                if(error.eq.10) return
!c                                       
!c ... find y_3
!c

                do i=1,n
                        y_3(i)=y_k(i)-DT_k*kRK_1(i)+two*DT_k*kRK_2(i)
                end do
!c

                call check_RKF_h(error_RKF,y_3,n,nasv,parms,nparms)
                if(error_RKF.eq.1) then 
		  error=3
		  return
		else
		  call rhs_h(y_3,n,nasv,parms,nparms,deps_np1,kRK_3,nfev,error)
                end if
                if(error.eq.10) return

!c                               
!c ... approx. solutions of 2nd (y_til) and 3rd (y_hat) order
!c
                do i=1,n        
                        y_til(i)=y_k(i)+DT_k*kRK_2(i)
                        y_hat(i)=y_k(i)+DT_k*&
               (onesixth*kRK_1(i)+twothirds*kRK_2(i)+onesixth*kRK_3(i))
                end do
!c
!c ... local error estimate
!c

                call norm_res_h(y_til,y_hat,n,nasv,norm_R)
!c				check if output y_hat can be used as an input into the next step
                call check_RKF_h(error_RKF,y_hat,n,nasv,parms,nparms)

                if (error_RKF.ne.0) then
!c                	error=1.d20
!c                	error_RKF=0
			error=3
			return
                end if
!c
!c ... time step size estimator according to Hull
!c       	
		if(norm_R .ne. 0) then
                	S_hull=ptnine*DT_k*(err_tol/norm_R)**onethird
                else
                	S_hull=1
                end if
!c

      if (norm_R.lt.err_tol) then                             
!c
!c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
!c
                 do i=1,n        
                        y_k(i)=y_hat(i)
                 end do
!c
                        T_k=T_k+DT_k
                        DT_k=min(four*DT_k,S_hull)
                        dtsub=DT_k*dtime
                        DT_k=min((one-T_k),DT_k)        
!c
      else
!c
!c ... substep is not accepted, recompute with new (smaller) substep size DT
!c
                 DT_k=max(DT_k/four,S_hull)
!c
!c ... check for minimum step size
!c
                 if(DT_k.lt.DTmin) then
                              write(1,*) 'substep size ',DT_k,&
                                  ' is too small, step rejected'
                              error=3
                              return
                 end if          
!c                                       
      end if                                                  
!c
!c ... bottom of while loop
!c
      end do
        
!c
!c ... recover final state
!c
      do i=1,n
                y(i)=y_k(i)
      end do
!c
      return
      end
!c

!c-----------------------------------------------------------------------------
      subroutine check_RKF_h(error_RKF,y,ny,nasv,parms,nparms)
!c-----------------------------------------------------------------------------
!c Checks is RKF23 solout vector y is OK for hypoplasticity
!c-----------------------------------------------------------------------------
      implicit none
!c
        integer error_RKF,ny,nasv,i,nparms,testnan,iopt
!c
        double precision y(ny),parms(nparms)
        double precision sig(6),pmean,sig_star(6)
        double precision xN1(3),xN2(3),xN3(3),S(3),P,Q,tmin
        double precision p_t,minstress
!c
        p_t    =parms(2)
	minstress=p_t/4.d0
        do i=1,6
                sig(i)=y(i)
        end do

        sig_star(1)=sig(1)-p_t
        sig_star(2)=sig(2)-p_t
        sig_star(3)=sig(3)-p_t
        sig_star(4)=sig(4)
!c	changed order due to prnsig convention different from abaqus
        sig_star(5)=sig(6)
        sig_star(6)=sig(5)
                
    	pmean=-(sig_star(1)+sig_star(2)+sig_star(3))/3
    	
!c       check for positive mean stress
        if(pmean .le. minstress) then
        	error_RKF=1
        end if
!c
!c		calculate minimum principal stress
!c
		iopt=0
        Call PrnSig_h(iopt, sig_star, xN1, xN2, xN3,&
     		S(1),S(2),S(3), P, Q)
        tmin     = 1.0d+20
        do i=1,3
                if(tmin .ge. -S(i)) then
                	 tmin=-S(i)
                endif	 
        enddo 
!c
!c		check for tension
!c
        if(tmin .le. minstress) then
                error_RKF=1
        end if
        
!c		check for NAN
 	  	testnan=0
        do i=1,ny
       	  call umatisnan_h(y(i),testnan)
        end do
        call umatisnan_h(tmin,testnan)
        
        if(testnan.eq.1) error_RKF=1
!c
!c
      return
      end
!c

!c-----------------------------------------------------------------------------
      subroutine solout_h(stress,ntens,asv,nasv,ddsdde,y,nydim,&
                       pore,depsv_np1,parms,nparms,DD)
!c-----------------------------------------------------------------------------
!c copy the vector of state variables to umat output
!c modified 7/2005 (Tamagnini, Sellari)
!c
!c NOTE: solid mechanics convention for stress and strain components
!c       pore is always positive in compression
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer nydim,nasv,nparms,ntens,i,j
!c
      double precision y(nydim),asv(nasv),stress(ntens)
        double precision ddsdde(ntens,ntens),DD(6,6)
        double precision parms(nparms),bulk_w,pore,depsv_np1 
!c
        bulk_w=parms(15)
!c
!c ... update excess pore pressure (if undrained conditions), compression positive
!c
        pore=pore-bulk_w*depsv_np1
!c
!c updated total stresses (effective stresses stored in y(1:6))
!c
      do i=1,ntens
                if (i.le.3) then
                        stress(i) = y(i)-pore
                else
                        stress(i) = y(i)
                end if
        enddo
!c
!c additional state variables (first 6 components are intergranular strains)
!c
      do i=1,nasv
                asv(i) = y(6+i)
      enddo
!c
!c consistent tangent stiffness
!c
      do j=1,ntens
        do i=1,ntens
          ddsdde(i,j) = DD(i,j)      
        enddo
      enddo
!c
      do j=1,3
        do i=1,3
          ddsdde(i,j) = ddsdde(i,j)+bulk_w        
        enddo
      enddo
      return
      end
!c-----------------------------------------------------------------------------
      subroutine wrista_h(mode,y,nydim,deps_np1,dtime,coords,statev,&
                nstatv,parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
!c-----------------------------------------------------------------------------
!c ... subroutine for managing output messages
!c
!c     mode
!c
!c     all = writes:             kstep, kinc, noel, npt
!c       2   = writes also:      error message,coords(3),parms(nparms),ndi,nshr,stress(nstress)
!c                                               deps(nstress),dtime,statev(nstatv)
!c     3   = writes also:        stress(nstress),deps(nstress),dtime,statev(nstatv)
!c-----------------------------------------------------------------------------
      implicit none
!c
      integer mode,nydim,nstatv,nparms,noel,npt,ndi,nshr,kstep,kinc,i    
!c
      double precision y(nydim),statev(nstatv),parms(nparms)
        double precision deps_np1(6),coords(3),dtime
!c
!c ... writes for mode = 2
!c
      if (mode.eq.2) then
        write(1,*) '==================================================='
        write(1,*) 'ERROR: abaqus job failed during call of UMAT'
        write(1,*) '==================================================='
        write(1,*) 'state dump:'
        write(1,*) 
      endif
!c
!c ... writes for all mode values
!c
      write(1,111) 'Step: ',kstep, 'increment: ',kinc,&
      'element: ', noel, 'Integration point: ',npt
      write(1,*) 
!c
!c ... writes for mode = 2
!c
      if (mode.eq.2) then
        write(1,*) 'Co-ordinates of material point:'
        write(1,104) 'x1 = ',coords(1),' x2 = ',coords(2),' x3 = ',&
         coords(3)
        write(1,*) 
        write(1,*) 'Material parameters:'
        write(1,*) 
        do i=1,nparms
          write(1,105) 'prop(',i,') = ',parms(i)
        enddo 
        write(1,*)
        write(1,102) 'No. of mean components:  ',ndi
        write(1,102) 'No. of shear components: ',nshr
        write(1,*)
      endif
!c
!c ... writes for mode = 2 or 3
!c
      if ((mode.eq.2).or.(mode.eq.3)) then
        write(1,*) 'Stresses:'
        write(1,*) 
        write(1,101) 'sigma(1) = ',y(1)
        write(1,101) 'sigma(2) = ',y(2)
        write(1,101) 'sigma(3) = ',y(3)
        write(1,101) 'sigma(4) = ',y(4)
        write(1,101) 'sigma(5) = ',y(5)
        write(1,101) 'sigma(6) = ',y(6)
        write(1,*) 
        write(1,*) 'Strain increment:'
        write(1,*) 
        write(1,101) 'deps_np1(1) = ',deps_np1(1)
        write(1,101) 'deps_np1(2) = ',deps_np1(2)
        write(1,101) 'deps_np1(3) = ',deps_np1(3)
        write(1,101) 'deps_np1(4) = ',deps_np1(4)
        write(1,101) 'deps_np1(5) = ',deps_np1(5)
        write(1,101) 'deps_np1(6) = ',deps_np1(6)
        write(1,*) 
        write(1,*) 'Time increment:'
        write(1,*) 
        write(1,108) 'dtime = ',dtime
        write(1,*) 
        write(1,*) 'Internal variables:'
        write(1,*) 
        write(1,109) 'del(1) = ',statev(1)
        write(1,109) 'del(2) = ',statev(2)
        write(1,109) 'del(3) = ',statev(3)
        write(1,109) 'del(4) = ',statev(4)
        write(1,109) 'del(5) = ',statev(5)
        write(1,109) 'del(6) = ',statev(6)
        write(1,109) 'void   = ',statev(7)
        write(1,*) 
        write(1,*) '==================================================='
!c
      endif
!c
101   format(1X,a15,e11.4)
102   format(1X,a25,i1)
103   format(1X,a7,i5)
104   format(1X,3(a5,f10.4,2X))
105   format(1X,a5,i2,a4,f20.3)
106   format(1X,3(a9,f12.4,2X))
107   format(1X,3(a10,f12.4,2X))
108   format(1X,a8,f12.4)
109   format(1X,a6,f10.4)
110   format(1X,a5,f10.4)
111   format(1X,a6,i4,2X,a11,i4,2X,a9,i10,2X,a19,i4)
!c       
      return
      end

      
!c-----------------------------------------------------------------------------
      subroutine calc_statev_h(stress,statev,parms,nparms,nasv,&
      nstatv,deps)
!c-----------------------------------------------------------------------------
!c
!c  computes additional state variables for postprocessing
!c
!c-----------------------------------------------------------------------------
        implicit none
!c 
        logical elprsw
!c
      integer ntens,jj,kk,i
      integer n,nasv,nparms,nfev,nstatv
        integer maxnint,error
!c
      double precision parms(nparms)!,dot_vect_h
        double precision stress(6),statev(nstatv)
        double precision deps(6),tmax,tmin
        double precision MM(6,6),HHtmp(nasv,6)
        double precision LL(6,6),NN(6)
        integer istrain
        double precision zero,two,four,iopt,three
        double precision I1,I2,I3,cos3t,pp,qq
        double precision sin2phi,sinphi,sig_star(6),p_t
        double precision norm_del,norm_del2,del(6)
!c
      parameter(zero=0.0d0,two=2.0d0,four=4.0d0,three=3.0d0)
!c

!c ... calc phimob (statev 11) from Matsuoka-Nakai YS

      p_t    =parms(2)
      do i=1,3
              sig_star(i)=stress(i)-p_t
      end do
      do i=4,6
              sig_star(i)=stress(i)
      end do
      call inv_sig_h(sig_star,pp,qq,cos3t,I1,I2,I3)
	  if(I3 .ne. 0) then
        sin2phi=(9.d0+I1*I2/I3)/(1.d0+I1*I2/I3)
      else 
      	sin2phi=0
      end if
	  if(sin2phi .lt. 0) then
        sin2phi=0
      end if 
	  if(sin2phi .gt. 1) then
        sin2phi=1
      end if 
      sinphi=sqrt(sin2phi)
      
      statev(11)= asin(sinphi)*&
        180.0d0/3.141592d0

!c ... calc norm. length of intergr. strain rho (statev 12)
      if(parms(10) .le. 0.5) then
          istrain=0 
      else 
          istrain=1
      end if

      if(istrain .eq. 1) then
        
      do i=1,6
          del(i)=statev(i)
      enddo       
        
      norm_del2=dot_vect_h(2,del,del,6)!dot_vect_h(2,del,del,6)
      norm_del=dsqrt(norm_del2)
      statev(12)=norm_del/parms(12)
     
      else
        statev(12)=0
      end if

      return
      end        
            
!c-----------------------------------------------------------------------------
      subroutine umatisnan_h(chcknum,testnan)
!c-----------------------------------------------------------------------------
!c
!c  checks whether number is NaN
!c
!c-----------------------------------------------------------------------------
        double precision chcknum
        integer testnan

	    if (.not.(chcknum .ge. 0. .OR. chcknum .lt. 0.)) testnan=1        
	    if (chcknum .gt. 1.d30) testnan=1        
	    if (chcknum .lt. -1.d30) testnan=1        
 	    if (chcknum .ne. chcknum) testnan=1        
       
        return
        end         
      
!c-----------------------------------------------------------------------------
        subroutine xit_h
!c-----------------------------------------------------------------------------
        stop
!c
        return
        end

!C***********************************************************************
      Subroutine PrnSig_h(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension S(*),xN1(*),xN2(*),xN3(*)

      If (iOpt.Eq.1) Then
        Call Eig_3_h(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
      Else
        Call Eig_3a_h(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
      End If
      Return
      End
!C***********************************************************************
      Subroutine Eig_3_h(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3),V(3,3),&
               xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
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
      Do While ( it.Lt.itMax .And.&
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
          If (abs(a(ip,iq)) .gt. Tol) Then
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

      Subroutine Eig_3a_h(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      Implicit Double Precision (A-H,O-Z)
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

          If (abs(a(ip,iq)) .gt. Tol) Then         ! ongelijk nul ?
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
      
!c-----------------------------------------------------------------------------
      subroutine calc_elasti_h(y,n,nasv,dtsub,err_tol,maxnint,DTmin,&
                             deps_np1,parms,nparms,nfev,elprsw,&
     				dtime,DDtan,youngel,nuel,error)
!c-----------------------------------------------------------------------------
!c
!c  numerical solution of y'=f(y)
!c  explicit, adapive RKF23 scheme with local time step extrapolation
!c
!c  Tamagnini, Sellari & Miriano 6/2005
!c
!c-----------------------------------------------------------------------------
        implicit none
!c
        logical elprsw
!c
      integer n,nasv,nparms,i,ksubst,kreject,nfev
      integer maxnint,error,error_RKF,tension,j
!c
      double precision y(n),parms(nparms),dtsub,err_tol,DTmin
        double precision zero,half,one,two,three,four,six
        double precision ptnine,onesixth,onethird,twothirds,temp
!c
        double precision deps_np1(6),y_k(n),y_2(n),y_3(n),y_til(n)
        double precision y_hat(n),DDtan(6,6)
        double precision T_k,DT_k,dtime,II(6,6),krondelta(6)
        double precision kRK_1(n),kRK_2(n),kRK_3(n)
        double precision norm_R,S_hull,youngel,nuel,F_sig(6)
!c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(four=4.0d0,six=6.0d0,half=0.5d0,ptnine=0.9d0)
!c
!c ... initialize y_k vector and other variables
!c
        do i=1,n
                y_k(i)=zero
        end do
!c
        onesixth=one/six
        onethird=one/three
        twothirds=two/three

!c
!c ... fourth order identity tensors in Voigt notation
!c
        do i = 1,6
          do j=1,6
            II(i,j)=zero
          end do
        end do
        
        II(1,1)=one
        II(2,2)=one
        II(3,3)=one
        II(4,4)=half
        II(5,5)=half
        II(6,6)=half
!c
        krondelta(1)=one
        krondelta(2)=one
        krondelta(3)=one
        krondelta(4)=zero
        krondelta(5)=zero
        krondelta(6)=zero
!c
!c ... Elastic stiffness tensor 
!c
	if(youngel.gt.0) then
        do i = 1,6
          do j=1,6
            DDtan(i,j)=(youngel/(1+nuel))*(II(i,j) + &
            	nuel/(1-2*nuel)*krondelta(i)*krondelta(j));
          end do
        end do
        end if
        
        call matmul_h(DDtan,deps_np1,F_sig,6,6,1)
        do i=1,6
                y(i)=y(i)+F_sig(i)
        end do

        return
        end
!c

      
      
 
      
      
      
      
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
        
        
        
        
        
        
        
        
        
        
        
        
        
        

end module ModExternalSoilModel
