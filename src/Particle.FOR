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
    !	Copyright (C) 2023  Members of the Anura3D MPM Research Community 
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
	  
	  
	  module ModParticle
      !**********************************************************************
      !
      !    Function:  Definition of a TYPE (class) ParticleType which forms the basis
      !               of the particle housekeeping. The data structure that stores
      !               the information of each individual particle contains items
      !               of this type.
      !
      !               This module should only be extended by variables that are missing in ParticleType.
      !               No routines should be added unless they solely involve manipulation of data stored
      !               in ParticleType.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      
      use ModGlobalConstants
      use ModReadMaterialData
      
      implicit none
      
        !> Stores the maximum assigned particle ID. This variable should only be modified through the routine InitParticle.
        integer(INTEGER_TYPE), private :: IDCounter = 0

        type ParticleType(tsize, vsize,LoadSystems)
          integer(INTEGER_TYPE), len :: tsize, vsize,LoadSystems 
          
          integer(INTEGER_TYPE) :: Kind = ANYPARTICLE
          !> Particle mass [kg] (corresponds to the gas density)(fixed during simulation - ensures conservation of mass)
          real(REAL_TYPE) :: MassGas = -1 ! Particle gas mass [kg]
          !> Particle mass [kg] (corresponds to the saturated density density)
          !> (fixed during simulation - ensures conservation of mass)
          real(REAL_TYPE) :: MassMixed = -1 ! Particle saturated mass [kg]
          real(REAL_TYPE) :: MaterialWeight = -1 ! Gamma (dry)
          real(REAL_TYPE) :: WaterWeight = -1 ! Gamma (water)
          real(REAL_TYPE) :: GasWeight = -1 ! Gamma (gas)
          real(REAL_TYPE) :: SatWeight = -1 ! Gamma (saturated)
          real(REAL_TYPE) :: UnsatWeight = -1 ! Gamma (unsaturated)
          real(REAL_TYPE) :: MixedWeight = -1 ! Gamma (saturated)
          real(REAL_TYPE) :: Conductivity = -1 ! Hydraulic Conductivity
          real(REAL_TYPE) :: ConductivityGas = -1 ! Hydraulic Conductivity of gas
          real(REAL_TYPE) :: Porosity = -1 ! Porosity
          real(REAL_TYPE) :: InitialPorosity = -1 ! Initial porosity
          real(REAL_TYPE) :: DegreeSaturation = -1 ! Degree of saturation (Sr)
          real(REAL_TYPE) :: BulkWater = -1 ! Bulk Modulus of water
          real(REAL_TYPE) :: BulkGas = -1 ! Bulk Modulus of gas
          real(REAL_TYPE) :: Density = -1 ! Soil dry Density (changes during simulation based on volume changes)
          real(REAL_TYPE) :: ConstDensity = -1 ! Initial density of solid and liquid
          real(REAL_TYPE) :: AirInWaterMassFraction = 0.0 ! Mass Fraction (Air in the liquid)
          real(REAL_TYPE) :: VapourInGasMassFraction = 0.0 ! Mass Fraction (Vapour in the gas)
          real(REAL_TYPE) :: DiffusionAirInWater = 0.0 ! Diffusion of the air in the liquid
          real(REAL_TYPE) :: DiffusionVapourInGas = 0.0 ! Diffusion of the vapour in the gas
          real(REAL_TYPE) :: WaterAdvectiveFlux = 0.0 ! Water Advective flux
          real(REAL_TYPE) :: AirAdvectiveFlux = 0.0 ! Air Advective flux
          !> Local 'volume' occupied by the particle (changes during simulation due to extension/compression of material)
          real(REAL_TYPE) :: IntegrationWeight = -1 ! Particle volume
          real(REAL_TYPE) :: ShearModulus = 0.0 ! Shear modulus of the particle
          real(REAL_TYPE) :: CohesionCosPhi = 0.0 ! Cohesion of the particle multiplied by Cos Phi (friction angle)
          real(REAL_TYPE) :: SFail = 0.0 ! CohesionCosPhi divided by Sin Phi
          real(REAL_TYPE) :: CohesionStSoft = 0.0 ! Cohesion of the particle (model Strain Softening MC)
          real(REAL_TYPE) :: PhiStSoft = 0.0 ! Friction angle of the particle (model Strain Softening MC)
          real(REAL_TYPE) :: PsiStSoft = 0.0 ! Dilatancy angle of the particle (model Strain Softening MC)
          
          integer(INTEGER_TYPE) :: IPL = 0 ! Plasticity state of the particle
          real(REAL_TYPE) :: EoedHP = 0.0 
          real(REAL_TYPE) :: Alpha = 0.0
          real(REAL_TYPE) :: Beta = 0.0
          
          ! One-layer liquid formulation
          real(REAL_TYPE) :: LiquidFreeSurface = 0.0
          real(REAL_TYPE) :: LiquidFreeSurfaceCumul = 0.0
          
          ! Two-layer formulation
          real(REAL_TYPE) :: ConcentrationRatioLiquidL = 0.0
          real(REAL_TYPE) :: ConcentrationRatioSolidL = 0.0
          real(REAL_TYPE) :: PorosityL = 0.0
          real(REAL_TYPE) :: EffConcentrationRatioSolid = 0.0
          real(REAL_TYPE) :: EffConcentrationRatioLiquid = 0.0
          real(REAL_TYPE) :: EffPorosity = 0.0
          real(REAL_TYPE) :: ConcentrationRatioLiquidS = 0.0
          real(REAL_TYPE) :: FillingRatio = 1.0
          
          integer(INTEGER_TYPE) :: PhaseStatus = PhaseStatusUNDEFINED
          integer(INTEGER_TYPE) :: ContainedMaterialTypes = ContainedMaterialTypeUNDEFINED
          
          ! Applied boundary conditions (loads, prescribed displacements)
          real(REAL_TYPE), dimension(vsize) :: FBody != 0.0 ! Body forces assigned to the particle (calculated from soil dry weight)
          real(REAL_TYPE), dimension(vsize) :: FBodyWater != 0.0 ! Body forces assigned to the particle (calculated from water weight)
          real(REAL_TYPE), dimension(vsize) :: FBodyGas != 0.0 ! Body forces assigned to the particle (calculated from gas weight)
          !> Body forces assigned to the particle (calculated from saturated weight)
          real(REAL_TYPE), dimension(vsize) :: FBodyMixed != 0.0
          !real(REAL_TYPE), dimension(vsize) :: FExt != 0.0 ! External loads assigned to the particle (corresponds to total stress)
          real(REAL_TYPE), dimension(vsize,LoadSystems) :: FExt != 0.0 ! External loads assigned to the particle (corresponds to total stress)
          !> External loads assigned to the particle (corresponds to water pressure)
          real(REAL_TYPE), dimension(vsize,LoadSystems) :: FExtWater != 0.0
          !> External loads assigned to the particle (corresponds to gas pressure)
          real(REAL_TYPE), dimension(vsize,LoadSystems) :: FExtGas != 0.0
          !prescribed velocity or displacement
          real(REAL_TYPE), dimension(vsize) :: PrescrDisp != 0.0 ! Prescribed particle displacements
          real(REAL_TYPE), dimension(vsize) :: PrescrVelo != 0.0 ! Prescribed particle displacements
          logical, dimension(vsize)  :: MPPrescribedVeloDir !true if a prescribed velocity is assigned
          

          ! Location
          real(REAL_TYPE), dimension(vsize) :: LocPos != 0.0 ! Local particle position xi = LocPos(1), eta = LocPos(2), zeta = LocPos(3)
          logical :: IsBoundaryParticle = .false. ! True, if the particle is located on the boundary of the material body
      
          logical :: IsHydrHeadBoundaryParticle = .false. ! True, if the particle is located on the boundary of the material body where also hydr head is prescribed

          ! Values specific to the Dynamic MPM
          real(REAL_TYPE), dimension(vsize) :: AccelerationWater != 0.0

          ! Displacements (water and gas)
          real(REAL_TYPE), dimension(vsize) :: UW != 0.0 ! Total displacements in x, y, z direction of the water
          real(REAL_TYPE), dimension(vsize) :: UG != 0.0 ! Total displacements in x, y, z direction of the gas
          real(REAL_TYPE), dimension(vsize) :: UStepWater != 0.0

          ! Strains
          !> Total strains (Exx, Eyy, Ezz, Gxy, Gyz, Gxz) of the particle
          real(REAL_TYPE), dimension(tsize) :: Eps ! = 0.0
          !> Total plastic strains (Exx, Eyy, Ezz, Gxy, Gyz, Gxz) of the particle
          real(REAL_TYPE), dimension(tsize) :: EpsP != 0.0
          !> Incremental strains (Exx, Eyy, Ezz, Gxy, Gyz, Gxz) of the particle during a load step
          real(REAL_TYPE), dimension(tsize) :: EpsStep != 0.0
          real(REAL_TYPE), dimension(tsize) :: EpsStepPrevious != 0.0
          !> Incremental strains (Exx, Eyy, Ezz, Gxy, Gyz, Gxz) of the particle during a load step
          !real(REAL_TYPE), dimension(:), allocatable :: EpsStepW 
          !> Strain Rates (dExx, dEyy, dEzz, dGxy, dGyz, dGxz) of the particle during a load step
          real(REAL_TYPE) :: WaterVolumetricStrain = 0.0 ! The volumetric strain (water phase)
          real(REAL_TYPE) :: GasVolumetricStrain = 0.0 ! The volumetric strain (gas phase)

          ! Stresses
          real(REAL_TYPE), dimension(tsize) :: SigmaEffStep != 0.0
          !> Sima principal (S1, S2, S3) of the particle during a load step
          real(REAL_TYPE), dimension(tsize) :: SigmaPrin != 0.0
          real(REAL_TYPE) :: WaterPressure = 0.0 ! Water pressure of the particle
          real(REAL_TYPE) :: WaterPressure0 = 0.0 ! Initial water pressure of the particle from the last load step
          real(REAL_TYPE) :: HydraulicHead = 0.0 ! 
          
          real(REAL_TYPE) :: GasPressure = 0.0 ! Gas pressure of the particle  (100kPa, Atmospheric pressure)
          !> Initial gas pressure of the particle from the last load step  (100kPa, Atmospheric pressure)
          real(REAL_TYPE) :: GasPressure0 = 0.0

          ! Heat parameters
          real(REAL_TYPE) :: Temperature = 20.0 ! Temperature of the particle
          real(REAL_TYPE) :: ThermalCapacity = 0.0 ! Thermal capacity of the particle

          ! Damping
          real(REAL_TYPE) :: Damping = 0.0 ! particle damping
          real(REAL_TYPE) :: DBulkViscousPressure = 0.0 ! bulk viscosity damping pressure


          ! HP Model
          real(REAL_TYPE), dimension(2) :: HPStateVariables          = 0.0 ! the state variables of HP model
          real(REAL_TYPE), dimension(2) :: ModifiedHPStateVariables  = 0.0 ! the state variables of the modified HP model
          real(REAL_TYPE), dimension(7) :: HPIGStateVariables        = 0.0 ! the state variables of HP intergranular strain model

          ! MCC model
          !real(REAL_TYPE) :: lambda = 0.0 ! slope of the CSL and NCL in e-ln(p') plane (model MCC)
          !real(REAL_TYPE) :: kappa = 0.0 ! slope of the unloading/reloading line in e-ln(p') plane (model MCC)
          !real(REAL_TYPE) :: M = 0.0 ! slope of the CSL in q-p' plane (model MCC)
          !real(REAL_TYPE) :: nu = 0.0 ! poisson ration (model MCC)
          real(REAL_TYPE) :: init_void_ratio = 0.0 ! initial void ration ration (model MCC)
          real(REAL_TYPE) :: pp = 0.0 ! the preconsolidation pressure (model MCC)

          ! User defined model
          !> equivalent unloading stiffness, needed for calculation of wave speed to determine critical time step size
          real(REAL_TYPE) :: ESM_UnloadingStiffness

          ! Arbitrary data (for testing)
          real(REAL_TYPE), dimension(vsize) :: VariableData != 0.0 ! Arbitrary data, dimension 3

        end type ParticleType

    contains ! Routines of this module
    
        integer(INTEGER_TYPE) function GetIDCounter()
        !**********************************************************************
        !
        !    Function:  Returns the counter which stores the currently assigned
        !               maximum particle ID. This value does not necessarily
        !               correspond to the number of particles.
        !
        ! O GetIDCounter : Counter of particle ID's
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          GetIDCounter = IDCounter
        
        end function GetIDCounter


        subroutine SetIDCounter(NewIDCounter)
        !**********************************************************************
        !
        !    Function:  Set the counter which stores the currently assigned
        !               maximum particle ID. This value does not necessarily
        !               correspond to the number of particles.
        !
        !    NewIDCounter : New value of the ID counter
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NewIDCounter

          IDCounter = NewIDCounter
        
        end subroutine SetIDCounter


         
        subroutine SetParticleStructureDefault(Particle)
        !**********************************************************************
        !
        !    SUBROUTINE:  SetParticleStructureDefault
        !>   Sets the particle structure to default values
        !
        !>   @param[inout] Particle: particle structure    
        !
        !**********************************************************************
        implicit none
        ! Local variables
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle

          Particle%MassGas = -1 ! Pasticle gas mass [kg]
          ! Particle mass [kg] (corresponds to the saturated density density)
          ! (fixed during simulation - ensures conservation of mass)
          Particle%MassMixed = -1 ! Particle saturated mass [kg]
          Particle%MaterialWeight = -1 ! Gamma (dry)
          Particle%WaterWeight = -1 ! Gamma (water)
          Particle%GasWeight = -1 ! Gamma (gas)
          Particle%SatWeight = -1 ! Gamma (saturated)
          Particle%UnsatWeight = -1 ! Gamma (unsaturated)
          Particle%MixedWeight = -1 ! Gamma (saturated)
          Particle%Conductivity = -1 ! Hydraulic Conductivity
          Particle%ConductivityGas = -1 ! Hydraulic Conductivity of gas
          Particle%Porosity = -1 ! Porosity
          Particle%ESM_UnloadingStiffness = 0.0 ! Eunloading
          Particle%InitialPorosity = -1 ! Initial porosity
          Particle%DegreeSaturation = -1 ! Degree of saturation (Sr)
          Particle%BulkWater = -1 ! Bulk Modulus of water
          Particle%BulkGas = -1 ! Bulk Modulus of gas
          Particle%Density = -1 ! Soil dry Density (changes during simulation based on volume changes)
          Particle%ConstDensity = -1 ! Initial density of solid and liquid
          Particle%AirInWaterMassFraction = 0.0 ! Mass Fraction (Air in the liquid)
          Particle%VapourInGasMassFraction = 0.0 ! Mass Fraction (Vapour in the gas)
          Particle%DiffusionAirInWater = 0.0 ! Diffusion of the air in the liquid
          Particle%DiffusionVapourInGas = 0.0 ! Diffusion of the vapour in the gas
          Particle%WaterAdvectiveFlux = 0.0 ! Water Advective flux
          Particle%AirAdvectiveFlux = 0.0 ! Air Advective flux
          ! Local 'volume' occupied by the particle (changes during simulation due to extension/compression of material)
          Particle%IntegrationWeight = -1 ! Particle volume
          Particle%ShearModulus = 0.0 ! Shear modulus of the particle
          Particle%CohesionCosPhi = 0.0 ! Cohesion of the particle multiplied by Cos Phi (friction angle)
          Particle%SFail = 0.0 ! CohesionCosPhi divided by Sin Phi
          Particle%CohesionStSoft = 0.0 ! Cohesion of the particle (model Strain Softening MC)
          Particle%PhiStSoft = 0.0 ! Friction angle of the particle (model Strain Softening MC)
          Particle%PsiStSoft = 0.0 ! Dilatancy angle of the particle (model Strain Softening MC)
 
          Particle%IPL = 0 ! Plasticity state of the particle
          Particle%EoedHP = 0.0 
          Particle%Alpha = 0.0
          Particle%Beta = 0.0
          
          ! One-layer liquid formulation
          Particle%LiquidFreeSurface = 0.0
          Particle%LiquidFreeSurfaceCumul = 0.0
          
          ! Two-layer formulation
          Particle%ConcentrationRatioLiquidL = 0.0
          Particle%ConcentrationRatioSolidL = 0.0
          Particle%PorosityL = 0.0
          Particle%EffConcentrationRatioSolid = 0.0
          Particle%EffConcentrationRatioLiquid = 0.0
          Particle%EffPorosity = 0.0
          Particle%ConcentrationRatioLiquidS = 0.0
          Particle%FillingRatio = 1.0
          
          Particle%PhaseStatus = PhaseStatusUNDEFINED
          Particle%ContainedMaterialTypes = ContainedMaterialTypeUNDEFINED
          
          ! Applied boundary conditions (loads, prescribed displacements)
          Particle%FBody = 0.0 ! Body forces assigned to the particle (calculated from soil dry weight)
          Particle%FBodyWater = 0.0 ! Body forces assigned to the particle (calculated from water weight)
          Particle%FBodyGas = 0.0 ! Body forces assigned to the particle (calculated from gas weight)
          !> Body forces assigned to the particle (calculated from saturated weight)
          Particle%FBodyMixed = 0.0
          Particle%FExt = 0.0 ! External loads assigned to the particle (corresponds to total stress)
          !> External loads assigned to the particle (corresponds to water pressure)
          Particle%FExtWater = 0.0
          !> External loads assigned to the particle (corresponds to gas pressure)
          Particle%FExtGas = 0.0
          Particle%PrescrDisp = 0.0 ! Prescribed particle displacements
          Particle%PrescrVelo = 0.0 ! Prescribed particle velocity
          Particle%MPPrescribedVeloDir = .false. !true if a prescribed velocity is assigned

          ! Location
          Particle%LocPos = 0.0 ! Local particle position xi = LocPos(1), eta = LocPos(2), zeta = LocPos(3)
          Particle%IsBoundaryParticle = .false. ! True, if the particle is located on the boundary of the material body
      
          Particle%IsHydrHeadBoundaryParticle = .false. !True, if the particle is located on the boundary of the material body 
          ! and belongs to a border where hydraulic head is prescribed

          ! Values specific to the Dynamic MPM
          Particle%AccelerationWater = 0.0

          ! Displacements (water and gas)
          Particle%UW = 0.0 ! Total displacements in x, y, z direction of the water
          Particle%UG = 0.0 ! Total displacements in x, y, z direction of the gas
          Particle%UStepWater = 0.0

          ! Strains
          ! Total strains (Exx, Eyy, Ezz, Gxy, Gyz, Gxz) of the particle
          Particle%Eps  = 0.0
          ! Total plastic strains (Exx, Eyy, Ezz, Gxy, Gyz, Gxz) of the particle
          Particle%EpsP = 0.0
          ! Incremental strains (Exx, Eyy, Ezz, Gxy, Gyz, Gxz) of the particle during a load step
          Particle%EpsStep = 0.0
          Particle%EpsStepPrevious = 0.0


          ! Strain Rates (dExx, dEyy, dEzz, dGxy, dGyz, dGxz) of the particle during a load step
          Particle%WaterVolumetricStrain = 0.0 ! The volumetric strain (water phase)
          Particle%GasVolumetricStrain = 0.0 ! The volumetric strain (gas phase)

          ! Stresses
          Particle%SigmaEffStep = 0.0
          !> Sima principal (S1, S2, S3) of the particle during a load step
          Particle%SigmaPrin = 0.0
          Particle%WaterPressure = 0.0 ! Water pressure of the particle
          Particle%WaterPressure0 = 0.0 ! Initial water pressure of the particle from the last load step
          Particle%HydraulicHead = 0.0 ! 
          
          Particle%GasPressure = 0.0 ! Gas pressure of the particle  (100kPa, Atmospheric pressure)
          !> Initial gas pressure of the particle from the last load step  (100kPa, Atmospheric pressure)
          Particle%GasPressure0 = 0.0

          ! Heat parameters
          Particle%Temperature = 20.0 ! Temperature of the particle
          Particle%ThermalCapacity = 0.0 ! Thermal capacity of the particle

          ! Damping
          Particle%Damping = 0.0 ! particle damping
          Particle%DBulkViscousPressure = 0.0 ! bulk viscosity damping pressure


          ! HP Model
          Particle%HPStateVariables          = 0.0 ! the state variables of HP model
          Particle%ModifiedHPStateVariables  = 0.0 ! the state variables of the modified HP model
          Particle%HPIGStateVariables        = 0.0 ! the state variables of HP intergranular strain model

          ! MCC model
          !Particle%lambda = 0.0 ! slope of the CSL and NCL in e-ln(p') plane (model MCC)
          !Particle%kappa = 0.0 ! slope of the unloading/reloading line in e-ln(p') plane (model MCC)
          !Particle%M = 0.0 ! slope of the CSL in q-p' plane (model MCC)
          !Particle%nu = 0.0 ! poisson ration (model MCC)
          Particle%init_void_ratio = 0.0 ! initial void ration ration (model MCC)
          Particle%pp = 0.0 ! the preconsolidation pressure (model MCC)



          ! Arbitrary data (for testing)
          Particle%VariableData = 0.0 ! Arbitrary data, dimension 3
          
        end subroutine SetParticleStructureDefault
        
        
        type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)) function InitParticle(MaterialID, &
                                                 MaterialPointType, &
                                                 IntegrationWeight, &
                                                 LocPos, &
                                                 Kind, &
                                                 VirtualParticleMass, &
                                                 GGrav, &
                                                 GravityVector, &
                                                 DampingFactor, &
                                                 ISSubmerged, &
                                                 ParticleMassValue, &
                                                 ParticleMassWater)
        !**********************************************************************
        !
        !     Function:  Allocates a particle and assigns a unique ID to the particle.
        !                Increases IDCounter by one.
        !
        !     ElementID : ID of the element where the particle is initially located
        !     EntityID : ID of the entity to which the particle belongs
        !     MaterialID : ID of the material of the particle
        !     MaterialPointType : MIXTURE, SOLID, LIQUID
        !     IntegrationWeight : Integration weight assigned to the particle
        !     MaterialWeight : Gamma(dry)
        !     WaterWeight : Gamma (water)
        !     GasWeight : Gamma (gas)
        !     UnsatWeight : Gamma (unsaturated)
        !     SatWeight : Gamma (saturated)
        !     MixedWeight : Gamma mixed (saturated or unsaturated)
        !     Conductivity : Hydraulic Conductivity of water
        !     ConductivityGas : Hydraulic Conductivity of gas
        !     Porosity : Porosity
        !     InitialPorosity : Initial porosity
        !     DegreeSaturation : Degree of saturation
        !     BulkWater : Bulk modulus of water
        !     BulkGas : Bulk modulus of gas
        !     AirInWaterMassFraction: mass fraction of the air in the liquid phase
        !     VapourInGasMassFraction: mass fraction of the vapour in the gas phase
        !     DiffusionAirInWater : Diffusion of the air in the liquid
        !     DiffusionVapourInGas : Diffusion of the vapour in the gas
        !     ShearModulus : Shear modulus of the particle
        !     RefCohesion : Cohesion of the particle
        !     SinPhi : Sinus Phi (friction angle)
        !     LocPos : Location of the particle inside the reference coordinate systems of the element
        !     Kind : Material or virtual particle
        !     VirtualParticleMass : Mass assigned to particle in case of virtual particles
        !     GGrav : Gravity acceleration
        !     GravityVector : Direction-vector of gravitational force
        !     DampingFactor : The damping factor of the element where the particle is initially located
        !     IElTyp : Number of node connectivities
        !
        ! O   InitParticle : Returns a new particle
        !
        !**********************************************************************
                
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: Kind
          integer(INTEGER_TYPE), intent(in) :: MaterialID
          integer(INTEGER_TYPE), intent(in) :: MaterialPointType
          real(REAL_TYPE), intent(in) :: DampingFactor
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: LocPos, GravityVector
          real(REAL_TYPE), intent(in) :: GGrav, IntegrationWeight
          logical :: ISSubmerged, IsUndrEffectiveStress, IsUndrTotalStress
          real(REAL_TYPE), intent(in) :: VirtualParticleMass
          real(REAL_TYPE), intent(inout) ::ParticleMassValue, ParticleMassWater
          
          ! Local variables
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)) :: NewParticle
          real(REAL_TYPE) :: MaterialWeight, WaterWeight, GasWeight, UnsatWeight, SatWeight, MixedWeight,  &
                              Conductivity, ConductivityGas, InitialPorosity, DegreeSaturation, BulkWater, BulkGas, &
                              ShearModulus, RefCohesion, SinPhi, CohesionStSoft, PhiStSoft, PsiStSoft, WaterWeightPorosity
          integer(INTEGER_TYPE) :: MPIDCounter
          
          call SetParticleStructureDefault(NewParticle)
          
          MPIDCounter = GetIDCounter() + 1
          call SetIDCounter(MPIDCounter)
          
          NewParticle%Kind = Kind
          
          IsUndrEffectiveStress = &
                !code version 2016 and previous
                ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialID)%MaterialType)=='2-phase')) .or. &
                !code version 2017.1 and following
                (trim(MatParams(MaterialID)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))

          !code version 2017.1 and following
          IsUndrTotalStress = trim(MatParams(MaterialID)%MaterialType)==SATURATED_SOIL_UNDRAINED_TOTAL

          if (IsUndrEffectiveStress.or.IsUndrTotalStress) then
            MaterialWeight = MatParams(MaterialID)%WeightMixture ! total weight of the mixture
          else
            MaterialWeight = MatParams(MaterialID)%DryWeight ! dry weight of the solid
          end if

          WaterWeight    = MatParams(MaterialID)%WeightLiquid ! weight liquid
          GasWeight      = MatParams(MaterialID)%WeightGas ! weight gas
          UnsatWeight    = MatParams(MaterialID)%WeightMixture ! unsaturated weight
          SatWeight      = MatParams(MaterialID)%WeightMixture ! saturated weight
          MixedWeight    = MatParams(MaterialID)%WeightMixture ! mixed weight (saturated or unsaturated)
          Conductivity   = MatParams(MaterialID)%HydraulicConductivityLiquid ! Hydraulic conductivity liquid
          ConductivityGas = MatParams(MaterialID)%HydraulicConductivityGas ! Hydraulic conductivity of gas
          InitialPorosity = MatParams(MaterialID)%InitialPorosity! initial porosity
          DegreeSaturation = MatParams(MaterialID)%InitialDegreeOfSaturation ! initial degree of saturation
          BulkWater      = MatParams(MaterialID)%BulkModulusLiquid ! Bulk modulus of water
          BulkGas        = MatParams(MaterialID)%BulkModulusGas ! Bulk modulus of gas
          ShearModulus   = MatParams(MaterialID)%ShearModulus ! Shear modulus
          RefCohesion    = MatParams(MaterialID)%Cohesion ! Reference cohesion
          SinPhi         = SIN(MatParams(MaterialID)%FrictionAngle*(Pi/180.0)) ! Sinus Phi (friction angle)
                
          CohesionStSoft = MatParams(MaterialID)%PeakCohesion ! cohesion (depending on the softening parameters)
          ! friction angle (depending on the softening parameters)
          PhiStSoft      = MatParams(MaterialID)%PeakFrictionAngle*Pi/180.0
          ! dilatancy angle (depending on the softening parameters)
          PsiStSoft      = MatParams(MaterialID)%PeakDilatancyAngle*Pi/180.0
          
          NewParticle%IntegrationWeight = IntegrationWeight
          NewParticle%Damping = DampingFactor 
          NewParticle%MaterialWeight = MaterialWeight 
          NewParticle%WaterWeight = WaterWeight
          NewParticle%GasWeight = GasWeight
          NewParticle%UnsatWeight = UnsatWeight 
          NewParticle%SatWeight = SatWeight 
          NewParticle%MixedWeight = MixedWeight 
          if (NFORMULATION==1) then
              NewParticle%Conductivity = Conductivity
              NewParticle%ConductivityGas = ConductivityGas
          else
               ! 2 Constituents and 1 Phase
          !    if(MaterialPointType=='SOLID') then
                  NewParticle%Conductivity = Conductivity
                  NewParticle%ConductivityGas = ConductivityGas
              !end if
              !if(MaterialPointType=='LIQUID') then
              !    NewParticle%Conductivity = 0.0
              !    NewParticle%ConductivityGas = 0.0
              !end if
          end if
          NewParticle%Porosity = InitialPorosity
          NewParticle%InitialPorosity = InitialPorosity
          NewParticle%DegreeSaturation = DegreeSaturation
          NewParticle%BulkWater = BulkWater
          NewParticle%BulkGas = BulkGas
          NewParticle%AirInWaterMassFraction = 0.0
          NewParticle%VapourInGasMassFraction = 0.0
          NewParticle%LocPos = LocPos
          if (NFORMULATION==1) then
              NewParticle%Density = NewParticle%MaterialWeight / GGrav !......Soil dry density
              NewParticle%ConstDensity = NewParticle%Density
          else
              ! 2 Constituents
              if(MaterialPointType==MaterialPointTypeSolid) then
                !......solid unit density
                NewParticle%ConstDensity = NewParticle%MaterialWeight / ((1 - NewParticle%InitialPorosity) * GGrav)
                NewParticle%Density = NewParticle%ConstDensity*(1 - NewParticle%InitialPorosity)
              end if
              if(MaterialPointType==MaterialPointTypeLiquid) then
                NewParticle%ConstDensity = NewParticle%WaterWeight /  GGrav !......liquid unit density
                NewParticle%Density = NewParticle%ConstDensity*NewParticle%InitialPorosity
              end if

          end if

          NewParticle%ShearModulus = ShearModulus
          NewParticle%CohesionCosPhi = RefCohesion * dsqrt(1D0 - SinPhi * SinPhi)
          if (SinPhi>1D-10) then
            NewParticle%SFail = NewParticle%CohesionCosPhi / SinPhi
          else
            NewParticle%SFail = 1D10
          end if

          ! Strain softening parameters MC
          NewParticle%CohesionStSoft = CohesionStSoft
          NewParticle%PhiStSoft = PhiStSoft
          NewParticle%PsiStSoft = PsiStSoft 

          ! MCC parameters
          NewParticle%Init_void_ratio = MatParams(MaterialID)%InitialVoidRatio
          NewParticle%pp = MatParams(MaterialID)%MCC_PC
          
          ! Determine the mass of the particle
          if (NFORMULATION==1) then
              ParticleMassValue = DeterminePartMass(MaterialWeight, GGrav, IntegrationWeight)
          else
              ! 2 Constituents
              if(MaterialPointType==MaterialPointTypeSolid) then
                  ParticleMassValue = DeterminePartMass(MaterialWeight, GGrav, IntegrationWeight)
              end if
              if(MaterialPointType==MaterialPointTypeLiquid) then
                  WaterWeightPorosity = InitialPorosity * WaterWeight
                  ! mass of liquid per unit volume of mixture
                  ParticleMassValue = DeterminePartMass(WaterWeightPorosity, GGrav, IntegrationWeight)
              end if
          end if

          if (Kind==NEWVIRTUALPARTICLE) then ! Use specified weight
            ParticleMassValue = VirtualParticleMass
            NewParticle%MaterialWeight = 0.0
          end if

          ParticleMassWater = DeterminePartMass(WaterWeight, GGrav, IntegrationWeight)
          NewParticle%MassGas = DeterminePartMass(GasWeight, GGrav, IntegrationWeight)
          NewParticle%MassMixed = DeterminePartMass(MixedWeight, GGrav, IntegrationWeight)

          if (NFORMULATION==1) then
              
            ! Determine dry weight
            NewParticle%FBody(:) = MaterialWeight * IntegrationWeight * GravityVector(:) 
            ! Determine water weight
            NewParticle%FBodyWater(:) = WaterWeight * IntegrationWeight * GravityVector(:) 
            ! Determine gas weight
            NewParticle%FBodyGas(:) = GasWeight * IntegrationWeight * GravityVector(:) 

          else ! 2 Constituents
            
            if (MaterialPointType==MaterialPointTypeSolid) then
                
              NewParticle%FBody(:) = (1.0 - NewParticle%Porosity) * NewParticle%ConstDensity * IntegrationWeight * CalParams%GravityData%GAccel * GravityVector(:)
              NewParticle%FBodyWater = 0.0 
              NewParticle%FBodyGas = 0.0
              
            end if
            
            if (MaterialPointType==MaterialPointTypeLiquid) then
                
              NewParticle%FBody(:) = NewParticle%Porosity * NewParticle%ConstDensity * IntegrationWeight * CalParams%GravityData%GAccel * GravityVector(:)
              NewParticle%FBodyWater = 0.0
              NewParticle%FBodyGas = 0.0
              
            end if
            
          end if

          ! water weight is zero if using submerged weight 
          if(ISSubmerged)then
            if (MixedWeight <= MaterialWeight) then
              NewParticle%FBodyWater = 0.0
            end if
          end if

          ! Determine mixed (saturated/unsaturated) weight
          NewParticle%FBodyMixed(:) = MixedWeight * IntegrationWeight * GravityVector(:)


          if (.not.(NFORMULATION==1)) then
            if (MaterialPointType==MaterialPointTypeLiquid) then
              NewParticle%BulkWater = BulkWater
              NewParticle%Porosity = InitialPorosity
              NewParticle%InitialPorosity = InitialPorosity

            else
              NewParticle%BulkWater = BulkWater
              NewParticle%Porosity = InitialPorosity
              NewParticle%InitialPorosity = InitialPorosity
            end if
            
            if (MatParams(MaterialID)%MaterialType=='1-phase-liquid'.or.MatParams(MaterialID)%MaterialPhases=='1-phase-liquid') then
              NewParticle%PhaseStatus = PhaseStatusLIQUID
              NewParticle%ContainedMaterialTypes = ContainedMaterialTypeLIQUID
              NewParticle%FillingRatio = 1.0
            else if(MatParams(MaterialID)%MaterialType=='1-phase-solid'.or.MatParams(MaterialID)%MaterialPhases=='1-phase-solid') then
              NewParticle%PhaseStatus = PhaseStatusSOLID
              NewParticle%ContainedMaterialTypes = ContainedMaterialTypeSOLID
              NewParticle%FillingRatio = 1.0 - InitialPorosity
            else if(MatParams(MaterialID)%MaterialType=='2-phase'.or.MatParams(MaterialID)%MaterialPhases=='2-phase') then
              NewParticle%PhaseStatus = PhaseStatusSOLID
              NewParticle%ContainedMaterialTypes = ContainedMaterialTypeSOLIDLIQUID
              NewParticle%FillingRatio = 1.0
            end if
            NewParticle%ConcentrationRatioSolidL = 1.0 - InitialPorosity
            NewParticle%EffConcentrationRatioSolid = 1.0 - InitialPorosity
            NewParticle%ConcentrationRatioLiquidS = InitialPorosity
            NewParticle%ConcentrationRatioLiquidL = InitialPorosity
            NewParticle%EffConcentrationRatioLiquid = InitialPorosity
            NewParticle%PorosityL = InitialPorosity
            NewParticle%EffPorosity = InitialPorosity
         end if

          InitParticle = NewParticle

        end function InitParticle


        type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)) function CopyParticle(Particle, ID)
        !**********************************************************************
        !
        !   Function:  Copies all data from Particle to a new particle except
        !              for the ID. A new ID will be assigned to the new particle.
        !
        !   Particle : Particle to be copied
        !   IElTyp : Number of nodes per element
        !
        ! O CopyParticle : Returns a new particle containing the data of Particle
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) ::  ID
          
          call SetParticleStructureDefault(CopyParticle)
          
          if (ID>0) then
     
            CopyParticle%Kind = Particle%Kind
            CopyParticle%MassGas = Particle%MassGas
            CopyParticle%MassMixed = Particle%MassMixed
            CopyParticle%Damping = Particle%Damping
            CopyParticle%DBulkViscousPressure = Particle%DBulkViscousPressure
            CopyParticle%MaterialWeight = Particle%MaterialWeight
            CopyParticle%WaterWeight = Particle%WaterWeight
            CopyParticle%GasWeight = Particle%GasWeight
            CopyParticle%SatWeight = Particle%SatWeight
            CopyParticle%UnsatWeight = Particle%UnsatWeight
            CopyParticle%MixedWeight = Particle%MixedWeight
            CopyParticle%Conductivity = Particle%Conductivity
            CopyParticle%ConductivityGas = Particle%ConductivityGas
            CopyParticle%Porosity = Particle%Porosity
            CopyParticle%ESM_UnloadingStiffness = Particle%ESM_UnloadingStiffness
            CopyParticle%InitialPorosity = Particle%InitialPorosity
            CopyParticle%DegreeSaturation = Particle%DegreeSaturation
            CopyParticle%BulkWater = Particle%BulkWater
            CopyParticle%BulkGas = Particle%BulkGas
            CopyParticle%AirInWaterMassFraction = Particle%AirInWaterMassFraction
            CopyParticle%VapourInGasMassFraction = Particle%VapourInGasMassFraction
            CopyParticle%DiffusionAirInWater = Particle%DiffusionAirInWater
            CopyParticle%DiffusionVapourInGas = Particle%DiffusionVapourInGas
            CopyParticle%WaterAdvectiveFlux = Particle%WaterAdvectiveFlux ! Water Advective flux
            CopyParticle%AirAdvectiveFlux = Particle%AirAdvectiveFlux
            CopyParticle%Density = Particle%Density
            CopyParticle%ConstDensity = Particle%ConstDensity
            CopyParticle%ShearModulus = Particle%ShearModulus
            CopyParticle%CohesionCosPhi = Particle%CohesionCosPhi
            CopyParticle%CohesionStSoft = Particle%CohesionStSoft
            CopyParticle%PhiStSoft = Particle%PhiStSoft
            CopyParticle%PsiStSoft = Particle%PsiStSoft 
            CopyParticle%pp = Particle%pp
            CopyParticle%SFail = Particle%SFail
            CopyParticle%IPL = Particle%IPL
            CopyParticle%IntegrationWeight = Particle%IntegrationWeight
            CopyParticle%EoedHP = Particle%EoedHP
            CopyParticle%Alpha = Particle%Alpha
            CopyParticle%Beta = Particle%Beta

            CopyParticle%LiquidFreeSurface = Particle%LiquidFreeSurface
            CopyParticle%LiquidFreeSurfaceCumul = Particle%LiquidFreeSurfaceCumul
            CopyParticle%PhaseStatus = Particle%PhaseStatus
            CopyParticle%ContainedMaterialTypes = Particle%ContainedMaterialTypes
            CopyParticle%ConcentrationRatioSolidL = Particle%ConcentrationRatioSolidL
            CopyParticle%ConcentrationRatioLiquidL = Particle%ConcentrationRatioLiquidL
            CopyParticle%FillingRatio = Particle%FillingRatio

            CopyParticle%PorosityL = Particle%PorosityL

            CopyParticle%ConcentrationRatioLiquidS = Particle%ConcentrationRatioLiquidS
            CopyParticle%EffConcentrationRatioSolid = Particle%EffConcentrationRatioSolid
            CopyParticle%EffConcentrationRatioLiquid = Particle%EffConcentrationRatioLiquid
            CopyParticle%EffPorosity = Particle%EffPorosity

            call SetFExt(CopyParticle, Particle%FExt)
            call SetFExtWater(CopyParticle, Particle%FExtWater)
            call SetFExtGas(CopyParticle, Particle%FExtGas)

            CopyParticle%FBody = Particle%FBody
            CopyParticle%FBodyWater = Particle%FBodyWater
            CopyParticle%FBodyGas = Particle%FBodyGas
            CopyParticle%FBodyMixed = Particle%FBodyMixed

            call SetPrescrDisp(CopyParticle, Particle%PrescrDisp)
            CopyParticle%PrescrVelo = Particle%PrescrVelo !Prescribed particle velocity
            CopyParticle%MPPrescribedVeloDir = Particle%MPPrescribedVeloDir
            CopyParticle%LocPos = Particle%LocPos
            CopyParticle%IsBoundaryParticle = Particle%IsBoundaryParticle
            CopyParticle%IsHydrHeadBoundaryParticle = Particle%IsHydrHeadBoundaryParticle

            CopyParticle%AccelerationWater = Particle%AccelerationWater
            CopyParticle%UStepWater = Particle%UStepWater
            CopyParticle%UW = Particle%UW
            CopyParticle%UG = Particle%UG
            CopyParticle%WaterPressure = Particle%WaterPressure
            CopyParticle%WaterPressure0 = Particle%WaterPressure0
            CopyParticle%HydraulicHead = Particle%HydraulicHead
            
            CopyParticle%GasPressure = Particle%GasPressure
            CopyParticle%GasPressure0 = Particle%GasPressure0

            CopyParticle%WaterVolumetricStrain = Particle%WaterVolumetricStrain
            CopyParticle%GasVolumetricStrain= Particle%GasVolumetricStrain
            
            CopyParticle%Temperature = Particle%Temperature
            CopyParticle%ThermalCapacity = Particle%ThermalCapacity
            
            call SetEps(CopyParticle, Particle%Eps)

            call SetEpsP(CopyParticle, Particle%EpsP)

            call SetEpsStep(CopyParticle, Particle%EpsStep)
            call SetEpsStepPrevious(CopyParticle, Particle%EpsStepPrevious)
            call SetSigmaEffStep(CopyParticle, Particle%SigmaEffStep)
            call SetSigmaPrin(CopyParticle, Particle%SigmaPrin)
            call SetVariableData(CopyParticle, Particle%VariableData)
            call SetHPStateVariables(CopyParticle, Particle%HPStateVariables)
            call SetHPIGStateVariables(CopyParticle, Particle%HPIGStateVariables)
            call SetModifiedHPStateVariables(CopyParticle, Particle%ModifiedHPStateVariables)
          end if
        
        end function CopyParticle

        function GetFExt(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the external load vector assigned to Particle.
        !               If no load vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose external force information is returned
        !
        ! O GetFExt : External force vector assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS) :: GetFExt

          GetFExt = Particle%FExt

        end function GetFExt

        function GetFExtWater(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the external water load vector assigned to Particle.
        !               If no load vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose external force information is returned
        !
        ! O GetFExtWater : External force vector assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS) :: GetFExtWater
        
          GetFExtWater = Particle%FExtWater
        
        end function GetFExtWater

        function GetFExtGas(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the external water load vector assigned to Particle.
        !               If no load vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose external force information is returned
        !
        ! O GetFExtWater : External force vector assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS) :: GetFExtGas
        
          GetFExtGas = Particle%FExtGas
        
        end function GetFExtGas

        function GetPrescrDisp(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the prescribed displacement vector assigned to Particle.
        !               If no vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose prescribed displacement information is returned
        !
        ! O GetPrescrDisp : Prescribed displacement vector assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NVECTOR) :: GetPrescrDisp
        
          GetPrescrDisp = Particle%PrescrDisp
        
        end function GetPrescrDisp




        real(REAL_TYPE) function GetEpsI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the total strain (1 to 6) assigned to Particle.
        !               If no strain information is assigned to Particle, 0.0 is returned.
        !
        !   Particle : Particle whose strain information is returned
        !   I : Strain component identifier (1 to 6)
        !
        ! O GetEpsI : Strain component I assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I

          GetEpsI = Particle%Eps(I)

        end function GetEpsI

        function GetEps(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the total strain vector assigned to Particle.
        !               If no vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose total strain information is returned
        !
        ! O GetEps : Vector containing the total strain components assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NTENSOR) :: GetEps

          GetEps = Particle%Eps

        end function GetEps
        
        real(REAL_TYPE) function GetEpsPI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the total plastic strain (1 to 6) assigned to Particle.
        !               If no strain information is assigned to Particle, 0.0 is returned.
        !
        !   Particle : Particle whose plastic strain information is returned
        !   I : Plastic Strain component identifier (1 to 6)
        !
        ! O GetEpsI : Plastic Strain component I assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I

          GetEpsPI = Particle%EpsP(I)

        end function GetEpsPI

        function GetEpsP(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the total plastic strain (1 to 6) assigned to Particle.
        !               If no vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose plastic strain information is returned
        !
        ! O GetEpsP : Vector containing the plastic strain components assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NTENSOR) :: GetEpsP

          GetEpsP = Particle%EpsP

        end function GetEpsP

        real(REAL_TYPE) function GetEpsStepI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the incremental strain (1 to 6) assigned to Particle.
        !               If no strain information is assigned to Particle, 0.0 is returned.
        !
        !   Particle : Particle whose strain information is returned
        !   I : Strain component identifier (1 to 6)
        !
        ! O GetEpsStepI : Strain component I assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I

          GetEpsStepI = Particle%EpsStep(I)

        end function GetEpsStepI

        real(REAL_TYPE) function GetEpsStepPreviousI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the incremental strain (1 to 6) assigned to Particle.
        !               If no strain information is assigned to Particle, 0.0 is returned.
        !
        !   Particle : Particle whose strain information is returned
        !   I : Strain component identifier (1 to 6)
        !
        ! O GetEpsStepI : Strain component I assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I

          GetEpsStepPreviousI = Particle%EpsStepPrevious(I)

        end function GetEpsStepPreviousI

        function GetEpsStep(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the incremental strain vector assigned to Particle.
        !               If no vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose incremental strain information is returned
        !
        ! O GetEpsStep : Vector containing the incremental strain components assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NTENSOR) :: GetEpsStep

          GetEpsStep = Particle%EpsStep

        end function GetEpsStep

        function GetSigmaEffStep(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the incremental stress vector assigned to Particle.
        !               If no vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose stress information is returned
        !
        ! O GetSigmaEffStep : Stress components assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NTENSOR) :: GetSigmaEffStep

          GetSigmaEffStep = Particle%SigmaEffStep

        end function GetSigmaEffStep
        
        real(REAL_TYPE) function GetSigmaPrinI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the principal stress (1 to 6) of the current step assigned to Particle.
        !               If no stress information is assigned to Particle, 0.0 is returned.
        !
        !   Particle : Particle whose stress information is returned
        !   I : Stress component identifier (1 to 6)
        !
        ! O GetSigmaEff0I : Stress component I assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          
          GetSigmaPrinI = Particle%SigmaPrin(I)
        
        end function GetSigmaPrinI

        function GetSigmaPrin(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the principal stress vector assigned to Particle.
        !               If no vector is assigned to Particle, a vector
        !               containing only 0.0 values is returned.
        !
        !   Particle : Particle whose stress information is returned
        !
        ! O GetSigmaEffStep : Stress components assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(NTENSOR):: GetSigmaPrin
        
          GetSigmaPrin = Particle%SigmaPrin
        
        end function GetSigmaPrin

        real(REAL_TYPE) function GetVariableDataI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the variable data (1 to 3) assigned to Particle.
        !               If no arbitrary data is assigned to Particle, 0.0 is returned.
        !
        !   Particle : Particle whose arbitrary information is returned
        !   I : Arbitrary array component identifier (1 to 3)
        !
        ! O GetVariableDataI : Data component I assigned to Particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          
          GetVariableDataI = Particle%VariableData(I)
        
        end function GetVariableDataI

        subroutine SetFExt(Particle, FExt)
        !**********************************************************************
        !
        !    Function:  Sets FExt of Particle to the provided FExt.
        !
        !   Particle : Particle whose external forces are set
        !   FExt : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS), intent(in) :: FExt

          Particle%FExt = FExt

        end subroutine SetFExt

        subroutine SetFExtWater(Particle, FExtWater)
        !**********************************************************************
        !
        !    Function:  Sets FExtWater of Particle to the provided FExtWater.
        !
        !   Particle : Particle whose external forces are set
        !   FExtWater : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS), intent(in) :: FExtWater
          
          Particle%FExtwater = FExtWater
        
        end subroutine SetFExtWater

        subroutine SetFExtGas(Particle, FExtGas)
        !**********************************************************************
        !
        !    Function:  Sets FExtGas of Particle to the provided FExtGas.
        !
        !   Particle : Particle whose external forces are set
        !   FExtGas : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS), intent(in) :: FExtGas
          
          Particle%FExtGas = FExtGas
        
        end subroutine SetFExtGas

        subroutine IncreaseFExt(Particle, FExtIncrease, ILoadSystem)
        !**********************************************************************
        !
        !    Function:  Increases the external load applied on Particle.
        !
        !   Particle : Particle whose load is increased
        !   FExtIncrease : Increase of the particle load
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: FExtIncrease
          integer(INTEGER_TYPE), intent(in) :: ILoadSystem

          Particle%FExt(:,ILoadSystem) = Particle%FExt(:,ILoadSystem) + FExtIncrease

        end subroutine IncreaseFExt
        
         subroutine IncreaseFExtWater(Particle, FExtWaterIncrease,ILoadSystem)
        !**********************************************************************
        !
        !    Function:  Increases the external water load applied on Particle.
        !
        !   Particle : Particle whose load is increased
        !   FExtWaterIncrease : Increase of the particle water load
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: FExtWaterIncrease
          integer(INTEGER_TYPE), intent(in) :: ILoadSystem

          Particle%FExtWater(:,ILoadSystem) = Particle%FExtWater(:,ILoadSystem) + FExtWaterIncrease

        end subroutine IncreaseFExtWater

         subroutine IncreaseFExtGas(Particle, FExtGasIncrease,ILoadSystem)
        !**********************************************************************
        !
        !    Function:  Increases the external gas load applied on Particle.
        !
        !   Particle : Particle whose load is increased
        !   FExtGasIncrease : Increase of the particle water load
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: FExtGasIncrease
          integer(INTEGER_TYPE), intent(in) :: ILoadSystem

          Particle%FExtGas(:,ILoadSystem) = Particle%FExtGas(:,ILoadSystem) + FExtGasIncrease

        end subroutine IncreaseFExtGas

        subroutine SetPrescrDisp(Particle, PrescrDisp)
        !**********************************************************************
        !
        !    Function:  Sets PrescrDisp of Particle to the provided PrescrDisp.
        !
        !   Particle : Particle whose prescribed displacements are set
        !   PrescrDisp : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: PrescrDisp
          
          Particle%PrescrDisp = PrescrDisp
        
        end subroutine SetPrescrDisp

        subroutine SetPrescrDispI(Particle, I, Value)
        !**********************************************************************
        !
        !    Function:  Sets component I of PrescrDisp of Particle to the provided value.
        !
        !   Particle : Particle whose prescribed displacement I is set
        !   I : Component of the prescribed displacement vector (1 .. NVECTOR)
        !   Value : Value assigned to the component I of particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          real(REAL_TYPE), intent(in) :: Value
          
          if ( (I>=1) .and. (I<=NVECTOR) ) then
            Particle%PrescrDisp(I) = Value
          end if

        end subroutine SetPrescrDispI


        subroutine SetEps(Particle, Eps)
        !**********************************************************************
        !
        !    Function:  Sets Eps of Particle to the provided Eps.
        !
        !   Particle : Particle whose total strain data is set
        !   Eps : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: Eps

          Particle%Eps = Eps

        end subroutine SetEps

        subroutine SetEpsI(Particle, I, Value)
        !**********************************************************************
        !
        !    Function:  Sets component I of Eps of Particle to the provided value.
        !
        !   Particle : Particle whose total strain I is set
        !   I : Component of the total strain vector (1 .. 6)
        !   Value : Value assigned to the component I of particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          real(REAL_TYPE), intent(in) :: Value

          Particle%Eps(I) = Value

        end subroutine SetEpsI
 
         subroutine SetEpsP(Particle, EpsP)
        !**********************************************************************
        !
        !    Function:  Sets Eps of Particle to the provided Eps.
        !
        !   Particle : Particle whose total strain data is set
        !   Eps : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: EpsP

          Particle%EpsP = EpsP

        end subroutine SetEpsP

        subroutine SetEpsPI(Particle, I, Value)
        !**********************************************************************
        !
        !    Function:  Sets component I of EpsP of Particle to the provided value.
        !
        !   Particle : Particle whose total plastic strain I is set
        !   I : Component of the total plastic strain vector (1 .. 6)
        !   Value : Value assigned to the component I of particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          real(REAL_TYPE), intent(in) :: Value

          Particle%EpsP(I) = Value

        end subroutine SetEpsPI
        
 
        subroutine IncreaseEps(Particle, EpsStep)
        !**********************************************************************
        !
        !    Function:  Adds EpsStep to Eps of Particle.
        !
        !   Particle : Particle whose total strain data is increased
        !   EpsStep : Values added to Eps of the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: EpsStep

          Particle%Eps = Particle%Eps + EpsStep

        end subroutine IncreaseEps

        subroutine IncreaseEpsP(Particle, DEpsP)
        !**********************************************************************
        !
        !    Function:  Adds DEpsP to EpsP of Particle.
        !
        !   Particle : Particle whose total plastic strain data is increased
        !   EpsStep : Values added to EpsP of the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: DEpsP

          Particle%EpsP = Particle%EpsP + DEpsP

        end subroutine IncreaseEpsP

        subroutine SetEpsStep(Particle, EpsStep)
        !**********************************************************************
        !
        !    Function:  Sets EpsStep of Particle to the provided EpsStep.
        !
        !   Particle : Particle whose incremental strain data is set
        !   EpsStep : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: EpsStep

          Particle%EpsStep = EpsStep

        end subroutine SetEpsStep

        subroutine SetEpsStepPrevious(Particle, EpsStepPrevious)
        !**********************************************************************
        !
        !    Function:  Sets EpsStepPrevious of Particle to the provided EpsStepPrevious.
        !
        !   Particle : Particle whose incremental strain data is set
        !   EpsStepPrevious : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: EpsStepPrevious

          Particle%EpsStepPrevious = EpsStepPrevious

        end subroutine SetEpsStepPrevious

        subroutine SetEpsStepPreviousI(Particle, I, Value)
        !**********************************************************************
        !
        !    Function:  Sets component I of EpsStepPrevious of Particle to the provided value.
        !
        !   Particle : Particle whose incremental strain I is set
        !   I : Component of the incremental strain vector (1 .. 6)
        !   Value : Value assigned to the component I of particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          real(REAL_TYPE), intent(in) :: Value

          Particle%EpsStepPrevious(I) = Value

        end subroutine SetEpsStepPreviousI


        subroutine SetSigmaEffStep(Particle, SigmaEffStep)
        !**********************************************************************
        !
        !    Function:  Sets SigmaEffStep of Particle to the provided SigmaEffStep.
        !
        !   Particle : Particle whose incremental stress data is set
        !   SigmaEffStep : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: SigmaEffStep

          Particle%SigmaEffStep = SigmaEffStep

        end subroutine SetSigmaEffStep
  
        subroutine SetSigmaPrin(Particle, SigmaPrin)
        !**********************************************************************
        !
        !    Function:  Sets principal stresses of Particle to the provided VariableData.
        !
        !   Particle : Particle whose arbitrary data is set (for testing purposes)
        !   VariableData : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: SigmaPrin

          Particle%SigmaPrin = SigmaPrin

        end subroutine SetSigmaPrin
 
        subroutine SetVariableData(Particle, VariableData)
        !**********************************************************************
        !
        !    Function:  Sets VariableData of Particle to the provided VariableData.
        !
        !   Particle : Particle whose arbitrary data is set (for testing purposes)
        !   VariableData : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: VariableData
          
          Particle%VariableData = VariableData
        
        end subroutine SetVariableData

        subroutine SetVariableDataI(Particle, I, Value)
        !**********************************************************************
        !
        !    Function:  Sets component I of VariableData of Particle to the provided value.
        !
        !   Particle : Particle whose arbitrary data I is set
        !   I : Component of the arbitrary data vector (1 .. 3)
        !   Value : Value assigned to the component I of particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          real(REAL_TYPE), intent(in) :: Value
          
          if ( (I>=1) .and. (I<=NVECTOR) ) then
            Particle%VariableData(I) = Value
          end if
        
        end subroutine SetVariableDataI

        real(REAL_TYPE) function DeterminePartMass(MaterialWeight, GGrav, IntegrationWeight)
        !**********************************************************************
        !
        !    Function:  Determines the (constant) mass of a particle.
        !
        !   MaterialWeight : Gamma
        !   GGrav : Gravity acceleration
        !   IntegrationWeight : Global integration weight
        !
        ! O DeterminePartMass : Returns the mass of a particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), intent(in) :: MaterialWeight, GGrav, IntegrationWeight

          DeterminePartMass = MaterialWeight / GGrav * IntegrationWeight
     
        end function DeterminePartMass
        



        subroutine SetHPStateVariablesI(Particle, I, Value)
        !**********************************************************************
        !
        !    Function:  Sets variable I of HPStateVariables of Particle to the provided value.
        !
        !   Particle : Particle whose state variable I is set
        !   I : the state variable
        !   Value : Value assigned to the state variable
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          real(REAL_TYPE), intent(in) :: Value
          
          if ( (I>=1).and.(I<=2) ) then
            Particle%HPStateVariables(I) = Value
          end if
        
        end subroutine SetHPStateVariablesI

        subroutine SetModifiedHPStateVariablesI(Particle, I, Value)
        !**********************************************************************
        !
        !    Function:  Sets variable I of ModifiedHPStateVariables of Particle to the provided value.
        !
        !   Particle : Particle whose state variable I is set
        !   I : the state variable
        !   Value : Value assigned to the state variable
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          real(REAL_TYPE), intent(in) :: Value
          
          if ( (I>=1).and.(I<=2) ) then
            Particle%ModifiedHPStateVariables(I) = Value
          end if
        
        end subroutine SetModifiedHPStateVariablesI
        
        real(REAL_TYPE) function GetHPStateVariablesI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the state variable 
        !
        !   Particle : Particle whose state variable is returned
        !   I : the state variable
        !
        ! O GetHPStateVariablesI : value assigned to the state variable
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          
          GetHPStateVariablesI = Particle%HPStateVariables(I)
        
        end function GetHPStateVariablesI
        
        real(REAL_TYPE) function GetModifiedHPStateVariablesI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the state variable 
        !
        !   Particle : Particle whose state variable is returned
        !   I : the state variable
        !
        ! O GetModifiedHPStateVariablesI : value assigned to the state variable
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          
          GetModifiedHPStateVariablesI = Particle%ModifiedHPStateVariables(I)
        
        end function GetModifiedHPStateVariablesI

        subroutine SetHPIGStateVariablesI(Particle, I, Value)
        !**********************************************************************
        !
        !    Function:  Sets variable I of HPIGStateVariables of Particle to the provided value.
        !
        !   Particle : Particle whose state variable I is set
        !   I : the state variable
        !   Value : Value assigned to the state variable
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          real(REAL_TYPE), intent(in) :: Value
          
          if ( (I>=1).and.(I<=7) ) then
            Particle%HPIGStateVariables(I) = Value
          end if

        end subroutine SetHPIGStateVariablesI

        real(REAL_TYPE) function GetHPIGStateVariablesI(Particle, I)
        !**********************************************************************
        !
        !    Function:  Returns the state variable 
        !
        !   Particle : Particle whose state variable is returned
        !   I : the state variable
        !
        ! O GetHPIGStateVariablesI : value assigned to the state variable
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          integer(INTEGER_TYPE), intent(in) :: I
          
          GetHPIGStateVariablesI = Particle%HPIGStateVariables(I)
        
        end function GetHPIGStateVariablesI
        

        
        function GetHPStateVariables(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the state variable array of the HP basic model
        !
        !   Particle : Particle whose state variable array is returned
        !
        ! 
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(2) :: GetHPStateVariables
        
          GetHPStateVariables = Particle%HPStateVariables
        
        end function GetHPStateVariables
        
        function GetModifiedHPStateVariables(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the state variable array of the HP basic model
        !
        !   Particle : Particle whose state variable array is returned
        !
        ! 
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(2) :: GetModifiedHPStateVariables
        
          GetModifiedHPStateVariables = Particle%ModifiedHPStateVariables
        
        end function GetModifiedHPStateVariables
        
        function GetHPIGStateVariables(Particle)
        !**********************************************************************
        !
        !    Function:  Returns the state variable array of the HP-IG basic model
        !
        !   Particle : Particle whose state variable array is returned
        !
        ! 
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(in) :: Particle
          real(REAL_TYPE), dimension(7) :: GetHPIGStateVariables
        
          GetHPIGStateVariables = Particle%HPIGStateVariables
        
        end function GetHPIGStateVariables
        
        
        subroutine SetHPStateVariables(Particle, HP)
        !**********************************************************************
        !
        !    Function:  Set the state variables of the HP basic model 
        !
        !   Particle : Particle whose external forces are set
        !   ICH : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(2), intent(in) :: HP
          
          Particle%HPStateVariables = HP
        
        end subroutine SetHPStateVariables

        subroutine SetModifiedHPStateVariables(Particle, HP)
        !**********************************************************************
        !
        !    Function:  Set the state variables of the HP basic model 
        !
        !   Particle : Particle whose external forces are set
        !   ICH : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(2), intent(in) :: HP
          ! Local variables

          Particle%ModifiedHPStateVariables = HP

        end subroutine SetModifiedHPStateVariables

        subroutine SetHPIGStateVariables(Particle, HPIG)
        !**********************************************************************
        !
        !    Function:  Set the state variables of the HP IG basic model 
        !
        !   Particle : Particle whose external forces are set
        !   ICH : Values assigned to the particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
          real(REAL_TYPE), dimension(7), intent(in) :: HPIG
          ! Local variables
          
          Particle%HPIGStateVariables = HPIG
        
        end subroutine SetHPIGStateVariables

      end module ModParticle