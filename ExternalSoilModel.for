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
    !	Copyright (C) 2021  Members of the Anura3D MPM Research Community 
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
    integer(INTEGER_TYPE), parameter :: nAddVar = 11
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
    real(REAL_TYPE) :: Kf
    real(REAL_TYPE) :: N ! Porosity         
    real(REAL_TYPE) :: DEpsVol ! Incremental volumetric strain (water)
    real(REAL_TYPE) :: DEpsVolW ! Incremental volumetric strain (water)
    real(REAL_TYPE) :: dT ! Change of Temperature at integration point
    real(REAL_TYPE) :: lambda !Used to compute water pressure
    real(REAL_TYPE) :: Sr, dSrdp(1) !Degree of saturation and derivative dSr/dp_w
    pointer (p, ESM)             
          
    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
    NameModel = MatParams(IDset)%MaterialModel ! name of constitutive model as specified in GOM-file
    ntens = MatParams(IDset)%UMATDimension     
          
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
    !if (CalParams%ApplyEffectiveStressAnalysis) then
    if (IsUndrEffectiveStress) then
        if (Particles(IDpt)%Porosity > 0.0) then
        Bulk = Particles(IDpt)%BulkWater / Particles(IDpt)%Porosity ! kN/m2
        DSigWP = Bulk * DEpsVol
        else
        DSigWP = 0.0
        end if
    end if ! effective stress analysis
          
    ! for 2-phase analysis
    if ((CalParams%NumberOfPhases==2).and.(NFORMULATION==1)) then
    if (Particles(IDpt)%WaterWeight > 0.0) then
      DEpsVolW = Particles(IDpt)%WaterVolumetricStrain ! Water phase
      N = Particles(IDpt)%Porosity
      Kf = Particles(IDpt)%BulkWater
        
      If (CalParams%ApplyPartialSaturation) then
        Sr = Particles(IDPt)%DegreeSaturation
        call CalculateDerivDegreeSaturation(IDPt,dSrdp,1)
        
        If ((Sr<1).and.(Sr>0)) then
          lambda = N/Kf - N/Sr * dSrdp(1)
          lambda = 1/lambda
        else
          lambda = Kf/N   
        end if
      else
        lambda = Kf/N
      end if
      DSigWP = lambda*( N * DEpsVolW + (1.0 - N) * DEpsVol)
        
        ! for submerged calculation
        if (CalParams%ApplySubmergedCalculation) then 
        if (CalParams%IStep <= CalParams%NumberSubmergedCalculation) then 
            DSigWP = 0.0 ! excess pore pressure is zero in gravity phase
        end if
        end if ! submerged calculation
              
    else 
        DSigWP = 0.0
    end if
    end if ! 2-phase analysis

    ! for 3-phase analysis (unsaturated soil)
    if (CalParams%NumberOfPhases==3) then ! solve mass balances and energy balance
    ! solving balance equations for the unsaturated soil
    ! the incremental water pressure, gas pressure and temperature (for one material point) are obtained
    call SolveBalanceEquations(IDpt, DEpsVol, DSigWP, DSigGP, dT)
    end if ! 3-phase analysis
          
    ! single phase quasi-static consolidation
    if (CalParams%ApplyImplicitQuasiStatic) then
        DsigWP = Particles(IDpt)%WaterPressure - Particles(IDpt)%WaterPressure0
    end if
          
    ! get stresses in integration/material point      
    do I = 1, NTENSOR
    Sig0(I) = SigmaEff0Array(IDpt, I) ! get initial stress of current step assigned to point 
    end do
    Stress=0.0
    do I=1, NTENSOR
        Stress(I) = Stress(I) + Sig0(I) !To use UMAT for 3D also for 2D
    enddo 
          
    ! initialise state variables (only for very first time and load step)
    if ((CalParams%IStep == 1).and.(CalParams%TimeStep == 1)) then
    StateVar = MatParams(IDset)%ESM_Statvar_in
    else 
    StateVar = ESMstatevArray(IDpt,:)
    end if 
   
    
    call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP) !Note that the subroutine checks Cavitation Threshold & Gas Pressure
          
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
          
    ! get name of DLL
    cmname = MatParams(IDSet)%SoilModelDLL
    ! get material properties  
    props = MatParams(IDSet)%ESM_Solid
         
    if (trim(NameModel)//char(0) == trim('linear_elasticity')//char(0)) then
    props(1) = Particles(IDpt)%ShearModulus ! shear modulus, G
    cmname = UMAT_LINEAR_ELASTICITY
    elseif (trim(NameModel)//char(0) == trim(ESM_MOHR_COULOMB_STANDARD)//char(0)) then
    props(1) = Particles(IDpt)%ShearModulus ! shear modulus, G
    props(2) = MatParams(IDSet)%PoissonRatio ! shear modulus, G
    props(3) = SIN(MatParams(IDSet)%FrictionAngle*(Pi/180.0)) ! shear modulus, G
    props(4) = Particles(IDpt)%CohesionCosPhi ! shear modulus, G
    props(5) = SIN(MatParams(IDSet)%DilatancyAngle*(Pi/180.0)) ! shear modulus, G
    props(6) = MatParams(IDSet)%TensileStrength ! shear modulus, G
    cmname = UMAT_MOHR_COULOMB_STANDARD
    endif          
    ! initialise UMAT
    p = GetProcAddress(MatParams(IDSet)%SoilModelDLLHandle, "ESM"C) ! TODO: Check to see if multiple materials works ?
    call ESM(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
    ! save unloading stiffness in Particles array  
    Particles(IDpt)%ESM_UnloadingStiffness = Eunloading
                 
    call SetIPL(IDpt, IDel, int(PlasticMultiplier))


    
    ! to using objective stress definition
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

end module ModExternalSoilModel
