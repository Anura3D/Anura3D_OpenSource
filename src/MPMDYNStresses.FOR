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
    
    
	  module ModMPMDYNStresses
      !**********************************************************************
      !
      !  Function: Contains the routines related to calculating the stresses at material points
      !            for partially filled elements and Gauss points for fully filled elements.
      !            This module is only used by the dynamic MPM.  
      !
      !            Note: 'Integration point' denotes either a Gauss point or material point. Except for
      !            the data that is used or modified, both types of points are treated identically.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      use ModGlobalConstants
      use ModCounters
      use ModReadCalculationData
      use ModReadMaterialData
      use ModElementEvaluation
      use ModMPMData
      use ModMPMInit
      use ModWriteTestData
      use ModMPMStresses
      use ModStrainSmoothing
      use ModMeshInfo
      use ModLiquid
      use ModTwoLayerFormulation
      use ModExternalSoilModel
      use user32
      use kernel32
      use ModMPMDYN2PhaseSP
       
      contains	  
	  
        subroutine MPMDYNGetSig()
        !**********************************************************************
        !
        !  Function:  Loops over the active elements and updates the stresses and
        !             plasticity state of Gauss points for fully filled elements and
        !             material points (particles) for partially filled elements by calling the appropriate
        !             constitutive model routine of the material assigned to a Gauss point
        !             or material point.
        !             Assign Unloading Elastic Modulus to material points (particles) for fully filled elements
        !
        !             Structure of      [ D1  D2  D2  o   o   o ]
        !             elastic D matrix  [ D2  D1  D2  o   o   o ]
        !                               [ D2  D2  D1  o   o   o ]
        !                               [ o   o   o  GG   o   o ]
        !                               [ o   o   o   o  GG   o ]
        !                               [ o   o   o   o   o  GG ]
        !
        !**********************************************************************

        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IntGlo,       & ! Global ID of Gauss point or particle 
                     IEl,          & ! Element ID &
                     IAEl,         & ! Active element ID &
                     NElemPart,    & ! Number of Gauss points or particles (material points) per element &
                     Int,          & ! Local integration point counter &
                     IMatSet,      & ! Counter on number of material sets &
                     IDof,         & ! Counter &
                     NMaterialSets,& ! Number of material sets &
                     IEntity,      &
                     I,            & ! Counter
                     IDMaterial,   &
                     NumParticles, &      ! ID of Material Set, number of particles (material points) inside the element
                     FirstParticleIndex
          real(REAL_TYPE) :: WTN,  FirstParticleEUnloading

          real(REAL_TYPE), dimension(Counters%N) :: DDisp, DDispLiquid, VelLiquid
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: BMatrix
          logical :: DoSkipStressComputation

          ! Initialise global variables
          CalParams%IntegrationPointData%NPlasticPoints = 0
          CalParams%IntegrationPointData%NNegativePlasticPoints = 0
          CalParams%IntegrationPointData%NApexPoints = 0
          CalParams%IntegrationPointData%NTensionCutOffPoints = 0
          CalParams%ConvergenceCheck%NInaccuratePlasticPoints = 0
          CalParams%ConvergenceCheck%SumLocalError = 0.0
          CalParams%ConvergenceCheck%SumIntegrationPointWeights = 0.0

          if (CalParams%ApplyBulkViscosityDamping) then
            RateVolStrain = 0.0
          end if
          
          if (CalParams%ApplyStrainSmoothing) then   !if strain smoothing is used
            ! initialise smoothened strain or strain rate data 
            NMaterialSets = Counters%NLayers
            do IMatSet = 1, NMaterialSets

              IEntity =getMaterialEntity(IMatSet)

              if (NFORMULATION==1) then                                ! one set of material points (1-point formulation)
                if (.not.CalParams%ApplyContactAlgorithm) then         ! no contact
                  do IDof  = 1, Counters%N                             ! loop all dof's
                    DDisp (IDof) = IncrementalDisplacementSoil(IDof,1) ! if no contact is used then the system always has one shared displacement
                  end do
                else if (CalParams%ApplyContactAlgorithm) then                ! contact is used
                  do IDof  = 1, Counters%N  
                    DDisp (IDof) = IncrementalDisplacementSoil(IDof,IEntity)  ! if contact is used then each entity has its own vector 
                  end do
                end if
              end if

              if (.not.(NFORMULATION==1)) then
                do IDof  = 1, Counters%N
                  DDisp (IDof) = IncrementalDisplacementSoil(IDof,1)
                end do
                call UpdateParticleStrains() !update strains (assign strains to all solid MPs)
              end if

              ElementStrain = 0.0
              if (CalParams%ApplyStrainSmoothing)  then
                call SmoothenStrains(DDisp, IMatSet)
              end if

            end do   ! NMaterialSets
          end if     ! strain smoothing
              
          if (CalParams%ApplyStrainSmoothingLiquidTwoLayer) then
            ! initialise smoothened strain or strain rate data 
            NMaterialSets = Counters%NLayers
            do IMatSet = 1, NMaterialSets

              if (.not.(NFORMULATION==1)) then
                do IDof  = 1, Counters%N
                  DDispLiquid (IDof) = IncrementalDisplacementWater(IDof,1)
                  VelLiquid (IDof)   = IncrementalDisplacementWater(IDof,1)/CalParams%TimeIncrement
                end do
                call UpdateParticleStrainsLiquidTwoLayer() !update strains (assign strains to all liquid MPs)
              end if

              ElementStrain = 0.0

              if ((.not.(CalParams%NumberOfPhases==1))) then
                call ComputeTwoLayerVolumetricStrainLiquid(IncrementalDisplacementWater, IncrementalDisplacementSoil)
                call SmoothenStrainsLiquid(DDispLiquid, IMatSet) ! strain smoothing for water
              else
                call SmoothenStrainsLiquid(DDispLiquid, IMatSet) ! strain smoothing for water
              end if
            end do   ! NMaterialSets
          end if     ! strain smoothing Liquid 2 layer
              
          if ((NFORMULATION==1).and.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and. &
                    CalParams%ApplyStrainSmoothing) then  
            ! initialise smoothened strain data 
            NMaterialSets = Counters%NLayers 
            do IMatSet = 1, NMaterialSets
              ElementStrain = 0.0
              call SmoothenStrainsWater(IMatSet)
            end do   ! NMaterialSets

            if (CalParams%NumberOfPhases==3) then
              do IMatSet = 1, NMaterialSets
                ElementStrain = 0.0
                call SmoothenStrainsGas(IMatSet)
              end do  ! NMaterialSets
            end if
          end if
          
          !!******STRAT calculation of fluid pressures (liquid and gas pressure)         
          if(NFORMULATION==1) then       ! one set of material points (1-point formulation)
          do IAEl = 1, Counters%NAEl     ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)            
            NElemPart = NumberOfIntegrationPoints(IEl)    ! Get the number of integration per element (Gaussian integration points or number of material points)
           
            ! Calculate stresses in integration/material points
            do Int = 1, NElemPart ! Loop over all integration/material points of the element
        
              IntGlo = GetParticleIndex(Int, IEl)   ! Determine global ID of integration point 
              IDMaterial = MaterialIDArray(IntGlo)  ! Material number stored in $$MATERIAL_INDEX in the GOM-file
              
                if (CalParams%ApplyContactAlgorithm) then
                  IEntity = EntityIDArray(IntGlo) 
                else
                  IEntity = 1 
                end if

                DoSkipStressComputation = EntityIDArray(IntGlo) == HARD_ENTITY
                do I = 1, NVECTOR
                  DoSkipStressComputation = DoSkipStressComputation .and. CalParams%ApplyPrescribedVelocity(I)
                end do

                if (.not.DoSkipStressComputation) then
                    call GetFluidPressure(IntGlo)     !Provide liquid and gas pressures for single-point formulation, 2phase and 3phase analysis
                else
                    Particles(IntGlo)%WaterPressure = 0.0
                    Particles(IntGlo)%GasPressure = 0.0
                end if
                        
            end do   ! Loop over all integration/material points of the element
          end do     ! Loop over all active elements for computation of stresses

          !Apply pressure increment smoothing
          if (CalParams%ApplySmootheningLiquidPressureIncrement) then
            NMaterialSets = Counters%NLayers 
            do IMatSet = 1, NMaterialSets
              call SmoothenLiquidPressureIncrement(IMatSet)
            end do ! NMaterialSets
          end if
         
          end if
  
          !!******STRAT calculation of effective stresses
          do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)
            
            if (CalParams%ApplyObjectiveStress.or.CalParams%ApplyImplicitQuasiStatic) then 
              call GetBMatrix(IEl, 1, NodalCoordinatesUpd, BMatrix, WtN)
            end if
            
            ! Get the number of integration per element (Gaussian integration points or number of material points)
            NElemPart = NumberOfIntegrationPoints(IEl) 
     
            ! Calculate stresses in integration/material points
            do Int = 1, NElemPart ! Loop over all integration/material points of the element
         
              IntGlo = GetParticleIndex(Int, IEl)   ! Determine global ID of integration point 
              IDMaterial = MaterialIDArray(IntGlo)  ! Material number stored in $$MATERIAL_INDEX in the GOM-file
                          
              if ( (MaterialPointTypeArray(IntGlo)==MaterialPointTypeMixture).or.(MaterialPointTypeArray(IntGlo)==MaterialPointTypeSolid) ) then
                if (CalParams%ApplyContactAlgorithm) then
                  IEntity = EntityIDArray(IntGlo) 
                else
                  IEntity = 1 
                end if
                
                if((NFORMULATION==1).or.(Particles(IntGlo)%PhaseStatus==PhaseStatusSOLID)) then 

                 DoSkipStressComputation = EntityIDArray(IntGlo) == HARD_ENTITY
                 do I = 1, NVECTOR
                   DoSkipStressComputation = DoSkipStressComputation .and. CalParams%ApplyPrescribedVelocity(I)
                 end do
    
		  	     if (MatParams(IDMaterial)%MaterialModel==ESM_RIGID_BODY) then
					DoSkipStressComputation = .true.
                 end if
                  
                 if (.not.DoSkipStressComputation) then  
                   call StressSolid(IntGlo, IEl, BMatrix, IEntity)   ! calculate stresses for mixture or solid material points
                 else
                   do I = 1, NVECTOR     ! only first three indices, both for 2D and 3D
                     SigmaEffArray(IntGlo,I) = 0.D0
                   end do
                   Particles(IntGlo)%WaterPressure = 0.0
                 end if
                
                else if(Particles(IntGlo)%PhaseStatus==PhaseStatusLIQUID) then
                  SigmaEffArray(IntGlo,:) = 0.D0
                end if 
                  
               else if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then
                  if (MatParams(IDMaterial)%MaterialModel==ESM_BINGHAM_LIQUID) then
                    call StressSolid(IntGlo, IEl, BMatrix, IEntity) 
                  else
                    call GetPressureLiquidMaterialPoint(IntGlo, IEl)  ! provide liquid pressure for material points
                  end if
               end if
               
            end do ! Loop over all integration/material points of the element
                        
            ! Assign EUnloading to all particles (material points) of each element
                      
            if (NElemPart==1) then ! only for fully filled elements               
             NumParticles = NPartEle(IEl)
             FirstParticleIndex = GetParticleIndex(1, IEl)
             FirstParticleEUnloading = Particles(FirstParticleIndex)%ESM_UnloadingStiffness
             do Int = 2, NumParticles
                IntGlo = GetParticleIndex(Int, IEl)
                Particles(IntGlo)%ESM_UnloadingStiffness = FirstParticleEUnloading
             enddo
            endif
            
          end do ! Loop over all active elements for computation of stresses

          
          do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)
            
            if (CalParams%ApplyObjectiveStress.or.CalParams%ApplyImplicitQuasiStatic) then 
              call GetBMatrix(IEl, 1, NodalCoordinatesUpd, BMatrix, WtN)
            end if
            
            ! Get the number of integration per element (Gaussian integration points or number of material points)
            NElemPart = NumberOfIntegrationPoints(IEl) 
          
            ! Calculate stresses in integration/material points
            do Int = 1, NElemPart ! Loop over all integration/material points of the element
         
             IntGlo = GetParticleIndex(Int, IEl)   ! Determine global ID of integration point
             IDMaterial = MaterialIDArray(IntGlo)  ! Material number stored in $$MATERIAL_INDEX in the GOM-file
              
             if (CalParams%ApplyContactAlgorithm) then
               IEntity = EntityIDArray(IntGlo)    ! Get the entity ID for the current particle          
             else
               IEntity = 1 
             end if
             if(.not.(CalParams%ApplyBinghamFluid.or.CalParams%ApplyFrictionalFluid & !v2016
              .or. (MatParams(IDMaterial)%MaterialModel==ESM_FRICTIONAL_LIQUID) &     !v2017.1 and following
              .or. (MatParams(IDMaterial)%MaterialModel==ESM_BINGHAM_LIQUID))) then
              !for Bingham and frictional fluid stresses are calculated elsewhere
               if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then
                 call GetStressesLiquidMaterialPoint(IntGlo, IEl) ! calculate stress for liquid material points
               end if
              end if 
             end do
          end do  ! Loop over all active elements for computation of stresses
               
        end subroutine MPMDYNGetSig

           
        subroutine GetStressesSoilMaterialPoint(IntGlo, IEl, BMatrix, IEntityID, WtN)
        !SUBROUTINE NOT USED. IT SHOULD BE REMOVED
        !*************************************************************************************
        !  Function:   Calculate Stress from constitutive model (lineal elastic)
        !              Calculate liquid and gas pressure increment from mass balance equations
        !
        !
        !*************************************************************************************        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IntGlo ! global ID of integration/material point
          integer(INTEGER_TYPE), intent(in) :: IEl ! current element ID

          ! B-matrix at the considered integration point (here only used if ApplyObjectiveStress=TRUE)
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: BMatrix
          integer(INTEGER_TYPE), intent(in) :: IEntityID ! entity ID (here only used if ApplyObjectiveStress=TRUE)

          ! global weight at considered integration point (here only used if ApplyImplicitQuasiStatic=TRUE)
          real(REAL_TYPE), intent(in) :: WtN
          
          ! local variables
          integer(INTEGER_TYPE) :: I ! counter
          integer(INTEGER_TYPE) :: ISet ! Material set ID of considered integration point
          integer(INTEGER_TYPE) :: ITens ! 1 in case of tension cut-off integration point
          integer(INTEGER_TYPE) :: IApex ! 1 in case of ... integration point
          real(REAL_TYPE) :: XNu ! Poisson's ratio of integration point
          real(REAL_TYPE) :: BFac ! Undrained material parameter of integration point
          real(REAL_TYPE) :: DEpsVol ! Incremental volumetric strain (used for undrained soil material)
          real(REAL_TYPE) :: DEpsVolW ! Incremental volumetric strain (water)
          real(REAL_TYPE) :: D2DDEpsV
          real(REAL_TYPE) :: Bulk ! Bulk modulus (?) for undrained behaviour if not user-defined models (?)
          real(REAL_TYPE) :: DSigWP ! Change of water pressure at integration point 
          real(REAL_TYPE) :: DSigGP ! Change of gas pressure at integration point 
          real(REAL_TYPE) :: Suction0 ! Suction
          real(REAL_TYPE) :: DSuction ! Suction increment
          real(REAL_TYPE) :: dT ! Change of Temperature at integration point 
          real(REAL_TYPE) :: SPhi ! Sinus Phi at integration point
          real(REAL_TYPE) :: SPsi ! Sinus Psi at integration point
          real(REAL_TYPE) :: GG ! Shear modulus at integration point
          real(REAL_TYPE) :: Cohec ! Cohesion at integration point
          real(REAL_TYPE) :: Tens ! Tension cut-off value at integration point
          real(REAL_TYPE) :: Fac, D1, D2 ! Local helper variables
          real(REAL_TYPE) :: TauMax ! ...
          !real(REAL_TYPE) :: cp, cr, phip, phir, psip, psir, factor, c, phi, psi !Softening parameters Strain Softening MC model
          !real(REAL_TYPE) :: lambda, kappa, M, init_void_ratio, pp, nu_u ! ModifiedCamClay: parameters
          real(REAL_TYPE) :: N ! Porosity
          real(REAL_TYPE) :: Kf, KSuc ! bulk modulus fluid, bulk modulus gas, elastic swelling index due to suction
          real(REAL_TYPE) :: ToleratedLocalError ! variables for ImplicitQuasiStatic calculation
          real(REAL_TYPE) :: LocalError ! variables for ImplicitQuasiStatic calculation
          real(REAL_TYPE), dimension(NTENSOR) :: Sig0 ! Local copy of initial stresses at integration point
          real(REAL_TYPE), dimension(NTENSOR) :: SigC ! Local copy of current stresses at integration point
          real(REAL_TYPE), dimension(NTENSOR) :: SigPrin ! Local copy of current principal stresses at integration point
          real(REAL_TYPE), dimension(NTENSOR) :: SigEQ
          real(REAL_TYPE), dimension(NTENSOR) :: DDSigE
          real(REAL_TYPE), dimension(NTENSOR, NTENSOR) :: D ! Local copy of material stiffness matrix and its inverse at integration point
          real(REAL_TYPE), dimension(NTENSOR, NTENSOR) :: DI ! Local copy of material stiffness matrix and its inverse at integration point

          ! Determined incremental strain at integration point  (Exx, Eyy, Ezz, Gxy, Gyz, Gzx)
          real(REAL_TYPE), dimension(NTENSOR) :: DEps

          ! Determined incremental elastic strain due to an increment of suction at integration point (Exx, Eyy, Ezz, 0, 0, 0)
          real(REAL_TYPE), dimension(NTENSOR) :: DEpsES
          real(REAL_TYPE), dimension(NTENSOR) :: DSig ! Determined incremental stress change at integration point
          real(REAL_TYPE), dimension(NTENSOR) :: DSigE ! Determined stress increment due to elastic strain increment
          real(REAL_TYPE), dimension(NTENSOR) :: SigE ! Total stress due to elastic strain SigE = Sig0 + DSigE
          real(REAL_TYPE), dimension(NTENSOR) :: DSigP ! Plastic stress increment at integration point DSigP = SigE - SigC
          real(REAL_TYPE), dimension(NTENSOR) :: DEpsP ! Plastic strain increment at integration point DEpsP = DI * DSigP
          real(REAL_TYPE), dimension(NTENSOR) :: DDEps
          ! strain increment of previous step (needed for implicit quasi-static integration)
          real(REAL_TYPE), dimension(NTENSOR) :: DEpsPrevious
          real(REAL_TYPE), dimension(NTENSOR) :: ErrorSig
          !real(REAL_TYPE), dimension(:), allocatable :: StVar ! state variables
          !real(REAL_TYPE), dimension(:), allocatable :: Props ! material properties
          character(len=64) :: SoilModel ! name of the constitutive model
          logical :: IsUndrEffectiveStress

     
          ! determine material data of integration/material point
          call GetMaterialData(IntGlo, ISet, XNu, BFac, SPhi, SPsi, GG, &
                               Cohec, Tens)  
          
          SoilModel = MatParams(ISet)%MaterialModel ! name of constitutive model as specified in GOM-file
          
          IsUndrEffectiveStress = &
              !code version 2016 and previous
              ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(ISet)%MaterialType)=='2-phase')) .or. &
              !code version 2017.1 and following
              (trim(MatParams(ISet)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))

          ! get stresses at last iteration step (Sig0=SigmaEff0), current stresses (SigC=SigmaEff), 
          ! and previous strain increment (DEpsPrevious=EpsStepPrevious)
          call AssignStressStrainToLocalArray(IntGlo, NTENSOR, Sig0, SigC, DEpsPrevious)

          ! initialise stress increment (DSig = difference between initial stress and current stress of load step
          DSig = 0.0
          
          ! get current strain increment (DEps=EpsStep)
          DEps = GetEpsStep(Particles(IntGlo))
          DEpsVol = DEps(1) + DEps(2) + DEps(3) ! volumetric strain, valid for 2D and 3D
                  
          ! initalise water pressure (only needed for undrained analyses)
          DSigWP = 0.0

          ! for effective stress analysis
          !if (CalParams%ApplyEffectiveStressAnalysis) then
          if (IsUndrEffectiveStress) then
              if (Particles(IntGlo)%Porosity > 0.0) then
                Bulk = Particles(IntGlo)%BulkWater / Particles(IntGlo)%Porosity ! kN/m2
                DSigWP = Bulk * DEpsVol
              else
                DSigWP = 0.0
              end if
          end if ! effective stress analysis
          
          ! for 2-phase analysis
          if ((CalParams%NumberOfPhases==2).and.(NFORMULATION==1)) then
            if (Particles(IntGlo)%WaterWeight > 0.0) then
              DEpsVolW = Particles(IntGlo)%WaterVolumetricStrain ! Water phase
              N = Particles(IntGlo)%Porosity
              Kf = Particles(IntGlo)%BulkWater
              DSigWP = Kf * DEpsVolW + ((1.0 - N) / N) * Kf * DEpsVol
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
            DSigWP = 0.0d0
            DSigGP = 0.0d0
            ! solving balance equations for the unsaturated soil
            ! the incremental water pressure, gas pressure and temperature (for one material point) are obtained
            call SolveBalanceEquations(IntGlo, DEpsVol, DSigWP, DSigGP, dT)
          end if ! 3-phase analysis
          
          ! single phase quasi-static consolidation
          if (CalParams%ApplyImplicitQuasiStatic) then
             DsigWP = Particles(IntGlo)%WaterPressure - Particles(IntGlo)%WaterPressure0
          end if
          
          call BuildDElastic(GG, XNu, D)  ! returns elastic stiffness matrix, D

          call BuildDElasticInverse(GG, XNu, DI) ! returns inverse of elastic stiffness matrix, DI
          
          ! calculate elastic stress increment (DSigE = elastic stiffness D * strain increment DEps)
          Fac = 2.0 * GG / (1.0 - 2.0 * XNu)
          D1 = FAC * (1.0 - XNU)
          D2 = FAC * XNU
          DSigE(1) = (D1 - D2) * DEps(1) + D2 * DEpsVol ! for 2D and 3D
          DSigE(2) = (D1 - D2) * DEps(2) + D2 * DEpsVol ! for 2D and 3D
          DSigE(3) = (D1 - D2) * DEps(3) + D2 * DEpsVol ! for 2D and 3D
          do I = 4, NTENSOR
            DSigE(I) = GG * DEps(I)
          end do
          
          ! for 3-phase analysis (unsaturated soil)
          ! DSigE = De * (DEps - DEpsES), where DEpsES is the elastic strain due to an increment of suction
          if (CalParams%NumberOfPhases==3) then
            Suction0 = Particles(IntGlo)%GasPressure - Particles(IntGlo)%WaterPressure  ! if is unsat --> Sution0 < 0
            DSuction = DSigGP - DSigWP 
            Suction0 = min(Suction0, 0.0d0)
            if (Suction0==0.0d0) then
              DSuction = 0.0d0
            end if
            KSuc = MatParams(ISet)%SwellingIndexSuc
            if (Suction0<0.0d0) then
              DEpsES = (KSuc/abs(Suction0+100))*DSuction
              DSigE(1) = DSigE(1) - (D(1,1) * DEpsES(1) + D(1,2) * DEpsES(2) + D(1,3) * DEpsES(3))
              DSigE(2) = DSigE(2) - (D(2,1) * DEpsES(1) + D(2,2) * DEpsES(2) + D(2,3) * DEpsES(3))
              DSigE(3) = DSigE(3) - (D(3,1) * DEpsES(1) + D(3,2) * DEpsES(2) + D(3,3) * DEpsES(3))
            end if
          end if ! 3-phase analysis
     
          ! elastic stress SigE = initial stress Sig0 + elastic stress increment DSigE
          SigE = Sig0 + DSigE
     
           
          ! for ImplicitQuasiStatic integratoin
          if (CalParams%ApplyImplicitQuasiStatic) then
            if (CalParams%ImplicitIntegration%Iteration > 1) then
              DDEps = DEps - DEpsPrevious
              D2DDEpsV = D2 * (DDEps(1) + DDEps(2) + DDEps(3) )
              DDSigE(1) = (D1 - D2) * DDEps(1) + D2DDEpsV
              DDSigE(2) = (D1 - D2) * DDEps(2) + D2DDEpsV
              DDSigE(3) = (D1 - D2) * DDEps(3) + D2DDEpsV
              do I = 4, NTENSOR
                DDSigE(I) = GG * DDEps(I)
              end do
              SigEQ = SigC + DDSigE
            else
              SigEQ = SigE 
            end if
          else ! ImplicitQuasiStatic calculation
          
          ! for ExplicitDynamic integration
              SigEQ = SigE ! equilibrium stress SigEQ in Mohr-Coulomb model
          end if 
          
          ! store elastic stress SigE in current stress SigC
          SigC = SigE
     
          ! begin stress correction to bring stress back on yield surface
          TauMax = 1
          ITens = 0
          IApex = 0
          select case (SoilModel)
            case (ESM_FRICTIONAL_LIQUID, ESM_BINGHAM_LIQUID)
              call SetIPL(IntGlo, IEl, 0)
            end select

          ! for using objective stress definition
          if (CalParams%ApplyObjectiveStress) then ! Consider large deformation terms
            call Hill(IEl, ELEMENTNODES, IncrementalDisplacementSoil(1:Counters%N, IEntityID),  &
                        ReducedDof, ElementConnectivities, BMatrix, Sig0, SigC, DEpsVol)
          end if ! objective stress  
            
          call CalculatePrincipalStresses(IntGlo, SigC, SigPrin)
          
          DSig = SigC - Sig0
          DSigP = SigE - SigC
          call MatVec(DI, NTENSOR, DSigP, NTENSOR, DEpsP) ! DEpsP = DI * DSigP

          ! Calculate stresses and incremental strains of current time step
          call AssignStressStrainToGlobalArray(IntGlo, NTENSOR, DSigWP, DSigGP, DSig, SigPrin, DEps)

          if (CalParams%ApplyBulkViscosityDamping) then
            RateVolStrain(IEl) = DEpsVol / CalParams%TimeIncrement
            call CalculateViscousDamping(IntGlo, IEl)
          end if

          ! determine plastic points
          if (GetIPL(IntGlo)/=0) then
            CalParams%IntegrationPointData%NPlasticPoints = CalParams%IntegrationPointData%NPlasticPoints + 1
          end if
          if (GetIPL(IntGlo)<0) then
            CalParams%IntegrationPointData%NNegativePlasticPoints = CalParams%IntegrationPointData%NNegativePlasticPoints + 1
          end if
          ! determine tension cut-off points
          if (ITens/=0) then 
            CalParams%IntegrationPointData%NTensionCutOffPoints = CalParams%IntegrationPointData%NTensionCutOffPoints + 1
          end if
          ! determine apex points
          if (IApex/=0) then
            CalParams%IntegrationPointData%NApexPoints = CalParams%IntegrationPointData%NApexPoints + 1
          endif

          ! convergence check for ImplicitQuasiStatic calculation
          if (CalParams%ApplyImplicitQuasiStatic.and.(GetIPL(IntGlo)/=0)) then
            ErrorSig = SigC - SigEQ
            LocalError = 0.0
            do I = 1, NTENSOR
              LocalError = LocalError + ErrorSig(I) * ErrorSig(I)
            end do
            CalParams%ConvergenceCheck%SumLocalError = CalParams%ConvergenceCheck%SumLocalError+0.5*dsqrt(LocalError)/TauMax*WtN
            CalParams%ConvergenceCheck%SumIntegrationPointWeights = CalParams%ConvergenceCheck%SumIntegrationPointWeights + WtN
            ToleratedLocalError = 4 * (CalParams%ToleratedErrorForce * TauMax) * (CalParams%ToleratedErrorForce * TauMax)
            if (LocalError > ToleratedLocalError) then
              CalParams%ConvergenceCheck%NInaccuratePlasticPoints = CalParams%ConvergenceCheck%NInaccuratePlasticPoints + 1
            end if
          end if ! convergence check for ImplicitQuasiStatic calculation
          
        end subroutine GetStressesSoilMaterialPoint !SUBROUTINE NOT USED. IT SHOULD BE REMOVED
        
        subroutine GetFluidPressure(IDpt)
        !*************************************************************************************
        !  Function: compute liquid and gas pressures. Solve mass banance equations
        !            Valid for single-point formulation, 2phase and 3phase analysis
        !
        !*************************************************************************************
        
        implicit none
        
        ! local variables
        integer(INTEGER_TYPE), intent(in) :: IDpt ! global integration/material point number

        integer(INTEGER_TYPE) :: I     ! counter
        integer(INTEGER_TYPE) :: IDset ! ID of material parameter set
        integer(INTEGER_TYPE) :: ntens ! Dimension of stress vector to pass to UDSM 

        real(REAL_TYPE), dimension(MatParams(MaterialIDArray(IDpt))%UMATDimension) :: StrainIncr ! strain increment in integration/material point
        real(REAL_TYPE), dimension(NTENSOR) :: TempStrainIncr                                    ! incremental strain vector assigned to point 

        real(REAL_TYPE) :: DSigWP    ! Change of water pressure at integration/material point 
        real(REAL_TYPE) :: DSigGP    ! Change of gas pressure at integration/material point 
        real(REAL_TYPE) :: Kf, WaterAdvectiveFlux, WaterDensity, g   !Bulk of water, Advective flow, water density, gravity acceleration
        real(REAL_TYPE) :: N         ! Porosity         
        real(REAL_TYPE) :: DEpsVol   ! Incremental volumetric strain
        real(REAL_TYPE) :: DEpsVolW  ! Incremental volumetric strain (water)
        real(REAL_TYPE) :: dT        ! Change of Temperature at integration point
        real(REAL_TYPE) :: lambda    ! Used to compute water pressure
        real(REAL_TYPE) :: Sr, dSrdp(1) !Degree of saturation, derivative dSr/dp_w
      
       if (.not.(((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1))) RETURN !Valid for single-point formulation, 2phase and 3phase analysis
    
       IDset = MaterialIDArray(IDpt) ! get constitutive model in integration/material point. It is the material number stored in $$MATERIAL_INDEX in the GOM-file
       
       ntens = MatParams(IDset)%UMATDimension     
          
       TempStrainIncr = GetEpsStep(Particles(IDpt)) ! get strain increments in integration/material point. It is the incremental strain vector assigned to point   
        
       StrainIncr = 0.0
       do I=1, NTENSOR
          StrainIncr(I) = StrainIncr(I) + TempStrainIncr(I)
       enddo 
        
       DEpsVol = StrainIncr(1) + StrainIncr(2) + StrainIncr(3) ! volumetric strain, valid for 2D and 3D      
             
       ! initalise water pressure
       DSigWP = 0.0d0
       DSigGP = 0.0d0

       ! for 2-phase analysis
       if ((CalParams%NumberOfPhases==2).and.(NFORMULATION==1)) then
        if (Particles(IDpt)%WaterWeight > 0.0) then
         DEpsVolW = Particles(IDpt)%WaterVolumetricStrain ! Water phase
         N = Particles(IDpt)%Porosity
         Kf = Particles(IDpt)%BulkWater
        
         if (CalParams%ApplyPartialSaturation) then
           Sr = Particles(IDPt)%DegreeSaturation
           WaterAdvectiveFlux = Particles(IDpt)%WaterAdvectiveFlux        
           g = CalParams%GravityData%GAccel                    !Gravity (m/s2)
           WaterDensity = Particles(IDpt)%WaterWeight/g
           call CalculateDerivDegreeSaturation(IDPt,dSrdp,1)
        
           if ((Sr<1).and.(Sr>0)) then              !(consider advective terms)
             lambda = N*Sr*WaterDensity/Kf - N*WaterDensity * dSrdp(1)
             lambda = 1/lambda
             DSigWP = lambda*((Sr*WaterDensity)*DEpsVol + WaterAdvectiveFlux)   
           else
             lambda = Kf/N
             DSigWP = lambda*( N * DEpsVolW + (1.0d0 - N) * DEpsVol)  
           end if
        
         else
           lambda = Kf/N
           DSigWP = lambda*( N * DEpsVolW + (1.0d0 - N) * DEpsVol)
         end if
        
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
       if (CalParams%NumberOfPhases==3) then ! the incremental water pressure, gas pressure and temperature (for one material point) are obtained
         call SolveBalanceEquations(IDpt, DEpsVol, DSigWP, DSigGP, dT) ! solve mass balances (for unsaturated porous media) and energy balance
       end if ! 3-phase analysis

        ! single phase quasi-static consolidation
        if (CalParams%ApplyImplicitQuasiStatic) then
          DsigWP = Particles(IDpt)%WaterPressure - Particles(IDpt)%WaterPressure0
        end if
    
        call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP)  !Note that the subroutine checks Cavitation Threshold & Gas Pressure
      
        end subroutine GetFluidPressure

    subroutine SmoothenLiquidPressureIncrement(IMatSet)
    !**********************************************************************
    !
    !    Function: Somoothen pressure increment of water phase
    !
    !    IMatSet : The number of the material set
    !
    !**********************************************************************
       
        implicit none
 
          integer(INTEGER_TYPE), intent(in) :: IMatSet
          ! Local variables  
          real(REAL_TYPE), dimension(Counters%NodTot,2) :: dPressureSmooth
          real(REAL_TYPE) :: dWP, Weight, dWPWeighted, EnhancedPressureIncrement
          integer(INTEGER_TYPE) :: IEl, IAEl, I, NElemPart, MaterialID, IParticle, ParticleIndex
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: LJ
        
          dPressureSmooth = 0.0

          !Calculate nodal incremental pressure mapping from MP to element nodes    
          do IAEl = 1, Counters%NAEl  ! active element
            IEl = ActiveElement(IAEl)
            LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
            if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
              NElemPart = NumberOfIntegrationPoints(IEl)
              do IParticle = 1, NElemPart
                ParticleIndex = GetParticleIndex(IParticle, IEl)
                MaterialID = MaterialIDArray(ParticleIndex)
                
                if (MaterialID==IMatSet)then ! This particle belongs to the considered entity (Material)
                  ! The volumetric strain in this element for the considered entity (Material)
                  dWP = Particles(ParticleIndex)%WaterPressure - Particles(ParticleIndex)%WaterPressure0
                  Weight = Particles(ParticleIndex)%IntegrationWeight
                  dWPWeighted = dWP * Weight
                  dPressureSmooth(LJ(:), 1) = dPressureSmooth(LJ(:), 1) + dWPWeighted
                  dPressureSmooth(LJ(:), 2) = dPressureSmooth(LJ(:), 2) + Weight
                end if
              end do ! particles

            end if ! entity
          end do ! active element     

          do I = 1, Counters%NodTot
            if(dPressureSmooth(I, 2)>0.0) then
              dPressureSmooth(I, 1) = dPressureSmooth(I, 1) / dPressureSmooth(I, 2) ! Only consider active nodes
            end if
          end do
       
          do IAEl = 1, Counters%NAEl  ! active element 
            IEl = ActiveElement(IAEl)
              if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
                LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
                EnhancedPressureIncrement = 0.0
                do I = 1, ELEMENTNODES
                  EnhancedPressureIncrement = EnhancedPressureIncrement + dPressureSmooth(LJ(I), 1) / dble(ELEMENTNODES)
                end do    

              ! Update particle pressure

                NElemPart = NumberOfIntegrationPoints(IEl)
                do IParticle = 1, NElemPart ! Loop over all particles of the element
                  ParticleIndex = GetParticleIndex(IParticle, IEl)
                  MaterialID = MaterialIDArray(ParticleIndex)
                    if (MaterialID==IMatSet)then ! this particle belongs to material IMatSet. Assign its strains
                      Particles(ParticleIndex)%WaterPressure = Particles(ParticleIndex)%WaterPressure0 + EnhancedPressureIncrement
                    end if
                end do ! particles
              end if ! material
    end  do ! active element
    
    end subroutine SmoothenLiquidPressureIncrement
    
        subroutine CalculateViscousDamping(ParticleID, IEl)
        !**********************************************************************
        !> Computes a pressure term introducing bulk viscosity damping to the equation of motion.
        !>
        !! \param[in] ParticleID ID of considered material point.
        !! \param[in] IEl ID of element of the considered material point.
        !! \param[in] DilationalWaveSpeed Current wave speed computed for the considered material point.
        !**********************************************************************
        
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

        end subroutine CalculateViscousDamping
        
        
        end module ModMPMDYNStresses
