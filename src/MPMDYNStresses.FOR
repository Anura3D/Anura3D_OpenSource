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
      !     $Revision: 9262 $
      !     $Date: 2021-04-29 05:31:48 +0200 (do, 29 apr. 2021) $
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

       
      contains
     
        subroutine MPMDYNGetSig()
        !**********************************************************************
        !
        !  Function:  Loops over the active elements and updates the stresses and
        !             plasticity state of Gauss points for fully filled elements and
        !             material points for partially filled elements by calling the appropriate
        !             constitutive model routine of the material assigned to a Gauss point
        !             or material point.
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
                     NElemPart,    & ! Number of Gauss points or particles per element &
                     Int,          & ! Local integration point counter &
                     IMatSet,      & ! Counter  &
                     IDof,         & ! Counter &
                     NMaterialSets,& ! Number of material sets &
                     IEntity,      &
                     I,     &        ! Counter
                     IDMaterial      ! ID of Material Set
          real(REAL_TYPE) :: WTN,         & ! Weight of Gauss point &
                              ShearStress, & ! shear stress in a material point    &
                              YieldStress    !yield stress for definition of plastic flow in bingham or frictional fluid    

          real(REAL_TYPE), dimension(Counters%N) :: DDisp, DDispLiquid, VelLiquid
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: BMatrix
          logical :: DoSkipStressComputation

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
          
          if (CalParams%ApplyStrainSmoothing) then
            ! initialise smoothened strain or strain rate data 
            NMaterialSets = Counters%NLayers
            do IMatSet = 1, NMaterialSets

              IEntity =getMaterialEntity(IMatSet)

              if (NFORMULATION==1) then
                if (.not.CalParams%ApplyContactAlgorithm) then ! no contact
                  do IDof  = 1, Counters%N
                    ! if no contact is used then the system always has one shared displacement
                    DDisp (IDof) = IncrementalDisplacementSoil(IDof,1)
                  end do
                else if (CalParams%ApplyContactAlgorithm) then ! contact is used
                  do IDof  = 1, Counters%N
                    ! if contact is used then each entity has its own vector 
                    DDisp (IDof) = IncrementalDisplacementSoil(IDof,IEntity)
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

            end do
          end if ! strain smoothing
              
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
                !call SmoothenStrainsWater2LayerForm()
              else
                call SmoothenStrainsLiquid(DDispLiquid, IMatSet) ! strain smoothing for water
              end if
            end do
          end if ! strain smoothing Liquid 2 layer
              
          if ((NFORMULATION==1).and.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and. &
                    CalParams%ApplyStrainSmoothing) then  
            ! initialise smoothened strain data -- later should be another flag for water
            NMaterialSets = Counters%NLayers 
            do IMatSet = 1, NMaterialSets
              ElementStrain = 0.0
              call SmoothenStrainsWater(IMatSet)
            end do ! NMaterialSets

            if (CalParams%NumberOfPhases==3) then
              do IMatSet = 1, NMaterialSets
                ElementStrain = 0.0
                call SmoothenStrainsGas(IMatSet)
              end do ! NMaterialSets
            end if
          end if

          !$omp parallel do private(IEl, IAEl, Int, IntGlo, NElemPart, BMatrix, WtN, IDMaterial, ShearStress, DoSkipStressComputation)
          do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)
            
            if (CalParams%ApplyObjectiveStress.or.CalParams%ApplyImplicitQuasiStatic) then 
              call GetBMatrix(IEl, 1, NodalCoordinatesUpd, BMatrix, WtN)
            end if
            
            ! Get the number of integration per element (Gaussian integration points or number of material points)
            NElemPart = NumberOfIntegrationPoints(IEl) 
            

          
            ! Calculate stresses in integration/material points
            do Int = 1, NElemPart ! Loop over all integration/material points of the element
         
              ! Determine global ID of integration point 
              IntGlo = GetParticleIndex(Int, IEl)
              IDMaterial = MaterialIDArray(IntGlo)
              
              
               if ( (MaterialPointTypeArray(IntGlo)==MaterialPointTypeMixture).or.(MaterialPointTypeArray(IntGlo)==MaterialPointTypeSolid) ) then
                if (CalParams%ApplyContactAlgorithm) then
                  IEntity = EntityIDArray(IntGlo) 
                else
                  IEntity = 1 ! ... which might not always be correct when considering different materials?
                end if
                
                if((NFORMULATION==1).or.(Particles(IntGlo)%PhaseStatus==PhaseStatusSOLID)) then 

                  DoSkipStressComputation = EntityIDArray(IntGlo) == PILE_ENTITY
                  do I = 1, NVECTOR
                    DoSkipStressComputation = DoSkipStressComputation .and. CalParams%ApplyPrescribedVelocity(I)
                  end do
				  
				  if (MatParams(IDMaterial)%MaterialModel==ESM_RIGID_BODY) then
						DoSkipStressComputation = .true.
                  end if
                  
                  if (.not.DoSkipStressComputation) then  
                    ! calculate stresses for mixture or solid material points
                    ! Note: needs to be checked (NGGS pilot - UMAT interface)
                      call StressSolid(IntGlo, IEl, BMatrix, IEntity)

                       
                  else
                    do I = 1, NVECTOR ! only first three indices, both for 2D and 3D
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
                    call GetPressureLiquidMaterialPoint(IntGlo, IEl) ! calculate pressure for liquid material points
                   end if
               end if
               
              
            end do
          end do

          !$omp end parallel do
          
          do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)
            
            if (CalParams%ApplyObjectiveStress.or.CalParams%ApplyImplicitQuasiStatic) then 
              call GetBMatrix(IEl, 1, NodalCoordinatesUpd, BMatrix, WtN)
            end if
            
            ! Get the number of integration per element (Gaussian integration points or number of material points)
            NElemPart = NumberOfIntegrationPoints(IEl) 
          
            ! Calculate stresses in integration/material points
            do Int = 1, NElemPart ! Loop over all integration/material points of the element
         
              ! Determine global ID of integration point 
              IntGlo = GetParticleIndex(Int, IEl)
              IDMaterial = MaterialIDArray(IntGlo)
              
               if (CalParams%ApplyContactAlgorithm) then
                IEntity = EntityIDArray(IntGlo)  ! Get the entity ID for the current particle          
               else
                IEntity = 1 ! ... which might not always be correct when considering different materials?
               end if
             if(.not.(CalParams%ApplyBinghamFluid.or.CalParams%ApplyFrictionalFluid & !v2016
              .or. (MatParams(IDMaterial)%MaterialModel==ESM_FRICTIONAL_LIQUID) & !v2017.1 and following
              .or. (MatParams(IDMaterial)%MaterialModel==ESM_BINGHAM_LIQUID))) then
              !for Bingham and frictional fluid stresses are calculated elsewhere
               if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then
                  call GetStressesLiquidMaterialPoint(IntGlo, IEl) ! calculate stress for liquid material points
               end if
              end if 
             end do
          end do
        
        end subroutine MPMDYNGetSig


        
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
