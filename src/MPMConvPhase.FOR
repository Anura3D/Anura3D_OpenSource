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
	

      module ModConvectivePhase
      !**********************************************************************
      !
      !    Function:  This module contains all routines which are used by both quasi-static 
      !               and dynamic MPM in the convective phase of the calculation process.
      !
      !               In order to keep the size of this source file reasonably small,
      !               this module only contains routines that are directly related to
      !               the updating of particle data that are identical for the quasi-static
      !               and dynamic MPM.
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 9767 $
      !     $Date: 2022-06-20 14:43:49 +0200 (ma, 20 jun 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluation
      use ModMPMData
      use ModMeshInfo
      use ModMPMMeshAdjustment
      use ModString
      use ModIsort
      use ModGlobalConstants
      use ModFileIO
      
      implicit none

      contains ! Routines of this module

        !----------------------------------------------------------------------

      
      subroutine GetWaveSpeedForTimeStepSize(IParticle, IEl, WaveSpeed)
      !**********************************************************************
      !
      !  Function: calculate wave speed for particle IParticle in element IEl
      !            calculation differs for 1-phase and 2-phase.
      !            mass scaling is taken into account.
      !
      !**********************************************************************
      
      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: IParticle, IEl
        real(REAL_TYPE), intent(inout) :: WaveSpeed
       
        ! local variables
        integer(INTEGER_TYPE) :: MaterialIndex
        real(REAL_TYPE) :: EUnloading
        real(REAL_TYPE) :: G, Nu, Density, Bulk
        real(REAL_TYPE) :: consolidationcoeff, WaveSpeed_cons, WaveSpeed_dyn

        real(REAL_TYPE) :: porosity, rho_solid, rho_liquid, rho_mixture, gravity,  &
                            hydraulic_conductivity, bulk_modulus, dx, Sr, dSrdp(1)
        real(REAL_TYPE) :: a_crit, b_crit, d_crit
        character(len=64) :: SoilModel ! name of the constitutive model
        logical :: IsUndrEffectiveStress = .false. 

        WaveSpeed = 0.0
        consolidationcoeff = 0.0
        WaveSpeed_dyn = 0.0
        WaveSpeed_cons = 0.0

        Sr = 1.0
        dSrdp(1) = 0.0

        MaterialIndex = MaterialIDArray_N(IParticle)
        SoilModel = MatParams(MaterialIndex)%MaterialModel
        ! name of constitutive model as specified in GOM-file
        
        IsUndrEffectiveStress = &
        !code version 2016 and previous
        ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
        !code version 2017.1 and following
        (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))


        ! calculate stiffness of soil
          G = MatParams(MaterialIndex)%ShearModulus !kPa
          Nu = MatParams(MaterialIndex)%PoissonRatio
          EUnloading = ( 2*G*(1-Nu))/ (1-2*Nu) !kPa  ! Eoed

          if (SoilModel==ESM_MODIFIED_CAM_CLAY) then ! MCC
            Eunloading = Particles(IParticle)%ESM_UnloadingStiffness
          end if                    
          
          if (SoilModel == ESM_HYPOPLASTICITY_SAND .or. &
              SoilModel==ESM_EXTERNAL_SOIL_MODEL) then ! user-defined soil model
            Eunloading = Particles(IParticle)%ESM_UnloadingStiffness
          end if


        !if (CalParams%ApplyEffectiveStressAnalysis) then
        if (IsUndrEffectiveStress) then
          if (Particles(IParticle)%Porosity>0.0) then
            Bulk = Particles(IParticle)%BulkWater / Particles(IParticle)%Porosity
            EUnloading = EUnloading + Bulk
          end if
        end if
               
        if ( (MatParams(MaterialIndex)%MaterialType=='1-phase-liquid') .or. (MatParams(MaterialIndex)%MaterialPhases=='1-phase-liquid') ) then ! if material point is free water
          Density = MatParams(MaterialIndex)%DensityLiquid / 1000.0
          EUnloading = MatParams(MaterialIndex)%BulkModulusLiquid
        else
          !if (CalParams%ApplyEffectiveStressAnalysis.or.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3))) then
           if (IsUndrEffectiveStress.or.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3))) then
            Density = MatParams(MaterialIndex)%DensityMixture / 1000.0
          else
            Density = (1 - MatParams(MaterialIndex)%InitialPorosity) * MatParams(MaterialIndex)%DensitySolid / 1000.0
          end if
        end if
        
        if (((CalParams%NumberOfPhases==2).and.(.not.(IsUndrEffectiveStress))) &
          .and.(.not.(MatParams(MaterialIndex)%MaterialType.eq.'dry_material')).or.(MatParams(MaterialIndex)%MaterialType.eq.'3-phase') &
          .or.(MatParams(MaterialIndex)%MaterialType.eq.'unsaturated_material_3phase_coupled')) then ! Excludes 'dry_material' for multi material analyses
          ! get necessary parameters
          porosity = MatParams(MaterialIndex)%InitialPorosity
          rho_solid = MatParams(MaterialIndex)%DensitySolid / 1000.0
          rho_liquid = MatParams(MaterialIndex)%DensityLiquid / 1000.0
          rho_mixture = MatParams(MaterialIndex)%DensityMixture / 1000.0
          gravity = CalParams%GravityData%GAccel
          !hydraulic_conductivity = MatParams(MaterialIndex)%HydraulicConductivityLiquid
          bulk_modulus = MatParams(MaterialIndex)%BulkModulusLiquid
          hydraulic_conductivity = Particles(IParticle)%Conductivity
          If (CalParams%ApplyPartialSaturation) then
            Sr = Particles(IParticle)%DegreeSaturation
            rho_mixture = (1-porosity)*rho_solid + Sr*porosity*rho_liquid
            call CalculateDerivDegreeSaturation(IParticle,dSrdp,1)
          end if
          
         ! consolidationcoeff = hydraulic_conductivity/(gravity* rho_liquid*(1/EUnloading + Sr*porosity/bulk_modulus-porosity*dSrdp(1)))
          consolidationcoeff = hydraulic_conductivity/(gravity* rho_liquid*(1/EUnloading + porosity/bulk_modulus))
          
          dx = ElementLMin(IEl)

          ! calculate a,b and d for wavespeed
          a_crit = porosity * rho_mixture * gravity / (1.0 - porosity) / rho_solid / hydraulic_conductivity
          b_crit = 4.0 * ( porosity * rho_mixture * bulk_modulus + (1.0 - 2.0 * porosity) * rho_liquid * bulk_modulus &
                   + porosity * rho_liquid * EUnloading) / porosity / (1 - porosity) / rho_solid / rho_liquid / dx**2
          d_crit = 16.0 * EUnloading * bulk_modulus / (1.0 - porosity) / rho_solid / rho_liquid / dx**4

          ! include mass scaling:
          a_crit = a_crit / CalParams%ScalingMassFactor
          b_crit = b_crit / CalParams%ScalingMassFactor
          d_crit = d_crit / (CalParams%ScalingMassFactor**2)
          
          ! calculate wavespeed
          WaveSpeed_dyn = dx / (- 2.0 * a_crit + sqrt(4.0 * a_crit**2 + 8.0 * (b_crit + sqrt(b_crit**2 - 4.0 * d_crit)))) &
                      * (b_crit + sqrt(b_crit**2 - 4.0 * d_crit))
          
 
          WaveSpeed_cons = 2*consolidationcoeff/dx

          
          WaveSpeed = max(WaveSpeed_dyn,WaveSpeed_cons)
        else
          ! calculate wavespeed^2
          WaveSpeed = EUnloading / Density
          
          ! include mass scaling:
          WaveSpeed = WaveSpeed / CalParams%ScalingMassFactor
          
          ! calculate wavespeed
          WaveSpeed = sqrt(WaveSpeed)
        end if

      end subroutine GetWaveSpeedForTimeStepSize
      
      subroutine GetWaveSpeed(IParticle, WaveSpeed)
      !**********************************************************************
      !
      !  Function: calculate wave speed for particle IParticle in element IEl
      !            calculation differs for 1-phase and 2-phase.
      !            mass scaling is taken into account.
      !
      !**********************************************************************
      
      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: IParticle
        real(REAL_TYPE), intent(inout) :: WaveSpeed
       
        ! local variables
        integer(INTEGER_TYPE) :: MaterialIndex
        real(REAL_TYPE) :: EUnloading
        real(REAL_TYPE) :: G, Nu, Density, Bulk
        
        character(len=64) :: SoilModel ! name of the constitutive model
        logical :: IsUndrEffectiveStress = .false. 

        WaveSpeed = 0.0

        MaterialIndex = MaterialIDArray_N(IParticle)
        SoilModel = MatParams(MaterialIndex)%MaterialModel ! name of constitutive model as specified in GOM-file
        
        IsUndrEffectiveStress = &
        !code version 2016 and previous
        ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
        !code version 2017.1 and following
        (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))


        ! calculate stiffness of soil
          G = MatParams(MaterialIndex)%ShearModulus !kPa
          Nu = MatParams(MaterialIndex)%PoissonRatio
          EUnloading = ( 2*G*(1-Nu))/ (1-2*Nu) !kPa  ! Eoed

          if (SoilModel==ESM_MODIFIED_CAM_CLAY) then ! MCC
            Eunloading = Particles(IParticle)%ESM_UnloadingStiffness
          end if          
          
          if (SoilModel==ESM_EXTERNAL_SOIL_MODEL .or. &
              SoilModel == ESM_HYPOPLASTICITY_SAND) then ! user-defined soil model
            Eunloading = Particles(IParticle)%ESM_UnloadingStiffness
          end if


        if (IsUndrEffectiveStress.or.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3))) then
          if (Particles(IParticle)%Porosity>0.0) then
            Bulk = Particles(IParticle)%BulkWater / Particles(IParticle)%Porosity
            EUnloading = EUnloading + Bulk
          end if
        end if

        if ( (MatParams(MaterialIndex)%MaterialType=='1-phase-liquid') .or. (MatParams(MaterialIndex)%MaterialPhases=='1-phase-liquid') ) then ! if material point is free water
          Density = MatParams(MaterialIndex)%DensityLiquid / 1000.0
          EUnloading = MatParams(MaterialIndex)%BulkModulusLiquid
        else

        if (IsUndrEffectiveStress.or.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3))) then
          Density = MatParams(MaterialIndex)%DensityMixture / 1000.0
        else
          Density = (1 - MatParams(MaterialIndex)%InitialPorosity) * MatParams(MaterialIndex)%DensitySolid / 1000.0
        end if
        end if

          ! calculate wavespeed^2
          WaveSpeed = EUnloading / Density
         
          ! calculate wavespeed
          WaveSpeed = sqrt(WaveSpeed)


      end subroutine GetWaveSpeed
      

        subroutine UpdateParticleStrains(IncrementalDisplacementSoil, DShapeValuesArray, &
            ActiveElement, NPartEle)!, NAEl)!, GetParticleIndex)!DShapeValues) ! all these variables need to be local 
        !**********************************************************************
        !
        !    Function:  Updates the total and incremental particle strains before mesh resetting
        !
        !**********************************************************************
       
          implicit none
        
          ! local variables
          real(REAL_TYPE), dimension(Counters%N) :: DUTotEnt
          real(REAL_TYPE), dimension(NTENSOR) :: Eps
          real(REAL_TYPE), dimension(NDIM, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac
          integer(INTEGER_TYPE) :: IElement, NElemPart, IParticle, ParticleIndex, I, iEntity, iEntityDefault, IActiveElement

          ! implicit scheme - shape values array inout 
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
          
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: IncrementalDisplacementSoil
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: TotalDisplacementSoil
          
          ! we input local variables into this subroutine
          integer(INTEGER_TYPE), dimension(Counters%NAEl), intent(in) :: ActiveElement !list of IDs of active elements
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle !number of particles in an element
          !integer(INTEGER_TYPE) :: NAEl ! number of active elements 
          !integer(INTEGER_TYPE), dimension(:,:), intent(in) :: GetParticleIndex ! we need this as an input
          
          ! set the nodal displacement increments to the first entity
          ! only need to change if the contact model is used and the next particle belongs to a different entity
          iEntityDefault = 1 ! set default to first entity
          do I = 1,Counters%N
            DUTotEnt(I) = IncrementalDisplacementSoil(I,iEntityDefault)
          end do

          do IActiveElement = 1, Counters%NAEl!Counters%NAEl ! loop over all active elements
            
              IElement = ActiveElement(IActiveElement) ! get element number of the active element
              NElemPart = NPartEle(IElement) ! get number of particles in the active element
              
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if  
              
              do IParticle = 1, NElemPart ! loop over all particles of the active element
              
                  ! there needs to be GetParticleIndex "local for iteration" and "actual"
                ParticleIndex = GetParticleIndex_N(IParticle, IElement) ! get global particle ID
                
                if( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then 

                  ! determine B matrix of deformed element at the particle local position
                  call BMatrix(Particles(ParticleIndex)%LocPos, &
                               ELEMENTNODES, Counters%NEl, &
                               Counters%NodTot, NDIM, &
                               IElement, ElementConnectivities, &
                               NodalCoordinatesUpd,  &
                               BMatrixDeformed, DetJac)!, & 
                              ! DShapeValuesArray)
                  
                  ! get the nodal displacement increments for the specific entity
                  !if (CalParams%ApplyContactAlgorithm) then ! CONTACT IS APPLIED
                  !  iEntity = EntityIDArray(ParticleIndex) ! get the entity ID for the current particle
                  !  if (iEntity /= iEntityDefault) then ! the particle belongs to a different entity from the previous entity
                  !    DUTotEnt(:) = IncrementalDisplacementSoil(1:Counters%N,iEntity)
                  !    iEntityDefault = iEntity !set the entity ID for the next particle loop
                  !  end if
                  !end if

                  ! Eps in 3D (Exx, Eyy, Ezz, Gxy, Gyz, Gzx), in 2D (Exx, Eyy, Ezz, Gxy)
                  call Get_Strain(IElement, IParticle, &
                                ElementConnectivities, & 
                                BMatrixDeformed, DUTotEnt, ReducedDof, & !--> DUTotEnt is where IncrementalDisplacementSoil
                                Eps) ! Eps = strain increment at material point
     
                  ! Update particle strains
                  call SetEpsStep(Particles(ParticleIndex), Eps)
                  ! we do not want to increase the total particle. We calculate the incremental strain and then use it to find the stress.  
                  ! skip increasing the total strain of the particle and leave it for later after we achieve the solution
                  !if (.not.CalParams%ImplicitIntegration%IsImplicitScheme) then !CalParams%ApplyImplicitQuasiStatic
                  !  call IncreaseEps(Particles(ParticleIndex), Eps)
                  !end if

                end if ! Material Point can be MIXTURE or SOLID

              end do ! particle loop

          end do ! active element loop

        end subroutine UpdateParticleStrains
        
        
        subroutine UpdateParticleStrainsLiquidTwoLayer()
        !**********************************************************************
        !
        !    Function:  Updates the total and incremental particle strains 
        !               before mesh resetting.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
          implicit none
        
          !Local variables
          !!CC - nodal displacement increment for a specific entity - passed to Get_Strain
          real(REAL_TYPE), dimension(Counters%N) :: DULiqEnt

          real(REAL_TYPE), dimension(NTENSOR) :: Eps
          real(REAL_TYPE), &
            dimension(NDIM, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac
          integer(INTEGER_TYPE) :: IElement, NElemPart, IParticle, &
                     ParticleIndex,iEntityDefault    !!CC

          !--------------------------------------------------------

          !========================================================             !!CC
          !set the nodal displacement increments to the first entity
          !only need to change if the contact model is used and 
          !the next particle belongs to a different entity
          iEntityDefault = 1    !set default to first entity
          DULiqEnt = IncrementalDisplacementWater(:,iEntityDefault)
          !=======================================================

         !----------------------------------------------------
          do IElement = 1, Counters%NEl     !loop over all elements
            if (IsActiveElement_N(IElement)) then  !active elements only
              
                NElemPart = NPartEle_N(IElement)    !number of particles in element
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if
              
              do IParticle = 1, NElemPart       !loop over all particles of the element
                ParticleIndex = GetParticleIndex_N(IParticle, IElement)   !particle global ID
                if((MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid).and. &
                         (.not.(NFORMULATION==1))) then !Material Point LIQUID

                  ! Determine B matrix of deformed element at the particle local position
                  call BMatrix(Particles(ParticleIndex)%LocPos, &
                             ELEMENTNODES, Counters%NEl,  &
                             Counters%NodTot, NDIM, &
                             IElement, ElementConnectivities, &
                             NodalCoordinatesUpd, &
                             BMatrixDeformed, DetJac)

                  ! Eps in 3D (Exx, Eyy, Ezz, Gxy, Gyz, Gzx), in 2D (Exx, Eyy, Ezz, Gxy)
                  call Get_Strain(IElement, IParticle, &
                                ElementConnectivities, & 
                                BMatrixDeformed, DULiqEnt, ReducedDof, &
                                Eps)

                  ! Update particle strains
                  call SetEpsStep(Particles(ParticleIndex), Eps)
                  call IncreaseEps(Particles(ParticleIndex), Eps)

                  ! Update Volumetric Strain Water
                  Particles(ParticleIndex)%WaterVolumetricStrain =  getEpsV(Eps)

                end if !Material Point LIQUID
              end do
            end if
          end  do

        end subroutine UpdateParticleStrainsLiquidTwoLayer


        subroutine UpdateParticlePos(ShapeValuesArray_N, DShapeValuesArray_N, &
                                     ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
                                     GlobPosArray_N, GlobPosArray_NPlus1, &
                                     UStepArray, UPhaseArray, UArray, &
                                     ElementIDArray, EleParticles, &
                                    IncrementalDisplacementSoil)
        !**********************************************************************
        !
        !    Function:  Updates the location of particles. Connected to this task are:
        !
        !               - calculate particle displacements from nodal displacements
        !               - calculate new particle global position
        !               - determine into which elements the particles moved
        !               - calculate new particle local position
        !               - update particle shape values
        !
        !   Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: I, IParticle, ElementID, NewElementID
          real(REAL_TYPE), dimension(NVECTOR) :: LocPos, GlopPos0
          logical :: Success, IsInside
          real(REAL_TYPE) :: diff, sumLocPos
           
          ! implicit scheme - shape values array inout 
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray_N 
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray_N
		 
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(out) :: ShapeValuesArray_NPlus1 
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(out) :: DShapeValuesArray_NPlus1
		 
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: IncrementalDisplacementSoil
          
          
          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(inout) :: GlobPosArray_N
                 
          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(inout) :: GlobPosArray_NPlus1

          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(inout) :: UStepArray 
          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(inout) :: UPhaseArray 
          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(inout) :: UArray

         integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(inout) :: ElementIDArray
          
         integer(kind = 8), dimension(Counters%NParticles), intent(inout) :: EleParticles !INTEGER_TYPE
          
          
           do IParticle = 1, Counters%NParticles ! Loop over particles

            
               !Particles(IParticle)%LocPos_NPlus1 = Particles(IParticle)%LocPos
               
               
               ElementID = ElementIDArray(IParticle) !Element ID determined from ElementIDArray
            
            

            ! Calculate particle displacements from nodal displacements
            call CalcParticleDisplacements(IParticle, ElementID, ShapeValuesArray_N, IncrementalDisplacementSoil, &
                UStepArray, UPhaseArray, UArray)

            GlopPos0 =  GlobPosArray_N(IParticle,:)  ! store the initial global position of the particle

            ! Calculate new particle global positions
            GlobPosArray_NPlus1(IParticle,:) = GlobPosArray_N(IParticle,:) + UStepArray(IParticle,:)

           if (CalParams%ApplyFixParticlesK0) then ! only 3D functionality
             !if (GlobPosArray(IParticle,1) >= (CalParams%Fix(1))) then ! allow sliding along y-axis only
             !    GlobPosArray(IParticle,1) = GlopPos0(1)
             !    GlobPosArray(IParticle,3) = GlopPos0(3)
             !end if
             !if (GlobPosArray(IParticle,2) <= (CalParams%Fix(2))) then
             !   GlobPosArray(IParticle,:) = GlopPos0 ! particles along the bottom are fixed
             !end if
           end if ! ApplyFixParticlesK0

          IsInside = IsInBoundingBox(GlobPosArray_N(IParticle,:))
          ! Check whether particle moved outside the mesh
          if (.not.IsInside ) then
             if (NDIM==3) call GiveError('  Particle '// trim(String(IParticle)) // & 
                            ' from element '// trim(String(ElementID)) // ' lies '// &
                            'outside the mesh. '// &
                            'X, Y, Z ='// trim(String(GlobPosArray_N(IParticle,1))) // trim(String(GlobPosArray_N(IParticle,2))) // trim(String(GlobPosArray_N(IParticle,3))) )
             if (NDIM==2) call GiveError('  Particle '// trim(String(IParticle)) // & 
                            ' from element '// trim(String(ElementID)) // ' lies '// &
                            'outside the mesh. '// &
                            'X, Y ='// trim(String(GlobPosArray_N(IParticle,1))) // trim(String(GlobPosArray_N(IParticle,2)))  )
            end if
            ! Determine new local particle position from global position
            call DetermineLocPosIterative(IParticle, &
                                 ElementID, GlopPos0, &
                                 Particles(IParticle)%LocPos, &
                                 NewElementID, LocPos, Success, &
                                 GlobPosArray_N)!, &
                                 !ShapeValuesArray)
        
            ! Correct if local particle position is outside element due to rounding
            if (LocPos(1)<0.0) then
                !LocPos(1) = 0.0
                !
                !if (CalParams%OutputDebugData) then
                !  call WriteInLogFile('Particle neg loc pos (x) '// trim(String(IParticle))    //' '// &
                !                                                trim(String(NewElementID)) //' '// &
                !                                                trim(String(LocPos(1))))
                !end if
            end if
            
            if (LocPos(2)<0.0) then
                !LocPos(2) = 0.0
                !
                !if (CalParams%OutputDebugData) then
                !  call WriteInLogFile('Particle neg loc pos (y) '// trim(String(IParticle))    //' '// &
                !                                                trim(String(NewElementID)) //' '// &
                !                                                trim(String(LocPos(2))))
                !end if
            end if
            
            if (NDIM == 3) then ! only for 3D case
              if (LocPos(3)<0.0) then
              !  LocPos(3) = 0.0
              !  
              !  if (CalParams%OutputDebugData) then
              !    call WriteInLogFile('Particle neg loc pos (z) '// trim(String(IParticle))    //' '// &
              !                                                  trim(String(NewElementID)) //' '// &
              !                                                  trim(String(LocPos(3))))
              !  end if
              end if  
            end if
 
           sumLocPos = 0.0
           do I = 1, NVECTOR
             sumLocPos = sumLocPos + LocPos(I)
           end do

           if (sumLocPos > 1.0) then
             diff = sumLocPos - 1.0
             do I = 1, NVECTOR
               LocPos(I) = LocPos(I) - diff/ sqrt(real(NVECTOR, REAL_TYPE))
             end do
           end if
            
            ! Update Particle%LocPos, Particle%ElementID and EleParticles array
            if (Success) then
              Particles(IParticle)%LocPos = LocPos ! LocPos --> we need to know which LocPos we are using
              call SetParticleElementID(IParticle, IParticle, NewElementID, ElementIDArray_NPlus1, EleParticles_NPlus1)

            else ! Error in determining the new location of the particle
              call GiveError('  Error in determining the local ' // &
                             'coordinate of particle '// trim(String(IParticle)) // &
                             ' in element '// trim(String(NewElementID)))
            end if

            ! Calculate new particle shape function values
            !call SetParticleShapeFunctionData(Particles(IParticle), IParticle, ShapeValuesArray, DShapeValuesArray)            
            call SetParticleShapeFunctionData(Particles(IParticle), IParticle, ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1) ! update shape functions


          end do

        end subroutine UpdateParticlePos

        subroutine CalcParticleDisplacements(ParticleIndex, ElementID, ShapeValuesArray, IncrementalDisplacementSoil, &
            UStepArray, UPhaseArray, UArray)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements
        !               for particle at ParticleIndex in Particles array
        !
        !     ParticleIndex : Index of the considered particle in the particle house-keeping arrays
        !     ElementID : ID of the element the particle is located in
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex, ElementID
          
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: IncrementalDisplacementSoil

          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement

          ! Shape function values input
          ! implicit scheme - shape values array inout   
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(inout) :: UArray
          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(inout) :: UPhaseArray
          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(inout) :: UStepArray
          
          !!get particle entity ID
            !if (CalParams%ApplyContactAlgorithm) then
            !  iEntity = EntityIDArray(ParticleIndex) 
            !else
              iEntity = 1
            !end if

            ParticleDisplacement = 0.0

            if ((NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid)) then
              if (.not.CalParams%SkipConvection) then
                do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
                  do INode = 1, ELEMENTNODES  ! Loop over all nodes in element
                    NodeID = iabs(ElementConnectivities(INode, ElementID)) ! this needs to be intent in if we are in multipatch implementation
                    DofID = ReducedDof(NodeID) + IDof
                      ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + &
                                                   ShapeValuesArray(ParticleIndex,INode) * &
                                                   IncrementalDisplacementSoil(DofID,iEntity)
                  enddo
                enddo
              else
                !do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
                !  do INode = 1, ELEMENTNODES  ! Loop over all nodes in element
                !    NodeID = iabs(ElementConnectivities(INode, ElementID))
                !    DofID = ReducedDof(NodeID) + IDof
                !    ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + &
                !                                 ShapeValuesArray(ParticleIndex,INode) * &
                !                                 AccumulatedIncDisplacementSoil(DofID,iEntity)
                !  enddo
                !enddo
              endif
            endif

            if ((MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid).and. &
                         (.not.(NFORMULATION==1))) then !Material Point LIQUID
              !do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              !  do INode = 1, ELEMENTNODES  ! Loop over all nodes in element
              !    NodeID = iabs(ElementConnectivities(INode, ElementID) ) 
              !    DofID = ReducedDof(NodeID) + IDof
              !
              !    ParticleDisplacement(IDof) = ParticleDisplacement(IDof) +  &
              !                                 ShapeValuesArray(ParticleIndex,INode) *  &
              !                                 IncrementalDisplacementWater(DofID,iEntity)
              !  end do
              !end do
           end if

            ! Incremental particle displacement
            UStepArray(ParticleIndex,:) = ParticleDisplacement

            ! Phase particle displacement
            UPhaseArray(ParticleIndex,:) = UPhaseArray(ParticleIndex,:) + ParticleDisplacement
            
            ! Total particle displacement
            UArray(ParticleIndex,:) =  UArray(ParticleIndex,:) + ParticleDisplacement
          
        end subroutine CalcParticleDisplacements
        
         subroutine DetermineLocPosIterative(ParticleID, OldElementID, GlobPos0, &
                       OldLocPos, NewElementID, NewLocPos, Success, &
                        GlobPosArray)!, &
                        !ShapeValuesArray)
        !**********************************************************************
        !
        !    Function:  Determine local coordinates and, eventually, the
        !               element, that belong to GlobPos.
        !
        !               Two recursive algorithms have been implemented:
        !               Version 1 checks whether GlobPos lies inside elements surrounding one of the
        !               corner nodes of NewElementID. If no element contains GlobPos the element
        !               whose centrepoint is closest to GlobPos is checked ... .
        !               Version 2 checks whether GlobPos lies inside NewElementID. If not, the element
        !               adjacent to the side, through which GlobPos has left NewElementID is checked ... .
        !               NOTE: Using Version 2 does not require to check whether GlobPos lies inside the
        !               original element of the considered particle!
        !               Version 2 seems to be faster as less elements are checked. Does it always work?
        !
        !     ParticleID : ID of the particle which moved to another element
        !     OldElementID : Current ID of the element the particle was located in
        !     OldLocPos : Current local position of the particle
        !
        ! O   NewElementID : ID of the element, the particle moved into
        ! O   NewLocPos : New local coordinates of the particle with GlobPos
        ! O   Success : True, if the local coordinate and element could be determined
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleID
          integer(INTEGER_TYPE), intent(in) :: OldElementID
          real(REAL_TYPE), dimension(:), intent(in) :: OldLocPos
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos0
          integer(INTEGER_TYPE), intent(out) :: NewElementID
          real(REAL_TYPE), dimension(NVECTOR), intent(out) :: NewLocPos
          logical, intent(out) :: Success
          ! Local variables
          logical :: OutsideElement, Fixed
          
          
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(inout) :: ShapeValuesArray 
          real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR), intent(in) :: GlobPosArray
          
          NewElementID = OldElementID
          NewLocPos = OldLocPos
          Success = .true.

          ! Try determining the local particle coordinates from the global ones
          call GetLocalCoordinates(GlobPosArray(ParticleID,:), &
                                   OldElementID, &
                                   ELEMENTNODES, Counters%NEl,  &
                                   Counters%NodTot, NDIM, &
                                   NodalCoordinates, &
                                   ElementConnectivities, &
                                   NewLocPos, &
                                   OutsideElement, Success)!, &
                                   !ShapeValuesArray)
          
          if (OutsideElement) then ! Particle moved to another element
           ! Determine new element and new local position of the particle
            call DetermineGlobPosElement(GlobPosArray(ParticleID,:), &
                                         OldElementID, &
                                         ELEMENTNODES,  &
                                         Counters%NEl,  &
                                         Counters%NodTot,  &
                                         NDIM, &
                                         NodalCoordinates,  &
                                         ElementConnectivities, &
                                         NewElementID, &
                                         NewLocPos, &
                                         Success)
            if (Success) then
              if (CalParams%OutputDebugData) then
              call WriteInLogFile('  Particle '   // trim(String(ParticleID))        // &
                                  ' moved from '  // trim(String(OldElementID))      // &
                                  ' to '          // trim(String(NewElementID))      // &
                                  ' at local pos:'// trim(String(NewLocPos, 1, 3))   // &
                                  'TimeStep '     // trim(String(Calparams%Timestep))  )
              end if
            else ! New element could not be found
              Fixed = .false.

              if (.not.Fixed) then
                call GiveError('  Failed finding new element for particle ' // trim(String(ParticleID)) // &
                               ' in element '  // trim(String(OldElementID)) // &
                               ' at position ' // trim(String(GlobPosArray(ParticleID,:),1,NVECTOR)) // &
                               ' in time step '// trim(String(CalParams%TimeStep)) )
              end if
            end if
          end if

        end subroutine DetermineLocPosIterative

        
       

        
        
        
        


        


        
        
         subroutine AssignStressesToParticles()
        !**********************************************************************
        !
        !    Function: Map data from first particle to other particles in the element using 
        !              interpolation functions.
        !              At the moment only stresses are mapped. This routine can be
        !              easily extended for arbitrary information passed on as
        !              additional parameters.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
 
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, IPart, FirstParticleIndex, NElemPart,  &
                     ParticleIndex, IAElement
          real(REAL_TYPE) :: FirstParticleWaterPressure
          real(REAL_TYPE) :: FirstParticleGasPressure
          real(REAL_TYPE) :: FirstDBulkViscousPressure
          real(REAL_TYPE), dimension(NTENSOR) :: FirstParticleStress
          
          FirstParticleWaterPressure = 0.0
          FirstParticleGasPressure = 0.0
          FirstParticleStress = 0.0
          FirstDBulkViscousPressure = 0.0
    
          do IAElement = 1, Counters%NAEl
            IElement = ActiveElement_N(IAElement)
            if ( (.not.IsParticleIntegration(IElement) ) ) then ! Loop over all active fully filled elements
     
             FirstParticleIndex = GetParticleIndex_N(1, IElement)
     
             FirstParticleStress = SigmaEffArray_N(FirstParticleIndex,:)
             FirstParticleWaterPressure = Particles(FirstParticleIndex)%WaterPressure
             FirstParticleGasPressure =  Particles(FirstParticleIndex)%GasPressure
             FirstDBulkViscousPressure = Particles(FirstParticleIndex)%DBulkViscousPressure
     
             NElemPart = NPartEle_N(IElement)  ! Number of particles in element
             
             do IPart = 2, NElemPart ! loop over particles. The stress of the first particle is already assigned.
               
               ParticleIndex = GetParticleIndex_N(IPart, IElement)
             
               SigmaEffArray_N(ParticleIndex,:) = FirstParticleStress
               Particles(ParticleIndex)%WaterPressure = FirstParticleWaterPressure
               Particles(ParticleIndex)%GasPressure = FirstParticleGasPressure
               Particles(ParticleIndex)%DBulkViscousPressure = FirstDBulkViscousPressure
             end do
                   
            end if
          end do
      
        end subroutine AssignStressesToParticles
        
          subroutine AssignStateParametersToParticles()
        !**********************************************************************
        !
        !    Function: Map the state variables from first particle to other particles in the element 
        !
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
 
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, IPart, FirstParticleIndex, NElemPart, &
                     ParticleIndex, IAElement
          real(REAL_TYPE) :: FirstParticleCohesionStSoft, &
                              FirstParticlePhiStSoft, &
                              FirstParticlePsiStSoft
          real(REAL_TYPE) ::  FirstParticleHPStateVariables (2), &
                               FirstParticleModifiedHPStateVariables(2), &
                               FirstParticleHPIGStateVariables (7), &
                               FirstParticleEpsP(NTENSOR), &
                               FirstParticleSigmaPrin(NTENSOR)
          real(REAL_TYPE) :: FirstParticlePP
          real(REAL_TYPE) :: FirstParticleESMstatev(NSTATEVAR)
     
         if (.not.IsMPMComputation()) RETURN ! FEM...no need ..already one particle
          
          FirstParticleHPStateVariables   = 0.0
          FirstParticleModifiedHPStateVariables = 0.0
          FirstParticleHPIGStateVariables = 0.0
          FirstParticleEpsP = 0.0
          FirstParticleSigmaPrin = 0.0
          FirstParticleCohesionStSoft = 0.0
          FirstParticlePhiStSoft = 0.0
          FirstParticlePsiStSoft = 0.0
          FirstParticlePP = 0.0
        
          do IAElement = 1, Counters%NAEl
            IElement = ActiveElement_N(IAElement)
            if ( (.not.IsParticleIntegration(IElement) ) ) then ! Loop over all active fully filled elements
     
             FirstParticleIndex = GetParticleIndex_N(1, IElement)
     

             FirstParticleHPStateVariables = GetHPStateVariables(Particles(FirstParticleIndex))
             FirstParticleHPIGStateVariables = GetHPIGStateVariables(Particles(FirstParticleIndex))
             FirstParticleModifiedHPStateVariables = GetModifiedHPStateVariables(Particles(FirstParticleIndex))

             FirstParticleEpsP = GetEpsP(Particles(FirstParticleIndex))
             FirstParticleSigmaPrin = GetSigmaPrin(Particles(FirstParticleIndex))
             FirstParticleCohesionStSoft = Particles(FirstParticleIndex)%CohesionStSoft
             FirstParticlePhiStSoft = Particles(FirstParticleIndex)%PhiStSoft
             FirstParticlePsiStSoft = Particles(FirstParticleIndex)%PsiStSoft
             FirstParticlePP = Particles(FirstParticleIndex)%PP
             
             !ESM
             FirstParticleESMstatev = ESMstatevArray(FirstParticleIndex,:)
     
             NElemPart = NPartEle_N(IElement)  ! Number of particles in element
             
             do IPart = 2, NElemPart ! loop ove particles. The state variables of the first particle are assigned.
               
               ParticleIndex = GetParticleIndex_N(IPart, IElement)
               
                
               call SetHPStateVariables(Particles(ParticleIndex), FirstParticleHPStateVariables)
                
               call SetHPIGStateVariables(Particles(ParticleIndex), FirstParticleHPIGStateVariables)
               call SetModifiedHPStateVariables(Particles(ParticleIndex), FirstParticleModifiedHPStateVariables)
            
               !MCStrainSoftening model
               call SetEpsP(Particles(ParticleIndex), FirstParticleEpsP)
               call SetSigmaPrin(Particles(ParticleIndex), FirstParticleSigmaPrin)
               Particles(ParticleIndex)%CohesionStSoft = FirstParticleCohesionStSoft
               Particles(ParticleIndex)%PhiStSoft = FirstParticlePhiStSoft
               Particles(ParticleIndex)%PsiStSoft = FirstParticlePsiStSoft
                
                ! MCC model    
                Particles(ParticleIndex)%PP = FirstParticlePP
                
                !ESM
                ESMstatevArray(ParticleIndex,:) = FirstParticleESMstatev
             end do
                   
            end if
          end do
      
         end subroutine AssignStateParametersToParticles
                
        

        subroutine SetInitialParticleStressForNextLoadStep(ActiveElement, NPartEle, SigmaEffArray, SigmaEff0Array) ! SigmaEffStep,
        !**********************************************************************
        !
        !    Function:  Stores particle stresses, water pressures and gas pressures
        !               determined in the Convective Phase
        !               as the initial stresses/pressures of the next load step.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IAElement, IElement, IParticle, ParticleIndex
          real(REAL_TYPE), dimension(NTENSOR) :: SigmaEffStep, SigmaEff !local varaibles

          
          ! implicit scheme variables vectorization 
          real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(in) :: SigmaEffArray !new stress computed in the convective phase
          real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(inout) :: SigmaEff0Array ! old stres from previous step
                    
          integer(INTEGER_TYPE), dimension(Counters%NAEl), intent(in) :: ActiveElement         
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          
           do IAElement = 1, Counters%NAEl
             IElement = ActiveElement(IAElement)
                do IParticle = 1, NPartEle(IElement) ! Loop over particles in IElement
                  ParticleIndex = GetParticleIndex_N(IParticle, IElement)
                  
                  SigmaEffStep = SigmaEffArray(ParticleIndex,:) - SigmaEff0Array(ParticleIndex,:)  ! local variable
                  call SetSigmaEffStep(Particles(ParticleIndex), SigmaEffStep) ! this just assigns the stresses
                  
                  SigmaEff = SigmaEffArray(ParticleIndex,:) ! local vaiable 
                  SigmaEff0Array(ParticleIndex,:) = SigmaEff !inout variable
                  
                  Particles(ParticleIndex)%WaterPressure0 = Particles(ParticleIndex)%WaterPressure
                  Particles(ParticleIndex)%GasPressure0 = Particles(ParticleIndex)%GasPressure
                  
                end do
          end do
     
        end subroutine SetInitialParticleStressForNextLoadStep

        subroutine MapDataFromNodesToParticles()
        !**********************************************************************
        !
        !    Function:  Updates the location of particles. Connected to this task are:
        !
        !               - calculate particle displacements from nodal displacements
        !
        !**********************************************************************
        
        implicit none        
          
          ! local variables        
          integer(INTEGER_TYPE) :: IParticle, ElementID
          
          if (.not.CalParams%ApplyFEMtoMPM) RETURN

          do IParticle = 1, Counters%NParticles ! loop over material points
            
            ElementID = ElementIDArray_N(IParticle)

            ! Calculate particle displacements from nodal displacements
            call MapTotalDisplacementsToParticles(IParticle, ElementID)

            ! Calculate particle velocities from nodal displacements
            call MapTotalVelocityToParticles(IParticle, ElementID)
          
          end do
          
        end subroutine MapDataFromNodesToParticles 


        subroutine MapTotalDisplacementsToParticles(ParticleIndex, &
                                                    ElementID)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements
        !               for particle at ParticleIndex in Particles array
        !
        !     ParticleIndex : Index of the considered particle in the particle house-keeping arrays
        !     ElementID : ID of the element the particle is located in
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex, ElementID

          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode
          integer(INTEGER_TYPE) :: NodeID, DofID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement

           ParticleDisplacement = 0.0

           do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
             do INode = 1, ELEMENTNODES  ! Loop over all nodes in element

               NodeID = iabs(ElementConnectivities(INode, ElementID) ) 
               DofID = ReducedDof(NodeID) + IDof
               ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray_N(ParticleIndex,INode) * TotalDisplacementSoil_N(DofID)
               
             end do
           end do

           ! Total particle displacement
           UArray(ParticleIndex,:) =  UArray(ParticleIndex,:) + ParticleDisplacement

        end subroutine MapTotalDisplacementsToParticles

        subroutine MapTotalVelocityToParticles(ParticleIndex, ElementID)
        !**********************************************************************
        !
        !    Function:  Calculate particle velocities from nodal displacements
        !               for particle at ParticleIndex in Particles array
        !
        !     ParticleIndex : Index of the considered particle in the particle house-keeping arrays
        !     ElementID : ID of the element the particle is located in
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex, ElementID
     
          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity

            !!get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(ParticleIndex) 
            else
              iEntity = 1
            end if 
              
            ParticleVelocity = 0.0
            
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                  
                NodeID = iabs(ElementConnectivities(INode, ElementID) ) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleVelocity(IDof) = ParticleVelocity(IDof) + ShapeValuesArray_N(ParticleIndex,INode) * TotalVelocitySoil_N(DofID,iEntity)
                
              end do                                             
            end do 
              
            ! Total particle velocity
             VelocityArray_N(ParticleIndex,:) = ParticleVelocity
          
        end subroutine MapTotalVelocityToParticles
                 
        
        !  subroutine Fix() ! 3D function
        !!**********************************************************************
        !!
        !!    Function:  Read some data related to pile driving
        !!    Note: 3D function  
        !!
        !!     FileName : Name of the PIL file to open
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  ! Local variables
        !  character(len = 255) :: CompleteFileName
        !  real(REAL_TYPE) :: x, y
        !  
        !  x = 0.0
        !  y = 0.0
        !  
        !  CompleteFileName = trim(CalParams%FileNames%ProjectName)//'.FIX'
        !
        !  if (FExist(CompleteFileName) ) then
        !  
        !   CalParams%ApplyFixParticlesK0 = .true.
        !
        !    call FileOpen(1, CompleteFileName)
        !
        !    read(1, *) x
        !    read(1, *) y
        !    
        !    CalParams%Fix (1) = x
        !    CalParams%Fix (2) = y
        !   
        !    
        !    close(1)
        !
        !  end if ! exists
        !
        !end subroutine Fix


        subroutine ResetDisplacements()
        !**********************************************************************
        !
        !    Function:  reset total displacements to zero
        !
        !*********************************************************************     

        implicit none
        
          if(.not.CalParams%ApplyResetDisplacements) RETURN
          
          TotalDisplacementSoil_N = 0.0
          PhaseDisplacementSoil = 0.0
          
          
          if (NFORMULATION==2) then !constituents
            TotalDisplacementWater = 0.0
            PhaseDisplacementWater = 0.0
          end if
          
        end subroutine ResetDisplacements


      end module ModConvectivePhase

      !subroutine UpdateParticleHouseKeepingMock(NPartEle, EleParticlesHelp, IsActiveElement)
      !!**********************************************************************
      !!
      !!    Function:  update housekeepting of material points
      !!
      !!*********************************************************************     
      !
      !!use ModConvectivePhase, only: UpdateParticleHouseKeeping
      !!use ModDynamicExplicit, only: UpdateParticleHouseKeeping
      !
      !implicit none
      !
      !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(inout) :: NPartEle
      !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(inout) :: EleParticlesHelp
      !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(inout) :: IsActiveElement
      !
      !!--HARDCODED removed this
      !call UpdateParticleHouseKeeping(NPartEle, EleParticlesHelp, IsActiveElement)
      !
      !end subroutine UpdateParticleHouseKeepingMock
      !
      !subroutine SetActiveElementMock()
      !!**********************************************************************
      !!
      !!    Function:  set active elements
      !!
      !!*********************************************************************  
      !use ModConvectivePhase, only: SetActiveElement
      !implicit none
      !call SetActiveElement()
      !end subroutine SetActiveElementMock