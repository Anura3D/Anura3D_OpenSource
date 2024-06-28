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
      !     $Revision: 10126 $
      !     $Date: 2024-04-08 12:14:40 +0200 (Mon, 08 Apr 2024) $
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

        MaterialIndex = MaterialIDArray(IParticle)
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

        MaterialIndex = MaterialIDArray(IParticle)
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
      

        subroutine UpdateParticleStrains()
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

          ! set the nodal displacement increments to the first entity
          ! only need to change if the contact model is used and the next particle belongs to a different entity
          iEntityDefault = 1 ! set default to first entity
          do I = 1,Counters%N
            DUTotEnt(I) = IncrementalDisplacementSoil(I,iEntityDefault)
          end do

          do IActiveElement = 1, Counters%NAEl ! loop over all active elements
            
              IElement = ActiveElement(IActiveElement) ! get element number of the active element
              NElemPart = NPartEle(IElement) ! get number of particles in the active element
              
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if  
              
              do IParticle = 1, NElemPart ! loop over all particles of the active element
              
                ParticleIndex = GetParticleIndex(IParticle, IElement) ! get global particle ID
                
                if( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then 

                  ! determine B matrix of deformed element at the particle local position
                  call BMatrix(Particles(ParticleIndex)%LocPos, &
                               ELEMENTNODES, Counters%NEl, &
                               Counters%NodTot, NDIM, &
                               IElement, ElementConnectivities, &
                               NodalCoordinatesUpd,  &
                               BMatrixDeformed, DetJac)
                  
                  ! get the nodal displacement increments for the specific entity
                  if (CalParams%ApplyContactAlgorithm) then ! CONTACT IS APPLIED
                    iEntity = EntityIDArray(ParticleIndex) ! get the entity ID for the current particle
                    if (iEntity /= iEntityDefault) then ! the particle belongs to a different entity from the previous entity
                      DUTotEnt(:) = IncrementalDisplacementSoil(1:Counters%N,iEntity)
                      iEntityDefault = iEntity !set the entity ID for the next particle loop
                    end if
                  end if

                  ! Eps in 3D (Exx, Eyy, Ezz, Gxy, Gyz, Gzx), in 2D (Exx, Eyy, Ezz, Gxy)
                  call Get_Strain(IElement, IParticle, &
                                ElementConnectivities, & 
                                BMatrixDeformed, DUTotEnt, ReducedDof, &
                                Eps) ! Eps = strain increment at material point
     
                  ! Update particle strains
                  call SetEpsStep(Particles(ParticleIndex), Eps)
                  if (.not.CalParams%ApplyImplicitQuasiStatic) then
                    call IncreaseEps(Particles(ParticleIndex), Eps)
                  end if

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
            if (IsActiveElement(IElement)) then  !active elements only
              
                NElemPart = NPartEle(IElement)    !number of particles in element
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if
              
              do IParticle = 1, NElemPart       !loop over all particles of the element
                ParticleIndex = GetParticleIndex(IParticle, IElement)   !particle global ID
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


        subroutine UpdateParticlePos()
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
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: I, IParticle, ElementID, NewElementID
          real(REAL_TYPE), dimension(NVECTOR) :: LocPos, GlopPos0
          logical :: Success, IsInside
          real(REAL_TYPE) :: diff, sumLocPos

           do IParticle = 1, Counters%NParticles ! Loop over particles

            ElementID = ElementIDArray(IParticle)

            ! Calculate particle displacements from nodal displacements
            call CalcParticleDisplacements(IParticle, ElementID)

            GlopPos0 =  GlobPosArray(IParticle,:)  ! store the initial global position of the particle

            ! Calculate new particle global positions
            GlobPosArray(IParticle,:) = GlobPosArray(IParticle,:) + UStepArray(IParticle,:)

           if (CalParams%ApplyFixParticlesK0) then ! only 3D functionality
             if (GlobPosArray(IParticle,1) >= (CalParams%Fix(1))) then ! allow sliding along y-axis only
                 GlobPosArray(IParticle,1) = GlopPos0(1)
                 GlobPosArray(IParticle,3) = GlopPos0(3)
             end if
             if (GlobPosArray(IParticle,2) <= (CalParams%Fix(2))) then
                GlobPosArray(IParticle,:) = GlopPos0 ! particles along the bottom are fixed
             end if
           end if ! ApplyFixParticlesK0

          IsInside = IsInBoundingBox(GlobPosArray(IParticle,:))
          ! Check whether particle moved outside the mesh
          if (.not.IsInside ) then
             if (NDIM==3) call GiveError('  Particle '// trim(String(IParticle)) // & 
                            ' from element '// trim(String(ElementID)) // ' lies '// &
                            'outside the mesh. '// &
                            'X, Y, Z ='// trim(String(GlobPosArray(IParticle,1))) // trim(String(GlobPosArray(IParticle,2))) // trim(String(GlobPosArray(IParticle,3))) )
             if (NDIM==2) call GiveError('  Particle '// trim(String(IParticle)) // & 
                            ' from element '// trim(String(ElementID)) // ' lies '// &
                            'outside the mesh. '// &
                            'X, Y ='// trim(String(GlobPosArray(IParticle,1))) // trim(String(GlobPosArray(IParticle,2)))  )
            end if
            ! Determine new local particle position from global position
            call DetermineLocPosIterative(IParticle, &
                                 ElementID, GlopPos0, &
                                 Particles(IParticle)%LocPos, &
                                 NewElementID, LocPos, Success)
        
            ! Correct if local particle position is outside element due to rounding
            if (LocPos(1)<0.0) then
                LocPos(1) = 0.0
                
                if (CalParams%OutputDebugData) then
                  call WriteInLogFile('Particle neg loc pos (x) '// trim(String(IParticle))    //' '// &
                                                                trim(String(NewElementID)) //' '// &
                                                                trim(String(LocPos(1))))
                end if
            end if
            
            if (LocPos(2)<0.0) then
                LocPos(2) = 0.0
                
                if (CalParams%OutputDebugData) then
                  call WriteInLogFile('Particle neg loc pos (y) '// trim(String(IParticle))    //' '// &
                                                                trim(String(NewElementID)) //' '// &
                                                                trim(String(LocPos(2))))
                end if
            end if
            
            if (NDIM == 3) then ! only for 3D case
              if (LocPos(3)<0.0) then
                LocPos(3) = 0.0
                
                if (CalParams%OutputDebugData) then
                  call WriteInLogFile('Particle neg loc pos (z) '// trim(String(IParticle))    //' '// &
                                                                trim(String(NewElementID)) //' '// &
                                                                trim(String(LocPos(3))))
                end if
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
              Particles(IParticle)%LocPos = LocPos
              call SetParticleElementID(IParticle, IParticle, NewElementID)

            else ! Error in determining the new location of the particle
              call GiveError('  Error in determining the local ' // &
                             'coordinate of particle '// trim(String(IParticle)) // &
                             ' in element '// trim(String(NewElementID)))
            end if

            ! Calculate new particle shape function values
            call SetParticleShapeFunctionData(Particles(IParticle), IParticle)

          end do

        end subroutine UpdateParticlePos

        subroutine CalcParticleDisplacements(ParticleIndex, ElementID)
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
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement

            !!get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(ParticleIndex) 
            else
              iEntity = 1
            end if

            ParticleDisplacement = 0.0

            if ((NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid)) then
              if (.not.CalParams%SkipConvection) then
                do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
                  do INode = 1, ELEMENTNODES  ! Loop over all nodes in element
                    NodeID = iabs(ElementConnectivities(INode, ElementID))
                    DofID = ReducedDof(NodeID) + IDof
                      ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + &
                                                   ShapeValuesArray(ParticleIndex,INode) * &
                                                   IncrementalDisplacementSoil(DofID,iEntity)
                  enddo
                enddo
              else
                do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
                  do INode = 1, ELEMENTNODES  ! Loop over all nodes in element
                    NodeID = iabs(ElementConnectivities(INode, ElementID))
                    DofID = ReducedDof(NodeID) + IDof
                    ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + &
                                                 ShapeValuesArray(ParticleIndex,INode) * &
                                                 AccumulatedIncDisplacementSoil(DofID,iEntity)
                  enddo
                enddo
              endif
            endif

            if ((MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid).and. &
                         (.not.(NFORMULATION==1))) then !Material Point LIQUID
              do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
                do INode = 1, ELEMENTNODES  ! Loop over all nodes in element
                  NodeID = iabs(ElementConnectivities(INode, ElementID) ) 
                  DofID = ReducedDof(NodeID) + IDof

                  ParticleDisplacement(IDof) = ParticleDisplacement(IDof) +  &
                                               ShapeValuesArray(ParticleIndex,INode) *  &
                                               IncrementalDisplacementWater(DofID,iEntity)
                end do
              end do
           end if

            ! Incremental particle displacement
            UStepArray(ParticleIndex,:) = ParticleDisplacement

            ! Phase particle displacement
            UPhaseArray(ParticleIndex,:) = UPhaseArray(ParticleIndex,:) + ParticleDisplacement
            
            ! Total particle displacement
            UArray(ParticleIndex,:) =  UArray(ParticleIndex,:) + ParticleDisplacement
          
        end subroutine CalcParticleDisplacements
        
         subroutine DetermineLocPosIterative(ParticleID, OldElementID, GlobPos0, &
                       OldLocPos, NewElementID, NewLocPos, Success)
        !**********************************************************************
        !
        !    Function:  Determine local coordinates and, eventually, the
        !               element in which particle is located.
        !
        !               This subroutine checks if the MP global coordinates can 
        !               be expressed as local coordinates of the element in which 
        !               the MP was in the previous time step. If it fails, it has 
        !               to find the new element for the MP.
        !               The search consists on an exhaustive algorithm around adjacent elements. If the 
        !               local coordinates can be determined in one of the adjacent 
        !               elements the search ends. 
        !               If the MP is not in the adjacent elements, a more time expensive 
        !               combinatorial algorithm capable of determining the location of the MP is called
        !               (subroutine DetermineLocPosIterative2).
        !
        ! I   ParticleID : ID of the particle which moved to another element
        ! I   OldElementID : Current ID of the element the particle was located in
        ! I   OldLocPos : Current local position of the particle
        !
        ! O   NewElementID : ID of the element, the particle moved into
        ! O   NewLocPos : New local coordinates of the particle with GlobPos
        ! O   Success : false, if the local coordinate and element could not be determined 
        !               in a maximum number of tries
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
          logical :: Fixed
          logical :: OutsideElement !true if LocPos can not be found in element 
          integer(INTEGER_TYPE) :: NAdjEl, I
          
          NewElementID = OldElementID
          NewLocPos = OldLocPos
          Success = .true.
          OutsideElement = .true.

          ! Try determining the local particle coordinates from the global ones in the previous time step element
          call GetLocalCoordinates(GlobPosArray(ParticleID,:), &
                                   OldElementID, &
                                   ELEMENTNODES, Counters%NEl,  &
                                   Counters%NodTot, NDIM, &
                                   NodalCoordinates, &
                                   ElementConnectivities, &
                                   NewLocPos, &
                                   OutsideElement, Success)
           if (OutsideElement) then ! Particle moved to another element
               NAdjEl = GetNElmOfElm(OldElementID) ! Get adjacent elements to the previous element
               
               do I = 1, NAdjEl      ! Try to determine the local position in adjacent elements
                   
               NewElementID = GetElmIOfElm(OldElementID, I)
               
                call GetLocalCoordinates(GlobPosArray(ParticleID,:), &
                                   NewElementID, &
                                   ELEMENTNODES, Counters%NEl,  &
                                   Counters%NodTot, NDIM, &
                                   NodalCoordinates, &
                                   ElementConnectivities, &
                                   NewLocPos, &
                                   OutsideElement, Success)
          
                    if (.not.OutsideElement) then       
                        if (CalParams%OutputDebugData) then ! If new element is found write debug data if required
                                                call WriteInLogFile('  Particle '   // trim(String(ParticleID))        // &
                                            ' moved from '  // trim(String(OldElementID))      // &
                                              ' to '          // trim(String(NewElementID))      // &
                                              ' at local pos:'// trim(String(NewLocPos, 1, 3))   // &
                                             'TimeStep '     // trim(String(Calparams%Timestep))  )
                        end if
                        exit ! Exit the adjacent elements loop if new element is found
                    end if
               end do
               !If new local position is not found in adjacent elements, try to find the new element with a combinatorial search algorithm
               if (OutsideElement) then     
                   call DetermineLocPosIterative2(ParticleID, OldElementID, GlobPos0, &
                       OldLocPos, NewElementID, NewLocPos, Success)
              end if
            end if

     end subroutine DetermineLocPosIterative
                       
                 subroutine DetermineLocPosIterative2(ParticleID, OldElementID, GlobPos0, &
                       OldLocPos, NewElementID, NewLocPos, Success)
        !**********************************************************************
        !
        !    Function:  Determine local coordinates and, eventually, the
        !               element, that belong to GlobPos.
        !
        !               Checks whether GlobPos lies inside NewElementID. If not, the element
        !               adjacent to the side, through which GlobPos has left NewElementID is checked.
        !
        ! I   ParticleID : ID of the particle which moved to another element
        ! I   OldElementID : Current ID of the element the particle was located in
        ! I   OldLocPos : Current local position of the particle
        ! I
        ! O   NewElementID : ID of the element, the particle moved into
        ! O   NewLocPos : New local coordinates of the particle with GlobPos
        ! O   Success : false, if the local coordinate and element could not be determined 
        !                in a maximum number of tries
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
                                   OutsideElement, Success)
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

        end subroutine DetermineLocPosIterative2


        subroutine UpdateParticleHouseKeeping()
        !**********************************************************************
        !
        !    Function:  Update the particle house-keeping lists (particle-element assignment)
        !               and element switches (activation of newly filled and deactivation of 
        !               empty element).
        !               NOTE: EleParticles and Particle%ElementID are already updated!!
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, IElement, FirstParticle

          NPartEle = 0
          call IHpSortKind8(EleParticles, Counters%NParticles)
          ! Count the particles in element, update NPartEle
          do IParticle = 1, Counters%NParticles ! Loop over all particles
            IElement = GetElementIDFromList(IParticle) ! Determine the element that IParticle is located in
            NPartEle(IElement) = NPartEle(IElement) + 1 ! Increase the particle counter for IElement
          end do


          ! Update the helper array and element switches
          FirstParticle = 1
          counters%NAEl = 0
          do IElement = 1, Counters%NEl ! Loop over all elements
            if (NPartEle(IElement)==0) then ! Empty element
              EleParticlesHelp(IElement) = -1

              IsActiveElement(IElement) = .false. ! Switch off the element
            else  ! Element contains particles
              counters%NAEl = counters%NAEl + 1
              EleParticlesHelp(IElement) = FirstParticle
              FirstParticle = FirstParticle + NPartEle(IElement)

              IsActiveElement(IElement) = .true. ! Switch on the element
            end if
          end do

        end subroutine UpdateParticleHouseKeeping
       

        subroutine SetUpEntityElements()
        !**********************************************************************
        !
        !  Function : For each entity a list of the elements containing particles belonging to this entity
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: IAEl, IEl, NElemPart, IPart, ParticleIndex, IEntity
         
          EntityElements = 0      !reset

          do IAEl = 1,Counters%NAel                                      !loop through all elements
            IEl = ActiveElement(IAEl)
            NElemPart = NPartEle(IEl)
            do IPart= 1, NElemPart                                 !loop through particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl)         !Get the particle ID
              IEntity = EntityIDArray(ParticleIndex)               !entity to which particle belong
              EntityElements(IEntity,IEl) = 1
            end do
          end do
         
        end subroutine SetUpEntityElements
        
        
        subroutine SetUpMaterialElements()
        !**********************************************************************
        !
        !  Function : For each material ID a list of the elements containing particles belonging to this material ID
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: IAEl, IEl, NElemPart, IPart, ParticleIndex, MaterialID
         
          MaterialElements = 0      !reset
          
                do IAEl = 1,Counters%NAel                                      !loop through all elements
            IEl = ActiveElement(IAEl)
            NElemPart = NPartEle(IEl)
            do IPart= 1, NElemPart                                 !loop through particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl)        ! Get the particle ID
              MaterialID = MaterialIDArray(ParticleIndex)       !material ID to which particle belong
              MaterialElements(MaterialID,IEl) = 1
            end do
          end do
                  
        end subroutine SetUpMaterialElements


        subroutine SetActiveElement()
        !**********************************************************************
        !
        !  Function : To set the active elements
        !
        !  ActiveElement: Active element vector
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE) :: IAEl, IEl, IError
        
          if (allocated(ActiveElement) ) then
            deallocate(ActiveElement, stat = IError)
          end if
          
          allocate(ActiveElement(counters%NAEl), stat = IError)
          ActiveElement = 0

          IAEl = 0
          do IEl=1, counters%NEl
            if (IsActiveElement(IEl)) then ! Active elements
              IAEl = IAEl + 1
              ActiveElement(IAEl) = IEl
            end if
          end do
 
          call DetermineActiveNodes()
          
        end subroutine SetActiveElement


        subroutine CheckFillingOfElements()
        !**********************************************************************
        !
        !    Function:  Determines whether an element is fully or partially filled.
        !               The global volume inside an element occupied by particles is
        !               calculated. If it lies below a threshold value and a
        !               boundary particle lies inside the element, the
        !               element is partially filled, else fully filled.
        !               The threshold value is defined as a percentage of the global element volume.
        !               Partially filled elements will be integrated by using mass points,
        !               fully filled elements will be integrated by using Gauss Point integration.
        !
        !               Also elements containing particles of different materials will be
        !               considered partially filled.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IElement, IAElement, IParticle, ParticleIndex, &
                                                           CheckMaterialID
          logical :: ContainsBoundaryParticle, IsMixedElement
          real(REAL_TYPE) :: SummedParticleWeights, &
                              RequiredFilledVolume

          if (.not.IsMPMWithMixedIntegration()) RETURN
            IsParticleIntegration = .false.

            do IAElement = 1, Counters%NAEl
             IElement = ActiveElement(IAElement)
              ContainsBoundaryParticle = .false.
              SummedParticleWeights = 0.0

              ParticleIndex = GetParticleIndex(1, IElement)
              CheckMaterialID = MaterialIDArray(ParticleIndex)
              IsMixedElement = .false.

              do IParticle = 1, NPartEle(IElement) ! Loop over particles in IElement
                ParticleIndex = GetParticleIndex(IParticle, IElement)
                ! Sum up volumes occupied by particles
                SummedParticleWeights = SummedParticleWeights + &
                  Particles(ParticleIndex)%IntegrationWeight
                ! Check for boundary particle inside the element
                if (Particles(ParticleIndex)%IsBoundaryParticle) then
                  ContainsBoundaryParticle = .true.
                end if

                if (MaterialIDArray(ParticleIndex)/= &
                    CheckMaterialID) then
                  IsMixedElement = .true.
                end if
              end do

              if (ContainsBoundaryParticle) then
                ! Only elements which contain at least one boundary particle can be considered as partially filled.
                ! If less than CalParams%RequiredDegreeOfFilling percent of the element volume is filled,
                ! the element is partially filled
                RequiredFilledVolume = ElementSpace(IElement) * CalParams%RequiredDegreeOfFilling
                if (SummedParticleWeights<=RequiredFilledVolume) then
                  IsParticleIntegration(IElement) = .true. ! Consider the element as partially filled
                end if
              end if

              if (IsMixedElement) then ! Elements with particles of different material
                IsParticleIntegration(IElement) = .true.
              end if
          end do

        end subroutine CheckFillingOfElements
        
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
            IElement = ActiveElement(IAElement)
            if ( (.not.IsParticleIntegration(IElement) ) ) then ! Loop over all active fully filled elements
     
             FirstParticleIndex = GetParticleIndex(1, IElement)
     
             FirstParticleStress = SigmaEffArray(FirstParticleIndex,:)
             FirstParticleWaterPressure = Particles(FirstParticleIndex)%WaterPressure
             FirstParticleGasPressure =  Particles(FirstParticleIndex)%GasPressure
             FirstDBulkViscousPressure = Particles(FirstParticleIndex)%DBulkViscousPressure
     
             NElemPart = NPartEle(IElement)  ! Number of particles in element
             
             do IPart = 2, NElemPart ! loop over particles. The stress of the first particle is already assigned.
               
               ParticleIndex = GetParticleIndex(IPart, IElement)
             
               SigmaEffArray(ParticleIndex,:) = FirstParticleStress
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
            IElement = ActiveElement(IAElement)
            if ( (.not.IsParticleIntegration(IElement) ) ) then ! Loop over all active fully filled elements
     
             FirstParticleIndex = GetParticleIndex(1, IElement)
     

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
     
             NElemPart = NPartEle(IElement)  ! Number of particles in element
             
             do IPart = 2, NElemPart ! loop ove particles. The state variables of the first particle are assigned.
               
               ParticleIndex = GetParticleIndex(IPart, IElement)
               
                
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
                
        subroutine SetInitialStressForNextLoadStep()
        !**********************************************************************
        !
        !    Function:  Call the routines in which the initial stresses of particles
        !               and Gauss points are stored
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Store intial particle stresses of next step
          call SetInitialParticleStressForNextLoadStep()

        end subroutine SetInitialStressForNextLoadStep

        subroutine SetInitialParticleStressForNextLoadStep()
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
          real(REAL_TYPE), dimension(NTENSOR) :: SigmaEffStep, SigmaEff

           do IAElement = 1, Counters%NAEl
             IElement = ActiveElement(IAElement)
                do IParticle = 1, NPartEle(IElement) ! Loop over particles in IElement
                  ParticleIndex = GetParticleIndex(IParticle, IElement)
                  
                  SigmaEffStep = SigmaEffArray(ParticleIndex,:) - SigmaEff0Array(ParticleIndex,:) 
                  call SetSigmaEffStep(Particles(ParticleIndex), SigmaEffStep)
                  
                  SigmaEff = SigmaEffArray(ParticleIndex,:)
                  SigmaEff0Array(ParticleIndex,:) = SigmaEff
                  
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
            
            ElementID = ElementIDArray(IParticle)

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
               ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(ParticleIndex,INode) * TotalDisplacementSoil(DofID)
               
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
                ParticleVelocity(IDof) = ParticleVelocity(IDof) + ShapeValuesArray(ParticleIndex,INode) * TotalVelocitySoil(DofID,iEntity)
                
              end do                                             
            end do 
              
            ! Total particle velocity
             VelocityArray(ParticleIndex,:) = ParticleVelocity
          
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
          
          TotalDisplacementSoil = 0.0
          PhaseDisplacementSoil = 0.0
          
          
          if (NFORMULATION==2) then !constituents
            TotalDisplacementWater = 0.0
            PhaseDisplacementWater = 0.0
          end if
          
        end subroutine ResetDisplacements


      end module ModConvectivePhase

      subroutine UpdateParticleHouseKeepingMock()
      !**********************************************************************
      !
      !    Function:  update housekeepting of material points
      !
      !*********************************************************************     
      use ModConvectivePhase, only: UpdateParticleHouseKeeping
      implicit none
      call UpdateParticleHouseKeeping()
      end subroutine UpdateParticleHouseKeepingMock
      
      subroutine SetActiveElementMock()
      !**********************************************************************
      !
      !    Function:  set active elements
      !
      !*********************************************************************  
      use ModConvectivePhase, only: SetActiveElement
      implicit none
      call SetActiveElement()
      end subroutine SetActiveElementMock