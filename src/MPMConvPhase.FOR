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
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr. 2022) $
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
      !use ModDYNConvectivePhase
      
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
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch_Temporary = 1

        WaveSpeed = 0.0
        consolidationcoeff = 0.0
        WaveSpeed_dyn = 0.0
        WaveSpeed_cons = 0.0

        Sr = 1.0
        dSrdp(1) = 0.0

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
          
          if (SoilModel == ESM_HYPOPLASTICITY_SAND .or. &
              SoilModel==ESM_EXTERNAL_SOIL_MODEL .or. & ! user-defined soil model
          SoilModel == ESM_ARB_Model) then 
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
          
          dx = ElementLMin(IEl,IPatch_Temporary)

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
          
          real(REAL_TYPE), dimension(NDIM,NDIM) :: RelativeDeformationGradient
          real(REAL_TYPE), dimension(NDIM,NDIM) :: RelativeDeformationGradientInverse
          
          real(REAL_TYPE), dimension(NDIM,NDIM) :: DeformationGradientPreviousInverse
          
          real(REAL_TYPE) :: RelativeDeformationGradientJacobian
          real(REAL_TYPE) :: RelativeDeformationGradientInverseJacobian
          
          
          real(REAL_TYPE) :: DeformationGradientPreviousJacobian
          real(REAL_TYPE) :: DeformationGradientPreviousInverseJacobian
          
          real(REAL_TYPE), dimension(NDIM, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac
          integer(INTEGER_TYPE) :: IElement, NElemPart, IParticle, ParticleIndex, I, iEntity, iEntityDefault, IActiveElement

          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          
          ! Volumetric locking variables 
          integer(INTEGER_TYPE) :: IDim
          real(REAL_TYPE), dimension(Maximum_NControlPoints,Counters%nEntity) :: &
              WeightedUpdatedDeformationGradientJacobian
          real(REAL_TYPE), dimension(Maximum_NControlPoints,Counters%nEntity) :: &
              UpdatedDeformationGradientJacobian
          real(REAL_TYPE), dimension(NDIM,NDIM) :: TensorialStrainIncrement
          
          ! set the nodal displacement increments to the first entity
          ! only need to change if the contact model is used and the next particle belongs to a different entity
          iEntityDefault = 1 ! set default to first entity
          do I = 1,Counters%N
            DUTotEnt(I) = IncrementalDisplacementSoil(I,iEntityDefault)
          end do

          
          !do IActiveElement = 1, Counters%NAEl ! loop over all active elements
           
          do IPatch = 1, Counters%NPatches ! loop over all patches
          do IActiveElement = 1, NAEl_NURBS(IPatch) !Counters%NAEl ! loop over all elements

              IElement = ActiveElement(IActiveElement, IPatch) ! get element number of the active element
              NElemPart = NPartEle(IElement,IPatch) ! get number of particles in the active element
              
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement, IPatch) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if  
              
              do IParticle = 1, NElemPart ! loop over all particles of the active element
              
                ParticleIndex = GetParticleIndex(IParticle, IElement, IPatch) ! get global particle ID
                
                if( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then 
                    
                  ! determine B matrix of deformed element at the particle local position
                  call BMatrix(Particles(ParticleIndex)%LocPos, &
                               ELEMENTNODES, nel_NURBS(IPatch), &!Counters%NEl, & !nel_NURBS(IPatch), &
                               NControlPoints(IPatch), NDIM, & !Counters%NodTot
                               IElement, ElementConnectivities(:,:,IPatch), &
                               NodalCoordinatesUpd(:,:,IPatch),  &
                               BMatrixDeformed, &
                               IPatch)!, DetJac) !update Jacob -> DetJac = | dx/dtildexi |
                  
                  ! get the nodal displacement increments for the specific entity
                  if (CalParams%ApplyContactAlgorithm) then ! CONTACT IS APPLIED
                    iEntity = EntityIDArray(ParticleIndex) ! get the entity ID for the current particle
                    if (iEntity /= iEntityDefault) then ! the particle belongs to a different entity from the previous entity
                      DUTotEnt(:) = IncrementalDisplacementSoil(1:Counters%N,iEntity)
                      iEntityDefault = iEntity !set the entity ID for the next particle loop
                    end if
                  end if

                  if (.not.CalParams%ApplyFBar) then 
                      
                      
                  ! Eps in 3D (Exx, Eyy, Ezz, Gxy, Gyz, Gzx), in 2D (Exx, Eyy, Ezz, Gxy)
                  call Get_Strain(IElement, ParticleIndex, &!IParticle, &
                                ElementConnectivities(:,:,IPatch), & 
                                BMatrixDeformed, DUTotEnt, ReducedDof, &
                                Eps, &
                                IPatch) ! Eps = strain increment at material point
                  
                  
                  
                  else
                  
                      ! Relative deformation gradient 
                      call Get_RelativeDeformationGradient(IElement, ParticleIndex, &
                                ElementConnectivities(:,:,IPatch), & 
                                BMatrixDeformed, DUTotEnt, ReducedDof, &
                                RelativeDeformationGradient, &
                                IPatch) ! RelativeDeformationGradient = relative deformation gradient at material point
                  
                      ! Relative deformation gradient jacobian
                      !call RJacInv_VolLockSmooth(NDIM, RelativeDeformationGradient, &!RelativeDeformationGradientInverse, &
                      !RelativeDeformationGradientJacobian)!, RelativeDeformationGradientInverseJacobian)
                      
                      call RJacInv(NDIM, RelativeDeformationGradient, RelativeDeformationGradientInverse, &
                      RelativeDeformationGradientJacobian, RelativeDeformationGradientInverseJacobian)
                        
                      ! Previous Relative deformation gradient jacobian  
                      !call RJacInv_VolLockSmooth(NDIM, Particles(ParticleIndex)%DeformationGradientPrevious, &!DeformationGradientPreviousInverse, &
                      !DeformationGradientPreviousJacobian)!, DeformationGradientPreviousInverseJacobian)
                  
                      call RJacInv(NDIM, Particles(ParticleIndex)%DeformationGradientPrevious, DeformationGradientPreviousInverse, &
                      DeformationGradientPreviousJacobian, DeformationGradientPreviousInverseJacobian)
                      
                      ! Set Relative deformation gradient in particles type
                      call SetRelDefGradient(Particles(ParticleIndex), RelativeDeformationGradient) 
                        
                      ! Set Relative deformation gradient jacobian in particles type
                      call SetRelDefGradientJacobian(Particles(ParticleIndex), RelativeDeformationGradientJacobian) 
                  
                      ! Set Previous  deformation gradient in particles type
                      call SetDefGradientJacobianPrevious(Particles(ParticleIndex), DeformationGradientPreviousJacobian) 
                  
                      
                  end if
                  
                  
                  ! Update particle strains
                  !if (.not.CalParams%ApplyBBar) then 
                  if (.not.CalParams%ApplyFBar) then 
                  call SetEpsStep(Particles(ParticleIndex), Eps) 
                  end if 
                  
                  ! if Applying BBar then add in the total strain later on
                  if (.not.CalParams%ApplyImplicitQuasiStatic .and. .not.CalParams%ApplyFBar) then !.and. .not.CalParams%ApplyBBar
                    call IncreaseEps(Particles(ParticleIndex), Eps)
                  end if

                end if ! Material Point can be MIXTURE or SOLID

              end do ! particle loop

          end do ! active element loop
          
          end do ! patch loop
          
          
          
 
          if (CalParams%ApplyFBar) then 
          
              call MapWeightedUpdatedDeformationGradientJacobianToNodes(WeightedUpdatedDeformationGradientJacobian) !CalculateUpdatedNodalDeformationGradientJacobian()
          
          
              call NormalizeWeightedUpdatedDeformationGradientJacobianAtNodes(UpdatedDeformationGradientJacobian, &
                                                                            WeightedUpdatedDeformationGradientJacobian)
                      
              call MapUpdatedDeformationGradientJacobianAtNodesToMPs(UpdatedDeformationGradientJacobian)
                      
          
              do IParticle = 1, Counters%NParticles
                  ! Update deformation gradients
                  call UpdateRelativeDeformationGradient(Particles(IParticle))
                  call UpdateDeformationGradientPrevious(Particles(IParticle)) 
                  
                  ! Update strains 
                  call CalculateTensorialInfinitesimalStrain(Particles(IParticle), TensorialStrainIncrement)
                  call CalculateEpsVoigt(Particles(IParticle), TensorialStrainIncrement, Eps)
                  call SetEpsStep(Particles(IParticle), Eps) 
                  call IncreaseEps(Particles(IParticle), Eps)
              end do 
          
          end if 
          
        end subroutine UpdateParticleStrains
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
         
       !-----------------------------------------------------------------------
       
       
       
       
       subroutine MapUpdatedDeformationGradientJacobianAtNodesToMPs(UpdatedDeformationGradientJacobian)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        ! hardcoded on the first patch
        real(REAL_TYPE), dimension(Maximum_NControlPoints,Counters%nEntity), intent(in) :: UpdatedDeformationGradientJacobian       !CC added nEntity
          
          
        ! Local variables
        integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
        integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: Nodes_Element
          
        real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity) :: NodUpdatedDefGradientJacobian !NodStrain NTENSOR, 
        real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch
       
        !!CC - changed to function call - needs to be set to zero here
        NoEn = Counters%nEntity
       
          
        do IPatch = 1, Counters%NPatches ! Loop over patches
          
            do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
            
                IEl = ActiveElement(IAEl, IPatch)
              
                Nodes_Element(1:ELEMENTNODES) = &
                        IEN(1:ELEMENTNODES, IEl, IPatch)
                
                !do I = 1, NTENSOR
              
                    
                        !ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
                        !IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       
                    ! start here 
                      NodUpdatedDefGradientJacobian(1:ELEMENTNODES, 1:Counters%nEntity) = &
                          UpdatedDeformationGradientJacobian(Nodes_Element(1:ELEMENTNODES), 1:NoEn) 
            
                !end do
           
                  
                do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                    if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 
                        !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
                        Particles(ParticleIndex)%DeformationGradientJacobian = 0.0
                        ParticleShape = ShapeValuesArray(ParticleIndex,:)
                 
                        if (CalParams%ApplyContactAlgorithm) then
                            iEntity = EntityIDArray(ParticleIndex)
                        else
                            iEntity = 1  
                        end if 
                 
                        do INode = 1, ELEMENTNODES
       
                            !do I = 1, NTENSOR!NVECTOR  
       
                                ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
                                ! Particle strain 
                                
                            !Particles(ParticleIndex)%EpsStep(I) = Particles(ParticleIndex)%EpsStep(I) + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
                            Particles(ParticleIndex)%DeformationGradientJacobian = Particles(ParticleIndex)%DeformationGradientJacobian + &
                                ( ParticleShape(INode) * NodUpdatedDefGradientJacobian(INode, iEntity) )
                            !end do    
                        
                        end do !Loop over nodes
                      
                    end if ! NumbOfLayers = 1 or SOLID MatPoint
                  
                    
                    !call IncreaseEps(Particles(ParticleIndex), Particles(ParticleIndex)%EpsStep)
                    
                end do !Loop over particles
              
            end do !elements    
          
        end do ! patches
        
       
       end subroutine MapUpdatedDeformationGradientJacobianAtNodesToMPs !MapVolumetricStrainToMPs_VolLockSmooth  
       
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        subroutine UpdateRelativeDeformationGradient(Particle) 
        !----------------------------------------------------------------
        ! To update the relative deformation gradient after smoothing the jacobians
        !----------------------------------------------------------------
        implicit none 
        
        type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
        
        integer(INTEGER_TYPE) :: IParticle
        
        real(REAL_TYPE) :: Exponent 
        
        Exponent = real(1)/real(NDIM)
        
        ! Calculate an assumed relative deformation gradient
        Particle%RelativeDeformationGradient = &
              (( Particle%DeformationGradientJacobian / (Particle%RelativeDeformationGradientJacobian * Particle%DeformationGradientPreviousJacobian) )**Exponent) &
              * Particle%RelativeDeformationGradient
        
        ! Set previous deformation gradient
        !Particle%DeformationGradientPrevious = matmul(Particle%RelativeDeformationGradient,&
        !                                    Particle%DeformationGradientPrevious)
        
        end subroutine UpdateRelativeDeformationGradient
        
       
       
       
       
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        
        
        
        subroutine UpdateDeformationGradientPrevious(Particle) 
        !----------------------------------------------------------------
        ! To update the relative deformation gradient after smoothing the jacobians
        !----------------------------------------------------------------
        implicit none 
        
        type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle

        integer(INTEGER_TYPE) :: IParticle
        
        ! Set previous deformation gradient
        Particle%DeformationGradientPrevious = matmul(Particle%RelativeDeformationGradient,&
                                            Particle%DeformationGradientPrevious)
        
        end subroutine UpdateDeformationGradientPrevious
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       
        subroutine CalculateTensorialInfinitesimalStrain(Particle, TensorialStrainIncrement)
        !----------------------------------------------------------------
        ! To calculate infinitesimal strains from relative deformation gradient
        !----------------------------------------------------------------
        implicit none 
        
        type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
        
        real(REAL_TYPE), dimension(NDIM, NDIM) :: RelativeDeformationGradientTransposed
        real(REAL_TYPE), dimension(NDIM, NDIM) :: IdentityMatrix
        real(REAL_TYPE), dimension(NDIM, NDIM), intent(out) :: TensorialStrainIncrement

        integer(INTEGER_TYPE) :: IParticle, IDim
        
        ! initialize 
        TensorialStrainIncrement = 0.0
        IdentityMatrix = 0.0
        
        do IDim = 1, NDIM
            IdentityMatrix(IDim, IDim) = 1.0
        end do
        
        
        RelativeDeformationGradientTransposed = transpose(Particle%RelativeDeformationGradient)
        
        TensorialStrainIncrement = 0.5 * ( Particle%RelativeDeformationGradient + RelativeDeformationGradientTransposed ) &
             - IdentityMatrix 
        
        
        
        end subroutine CalculateTensorialInfinitesimalStrain
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       
        subroutine CalculateEpsVoigt(Particle, TensorialStrainIncrement, Eps)
        !----------------------------------------------------------------
        ! To calculate infinitesimal strains from relative deformation gradient
        !----------------------------------------------------------------
        implicit none 
        
        type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)), intent(inout) :: Particle
        
        real(REAL_TYPE), dimension(NDIM, NDIM) :: RelativeDeformationGradientTransposed
        !real(REAL_TYPE), dimension(NDIM, NDIM) :: IdentityMatrix
        real(REAL_TYPE), dimension(NDIM, NDIM), intent(in) :: TensorialStrainIncrement
        real(REAL_TYPE), dimension(NTENSOR), intent(out) :: Eps

        integer(INTEGER_TYPE) :: IParticle, IDim
        
        ! initialize 
        Eps = 0.0
        
        
        !TensorialStrainIncrement = 0.5 * ( Particle%RelativeDeformationGradient + RelativeDeformationGradientTransposed ) &
        !     - IdentityMatrix 
        
        if (NDIM == 2) then 
            
            Eps(1) = TensorialStrainIncrement(1,1)
            Eps(2) = TensorialStrainIncrement(2,2)
            Eps(3) = 0.0
            Eps(4) = 0.5*(TensorialStrainIncrement(1,2) + TensorialStrainIncrement(2,1))
            
        elseif (NDIM == 3) then 
            
            Eps(1) = TensorialStrainIncrement(1,1)
            Eps(2) = TensorialStrainIncrement(2,2)
            Eps(3) = TensorialStrainIncrement(3,3)
            Eps(4) = 0.5*(TensorialStrainIncrement(1,2) + TensorialStrainIncrement(2,1))
            Eps(5) = 0.5*(TensorialStrainIncrement(2,3) + TensorialStrainIncrement(3,2))
            Eps(6) = 0.5*(TensorialStrainIncrement(3,1) + TensorialStrainIncrement(1,3))
        
        end if 
        
        
        end subroutine CalculateEpsVoigt
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine MapWeightedUpdatedDeformationGradientJacobianToNodes(WeightedUpdatedDeformationGradientJacobian)
        !ELEMENTNODES, 
        !ShapeValuesArray, 
        !IEN, 
        !) !(Momentum)
        !**********************************************************************
        !
        !    Function:  To map weighted deformation gradient to the nodes. 
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
        implicit none
        
        !integer(INTEGER_TYPE), intent(in) :: IPatch
        ! this is a scalar for each node
        real(REAL_TYPE), dimension(Maximum_NControlPoints,Counters%nEntity), intent(out) :: WeightedUpdatedDeformationGradientJacobian       !CC added nEntity
        
        ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn ! INode,
          !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: INode !NVECTOR, 
          
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) :: ParticleVolStrain
          real(REAL_TYPE) :: Time, PartilceMass
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: aa 
          
          !!CC - changed to function call - needs to be set to zero here
          WeightedUpdatedDeformationGradientJacobian = 0.0
          NoEn = Counters%nEntity
       
          do IPatch = 1, Counters%NPatches ! Loop over patches
          
              do IAEl = 1, nael_NURBS(IPatch) ! Loop over all elements 
                  
                  ! the background grid should have the same exact active element discretization 
                  IEl = ActiveElement(IAEl, IPatch)
            
                  !IDof(I, 1:ELEMENTNODES) = &
                  !    ReducedDof(Multipatch_Connecting_Local_To_Global_ControlPoints(ElementConnectivities(1:ELEMENTNODES, IEl, IPatch), IPatch) ) + I

                      ! degree of freedom for active element 
                      INode(1:ELEMENTNODES) = &
                          IEN(1:ELEMENTNODES, IEl, IPatch) !ReducedDof_VolLockSmooth( ) + I
              
                      
                      
                  !end do
            
                  do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                      ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                      if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       
                          PartilceMass = MassArray(ParticleIndex) ! particle mass
                          ParticleShape = ShapeValuesArray(ParticleIndex,:) ! particle shape values array on the lower order grid
                 
                          ! reset volumetric strain 
                          !Particles(ParticleIndex)%VolumetricStrain = 0.0 
                 
                          if (CalParams%ApplyContactAlgorithm) then
                             iEntity = EntityIDArray(ParticleIndex) 
                          else
                             iEntity = 1
                          end if 
                 
                          
                          
                          !do I = 1, NTENSOR ! nodal i-momentum
                              WeightedUpdatedDeformationGradientJacobian(INode(1:ELEMENTNODES), iEntity) = &
                                  WeightedUpdatedDeformationGradientJacobian(INode(1:ELEMENTNODES), iEntity ) &
                                  + ( PartilceMass * ParticleShape * &
                                  Particles(ParticleIndex)%RelativeDeformationGradientJacobian * &
                                  Particles(ParticleIndex)%DeformationGradientPreviousJacobian )
                          !end do     
       
                      end if ! NumbOfLayers = 1 or SOLID MatPoint
            
                  end do !Loop over particles
              
              end do !elements    
          
          end do ! patches
       
        end subroutine MapWeightedUpdatedDeformationGradientJacobianToNodes  
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        !-----------------------------------------------------------------
        
        
        
        
        subroutine NormalizeWeightedUpdatedDeformationGradientJacobianAtNodes(UpdatedDeformationGradientJacobian, &
            WeightedUpdatedDeformationGradientJacobian)
        !**********************************************************************
        !
        !    Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        
        integer(INTEGER_TYPE) :: IDOF, J, IPatch_Temporary
        integer(INTEGER_TYPE) :: INode, ITens
        
        real(REAL_TYPE), dimension(Maximum_NControlPoints, Counters%nEntity), intent(out) :: UpdatedDeformationGradientJacobian
        real(REAL_TYPE), dimension(Maximum_NControlPoints, Counters%nEntity), intent(in) :: WeightedUpdatedDeformationGradientJacobian
        
        
        IPatch_Temporary = 1
        ! initialisation --> not necessary here
        UpdatedDeformationGradientJacobian = 0.0
        
        do INode = 1, NControlPoints(IPatch_Temporary) !Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
            
                IDOF = ReducedDof(INode)+1 ! obtain the x dof so that we can find the mass. 
                                           ! It does not matter what dof you choose (here x-direction), 
                                           ! because mass is scalar. 
                
                !do ITens = 1, NTENSOR
                
                    if (LumpedMassDry(IDOF,J)/=0) then
                        UpdatedDeformationGradientJacobian(INode,J) = ( WeightedUpdatedDeformationGradientJacobian(INode,J) &
                                                                            / LumpedMassDry(IDOF,J) )
                    else
                        UpdatedDeformationGradientJacobian(INode,J) = 0.0 
                    end if
                
                !end do 
            
            end do
          
        end do
        
        
            end subroutine NormalizeWeightedUpdatedDeformationGradientJacobianAtNodes
            
            
            
            
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !subroutine UpdateParticleStrains_VolLockSmooth()
        !!**********************************************************************
        !!
        !!    Function:  Updates the total and incremental particle strains before mesh resetting
        !!
        !!**********************************************************************
        !
        !  implicit none
        !
        !  ! local variables
        !  real(REAL_TYPE), dimension(Counters%N_VolLockSmooth) :: DUTotEnt
        !  real(REAL_TYPE), dimension(NTENSOR) :: Eps
        !  real(REAL_TYPE), dimension(NDIM, ELEMENTNODES_VolLockSmooth) :: BMatrixDeformed
        !  real(REAL_TYPE) :: DetJac
        !  integer(INTEGER_TYPE) :: IElement, NElemPart, IParticle, ParticleIndex, I, iEntity, iEntityDefault, IActiveElement
        !
        !  
        !  ! Multipatch variables 
        !  integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        !  integer(INTEGER_TYPE) :: IPatch
        !  
        !  ! Volumetric locking variables 
        !  integer(INTEGER_TYPE) :: IDim 
        !  
        !  ! set the nodal displacement increments to the first entity
        !  ! only need to change if the contact model is used and the next particle belongs to a different entity
        !  iEntityDefault = 1 ! set default to first entity
        !  do I = 1,Counters%N_VolLockSmooth
        !    DUTotEnt(I) = IncrementalDisplacementSoil_VolLockSmooth(I,iEntityDefault)
        !  end do
        !
        !  
        !  !do IActiveElement = 1, Counters%NAEl ! loop over all active elements
        !   
        !  do IPatch = 1, Counters%NPatches ! loop over all patches
        !  do IActiveElement = 1, NAEl_NURBS(IPatch) !Counters%NAEl ! loop over all elements
        !
        !      IElement = ActiveElement(IActiveElement, IPatch) ! get element number of the active element
        !      NElemPart = NPartEle(IElement,IPatch) ! get number of particles in the active element
        !      
        !      if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement, IPatch) ) then
        !        NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
        !      end if  
        !      
        !      do IParticle = 1, NElemPart ! loop over all particles of the active element
        !      
        !        ParticleIndex = GetParticleIndex(IParticle, IElement, IPatch) ! get global particle ID
        !        
        !        if( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then 
        !            
        !          ! determine B matrix of deformed element at the particle local position
        !          call BMatrix_VolLockSmooth(Particles(ParticleIndex)%LocPos, &
        !                       ELEMENTNODES_VolLockSmooth, nel_NURBS(IPatch), &!Counters%NEl, & !nel_NURBS(IPatch), &
        !                       NControlPoints_VolLockSmooth(IPatch), NDIM, & !Counters%NodTot
        !                       IElement, IEN_VolLockSmooth(:,:,IPatch), &
        !                       NodalCoordinates_VolLockSmooth(:,:,IPatch),  & !            NodalCoordinatesUpd(:,:,IPatch),  &
        !                       BMatrixDeformed, &
        !                       IPatch)     !, DetJac) !update Jacob -> DetJac = | dx/dtildexi |
        !          
        !          ! only in UL FEM do we update the mesh --> NodalCoordinatesUpd
        !          
        !          ! get the nodal displacement increments for the specific entity
        !          !if (CalParams%ApplyContactAlgorithm) then ! CONTACT IS APPLIED
        !          !  iEntity = EntityIDArray(ParticleIndex) ! get the entity ID for the current particle
        !          !  if (iEntity /= iEntityDefault) then ! the particle belongs to a different entity from the previous entity
        !          !    DUTotEnt(:) = IncrementalDisplacementSoil(1:Counters%N,iEntity)
        !          !    iEntityDefault = iEntity !set the entity ID for the next particle loop
        !          !  end if
        !          !end if
        !            
        !          ! Eps in 3D (Exx, Eyy, Ezz, Gxy, Gyz, Gzx), in 2D (Exx, Eyy, Ezz, Gxy)
        !          call Get_Strain_VolLockSmooth(IElement, ParticleIndex, &!IParticle, &
        !                        IEN_VolLockSmooth(:,:,IPatch), & 
        !                        BMatrixDeformed, DUTotEnt, ReducedDof_VolLockSmooth, &
        !                        Eps, &
        !                        IPatch) ! Eps = strain increment at material point
        !          
        !          
        !          !call Get_Strain(IElement, ParticleIndex, &!IParticle, &
        !          !              ElementConnectivities(:,:,IPatch), & 
        !          !              BMatrixDeformed, DUTotEnt, ReducedDof, &
        !          !              Eps, &
        !          !              IPatch) ! Eps = strain increment at material point
        !          
        !          
        !          !! we have the strain increment Eps --> VolLockSmooth_BOOKMARK
        !          !! - get deviatoric strain increment 
        !          !! - get volumetric strain increment 
        !          !! - project the volumetric strain into the lower order shape functions 
        !          !! - project the volumetric strain back to the material points using the lower order shape functions 
        !          !! - sum the projected volumetric strain with the original deviatoric strain to get the total strain 
        !          !
        !          !if (NDIM == 2) then 
        !          !    
        !          !    ! ( sigma_xx + sigma_yy )/2
        !          !    MeanStrainIncrement = ( Eps(1) + Eps(2) ) * (1.0/2.0)
        !          !    
        !          !    ! volumetric strain = mean effective stress on the diagonals
        !          !    VolumetricStrainIncrement = 0.0
        !          !    VolumetricStrainIncrement(1) = MeanStrainIncrement
        !          !    VolumetricStrainIncrement(2) = MeanStrainIncrement
        !          !    
        !          !    ! deviatoric strain = total strain - volumetric strain
        !          !    DevStrainIncrement = Eps - VolumetricStrainIncrement
        !          !    
        !          !    ! 
        !          !    call MapVolStrainToNodesLowerOrderShapeFunctions
        !          !     
        !          !elseif (NDIM == 3) then 
        !          !    
        !          !    ! ( sigma_xx + sigma_yy + sigma_zz )/3
        !          !    HydrostaticStressIncrement = ( Eps(1) + Eps(2) + Eps(3) ) * (1.0/3.0)
        !          !
        !          !end if 
        !          
        !          
        !          ! Update particle strains
        !          call SetEpsStep_VolLockSmooth(Particles(ParticleIndex), Eps)
        !          
        !          
        !          ! calculate volumetric strain 
        !          Particles(ParticleIndex)%VolumetricStrainIncrement_VolLockSmooth = 0.0
        !          do IDim = 1, NVECTOR 
        !              Particles(ParticleIndex)%VolumetricStrainIncrement_VolLockSmooth = &
        !                  Particles(ParticleIndex)%VolumetricStrainIncrement_VolLockSmooth + & 
        !                  Particles(ParticleIndex)%EpsStep_VolLockSmooth(IDim)
        !          end do 
        !          
        !          !if (.not.CalParams%ApplyImplicitQuasiStatic .and. .not.CalParams%ApplyBBar) then
        !          !  call IncreaseEps(Particles(ParticleIndex), Eps)
        !          !end if
        !
        !        end if ! Material Point can be MIXTURE or SOLID
        !
        !      end do ! particle loop
        !
        !  end do ! active element loop
        !  
        !  end do ! patch loop
        !
        !end subroutine UpdateParticleStrains_VolLockSmooth
        
        
        
        
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !subroutine UpdateNodalStressGradient()
        !!**********************************************************************
        !!
        !!    Function:  Updates the total and incremental particle strains before mesh resetting
        !!
        !!**********************************************************************
        !
        !  implicit none
        !
        !  ! local variables
        !  !real(REAL_TYPE), dimension(Counters%N) :: DUTotEnt
        !  real(REAL_TYPE), dimension(NVECTOR) :: MPStressGradient !Eps
        !  real(REAL_TYPE), dimension(NDIM, ELEMENTNODES) :: BMatrixDeformed
        !  real(REAL_TYPE) :: DetJac
        !  integer(INTEGER_TYPE) :: IElement, NElemPart, IParticle, ParticleIndex, I, iEntity, iEntityDefault, IActiveElement
        !
        !  ! Multipatch variables 
        !  integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        !  integer(INTEGER_TYPE) :: IPatch
        !  
        !  ! set the nodal displacement increments to the first entity
        !  ! only need to change if the contact model is used and the next particle belongs to a different entity
        !  iEntityDefault = 1 ! set default to first entity
        !   
        !  do IPatch = 1, Counters%NPatches ! loop over all patches
        !  do IActiveElement = 1, NAEl_NURBS(IPatch) !Counters%NAEl ! loop over all elements
        !
        !      IElement = ActiveElement(IActiveElement, IPatch) ! get element number of the active element
        !      NElemPart = NPartEle(IElement,IPatch) ! get number of particles in the active element
        !      
        !      if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement, IPatch) ) then
        !        NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
        !      end if  
        !      
        !      do IParticle = 1, NElemPart ! loop over all particles of the active element
        !      
        !        ParticleIndex = GetParticleIndex(IParticle, IElement, IPatch) ! get global particle ID
        !        
        !        if( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then 
        !            
        !          ! determine B matrix of deformed element at the particle local position
        !          call BMatrix(Particles(ParticleIndex)%LocPos, &
        !                       ELEMENTNODES, nel_NURBS(IPatch), &!Counters%NEl, & !nel_NURBS(IPatch), &
        !                       NControlPoints(IPatch), NDIM, & !Counters%NodTot
        !                       IElement, ElementConnectivities(:,:,IPatch), &
        !                       NodalCoordinatesUpd(:,:,IPatch),  &
        !                       BMatrixDeformed, &
        !                       IPatch)!, DetJac) !update Jacob -> DetJac = | dx/dtildexi |
        !          
        !          !subroutine Get_StressGradient(IEl, IPoint, ICon, B,  NDofEx, IPatch) !Eps, !Disp,
        !          ! Eps in 3D (Exx, Eyy, Ezz, Gxy, Gyz, Gzx), in 2D (Exx, Eyy, Ezz, Gxy)
        !          call Get_StressGradient(IElement, ParticleIndex, &!IParticle, &
        !                        ElementConnectivities(:,:,IPatch), & 
        !                        BMatrixDeformed,  ReducedDof, & !DUTotEnt,
        !                        IPatch, &
        !                        MPStressGradient) ! Eps = strain increment at material point
        !          
        !          call SetStressGradient(Particles(ParticleIndex), MPStressGradient)
        !
        !        end if ! Material Point can be MIXTURE or SOLID
        !
        !      end do ! particle loop
        !
        !  end do ! active element loop
        !  
        !  end do ! patch loop
        !
        !end subroutine UpdateNodalStressGradient
        
        
        
        
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

          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          
          !--------------------------------------------------------

          !========================================================             !!CC
          !set the nodal displacement increments to the first entity
          !only need to change if the contact model is used and 
          !the next particle belongs to a different entity
          iEntityDefault = 1    !set default to first entity
          DULiqEnt = IncrementalDisplacementWater(:,iEntityDefault)
          !=======================================================

         !----------------------------------------------------
          do IPatch = 1, Counters%NPatches
          do IElement = 1, nel_NURBS(IPatch)!Counters%NEl     !loop over all elements
            if (IsActiveElement(IElement, IPatch)) then  !active elements only
              
                NElemPart = NPartEle(IElement,IPatch)    !number of particles in element
              if ( ISAXISYMMETRIC .and. .not. IsParticleIntegration(IElement, IPatch) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if
              
              do IParticle = 1, NElemPart       !loop over all particles of the element
                ParticleIndex = GetParticleIndex(IParticle, IElement, IPatch)   !particle global ID
                if((MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid).and. &
                         (.not.(NFORMULATION==1))) then !Material Point LIQUID

                  ! Determine B matrix of deformed element at the particle local position
                  call BMatrix(Particles(ParticleIndex)%LocPos, &
                             ELEMENTNODES, nel_NURBS(IPatch), &!Counters%NEl,  &
                             NControlPoints(IPatch), NDIM, & ! Counters%NodTot
                             IElement, ElementConnectivities(:,:,IPatch), &
                             NodalCoordinatesUpd(:,:,IPatch), &!NodalCoordinatesUpd, &
                             BMatrixDeformed, &
                             IPatch)!, DetJac)

                  ! Eps in 3D (Exx, Eyy, Ezz, Gxy, Gyz, Gzx), in 2D (Exx, Eyy, Ezz, Gxy)
                  call Get_Strain(IElement, IParticle, &
                                ElementConnectivities(:,:,IPatch), & 
                                BMatrixDeformed, DULiqEnt, ReducedDof, &
                                Eps, &
                                IPatch)

                  ! Update particle strains
                  call SetEpsStep(Particles(ParticleIndex), Eps)
                  call IncreaseEps(Particles(ParticleIndex), Eps)

                  ! Update Volumetric Strain Water
                  Particles(ParticleIndex)%WaterVolumetricStrain =  getEpsV(Eps)

                end if !Material Point LIQUID
              end do
            end if
          end  do ! loop over elements
          end do ! loop over patches

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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: PatchID
          integer(INTEGER_TYPE) :: ParticleID
          
          GetParticleIndexCounter = Counters%NParticles

           do IParticle = 1, Counters%NParticles ! Loop over particles

               
               
               
            PatchID = PatchIDArray(IParticle)    ! look up patch ID 
            ElementID = ElementIDArray(IParticle)! loop up local element ID in the PatchID
            
            
            !ParticleID = GetParticleIndex(IParticle, ElementID, PatchID)
            

            ! Calculate particle displacements from nodal displacements
            call CalcParticleDisplacements(IParticle, ElementID, PatchID) !IParticle

            GlopPos0 =  GlobPosArray(IParticle,:)  ! store the initial global position of the particle !IParticle

            ! Calculate new particle global positions
            GlobPosArray(IParticle,:) = GlobPosArray(IParticle,:) + UStepArray(IParticle,:) ! IParticle

           if (CalParams%ApplyFixParticlesK0) then ! only 3D functionality
             if (GlobPosArray(IParticle,1) >= (CalParams%Fix(1))) then ! allow sliding along y-axis only ! IParticle
                 GlobPosArray(IParticle,1) = GlopPos0(1) ! IParticle
                 GlobPosArray(IParticle,3) = GlopPos0(3) ! IParticle
             end if
             if (GlobPosArray(IParticle,2) <= (CalParams%Fix(2))) then ! IParticle
                GlobPosArray(IParticle,:) = GlopPos0 ! particles along the bottom are fixed ! IParticle
             end if
           end if ! ApplyFixParticlesK0

            ! Determine new local particle position from global position
            call DetermineLocPosIterative (IParticle, & !IParticle
                                 ElementID, GlopPos0, &
                                 Particles(IParticle)%LocPos, & ! IParticle
                                 NewElementID, LocPos, Success, &
                                 PatchID)
        
            
            ! Update Particle%LocPos, Particle%ElementID and EleParticles array
            if (Success) then
              Particles(IParticle)%LocPos = LocPos ! IParticle
              call SetParticleElementID(IParticle, IParticle, NewElementID) ! IParticle
              !! I added this to see how this is going to go... 
              !GetParticleIndex(IParticle, ElementIDArray(IParticle), PatchIDArray(IParticle)) = IParticle

            else ! Error in determining the new location of the particle
              call GiveError('  Error in determining the local ' // &
                             'coordinate of particle '// trim(String(IParticle)) // & ! IParticle
                             ' in element '// trim(String(NewElementID)) // &
                             ' in patch ' // trim(String(PatchID)))
            end if

            ! Calculate new particle shape function values -> after updating position 
            call SetParticleShapeFunctionData(Particles(IParticle), IParticle, ElementIDArray(IParticle), ElementConnectivities(:,:,PatchIDArray(IParticle)), NodalCoordinates(:,:,PatchIDArray(IParticle)), PatchIDArray(IParticle)) ! -> this needs to be for singular MP
            ! Calculate new particle shape function values -> after updating position --> Volumetric locking smoothing
            !if (CalParams%ApplyBBar) then
            !    call SetParticleShapeFunctionData_VolLockSmooth(Particles(IParticle), IParticle, ElementIDArray(IParticle), IEN_VolLockSmooth(:,:,PatchIDArray(IParticle)), NodalCoordinates_VolLockSmooth(:,:,PatchIDArray(IParticle)), PatchIDArray(IParticle)) ! -> this needs to be for singular MP
            !end if 
            
                ! IParticle
           end do
           
           

        end subroutine UpdateParticlePos

        subroutine CalcParticleDisplacements(ParticleIndex, ElementID, IPatch)
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

          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
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
                  do INode = 1, ELEMENTNODES!(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)!ELEMENTNODES  ! Loop over all nodes in element
                    NodeID = &
                        Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID,IPatch)), IPatch)
                    DofID = ReducedDof(NodeID) + IDof
                      ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + &
                                                   ( ShapeValuesArray(ParticleIndex,INode) * &
                                                   IncrementalDisplacementSoil(DofID,iEntity) )
                  enddo
                enddo
              else
                do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
                  do INode = 1, ELEMENTNODES!(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)!ELEMENTNODES  ! Loop over all nodes in element
                    NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID, IPatch)), IPatch)
                    DofID = ReducedDof(NodeID) + IDof
                    ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + &
                                                 ( ShapeValuesArray(ParticleIndex,INode) * &
                                                 AccumulatedIncDisplacementSoil(DofID,iEntity) )
                  enddo
                enddo
              endif
            endif

            if ((MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid).and. &
                         (.not.(NFORMULATION==1))) then !Material Point LIQUID
              do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
                do INode = 1, ELEMENTNODES  ! Loop over all nodes in element
                  NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID, IPatch) ), IPatch )
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
                       OldLocPos, NewElementID, NewLocPos, Success, IPatch)
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
          integer(INTEGER_TYPE):: ElemTrial
          ! Local variables
          logical :: OutsideElement, Fixed
          
          integer(INTEGER_TYPE) :: while_counter 
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(inout) :: IPatch
          integer(INTEGER_TYPE) :: PatchTrial
          integer(INTEGER_TYPE) :: IParticle
          
          NewElementID = OldElementID
          NewLocPos = OldLocPos
          Success = .true.

          
          if (NDIM ==2 ) then 
              
              ! TODO: I commented this temporarily 5/13/2024
              !call GetLocalCoordinates_2D(GlobPosArray(ParticleID,:), &
              !                     UStepArray(ParticleID,:), &
              !                     OldElementID, &
              !                     ELEMENTNODES, &
              !                     nel_NURBS(IPatch),  & 
              !                     NControlPoints(IPatch), NDIM, & 
              !                     NodalCoordinates(:,:,IPatch), & 
              !                     ElementConnectivities(1:ELEMENTNODES,1:nel_NURBS(IPatch),IPatch), &
              !                     NewLocPos, &
              !                     OutsideElement, Success, ElementAdjacencies(1:nel_NURBS(IPatch),1:ELEMENTSIDES, IPatch), IPatch)!, &
                  
                  call GetLocalCoordinates_StructuredMesh(GlobPosArray(ParticleID,:),  OldElementID, ELEMENTNODES, nel_NURBS(IPatch),  NControlPoints(IPatch), NDIM, NodalCoordinates(:,:,IPatch), ElementConnectivities(1:ELEMENTNODES,1:nel_NURBS(IPatch),IPatch), &
                                     NewLocPos, OutsideElement, Success, IPatch)
                  
          elseif (NDIM ==3 ) then 
              
          
          call GetLocalCoordinates(GlobPosArray(ParticleID,:), &
                                   UStepArray(ParticleID,:), &
                                   OldElementID, &
                                   ELEMENTNODES, &
                                   nel_NURBS(IPatch),  & 
                                   NControlPoints(IPatch), NDIM, & 
                                   NodalCoordinates(:,1:NDIM,IPatch), & 
                                   ElementConnectivities(:,:,IPatch), & 
                                   NewLocPos, &
                                   OutsideElement, Success, ElementAdjacencies(1:nel_NURBS(IPatch),1:ELEMENTSIDES, IPatch), IPatch)
          end if 
          
          !--------------------------------------------------

           !if (OutsideElement) then ! Particle moved to another element
           !
           !     
           !     do PatchTrial = 1, Counters%NPatches ! loop over patches
           !     do ElemTrial = 1, nel_NURBS(PatchTrial)!CalParams%NElem ! loop over elements
           !     
           !
           !         if (NDIM ==2 ) then 
           !             
           !             call GetLocalCoordinates_2D(GlobPosArray(ParticleID,:), &
           !                        UStepArray(ParticleID,:), &
           !                        ElemTrial, &
           !                        ELEMENTNODES, nel_NURBS(PatchTrial),  & !Counters%NEl
           !                        NControlPoints(PatchTrial), NDIM, & !Counters%NodTot
           !                        NodalCoordinates(:,:,PatchTrial), &
           !                        ElementConnectivities(:,:,PatchTrial), &
           !                        NewLocPos, &
           !                        OutsideElement, Success, ElementAdjacencies(:,:, PatchTrial), PatchTrial)!, &
           !             
           !             
           !         elseif (NDIM ==3) then 
           !             
           !         call GetLocalCoordinates(GlobPosArray(ParticleID,:), &
           !                        UStepArray(ParticleID,:), &
           !                        ElemTrial, &
           !                        ELEMENTNODES, nel_NURBS(PatchTrial),  & !Counters%NEl
           !                        NControlPoints(PatchTrial), NDIM, & !Counters%NodTot
           !                        NodalCoordinates(:,:,PatchTrial), &
           !                        ElementConnectivities(:,:,PatchTrial), &
           !                        NewLocPos, &
           !                        OutsideElement, Success, ElementAdjacencies(:,:, PatchTrial), PatchTrial)!, &
           !             
           !         end if 
           !         
           !             if (success) then 
           !                  ElementIDArray(ParticleID)= ElemTrial!NewElementID
           !                 PatchIDArray(ParticleID)= PatchTrial
           !                 NewElementID = ElemTrial
           !                 ! remove it from old patch
           !                 do IParticle = 1, Counters%NParticles
           !                     if (GetParticleIndex(IParticle,OldElementID,IPatch) == ParticleID) then 
           !                         GetParticleIndex(IParticle, OldElementID, IPatch) = 0
           !                         !GetParticleIndex(Counters%NParticles-IParticle, ElementIDArray(ParticleID), PatchIDArray(ParticleID)) = ParticleID
           !                         exit
           !                     end if 
           !                 end do 
           !                 
           !                 ! add it in new patch  
           !                 GetParticleIndex(GetParticleIndexCounter, ElementIDArray(ParticleID), PatchIDArray(ParticleID)) = ParticleID    
           !                 GetParticleIndexCounter = GetParticleIndexCounter-1
           !                  
           !                 call WriteInLogFile('  Particle '   // trim(String(ParticleID))        // &
           !                       ' moved from element'  // trim(String(OldElementID))      // &
           !                       ' in patch '  // trim(String(IPatch)) // &   
           !                       ' to '          // trim(String(ElemTrial))      // &
           !                       ' in patch '  // trim(String(PatchTrial)) // &   
           !                       ' at local pos:'// trim(String(NewLocPos, 1, 3))   // &
           !                       'TimeStep '     // trim(String(Calparams%Timestep))  )
           !                 
           !                 
           !                 exit
           !                 
           !             end if 
           !             
           !     !--------------------------------------------------------------------------
           !     end do ! loop over elements
           !     
           !     if (success) then 
           !                 !ElementIDArray(ParticleID)= ElemTrial!NewElementID
           !                 !PatchIDArray(ParticleID)= PatchTrial
           !                 !NewElementID = ElemTrial
           !                 !
           !                 !call WriteInLogFile('  Particle '   // trim(String(ParticleID))        // &
           !                 !      ' moved from element'  // trim(String(OldElementID))      // &
           !                 !      ' in patch '  // trim(String(IPatch)) // &   
           !                 !      ' to '          // trim(String(ElemTrial))      // &
           !                 !      ' in patch '  // trim(String(PatchTrial)) // &   
           !                 !      ' at local pos:'// trim(String(NewLocPos, 1, 3))   // &
           !                 !      'TimeStep '     // trim(String(Calparams%Timestep))  )
           !                 ! I added this to see how this is going to go... 
           !     
           !         
           !                 
           !                 exit
           !                 
           !             end if 
           !     
           !     end do ! loop over patches
           !     
           !     
           !     if (.not.Success) then
           !   Fixed = .false.
           !
           !   if (.not.Fixed) then
           !     call GiveError('  Failed finding new element for particle ' // trim(String(ParticleID)) // &
           !                    ' in element '  // trim(String(OldElementID)) // &
           !                    ' in patch '  // trim(String(IPatch)) // & 
           !                    ' at position ' // trim(String(GlobPosArray(ParticleID,:),1,NVECTOR)) // &
           !                    ' in time step '// trim(String(CalParams%TimeStep)) )
           !   end if
           !     end if
           !     
           ! end if
          !end if
          
          !--------------------------------------------------
          
          if (OutsideElement) then ! Particle moved to another element
           
               !while_counter = 0 
               !do while ( (Success == .false.) .and. (while_counter <5))  
               ! Determine new element and new local position of the particle
            call DetermineGlobPosElement(GlobPosArray(ParticleID,:), &
                                         UStepArray(ParticleID,:), &
                                         OldElementID, &
                                         ELEMENTNODES,  &
                                         nel_NURBS(PatchIDArray(ParticleID)), &!Counters%NEl,  &
                                         NControlPoints(PatchIDArray(ParticleID)), &!Counters%NodTot,  &
                                         NDIM, &
                                         NodalCoordinates(1:NControlPoints(IPatch),1:NDIM,IPatch),  &
                                         ElementConnectivities(1:ELEMENTNODES,1:nel_NURBS(IPatch),IPatch), &
                                         NewElementID, &
                                         NewLocPos, &
                                         Success, &
                                         IPatch)
            
            ElementIDArray(ParticleID)= NewElementID
            PatchIDArray(ParticleID)= IPatch!NewPatchID
            
            !while_counter = while_counter + 1
            
               !end do 
               
            if (Success) then
                
                ! I added this to see how this is going to go... 
                !GetParticleIndex(ParticleID, OldElementID, PatchIDArray(ParticleID)) = 0
                !GetParticleIndex(ParticleID, ElementIDArray(ParticleID), PatchIDArray(ParticleID)) = ParticleID
                !GetParticleIndex(IParticle, ElementIDArray(IParticle), PatchIDArray(IParticle)) = IParticle
                
                ! remove it from old patch
                            do IParticle = 1, Counters%NParticles
                                if (GetParticleIndex(IParticle,OldElementID,IPatch) == ParticleID) then 
                                    GetParticleIndex(IParticle, OldElementID, IPatch) = 0
                                    !GetParticleIndex(Counters%NParticles-IParticle, ElementIDArray(ParticleID), PatchIDArray(ParticleID)) = ParticleID
                                    exit
                                end if 
                            end do 
                            ! add it in new patch
                            !do IParticle = 1, Counters%NParticles
                            !    if (GetParticleIndex(IParticle,ElementIDArray(ParticleID),PatchIDArray(ParticleID)) == ParticleID) then 
                                    GetParticleIndex(GetParticleIndexCounter, ElementIDArray(ParticleID), PatchIDArray(ParticleID)) = ParticleID
                            !        exit
                            !    end if 
                            !end do 
                                    GetParticleIndexCounter = GetParticleIndexCounter-1

                
              if (CalParams%OutputDebugData) then
              call WriteInLogFile('  Particle '   // trim(String(ParticleID))        // &
                                  ' moved from element '  // trim(String(OldElementID))      // &
                                  ' in patch '  // trim(String(IPatch)) // &   
                                  ' to element '          // trim(String(NewElementID))      // &
                                  ' in patch '  // trim(String(IPatch)) // &   
                                  ' at local pos:'// trim(String(NewLocPos, 1, 3))   // &
                                  'TimeStep '     // trim(String(Calparams%Timestep))  )
              
              
              end if
            else ! New element could not be found
                ! loop accross all elements and check if it would fit into an element 
                
                do PatchTrial = 1, Counters%NPatches ! loop over patches
                do ElemTrial = 1, nel_NURBS(PatchTrial)!CalParams%NElem ! loop over elements
                
                !--------------------------------------------------------------------
                ! Abdel: repeat again to see if we can catch it in the new element
                !call DetermineGlobPosElement(GlobPosArray(ParticleID,:), &
                !                         UStepArray(ParticleID,:), &
                !                         ElemTrial, &
                !                         ELEMENTNODES,  &
                !                         Counters%NEl,  &
                !                         Counters%NodTot,  &
                !                         NDIM, &
                !                         NodalCoordinates,  &
                !                         ElementConnectivities, &
                !                         NewElementID, &
                !                         NewLocPos, &
                !                         Success)
                !ElementIDArray(ParticleID)= NewElementID
                    if (NDIM ==2 ) then 
                        
                        !call GetLocalCoordinates_2D(GlobPosArray(ParticleID,:), &
                        !           UStepArray(ParticleID,:), &
                        !           ElemTrial, &
                        !           ELEMENTNODES, nel_NURBS(PatchTrial),  & !Counters%NEl
                        !           NControlPoints(PatchTrial), NDIM, & !Counters%NodTot
                        !           NodalCoordinates(:,:,PatchTrial), &
                        !           ElementConnectivities(:,:,PatchTrial), &
                        !           NewLocPos, &
                        !           OutsideElement, Success, ElementAdjacencies(:,:, PatchTrial), PatchTrial)!, &
                        
                            call GetLocalCoordinates_StructuredMesh(GlobPosArray(ParticleID,:),  ElemTrial, ELEMENTNODES, &
                                nel_NURBS(PatchTrial),  NControlPoints(PatchTrial), NDIM, NodalCoordinates(:,:,PatchTrial),&
                                ElementConnectivities(1:ELEMENTNODES,1:nel_NURBS(PatchTrial),PatchTrial), &
                                     NewLocPos, OutsideElement, Success, IPatch)
                        
                    elseif (NDIM ==3) then 
                        
                    call GetLocalCoordinates(GlobPosArray(ParticleID,:), &
                                   UStepArray(ParticleID,:), &
                                   ElemTrial, &
                                   ELEMENTNODES, nel_NURBS(PatchTrial),  & !Counters%NEl
                                   NControlPoints(PatchTrial), NDIM, & !Counters%NodTot
                                   NodalCoordinates(:,:,PatchTrial), &
                                   ElementConnectivities(:,:,PatchTrial), &
                                   NewLocPos, &
                                   OutsideElement, Success, ElementAdjacencies(:,:, PatchTrial), PatchTrial)!, &
                        
                    end if 
                    
                        if (success) then 
                            ElementIDArray(ParticleID)= ElemTrial!NewElementID
                            PatchIDArray(ParticleID)= PatchTrial
                            NewElementID = ElemTrial
                            ! remove it from old patch
                            do IParticle = 1, Counters%NParticles
                                if (GetParticleIndex(IParticle,OldElementID,IPatch) == ParticleID) then 
                                    GetParticleIndex(IParticle, OldElementID, IPatch) = 0
                                    !GetParticleIndex(Counters%NParticles-IParticle, ElementIDArray(ParticleID), PatchIDArray(ParticleID)) = ParticleID
                                    exit
                                end if 
                            end do 
                            ! add it in new patch
                            !do IParticle = 1, Counters%NParticles
                            !    if (GetParticleIndex(IParticle,ElementIDArray(ParticleID),PatchIDArray(ParticleID)) == ParticleID) then 
                                    GetParticleIndex(GetParticleIndexCounter, ElementIDArray(ParticleID), PatchIDArray(ParticleID)) = ParticleID
                            !        exit
                            !    end if 
                            !end do 
                                    GetParticleIndexCounter = GetParticleIndexCounter-1
                             
                                    
                            
                            !GetParticleIndex(ParticleID, ElementIDArray(ParticleID), PatchIDArray(ParticleID)) = ParticleID
                            
                            call WriteInLogFile('  Particle '   // trim(String(ParticleID))        // &
                                  ' moved from element'  // trim(String(OldElementID))      // &
                                  ' in patch '  // trim(String(IPatch)) // &   
                                  ' to '          // trim(String(ElemTrial))      // &
                                  ' in patch '  // trim(String(PatchTrial)) // &   
                                  ' at local pos:'// trim(String(NewLocPos, 1, 3))   // &
                                  'TimeStep '     // trim(String(Calparams%Timestep))  )
                            
                            
                            exit
                            
                        end if 
                        
                !--------------------------------------------------------------------------
                end do ! loop over elements
                
                if (success) then 
                            !ElementIDArray(ParticleID)= ElemTrial!NewElementID
                            !PatchIDArray(ParticleID)= PatchTrial
                            !NewElementID = ElemTrial
                            !
                            !call WriteInLogFile('  Particle '   // trim(String(ParticleID))        // &
                            !      ' moved from element'  // trim(String(OldElementID))      // &
                            !      ' in patch '  // trim(String(IPatch)) // &   
                            !      ' to '          // trim(String(ElemTrial))      // &
                            !      ' in patch '  // trim(String(PatchTrial)) // &   
                            !      ' at local pos:'// trim(String(NewLocPos, 1, 3))   // &
                            !      'TimeStep '     // trim(String(Calparams%Timestep))  )
                            ! I added this to see how this is going to go... 
                
                    
                            
                            exit
                            
                        end if 
                
                end do ! loop over patches
                
                
                if (.not.Success) then
              Fixed = .false.
          
              if (.not.Fixed) then
                call GiveError('  Failed finding new element for particle ' // trim(String(ParticleID)) // &
                               ' in element '  // trim(String(OldElementID)) // &
                               ' in patch '  // trim(String(IPatch)) // & 
                               ' at position ' // trim(String(GlobPosArray(ParticleID,:),1,NVECTOR)) // &
                               ' in time step '// trim(String(CalParams%TimeStep)) )
              end if
                end if
                
            end if
          end if
          
          
          !--------------------------------------------------
          
          
          
          

        end subroutine DetermineLocPosIterative

        subroutine UpdateParticleHouseKeeping_MPInitFile()
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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch

          
          NPartEle = 0
          !I commented this because it was sorting things in a manner that was breaking the code 
          !call IHpSortKind8(EleParticles, Counters%NParticles)
          ! Count the particles in element, update NPartEle
          do IParticle = 1, Counters%NParticles ! Loop over all particles
              
            IPatch = PatchIDArray(IParticle)  
            IElement = ElementIDArray(IParticle)!GetElementIDFromList(IParticle) ! Determine the element that IParticle is located in
            NPartEle(IElement,IPatch) = NPartEle(IElement,IPatch) + 1 ! Increase the particle counter for IElement
          end do


          ! Update the helper array and element switches
          
          do IPatch = 1, Counters%NPatches ! Loop over patches
              
              !FirstParticle = 1!GetParticleIndex(1, )!1
              counters%NAEl = 0
              nael_NURBS(IPatch) = 0
              
          do IElement = 1, nel_NURBS(IPatch)!Counters%NEl ! Loop over all elements
            
              FirstParticle = GetParticleIndex(1, IElement, IPatch)!1
              
              if (NPartEle(IElement,IPatch)==0) then ! Empty element
              EleParticlesHelp(IElement, IPatch) = -1

              IsActiveElement(IElement, IPatch) = .false. ! Switch off the element
            else  ! Element contains particles
              !counters%NAEl = counters%NAEl + 1
              nael_NURBS(IPatch) = nael_NURBS(IPatch) + 1 ! <-- might need to find the sum of this to know the total number of active elements
              EleParticlesHelp(IElement, IPatch) = FirstParticle
              FirstParticle = FirstParticle + NPartEle(IElement,IPatch)

              IsActiveElement(IElement, IPatch) = .true. ! Switch on the element
            end if
          end do ! loop over elements
          end do ! loop over patches

        end subroutine UpdateParticleHouseKeeping_MPInitFile
        
        
        
        subroutine UpdateParticleHouseKeeping_MPInitGP()
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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch

          
          NPartEle = 0
          !I commented this because it was sorting things in a manner that was breaking the code 
          call IHpSortKind8(EleParticles, Counters%NParticles)
          ! Count the particles in element, update NPartEle
          do IParticle = 1, Counters%NParticles ! Loop over all particles
              
            IPatch = PatchIDArray(IParticle)  
            IElement = ElementIDArray(IParticle)!GetElementIDFromList(IParticle) ! Determine the element that IParticle is located in
            NPartEle(IElement,IPatch) = NPartEle(IElement,IPatch) + 1 ! Increase the particle counter for IElement
          end do


          ! Update the helper array and element switches
          
          do IPatch = 1, Counters%NPatches ! Loop over patches
              
              FirstParticle = 1!GetParticleIndex(1, )!1
              counters%NAEl = 0
              nael_NURBS(IPatch) = 0
              
          do IElement = 1, nel_NURBS(IPatch)!Counters%NEl ! Loop over all elements
              if (NPartEle(IElement,IPatch)==0) then ! Empty element
              EleParticlesHelp(IElement, IPatch) = -1

              IsActiveElement(IElement, IPatch) = .false. ! Switch off the element
            else  ! Element contains particles
              !counters%NAEl = counters%NAEl + 1
              nael_NURBS(IPatch) = nael_NURBS(IPatch) + 1 ! <-- might need to find the sum of this to know the total number of active elements
              EleParticlesHelp(IElement, IPatch) = FirstParticle
              FirstParticle = FirstParticle + NPartEle(IElement,IPatch)

              IsActiveElement(IElement, IPatch) = .true. ! Switch on the element
            end if
          end do ! loop over elements
          end do ! loop over patches

        end subroutine UpdateParticleHouseKeeping_MPInitGP
       

        subroutine SetUpEntityElements()
        !**********************************************************************
        !
        !  Function : For each entity a list of the elements containing particles belonging to this entity
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: IAEl, IEl, NElemPart, IPart, ParticleIndex, IEntity
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
         
          EntityElements = 0      !reset

          !do IAEl = 1,Counters%NAel                                      !loop through all elements
      
              
          do IPatch = 1, Counters%NPatches ! Loop over patches
          do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements
            IEl = ActiveElement(IAEl, IPatch)
            NElemPart = NPartEle(IEl,IPatch)
            do IPart= 1, NElemPart                                 !loop through particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl, IPatch)         !Get the particle ID
              IEntity = EntityIDArray(ParticleIndex)               !entity to which particle belong
              EntityElements(IEntity,IEl,IPatch) = 1
            end do
          end do ! loop over elements
          end do ! loop over patches
          
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
         
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          
          MaterialElements = 0      !reset
          
                
          !do IAEl = 1,Counters%NAel                                      !loop through all elements
      
          do IPatch = 1, Counters%NPatches ! Loop over patches
          do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements
            IEl = ActiveElement(IAEl, IPatch)
            !if (IEl == 0) cycle ! if no active element then cycle to the next patch
            NElemPart = NPartEle(IEl,IPatch)
            do IPart= 1, NElemPart                                 !loop through particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl, IPatch)        ! Get the particle ID
              MaterialID = MaterialIDArray(ParticleIndex)       !material ID to which particle belong
              MaterialElements(MaterialID,IEl, IPatch) = 1
            end do
          end do ! loop over elements
          end do ! loop over patches
                  
        end subroutine SetUpMaterialElements


        subroutine SetActiveElement()
        !**********************************************************************
        !
        !  Function : To set the active elements
        !
        !  ActiveElement: Active element vector
        !
        !**********************************************************************
        
        ! Abdelrahman Alsardi: Stamp of approval for multipatch generlization
        
        implicit none
        
          integer(INTEGER_TYPE) :: IAEl, IEl, IError
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
        
          if (allocated(ActiveElement) ) then
            deallocate(ActiveElement, stat = IError)
          end if
          
          sum_nael_NURBS = 0
          do IPatch = 1, Counters%NPatches
          sum_nael_NURBS = sum_nael_NURBS + nael_NURBS(IPatch)
          end do 
          
          !allocate(ActiveElement(counters%NAEl, Counters%NPatches), stat = IError) ! this needs to be the sum of all active elements
          allocate(ActiveElement(sum_nael_NURBS, Counters%NPatches), stat = IError) ! this needs to be the sum of all active elements
          ActiveElement = 0

          IAEl = 0
          do IPatch = 1, Counters%NPatches ! loop over patches
              IAEl = 0
          do IEl=1, nel_NURBS(IPatch) !counters%NEl ! loop over elements
            if (IsActiveElement(IEl, IPatch)) then ! Active elements
              IAEl = IAEl + 1
              ActiveElement(IAEl, IPatch) = IEl
            end if
          end do ! loop over elements
          end do ! loop over patches
 
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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch

          if (.not.IsMPMWithMixedIntegration()) RETURN
            IsParticleIntegration = .false.

            !do IAElement = 1, Counters%NAEl
            do IPatch = 1, Counters%NPatches ! Loop over patches
          
                do IAElement = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
   
              IElement = ActiveElement(IAElement, IPatch)
              ContainsBoundaryParticle = .false.
              SummedParticleWeights = 0.0

              ParticleIndex = GetParticleIndex(1, IElement, IPatch)
              CheckMaterialID = MaterialIDArray(ParticleIndex)
              IsMixedElement = .false.

              do IParticle = 1, NPartEle(IElement,IPatch) ! Loop over particles in IElement
                ParticleIndex = GetParticleIndex(IParticle, IElement, IPatch)
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
                RequiredFilledVolume = ElementSpace(IElement,IPatch) * CalParams%RequiredDegreeOfFilling
                if (SummedParticleWeights<=RequiredFilledVolume) then
                  IsParticleIntegration(IElement, IPatch) = .true. ! Consider the element as partially filled
                end if
              end if

              if (IsMixedElement) then ! Elements with particles of different material
                IsParticleIntegration(IElement, IPatch) = .true.
              end if
                end do
                
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
          
          integer(INTEGER_TYPE) :: IPatch
          
          FirstParticleWaterPressure = 0.0
          FirstParticleGasPressure = 0.0
          FirstParticleStress = 0.0
          FirstDBulkViscousPressure = 0.0
    
          !do IAElement = 1, Counters%NAEl
              do IPatch = 1, Counters%NPatches ! Loop over patches
          
                do IAElement = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
            IElement = ActiveElement(IAElement, IPatch)
            if ( (.not.IsParticleIntegration(IElement, IPatch) ) ) then ! Loop over all active fully filled elements
     
             FirstParticleIndex = GetParticleIndex(1, IElement, IPatch)
     
             FirstParticleStress = SigmaEffArray(FirstParticleIndex,:)
             FirstParticleWaterPressure = Particles(FirstParticleIndex)%WaterPressure
             FirstParticleGasPressure =  Particles(FirstParticleIndex)%GasPressure
             FirstDBulkViscousPressure = Particles(FirstParticleIndex)%DBulkViscousPressure
     
             NElemPart = NPartEle(IElement,IPatch)  ! Number of particles in element
             
             do IPart = 2, NElemPart ! loop over particles. The stress of the first particle is already assigned.
               
               ParticleIndex = GetParticleIndex(IPart, IElement, IPatch)
             
               SigmaEffArray(ParticleIndex,:) = FirstParticleStress
               Particles(ParticleIndex)%WaterPressure = FirstParticleWaterPressure
               Particles(ParticleIndex)%GasPressure = FirstParticleGasPressure
               Particles(ParticleIndex)%DBulkViscousPressure = FirstDBulkViscousPressure
             end do
                   
            end if
                end do
                
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
          
          integer(INTEGER_TYPE) :: IPatch
     
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
        
          !do IAElement = 1, Counters%NAEl
              
              do IPatch = 1, Counters%NPatches ! Loop over patches
          
                do IAElement = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
                    
            IElement = ActiveElement(IAElement, IPatch)
            if ( (.not.IsParticleIntegration(IElement, IPatch) ) ) then ! Loop over all active fully filled elements
     
             FirstParticleIndex = GetParticleIndex(1, IElement, IPatch)
     

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
     
             NElemPart = NPartEle(IElement,IPatch)  ! Number of particles in element
             
             do IPart = 2, NElemPart ! loop ove particles. The state variables of the first particle are assigned.
               
               ParticleIndex = GetParticleIndex(IPart, IElement, IPatch)
               
                
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
                
              end do ! loop over patches
              
      
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
          
          integer(INTEGER_TYPE) :: IPatch

           !do IAElement = 1, Counters%NAEl
               
               do IPatch = 1, Counters%NPatches ! Loop over patches
          
                do IAElement = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
                    
             
                    IElement = ActiveElement(IAElement, IPatch)
      
                    do IParticle = 1, NPartEle(IElement,IPatch) ! Loop over particles in IElement
                  ParticleIndex = GetParticleIndex(IParticle, IElement, IPatch)
                  
                  SigmaEffStep = SigmaEffArray(ParticleIndex,:) - SigmaEff0Array(ParticleIndex,:) 
                  call SetSigmaEffStep(Particles(ParticleIndex), SigmaEffStep)
                  
                  SigmaEff = SigmaEffArray(ParticleIndex,:)
                  SigmaEff0Array(ParticleIndex,:) = SigmaEff
                  
                  Particles(ParticleIndex)%WaterPressure0 = Particles(ParticleIndex)%WaterPressure
                  Particles(ParticleIndex)%GasPressure0 = Particles(ParticleIndex)%GasPressure
                  
                end do
                end do ! loop over active elements
                end do ! loop over patches
     
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
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: PatchID
          
          if (.not.CalParams%ApplyFEMtoMPM) RETURN

          do IParticle = 1, Counters%NParticles ! loop over material points
            
            ElementID = ElementIDArray(IParticle)
            PatchID = PatchIDArray(IParticle)

            ! Calculate particle displacements from nodal displacements
            call MapTotalDisplacementsToParticles(IParticle, ElementID, PatchID)

            ! Calculate particle velocities from nodal displacements
            call MapTotalVelocityToParticles(IParticle, ElementID, PatchID)
          
          end do
          
        end subroutine MapDataFromNodesToParticles 


        subroutine MapTotalDisplacementsToParticles(ParticleIndex, &
                                                    ElementID, &
                                                    PatchID)
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
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: PatchID

           ParticleDisplacement = 0.0

           do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
             do INode = 1, ELEMENTNODES  ! Loop over all nodes in element

               NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID, PatchID) ) , PatchID)
               DofID = ReducedDof(NodeID) + IDof
               ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(ParticleIndex,INode) * TotalDisplacementSoil(DofID)
               
             end do
           end do

           ! Total particle displacement
           UArray(ParticleIndex,:) =  UArray(ParticleIndex,:) + ParticleDisplacement

        end subroutine MapTotalDisplacementsToParticles

        subroutine MapTotalVelocityToParticles(ParticleIndex, ElementID, IPatch)
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
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary 
          integer(INTEGER_TYPE), intent(in) :: IPatch
          

            !!get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(ParticleIndex) 
            else
              iEntity = 1
            end if 
              
            ParticleVelocity = 0.0
            
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                  
                NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID, IPatch) ) , IPatch)
                DofID = ReducedDof(NodeID) + IDof
                ParticleVelocity(IDof) = ParticleVelocity(IDof) + ShapeValuesArray(ParticleIndex,INode) * TotalVelocitySoil(DofID,iEntity)
                
              end do                                             
            end do 
              
            ! Total particle velocity
             VelocityArray(ParticleIndex,:) = ParticleVelocity
          
        end subroutine MapTotalVelocityToParticles
                 
        


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
      use ModConvectivePhase, only: UpdateParticleHouseKeeping_MPInitGP!UpdateParticleHouseKeeping
      implicit none
      !call UpdateParticleHouseKeeping()
      !if (MP_INIT_FILE) then 
      !        call UpdateParticleHouseKeeping_MPInitFile()
      !    else ! MP_INIT_GP
              
      call UpdateParticleHouseKeeping_MPInitGP()
      !    end if
      
      end subroutine UpdateParticleHouseKeepingMock
      
      subroutine SetActiveElementMock()
      !**********************************************************************
      !
      !    Function:  set active elements
      !
      !*********************************************************************  
      use ModConvectivePhase, only: SetActiveElement
      implicit none
      ! hardcoded to remove active element error
      !call SetActiveElement()
      end subroutine SetActiveElementMock