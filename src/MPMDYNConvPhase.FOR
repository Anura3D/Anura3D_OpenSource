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
    !   and soilâ€“waterâ€“structure interaction using the material point method (MPM)
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
	  
	  
	  module ModDYNConvectivePhase
      !**********************************************************************
      !
      !    Function:  Contains the routines related to updating the mesh and
      !               particle data at the end of the Lagrangian Phase.
      !               This module is only used by the Dynamic MPM.
      !               In order to keep the size of this source file reasonably small,
      !               this module only contains routines that are directly related to
      !               the updating of particle data.
      !
      !     $Revision: 9793 $
      !     $Date: 2022-09-20 15:19:37 +0200 (di, 20 sep 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluation
      use ModMPMData
      use ModWriteTestData
      use ModLagrangianPhase
      use ModConvectivePhase
      use ModMPMDYNStresses
      use ModMeshAdjacencies
      use ModMPMMeshAdjustment
      use ModMeshInfo
      use ModRotBoundCond
      use ModEmptyElements
      use ModAdjustParticleDiscretisation
      use ModGlobalConstants
      use ModTwoLayerFormulation
      use ModString
      use ModRigidBody
      use ModLiquid
      use ModReadGeometryData
      
      implicit none

        real(REAL_TYPE), dimension(:),  &
          allocatable :: TemporaryMappingVector

      contains ! Routines of this module

      subroutine DYNConvectivePhase()
        !**********************************************************************
        !
        !    Function:  Calls the different subroutines required for updating the
        !               particle data. The basic steps are:
        !               
        !               - update particle velocity
        !               - update particle displacements and global position
        !               - map the new particles ve  locities to the nodes
        !               - update the nodal coordinates from the new nodal velocities
        !               - calculate particle strains 
        !               - calculate stresses for integration points
        !               - map stresses from Gauss points to particles for fully filled elements
        !               - reset the mesh
        !               - determine the elements that particles moved into
        !               - determine the new local particle coordinates
        !               - smoothen particle stresses within each element
        !               - update particle shape values
        !               - update the particle house-keeping data structure
        !               - ... further checks ...
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
         real(REAL_TYPE), dimension(Counters%N) :: NodalIncDisplacementMeshAdjust
         ! Local variables
         integer(INTEGER_TYPE) :: iOpt, I, EntityUsed

         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: MomentumW
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: MomentumG

         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: NodalIncDisplacementG

         real(REAL_TYPE), dimension(Counters%N) :: DummyRotationVector1

         integer(INTEGER_TYPE) :: IEntity, StructureMaterialEntity, LoadedSides, IDOF

          EntityUsed = 1   !!CC - each node has nEntity displacements, choose entity used for nodal update
          
          call GetRigidBodyAverageAcceleration()
         
          ! Update particles total velocities and Accelerations (accelerations added by CC)    
        
         call ApplyMPPrescribedVelocity(TotalVelocitySys)
         call ApplyNodalPrescribedVelocity(DummyRotationVector1)

         if (IsMPMComputation()) then ! if MPM update particle velocity....
           call UpdateParticleVelocityAndMapMomentum(Momentum) ! From accelerations in global coordinate system
           if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3).or.(.not.(NFORMULATION==1))) then
             call UpdateParticleWaterVelocityAndMapMomentumW(MomentumW) ! From accelerations in global coordinate system
           end if
           if (CalParams%NumberOfPhases==3) then
             call UpdateParticleGasVelocityAndMapMomentumG(MomentumG)
           end if
         else !FEM  get velocity directly at the nodes
           ! Accelerations are in global coordinate system but TotalVelocitySoil, too?
           ! TotalVelocitySoil is global if contact formulation used and local without contact
           ! TotalVelocityWater is local with and without contact
           ! Rotate TotalVelocitySoil in case of no contact to global coordinate system
           ! Rotate TotalVelocityWater to global coordinate system
           if (IS3DCYLINDRIC) then
             do IEntity = 1, Counters%nEntity
               call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocityWater(:, IEntity), TotalVelocityWater(:, IEntity))
               if (.not.CalParams%ApplyContactAlgorithm.and.(.not.IsMPMSkipConvection())) then
                 call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocitySoil(:, IEntity), TotalVelocitySoil(:, IEntity))
               endif
             end do
           end if

           do IEntity = 1, Counters%nEntity
             do I = 1, Counters%N
               TotalVelocitySoil(I, IEntity) = TotalVelocitySoil(I, IEntity) + AccelerationSoil(I, IEntity) * CalParams%TimeIncrement  ! solid
               if (CalParams%NumberOfPhases==3) then
                 TotalVelocityGas (I, IEntity)= TotalVelocityGas(I, IEntity)+ AccelerationGas(I, IEntity) * CalParams%TimeIncrement
               end if
             end do
           end do

         endif
      
          ! Rotate Momentum from global to local coordinate system
          if (IS3DCYLINDRIC) then
            do IEntity = 1, Counters%nEntity
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, Momentum(:, IEntity), Momentum(:, IEntity))
            end do
          end if

           if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1)) then
            ! Rotate water momentum from global to local coordinate system
            if (IS3DCYLINDRIC) then
              do IEntity = 1, Counters%nEntity
                call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, MomentumW(:, IEntity), MomentumW(:, IEntity))
              end do
            end if ! rotation
          end if ! consolidation
      
          if (CalParams%NumberOfPhases==3) then
            if (IS3DCYLINDRIC) then ! rotation is needed
              do IEntity = 1, Counters%nEntity 
                call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, MomentumG(:, IEntity), MomentumG(:, IEntity))
              end do
            end if ! rotation
          end if ! unsat consolidation

         ! DoSystem = .false.  !not needed to map system velocity here
         ! Calculate new nodal velocities from the new momentums
         if (IsMPMComputation()) then !MPM
           call GetNodalVelocityFromNodalMomentumConv(Momentum) ! In local coordinate system
           if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3).or.(.not.(NFORMULATION==1))) then
             call GetNodalWaterVelocityFromNodalWaterMomentum(MomentumW) ! In local coordinate system
             if (CalParams%NumberOfPhases==3) then
               call GetNodalGasVelocityFromNodalGasMomentum(MomentumG)
             end if
           end if
         end if      

          ! Rotate Momentum and TotalVelocitySoil from local to global coordinate system (but not for FEM which is already global)
          if (IS3DCYLINDRIC) then
            do IEntity = 1, Counters%nEntity
              call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, Momentum(:, IEntity), Momentum(:, IEntity))
              call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocitySoil(:, IEntity), TotalVelocitySoil(:, IEntity))
            end do
          end if

           ! Rotate momentum and velocity from local to global coordinate system
           if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1)) then
             if (IS3DCYLINDRIC) then
               do IEntity = 1, Counters%nEntity
                 call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, MomentumW(:, IEntity), MomentumW(:, IEntity))
                 if (IsMPMComputation()) then
                   call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocityWater(:, IEntity), TotalVelocityWater(:, IEntity))
                 end if
               end do
             end if
           end if ! consolidation
          
          ! Rotate momentum and velocity from local to global coordinate system
          if (CalParams%NumberOfPhases==3) then ! rotation is needed
            if (IS3DCYLINDRIC) then ! rotation is needed
              do IEntity = 1, Counters%nEntity
                call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, MomentumG(:, IEntity), MomentumG(:, IEntity))
                call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocityGas(:, IEntity), TotalVelocityGas(:, IEntity))
              end do
            end if
          end if ! Unsat calculation 
          
          ! Overwrite nodal velocity and particle velocity if prescribed
           call ApplyMPPrescribedVelocity(TotalVelocitySys)
           call ApplyNodalPrescribedVelocity(DummyRotationVector1)    
        
           if (CalParams%BoundaryConditions%ApplyInfiltrationRate.or.CalParams%BoundaryConditions%ApplySeepageFace) then  
              call ApplyNodalInfiltrationRate(TotalVelocityWater,TotalVelocitySoil,LumpedNodalPorosityDegSat,NodalUnitMassGradient)
           end if
          
          call GetNodalIncrementalDisplacement(IncrementalDisplacementSoil, TotalVelocitySoil)

          if (IsMPMSkipConvection()) then
            AccumulatedIncDisplacementSoil = AccumulatedIncDisplacementSoil + IncrementalDisplacementSoil
          end if

          if ((.not.(NFORMULATION==1))) then
            call GetNodalIncrementalDisplacement(IncrementalDisplacementWater, TotalVelocityWater) ! Global coordinate system
          
            if (.not.IsMPMComputation()) then !FEM
             do I = 1, Counters%N
               TotalDisplacementWater(I)  = TotalDisplacementWater(I) + IncrementalDisplacementWater(I,1) ! Global coordinate system
             end do 
            end if
          end if

          if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1)) then
            call GetNodalIncrementalDisplacement(IncrementalDisplacementWater, TotalVelocityWater) ! Global coordinate system
     
            if ( CalParams%ApplyAbsorbingBoundary.and.IsMPMComputation()) then
              !MPM ...need to get the accumulated displacement at the nodes
              call GetNodalAccumulatedDisplacementsWater(IncrementalDisplacementWater) ! Global coordinate system
            end if

            if ( CalParams%ApplyAbsorbingBoundary) then
              call CalcParticleDisplacementsWater()
              call UpdateParticleWaterAcceleration()
            end if

            if (.not.IsMPMComputation()) then !FEM
              do I = 1, Counters%N
                TotalDisplacementWater(I)  = TotalDisplacementWater(I) + IncrementalDisplacementWater(I,1) ! Global coordinate system
             end do 
            end if
          
          end if ! consolidation
      
          if (CalParams%NumberOfPhases==3) then
            call GetNodalIncrementalDisplacement(NodalIncDisplacementG, TotalVelocityGas)
     
            if ( CalParams%ApplyAbsorbingBoundary.and. IsMPMComputation()) then !MPM ...need to get the accumulated displacement at the nodes
              call GetNodalAccumulatedDisplacementsGas(NodalIncDisplacementG)
            end if 
                
            if ( CalParams%ApplyAbsorbingBoundary) then 
              call GetIncrementalDisplacementGas(NodalIncDisplacementG)
            end if 
 
            if (.not.IsMPMComputation()) then !FEM
              do I = 1, Counters%N
                TotalDisplacementGas(I)  = TotalDisplacementGas(I) + NodalIncDisplacementG(I,1)
              end do 
            end if
          
          end if ! Unsat Calculation
      
          ! Update the nodal coordinates in case of updated mesh analysis
          if (IsULFEMComputation()) then
            call UpdateNodes(EntityUsed)
            do I = 1, Counters%NEl
              call CoordLocalToGlobal(I, NodalCoordinatesUpd)
            end do
          end if
         
         !call ComputeStructureMoments()
         
         call UpdateNodalTotalDisplacement(EntityUsed)
         
         ! Determine nodal density field
          if(NFORMULATION==1) then ! 1 Constituent
              call DetermineLiquidDensityField()
          else ! 2 Constituents
              if(CalParams%NumberOfPhases==1) then ! 1 Phase
                call DetermineLiquidDensityField() !TEMPORARY
              end if
          end if
    
        if (.not.CalParams%ApplyStrainSmoothing) then 
           ! Update particle strain data
            call UpdateParticleStrains() ! update particle strains for solid and mixture
        end if
        
        if ((.not.CalParams%ApplyStrainSmoothingLiquidTwoLayer)) then 
            call UpdateParticleStrainsLiquidTwoLayer() ! update particle strains for liquid MPs in TwoLayerForm
            if ((.not.(CalParams%NumberOfPhases==1)).and.(.not.(NFORMULATION==1))) then
              call ComputeTwoLayerVolumetricStrainLiquid(IncrementalDisplacementWater, IncrementalDisplacementSoil)
      end if
        end if
          
          
         if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1)) then
           call CalculateWaterVolumetricStrain(IncrementalDisplacementWater)
             if (IsMPMComputation()) then  ! MPM
               !update particle displacement (water)
               call UpdateParticleDisplacementsWater(IncrementalDisplacementWater)
            else !(.not.IsMPMComputation()) then !FEM
              call UpdateParticleVelocityWater(TotalVelocityWater)
             end if
         end if    
         
         if (CalParams%NumberOfPhases==3) then
            call CalculateGasVolumetricStrain(NodalIncDisplacementG)
            if (IsMPMComputation()) then  ! MPM
              !update particle displacement (gas)
              call UpdateParticleDisplacementsGas(NodalIncDisplacementG)
            else ! (.not.IsMPMComputation()) then !FEM
              call UpdateParticleVelocityGas(TotalVelocityGas)
           end if
         end if

         if ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation).and.(NFORMULATION==1)) then
            call CalculateParticleWaterAdvectiveFlux(IncrementalDisplacementWater)
         else if ((CalParams%NumberOfPhases==3).and.(NFORMULATION==1)) then
            call CalculateParticleWaterAdvectiveFlux(IncrementalDisplacementWater)
            call CalculateParticleAirAdvectiveFlux(NodalIncDisplacementG)
         end if

         iOpt = 0

         !Detect Liquid Free Surface
         call DetectLiquidFreeSurface();
         
         if(.not.(NFORMULATION==1)) then
            call TwoLayerData%DetermineTwoLayerStatus() 
         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The balance equations are solved--> Incremental Water Pressure, Gas Pressure and Temperature are calculated, 
!and the new stresses are obtained applying the proper Constitutive equation
!Global coordinate system 
        ! The stresses, the water pressure and the gass pressure are updated (The mass balance equations are solved)
        call MPMDYNGetSig()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Update particle stress data for fully filled elements
       if (IsMPMWithMixedIntegration() .and. .not.IsMPMSkipConvection()) then ! Not needed for pure MPM and FEM
         call AssignStressesToParticles() ! Need to check for parallelization
         call AssignStateParametersToParticles()  
       end if


!\\\\\\\ UPDATE PARTICLE PROPERTIES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      


       ! Update particle weights
      if (IsMPMComputation()) then
       call DynUpdateParticleWeights( )
      end if

      !Update porosity (Directly depends on: DEpsVol)
      if (CalParams%ApplyPorosityUpdate) then !Update porosity (Directly depends on: DEpsVol)
        call DynUpdateParticlePorosity( )
      end if

!\\\\\\\ END UPDATE PARTICLE PROPERTIES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      

       ! VAHID: This does not seem to be true! isMPMSkip..() must be swaped with %SkipConvection
       if (.not.IsMPMSkipConvection()) then
         if (CalParams%ApplyMeshSmoothing) then
           NodalIncDisplacementMeshAdjust = 0.0
           if (.not.CalParams%SkipConvection) then
             if (CalParams%ApplyContactAlgorithm) then ! contact...use pile entity
                 StructureMaterialEntity = getMaterialEntity(CalParams%MovingMesh%MovingMaterialID)
               do I = 1, Counters%N
                 NodalIncDisplacementMeshAdjust(I) = IncrementalDisplacementSoil(I,StructureMaterialEntity) !2 is the entity of the pile
               end do
             else
               ! no contact it is just one entity
               do I = 1, Counters%N
                 NodalIncDisplacementMeshAdjust(I) = IncrementalDisplacementSoil(I,1)
               end do
             endif
           else
             do I = 1, Counters%N
               NodalIncDisplacementMeshAdjust(I) = IncrementalDisplacementSoil(I,1)
             end do
           endif

           call UpdateNodesMeshAdjust(NodalIncDisplacementMeshAdjust)

           ! Perform mesh adjustment by modifying the initial nodal coordinates
           call MeshAdjustment()

           if (IsULFEMComputation()) then
             call UpdateMeshAdjacencyInformation(NodalCoordinatesUpd)
           else
             call UpdateMeshAdjacencyInformation(NodalCoordinates)
           end if
         end if

         ! Update particle global positions, particle-element assignment, particle local positions and shape function values
         ! NOTE: EleParticles and Particle%ElementID get already updated while
         ! the remaining house-keeping data is updated in the succeeding routine!!
         if (IsMPMComputation()) then ! MPM
           call UpdateParticlePos()
         end if

         if ( CalParams%ApplyMeshSmoothing .and. ISAXISYMMETRIC ) then
            ! updating global position of GPs is only needed in case of moving mesh
            call SetGPGlobalPositionElement(GPGlobalPositionElement)
         end if
         
         ! Update particle house-keeping lists and element switches
         if (IsMPMComputation()) then ! MPM
            call UpdateParticleHouseKeeping()
         end if
         ! Reset mesh (nodal coordinates)
         if (.not.IsULFEMComputation()) then
           NodalCoordinatesUpd = NodalCoordinates
         end if

         call SetActiveElement()
         call SetParticleIndex()
         
         if (CalParams%ApplyEmptyElements) then
           ! Check whether elements are deactivated that should be considered fully filled
           call CheckEmptyElements()
           call AdjustParticleDiscretisation()
         end if

         call SetActiveElement()
         call SetParticleIndex()

         call CheckFillingOfElements()

         call SetUpEntityElements ()
         call SetUpMaterialElements ()
       endif

       ! Smoothing of particle stresses for fully filled elements
       if (IsMPMWithMixedIntegration() .and. .not.IsMPMSkipConvection()) then ! Not needed for pure MPM and FEM
         call StressAndPorePressureSmoothening()
         call StateParametersSmoothening ()
       else
         ! Store initial stresses of next step
         call SetInitialStressForNextLoadStep()
       end if
      
   !UPDATE particle properties that depends on stresses (NB: done after stress and pore pressure smoothening)    
       if ((CalParams%NumberOfPhases==3).or. &
          ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation))) then
      
        call DynUpdateParticleDegreeOfSaturation( ) !Update Degree of Saturation (Retention Curve) (Directly depends on: Pg, Pl) (Undirectly depends on: Smin, Smax,T,...)
        call DynUpdateParticleMixedWeight()!Update Mixed Weight according to degree of saturation  
        call DynUpdateParticleHydraulicConductivity( )!Update hydraulic condictivity according to degree of saturation
       end if
          


       if ((.not.(NFORMULATION==1)).and.(.not.(CalParams%NumberOfPhases==1))) then
         call TwoLayerData%DetermineContainedMaterialTypes()
         call DetermineDensityField(TwoLayerData%Nodes%DensityLiquidL, MaterialPointTypeLiquid, .false., .true.)
         call DetermineDensityField(TwoLayerData%Nodes%DensitySolidL, MaterialPointTypeSolid, .false., .true.)
         call DetermineDensityField(TwoLayerData%Nodes%DensityLiquidS, MaterialPointTypeLiquid, .false., .false.)
         call DetermineDensityField(TwoLayerData%Nodes%DensitySolidS, MaterialPointTypeSolid, .false., .false.)


         call DetermineConcentrationRatio(TwoLayerData%Nodes%DensitySolidL, TwoLayerData%Nodes%DensityLiquidL, MaterialPointTypeLiquid,.true.)
         call DetermineConcentrationRatio(TwoLayerData%Nodes%DensitySolidL, TwoLayerData%Nodes%DensityLiquidL, MaterialPointTypeSolid,.true.)
         call DetermineConcentrationRatio(TwoLayerData%Nodes%DensitySolidS, TwoLayerData%Nodes%DensityLiquidS, MaterialPointTypeSolid,.false.)

         call TwoLayerData%DetermineConcentrationRatioElementFND()

         call DetermineBoundaryElementSolidDomain2LayerForm()
         call DetermineBoundaryElementLiquidDomain2LayerForm()

         call DetermineFillingRatioField()
         call DetermineFillingRatioLiquid(TwoLayerData)

       end if

       ! calculate system initial kinetic energy
       if (CalParams%ConvergenceCheck%KineticEnergy0 == -1.0) then
         if (IsMPMComputation()) then ! MPM
           if (NFORMULATION==1) then
             call CalculateKineticEnergy(TotalVelocitySoilPrevious, TotalVelocityWaterPrevious)
           else
             call CalculateKineticEnergy2LayForm(TotalVelocitySoilPrevious, TotalVelocityWaterPrevious)
           end if
         else ! FEM
           call CalculateKineticEnergyFEM(TotalVelocitySoilPrevious, TotalVelocityWaterPrevious)
         end if

         CalParams%ConvergenceCheck%KineticEnergy0 = CalParams%ConvergenceCheck%KineticEnergy
         CalParams%ConvergenceCheck%KineticEnergySoil0 = CalParams%ConvergenceCheck%KineticEnergySoil
         CalParams%ConvergenceCheck%KineticEnergyWater0 = CalParams%ConvergenceCheck%KineticEnergyWater
       end if

       ! calculate system kinetic energy
       if (.not.IsMPMSkipConvection()) then
         if (IsMPMComputation()) then
           ! MPM
           if (NFORMULATION==1) then
             call CalculateKineticEnergy(TotalVelocitySoil, TotalVelocityWater)
           else
             call CalculateKineticEnergy2LayForm(TotalVelocitySoil, TotalVelocityWater)
           end if
         else
           ! FEM
           call CalculateKineticEnergyFEM(TotalVelocitySoil, TotalVelocityWater)
         end if
         ! calculate system internal and external works
       endif

       if (NFORMULATION==1) then
         call CalculateIntAndExtWorks(IncrementalDisplacementWater)
       else
         call CalculateIntAndExtWorks2LayForm(IncrementalDisplacementWater)
       end if
      
       if ((CalParams%ApplyFixedSolidSkeleton == .false.).and.(CalParams%PrescribedHead%HydraulicHead == .true.)) then
           HydraulicHeadLoadedElemID = .false.
           call GetAndUpdateHydraulicHeadLoad(HydraulicHeadLoadedElemID, LoadedSides, HydraulicHeadNodesConnectivities, HydraulicHeadLoad)
           Counters%HydraulicHeadSides = LoadedSides
           HydraulicHeadVector = 0.0
           if (NDIM == 3) then
               do I = 1, Counters%HydraulicHeadSides
                   call Load3D(HydraulicHeadVector, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, HydraulicHeadNodesConnectivities(1, I), HydraulicHeadLoad(1, 1, I), LOADTYPE_LIQUID)
               end do
           elseif (NDIM == 2) then
               do I = 1, Counters%HydraulicHeadSides
                   call Load2D(HydraulicHeadVector, ReducedDof, NodalCoordinates, 1, HydraulicHeadNodesConnectivities(1, I), HydraulicHeadLoad(1, 1, I), LOADTYPE_LIQUID)
               end do
           end if
           HydraulicHeadLoadTotal = 0.0
           do IEntity=1,Counters%NEntity
               do IDof=1,Counters%N
                 HydraulicHeadLoadTotal(IDof, IEntity) = HydraulicHeadVector(IDof)
               end do
            end do
       end if
       
       if (CalParams%SkipConvection) then
         if (.not.IsMPMSkipConvection()) then
           AccumulatedIncDisplacementSoil = 0.0
           call setElementDeterminant()
           call setIsDistorted(.false.)
         else
           call setIsDistorted()
         endif
       end if

       end subroutine DYNConvectivePhase

       subroutine ComputeForceError()
        !**********************************************************************
        !
        !    Function:  Computes the out-of-balance forces error, for the
        !               mixture and for the water phase.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
       implicit none

         ! Local variables
         integer(INTEGER_TYPE) :: I, J, K, INode, IDof
         real(REAL_TYPE), dimension(Counters%N) :: ExternalLoadDueToVelocity
         real(REAL_TYPE) :: OutOfBalanceNorm, ExternalLoadNorm
         real(REAL_TYPE) :: OutOfBalanceNormWater, ExternalLoadNormWater
         logical :: IsPrescribedVelocity, IsPrescribedSurfaceVelocity

         ! Compute ForceError of the mixture

          ExternalLoadDueToVelocity = 0.0
          TemporaryMappingVector = 0.0

           IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do

          if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity
            do INode = 1, Counters%Nodtot
              do J = 1, NVECTOR
                IDof = (INode-1) * NDOFL + J      ! global Dof
                do K = 1, Counters%NEntity
                  ExternalLoadDueToVelocity(IDof) = ExternalLoadDueToVelocity(IDof) + FReaction(IDoF, K)
                end do
              end do
            end do 
         end if
         
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
             if((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then    ! consolidation
               if(CalParams%ApplyAbsorbingBoundary) then  ! absorbing boundary
                 TemporaryMappingVector(I) = TemporaryMappingVector(I)+  &
                                         GravityLoadMixture(I, J) * PBoundary(I) + &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         DummyIntLoadWaterPrevious(I, J) * PBoundaryWater(I) &
                                         - VisDampForceSld(I, J) &
                                         - VisDampForceWat(I, J)
               else
                 TemporaryMappingVector(I) = TemporaryMappingVector(I)+   & ! no absorbing boundary
                                         GravityLoadMixture(I, J) * PBoundary(I) + &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         DummyIntLoadWaterPrevious(I, J) * PBoundaryWater(I)
               end if
             else     ! no consolidation
               if(CalParams%ApplyAbsorbingBoundary) then ! absorbing boundary
                 TemporaryMappingVector(I) = TemporaryMappingVector(I)+  &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I) &
                                         - VisDampForceSld(I, J)
               else
                 TemporaryMappingVector(I) = TemporaryMappingVector(I) + &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I)
               end if
             end if
           end do
         end do
         
         
         call Norm(Counters%NodTot, NDOFL, TemporaryMappingVector, ReducedDof, OutOfBalanceNorm)
     
         TemporaryMappingVector = 0.0
          do J = 1, Counters%NEntity
             do I = 1, Counters%N
             if((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoadMixture(I, J) * PBoundary(I) + &
                                         ExtLoad(I, J) * PBoundary(I)
             else
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I)
             end if
           end do
         end do

           IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do

         if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity

         TemporaryMappingVector = 0.0
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) + ExternalLoadDueToVelocity (I)
           end do
         end do
        end if
    
         call Norm(Counters%NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof, &
                   ExternalLoadNorm)
     
        if (ExternalLoadNorm==0.0) then
          CalParams%ConvergenceCheck%ForceError = CONVERGENCE_ERROR_NOT_USED
        else
          CalParams%ConvergenceCheck%ForceError = OutOfBalanceNorm / ExternalLoadNorm
        end if       
        
        if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then 


        CalParams%ConvergenceCheck%ForceErrorSoil = CONVERGENCE_ERROR_NOT_USED
       
            
         ! Compute ForceError of the water   
            
         TemporaryMappingVector = 0.0
          
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
             if(CalParams%ApplyAbsorbingBoundary) then  ! absorbing boundary
               TemporaryMappingVector(I) = TemporaryMappingVector(I)+  &
                                         GravityLoadWater(I, J) * PBoundaryWater(I) + &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         IntLoadWaterPrevious(I, J) * PBoundaryWater(I) + &
                                         QVW(I,J) - &
                                         VisDampForceWat(I, J)
             else
                TemporaryMappingVector(I) = TemporaryMappingVector(I)+   & ! no absorbing boundary 
                                         GravityLoadWater(I, J) * PBoundaryWater(I)+ &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         IntLoadWaterPrevious(I, J) * PBoundaryWater(I) + &
                                         QVW(I,J)
             end if
           end do
         end do          
         
         call Norm(Counters%NodTot, NDOFL, TemporaryMappingVector, ReducedDof, OutOfBalanceNormWater)
     
         TemporaryMappingVector = 0.0
         
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoadWater(I, J) * PBoundaryWater(I) + &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) +   &
                                         QVW(I,J)
 
           end do
         end do

         call Norm(Counters%NodTot, NDOFL, TemporaryMappingVector, ReducedDof, ExternalLoadNormWater)
     
        if (ExternalLoadNormWater==0.0) then
          CalParams%ConvergenceCheck%ForceErrorWater = CONVERGENCE_ERROR_NOT_USED
        else
          CalParams%ConvergenceCheck%ForceErrorWater = OutOfBalanceNormWater / ExternalLoadNormWater
        end if       
        
        end if
        
        end subroutine ComputeForceError
        
        
          subroutine ComputeForceError2LayForm()
        !**********************************************************************
        !
        !    Function:  For 2 Layer Formulation :: Computes the out-of-balance forces error, 
        !               for the solid and for the water phase.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
       implicit none

         ! Local variables
         integer(INTEGER_TYPE) :: I, J, K, INode, IDof
         real(REAL_TYPE) :: OutOfBalanceNormSolid, ExternalLoadNormSolid
         real(REAL_TYPE), dimension(Counters%N) :: ExternalLoadDueToVelocity
         real(REAL_TYPE) :: OutOfBalanceNormWater, ExternalLoadNormWater
         logical :: IsPrescribedVelocity, IsPrescribedSurfaceVelocity
 
        ! Compute ForceError of the SOLID
                  
          ExternalLoadDueToVelocity = 0.0
          TemporaryMappingVector = 0.0
          
           IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do
           
           if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity

            do INode = 1, Counters%Nodtot
              do J = 1, NVECTOR
                IDof = (INode-1) * NDOFL      ! global Dof
                do K = 1, Counters%NEntity
                  ExternalLoadDueToVelocity(IDof) = ExternalLoadDueToVelocity(IDof) + FReaction(IDoF, K)
                end do
              end do
            end do 
          end if 
         
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
                 TemporaryMappingVector(I) = TemporaryMappingVector(I) + &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I) - QVW(I, J) + &
                                         TwoLayerData%InteractionForceSolid(I, J)
           end do
         end do        
         
         call Norm(Counters%NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof,  &
                   OutOfBalanceNormSolid)
     
         TemporaryMappingVector = 0.0
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I) - QVW(I, J) +  &
                                         TwoLayerData%InteractionForceSolid(I, J)
           end do
         end do

            IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do

          if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity
         TemporaryMappingVector = 0.0
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                    ExternalLoadDueToVelocity (I)
           end do
         end do
        end if
    
         call Norm(Counters%NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof, &
                   ExternalLoadNormSolid)
     
        if (ExternalLoadNormSolid==0.0) then
          CalParams%ConvergenceCheck%ForceErrorSoil = CONVERGENCE_ERROR_NOT_USED
        else
          CalParams%ConvergenceCheck%ForceErrorSoil = OutOfBalanceNormSolid / ExternalLoadNormSolid
        end if       
       
        ! Compute ForceError of the LIQUID  

         TemporaryMappingVector = 0.0
          
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
                TemporaryMappingVector(I) = TemporaryMappingVector(I)+   & ! no absorbing boundary 
                                         GravityLoadWater(I, J) * PBoundaryWater(I)+ &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         IntLoadWaterPrevious(I, J) * PBoundaryWater(I) + &
                                         QVW(I,J) +  &
                                         TwoLayerData%InteractionForceLiquid(I, J)
           end do
         end do       
         
         call Norm(Counters%NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof,  &
                   OutOfBalanceNormWater)
     
         TemporaryMappingVector = 0.0
         
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoadWater(I, J) * PBoundaryWater(I) + &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) +   &
                                         QVW(I,J) +  &
                                         TwoLayerData%InteractionForceLiquid(I, J)
 
           end do
         end do        

         call Norm(Counters%NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof, &
                   ExternalLoadNormWater)
     
        if (ExternalLoadNormWater==0.0) then
          CalParams%ConvergenceCheck%ForceErrorWater = CONVERGENCE_ERROR_NOT_USED
        else
          CalParams%ConvergenceCheck%ForceErrorWater = OutOfBalanceNormWater / ExternalLoadNormWater
        end if       
        
        
       end subroutine ComputeForceError2LayForm
     
       subroutine UpdateParticleVelocityAndMapMomentum(Momentum)
                              
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum       !CC added nEntity
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: IDof
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleIncrementalVelocity
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: NodAcc
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) :: Time, PartilceMass
          
          !!CC - changed to function call - needs to be set to zero here
          Momentum = 0.0 
          Time = CalParams%TimeIncrement
          NoEn = Counters%nEntity

          do IAEl = 1, Counters%NAEl   ! Loop over all elements
            IEl = ActiveElement(IAEl)
            do I = 1, NVECTOR
              IDof(I, 1:ELEMENTNODES) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
              NodAcc(I, 1:ELEMENTNODES, 1:Counters%nEntity) = AccelerationSoil(IDof(I, 1:ELEMENTNODES), 1:NoEn)
            end do
            
            do IPart = 1, NPartEle(IEl)   ! Loop over all particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 ParticleVelocity = VelocityArray(ParticleIndex,:)
                 PartilceMass = MassArray(ParticleIndex)
                 ParticleShape = ShapeValuesArray(ParticleIndex,:)
                 ParticleIncrementalVelocity = 0.0
                 ParticleAcceleration = 0.0
                 if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex) 
                 else
                  iEntity = 1
                 end if 
                  
                 do INode = 1, ELEMENTNODES  ! loop over element nodes
                      
                  do I = 1, NVECTOR  
                    ! Particle accelerations
                    ParticleAcceleration(I) = ParticleAcceleration(I) + ParticleShape(INode) * NodAcc(I, INode, iEntity)
                    ! Particle x-velocity
                    ParticleIncrementalVelocity(I) = ParticleIncrementalVelocity(I) + Time * ParticleShape(INode) * NodAcc(I, INode, iEntity)
                  end do    
                end do !Loop over nodes
               
                ParticleVelocity = ParticleVelocity + ParticleIncrementalVelocity

                do I = 1, NVECTOR ! nodal i-momentum
                  Momentum(IDof(I,1:ELEMENTNODES), iEntity) = Momentum(IDof(I,1:ELEMENTNODES), iEntity) + PartilceMass * ParticleShape * ParticleVelocity(I)
                end do     

                VelocityArray(ParticleIndex,:) = ParticleVelocity
                AccelerationArray(ParticleIndex,:) =  ParticleAcceleration
                 
              end if ! NumbOfLayers = 1 or SOLID MatPoint
            end do !Loop over particles
          end do !elements    
          
         end subroutine UpdateParticleVelocityAndMapMomentum  
         
         subroutine GetNodalIncrementalDisplacement(NodalIncDisplacement, NodalVelocities)
                                                   
        !**********************************************************************
        !
        !    Function:  To calculate the nodal incremental displacements from the nodal velocities
        !
        !    NodalIncDisplacement: The nodal incremental displacements vector
        !    NodalVelocities : Nodal velocities vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: NodalIncDisplacement
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: NodalVelocities
     
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          
          do IDOF = 1, Counters%N       !Loop over all degrees of freedom
            do J = 1,Counters%nEntity   !loop over all entities
              NodalIncDisplacement(IDOF,J) = NodalVelocities(IDOF,J) * CalParams%TimeIncrement
            end do
          end do
                   
         end subroutine GetNodalIncrementalDisplacement    
         
         subroutine UpdateNodalTotalDisplacement(EntityUsed)
                                                   
        !**********************************************************************
        !
        !    Function:  To add the nodal incremental displacements to the total displacemnts
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, EntityUsed
          
          do IDOF = 1, Counters%N    !Loop over all degrees of freedom
            TotalDisplacementSoil(IDOF)= TotalDisplacementSoil(IDOF)+ IncrementalDisplacementSoil(IDOF,EntityUsed)
          end do
          
          if(.not.(NFORMULATION==1)) then
          do IDOF = 1, Counters%N    !Loop over all degrees of freedom
            TotalDisplacementWater(IDOF)= TotalDisplacementWater(IDOF)+ IncrementalDisplacementWater(IDOF,EntityUsed)
          end do
          end if
          
        end subroutine UpdateNodalTotalDisplacement    
        
         
        subroutine CalculateKineticEnergy(TotalVelocitySoilLoc, TotalVelocityWaterLoc)
        !**********************************************************************
        !
        !   Function:  Calculate system kinetic energy
        !
        !   Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: TotalVelocitySoilLoc, TotalVelocityWaterLoc
          ! local variables
          integer(INTEGER_TYPE) :: I, INode, IDof
          real(REAL_TYPE) :: SquaredVelocityS, KES, KESn, SquaredVelocityW, KEW, KEWn
          
          ! Calculation of kinetic energy of the system using the velocity fields given by nodal values
          KES = 0.0
          KEW = 0.0  
          KESn = 0.0
          KEWn = 0.0
          
          do INode = 1, Counters%NodTot ! loop over all nodes
          
            IDof = ReducedDof(INode)

            SquaredVelocityS = 0.0
            do I = 1, NVECTOR
              SquaredVelocityS = SquaredVelocityS + (TotalVelocitySoilLoc(IDof+I, 1)) * (TotalVelocitySoilLoc(IDof+I, 1))
            end do
            KESn = 0.5 * LumpedMassDry(IDof+1, 1) * SquaredVelocityS
            
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then ! if consolidation
              SquaredVelocityW = 0.0
              do I = 1, NVECTOR
                SquaredVelocityW = SquaredVelocityW + (TotalVelocityWaterLoc(IDof+I, 1)) * (TotalVelocityWaterLoc(IDof+I, 1))
              end do
              KEWn = 0.5 * LumpedMassNWater(IDof+1, 1) * SquaredVelocityW
            end if ! end if consolidation
            KES = KES + KESn
            KEW = KEW + KEWn

          end do ! end loop over nodes
          
          CalParams%ConvergenceCheck%KineticEnergy = KES + KEW
          CalParams%ConvergenceCheck%KineticEnergySoil = KES
          CalParams%ConvergenceCheck%KineticEnergyWater = KEW
          
        end subroutine CalculateKineticEnergy
        
        
        subroutine CalculateKineticEnergy2LayForm(TotalVelocitySoilLoc, TotalVelocityWaterLoc)
        !**********************************************************************
        !
        !   Function:  Calculate system kinetic energy in 2 Layer formulation
        !
        !   Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: TotalVelocitySoilLoc, TotalVelocityWaterLoc    
          ! local variables
          integer(INTEGER_TYPE) :: I, INode, IDof
          real(REAL_TYPE) :: SquaredVelocityS, KES, KESn, SquaredVelocityW, KEW, KEWn
          
          ! Calculation of kinetic energy of the system
          KES = 0.0
          KEW = 0.0 
          KESn = 0.0
          KEWn = 0.0
          
          do INode = 1, Counters%NodTot ! loop over all nodes
          
            IDof = ReducedDof(INode)
            
            ! Kinetic energy of solid phase
            SquaredVelocityS = 0.0
            do I = 1, NVECTOR
              SquaredVelocityS = SquaredVelocityS + (TotalVelocitySoilLoc(IDof+I, 1)) * (TotalVelocitySoilLoc(IDof+I, 1))
            end do
            KESn = 0.5 * LumpedMassDry(IDof+1, 1) * SquaredVelocityS
            
            ! Kinetic energy of water phase
            SquaredVelocityW = 0.0
            do I = 1, NVECTOR
                SquaredVelocityW = SquaredVelocityW + (TotalVelocityWaterLoc(IDof+I, 1)) * (TotalVelocityWaterLoc(IDof+I, 1))
            end do
            KEWn = 0.5 * LumpedMassWater(IDof+1, 1) * SquaredVelocityW
            KES = KES + KESn
            KEW = KEW + KEWn

          end do ! end loop over nodes

          
          CalParams%ConvergenceCheck%KineticEnergy = KES + KEW       
          CalParams%ConvergenceCheck%KineticEnergySoil = KES          
          CalParams%ConvergenceCheck%KineticEnergyWater = KEW
            
        end subroutine CalculateKineticEnergy2LayForm
        
        
        subroutine CalculateKineticEnergyFEM(TotalVelocitySoilLoc, TotalVelocityWaterLoc)
        !**********************************************************************
        !
        !   Function:  Calculate system kinetic energy in FEM
        !
        !**********************************************************************
        implicit none 
        
          ! arguments
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: TotalVelocitySoilLoc, TotalVelocityWaterLoc
          
          !local variable
          integer(INTEGER_TYPE) :: I, INode, IDof
          real(REAL_TYPE) :: SquaredVelocityS, KES, SquaredVelocityW, KEW
          
          KES = 0.0
          KEW = 0.0          
          
          ! Calculation of kinetic energy of the system

          do INode = 1, Counters%NodTot ! loop over all nodes
              
            IDof = ReducedDof(INode)
            
            SquaredVelocityS = 0.0
            do I = 1, NVECTOR
              SquaredVelocityS = SquaredVelocityS + (TotalVelocitySoilLoc(IDof+I, 1)) * (TotalVelocitySoilLoc(IDof+I, 1))
            end do
            KES  = KES + 0.5 * LumpedMassDry(IDof+1, 1) * SquaredVelocityS
          
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then ! if consolidation
              SquaredVelocityW = 0.0
              do I = 1, NVECTOR
                SquaredVelocityW = SquaredVelocityW + (TotalVelocityWaterLoc(IDof+I, 1)) * (TotalVelocityWaterLoc(IDof+I, 1))
              end do
              KEW  = KEW + 0.5 * LumpedMassNWater(IDof+1, 1) * SquaredVelocityW
            end if ! end if consolidation
          
          end do ! end loop over nodes
         
          CalParams%ConvergenceCheck%KineticEnergy = KES + KEW       
          CalParams%ConvergenceCheck%KineticEnergySoil = KES          
          CalParams%ConvergenceCheck%KineticEnergyWater = KEW
        
        end subroutine CalculateKineticEnergyFEM

        subroutine CalculateIntAndExtWorks(IncrementalDisplacementWater)
        !**********************************************************************
        !
        !    Function:  Calculate internal and external works
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IncrementalDisplacementWater
          
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, I, J, K, INode
          real(REAL_TYPE), dimension(Counters%N) :: ExternalLoadDueToVelocity
          real(REAL_TYPE) :: IncrementalExtWorkSoil
          real(REAL_TYPE) :: IncrementalIntWorkSoil
          real(REAL_TYPE) :: IncrementalExtWorkWater
          real(REAL_TYPE) :: IncrementalIntWorkWater
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: DummyIntLoadWater!, DummyIntLoadWaterPrevious
          logical :: IsPrescribedVelocity, IsPrescribedSurfaceVelocity

          if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then ! if consolidation

            IntLoadPrevious = IntLoad
            IntLoadWaterPrevious = IntLoadWater
            IntLoadWaterPorosityPrevious = IntLoadWaterPorosity
            QVWPorosityPrevious = QVWPorosity
            
            if(CalParams%ApplyPartialSaturation) then
             DummyIntLoadWaterPrevious =  BishopIntLoad  !account for partial saturation to be used in mixture equation
            else
             DummyIntLoadWaterPrevious = IntLoadWater
            end if
            
            call GetNodalIntForces()
            call ConsolidationIntForces(IntLoadWater)
            call ConsolidationIntForcesPorosity(IntLoadWaterPorosity)
            call GetQVWArrayPorosity(QVWPorosity)
            
            if((CalParams%NumberOfPhases==3).or.(CalParams%ApplyPartialSaturation)) then
              call ConsolidationForcesBishop(BishopIntLoad)  
              DummyIntLoadWater =  BishopIntLoad  !account for partial saturation
            else
              DummyIntLoadWater = IntLoadWater
            end if
              
            IncrementalIntWorkSoil = 0.0
            IncrementalExtWorkSoil = 0.0
            IncrementalIntWorkWater = 0.0
            IncrementalExtWorkWater = 0.0
            
            do IDOF = 1, Counters%N ! Loop over all degrees of freedom
                do J = 1, Counters%nEntity
                    
                    ! calculate incremental internal and external work of soil
                    IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                             0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J) + DummyIntLoadWater(IDOF,J) + &
                                             DummyIntLoadWaterPrevious(IDOF,J) - IntLoadWaterPorosity(IDOF,J) - &
                                             IntLoadWaterPorosityPrevious(IDOF,J)) * &
                                             IncrementalDisplacementSoil(IDOF,J)                        ! mid point rule
                    IncrementalExtWorkSoil = IncrementalExtWorkSoil + &
                                             (ExtLoad(IDOF,J) - ExtLoadWaterPorosity(IDOF,J) + &
                                             GravityLoadMixture(IDOF,J) - GravityLoadWaterPorosity(IDOF,J) - &
                                             0.5 * (QVWPorosity(IDOF,J) + QVWPorosityPrevious(IDOF,J))) *  &
                                             IncrementalDisplacementSoil(IDOF,J)                        ! left/mid point rule
                    
                    ! calculate incremental internal and external work of water
                    IncrementalIntWorkWater = IncrementalIntWorkWater +  &
                                              0.5 * (IntLoadWaterPorosity(IDOF,J) + IntLoadWaterPorosityPrevious(IDOF,J)) *   &
                                              IncrementalDisplacementWater(IDOF,J)                      ! mid point rule
                    IncrementalExtWorkWater = IncrementalExtWorkWater +  &
                                              (ExtLoadWaterPorosity(IDOF,J) + &
                                              GravityLoadWaterPorosity(IDOF,J) + &
                                              0.5 * (QVWPorosity(IDOF,J) + QVWPorosityPrevious(IDOF,J))) *  &
                                              IncrementalDisplacementWater(IDOF,J)                      ! left/mid point rule
                end do
            end do
          
            CalParams%ConvergenceCheck%InternalWorkSoil = CalParams%ConvergenceCheck%InternalWorkSoil + IncrementalIntWorkSoil
            CalParams%ConvergenceCheck%ExternalWorkSoil = CalParams%ConvergenceCheck%ExternalWorkSoil + IncrementalExtWorkSoil
            
            CalParams%ConvergenceCheck%InternalWorkWater = CalParams%ConvergenceCheck%InternalWorkWater + IncrementalIntWorkWater 
            CalParams%ConvergenceCheck%ExternalWorkWater = CalParams%ConvergenceCheck%ExternalWorkWater + IncrementalExtWorkWater
            
            CalParams%ConvergenceCheck%InternalWork = CalParams%ConvergenceCheck%InternalWorkSoil + CalParams%ConvergenceCheck%InternalWorkWater
            CalParams%ConvergenceCheck%ExternalWork = CalParams%ConvergenceCheck%ExternalWorkSoil + CalParams%ConvergenceCheck%ExternalWorkWater
           
          else ! 1-phase calculation
          
            IntLoadPrevious = IntLoad
            call GetNodalIntForces()  
              
            IncrementalIntWorkSoil = 0.0
            IncrementalExtWorkSoil = 0.0
            ExternalLoadDueToVelocity = 0.0
          
            do INode = 1, Counters%Nodtot
              do J = 1, NVECTOR
                IDof = (INode-1) * NDOFL + J      ! global Dof
                do K = 1, Counters%NEntity
                    ExternalLoadDueToVelocity(IDof) = ExternalLoadDueToVelocity(IDof) + FReaction(IDoF, K)
                end do
              end do
            end do 

           IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do

           if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity 
              do IDOF = 1, Counters%N ! Loop over all degrees of freedom
                do J = 1, Counters%nEntity   
                  IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                           0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J)) *  &
                                           IncrementalDisplacementSoil(IDOF,J)                  ! mid point rule
                  IncrementalExtWorkSoil = IncrementalExtWorkSoil +  &
                                           ExternalLoadDueToVelocity(IDOF) *  &
                                           TotalVelocitySoil(IDOF,J) * CalParams%TimeIncrement  ! left point rule
                end do
              end do
            else
            ! load
              do IDOF = 1, Counters%N ! Loop over all degrees of freedom
                do J = 1, Counters%nEntity
                  IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                           0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J)) *  &
                                           IncrementalDisplacementSoil(IDOF,J)                  ! mid point rule
                  IncrementalExtWorkSoil = IncrementalExtWorkSoil +  &
                                           (ExtLoad(IDOF,J) + GravityLoad(IDOF,J)) *  &
                                           IncrementalDisplacementSoil(IDOF,J)                  ! left point rule
                end do
              end do
            endif

            CalParams%ConvergenceCheck%InternalWork = CalParams%ConvergenceCheck%InternalWork + IncrementalIntWorkSoil 
            CalParams%ConvergenceCheck%ExternalWork = CalParams%ConvergenceCheck%ExternalWork + IncrementalExtWorkSoil

          end if
          
        end subroutine CalculateIntAndExtWorks
        
        subroutine CalculateIntAndExtWorks2LayForm(IncrementalDisplacementWater)
        !**********************************************************************
        !
        !    Function:  Calculate internal and external works
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none       
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IncrementalDisplacementWater
          
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, I, J, K, INode
          real(REAL_TYPE), dimension(Counters%N) :: ExternalLoadDueToVelocity
          real(REAL_TYPE) :: IncrementalExtWorkSoil
          real(REAL_TYPE) :: IncrementalIntWorkSoil
          real(REAL_TYPE) :: IncrementalExtWorkWater
          real(REAL_TYPE) :: IncrementalIntWorkWater
          logical :: IsPrescribedVelocity, IsPrescribedSurfaceVelocity
   
          ! Calculate Internal and External Work of the SOLID and WATER

          IntLoadPrevious = IntLoad
          IntLoadWaterPrevious = IntLoadWater
            
          call GetNodalIntForces()
          call ConsolidationIntForces(IntLoadWater)
            
          IncrementalIntWorkSoil = 0.0
          IncrementalExtWorkSoil = 0.0
          IncrementalIntWorkWater = 0.0
          IncrementalExtWorkWater = 0.0
          ExternalLoadDueToVelocity = 0.0
          
          do INode = 1, Counters%Nodtot
            do J = 1, NVECTOR
              IDof = ReducedDof(INode)+J!(INode-1) * NDOFL      ! global Dof
              do K = 1, Counters%NEntity
                ExternalLoadDueToVelocity(IDof) =  &
                  ExternalLoadDueToVelocity(IDof) +  &
                  FReaction(IDoF, K)
              end do
            end do
          end do 
          
           do IDOF = 1, Counters%N ! Loop over all degrees of freedom
            do J = 1, Counters%nEntity   
 
              IsPrescribedVelocity = .false.
              IsPrescribedSurfaceVelocity = .false.
              do I = 1, NVECTOR
                IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
                IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
              end do

              if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity 
            
                IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                         0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J)) *  &
                                         IncrementalDisplacementSoil(IDOF,J)                  ! mid point rule
                IncrementalExtWorkSoil = IncrementalExtWorkSoil +  &
                                         ExternalLoadDueToVelocity(IDOF) *  &
                                         TotalVelocitySoil(IDOF,J) * CalParams%TimeIncrement

                IncrementalIntWorkWater = IncrementalIntWorkWater +  &
                                          0.5 * (IntLoadWater(IDOF,J) + IntLoadWaterPrevious(IDOF,J)) *  &
                                          IncrementalDisplacementWater(IDOF,J)                      ! mid point rule
                IncrementalExtWorkWater = IncrementalExtWorkWater +  &
                                          (ExtLoadWater(IDOF,J) + GravityLoadWater(IDOF,J) +  &
                                          QVW(IDOF,J) + TwoLayerData%InteractionForceLiquid(IDOF, J) ) *  &
                                          IncrementalDisplacementWater(IDOF,J)                      ! left point rule
     
             else ! load
                
                IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                         0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J)) *  &
                                         IncrementalDisplacementSoil(IDOF,J)                        ! mid point rule
                IncrementalExtWorkSoil = IncrementalExtWorkSoil +  &
                                         (ExtLoad(IDOF,J) + GravityLoad(IDOF,J) - &
                                         QVW(IDOF, J) + TwoLayerData%InteractionForceSolid(IDOF, J)) *  &
                                         IncrementalDisplacementSoil(IDOF,J)                        ! left point rule
                    
                IncrementalIntWorkWater = IncrementalIntWorkWater +  &
                                          0.5 * (IntLoadWater(IDOF,J) + IntLoadWaterPrevious(IDOF,J)) *  &
                                          IncrementalDisplacementWater(IDOF,J)                      ! mid point rule 
                IncrementalExtWorkWater = IncrementalExtWorkWater +  &
                                          (ExtLoadWater(IDOF,J) + GravityLoadWater(IDOF,J) +  &
                                          QVW(IDOF, J) + TwoLayerData%InteractionForceLiquid(IDOF, J)) *  &
                                          IncrementalDisplacementWater(IDOF,J)                      ! left point rule
     
             end if
            end do   
          end do
          
          CalParams%ConvergenceCheck%InternalWorkSoil = CalParams%ConvergenceCheck%InternalWorkSoil + IncrementalIntWorkSoil 
          CalParams%ConvergenceCheck%ExternalWorkSoil = CalParams%ConvergenceCheck%ExternalWorkSoil + IncrementalExtWorkSoil

          CalParams%ConvergenceCheck%InternalWorkWater = CalParams%ConvergenceCheck%InternalWorkWater + IncrementalIntWorkWater 
          CalParams%ConvergenceCheck%ExternalWorkWater = CalParams%ConvergenceCheck%ExternalWorkWater + IncrementalExtWorkWater
            
        end subroutine CalculateIntAndExtWorks2LayForm

        subroutine  GetAndUpdateHydraulicHeadLoad(HHElemID, NumberOfSides, HHNodesCon, HHDummyLoad)
        !HH= HydraulicHead  
 
        use ModMeshAdjacencies
        use ModReadCalculationData
        use ModMeshInfo
                
        logical, dimension(:), allocatable, intent(out) :: HHElemID
        integer(INTEGER_TYPE), dimension(:,:), allocatable, intent(out) :: HHNodesCon
        integer(INTEGER_TYPE),intent(out) ::NumberOfSides
        real(REAL_TYPE), dimension(:, :, :), allocatable,intent(out)  :: HHDummyLoad
        !local
        integer(INTEGER_TYPE) :: IElement, IError,ISide, IsBoundarySide
        integer(INTEGER_TYPE) :: NodeID, J, I, K
            
        
        integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: LocalNodes
        logical, dimension (ELEMENTBOUNDARYNODES) :: NodePair
        logical, dimension (Counters%NEl, ELEMENTSIDES) :: HHElemSide
           
       ! initialization
        NumberOfSides = 0
        allocate (HHElemID(Counters%NEl), stat = IError)
        HHElemID = .false.
        
        HHElemSide = .false.
       
        do IElement = 1, Counters%NEl ! loop all elements
            if (IsActiveElement(IElement)) then
                do ISide =1,ELEMENTSIDES !loop over side nodes
                    IsBoundarySide = BoundaryElementSurface(IElement,ISide,IsActiveElement, Counters%NEl) !Give 1 if ISide is adiacent to inactive element
                    if (IsBoundarySide==1) then !side is on a free surface
                        call DetermineSideNodes(ISide,LocalNodes) !Local ID (1, 2, 3...) of nodes at the boundary of ISide, TRI3,TRI6 = 2 nodes, TETRA = 3 nodes
                        do J=1,ELEMENTBOUNDARYNODES

                            NodeID = ElementConnectivities(LocalNodes(J),IElement) !Global name of boundary node
                            call DetectHydraulicHeadNode(IsHydraulicHeadNode)
                            NodePair(J)= IsHydraulicHeadNode(NodeID)
                        end do
                        
                        if (count(NodePair) == ELEMENTBOUNDARYNODES) then !If both nodes of the ISide are inside infiltration area then store name of element and side

                            NumberOfSides = NumberOfSides + 1      !count sides
                            HHElemID(IElement) = .true.
                            HHElemSide(IElement,ISide) = .true. !matrix which store loaded side for the elements
                        end if
                    end if
                end do
            end if
        end do
        
        NumberOfSides = count(HHElemSide)
              
        if (allocated(HHNodesCon)) then
            deallocate(HHNodesCon, stat = IError)
        end if
                
        if (NDIM == 3) then
            allocate(HHNodesCon(N_BOUNDARY_NODES_HOE,NumberOfSides), stat = IError)
        elseif (NDIM == 2) then
            allocate(HHNodesCon(ELEMENTBOUNDARYNODES,NumberOfSides), stat = IError)
        end if
      
        K = 1
        do IElement = 1, Counters%NEl
            do ISide = 1,ELEMENTSIDES
                if (HHElemSide(IElement,ISide)) then
                    call DetermineSideNodes(ISide,LocalNodes)
                    do I=1,ELEMENTBOUNDARYNODES
                        NodeID = ElementConnectivities(LocalNodes(I), IElement) !Global name of boundary node
                        HHNodesCon(I, K) =  NodeID
                    end do
                    K = K+1
                end if
            end do
        end do
                 
        if (allocated(HHDummyLoad)) then
            deallocate(HHDummyLoad, stat = IError)
        end if
        
         if (NDIM == 3) then
            allocate(HHDummyLoad(N_BOUNDARY_NODES_HOE, NDOFL,NumberOfSides), stat = IError)
        elseif (NDIM == 2) then
            allocate(HHDummyLoad(ELEMENTBOUNDARYNODES, NDOFL,NumberOfSides), stat = IError)
        end if
        HHDummyLoad = 1.0
                 
        end subroutine GetAndUpdateHydraulicHeadLoad       
        
        logical function ConvergenceCheck() 
        !**********************************************************************
        !
        !    Function:  Check Convergence
        !
        !**********************************************************************
        
        implicit none        
        
          logical :: ConvergenceSoil, ConvergenceWater, Convergence
        
          if(NFORMULATION==1) then
              call ComputeForceError()
          else
              call ComputeForceError2LayForm()
          end if
 
          if (CalParams%ConvergenceCheck%ExternalWork==0.0) then
            CalParams%ConvergenceCheck%KineticError = CONVERGENCE_ERROR_NOT_USED
          else
            CalParams%ConvergenceCheck%KineticError =  &
              abs(CalParams%ConvergenceCheck%KineticEnergy / CalParams%ConvergenceCheck%ExternalWork)
          end if
          
          if (CalParams%ConvergenceCheck%ExternalWorkSoil==0.0) then
            CalParams%ConvergenceCheck%KineticErrorSoil = CONVERGENCE_ERROR_NOT_USED
          else
            CalParams%ConvergenceCheck%KineticErrorSoil =  &
              abs(CalParams%ConvergenceCheck%KineticEnergySoil / CalParams%ConvergenceCheck%ExternalWorkSoil)
          end if

          if (CalParams%ConvergenceCheck%ExternalWorkWater==0.0) then
            CalParams%ConvergenceCheck%KineticErrorWater = CONVERGENCE_ERROR_NOT_USED
          else
            CalParams%ConvergenceCheck%KineticErrorWater =  &
              abs(CalParams%ConvergenceCheck%KineticEnergyWater / CalParams%ConvergenceCheck%ExternalWorkWater)            
          end if
          
          if (CalParams%ApplyQuasiStatic) then
            Convergence = &
              (CalParams%ConvergenceCheck%ForceError<=CalParams%ToleratedErrorForce.or. &
               CalParams%ConvergenceCheck%ForceError==CONVERGENCE_ERROR_NOT_USED).and. &
              (CalParams%ConvergenceCheck%KineticError<=CalParams%ToleratedErrorEnergy.or. &
               CalParams%ConvergenceCheck%KineticError==CONVERGENCE_ERROR_NOT_USED)
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
              ConvergenceSoil =  &
                (CalParams%ConvergenceCheck%ForceErrorSoil<=CalParams%ToleratedErrorForce.or. &
                 CalParams%ConvergenceCheck%ForceErrorSoil==CONVERGENCE_ERROR_NOT_USED).and. &
                (CalParams%ConvergenceCheck%KineticErrorSoil<=CalParams%ToleratedErrorEnergy.or. &
                 CalParams%ConvergenceCheck%KineticErrorSoil==CONVERGENCE_ERROR_NOT_USED)
              ConvergenceWater =  &
                (CalParams%ConvergenceCheck%ForceErrorWater<=CalParams%ToleratedErrorForceWater.or. &
                 CalParams%ConvergenceCheck%ForceErrorWater==CONVERGENCE_ERROR_NOT_USED).and. &
                (CalParams%ConvergenceCheck%KineticErrorWater<=CalParams%ToleratedErrorEnergyWater.or. &
                 CalParams%ConvergenceCheck%KineticErrorWater==CONVERGENCE_ERROR_NOT_USED)
            else
              ConvergenceSoil = .true.
              ConvergenceWater = .true.
            end if
            ConvergenceCheck = ((Convergence.and.ConvergenceSoil.and.ConvergenceWater).or. &
                                 (CalParams%TimeStep>CalParams%MaxTimeSteps))
          else 
            ! Use the real time instead of number of increments
            ConvergenceCheck = (CalParams%TotalRealTime>=CalParams%TotalTime)
          end if
                    
        end function ConvergenceCheck
        
        logical function DivergenceCheck() 
        !**********************************************************************
        !
        !    Function:  Check Divergence
        !
        !    Divergence in case of negative energy dissipation
        !    Divergence tolerance is -ToleratedDivergence (in J)
        !
        !**********************************************************************
        
        implicit none        
        
          logical :: DivergenceSoil, DivergenceWater, Divergence
 
          CalParams%ConvergenceCheck%EnergyDissipation = CalParams%ConvergenceCheck%ExternalWork -  &
                                                         CalParams%ConvergenceCheck%InternalWork - &
                                                         CalParams%ConvergenceCheck%KineticEnergy + &
                                                         CalParams%ConvergenceCheck%KineticEnergy0
          
          Divergence = (CalParams%ConvergenceCheck%EnergyDissipation<=(-1 * CalParams%ToleratedDivergence))
          
          if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
            CalParams%ConvergenceCheck%EnergyDissipationSoil = CalParams%ConvergenceCheck%ExternalWorkSoil -  &
                                                               CalParams%ConvergenceCheck%InternalWorkSoil - &
                                                               CalParams%ConvergenceCheck%KineticEnergySoil + &
                                                               CalParams%ConvergenceCheck%KineticEnergySoil0
            
            CalParams%ConvergenceCheck%EnergyDissipationWater = CalParams%ConvergenceCheck%ExternalWorkWater -  &
                                                                CalParams%ConvergenceCheck%InternalWorkWater - &
                                                                CalParams%ConvergenceCheck%KineticEnergyWater + &
                                                                CalParams%ConvergenceCheck%KineticEnergyWater0
            
            DivergenceSoil = (CalParams%ConvergenceCheck%EnergyDissipationSoil<=(-1 * CalParams%ToleratedDivergence))
            
            DivergenceWater = (CalParams%ConvergenceCheck%EnergyDissipationWater<=(-1 * CalParams%ToleratedDivergence))
          else
            DivergenceSoil = .false.
            DivergenceWater = .false.
          end if
          
          DivergenceCheck = (Divergence.or.DivergenceSoil.or.DivergenceWater)
          
        end function DivergenceCheck
                
         subroutine DynUpdateParticleWeights( )
        !**********************************************************************
        !
        !    Function:  Update the weights of the material points.
        ! 
        !
        !**********************************************************************
 
        implicit none
       
          ! Local variables
          integer(INTEGER_TYPE) :: I, IAEl, IEl, Int, IntGlo, MaterialID, J
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
          real(REAL_TYPE) :: VolumetricStrain, ConstDensitySolid, ConstDensityLiquid
          real(REAL_TYPE) :: LiquidPressure, Density, ConcRatioSolid
          real(REAL_TYPE) :: RatioDensity
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ElementNodalDensity, ElementNodalConcRatio
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          logical :: PressureGreaterThreshold, PressureLowerThreshold, FullyFilled, FreeSurface 
          logical :: UpdDensityLiquid
        
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! 1-layer formulation or 2-layer form with 1 Phase
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if((NFORMULATION==1).or. &
               ((NFORMULATION==2).and.(CalParams%NumberOfPhases==1))) then
                                
          do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)   
            
            do Int = 1, NPartEle(IEl) ! Loop over all material points of the element
                
              ! Determine global ID of Material Point 
              IntGlo = GetParticleIndex(Int, IEl)
              
              ! Determine volumetric strain of Material Point 
              VolumetricStrain = 0.0
              do I = 1, NVECTOR ! for 2D and 3D all three components have to be taken into account as eps_tt in axisymmetric \=0
                 VolumetricStrain = VolumetricStrain + GetEpsStepI(Particles(IntGlo),I)
              end do
          
             if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then
                 
                 ! Consider stress state 
                 LiquidPressure = (SigmaEffArray(IntGlo,1) + SigmaEffArray(IntGlo,2) + SigmaEffArray(IntGlo,3))/3.0 ! this is valid for 2D and 3D
     
                 MaterialID = MaterialIDArray(IntGlo)
     
                 ! Determine interpolated density
                 call CalculateRatioDensity(IEl,IntGlo,RatioDensity)
      
                 PressureGreaterThreshold = (Particles(IntGlo)%WaterPressure>=CalParams%LiquidPressureCavitationThreshold)
                 PressureLowerThreshold = (Particles(IntGlo)%WaterPressure<CalParams%LiquidPressureCavitationThreshold)
                 FreeSurface = (Particles(IntGlo)%LiquidFreeSurface==1.0)
                 FullyFilled = (RatioDensity>1.0)
               

                 UpdDensityLiquid = (((.not.FreeSurface).and.PressureLowerThreshold).or. &
                    ((PressureGreaterThreshold.or.FreeSurface).and.FullyFilled))
                 
                 if (.not.UpdDensityLiquid) then
                    Particles(IntGlo)%Density = MatParams(MaterialID)%FluidThresholdDensity/1000         
                 else 
                    Particles(IntGlo)%LiquidFreeSurface = 0.0 ! Particle doesn't belong to the free surface
                    Particles(IntGlo)%Density = Particles(IntGlo)%Density /  (VolumetricStrain + 1.0)
                 end if
                 
                 Particles(IntGlo)%IntegrationWeight = MassArray(IntGlo) / Particles(IntGlo)%Density

             else
                 ! for SOLID or MIXTURE Material Point               
                 Particles(IntGlo)%IntegrationWeight = Particles(IntGlo)%IntegrationWeight * (1.0 + VolumetricStrain)
             end if
             end do !loop over material points
            end do ! loop over elements
            
          end if ! 1-layer formulation or 2-layer form with 1 Phase
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! 2-layer formulation with 2 Phases
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if((NFORMULATION==2).and.(CalParams%NumberOfPhases==2)) then
              
          !--- Calculate  density  ----
              ConstDensityLiquid = 0.0
              ConstDensitySolid = 0.0
                
              do J = 1, Counters%NLayers ! loop over all material points in element
                if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-solid'.or.MatParams(j)%MaterialPhases=='1-phase-solid') then
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)    
                end if
              end do 
          
          do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)   
            
            do Int = 1, NPartEle(IEl) ! Loop over all material points of the element
                
              ! Determine global ID of Material Point 
              IntGlo = GetParticleIndex(Int, IEl)
              
              ! Determine volumetric strain of Material Point 
              VolumetricStrain = 0.0
              do I = 1, NVECTOR ! for 2D and 3D all three components have to be taken into account as eps_tt in axisymmetric \=0
                 VolumetricStrain = VolumetricStrain + GetEpsStepI(Particles(IntGlo),I)
              end do
          
             if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then
     
              PressureGreaterThreshold = (Particles(IntGlo)%WaterPressure>=CalParams%LiquidPressureCavitationThreshold)
              PressureLowerThreshold = (Particles(IntGlo)%WaterPressure<CalParams%LiquidPressureCavitationThreshold)
              FreeSurface = (Particles(IntGlo)%LiquidFreeSurface==1.0)
              if(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0) then
               FullyFilled = (Particles(IntGlo)%FillingRatio>0.98)
              else
               FullyFilled = (Particles(IntGlo)%FillingRatio>1.0)
              end if

              UpdDensityLiquid = (((.not.FreeSurface).and.PressureLowerThreshold).or. &
              ((PressureGreaterThreshold.or.FreeSurface).and.FullyFilled))
                 
                 if (.not.UpdDensityLiquid) then
                  if(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0) then
                    ! Compute density from nodal density field
                    NodeIDs = ElementConnectivities(1:ELEMENTNODES, IEl)    
                    ElementNodalConcRatio = TwoLayerData%Nodes(NodeIDs)%DensitySolidS / ConstDensitySolid
                 
                    ConcRatioSolid = 0.0
                    ParticleShape = ShapeValuesArray(IntGlo,:)
                    Particles(IntGlo)%Density = 0.0
                    do I = 1, ELEMENTNODES
                      ConcRatioSolid = ConcRatioSolid + ParticleShape(I) * ElementNodalConcRatio(I) 
                    end do
                  else
                    ConcRatioSolid = 0.0
                  end if
                    Particles(IntGlo)%Density = Particles(IntGlo)%ConstDensity * (1.0 - ConcRatioSolid)        
                 else ! For output purposes, possible future use, update density
                    Particles(IntGlo)%LiquidFreeSurface = 0.0d0 ! Particle doesn't belong to the free surface
                    Particles(IntGlo)%Density = Particles(IntGlo)%Density / (VolumetricStrain + 1.0)
                 end if
                 
                 Particles(IntGlo)%IntegrationWeight = MassArray(IntGlo) / Particles(IntGlo)%Density
                 Particles(IntGlo)%EffConcentrationRatioLiquid = Particles(IntGlo)%Density / Particles(IntGlo)%ConstDensity
                 
              end if ! LIQUID Material Point
                 
             if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeSolid) then
                 
                 if(Particles(IntGlo)%PhaseStatus==PhaseStatusLIQUID) then 
                    ! Compute density from nodal density field
                    NodeIDs = ElementConnectivities(1:ELEMENTNODES, IEl)
                    ElementNodalDensity = TwoLayerData%Nodes(NodeIDs)%DensitySolidS
                                        
                    Density = 0.0
                    ParticleShape = ShapeValuesArray(IntGlo,:)
                    Particles(IntGlo)%Density = 0.0
                    do I = 1, ELEMENTNODES
                      Density = Density + ParticleShape(I) * ElementNodalDensity(I) 
                    end do
                    Particles(IntGlo)%Density = Density
                    
                    ! Calculate new effective porosity
                    Particles(IntGlo)%EffConcentrationRatioSolid = Particles(IntGlo)%Density / Particles(IntGlo)%ConstDensity
                    Particles(IntGlo)%EffPorosity = 1.0 - Particles(IntGlo)%EffConcentrationRatioSolid
                    
                    if(Particles(IntGlo)%EffPorosity<CalParams%LimitPorosity) then
                      Particles(IntGlo)%EffPorosity = CalParams%LimitPorosity        
                      Particles(IntGlo)%EffConcentrationRatioSolid = 1.0 - Particles(IntGlo)%EffPorosity
                      Particles(IntGlo)%Density = Particles(IntGlo)%EffConcentrationRatioSolid * Particles(IntGlo)%ConstDensity
                    end if
                    
                 else  ! Not Fluidized soil
                    Particles(IntGlo)%Density = Particles(IntGlo)%Density / (VolumetricStrain + 1.0)
                 end if
                 
                 Particles(IntGlo)%IntegrationWeight = MassArray(IntGlo) / Particles(IntGlo)%Density
                 Particles(IntGlo)%EffConcentrationRatioSolid = Particles(IntGlo)%Density / Particles(IntGlo)%ConstDensity
                 Particles(IntGlo)%EffPorosity = 1.0 - Particles(IntGlo)%EffConcentrationRatioSolid
                 
              end if   ! SOLID Material Point              
                 
             end do !loop over material points
            end do ! loop over elements
            
          end if !2-layer formulation with 2 Phases
         
        end subroutine DynUpdateParticleWeights
        

        
        subroutine DynUpdateParticlePorosity( )
        !**********************************************************************
        !
        !    Function:  Updates particle porosity
        !
        !    - The "mass balance of the solid" has been considered. 
        !    - The variation of the solid density is neglected.          
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 
        implicit none
       
          ! Local variables
          integer(INTEGER_TYPE) :: I, IParticle, ParticleIndex
          real(REAL_TYPE) :: DEpsVol

          do IParticle = 1, Counters%NParticles  ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            DEpsVol = 0.0
             do I = 1, NVECTOR
              DEpsVol = DEpsVol + GetEpsStepI(Particles(ParticleIndex),I)
             end do

            Particles(ParticleIndex)%Porosity = Particles(ParticleIndex)%Porosity + (1.0d0 - Particles(ParticleIndex)%Porosity) * DEpsVol
            
            if (Particles(ParticleIndex)%Porosity<0.0) then
                Particles(ParticleIndex)%Porosity = 0.0
            end if
            
          end do ! loop over particles
         
        end subroutine DynUpdateParticlePorosity

        
        subroutine DynUpdateParticleDegreeOfSaturation( )
        !**********************************************************************
        !
        !    Function:  Updates Degree of Saturation (Sr) of the particle 
        !               The expression of the Degree of Saturation is the RETENTION CURVE
        !
        !**********************************************************************
        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: ParticleIndex, IParticle, ISet
        real(REAL_TYPE) :: Sr
        
          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            ISet = MaterialIDArray(ParticleIndex)
        
              if (MatParams(ISet)%RetentionCurve== SWRC_VANGENUCHTEN) then
                  call DynUpdateParticleDegreeOfSaturationVanGenuchten(ParticleIndex,ISet,Sr)
              else if (MatParams(ISet)%RetentionCurve==SWRC_LINEAR) then
                  call DynUpdateParticleDegreeOfSaturationLinear(ParticleIndex,ISet,Sr)
              end if
          
          end do
        
        end subroutine DynUpdateParticleDegreeOfSaturation
        
        
        subroutine DynUpdateParticleHydraulicConductivity( )
         !**********************************************************************
        !
        !    Function:  Updates hydraulic conductivity (K) of the particle 
        !    as a function of Degree of Saturation (Sr)
        !               
        !
        !**********************************************************************
        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: ParticleIndex, IParticle, ISet
        
          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            ISet = MaterialIDArray(ParticleIndex)
            
            if (MatParams(ISet)%HydraulicConductivityCurve==HCC_CONSTANT) then
            ! do nothing, Kw doesn't vary
            else if (MatParams(ISet)%HydraulicConductivityCurve==HCC_HILLEL) then
                  call DynUpdateParticleHydraulicConductivityHillel(ParticleIndex,ISet)
            else if (MatParams(ISet)%HydraulicConductivityCurve==HCC_MUALEM) then
                  call DynUpdateParticleHydraulicConductivityMualem(ParticleIndex,ISet)  
            end if
          
          end do
        
        end subroutine DynUpdateParticleHydraulicConductivity
        
        subroutine DynUpdateParticleHydraulicConductivityHillel(ParticleIndex,ISet)
        !**********************************************************************
        !
        !    Function:  Updates hydraulic conductivity (K) of the particle 
        !    as a function of Degree of Saturation (Sr)
        !   The expression of the Degree of Saturation is the RETENTION CURVE, 
        !   while the expression of K is based on   Hillel (1971)
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 

        implicit none
        ! Local variables
        integer(INTEGER_TYPE), intent (in):: ParticleIndex, ISet
        real(REAL_TYPE) :: Sr, r, ksat
        
        Sr = Particles(ParticleIndex)%DegreeSaturation
        if (Sr >=1.0)  RETURN

        r = MatParams(ISet)%rexponentHillel_HCC
        ksat = MatParams(ISet)%HydraulicConductivityLiquid !% this is the saturated permeability
        
        Particles(ParticleIndex)%Conductivity = ksat * (Sr**(r))
        
        end subroutine DynUpdateParticleHydraulicConductivityHillel
        
        subroutine DynUpdateParticleHydraulicConductivityMualem(ParticleIndex,ISet)
        !**********************************************************************
        !
        !    Function:  Updates hydraulic conductivity (K) of the particle 
        !    as a function of Degree of Saturation (Sr)
        !   The expression of the Degree of Saturation is the RETENTION CURVE, 
        !   while the expression of K is based on   Mualem (1976)
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 
        implicit none
        ! Local variables
        integer(INTEGER_TYPE), intent (in):: ParticleIndex, ISet
        real(REAL_TYPE) :: Sr, ksat, krel
        real(REAL_TYPE) :: Smin, Smax, L
        real(REAL_TYPE) :: n00, w1, w2, w3, w4, w5, w6
        
        Sr = Particles(ParticleIndex)%DegreeSaturation
         if (Sr >=1.0)  return
        !
        ksat = MatParams(ISet)%HydraulicConductivityLiquid
        Smin = MatParams(ISet)%Smin_SWRC
        Smax = MatParams(ISet)%Smax_SWRC
        L = MatParams(ISet)%Lambda_SWRC
        
        n00 = 1.0d0
        w1 = sqrt(Sr)
        w3 = n00/L
        w2 = Sr**(w3)
        w4 = (1-w2)**(L)
        w5 = 1-w4
        w6 = w5*w5
        
        krel = w1 * w6
        Particles(ParticleIndex)%Conductivity = ksat * krel
        
        end subroutine DynUpdateParticleHydraulicConductivityMualem
       
        
        subroutine DynUpdateParticleLiquidWeight()
        !**********************************************************************
        !
        !    Function: updates the liquid density (WD)
        !               
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: IParticle, ParticleIndex
        real(REAL_TYPE) ::  Pw, T, Kw, g
        real(REAL_TYPE) ::  Alpha, Beta, Pw0, WD0

          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)

            Pw = -Particles(ParticleIndex)%WaterPressure
            T = Particles(ParticleIndex)%Temperature   !Temperature (ºC)
            Kw = Particles(ParticleIndex)%BulkWater
            g = CalParams%GravityData%GAccel  !Gravity (m/s2)
            
            Alpha = -0.00034d0  !Default=-0.00034, volumetric thermal expansion coefficient for water (1/ºC)
            Beta = 1/Kw         !Default=0.00045, water compressibility (1/MPa) ---> Inverse of the bulk modulus
            Pw0 =  0.0d0        !Default=0.1, Reference water pressure (MPa) 
            WD0 = 1.0026d0      !Default=1002.6, Reference water density (kg/m3)
        
            Particles(ParticleIndex)%WaterWeight = WD0*exp(Beta*(Pw-Pw0) + Alpha*T)*g
     
          end do ! loop over particles

        end subroutine DynUpdateParticleLiquidWeight


        subroutine DynUpdateParticleGasWeight()
        !**********************************************************************
        !
        !    Function: Updates the gas density (GD)
        !               with respect de Pw, Pg, and T
        !               The expression of the Gas Density is based
        !               on the IDEAL GASES LAW (PV=nRT)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: IParticle, ParticleIndex
        real(REAL_TYPE) ::  Pw, Pg, T, WD
        real(REAL_TYPE) ::  g, Suc, Ma, Mw, R, Pv0, F, Pv
        real(REAL_TYPE) ::  n1, n2

          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)

            g = CalParams%GravityData%GAccel  !Gravity (m/s2)
            
            Pw = Particles(ParticleIndex)%WaterPressure
            Pg = Particles(ParticleIndex)%GasPressure
            T = Particles(ParticleIndex)%Temperature           !Temperature (ºC)
            WD = Particles(ParticleIndex)%WaterWeight/g        !Water density of the particle
      
            Suc = Pw-Pg     !Suction

            if (Suc<=0.0) RETURN

            T = T+273.15d0              !Temperature (ºK)
            Ma = 0.028d0                !Default=0.02895, Molecular gas of dry air (kg/mol)
            Mw = 0.018d0                !Default=0.01801528, Molecular mass of water (=vapour) (kg/mol)
            R = 0.008314d0              !Default=8.3144621, Gas constant value (J/mol*K)

            !Calculation of the vapour Pressure (Pg = Pv + Pa) using the psychrometric law 
            n1 = 136075.0d0
            n2 = -5239.7d0
            Pv0 = n1*exp(n2/T)
            F = exp(-(Mw*Suc)/(R*T*WD))     !Psychrometric Law
            Pv = Pv0*F                      !Vapour Pressure

            Particles(ParticleIndex)%GasWeight = (Pv*(Mw-Ma) + Pg*Ma)/(R*T)*g
        
          end do ! loop over particles

        end subroutine DynUpdateParticleGasWeight


        subroutine DynUpdateParticleMixedWeight
        !**********************************************************************
        !
        !    Function: Updates the density of the mixture
        !               with respect de new densities, porosity and degree of saturation
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: IParticle, ParticleIndex
        real(REAL_TYPE) ::  MW, WW, GW, Sr, Sg, N, MixedWeight
        real(REAL_TYPE), dimension(NVECTOR) :: GravityVector
        
        GravityVector = CalParams%GravityData%GravityVector

          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
 
          MW = Particles(ParticleIndex)%MaterialWeight      ! Material weight (dry) of the particle       
          WW = Particles(ParticleIndex)%WaterWeight         ! Water weight of the particle
          GW = Particles(ParticleIndex)%GasWeight           ! Gas weight of the particle
          Sr = Particles(ParticleIndex)%DegreeSaturation    ! Degree of Saturation (liquid) of the particle
          Sg = 1.0d0 - Sr                                   ! Degree of Satutation (gas) of the particle
          N = Particles(ParticleIndex)%Porosity             ! Porosity of the particle
          MixedWeight = MW + N*Sr*WW + N*Sg*GW
          

          Particles(ParticleIndex)%MixedWeight = MixedWeight
          Particles(ParticleIndex)%FbodyMixed(:) = MixedWeight*GravityVector(:)*Particles(ParticleIndex)%IntegrationWeight
    
          end do ! loop over particles
        
        end subroutine DynUpdateParticleMixedWeight
        
        

        subroutine DynUpdateParticleAirInWaterMassFraction( )
        !**********************************************************************
        !
        !    Function:  Update Mass Fraction (Henry's Law)
        !
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 
        implicit none
       
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex
          real(REAL_TYPE) :: H, Mw, Ma, HM     ! Henry's constant, molecular water mass, molecular air mass   
     
          H = 100.0              ! Taking Henry's coefficient as a constant. Default=10000000 (KPa)
          Ma = 0.028d0                ! Default=0.02895, Molecular gas of dry air (kg/mol)
          Mw = 0.018d0                ! Default=0.01801528, Molecular mass of water (kg/mol)
          HM = Ma/(H*Mw)
            
          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            
            Particles(ParticleIndex)%AirInWaterMassFraction = Particles(ParticleIndex)%GasPressure*HM
     
           end do ! loop over particles
         
        end subroutine DynUpdateParticleAirInWaterMassFraction
        
        
        
        subroutine DynUpdateParticleVapourInGasMassFraction( )
        !**********************************************************************
        !
        !    Function:  Update Mass Fraction (Psychrometric's Law)
        !
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 
        implicit none
       
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex
          real(REAL_TYPE) ::  Pw, Pg, T, wVG, WD
          real(REAL_TYPE) ::  g, Mw, R, Suc

          g = CalParams%GravityData%GAccel  !Gravity (m/s2)
          Mw = 0.018d0                !Default=0.01801528, Molecular mass of water (kg/mol)
          R = 0.008314d0              !Default= 8.3144621, Gas constant value (J/mol*K)

          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            
            Pw = Particles(ParticleIndex)%WaterPressure
            Pg = Particles(ParticleIndex)%GasPressure
            T = Particles(ParticleIndex)%Temperature           !Temperature (ºC)
            wVG = Particles(ParticleIndex)%VapourInGasMassFraction
            WD = Particles(ParticleIndex)%WaterWeight/g        !Water density of the particle
      
            T = T+273.15d0  !Temperature (ºK)
            Suc = Pw-Pg     !Suction

           Particles(ParticleIndex)%VapourInGasMassFraction = exp(-Mw*Suc/(R*T*WD))
     
           end do ! loop over particles
         
        end subroutine DynUpdateParticleVapourInGasMassFraction



        subroutine UpdateNodes(EntityUsed)
        !**********************************************************************
        !
        !    Function:  Updating of nodal coordinates from displacements for
        !               updated mesh analysis.
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE) :: EntityUsed
          
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, J0
          
          do I = 1, Counters%NodTot
            J0 = ReducedDof(I)
            do J = 1, NVECTOR
              
              NodalCoordinatesUpd(I, J) = NodalCoordinatesUpd(I, J) +  IncrementalDisplacementSoil(J0 + J,EntityUsed)
                
            end do
          end do
          
        end subroutine UpdateNodes

        
        subroutine UpdateNodesMeshAdjust(IncUTot)
        !**********************************************************************
        !
        !    Function:  Updating of nodal coordinates from displacements for
        !               updated mesh analysis.
        ! 
        !     IncUTot : Incremental nodal displacements of load step
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N),intent(in) :: IncUTot  
     
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, J0
          
          do I = 1, Counters%NodTot
            J0 = ReducedDof(I)
            do J = 1, NVECTOR
              
              NodalCoordinatesUpd(I, J) = NodalCoordinates(I, J) + IncUTot(J0 + J)  
             end do
          end do
                   
        end subroutine UpdateNodesMeshAdjust
    

        subroutine StressAndPorePressureSmoothening()
        !**********************************************************************
        !
        !    Function: Smoothing of stress and pore pressure in case of GP integration.
        !
        ! Implemented in the frame of the MPM project.
        !
        ! Note : This subroutine works only if NumberOfLayers = 1 !!! Temprary Solution
        !
        !**********************************************************************
        
        implicit none
      
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, IEl, NElemPart, IPart, ParticleIndex
          integer(INTEGER_TYPE) :: CheckMaterialID, IAEl
          real(REAL_TYPE) :: ParticleStress
          real(REAL_TYPE) :: WPAverage, DPBVDAverage
          real(REAL_TYPE) :: GPAverage
          real(REAL_TYPE), dimension(NTENSOR) :: StressAverage
          logical :: IsMixedElement
          real(REAL_TYPE) :: ParticleVol, ParticlesVols, AverageWeight
          real(REAL_TYPE), dimension(NTENSOR) :: SigmaEffStep, SigmaEff

          if(.not.(NFORMULATION==1)) RETURN

          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
              NElemPart = NPartEle(IEl) 
              StressAverage = 0.0 ! Initialise average stress of element IEl
              WPAverage = 0.0 ! Initialise average water pressure of element IEl
              GPAverage = 0.0 ! Initialise average gas pressure of element IEl
              DPBVDAverage = 0.0
              
              ParticleIndex = GetParticleIndex(1, IEl)
              CheckMaterialID = MaterialIDArray(ParticleIndex)
              IsMixedElement = .false.
              
              ParticlesVols = 0.0
              ParticleVol = 0.0

              if (CalParams%ApplyEmptyElements) then
                AverageWeight = 0.0
                if (NMaterialParticlesEle(IEl)>0) then
                  do IPart = 1, NElemPart
                    ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                    if (IsMaterialParticle(ParticleIndex)) then
                      AverageWeight = AverageWeight + Particles(ParticleIndex)%IntegrationWeight
                    end if
                  end do
                  AverageWeight = AverageWeight / NMaterialParticlesEle(IEl)
                else ! Only virtual particle(s)
                  AverageWeight = 1.0
                end if
              end if

              do IPart = 1, NElemPart ! Loop over particles
                ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                
                if (MaterialIDArray(ParticleIndex)/= CheckMaterialID) then
                  IsMixedElement = .true.
                end if 
                
                if (CalParams%ApplyEmptyElements .and. (.not.IsMaterialParticle(ParticleIndex))) then
                  ParticleVol = AverageWeight ! Use average weight
                else  
                  ParticleVol = Particles(ParticleIndex)%IntegrationWeight
                end if
                ! Sum up stress components, water pressure and gas pressure
                do I = 1, NTENSOR
                  ParticleStress = SigmaEffArray(ParticleIndex,I) 
                  StressAverage(I) = StressAverage(I) + ParticleStress * ParticleVol
                end do
                
                WPAverage = WPAverage + Particles(ParticleIndex)%WaterPressure * ParticleVol
                DPBVDAverage = DPBVDAverage + Particles(ParticleIndex)%DBulkViscousPressure * ParticleVol
                GPAverage = GPAverage + Particles(ParticleIndex)%GasPressure * ParticleVol
                ParticlesVols = ParticlesVols + ParticleVol
              end do

              ! Take the average of each stress component, the water pressure and the gas pressure
              do J = 1, NTENSOR
                StressAverage(J) = StressAverage(J) / ParticlesVols
              end do
              WPAverage = WPAverage / ParticlesVols
              GPAverage = GPAverage / ParticlesVols
              DPBVDAverage = DPBVDAverage / ParticlesVols

              ! Assign the average stress, water pressure and gas pressure to particles
              do IPart = 1, NElemPart
                ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                if (.not.IsMixedElement) then
                    SigmaEffArray(ParticleIndex, :) = StressAverage(:)
                    Particles(ParticleIndex)%WaterPressure = WPAverage
                    Particles(ParticleIndex)%DBulkViscousPressure = DPBVDAverage
                    Particles(ParticleIndex)%GasPressure = GPAverage
                end if
              end do

              do IPart = 1, NPartEle(IEl) ! Loop over particles in IElement
                ParticleIndex = GetParticleIndex(IPart, IEl)
 
                SigmaEffStep = SigmaEffArray(ParticleIndex,:) - SigmaEff0Array(ParticleIndex,:)
                call SetSigmaEffStep(Particles(ParticleIndex), SigmaEffStep)

                SigmaEff =SigmaEffArray(ParticleIndex,:)
                SigmaEff0Array(ParticleIndex,:)=SigmaEff

                Particles(ParticleIndex)%WaterPressure0 = Particles(ParticleIndex)%WaterPressure

                Particles(ParticleIndex)%GasPressure0 = Particles(ParticleIndex)%GasPressure

              end do
          end do

        end subroutine StressAndPorePressureSmoothening
        
        subroutine StateParametersSmoothening()
        !**********************************************************************
        !
        !    Function: Smoothing of state parameters in case of GP integration.
        !
        ! Implemented in the frame of the MPM project.
        !
        ! Note : This subroutine works only if NumberOfLayers = 1 !! Temprary Solution
        !
        !**********************************************************************
        implicit none
      
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, NElemPart, IPart, ParticleIndex
          integer(INTEGER_TYPE) :: CheckMaterialID, IAEl
          real(REAL_TYPE) :: ParticleStateVariable
          real(REAL_TYPE) :: AverageHPStateVariables (2), &
                              AverageModifiedHPStateVariables(2), &
                              AverageHPIGStateVariables (7), &
                              AverageEpsP(NTENSOR), &
                              AverageSigmaPrin(NTENSOR)
          real(REAL_TYPE) :: AverageCohesionStSoft, &
                              AveragePhiStSoft, &
                              AveragePsiStSoft, &
                              AveragePP
          real(REAL_TYPE) :: AverageESMstatev(NSTATEVAR) !user defined model
          logical :: IsMixedElement
          real(REAL_TYPE) :: ParticleVol, ParticlesVols, AverageWeight
     
         if(.not.(NFORMULATION==1)) RETURN
          
         if (.not.IsMPMComputation()) RETURN ! if FEM no need for mapping...already one particle
         
          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
              NElemPart = NPartEle(IEl)
              AverageHPStateVariables         = 0.0
              AverageModifiedHPStateVariables = 0.0
              AverageHPIGStateVariables       = 0.0
              AverageEpsP = 0.0
              AverageSigmaPrin = 0.0
              AverageCohesionStSoft = 0.0
              AveragePhiStSoft = 0.0
              AveragePsiStSoft = 0.0
              AveragePP = 0.0
              
              AverageESMstatev = 0.0
              
              ParticleIndex = GetParticleIndex(1, IEl)
              CheckMaterialID = MaterialIDArray(ParticleIndex) ! the same as material ID
              IsMixedElement = .false.
              
              ParticlesVols = 0.0
              ParticleVol = 0.0

              if (CalParams%ApplyEmptyElements) then
                AverageWeight = 0.0
                if (NMaterialParticlesEle(IEl)>0) then
                  do IPart = 1, NElemPart
                    ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                    if (IsMaterialParticle(ParticleIndex)) then
                      AverageWeight = AverageWeight + Particles(ParticleIndex)%IntegrationWeight
                    end if
                  end do
                  AverageWeight = AverageWeight / NMaterialParticlesEle(IEl)
                else ! Only virtual particle(s)
                  AverageWeight = 1.0
                end if
              end if

              do IPart = 1, NElemPart ! Loop over particles
                ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                
                if (MaterialIDArray(ParticleIndex) /= CheckMaterialID) then
                  IsMixedElement = .true.
                  goto 100 ! skip this element 
                end if 
              
                if (CalParams%ApplyEmptyElements .and. (.not.IsMaterialParticle(ParticleIndex))) then
                  ParticleVol = AverageWeight ! Use average weight
                else  
                  ParticleVol = Particles(ParticleIndex)%IntegrationWeight
                end if

               ! Sum up state variables for MCStrainSoftening model
                do I = 1, NTENSOR
                  ParticleStateVariable = GetEpsPI(Particles(ParticleIndex), I)  
                  AverageEpsP(I) = AverageEpsP(I) + ParticleStateVariable * ParticleVol
                  ParticleStateVariable = GetSigmaPrinI(Particles(ParticleIndex), I)  
                  AverageSigmaPrin(I) = AverageSigmaPrin(I) + ParticleStateVariable * ParticleVol
                end do
                AverageCohesionStSoft = AverageCohesionStSoft + Particles(ParticleIndex)%CohesionStSoft * ParticleVol
                AveragePhiStSoft = AveragePhiStSoft + Particles(ParticleIndex)%PhiStSoft * ParticleVol
                AveragePsiStSoft = AveragePsiStSoft + Particles(ParticleIndex)%PsiStSoft * ParticleVol
                
              ! Sum up state variables for MCC model
                AveragePP = AveragePP + Particles(ParticleIndex)%PP * ParticleVol
                

                do I = 1, 2  ! HP model 
                  ParticleStateVariable = GetHPStateVariablesI(Particles(ParticleIndex), I)
                  AverageHPStateVariables(I) = AverageHPStateVariables(I) + ParticleStateVariable * ParticleVol
                end do

                do I = 1, 7  ! HPIG model 
                  ParticleStateVariable = GetHPIGStateVariablesI(Particles(ParticleIndex), I)
                  AverageHPIGStateVariables(I) = AverageHPIGStateVariables(I) + ParticleStateVariable * ParticleVol
                end do

                ! Modified HP model
                do I = 1, 2
                  ParticleStateVariable = GetModifiedHPStateVariablesI(Particles(ParticleIndex), I)
                  AverageModifiedHPStateVariables(I) = AverageModifiedHPStateVariables(I) + ParticleStateVariable * ParticleVol
                end do

                AverageESMstatev=  AverageESMstatev +  ESMstatevArray(ParticleIndex,:) * ParticleVol
   
                ParticlesVols = ParticlesVols + ParticleVol
              end do ! loop over particles


              ! Take the average of each state parameter

              ! HP model
              AverageHPStateVariables = AverageHPStateVariables / ParticlesVols

              ! HPIG model
              AverageHPIGStateVariables = AverageHPIGStateVariables / ParticlesVols
          

              ! modified HP model
              AverageModifiedHPStateVariables = AverageModifiedHPStateVariables / ParticlesVols

              AverageESMstatev =  AverageESMstatev / ParticlesVols
 
              ! Strain softening model
              AverageEpsP = AverageEpsP / ParticlesVols
              AverageSigmaPrin = AverageSigmaPrin / ParticlesVols

              AverageCohesionStSoft = AverageCohesionStSoft / ParticlesVols
              AveragePhiStSoft = AveragePhiStSoft / ParticlesVols
              AveragePsiStSoft = AveragePsiStSoft / ParticlesVols
              
              ! MCC model
              AveragePP = AveragePP / ParticlesVols
              
              ! Assign the average state parameters to particles
              if (.not.IsMixedElement) then ! double check for mixed element
                do IPart = 1, NElemPart
                  ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
            
                  call SetHPStateVariables(Particles(ParticleIndex), AverageHPStateVariables)
                  call SetHPIGStateVariables(Particles(ParticleIndex), AverageHPIGStateVariables)
                  call SetModifiedHPStateVariables(Particles(ParticleIndex), AverageModifiedHPStateVariables)
     
                   ! Strain softening model
                   call SetEpsP(Particles(ParticleIndex), AverageEpsP)
                   call SetSigmaPrin(Particles(ParticleIndex), AverageSigmaPrin)
                   Particles(ParticleIndex)%CohesionStSoft = AverageCohesionStSoft
                   Particles(ParticleIndex)%PhiStSoft = AveragePhiStSoft
                   Particles(ParticleIndex)%PsiStSoft = AveragePsiStSoft
                   
                   ! MCC model
                   Particles(ParticleIndex)%PP = AveragePP
                   !user-defined model
                   ESMstatevArray(ParticleIndex,:) = AverageESMstatev
                end do
              end if
                         
100        end do ! elements

          end subroutine StateParametersSmoothening
                
        subroutine MapVeloToParticle(NodalVelocities)
                              
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        !    NodalVelocities : Nodal accelerations vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: NodalVelocities
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IPart, INode, IDof, ParticleIndex, NodeID
          real(REAL_TYPE), dimension (NVECTOR) :: ParticleIncrementalVelocity
          real(REAL_TYPE), dimension (NVECTOR) :: ParIncVel

          do IEl = 1, Counters%NEl   ! Loop over all elements
              do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
                
                ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                ParticleIncrementalVelocity = 0.0
                ParIncVel = 0.0
                
                  do INode = 1, ELEMENTNODES

                    NodeID = ElementConnectivities(INode, IEl)    ! Global node ID
                    IDof = ReducedDof(NodeID)

                    do I = 1, NVECTOR
                      ParIncVel(I) = ShapeValuesArray(ParticleIndex,INode) * NodalVelocities(IDof+1, 1)
                      ParticleIncrementalVelocity(I) = ParticleIncrementalVelocity(I) + ParIncVel(I)
                    end do  
               
               end do !Loop over nodes
                   
                 VelocityArray(ParticleIndex,:) = ParticleIncrementalVelocity
     
            end do !Loop over particles
          end do !elements
          
         end subroutine MapVeloToParticle

        subroutine GetNodalVelocityFromNodalMomentumConv(Momentum)
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
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Momentum
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J

         do IDOF = 1, Counters%N                                       !Loop over all degrees of freedom
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
                if(LumpedMassDry(IDOF,J)/=0) then
                    TotalVelocitySoil(IDOF,J) = ( Momentum(IDOF,J) /  &
                                        LumpedMassDry(IDOF,J) )*  &
                                        PBoundary(IDOF)
                else
                    TotalVelocitySoil(IDOF,J) = 0.0
                end if
            end do
          end do
          
         end subroutine GetNodalVelocityFromNodalMomentumConv

         subroutine UpdateParticleDisplacementsWater(DUTot)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements (water)
        !
        !     DUTot : Obtained incremental displacements (water)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTot
          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID, IParticle
          integer(INTEGER_TYPE) :: ElementID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement
          
          do IParticle = 1, Counters%NParticles ! Loop over particles
            
            ElementID = ElementIDArray(IParticle)
          
            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle) 
            else
              iEntity = 1
            end if 
              
            ParticleDisplacement = 0.0
                    
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                NodeID = iabs(ElementConnectivities(INode, ElementID)) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * DUTot(DofID,iEntity)                             
              end do                                             
            end do 
                          
            ! Total particle displacement
            Particles(IParticle)%UW = Particles(IParticle)%UW + ParticleDisplacement
        
         end do ! particles
          
         end subroutine UpdateParticleDisplacementsWater
         
         subroutine UpdateParticleDisplacementsGas(DUTot)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements (gas)
        !
        !     DUTot : Obtained incremental displacements (gas)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTot
          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID, IParticle
          integer(INTEGER_TYPE) :: ElementID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement
          
           do IParticle = 1, Counters%NParticles ! Loop over particles

            ElementID = ElementIDArray(IParticle)
          
            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle)
            else
              iEntity = 1
            end if 
              
            ParticleDisplacement = 0.0
                    
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                NodeID = iabs(ElementConnectivities(INode, ElementID)) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * DUTot(DofID,iEntity)                             
              end do                                             
            end do 
                          
            ! Total particle displacement
            Particles(IParticle)%UG = Particles(IParticle)%UG + ParticleDisplacement
        
         end do ! particles
          
         end subroutine UpdateParticleDisplacementsGas
         
          subroutine UpdateParticleVelocityWater(DUTot)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements (water)
        !
        !     DUTot : Obtained incremental displacements (water)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        
     
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTot          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID, IParticle
          integer(INTEGER_TYPE) :: ElementID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement
          
           do IParticle = 1, Counters%NParticles ! Loop over particles

            ElementID = ElementIDArray(IParticle)
          
            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle) 
            else
              iEntity = 1
            end if 
              
            ParticleDisplacement = 0.0
                    
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                NodeID = iabs(ElementConnectivities(INode, ElementID)) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * DUTot(DofID,iEntity)                             
              end do                                             
            end do 
                          
            ! Total particle velocity
            
            VelocityWaterArray(IParticle,:) = ParticleDisplacement 
            
                   
         end do ! particles
          
         end subroutine UpdateParticleVelocityWater
         
         subroutine UpdateParticleVelocityGas(DUTot)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements (gas)
        !
        !     DUTot : Obtained incremental displacements (gas)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        
     
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTot
          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID, IParticle
          integer(INTEGER_TYPE) :: ElementID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement
          
           do IParticle = 1, Counters%NParticles ! Loop over particles

            ElementID = ElementIDArray(IParticle)
          
            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle) 
            else
              iEntity = 1
            end if 
              
            ParticleDisplacement = 0.0
                    
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                NodeID = iabs(ElementConnectivities(INode, ElementID)) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * DUTot(DofID,iEntity)                             
              end do                                             
            end do 
                          
            ! Total particle velocity
            
            VelocityGasArray(IParticle,:) = ParticleDisplacement
                   
         end do ! particles
          
         end subroutine UpdateParticleVelocityGas
         
         
 
         subroutine InitialiseConvectivePhaseData()   
        !**********************************************************************
        !    SUBROUTINE: InitialiseConvectivePhaseData
        ! 
        !    DESCRIPTION: 
        !>   Initialises the arrays globally defined in this module
        !
        !**********************************************************************     
           implicit none
         
           ! Local variable
           integer(INTEGER_TYPE) :: IError
           
           call DestroyConvectivePhaseData()
           
           allocate(TemporaryMappingVector(Counters%N), stat = IError)
           TemporaryMappingVector = 0.0
         
         end subroutine InitialiseConvectivePhaseData
         
   
         subroutine DestroyConvectivePhaseData()    
         !**********************************************************************
         !    SUBROUTINE: DestroyConvectivePhaseData
         ! 
         !    DESCRIPTION: 
         !>   Destroys the arrays globally defined in this module
         !
         !**********************************************************************   
           implicit none
         
           ! Local variable
           integer(INTEGER_TYPE) :: IError
           
           if (allocated(TemporaryMappingVector) ) then
             deallocate(TemporaryMappingVector, stat = IError)
           end if
         
    end subroutine DestroyConvectivePhaseData

 
      subroutine CalculateCriticalTimeStep()
      !**********************************************************************
      !
      !  Function: calculate critical time step size
      !            1-phase: based on CFL condition
      !            2-phase: based on MSc Thesis Mieremet (TU Delft, 2015)
      !            mass scaling is taken into account
      !
      !**********************************************************************
      
      implicit none
      
        !local variable
        integer(INTEGER_TYPE) :: IAEl, IEl, NParticles, ParticleIndex, IParticle, J, MaterialIndex
        real(REAL_TYPE) :: WaveSpeed, MaxWaveSpeed, TimeStep, MinTimeStep, CourantFactor
        real(REAL_TYPE) :: ConstDensityLiquid, ConstDensitySolid, BulkW, Nu, G
        real(REAL_TYPE) :: WaveSpeed_c1, WaveSpeed_c2, WaveSpeed_c3
        real(REAL_TYPE) :: rho_tilde, Timestep_3
        real(REAL_TYPE) :: TimeIncrementNew
        real(REAL_TYPE) :: Chi = 0.0
        
        !--------------------------------------------------------------------
        !-------------------- 1 layer formulation ---------------------------
        !--------------------------------------------------------------------
        
        if((NFORMULATION==1).or. &
          ((NFORMULATION==2).and.(CalParams%NumberOfPhases==1)))then
          TimeStep = 0.0
          MinTimeStep = huge(MinTimeStep)
          CourantFactor = CalParams%CourantNumber
         
          do IAEl = 1, Counters%NAEl ! loop over all active elements
            IEl = ActiveElement(IAEl)
            NParticles = NPartEle(IEl) 
            
            Chi = 0.0
            WaveSpeed = 0.0
            MaxWaveSpeed = 0.0
            
            do IParticle = 1, NParticles ! loop over all particles in element
              ParticleIndex = GetParticleIndex(IParticle, IEl)
			  MaterialIndex = MaterialIDArray(ParticleIndex)
			  if (.not.MatParams(MaterialIndex)%MaterialModel==ESM_RIGID_BODY) then ! do not compute wavespeed for a rigid body
              call GetWaveSpeedForTimeStepSize(ParticleIndex, IEl, WaveSpeed)
              
              if (IsNan(WaveSpeed).or.(abs(WaveSpeed) > huge(WaveSpeed))) then
                call GiveError('WaveSpeed is not a number or infinity')
              end if
              
              if (WaveSpeed>MaxWaveSpeed) MaxWaveSpeed = WaveSpeed
			  end if
            end do ! end loop over all particles in element
            
            if (MaxWaveSpeed>(0.0)) then
              TimeStep = ElementLMin(IEl) / MaxWaveSpeed
              
              if (CalParams%ApplyBulkViscosityDamping) then
                Chi = CalParams%BulkViscosityDamping1 - CalParams%BulkViscosityDamping2**2 * TimeStep * RateVolStrain(IEl)
                TimeStep = (sqrt(1 + Chi**2) - Chi) * TimeStep
              end if
              
			end if
			
			if (TimeStep>0.0) then
				if (TimeStep<MinTimeStep) MinTimeStep = TimeStep
			end if
               
          end do ! end loop over all active elements
          
        end if ! end 1 layer formulation
        
        !------------------------------------------------------------------------- 
        !--------------- 2 layer formulation and 2 phases ------------------------
        !-------------------------------------------------------------------------
        
        if ((NFORMULATION==2).and.(CalParams%NumberOfPhases==2)) then
            
          !--- Calculate  Soil and Water parameters  ----
          ConstDensityLiquid = 0.0
          ConstDensitySolid = 0.0
          BulkW  = 0.0
          Nu = 0.0
          G = 0.0
                
          do J = 1, Counters%NLayers ! loop over all materials
            if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
              ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
              ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
              BulkW  = MatParams(J)%BulkModulusLiquid
              Nu  = MatParams(J)%PoissonRatio
              G = MatParams(J)%ShearModulus
            else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
              ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
              BulkW  = MatParams(J)%BulkModulusLiquid 
            else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-solid'.or.MatParams(j)%MaterialPhases=='1-phase-solid') then
              ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
              Nu  = MatParams(J)%PoissonRatio
              G = MatParams(J)%ShearModulus
            end if
          end do 
             
          TimeStep = 0.0
          MinTimeStep = huge(MinTimeStep)
          CourantFactor = CalParams%CourantNumber

          do IAEl = 1, Counters%NAEl ! loop over all active elements
            IEl = ActiveElement(IAEl)
            NParticles = NPartEle(IEl)
            
            WaveSpeed_c1 = 0.0
            WaveSpeed_c2 = 0.0
            WaveSpeed_c3 = 0.0
            MaxWaveSpeed = 0.0
            rho_tilde = 0.0
            Chi = 0.0

            do IParticle = 1, NParticles ! loop over all particles in element
              ParticleIndex = GetParticleIndex(IParticle, IEl)
              call GetWaveSpeed2LayerForm (ParticleIndex, IEl, ConstDensityLiquid, ConstDensitySolid, &
                                           BulkW, Nu, G, WaveSpeed_c1, WaveSpeed_c2)
              if (IsNan(WaveSpeed_c1).or.(abs(WaveSpeed_c1) > huge(WaveSpeed_c1))) then
                call GiveError('WaveSpeed_c1 is not a number or infinity'       //NEW_LINE('A')// &
                               'ParticleIndex: ' // trim(String(ParticleIndex)) //NEW_LINE('A')// &
                               'Element      : ' // trim(String(IEl)))
              end if

              if (IsNan(WaveSpeed_c2).or.(abs(WaveSpeed_c2) > huge(WaveSpeed_c2))) then
                call GiveError('WaveSpeed_c2 is not a number or infinity'       //NEW_LINE('A')// &
                               'ParticleIndex: ' // trim(String(ParticleIndex)) //NEW_LINE('A')// &
                               'Element      : ' // trim(String(IEl)))
              end if

              if(WaveSpeed_c1>MaxWaveSpeed) MaxWaveSpeed = WaveSpeed_c1
              if(WaveSpeed_c2>MaxWaveSpeed) MaxWaveSpeed = WaveSpeed_c2

              ! ---------- Permeability criterion ------------!
              if((TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLIDLIQUID).and. &
                 (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid))then

                rho_tilde = (Particles(ParticleIndex)%EffConcentrationRatioSolid * ConstDensitySolid +  &
                            Particles(ParticleIndex)%EffPorosity * ConstDensityLiquid + &
                            (1.0 / Particles(ParticleIndex)%EffPorosity - 2.0) * ConstDensityLiquid)*CalParams%ScalingMassFactor

                Timestep_3 = 2.0 * rho_tilde * Particles(ParticleIndex)%Conductivity /  &
                             (ConstDensityLiquid * CalParams%GravityData%GAccel)

                WaveSpeed_c3 = ElementLMin(IEl) / Timestep_3

                if (IsNan(WaveSpeed_c3).or.(abs(WaveSpeed_c3) > huge(WaveSpeed_c3))) then
                  call GiveError('WaveSpeed_c3 is not a number or infinity'       //NEW_LINE('A')// &
                                 'ParticleIndex: ' // trim(String(ParticleIndex)) //NEW_LINE('A')// &
                                 'Element      : ' // trim(String(IEl)))

                end if

                if(WaveSpeed_c3>MaxWaveSpeed) MaxWaveSpeed = WaveSpeed_c3

              end if
            end do ! end loop over all particles in element
            
            if (MaxWaveSpeed>(0.0)) then
               TimeStep = ElementLMin(IEl) / MaxWaveSpeed
               if (CalParams%ApplyBulkViscosityDamping) then
                 Chi = CalParams%BulkViscosityDamping1 - CalParams%BulkViscosityDamping2**2 * TimeStep * RateVolStrain(IEl)
                 TimeStep = (sqrt(1 + Chi**2) - Chi) * TimeStep
               end if
            end if

            if (TimeStep<MinTimeStep) MinTimeStep = TimeStep

         end do ! end loop over all active elements

        end if ! 2 layer formulation and 2 phases

        TimeIncrementNew = MinTimeStep * CourantFactor

        if (.not.CalParams%ApplyQuasiStatic) then
          if ((CalParams%TotalRealTime + TimeIncrementNew)>CalParams%TotalTime) then
            TimeIncrementNew = CalParams%TotalTime - CalParams%TotalRealTime
          end if
        end if

        ! hardcoded minimum time step size, to prevent time steps smaller than 1d-10 s
        if (TimeIncrementNew <= 1d-10) then
          TimeIncrementNew = 1d-10
        end if

        if (CalParams%TimeStep == 1 .or. (CalParams%IStep == 1.and.CalParams%TimeStep == 2)) then
          ! only for the first time step
          if (abs(TimeIncrementNew - CalParams%TimeIncrement) > SMALL) then
            call GiveMessage(' MaxWaveSpeed : ' // trim(String(MaxWaveSpeed)))
            call GiveMessage(' MinTimeStep  : ' // trim(String(MinTimeStep)))
            call GiveMessage(' TimeIncrement: ' // trim(String(TimeIncrementNew)))
          endif
        endif

        CalParams%TimeIncrement = TimeIncrementNew

      end subroutine CalculateCriticalTimeStep




      subroutine GetWaveSpeed2LayerForm(IParticle, IEl, ConstDensityLiquid, ConstDensitySolid, &
                                        BulkW, Nu, G, WaveSpeed_c1, WaveSpeed_c2)
      !**********************************************************************
      !
      !  Function: calculate wave speed for particle IParticle in element IEl
      !
      !********************************************************************** 
      
      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: IParticle, IEl
        real(REAL_TYPE), intent(in) :: ConstDensityLiquid, ConstDensitySolid, BulkW, Nu, G
        real(REAL_TYPE), intent(inout) :: WaveSpeed_c1, WaveSpeed_c2

        ! local variables
        real(REAL_TYPE) :: Eoed_Skel, Eoed, rho_sat, Beta_s

        WaveSpeed_c1 = 0.0
        WaveSpeed_c2 = 0.0
        rho_sat = 0.0
        Beta_s = 0.0

        !!=========== ELEMENT CONTAINS = SOLID+LIQUID =====================!!
        if (TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLIDLIQUID) then

          if (MaterialPointTypeArray(IParticle)==MaterialPointTypeSolid) then !this takes a lot of time!

            if (Particles(IParticle)%PhaseStatus==PhaseStatusSOLID) then

              Eoed_Skel = (2 * G * (1 - Nu)) / (1 - 2 * Nu) ! kPa  ! Eoed

              ! Assumed to be fully saturated :: %EffPorosity = %ConcRatioLiquid
              Eoed = Eoed_Skel + BulkW / Particles(IParticle)%EffPorosity
              rho_sat = Particles(IParticle)%EffConcentrationRatioSolid * ConstDensitySolid +  &
                        Particles(IParticle)%EffPorosity * ConstDensityLiquid

              Beta_s = sqrt((Particles(IParticle)%EffPorosity * Eoed_Skel / BulkW) /  &
                            (Particles(IParticle)%EffConcentrationRatioSolid +  &
                             Particles(IParticle)%EffPorosity * Eoed_Skel / BulkW))
            end if

            if (Particles(IParticle)%PhaseStatus==PhaseStatusLIQUID) then
              Beta_s = 1.0
            end if

          end if ! Solid Material point

          if (MaterialPointTypeArray(IParticle)==MaterialPointTypeLiquid) then 
            Beta_s = 1.0
          end if ! Liquid Material point
        end if
           
        !!=========== ELEMENT CONTAINS = LIQUID =====================!!
        if (TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) then
          Beta_s = 1.0
        end if

        !!=========== ELEMENT CONTAINS = SOLID =====================!!
        if (TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLID) then
          Eoed_Skel = ( 2 * G * (1 - Nu))/ (1 - 2 * Nu) ! kPa  ! Eoed
          Eoed = Eoed_Skel 
          rho_sat = Particles(IParticle)%EffConcentrationRatioSolid * ConstDensitySolid
          Beta_s = 0.0
        end if

        !--------- Calculate wave speed --------------!
        if(rho_sat>0.0) then
          WaveSpeed_c1 = sqrt((Eoed / rho_sat)/CalParams%ScalingMassFactor)   
        else
          WaveSpeed_c1 = 0.0
        end if
   
        if(ConstDensityLiquid>0.0) then
          WaveSpeed_c2 = Beta_s * sqrt((BulkW / ConstDensityLiquid)/CalParams%ScalingMassFactor)
        else
          WaveSpeed_c2 = 0.0
        end if
         
      end subroutine GetWaveSpeed2LayerForm

        
      end module ModDYNConvectivePhase
