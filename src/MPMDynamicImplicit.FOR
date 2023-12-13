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
	  
	  
	  module ModDynamicImplicit
      !**********************************************************************
      !
      !    Function:  This module contains the routine for running a dynamic
      !               load step using explicit time integration.
      !
      !     $Revision: 9794 $
      !     $Date: 2022-09-20 15:20:25 +0200 (di, 20 sep 2022) $
      !
      !**********************************************************************

      use ModReadCalculationData
      use ModReadMaterialData
      use ModWriteTestData
      use ModMPMData
      use ModMeshInfo
      use ModRotBoundCond
      use ModMPMMeshAdjustment
      use ModReadMPMData
      use ModElementEvaluation
      use ModDYNConvectivePhase
      use ModMPMDynContact
      use ModCounters
      use ModEmptyElements
      use ModLiquid
      use ModWriteVTKOutput
      use ModParticle
      use ModWriteResultData
      use ModTwoLayerFormulation
      use ModGlobalConstants
      
      implicit none

      contains ! Routines of this module

        subroutine RunImplicitDynamicLoadStep()
        !**********************************************************************
        !
        !  Function:  Routine called from the main routine for performing a
        !             dynamic load step using explicit time integration.
        !             Contains the 'time step loop' of explicit integration.
        !
        !**********************************************************************

        implicit none
        
        ! Local variables
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: &
            Momentum
        
        ! Integer variable 
        integer(INTEGER_TYPE) :: IIteration_Newton, EntityUsed, ParticleIndex
        
        ! PBoundary variables 
        integer(INTEGER_TYPE) :: IDOF, J
        
        
        
        EntityUsed = 1 

          !********** 4b - TIME STEP LOOP ******************************
          do while( (.not.CalParams%ConvergenceCheck%DoesConverge) .and. (.not.CalParams%ConvergenceCheck%DoesDiverge) )
            
            ! Increase time
            CalParams%TimeStep = CalParams%TimeStep + 1
            
            !if (CalParams%IStep .ge. 99) then 
            !    
            !    CalParams%ImplicitIntegration%HH_DirDerivParameter_SK04 = 1e-5
            !    
            !end if 
            
            CalParams%TotalRealTime = CalParams%TotalRealTime + CalParams%TimeIncrement
            CalParams%OverallRealTime = CalParams%OverallRealTime + CalParams%TimeIncrement
            call WriteInLogFile('TimeStep: ' // trim(String(CalParams%TimeStep))//' '//trim(String(CalParams%TimeIncrement)) &
                                                                                //' '//trim(String(CalParams%TotalRealTime)))
            call WriteInLogFile('  Skipconvection? '//trim(String(MinimumDeterminantRatioReached))//' '// &
                                                      trim(String(IsMPMSkipConvection())))
            
            call UpdateMultipliersForTimeDepencency() 
           
            
            
            ! INITIALISATION ------------------------------------------------------------------------------------
            ! SUBROUTINE #1: Map the mass from the MPs to the nodes using shape functions ShapeValuesArray.  
            call MapMassMP2NodesInitialization_ImplicitIter(LumpedMassDry) 
            
            ! SUBROUTINE #2: Map the momentum (MassArray * VelocityArray) of the MPs to the nodes in a local variable, Momentum.  
            !                Momentum is a local variable as we only care about velocities.             
            call MapMomentumMP2NodesInitialization_ImplicitIter(Momentum) 
            
            ! SUBROUTINE #3: Sulsky et al. (1995) implementation where we divide by Momentum by LumpedMassDry to find the nodal velocity TotalVelocitySoil. 
            !                Initial velocity TotalVelocitySoil can be a triggering mechanism.  
            call GetNodalVelocityFromNodalMomentum_ImplicitIter(Momentum, LumpedMassDry, TotalVelocitySoil)
            
            ! SUBROUTINE #4: Compute IntLoad based on the initial stress of the MPs SigmaEffArray. 
            call GetNodalIntForcesBefore_ImplicitIter(IntLoad, SigmaEffArray) 
            
            ! SUBROUTINE #5: Compute newly introduced traction forces, ExtLoad_ImplicitIter, and gravity forces, GravityLoad_ImplicitIter.
            !                These newly introduced external forces can be triggering mechanism.
            call GetNodalExternalForces(ExtLoad, GravityLoad)
            !call GetNodalExternalForces(ExtLoad_ImplicitIter, GravityLoad_ImplicitIter) ! 
            
            ! SUBROUTINE #6: Guess the new velocity, TotalVelocitySoil_ImplicitIter, based on the initial velocity, TotalVelocitySoil.  
            call GetFirstGuessOfVelocity_ImplicitIter(TotalVelocitySoil_ImplicitIter, TotalVelocitySoil)
            
            ! SUBROUTINE #7: Prescribe velocity at the nodes. This step affects TotalVelocitySoil_ImplicitIter. 
            call ApplyMPPrescribedVelocity_ImplicitIter(TotalVelocitySys) ! prescribe and map MP velocity to the nodes. 
            call ApplyNodalPrescribedVelocity_ImplicitIter(TotalVelocitySys) ! prescribe velocity to the nodes. 
            
            ! NEWTON LOOP ------------------------------------------------------------------------------------
            
            do IIteration_Newton = 1, CalParams%ImplicitIntegration%MaxNewtonIterations ! Newton iterations
            
            ! SUBROUTINE #8: Compute IncrementalDisplacementSoil, the strain increment, Eps, and store in 
            !                Particles(ParticleIndex)%EpsStep_ImplicitIter. 
            call GetMPStrainFromNodeVelocity_ImplicitIter(TotalVelocitySoil_ImplicitIter)
            
            ! SUBROUTINE #9: Compute SigmaEffArray_ImplicitIter using the constitutive model where DeltaSigma = DeltaEps * D
            !                (the water pressure and the gass pressure are updated - mass balance equations are solved). 
            call MPMDYNGetSig_ImplicitIter()
           
            ! SUBROUTINE #10: Update the weight of the point, Particles(IntGlo)%IntegrationWeight_ImplicitIter, in the physical space based 
            !                 on the volumetric strain.  
            call DynUpdateParticleWeights_ImplicitIter()
            
            ! SUBROUTINE #11: Update the particle position, GlobPosArray_ImplicitIter, local positions, Particles(IParticle)%LocPos_ImplicitIter, 
            !                 and shape functions, ShapeValuesArray_ImplicitIter. 
            call UpdateParticlePos_ImplicitIter()
           
            ! SUBROUTINE #12: Update particle housekeeping information.   
            call UpdateParticleHouseKeeping_ImplicitIter() !-> uses EleParticles to updated the rest of the housekeeping parameters
            call SetActiveElement_ImplicitIter()
            call SetParticleIndex_ImplicitIter()
            call CheckFillingOfElements()
            call SetUpEntityElements_ImplicitIter()
            call SetUpMaterialElements_ImplicitIter()
           
            ! SUBROUTINE #13: Compute new internal force, IntLoad_ImplicitIter, based on the newly computed stress SigmaEffArray_ImplicitIter
            call GetNodalIntForcesAfter_ImplicitIter(IntLoad_ImplicitIter, SigmaEffArray_ImplicitIter)  
            
            ! SUBROUTINE #: Map mass to the nodes after updating the position of the MPs using ShapeValuesArray_ImplicitIter
            call MapMassMP2Nodes_ImplicitIter(LumpedMassDry_ImplicitIter)
            
            call GetNodalExternalForces_ImplicitIter(ExtLoad_ImplicitIter, GravityLoad_ImplicitIter)
           
            ! SUBROUTINE #14: Compute internal force, IntLoad_ImplicitIter, based on prescribe velocity at the nodes. 
            !                 This step affects IntLoad_ImplicitIter. 
            call ApplyMPPrescribedVelocity_IntLoad(IntLoad_ImplicitIter)
            call ApplyNodalPrescribedVelocity_IntLoad(IntLoad_ImplicitIter)
            
            
            
            ! SUBROUTINE #15: Calculate the implicit residual based on the momentum balance equation. 
            call GetImplicitResidualEquation_SK04(LumpedMassDry_ImplicitIter, LumpedMassDry_ImplicitIter, &  !LumpedMassDry_ImplicitIter, LumpedMassDry, &  
                                                        ExtLoad_ImplicitIter, ExtLoad, & 
                                                    GravityLoad_ImplicitIter, GravityLoad, &  
                                                        IntLoad_ImplicitIter, IntLoad, & 
                                              TotalVelocitySoil_ImplicitIter, TotalVelocitySoil, &             
                                                    ImplicitResidual_ImplicitIter)
                
            ! SUBROUTINE #16: The below subroutine is not necessary. Make sure that ImplicitResidual_ImplicitIter is zero 
            !                 at the nodes where velocity is presrcibed. 
            !                 We need to apply the prescribed velocity boundary condition in a way that we 
            !                 maintain the prescribed velocity to be that value 
            !                 if we are applying that prescribed velocity there then the residual implicit 
            !                 there should be zero after one step. BECAUSE IT SHOULDN'T CHANGE IF CONSTANT VELOCITY.
            !                 Directional derivative should also be zero at those nodes where we are prescribing velocity. 
            !                 The rest should be fine. 
            !                 Then we just map the velocities appropriately to the nodes    
            call ApplyNodalPrescribedVelocity_ZeroVector(ImplicitResidual_ImplicitIter)

            ! SUBROUTINE #17: Calculate the normal of the ImplicitResidual_ImplicitIter 
            ImplicitResidualMagnitude_ImplicitIter = norm2(ImplicitResidual_ImplicitIter)
                
            ! SUBROUTINE #18: Check error tolerance is below  
            if (ImplicitResidualMagnitude_ImplicitIter < CalParams%ImplicitIntegration%NewtonErrorThreshold) then 
                !... we reached the solution
                CalParams%ImplicitIntegration%IsEnterGMRES = .false. 
                EXIT
            else 
                CalParams%ImplicitIntegration%IsEnterGMRES = .true.
            end if     
            
            ! SUBROUTINE #19: GMRES algorithm
            if (CalParams%ImplicitIntegration%IsEnterGMRES == .true.) then 
                
                ! SUBROUTINE #20: Caluclate TotalVelocitySoil_ImplicitIter and IntLoad_ImplicitIter
                call GMRES_SK04()
            
                ! SUBROUTINE #21: Make sure that the fixities are applied according to PBoundary
                do IDOF = 1, Counters%N ! loop over all degrees of freedom
                    do J = 1, Counters%nEntity ! loop over all entities
                        TotalVelocitySoil_ImplicitIter(IDOF,J) = TotalVelocitySoil_ImplicitIter(IDOF,J) * PBoundary(IDOF)
                    end do ! loop over all entities
                end do ! loop over all degrees of freedom
            
            end if            
            
            end do !CalParams%ImplicitIntegration%MaxNewtonIterations
            
            
            ! SUBROUTINE #22: Map velocity, TotalVelocityArray_ImplicitIter, from nodes to material points, VelocityArray, using 
            !                 the momentum formulation. 
            call UpdateParticleVelocityAndMapMomentum_ImplicitIter(Momentum,TotalVelocitySoil_ImplicitIter)
            
            ! SUBROUTINE #23: Increase particle strain, Particles()%Eps, by the implicit strain  increment, Particle%EpsStep_ImplicitIter. 
            call IncreaseParticleStrains_ImplicitIter()
            
            !call UpdateParticlePos_ImplicitIter()
            
            ! SUBROUTINE #24: Transfer variables from ImplicitIter to regular variables. This is important if 0<mu<1 and 
            !                  for VTK outputs. 
            call ImplicitVariablesTransfer_ImplicitIter_()
            
            ! SUBROUTINE #25: Update Sigma0EffArray = SigmaEffArray and other stress steps.
            call SetInitialStressForNextLoadStep()

            ! Optional writing of time step results
            ! write intermediate time step results and if ApplyQuickCheckOutput
            call WriteTimeStepResults(.false.)

            ! Convergence check
            CalParams%ConvergenceCheck%DoesConverge = ConvergenceCheck()

            ! Divergence check
            if (CalParams%ConvergenceCheck%ApplyDivergenceCheck) then
              CalParams%ConvergenceCheck%DoesDiverge = DivergenceCheck()

              if (CalParams%ConvergenceCheck%DoesDiverge) then
                call GiveWarning('Calculation is diverging...')
              end if
            end if

              if (CalParams%OutputDebugData) then
            ! Time step output
            if (NFORMULATION==1) then
              call EnergyOutput()
            else
              call EnergyOutput2LayForm()
            end if

            call MaterialPointOutput() ! only if OutputNumberParticles>0
            
            !compute and write on file the sum of nodal reactions on selected surfaces
            call SurfaceReactionOutput() 
            
             end if

            ! update ElementLMin in case of ULFEM
            if (IsULFEMComputation()) then
              call DetermineElementLMin() ! calulate minimum element altitude
            end if

            call CalculateCriticalTimeStep()
            call WriteFEMNodeData() ! only if .not.IsMPMComputation

          end do ! ----- end TIME STEP LOOP -----

        end subroutine RunImplicitDynamicLoadStep
        
        
        
        
        
        
        
            !-----------------------------------------------------------------------------------------------------------------
            subroutine ImplicitVariablesTransfer_ImplicitIter_() 
            !-------------------------------------------------------
            ! Function: to transfer variables from ImplicitIter variables 
            !           to regular variables.  
            ! 
            !-------------------------------------------------------
            
            
            implicit none 
            
            integer(INTEGER_TYPE) :: ParticleIndex
            
            TotalVelocitySoil = TotalVelocitySoil_ImplicitIter
            IntLoad = IntLoad_ImplicitIter
            ExtLoad = ExtLoad_ImplicitIter
            GravityLoad = GravityLoad_ImplicitIter
            
            SigmaEffArray = SigmaEffArray_ImplicitIter
            
            ESMStatevArray = ESMStatevArray_ImplicitIter
            GlobPosArray = GlobPosArray_ImplicitIter
            ElementIDArray = ElementIDArray_ImplicitIter
            EleParticles = EleParticles_ImplicitIter
            ShapeValuesArray = ShapeValuesArray_ImplicitIter
            DShapeValuesArray = DShapeValuesArray_ImplicitIter
            
            NPartEle = NPartEle_ImplicitIter
            EleParticlesHelp = EleParticlesHelp_ImplicitIter
            IsActiveElement = IsActiveElement_ImplicitIter
            ActiveElement = ActiveElement_ImplicitIter
            ActiveNode = ActiveNode_ImplicitIter
            GetParticleIndex = GetParticleIndex_ImplicitIter
            EntityElements = EntityElements_ImplicitIter
            MaterialElements = MaterialElements_ImplicitIter
            
            Counters%NAEl = Counters%NAEl_ImplicitIter
            TotalVelocitySoil = TotalVelocitySoil_ImplicitIter
            
            ! need to loop over particles 
            do ParticleIndex = 1, Counters%NPARTICLES
            Particles(ParticleIndex)%EpsStep = Particles(ParticleIndex)%EpsStep_ImplicitIter
            Particles(ParticleIndex)%IntegrationWeight = Particles(ParticleIndex)%IntegrationWeight_ImplicitIter
            Particles(ParticleIndex)%LocPos = Particles(ParticleIndex)%LocPos_ImplicitIter
            
            UArray(ParticleIndex,:) = UArray(ParticleIndex,:) + UStepArray(ParticleIndex,:)
            
            end do
            
            
            
            end subroutine ImplicitVariablesTransfer_ImplicitIter_
            !-----------------------------------------------------------------------------------------------------------------
            
        
        
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        ! SUBROUTINE #1
        subroutine MapMassMP2NodesInitialization_ImplicitIter(LumpedMass)
        !**********************************************************************
        !
        ! Function:  to form the lumped mass vector and map momentum from material points to nodes
        !
        ! O  Momentum : nodal momentum vector, the output of this subroutine
        ! O  LumpedMass : lumped soil mass
        !
        !**********************************************************************

        implicit none
        
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: Momentum
          
          
          
          ! Local variables
          !real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(ELEMENTNODES) :: PartShape
          real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity) :: LumMassX
          real(REAL_TYPE) :: ParticleMass
          
          integer(INTEGER_TYPE) :: I, IEl, IPart, IAEl, ParticleIndex, iEntity, NoEn

          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: Ni ! number of columns is equal to the number of nodes of the element
          !real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: Mom
          
			
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES) :: ShapeValuesArray !, intent(in)
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: LumpedMass
          !real(REAL_TYPE), dimension(Counters%NParticles), intent(in) :: MassArray
         
          !integer(INTEGER_TYPE), dimension(Counters%NEl) :: ActiveElement !, intent(in)
          
          LumpedMass = 0.0
          !Momentum = 0.0
          NoEn = Counters%nEntity
          
          do IAEl = 1, Counters%NAEl ! loop over all elements

              IEl = ActiveElement(IAEl)
                
              do I = 1, NVECTOR
                  Ni(I, 1:ELEMENTNODES) = ReducedDof( ElementConnectivities(1:ELEMENTNODES, IEl) ) + I
                  !Mom(I, 1:ELEMENTNODES, 1:Counters%nEntity) = Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn)
              end do
              LumMassX = LumpedMass( Ni(1, 1:ELEMENTNODES), 1:NoEn)
                
              do IPart = 1, NPartEle(IEl) ! loop over all material points in element
                  ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID

                  if((NFORMULATION==1).or.(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid)) then ! Only if NumberOfLayers = 1 or SOLID Material Point

                      !ParticleVelocity = VelocityArray(ParticleIndex,:) !get particle velocity vector
                      ParticleMass = MassArray(ParticleIndex) !get particle's mass
                      PartShape = ShapeValuesArray(ParticleIndex,:) !get particle shape fucntions

                      !if (CalParams%ApplyContactAlgorithm) then
                      !    iEntity = EntityIDArray(ParticleIndex)
                      !else
                          iEntity = 1
                      !end if

                      ! nodal x-mass
                      LumMassX(:, iEntity) = LumMassX(:, iEntity) + ParticleMass * PartShape(:)

                      ! nodal momentum in directions from 1 to NDIM
                      !do I=1,NVECTOR
                      !    Mom(I, :, iEntity) = Mom(I, :, iEntity) + ParticleMass * PartShape(:) * ParticleVelocity(I)
                      !end do

                  end if ! Only if NumberOfLayers = 1 or SOLID Material Point

              end do ! loop over material points in element

              do I = 1, NVECTOR
                  LumpedMass( Ni(I, 1:ELEMENTNODES), 1:NoEn ) = LumMassX( 1:ELEMENTNODES, 1:NoEn )
                  !Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn ) = Mom( I, 1:ELEMENTNODES, 1:NoEn )
              end do

          end do ! loop over Eelements
          
        
        
        end subroutine MapMassMP2NodesInitialization_ImplicitIter
        
        !-------------------------------------------------------------------------------------------
        
        
        subroutine MapMassMP2Nodes_ImplicitIter(LumpedMass)
        !**********************************************************************
        !
        ! Function:  to form the lumped mass vector and map momentum from material points to nodes
        !
        ! O  Momentum : nodal momentum vector, the output of this subroutine
        ! O  LumpedMass : lumped soil mass
        !
        !**********************************************************************

        implicit none
        
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: Momentum
          
          
          
          ! Local variables
          !real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(ELEMENTNODES) :: PartShape
          real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity) :: LumMassX
          real(REAL_TYPE) :: ParticleMass
          
          integer(INTEGER_TYPE) :: I, IEl, IPart, IAEl, ParticleIndex, iEntity, NoEn

          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: Ni ! number of columns is equal to the number of nodes of the element
          !real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: Mom
          
			
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES) :: ShapeValuesArray !, intent(in)
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: LumpedMass
          !real(REAL_TYPE), dimension(Counters%NParticles), intent(in) :: MassArray
         
          !integer(INTEGER_TYPE), dimension(Counters%NEl) :: ActiveElement !, intent(in)
          
          LumpedMass = 0.0
          !Momentum = 0.0
          NoEn = Counters%nEntity
          
          do IAEl = 1, Counters%NAEl_ImplicitIter ! loop over all elements
            
              IEl = ActiveElement_ImplicitIter(IAEl)
                
              do I = 1, NVECTOR
                  Ni(I, 1:ELEMENTNODES) = ReducedDof( ElementConnectivities(1:ELEMENTNODES, IEl) ) + I
                  !Mom(I, 1:ELEMENTNODES, 1:Counters%nEntity) = Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn)
              end do
              LumMassX = LumpedMass( Ni(1, 1:ELEMENTNODES), 1:NoEn)
                
              do IPart = 1, NPartEle_ImplicitIter(IEl) ! loop over all material points in element
                  ParticleIndex = GetParticleIndex_ImplicitIter(IPart, IEl) ! Get the particle ID

                  if((NFORMULATION==1).or.(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid)) then ! Only if NumberOfLayers = 1 or SOLID Material Point

                      !ParticleVelocity = VelocityArray(ParticleIndex,:) !get particle velocity vector
                      ParticleMass = MassArray(ParticleIndex) !get particle's mass
                      PartShape = ShapeValuesArray_ImplicitIter(ParticleIndex,:) !get particle shape fucntions
                        
                      !if (CalParams%ApplyContactAlgorithm) then
                      !    iEntity = EntityIDArray(ParticleIndex)
                      !else
                          iEntity = 1
                      !end if

                      ! nodal x-mass
                      LumMassX(:, iEntity) = LumMassX(:, iEntity) + ParticleMass * PartShape(:)

                      ! nodal momentum in directions from 1 to NDIM
                      !do I=1,NVECTOR
                      !    Mom(I, :, iEntity) = Mom(I, :, iEntity) + ParticleMass * PartShape(:) * ParticleVelocity(I)
                      !end do

                  end if ! Only if NumberOfLayers = 1 or SOLID Material Point

              end do ! loop over material points in element

              do I = 1, NVECTOR
                  LumpedMass( Ni(I, 1:ELEMENTNODES), 1:NoEn ) = LumMassX( 1:ELEMENTNODES, 1:NoEn )
                  !Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn ) = Mom( I, 1:ELEMENTNODES, 1:NoEn )
              end do
                
          end do ! loop over Eelements
          
        
        
        end subroutine MapMassMP2Nodes_ImplicitIter
        !-------------------------------------------------------------------------------------------
        
        
        subroutine MapMassMP2Nodes_DirDeriv_ImplicitIter(LumpedMass)
        !**********************************************************************
        !
        ! Function:  to form the lumped mass vector and map momentum from material points to nodes
        !
        ! O  Momentum : nodal momentum vector, the output of this subroutine
        ! O  LumpedMass : lumped soil mass
        !
        !**********************************************************************

        implicit none
        
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: Momentum
          
          
          
          ! Local variables
          !real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(ELEMENTNODES) :: PartShape
          real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity) :: LumMassX
          real(REAL_TYPE) :: ParticleMass
          
          integer(INTEGER_TYPE) :: I, IEl, IPart, IAEl, ParticleIndex, iEntity, NoEn

          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: Ni ! number of columns is equal to the number of nodes of the element
          !real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: Mom
          
			
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES) :: ShapeValuesArray !, intent(in)
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: LumpedMass
          !real(REAL_TYPE), dimension(Counters%NParticles), intent(in) :: MassArray
         
          !integer(INTEGER_TYPE), dimension(Counters%NEl) :: ActiveElement !, intent(in)
          
          LumpedMass = 0.0
          !Momentum = 0.0
          NoEn = Counters%nEntity
          
          do IAEl = 1, Counters%NAEl_DirDeriv_ImplicitIter ! loop over all elements
        
              IEl = ActiveElement_DirDeriv_ImplicitIter(IAEl)

              do I = 1, NVECTOR
                  Ni(I, 1:ELEMENTNODES) = ReducedDof( ElementConnectivities(1:ELEMENTNODES, IEl) ) + I
                  !Mom(I, 1:ELEMENTNODES, 1:Counters%nEntity) = Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn)
              end do
              
              LumMassX = LumpedMass( Ni(1, 1:ELEMENTNODES), 1:NoEn)
                
              do IPart = 1, NPartEle_DirDeriv_ImplicitIter(IEl) ! loop over all material points in element
                  ParticleIndex = GetParticleIndex_DirDeriv_ImplicitIter(IPart, IEl) ! Get the particle ID

                  if((NFORMULATION==1).or.(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid)) then ! Only if NumberOfLayers = 1 or SOLID Material Point

                      !ParticleVelocity = VelocityArray(ParticleIndex,:) !get particle velocity vector
                      ParticleMass = MassArray(ParticleIndex) !get particle's mass
                      PartShape = ShapeValuesArray_DirDeriv_ImplicitIter(ParticleIndex,:) !get particle shape fucntions
                        
                      !if (CalParams%ApplyContactAlgorithm) then
                      !    iEntity = EntityIDArray(ParticleIndex)
                      !else
                          iEntity = 1
                      !end if

                      ! nodal x-mass
                      LumMassX(:, iEntity) = LumMassX(:, iEntity) + ParticleMass * PartShape(:)

                      ! nodal momentum in directions from 1 to NDIM
                      !do I=1,NVECTOR
                      !    Mom(I, :, iEntity) = Mom(I, :, iEntity) + ParticleMass * PartShape(:) * ParticleVelocity(I)
                      !end do

                  end if ! Only if NumberOfLayers = 1 or SOLID Material Point

              end do ! loop over material points in element

              do I = 1, NVECTOR
                  LumpedMass( Ni(I, 1:ELEMENTNODES), 1:NoEn ) = LumMassX( 1:ELEMENTNODES, 1:NoEn )
                  !Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn ) = Mom( I, 1:ELEMENTNODES, 1:NoEn )
              end do
                
          end do ! loop over Eelements
          
        
        
        end subroutine MapMassMP2Nodes_DirDeriv_ImplicitIter
        
        
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        
        
        
        
        
        
        
        
        ! SUBROUTINE #2
        
        subroutine MapMomentumMP2NodesInitialization_ImplicitIter(Momentum)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: Momentum       !CC added nEntity
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: IDof
          !real(REAL_TYPE), dimension(NVECTOR) :: ParticleIncrementalVelocity
          !real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          !real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: NodAcc
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) :: Time, PartilceMass
          
            
          ! implicit variables - implicit bookmark
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_NPlus1

          ! implicit variables - implicit bookmark 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          !real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(in) :: VelocityArray
          !real(REAL_TYPE), dimension(Counters%NParticles), intent(in) :: MassArray
          !integer(INTEGER_TYPE), dimension(Counters%NEl, ELEMENTNODES), intent(in) :: ElementConnectivities
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          
          !!CC - changed to function call - needs to be set to zero here
          Momentum = 0.0 
          Time = CalParams%TimeIncrement
            
          do IAEl = 1, Counters%NAEl   ! Loop over all elements --> NAEl N and NPlus1
            IEl = ActiveElement(IAEl)
            do I = 1, NVECTOR
              IDof(I, 1:ELEMENTNODES) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
            end do
            
            do IPart = 1, NPartEle(IEl)   ! Loop over all particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 ParticleVelocity = VelocityArray(ParticleIndex,:)
                 PartilceMass = MassArray(ParticleIndex)
                 ParticleShape = ShapeValuesArray(ParticleIndex,:)

                  iEntity = 1


                do I = 1, NVECTOR ! nodal i-momentum
                  Momentum(IDof(I,1:ELEMENTNODES), iEntity) = Momentum(IDof(I,1:ELEMENTNODES), iEntity) + PartilceMass * ParticleShape * ParticleVelocity(I)
                end do     

                 
              end if ! NumbOfLayers = 1 or SOLID MatPoint
            end do !Loop over particles
          end do !elements    
          
        end subroutine MapMomentumMP2NodesInitialization_ImplicitIter  
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
         ! SUBROUTINE #3
        
        subroutine GetNodalVelocityFromNodalMomentum_ImplicitIter(Momentum, LumpedMassDry, TotalVelocitySoil)!, DoSystem)
        !**********************************************************************
        !
        ! Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Momentum !input
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry !input
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: TotalVelocitySoil !input
          
          !real(REAL_TYPE), dimension(Counters%N), intent(in) :: PBoundary
          
          logical DoSystem !inout 
                   
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF,J
          real(REAL_TYPE), dimension(Counters%N) :: MomentumSystem
          

          do IDOF = 1, Counters%N ! loop over all degrees of freedom
            do J = 1, Counters%nEntity ! loop over all entities
                
                if (LumpedMassDry(IDOF,J)/=0) then
                  TotalVelocitySoil(IDOF,J) = ( Momentum(IDOF,J) / LumpedMassDry(IDOF,J) ) * PBoundary(IDOF)
                else
                  TotalVelocitySoil(IDOF,J) = 0.0
                end if
                              
            end do ! loop over all entities
          end do ! loop over all degrees of freedom
          
        end subroutine GetNodalVelocityFromNodalMomentum_ImplicitIter
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #4 - BOOKMARK
        subroutine GetNodalExternalForces(ExtLoad, GravityLoad)
        !**********************************************************************
        !
        !    Function:  To extrapolate loads from material points to nodes using 
        !               the shape function values evaluated at the material points local position.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: IEntity
          integer(INTEGER_TYPE) :: IDim, ILoadSystem

          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: ExtLoad
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: GravityLoad
          
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          !integer(INTEGER_TYPE), dimension(counters%NEl), intent(inout) :: ActiveElement
          
          IDim = NVECTOR
           
          call MPMDYNLoad(ExtLoad, GravityLoad)

           if ((Counters%NLoadedElementSidesSolidNodes+Counters%NLoadedElementSidesSolidNodesB) > 0) then
             do ILoadSystem =1, Counters%NSoilLoadSystems
                 ExtLoad = ExtLoad + ExtLoadTotal(:,:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
             end do
           end if
           

        end subroutine GetNodalExternalForces
     
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
        ! SUBROUTINE #4 - BOOKMARK
        subroutine GetNodalExternalForces_ImplicitIter(ExtLoad, GravityLoad)
        !**********************************************************************
        !
        !    Function:  To extrapolate loads from material points to nodes using 
        !               the shape function values evaluated at the material points local position.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: IEntity
          integer(INTEGER_TYPE) :: IDim, ILoadSystem

          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: ExtLoad
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: GravityLoad
          
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          !integer(INTEGER_TYPE), dimension(counters%NEl), intent(inout) :: ActiveElement
          
          IDim = NVECTOR
           
          call MPMDYNLoad_ImplicitIter(ExtLoad, GravityLoad)

           if ((Counters%NLoadedElementSidesSolidNodes+Counters%NLoadedElementSidesSolidNodesB) > 0) then
             do ILoadSystem =1, Counters%NSoilLoadSystems
                 ExtLoad = ExtLoad + ExtLoadTotal(:,:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
             end do
           end if
           

        end subroutine GetNodalExternalForces_ImplicitIter
        
        !--------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------
        
        ! SUBROUTINE #4 - BOOKMARK
        subroutine GetNodalExternalForces_DirDeriv_ImplicitIter(ExtLoad, GravityLoad)
        !**********************************************************************
        !
        !    Function:  To extrapolate loads from material points to nodes using 
        !               the shape function values evaluated at the material points local position.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: IEntity
          integer(INTEGER_TYPE) :: IDim, ILoadSystem

          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: ExtLoad
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: GravityLoad
          
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          !integer(INTEGER_TYPE), dimension(counters%NEl), intent(inout) :: ActiveElement
          
          IDim = NVECTOR
           
          call MPMDYNLoad_DirDeriv_ImplicitIter(ExtLoad, GravityLoad)

           if ((Counters%NLoadedElementSidesSolidNodes+Counters%NLoadedElementSidesSolidNodesB) > 0) then
             do ILoadSystem =1, Counters%NSoilLoadSystems
                 ExtLoad = ExtLoad + ExtLoadTotal(:,:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
             end do
           end if
           

        end subroutine GetNodalExternalForces_DirDeriv_ImplicitIter
        
        !--------------------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------
        
        
        ! SUBROUTINE #6 - BOOKMARK
        
        subroutine GetFirstGuessOfVelocity_ImplicitIter(TotalVelocitySoil_ImplicitIter, TotalVelocitySoil)
        ! purpose: obtain a first guess for velocity
        !
        implicit none 
        
        real(REAL_TYPE), dimension(Counters%N), intent(inout) :: TotalVelocitySoil
        real(REAL_TYPE), dimension(Counters%N), intent(inout) :: TotalVelocitySoil_ImplicitIter
            
        ! this first precition can be modified depending on what we are trying to do
        TotalVelocitySoil_ImplicitIter = 1*TotalVelocitySoil ! hardcoded increase by 20 percent
        
        end subroutine GetFirstGuessOfVelocity_ImplicitIter
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #7 - BOOKMARK
         subroutine GetMPStrainFromNodeVelocity_ImplicitIter(TotalVelocitySoil)
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

         ! implicit bookmark
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum !, intent(in)... we do not need this intent(in)
         
         integer(INTEGER_TYPE) :: IEntity, StructureMaterialEntity, LoadedSides, IDOF
         
         !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: IntLoad
         
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalDisplacementSoil !local
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IncrementalDisplacementSoil !local
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccumulatedIncDisplacementSoil !local 
         
         ! implicit scheme - shape values array inout 
         !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
         !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
         integer(INTEGER_TYPE), dimension(Counters%NEl) :: ActiveElement
         
         ! local stresses 
         !real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR) :: SigmaEffArray_Local
              
         EntityUsed = 1   !!CC - each node has nEntity displacements, choose entity used for nodal update
          
          ! SUBROUTINE #1
          call GetNodalIncrementalDisplacement(IncrementalDisplacementSoil_ImplicitIter, TotalVelocitySoil)!_ImplicitIter)
          
          !AccumulatedIncDisplacementSoil = AccumulatedIncDisplacementSoil + IncrementalDisplacementSoil !--> I do not care about this accumulation!
          
          ! SUBROUTINE #2: calculate IncrementalDisplacementSoil
          !call UpdateNodalTotalDisplacementImplicit(EntityUsed, IncrementalDisplacementSoil, TotalDisplacementSoil)
    
          ! SUBROUTINE #3: Update particle strain data
          call UpdateParticleStrains_ImplicitIter(IncrementalDisplacementSoil_ImplicitIter)!, DShapeValuesArray, ActiveElement, NPartEle)!, GetParticleIndex_Local) ! update particle strains for solid and mixture
          
         end subroutine GetMPStrainFromNodeVelocity_ImplicitIter
         
         !-----------------------------------------------------------------------------------------------------
         
         
         !-----------------------------------------------------------------------------------------------------
         
         !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #7 - BOOKMARK
         subroutine GetMPStrainFromNodeVelocity_DirDeriv_ImplicitIter(TotalVelocitySoil)
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

         ! implicit bookmark
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum !, intent(in)... we do not need this intent(in)
         
         integer(INTEGER_TYPE) :: IEntity, StructureMaterialEntity, LoadedSides, IDOF
         
         !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: IntLoad
         
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalDisplacementSoil !local
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IncrementalDisplacementSoil !local
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccumulatedIncDisplacementSoil !local 
         
         ! implicit scheme - shape values array inout 
         !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
         !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
         integer(INTEGER_TYPE), dimension(Counters%NEl) :: ActiveElement
         
         ! local stresses 
         !real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR) :: SigmaEffArray_Local
              
         EntityUsed = 1   !!CC - each node has nEntity displacements, choose entity used for nodal update
          
          ! SUBROUTINE #1
          call GetNodalIncrementalDisplacement(IncrementalDisplacementSoil_DirDeriv_ImplicitIter, TotalVelocitySoil_DirDeriv_ImplicitIter)
          
          !AccumulatedIncDisplacementSoil = AccumulatedIncDisplacementSoil + IncrementalDisplacementSoil !--> I do not care about this accumulation!
          
          ! SUBROUTINE #2: calculate IncrementalDisplacementSoil
          !call UpdateNodalTotalDisplacementImplicit(EntityUsed, IncrementalDisplacementSoil, TotalDisplacementSoil)
    
          ! SUBROUTINE #3: Update particle strain data
          call UpdateParticleStrains_DirDeriv_ImplicitIter(IncrementalDisplacementSoil_DirDeriv_ImplicitIter)!, DShapeValuesArray, ActiveElement, NPartEle)!, GetParticleIndex_Local) ! update particle strains for solid and mixture
          !call UpdateParticleStrains(IncrementalDisplacementSoil_DirDeriv_ImplicitIter)!, DShapeValuesArray, ActiveElement, NPartEle)!, GetParticleIndex_Local) ! update particle strains for solid and mixture

          
         end subroutine GetMPStrainFromNodeVelocity_DirDeriv_ImplicitIter
         
         !-----------------------------------------------------------------------------------------------------
         !-----------------------------------------------------------------------------------------------------
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
         
         subroutine UpdateNodalTotalDisplacementImplicit(EntityUsed, IncrementalDisplacementSoil, TotalDisplacementSoil)
        !**********************************************************************
        !
        !    Function:  To add the nodal incremental displacements to the total displacemnts
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF 
          integer(INTEGER_TYPE), intent(in) :: EntityUsed
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: IncrementalDisplacementSoil
          ! why is TotalDisplacementSoil not sperated for different entities
          real(REAL_TYPE), dimension(Counters%N), intent(out) :: TotalDisplacementSoil
          
          do IDOF = 1, Counters%N    !Loop over all degrees of freedom
            TotalDisplacementSoil(IDOF)= TotalDisplacementSoil(IDOF)+ IncrementalDisplacementSoil(IDOF,EntityUsed)
          end do
          
          if(.not.(NFORMULATION==1)) then
          do IDOF = 1, Counters%N    !Loop over all degrees of freedom
            TotalDisplacementWater(IDOF)= TotalDisplacementWater(IDOF)+ IncrementalDisplacementWater(IDOF,EntityUsed)
          end do
          end if
          
         end subroutine UpdateNodalTotalDisplacementImplicit  
         
         !----------------------------------------------------------------------------------------------
         !----------------------------------------------------------------------------------------------
         
         
         subroutine UpdateParticleStrainsImplicit(IncrementalDisplacementSoil)!, DShapeValuesArray, &
            !ActiveElement, NPartEle)!, NAEl)!, GetParticleIndex)!DShapeValues) ! all these variables need to be local 
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
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
          
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: IncrementalDisplacementSoil
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: TotalDisplacementSoil
          
          ! we input local variables into this subroutine
          !integer(INTEGER_TYPE), dimension(Counters%NAEl), intent(in) :: ActiveElement !list of IDs of active elements
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle !number of particles in an element
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
                ParticleIndex = GetParticleIndex(IParticle, IElement) ! get global particle ID
                
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

         end subroutine UpdateParticleStrainsImplicit
         
         
         !----------------------------------------------------------------------------------------------
         !----------------------------------------------------------------------------------------------
         
         
         
         !----------------------------------------------------------------------------------------------
         
         
         subroutine UpdateParticleStrains_ImplicitIter(IncrementalDisplacementSoil)!, DShapeValuesArray, &
            !ActiveElement, NPartEle)!, NAEl)!, GetParticleIndex)!DShapeValues) ! all these variables need to be local 
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
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
          
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: IncrementalDisplacementSoil
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: TotalDisplacementSoil
          
          ! we input local variables into this subroutine
          !integer(INTEGER_TYPE), dimension(Counters%NAEl), intent(in) :: ActiveElement !list of IDs of active elements
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle !number of particles in an element
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
                ParticleIndex = GetParticleIndex(IParticle, IElement) ! get global particle ID
                
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
                  call SetEpsStep_ImplicitIter(Particles(ParticleIndex), Eps)
                  ! we do not want to increase the total particle. We calculate the incremental strain and then use it to find the stress.  
                  ! skip increasing the total strain of the particle and leave it for later after we achieve the solution
                  !if (.not.CalParams%ImplicitIntegration%IsImplicitScheme) then !CalParams%ApplyImplicitQuasiStatic
                    !call IncreaseEps(Particles(ParticleIndex), Eps)
                  !end if
                
                end if ! Material Point can be MIXTURE or SOLID

              end do ! particle loop

          end do ! active element loop

         end subroutine UpdateParticleStrains_ImplicitIter
         
         
         !----------------------------------------------------------------------------------------------
         
         
         
         
         
         
         
         
         
         
         
         !----------------------------------------------------------------------------------------------
         
         
         subroutine IncreaseParticleStrains_ImplicitIter!(IncrementalDisplacementSoil)!, DShapeValuesArray, &
            !ActiveElement, NPartEle)!, NAEl)!, GetParticleIndex)!DShapeValues) ! all these variables need to be local 
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
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
          
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: IncrementalDisplacementSoil
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: TotalDisplacementSoil
          
          ! we input local variables into this subroutine
          !integer(INTEGER_TYPE), dimension(Counters%NAEl), intent(in) :: ActiveElement !list of IDs of active elements
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle !number of particles in an element
          !integer(INTEGER_TYPE) :: NAEl ! number of active elements 
          !integer(INTEGER_TYPE), dimension(:,:), intent(in) :: GetParticleIndex ! we need this as an input
          
          ! set the nodal displacement increments to the first entity
          ! only need to change if the contact model is used and the next particle belongs to a different entity
          iEntityDefault = 1 ! set default to first entity
          !do I = 1,Counters%N
          !  DUTotEnt(I) = IncrementalDisplacementSoil(I,iEntityDefault)
          !end do
            
          do IActiveElement = 1, Counters%NAEl_ImplicitIter!Counters%NAEl ! loop over all active elements
                
              IElement = ActiveElement_ImplicitIter(IActiveElement) ! get element number of the active element
              NElemPart = NPartEle_ImplicitIter(IElement) ! get number of particles in the active element
              
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if  
              
              do IParticle = 1, NElemPart ! loop over all particles of the active element
              
                  ! there needs to be GetParticleIndex "local for iteration" and "actual"
                ParticleIndex = GetParticleIndex_ImplicitIter(IParticle, IElement) ! get global particle ID
                
                if( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then 

                  ! determine B matrix of deformed element at the particle local position
                  !call BMatrix(Particles(ParticleIndex)%LocPos, &
                  !             ELEMENTNODES, Counters%NEl, &
                  !             Counters%NodTot, NDIM, &
                  !             IElement, ElementConnectivities, &
                  !             NodalCoordinatesUpd,  &
                  !             BMatrixDeformed, DetJac)!, & 
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
                  !call Get_Strain(IElement, IParticle, &
                  !              ElementConnectivities, & 
                  !              BMatrixDeformed, DUTotEnt, ReducedDof, & !--> DUTotEnt is where IncrementalDisplacementSoil
                  !              Eps) ! Eps = strain increment at material point
                    
                  ! Update particle strains
                  !call SetEpsStep_ImplicitIter(Particles(ParticleIndex), Eps)
                  ! we do not want to increase the total particle. We calculate the incremental strain and then use it to find the stress.  
                  ! skip increasing the total strain of the particle and leave it for later after we achieve the solution
                  !if (.not.CalParams%ImplicitIntegration%IsImplicitScheme) then !CalParams%ApplyImplicitQuasiStatic
                    call IncreaseEps_ImplicitIter(Particles(ParticleIndex), Eps)
                  !end if
                
                end if ! Material Point can be MIXTURE or SOLID

              end do ! particle loop

          end do ! active element loop

         end subroutine IncreaseParticleStrains_ImplicitIter
         
         
         !----------------------------------------------------------------------------------------------
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         subroutine UpdateParticleStrains_DirDeriv_ImplicitIter(IncrementalDisplacementSoil)!, DShapeValuesArray, &
            !ActiveElement, NPartEle)!, NAEl)!, GetParticleIndex)!DShapeValues) ! all these variables need to be local 
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
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
          
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: IncrementalDisplacementSoil
          !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: TotalDisplacementSoil
          
          ! we input local variables into this subroutine
          !integer(INTEGER_TYPE), dimension(Counters%NAEl), intent(in) :: ActiveElement !list of IDs of active elements
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle !number of particles in an element
          !integer(INTEGER_TYPE) :: NAEl ! number of active elements 
          !integer(INTEGER_TYPE), dimension(:,:), intent(in) :: GetParticleIndex ! we need this as an input
          
          ! set the nodal displacement increments to the first entity
          ! only need to change if the contact model is used and the next particle belongs to a different entity
          iEntityDefault = 1 ! set default to first entity
          do I = 1,Counters%N
            DUTotEnt(I) = IncrementalDisplacementSoil(I,iEntityDefault)
          end do
            
          do IActiveElement = 1, Counters%NAEl_ImplicitIter!Counters%NAEl_ImplicitIter!Counters%NAEl ! loop over all active elements
                
              IElement = ActiveElement_ImplicitIter(IActiveElement) !ActiveElement_ImplicitIter(IActiveElement) ! get element number of the active element
              NElemPart = NPartEle_ImplicitIter(IElement)!NPartEle_ImplicitIter(IElement) ! get number of particles in the active element
              
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if  
              
              do IParticle = 1, NElemPart ! loop over all particles of the active element
              
                  ! there needs to be GetParticleIndex "local for iteration" and "actual"
                ParticleIndex = GetParticleIndex_ImplicitIter(IParticle, IElement)!GetParticleIndex_ImplicitIter(IParticle, IElement) ! get global particle ID
                
                if( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then 
                
                  ! determine B matrix of deformed element at the particle local position
                  call BMatrix(Particles(ParticleIndex)%LocPos_ImplicitIter, &!_ImplicitIter, &
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
                  !call SetEpsStep_DirDeriv_ImplicitIter(Particles(ParticleIndex), Eps)
                  call SetEpsStep_DirDeriv_ImplicitIter(Particles(ParticleIndex), Eps)
                  ! we do not want to increase the total particle. We calculate the incremental strain and then use it to find the stress.  
                  ! skip increasing the total strain of the particle and leave it for later after we achieve the solution
                  !if (.not.CalParams%ImplicitIntegration%IsImplicitScheme) then !CalParams%ApplyImplicitQuasiStatic
                  !  call IncreaseEps(Particles(ParticleIndex), Eps)
                  !end if
                
                end if ! Material Point can be MIXTURE or SOLID

              end do ! particle loop

          end do ! active element loop

        end subroutine UpdateParticleStrains_DirDeriv_ImplicitIter
         
         !----------------------------------------------------------------------------------------------
         
         
         
         
         !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #8 - BOOKMARK
        subroutine MPMDYNGetSig_ImplicitIter()!SigmaEffArray)!(IntLoad_NPlus1) !MPMDYNGetSig ! NAEl, 
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
          
          ! implicit scheme variable - implicit bookmark
          !real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(inout) :: SigmaEffArray
          !real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(in) :: SigmaEff0Array
          
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement         
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          !
          !real(REAL_TYPE), dimension(Counters%NParticles, NSTATEVAR), intent(inout) :: ESMStatevArray
          
          !integer(INTEGER_TYPE), dimension(NElements), intent(in) :: NumberOfIntegrationPoints
          
          !integer(INTEGER_TYPE), intent(in) :: NAEl
          !integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(in) :: MaterialIDArray
          
          IEntity = 1 
          
          !! implicit scheme variable - implicit bookmark 
          !real(REAL_TYPE), dimension(Counters%N,Counters%NEntity), intent(out) :: IntLoad_NPlus1

          ! Initialise global variables
          CalParams%IntegrationPointData%NPlasticPoints = 0
          CalParams%IntegrationPointData%NNegativePlasticPoints = 0
          CalParams%IntegrationPointData%NApexPoints = 0
          CalParams%IntegrationPointData%NTensionCutOffPoints = 0
          CalParams%ConvergenceCheck%NInaccuratePlasticPoints = 0
          CalParams%ConvergenceCheck%SumLocalError = 0.0
          CalParams%ConvergenceCheck%SumIntegrationPointWeights = 0.0
            
          !!******START calculation of effective stresses
          do IAEl = 1, Counters%NAEl!Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)
            
            ! Get the number of integration per element (Gaussian integration points or number of material points)
            NElemPart = NumberOfIntegrationPoints(IEl) 
     
            ! Calculate stresses in integration/material points
            do Int = 1, NElemPart ! Loop over all integration/material points of the element
         
              IntGlo = GetParticleIndex(Int, IEl)   ! Determine global ID of integration point 
              IDMaterial = MaterialIDArray(IntGlo)  ! Material number stored in $$MATERIAL_INDEX in the GOM-file
                          
              if ( (MaterialPointTypeArray(IntGlo)==MaterialPointTypeMixture).or.(MaterialPointTypeArray(IntGlo)==MaterialPointTypeSolid) ) then

                if((NFORMULATION==1).or.(Particles(IntGlo)%PhaseStatus==PhaseStatusSOLID)) then 

                 DoSkipStressComputation = EntityIDArray(IntGlo) == HARD_ENTITY
                 do I = 1, NVECTOR
                   DoSkipStressComputation = DoSkipStressComputation .and. CalParams%ApplyPrescribedVelocity(I)
                 end do
                    
                 call StressSolid_ImplicitIter(IntGlo, IEl, BMatrix, IEntity, SigmaEff0Array, SigmaEffArray_ImplicitIter,  &
                     ESMStatevArray, ESMStatevArray_ImplicitIter)   ! calculate stresses for mixture or solid material points
                     ! ESMStatevArray remains unchanged here
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

        end subroutine MPMDYNGetSig_ImplicitIter
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        subroutine MPMDYNGetSig_DirDeriv_ImplicitIter()!SigmaEffArray)!(IntLoad_NPlus1) !MPMDYNGetSig ! NAEl, 
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
          
          ! implicit scheme variable - implicit bookmark
          !real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(inout) :: SigmaEffArray
          !real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(in) :: SigmaEff0Array
          
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement         
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          !
          !real(REAL_TYPE), dimension(Counters%NParticles, NSTATEVAR), intent(inout) :: ESMStatevArray
          
          !integer(INTEGER_TYPE), dimension(NElements), intent(in) :: NumberOfIntegrationPoints
          
          !integer(INTEGER_TYPE), intent(in) :: NAEl
          !integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(in) :: MaterialIDArray
          
          IEntity = 1 
          
          !! implicit scheme variable - implicit bookmark 
          !real(REAL_TYPE), dimension(Counters%N,Counters%NEntity), intent(out) :: IntLoad_NPlus1

          ! Initialise global variables
          CalParams%IntegrationPointData%NPlasticPoints = 0
          CalParams%IntegrationPointData%NNegativePlasticPoints = 0
          CalParams%IntegrationPointData%NApexPoints = 0
          CalParams%IntegrationPointData%NTensionCutOffPoints = 0
          CalParams%ConvergenceCheck%NInaccuratePlasticPoints = 0
          CalParams%ConvergenceCheck%SumLocalError = 0.0
          CalParams%ConvergenceCheck%SumIntegrationPointWeights = 0.0
            
          !!******START calculation of effective stresses
          do IAEl = 1, Counters%NAEl_ImplicitIter!Counters%NAEl_ImplicitIter!Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement_ImplicitIter(IAEl) !ActiveElement_ImplicitIter(IAEl)
            
            ! Get the number of integration per element (Gaussian integration points or number of material points)
            NElemPart = NumberOfIntegrationPoints(IEl) 
     
            ! Calculate stresses in integration/material points
            do Int = 1, NElemPart ! Loop over all integration/material points of the element
         
              IntGlo = GetParticleIndex_ImplicitIter(Int, IEl)!GetParticleIndex_ImplicitIter(Int, IEl)   ! Determine global ID of integration point 
              IDMaterial = MaterialIDArray(IntGlo)  ! Material number stored in $$MATERIAL_INDEX in the GOM-file
                          
              if ( (MaterialPointTypeArray(IntGlo)==MaterialPointTypeMixture).or.(MaterialPointTypeArray(IntGlo)==MaterialPointTypeSolid) ) then

                if((NFORMULATION==1).or.(Particles(IntGlo)%PhaseStatus==PhaseStatusSOLID)) then 

                 DoSkipStressComputation = EntityIDArray(IntGlo) == HARD_ENTITY
                 do I = 1, NVECTOR
                   DoSkipStressComputation = DoSkipStressComputation .and. CalParams%ApplyPrescribedVelocity(I)
                 end do
                    
                 !call StressSolid_DirDeriv_ImplicitIter(IntGlo, IEl, BMatrix, IEntity, SigmaEff0Array, SigmaEffArray_DirDeriv_ImplicitIter,  &
                 !    ESMStatevArray_ImplicitIter, ESMStatevArray_DirDeriv_ImplicitIter)   ! calculate stresses for mixture or solid material points
                 call StressSolid_DirDeriv_ImplicitIter(IntGlo, IEl, BMatrix, IEntity, SigmaEffArray_ImplicitIter, SigmaEffArray_DirDeriv_ImplicitIter,  &
                     ESMStatevArray_ImplicitIter, ESMStatevArray_DirDeriv_ImplicitIter)   ! calculate stresses for mixture or solid material points
                     ! ESMStatevArray remains unchanged here
                end if 
                  
               end if
               
            end do ! Loop over all integration/material points of the element
                        
            ! Assign EUnloading to all particles (material points) of each element
            if (NElemPart==1) then ! only for fully filled elements               
             NumParticles = NPartEle_ImplicitIter(IEl)!NPartEle_ImplicitIter(IEl)
             FirstParticleIndex = GetParticleIndex_ImplicitIter(1, IEl) !GetParticleIndex_ImplicitIter(1, IEl)
             FirstParticleEUnloading = Particles(FirstParticleIndex)%ESM_UnloadingStiffness
             do Int = 2, NumParticles
                IntGlo = GetParticleIndex(Int, IEl)!GetParticleIndex_ImplicitIter(Int, IEl)
                Particles(IntGlo)%ESM_UnloadingStiffness = FirstParticleEUnloading
             enddo
            endif
            
          end do ! Loop over all active elements for computation of stresses

        end subroutine MPMDYNGetSig_DirDeriv_ImplicitIter
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
         subroutine StressSolid_ImplicitIter(IDpt, IDel, BMatrix, IEntityID, SigmaEff0Array, SigmaEffArray,  ESMstatevArray, ESMstatevArray_ImplicitIter) !MaterialIDArray, 
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


    ! intent in 
    real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(inout) :: SigmaEffArray
    real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(in) :: SigmaEff0Array
    !integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(in) :: MaterialIDArray
    real(REAL_TYPE), dimension(Counters%NParticles, NSTATEVAR), intent(in) :: ESMstatevArray
    real(REAL_TYPE), dimension(Counters%NParticles, NSTATEVAR), intent(out) :: ESMstatevArray_ImplicitIter
    
    
    pointer (p, ESM)             
    
        
    StateVar = 0

    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
    NameModel = MatParams(IDset)%MaterialModel ! name of constitutive model as specified in GOM-file
    ntens = MatParams(IDset)%UMATDimension  ! 2D or 3D formulation of the External soil model   
          
    ! get strain increments in integration/material point
    TempStrainIncr = GetEpsStep_ImplicitIter(Particles(IDpt)) ! incremental strain vector assigned to point
    
    !if (CalParams%ApplyImplicitQuasiStatic) then
    !    if (CalParams%ImplicitIntegration%Iteration > 1) then
    !        do I = 1, NTENSOR
    !            TempStrainIncrPrevious(I) = GetEpsStepPreviousI(Particles(IDpt), I)
    !        end do
    !        
    !        TempStrainIncr = TempStrainIncr - TempStrainIncrPrevious
    !        
    !    end if
    !end if
        
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
    !if (IsUndrEffectiveStress) then
    !    if (Particles(IDpt)%Porosity > 0.0) then
    !    Bulk = Particles(IDpt)%BulkWater / Particles(IDpt)%Porosity ! kN/m2
    !    DSigWP = Bulk * DEpsVol
    !    else
    !    DSigWP = 0.0
    !    end if
    !    call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP)
    !end if ! effective stress analysis
          
    ! get stresses in integration/material point      
    do I = 1, NTENSOR
    Sig0(I) = SigmaEff0Array(IDpt, I) ! get initial stress of current step assigned to point 
    end do
    Stress=0.0
    do I=1, NTENSOR
        Stress(I) = Stress(I) + Sig0(I) !To use 3D UMAT also for 2D formulations
    enddo 
          
    ! initialise state variables (only for very first time and load step)
    !if ((CalParams%IStep == 1).and.(CalParams%TimeStep == 1)) then
    !StateVar = MatParams(IDset)%ESM_Statvar_in
    !else 
    StateVar = ESMstatevArray(IDpt,:) ! temporary
    !end if 
          
    
    
    !call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP) !Note that the subroutine checks Cavitation Threshold & Gas Pressure
          
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
    !props = ESMpropsArray(IDpt,:)
         
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
    call ESM(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
    ! save unloading stiffness in Particles array  
    Particles(IDpt)%ESM_UnloadingStiffness = Eunloading
                 
    if (IsUndrEffectiveStress) then
        Particles(IDpt)%BulkWater = AdditionalVar(12)
    end if
    
    call SetIPL(IDpt, IDel, int(PlasticMultiplier))


    
    ! to use objective stress definition
    !if (CalParams%ApplyObjectiveStress) then ! Consider large deformation terms
    !call Hill(IdEl, ELEMENTNODES, IncrementalDisplacementSoil(1:Counters%N, IEntityID),  &
    !                 ReducedDof, ElementConnectivities, BMatrix, Sig0(1:NTENSOR), Stress(1:NTENSOR), DEpsVol)
    !end if ! objective stress            
            
    ! write new stresses to global array
    do I=1, NTENSOR
        StressIncr(I) = Stress(I) - Sig0(I)
    enddo             
                               
    ! save updated state variables and in Particles array
    ESMstatevArray_ImplicitIter(IDpt,:) = StateVar ! temporary
          
    call CalculatePrincipalStresses(IDpt, Stress(1:NTENSOR), StressPrinc)
    call AssignStressStrainToGlobalArrayESM(IDpt, NTENSOR, StressIncr, StressPrinc, StrainIncr, SigmaEffArray_ImplicitIter )
    
    ! write plasticity state to global array
    !  call SetIPL(IDpt, IDel, int(StateVar(50)))
    !if (CalParams%ApplyBulkViscosityDamping) then
    !RateVolStrain(IDEl) = DEpsVol / CalParams%TimeIncrement
    !call CalculateViscousDamping_interface(IDpt, IDEl)
    !end if  
    
    end subroutine StressSolid_ImplicitIter
         
         !--------------------------------------------------------------------------------------------------------
         !--------------------------------------------------------------------------------------------------------
         !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
         subroutine StressSolid_DirDeriv_ImplicitIter(IDpt, IDel, BMatrix, IEntityID, SigmaEff0Array, SigmaEffArray,  ESMstatevArray, ESMstatevArray_ImplicitIter) !MaterialIDArray, 
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


    ! intent in 
    real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(inout) :: SigmaEffArray
    real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(in) :: SigmaEff0Array
    !integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(in) :: MaterialIDArray
    real(REAL_TYPE), dimension(Counters%NParticles, NSTATEVAR), intent(in) :: ESMstatevArray
    real(REAL_TYPE), dimension(Counters%NParticles, NSTATEVAR), intent(out) :: ESMstatevArray_ImplicitIter
    
    
    pointer (p, ESM)             
    
        
    StateVar = 0

    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
    NameModel = MatParams(IDset)%MaterialModel ! name of constitutive model as specified in GOM-file
    ntens = MatParams(IDset)%UMATDimension  ! 2D or 3D formulation of the External soil model   
          
    ! get strain increments in integration/material point
    TempStrainIncr = GetEpsStep_DirDeriv_ImplicitIter(Particles(IDpt)) ! incremental strain vector assigned to point
    
    !if (CalParams%ApplyImplicitQuasiStatic) then
    !    if (CalParams%ImplicitIntegration%Iteration > 1) then
    !        do I = 1, NTENSOR
    !            TempStrainIncrPrevious(I) = GetEpsStepPreviousI(Particles(IDpt), I)
    !        end do
    !        
    !        TempStrainIncr = TempStrainIncr - TempStrainIncrPrevious
    !        
    !    end if
    !end if
        
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
    !if (IsUndrEffectiveStress) then
    !    if (Particles(IDpt)%Porosity > 0.0) then
    !    Bulk = Particles(IDpt)%BulkWater / Particles(IDpt)%Porosity ! kN/m2
    !    DSigWP = Bulk * DEpsVol
    !    else
    !    DSigWP = 0.0
    !    end if
    !    call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP)
    !end if ! effective stress analysis
          
    ! get stresses in integration/material point      
    do I = 1, NTENSOR
    Sig0(I) = SigmaEff0Array(IDpt, I) ! get initial stress of current step assigned to point 
    end do
    Stress=0.0
    do I=1, NTENSOR
        Stress(I) = Stress(I) + Sig0(I) !To use 3D UMAT also for 2D formulations
    enddo 
          
    ! initialise state variables (only for very first time and load step)
    !if ((CalParams%IStep == 1).and.(CalParams%TimeStep == 1)) then
    !StateVar = MatParams(IDset)%ESM_Statvar_in
    !else 
    StateVar = ESMstatevArray_ImplicitIter(IDpt,:) ! temporary
    !end if 
          
    
    
    !call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP) !Note that the subroutine checks Cavitation Threshold & Gas Pressure
          
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
    !props = ESMpropsArray(IDpt,:)
         
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
    call ESM(IDpt, IDel, IDset, Stress, Eunloading, PlasticMultiplier, StrainIncr, NSTATEVAR, StateVar, nAddVar, AdditionalVar,cmname, NPROPERTIES, props, CalParams%NumberOfPhases, ntens)
    ! save unloading stiffness in Particles array  
    Particles(IDpt)%ESM_UnloadingStiffness = Eunloading
                 
    if (IsUndrEffectiveStress) then
        Particles(IDpt)%BulkWater = AdditionalVar(12)
    end if
    
    call SetIPL(IDpt, IDel, int(PlasticMultiplier))


    
    ! to use objective stress definition
    !if (CalParams%ApplyObjectiveStress) then ! Consider large deformation terms
    !call Hill(IdEl, ELEMENTNODES, IncrementalDisplacementSoil(1:Counters%N, IEntityID),  &
    !                 ReducedDof, ElementConnectivities, BMatrix, Sig0(1:NTENSOR), Stress(1:NTENSOR), DEpsVol)
    !end if ! objective stress            
            
    ! write new stresses to global array
    do I=1, NTENSOR
        StressIncr(I) = Stress(I) - Sig0(I)
    enddo             
                               
    ! save updated state variables and in Particles array
    ESMstatevArray_ImplicitIter(IDpt,:) = StateVar ! temporary
          
    call CalculatePrincipalStresses(IDpt, Stress(1:NTENSOR), StressPrinc)
    call AssignStressStrainToGlobalArrayESM(IDpt, NTENSOR, StressIncr, StressPrinc, StrainIncr, SigmaEffArray )
    
    ! write plasticity state to global array
    !  call SetIPL(IDpt, IDel, int(StateVar(50)))
    !if (CalParams%ApplyBulkViscosityDamping) then
    !RateVolStrain(IDEl) = DEpsVol / CalParams%TimeIncrement
    !call CalculateViscousDamping_interface(IDpt, IDEl)
    !end if  
    
    end subroutine StressSolid_DirDeriv_ImplicitIter
         
         !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
        subroutine DynUpdateParticleWeights_ImplicitIter()!( ShapeValuesArray, ActiveElement, NPartEle )
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
        
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! 1-layer formulation or 2-layer form with 1 Phase
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if((NFORMULATION==1).or. &
               ((NFORMULATION==2).and.(CalParams%NumberOfPhases==1))) then
                                
          do IAEl = 1, Counters%NAEl!_ImplicitIter ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)!_ImplicitIter(IAEl)   
            
            do Int = 1, NPartEle(IEl)!_ImplicitIter(IEl) ! Loop over all material points of the element
                
              ! Determine global ID of Material Point 
              IntGlo = GetParticleIndex(Int, IEl)!_ImplicitIter(Int, IEl)
              
              ! we need to keep track of this...
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight

              
              ! Determine volumetric strain of Material Point 
              VolumetricStrain = 0.0
              do I = 1, NVECTOR ! for 2D and 3D all three components have to be taken into account as eps_tt in axisymmetric \=0
                 VolumetricStrain = VolumetricStrain + GetEpsStepI_ImplicitIter(Particles(IntGlo),I)
              end do
          
             
              ! for SOLID or MIXTURE Material Point               
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight_Local * (1.0 + VolumetricStrain)
              ! we calculate the dir deriv based on the previous result so the %Integration Weight does not need updating
              Particles(IntGlo)%IntegrationWeight_ImplicitIter = Particles(IntGlo)%IntegrationWeight * (1.0 + VolumetricStrain)
            
             end do !loop over material points
            end do ! loop over elements
            
          end if ! 1-layer formulation or 2-layer form with 1 Phase
          
         
        end subroutine DynUpdateParticleWeights_ImplicitIter
        
        !----------------------------------------------------------------------------------------------------
        
        
        
        subroutine DynUpdateParticleWeights_DirDeriv_ImplicitIter()!( ShapeValuesArray, ActiveElement, NPartEle )
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
        
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! 1-layer formulation or 2-layer form with 1 Phase
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if((NFORMULATION==1).or. &
               ((NFORMULATION==2).and.(CalParams%NumberOfPhases==1))) then
                                
          do IAEl = 1, Counters%NAEl_ImplicitIter!_ImplicitIter!Counters%NAEl_ImplicitIter!_ImplicitIter ! Loop over all active elements for computation of stresses
            IEl = ActiveElement_ImplicitIter(IAEl)!ActiveElement_ImplicitIter(IAEl)!_ImplicitIter(IAEl)   
            
            do Int = 1, NPartEle_ImplicitIter(IEl)!NPartEle_ImplicitIter(IEl)!_ImplicitIter(IEl) ! Loop over all material points of the element
                
              ! Determine global ID of Material Point 
              IntGlo = GetParticleIndex_ImplicitIter(Int, IEl) !GetParticleIndex_ImplicitIter(Int, IEl)!_ImplicitIter(Int, IEl)
              
              ! we need to keep track of this...
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight

              
              ! Determine volumetric strain of Material Point 
              VolumetricStrain = 0.0
              do I = 1, NVECTOR ! for 2D and 3D all three components have to be taken into account as eps_tt in axisymmetric \=0
                 VolumetricStrain = VolumetricStrain + GetEpsStepI_DirDeriv_ImplicitIter(Particles(IntGlo),I)
              end do
          
             
              ! for SOLID or MIXTURE Material Point               
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight_Local * (1.0 + VolumetricStrain)
              ! we calculate the dir deriv based on the previous result so the %Integration Weight does not need updating
              Particles(IntGlo)%IntegrationWeight_DirDeriv_ImplicitIter = Particles(IntGlo)%IntegrationWeight_ImplicitIter * (1.0 + VolumetricStrain)!_ImplicitIter
            
             end do !loop over material points
            end do ! loop over elements
            
          end if ! 1-layer formulation or 2-layer form with 1 Phase
          
         
        end subroutine DynUpdateParticleWeights_DirDeriv_ImplicitIter
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #18 - BOOKMARK
        
        subroutine GetImplicitResidualEquation_SK04(LumpedMassDry_NPlus1, LumpedMassDry, &  
                                                            ExtLoad_NPlus1, ExtLoad, & 
                                                        GravityLoad_NPlus1, GravityLoad, &  
                                                            IntLoad_NPlus1, IntLoad, & 
                                                  TotalVelocitySoil_NPlus1, TotalVelocitySoil, &
                                                          ImplicitResidual) !ImplicitGaugeParameter, &
          !------------------------------------------------------
          ! purpose: implicit residual function
          
         implicit none 
         
         ! mass --> input 
         real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: LumpedMassDry_NPlus1 
         real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: LumpedMassDry
         
          ! internal force --> input 
          real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: ExtLoad_NPlus1
          real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: ExtLoad
          
          ! external force --> input 
          real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: IntLoad_NPlus1
          real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: IntLoad
          
          ! gravity force --> input 
          real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: GravityLoad_NPlus1
          real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: GravityLoad
          
          ! NPlus1 velocity --> input 
          real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: TotalVelocitySoil_NPlus1
          real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: TotalVelocitySoil
          
          ! N velocity --> input 
          !real(REAL_TYPE), intent(in), dimension(Counters%N, Counters%nEntity) :: TotalVelocitySoil
          
		  ! residual --> output 
		  real(REAL_TYPE), intent(out), dimension(Counters%N, Counters%nEntity) :: ImplicitResidual
          
          !real(REAL_TYPE), intent(in) :: ImplicitGaugeParameter
          
          ! local terms 
          integer(INTEGER_TYPE) :: IDof
          integer(INTEGER_TYPE) :: iEntity
		  
          !iEntity_Temporary =1
          
         do IEntity = 1, Counters%nEntity ! IEntity loop
         
             do IDof = 1, Counters%N ! IDof loop 
             
                 ImplicitResidual(IDof, IEntity) = ( (LumpedMassDry_NPlus1(IDof,IEntity) * TotalVelocitySoil_NPlus1(IDof,IEntity)) + &
                                               - (LumpedMassDry(IDof,IEntity) * TotalVelocitySoil(IDof,IEntity)) + &
                                               - (CalParams%TimeIncrement * CalParams%ImplicitIntegration%ImplicitGaugeParameter * &
                                                   (ExtLoad_NPlus1(IDof,IEntity) + GravityLoad_NPlus1(IDof,IEntity) - IntLoad_NPlus1(IDof,IEntity)) ) + &
                                               - (CalParams%TimeIncrement * (1-CalParams%ImplicitIntegration%ImplicitGaugeParameter) * &
                                                   (ExtLoad(IDof,IEntity) + GravityLoad(IDof,IEntity) - IntLoad(IDof,IEntity)) ) ) * PBoundary(IDof)!, IEntity)
         
             end do ! IDof loop 
         end do ! IEntity loop 
         
         
                                                          end subroutine GetImplicitResidualEquation_SK04  
                                                          
        !------------------------------------------------------------------------------------ 
        !------------------------------------------------------------------------------------  
                                                
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #19 - BOOKMARK
        
        
        !%---------------------------------------------------------------------%
        !%                                GMRES                                %
        !%---------------------------------------------------------------------%
        !---this is the mother subroutine where we use the directional derivative in the GMRES algorithm 
        subroutine GMRES_SK04()  

          implicit none 
          ! -- Input -------------------------------------------------------------------------
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_NPlus1 !, intent(inout)
          real(REAL_TYPE), dimension(Counters%N) :: DeltaTotalVelocitySoil !,Counters%nEntity
          real(REAL_TYPE) :: DeltaTotalVelocitySoilMagnitude
          real(REAL_TYPE) :: TotalVelocitySoilMagnitude_DirDeriv_ImplicitIter
          real(REAL_TYPE) :: TotalVelocitySoilMagnitude_ImplicitIter
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: LumpedMassDry_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: ExtLoad_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: GravityLoad_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: ImplicitResidual_DirDeriv
          real(REAL_TYPE), dimension(Counters%N) :: ProdJacobDeltaVeloDirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: SS
          
          
          			
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum !input
            
          ! -----------------------------------------------------------------------------------
          ! -- Local --------------------------------------------------------------------------
          ! -----------------------------------------------------------------------------------
          
          
          
          !---inout variables
          ! should this be inout or just output??
          integer(INTEGER_TYPE) :: mm_max_iteration !--> this needs to be local !, intent(inout)
          !---local variable
          integer(INTEGER_TYPE) :: ii, jj, kk, nn
          !----local allocatable GMRES variables 
          real(REAL_TYPE), allocatable, dimension(:) :: sn ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: cs ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: e1 ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: beta ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: HH_Arnoldi ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: ee_saved ! error saved 
          real(REAL_TYPE), allocatable, dimension(:) :: YY
          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ ! Rank 2 matrix?   
          real(REAL_TYPE), allocatable, dimension(:,:) :: HH_Arnoldi_Full
          !real(REAL_TYPE), allocatable, dimension(:,:) :: QQ_Arnoldi ! vectors? --> this needs to be assembled in QQ
          real(REAL_TYPE), dimension(Counters%N,1) :: QQ_Arnoldi
          real(REAL_TYPE), allocatable, dimension(:,:) :: HH_Arnoldi_New ! vectors?
          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ_Full
          real(REAL_TYPE), allocatable, dimension(:,:) :: HH_Full
          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ_Current
          real(REAL_TYPE), allocatable, dimension(:,:) :: HH_Current
          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ_Stored
          real(REAL_TYPE), allocatable, dimension(:,:) :: HH_Stored
          real(REAL_TYPE), allocatable, dimension(:,:) :: HH_inverted
          !----local GMRES variables
          real(REAL_TYPE) :: bb_norm 
          !real(REAL_TYPE) :: ImplicitResidualMagnitude
          real(REAL_TYPE) :: error 
          real(REAL_TYPE) :: ee, VerySmallNumber
          integer(INTEGER_TYPE) :: INFO
          
          real(REAL_TYPE) :: h_DirDerivParameter_SK04
          
          integer(INTEGER_TYPE) :: IIteration
          
          logical :: restart_gmres 
          
          restart_gmres = .false.
          VerySmallNumber = 1e-10
          !---the residual is already evaluated so we do not need to evaluate this again here 
          !---using norm2 function from FOTRAN90. Inputing rank as 1, but this is not necessary (i.e., optional)
          
          !---we define the error so many times. Here is it a function of the magnitude but we have to settle on a particular option. 
          error = ImplicitResidualMagnitude_ImplicitIter !/ !bb_norm ! real scalar in the numerator and denominator
          
          !---initialize the 1D vectors
          allocate(sn(CalParams%ImplicitIntegration%MaxGMRESIterations))
          allocate(cs(CalParams%ImplicitIntegration%MaxGMRESIterations))
          allocate(e1(CalParams%ImplicitIntegration%MaxGMRESIterations+1))
          allocate(beta(CalParams%ImplicitIntegration%MaxGMRESIterations+1))
          allocate(QQ(Counters%N,Counters%nEntity))
          
          cs = 0 
          sn = 0 
          e1 = 0 
          beta = 0
          QQ = 0 
          
          e1(1) = 1 ! assigning the first element of ee to be 1
          ee = error ! ee is a scalar. Set it equal to the value of error
          
          !---this would be p_1 in SK04 Algorithm 4 step 1
          QQ(:,1) = (-1/ImplicitResidualMagnitude_ImplicitIter) * ImplicitResidual_ImplicitIter(:,1) ! assigning this for the first iteration 
          
          ! equate DeltaTotalVelocitySoil to QQ(:,1) initially
          DeltaTotalVelocitySoil = QQ(:,1)
          DeltaTotalVelocitySoilMagnitude = norm2(DeltaTotalVelocitySoil)
          
          beta = ImplicitResidualMagnitude_ImplicitIter * e1 ! e1 and beta are the same size vectors 
          
          nn = 1
          
          IIteration = 1
          HH_Arnoldi = 0 
          QQ_Arnoldi = 0
          HH_Arnoldi_New = 0 
          
          allocate(QQ_Stored(Counters%N,IIteration))  !kk+1
          allocate(HH_Stored(IIteration+1,IIteration))
          
          QQ_Stored = 0
          HH_Stored = 0
          
          TotalVelocitySoilMagnitude_ImplicitIter = norm2(TotalVelocitySoil_ImplicitIter) ! v in SK04
          
          ! start of GMRES loop 
          do IIteration = 1, CalParams%ImplicitIntegration%MaxGMRESIterations ! loop accross the number of maximum iterations !kk
              
              ! if allocated we need to deallocate here
              if (allocated(QQ_Full)) then 
                  deallocate(QQ_Full)
              end if 
              
              ! if allocated we need to deallocate here
              if (allocated(HH_Full)) then 
                  deallocate(HH_Full)
              end if
              
              allocate(QQ_Full(Counters%N,IIteration+1)) !2
              allocate(HH_Full(IIteration+1,IIteration))
              
              
              QQ_Full = 0
              HH_Full = 0
              
              QQ_Stored(:,1) = QQ(:,1) ! this is a bit redundant in steps above k=1
              QQ_Full(:,1:IIteration) = QQ_Stored(:,1:IIteration) ! QQ_Full in column kk+1 should be left empty
              
              ! if larger than 2 then we need to store HH
              if (IIteration .ge. 2) then 
                  HH_Full(1:IIteration,1:IIteration-1) = HH_Stored(:,:)
              end if
              
              ! equate DeltaTotalVelocitySoil to QQ(:,1) initially -> update every GMRES step
              DeltaTotalVelocitySoil = QQ_Full(:,IIteration)
              DeltaTotalVelocitySoilMagnitude = norm2(DeltaTotalVelocitySoil)
              
              ! SUBROUTINE #1: Checked
              call GetDirDerivVelocity_SK04_ImplicitIter(TotalVelocitySoil_ImplicitIter, TotalVelocitySoilMagnitude_ImplicitIter, & !input 
                                                        DeltaTotalVelocitySoil, DeltaTotalVelocitySoilMagnitude)!, & !input
              
              call ApplyMPPrescribedVelocity_DirDeriv_ImplicitIter(TotalVelocitySys)
              call ApplyNodalPrescribedVelocity_DirDeriv_ImplicitIter(TotalVelocitySys)
              
              ! SUBROUTINE #2: Checked
              ! ------ I STOPPED HERE ------------
              call GetMPStrainFromNodeVelocity_DirDeriv_ImplicitIter(TotalVelocitySoil_DirDeriv_ImplicitIter)!, ShapeValuesArray_ImplicitIter, DShapeValuesArray_ImplicitIter, ActiveElement)                                           
              
              ! SUBROUTINE #3:
              call MPMDYNGetSig_DirDeriv_ImplicitIter()
                                !SigmaEffArray_N, SigmaEff0Array, ActiveElement_N, &!NumberOfIntegrationPoints, &
                                !MaterialIDArray_N, NPartEle_N, ESMStatevArray_N) !NAEl
              
              ! SUBROUTINE #4:
              call DynUpdateParticleWeights_DirDeriv_ImplicitIter() !DynUpdateParticleWeightsImplicit() !ShapeValuesArray_N, ActiveElement_N, NPartEle_N )
              
              ! SUBROUTINE #
              call UpdateParticlePos_DirDeriv_ImplicitIter()
              
              ! SUBROUTINE #5:
              call UpdateParticleHouseKeeping_DirDeriv_ImplicitIter()!NPartEle_N, EleParticlesHelp_N, IsActiveElement_N, EleParticles_N) 
              ! --> we update NAEl_ImplicitIter
              call SetActiveElement_DirDeriv_ImplicitIter()!ActiveElement_N, ActiveNode_N, IsActiveElement_N)
              call SetParticleIndex_DirDeriv_ImplicitIter()!(ActiveElement_N,NPartEle_N, EleParticlesHelp_N, EleParticles_N, GetParticleIndex_N)
              call CheckFillingOfElements()!(ActiveElement_N, MaterialIDArray_N, NPartEle_N)     !NPartEle, EleParticlesHelp, IsActiveElement)
              call SetUpEntityElements_DirDeriv_ImplicitIter()!ActiveElement_N, NPartEle_N, GetParticleIndex_N)
              call SetUpMaterialElements_DirDeriv_ImplicitIter()!(ActiveElement_N, MaterialIDArray_N, NPartEle_N) !--> why are we even doing this???
              
              ! SUBROUTINE #12:
              call GetNodalIntForcesAfter_DirDeriv_ImplicitIter(IntLoad_DirDeriv_ImplicitIter, SigmaEffArray_DirDeriv_ImplicitIter)  
              !call CalculateIntAndExtWorks_DirDeriv_ImplicitIter(IncrementalDisplacementWater)!, IntLoad, IntLoad_NPlus1, SigmaEffArray_N)
              
              ! SUBROUTINE #14
              call MapMassMP2Nodes_DirDeriv_ImplicitIter(LumpedMassDry_DirDeriv_ImplicitIter)!, ShapeValuesArray_NPlus1, ActiveElement_N)
              
              ! SUBROUTINE #13:
              call GetNodalExternalForces_DirDeriv_ImplicitIter(ExtLoad_DirDeriv_ImplicitIter, GravityLoad_DirDeriv_ImplicitIter)!, ActiveElement_N) ShapeValuesArray_NPlus1
              
              call ApplyMPPrescribedVelocity_IntLoad(IntLoad_DirDeriv_ImplicitIter)
              call ApplyNodalPrescribedVelocity_IntLoad(IntLoad_DirDeriv_ImplicitIter)
              
              
              call GetImplicitResidualEquation_SK04(LumpedMassDry_DirDeriv_ImplicitIter, LumpedMassDry_DirDeriv_ImplicitIter, &!LumpedMassDry_DirDeriv_ImplicitIter, LumpedMassDry_DirDeriv_ImplicitIter, & !input !LumpedMassDry_DirDeriv_ImplicitIter, LumpedMassDry, & !input 
                                                          ExtLoad_DirDeriv_ImplicitIter, ExtLoad, &  !input
                                                      GravityLoad_DirDeriv_ImplicitIter, GravityLoad, &  !input
                                                          IntLoad_DirDeriv_ImplicitIter, IntLoad, &  !input
                                                TotalVelocitySoil_DirDeriv_ImplicitIter, TotalVelocitySoil, &  !input
                                                 ImplicitResidual_DirDeriv_ImplicitIter)  !output
              
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        allocate(HH_Arnoldi(IIteration)) ! 1 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        QQ_Arnoldi = 0 
        HH_Arnoldi = 0 
        
              ! SUBROUTINE #18: 
              call GetProdJacobDeltaVeloDirDeriv(ProdJacobDeltaVeloDirDeriv, DeltaTotalVelocitySoilMagnitude, &
                                                       TotalVelocitySoilMagnitude_ImplicitIter, & !TotalVelocitySoil, 
                                                        ImplicitResidual_DirDeriv_ImplicitIter, ImplicitResidual_ImplicitIter)!, &
                                                        !h_DirDerivParameter_SK04)
              
              ! SUBROUTINE #19: p_k+1 = DirectionalDerivative(AccelerationSoil_NPlus1, OldImplicitResidual))  
              QQ_Arnoldi(:,1) = ProdJacobDeltaVeloDirDeriv!DirectionalDerivative
              
              ! SUBROUTINE #20: THIS SHOULD BE QQ(:,kk) --> It should depend on the timestep 
              !---step 2b h_k,i = p_k+1 * p_i 
              do ii = 1, IIteration !% Modified Gram-Schmidt, keeping the Hessenberg matrix            
                  do jj = 1,Counters%N !2 
                      HH_Arnoldi(ii) = HH_Arnoldi(ii) + ( QQ_Arnoldi(jj,1) * QQ_Full(jj,ii) ) !QQ_Arnoldi is the p_k+1 term. QQ is the p_i term 
                  end do 
              end do 
              
              
              ! SUBROUTINE #21: !---step 2b in algorithm 4 in SK04
              !----Equation p_k+1 = p_k+1 - (h_k,i * p_i)
              do ii = 1, IIteration !% Modified Gram-Schmidt, keeping the Hessenberg matrix
                  QQ_Arnoldi(:,1) = QQ_Arnoldi(:,1) - ( HH_Arnoldi(ii) * QQ_Full(:,ii) ) !QQ_Arnoldi is the p_k+1 term. QQ is the p_i term 
              end do 
              
              ! SUBROUTINE #22: this needs to be stored outside this loop within the entire module as a global variable 
              if (allocated(HH_Arnoldi_New)) then 
                  deallocate(HH_Arnoldi_New)
              end if 
              allocate( HH_Arnoldi_New (IIteration+1, 1) )
              HH_Arnoldi_New = 0
              HH_Arnoldi_New(1:IIteration, 1) = HH_Arnoldi(1:IIteration) 
              
              if (allocated(HH_Arnoldi)) then 
                  deallocate(HH_Arnoldi)
              end if 
              
              ! SUBROUTINE #23: 
              !---step 1c h_k+1,k = |p_k+1|
              HH_Arnoldi_New(IIteration+1, 1) =  norm2(QQ_Arnoldi(:,1),1)
          
              
              if (abs(HH_Arnoldi_New(IIteration+1, 1)) == 0.0) then
                  !IIteration_PreventOrtho = Iteration -1 
                  restart_gmres = .true.
                  exit 
              else 
                  !---step 1d p_k+1 = p_k+1/|p_k+1| 
                  ! if HH_Arnoldi_New(IIteration+1, 1) is zero it means that we have reached orthogonality and the solution cannot be further improved
                  QQ_Arnoldi = QQ_Arnoldi * (1/HH_Arnoldi_New(IIteration+1, 1) ) ! divide by a scalar
          
              end if
              
              !---pouring q (or called as p_k+1 in SK04) back to QQ_Full to store    
              QQ_Full(:,IIteration+1) = QQ_Arnoldi(:,1)
              
              ! SUBROUTINE #24: eliminate the last element in H ith row and update the rotation matrix 
              call apply_givens_rotation( HH_Arnoldi_New, cs, sn, IIteration, CalParams%ImplicitIntegration%MaxGMRESIterations ) ! HH_Arnoldi_New compiles both old and new
              
              ! SUBROUTINE #25: 
              if (allocated(HH_Arnoldi_Full)) then 
                  deallocate(HH_Arnoldi_Full)
              end if 
              allocate(HH_Arnoldi_Full(nn,IIteration))
              HH_Arnoldi_Full = 0
              HH_Arnoldi_Full = HH_Arnoldi_New
              
              ! SUBROUTINE #26: update the residual vector 
              beta(IIteration+1) = -sn(IIteration) * beta(IIteration)
              beta(IIteration) = cs(IIteration) * beta(IIteration)
              
              ! SUBROUTINE #27: evaluate error
              error = abs(beta(IIteration+1)) !/ bb_norm
              
              ! SUBROUTINE #28: allocate/deallocate QQ and HH
              deallocate(QQ_Stored)  !kk+1
              deallocate(HH_Stored)
              allocate(QQ_Stored(Counters%N,IIteration+1))  !kk+1
              allocate(HH_Stored(IIteration+1,IIteration))
              QQ_Stored = QQ_Full
              
              
              ! SUBROUTINE #29: 
              if (IIteration .ge. 2) then 
                  HH_Stored = 0
                  HH_Stored(:, 1:IIteration-1) = HH_Full(:, 1:IIteration-1)
                  HH_Stored(:, IIteration) = HH_Arnoldi_New(:,1)
              else 
                  HH_Stored = 0          
                  HH_Stored = HH_Arnoldi_New
              end if 
              
              ! SUBROUTINE #30: check error against threshold
              if (error .le. CalParams%ImplicitIntegration%GMRESErrorThreshold ) then 
                  exit
              end if 
          
          end do ! loop over mm maximum iterations
          
          
          if (restart_gmres == .true.) then 
              IIteration = IIteration -1
          end if
          if (IIteration == CalParams%ImplicitIntegration%MaxGMRESIterations+1) then 
          
              IIteration = CalParams%ImplicitIntegration%MaxGMRESIterations!IIteration -1 
              
          end if 
          
          
          ! SUBROUTINE #31: allocate HH based on the number of iterations kk
          allocate(HH_inverted(IIteration,IIteration))
          HH_inverted = 0.0
          
          ! SUBROUTINE #32: find the inverse of upper triangular H matrix
          call L_inv(HH_Stored(1:IIteration, 1:IIteration), & ! in variable, 
                     HH_inverted, & ! out variable 
                     IIteration)  ! in variable  
          
          ! SUBROUTINE #33: GMRES Step 3 to solve the system Hy=beta for y
          allocate(YY(IIteration)) 
          YY = 0 
          
          ! SUBROUTINE #34: conduct the multiplication to find YY 
          ! loop over rows  
          do ii = 1, IIteration ! Modified Gram-Schmidt, keeping the Hessenberg matrix
              ! loop over columns 
              do jj = 1, IIteration 
                  YY(ii) = YY(ii) + HH_inverted(ii, jj) * beta(jj)   !HH_Stored
              end do 
          end do 
          
          ! SUBROUTINE #35: GMRES Step 4 to Calculate the solution s = s_0 + y1*p_1 + ... + y_k*p_k
          ! xx = xx + (QQ_Full(:,1:kk) * YY)
          ! loop over rows 
          !do ii = 1, Counters%N ! rank of matrix 
          SS = 0
          do jj = 1, IIteration
          
              SS(:,1) = SS(:,1) + QQ_Full(:,jj)*YY(jj)    
          end do
          !end do
          
          ! project total velocity soil
          TotalVelocitySoil_ImplicitIter(:,Counters%nEntity) = TotalVelocitySoil_ImplicitIter(:,Counters%nEntity) & !TotalVelocitySoil(:,Counters%nEntity) &
              + SS(:,1) ! this becomes the new velocity
          
          
          !TotalVelocitySoil = TotalVelocitySoil_ImplicitIter
          !do ii = 1, Counters%NDof
          !    
          !    
          !    
          !end do 
          !do ii = 1, Counters%N ! loop over all degrees of freedom
          !          do jj = 1, Counters%nEntity ! loop over all entities
          !      
          !      !if (LumpedMassDry(IDOF,J)/=0) then
          !      !  TotalVelocitySoil(IDOF,J) = ( Momentum(IDOF,J) / LumpedMassDry(IDOF,J) ) * PBoundary(IDOF)
          !      !else
          !      !  TotalVelocitySoil(IDOF,J) = 0.0
          !      !end if
          !        
          !      TotalVelocitySoil_ImplicitIter(ii,jj) = TotalVelocitySoil_ImplicitIter(ii,jj) * PBoundary(ii)
          ! 
          !          end do ! loop over all entities
          !  end do ! loop over all degrees of freedom
          
          ! Apply material point prescribed velocity
            !call ApplyMPPrescribedVelocity_ImplicitIter(TotalVelocitySys)
          
          end subroutine GMRES_SK04                                                  
                      
          !--------------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------
          
          ! SUBROUTINE #21 - BOOKMARK
        
        
            subroutine GetDirDerivVelocity_SK04_ImplicitIter(TotalVelocitySoil, TotalVelocitySoilMagnitude, &
                                            DeltaTotalVelocitySoil, DeltaTotalVelocitySoilMagnitude)!, &
                                                
            
                                                  !h_DirDerivParameter_SK04, &
                                                  !TotalVelocitySoil_DirDeriv)
              ! purpose: get the velocity to be used in the directional derivative equation 
              implicit none 
              !---inputs
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: TotalVelocitySoil
              real(REAL_TYPE), intent(in) :: TotalVelocitySoilMagnitude
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: DeltaTotalVelocitySoil
              real(REAL_TYPE), intent(in) :: DeltaTotalVelocitySoilMagnitude
              !real(REAL_TYPE), intent(in) :: h_DirDerivParameter_SK04
              !---outputs
              !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: TotalVelocitySoil_DirDeriv
              !---local
              real(REAL_TYPE) :: VerySmallNumber
              
              VerySmallNumber = 1e-10!10
             
              ! this subroutine find out what is the value of the velocity on the left hand side of the the 
              !     directional derivative approximation. It is always TotalVelocitySoil + ... with value of
              !     TotalVelocitySoil okay being zero
              
              ! SUBROUTINE #18: 
              !---case 1
              if (DeltaTotalVelocitySoilMagnitude .lt. VerySmallNumber) then ! equals to zero     
                  TotalVelocitySoil_DirDeriv_ImplicitIter = TotalVelocitySoil ! this was not written in SK04
                  
              !---case 2    
              elseif ( (DeltaTotalVelocitySoilMagnitude .gt. VerySmallNumber) .and. &
              (TotalVelocitySoilMagnitude .gt. VerySmallNumber) ) then ! greater than zero 
                  TotalVelocitySoil_DirDeriv_ImplicitIter = TotalVelocitySoil + &
                      ( CalParams%ImplicitIntegration%HH_DirDerivParameter_SK04 * TotalVelocitySoilMagnitude * DeltaTotalVelocitySoil * (1/DeltaTotalVelocitySoilMagnitude) )
                  
              !---case 3
              elseif ( (DeltaTotalVelocitySoilMagnitude .gt. VerySmallNumber) .and. &
              (TotalVelocitySoilMagnitude .lt. VerySmallNumber) ) then 
                  TotalVelocitySoil_DirDeriv_ImplicitIter = &!TotalVelocitySoil + & ! this should be TotalVelocitySoil but since it is small = approx zero
                      ( CalParams%ImplicitIntegration%HH_DirDerivParameter_SK04 * DeltaTotalVelocitySoil * (1/DeltaTotalVelocitySoilMagnitude) )
                  
              end if 
              
              
                                            
              end subroutine GetDirDerivVelocity_SK04_ImplicitIter
        
        !--------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------
                                                  
                                                  
                                                  
                                                  
                                                  
         !--------------------------------------------------------------------------------------------------
         !--------------------------------------------------------------------------------------------------                           
         !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #20 - BOOKMARK
                      
                      
                      
        subroutine GetProdJacobDeltaVeloDirDeriv(ProdJacobDeltaVeloDirDeriv, DeltaTotalVelocitySoilMagnitude, &
                                                       TotalVelocitySoilMagnitude, & !TotalVelocitySoil, 
                                                        ImplicitResidual_DirDeriv, ImplicitResidual)!, &
                                                        !h_DirDerivParameter_SK04)
              ! purpose: to evaluate ProdJacobDeltaVeloDirDeriv
              !
              implicit none 
              
              
              !---Input
              real(REAL_TYPE), intent(in) :: DeltaTotalVelocitySoilMagnitude
              real(REAL_TYPE), intent(in) :: TotalVelocitySoilMagnitude
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: ImplicitResidual_DirDeriv
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: ImplicitResidual
              !real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: TotalVelocitySoil
              !real(REAL_TYPE), intent(in) :: h_DirDerivParameter_SK04
              !---Output
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: ProdJacobDeltaVeloDirDeriv
              !---Local
              real(REAL_TYPE) :: VerySmallNumber
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: Zeros_Vector !, intent(in)
              integer(INTEGER_TYPE) :: ii
              
              ! setting up local variables 
              VerySmallNumber = 1e-10!10
              Zeros_Vector = 0.0
              
              !---case 1
              if (DeltaTotalVelocitySoilMagnitude .lt. VerySmallNumber) then ! equals to zero     
                  ProdJacobDeltaVeloDirDeriv = Zeros_Vector !0.0
                
              !---case 2    
              elseif ( (DeltaTotalVelocitySoilMagnitude .gt. VerySmallNumber) .and. (TotalVelocitySoilMagnitude .gt. VerySmallNumber) ) then ! greater than zero 
                  !---Directional derivative calculation
                  do ii = 1,Counters%N
                      ProdJacobDeltaVeloDirDeriv(ii,1) = DeltaTotalVelocitySoilMagnitude * & !<-- this is p_k
                      (ImplicitResidual_DirDeriv(ii,1) - ImplicitResidual(ii,1)) * &
                          (1/CalParams%ImplicitIntegration%HH_DirDerivParameter_SK04) * (1/TotalVelocitySoilMagnitude)  * PBoundary(ii) 
                  end do
                  
              !---case 3
              elseif ( (DeltaTotalVelocitySoilMagnitude .gt. VerySmallNumber) .and. (TotalVelocitySoilMagnitude .lt. VerySmallNumber) ) then 
              
                  !---Directional derivative calculation
                  !----the reason there is a 1 is because the directional derivative is a vector (rank 1). But Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside is rank 2.  
                  do ii = 1,Counters%N
                      ProdJacobDeltaVeloDirDeriv(ii,1) = DeltaTotalVelocitySoilMagnitude * &
                          (ImplicitResidual_DirDeriv(ii,1) - ImplicitResidual(ii,1)) * & ! I removed generalizedalpha naming from here...
                          (1/CalParams%ImplicitIntegration%HH_DirDerivParameter_SK04)  * PBoundary(ii)
                  end do
              end if 
              
                                                        end subroutine GetProdJacobDeltaVeloDirDeriv                   
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------                     
                                                          
                                                        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #22 - BOOKMARK
        !%---------------------------------------------------------------------%
        !%                  Applying Givens Rotation to H col                  %
        !%---------------------------------------------------------------------%
        subroutine apply_givens_rotation(HH_Arnoldi, cs, sn, kk, mm_max_iteration)
        
        integer(INTEGER_TYPE), intent(in) :: mm_max_iteration
        integer(INTEGER_TYPE), intent(in) :: kk ! current number of iteration
        
        integer(INTEGER_TYPE) :: ii
        
        real(REAL_TYPE),  intent(inout), dimension(mm_max_iteration) :: cs ! vector (i.e., rank 1) !allocatable,
        real(REAL_TYPE),  intent(inout), dimension(mm_max_iteration) :: sn ! vector (i.e., rank 1) !allocatable,
        real(REAL_TYPE),  intent(inout), dimension(kk+1) :: HH_Arnoldi ! vector (i.e., rank 1) !allocatable,
        
        
        real(REAL_TYPE) :: temp
              
              
        
        ! apply for ith column --> apply previous Givens rotations
        do ii = 1,kk-1
            temp          = (  cs(ii) * HH_Arnoldi(ii) ) + ( sn(ii) * HH_Arnoldi(ii+1) )
            HH_Arnoldi(ii+1)  = ( -sn(ii) * HH_Arnoldi(ii) ) + ( cs(ii) * HH_Arnoldi(ii+1) )
            HH_Arnoldi(ii)    = temp    
        end do 
        
        ! update the next sin cos values for rotation --> form current givens rotation
        ![cs_k, sn_k] = givens_rotation(h(k), h(k + 1));
        call givens_rotation(HH_Arnoldi(kk), HH_Arnoldi(kk + 1), cs(kk), sn(kk))  
        
        ! eliminate HH_Arnoldi(ii+1, ii) 
        HH_Arnoldi(kk) = ( cs(kk) * HH_Arnoldi(kk) ) + ( sn(kk) * HH_Arnoldi(kk+1) ) ! rewriting the value of HH_Arnoldi(kk)
        HH_Arnoldi(kk+1) = 0.0
        
        end subroutine apply_givens_rotation
        
        
        
        
        
        
        
        
        
        
        !%%----Calculate the Givens rotation matrix----%%
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #23 - BOOKMARK
        
        
        
        subroutine givens_rotation(v1, v2, cs, sn)
        
        !subroutine givens_rotation(HH_Arnoldi_kk, HH_Arnoldi_kkPlus1, cs, sn)
        
        implicit none 
        
        real(REAL_TYPE) :: tt
        real(REAL_TYPE), intent(in) :: v1
        real(REAL_TYPE), intent(in) :: v2
              
        real(REAL_TYPE), intent(out) :: cs
        real(REAL_TYPE), intent(out) :: sn
        
        tt = (v1**2 + v2**2)**0.5 ! finding the norm of the input vector 
        
        cs = v1/tt !% see http://www.netlib.org/eispack/comqr.f
        
        sn = v2/tt
        
        end subroutine givens_rotation
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #24 - BOOKMARK
        
        
        
        
        
        
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%---------------------------------------------------------------------%
        !%                      Inverting lower triangular matrix              %
        !%---------------------------------------------------------------------%
        
        subroutine L_inv(A, & ! in variable, 
                         L, & ! out variable 
                         kk)  ! in variable  
        
        implicit none 
        
        real(REAL_TYPE), dimension(kk,kk), intent(in) :: A
        real(REAL_TYPE), dimension(kk,kk) :: A_Transposed
        real(REAL_TYPE), dimension(kk,kk), intent(out) :: L
        integer(INTEGER_TYPE) , intent(in) :: kk
        integer(INTEGER_TYPE) :: n, i, j, oo
            
        n = kk
        L = 0.0 
           
        
        A_Transposed = transpose(A)
        
        do i=1,kk
            L(i,i) = 1/A(i,i)
            do j=i+1,kk
                 !n = kk
                do oo = i,j-1 
                L(j,i)=L(j,i)-A_Transposed(j, oo)*L(oo,i)!*(1/A_Transposed(j,j))
                !L(j,i)=-A_Transposed(j, i:j-1)*L(i:j-1,i)*(1/A_Transposed(j,j))
                end do 
                L(j,i)= L(j,i) * (1/A_Transposed(j,j))
            end do       
        end do            
        
        L = transpose(L)
   
        end subroutine L_inv
        
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
               
                         
                                                                   
                                                  
                                                  
                                                  
                                                          
                                                          
      end module ModDynamicImplicit
