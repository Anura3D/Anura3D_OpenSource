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
	  
	  
	  module ModDynamicExplicit
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
        
        real(REAL_TYPE) :: ImplicitResidual_Magnitude
        logical :: IsEnterGMRES
        
        real(REAL_TYPE) :: ImplicitGMRESThreshold
        integer(INTEGER_TYPE) :: mm_max_iteration, ii
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum       !CC added nEntity !, intent(out)
        real(REAL_TYPE) :: ImplicitResidualMagnitude
        real(REAL_TYPE) :: ImplicitGaugeParameter
        
        ImplicitGaugeParameter = 1
        
        ImplicitGMRESThreshold = 1e-7
        
        mm_max_iteration = MAX_ITERATIONS_IMPLICIT 
        

          !********** 4b - TIME STEP LOOP ******************************
          !do while( (.not.CalParams%ConvergenceCheck%DoesConverge) .and. (.not.CalParams%ConvergenceCheck%DoesDiverge) )

            ! Increase time -> explicit style 
            CalParams%TimeStep = CalParams%TimeStep + 1
            CalParams%TotalRealTime = CalParams%TotalRealTime + CalParams%TimeIncrement
            CalParams%OverallRealTime = CalParams%OverallRealTime + CalParams%TimeIncrement
            call WriteInLogFile('TimeStep: ' // trim(String(CalParams%TimeStep))//' '//trim(String(CalParams%TimeIncrement)) &
                                                                                //' '//trim(String(CalParams%TotalRealTime)))
            call WriteInLogFile('  Skipconvection? '//trim(String(MinimumDeterminantRatioReached))//' '// &
                                                      trim(String(IsMPMSkipConvection())))

            !call UpdateMultipliersForTimeDepencency() 
            
            ! Sulsky and Kaul 2004 - implicit scheme formulation 
            !if ((CalParams%ImplicitIntegration%IsImplicitFormulationSK04 == .true.) .and. (CalParams%ImplicitIntegration%IsImplicitScheme == .true.)) then     
            
            ! SUBROUTINE #1
            call MapMassMP2Nodes(MassArray, LumpedMassDry_N, ShapeValuesArray_N, &
                                  ActiveElement_N) 
                                                !--> O: LumpedMassDry_N
                                                !--> I: MassDry
                                                !--> I: ShapeValuesArray_N
            ! SUBROUTINE #2
            call MapMomentumMP2Nodes(MassArray, Momentum_N, VelocityArray_N, ShapeValuesArray_N, ElementConnectivities, &
                                    ActiveElement_N)            
                                                !--> O: Momentum_N
                                                !--> I: MassDry
                                                !--> I: VelocityArray
                                                !--> I: ShapeValuesArray                       
            ! SUBROUTINE #3
            call GetNodalVelocityFromNodalMomentum(Momentum_N, LumpedMassDry_N, TotalVelocitySoil_N, PBoundary)
                                                !--> O: TotalVelocitySoil_N
                                                !--> I: Momentum_N
                                                !--> I: LumpedMassDry_N
                                                !--> I: PBoundary
            ! SUBROUTINE #4
            call GetNodalExternalForces(ShapeValuesArray_N, ExtLoad_N, GravityLoad_N, ActiveElement_N)
                                                !--> O: ExtLoad_N
                                                !--> O: GravityLoad_N
                                                !--> I: ShapeValuesArray_N            
            ! SUBROUTINE #5
            call GetNodalIntForces(IntLoad_N, SigmaEffArray_N, NPartEle_N, ActiveElement_N) ! SigmaEffArray comes from existing stress
            
            ! SUBROUTINE #6
            call GetFirstGuessOfVelocity(TotalVelocitySoil_NPlus1)
                                                !--> O: TotalVelocitySoil_NPlus1
            ! SUBROUTINE #7
            call GetMPStrainFromNodeVelocity(TotalVelocitySoil_NPlus1, ShapeValuesArray_N, DShapeValuesArray_N, ActiveElement_N)!, &
                                            !ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1)                                             
                                                !--> O: ShapeValuesArray_NPlus1
                                                !--> O: DShapeValuesArray_NPlus1
                                                !--> O: STRAIN?????
                                                !--> I: TotalVelocitySoil_NPlus1            
                                                !--> I: ShapeValuesArray_N     
                                                !--> I: DShapeValuesArray_N     
            
           ! SUBROUTINE #8: stresses, the water pressure and the gass pressure are updated (The mass balance equations are solved)
           call MPMDYNGetSig(SigmaEffArray_NPlus1, SigmaEff0Array, ActiveElement_N, & ! we are calculating new stress
                          MaterialIDArray_N, NPartEle_N)!(IntLoad_NPlus1) !GetParticleIndex,          !NAEl,
           
           ! SUBROUTINE #9: 
           call DynUpdateParticleWeights( ShapeValuesArray_N, ActiveElement_N, NPartEle_N)
                                
           call UpdateParticlePos(ShapeValuesArray_N, DShapeValuesArray_N, &
                                     ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
                                     GlobPosArray_N, GlobPosArray_NPlus1, &
                                     UStepArray, UPhaseArray, UArray, &
                                     ElementIDArray_N, EleParticles_N, &
                                    IncrementalDisplacementSoil_N)!ShapeValuesArray, DShapeValuesArray, ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
                                   !IncrementalDisplacementSoil, GlobPosArray, UStepArray, UPhaseArray, UArray, &
                                   !   ElementIDArray, EleParticles)
           
           ! SUBROUTINE #10:
                !Variables populated in UpdateParticleHouseKeeping:
			    !NPartEle
			    !EleParticleHelp
			    !IsActiveElement
           call UpdateParticleHouseKeeping(NPartEle_NPlus1, EleParticlesHelp_NPlus1, IsActiveElement_NPlus1, EleParticles_NPlus1)
                                                               
           ! SUBROUTINE #11:
                !Variables populated in SetActiveElement: 
                !ActiveElement
			    !ActiveNode
           call SetActiveElement(ActiveElement_NPlus1, ActiveNode_NPlus1)!(ActiveElement, ActiveNode)
                                                        
           ! SUBROUTINE #12: 
                !Variables populated in SetParticleIndex: 
			    !GetParticleIndex
           call SetParticleIndex(ActiveElement_NPlus1, NPartEle_NPlus1, EleParticlesHelp_NPlus1, EleParticles_NPlus1)
           
           ! SUBROUTINE #13: 
           call CheckFillingOfElements(ActiveElement_NPlus1, MaterialIDArray_NPlus1, NPartEle_NPlus1)     !NPartEle, EleParticlesHelp, IsActiveElement)
           
           ! SUBROUTINE #14: 
           call SetUpEntityElements(ActiveElement_NPlus1, NPartEle_NPlus1)
           
           ! SUBROUTINE #15: 
           call SetUpMaterialElements (ActiveElement_NPlus1, MaterialIDArray_NPlus1, NPartEle_NPlus1) !--> why are we even doing this???
           
           ! SUBROUTINE #16: 
           call SetInitialStressForNextLoadStep(ActiveElement_NPlus1, NPartEle_NPlus1, SigmaEffArray_NPlus1, SigmaEff0Array)
           
           ! SUBROUTINE #17: 
           !call CalculateIntAndExtWorks(IncrementalDisplacementWater_NPlus1, IntLoad_N, IntLoad_NPlus1, SigmaEffArray_N)
           call CalculateIntAndExtWorks(IncrementalDisplacementWater, IntLoad_N, IntLoad_NPlus1, SigmaEffArray_N)

           
            ! SUBROUTINE #8
            call GetNodalExternalForces(ShapeValuesArray_NPlus1, ExtLoad_NPlus1, GravityLoad_NPlus1, ActiveElement_N)
            
            ! SUBROUTINE #9
            call MapMassMP2Nodes(MassArray, LumpedMassDry_NPlus1, ShapeValuesArray_NPlus1, &
                                  ActiveElement_N)
            
            ! SUBROUTINE #10
            call GetImplicitResidualEquation_SK04(LumpedMassDry_NPlus1, LumpedMassDry_N, &  
                                                            ExtLoad_NPlus1, ExtLoad_N, & 
                                                        GravityLoad_NPlus1, GravityLoad_N, &  
                                                            IntLoad_NPlus1, IntLoad_N, & 
                                                  TotalVelocitySoil_NPlus1, TotalVelocitySoil_N, &
                                                    ImplicitGaugeParameter, &
                                                          ImplicitResidual_NPlus1)
            
            ! SUBROUTINE #11
            !ImplicitResidualMagnitude_NPlus1 = norm2(ImplicitResidual_NPlus1)
            !ImplicitResidualMagnitude_NPlus1 = norm2(ImplicitResidual_NPlus1)
            
            ! SUBROUTINE #12
            if (ImplicitResidual_Magnitude<ImplicitGMRESThreshold) then 
                !... we reached the solution
                CalParams%ImplicitIntegration%IsEnterGMRES = .false. 
            else 
                CalParams%ImplicitIntegration%IsEnterGMRES = .true.     
            end if 
                
            if (CalParams%ImplicitIntegration%IsEnterGMRES == .true.) then 
            
                !call GMRES_SK04()
                
                !call GMRES_DirectionalDerivative_SK04(& !AA, bb, xx, threshold, nn_length_AA, mm_max_iteration)
                !      LumpedMassDry, & !LumpedMassDry_NPlus1, & ! in variables 
                !      AccelerationSoil, & !Explicit velocity
                !      TotalVelocitySoil, &
                !      IncrementalDisplacementSoil, &
                !      AccumulatedIncDisplacementSoil, &
                !      TotalDisplacementSoil, &
                !      IntLoad, & !IntLoad_NPlus1, & ! in variables 
                !      ImplicitResidual, & ! in variables 
                !      AccelerationSoil_Increment(:,1), & ! in variable
                !      ImplicitGMRESThreshold, & ! in variable
                !      !nn_Length_AA, & ! in variable
                !      mm_max_iteration)!, & ! this is s_0 in SK04)
             
            end if 
                
                
                
                
                
            !end if 
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
            ! First stage effects 
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Lagrangian phase
            !call LagrangianPhase(LumpedMassDry, ShapeValuesArray, ExtLoad, GravityLoad) !--> gives us the base items of the residual equation
            
            ! obtain first guess of acceleration 
            !AccelerationSoil_NPlus1 = AccelerationSoil + AccelerationSoil_Increment ! Acceleration increment is projected to be zero to start with
            
            
            !call ApplyMPMDynamicContact() ! only if ApplyContactAlgorithm
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Part 1 assuming an acceleration increment of zero 
            ! Convective phase
            !call DYNConvectivePhase(LumpedMassDry, & ! input
            !                        IntLoad_NPlus1,  & ! output
            !                        AccelerationSoil_NPlus1, & ! input
            !                        TotalVelocitySoil_NPlus1, & ! input 
            !                        IncrementalDisplacementSoil_NPlus1, & ! input 
            !                        AccumulatedIncDisplacementSoil_NPlus1, &
            !                        TotalDisplacementSoil_NPlus1, & !input
            !                        ShapeValuesArray, DShapeValuesArray, & !input
            !                        ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1) ! out --> gives us the NPlus1 terms to be used to find the incremental terms
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                 
            ! Part 2 to update mass and external forces                        
            !call LagrangianPhase(LumpedMassDry_NPlus1, ShapeValuesArray_NPlus1, ExtLoad_NPlus1, GravityLoad_NPlus1) !--> gives us the base items of the residual equation
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            !CalculateIncrementalResidualTermsImplicitScheme
            !%step 4: calculate implicit scheme incremental terms
            !call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_NPlus1, LumpedMassDry_Increment, &
            !                                                     ExtLoad, ExtLoad_NPlus1, ExtLoad_Increment, &
            !                                                     IntLoad, IntLoad_NPlus1, IntLoad_Increment, &
            !                                                     GravityLoad, GravityLoad_NPlus1, GravityLoad_Increment, &
            !                                                     AccelerationSoil, AccelerationSoil_NPlus1, AccelerationSoil_Increment)
            ! calculating incremental terms for use in the residual equation based on the assumption of AccelerationSoil_Increment
            
            !%step 5: at the end of the lagrangian phase, we evaluate the residual assuming that everything has initialized
            ! --> this needs to happen after updating of the new terms 
            !call CalculateResidualImplicitScheme(LumpedMassDry, LumpedMassDry_Increment, & !LumpedMassDry_NPlus1, 
            !                                           ExtLoad, ExtLoad_Increment, & ! ExtLoad_NPlus1,
            !                                       GravityLoad, GravityLoad_Increment, & ! GravityLoad_NPlus1, 
            !                                           IntLoad, IntLoad_Increment,& ! IntLoad_NPlus1, 
            !                                  ImplicitResidual, &
            !                                  AccelerationSoil, AccelerationSoil_Increment)
            !
            !do ii = 1, Counters%N 
            !  ImplicitResidual(ii,1) = ImplicitResidual(ii,1) * PBoundary(ii) !--> applying boundary conditions
            !end do

        
            call FindVectorMagnitude(ImplicitResidual_NPlus1, Counters%N, ImplicitResidual_Magnitude) 
            
            if (ImplicitResidual_Magnitude<ImplicitGMRESThreshold) then 
                !... we reached the solution
                IsEnterGMRES = .false. 
            else 
                IsEnterGMRES = .true. 
            end if 
            
            if (IsEnterGMRES == .true.) then 
                
                    
                !call GMRES_DirectionalDerivative_AccelerationFormulation(& !AA, bb, xx, threshold, nn_length_AA, mm_max_iteration)
                !      LumpedMassDry, & !LumpedMassDry_NPlus1, & ! in variables 
                !      AccelerationSoil, & !Explicit velocity
                !      TotalVelocitySoil, &
                !      IncrementalDisplacementSoil, &
                !      AccumulatedIncDisplacementSoil, &
                !      TotalDisplacementSoil, &
                !      IntLoad, & !IntLoad_NPlus1, & ! in variables 
                !      ImplicitResidual, & ! in variables 
                !      AccelerationSoil_Increment(:,1), & ! in variable
                !      ImplicitGMRESThreshold, & ! in variable
                !      !nn_Length_AA, & ! in variable
                !      mm_max_iteration)!, & ! this is s_0 in SK04)
                    
                    
                
            
                AccelerationSoil = AccelerationSoil + AccelerationSoil_Increment
                AccelerationSoil_NPlus1 = AccelerationSoil
                
            end if 
            
            
            
            
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
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

            !  if (CalParams%OutputDebugData) then
            ! Time step output
            if (NFORMULATION==1) then
              call EnergyOutput()
            else
              call EnergyOutput2LayForm()
            end if

            call MaterialPointOutput() ! only if OutputNumberParticles>0
            
            !compute and write on file the sum of nodal reactions on selected surfaces
            call SurfaceReactionOutput() 
            
            ! end if

            ! update ElementLMin in case of ULFEM
            if (IsULFEMComputation()) then
              call DetermineElementLMin() ! calulate minimum element altitude
            end if

            call CalculateCriticalTimeStep()
            call WriteFEMNodeData() ! only if .not.IsMPMComputation

          !end do ! ----- end TIME STEP LOOP -----

        end subroutine RunImplicitDynamicLoadStep
        
        
        
        
        
        !--------------------------------------------------------------------------------------------        
        !--------------------------------------------------------------------------------------------        
        !--------------------------------------------------------------------------------------------      
        !--------------------------------------------------------------------------------------------
       
          !--------------------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------------
          

        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                               
                

        
        
        
        
        
        
        
        
        
        
        !-----------------------------------------------------------------------------------------------------
            
            !-----------------------------------------------------------------------------------------------------
        
        
        
        
        
        
        !---------------------------------------------------------------------------------------------------
        
        
        !-----------------------------------------------------------------------------         
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        
        
        
        
        
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        ! SUBROUTINE #1
        subroutine MapMassMP2Nodes(MassArray, LumpedMass, ShapeValuesArray, &
                                  ActiveElement)
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
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(ELEMENTNODES) :: PartShape
          real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity) :: LumMassX
          real(REAL_TYPE) :: ParticleMass
          
          integer(INTEGER_TYPE) :: I, IEl, IPart, IAEl, ParticleIndex, iEntity, NoEn

          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: Ni ! number of columns is equal to the number of nodes of the element
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: Mom
          
			
          ! implicit scheme - shape values array inout 
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: LumpedMass
          real(REAL_TYPE), dimension(Counters%NParticles), intent(in) :: MassArray
         
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          
          LumpedMass = 0.0
          !Momentum = 0.0
          NoEn = Counters%nEntity
          
          do IAEl = 1, Counters%NAEl ! loop over all elements

              IEl = ActiveElement_N(IAEl)

              do I = 1, NVECTOR
                  Ni(I, 1:ELEMENTNODES) = ReducedDof( ElementConnectivities(1:ELEMENTNODES, IEl) ) + I
                  !Mom(I, 1:ELEMENTNODES, 1:Counters%nEntity) = Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn)
              end do
              LumMassX = LumpedMass( Ni(1, 1:ELEMENTNODES), 1:NoEn)

              do IPart = 1, NPartEle_N(IEl) ! loop over all material points in element
                  ParticleIndex = GetParticleIndex_N(IPart, IEl) ! Get the particle ID

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
          
        
        
        end subroutine MapMassMP2Nodes
        
        
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        ! SUBROUTINE #2
        
        subroutine MapMomentumMP2Nodes(MassArray, Momentum, VelocityArray, ShapeValuesArray, ElementConnectivities, &
            ActiveElement)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
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
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          real(REAL_TYPE), dimension(Counters%NParticles, NVECTOR), intent(in) :: VelocityArray
          real(REAL_TYPE), dimension(Counters%NParticles), intent(in) :: MassArray
          integer(INTEGER_TYPE), dimension(Counters%NEl, ELEMENTNODES), intent(in) :: ElementConnectivities
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          
          !!CC - changed to function call - needs to be set to zero here
          Momentum = 0.0 
          Time = CalParams%TimeIncrement

          do IAEl = 1, Counters%NAEl   ! Loop over all elements --> NAEl N and NPlus1
            IEl = ActiveElement(IAEl)
            do I = 1, NVECTOR
              IDof(I, 1:ELEMENTNODES) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
            end do
            
            do IPart = 1, NPartEle_N(IEl)   ! Loop over all particles in element
              ParticleIndex = GetParticleIndex_N(IPart, IEl) ! Get the particle ID
              if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 ParticleVelocity = VelocityArray_N(ParticleIndex,:)
                 PartilceMass = MassArray(ParticleIndex)
                 ParticleShape = ShapeValuesArray(ParticleIndex,:)

                  iEntity = 1


                do I = 1, NVECTOR ! nodal i-momentum
                  Momentum(IDof(I,1:ELEMENTNODES), iEntity) = Momentum(IDof(I,1:ELEMENTNODES), iEntity) + PartilceMass * ParticleShape * ParticleVelocity(I)
                end do     

                 
              end if ! NumbOfLayers = 1 or SOLID MatPoint
            end do !Loop over particles
          end do !elements    
          
        end subroutine MapMomentumMP2Nodes  
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #3
        
        subroutine GetNodalVelocityFromNodalMomentum(Momentum, LumpedMassDry, TotalVelocitySoil, PBoundary)!, DoSystem)
        !**********************************************************************
        !
        ! Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum !input
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: LumpedMassDry !input
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: TotalVelocitySoil !input
          
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: PBoundary
          
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
          
        end subroutine GetNodalVelocityFromNodalMomentum
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #4 - BOOKMARK
        subroutine GetNodalExternalForces(ShapeValuesArray, ExtLoad, GravityLoad, ActiveElement)
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
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          integer(INTEGER_TYPE), dimension(counters%NEl), intent(inout) :: ActiveElement
          
          IDim = NVECTOR
           
          call MPMDYNLoad(ExtLoad, GravityLoad, ShapeValuesArray, ActiveElement)

           if ((Counters%NLoadedElementSidesSolidNodes+Counters%NLoadedElementSidesSolidNodesB) > 0) then
             do ILoadSystem =1, Counters%NSoilLoadSystems
                 ExtLoad = ExtLoad + ExtLoadTotal(:,:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
             end do
           end if
           

        end subroutine GetNodalExternalForces
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #5 - BOOKMARK
        subroutine GetNodalIntForces(IntLoad, SigmaEffArray, NPartEle, ActiveElement)
        !**********************************************************************
        !
        !    Function:  get nodal internal forces
        !
        !*********************************************************************     
        implicit none

          integer(INTEGER_TYPE) :: IEntity, I, IDoF, J
          logical :: DoConsiderReactionForces, IsPrescribedVelocity

          ! implicit scheme variables - implicit bookmark
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity), intent(out) :: IntLoad
          
          ! intent in 
          real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(in) :: SigmaEffArray
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          
           !call MPMDYNBTSigOnly(IntLoad, FReaction, BulkViscLoad)
           call MPMDYNBTSigOnly(IntLoad, SigmaEffArray, NPartEle, ActiveElement) !FReaction, BulkViscLoad,


        end subroutine GetNodalIntForces
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #6 - BOOKMARK
        
        subroutine GetFirstGuessOfVelocity(TotalVelocitySoil)
        ! purpose: obtain a first guess for velocity
        !
        implicit none 
        
        real(REAL_TYPE), dimension(Counters%N), intent(inout) :: TotalVelocitySoil
            
        
        TotalVelocitySoil = 0
        
        end subroutine GetFirstGuessOfVelocity
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #7 - BOOKMARK
         subroutine GetMPStrainFromNodeVelocity(TotalVelocitySoil, &
                                               ShapeValuesArray, DShapeValuesArray, &
                                               ActiveElement)
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
         real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
         real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
         integer(INTEGER_TYPE), dimension(Counters%NEl) :: ActiveElement
         
         ! local stresses 
         !real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR) :: SigmaEffArray_Local
              
         EntityUsed = 1   !!CC - each node has nEntity displacements, choose entity used for nodal update
          
          ! SUBROUTINE #1
          call GetNodalIncrementalDisplacement(IncrementalDisplacementSoil, TotalVelocitySoil)

          AccumulatedIncDisplacementSoil = AccumulatedIncDisplacementSoil + IncrementalDisplacementSoil

          ! SUBROUTINE #2: calculate IncrementalDisplacementSoil
          call UpdateNodalTotalDisplacement(EntityUsed, IncrementalDisplacementSoil, TotalDisplacementSoil)
    
          ! SUBROUTINE #3: Update particle strain data
          call UpdateParticleStrains(IncrementalDisplacementSoil, DShapeValuesArray, ActiveElement, NPartEle_N)!, GetParticleIndex_Local) ! update particle strains for solid and mixture
                                                        
          end subroutine GetMPStrainFromNodeVelocity
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #8 - BOOKMARK
        subroutine MPMDYNGetSig(SigmaEffArray, SigmaEff0Array, ActiveElement, &
                                 MaterialIDArray, NPartEle)!(IntLoad_NPlus1) !MPMDYNGetSig ! NAEl, 
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
          real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(inout) :: SigmaEffArray
          real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(in) :: SigmaEff0Array
          
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement         
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          
          !integer(INTEGER_TYPE), dimension(NElements), intent(in) :: NumberOfIntegrationPoints
          
          !integer(INTEGER_TYPE), intent(in) :: NAEl
          integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(in) :: MaterialIDArray

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
         
              IntGlo = GetParticleIndex_N(Int, IEl)   ! Determine global ID of integration point 
              IDMaterial = MaterialIDArray(IntGlo)  ! Material number stored in $$MATERIAL_INDEX in the GOM-file
                          
              if ( (MaterialPointTypeArray(IntGlo)==MaterialPointTypeMixture).or.(MaterialPointTypeArray(IntGlo)==MaterialPointTypeSolid) ) then

                if((NFORMULATION==1).or.(Particles(IntGlo)%PhaseStatus==PhaseStatusSOLID)) then 

                 DoSkipStressComputation = EntityIDArray(IntGlo) == HARD_ENTITY
                 do I = 1, NVECTOR
                   DoSkipStressComputation = DoSkipStressComputation .and. CalParams%ApplyPrescribedVelocity(I)
                 end do
    
                 call StressSolid(IntGlo, IEl, BMatrix, IEntity, MaterialIDArray, SigmaEff0Array, SigmaEffArray,  ESMStatevArray)   ! calculate stresses for mixture or solid material points
                
                end if 
                  
               end if
               
            end do ! Loop over all integration/material points of the element
                        
            ! Assign EUnloading to all particles (material points) of each element
            if (NElemPart==1) then ! only for fully filled elements               
             NumParticles = NPartEle(IEl)
             FirstParticleIndex = GetParticleIndex_N(1, IEl)
             FirstParticleEUnloading = Particles(FirstParticleIndex)%ESM_UnloadingStiffness
             do Int = 2, NumParticles
                IntGlo = GetParticleIndex_N(Int, IEl)
                Particles(IntGlo)%ESM_UnloadingStiffness = FirstParticleEUnloading
             enddo
            endif
            
          end do ! Loop over all active elements for computation of stresses

        end subroutine MPMDYNGetSig
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #9 - BOOKMARK
        subroutine DynUpdateParticleWeights( ShapeValuesArray, ActiveElement, NPartEle )
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
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(in) :: ShapeValuesArray 
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! 1-layer formulation or 2-layer form with 1 Phase
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if((NFORMULATION==1).or. &
               ((NFORMULATION==2).and.(CalParams%NumberOfPhases==1))) then
                                
          do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)   
            
            do Int = 1, NPartEle(IEl) ! Loop over all material points of the element
                
              ! Determine global ID of Material Point 
              IntGlo = GetParticleIndex_N(Int, IEl)
              
              ! we need to keep track of this...
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight

              
              ! Determine volumetric strain of Material Point 
              VolumetricStrain = 0.0
              do I = 1, NVECTOR ! for 2D and 3D all three components have to be taken into account as eps_tt in axisymmetric \=0
                 VolumetricStrain = VolumetricStrain + GetEpsStepI(Particles(IntGlo),I)
              end do
          
             
              ! for SOLID or MIXTURE Material Point               
              !Particles(IntGlo)%IntegrationWeight_Local = Particles(IntGlo)%IntegrationWeight_Local * (1.0 + VolumetricStrain)                
              Particles(IntGlo)%IntegrationWeight = Particles(IntGlo)%IntegrationWeight * (1.0 + VolumetricStrain)
            
             end do !loop over material points
            end do ! loop over elements
            
          end if ! 1-layer formulation or 2-layer form with 1 Phase
          
         
        end subroutine DynUpdateParticleWeights
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #10 - BOOKMARK
        subroutine UpdateParticleHouseKeeping(NPartEle, EleParticlesHelp, IsActiveElement, EleParticles)
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
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(out) :: NPartEle
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(out) :: EleParticlesHelp
          logical, dimension(Counters%NEl), intent(out) :: IsActiveElement
          integer(kind = 8), dimension(Counters%NEl), intent(inout) :: EleParticles

          NPartEle = 0
          call IHpSortKind8(EleParticles, Counters%NParticles)
          ! Count the particles in element, update NPartEle
          do IParticle = 1, Counters%NParticles ! Loop over all particles
            IElement = GetElementIDFromList(IParticle, EleParticles) ! Determine the element that IParticle is located in
            NPartEle(IElement) = NPartEle(IElement) + 1 ! Increase the particle counter for IElement
          end do


          ! Update the helper array and element switches
          FirstParticle = 1
          counters%NAEl = 0 !--> this needs to be an output
          do IElement = 1, Counters%NEl ! Loop over all elements
            if (NPartEle(IElement)==0) then ! Empty element
              EleParticlesHelp(IElement) = -1

              IsActiveElement(IElement) = .false. ! Switch off the element
            else  ! Element contains particles
              counters%NAEl = counters%NAEl + 1 !this needs to be an input
              EleParticlesHelp(IElement) = FirstParticle
              FirstParticle = FirstParticle + NPartEle(IElement)

              IsActiveElement(IElement) = .true. ! Switch on the element
            end if
          end do

        end subroutine UpdateParticleHouseKeeping
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #11 - BOOKMARK
        subroutine SetActiveElement(ActiveElement, ActiveNode)!(ActiveElement, ActiveNode)
        !**********************************************************************
        !
        !  Function : To set the active elements
        !
        !  ActiveElement: Active element vector
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE) :: IAEl, IEl, IError
          
          ! local variables within an iteration 
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: ActiveElement !NAEl or NAEl_Plus1??? !Counters%NAEl
          logical, dimension(Counters%NodTot), intent(inout) :: ActiveNode
        
          if (allocated(ActiveElement) ) then
            deallocate(ActiveElement, stat = IError)
          end if
          
          allocate(ActiveElement(counters%NAEl), stat = IError)
          ActiveElement = 0

          IAEl = 0
          do IEl=1, counters%NEl
            if (IsActiveElement_N(IEl)) then ! Active elements
              IAEl = IAEl + 1
              ActiveElement(IAEl) = IEl
            end if
          end do
 
          call DetermineActiveNodes(ActiveElement, ActiveNode)
          
        end subroutine SetActiveElement
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #12 - BOOKMARK
        
        subroutine SetParticleIndex(ActiveElement,NPartEle, &
            EleParticlesHelp, EleParticles)!,GetParticleIndex)
        !**********************************************************************
        !
        !    Function:  Sets/updates GetParticleIndex.
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        integer(INTEGER_TYPE) :: maxParticle, iError, IAEl, iEl, int, maxEl

        integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
        integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
        !integer(INTEGER_TYPE), dimension(maxParticle,NElements), intent(in) :: GetParticleIndex !maxParticle is not defined...
        
        
        ! Input variables 
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: EleParticlesHelp
          integer(kind = 8), dimension(Counters%NParticles), intent(in) :: EleParticles
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle

        
        maxParticle = 1
        maxEl       = 1
        do IAEl = 1, Counters%NAEl
          iEl = ActiveElement(IAEl)
          maxEl = max(maxEl, iEl)
          maxParticle = max(maxParticle, NPartEle(iEl))
        enddo

        iError = 0
        if (allocated(GetParticleIndex_N)) then
          deallocate(GetParticleIndex_N, stat = IError)
          call DeAllocationError(iError, 'GetParticleIndex', 'SetParticleIndex')
        endif

        allocate(GetParticleIndex_N(maxParticle, maxEl), stat = IError)
        call AllocationError(iError, 'GetParticleIndex', 'SetParticleIndex')
        GetParticleIndex_N = -huge(GetParticleIndex_N)

        do IAEl = 1, Counters%NAEl
          iEl = ActiveElement(IAEl)
          do int = 1, NPartEle(iEl)
            GetParticleIndex_N(int, iEl) = GetParticleIndexFunction(int, iEl)!, &
                !EleParticlesHelp, EleParticles, NPartEle)
          enddo
        enddo
        
        end subroutine SetParticleIndex
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #13 - BOOKMARK
        
        subroutine CheckFillingOfElements(ActiveElement, MaterialIDArray, NPartEle)
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
          
          ! implicit scheme vectorization of subroutine
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(in) :: MaterialIDArray
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle

          if (.not.IsMPMWithMixedIntegration()) RETURN
            IsParticleIntegration = .false.

            do IAElement = 1, Counters%NAEl
             IElement = ActiveElement(IAElement)
              ContainsBoundaryParticle = .false.
              SummedParticleWeights = 0.0

              ParticleIndex = GetParticleIndex_N(1, IElement)
              CheckMaterialID = MaterialIDArray(ParticleIndex)
              IsMixedElement = .false.

              do IParticle = 1, NPartEle(IElement) ! Loop over particles in IElement
                ParticleIndex = GetParticleIndex_N(IParticle, IElement)
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
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #14 - BOOKMARK
        
        
        subroutine SetUpEntityElements(ActiveElement, NPartEle)
        !**********************************************************************
        !
        !  Function : For each entity a list of the elements containing particles belonging to this entity
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: IAEl, IEl, NElemPart, IPart, ParticleIndex, IEntity
          
          ! implicit scheme vectorization of subroutine
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          
          
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle_NPlus1
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement_NPlus1
          
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement_N
          !ActiveElement_NPlus1
          
          EntityElements = 0      !reset

          do IAEl = 1,Counters%NAel                                      !loop through all elements
            IEl = ActiveElement(IAEl)
            NElemPart = NPartEle(IEl)
            do IPart= 1, NElemPart                                 !loop through particles in element
              ParticleIndex = GetParticleIndex_N(IPart, IEl)         !Get the particle ID
              IEntity = EntityIDArray(ParticleIndex)               !entity to which particle belong
              EntityElements(IEntity,IEl) = 1
            end do
          end do
         
        end subroutine SetUpEntityElements
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #15 - BOOKMARK
        subroutine SetUpMaterialElements(ActiveElement, MaterialIDArray, NPartEle)
        !**********************************************************************
        !
        !  Function : For each material ID a list of the elements containing particles belonging to this material ID
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: IAEl, IEl, NElemPart, IPart, ParticleIndex, MaterialID
          
          ! implicit scheme vectorization of subroutine
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(in) :: MaterialIDArray
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle

         
          MaterialElements = 0      !reset
          
                do IAEl = 1,Counters%NAel                                      !loop through all elements
            IEl = ActiveElement(IAEl)
            NElemPart = NPartEle(IEl)
            do IPart= 1, NElemPart                                 !loop through particles in element
              ParticleIndex = GetParticleIndex_N(IPart, IEl)        ! Get the particle ID
              MaterialID = MaterialIDArray(ParticleIndex)       !material ID to which particle belong
              MaterialElements(MaterialID,IEl) = 1
            end do
          end do
                  
        end subroutine SetUpMaterialElements
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #16 - BOOKMARK
        
        subroutine SetInitialStressForNextLoadStep(ActiveElement, NPartEle, SigmaEffArray, SigmaEff0Array)
        !**********************************************************************
        !
        !    Function:  Call the routines in which the initial stresses of particles
        !               and Gauss points are stored
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
        
        ! implicit scheme variables vectorization 
          real(REAL_TYPE), dimension(NTENSOR), intent(in) :: SigmaEffArray !new stress computed in the convective phase
          real(REAL_TYPE), dimension(NTENSOR), intent(inout) :: SigmaEff0Array ! old stres from previous step
                    
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement         
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
        

          ! Store intial particle stresses of next step
          call SetInitialParticleStressForNextLoadStep(ActiveElement, NPartEle, SigmaEffArray, SigmaEff0Array)

        end subroutine SetInitialStressForNextLoadStep
        
        
        
        
        
        
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #17 - BOOKMARK
        
        subroutine CalculateIntAndExtWorks(IncrementalDisplacementWater, IntLoad, IntLoad_NPlus1, SigmaEffArray)
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
          
          ! implicit scheme variable - implicit bookmark 
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(out) :: IntLoad_NPlus1
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: SigmaEffArray
          
          if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then ! if consolidation

          else ! 1-phase calculation
          
            !IntLoadPrevious = IntLoad ! I commented this on 12 September 2023
            call GetNodalIntForces(IntLoad_NPlus1, SigmaEffArray, NPartEle_N, ActiveElement_N)  
              
            IncrementalIntWorkSoil = 0.0
            IncrementalExtWorkSoil = 0.0
            ExternalLoadDueToVelocity = 0.0
          
          end if
          
        end subroutine CalculateIntAndExtWorks
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #18 - BOOKMARK
        
        subroutine GetImplicitResidualEquation_SK04(LumpedMassDry_NPlus1, LumpedMassDry, &  
                                                            ExtLoad_NPlus1, ExtLoad, & 
                                                        GravityLoad_NPlus1, GravityLoad, &  
                                                            IntLoad_NPlus1, IntLoad, & 
                                                  TotalVelocitySoil_NPlus1, TotalVelocitySoil, &
                                                    ImplicitGaugeParameter, &
                                                          ImplicitResidual)
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
          
          real(REAL_TYPE), intent(in) :: ImplicitGaugeParameter
          
          
          ! local terms 
          integer(INTEGER_TYPE) :: IDof
          integer(INTEGER_TYPE) :: iEntity
		  
          !iEntity_Temporary =1
          
         do IEntity = 1, Counters%nEntity ! IEntity loop
         
             do IDof = 1, Counters%N ! IDof loop 
             
                 ImplicitResidual(IDof, IEntity) = (LumpedMassDry(IDof,IEntity) * TotalVelocitySoil_NPlus1(IDof,IEntity)) + &
                                               - (LumpedMassDry(IDof,IEntity) * TotalVelocitySoil(IDof,IEntity)) + &
                                               - (CalParams%TimeIncrement * ImplicitGaugeParameter * &
                                                   (ExtLoad_NPlus1(IDof,IEntity) + GravityLoad_NPlus1(IDof,IEntity) - IntLoad_NPlus1(IDof,IEntity)) ) + &
                                               - (CalParams%TimeIncrement * (1-ImplicitGaugeParameter) * &
                                                   (ExtLoad(IDof,IEntity) + GravityLoad(IDof,IEntity) - IntLoad(IDof,IEntity)) )
         
             end do ! IDof loop 
         end do ! IEntity loop 
         
         
        end subroutine GetImplicitResidualEquation_SK04  
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #19 - BOOKMARK
        
        
        !%---------------------------------------------------------------------%
        !%                                GMRES                                %
        !%---------------------------------------------------------------------%
        !---this is the mother subroutine where we use the directional derivative in the GMRES algorithm 
        subroutine GMRES_SK04( & !AA, bb, xx, threshold, nn_length_AA, mm_max_iteration)
                      LumpedMassDry, & !LumpedMassDry_NPlus1, & ! in variables 
                      TotalVelocitySoil, &
                      IntLoad, & !IntLoad_NPlus1, & ! in variables 
                      ImplicitResidual, & ! in variables 
                      Threshold, & ! in variable
                      !nn_Length_AA, & ! in variable
                      mm_max_iteration, &
                      ImplicitGaugeParameter)!, & ! this is s_0 in SK04
                      !AccelerationSoil, & !Explicit velocity
                      !IncrementalDisplacementSoil, &
                      !AccumulatedIncDisplacementSoil, &
                      !TotalDisplacementSoil, &
                      !AccelerationSoil_Incremental, & ! in variable

          implicit none 
          ! -- Input -------------------------------------------------------------------------
          !-----------Mass
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: LumpedMassDry_NPlus1!, intent(in)
          !-----------Velocity
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil
          !-----------Internal load
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity), intent(in) :: IntLoad
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: IntLoad_NPlus1 !, intent(in)
          !-----------External load
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: ExtLoad !, intent(in)
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: ExtLoad_NPlus1 !, intent(in)
          !-----------Gravity load
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: GravityLoad !, intent(in)
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: GravityLoad_NPlus1 !, intent(in)
          !-----------Implicit residual
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity), intent(in) :: ImplicitResidual
          real(REAL_TYPE) :: ImplicitResidualMagnitude !,                                         intent(in)
          ! -----------------------------------------------------------------------------------
          !real(REAL_TYPE), dimension(Counters%N,Counters%NEntity), intent(in) :: 
          ! -- Output -------------------------------------------------------------------------
          ! -----------------------------------------------------------------------------------
          ! -- Inout --------------------------------------------------------------------------
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_NPlus1 !, intent(inout)
          real(REAL_TYPE), dimension(Counters%N) :: DeltaTotalVelocitySoil !,Counters%nEntity
          real(REAL_TYPE) :: DeltaTotalVelocitySoilMagnitude
          real(REAL_TYPE) :: TotalVelocitySoilMagnitude
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: LumpedMassDry_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: ExtLoad_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: GravityLoad_DirDeriv
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: ImplicitResidual_DirDeriv
          real(REAL_TYPE), dimension(Counters%N) :: ProdJacobDeltaVeloDirDeriv

          ! -----------------------------------------------------------------------------------
          ! -- Local --------------------------------------------------------------------------
          ! -----------------------------------------------------------------------------------
          
          
          
          real(REAL_TYPE), intent(in) :: Threshold
          !integer(INTEGER_TYPE), intent(in) :: nn_length_AA
          !---inout variables
          ! should this be inout or just output??
          integer(INTEGER_TYPE), intent(inout) :: mm_max_iteration 
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
          real(REAL_TYPE) :: ee
          integer(INTEGER_TYPE) :: INFO
          
          real(REAL_TYPE) :: h_DirDerivParameter_SK04
          real(REAL_TYPE), intent(in) :: ImplicitGaugeParameter
        
          !---the residual is already evaluated so we do not need to evaluate this again here 
          !---using norm2 function from FOTRAN90. Inputing rank as 1, but this is not necessary (i.e., optional)
          
          !do ii = 1, Counters%N 
          !    ImplicitResidual(ii,1) = ImplicitResidual(ii) * PBoundary(ii,1) !--> applying boundary conditions
          !end do 
          mm_max_iteration = 10 ! # of iterattions
          h_DirDerivParameter_SK04 = 1e-5
          
          !ImplicitResidualMagnitude = norm2(ImplicitResidual, 1) !rr_norm
          
          !---we define the error so many times. Here is it a function of the magnitude but we have to settle on a particular option. 
          error = ImplicitResidualMagnitude!/ !bb_norm ! real scalar in the numerator and denominator
          
          !---initialize the 1D vectors
          allocate(sn(mm_max_iteration))
          allocate(cs(mm_max_iteration))
          allocate(e1(mm_max_iteration+1))
          allocate(beta(mm_max_iteration+1))
          allocate(QQ(Counters%N,Counters%nEntity))
          !allocate(QQ_Full(nn,nn))
          
          cs = 0 
          sn = 0 
          e1 = 0 
          beta = 0
          QQ = 0 
          
          e1(1) = 1 ! assigning the first element of ee to be 1
          ee = error ! ee is a scalar. Set it equal to the value of error
          
          !---this would be p_1 in SK04 Algorithm 4 step 1
          QQ(:,1) = (1/ImplicitResidualMagnitude) * ImplicitResidual(:,1) ! assigning this for the first iteration 
          
          beta = ImplicitResidualMagnitude * e1 ! e1 and beta are the same size vectors 
          
          nn = 1
          kk = 1
          HH_Arnoldi = 0 
          QQ_Arnoldi = 0
          HH_Arnoldi_New = 0 
          !QQ_Arnoldi_New = 0
          
          allocate(QQ_Stored(Counters%N,kk))  !kk+1
          allocate(HH_Stored(kk+1,kk))
          
          QQ_Stored = 0
          HH_Stored = 0
          
          do kk = 1, mm_max_iteration ! loop accross the number of maximum iterations
              
              ! if allocated we need to deallocate here
              if (allocated(QQ_Full)) then 
                  deallocate(QQ_Full)
              end if 
              
              ! if allocated we need to deallocate here
              if (allocated(HH_Full)) then 
                  deallocate(HH_Full)
              end if
              
              allocate(QQ_Full(Counters%N,kk+1)) !2
              allocate(HH_Full(kk+1,kk))
              
              
              QQ_Full = 0
              HH_Full = 0
              !QQ_Full(:,1) = QQ(:,1)
              QQ_Stored(:,1) = QQ(:,1) ! this is a bit redundant in steps above k=1
              QQ_Full(:,1:kk) = QQ_Stored(:,1:kk) ! QQ_Full in column kk+1 should be left empty
              
              ! if larger than 2 then we need to store HH
              if (kk .ge. 2) then 
                  HH_Full(1:kk,1:kk-1) = HH_Stored(:,:)
              end if
              
              ! SUBROUTINE #1:
              call GetDirDerivVelocity_SK04(TotalVelocitySoil, TotalVelocitySoilMagnitude, & !input 
                                                  DeltaTotalVelocitySoil, DeltaTotalVelocitySoilMagnitude, & !input
                                                  h_DirDerivParameter_SK04, & !input
                                                  TotalVelocitySoil_DirDeriv) !output
              
              ! SUBROUTINE #2:
              call GetMPStrainFromNodeVelocity(TotalVelocitySoil_N, ShapeValuesArray_N, DShapeValuesArray_N, ActiveElement_N)                                           
              
              ! SUBROUTINE #3:
              call MPMDYNGetSig(SigmaEffArray_N, SigmaEff0Array, ActiveElement_N, &!NumberOfIntegrationPoints, &
                                MaterialIDArray_N, NPartEle_N) !NAEl
              
              ! SUBROUTINE #4:
              call DynUpdateParticleWeights( ShapeValuesArray_N, ActiveElement_N, NPartEle_N )
              
              ! SUBROUTINE #5:
              call UpdateParticleHouseKeeping(NPartEle_N, EleParticlesHelp_N, IsActiveElement_N, EleParticles_N)
              
              ! SUBROUTINE #6:
              call SetActiveElement(ActiveElement_N, ActiveNode_N)
              
              ! SUBROUTINE #7: 
              call SetParticleIndex(ActiveElement_N,NPartEle_N, EleParticlesHelp_N, EleParticles_N)
              
              ! SUBROUTINE #8:
              call CheckFillingOfElements(ActiveElement_N, MaterialIDArray_N, NPartEle_N)     !NPartEle, EleParticlesHelp, IsActiveElement)
              
              ! SUBROUTINE #9:
              call SetUpEntityElements(ActiveElement_N, NPartEle_N)
              
              ! SUBROUTINE #10:
              call SetUpMaterialElements (ActiveElement_N, MaterialIDArray_N, NPartEle_N) !--> why are we even doing this???
              
              ! SUBROUTINE #11:
              !call SetInitialStressForNextLoadStep(ActiveElement, NPartEle, SigmaEffArray, SigmaEff0Array)
              
              ! SUBROUTINE #12:
              call CalculateIntAndExtWorks(IncrementalDisplacementWater, IntLoad, IntLoad_NPlus1, SigmaEffArray_N)
              
              ! SUBROUTINE #13:
              call GetNodalExternalForces(ShapeValuesArray_NPlus1, ExtLoad_NPlus1, GravityLoad_NPlus1, ActiveElement_N)
              
              ! SUBROUTINE #14
              call MapMassMP2Nodes(MassArray, LumpedMassDry_NPlus1, ShapeValuesArray_NPlus1, ActiveElement_N)
              
              ! SUBROUTINE #15
              call GetImplicitResidualEquation_SK04(LumpedMassDry_DirDeriv, LumpedMassDry, & !input 
                                                          ExtLoad_DirDeriv, ExtLoad, &  !input
                                                      GravityLoad_DirDeriv, GravityLoad, &  !input
                                                          IntLoad_DirDeriv, IntLoad, &  !input
                                                TotalVelocitySoil_DirDeriv, TotalVelocitySoil, &  !input
                                                    ImplicitGaugeParameter, &  !input
                                                 ImplicitResidual_DirDeriv)  !output
              
              ! SUBROUTINE #: 
              !call GetDirDerivVelocity_SK04(TotalVelocitySoil, TotalVelocitySoilMagnitude, & !input 
              !                                    DeltaTotalVelocitySoil, DeltaTotalVelocitySoilMagnitude, & !input
              !                                    HH, & !input
              !                                    TotalVelocitySoil_DirDeriv) !output
              
              
              ! SUBROUTINE #16
              DeltaTotalVelocitySoil = TotalVelocitySoil_DirDeriv(:,1) - TotalVelocitySoil(:,1)
              
              ! SUBROUTINE #17: finding the magnitudes for directional derivative evaluation
              DeltaTotalVelocitySoilMagnitude = norm2(DeltaTotalVelocitySoil) ! -> s in SK04 
                                                ! 1 is for looking at only the first entity
              TotalVelocitySoilMagnitude = norm2(TotalVelocitySoil) ! v in SK04
                                                ! 1 is for looking at only the first entity
              
              ! SUBROUTINE #18: 
              call GetProdJacobDeltaVeloDirDeriv(ProdJacobDeltaVeloDirDeriv, DeltaTotalVelocitySoilMagnitude, &
                                                       TotalVelocitySoil, TotalVelocitySoilMagnitude, &
                                                        ImplicitResidual_DirDeriv, ImplicitResidual, &
                                                        h_DirDerivParameter_SK04)
              
              ! SUBROUTINE #19: p_k+1 = DirectionalDerivative(AccelerationSoil_NPlus1, OldImplicitResidual))  
              QQ_Arnoldi(:,1) = ProdJacobDeltaVeloDirDeriv!DirectionalDerivative
              
              ! SUBROUTINE #20: THIS SHOULD BE QQ(:,kk) --> It should depend on the timestep 
              !---step 2b h_k,i = p_k+1 * p_i 
              do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix            
                  do jj = 1,Counters%N !2 
                      HH_Arnoldi(ii) = HH_Arnoldi(ii) + ( QQ_Arnoldi(jj,1) * QQ_Full(jj,ii) ) !QQ_Arnoldi is the p_k+1 term. QQ is the p_i term 
                  end do 
              end do 
              
              ! SUBROUTINE #21: !---step 2b in algorithm 4 in SK04
              !----Equation p_k+1 = p_k+1 - (h_k,i * p_i)
              do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix
                  QQ_Arnoldi(:,1) = QQ_Arnoldi(:,1) - ( HH_Arnoldi(ii) * QQ_Full(:,ii)) !QQ_Arnoldi is the p_k+1 term. QQ is the p_i term 
              end do 
              
              ! SUBROUTINE #22: this needs to be stored outside this loop within the entire module as a global variable 
              if (allocated(HH_Arnoldi_New)) then 
                  deallocate(HH_Arnoldi_New)
              end if 
              allocate(HH_Arnoldi_New (kk+1, 1) ) 
              HH_Arnoldi_New = 0
              HH_Arnoldi_New(1:kk, 1) = HH_Arnoldi(1:kk) 
              
              ! SUBROUTINE #23: 
              !---step 1c h_k+1,k = |p_k+1|
              HH_Arnoldi_New(kk+1, 1) =  norm2(QQ_Arnoldi(:,1),1)
        
              !---step 1d p_k+1 = p_k+1/|p_k+1| 
              QQ_Arnoldi = QQ_Arnoldi * (1/HH_Arnoldi_New(kk+1, 1) ) ! divide by a scalar

              !---pouring q (or called as p_k+1 in SK04) back to QQ_Full to store    
              QQ_Full(:,kk+1) = QQ_Arnoldi(:,1)
              
              ! SUBROUTINE #24: eliminate the last element in H ith row and update the rotation matrix 
              call apply_givens_rotation( HH_Arnoldi_New, cs, sn, kk, mm_max_iteration )
              
              ! SUBROUTINE #25: 
              if (allocated(HH_Arnoldi_Full)) then 
                  deallocate(HH_Arnoldi_Full)
              end if 
              allocate(HH_Arnoldi_Full(nn,kk))
              HH_Arnoldi_Full = 0
              HH_Arnoldi_Full = HH_Arnoldi_New
              
              ! SUBROUTINE #26: update the residual vector 
              beta(kk+1) = -sn(kk) * beta(kk)
              beta(kk) = cs(kk) * beta(kk)
              
              ! SUBROUTINE #27: evaluate error
              error = abs(beta(kk+1)) !/ bb_norm
              
              ! SUBROUTINE #28: allocate/deallocate QQ and HH
              deallocate(QQ_Stored)  !kk+1
              deallocate(HH_Stored)
              allocate(QQ_Stored(Counters%N,kk+1))  !kk+1
              allocate(HH_Stored(kk+1,kk))
              QQ_Stored = QQ_Full
              
              ! SUBROUTINE #29: 
              if (kk .ge. 2) then 
                  HH_Stored = 0
                  HH_Stored(:, 1:kk-1) = HH_Full(:, 1:kk-1)
                  HH_Stored(:, kk) = HH_Arnoldi_New(:,1)
              else 
                  HH_Stored = 0          
                  HH_Stored = HH_Arnoldi_New
              end if 
              
              ! SUBROUTINE #30: check error against threshold
              if (error .le. threshold ) then 
                  exit
              end if 

          end do ! loop over mm maximum iterations
          
          ! SUBROUTINE #31: allocate HH based on the number of iterations kk
          allocate(HH_inverted(kk,kk))
          HH_inverted = 0.0
          
          ! SUBROUTINE #32: find the inverse of upper triangular H matrix
          call L_inv(HH_Stored(1:kk, 1:kk), & ! in variable, 
                     HH_inverted, & ! out variable 
                     kk)  ! in variable  
          
          ! SUBROUTINE #33: GMRES Step 3 to solve the system Hy=beta for y
          allocate(YY(kk)) 
          YY = 0 
          
          ! SUBROUTINE #34: conduct the multiplication to find YY 
          ! loop over rows  
          do ii = 1, kk ! Modified Gram-Schmidt, keeping the Hessenberg matrix
              ! loop over columns 
              do jj = 1, kk 
                  YY(ii) = YY(ii) + HH_Stored(ii, jj) * beta(jj)  
              end do 
          end do 
          
          ! SUBROUTINE #35: GMRES Step 4 to Calculate the solution s = s_0 + y1*p_1 + ... + y_k*p_k
          ! xx = xx + (QQ_Full(:,1:kk) * YY)
          ! loop over rows 
          do ii = 1, Counters%N ! rank of matrix 
              do jj = 1, kk
                  DeltaTotalVelocitySoil(ii) = DeltaTotalVelocitySoil(ii) + QQ_Full(ii,jj)*YY(jj)
              end do
          end do
          
          !--------------------------------------------------------------------------------------
          !--------------------------------------------------------------------------------------
          
           end subroutine GMRES_SK04
        
        
        
                                                          
                                                          
                                                          
                                                          
                                                          
                                                          
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #20 - BOOKMARK
                      
                      
                      
        subroutine GetProdJacobDeltaVeloDirDeriv(ProdJacobDeltaVeloDirDeriv, DeltaTotalVelocitySoilMagnitude, &
                                                       TotalVelocitySoil, TotalVelocitySoilMagnitude, &
                                                        ImplicitResidual_DirDeriv, ImplicitResidual, &
                                                        h_DirDerivParameter_SK04)
              ! purpose: to evaluate ProdJacobDeltaVeloDirDeriv
              !
              implicit none 
              
              
              !---Input
              real(REAL_TYPE), intent(in) :: DeltaTotalVelocitySoilMagnitude
              real(REAL_TYPE), intent(in) :: TotalVelocitySoilMagnitude
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: ImplicitResidual_DirDeriv
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: ImplicitResidual
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: TotalVelocitySoil
              real(REAL_TYPE), intent(in) :: h_DirDerivParameter_SK04
              !---Output
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: ProdJacobDeltaVeloDirDeriv
              !---Local
              real(REAL_TYPE) :: VerySmallNumber
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: Zeros_Vector !, intent(in)
              integer(INTEGER_TYPE) :: ii
              
              ! setting up local variables 
              VerySmallNumber = 1e-3
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
                          (1/h_DirDerivParameter_SK04) * (1/TotalVelocitySoilMagnitude)  
                  end do
                  
              !---case 3
              elseif ( (DeltaTotalVelocitySoilMagnitude .gt. VerySmallNumber) .and. (TotalVelocitySoilMagnitude .lt. VerySmallNumber) ) then 
              
                  !---Directional derivative calculation
                  !----the reason there is a 1 is because the directional derivative is a vector (rank 1). But Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside is rank 2.  
                  do ii = 1,Counters%N
                      ProdJacobDeltaVeloDirDeriv(ii,1) = DeltaTotalVelocitySoilMagnitude * &
                          (ImplicitResidual_DirDeriv(ii,1) - ImplicitResidual(ii,1)) * & ! I removed generalizedalpha naming from here...
                          (1/h_DirDerivParameter_SK04) !* PBoundary(ii)
                  end do
              end if 
              
              end subroutine GetProdJacobDeltaVeloDirDeriv                   
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #21 - BOOKMARK
        
        
        subroutine GetDirDerivVelocity_SK04(TotalVelocitySoil, TotalVelocitySoilMagnitude, &
                                                  DeltaTotalVelocitySoil, DeltaTotalVelocitySoilMagnitude, &
                                                  h_DirDerivParameter_SK04, &
                                                  TotalVelocitySoil_DirDeriv)
              ! purpose: get the velocity to be used in the directional derivative equation 
              implicit none 
              !---inputs
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: TotalVelocitySoil
              real(REAL_TYPE), intent(in) :: TotalVelocitySoilMagnitude
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: DeltaTotalVelocitySoil
              real(REAL_TYPE), intent(in) :: DeltaTotalVelocitySoilMagnitude
              real(REAL_TYPE), intent(in) :: h_DirDerivParameter_SK04
              !---outputs
              real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: TotalVelocitySoil_DirDeriv
              !---local
              real(REAL_TYPE) :: VerySmallNumber
              
              VerySmallNumber = 1e-3
              
              ! SUBROUTINE #18: 
              !---case 1
              if (DeltaTotalVelocitySoilMagnitude .lt. VerySmallNumber) then ! equals to zero     
                  TotalVelocitySoil_DirDeriv = 0
                  
              !---case 2    
              elseif ( (DeltaTotalVelocitySoilMagnitude .gt. VerySmallNumber) .and. &
              (TotalVelocitySoilMagnitude .gt. VerySmallNumber) ) then ! greater than zero 
                  TotalVelocitySoil_DirDeriv = TotalVelocitySoil + &
                      ( h_DirDerivParameter_SK04 * TotalVelocitySoilMagnitude * DeltaTotalVelocitySoil * (1/DeltaTotalVelocitySoilMagnitude) )
                  
              !---case 3
              elseif ( (DeltaTotalVelocitySoilMagnitude .gt. VerySmallNumber) .and. &
              (TotalVelocitySoilMagnitude .lt. VerySmallNumber) ) then 
                  TotalVelocitySoil_DirDeriv = 0 + & ! this should be TotalVelocitySoil but since it is small = approx zero
                      ( h_DirDerivParameter_SK04 * DeltaTotalVelocitySoil * (1/DeltaTotalVelocitySoilMagnitude) )
                  
              end if 
              
              
           end subroutine GetDirDerivVelocity_SK04
        
        
        
        
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
              
              
        
        ! apply for ith column 
        do ii = 1,kk-1
            temp          = (  cs(ii) * HH_Arnoldi(ii) ) + ( sn(ii) * HH_Arnoldi(ii+1) )
            HH_Arnoldi(ii+1)  = ( -sn(ii) * HH_Arnoldi(ii) ) + ( cs(ii) * HH_Arnoldi(ii+1) )
            HH_Arnoldi(ii)    = temp    
        end do 
        
        ! update the next sin cos values for rotation 
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
        ! SUBROUTINE #25 - BOOKMARK
        
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subroutine EvaluateResidualMagnitude(Residual_Implicit, Residual_Implicit_Magnitude) 
         ! Part of Step 1 in Algorithm 4
         
         implicit none 
         
         !output 
         real(REAL_TYPE), intent(out) :: Residual_Implicit_Magnitude
         !input 
         real(REAL_TYPE), intent(in), dimension(Counters%N) :: Residual_Implicit
         !local
         real(REAL_TYPE) :: Residual_Implicit_Magnitude_Squared
         integer(INTEGER_TYPE) :: ii 
         
         Residual_Implicit_Magnitude_Squared = 0
         
         do ii = 1, Counters%N

             Residual_Implicit_Magnitude_Squared = Residual_Implicit_Magnitude_Squared + (Residual_Implicit(ii))**2
             
             !(Residual_Implicit(1)**2 + Residual_Implicit(2)**2)**0.5
         
         end do 
         
         
         Residual_Implicit_Magnitude = Residual_Implicit_Magnitude_Squared**0.5
         
         
        end subroutine            
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #26 - BOOKMARK
        
        subroutine InitialiseMaterialPointHousekeeping(&
            IsActiveElement_N, &            
            IsActiveElement_NPlus1, &
            NSolidEle_N, &
            NSolidEle_NPlus1, &
            NLiquidEle_N, &
            NLiquidEle_NPlus1, &
            ElementIDArray_N, &
            ElementIDArray_NPlus1, &
            ShapeValuesArray_N, &
            ShapeValuesArray_NPlus1, &
            DShapeValuesArray_N, &
            DShapeValuesArray_NPlus1, &
            MaterialIDArray_N, &
            MaterialIDArray_NPlus1, &
            EleParticles_N, &
            EleParticles_NPlus1, &
            EleParticlesHelp_N, &
            EleParticlesHelp_NPlus1)
        
        
        !(ShapeValuesArray, DShapeValuesArray)
        !**********************************************************************
        !
        !  Function : Initialises the list of material points
        !
        !             Material point information is initialised in case of the calculation of an initial phase.
        !             If a 'follow-up' phase is calculated, material point information is reinitialised
        !             from data from the BRF file of the previous load phase. 
        !
        !**********************************************************************

        implicit none
        
        !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(inout) :: ShapeValuesArray
        !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NDIM), intent(inout) :: DShapeValuesArray

        ! implicit scheme variables - implicit bookmark
          logical, dimension(Counters%NEl) :: IsActiveElement_N !integer(INTEGER_TYPE)
          logical, dimension(Counters%NEl) :: IsActiveElement_NPlus1 !integer(INTEGER_TYPE)
          
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: NSolidEle_N
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: NSolidEle_NPlus1
          
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: NLiquidEle_N
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: NLiquidEle_NPlus1

          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: ElementIDArray_N !dimension(Counters%NParticles)
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: ElementIDArray_NPlus1
          
          real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: ShapeValuesArray_N !dimension(Counters%NParticles, ELEMENTNODES),
          real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: ShapeValuesArray_NPlus1 !dimension(Counters%NParticles, ELEMENTNODES),
          
          real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: DShapeValuesArray_N !dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), 
          real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: DShapeValuesArray_NPlus1 !dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), 
          
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: MaterialIDArray_N
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: MaterialIDArray_NPlus1
          
          integer(kind = 8), allocatable, dimension(:), intent(inout) :: EleParticles_N
          integer(kind = 8), allocatable, dimension(:), intent(inout) :: EleParticles_NPlus1
          
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: EleParticlesHelp_N
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: EleParticlesHelp_NPlus1
          
          
          
           !I COMMENTED THIS BUT I NEED TO UNCOMMENT THIS!
          call DestroyHouseKeeping(&
              NPartEle_N, &
              NPartEle_NPlus1, &
              EleParticlesHelp_N, &
              EleParticlesHelp_NPlus1, &
              EleParticles_N, &
              EleParticles_NPlus1, &
              ElementIDArray_N, &
              ElementIDArray_NPlus1, &
              MaterialIDArray_N, &
              MaterialIDArray_NPlus1, &
              ShapeValuesArray_N, &
              ShapeValuesArray_NPlus1, &
              DShapeValuesArray_N, &
              DShapeValuesArray_NPlus1, &
              GlobPosArray_N, &
              GlobPosArray_NPlus1, &
              SigmaEffArray_N, &
              SigmaEffArray_NPlus1)

          if (.not.IsFollowUpPhase()) then ! Initial phase: initialise material point data

              ! I NEED TO UNCOMMENT THIS ! 
            !call InitialiseMaterialPoints(&
            !IsActiveElement_N, NSolidEle, NLiquidEle, &
            !IsActiveElement_NPlus1, &
            !ElementIDArray_N, ElementIDArray_NPlus1, &
            !ShapeValuesArray_N, DShapeValuesArray_N, & !-> this where we define the shape functions
            !ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
            !MaterialIDArray_N, MaterialIDArray_NPlus1, &
            !EleParticles_N, EleParticlesHelp_N, &
            !EleParticles_NPlus1, EleParticlesHelp_NPlus1) ! Initialise material point data structure

          else ! Follow-up phase: reinitialise material point data from the BRF file of the previous load phase

            call ReinitialiseMaterialPoints(ShapeValuesArray_N, DShapeValuesArray_N) ! Reinitialise material point data structure
            CalParams%ApplyPrescribedDisplacements = ReadPPD() ! Initialise prescribed displacements if the PPD file exists

          end if

          call SetActiveElement(ActiveElement_N, ActiveNode_N) ! allocate and assign array ActiveElement
          call SetParticleIndex(ActiveElement_N,NPartEle_N, &
            EleParticlesHelp_N, EleParticles_N)

        end subroutine InitialiseMaterialPointHousekeeping
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #27 - BOOKMARK
        
        subroutine AdjustParticleHouseKeeping()
        !**********************************************************************
        !
        !    Function:  After the number of particles has been changed, the
        !               remaining particle house keeping data structure is updated.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
          
          ! Local variables
          integer(kind = 8) :: ParticleIndex,ElementID

          ! Update EleParticles array
          do ParticleIndex = 1, Counters%NParticles ! Loop over particles
             ElementID = ElementIDArray_N(ParticleIndex)
             
            EleParticles_N(ParticleIndex) =  &
              ElementID * GetPartFactor() + ParticleIndex
          end do

          ! Sort EleParticles array, update EleParticlesHelp and NPartEle arrays
          call UpdateParticleHouseKeeping(NPartEle_N, EleParticlesHelp_N, IsActiveElement_N, EleParticles_N)

        end subroutine AdjustParticleHouseKeeping
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #28 - BOOKMARK
        
        subroutine ReinitialiseMaterialPoints(ShapeValuesArray, DShapeValuesArray)
        !**********************************************************************
        !
        !    Function:  Reinitialises the list of particles and the house-keeping arrays
        !               from particle data from the BRF file of the previous load phase.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IError
          character(len = MAX_FILENAME_LENGTH) :: FileName
          
          ! implicit scheme - shape values array inout 
         real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(inout) :: ShapeValuesArray 
         real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(inout) :: DShapeValuesArray
		 
        
          FileName = trim(CalParams%FileNames%ProjectName)//BRF_FILE_EXTENSION//trim(CalParams%FileNames%PreviousStepExt)

          ! Read number of particles from the BRF file of the previous load phase
          call ReadHouseKeepingCountersFromFile(FileName)

          ! Initialise particle house-keeping arrays with the counters from the BRF file of the previous load phase
          ! I COMMENTED THIS BUT I NEED TO UNCOMMENT THIS!
          !call InitialiseHouseKeepingArrays(&
          !  EleParticlesHelp_N, EleParticles_N, &
          !  EleParticlesHelp_NPlus1, EleParticles_NPlus1, &
          !  ElementIDArray_N, MaterialIDArray_N, &
          !  ElementIDArray_NPlus1, MaterialIDArray_NPlus1, &
          !  ShapeValuesArray_N, DShapeValuesArray_N, &
          !  ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
          !  GlobPosArray_N, GlobPosArray_NPlus1, &
          !  SigmaEffArray_N, SigmaEffArray_NPlus1)
          
          if (CalParams%ApplyEmptyElements.or.CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix) then
            call InitialiseModEmptyElementsData()
          end if
          
          allocate(NPartEle_N(Counters%NEl), stat = IError)
          NPartEle_N = 0 ! Initialise to zero
          
          allocate(NPartEle_NPlus1(Counters%NEl), stat = IError)
          NPartEle_NPlus1 = 0 ! Initialise to zero

          ! Fill particle house-keeping arrays with data from the BRF file of the previous load phase
          call ReadHouseKeepingDataFromFile(FileName)
          if (IsMPMWithMPIntegration()) then !Correction is required if change from MIXED- to MP-integration is considered in latest CPS
            IsParticleIntegration = .true.
          end if      
          
          ! Fill particle list with particle data from the BRF file of the previous load phase
          call ReadParticleDataFromFile(FileName)

          ! Fill particle list with stress strain parameters from the BRF file of the previous load phase
          call ReadParticleStateParametersFromFile(FileName)

          ! Fill particle list with dynamic related information from the BRF file of the previous load phase
          call ReadParticleDynamicDataFromFile(FileName)

          ! Fill the particle information related to consolidation from the BRF file of the previous load phase
          call ReadConsolidationDataFromFile(FileName)

          if(CalParams%NumberOfPhases==3) then
            ! Fill the particle information related to unsaturated material from the BRF file of the previous load phase
            call ReadConsolidationDataGasFromFile(FileName)
          end if

          ! rewrite the particle material parameters from new GOM file 
          ! if it is not the first phase and ApplyMaterialUpdate is true
          if (CalParams%ApplyMaterialUpdate) then
            call RewriteParticleMaterialParameters()
          end if
          
          ! Set particle shape function values and derivatives
          call InitialiseParticleShapeFunctionData(ShapeValuesArray, DShapeValuesArray)

          ! Adapt IsActiveElement to list of elements which contain particles
          call AdaptIsElemToParticleData()

          ! Set the particle IDCounter to the ID of the last read particle
          call SetIDCounter(MAXVAL(IDArray))

          if (.not.allocated(ActiveElement_N)) then
            call SetActiveElement(ActiveElement_N, ActiveNode_N)
          endif

          call SetParticleIndex(ActiveElement_N,NPartEle_N, &
            EleParticlesHelp_N, EleParticles_N)
          if (CalParams%ApplyEmptyElements.or.CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix) then
            ! Check whether elements are deactivated that should be considered fully filled
            call CheckEmptyElements()
          end if
          if (CalParams%ApplyEmptyElements) then
            call AdjustParticleDiscretisation()
          end if
          
          if (CalParams%ApplyContactAlgorithm) then
            call ReAssignEntityID()
          end if
        
        end subroutine ReinitialiseMaterialPoints
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #29 - BOOKMARK
        
        subroutine AdjustParticleDiscretisation()
        !**********************************************************************
        !
        !    Function:  Adds and removes particles from the discretisation.
        !
        !     FileName : Name of the project file (without extension)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: NUpdateParticles, NRemovedParticles, &
                     NAddedParticles, I
          real(REAL_TYPE) :: TotalVirtualParticleMass, &
                              TotalParticleMass

          ! Determine number of particles removed from the discretisation
          NRemovedParticles = CheckRemovalParticles()
          
          ! Determine number of particles added to the discretisation
          NAddedParticles = CalParams%NVirtualParticles * NEmptyElements
               
          ! Determine change of particle discretisation
          NUpdateParticles = NAddedParticles - NRemovedParticles

          if ((NRemovedParticles > 0) .or. (NAddedParticles > 0)) then

            if (CalParams%OutputDebugData) then
              call WriteInLogFile(' Removing ' // trim(String(NRemovedParticles)) // &
                                  ' particles from the discretisation.')

              call WriteInLogFile(' Adding '// trim(String(NAddedParticles)) // &
                                  ' particles to the discretisation.')
            end if

            ! Update array Particles, reset incorrect and redundant particle house keeping information
            call AdjustParticlesArrays(NUpdateParticles)
            ! Add particles to array Particles
            if (NAddedParticles > 0) then
              call FillParticlesArray(NUpdateParticles, NRemovedParticles)
            end if
            
             ! Adjust the remaining particle housekeeping data structure to the new discretisation
            call AdjustParticleHouseKeeping()
     
             ! Update remaining particle data which might depend on surrounding elements/particles
            call CompleteAddedParticlesData()

            TotalVirtualParticleMass = 0.0
            TotalParticleMass = 0.0
            do I = 1, Counters%NParticles
              if (IsNewVirtualParticle(I)) then
                Particles(I)%Kind = VIRTUALPARTICLE
              end if
              TotalParticleMass = TotalParticleMass + MassArray(I)
              if (IsVirtualParticle(I)) then
                TotalVirtualParticleMass =TotalVirtualParticleMass + MassArray(I)
              end if
            end do
            write(TVMunit, '(I12, G12.4, G12.4, I12, I12, I12)')  &
                                CalParams%TimeStep, &
                                TotalParticleMass, &
                                TotalVirtualParticleMass, &
                                Counters%NParticles, &
                                NAddedParticles, &
                                NRemovedParticles

          end if
                    
          ! Update auxiliary data
          IsEmptyElement = .false. ! All empty elements removed
          NEmptyElements = 0
    
        end subroutine AdjustParticleDiscretisation
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #30 - BOOKMARK
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #31 - BOOKMARK
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #32 - BOOKMARK
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #33 - BOOKMARK
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #34 - BOOKMARK
        
        
        
        
        
        
        
        
        
        
        
        
                                                 
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            
                                       
                                                 
                                                 
                                                 
                                                 
                                                 
                                                 
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            
            
            
        
            
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
              
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
        
                         
                         
                         
                         
                         
                         
                         
                         
        
        
      end module ModDynamicExplicit
