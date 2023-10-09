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

          !********** 4b - TIME STEP LOOP ******************************
          do while( (.not.CalParams%ConvergenceCheck%DoesConverge) .and. (.not.CalParams%ConvergenceCheck%DoesDiverge) )

            ! Increase time
            CalParams%TimeStep = CalParams%TimeStep + 1
            CalParams%TotalRealTime = CalParams%TotalRealTime + CalParams%TimeIncrement
            CalParams%OverallRealTime = CalParams%OverallRealTime + CalParams%TimeIncrement
            !call WriteInLogFile('TimeStep: ' // trim(String(CalParams%TimeStep))//' '//trim(String(CalParams%TimeIncrement)) &
            !                                                                    //' '//trim(String(CalParams%TotalRealTime)))
            !call WriteInLogFile('  Skipconvection? '//trim(String(MinimumDeterminantRatioReached))//' '// &
            !                                          trim(String(IsMPMSkipConvection())))

            call UpdateMultipliersForTimeDepencency() 
           
            ! SUBROUTINE #1
            call MapMassMP2Nodes(LumpedMassDry) 
            
            ! SUBROUTINE #2
            call MapMomentumMP2Nodes(Momentum)
            
            ! SUBROUTINE #3
            call GetNodalVelocityFromNodalMomentumImplicit(Momentum, LumpedMassDry, TotalVelocitySoil)
            
            ! SUBROUTINE #4
            call GetNodalExternalForces(ExtLoad, GravityLoad)

            ! SUBROUTINE #5
            call GetNodalIntForcesImplicit(IntLoad, SigmaEffArray) ! SigmaEffArray comes from existing stress

            ! SUBROUTINE #6
            call GetFirstGuessOfVelocity(TotalVelocitySoil_ImplicitIter)
            
            ! SUBROUTINE #7
            call GetMPStrainFromNodeVelocity(TotalVelocitySoil_ImplicitIter)!, ShapeValuesArray_N, DShapeValuesArray_N, ActiveElement_N)!, &
            
           ! SUBROUTINE #8: stresses, the water pressure and the gass pressure are updated (The mass balance equations are solved)
           call MPMDYNGetSigImplicit()!(SigmaEffArray) ! we are calculating new stress
                          
            ! SUBROUTINE #9: 
           call DynUpdateParticleWeights()
           
           call UpdateParticlePos()
           !ShapeValuesArray_N, DShapeValuesArray_N, &
           !                          ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
           !                          GlobPosArray_N, GlobPosArray_NPlus1, &
           !                          UStepArray, UPhaseArray, UArray, &
           !                          ElementIDArray_N, EleParticles_N, &
           !                         IncrementalDisplacementSoil_N)!ShapeValuesArray, DShapeValuesArray, ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
           
           ! SUBROUTINE #10:
                !Variables populated in UpdateParticleHouseKeeping:
			    !NPartEle
			    !EleParticleHelp
			    !IsActiveElement
                
                
                
            if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
                call ComputeWaterSpecificWeight() ! compute gamma_w = rhow*g(gravity_acc*current multiplier)
            end if
            
            ! Lagrangian phase
            call LagrangianPhase()
            call ApplyMPMDynamicContact() ! only if ApplyContactAlgorithm

            ! Convective phase
            call DYNConvectivePhase()

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

          end do ! ----- end TIME STEP LOOP -----

        end subroutine RunImplicitDynamicLoadStep
        
        
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        ! SUBROUTINE #1
        subroutine MapMassMP2Nodes(LumpedMass)
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
          real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES) :: ShapeValuesArray !, intent(in)
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(out) :: LumpedMass
          !real(REAL_TYPE), dimension(Counters%NParticles), intent(in) :: MassArray
         
          integer(INTEGER_TYPE), dimension(Counters%NEl) :: ActiveElement !, intent(in)
          
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
          
        
        
        end subroutine MapMassMP2Nodes
        
        
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        ! SUBROUTINE #2
        
        subroutine MapMomentumMP2Nodes(Momentum)
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
          
        end subroutine MapMomentumMP2Nodes  
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
         ! SUBROUTINE #3
        
        subroutine GetNodalVelocityFromNodalMomentumImplicit(Momentum, LumpedMassDry, TotalVelocitySoil)!, DoSystem)
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
          
        end subroutine GetNodalVelocityFromNodalMomentumImplicit
        
        
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
        ! SUBROUTINE #5 - BOOKMARK
        subroutine GetNodalIntForcesImplicit(IntLoad, SigmaEffArray)
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
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          
           !call MPMDYNBTSigOnly(IntLoad, FReaction, BulkViscLoad)
           call MPMDYNBTSigOnlyImplicit(IntLoad, SigmaEffArray) !FReaction, BulkViscLoad,


        end subroutine GetNodalIntForcesImplicit
        
        
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
        
        subroutine MPMDYNBTSigOnlyImplicit(InternalLD, SigmaEffArray) !Reactions, BulkViscousLD,
        !**********************************************************************
        !
        ! Function:  Calculation of the equivalent nodal forces due to
        !            a given stress InternalLD = Integral {BT*Sigma}.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(out) :: InternalLD
          real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR), intent(in) :: SigmaEffArray
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: Reactions
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: BulkViscousLD
          
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: ActiveElement
          
          
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
          real(REAL_TYPE), dimension(NTENSOR) :: S
          real(REAL_TYPE) :: Position
          real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)
          integer(INTEGER_TYPE) :: I, IntGlo, IEl, Int, J, NN, IAEl, NElemPart, iEntityID, iNode, IDof, MaterialIndex
          real(REAL_TYPE) :: WtN, Det, WPP, VPressure = 0.0
          logical :: DoConsiderReactionForces, IsUndrEffectiveStress
 
          
          
          
          IntGlo = 0
          S = 0.0
          InternalLD = 0.0
          !Reactions = 0.0
          !BulkViscousLD = 0.0

          do IAEl = 1, Counters%NAEl
           IEl = ActiveElement(IAEl)

           if(NFORMULATION==1) then
               
            if (IsParticleIntegration(IEl) ) then
              NElemPart = NPartEle(IEl)
            else
              NElemPart = ELEMENTGAUSSPOINTS
            end if
             
            call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN)
        
            do Int = 1, NElemPart
              
              IntGlo = GetParticleIndex(Int, IEl)
              MaterialIndex = MaterialIDArray(IntGlo)
              
              IsUndrEffectiveStress = &
              !code version 2016 and previous
              ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
              !code version 2017.1 and following
              (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
              
              if (IsParticleIntegration(IEl) ) then
                !if ( ISAXISYMMETRIC ) then
                !  ! the integration weight of the MP is not corrected due to the axisymmetry
                !  Position = GlobPosArray(IntGlo, 1) ! INDEX_X = 1
                !  ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                !end if
                WTN = Particles(IntGlo)%IntegrationWeight
              else
                !if ( ISAXISYMMETRIC ) then
                !  Position = GPGlobalPositionElement(1, Int, IEl) !(GaussPointGlobalCoord(INDEX_X)
                !  ShapeValues(:) = GPShapeFunction(Int, :)
                !  WTN = WTN * Position ! the volume of the element is corrected due to the axisymmetry
                !end if
              end if
              
              if (IsUndrEffectiveStress) then
                  WPP = Particles(IntGlo)%WaterPressure * WTN
              end if    
              do J = 1, 3 ! always first three stress components, both for 2D and 3D
                !if (IsUndrEffectiveStress.or.CalParams%ApplyImplicitQuasiStatic) then
                !  WPP = Particles(IntGlo)%WaterPressure * WTN 
                !  S(J) = SigmaEffArray(IntGlo, J) * WTN + WPP
                !else
                  S(J) = SigmaEffArray(IntGlo, J) * WTN
                !end if
              end do
              do J = 4, NTENSOR
                S(J) = SigmaEffArray(IntGlo,J)  * WTN
              end do
                         
              !if (CalParams%ApplyBulkViscosityDamping) then
              !  VPressure = Particles(IntGlo)%DBulkViscousPressure * WTN
              !end if

              !if (CalParams%ApplyContactAlgorithm) then
              !  iEntityID = EntityIDArray(IntGlo) 
              !else 
                iEntityID = 1
              !end if        
                
              do iNode=1,ELEMENTNODES
                nn=ElementConnectivities(iNode,iel)
                IDof = ReducedDof(nn)
               
                if ( NVECTOR == 3 ) then ! 3D
                  InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4) + B(3,iNode)*S(6)
                  InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2) + B(3,iNode)*S(5)
                  InternalLD(IDof+3, iEntityID) = InternalLD(IDof+3, iEntityID) + B(1,iNode)*S(6) + B(2,iNode)*S(5) + B(3,iNode)*S(3)
                else if ( NVECTOR == 2 ) then ! 2D  
                  InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4) 
                  InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2)
                  !if ( ISAXISYMMETRIC ) then
                  !  ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912 
                  !  InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + S(3) * ShapeValues(iNode) / Position
                  !end if
                else
                  call GiveError('Dimension is not correct. It must be 2 or 3.')
                end if
                !
                !if (CalParams%ApplyBulkViscosityDamping) then
                !  do I = 1, NVECTOR  
                !    BulkViscousLD(IDof+I, IEntityID) = BulkViscousLD(IDof+I, IEntityID) + B(I, INode) * VPressure
                !    if ( ISAXISYMMETRIC ) then
                !      ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                !      BulkViscousLD(IDof+I, IEntityID) = BulkViscousLD(IDof+I, IEntityID) + VPressure * ShapeValues(iNode) / Position
                !    end if
                !  end do  
                !end if

                !if ( NVECTOR == 3 ) then ! 3D
                !  DoConsiderReactionForces = CalParams%PrescribedVelo%ApplyPrescribedVelo
                !else if ( NVECTOR == 2 ) then ! 2D  
                !  DoConsiderReactionForces =  CalParams%PrescribedVelo%ApplyPrescribedVelo            
                !end if

                !if (DoConsiderReactionForces) then
                !  if ( NVECTOR == 3 ) then ! 3D
                !    Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + B(1, INode) * S(1) + B(2, INode) * S(4) + B(3, INode) * S(6) + B(1, INode) * VPressure
                !    Reactions(IDof+2, IEntityID) = Reactions(IDof+2, IEntityID) + B(1, INode) * S(4) + B(2, INode) * S(2) + B(3, INode) * S(5) + B(2, INode) * VPressure
                !    Reactions(IDof+3, IEntityID) = Reactions(IDof+3, IEntityID) + B(1, INode) * S(6) + B(2, INode) * S(5) + B(3, INode) * S(3) + B(3, INode) * VPressure
                !  else if ( NVECTOR == 2 ) then ! 2D  
                !    Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + B(1, INode) * S(1) + B(2, INode) * S(4) + B(1, INode) * VPressure
                !    Reactions(IDof+2, IEntityID) = Reactions(IDof+2, IEntityID) + B(1, INode) * S(4) + B(2, INode) * S(2) + B(2, INode) * VPressure
                !    if ( ISAXISYMMETRIC ) then
                !      ! Equation below is to be double checked
                !      Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + S(3) * ShapeValues(iNode) / Position
                !    end if
                !  else
                !    call GiveError('Dimension is not correct. It must be 2 or 3.')
                !  end if
                !end if
                
              end do
            end do
           end if
          end do


        end subroutine MPMDYNBTSigOnlyImplicit
        
        
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
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #7 - BOOKMARK
         subroutine GetMPStrainFromNodeVelocity(TotalVelocitySoil)
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
          call GetNodalIncrementalDisplacement(IncrementalDisplacementSoil, TotalVelocitySoil)

          AccumulatedIncDisplacementSoil = AccumulatedIncDisplacementSoil + IncrementalDisplacementSoil

          ! SUBROUTINE #2: calculate IncrementalDisplacementSoil
          call UpdateNodalTotalDisplacementImplicit(EntityUsed, IncrementalDisplacementSoil, TotalDisplacementSoil)
    
          ! SUBROUTINE #3: Update particle strain data
          call UpdateParticleStrainsImplicit(IncrementalDisplacementSoil)!, DShapeValuesArray, ActiveElement, NPartEle)!, GetParticleIndex_Local) ! update particle strains for solid and mixture
                                                        
         end subroutine GetMPStrainFromNodeVelocity
        
        
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
         !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #8 - BOOKMARK
        subroutine MPMDYNGetSigImplicit()!SigmaEffArray)!(IntLoad_NPlus1) !MPMDYNGetSig ! NAEl, 
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
    
                 call StressSolidImplicit(IntGlo, IEl, BMatrix, IEntity, SigmaEff0Array, SigmaEffArray,  ESMStatevArray)   ! calculate stresses for mixture or solid material points
                
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

        end subroutine MPMDYNGetSigImplicit
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
         subroutine StressSolidImplicit(IDpt, IDel, BMatrix, IEntityID, SigmaEff0Array, SigmaEffArray,  ESMstatevArray) !MaterialIDArray, 
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
    real(REAL_TYPE), dimension(Counters%NParticles, NSTATEVAR), intent(inout) :: ESMstatevArray
    
    
    
    pointer (p, ESM)             
          
    ! get constitutive model in integration/material point
    IDset = MaterialIDArray(IDpt) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
    NameModel = MatParams(IDset)%MaterialModel ! name of constitutive model as specified in GOM-file
    ntens = MatParams(IDset)%UMATDimension  ! 2D or 3D formulation of the External soil model   
          
    ! get strain increments in integration/material point
    TempStrainIncr = GetEpsStep(Particles(IDpt)) ! incremental strain vector assigned to point
    
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
    ESMstatevArray(IDpt,:) = StateVar ! temporary
          
    call CalculatePrincipalStresses(IDpt, Stress(1:NTENSOR), StressPrinc)
    call AssignStressStrainToGlobalArrayESM(IDpt, NTENSOR, StressIncr, StressPrinc, StrainIncr)!, SigmaEffArray )
    
    ! write plasticity state to global array
    !  call SetIPL(IDpt, IDel, int(StateVar(50)))
    !if (CalParams%ApplyBulkViscosityDamping) then
    !RateVolStrain(IDEl) = DEpsVol / CalParams%TimeIncrement
    !call CalculateViscousDamping_interface(IDpt, IDEl)
    !end if  
         end subroutine StressSolidImplicit
         
         !--------------------------------------------------------------------------------------------------------
         !--------------------------------------------------------------------------------------------------------
         
         
         !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        ! SUBROUTINE #9 - BOOKMARK
        subroutine DynUpdateParticleWeightsImplicit()!( ShapeValuesArray, ActiveElement, NPartEle )
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
                                
          do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            IEl = ActiveElement(IAEl)   
            
            do Int = 1, NPartEle(IEl) ! Loop over all material points of the element
                
              ! Determine global ID of Material Point 
              IntGlo = GetParticleIndex(Int, IEl)
              
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
          
         
        end subroutine DynUpdateParticleWeightsImplicit
        
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        
         
         
      end module ModDynamicImplicit
