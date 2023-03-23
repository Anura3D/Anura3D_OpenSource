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
	  
	  
	  module ModDynamicExplicit
      !**********************************************************************
      !
      !    Function:  This module contains the routine for running a dynamic
      !               load step using explicit time integration.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr. 2022) $
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
      use ModDynamicImplicitGeneralizedAlphaScheme
      
      implicit none

      contains ! Routines of this module

        !subroutine RunExplicitDynamicLoadStep()
        !!**********************************************************************
        !!
        !!  Function:  Routine called from the main routine for performing a
        !!             dynamic load step using explicit time integration.
        !!             Contains the 'time step loop' of explicit integration.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum_NPlus1
        !
        !  !********** 4b - TIME STEP LOOP ******************************
        !  do while( (.not.CalParams%ConvergenceCheck%DoesConverge) .and. (.not.CalParams%ConvergenceCheck%DoesDiverge) )
        !
        !    ! Increase time
        !    CalParams%TimeStep = CalParams%TimeStep + 1
        !    CalParams%TotalRealTime = CalParams%TotalRealTime + CalParams%TimeIncrement
        !    CalParams%OverallRealTime = CalParams%OverallRealTime + CalParams%TimeIncrement
        !    call WriteInLogFile('TimeStep: ' // trim(String(CalParams%TimeStep))//' '//trim(String(CalParams%TimeIncrement)) &
        !                                                                        //' '//trim(String(CalParams%TotalRealTime)))
        !    call WriteInLogFile('  Skipconvection? '//trim(String(MinimumDeterminantRatioReached))//' '// &
        !                                              trim(String(IsMPMSkipConvection())))
        !
        !    call UpdateMultipliersForTimeDepencency() 
        !
        !    if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
        !        call ComputeWaterSpecificWeight() ! compute gamma_w = rhow*g(gravity_acc*current multiplier)
        !    end if
        !    
        !    ! Generalized alpha variable initialization 
        !    call InitializeGeneralizedAlphaTimeIntegrationVariables()!(Gamma_GAlpha, Beta_GAlpha, AlphaM_GAlpha, AlphaF_GAlpha)
        !
        !    ! Predictor phase 
        !    call PredictorPhase() ! -> Predict NPlus1 values  
        !
        !    !! Projection phase
        !    !call ProjectionPhase()
        !    
        !    
        !    
        !    
        !    ! In a way, the first explicit 'cycle' is kind of an initialization process 
        !    !-----------------------------------------------------------------------
        !    ! Lagrangian phase
        !    call LagrangianPhase() ! -> particles to nodes -> update AccelerationSoil 
        !    call ApplyMPMDynamicContact() ! only if ApplyContactAlgorithm
        !
        !    ! Convective phase
        !    call DYNConvectivePhase() ! -> nodes to particles -> update UArray, VeloArray
        !    !-----------------------------------------------------------------------
        !    
        !    call MapMomentumAndMassP2N(Momentum_NPlus1, LumpedMassDry_NPlus1) ! we need the "N+1" mass from here 
        !    
        !    ! Thinking of the tensile bar problem --> we would have AccelerationSoil to be 
        !    ! all zeros except the boundary condition nodes 
        !    ! so far we should be okay here....
        !    
        !    
        !    if (CalParams%ImplicitGeneralizedAlphaTimeScheme%Switch == .TRUE. ) then 
        !        ! we have the first estimate of the values 
        !        call DynamicImplicitGeneralizedAlphaScheme_() ! iterative implicit scheme currently being coded
        !    end if
        !    
        !    !call GetNodalExtAndIntForces_NPlus1() ! rotated to local coordinate system
        !    
        !    
        !    
        !    ! Optional writing of time step results
        !    ! write intermediate time step results and if ApplyQuickCheckOutput
        !    call WriteTimeStepResults(.false.)
        !
        !    ! Convergence check
        !    CalParams%ConvergenceCheck%DoesConverge = ConvergenceCheck()
        !
        !    ! Divergence check
        !    if (CalParams%ConvergenceCheck%ApplyDivergenceCheck) then
        !      CalParams%ConvergenceCheck%DoesDiverge = DivergenceCheck()
        !
        !      if (CalParams%ConvergenceCheck%DoesDiverge) then
        !        call GiveWarning('Calculation is diverging...')
        !      end if
        !    end if
        !
        !    !  if (CalParams%OutputDebugData) then
        !    ! Time step output
        !    if (NFORMULATION==1) then
        !      call EnergyOutput()
        !    else
        !      call EnergyOutput2LayForm()
        !    end if
        !
        !    call MaterialPointOutput() ! only if OutputNumberParticles>0
        !    
        !    !compute and write on file the sum of nodal reactions on selected surfaces
        !    call SurfaceReactionOutput() 
        !    
        !    ! end if
        !
        !    ! update ElementLMin in case of ULFEM
        !    if (IsULFEMComputation()) then
        !      call DetermineElementLMin() ! calulate minimum element altitude
        !    end if
        !
        !    call CalculateCriticalTimeStep()
        !    call WriteFEMNodeData() ! only if .not.IsMPMComputation
        !
        !  end do ! ----- end TIME STEP LOOP -----
        !
        !end subroutine RunExplicitDynamicLoadStep
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! IMPLICIT GENERALIZED ALPHA 
        !
        !
        !
        !
        
        
        
        
        
        
        
        subroutine RunImplicitDynamicLoadStep()
        !**********************************************************************
        !
        !  Function:  Routine called from the main routine for performing a
        !             dynamic load step using explicit time integration.
        !             Contains the 'time step loop' of explicit integration.
        !
        !**********************************************************************

        implicit none
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum_NPlus1

          !********** 4b - TIME STEP LOOP ******************************
          do while( (.not.CalParams%ConvergenceCheck%DoesConverge) .and. (.not.CalParams%ConvergenceCheck%DoesDiverge) )

            ! Increase time
            CalParams%TimeStep = CalParams%TimeStep + 1 ! CalParams%TimeStep = 0 initially 
            ! CalParams%TimeIncrement need to be user defined in the implicit scheme
            CalParams%TotalRealTime = CalParams%TotalRealTime + CalParams%TimeIncrement
            CalParams%OverallRealTime = CalParams%OverallRealTime + CalParams%TimeIncrement
            call WriteInLogFile('TimeStep: ' // trim(String(CalParams%TimeStep))//' '//trim(String(CalParams%TimeIncrement)) &
                                                                                //' '//trim(String(CalParams%TotalRealTime)))
            call WriteInLogFile('  Skipconvection? '//trim(String(MinimumDeterminantRatioReached))//' '// &
                                                      trim(String(IsMPMSkipConvection())))

            call UpdateMultipliersForTimeDepencency() 

            if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
                call ComputeWaterSpecificWeight() ! compute gamma_w = rhow*g(gravity_acc*current multiplier)
            end if
            
            ! Generalized alpha variable initialization 
            !call InitializeGeneralizedAlphaTimeIntegrationVariables() !(Gamma_GAlpha, Beta_GAlpha, AlphaM_GAlpha, AlphaF_GAlpha)

            ! Set initial conditions here 
            !call SetInitialConditionsOnPreviousVectors() 
            
            ! Predictor phase 
            !call PredictorPhase() ! -> Predict NPlus1 values  
            
            
            !AccelerationSoil(2,1) = 9.81
            !AccelerationSoil(4,1) = 9.81
            
            
            ! Do I even need velocity and displacement???
            ! WE DO NOT NEED TO GO THROUGH THE EXPLICIT CYCLE!
            !TotalVelocitySoil(2,1) = TotalVelocitySoilPrevious(2,1) + (AccelerationSoil(2,1)*CalParams%TimeIncrement)
            !TotalVelocitySoil(4,1) = TotalVelocitySoilPrevious(4,1) + (AccelerationSoil(4,1)*CalParams%TimeIncrement)
            
            
            ! Projection phase
            !call ProjectionPhase()
            
            
            

            
            ! In a way, the first explicit 'cycle' is kind of an initialization process 
            !-----------------------------------------------------------------------
            ! Lagrangian phase
            !call LagrangianPhase() ! -> particles to nodes -> update AccelerationSoil 
            !call ApplyMPMDynamicContact() ! only if ApplyContactAlgorithm

            ! Convective phase
            !call DYNConvectivePhase_InputAccelerationSoil(AccelerationSoil, ForceInt) ! -> nodes to particles -> update UArray, VeloArray
            !-----------------------------------------------------------------------
            ! updating nodal values of mass using MapMomentumAndMassP2N... we do not need Momemtum value 
            !call MapMomentumAndMassP2N(Momentum_NPlus1, LumpedMassDry_NPlus1) ! we need the "N+1" mass from here 
            
            ! Thinking of the tensile bar problem --> we would have AccelerationSoil to be 
            ! all zeros except the boundary condition nodes 
            ! so far we should be okay here....
            
            
            ! Before we correct acceleration, we check the tolerance values here         
            !call CalculateGeneralizedAlphaResidualEquation()!Residual_Implicit)
            
            ! Step 2 in Sulsky and Kaul (2004)
            !call CheckResidualTolerance()!ImplicitConverged, Tolerance_ImplicitConvergence)
            
            ! Step 3 and 4 combined --> evaluation of the directional derivative 
            !call EstimateDirectionalDerivativeValue()
            
            !CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector
            
            
            if (CalParams%ImplicitGeneralizedAlphaTimeScheme%Switch == .TRUE. ) then 
                ! we have the first estimate of the values 
                call DynamicImplicitScheme_() ! iterative implicit scheme currently being coded
            end if
            
            !call GetNodalExtAndIntForces_NPlus1() ! rotated to local coordinate system
            
            
            
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine SetInitialConditionsOnPreviousVectors()
        ! set initial values 
        implicit none 
        
        
        AccelerationSoilPrevious = 0 
        TotalVelocitySoilPrevious = 0 
        TotalDisplacementSoilPrevious = 0 
        IncrementalDisplacementSoilPrevious = 0 
        
        
        end subroutine SetInitialConditionsOnPreviousVectors
        
        
        
        
        
        
        
        
        
        subroutine ProjectionPhase()
        
        implicit none 
        
        AccelerationSoilAlphaM = ((1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m)*AccelerationSoil) + &
            (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m*AccelerationSoilPrevious)
        
        
        ! we do not need this velocity value in the residual evaluation
        !TotalVelocitySoilAlphaF = ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f)*TotalVelocitySoil) + &
        !    (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f*TotalVelocitySoilPrevious) ! these need to be updated based on the 9.81 m/s2 value 
        !
        !! we do not need this displacement value in the residual evaluation
        !TotalDisplacementSoilAlphaF = ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f)*TotalDisplacementSoil) + &
        !    (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f*TotalDisplacementSoilPrevious)
        
        !
        !AccelerationSoil + (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m*(AccelerationSoil - AccelerationSoilPrevious) )
        !
        !TotalVelocitySoilAlphaF = TotalVelocitySoilPrevious + (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f*(TotalVelocitySoil - TotalVelocitySoilPrevious) )
        !
        !TotalDisplacementSoilAlphaF = TotalDisplacementSoilPrevious + (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f*(TotalDisplacementSoil - TotalDisplacementSoilPrevious) )
        !
        
        
        end subroutine ProjectionPhase
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !subroutine EstimateDirectionalDerivativeValue
        !
        !implicit none 
        !
        !real(REAL_TYPE) :: DirectionalDerivative 
        !real(REAL_TYPE), dimension(Counters%N,1) :: IncrementalAcceleration_SS ! this should be the same size as acceleration soil 
        !real(REAL_TYPE), dimension(Counters%N,1) :: AdditionalTermForResidualEvaluation
        !real(REAL_TYPE), dimension(Counters%N,1) :: AccelerationSoil_DirectionalDerivative
        !real(REAL_TYPE) :: SS_Mod
        !real(REAL_TYPE) :: HH
        !real(REAL_TYPE) :: AccelerationSoil_Mod
        !
        !!integer(INTEGER_TYPE) :: ii
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !! input
        !real(REAL_TYPE) :: Mass_I_NPlus1 
        !real(REAL_TYPE) :: Mass_I_N
        !
        !
        !real(REAL_TYPE) :: Acceleration_I_NPlus1
        !real(REAL_TYPE) :: Acceleration_I_N
        !
        !
        !real(REAL_TYPE) :: ForceExt_I_NPlus1 
        !real(REAL_TYPE) :: ForceExt_I_N
        !
        !
        !real(REAL_TYPE) :: ForceInt_I_NPlus1
        !real(REAL_TYPE) :: ForceInt_I_N
        !
        !
        !real(REAL_TYPE) :: AlphaM_GAlpha
        !real(REAL_TYPE) :: AlphaF_GAlpha
        !
        !
        !! local 
        !real(REAL_TYPE), dimension(Counters%N) :: Residual_LHS
        !real(REAL_TYPE), dimension(Counters%N) :: Residual_RHS
        !real(REAL_TYPE), dimension(Counters%N) :: Residual_Implicit!, intent(out) 
        !
        !real(REAL_TYPE), dimension(Counters%N) :: DirectionalDerivativeValue
        !
        !real(REAL_TYPE) :: Residual_Implicit_X
        !real(REAL_TYPE) :: Residual_Implicit_Y
        !
        !integer(INTEGER_TYPE) :: ii, jj, kk
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !! This is SS in Sulsky and Kaul (2004)
        !!IncrementalAcceleration_SS(1,1) = 1
        !!IncrementalAcceleration_SS(2,1) = 1
        !!IncrementalAcceleration_SS(3,1) = 1
        !!IncrementalAcceleration_SS(4,1) = 1
        !!IncrementalAcceleration_SS(5,1) = 1
        !!IncrementalAcceleration_SS(6,1) = 1
        !!IncrementalAcceleration_SS(7,1) = 1
        !!IncrementalAcceleration_SS(8,1) = 1
        !
        !AccelerationSoil_Mod = 0
        !
        !AccelerationSoil_DirectionalDerivative = 0
        !
        !SS_Mod = 0
        !
        !HH = 1e-5
        !
        !do ii = 1, Counters%N
        !
        !    AccelerationSoil_Mod = AccelerationSoil_Mod + AccelerationSoil(ii,1)**2
        !
        !    SS_Mod = SS_Mod + IncrementalAcceleration_SS(ii,1)**2
        !    
        !end do 
        !
        !AccelerationSoil_Mod = sqrt(AccelerationSoil_Mod)
        !
        !SS_Mod = sqrt(SS_Mod)
        !
        !
        !
        !AdditionalTermForResidualEvaluation = ((HH*AccelerationSoil_Mod)/SS_Mod ) * IncrementalAcceleration_SS
        !
        !
        !AccelerationSoil_DirectionalDerivative = AdditionalTermForResidualEvaluation + AccelerationSoil 
        !
        !! I need to evaluate the residual based on the new AccelerationSoil_DirectionalDerivative value 
        !
        !
        !
        !
        !
        !! This is called D_h * F
        !!DirectionalDerivative_Dh_F =  
        !
        !!CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector
        !
        !
        !
        !!if 
        !!
        !!CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector ! this is equivalent to F(v) in Sulsky and Kaul (2004)
        !!
        !!
        !!DirectionalDerivative = 
        !!end if
        !
        !
        !
        !
        !!Note that you need to update the mass corresponding to that acceleration AccelerationSoil_DirectionalDerivative
        !
        !
        !do ii = 1, Counters%N
        !
        !    Residual_LHS(ii) = ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m) * LumpedMassDry_NPlus1(ii,1) * AccelerationSoil_DirectionalDerivative(ii,1) ) + ( CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m * LumpedMassDry(ii,1) * AccelerationSoilPrevious(ii,1) )
        !    
        !    Residual_RHS(ii) = ( ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f) * ExtLoad(ii,1) ) + ( CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f * ExtLoad(ii,1) ) )      -       &
        !             ( ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f) * IntLoad(ii,1) ) + ( CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f * IntLoadPrevious(ii,1) ) )
        !
        !end do 
        !
        !
        !
        !
        !
        !Residual_Implicit = Residual_LHS - Residual_RHS
        !
        !ResidualLoadGeneralizedAlpha(:,1) = Residual_Implicit
        !
        !Residual_Implicit_X = 0
        !
        !Residual_Implicit_Y = 0
        !
        !kk = 1
        !do jj = 1, Counters%N/2
        !!ForceExt
        !    ! looking at elements 1, 3, 5, 7
        !    Residual_Implicit_X = Residual_Implicit_X + ResidualLoadGeneralizedAlpha(kk,1)
        !    
        !    
        !    !Residual_Implicit_Y(jj)
        !
        !    kk = kk+2
        !end do 
        !
        !kk = 2
        !do jj = 1, Counters%N/2
        !!ForceExt
        !    ! looking at elements 1, 3, 5, 7
        !    Residual_Implicit_Y = Residual_Implicit_Y + ResidualLoadGeneralizedAlpha(kk,1)
        !    
        !    
        !    !Residual_Implicit_Y(jj)
        !
        !    kk = kk+2
        !end do 
        !
        !
        !CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector_DirectionalDerivative(1) = Residual_Implicit_X
        !
        !CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector_DirectionalDerivative(2) = Residual_Implicit_Y
        !
        !
        !
        !
        !DirectionalDerivativeValue = (SS_Mod/ (HH*AccelerationSoil_Mod) ) * (CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector_DirectionalDerivative - CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector)
        !
        !
        !end subroutine EstimateDirectionalDerivativeValue
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
      
      end module ModDynamicExplicit
