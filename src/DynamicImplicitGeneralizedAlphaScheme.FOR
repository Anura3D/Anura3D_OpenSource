      !*****************************************************************************
	  !                                       ____  _____  
      !           /\                         |___ \|  __ \ 
      !          /  \   _ __  _   _ _ __ __ _  __) | |  | |
      !         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
      !        / ____ \| | | | |_| | | | (_| |___) | |__| |
      !       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/ 
      !
	  !
	  !	  Anura3D - Numerical modelling and simulation of large deformations 
	  !   and soil–water–structure interaction using the material point method (MPM)
      !
	  !	  Copyright (C) 2022  Members of the Anura3D MPM Research Community 
	  !   (See Contributors file "Contributors.txt")
	  !
      !	  This program is free software: you can redistribute it and/or modify
      !	  it under the terms of the GNU Lesser General Public License as published by
      !	  the Free Software Foundation, either version 3 of the License, or
      !	  (at your option) any later version.
	  !
      !	  This program is distributed in the hope that it will be useful,
      !	  but WITHOUT ANY WARRANTY; without even the implied warranty of
      !	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      !	  GNU Lesser General Public License for more details.
	  !
      !	  You should have received a copy of the GNU Lesser General Public License
      !	  along with this program.  If not, see <https://www.gnu.org/licenses/>.
	  !
	  !*****************************************************************************
	  
	  
	  module ModDynamicImplicitGeneralizedAlphaScheme
      !**********************************************************************
      !
      !    Function:  
      !
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 9748 $
      !     $Date: 2022-04-21 14:53:15 +0200 (do, 21 apr. 2022) $
      !
      !**********************************************************************
        
      ! check which of these are needed in this module
      use ModCounters 
      use ModReadCalculationData
      use ModMPMData
      use ModParticle
      use ModMeshInfo
      use ModWriteTestData
      use ModEmptyElements
      use ModConvectivePhase
      use ModWriteMPMData
      use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
      
      use ModLagrangianPhase 
      use ModDYNConvectivePhase
      
      use ModConvectivePhase
          
      
      

      implicit none

      
      real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Arnoldi_Compiled
    
    
    
        contains ! Routines of this module

    
    
    
        ! ------------ main implicit kernel    
    
        
        subroutine DynamicImplicitGeneralizedAlphaScheme_()
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
        
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum_NPlus1
        
        Logical :: ImplicitConverged ! 1/0 switch 
        
        real(REAL_TYPE) :: Residual_Implicit_Magnitude
      
        ! Acceleration increment
        !real(REAL_TYPE), dimension(Counters%N) :: DeltaAcceleration
        real(REAL_TYPE), dimension(Counters%N) :: DeltaVelocity
        
        real(REAL_TYPE), dimension(Counters%N) :: Residual_Zero
        
        ! Note that this needs to be 
        real(REAL_TYPE), dimension(Counters%N) :: PP_1
        real(REAL_TYPE), dimension(Counters%N) :: PP_2
        
        real(REAL_TYPE), dimension(Counters%N) :: Beta
        
        real(REAL_TYPE) :: Rho_Zero
        
        real(REAL_TYPE) :: Iteration_kk 
        
        real(REAL_TYPE) :: Residual_Zero_Magnitude
        
        real(REAL_TYPE) :: TOLERANCE_EVALUATION
        
        real(REAL_TYPE) :: TOLERANCE_GMRES
        
        real(REAL_TYPE), dimension(Counters%N) :: DirectionalDerivative 
        
        integer(INTEGER_TYPE) :: ii, jj
                 
        ! This should not be kept as dimension(1)
        real(REAL_TYPE), dimension(1) :: h_kk ! this should be an allocatable and should be destroyed and expanded with each iteration

        ! I am hardcoding the iteration number (kk) and ii is basically from 1 to kk  
        real(REAL_TYPE) :: h_kPlus1_2_i_1
        
        real(REAL_TYPE) :: h_k_1_i_1
        
        integer(INTEGER_TYPE) :: nn_iteration_number
        
        real(REAL_TYPE) :: TOLERANCE
        
        ! GMRES example --------------------------
        !real(REAL_TYPE), dimension(2,2) :: AA
        !real(REAL_TYPE), dimension(2) :: bb
        !real(REAL_TYPE), dimension(2) :: xx  
        
        
          real(REAL_TYPE), allocatable, dimension(:,:) :: AA
          real(REAL_TYPE), allocatable, dimension(:) :: bb
          real(REAL_TYPE), allocatable, dimension(:) :: xx
          
          integer(INTEGER_TYPE) :: max_iteration 
          
          integer(INTEGER_TYPE) :: nn_length_AA
          integer(INTEGER_TYPE) :: mm_max_iteration 
          real(REAL_TYPE) :: threshold
        !--------------------------------------------------
        
        
        

          threshold = 1e-6!10e-3
           DeltaVelocity = 0
          !DeltaAcceleration = 0
        
          ! We need the following varaibles to use in the generalized alpha time scheme: 
          !     1) CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m
          !     2) CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f
          !     3) CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha
          !     4) CalParams%ImplicitGeneralizedAlphaTimeScheme%Delta
        
            
          !  Step 1: Initialize generalized alpha variables 
          !call InitializeGeneralizedAlphaTimeIntegrationVariables()!(Gamma_GAlpha, Beta_GAlpha, AlphaM_GAlpha, AlphaF_GAlpha)
        
          ! Choose initial value for acceleration
          AccelerationSoil = 0 
          TotalVelocitySoil = 0 
          nn_iteration_number = 0 
        
          
          ! Step 2: Initialize initial conditions/state variables for nodes and material points. This includes:
          !                           Nodes                       Material points
          !     1) Acceleration     - AccelerationSoil          - AccelerationArray
          !     2) Velocity         - TotalVelocitySoil         - VelocityArray
          !     3) Displacement     - TotalDisplacementSoil     - UArray 
          !     4) Stress           - N/A                       - 
          !     5) InternalForce    - IntLoad                   - N/A
          !     6) ExternalForce    - ExtLoad                   - N/A
          !     7) Mass             - LumpedMassDry             - N/A
          !     8) Global position  - N/A                       - GlobPosArray
          
          
          ! Let's set initial acceleration and velocity to zero
          !AccelerationSoil = 0 
          !TotalVelocitySoil = 0 
          !nn_iteration_number = 0 
          
          
          call LagrangianPhase_GeneralizedAlphaInitialization(AccelerationSoil, RateofMomentum, ShapeValuesArray, DShapeValuesArray, IntLoad, ExtLoad) !-> AccelerationSoil, RateOfMomentum are nodal outputs here
                                                                                                !-> No outputs here
          !We don't need DShapeValuesArray as an input
          !In the first time step we assume that all the material points are the gauss points
          
          !MM_exp = LumpedMassDry
          !
          !IntLoad
          
          IntLoadPrevious = IntLoad
          
          call DYNConvectivePhase_GeneralizedAlphaInitialization(AccelerationSoil, & ! nodal input 
                    AccelerationArray, VelocityArray, UArray, &                      ! MP output  
                    TotalVelocitySoil, IncrementalDisplacementSoil, &
                    ShapeValuesArray, DShapeValuesArray, & 
                    ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
                    IntLoad)                       
          
                    !-> nodal output 
          !Velocity, Displacement, Global position, Stress, IntForce 
          
          ! Step 3: Update LumpedMassDry based on the initial acceleration 
          call UpdateMassAfterDYNConvectivePhase(LumpedMassDry_NPlus1, ShapeValuesArray_NPlus1)
          
          ! Step 4: Predictor phase for new acceleration 
          !call PredictorPhase()
          !call DYNConvectivePhase_GeneralizedAlphaInitialization(AccelerationSoil_NPlus1, & ! nodal input 
          !          AccelerationArray_NPlus1, VelocityArray_NPlus1, UArray_NPlus1, &                      ! MP output  
          !          TotalVelocitySoil_NPlus1, IncrementalDisplacementSoil_NPlus1, &
          !          ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, & !these two should be identical because we are inputting them and then overwriting them 
          !          ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
          !          IntLoad_NPlus1)     
          
          !call DYNConvectivePhase_GeneralizedAlphaInitialization(AccelerationSoil_NPlus1 ) !-> Velocity, Displacement, Global position, Stress, IntForce 
          ! Step 5: Update LumpedMassDry based on the initial acceleration 
          !call UpdateMassAfterDYNConvectivePhase(LumpedMassDry_NPlus1, ShapeValuesArray_NPlus1 )
          
          
          ! Step 7: Evaluate the residual of the linear momentum equation 
          !call EvaluateResidualValue(AccelerationSoil_NPlus_1MinusAlphaM_1, LumpedMassDry_NPlus_1MinusAlphaM_1, ExtLoad_NPlus_1MinusAlphaF_1, IntLoad_NPlus_1MinusAlphaF_1, Residual_GeneralizedAlpha)
          !call EvaluateResidualValueVelocity(AccelerationSoil_NPlus_1MinusAlphaM_1, LumpedMassDry_NPlus_1MinusAlphaM_1, ExtLoad_NPlus_1MinusAlphaF_1, IntLoad_NPlus_1MinusAlphaF_1, Residual_GeneralizedAlpha)
          call EvaluateResidualValueVelocity(TotalVelocitySoil_NPlus1, LumpedMassDry_NPlus1, &
                                                   IntLoad_NPlus1, IntLoad, &
                                                   TotalVelocitySoil, LumpedMassDry, &
                                                   Residual_GeneralizedAlpha)
          
          !subroutine EvaluateResidualValueVelocity(TotalVelocitySoil_NPlus1, LumpedMassDry_NPlus1, &
          !                                         IntLoad_NPlus1, IntLoad, &
          !                                         TotalVelocitySoil, LumpedMassDry, &
          !                                         Residual)
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          

          !From velocity standpoint, I need the following: 
          ! 1) previous mass, MM_exp 
          ! 2) current mass, MM_imp 
          ! 3) explicit velocity, vv_exp
          ! 4) implicit velocity, vv_imp
          ! 5) previous internal force, FF_int_exp
          ! 6) current internal force, FF_int_imp
          ! 7) previous external force, FF_ext_exp
          ! 8) current external force, FF_ext_imp
          
          
          
          
          
          !Residual = (MM_imp * vv_imp) - (Implicit_mu * CalParams%TimeIncrement * (FF_int_imp - FF_int_exp) ) &
          !                             - (Implicit_mu * CalParams%TimeIncrement * (FF_ext_imp - FF_ext_exp) ) &
          !                             - (MM_exp * vv_imp)
          !
          !
          !
          !
          !
          !
          !call EvaluateResidualValueVelocity(AccelerationSoil_NPlus_1MinusAlphaM_1, LumpedMassDry_NPlus_1MinusAlphaM_1, ExtLoad_NPlus_1MinusAlphaF_1, IntLoad_NPlus_1MinusAlphaF_1, Residual_GeneralizedAlpha)
          !
          !TotalVelocitySoil
          
          ! Step 6: Evaluate terms for residual of the linear momentum balance equation 
          !call EvaluatePartialValuesForResidualEquation(AccelerationSoil_NPlus_1MinusAlphaM_1, AccelerationSoil, AccelerationSoil_NPlus1, &
          !                                                  LumpedMassDry_NPlus_1MinusAlphaM_1, LumpedMassDry, LumpedMassDry_NPlus1, &
          !                                                  ExtLoad_NPlus_1MinusAlphaF_1, ExtLoad, ExtLoad_NPlus1, &
          !                                                  IntLoad_NPlus_1MinusAlphaF_1, IntLoad, IntLoad_NPlus1)
          
          
          
          
          
          ! Step 8: Evaluate residual magnitude 
          !call EvaluateResidualMagnitude(Residual_GeneralizedAlpha, Residual_Implicit_Magnitude) 
          
          ! This tolerance need to be dynamic 
          TOLERANCE = 1e-10
          
          if (Residual_Implicit_Magnitude < TOLERANCE) then 
              
              RETURN ! terminate the iteration 
              
          end if
          
          
          ! Step 9: Algorithm 4 regarding the acceleration increment 
          ! This is SS_0 in Algorithm 4 in Sulsky and Kaul (2004)
          ! Step 1 a) 
          !call SetAccelerationVectorValueToZero(DeltaAcceleration) ! not sure what the initial projection should be ! This subroutine is not necessary
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!! STEP 1 GMRES
          ! Step 1 b) 
          call SetResidual_Zero_To_NegativeResidualValue(Residual_GeneralizedAlpha, Residual_Zero)
          
          
          ! Step 1 c) 
          call EvaluateResidualMagnitude(Residual_Zero, Residual_Zero_Magnitude) 
          Rho_Zero = Residual_Zero_Magnitude**2
          
          ! Step 1 d) set kk = 1
          Iteration_kk = 1
          
          ! Step 1 e) set p_1 
          call Calculate_PP_1(Residual_Zero, Residual_Zero_Magnitude, PP_1)
          
          ! Step 1 f) 
          call Calculate_Beta(Residual_Zero_Magnitude, Beta)
          
          ! Step 1 g)
          
          
          call EvaluateErrorToleranceGMRES 
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!! STEP 2 GMRES          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Repeat Step 2 until the error tolerance is met 


                    
          
          !call GMRES( AA, bb, xx, threshold, nn_length_AA, mm_max_iteration) ![x, e] = gmres( A, b, x, max_iterations, threshold); 

          nn_length_AA = Counters%N
          
          
          call GMRES_DirectionalDerivative( &
                      LumpedMassDry, LumpedMassDry_NPlus1, &
                      TotalVelocitySoil, TotalVelocitySoil_NPlus1, &
                      IntLoad, IntLoad_NPlus1, &
                      Residual_GeneralizedAlpha, &
                      DeltaVelocity, &
                      Threshold, &
                      nn_Length_AA, &
                      mm_max_iteration )
                      !Residual_GeneralizedAlpha) 
          !TrialSolutionValue, & ! this is s_0 in SK04
          
          
          !call EvaluateResidualValueVelocity(TotalVelocitySoil, LumpedMassDry_NPlus1, &
          !                                         
          !                                         TotalVelocitySoilPrevious, LumpedMassDry, &
          !                                         Residual_GeneralizedAlpha)
          
          
          
          !AA, bb, xx, threshold, nn_length_AA, mm_max_iteration) ![x, e] = gmres( A, b, x, max_iterations, threshold); 

          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
        
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          !if (TOLERANCE_EVALUATION < TOLERANCE_GMRES) then 
          !
          !    ImplicitConverged = .TRUE.
          !    
          !else 
          !    
          !    ImplicitConverged = .FALSE. 
          !    
          !end if 
          
          ! Step 2 a) evaluate the directional derivative using acceleration as an input 
          !call EvaluateDirectionalDerivativeUsingAccelerationAsAnInput(Residual_GeneralizedAlpha, DeltaAcceleration, AccelerationSoil, DirectionalDerivative)
          !call EvaluateDirectionalDerivativeUsingAccelerationAsAnInput(Residual_GeneralizedAlpha, PP_1, AccelerationSoil, DirectionalDerivative)
          
          !subroutine EvaluateDirectionalDerivativeUsingAccelerationAsAnInput(Old_Residual, DeltaAcceleration, AccelerationSoil)
          !real(REAL_TYPE), dimension(Counters%N), intent(out) :: DirectionalDerivative 
          
          ! PP_2 is equal to the directional derivative 
          !PP_2 = DirectionalDerivative
          
          !h_k_1_i_1 = 0.0
          
          ! Step 2 b) for i = 1,2,...,k set h_k,i = p_k+1 * p_i; p_k+1 = p_k+1 - h_k,i * p_i
          do ii = 1,Iteration_kk ! need to change this to a do while loop 
             
              !do jj = 1, Counters%N
             
                  !h_k_1_i_1 = h_k_1_i_1 +  PP_2(jj)*PP_1(jj) ! This variable needs to be an allocatable as we need to sift through it 
        
                  !h_k_1_i_1(ii) = h_k_1_i_1(ii) +  PP_2(jj)*PP_1(jj) ! This variable needs to be an allocatable as we need to sift through it 
                  !h_kk           !h_kk
                  
                  ! This h_kk is a scalar in the first iteration 
                  
                  
                  ! ii is fine 
                  ! jj needs to be modified, along with PP values 
                  
                  
              !end do
              
          end do 
          
          ! This needs to occur within the loop iteration 
          !PP_2 = PP_2 - (h_k_1_i_1 * PP_1)

          
          
          ! Step 2 c) find hh_k+1,k
          !EvaluateResidualMagnitude(Residual_Implicit, Residual_Implicit_Magnitude) 
          !call EvaluateResidualMagnitude(PP_2, h_kPlus1_2_i_1)    
          
          ! Step 2 d) PP_2 = PP_2 * (1/h_k_2_i_1)
          !PP_2 = PP_2 * (1/h_kPlus1_2_i_1)
          
          !!!!!!!!!!!!!!!!! NEED TO FIX THE GIVENS ROTATION PART
          ! Step 2 e) Apply all previous Givens rotation to (h_1,k , ..., h_k+1,k)
          
          
          ! Step 2 f) Construct next Givens rotation on h_k,k and h_k+1,k to eliminate 
          
          
          ! Step 2 g) Apply new Givens rotation to bold{Beta}, Beta_k+1 is the new norm estimate of the residual 
          !!!!!!!!!!!!!!!!! NEED TO FIX THE GIVENS ROTATION PART
          
          
          ! Step 2 h) replace k with k+1 
          Iteration_kk = Iteration_kk + 1
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Repeat Step 2 until the error tolerance is met 
          
          ! Step 3 Solve the system Hy = Beta --
          
          
          
          
          
          
          
          
          
          !subroutine EvaluateNewResidual(AccelerationSoil, DeltaAcceleration, HH, & ! Input 
          !                                     New_Residual) !Output
          !
          !implicit none
          !
          !
          !
          !
          !
          !
          !
          !
          !end subroutine EvaluateNewResidual
          
          
          
          
          !call SetRho_Zero_To_MagnitudeSquare()
          
          
          
          
           != 0 ! set to zero 
          
          
          !!
          !! CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m
          !! CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f
          !! CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha
          !! CalParams%ImplicitGeneralizedAlphaTimeScheme%Delta
          !! 
          !
          !
          !
          !!!!!!!!!!!!!!!!! STEP 1: MAP FROM MP TO NODES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!  All the variables are initially at the material points after a time step. So first this here to map acceleration 
          !
          !
          !
          !! the internal forces for N+1 need to be initialized properly 
          !! We updated IntLoad_NPlus1 here 
          !
          !
          !
          !
          !
          !
          !!call MapAccelerationFromNodesToMPs()
          !!call MapVelocityFromNodesToMPs()
          !!call MapDisplacementFromNodesToMPs()
          !
          !!AccelerationSoil ! size = number of CPs x NDIM
          !!TotalVelocitySoil ! size = number of CPs x NDIM
          !!TotalDisplacementSoil ! size = number of CPs x NDIM
          !
          !!  Get mass soil at the nodes 
          !!Mass_I_N = LumpedMassDry ! size = number of CPs x NDIM
          !
          !!  Get kinematic variables at the nodes 
          !!  Here I relabel them so that it would make more sense
          !!Acceleration_I_N = AccelerationSoil ! size = number of CPs x NDIM
          !!Velocity_I_N = TotalVelocitySoil ! size = number of CPs x NDIM
          !!Displacement_I_N = TotalDisplacementSoil ! size = number of CPs x NDIM
          !
          !
          !
          !
          !
          !
          !!!!!!!!!!!!!!!!! STEP 2: PREDICTOR STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!  Predictor stage evaluation of the kinematic vector variables 
          !!ImplicitIterations = 0 
          !!Velocity_I_NPlus1 = Velocity_I_N ! size = number of CPs x NDIM
          !TotalVelocitySoil_NPlus1 = TotalVelocitySoil
          !
          !!Acceleration_I_NPlus1 = ( (Gamma_GAlpha-1)/Gamma_GAlpha ) * Acceleration_I_N ! size = number of CPs x NDIM
          !AccelerationSoil_NPlus1 = ( (CalParams%ImplicitGeneralizedAlphaTimeScheme%Delta-1)/CalParams%ImplicitGeneralizedAlphaTimeScheme%Delta ) * &
          !                                                                                              AccelerationSoil ! size = number of CPs x NDIM
          !
          !
          !! There is a discrepancy between the sizes of TotalDisplacementSoil and TotalVelocitySoil and AccelerationSoil
          !TotalDisplacementSoil_NPlus1 = TotalDisplacementSoil + &
          !    ( CalParams%TimeIncrement * TotalVelocitySoil(:,1)  ) + &
          !    ( 0.5 * (CalParams%TimeIncrement**2) * ( ((1 - 2*CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha)*AccelerationSoil(:,1)  ) + &
          !    (2*CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha*AccelerationSoil_NPlus1(:,1) ) ) )
          !
          !
          !
          !!IncrementalDisplacementSoil = TotalDisplacementSoil_NPlus1 - TotalDisplacementSoil
          !
          !
          !! We really want convective phase here with NPlus1 Values 
          !
          !
          !! size = number of CPs x NDIM
          !
          !!  Update positions of the material points based on the new values
          !!
          !!call UpdateParticleDisplacementImplicitGeneralizedAlpha(GlobPos_I_NPlus1)
          !!
          !!call UpdateParticleStrains()
          !!!
          !!call MPMDYNGetSig() !call UpdateParticleStress
          !!
          !!! Update particle weights
          !!if (IsMPMComputation()) then
          !!    call DynUpdateParticleWeights( )
          !!end if
          !!
          !!!Update porosity (Directly depends on: DEpsVol)
          !!if (CalParams%ApplyPorosityUpdate) then !Update porosity (Directly depends on: DEpsVol)
          !!    call DynUpdateParticlePorosity( )
          !!end if
          !!
          !!! Update particle global positions, particle-element assignment, particle local positions and shape function values
          !!! NOTE: EleParticles and Particle%ElementID get already updated while
          !!! the remaining house-keeping data is updated in the succeeding routine!!
          !!if (IsMPMComputation()) then ! MPM
          !!    call UpdateParticlePos()
          !!end if
          !!
          !
          !! Obtain "updated" momentum and lumped mass values  
          !call MapMomentumAndMassP2N(Momentum_NPlus1, LumpedMassDry_NPlus1) ! we need the "N+1" mass from here 
          !
          !! Note that AccelerationSoil is already on the "N+1" updated value 
          !! AccelerationSoil_NPlus1 = AccelerationSoil 
          !! We updated Momentum_NPlus1, LumpedMassDry_NPlus1
          !
          !
          !
          !! we need to recalculate the internal and external forces 
          !
          !
          !! Train of thought --> we have updated stresses --> and now we need to obtain internal forces 
          !         
          !! I got this from the Lagrangian Phase 
          !!call GetNodalExtAndIntForces_NPlus1() ! rotated to local coordinate system
          !
          !
          !!
          !!call UpdateParticlePositionsImplicitGeneralizedAlpha(GlobPos_I_NPlus1)
          !!
          !!call UpdateNodalMassImplicitGeneralizedAlpha(Mass_I_NPlus1)  ! --> Based on the new positions we can calculate Mass_I_NPlus1
          !!
          !!call UpdateNodalInternalForceImplicitGeneralizedAlpha(ForceInt_I_NPlus1) ! -> using new particle stress
          !!
          !!call UpdateNodalExternalForceImplicitGeneralizedAlpha(ForceExt_I_NPlus1)
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!!!!!!!!!!!!!!!!! STEP 3: CHECK CONVERGENCE BASED ON RESIDUAL EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!
          !!!  Calculate residual of the equation 
          !!
          !!! plug in new mass values, new acceleration values and evaluate residual Residual_Implicit
          !!
          !call CalculateGeneralizedAlphaResidualEquation()!Residual_Implicit)
          !!
          !!! Check residual based on tolerance --> need to be smaller than 10**-12 
          !call CheckResidualTolerance()!ImplicitConverged, Tolerance_ImplicitConvergence)
          !!
          !!
          !!if (ImplicitConverged == .TRUE.) then 
          !!  
          !!    ! if converged then exit here 
          !!    RETURN 
          !!  
          !!end if 
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!!!!!!!!!!!!!!!!! STEP 4: FINDING INCREMENT STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! In this stage we solve for acceleration increment that would achieve the tolerance level 
          !!
          !!
          !!if (ImplicitConverged == .NO.) then 
          !!
          !!    call GMRES_GeneralizedAlphaCalculation(Tolerance_ImplicitConvergence, IncrementAcceleration)
          !!
          !!    ! --> output here is IncrementAcceleration
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!!!!!!!!!!!!!!!!! STEP 5: CORRECTOR STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! In this stage we use IncrementAcceleration to find N+1 values 
          !!
          !!
          !!Acceleration_I_NPlus1 = AccelerationSoil + IncrementAcceleration
          !!Velocity_I_NPlus1 = TotalVelocitySoil + ( IncrementAcceleration*Gamma_GAlpha*CalParams%TimeIncrement )
          !!Displacement_I_NPlus1 = TotalDisplacementSoil + ( IncrementAcceleration*Beta_GAlpha* (CalParams%TimeIncrement)**2 )
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!
          !!end if 
          !!
          !!
          !!
          !!          
          !!!!!!!!!!!!!!!!!! STEP 6: ERROR TESTING STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!
          !!Factor = Residual_PreviousIteration/Residual_CurrentIteration
          !!
          !!if (Factor < Tolerance_ImplicitConvergence) then 
          !!    
          !!    ImplicitConverged = .YES. 
          !!    
          !!end if 
          !!
          !!
          !!
          !!
          !
          
          
          
          
          
          
          
          
          
          
          
          
          !mm_max_iteration = 100 ! user set limit 
          !
          !nn_length_AA = 5!2
          !!mm_max_iteration = max_iteration
          !
          !! Define example of AA 
          !! Attempt 1 
          !!AA(1,1) = 1
          !!AA(1,2) = 20
          !!AA(2,1) = 5
          !!AA(2,2) = 7
          !!real(REAL_TYPE), allocatable, dimension(:,:) :: AA
          !!real(REAL_TYPE), allocatable, dimension(:) :: bb
          !!real(REAL_TYPE), allocatable, dimension(:) :: xx
          !
          !allocate(AA(nn_length_AA, nn_length_AA))
          !allocate(bb(nn_length_AA))
          !allocate(xx(nn_length_AA))
          !
          !
          !
          !! Attempt 2 
          !AA(1,1) = 10
          !AA(1,2) = 50
          !AA(1,3) = 20
          !AA(1,4) = 5
          !AA(1,5) = 55
          !AA(2,1) = 5
          !AA(2,2) = 9
          !AA(2,3) = 1
          !AA(2,4) = 20
          !AA(2,5) = 7
          !AA(3,1) = 6
          !AA(3,2) = 9
          !AA(3,3) = 3
          !AA(3,4) = 2
          !AA(3,5) = 1
          !AA(4,1) = 5
          !AA(4,2) = 6
          !AA(4,3) = 1
          !AA(4,4) = 4
          !AA(4,5) = 9
          !AA(5,1) = 4
          !AA(5,2) = 8
          !AA(5,3) = 4
          !AA(5,4) = 7
          !AA(5,5) = 6
          !
          !! Define example of bb 
          !bb(1) = 1
          !bb(2) = 2 
          !bb(3) = 9
          !bb(4) = 20
          !bb(5) = 4
          !
          !! Define threshold 
          !
          !! Define example of xx 
          !!xx(1) = 2
          !!xx(2) = 3
          !!xx(1) = 1
          !!xx(2) = 5
          !!xx(1) = 20
          !!xx(2) = 7
          !
          !xx(1) = 2
          !xx(2) = 0 !3
          !xx(3) = 5 
          !xx(4) = 6
          !xx(5) = 4
          !
          !! call GMRES algorithm 
          !
          !
          !! we need to know what needs to be input and what is output from the GMRES 
          !! Input: 
          !! 1) we need to evaluate p_k+1 using the directional derivative --> vector pp_k+1 which is q
          !! 2) we need to input the residual into the GMRES subroutine 
          !
          !! Output: 
          !! 1) s (i.e., velocity increment)
          !! 2) that's it.... 
          !call GMRES( AA, bb, xx, threshold, nn_length_AA, mm_max_iteration) ![x, e] = gmres( A, b, x, max_iterations, threshold); 
          
          
          
          
          
          
          
          
          
          
    
        end subroutine DynamicImplicitGeneralizedAlphaScheme_

        
        
        
        !%---------------------------------------------------------------------%
        !%                                GMRES                                %
        !%---------------------------------------------------------------------%
        subroutine GMRES_DirectionalDerivative( & !AA, bb, xx, threshold, nn_length_AA, mm_max_iteration)
                      LumpedMassDry, LumpedMassDry_NPlus1, &
                      TotalVelocitySoil, TotalVelocitySoil_NPlus1, &
                      IntLoad, IntLoad_NPlus1, &
                      Residual_GeneralizedAlpha, &
                      DeltaVelocity, &
                      Threshold, &
                      nn_Length_AA, &
                      mm_max_iteration)!, & ! this is s_0 in SK04
                      !Residual_GeneralizedAlpha) 
          implicit none 
          
          
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus1
          
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil ! this is the initial velocity          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NPlus1 ! this is the initial velocity

          
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad_NPlus1
          
          
          ! need to generalize using Counters%NEntity if we were to use this contact algorithm
          !real(REAL_TYPE), dimension(Counters%N), intent(inout) :: TrialSolutionValue !this is s_0 in SK04
          
          
          
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: InternalForceIncrement
          
          
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: Residual_GeneralizedAlpha

          !Residual_GeneralizedAlpha
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NPlus1
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          !real(REAL_TYPE), dimension(nn_length_AA,nn_length_AA), intent(in) :: AA
          !real(REAL_TYPE), dimension(nn_length_AA), intent(in) :: bb 
          real(REAL_TYPE), dimension(nn_length_AA), intent(inout) :: DeltaVelocity
          
          real(REAL_TYPE), intent(in) :: threshold
          integer(INTEGER_TYPE), intent(in) :: nn_length_AA

          
          integer(INTEGER_TYPE), intent(in) :: mm_max_iteration 
          
          integer(INTEGER_TYPE) :: ii, jj, kk, nn
          
          real(REAL_TYPE), dimension(nn_length_AA)  :: rr !           intent(in)?
          
          real(REAL_TYPE) :: bb_norm 
          real(REAL_TYPE) :: rr_norm 
          real(REAL_TYPE) :: error 
          real(REAL_TYPE) :: ee
          
          real(REAL_TYPE), allocatable, dimension(:) :: sn ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: cs ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: e1 ! vectors?
          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ ! Rank 2 matrix?
          
          real(REAL_TYPE), allocatable, dimension(:) :: beta ! vectors?
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Arnoldi_Full

          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ_Arnoldi ! vectors? --> this needs to be assembled in QQ
          real(REAL_TYPE), allocatable, dimension(:) :: HH_Arnoldi ! vectors?
          
          !real(REAL_TYPE), allocatable, dimension(:) :: QQ_Arnoldi_New ! vectors? --> this needs to be assembled in QQ
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Arnoldi_New ! vectors?

          
          real(REAL_TYPE), allocatable, dimension(:) :: ee_saved ! error saved 
          
          
          
          real(REAL_TYPE), allocatable, dimension(:, :) :: QQ_Full
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Full
          
          real(REAL_TYPE), allocatable, dimension(:, :) :: QQ_Current
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Current
          
          real(REAL_TYPE), allocatable, dimension(:, :) :: QQ_Stored
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Stored
          
          real(REAL_TYPE), allocatable, dimension(:) :: YY
          
          integer(INTEGER_TYPE) :: INFO
          !nn_length_AA --> nn in Matlab code
          !mm_max_iteration --> mm in Matlab code
          
          ! use xx as the initial vector 
          
          !do ii = 1,2
          !    do jj = 1,2
          !        rr(ii,1) = bb(ii,1) - (AA * xx)
          !    end do 
          !end do
          
          !nn_length_AA = Counters%N
          
          rr = 0.0
          
          !do ii = 1,nn_length_AA!2
          !   do jj = 1,nn_length_AA!2
          !       
          !       rr(ii) = rr(ii) + ( bb(ii) - ( AA(ii,jj)*xx(jj) ) )        
          !       
          !   end do 
          !end do 
          
          ! the residual is already evaluated so we do not need to evaluate this again here 
          rr = Residual_GeneralizedAlpha
          !do ii = 1,nn_length_AA!2
          !   do jj = 1,nn_length_AA!2
          !       
          !       rr(ii) = rr(ii) + ( ( AA(ii,jj)*xx(jj) ) )        
          !       
          !   end do 
          !end do 
          !
          !
          !!do ii = 1,nn_length_AA!2
          !!
          !!    rr(ii) = rr(ii) - bb(ii)
          !!    
          !!end do 
          !
          !!rr = rr - (nn_length_AA*bb)
          !
          !rr = bb - rr

          
          !
          !
          !rr(1,1) = bb(1,1) - (AA(1,1)*xx(1,1) + AA(1,2)*xx(2,1))
          !rr(2,1) = bb(2,1) - (AA(2,1)*xx(1,1) + AA(2,2)*xx(2,1))
          !    
          !    
          !rr(ii) = bb(ii) - (AA(ii,ii)*xx(ii) + AA(ii,ii+1)*xx(ii+1))
          !    
              
          
          
          ! is this going to work? if not, then you need to do 2 for loops 
          
          !bb_norm = 1!norm2(bb, 1) 
          ! using norm2 function from FOTRAN90. Inputing rank as 1, but this is not necessary (i.e., optional)
          
          rr_norm = norm2(rr, 1)
          ! using norm2 function from FOTRAN90. Inputing rank as 1, but this is not necessary (i.e., optional)
          
          error = rr_norm!/ !bb_norm ! real scalar in the numerator and denominator
          
          ! initialize the 1D vectors
          allocate(sn(mm_max_iteration))
          allocate(cs(mm_max_iteration))
          allocate(e1(mm_max_iteration+1))
          allocate(beta(mm_max_iteration+1))
          allocate(QQ(nn_length_AA,1))
          !allocate(QQ_Full(nn,nn))
          
          cs = 0 
          sn = 0 
          e1 = 0 
          beta = 0
          QQ = 0 
          
          e1(1) = 1 ! assigning the first element of ee to be 1
          ee = error ! ee is a scalar. Set it equal to the value of error
          
          QQ(:,1) = (1/rr_norm) * rr ! assigning this for the first iteration 
          
          beta = rr_norm * e1 ! e1 and beta are the same size vectors 
          
          nn = 1
          kk = 1
          HH_Arnoldi = 0 
          QQ_Arnoldi = 0
          HH_Arnoldi_New = 0 
          !QQ_Arnoldi_New = 0
          
          allocate(QQ_Stored(nn_length_AA,kk))  !kk+1
          allocate(HH_Stored(kk+1,kk))
          
          QQ_Stored = 0
          HH_Stored = 0
          
          do kk = 1, mm_max_iteration ! loop accross the number of maximum iterations
              !do while (error .ge. threshold ) 
              ! run Arnoldi 
              
              if (allocated(QQ_Full)) then 
                  deallocate(QQ_Full)
              end if 
              
              if (allocated(HH_Full)) then 
                  deallocate(HH_Full)
              end if
              
              allocate(QQ_Full(nn_length_AA,kk+1)) !2
              allocate(HH_Full(kk+1,kk))
              
              
              QQ_Full = 0
              HH_Full = 0
              !QQ_Full(:,1) = QQ(:,1)
              QQ_Stored(:,1) = QQ(:,1) ! this is a bit redundant in steps above k=1
              QQ_Full(:,1:kk) = QQ_Stored(:,1:kk) ! QQ_Full in column kk+1 should be left empty
              if (kk .ge. 2) then 
              HH_Full(1:kk,1:kk-1) = HH_Stored(:,:)
              end if
              
              !PP_1 = QQ(:,1)
              !call Arnoldi_DirectionalDerivative(AA, QQ_Full, kk, HH_Arnoldi_New, QQ_Arnoldi, nn_length_AA)!, HH_Arnoldi_New, QQ_Arnoldi_New)
              
              call Arnoldi_DirectionalDerivative(LumpedMassDry, LumpedMassDry_NPlus1, &
                                                 TotalVelocitySoil, TotalVelocitySoil_NPlus1, &
                                                 IntLoad, IntLoad_NPlus1, & 
                                                 QQ, kk, HH_Arnoldi_New, QQ_Arnoldi, nn_length_AA, &
                                                 rr)!, &
              
              nn = nn + 1
              
              !HH_Full (1:kk+1,kk) = HH_Arnoldi_New
              !QQ_Full = QQ
              
              
              
              
              
              
              
              !QQ_Full(:,1) = QQ(:,1)
              !QQ_Full(:,2) = QQ_Arnoldi(:,kk)
              
              !deallocate(QQ)
              !allocate(QQ(nn,nn))
              !
              !!QQ(:,1) = QQ_Full(:,2) 
              !QQ = QQ_Full
              
              
              !We need to append QQ and QQ_Arnoldi 
              ! QQ_Full here is same as Q in matlab
              ! H here is same as HH_Arnoldi_New
              
          !    !call Arnoldi(A, Q, kk, HH_Arnoldi(1:k+1, k), QQ_Arnoldi(:, k+1))
          !    ! the output from Arnoldi should be a vectors --> HH and QQ
          !    HH(1:kk+1,kk) = HH_Arnoldi(:)
          !    QQ(1:kk+1,kk) = QQ_Arnoldi(:)
          !    
          !    ! assign values output from the Arnoldi 
          !    
          !    ! eliminate the last element in H ith row and update the rotation matrix 
          !    ![H(1:k+1, k), cs(k), sn(k)] = apply_givens_rotation(H(1:k+1,k), cs, sn, k);
          !    call apply_givens_rotation(HH(1:kk+1,kk), cs, sn, kk)

              ! eliminate the last element in H ith row and update the rotation matrix 
              call apply_givens_rotation( HH_Arnoldi_New, cs, sn, kk, mm_max_iteration )
              
              if (allocated(HH_Arnoldi_Full)) then 
                  
                  deallocate(HH_Arnoldi_Full)
                  
              end if 
              
              allocate(HH_Arnoldi_Full(nn,kk))
              
              HH_Arnoldi_Full = 0 
              
              HH_Arnoldi_Full = HH_Arnoldi_New
              
              ! update the residual vector 
              beta(kk+1) = -sn(kk) * beta(kk)
              beta(kk) = cs(kk) * beta(kk)
              
              error = abs(beta(kk+1)) !/ bb_norm
              
              !allocate( ee_saved(nn) )
              ! save the error 
              
              ! we need to loop to save the error in the error_saved array 
              !ee_saved(kk) = ee
              !ee_saved(kk+1) = error
              
              
              ! we need to save what we have so far from H and Q before we deallocate ok? 
              !allocate(QQ_Current(2,1))
              !allocate(HH_Current(kk+1,1))
              
              

              deallocate(QQ_Stored)  !kk+1
              deallocate(HH_Stored)
              
              allocate(QQ_Stored(nn_length_AA,kk+1))  !kk+1
              allocate(HH_Stored(kk+1,kk))
              
              QQ_Stored = QQ_Full
              
              if (kk .ge. 2) then 
              HH_Stored = 0
              HH_Stored(:, 1:kk-1) = HH_Full(:, 1:kk-1)
              HH_Stored(:, kk) = HH_Arnoldi_New(:,1)
              
              else 
                  
                  HH_Stored = 0          
                  HH_Stored = HH_Arnoldi_New

              end if 
              !if (kk>2) then 
              !
              !    HH_Stored(1:kk-1,  ) = 
              !    
              !end if 
              
              
              !HH_Full(1:kk,) 
              if (error .le. threshold ) then 
                  
                  exit
                  
              end if 
              
                  
                  
              
              
              
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
          end do 
          
          
          call dtrtri( 'U', 'N', nn_length_AA, HH_Stored(1:kk, 1:kk), nn_length_AA, INFO  )
          
          !end do 
          
          ! There is probably more stuff here....
              
          ! I need to multiply the HH_Stored(1:kk, 1:kk) and beta(1:kk)
          
          !HH_Stored(1:kk, 1:kk) * beta(1:kk)
          
          allocate(YY(kk)) 
          YY = 0 
          
          ! conduct the multiplication to find YY 
          ! loop over rows  
          do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix
              
              ! loop over columns 
              do jj = 1, kk 
                  YY(ii) = YY(ii) + HH_Stored(ii, jj) * beta(jj)  
              end do 
              
          end do 
          
          !xx = xx + (QQ_Full(:,1:kk) * YY)
          ! loop over rows 
          do ii = 1, nn_length_AA ! rank of matrix 
              do jj = 1, kk
                  DeltaVelocity(ii) = DeltaVelocity(ii) + QQ_Full(ii,jj)*YY(jj)
              end do 
          end do 
          
          
          
        end subroutine GMRES_DirectionalDerivative
        
        
        
        
        
        
        
        
        
        
        
        subroutine EvaluateResidualValueVelocity(TotalVelocitySoil_NPlus1, LumpedMassDry_NPlus1, &
                                                   IntLoad_NPlus1, IntLoad, &
                                                   TotalVelocitySoil, LumpedMassDry, &
                                                   Residual)
          
          
          !AccelerationSoil_NPlus_1MinusAlphaM_1, LumpedMassDry_NPlus_1MinusAlphaM_1, ExtLoad_NPlus_1MinusAlphaF_1, IntLoad_NPlus_1MinusAlphaF_1, &
          !                                 Residual_GeneralizedAlpha)
          
          implicit none 
          
          
          

          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil ! this is the initial velocity
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus1
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad_NPlus1
          
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: InternalForceIncrement
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NPlus1
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry
          
          real(REAL_TYPE) :: Implicit_mu
          
          !---------------------------------------------------

          
          real(REAL_TYPE), dimension(Counters%N), intent(out) :: Residual
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_NPlus_1MinusAlphaM_1
          !
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus_1MinusAlphaM_1
          !
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad_NPlus_1MinusAlphaF_1
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: ExtLoad_NPlus_1MinusAlphaF_1
          !
          !! local 
          integer(INTEGER_TYPE) :: ii
          !real(REAL_TYPE) :: Implicit_mu
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NN ! this is the initial velocity
          !
          !
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: InternalForceIncrement
          
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: LumpedMassDry
          
          Implicit_mu = 1 ! gauge set to fully implicit
          
          ! we need to update IntLoad_NPlus1 after every increment
          InternalForceIncrement = IntLoad_NPlus1 - IntLoad 
          !                        Implicit - Explicit
          do ii = 1, Counters%N
          
            !Residual_GeneralizedAlpha(ii) = AccelerationSoil_NPlus_1MinusAlphaM_1(ii, 1)*LumpedMassDry_NPlus_1MinusAlphaM_1(ii, 1) - ExtLoad_NPlus_1MinusAlphaF_1(ii, 1) + IntLoad_NPlus_1MinusAlphaF_1(ii, 1)
            ! Entity = 1 
            
            !InternalForceIncrement = 
            
            ! Implicit terms - explicit terms --> lead to this expression 
            Residual(ii) = (TotalVelocitySoil_NPlus1(ii, 1)*LumpedMassDry_NPlus1(ii, 1)) &
                         - (Implicit_mu*CalParams%TimeIncrement*InternalForceIncrement(ii, 1)) &
                         - (TotalVelocitySoil(ii, 1)*LumpedMassDry(ii, 1))
                                                    ! M*v^imp                                 -                 mu*DeltaT*Delta_F                                  -              Mv^exp
            
            ! I need to use the and the new internal force to be used here
            
            
            !- ExtLoad_NPlus_1MinusAlphaF_1(ii, 1) + IntLoad_NPlus_1MinusAlphaF_1(ii, 1)
          
            !TotalVelocitySoilPrevious
            
            
          end do 
          
                                           
        end subroutine EvaluateResidualValueVelocity
                                                   
        !-------------------------------------------------------------------------------------------
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !%---------------------------------------------------------------------%
        !%                                GMRES                                %
        !%---------------------------------------------------------------------%
        subroutine GMRES( AA, bb, xx, threshold, nn_length_AA, mm_max_iteration)
          
          implicit none 
          
          real(REAL_TYPE), dimension(nn_length_AA,nn_length_AA), intent(in) :: AA
          real(REAL_TYPE), dimension(nn_length_AA), intent(in) :: bb 
          real(REAL_TYPE), dimension(nn_length_AA), intent(inout) :: xx
          
          real(REAL_TYPE), intent(in) :: threshold
          integer(INTEGER_TYPE), intent(in) :: nn_length_AA

          
          integer(INTEGER_TYPE), intent(in) :: mm_max_iteration 
          
          integer(INTEGER_TYPE) :: ii, jj, kk, nn
          
          real(REAL_TYPE), dimension(nn_length_AA)  :: rr !           intent(in)?
          
          real(REAL_TYPE) :: bb_norm 
          real(REAL_TYPE) :: rr_norm 
          real(REAL_TYPE) :: error 
          real(REAL_TYPE) :: ee
          
          real(REAL_TYPE), allocatable, dimension(:) :: sn ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: cs ! vectors?
          real(REAL_TYPE), allocatable, dimension(:) :: e1 ! vectors?
          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ ! Rank 2 matrix?
          
          real(REAL_TYPE), allocatable, dimension(:) :: beta ! vectors?
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Arnoldi_Full

          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ_Arnoldi ! vectors? --> this needs to be assembled in QQ
          real(REAL_TYPE), allocatable, dimension(:) :: HH_Arnoldi ! vectors?
          
          !real(REAL_TYPE), allocatable, dimension(:) :: QQ_Arnoldi_New ! vectors? --> this needs to be assembled in QQ
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Arnoldi_New ! vectors?

          
          real(REAL_TYPE), allocatable, dimension(:) :: ee_saved ! error saved 
          
          
          
          real(REAL_TYPE), allocatable, dimension(:, :) :: QQ_Full
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Full
          
          real(REAL_TYPE), allocatable, dimension(:, :) :: QQ_Current
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Current
          
          real(REAL_TYPE), allocatable, dimension(:, :) :: QQ_Stored
          real(REAL_TYPE), allocatable, dimension(:, :) :: HH_Stored
          
          real(REAL_TYPE), allocatable, dimension(:) :: YY
          
          integer(INTEGER_TYPE) :: INFO
          !nn_length_AA --> nn in Matlab code
          !mm_max_iteration --> mm in Matlab code
          
          ! use xx as the initial vector 
          
          !do ii = 1,2
          !    do jj = 1,2
          !        rr(ii,1) = bb(ii,1) - (AA * xx)
          !    end do 
          !end do
          
          rr = 0.0
          
          !do ii = 1,nn_length_AA!2
          !   do jj = 1,nn_length_AA!2
          !       
          !       rr(ii) = rr(ii) + ( bb(ii) - ( AA(ii,jj)*xx(jj) ) )        
          !       
          !   end do 
          !end do 
          
          do ii = 1,nn_length_AA!2
             do jj = 1,nn_length_AA!2
                 
                 rr(ii) = rr(ii) + ( ( AA(ii,jj)*xx(jj) ) )        
                 
             end do 
          end do 
          
          
          !do ii = 1,nn_length_AA!2
          !
          !    rr(ii) = rr(ii) - bb(ii)
          !    
          !end do 
          
          !rr = rr - (nn_length_AA*bb)

          rr = bb - rr

          
          !
          !
          !rr(1,1) = bb(1,1) - (AA(1,1)*xx(1,1) + AA(1,2)*xx(2,1))
          !rr(2,1) = bb(2,1) - (AA(2,1)*xx(1,1) + AA(2,2)*xx(2,1))
          !    
          !    
          !rr(ii) = bb(ii) - (AA(ii,ii)*xx(ii) + AA(ii,ii+1)*xx(ii+1))
          !    
              
          
          
          ! is this going to work? if not, then you need to do 2 for loops 
          
          bb_norm = norm2(bb, 1) 
          ! using norm2 function from FOTRAN90. Inputing rank as 1, but this is not necessary (i.e., optional)
          
          rr_norm = norm2(rr, 1)
          ! using norm2 function from FOTRAN90. Inputing rank as 1, but this is not necessary (i.e., optional)
          
          error = rr_norm/ bb_norm ! real scalar in the numerator and denominator
          
          ! initialize the 1D vectors
          allocate(sn(mm_max_iteration))
          allocate(cs(mm_max_iteration))
          allocate(e1(mm_max_iteration+1))
          allocate(beta(mm_max_iteration+1))
          allocate(QQ(nn_length_AA,1))
          !allocate(QQ_Full(nn,nn))
          
          cs = 0 
          sn = 0 
          e1 = 0 
          beta = 0
          QQ = 0 
          
          e1(1) = 1 ! assigning the first element of ee to be 1
          ee = error ! ee is a scalar. Set it equal to the value of error
          
          QQ(:,1) = (1/rr_norm) * rr ! assigning this for the first iteration 
          
          beta = rr_norm * e1 ! e1 and beta are the same size vectors 
          
          nn = 1
          kk = 1
          HH_Arnoldi = 0 
          QQ_Arnoldi = 0
          HH_Arnoldi_New = 0 
          !QQ_Arnoldi_New = 0
          
          allocate(QQ_Stored(nn_length_AA,kk))  !kk+1
          allocate(HH_Stored(kk+1,kk))
          
          QQ_Stored = 0
          HH_Stored = 0
          
          do kk = 1, mm_max_iteration ! loop accross the number of maximum iterations
              !do while (error .ge. threshold ) 
              ! run Arnoldi 
              
              if (allocated(QQ_Full)) then 
                  deallocate(QQ_Full)
              end if 
              
              if (allocated(HH_Full)) then 
                  deallocate(HH_Full)
              end if
              
              allocate(QQ_Full(nn_length_AA,kk+1)) !2
              allocate(HH_Full(kk+1,kk))
              
              
              QQ_Full = 0
              HH_Full = 0
              !QQ_Full(:,1) = QQ(:,1)
              QQ_Stored(:,1) = QQ(:,1) ! this is a bit redundant in steps above k=1
              QQ_Full(:,1:kk) = QQ_Stored(:,1:kk) ! QQ_Full in column kk+1 should be left empty
              if (kk .ge. 2) then 
              HH_Full(1:kk,1:kk-1) = HH_Stored(:,:)
              end if
              call Arnoldi(AA, QQ_Full, kk, HH_Arnoldi_New, QQ_Arnoldi, nn_length_AA)!, HH_Arnoldi_New, QQ_Arnoldi_New)
              
              nn = nn + 1
              
              !HH_Full (1:kk+1,kk) = HH_Arnoldi_New
              !QQ_Full = QQ
              
              
              
              
              
              
              
              !QQ_Full(:,1) = QQ(:,1)
              !QQ_Full(:,2) = QQ_Arnoldi(:,kk)
              
              !deallocate(QQ)
              !allocate(QQ(nn,nn))
              !
              !!QQ(:,1) = QQ_Full(:,2) 
              !QQ = QQ_Full
              
              
              !We need to append QQ and QQ_Arnoldi 
              ! QQ_Full here is same as Q in matlab
              ! H here is same as HH_Arnoldi_New
              
          !    !call Arnoldi(A, Q, kk, HH_Arnoldi(1:k+1, k), QQ_Arnoldi(:, k+1))
          !    ! the output from Arnoldi should be a vectors --> HH and QQ
          !    HH(1:kk+1,kk) = HH_Arnoldi(:)
          !    QQ(1:kk+1,kk) = QQ_Arnoldi(:)
          !    
          !    ! assign values output from the Arnoldi 
          !    
          !    ! eliminate the last element in H ith row and update the rotation matrix 
          !    ![H(1:k+1, k), cs(k), sn(k)] = apply_givens_rotation(H(1:k+1,k), cs, sn, k);
          !    call apply_givens_rotation(HH(1:kk+1,kk), cs, sn, kk)

              ! eliminate the last element in H ith row and update the rotation matrix 
              call apply_givens_rotation( HH_Arnoldi_New, cs, sn, kk, mm_max_iteration )
              
              if (allocated(HH_Arnoldi_Full)) then 
                  
                  deallocate(HH_Arnoldi_Full)
                  
              end if 
              
              allocate(HH_Arnoldi_Full(nn,kk))
              
              HH_Arnoldi_Full = 0 
              
              HH_Arnoldi_Full = HH_Arnoldi_New
              
              ! update the residual vector 
              beta(kk+1) = -sn(kk) * beta(kk)
              beta(kk) = cs(kk) * beta(kk)
              
              error = abs(beta(kk+1)) / bb_norm
              
              !allocate( ee_saved(nn) )
              ! save the error 
              
              ! we need to loop to save the error in the error_saved array 
              !ee_saved(kk) = ee
              !ee_saved(kk+1) = error
              
              
              ! we need to save what we have so far from H and Q before we deallocate ok? 
              !allocate(QQ_Current(2,1))
              !allocate(HH_Current(kk+1,1))
              
              

              deallocate(QQ_Stored)  !kk+1
              deallocate(HH_Stored)
              
              allocate(QQ_Stored(nn_length_AA,kk+1))  !kk+1
              allocate(HH_Stored(kk+1,kk))
              
              QQ_Stored = QQ_Full
              
              if (kk .ge. 2) then 
              HH_Stored = 0
              HH_Stored(:, 1:kk-1) = HH_Full(:, 1:kk-1)
              HH_Stored(:, kk) = HH_Arnoldi_New(:,1)
              
              else 
                  
                  HH_Stored = 0          
                  HH_Stored = HH_Arnoldi_New

              end if 
              !if (kk>2) then 
              !
              !    HH_Stored(1:kk-1,  ) = 
              !    
              !end if 
              
              
              !HH_Full(1:kk,) 
              if (error .le. threshold ) then 
                  
                  exit
                  
              end if 
              
                  
                  
              
              
              
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
          end do 
          
          
          call dtrtri( 'U', 'N', nn_length_AA, HH_Stored(1:kk, 1:kk), nn_length_AA, INFO  )
          
          !end do 
          
          ! There is probably more stuff here....
              
          ! I need to multiply the HH_Stored(1:kk, 1:kk) and beta(1:kk)
          
          !HH_Stored(1:kk, 1:kk) * beta(1:kk)
          
          allocate(YY(kk)) 
          YY = 0 
          
          ! conduct the multiplication to find YY 
          ! loop over rows  
          do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix
              
              ! loop over columns 
              do jj = 1, kk 
                  YY(ii) = YY(ii) + HH_Stored(ii, jj) * beta(jj)  
              end do 
              
          end do 
          
          !xx = xx + (QQ_Full(:,1:kk) * YY)
          ! loop over rows 
          do ii = 1, nn_length_AA ! rank of matrix 
              do jj = 1, kk
                  xx(ii) = xx(ii) + QQ_Full(ii,jj)*YY(jj)
              end do 
          end do 
          
          
          
        end subroutine GMRES 
        
        

        !%----------------------------------------------------%
        !%                  Arnoldi Function                  %
        !%----------------------------------------------------%
        subroutine Arnoldi(AA, QQ, kk, HH_Arnoldi_New, QQ_Arnoldi, nn_length_AA)
        !% inputs 
        !% - AA is the matrix on the left hand side to be inverted 
        !% - QQ is initially the normalized residual vector 
        !% - kk is the number of iterations 
        !% outputs 
        !% - HH_Arnoldi is the step vector to be placed in H
        !% - QQ_Arnoldi is the new Q step vector used to expand Q 
        
        implicit none 
                  
        ! input 
        real(REAL_TYPE), dimension(nn_length_AA,nn_length_AA), intent(in) :: AA    
        real(REAL_TYPE), dimension(:, :), intent(inout) :: QQ ! Rank 2 matrix? !allocatable, 
        
        
        
        integer(INTEGER_TYPE), intent(in) :: kk
        integer(INTEGER_TYPE), intent(in) :: nn_length_AA
        
        ! output 
        !real(REAL_TYPE), dimension(kk), intent(inout) :: HH_Arnoldi 
        !real(REAL_TYPE), dimension(2), intent(inout) :: QQ_Arnoldi
        real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: QQ_Arnoldi ! vectors? --> this needs to be assembled in QQ
        real(REAL_TYPE), allocatable, dimension(:) :: HH_Arnoldi ! vectors?

        real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HH_Arnoldi_New ! vectors?
        
        ! local 
        integer(INTEGER_TYPE) :: ii, jj
        
        
        ! consider moving this outside
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !allocate(HH_Arnoldi(kk))
        allocate(HH_Arnoldi(kk)) ! 1 
        !allocate(QQ_Arnoldi(kk+1)) ! 2
        
        if (.not.allocated(QQ_Arnoldi)) then
            allocate(QQ_Arnoldi(nn_length_AA,1)) ! 2
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        QQ_Arnoldi = 0 
        HH_Arnoldi = 0 
        do ii = 1,nn_length_AA!2
            do jj = 1,nn_length_AA!2
                QQ_Arnoldi(ii,1) = QQ_Arnoldi(ii,1) + (AA(ii,jj) * QQ(jj,kk)) !% Krylov Vector
            end do
        end do 
        
        ! need to call directional derivative here 
        !call EvaluateDirectionalDerivativeUsingAccelerationAsAnInput(Old_Residual, DeltaAcceleration, AccelerationSoil, DirectionalDerivative)
        
        do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix
            
            do jj = 1,nn_length_AA !2 
                HH_Arnoldi(ii) = HH_Arnoldi(ii) + ( QQ_Arnoldi(jj,1) * QQ(jj,ii) )
            end do 
            
            !QQ_Arnoldi(:,1) = QQ_Arnoldi(:,1) - ( HH_Arnoldi(ii) * QQ(:,ii))
            
            !deallocate(HH_Arnoldi)
            !allocate(HH_Arnoldi(ii)
            
        end do 
        
        do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix
            
            QQ_Arnoldi(:,1) = QQ_Arnoldi(:,1) - ( HH_Arnoldi(ii) * QQ(:,ii))
            
            !if ( abs(QQ_Arnoldi(:,1)) < 1e-14) then
            !    
            !    QQ_Arnoldi(:,1)  = 0.0
            !
            !
            !
            !end if 
            
            
            
        end do 
        
        
        ! need to store these 
        ! this needs to be stored outside this loop within the entire module as a global variable 
        !if (kk>1) then 
        !allocate(HH_Arnoldi_Full (kk+1, kk) )
        !HH_Arnoldi_Full(1:kk, kk-1) = HH_Arnoldi_New
        !end if
        
        
        
        if (allocated(HH_Arnoldi_New)) then 
            deallocate(HH_Arnoldi_New)
        end if 
        
        !The whole point is to amalgamate HH_Arnoldi_New and HH_Arnoldi, keeping that last row zero!
        
        allocate(HH_Arnoldi_New (kk+1, 1) ) 
        
        ! HH_Arnoldi_New needs to be rank 2 and we need to store a column for each step 
        ! so HH_Arnoldi_New should be two columns 
        ! HH_Arnoldi_New should be 3 rows and 2 columns
        HH_Arnoldi_New = 0
        
        HH_Arnoldi_New(1:kk, 1) = HH_Arnoldi(1:kk) 
        
        HH_Arnoldi_New(kk+1, 1) =  norm2(QQ_Arnoldi(:,1),1)
        
        
        ! we need to algamate QQ too here to it is 3 columns for kk =2 
        
        !allocate()
        
        ! rank 2
        do ii = 1, nn_length_AA!2 !% Modified Gram-Schmidt, keeping the Hessenberg matrix            
            if ( abs(QQ_Arnoldi(ii,1)) < 1e-14) then
                
                QQ_Arnoldi(ii,1)  = 1e-15
            
            
            
            end if 
            
        end do 
        
        
        QQ_Arnoldi = QQ_Arnoldi * (1/HH_Arnoldi_New(kk+1, 1) ) ! divide by a scalar
        
        !H size is (number of iterations +1, number of iterations)
        !Q size is (2 i.e., rank of matrix or how many rows in matrix, number of iterations +1)
        !
        !I NEED TO UNDERSTAND THOSE SIZES!!!!!
        
        
        QQ(:,kk+1) = QQ_Arnoldi(:,1)
        
        
        end subroutine Arnoldi 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine Arnoldi_DirectionalDerivative(LumpedMassDry, LumpedMassDry_NPlus1, &
                                                 TotalVelocitySoil, TotalVelocitySoil_NPlus1, &
                                                 IntLoad, IntLoad_NPlus1, & 
                                                 QQ, kk, HH_Arnoldi_New, QQ_Arnoldi, nn_length_AA, &
                                                 Old_Residual)!, &
                                                 !TrialSolution ) !AA, 
        
        !call EvaluateResidualValueVelocity(TotalVelocitySoil_NPlus1, LumpedMassDry_NPlus1, &
        !                                           IntLoad_NPlus1, IntLoad, &
        !                                           TotalVelocitySoil, LumpedMassDry, &
        !                                           Residual_GeneralizedAlpha)
        
        
        
        
        
        !AA, QQ, kk, HH_Arnoldi_New, QQ_Arnoldi, nn_length_AA)
        
        
        !subroutine Arnoldi_DirectionalDerivative(AA, QQ, kk, HH_Arnoldi_New, QQ_Arnoldi, nn_length_AA)
        
        !% inputs 
        !% - AA is the matrix on the left hand side to be inverted 
        !% - QQ is initially the normalized residual vector 
        !% - kk is the number of iterations 
        !% outputs 
        !% - HH_Arnoldi is the step vector to be placed in H
        !% - QQ_Arnoldi is the new Q step vector used to expand Q 
        
        implicit none 
                  
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus1
  
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil ! this is the initial velocity
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NPlus1

        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: IntLoad
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: IntLoad_NPlus1
        
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: DeltaVelocity
        
        
        ! input 
        !real(REAL_TYPE), dimension(nn_length_AA,nn_length_AA), intent(in) :: AA    
        real(REAL_TYPE), dimension(:, :), intent(inout) :: QQ ! Rank 2 matrix? !allocatable, 
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Old_Residual
        
        integer(INTEGER_TYPE), intent(in) :: kk
        integer(INTEGER_TYPE), intent(in) :: nn_length_AA
        
        ! output 
        !real(REAL_TYPE), dimension(kk), intent(inout) :: HH_Arnoldi 
        !real(REAL_TYPE), dimension(2), intent(inout) :: QQ_Arnoldi
        real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: QQ_Arnoldi ! vectors? --> this needs to be assembled in QQ
        real(REAL_TYPE), allocatable, dimension(:) :: HH_Arnoldi ! vectors?

        real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HH_Arnoldi_New ! vectors?
        
        ! local 
        integer(INTEGER_TYPE) :: ii, jj
        
        
        
        
          
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad_NPlus1
          
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: InternalForceIncrement
          
          !Old_Residual
        
        ! consider moving this outside
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !allocate(HH_Arnoldi(kk))
        allocate(HH_Arnoldi(kk)) ! 1 
        !allocate(QQ_Arnoldi(kk+1)) ! 2
        
        if (.not.allocated(QQ_Arnoldi)) then
            allocate(QQ_Arnoldi(nn_length_AA,1)) ! 2
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        QQ_Arnoldi = 0 
        HH_Arnoldi = 0 
        
        ! we do not need to evaluate this here as we need to calculate the directional derivative here 
        
        
        
        !do ii = 1,nn_length_AA!2
        !    do jj = 1,nn_length_AA!2
        !        QQ_Arnoldi(ii,1) = QQ_Arnoldi(ii,1) + (AA(ii,jj) * QQ(jj,kk)) !% Krylov Vector --> this is p_k+1
        !    end do
        !end do 
        
        
        DeltaVelocity = TotalVelocitySoil_NPlus1 - TotalVelocitySoil
        
        
        
        ! need to call directional derivative here 
        call EvaluateDirectionalDerivativeUsingAccelerationAsAnInput(Old_Residual, & ! F(v) 
                                                                     DeltaVelocity, &
                                                                     TotalVelocitySoil, &
                                                                     QQ_Arnoldi(ii,1), & !DirectionalDerivative, &
                                                                     QQ(:,kk)) ! this is PP_1 
        ! EDIT: THIS SHOULD BE QQ(:,kk) --> It should depend on the timestep 
        
        !QQ_Arnoldi(ii,1) = DirectionalDerivative
        
        do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix
            
            do jj = 1,nn_length_AA !2 
                HH_Arnoldi(ii) = HH_Arnoldi(ii) + ( QQ_Arnoldi(jj,1) * QQ(jj,ii) )
            end do 
            
            !QQ_Arnoldi(:,1) = QQ_Arnoldi(:,1) - ( HH_Arnoldi(ii) * QQ(:,ii))
            
            !deallocate(HH_Arnoldi)
            !allocate(HH_Arnoldi(ii)
            
        end do 
        
        do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix
            
            QQ_Arnoldi(:,1) = QQ_Arnoldi(:,1) - ( HH_Arnoldi(ii) * QQ(:,ii))
            
            !if ( abs(QQ_Arnoldi(:,1)) < 1e-14) then
            !    
            !    QQ_Arnoldi(:,1)  = 0.0
            !
            !
            !
            !end if 
            
            
            
        end do 
        
        
        ! need to store these 
        ! this needs to be stored outside this loop within the entire module as a global variable 
        !if (kk>1) then 
        !allocate(HH_Arnoldi_Full (kk+1, kk) )
        !HH_Arnoldi_Full(1:kk, kk-1) = HH_Arnoldi_New
        !end if
        
        
        
        if (allocated(HH_Arnoldi_New)) then 
            deallocate(HH_Arnoldi_New)
        end if 
        
        !The whole point is to amalgamate HH_Arnoldi_New and HH_Arnoldi, keeping that last row zero!
        
        allocate(HH_Arnoldi_New (kk+1, 1) ) 
        
        ! HH_Arnoldi_New needs to be rank 2 and we need to store a column for each step 
        ! so HH_Arnoldi_New should be two columns 
        ! HH_Arnoldi_New should be 3 rows and 2 columns
        HH_Arnoldi_New = 0
        
        HH_Arnoldi_New(1:kk, 1) = HH_Arnoldi(1:kk) 
        
        HH_Arnoldi_New(kk+1, 1) =  norm2(QQ_Arnoldi(:,1),1)
        
        
        ! we need to algamate QQ too here to it is 3 columns for kk =2 
        
        !allocate()
        
        ! rank 2
        do ii = 1, nn_length_AA!2 !% Modified Gram-Schmidt, keeping the Hessenberg matrix            
            if ( abs(QQ_Arnoldi(ii,1)) < 1e-14) then
                
                QQ_Arnoldi(ii,1)  = 1e-15          
            
            end if 
            
        end do 
        
        
        QQ_Arnoldi = QQ_Arnoldi * (1/HH_Arnoldi_New(kk+1, 1) ) ! divide by a scalar
        
        !H size is (number of iterations +1, number of iterations)
        !Q size is (2 i.e., rank of matrix or how many rows in matrix, number of iterations +1)
        !
        !I NEED TO UNDERSTAND THOSE SIZES!!!!!
        
        
        QQ(:,kk+1) = QQ_Arnoldi(:,1)
        
        
        end subroutine Arnoldi_DirectionalDerivative 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

              
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine EvaluateDirectionalDerivativeUsingAccelerationAsAnInput(Old_Residual, DeltaVelocity, TotalVelocitySoil, DirectionalDerivative, PP)
          
          ! we don't need to set DeltaAcceleration to zero
          
          
          implicit none 
          
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: Old_Residual ! F(a)
          
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: DeltaVelocity ! acceleration increment
          
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: TotalVelocitySoil ! acceleration increment
          
          real(REAL_TYPE), dimension(Counters%N) :: New_Residual ! F(a + h )
          
          real(REAL_TYPE), dimension(Counters%N), intent(out) :: DirectionalDerivative 
          

          real(REAL_TYPE), dimension(Counters%N), intent(in) :: PP ! F(a)

          
          integer(INTEGER_TYPE) :: ii
          
          real(REAL_TYPE) :: HH
          
          real(REAL_TYPE) :: TotalVelocitySoil_Magnitude
          
          real(REAL_TYPE) :: DeltaVelocity_Magnitude
          
          real(REAL_TYPE), dimension(Counters%N) :: TotalVelocitySoil_NPlus1_DirectionalDerivative
          
          !DeltaVelocity = PP
          
          HH = 1E-5 ! based on Sulsky and Kaul (2004)
          
          TotalVelocitySoil_NPlus1_DirectionalDerivative = 0
          
          DirectionalDerivative = 0
          
          call EvaluateResidualMagnitude(TotalVelocitySoil, TotalVelocitySoil_Magnitude) 
          
          call EvaluateResidualMagnitude(DeltaVelocity, DeltaVelocity_Magnitude)
          
          
          do ii = 1, Counters%N
              
              if (DeltaVelocity(ii) .eq. 0.0) then ! equals to zero 
                  
                  
                  DirectionalDerivative(ii) = 0
                  
                  
                  
              elseif ( (DeltaVelocity(ii) .ne. 0.0) .and. (TotalVelocitySoil(ii) .ne. 0.0 ) ) then ! greater than zero 
                  
                  
                  !do jj = 1, Counters%N
                  TotalVelocitySoil_NPlus1_DirectionalDerivative(ii) = TotalVelocitySoil(ii) + (HH * TotalVelocitySoil_Magnitude * (1/DeltaVelocity_Magnitude) * DeltaVelocity(ii)  )
                  !end do 
                  
                  call EvaluateResidualValue_Scalar(TotalVelocitySoil_NPlus1_DirectionalDerivative(ii), LumpedMassDry_NPlus_1MinusAlphaM_1(ii,1), ExtLoad_NPlus_1MinusAlphaF_1(ii,1), IntLoad_NPlus_1MinusAlphaF_1(ii,1), &
                                           New_Residual(ii))
                  
                  DirectionalDerivative = DeltaVelocity_Magnitude * ( New_Residual(ii) - Old_Residual(ii) )/(HH * TotalVelocitySoil_Magnitude)
                  
          !          subroutine EvaluateResidualValue(AccelerationSoil_NPlus_1MinusAlphaM_1, LumpedMassDry_NPlus_1MinusAlphaM_1, ExtLoad_NPlus_1MinusAlphaF_1, IntLoad_NPlus_1MinusAlphaF_1, &
          !                                 Residual_GeneralizedAlpha)
          !
          !implicit none 
          !
          !
          !real(REAL_TYPE), dimension(Counters%N), intent(out) :: Residual_GeneralizedAlpha
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_NPlus_1MinusAlphaM_1
          !
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus_1MinusAlphaM_1
          !
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad_NPlus_1MinusAlphaF_1
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: ExtLoad_NPlus_1MinusAlphaF_1
          !
          !! local 
          !integer(INTEGER_TYPE) :: ii
          !!real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: LumpedMassDry
          !
          !
          !
          !do ii = 1, Counters%N
          !
          !  Residual_GeneralizedAlpha(ii) = AccelerationSoil_NPlus_1MinusAlphaM_1(ii, 1)*LumpedMassDry_NPlus_1MinusAlphaM_1(ii, 1) - ExtLoad_NPlus_1MinusAlphaF_1(ii, 1) + IntLoad_NPlus_1MinusAlphaF_1(ii, 1)
          !  ! Entity = 1 
          !  
          !
          !end do 
          !
          !end subroutine EvaluateResidualValue
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                      !call EvaluateNewResidual(AccelerationSoil, DeltaAcceleration, HH, & ! Input 
                      !                         New_Residual) !Output
                      !
                      !
                      !
                      !
                      !
                      !
                      !call EvaluateAccelerationMagnitude(AccelerationSoilMagnitude, &)
                      !
                      !call EvaluateAccelerationIncrementMagnitude(AccelerationIncrementMagnitude, &)
                      !
                      !call EvaluateDirectionalDerivative(AccelerationSoil, Old_Residual, New_Residual, HH, AccelerationSoilMagnitude, AccelerationIncrementMagnitude &
                      !    )
                          
                          
                          
              elseif ( (DeltaVelocity(ii) .ne. 0.0) .and. (TotalVelocitySoil(ii) .eq. 0.0) ) then 
                  
                    
                  
                  !do jj = 1, Counters%N
                  TotalVelocitySoil_NPlus1_DirectionalDerivative(ii) = TotalVelocitySoil(ii) + (HH * (1/DeltaVelocity_Magnitude) * DeltaVelocity(ii)  )
                  !end do 
                  
                  call EvaluateResidualValue_Scalar(TotalVelocitySoil_NPlus1_DirectionalDerivative(ii), LumpedMassDry_NPlus_1MinusAlphaM_1(ii,1), ExtLoad_NPlus_1MinusAlphaF_1(ii,1), IntLoad_NPlus_1MinusAlphaF_1(ii,1), &
                                           New_Residual(ii))
                  
                  DirectionalDerivative = DeltaVelocity_Magnitude * ( New_Residual(ii) - Old_Residual(ii) )/(HH)
                  
                    
                  
                      !call EvaluateNewResidual(AccelerationSoil, DeltaAcceleration, HH, & ! Input 
                      !                         New_Residual) !Output
                      !
                      !call EvaluateAccelerationIncrementMagnitude(AccelerationIncrementMagnitude, &)
                  
                      
                  
              
              
              
              
              
              end if 
          
          
          
          
          
              
              
          end do 
          
          
          
          
          
          
          
          
          
          end subroutine 
        
        
        
        
        
        subroutine Calculate_Beta(Residual_Zero_Magnitude, Beta)
          
          implicit none 
          
          real(REAL_TYPE), intent(in) :: Residual_Zero_Magnitude
          
          real(REAL_TYPE), dimension(Counters%N), intent(out) :: Beta 
          real(REAL_TYPE), dimension(Counters%N) :: EE_1
          
          Beta = 0
          EE_1 = 0
          
          EE_1(1) = 1
          
          Beta = Residual_Zero_Magnitude*EE_1
          
          
        end subroutine Calculate_Beta
        
        
        
        
        
        
        subroutine EvaluateErrorToleranceGMRES 
        
        implicit none 
        
        ! empty template for GMRES tolerance evaluation
        
        real(REAL_TYPE) :: TOLERANCE 
        
        
        
        
        end subroutine
        
        
        
        
        
        subroutine Calculate_PP_1(Residual_Zero, Residual_Zero_Magnitude, PP_1)
        
        implicit none 
          
          real(REAL_TYPE), intent(in) :: Residual_Zero_Magnitude
          
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: Residual_Zero
          
          real(REAL_TYPE), dimension(Counters%N), intent(out) :: PP_1
          
          
          PP_1 =  Residual_Zero*(1/Residual_Zero_Magnitude)
          
          
          
        end subroutine Calculate_PP_1
        
        
        
        
        
        
        
        
        subroutine SetResidual_Zero_To_NegativeResidualValue(Residual_GeneralizedAlpha, Residual_Zero)
        
        implicit none 
        
        real(REAL_TYPE), dimension(Counters%N) :: Residual_GeneralizedAlpha
        
        real(REAL_TYPE), dimension(Counters%N) :: Residual_Zero
        
        
        
        Residual_Zero = -Residual_GeneralizedAlpha
        
        
        
        
        
        end subroutine 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine SetAccelerationVectorValueToZero(DeltaAcceleration)
        
        implicit none 
        
        real(REAL_TYPE), dimension(Counters%N), intent(inout) :: DeltaAcceleration

        
        DeltaAcceleration = 0
        
        
        end subroutine SetAccelerationVectorValueToZero
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine EvaluateResidualValue(AccelerationSoil_NPlus_1MinusAlphaM_1, LumpedMassDry_NPlus_1MinusAlphaM_1, ExtLoad_NPlus_1MinusAlphaF_1, IntLoad_NPlus_1MinusAlphaF_1, &
                                           Residual_GeneralizedAlpha)
          
          implicit none 
          
          
          real(REAL_TYPE), dimension(Counters%N), intent(out) :: Residual_GeneralizedAlpha
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_NPlus_1MinusAlphaM_1
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus_1MinusAlphaM_1
          
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad_NPlus_1MinusAlphaF_1
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: ExtLoad_NPlus_1MinusAlphaF_1
          
          ! local 
          integer(INTEGER_TYPE) :: ii
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: LumpedMassDry

          
          
          do ii = 1, Counters%N
          
            Residual_GeneralizedAlpha(ii) = AccelerationSoil_NPlus_1MinusAlphaM_1(ii, 1)*LumpedMassDry_NPlus_1MinusAlphaM_1(ii, 1) - ExtLoad_NPlus_1MinusAlphaF_1(ii, 1) + IntLoad_NPlus_1MinusAlphaF_1(ii, 1)
            ! Entity = 1 
            
          
          end do 
          
                                           
                                           end subroutine EvaluateResidualValue
        
                                           
                                           
                                           
                                           
                                           
          !subroutine EvaluateResidualValueVelocity(AccelerationSoil_NPlus_1MinusAlphaM_1, LumpedMassDry_NPlus_1MinusAlphaM_1, ExtLoad_NPlus_1MinusAlphaF_1, IntLoad_NPlus_1MinusAlphaF_1, &
          !                                 Residual_GeneralizedAlpha)
          !
          !implicit none 
          !
          !
          !real(REAL_TYPE), dimension(Counters%N), intent(out) :: Residual_GeneralizedAlpha
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_NPlus_1MinusAlphaM_1
          !
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus_1MinusAlphaM_1
          !
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad_NPlus_1MinusAlphaF_1
          !real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: ExtLoad_NPlus_1MinusAlphaF_1
          !
          !! local 
          !integer(INTEGER_TYPE) :: ii
          !real(REAL_TYPE) :: Implicit_mu
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NN ! this is the initial velocity
          !
          !
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: InternalForceIncrement
          !
          !!real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: LumpedMassDry
          !
          !Implicit_mu = 1 ! gauge set to fully implicit
          !
          !
          !
          !
          !do ii = 1, Counters%N
          !
          !  Residual_GeneralizedAlpha(ii) = AccelerationSoil_NPlus_1MinusAlphaM_1(ii, 1)*LumpedMassDry_NPlus_1MinusAlphaM_1(ii, 1) - ExtLoad_NPlus_1MinusAlphaF_1(ii, 1) + IntLoad_NPlus_1MinusAlphaF_1(ii, 1)
          !  ! Entity = 1 
          !  Residual_GeneralizedAlpha(ii) = (TotalVelocitySoil(ii, 1)*LumpedMassDry(ii, 1)) - (Implicit_mu*CalParams%TimeIncrement*InternalForceIncrement(ii, 1)) - (TotalVelocitySoil(ii, 1)*LumpedMassDry(ii, 1))
          !                                          ! M*v^imp                                 -                 mu*DeltaT*Delta_F                                  -              Mv^exp
          !  
          !  ! I need to use the and the new internal force to be used here
          !  
          !  
          !  !- ExtLoad_NPlus_1MinusAlphaF_1(ii, 1) + IntLoad_NPlus_1MinusAlphaF_1(ii, 1)
          !
          !end do 
          !
          !                                 
          !end subroutine EvaluateResidualValue
                                           
                                           
                                           
                
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
          subroutine EvaluateResidualValue_Scalar(AccelerationSoil_NPlus_1MinusAlphaM_1, LumpedMassDry_NPlus_1MinusAlphaM_1, ExtLoad_NPlus_1MinusAlphaF_1, IntLoad_NPlus_1MinusAlphaF_1, &
                                           Residual_GeneralizedAlpha)
          
          implicit none 
          
          
          real(REAL_TYPE), intent(out) :: Residual_GeneralizedAlpha
          real(REAL_TYPE), intent(in) :: AccelerationSoil_NPlus_1MinusAlphaM_1
          
          real(REAL_TYPE), intent(in) :: LumpedMassDry_NPlus_1MinusAlphaM_1
          
          real(REAL_TYPE), intent(in) :: IntLoad_NPlus_1MinusAlphaF_1
          real(REAL_TYPE), intent(in) :: ExtLoad_NPlus_1MinusAlphaF_1
          
          ! local 
          integer(INTEGER_TYPE) :: ii
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: LumpedMassDry

          
          
          !do ii = 1, Counters%N
            Residual_GeneralizedAlpha = AccelerationSoil_NPlus_1MinusAlphaM_1*LumpedMassDry_NPlus_1MinusAlphaM_1 - ExtLoad_NPlus_1MinusAlphaF_1 + IntLoad_NPlus_1MinusAlphaF_1
            ! Entity = 1 
          !end do 
          
          end subroutine EvaluateResidualValue_Scalar
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
        
        
                                           
                                           
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine EvaluatePartialValuesForResidualEquation(AccelerationSoil_NPlus_1MinusAlphaM_1, AccelerationSoil, AccelerationSoil_NPlus1, &
                                                            LumpedMassDry_NPlus_1MinusAlphaM_1, LumpedMassDry, LumpedMassDry_NPlus1, &
                                                            ExtLoad_NPlus_1MinusAlphaF_1, ExtLoad, ExtLoad_NPlus1, &
                                                            IntLoad_NPlus_1MinusAlphaF_1, IntLoad, IntLoad_NPlus1)
        
        ! Evaluating vectors required for residual evaluation
        
        implicit none 
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_NPlus1
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus1
        
        real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: ExtLoad
        real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: ExtLoad_NPlus1
        
        real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad
        real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IntLoad_NPlus1
        
        
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: AccelerationSoil_NPlus_1MinusAlphaM_1
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(out) :: LumpedMassDry_NPlus_1MinusAlphaM_1
        
        real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(out) :: ExtLoad_NPlus_1MinusAlphaF_1
        real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(out) :: IntLoad_NPlus_1MinusAlphaF_1
        
                
        AccelerationSoil_NPlus_1MinusAlphaM_1 = (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m*AccelerationSoil) + ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m)*AccelerationSoil_NPlus1)
        
        LumpedMassDry_NPlus_1MinusAlphaM_1 = (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m*LumpedMassDry) + ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m)*LumpedMassDry_NPlus1)
        
        ExtLoad_NPlus_1MinusAlphaF_1 = (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m*ExtLoad) + ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m)*ExtLoad_NPlus1)
        
        IntLoad_NPlus_1MinusAlphaF_1 = (CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m*IntLoad) + ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m)*IntLoad_NPlus1)
        
        
        
        
        
        
        
        
        end subroutine EvaluatePartialValuesForResidualEquation
        
        ! ---------
        
        
        ! ------------ Updating the mass to N+1 value 
        !subroutine UpdateNodalMassImplicitGeneralizedAlpha
        !! output
        !! input
        !! local 
        !
        !
        !
        !
        !end subroutine UpdateNodalMassImplicitGeneralizedAlpha
        
        ! ------------
        
        ! ------------ Updating the internal force to N+1 value
        
        !subroutine UpdateNodalInternalForceImplicitGeneralizedAlpha
        !! output
        !! input
        !! local 
        !
        !
        !end subroutine UpdateNodalInternalForceImplicitGeneralizedAlpha
        
        
        
        ! ------------ residual equation 
        
        subroutine CalculateGeneralizedAlphaResidualEquation()
        !Residual_Implicit, & !output 
        !Mass_I_NPlus1, Mass_I_N, & !mass 
        !Acceleration_I_NPlus1, Acceleration_I_N, & !acceleration 
        !ForceExt_I_NPlus1, ForceExt_I_N, & !external
        !ForceInt_I_NPlus1, ForceInt_I_N, & !internal 
        !AlphaM_GAlpha, AlphaF_GAlpha ) !
        
        implicit none 
        ! output
        !Residual_Implicit
        
        ! input
        real(REAL_TYPE) :: Mass_I_NPlus1 
        real(REAL_TYPE) :: Mass_I_N
        
        
        real(REAL_TYPE) :: Acceleration_I_NPlus1
        real(REAL_TYPE) :: Acceleration_I_N
        
        
        real(REAL_TYPE) :: ForceExt_I_NPlus1 
        real(REAL_TYPE) :: ForceExt_I_N
        
        
        real(REAL_TYPE) :: ForceInt_I_NPlus1
        real(REAL_TYPE) :: ForceInt_I_N
        
        
        real(REAL_TYPE) :: AlphaM_GAlpha
        real(REAL_TYPE) :: AlphaF_GAlpha
        
        
        ! local 
        real(REAL_TYPE), dimension(Counters%N) :: Residual_LHS
        real(REAL_TYPE), dimension(Counters%N) :: Residual_RHS
        real(REAL_TYPE), dimension(Counters%N) :: Residual_Implicit!, intent(out) 
        
        real(REAL_TYPE) :: Residual_Implicit_X
        real(REAL_TYPE) :: Residual_Implicit_Y
        
        integer(INTEGER_TYPE) :: ii, jj, kk
        
        Residual_LHS = 0
        Residual_RHS = 0
        
        ! AccelerationSoil size --> Counters%N, Counters%NEntity
        ! Note that I am assuming Counters%NEntity=1
        
        !Residual_LHS = ( (1-AlphaM_GAlpha) * Mass_I_NPlus1 * Acceleration_I_NPlus1 ) + ( AlphaM_GAlpha * Mass_I_N * Acceleration_I_N )
        
        do ii = 1, Counters%N
        
            Residual_LHS(ii) = ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m) * LumpedMassDry_NPlus1(ii,1) * AccelerationSoil(ii,1) ) + ( CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m * LumpedMassDry(ii,1) * AccelerationSoilPrevious(ii,1) )
            
            Residual_RHS(ii) = ( ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f) * ExtLoad(ii,1) ) + ( CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f * ExtLoad(ii,1) ) )      -       &
                     ( ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f) * IntLoad(ii,1) ) + ( CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f * IntLoadPrevious(ii,1) ) )
        
        end do 
        
        !LumpedMassDry_NPlus1
        
        
        !Residual_RHS = ( ( (1-AlphaF_GAlpha) * ForceExt_I_NPlus1 ) + ( AlphaF_GAlpha * ForceExt_I_N ) )      -       &
        !             ( ( (1-AlphaF_GAlpha) * ForceInt_I_NPlus1 ) + ( AlphaF_GAlpha * ForceInt_I_N ) )
        
        !Residual_RHS = ( ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f) * ExtLoad ) + ( CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f * ExtLoad ) )      -       &
        !             ( ( (1-CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f) * IntLoad ) + ( CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f * IntLoadPrevious ) )
        !
        Residual_Implicit = Residual_LHS - Residual_RHS
        
        ResidualLoadGeneralizedAlpha(:,1) = Residual_Implicit
        
        Residual_Implicit_X = 0
        
        Residual_Implicit_Y = 0
        
        kk = 1
        do jj = 1, Counters%N/2
        !ForceExt
            ! looking at elements 1, 3, 5, 7
            Residual_Implicit_X = Residual_Implicit_X + ResidualLoadGeneralizedAlpha(kk,1)
            
            
            !Residual_Implicit_Y(jj)
        
            kk = kk+2
        end do 
        
        kk = 2
        do jj = 1, Counters%N/2
        !ForceExt
            ! looking at elements 1, 3, 5, 7
            Residual_Implicit_Y = Residual_Implicit_Y + ResidualLoadGeneralizedAlpha(kk,1)
            
            
            !Residual_Implicit_Y(jj)
        
            kk = kk+2
        end do 
        
        CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector(1) = Residual_Implicit_X
        
        CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector(2) = Residual_Implicit_Y
        
        end subroutine CalculateGeneralizedAlphaResidualEquation
        
        
        
        
        
        
        
        subroutine CheckResidualTolerance()!ImplicitConverged, Tolerance_ImplicitConvergence)
        
        implicit none 
        
        ! output 
        logical :: ImplicitConverged ! 1/0 variable 
        
        ! input
        real(REAL_TYPE) :: Residual_Implicit
        real(REAL_TYPE) :: Tolerance_ImplicitConvergence 
        
        real(REAL_TYPE) :: Residual_Implicit_SquareAndSumElements
        
        integer(INTEGER_TYPE) :: ii 
        ! local
        
        
        real(REAL_TYPE) :: Residual_Implicit_SquareAndSumElements_SqrRoot
        
        Residual_Implicit_SquareAndSumElements = 0
        
        Tolerance_ImplicitConvergence = 1e-2

        
        ! find the modulus of the residual implicit vector here 
        do ii = 1, NDIM !Counters%N
 
            !Residual_Implicit_SquareAndSumElements = Residual_Implicit_SquareAndSumElements + (ResidualLoadGeneralizedAlpha(ii,1))**2 ! assuming entity 1
        
            Residual_Implicit_SquareAndSumElements = Residual_Implicit_SquareAndSumElements + (CalParams%ImplicitGeneralizedAlphaTimeScheme%ResidualVector(ii)**2)
            
        end do 
        
        
        Residual_Implicit_SquareAndSumElements_SqrRoot = sqrt(Residual_Implicit_SquareAndSumElements)
        
        if (Residual_Implicit_SquareAndSumElements_SqrRoot < Tolerance_ImplicitConvergence) then
            !ImplicitConverged = .YES. 
            CalParams%ImplicitGeneralizedAlphaTimeScheme%GeneralizedAlphaImplicitTimeSchemeConverged = .TRUE. 
        else    
            !ImplicitConverged = .NO.
            CalParams%ImplicitGeneralizedAlphaTimeScheme%GeneralizedAlphaImplicitTimeSchemeConverged = .FALSE. 
        end if 
        
        
        
        
        
        
        end subroutine CheckResidualTolerance

        
        
        
        
        
        
        ! ---------------------------------GMRES CALCULATION-----------------------------------------------
        ! ---------------------------------GMRES CALCULATION-----------------------------------------------
        ! ---------------------------------GMRES CALCULATION-----------------------------------------------
        
        
        !subroutine GMRES_GeneralizedAlphaCalculation(Tolerance_ImplicitConvergence, IncrementAcceleration)
        !! This should follow Algorithm 4 in Sulsky and Kaul (2004) 
        !!
        !! output 
        !implicit none 
        !
        !real(REAL_TYPE) :: IncrementAcceleration
        !
        !
        !! input 
        !real(REAL_TYPE) :: Residual_Implicit 
        !real(REAL_TYPE) :: Tolerance_ImplicitConvergence
        !
        !! local
        !integer(INTEGER_TYPE) :: NumberOfIteration
        !
        !! inout 
        !real(REAL_TYPE), intent(inout) :: NumberOfIteration_kk
        !
        !!----------- STEP 1 -------------------------------------------------------
        !! Step 1 in Algorithm 4 in Sulsky and Kaul (2004): 
        !! "set s_0 = 0"
        !! Initialize the initial increment of acceleration equal to zero 
        !IncrementAcceleration = 0 ! DeltaAcceleration or DeltaS set to zero 
        !NumberOfIteration_kk = 0 ! kk
        !
        !! "calculate r_0 = -F(a_n) ..."
        !call CalculateGeneralizedAlphaResidualEquation(Residual_Implicit, & !output 
        !Mass_I_NPlus1, Mass_I_N, & !mass 
        !Acceleration_I_NPlus1, Acceleration_I_N, & !acceleration 
        !ForceExt_I_NPlus1, ForceExt_I_N, & !external
        !ForceInt_I_NPlus1, ForceInt_I_N, & !internal 
        !AlphaM_GAlpha, AlphaF_GAlpha ) 
        !
        !! "and rho_0 = ||r_0||^2."
        !call EvaluateResidualMagnitude(Residual_Implicit, Residual_Implicit_Magnitude)
        !rho_0 = Residual_Implicit_Magnitude**2
        !
        !! "Set k = 1" 
        !NumberOfIteration_kk = 1
        !
        !! calculate pp_1
        !pp_1 = (1/(rho**2) ) * r_0 
        !
        !! define e_1 = ( 1, 0, ..., 0 )^T --> transposed so vector form. 
        !EE_1 =  [1; 0] ! --> 
        !
        !! calculate beta
        !beta = Residual_Implicit_Magnitude * EE_1
        !
        !! error tolerance was already defined before = 10**-12 
        !
        !
        !pp_ii(1) = pp_1 ! equate pp_ii to pp_1 initially 
        !
        !
        !!----------- STEP 2 -------------------------------------------------------
        !do while (Residual_Implicit < Tolerance_ImplicitConvergence)
        !    
        !    ! deallocate pp_ii
        !    ! allocate pp_ii based on NumberOfIteration_kk
        !    ! need to store pp in each iteration so rewrite in newly allocated pp
        !    
        !    
        !    ! Step (a): "p_k+1 = D_h F(v_n, p_k)"
        !    call CalculateDirectionalDerivative_Jacob_Multiplied_DeltaAccl(pp_ii(NumberOfIteration_kk+1) )
        !    
        !    ! Step (b): " for i = 1,2,...,k set h_k,i = p_k+1^T * p_i , p_k+1 = p_k+1 - h_k,i * p_i " 
        !    ! consider putting the below in a subroutine 
        !    do ii_IterationLoop = 1, NumberOfIteration_kk
        !       
        !        ! Step (b)
        !        hh_(NumberOfIteration_kk, ii_IterationLoop) =  pp_ii(ii_IterationLoop) * pp_ii(ii_IterationLoop+1) ! I think this is a scalar
        !        pp_ii(ii_IterationLoop+1) = pp_ii(ii_IterationLoop+1) - (h_kk_ii * pp_ii(ii_IterationLoop) )
        !        
        !        ! Step (c)
        !        call EvaluateResidualMagnitude(pp_ii(ii_IterationLoop+1), pp_iiPlus1_Magnitude)    
        !        hh_(NumberOfIteration_kk+1, ii_IterationLoop) = pp_iiPlus1_Magnitude
        !        
        !        ! Step (d)
        !        pp_ii(ii_IterationLoop+1) = pp_ii(ii_IterationLoop+1) * (1/hh_(NumberOfIteration_kk+1, ii_IterationLoop))
        !       
        !        
        !        ! Step (e): Apply all previous Givens rotations to (h_1,k,..., h_k+1,k)^T
        !        if (NumberOfIteration_kk .ne. 1) then 
        !        
        !            call ApplyAllPreviousGivenRotationTohh_(NumberOfIteration_kk) ! if kk=1 then return 
        !        
        !        end if
        !        
        !        
        !        ! Step (f): Construct Givens rotations on h_k,k and h_k+1,k to eliminate h_k+1,k
        !        call ConstructNewGivensRotation() !--> only construct for the last new one
        !        
        !        
        !        ! Step (g): Apply new Given rotations to Beta; Beta_k+1 is the new norm estimate of the residual.
        !        ! NB: Beta is a rank 1 tensor (i.e., vector)
        !        call ApplyNewGivensRotationToBeta !--> only one application of Givens rotation. Before, we applied all previous. 
        !        
        !        
        !        
        !        
        !        
        !        
        !        
        !        
        !        
        !    end do 
        !    
        !    NumberOfIteration_kk = NumberOfIteration_kk + 1 
        !    
        !    
        !    
        !    !  Calculate residual of the equation 
        !    call CalculateGeneralizedAlphaResidualEquation(Residual_Implicit)
        !  
        !
        !    ! Check residual based on tolerance --> need to be smaller than 10**-12 
        !    call CheckResidualTolerance(ImplicitConverged, Tolerance_ImplicitConvergence)
        !  
        !    
        !    
        !    
        !
        !end do ! end for do while 
        !
        !
        !! Step 3: Solve the system Hy=beta where H is upper triangular matrix
        !
        !
        !
        !
        !! Step 4: Calculate the solution s = s_0 + y_1*p_1 + ... + y_k*p_k
        !
        !
        !
        !
        !
        !
        !end subroutine GMRES_GeneralizedAlphaCalculation
        
        ! ---------------------------------GMRES CALCULATION-----------------------------------------------
        ! ---------------------------------GMRES CALCULATION-----------------------------------------------
        ! ---------------------------------GMRES CALCULATION-----------------------------------------------
        
        !
        !subroutine CalculateDirectionalDerivative_Jacob_Multiplied_DeltaAccl
        !
        !! output 
        !
        !! input 
        !
        !! local
        !
        !Acceleration_N 
        !Acceleration_I_NPlus1 
        !
        ! if ( DeltaAccl = 0 ) then 
        !                   Jacob_Multiplied_DeltaAccl = 0     
        !           -->    elseif ( (accl_i_n+AlphaM .ne. 0) & (DeltaAccl .ne. 0)) then 
        !                   Jacob_Multiplied_DeltaAccl = Modulus_DeltaAccl * &
        !                   (  RESIDUAL_FUNCTION( Accl_i_n+AlphaM  +  HH_DIRECTIONAL_DERIVATIVE_COSTANT * Modulus_Accl_i_n+AlphaM * DeltaAccl/Modulus_DeltaAccl ) - RESIDUAL_FUNCTION(Accl_i_n+AlphaM) ) *
        !                   (   1/HH_DIRECTIONAL_DERIVATIVE_COSTANT    ) * (    1/Modulus_Accl_i_n+AlphaM    )
        !           -->     elseif (accl_i_n+AlphaM == 0) & (DeltaAccl .ne. 0)) then
        !                   Jacob_Multiplied_DeltaAccl = Modulus_DeltaAccl * &
        !                   (  RESIDUAL_FUNCTION( HH_DIRECTIONAL_DERIVATIVE_COSTANT * DeltaAccl/Modulus_DeltaAccl ) - RESIDUAL_FUNCTION(Acceleration = 0) ) * &
        !                   (   1/HH_DIRECTIONAL_DERIVATIVE_COSTANT    ) 
        !
        !       Jacob_Multiplied_DeltaAccl --> vector 
        !       
        !       
        !       
        !               
        !       !if ( DeltaAccl = 0 ) then 
        !       !            Jacob_Multiplied_DeltaAccl = 0     
        !       !    -->    elseif ( (accl_i_n+AlphaM .ne. 0) & (DeltaAccl .ne. 0)) then 
        !       !            Jacob_Multiplied_DeltaAccl = Modulus_DeltaAccl * &
        !       !            (  RESIDUAL_FUNCTION( Accl_i_n+AlphaM  +  HH_DIRECTIONAL_DERIVATIVE_COSTANT * Modulus_Accl_i_n+AlphaM * DeltaAccl/Modulus_DeltaAccl ) - RESIDUAL_FUNCTION(Accl_i_n+AlphaM) ) *
        !       !            (   1/HH_DIRECTIONAL_DERIVATIVE_COSTANT    ) * (    1/Modulus_Accl_i_n+AlphaM    )
        !       !    -->     elseif (accl_i_n+AlphaM == 0) & (DeltaAccl .ne. 0)) then
        !       !            Jacob_Multiplied_DeltaAccl = Modulus_DeltaAccl * &
        !       !            (  RESIDUAL_FUNCTION( HH_DIRECTIONAL_DERIVATIVE_COSTANT * DeltaAccl/Modulus_DeltaAccl ) - RESIDUAL_FUNCTION(Acceleration = 0) ) * &
        !       !            (   1/HH_DIRECTIONAL_DERIVATIVE_COSTANT    ) 
        !       !
        !       !Jacob_Multiplied_DeltaAccl --> vector 
        !
        !
        ! end subroutine CalculateDirectionalDerivative_Jacob_Multiplied_DeltaAccl
        
         
         
        
        subroutine SetRho_Zero_To_MagnitudeSquare(Residual_Implicit, Residual_Implicit_Magnitude_Squared) 
         ! Part of Step 1 in Algorithm 4
         
         implicit none 
         
         !output 
         real(REAL_TYPE), intent(out) :: Residual_Implicit_Magnitude_Squared
         !input 
         real(REAL_TYPE), intent(in), dimension(Counters%N) :: Residual_Implicit
         !local
         !real(REAL_TYPE) :: Residual_Implicit_Magnitude_Squared
         integer(INTEGER_TYPE) :: ii 
         
         Residual_Implicit_Magnitude_Squared = 0
         
         do ii = 1, Counters%N

             Residual_Implicit_Magnitude_Squared = Residual_Implicit_Magnitude_Squared + (Residual_Implicit(ii))**2
             
             !(Residual_Implicit(1)**2 + Residual_Implicit(2)**2)**0.5
         
         end do 
         
         
         !Residual_Implicit_Magnitude = Residual_Implicit_Magnitude_Squared**0.5
         
         
         end subroutine 
        
        
        
        
        
         
         
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
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         subroutine InitializeGeneralizedAlphaTimeIntegrationVariables()!Gamma_GAlpha, Beta_GAlpha, AlphaM_GAlpha, AlphaF_GAlpha)
         !
         ! Goal: calculate numerical time scheme parameters for the implicit generalized alpha time scheme 
         !
         implicit none 
         
         
         ! based on the input spectral ratio value 
         ! we can calculate all the parameters based on Kontoe et al. (2008)
         ! Alpha_m
         CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m = ( (2*CalParams%ImplicitGeneralizedAlphaTimeScheme%Rho_Infinity) -1)/(CalParams%ImplicitGeneralizedAlphaTimeScheme%Rho_Infinity +1) !Eq 3.71 in Kontoe (2006)
         
         ! Alpha_f
         CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f = CalParams%ImplicitGeneralizedAlphaTimeScheme%Rho_Infinity/(CalParams%ImplicitGeneralizedAlphaTimeScheme%Rho_Infinity+1) !Eq 3.71 in Kontoe (2006)
         
         ! Alpha --> where do we use this?
         CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha = 0.25 * (1 - CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m - CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f)**2 !Eq 3.71 in Kontoe (2006)
         
         ! Delta
         CalParams%ImplicitGeneralizedAlphaTimeScheme%Delta = 0.5 - &
             CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m + &
             CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f !Eq 3.70 in Kontoe (2006)
         
         
         
         end subroutine InitializeGeneralizedAlphaTimeIntegrationVariables
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
        ! subroutine GetNodalExtAndIntForces_NPlus1(DShapeValuesArray)
        !!**********************************************************************
        !!
        !!    Function:  To extrapolate loads from material points to nodes using 
        !!               the shape function values evaluated at the material points local position.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !
        !  integer(INTEGER_TYPE) :: IEntity, I, IDoF, J
        !  integer(INTEGER_TYPE) :: IDim, ILoadSystem
        !  
        !  logical :: DoConsiderReactionForces, IsPrescribedVelocity
        !  
        !
        !  real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
        !
        !  
        !  IDim = NVECTOR
        !  
        !   ! Internal and External forces are calculated inside "BtSig"
        !   ! Calculate internal nodal forces
        !  
        !  ! NB: I have not defined FReaction_NPlus1, FReactionWater_NPlus1, BulkViscLoad_NPlus1
        !  !     I will keep them without NPlus1
        !  
        !   call MPMDYNBTSig(ExtLoad_NPlus1, IntLoad_NPlus1, GravityLoad_NPlus1, FReaction_NPlus1, FReactionWater_NPlus1, BulkViscLoad_NPlus1, DShapeValuesArray)
        !   
        !   !if ((Counters%NLoadedElementSidesSolidNodes+Counters%NLoadedElementSidesSolidNodesB) > 0) then
        !   ! do ILoadSystem = 1, Counters%NSoilLoadSystems    
        !   !   ExtLoad_NPlus1 = ExtLoad_NPlus1 + ExtLoadTotal_NPlus1(:,:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
        !   ! end do
        !   !
        !   !end if
        !   !
        !   ! if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
        !   !ExtLoad = ExtLoad + SpaceTimeExTLoadWater
        !   !end if
        !   !
        !   !IsPrescribedVelocity = CalParams%PrescribedVelo%ApplyPrescribedVelo
        !   !
        !   !DoConsiderReactionForces =  IsPrescribedVelocity  
        !   !
        !   !if (DoConsiderReactionForces) then ! Compute in global coordinate system for output
        !   !  do I = 1, Counters%NodTot
        !   !    IDoF = ReducedDof(I)
        !   !    do J = 1, IDim
        !   !      FReaction(IDoF + J, 1:Counters%NEntity) =  &
        !   !         FReaction(IDoF + J, 1:Counters%NEntity) - GravityLoad(IDoF + J, 1:Counters%NEntity)
        !   !       if (CalParams%NumberOfPhases==2) then
        !   !      FReactionWater(IDoF + J, 1:Counters%NEntity) = &
        !   !         FReactionWater(IDoF + J, 1:Counters%NEntity) - GravityLoadWater(IDoF + J, 1:Counters%NEntity)
        !   !       else if (CalParams%NumberOfPhases==3) then
        !   !      FReactionGas(IDoF + J, 1:Counters%NEntity) = &
        !   !         FReactionGas(IDoF + J, 1:Counters%NEntity) - GravityLoadGas(IDoF + J, 1:Counters%NEntity)
        !   !       end if
        !   !    end do
        !   !  end do 
        !   !end if ! only for ConsiderReactionForeces
        !   !
        !   !if (IS3DCYLINDRIC) then ! Rotate vectors from global to local coordinate system ! 3D function
        !   !  do IEntity = 1, Counters%nEntity 
        !   !    call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoad(:, IEntity),     ExtLoad(:, IEntity))
        !   !    call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, GravityLoad(:, IEntity), GravityLoad(:, IEntity))
        !   !    call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, IntLoad(:, IEntity),     IntLoad(:, IEntity))
        !   !  end do
        !   !end if ! only for RotBoundCond
        !  
        !end subroutine GetNodalExtAndIntForces_NPlus1
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
!            SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!!*
!!*  -- LAPACK routine (version 3.1) --
!!*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!!*     November 2006
!!*
!!*     .. Scalar Arguments ..
!!      CHARACTER          DIAG, UPLO
!!      INTEGER            INFO, LDA, N
!!*     ..
!!*     .. Array Arguments ..
!!      DOUBLE PRECISION   A( LDA, * )
!!*     ..
!!*
!!*  Purpose
!!*  =======
!!*
!!*  DTRTRI computes the inverse of a real upper or lower triangular
!!*  matrix A.
!!*
!!*  This is the Level 3 BLAS version of the algorithm.
!!*
!!*  Arguments
!!*  =========
!!*
!!*  UPLO    (input) CHARACTER*1
!!*          = 'U':  A is upper triangular;
!!*          = 'L':  A is lower triangular.
!!*
!!*  DIAG    (input) CHARACTER*1
!!*          = 'N':  A is non-unit triangular;
!!*          = 'U':  A is unit triangular.
!!*
!!*  N       (input) INTEGER
!!*          The order of the matrix A.  N >= 0.
!!*
!!*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!!*          On entry, the triangular matrix A.  If UPLO = 'U', the
!!*          leading N-by-N upper triangular part of the array A contains
!!*          the upper triangular matrix, and the strictly lower
!!*          triangular part of A is not referenced.  If UPLO = 'L', the
!!*          leading N-by-N lower triangular part of the array A contains
!!*          the lower triangular matrix, and the strictly upper
!!*          triangular part of A is not referenced.  If DIAG = 'U', the
!!*          diagonal elements of A are also not referenced and are
!!*          assumed to be 1.
!!*          On exit, the (triangular) inverse of the original matrix, in
!!*          the same storage format.
!!*
!!*  LDA     (input) INTEGER
!!*          The leading dimension of the array A.  LDA >= max(1,N).
!!*
!!*  INFO    (output) INTEGER
!!*          = 0: successful exit
!!*          < 0: if INFO = -i, the i-th argument had an illegal value
!!*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!!*               matrix is singular and its inverse can not be computed.
!!*
!!*  =====================================================================
!!*
!!*     .. Parameters ..
!      DOUBLE PRECISION   ONE, ZERO
!      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!!*     ..
!!*     .. Local Scalars ..
!      LOGICAL            NOUNIT, UPPER
!      INTEGER            J, JB, NB, NN
!!*     ..
!!*     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            ILAENV
!      EXTERNAL           LSAME, ILAENV
!!*     ..
!!*     .. External Subroutines ..
!      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
!!*     ..
!!*     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!!*     ..
!!*     .. Executable Statements ..
!!*
!!*     Test the input parameters.
!!*
!      
!      CHARACTER(LEN=1) :: UPLO, DIAG
!      INTEGER :: N, LDA, INFO
!      REAL(8), DIMENSION(:,:) :: A
!      
!      INFO = 0
!      UPPER = LSAME( UPLO, 'U' )
!      NOUNIT = LSAME( DIAG, 'N' )
!      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
!         INFO = -1
!      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
!         INFO = -2
!      ELSE IF( N.LT.0 ) THEN
!         INFO = -3
!      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
!         INFO = -5
!      END IF
!      IF( INFO.NE.0 ) THEN
!         CALL XERBLA( 'DTRTRI', -INFO )
!         RETURN
!      END IF
!!*
!!*     Quick return if possible
!!*
!      IF( N.EQ.0 ) then
!        RETURN
!!*
!!*     Check for singularity if non-unit.
!!*
!      !IF( NOUNIT ) THEN
!      !   DO  INFO = 1, N
!      !      IF( A( INFO, INFO ).EQ.ZERO ) then
!      !        RETURN
!      ! CONTINUE
!      !   INFO = 0
!      !END IF
!!*
!!*     Determine the block size for this environment.
!!*
!      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!!*
!!*        Use unblocked code
!!*
!         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!      ELSE
!!*
!!*        Use blocked code
!!*
!         IF( UPPER ) THEN
!!*
!!*           Compute inverse of upper triangular matrix
!!*
!            DO 20 J = 1, N, NB
!               JB = MIN( NB, N-J+1 )
!!*
!!*              Compute rows 1:j-1 of current block column
!!*
!               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1, &
!                          JB, ONE, A, LDA, A( 1, J ), LDA )
!               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1, &
!                          JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
!!*
!!*              Compute inverse of current diagonal block
!!*
!               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
!   20       CONTINUE
!         ELSE
!!*
!!*           Compute inverse of lower triangular matrix
!!*
!            NN = ( ( N-1 ) / NB )*NB + 1
!            DO 30 J = NN, 1, -NB
!               JB = MIN( NB, N-J+1 )
!               IF( J+JB.LE.N ) THEN
!!*
!!*                 Compute rows j+jb:n of current block column
!!*
!                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG, &
!                             N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, &
!                             A( J+JB, J ), LDA )
!                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG, &
!                             N-J-JB+1, JB, -ONE, A( J, J ), LDA, &
!                             A( J+JB, J ), LDA )
!               END IF
!!*
!!*              Compute inverse of current diagonal block
!!*
!               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
!          CONTINUE
!         END IF
!      END IF
!!*
!      RETURN
!!*
!!*     End of DTRTRI
!!*
!      END

         
         
         
         
         
         
         
         
         
         
         
         
         
        
        
      end module ModDynamicImplicitGeneralizedAlphaScheme
