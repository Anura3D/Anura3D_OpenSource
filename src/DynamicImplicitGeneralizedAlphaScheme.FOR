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
      implicit none

    contains ! Routines of this module

    
        ! ------------ main implicit kernel    
    
        subroutine DynamicImplicitGeneralizedAlphaScheme()
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
        
        
        
        ! output
        ! input
        
        
        
        Logical :: ImplicitConverged ! 1/0 switch 
        ! local 
        Tolerance_ImplicitConvergence = 10**-12
        
        
          !! Local variables
          !integer(INTEGER_TYPE) :: NUpdateParticles, NRemovedParticles, &
          !           NAddedParticles, I
          !real(REAL_TYPE) :: TotalVirtualParticleMass, &
          !                    TotalParticleMass
        
            
          !  Initialize generalized alpha variables 
          call InitializeGeneralizedAlphaTimeIntegrationVariables(Gamma_GAlpha, Beta_GAlpha, AlphaM_GAlpha, AlphaF_GAlpha)
        
          
          !!!!!!!!!!!!!!!! STEP 1: MAP FROM MP TO NODES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  All the variables are initially at the material points after a time step. So first this here to map acceleration 
          AccelerationSoil ! size = number of CPs x NDIM
          TotalVelocitySoil ! size = number of CPs x NDIM
          TotalDisplacementSoil ! size = number of CPs x NDIM
        
          !  Get mass soil at the nodes 
          Mass_I_N = LumpedMassDry ! size = number of CPs x NDIM
          
          !  Get kinematic variables at the nodes 
          !  Here I relabel them so that it would make more sense
          Acceleration_I_N = AccelerationSoil ! size = number of CPs x NDIM
          Velocity_I_N = TotalVelocitySoil ! size = number of CPs x NDIM
          Displacement_I_N = TotalDisplacementSoil ! size = number of CPs x NDIM
        
          
          
          
          

          !!!!!!!!!!!!!!!! STEP 2: PREDICTOR STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  Predictor stage evaluation of the kinematic vector variables 
          ImplicitIterations = 0 
          Velocity_I_NPlus1 = Velocity_I_N ! size = number of CPs x NDIM
          Acceleration_I_NPlus1 = ( (Gamma_GAlpha-1)/Gamma_GAlpha ) * Acceleration_I_N ! size = number of CPs x NDIM
          Displacement_I_NPlus1 = TotalDisplacementSoil + &
              ( CalParams%TimeStep * Velocity_I_N) + &
              ( 0.5 * (CalParams%TimeStep**2) * ( ((1 - 2*Beta_GAlpha)*Acceleration_I_N) + (2*Beta_GAlpha*Acceleration_I_NPlus1 ) ) )
          ! size = number of CPs x NDIM
        
          !  Update positions of the material points based on the new values
        
          call UpdateParticleDisplacementImplicitGeneralizedAlpha(GlobPos_I_NPlus1)
          
          call UpdateParticleStrains
          
          call UpdateParticleStress
          
          call UpdateParticlePositionsImplicitGeneralizedAlpha(GlobPos_I_NPlus1)
        
          call UpdateNodalMassImplicitGeneralizedAlpha(Mass_I_NPlus1)  ! --> Based on the new positions we can calculate Mass_I_NPlus1
          
          call UpdateNodalInternalForceImplicitGeneralizedAlpha(ForceInt_I_NPlus1) ! -> using new particle stress
          
          call UpdateNodalExternalForceImplicitGeneralizedAlpha(ForceExt_I_NPlus1)
          
          
          
          
          
          
          
          
          
          
          !!!!!!!!!!!!!!!! STEP 3: CHECK CONVERGENCE BASED ON RESIDUAL EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          !  Calculate residual of the equation 
          
          ! plug in new mass values, new acceleration values and evaluate residual Residual_Implicit
          
          call CalculateGeneralizedAlphaResidualEquation(Residual_Implicit)
          
          ! Check residual based on tolerance --> need to be smaller than 10**-12 
          call CheckResidualTolerance(ImplicitConverged, Tolerance_ImplicitConvergence)
          
          
          if (ImplicitConverged == .YES.) then 
            
              ! if converged then exit here 
              RETURN 
            
          end if 
          
          
          
          
          
          
          
          
          
          
          
          
          !!!!!!!!!!!!!!!! STEP 4: FINDING INCREMENT STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! In this stage we solve for acceleration increment that would achieve the tolerance level 
          

          if (ImplicitConverged == .NO.) then 
          
              call GMRES_GeneralizedAlphaCalculation(Tolerance_ImplicitConvergence, IncrementAcceleration)
          
              ! --> output here is IncrementAcceleration
          
          
          
          
          
          
          
          
          
          
          
          
          
          !!!!!!!!!!!!!!!! STEP 5: CORRECTOR STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! In this stage we use IncrementAcceleration to find N+1 values 
          
          
          Acceleration_I_NPlus1 = AccelerationSoil + IncrementAcceleration
          Velocity_I_NPlus1 = TotalVelocitySoil + ( IncrementAcceleration*Gamma_GAlpha*CalParams%TimeIncrement )
          Displacement_I_NPlus1 = TotalDisplacementSoil + ( IncrementAcceleration*Beta_GAlpha* (CalParams%TimeIncrement)**2 )
          
          
          
          
          
          
          
          
          
          
          
          
          end if 
          
          
          
                    
          !!!!!!!!!!!!!!!! STEP 6: ERROR TESTING STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          Factor = Residual_PreviousIteration/Residual_CurrentIteration
          
          if (Factor < Tolerance_ImplicitConvergence) then 
              
              ImplicitConverged = .YES. 
              
          end if 
          
          
          
          
          
          
          
          
          
          
          
    
        end subroutine DynamicImplicitGeneralizedAlphaScheme

        ! ------------
        
        ! ------------ Updating the mass to N+1 value 
        subroutine UpdateNodalMassImplicitGeneralizedAlpha
        ! output
        ! input
        ! local 
        
        
        
        
        end subroutine UpdateNodalMassImplicitGeneralizedAlpha
        
        ! ------------
        
        ! ------------ Updating the internal force to N+1 value
        
        subroutine UpdateNodalInternalForceImplicitGeneralizedAlpha
        ! output
        ! input
        ! local 
        
        
        end subroutine UpdateNodalInternalForceImplicitGeneralizedAlpha
        
        
        
        ! ------------ residual equation 
        
        subroutine CalculateGeneralizedAlphaResidualEquation(Residual_Implicit, & !output 
        Mass_I_NPlus1, Mass_I_N, & !mass 
        Acceleration_I_NPlus1, Acceleration_I_N, & !acceleration 
        ForceExt_I_NPlus1, ForceExt_I_N, & !external
        ForceInt_I_NPlus1, ForceInt_I_N, & !internal 
        AlphaM_GAlpha, AlphaF_GAlpha ) !
        
        implicit none 
        ! output
        Residual_Implicit
        
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
        real(REAL_TYPE) :: Residual_LHS
        real(REAL_TYPE) :: Residual_RHS 
        
        
        Residual_LHS = ( (1-AlphaM_GAlpha) * Mass_I_NPlus1 * Acceleration_I_NPlus1 ) + ( AlphaM_GAlpha * Mass_I_N * Acceleration_I_N )
        
        Residual_RHS = ( ( (1-AlphaF_GAlpha) * ForceExt_I_NPlus1 ) + ( AlphaF_GAlpha * ForceExt_I_N ) )      -       &
                     ( ( (1-AlphaF_GAlpha) * ForceInt_I_NPlus1 ) + ( AlphaF_GAlpha * ForceInt_I_N ) )
        
        Residual_Implicit = Residual_LHS - Residual_RHS
        
        
        end subroutine CalculateGeneralizedAlphaResidualEquation
        
        
        
        
        
        
        
        subroutine CheckResidualTolerance(ImplicitConverged, Tolerance_ImplicitConvergence)
        
        ! output 
        logical :: ImplicitConverged ! 1/0 variable 
        
        ! input
        real(REAL_TYPE) :: Residual_Implicit
        real(REAL_TYPE) :: Tolerance_ImplicitConvergence 
        
        ! local
        
        
        
        
        
        
        
        if (Residual_Implicit < Tolerance_ImplicitConvergence) then
            ImplicitConverged = .YES. 
            
        else    
            ImplicitConverged = .NO.
            
        end if 
        
        
        
        
        
        
        end subroutine CheckResidualTolerance

        
        
        
        
        
        
        
        
        subroutine GMRES_GeneralizedAlphaCalculation(Tolerance_ImplicitConvergence, IncrementAcceleration)
        ! This should follow Algorithm 4 in Sulsky and Kaul (2004) 
        !
        ! output 
        real(REAL_TYPE) :: IncrementAcceleration
        
        
        ! input 
        real(REAL_TYPE) :: Residual_Implicit 
        real(REAL_TYPE) :: Tolerance_ImplicitConvergence
        
        ! local
        integer(INTEGER_TYPE) :: NumberOfIteration
        
        ! inout 
        real(REAL_TYPE), intent(inout) :: NumberOfIteration_kk
        
        !----------- STEP 1 -------------------------------------------------------
        ! Step 1 in Algorithm 4 in Sulsky and Kaul (2004): 
        ! "set s_0 = 0"
        ! Initialize the initial increment of acceleration equal to zero 
        IncrementAcceleration = 0 ! DeltaAcceleration or DeltaS set to zero 
        NumberOfIteration_kk = 0 ! kk
        
        ! "calculate r_0 = -F(a_n) ..."
        call CalculateGeneralizedAlphaResidualEquation(Residual_Implicit, & !output 
        Mass_I_NPlus1, Mass_I_N, & !mass 
        Acceleration_I_NPlus1, Acceleration_I_N, & !acceleration 
        ForceExt_I_NPlus1, ForceExt_I_N, & !external
        ForceInt_I_NPlus1, ForceInt_I_N, & !internal 
        AlphaM_GAlpha, AlphaF_GAlpha ) 
        
        ! "and rho_0 = ||r_0||^2."
        call EvaluateResidualMagnitude(Residual_Implicit, Residual_Implicit_Magnitude)
        rho_0 = Residual_Implicit_Magnitude**2
        
        ! "Set k = 1" 
        NumberOfIteration_kk = 1
        
        ! calculate pp_1
        pp_1 = (1/(rho**2) ) * r_0 
        
        ! define e_1 = ( 1, 0, ..., 0 )^T --> transposed so vector form. 
        EE_1 =  [1; 0] ! --> 
        
        ! calculate beta
        beta = Residual_Implicit_Magnitude * EE_1
        
        ! error tolerance was already defined before = 10**-12 
        
        
        pp_ii(1) = pp_1 ! equate pp_ii to pp_1 initially 
        
        
        !----------- STEP 2 -------------------------------------------------------
        do while (Residual_Implicit < Tolerance_ImplicitConvergence)
            
            ! deallocate pp_ii
            ! allocate pp_ii based on NumberOfIteration_kk
            ! need to store pp in each iteration so rewrite in newly allocated pp
            
            
            ! Step (a): "p_k+1 = D_h F(v_n, p_k)"
            call CalculateDirectionalDerivative_Jacob_Multiplied_DeltaAccl(pp_ii(NumberOfIteration_kk+1) )
            
            ! Step (b): " for i = 1,2,...,k set h_k,i = p_k+1^T * p_i , p_k+1 = p_k+1 - h_k,i * p_i " 
            ! consider putting the below in a subroutine 
            do ii_IterationLoop = 1, NumberOfIteration_kk
               
                ! Step (b)
                hh_(NumberOfIteration_kk, ii_IterationLoop) =  pp_ii(ii_IterationLoop) * pp_ii(ii_IterationLoop+1) ! I think this is a scalar
                pp_ii(ii_IterationLoop+1) = pp_ii(ii_IterationLoop+1) - (h_kk_ii * pp_ii(ii_IterationLoop) )
                
                ! Step (c)
                call EvaluateResidualMagnitude(pp_ii(ii_IterationLoop+1), pp_iiPlus1_Magnitude)    
                hh_(NumberOfIteration_kk+1, ii_IterationLoop) = pp_iiPlus1_Magnitude
                
                ! Step (d)
                pp_ii(ii_IterationLoop+1) = pp_ii(ii_IterationLoop+1) * (1/hh_(NumberOfIteration_kk+1, ii_IterationLoop))
               
                
                ! Step (e): Apply all previous Givens rotations to (h_1,k,..., h_k+1,k)^T
                if (NumberOfIteration_kk .ne. 1) then 
                
                    call ApplyAllPreviousGivenRotationTohh_(NumberOfIteration_kk) ! if kk=1 then return 
                
                end if
                
                
                ! Step (f): Construct Givens rotations on h_k,k and h_k+1,k to eliminate h_k+1,k
                call ConstructNewGivensRotation() !--> only construct for the last new one
                
                
                ! Step (g): Apply new Given rotations to Beta; Beta_k+1 is the new norm estimate of the residual.
                ! NB: Beta is a rank 1 tensor (i.e., vector)
                call ApplyNewGivensRotationToBeta !--> only one application of Givens rotation. Before, we applied all previous. 
                
                
                
                
                
                
                
                
                
            end do 
            
            NumberOfIteration_kk = NumberOfIteration_kk + 1 
            
            
            
            !  Calculate residual of the equation 
            call CalculateGeneralizedAlphaResidualEquation(Residual_Implicit)
          
 
            ! Check residual based on tolerance --> need to be smaller than 10**-12 
            call CheckResidualTolerance(ImplicitConverged, Tolerance_ImplicitConvergence)
          
            
            
            
        
        end do ! end for do while 
        
        
        ! Step 3: Solve the system Hy=beta where H is upper triangular matrix
        
        
        
        
        ! Step 4: Calculate the solution s = s_0 + y_1*p_1 + ... + y_k*p_k
        
        
        
        
        
        
        end subroutine GMRES_GeneralizedAlphaCalculation
        
        
        
        subroutine CalculateDirectionalDerivative_Jacob_Multiplied_DeltaAccl
        
        ! output 
        
        ! input 
        
        ! local
        
        Acceleration_N 
        Acceleration_I_NPlus1 
        
         if ( DeltaAccl = 0 ) then 
                           Jacob_Multiplied_DeltaAccl = 0     
                   -->    elseif ( (accl_i_n+AlphaM .ne. 0) & (DeltaAccl .ne. 0)) then 
                           Jacob_Multiplied_DeltaAccl = Modulus_DeltaAccl * &
                           (  RESIDUAL_FUNCTION( Accl_i_n+AlphaM  +  HH_DIRECTIONAL_DERIVATIVE_COSTANT * Modulus_Accl_i_n+AlphaM * DeltaAccl/Modulus_DeltaAccl ) - RESIDUAL_FUNCTION(Accl_i_n+AlphaM) ) *
                           (   1/HH_DIRECTIONAL_DERIVATIVE_COSTANT    ) * (    1/Modulus_Accl_i_n+AlphaM    )
                   -->     elseif (accl_i_n+AlphaM == 0) & (DeltaAccl .ne. 0)) then
                           Jacob_Multiplied_DeltaAccl = Modulus_DeltaAccl * &
                           (  RESIDUAL_FUNCTION( HH_DIRECTIONAL_DERIVATIVE_COSTANT * DeltaAccl/Modulus_DeltaAccl ) - RESIDUAL_FUNCTION(Acceleration = 0) ) * &
                           (   1/HH_DIRECTIONAL_DERIVATIVE_COSTANT    ) 
    
               Jacob_Multiplied_DeltaAccl --> vector 
               
               
               
                       
               !if ( DeltaAccl = 0 ) then 
               !            Jacob_Multiplied_DeltaAccl = 0     
               !    -->    elseif ( (accl_i_n+AlphaM .ne. 0) & (DeltaAccl .ne. 0)) then 
               !            Jacob_Multiplied_DeltaAccl = Modulus_DeltaAccl * &
               !            (  RESIDUAL_FUNCTION( Accl_i_n+AlphaM  +  HH_DIRECTIONAL_DERIVATIVE_COSTANT * Modulus_Accl_i_n+AlphaM * DeltaAccl/Modulus_DeltaAccl ) - RESIDUAL_FUNCTION(Accl_i_n+AlphaM) ) *
               !            (   1/HH_DIRECTIONAL_DERIVATIVE_COSTANT    ) * (    1/Modulus_Accl_i_n+AlphaM    )
               !    -->     elseif (accl_i_n+AlphaM == 0) & (DeltaAccl .ne. 0)) then
               !            Jacob_Multiplied_DeltaAccl = Modulus_DeltaAccl * &
               !            (  RESIDUAL_FUNCTION( HH_DIRECTIONAL_DERIVATIVE_COSTANT * DeltaAccl/Modulus_DeltaAccl ) - RESIDUAL_FUNCTION(Acceleration = 0) ) * &
               !            (   1/HH_DIRECTIONAL_DERIVATIVE_COSTANT    ) 
               !
               !Jacob_Multiplied_DeltaAccl --> vector 
        
        
         end subroutine CalculateDirectionalDerivative_Jacob_Multiplied_DeltaAccl
        
         
         
         
         
         
         subroutine EvaluateResidualMagnitude(Residual_Implicit, Residual_Implicit_Magnitude) 
         ! Part of Step 1 in Algorithm 4
         
         implicit none 
         
         !output 
         real(REAL_TYPE), intent(out) :: Residual_Implicit_Magnitude
         !input 
         real(REAL_TYPE), intent(in) :: Residual_Implicit
         !local
         
         
         Residual_Implicit_Magnitude = (Residual_Implicit(1)**2 + Residual_Implicit(2)**2)**0.5
         
         
         end subroutine 
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         subroutine InitializeGeneralizedAlphaTimeIntegrationVariables()!Gamma_GAlpha, Beta_GAlpha, AlphaM_GAlpha, AlphaF_GAlpha)
         
         ! Goal: calculate numerical time scheme parameters for the implicit generalized alpha time scheme 
         
         implicit none 
         
         ! local
         
         
         ! input 
         !real(REAL_TYPE), intent(in) :: Rho_Infinity 
         
         ! output 
         
         
         ! These variables have been initialized in ReadCalculcationData.FOR 
         ! under type(ImplicitImplicitGeneralizedAlphaTimeSchemeType)
         
         
         ! based on the input spectral ratio value 
         ! we can calculate all the parameters based on Kontoe et al. (2008)
         ! Alpha_m
         CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_m = ( (2*CalParams%GeneralizedAlphaTimeScheme%Rho_Infinity) -1)/(CalParams%GeneralizedAlphaTimeScheme%Rho_Infinity +1) !Eq 3.71 in Kontoe (2006)
         
         ! Alpha_f
         CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha_f = CalParams%GeneralizedAlphaTimeScheme%Rho_Infinity/(CalParams%GeneralizedAlphaTimeScheme%Rho_Infinity+1) !Eq 3.71 in Kontoe (2006)
         
         ! Alpha --> where do we use this?
         CalParams%ImplicitGeneralizedAlphaTimeScheme%GAlpha = 0.25 * (1 - CalParams%GeneralizedAlphaTimeScheme%GAlpha_m - CalParams%GeneralizedAlphaTimeScheme%GAlpha_f)**2 !Eq 3.71 in Kontoe (2006)
         
         ! Delta
         CalParams%ImplicitGeneralizedAlphaTimeScheme%Delta = 0.5 - &
             CalParams%GeneralizedAlphaTimeScheme%GAlpha_m + &
             CalParams%GeneralizedAlphaTimeScheme%GAlpha_f !Eq 3.70 in Kontoe (2006)
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         end subroutine InitializeGeneralizedAlphaTimeIntegrationVariables
         
         
         
         
         
         
         
        
        
      end module ModDynamicImplicitGeneralizedAlphaScheme
