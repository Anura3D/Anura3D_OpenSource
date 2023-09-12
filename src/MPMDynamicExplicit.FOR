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

        subroutine RunExplicitDynamicLoadStep()
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
        integer(INTEGER_TYPE) :: mm_max_iteration
        
        ImplicitGMRESThreshold = 1e-7
        
        mm_max_iteration = MAX_ITERATIONS_IMPLICIT 
        

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

            if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
                call ComputeWaterSpecificWeight() ! compute gamma_w = rhow*g(gravity_acc*current multiplier)
            end if
        
            ! First stage effects 
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Lagrangian phase
            call LagrangianPhase(LumpedMassDry, ShapeValuesArray, ExtLoad, GravityLoad) !--> gives us the base items of the residual equation
            call ApplyMPMDynamicContact() ! only if ApplyContactAlgorithm
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ! Part 1 assuming an acceleration increment of zero 
            ! Convective phase
            call DYNConvectivePhase(LumpedMassDry, &
                                    IntLoad_NPlus1,  & ! intent out
                                    AccelerationSoil_NPlus1, &
                                    ShapeValuesArray, DShapeValuesArray, &
                                    ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1) !--> gives us the NPlus1 terms to be used to find the incremental terms

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                                    
            ! Part 2 to update mass and external forces                        
            call LagrangianPhase(LumpedMassDry_NPlus1, ShapeValuesArray_NPlus1, ExtLoad_NPlus1, GravityLoad_NPlus1) !--> gives us the base items of the residual equation
            
            
            
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            !CalculateIncrementalResidualTermsImplicitScheme
            !%step 4: calculate implicit scheme incremental terms
            call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_NPlus1, LumpedMassDry_Increment, &
                                                                    ExtLoad, ExtLoad_NPlus1, ExtLoad_Increment, &
                                                                    IntLoad, IntLoad_NPlus1, IntLoad_Increment, &
                                                                    GravityLoad, GravityLoad_NPlus1, GravityLoad_Increment, &
                                                                    AccelerationSoil, AccelerationSoil_NPlus1, AccelerationSoil_Increment)
            ! calculating incremental terms for use in the residual equation based on the assumption of AccelerationSoil_Increment
         
    
            !%step 5: at the end of the lagrangian phase, we evaluate the residual assuming that everything has initialized
            ! --> this needs to happen after updating of the new terms 
            call CalculateResidualImplicitScheme(LumpedMassDry, LumpedMassDry_Increment, & !LumpedMassDry_NPlus1, 
                                                          ExtLoad, ExtLoad_Increment, & ! ExtLoad_NPlus1,
                                                      GravityLoad, GravityLoad_Increment, & ! GravityLoad_NPlus1, 
                                                          IntLoad, IntLoad_Increment,& ! IntLoad_NPlus1, 
                                                                   ImplicitResidual, &
                                                    AccelerationSoil, AccelerationSoil_Increment)
            
        
            call FindVectorMagnitude(ImplicitResidual, Counters%N, ImplicitResidual_Magnitude) 
            
            if (ImplicitResidual_Magnitude<ImplicitGMRESThreshold) then 
                !... we reached the solution
                IsEnterGMRES = .false. 
            else 
                IsEnterGMRES = .true. 
            end if 
            
            if (IsEnterGMRES == .true.) then 
                
                !do while (ImplicitResidual_Magnitude>ImplicitGMRESThreshold)
                    
                call GMRES_DirectionalDerivative_AccelerationFormulation(& !AA, bb, xx, threshold, nn_length_AA, mm_max_iteration)
                      LumpedMassDry, & !LumpedMassDry_NPlus1, & ! in variables 
                      AccelerationSoil, & !Explicit velocity
                     !AccelerationSoil, !AccelerationSoil_NPlus1, & ! in variables 
                      IntLoad, & !IntLoad_NPlus1, & ! in variables 
                      ImplicitResidual, & ! in variables 
                      AccelerationSoil_Increment(:,1), & ! in variable
                      ImplicitGMRESThreshold, & ! in variable
                      !nn_Length_AA, & ! in variable
                      mm_max_iteration)!, & ! this is s_0 in SK04)
                    
                    
                    
                    
                !end do 
                
            
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

          end do ! ----- end TIME STEP LOOP -----

        end subroutine RunExplicitDynamicLoadStep
      
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ! implicit subroutines 
        
        !%---------------------------------------------------------------------%
        !%                                GMRES                                %
        !%---------------------------------------------------------------------%
        !---this is the mother subroutine where we use the directional derivative in the GMRES algorithm 
        subroutine GMRES_DirectionalDerivative_AccelerationFormulation( & !AA, bb, xx, threshold, nn_length_AA, mm_max_iteration)
                      LumpedMassDry, & !LumpedMassDry_NPlus1, & ! in variables 
                      AccelerationSoil, & !Explicit velocity
                      IntLoad, & !IntLoad_NPlus1, & ! in variables 
                      ImplicitResidual, & ! in variables 
                      AccelerationSoil_Incremental, & ! in variable
                      Threshold, & ! in variable
                      !nn_Length_AA, & ! in variable
                      mm_max_iteration)!, & ! this is s_0 in SK04
                      
          implicit none 
          
          !---input variables 
          !----mass 
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus1
          !----acceleration
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil ! explicit trial
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil ! this is the initial velocity          
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_NPlus1 ! this is the initial velocity
          !----internal force
          real(REAL_TYPE), dimension(Counters%N,Counters%NEntity), intent(in) :: IntLoad
          !real(REAL_TYPE), dimension(Counters%N,Counters%NEntity), intent(in) :: IntLoad_NPlus1
          !----implicit residual 
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: ImplicitResidual
          real(REAL_TYPE), intent(in) :: Threshold
          !integer(INTEGER_TYPE), intent(in) :: nn_length_AA
          !---inout variables
          real(REAL_TYPE), dimension(Counters%N), intent(inout) :: AccelerationSoil_Incremental
          integer(INTEGER_TYPE), intent(inout) :: mm_max_iteration 
          !---local variable
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: InternalForceIncrement
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
          real(REAL_TYPE), allocatable, dimension(:,:) :: QQ_Arnoldi ! vectors? --> this needs to be assembled in QQ
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
          real(REAL_TYPE) :: ImplicitResidualMagnitude
          real(REAL_TYPE) :: error 
          real(REAL_TYPE) :: ee
          integer(INTEGER_TYPE) :: INFO
        
          !---the residual is already evaluated so we do not need to evaluate this again here 
          !---using norm2 function from FOTRAN90. Inputing rank as 1, but this is not necessary (i.e., optional)
          ImplicitResidualMagnitude = norm2(ImplicitResidual, 1) !rr_norm
          
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
          QQ(:,1) = (1/ImplicitResidualMagnitude) * ImplicitResidual ! assigning this for the first iteration 
          
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
              !do while (error .ge. threshold ) 
              ! run Arnoldi 
              
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
              
              !PP_1 = QQ(:,1)
              !call Arnoldi_DirectionalDerivative(AA, QQ_Full, kk, HH_Arnoldi_New, QQ_Arnoldi, nn_length_AA)!, HH_Arnoldi_New, QQ_Arnoldi_New)
              !---the below subroutine pours directional derivative results into QQ_Full (p in SK04) and h into HH_Full 
              !---the below subroutine evaluates p in step 2(a)
              call Arnoldi_DirectionalDerivative_AccelerationFormulation(LumpedMassDry, & ! in variables
                                                 AccelerationSoil, & ! in variables 
                                                 IntLoad, & ! in variables 
                                                 kk, & ! in variable                                    
                                                 QQ_Full, & ! inout variable ! this has the p_k value for directional derivative
                                                 HH_Arnoldi_New, QQ_Arnoldi, & ! inout variable !nn_length_AA, &
                                                 ImplicitResidual) ! in variable
              
              !subroutine Arnoldi_DirectionalDerivative_VelocityFormulation(LumpedMassDry, & ! in variables 
              !                                   TotalVelocitySoil_NPlus1, & ! in variables 
              !                                   IntLoad_NPlus1, & ! in variables
              !                                   kk, &! in variable                                          
              !                                   QQ_Full, &! inout variable 
              !                                   HH_Arnoldi_New, QQ_Arnoldi, & ! inout variable                 !nn_length_AA, &
              !                                   Old_ImplicitEquationResidual) ! in variable 
                  
              
              !---not sure what this nn paramter is doing. Consider removing. 
              nn = nn + 1
              

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
              
              allocate(QQ_Stored(Counters%N,kk+1))  !kk+1
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
             
              if (error .le. threshold ) then 
                  exit
              end if 
              
          end do 
          
          
          allocate(HH_inverted(kk,kk))
          HH_inverted = 0.0
          
          ! find the inverse of upper triangular H matrix
          call L_inv(HH_Stored(1:kk, 1:kk), & ! in variable, 
                     HH_inverted, & ! out variable 
                     kk)  ! in variable  
          
          
          !call dtrtri( 'U', 'N', Counters%N, HH_inverted, Counters%N, INFO  )        
          !call dtrtri( 'U', 'N', Counters%N, HH_Stored(1:kk, 1:kk), Counters%N, INFO  )

          
          !!!!!!!!!!!!!!!!!!!!!!!!! GMRES Step 3: solve the system Hy=beta for y
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
          
          
          
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!! GMRES Step 4: Calculate the solution s = s_0 + y1*p_1 + ... + y_k*p_k
          !xx = xx + (QQ_Full(:,1:kk) * YY)
          ! loop over rows 
          do ii = 1, Counters%N ! rank of matrix 
              do jj = 1, kk
                  AccelerationSoil_Incremental(ii) = AccelerationSoil_Incremental(ii) + QQ_Full(ii,jj)*YY(jj)
              end do 
          end do 
          
          
          
        end subroutine GMRES_DirectionalDerivative_AccelerationFormulation
        
                      
                      
        
                      
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        subroutine Arnoldi_DirectionalDerivative_AccelerationFormulation(LumpedMassDry, & ! in variables 
                                                 AccelerationSoil, &!TotalVelocitySoil_NPlus1, & ! in variables 
                                                 IntLoad, & ! in variables
                                                 kk, &! in variable                                          
                                                 QQ_Full, &! inout variable 
                                                 HH_Arnoldi_New, QQ_Arnoldi, & ! inout variable                 !nn_length_AA, &
                                                 Old_ImplicitEquationResidual) ! in variable 
        !---description: this subroutine covers steps 2 a,b,c,d in Algorithm 4 in SK04
        
        !% inputs 
        !% - AA is the matrix on the left hand side to be inverted 
        !% - QQ is initially the normalized residual vector 
        !% - kk is the number of iterations 
        !% outputs 
        !% - HH_Arnoldi is the step vector to be placed in H
        !% - QQ_Arnoldi is the new Q step vector used to expand Q 
        
        implicit none 
                  
        !---input variables 
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry
        !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NPlus1
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil
        !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: LumpedMassDry_NPlus1
        !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil ! this is the initial velocity
        !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_NPlus1
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: IntLoad
        !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: IntLoad_NPlus1
        integer(INTEGER_TYPE), intent(in) :: kk
        real(REAL_TYPE), dimension(Counters%N), intent(in) :: Old_ImplicitEquationResidual
        !---local variables 
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: DeltaAcceleration
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: DirectionalDerivative
        real(REAL_TYPE), allocatable, dimension(:) :: HH_Arnoldi ! vectors?
        integer(INTEGER_TYPE) :: ii, jj
        real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: InternalForceIncrement
        !---inout variables 
        real(REAL_TYPE), dimension(:,:), intent(inout) :: QQ_Full ! Rank 2 matrix? !allocatable, 
        real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: QQ_Arnoldi ! vectors? --> this needs to be assembled in QQ
        real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HH_Arnoldi_New ! vectors?
        
        
        
        ! consider moving this outside
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !allocate(HH_Arnoldi(kk))
        allocate(HH_Arnoldi(kk)) ! 1 
        !allocate(QQ_Arnoldi(kk+1)) ! 2
        if (.not.allocated(QQ_Arnoldi)) then
            allocate(QQ_Arnoldi(Counters%N,1)) ! 2
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        QQ_Arnoldi = 0 
        HH_Arnoldi = 0 

        
        ! need to call directional derivative here 
        !call EvaluateDirectionalDerivativeUsingVelocityAsAnInput(& !Old_ImplicitEquationResidual, & ! F(v) 
        !                                                             !AccelerationSoil, &
        !                                                             QQ_Full, &
        !                                                             kk, &
        !                                                             TotalVelocitySoil_NPlus1, &
        !                                                             DirectionalDerivative)
        
        
        call EvaluateDirectionalDerivativeUsingAccelerationAsAnInput(Old_ImplicitEquationResidual, & ! F(v) 
                                                                     !AccelerationSoil, &
                                                                     QQ_Full, &
                                                                     kk, &
                                                                     AccelerationSoil, &
                                                                     DirectionalDerivative)
        
        
        
        ! p_k+1 = DirectionalDerivative(AccelerationSoil_NPlus1, OldImplicitResidual))
        QQ_Arnoldi = DirectionalDerivative
        
        !---EDIT: THIS SHOULD BE QQ(:,kk) --> It should depend on the timestep 
        !---step 2b h_k,i = p_k+1 * p_i 
        do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix            
            do jj = 1,Counters%N !2 
                HH_Arnoldi(ii) = HH_Arnoldi(ii) + ( QQ_Arnoldi(jj,1) * QQ_Full(jj,ii) ) !QQ_Arnoldi is the p_k+1 term. QQ is the p_i term 
            end do 
        end do 
        
        !---step 2b in algorithm 4 in SK04
        !----Equation p_k+1 = p_k+1 - (h_k,i * p_i)
        do ii = 1, kk !% Modified Gram-Schmidt, keeping the Hessenberg matrix
            QQ_Arnoldi(:,1) = QQ_Arnoldi(:,1) - ( HH_Arnoldi(ii) * QQ_Full(:,ii)) !QQ_Arnoldi is the p_k+1 term. QQ is the p_i term 
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
        
        !---step 1c h_k+1,k = |p_k+1|
        HH_Arnoldi_New(kk+1, 1) =  norm2(QQ_Arnoldi(:,1),1)
        
        !---step 1d p_k+1 = p_k+1/|p_k+1| 
        QQ_Arnoldi = QQ_Arnoldi * (1/HH_Arnoldi_New(kk+1, 1) ) ! divide by a scalar
        
        !H size is (number of iterations +1, number of iterations)
        !Q size is (2 i.e., rank of matrix or how many rows in matrix, number of iterations +1)
        
        
        !QQ(:,kk+1) = QQ_Arnoldi(:,1)
        
        !---pouring q (or called as p_k+1 in SK04) back to QQ_Full to store    
        QQ_Full(:,kk+1) = QQ_Arnoldi(:,1)
        
                                       

        end subroutine Arnoldi_DirectionalDerivative_AccelerationFormulation 
        
        
        
        
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
                                                 
                                               
         
         !subroutine EvaluateDirectionalDerivativeUsingVelocityAsAnInput &
         !    !(Residual_N, &
         !    !(AccelerationSoil, 
         !   (QQ_Full, &
         !   kk, &
         !   TotalVelocitySoil_NPlus1, & !in variables 
         !   DirectionalDerivative) !out variable
         ! !---step 1 a in Algorithm 4 in SK04
         ! !---This function evaluates the directional derivative in order to be used with 
         ! !---the GMRES algorithm 
         !   
         ! 
         ! implicit none 
         ! 
         ! !---input
         ! !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NPlus1
         ! real(REAL_TYPE), dimension(:,:), intent(in) :: QQ_Full ! p_k in step 2a 
         ! integer(INTEGER_TYPE), intent(in) :: kk
         ! !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Residual_N
         ! !---output 
         ! real(REAL_TYPE), dimension(Counters%N), intent(out) :: DirectionalDerivative 
         ! !--local 
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Residual_NPlus1
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: DeltaTotalVelocitySoil
         ! real(REAL_TYPE), dimension(Counters%N) :: Zeros_Vector 
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_DirectionalDerivative
         ! real(REAL_TYPE) :: HH
         ! real(REAL_TYPE) :: VerySmallNumber
         ! real(REAL_TYPE) :: TotalVelocitySoil_NPlus1_Magnitude, DeltaTotalVelocitySoil_Magnitude
         ! real(REAL_TYPE) :: DeltaVelocitySoil_Magnitude
         ! integer(INTEGER_TYPE) :: ii
         ! 
         ! 
         ! ! left-hand side 
         ! !---directional derivative variables initialization here
         ! !----particle
         ! real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: AccelerationArray_DirectionalDerivative
         ! real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: LumpedMassDry_DirectionalDerivative
         ! real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: VelocityArray_DirectionalDerivative
         ! real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: UArray_DirectionalDerivative
         ! !----nodes
         ! !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccelerationSoil_DirectionalDerivative
         ! !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_DirectionalDerivative
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IncrementalDisplacementSoil_DirectionalDerivative
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_NPlus1_DirectionalDerivative
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: ExtLoad_DirectionalDerivative
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: GravityLoad_DirectionalDerivative
         ! !----shape functions
         ! real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES) :: ShapeValuesArray_NPlus1_DirectionalDerivative
         ! real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES,NVECTOR) :: DShapeValuesArray_NPlus1_DirectionalDerivative
         ! !----residual of new acceleration directional derivative 
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_DirectionalDerivative
         !
         ! ! right-hand side
         ! !---directional derivative variables initialization here
         ! !----particle
         ! real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: AccelerationArray_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: LumpedMassDry_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: VelocityArray_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: UArray_DirectionalDerivative_Zero
         ! !----nodes
         ! !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccelerationSoil_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IncrementalDisplacementSoil_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: ExtLoad_NPlus1_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_NPlus1_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_DirectionalDerivative_Zero
         ! !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: GravityLoad_NPlus1_DirectionalDerivative_Zero
         ! !----shape functions
         ! real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES) :: ShapeValuesArray_NPlus1_DirectionalDerivative_Zero
         ! real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES,NVECTOR) :: DShapeValuesArray_NPlus1_DirectionalDerivative_Zero
         ! !----residual of new acceleration directional derivative 
         ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside
         ! 
         ! 
         ! 
         ! !---define parameters from SK04
         ! !---the value of h does not matter too much because it gets normalized when considering p (see step 1d in Algorithm 4 in SK04)
         ! HH = 1E-5!1E-5 !this should be a small number
         ! VerySmallNumber=1E-10
         ! !---initialization of variables
         ! Residual_NPlus1 = 0.0 !vector
         ! DirectionalDerivative = 0.0 !vector
         ! DeltaTotalVelocitySoil = 0.0 !vector
         ! TotalVelocitySoil_NPlus1_Magnitude = 0.0 !scalar
         ! DeltaTotalVelocitySoil_Magnitude = 0.0 !scalar
         ! !---initialization of directional derivative vectors --> left-hand side parameter. In case 2, that is all you need.  
         ! AccelerationArray_DirectionalDerivative = 0.0
         ! VelocityArray_DirectionalDerivative = 0.0
         ! UArray_DirectionalDerivative = 0.0
         ! TotalVelocitySoil_DirectionalDerivative = 0.0
         ! IncrementalDisplacementSoil_DirectionalDerivative = 0.0
         ! IntLoad_NPlus1_DirectionalDerivative = 0.0
         ! !IntLoad_DirectionalDerivative = 0.0
         ! ShapeValuesArray_NPlus1_DirectionalDerivative = 0.0
         ! DShapeValuesArray_NPlus1_DirectionalDerivative = 0.0
         ! LumpedMassDry_DirectionalDerivative = 0.0
         ! Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside = 0.0
         ! !AccelerationSoil_DirectionalDerivative = 0.0
         ! !---initialization of directional derivative vectors --> right-hand side parameter. In case 3, you need this. 
         ! AccelerationArray_DirectionalDerivative_Zero = 0.0
         ! VelocityArray_DirectionalDerivative_Zero = 0.0
         ! UArray_DirectionalDerivative_Zero = 0.0
         ! TotalVelocitySoil_DirectionalDerivative_Zero = 0.0
         ! IncrementalDisplacementSoil_DirectionalDerivative_Zero = 0.0
         ! IntLoad_NPlus1_DirectionalDerivative_Zero = 0.0
         ! ExtLoad_NPlus1_DirectionalDerivative_Zero = 0.0
         ! ShapeValuesArray_NPlus1_DirectionalDerivative_Zero = 0.0
         ! DShapeValuesArray_NPlus1_DirectionalDerivative_Zero = 0.0
         ! LumpedMassDry_DirectionalDerivative_Zero = 0.0
         ! Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside = 0.0
         ! !AccelerationSoil_DirectionalDerivative_Zero = 0.0
         ! 
         ! 
         ! 
         ! 
         ! !---comparison vectors 
         ! !DeltaAcceleration = AccelerationSoil_NPlus1 - AccelerationSoil
         ! DeltaTotalVelocitySoil(:,1) = QQ_Full(:,kk) ! p_k
         ! Zeros_Vector = 0.0
         ! 
         ! !---evaluate the magnitude 
         ! call EvaluateResidualMagnitude(TotalVelocitySoil_NPlus1, TotalVelocitySoil_NPlus1_Magnitude) 
         ! call EvaluateResidualMagnitude(DeltaTotalVelocitySoil, DeltaTotalVelocitySoil_Magnitude)
         !
         !
         ! if (DeltaTotalVelocitySoil_Magnitude .lt. VerySmallNumber) then ! equals to zero 
         !     DirectionalDerivative = Zeros_Vector !0.0
         !             
         ! elseif ( (DeltaTotalVelocitySoil_Magnitude .gt. VerySmallNumber) .and. (TotalVelocitySoil_NPlus1_Magnitude .gt. VerySmallNumber) ) then ! greater than zero 
         !     
         !     !----------------- Acceleration left term SK04 F(h * s/|s|)
         !     !AccelerationSoil_DirectionalDerivative = AccelerationSoil_NPlus1 + &
         !     !    (HH * AccelerationSoil_NPlus1_Magnitude * DeltaAcceleration * (1/DeltaAcceleration_Magnitude) )
         !     
         !     
         !     TotalVelocitySoil_DirectionalDerivative = TotalVelocitySoil_NPlus1 + & ! this would start as the explicit velocity 
         !         (HH * TotalVelocitySoil_NPlus1_Magnitude * DeltaTotalVelocitySoil * (1/DeltaTotalVelocitySoil_Magnitude) )
         !     
         !     ! find its corresponding internal force               
         !     call EvaluateForceFromVelocity(TotalVelocitySoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ShapeValuesArray, DShapeValuesArray)                       
         !         
         !     
         !     !---update LumpedMassDry based on the initial acceleration 
         !     !call UpdateMassAfterConvectivePhase(LumpedMassDry_DirectionalDerivative, & ! out variable
         !     !        ShapeValuesArray_NPlus1_DirectionalDerivative) ! in variable
         !
         !     !---note below that I have put ExtLoad for both ExtLoad and ExtLoad_NPlus1. This would be different if we were to change the ExtLoad_NPlus1 in a time dependent manner. 
         !     !call EvaluateImplicitResidualFunction_SK04_VelocityFormulation & ! Directional derivative
         !     !        (TotalVelocitySoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative, &   ! NPlus1 terms (inputs)
         !     !    IntLoad, ExtLoad, GravityLoad, &                                               ! N terms      (inputs) 
         !     !    Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside, &                                                   ! Residual of implicit equation (outputs)
         !     !    TotalVelocitySoil, LumpedMassDry)                                                ! N Acceleration terms (inputs)
         ! 
         ! 
         !     Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside(:,1) = Residual_GeneralizedAlpha
         !     !---Directional derivative calculation
         !     !----the reason there is a 1 is because the directional derivative is a vector (rank 1). But Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside is rank 2. 
         !     DirectionalDerivative = DeltaTotalVelocitySoil_Magnitude * &
         !         (Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside(:,1) - Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside(:,1)) * &
         !         (1/HH) * (1/TotalVelocitySoil_NPlus1_Magnitude)
         !     
         !     
         ! elseif ( (DeltaTotalVelocitySoil_Magnitude .gt. VerySmallNumber) .and. (TotalVelocitySoil_NPlus1_Magnitude .lt. VerySmallNumber) ) then 
         !     
         !     
         !     !----------------- Acceleration left term SK04 F(h * s/|s|)
         !     TotalVelocitySoil_DirectionalDerivative = HH * DeltaTotalVelocitySoil * (1/DeltaTotalVelocitySoil_Magnitude)
         !     
         !     call EvaluateForceFromVelocity(TotalVelocitySoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ShapeValuesArray, DShapeValuesArray)                       
         !
         !     !---update LumpedMassDry based on the initial acceleration 
         !     !call UpdateMassAfterConvectivePhase(LumpedMassDry_DirectionalDerivative, & ! out variable
         !     !        ShapeValuesArray_NPlus1_DirectionalDerivative) ! in variable
         !
         !     !---note below that I have put ExtLoad for both ExtLoad and ExtLoad_NPlus1. This would be different if we were to change the ExtLoad_NPlus1 in a time dependent manner.                                                                                                         
         !     call EvaluateImplicitResidualFunction_SK04_VelocityFormulation & ! Directional derivative
         !             (TotalVelocitySoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative, &   ! NPlus1 terms (inputs)
         !               IntLoad, ExtLoad, GravityLoad, &                                               ! N terms      (inputs) 
         !               Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside, &                                                   ! Residual of implicit equation (outputs)
         !               TotalVelocitySoil, LumpedMassDry)                                                ! N Acceleration terms (inputs)
         ! 
         !     
         !     !----------------- Acceleration right term SK04 F(0) (i.e., assuming acceleration is equal to zero).
         !     !---it might not be worth saving all these variables
         !     TotalVelocitySoil_DirectionalDerivative_Zero = 0.0
         !     
         !     !!---update LumpedMassDry based on the initial acceleration 
         !     !call UpdateMassAfterConvectivePhase(LumpedMassDry_DirectionalDerivative_Zero, & ! out variable
         !     !        ShapeValuesArray_NPlus1_DirectionalDerivative_Zero) ! in variable
         !     !---note below that I have put ExtLoad for both ExtLoad and ExtLoad_NPlus1. This would be different if we were to change the ExtLoad_NPlus1 in a time dependent manner. 
         !     call EvaluateImplicitResidualFunction_SK04_VelocityFormulation & ! Directional derivative
         !             (TotalVelocitySoil_DirectionalDerivative_Zero, IntLoad_NPlus1_DirectionalDerivative_Zero, ExtLoad_NPlus1_DirectionalDerivative_Zero, GravityLoad_NPlus1_DirectionalDerivative_Zero, &   ! NPlus1 terms (inputs) --> ExtLoad_NPlus1
         !                 IntLoad, ExtLoad, GravityLoad, &                                                                                                                                                    ! N terms      (inputs) 
         !                 Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside, &                                                                                                                    ! Residual of implicit equation (outputs)
         !                 TotalVelocitySoil, LumpedMassDry)                                                                                                                                                   ! N Acceleration terms (inputs)
         !         
         !     !---Directional derivative calculation
         !     !----the reason there is a 1 is because the directional derivative is a vector (rank 1). But Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside is rank 2. 
         !     DirectionalDerivative = DeltaTotalVelocitySoil_Magnitude * &
         !         (Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside(:,1) - Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside(:,1)) * &
         !         (1/HH)
         !         
         !      
         ! end if 
         ! 
         ! 
         ! 
         !   end subroutine EvaluateDirectionalDerivativeUsingVelocityAsAnInput
        
        
            
            
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            
            
            subroutine EvaluateDirectionalDerivativeUsingAccelerationAsAnInput &
             !(Residual_N, &
             !(AccelerationSoil, 
            (Old_ImplicitEquationResidual, &
            QQ_Full, &
            kk, &
            AccelerationSoil, &!TotalVelocitySoil_NPlus1, & !in variables 
            DirectionalDerivative) !out variable
          !---step 1 a in Algorithm 4 in SK04
          !---This function evaluates the directional derivative in order to be used with 
          !---the GMRES algorithm 
            
          
          implicit none 
          
          !---input
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: Old_ImplicitEquationResidual !,Counters%nEntity
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: TotalVelocitySoil_NPlus1
          real(REAL_TYPE), dimension(:,:), intent(in) :: QQ_Full ! p_k in step 2a 
          integer(INTEGER_TYPE), intent(in) :: kk
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Residual_N
          !---output 
          real(REAL_TYPE), dimension(Counters%N), intent(out) :: DirectionalDerivative 
          !--local 
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Residual_NPlus1
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: DeltaTotalVelocitySoil
          real(REAL_TYPE), dimension(Counters%N) :: Zeros_Vector 
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_DirectionalDerivative
          real(REAL_TYPE) :: HH
          real(REAL_TYPE) :: VerySmallNumber
          real(REAL_TYPE) :: TotalVelocitySoil_NPlus1_Magnitude, DeltaTotalVelocitySoil_Magnitude
          real(REAL_TYPE) :: DeltaVelocitySoil_Magnitude
          integer(INTEGER_TYPE) :: ii
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: DeltaAccelerationSoil
          real(REAL_TYPE) :: DeltaAccelerationSoil_Magnitude, AccelerationSoil_Magnitude
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccelerationSoil_DirectionalDerivative
          ! left-hand side 
          !---directional derivative variables initialization here
          !----particle
          real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: AccelerationArray_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: VelocityArray_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: UArray_DirectionalDerivative
          !----nodes
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccelerationSoil_DirectionalDerivative
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IncrementalDisplacementSoil_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_NPlus1_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: ExtLoad_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: GravityLoad_DirectionalDerivative

          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: LumpedMassDry_DirectionalDerivative

          !----shape functions
          real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES) :: ShapeValuesArray_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES,NVECTOR) :: DShapeValuesArray_DirectionalDerivative
          !----residual of new acceleration directional derivative 
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_DirectionalDerivative

          ! right-hand side
          !---directional derivative variables initialization here
          !----particle
          real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: AccelerationArray_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: LumpedMassDry_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: VelocityArray_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%NParticles,NVECTOR) :: UArray_DirectionalDerivative_Zero
          !----nodes
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccelerationSoil_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: TotalVelocitySoil_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IncrementalDisplacementSoil_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: ExtLoad_NPlus1_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_NPlus1_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_DirectionalDerivative_Zero
          !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: GravityLoad_NPlus1_DirectionalDerivative_Zero
          !----shape functions
          real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES) :: ShapeValuesArray_NPlus1_DirectionalDerivative_Zero
          real(REAL_TYPE), dimension(Counters%NParticles,ELEMENTNODES,NVECTOR) :: DShapeValuesArray_NPlus1_DirectionalDerivative_Zero
          !----residual of new acceleration directional derivative 
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Residual_DirectionalDerivative_Lefthandside
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Residual_DirectionalDerivative_Righthandside
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: LumpedMassDry_Increment
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: ExtLoad_Increment
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: GravityLoad_Increment
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_Increment
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccelerationSoil_DirectionalDerivative_Zero
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          ! implicit scheme - shape values array inout 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES) :: ShapeValuesArray_DirectionalDerivative 
          !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR) :: DShapeValuesArray_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: LumpedMassDry_Increment_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: ExtLoad_Increment_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: IntLoad_Increment_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: GravityLoad_Increment_DirectionalDerivative
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: AccelerationSoil_Increment_DirectionalDerivative
          
          !---define parameters from SK04
          !---the value of h does not matter too much because it gets normalized when considering p (see step 1d in Algorithm 4 in SK04)
          HH = 1E-5!1E-5 !this should be a small number
          VerySmallNumber=1E-10
          !---initialization of variables
          Residual_NPlus1 = 0.0 !vector
          DirectionalDerivative = 0.0 !vector
          DeltaTotalVelocitySoil = 0.0 !vector
          TotalVelocitySoil_NPlus1_Magnitude = 0.0 !scalar
          DeltaTotalVelocitySoil_Magnitude = 0.0 !scalar
          !---initialization of directional derivative vectors --> left-hand side parameter. In case 2, that is all you need.  
          AccelerationArray_DirectionalDerivative = 0.0
          VelocityArray_DirectionalDerivative = 0.0
          UArray_DirectionalDerivative = 0.0
          TotalVelocitySoil_DirectionalDerivative = 0.0
          IncrementalDisplacementSoil_DirectionalDerivative = 0.0
          IntLoad_NPlus1_DirectionalDerivative = 0.0
          !IntLoad_DirectionalDerivative = 0.0
          !ShapeValuesArray_NPlus1_DirectionalDerivative = 0.0
          !DShapeValuesArray_NPlus1_DirectionalDerivative = 0.0
          LumpedMassDry_DirectionalDerivative = 0.0
          Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside = 0.0
          !AccelerationSoil_DirectionalDerivative = 0.0
          !---initialization of directional derivative vectors --> right-hand side parameter. In case 3, you need this. 
          AccelerationArray_DirectionalDerivative_Zero = 0.0
          VelocityArray_DirectionalDerivative_Zero = 0.0
          UArray_DirectionalDerivative_Zero = 0.0
          TotalVelocitySoil_DirectionalDerivative_Zero = 0.0
          IncrementalDisplacementSoil_DirectionalDerivative_Zero = 0.0
          IntLoad_NPlus1_DirectionalDerivative_Zero = 0.0
          ExtLoad_NPlus1_DirectionalDerivative_Zero = 0.0
          !ShapeValuesArray_NPlus1_DirectionalDerivative_Zero = 0.0
          !DShapeValuesArray_NPlus1_DirectionalDerivative_Zero = 0.0
          LumpedMassDry_DirectionalDerivative_Zero = 0.0
          Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside = 0.0
          !AccelerationSoil_DirectionalDerivative_Zero = 0.0
          
          ! internal force initialization
          IntLoad_NPlus1 = 0 
          
          
          !---comparison vectors 
          !DeltaAcceleration = AccelerationSoil_NPlus1 - AccelerationSoil
          ! we obtain DeltaAcceleration from the p_k evaluation poured in QQ_Full
          DeltaAccelerationSoil(:,1) = QQ_Full(:,kk) ! p_k
          Zeros_Vector = 0.0
          
          !---evaluate the magnitude 
          !call EvaluateResidualMagnitude(AccelerationSoil, AccelerationSoil_Magnitude) 
          call EvaluateResidualMagnitude(AccelerationSoil, AccelerationSoil_Magnitude)
          call EvaluateResidualMagnitude(DeltaAccelerationSoil, DeltaAccelerationSoil_Magnitude)

  
          if (DeltaAccelerationSoil_Magnitude .lt. VerySmallNumber) then ! equals to zero 
              DirectionalDerivative = Zeros_Vector !0.0
                      
          elseif ( (DeltaAccelerationSoil_Magnitude .gt. VerySmallNumber) .and. (AccelerationSoil_Magnitude .gt. VerySmallNumber) ) then ! greater than zero 
              
              !----------------- Acceleration left term SK04 F(h * s/|s|)
              !AccelerationSoil_DirectionalDerivative = AccelerationSoil_NPlus1 + &
              !    (HH * AccelerationSoil_NPlus1_Magnitude * DeltaAcceleration * (1/DeltaAcceleration_Magnitude) )
              
              !call LagrangianPhase(LumpedMassDry_DirectionalDerivative) !--> gives us the base items of the residual equation
              
              AccelerationSoil_DirectionalDerivative = AccelerationSoil + & ! this would start as the explicit velocity 
                  (HH * AccelerationSoil_Magnitude * DeltaAccelerationSoil * (1/DeltaAccelerationSoil_Magnitude) )
              
              !-->we need a to do a convective phase to find DeltaIntForce, new shape functions, 
              !--> we need to do a lagrangian phase to find DeltaMass
              
              !call DYNConvectivePhase(LumpedMassDry, IntLoad_NPlus1, ExtLoad_NPlus1, GravityLoad_NPlus1, & ! intent out
              !                      AccelerationSoil_DirectionalDerivative) !--> gives us the NPlus1 terms to be used to find the incremental terms
             call DYNConvectivePhase(LumpedMassDry, &!_DirectionalDerivative, &
                                     IntLoad_DirectionalDerivative, &
                                     !ExtLoad_DirectionalDerivative, &
                                     !GravityLoad_DirectionalDerivative, & ! intent out    
                                     AccelerationSoil_DirectionalDerivative, &
                                     ShapeValuesArray, &
                                     DShapeValuesArray, &
                                     ShapeValuesArray_DirectionalDerivative, & 
                                     DShapeValuesArray_DirectionalDerivative) !--> gives us the NPlus1 terms to be used to find the incremental terms
             
              
              
              ! Lagrangian phase
              ! calculating the shape functions based on the new shape values
              call LagrangianPhase(LumpedMassDry_DirectionalDerivative, ShapeValuesArray_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative) !--> gives us the base items of the residual equation
            !call ApplyMPMDynamicContact() ! only if ApplyContactAlgorithm
            !
            !! Convective phase
            !call DYNConvectivePhase(LumpedMassDry_NPlus1, IntLoad_NPlus1, ExtLoad_NPlus1, GravityLoad_NPlus1, & ! intent out
            !                        AccelerationSoil_NPlus1) !--> gives us the NPlus1 terms to be used to find the incremental terms
            !
            !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
            !call LagrangianPhase(LumpedMassDry_NPlus1) !--> gives us the base items of the residual equation
            !
            !!CalculateIncrementalResidualTermsImplicitScheme
            !!%step 4: calculate implicit scheme incremental terms
            !call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_NPlus1, LumpedMassDry_Increment, &
            !                                                        ExtLoad, ExtLoad_NPlus1, ExtLoad_Increment, &
            !                                                        IntLoad, IntLoad_NPlus1, IntLoad_Increment, &
            !                                                        GravityLoad, GravityLoad_NPlus1, GravityLoad_Increment)
              
              
              ! find its corresponding internal force               
              ! I commented this on 9/11/2023
              !call EvaluateForceFromAcceleration(AccelerationSoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ShapeValuesArray, DShapeValuesArray)                       
                  
              
              !---update LumpedMassDry based on the initial acceleration 
              !call UpdateMassAfterConvectivePhase(LumpedMassDry_DirectionalDerivative, & ! out variable
              !        ShapeValuesArray_NPlus1_DirectionalDerivative) ! in variable

              !---note below that I have put ExtLoad for both ExtLoad and ExtLoad_NPlus1. This would be different if we were to change the ExtLoad_NPlus1 in a time dependent manner. 
              ! THIS NEEDS TO BE CHANGED TO ACCELERATION FORMULATION
              !call EvaluateImplicitResidualFunction_SK04_VelocityFormulation & ! Directional derivative
              !        (TotalVelocitySoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative, &   ! NPlus1 terms (inputs)
              !    IntLoad, ExtLoad, GravityLoad, &                                               ! N terms      (inputs) 
              !    Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside, &                                                   ! Residual of implicit equation (outputs)
              !    TotalVelocitySoil, LumpedMassDry)                                                ! N Acceleration terms (inputs)
          
              
              
              !call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_NPlus1, LumpedMassDry_Increment, &
              !                                                      ExtLoad, ExtLoad_NPlus1, ExtLoad_Increment, &
              !                                                      IntLoad, IntLoad_NPlus1, IntLoad_Increment, &
              !                                                      GravityLoad, GravityLoad_NPlus1, GravityLoad_Increment)
              call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_DirectionalDerivative, LumpedMassDry_Increment_DirectionalDerivative, &
                                                                    ExtLoad, ExtLoad_DirectionalDerivative, ExtLoad_Increment_DirectionalDerivative, &
                                                                    IntLoad, IntLoad_DirectionalDerivative, IntLoad_Increment_DirectionalDerivative, &
                                                                    GravityLoad, GravityLoad_DirectionalDerivative, GravityLoad_Increment_DirectionalDerivative, &
                                                                    AccelerationSoil, AccelerationSoil_DirectionalDerivative, AccelerationSoil_Increment_DirectionalDerivative)
            ! calculating incremental terms for use in the residual equation based on the assumption of AccelerationSoil_Increment
         
    
            !%step 5: at the end of the lagrangian phase, we evaluate the residual assuming that everything has initialized
            ! --> this needs to happen after updating of the new terms 
            call CalculateResidualImplicitScheme(LumpedMassDry, LumpedMassDry_Increment_DirectionalDerivative, & !LumpedMassDry_NPlus1, 
                                                          ExtLoad, ExtLoad_Increment_DirectionalDerivative, & ! ExtLoad_NPlus1,
                                                      GravityLoad, GravityLoad_Increment_DirectionalDerivative, & ! GravityLoad_NPlus1, 
                                                          IntLoad, IntLoad_Increment_DirectionalDerivative,& ! IntLoad_NPlus1, 
                                                                   Residual_DirectionalDerivative_Lefthandside(:,1), &
                                                    AccelerationSoil, AccelerationSoil_Increment_DirectionalDerivative)
              
              
              
          
              Residual_DirectionalDerivative_Righthandside(:,1) = Old_ImplicitEquationResidual
              !---Directional derivative calculation
              !----the reason there is a 1 is because the directional derivative is a vector (rank 1). But Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside is rank 2. 
              DirectionalDerivative = DeltaAccelerationSoil_Magnitude * & !<-- this is p_k
                  (Residual_DirectionalDerivative_Lefthandside(:,1) - Residual_DirectionalDerivative_Righthandside(:,1)) * &
                  (1/HH) * (1/AccelerationSoil_Magnitude)
              
              
          elseif ( (DeltaAccelerationSoil_Magnitude .gt. VerySmallNumber) .and. (AccelerationSoil_Magnitude .lt. VerySmallNumber) ) then 
              
              
              !----------------- Acceleration left term SK04 F(h * s/|s|)
              AccelerationSoil_DirectionalDerivative = HH * DeltaAccelerationSoil * (1/DeltaAccelerationSoil_Magnitude)
              
              
              
              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
              
              !-->we need a to do a convective phase to find DeltaIntForce, new shape functions, 
              !--> we need to do a lagrangian phase to find DeltaMass
              
              call DYNConvectivePhase(LumpedMassDry, &
                                      IntLoad_NPlus1, & ! intent out                      
                                      AccelerationSoil_DirectionalDerivative, &
                                      ShapeValuesArray, DShapeValuesArray, &
                                      ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1) !--> gives us the NPlus1 terms to be used to find the incremental terms
             
              
              
              ! Lagrangian phase
              ! calculating the shape functions based on the new shape values
              call LagrangianPhase(LumpedMassDry_DirectionalDerivative, ShapeValuesArray_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative) !--> gives us the base items of the residual equation
            !call ApplyMPMDynamicContact() ! only if ApplyContactAlgorithm
            !
            !! Convective phase
            !call DYNConvectivePhase(LumpedMassDry_NPlus1, IntLoad_NPlus1, ExtLoad_NPlus1, GravityLoad_NPlus1, & ! intent out
            !                        AccelerationSoil_NPlus1) !--> gives us the NPlus1 terms to be used to find the incremental terms
            !
            !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
            !call LagrangianPhase(LumpedMassDry_NPlus1) !--> gives us the base items of the residual equation
            !
            !!CalculateIncrementalResidualTermsImplicitScheme
            !!%step 4: calculate implicit scheme incremental terms
            !call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_NPlus1, LumpedMassDry_Increment, &
            !                                                        ExtLoad, ExtLoad_NPlus1, ExtLoad_Increment, &
            !                                                        IntLoad, IntLoad_NPlus1, IntLoad_Increment, &
            !                                                        GravityLoad, GravityLoad_NPlus1, GravityLoad_Increment)
              
              
              ! find its corresponding internal force               
              ! I commented this on 9/11/2023
              !call EvaluateForceFromAcceleration(AccelerationSoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ShapeValuesArray, DShapeValuesArray)                       
                  
              
              !---update LumpedMassDry based on the initial acceleration 
              !call UpdateMassAfterConvectivePhase(LumpedMassDry_DirectionalDerivative, & ! out variable
              !        ShapeValuesArray_NPlus1_DirectionalDerivative) ! in variable

              !---note below that I have put ExtLoad for both ExtLoad and ExtLoad_NPlus1. This would be different if we were to change the ExtLoad_NPlus1 in a time dependent manner. 
              ! THIS NEEDS TO BE CHANGED TO ACCELERATION FORMULATION
              !call EvaluateImplicitResidualFunction_SK04_VelocityFormulation & ! Directional derivative
              !        (TotalVelocitySoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative, &   ! NPlus1 terms (inputs)
              !    IntLoad, ExtLoad, GravityLoad, &                                               ! N terms      (inputs) 
              !    Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside, &                                                   ! Residual of implicit equation (outputs)
              !    TotalVelocitySoil, LumpedMassDry)                                                ! N Acceleration terms (inputs)
          
              
              
              call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_NPlus1, LumpedMassDry_Increment, &
                                                                    ExtLoad, ExtLoad_NPlus1, ExtLoad_Increment, &
                                                                    IntLoad, IntLoad_NPlus1, IntLoad_Increment, &
                                                                    GravityLoad, GravityLoad_NPlus1, GravityLoad_Increment, &
                                                                    AccelerationSoil, AccelerationSoil_NPlus1, AccelerationSoil_Increment)
            ! calculating incremental terms for use in the residual equation based on the assumption of AccelerationSoil_Increment
         
              
              call CalculateResidualImplicitScheme(LumpedMassDry, LumpedMassDry_Increment, & !LumpedMassDry_NPlus1, 
                                                         ExtLoad, ExtLoad_Increment, & ! ExtLoad_NPlus1,
                                                     GravityLoad, GravityLoad_Increment, & ! GravityLoad_NPlus1, 
                                                         IntLoad, IntLoad_Increment,& ! IntLoad_NPlus1, 
                                                    Residual_DirectionalDerivative_Lefthandside(:,1), &
                                                AccelerationSoil, AccelerationSoil_DirectionalDerivative)
              
              
              
              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
              
              
              
              !call EvaluateForceFromAcceleration(AccelerationSoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ShapeValuesArray, DShapeValuesArray)                       

              !---update LumpedMassDry based on the initial acceleration 
              !call UpdateMassAfterConvectivePhase(LumpedMassDry_DirectionalDerivative, & ! out variable
              !        ShapeValuesArray_NPlus1_DirectionalDerivative) ! in variable

              !---note below that I have put ExtLoad for both ExtLoad and ExtLoad_NPlus1. This would be different if we were to change the ExtLoad_NPlus1 in a time dependent manner.                                                                                                         
              ! THIS NEEDS TO BE CHANGED TO ACCELERATION FORMULATION
              !call EvaluateImplicitResidualFunction_SK04_VelocityFormulation & ! Directional derivative
              !        (AccelerationSoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative, &   ! NPlus1 terms (inputs)
              !          IntLoad, ExtLoad, GravityLoad, &                                               ! N terms      (inputs) 
              !          Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside, &                                                   ! Residual of implicit equation (outputs)
              !          TotalVelocitySoil, LumpedMassDry)                                                ! N Acceleration terms (inputs)
              !%step 5: at the end of the lagrangian phase, we evaluate the residual assuming that everything has initialized
              ! --> this needs to happen after updating of the new terms 
              
              
              
              
              
              !----------------- Acceleration right term SK04 F(0) (i.e., assuming acceleration is equal to zero).
              !---it might not be worth saving all these variables
              !TotalVelocitySoil_DirectionalDerivative_Zero = 0.0
              AccelerationSoil_DirectionalDerivative_Zero = 0.0

              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
              
              !-->we need a to do a convective phase to find DeltaIntForce, new shape functions, 
              !--> we need to do a lagrangian phase to find DeltaMass
              
              call DYNConvectivePhase(LumpedMassDry, &
                                    IntLoad_NPlus1, &
                                    ShapeValuesArray, DShapeValuesArray, &
                                    ShapeValuesArray_DirectionalDerivative, DShapeValuesArray_DirectionalDerivative, &
                                    AccelerationSoil_DirectionalDerivative_Zero) !--> gives us the NPlus1 terms to be used to find the incremental terms
             !ExtLoad_NPlus1, GravityLoad_NPlus1, & ! intent out
              
              
              ! Lagrangian phase
              ! calculating the shape functions based on the new shape values
              call LagrangianPhase(LumpedMassDry_DirectionalDerivative, ShapeValuesArray_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative) !--> gives us the base items of the residual equation
            !call ApplyMPMDynamicContact() ! only if ApplyContactAlgorithm
            !
            !! Convective phase
            !call DYNConvectivePhase(LumpedMassDry_NPlus1, IntLoad_NPlus1, ExtLoad_NPlus1, GravityLoad_NPlus1, & ! intent out
            !                        AccelerationSoil_NPlus1) !--> gives us the NPlus1 terms to be used to find the incremental terms
            !
            !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !
            !call LagrangianPhase(LumpedMassDry_NPlus1) !--> gives us the base items of the residual equation
            !
            !!CalculateIncrementalResidualTermsImplicitScheme
            !!%step 4: calculate implicit scheme incremental terms
            !call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_NPlus1, LumpedMassDry_Increment, &
            !                                                        ExtLoad, ExtLoad_NPlus1, ExtLoad_Increment, &
            !                                                        IntLoad, IntLoad_NPlus1, IntLoad_Increment, &
            !                                                        GravityLoad, GravityLoad_NPlus1, GravityLoad_Increment)
              
              
              ! find its corresponding internal force               
              ! I commented this on 9/11/2023
              !call EvaluateForceFromAcceleration(AccelerationSoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ShapeValuesArray, DShapeValuesArray)                       
                  
              
              !---update LumpedMassDry based on the initial acceleration 
              !call UpdateMassAfterConvectivePhase(LumpedMassDry_DirectionalDerivative, & ! out variable
              !        ShapeValuesArray_NPlus1_DirectionalDerivative) ! in variable

              !---note below that I have put ExtLoad for both ExtLoad and ExtLoad_NPlus1. This would be different if we were to change the ExtLoad_NPlus1 in a time dependent manner. 
              ! THIS NEEDS TO BE CHANGED TO ACCELERATION FORMULATION
              !call EvaluateImplicitResidualFunction_SK04_VelocityFormulation & ! Directional derivative
              !        (TotalVelocitySoil_DirectionalDerivative, IntLoad_DirectionalDerivative, ExtLoad_DirectionalDerivative, GravityLoad_DirectionalDerivative, &   ! NPlus1 terms (inputs)
              !    IntLoad, ExtLoad, GravityLoad, &                                               ! N terms      (inputs) 
              !    Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside, &                                                   ! Residual of implicit equation (outputs)
              !    TotalVelocitySoil, LumpedMassDry)                                                ! N Acceleration terms (inputs)
          
              
              
              call CalculateIncrementalResidualTermsImplicitScheme(LumpedMassDry, LumpedMassDry_NPlus1, LumpedMassDry_Increment, &
                                                                    ExtLoad, ExtLoad_NPlus1, ExtLoad_Increment, &
                                                                    IntLoad, IntLoad_NPlus1, IntLoad_Increment, &
                                                                    GravityLoad, GravityLoad_NPlus1, GravityLoad_Increment, &
                                                                    AccelerationSoil, AccelerationSoil_NPlus1, AccelerationSoil_Increment)
            
              ! calculating incremental terms for use in the residual equation based on the assumption of AccelerationSoil_Increment
              call CalculateResidualImplicitScheme(LumpedMassDry, LumpedMassDry_Increment, & !LumpedMassDry_NPlus1, 
                                                          ExtLoad, ExtLoad_Increment, & ! ExtLoad_NPlus1,
                                                      GravityLoad, GravityLoad_Increment, & ! GravityLoad_NPlus1, 
                                                          IntLoad, IntLoad_Increment,& ! IntLoad_NPlus1, 
                                                                   Residual_DirectionalDerivative_Righthandside(:,1), &
                                                    AccelerationSoil, AccelerationSoil_Increment)
              
              
              
              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              !!---update LumpedMassDry based on the initial acceleration 
              !call UpdateMassAfterConvectivePhase(LumpedMassDry_DirectionalDerivative_Zero, & ! out variable
              !        ShapeValuesArray_NPlus1_DirectionalDerivative_Zero) ! in variable
              !---note below that I have put ExtLoad for both ExtLoad and ExtLoad_NPlus1. This would be different if we were to change the ExtLoad_NPlus1 in a time dependent manner. 
              ! THIS NEEDS TO BE ACCELERATION FORMULATION
              !call EvaluateImplicitResidualFunction_SK04_VelocityFormulation & ! Directional derivative
              !        (TotalVelocitySoil_DirectionalDerivative_Zero, IntLoad_NPlus1_DirectionalDerivative_Zero, ExtLoad_NPlus1_DirectionalDerivative_Zero, GravityLoad_NPlus1_DirectionalDerivative_Zero, &   ! NPlus1 terms (inputs) --> ExtLoad_NPlus1
              !            IntLoad, ExtLoad, GravityLoad, &                                                                                                                                                    ! N terms      (inputs) 
              !            Residual_GeneralizedAlpha_DirectionalDerivative_Righthandside, &                                                                                                                    ! Residual of implicit equation (outputs)
              !            TotalVelocitySoil, LumpedMassDry)                                                                                                                                                   ! N Acceleration terms (inputs)
                  
              !---Directional derivative calculation
              !----the reason there is a 1 is because the directional derivative is a vector (rank 1). But Residual_GeneralizedAlpha_DirectionalDerivative_Lefthandside is rank 2. 
              DirectionalDerivative = DeltaAccelerationSoil_Magnitude * &
                  (Residual_DirectionalDerivative_Lefthandside(:,1) - Residual_DirectionalDerivative_Righthandside(:,1)) * & ! I removed generalizedalpha naming from here...
                  (1/HH)
                  
               
          end if 
          
          
          
            end subroutine EvaluateDirectionalDerivativeUsingAccelerationAsAnInput
            
            
            
        
            
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
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
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
                         
        
                         
                         
                         
                         
                         
                         
                         
                         
        
        
      end module ModDynamicExplicit
