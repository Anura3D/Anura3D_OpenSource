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

        end subroutine RunExplicitDynamicLoadStep
      
      end module ModDynamicExplicit
