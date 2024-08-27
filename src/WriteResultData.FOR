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
	  
	  
	  module ModWriteResultData
      !**********************************************************************
      !
      !  Function : Contains routines related to writing results data to files
      !
      !             This module should only contain routines that are directly related to
      !             the output of result data of load steps and time steps.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModMPMData
      use ModCounters
      use ModReadCalculationData
      use ModMeshInfo
      use ModMPMDYN2PhaseSP 
      use ModWriteVTKTwoLayer
      use ModWriteNodalData
      use ModWriteMPMData 
      use ModFeedback
      use WriteOutput_GiD
      use ModString
      use kernel32
      use ModFileIO
      
#ifdef __INTEL_COMPILER
      use IFPORT
#endif 

      implicit none

      contains
      
        subroutine WriteTimeStepResults(AnErrorHappend)
        !**********************************************************************
        !
        !  Function : Write result data for intermediate time steps or quick check output
        !     
        !**********************************************************************
            
        implicit none
          
          logical, intent(in) :: AnErrorHappend

          ! Local variables
          character(len = 1023) :: Filename
          logical :: WriteQuickOutput, WriteTimeStep
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: OutputTime(MAXOUTPUTSUBSTEPS) = -1.0
          
          WriteQuickOutput = (CalParams%ApplyQuickCheckOutput).and.(CalParams%IStep==1).and.(CalParams%TimeStep==1)
          
          WriteTimeStep = (CalParams%OutputNumberOfSubsteps>0).or. &
                          (CalParams%OutputNumberOfTimeSteps>0).or. &
                          (CalParams%OutputNumberOfRealTime>0)

          if (.not.AnErrorHappend) then
            ! in case of error, the last time step result should be written
            if ((.not.WriteTimeStep).and.(.not.WriteQuickOutput)) RETURN
            if (.not.CalParams%OutputEachLoadStep) RETURN
          endif

          call DefineVariableData()

          call GetStepExt(CalParams%IStep, CalParams%FileNames%LoadStepExt) ! current load step stored as string
          call GetTimeStep(CalParams%TimeStep, CalParams%FileNames%TimeStepExt) ! current time step stored as string

          call CreateResultFileHeader(CalParams%FileNames%LoadStepExt)

          Filename = trim(CalParams%FileNames%ProjectName)// '_' //trim(CalParams%FileNames%TimeStepExt) &
                                                          // BRF_FILE_EXTENSION //trim(CalParams%FileNames%LoadStepExt) 

          AccumulatedDisplacementSoil(1:Counters%N) = AccumulatedDisplacementSoil(1:Counters%N) +  &
                                                         IncrementalDisplacementSoil(1:Counters%N, 1)
          PhaseDisplacementSoil = TotalDisplacementSoil

          if (WriteQuickOutput) then ! QuickCheckOutput: time step result at first time step of first load step
              if (CalParams%Visualization=='GiD-ASCII') then
                  call WriteGiDMeshASCII
                  call WriteGiDResultsASCII
              else if (CalParams%Visualization=='GiD-Binary') then
                  call WriteGiDResultsBIN()
              elseif (CalParams%Visualization=='Paraview-GiD') then
                  call WriteResultsOfOneStep(Filename)
                  call WriteGiDResultsBIN()
              else
                  call WriteResultsOfOneStep(Filename)
              endif
              
                
          endif
 
          if (CalParams%OutputNumberOfTimeSteps>0) then ! time step results at specified time steps
            do I = 1, CalParams%OutputNumberOfTimeSteps
              if (CalParams%TimeStep==CalParams%OutputTimeStepID(I)) then
                call WriteResultsOfOneStep(Filename)
              end if
            end do
          end if

          if (AnErrorHappend) then ! time step results when an error happens
            call WriteResultsOfOneStep(Filename)
          end if

          if (CalParams%OutputNumberOfRealTime>0) then ! time step results at specified real times
            do I = 1, CalParams%OutputNumberOfRealTime
              if (CalParams%TotalRealTime>=CalParams%OutputRealTimeID(I)) then
                call WriteResultsOfOneStep(Filename)
              end if
            end do
          end if
          
          if (CalParams%OutputNumberOfSubsteps>0) then ! time step results for specified number of substeps
            do I = 1, CalParams%OutputNumberOfSubsteps-1
              OutputTime(I) = I * CalParams%TotalTime / CalParams%OutputNumberOfSubsteps
              if (CalParams%TotalRealTime>=OutputTime(I)) then
                call WriteResultsOfOneStep(Filename)
              end if
            end do
          end if
          
        end subroutine WriteTimeStepResults
      
        subroutine WriteResultsOfOneStep(Filename)
        !**********************************************************************
        !
        !  Function : Writes load phase result data
        !     
        !**********************************************************************
        implicit none
        character*(*), intent(in) :: Filename

#ifdef __INTEL_COMPILER        
        call CheckFreeSpace(Filename)
#endif        
        call WriteResultData(Filename)

        if (NFORMULATION==1) then
           call WriteVTKOutput() ! VTK output
        else
           call WriteVTKOutput_2LayForm_Solid() ! VTK output Solid Material Point
           call WriteVTKOutput_2LayForm_Liquid() ! VTK output Liquid Material Point
        end if

        end subroutine WriteResultsOfOneStep

        
        subroutine CheckFreeSpace(FileName)
        !**********************************************************************
        !
        !  Function : Checks free space on the drive where FileName is going to be saved.
        !
        !**********************************************************************
        implicit none
        integer(INTEGER_TYPE), parameter :: MINIMUM_FREE_SPACE_KB = 100000 ! 100 MB
        character*(*), intent(in) :: FileName
        integer(INTEGER_TYPE) :: FreeSpaceKB

        FreeSpaceKB = kBSpace(FileName)
        if (FreeSpaceKB < MINIMUM_FREE_SPACE_KB .and. FreeSpaceKB > 0) then
          call GiveMessage('There is not free space on project drive! Available space (KB): '// trim(String(FreeSpaceKB)))
          call GiveMessage('Free up disk space and press any key to continue!')
          read(*,*)
        endif

        end subroutine CheckFreeSpace

        subroutine WriteLoadPhaseResults()
        !**********************************************************************
        !
        !  Function : Writes load phase result data
        !     
        !**********************************************************************
        
        implicit none

          ! Local variables
          character(len = MAX_FILENAME_LENGTH) :: Filename

          if (((.not.CalParams%OutputEachLoadStep).and.(CalParams%IStep/=CalParams%NLoadSteps)).or. &
              (CalParams%ApplyImplicitQuasiStatic.and.((CalParams%OutputNumberOfLoadSteps>0).and. &
               (mod(CalParams%IStep, CalParams%OutputNumberOfLoadSteps)/=0)))) RETURN 

          call GetTimeStep(CalParams%TimeStep, CalParams%FileNames%TimeStepExt) ! current time step stored as string

          PhaseDisplacementSoil(1:Counters%N) = PhaseDisplacementSoil(1:Counters%N) + IncrementalDisplacementSoil(1:Counters%N, 1)
         
          if (CalParams%ApplyImplicitQuasiStatic.and.CalParams%ImplicitIntegration%IsZeroLoadStep) then
            Filename = trim(CalParams%FileNames%ProjectName)//BRF_FILE_EXTENSION//trim(CalParams%FileNames%LoadStepExt)//'ZLS'
          else
            Filename = trim(CalParams%FileNames%ProjectName)//BRF_FILE_EXTENSION//trim(CalParams%FileNames%LoadStepExt)
          end if
          
          call CreateResultFileHeader(CalParams%FileNames%LoadStepExt)
          call WriteResultData(Filename)

          if(NFORMULATION==1) then
              if (CalParams%Visualization=='GiD-ASCII') then
                  call WriteGiDMeshASCII
                  call WriteGiDResultsASCII
              else if (CalParams%Visualization=='GiD-Binary') then
                  call WriteGiDResultsBIN()
              elseif (CalParams%Visualization=='Paraview-GiD') then
                  call WriteResultsOfOneStep(Filename)
                  call WriteGiDResultsBIN()
              else
                  call WriteResultsOfOneStep(Filename)
              endif
          else
              call WriteVTKOutput_2LayForm_Solid() ! VTK output Solid Material Point
              call WriteVTKOutput_2LayForm_Liquid() ! VTK output Liquid Material Point
          end if
         
        end subroutine WriteLoadPhaseResults


        subroutine WriteResultData(BRFfilename)
        !**********************************************************************
        !
        !    Function:  Starts the binary result file (BRF).
        !               Appends the binary result file by blocks 
        !               which contain the nodal data.
        !               Appends the binary result file by blocks
        !               which contain the particle data.
        !
        !     BRFfilename : Name of the result file
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          ! arguments
          character(len = MAX_FILENAME_LENGTH), intent(in) :: BRFfilename

          call EraseFile(BRFfilename)

          call FileOpenAction(BRFunit, BRFfilename,'W')

          write(BRFunit)  CalParams%FileNames%ResultFileHeader

          call AppendNodalData(BRFunit)
          call AppendMPMData(BRFunit)

          write(BRFunit)'$$ENDOFFII$$', 0

          endfile BRFunit

          close(BRFunit)

        end subroutine WriteResultData


        subroutine WriteFEMNodeData()
        !**********************************************************************
        !
        !  Function : Writes nodal displacement data do TST file in case of FEM calculation
        !     
        !**********************************************************************
              
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: I, J, IDoF, NColumns
          real(REAL_TYPE) :: Values(2*MAXOUTPUTNODES * NVECTOR)

          if (IsMPMComputation()) RETURN
          
          
          Values = 0.0
        
          NColumns = 0
          
          
          do I = 1, MAXOUTPUTNODES
            if (CalParams%OutputNodes(I)==-1) EXIT

            IDoF = ReducedDof(CalParams%OutputNodes(I))
            do J = 1, NVECTOR
              Values(NColumns + J) = TotalDisplacementSoil(IDoF + J)

            end do
            NColumns = NColumns + NVECTOR
            do J = 1, NVECTOR
              Values(NColumns + J) = TotalVelocitySoil(IDoF + J,1)
            end do
            NColumns = NColumns + NVECTOR
          end do

        !TST file writes as follows UX(Node1) Uy(Node1) Uz(node1) Vx(Node1) Vy(Node1) Vz(node1) Ux(Node2) Uy(Node2) Uz(Node2)
          write(TSTUnit, '(I12, I12, 60G12.4)') CalParams%IStep, CalParams%TimeStep, (Values(I), I = 1, NColumns)
        
        end subroutine WriteFEMNodeData

      !-----------------------------------------------------------------
      integer(INTEGER_TYPE) function kBSpace(S)
      ! Function: Returns free disk space (in kB) of drive
      ! Input   : S = Directory
     
      implicit none
      character S*(*)
#ifdef __INTEL_COMPILER      
      character($MAXPATH) dir
#else
      character(260) dir
#endif
      integer(4) length
      Integer(8) I8_Total, I8_Avail

      if (S(2:2)/=':') then
        !  Get current directory
#ifdef __INTEL_COMPILER        
        dir = FILE$CURDRIVE
        length = GETDRIVEDIRQQ(dir)
        if (length > 0) then
          call WriteInLogFile('Current directory is: ' // trim(dir))
        else
          call WriteInLogFile('Failed to get current directory')
        end if
#else
	status = getcwd(dir)
        if (status > 0) then
          call WriteInLogFile('Current directory is: ' // trim(dir))
        else
          call WriteInLogFile('Failed to get current directory')
        end if
#endif	

      else
        Dir = trim(S)
      end if

      kBSpace = 1
#ifdef __INTEL_COMPILER      
      if (GetDriveSizeQQ( dir, I8_Total, I8_Avail )) then
        call WriteInLogFile('Capacity, available  Byte:'// trim(String(I8_Total)) // trim(String(I8_Avail)))
        call WriteInLogFile('Capacity, available KByte:'// trim(String(I8_Total/1024d0**1)) // trim(String(I8_Avail/1024d0**1)))
        call WriteInLogFile('Capacity, available GByte:'// trim(String(I8_Total/1024d0**3)) // trim(String(I8_Avail/1024d0**3)))
        kBSpace = nint(I8_Total/1024d0**1)
        RETURN
      else
        kBSpace = -1
      end if
#endif
      
      end function kBSpace

      end module ModWriteResultData
