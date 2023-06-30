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


      module ModWriteTestData
      !**********************************************************************
      !
      !  Function : Contains routines related to writing data to files for testing the MPM code.
      !
      !             In order to keep the size of this source file reasonably small,
      !             this module only contains routines that are directly related to
      !             the output of (material point) data to text files for testing.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModParticle
      use ModReadCalculationData
      use ModMPMData
      use ModMPMDYN2PhaseSP
      use ModMPMDYN3PhaseSP
      use ModDynViscousBoudary
      use ModMPMDynContact

      implicit none

      contains
    

        subroutine OpenTextOutputFiles()
        !**********************************************************************
        !
        !  Function : Open text files for output of test data:
        !             OUT: Log file for FEM routine calls.
        !             TST: Arbitrary output data for testing purposes
        !             MLG: Log file for MPM routine calls (for initialisation and each load step)
        !             BMR: Text file for writing output results for benchmarking 
        !             BMS: text file for writing output results for benchmarking stresses and strains
        !
        !**********************************************************************
        
        implicit none

          character(len = MAX_FILENAME_LENGTH) :: ProjectName = ' ' ! name of project file
          
          call getarg(1, ProjectName)
          
          ! Open problem.OUT file (log information)
          call FileOpen(OUTUnit, trim(ProjectName)//'.OUT')

          ! Open problem.TST file for testing
          call FileOpen(TSTUnit, trim(ProjectName)//'.TST')

          ! Open problem.MLG log file for material point initialisation
          call FileOpen(LOGUnit, trim(ProjectName)//'_000.MLG')

          ! Open problem.BMR file for writing benchmark result data
          call FileOpen(BMRUnit, trim(ProjectName)//'.BMR')
          
          ! Open problem.BMR file for writing benchmark result data
          call FileOpen(BMSUnit, trim(ProjectName)//'.BMS')

        end subroutine OpenTextOutputFiles

        
        subroutine InitialiseTextOutputFiles()
        !**********************************************************************
        !
        !  Function : Initialise text files for output of test data:
        !             OUT: Log file for FEM routine calls.
        !             TST: Arbitrary output data for testing purposes
        !             MLG: Log file for MPM routine calls (for initialisation and each load step)
        !             BMR: Text file for writing output results for benchmarking 
        !             BMS: text file for writing output results for benhcmarking stresses ands strains
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          character(len = 15) :: XXX

          call GetStepExt(CalParams%FileCounter, XXX)

          if (CalParams%IStep>1) then
            CalParams%FileCounter = CalParams%FileCounter + 1
          end if
        
          call FileOpen(RXunit, trim(CalParams%FileNames%ProjectName)//'_'//trim(XXX)//'.RX')
     
        end subroutine InitialiseTextOutputFiles

        !**********************************************************************
        !
        !    FUNCTION:  InitialiseMaterialPointOutputFiles
        ! 
        !    DESCRIPTION:
        !>   Initialise text files for output of material point data
        !
        !**********************************************************************
        subroutine InitialiseMaterialPointOutputFiles()
          implicit none
        
          ! Local variables
          character(len = 8) :: XXX
          character(len = 8) :: XXXCounter
          integer(INTEGER_TYPE) :: I
          character(len = 255) :: CompleteFileName
        
          if ( (CalParams%OutputNumberParticles > 0) .and. (CalParams%OutputNumberParticles <= MAXOUTPUTPARTICLES) ) then
            ! Open files for output of particle data
            CalParams%ParticleFileCounter = CalParams%ParticleFileCounter + 1
            call GetStepExt(CalParams%ParticleFileCounter, XXXCounter)
            do I = 1, CalParams%OutputNumberParticles
              if (CalParams%OutputParticles(I)>0) then
                call GetStepExt(IDArray(CalParams%OutputParticles(I)), XXX)
                CompleteFileName = trim(CalParams%FileNames%ProjectName)// '.PAR_' // XXX 
                if ( FExist(trim(CompleteFileName)) ) then
                  call FileOpenAppend(PARUnit + I, CompleteFileName) ! adapt existing file
                else
                  call FileOpen(PARUnit + I, CompleteFileName) ! create new file
                  
                  select case(NDIM)
                  
                    case(3)
                      write(PARUnit + I, '(28A12)')  &
                            'LoadStep ', 'TimeStep ', 'Time', 'PA(1) ','PA(2) ', 'PGravity ', 'ID_MP ', 'X ', 'Y ', 'Z ', &
                            'Ux ', 'Uy ', 'Uz ', 'SigmaXX ', 'SigmaYY ', 'SigmaZZ ', 'SigmaXY ', 'SigmaYZ ', 'SigmaZX ', &
                            'WPressure', 'EpsilonXX ', 'EpsilonYY ', 'EpsilonZZ ', 'GammaXY ', 'GammaYZ ', 'GammaZX ', &
                            'IntWeight ', 'MatID '
                    case(2)
                      write(PARUnit + I, '(22A12)')  &
                            'LoadStep ', 'TimeStep ', 'Time', 'PA(1) ','PA(2) ', 'PGravity ', 'ID_MP ', 'X ', 'Y ', &
                            'Ux ', 'Uy ', 'SigmaXX ', 'SigmaYY ', 'SigmaZZ ','SigmaXY ', &
                            'WPressure', 'EpsilonXX ', 'EpsilonYY ', 'EpsilonZZ ','GammaXY ', &
                            'IntWeight ', 'MatID '
                  end select

                end if  
              end if
            end do
          end if

        end subroutine InitialiseMaterialPointOutputFiles
        
        subroutine InitialiseSurfaceReactionOutputFiles()
        !**********************************************************************
        !
        !  Function : Initialise text files for output of surface reaction forces
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          character(len = 8) :: XXX
          integer(INTEGER_TYPE) :: ISurface, MaterialIndex
          logical :: IsUndrEffectiveStress
          character(len = 255) :: CompleteFileName
          
          if(Counters%NReactionSurfaceOutput<1) RETURN
          IsUndrEffectiveStress = .false.
          
            do MaterialIndex=1, CalParams%NumberOfMaterials
              if (.not.IsUndrEffectiveStress) then
                IsUndrEffectiveStress = &
            (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE)
              end if
             end do
  

            do ISurface = 1, Counters%NReactionSurfaceOutput
                call GetStepExt(ISurface, XXX)
                CompleteFileName = trim(CalParams%FileNames%ProjectName)// '.RSurf_' // XXX 
                if ( (FExist(trim(CompleteFileName))).and.(CalParams%IStep>1) ) then
                  call FileOpenAppend(SURFReacUnit + ISurface, CompleteFileName) ! adapt existing file
                else
                  call FileOpen(SURFReacUnit + ISurface, CompleteFileName) ! create new file
                  write(SURFReacUnit + ISurface, '(2A)') 'Surface name:', OutputSurfaceName(ISurface)
                  if (CalParams%NumberOfMaterials==1) then
                    if (CalParams%NumberOfPhases==3) then
                      if (NDIM == 3) then ! 3D case 
                        write(SURFReacUnit + ISurface, '(12A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', 'SumReactionSolidZ_1 ', &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1', 'SumReactionWaterZ_1 ', &
                                   'SumReactionGasX_1', 'SumReactionGasY_1', 'SumReactionGasZ_1'
                      elseif (NDIM == 2) then
                          write(SURFReacUnit + ISurface, '(12A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1', &
                                   'SumReactionGasX_1', 'SumReactionGasY_1'
                          
                      end if 
                    else if ((CalParams%NumberOfPhases==2).or.IsUndrEffectiveStress) then
                      if (NDIM == 3) then ! 3D case   
                        write(SURFReacUnit + ISurface, '(9A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', 'SumReactionSolidZ_1 ', &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1', 'SumReactionWaterZ_1 '
                      elseif (NDIM == 2) then
                        write(SURFReacUnit + ISurface, '(12A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                              'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1'
                      end if ! add 2D in sprint#2  
                    else
                      if (NDIM == 3) then ! 3D case   
                        write(SURFReacUnit + ISurface, '(6A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', 'SumReactionSolidZ_1 '
                      elseif (NDIM == 2) then
                        write(SURFReacUnit + ISurface, '(5A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                            'SumReactionSolidX_1 ', 'SumReactionSolidY_1 '
                      end if ! add 2D in sprint#2  
                    end if !number of phases
                    else if (CalParams%NumberOfMaterials==2) then
                    if (CalParams%NumberOfPhases==3) then
                      if (NDIM == 3) then ! 3D case   
                        write(SURFReacUnit + ISurface, '(21A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', 'SumReactionSolidZ_1 ', &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ', 'SumReactionSolidZ_2 ', &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1', 'SumReactionWaterZ_1 ', &
                                   'SumReactionWaterX_2 ', 'SumReactionWaterY_2', 'SumReactionWaterZ_2 ', &
                                   'SumReactionGasX_1', 'SumReactionGasY_1', 'SumReactionGasZ_1', &
                                   'SumReactionGasX_2', 'SumReactionGasY_2', 'SumReactionGasZ_2'
                      elseif (NDIM == 2) then
                        write(SURFReacUnit + ISurface, '(15A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ',  &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1', &
                                   'SumReactionWaterX_2 ', 'SumReactionWaterY_2', &
                                   'SumReactionGasX_1', 'SumReactionGasY_1', &
                                   'SumReactionGasX_2', 'SumReactionGasY_2'
                      end if  ! add 2D in sprint#2  
                    else if ((CalParams%NumberOfPhases==2).or.IsUndrEffectiveStress) then
                      if (NDIM == 3) then ! 3D case   
                        write(SURFReacUnit + ISurface, '(15A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1', 'SumReactionSolidY_1', 'SumReactionSolidZ_1', &
                                   'SumReactionSolidX_2', 'SumReactionSolidY_2', 'SumReactionSolidZ_2', &
                                   'SumReactionWaterX_1', 'SumReactionWaterY_1', 'SumReactionWaterZ_1', &                                   
                                   'SumReactionWaterX_2', 'SumReactionWaterY_2', 'SumReactionWaterZ_2'
                      elseif (NDIM == 2) then
                        write(SURFReacUnit + ISurface, '(11A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1', 'SumReactionSolidY_1', &
                                   'SumReactionSolidX_2', 'SumReactionSolidY_2', &
                                   'SumReactionWaterX_1', 'SumReactionWaterY_1', &                                   
                                   'SumReactionWaterX_2', 'SumReactionWaterY_2'
                      end if  ! add 2D in sprint#2  
                    else
                      if (NDIM == 3) then ! 3D case   
                        write(SURFReacUnit + ISurface, '(9A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', 'SumReactionSolidZ_1 ', &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ', 'SumReactionSolidZ_2 '
                      elseif (NDIM == 2) then
                        write(SURFReacUnit + ISurface, '(7A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                            'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', &
                             'SumReactionSolidX_2 ', 'SumReactionSolidY_2 '
                      end if  ! add 2D in sprint#2  
                    end if   
                        
                    else if (CalParams%NumberOfMaterials==3) then
                    if (CalParams%NumberOfPhases==3) then
                      if (NDIM == 3) then ! 3D case    
                        write(SURFReacUnit + ISurface, '(30A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', 'SumReactionSolidZ_1 ', &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ', 'SumReactionSolidZ_2 ', &
                                   'SumReactionSolidX_3 ', 'SumReactionSolidY_3', 'SumReactionSolidZ_3 ', &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1', 'SumReactionWaterZ_1 ', &
                                   'SumReactionWaterX_2 ', 'SumReactionWaterY_2', 'SumReactionWaterZ_2 ', &
                                   'SumReactionWaterX_3 ', 'SumReactionWaterY_3', 'SumReactionWaterZ_3 ', &
                                   'SumReactionGasX_1', 'SumReactionGasY_1', 'SumReactionGasZ_1', &
                                   'SumReactionGasX_2', 'SumReactionGasY_2', 'SumReactionGasZ_2', &
                                   'SumReactionGasX_3', 'SumReactionGasY_3', 'SumReactionGasZ_3'
                      else if (NDIM == 2) then
                        write(SURFReacUnit + ISurface, '(21A21)') 'Step ', 'TimeStep ','TotalTime ',  & 
                                  'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ',  &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ', &
                                   'SumReactionSolidX_3 ', 'SumReactionSolidY_3',  &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1', &
                                   'SumReactionWaterX_2 ', 'SumReactionWaterY_2', &
                                   'SumReactionWaterX_3 ', 'SumReactionWaterY_3',  &
                                   'SumReactionGasX_1', 'SumReactionGasY_1',  &
                                   'SumReactionGasX_2', 'SumReactionGasY_2',  &
                                   'SumReactionGasX_3', 'SumReactionGasY_3'
                      end if ! add 2D in sprint#2    
                    else if ((CalParams%NumberOfPhases==2).or.IsUndrEffectiveStress) then
                      if (NDIM == 3) then ! 3D case    
                        write(SURFReacUnit + ISurface, '(21A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', 'SumReactionSolidZ_1 ', &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ', 'SumReactionSolidZ_2 ', &
                                   'SumReactionSolidX_3 ', 'SumReactionSolidY_3', 'SumReactionSolidZ_3 ', &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1', 'SumReactionWaterZ_1 ', &
                                   'SumReactionWaterX_2 ', 'SumReactionWaterY_2', 'SumReactionWaterZ_2 ', &
                                   'SumReactionWaterX_3 ', 'SumReactionWaterY_3', 'SumReactionWaterZ_3 '
                      elseif (NDIM == 2) then
                        write(SURFReacUnit + ISurface, '(21A21)') 'Step ', 'TimeStep ','TotalTime ',  & 
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ',  &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ',  &
                                   'SumReactionSolidX_3 ', 'SumReactionSolidY_3',  &
                                   'SumReactionWaterX_1 ', 'SumReactionWaterY_1',  &
                                   'SumReactionWaterX_2 ', 'SumReactionWaterY_2',  &
                                   'SumReactionWaterX_3 ', 'SumReactionWaterY_3'
                      end if     
                    else
                      if (NDIM == 3) then ! 3D case    
                        write(SURFReacUnit + ISurface, '(12A21)') 'Step ', 'TimeStep ','TotalTime ',  &
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', 'SumReactionSolidZ_1 ', &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ', 'SumReactionSolidZ_2 ', &
                                   'SumReactionSolidX_3 ', 'SumReactionSolidY_3', 'SumReactionSolidZ_3 '
                      elseif (NDIM == 2) then
                        write(SURFReacUnit + ISurface, '(21A21)') 'Step ', 'TimeStep ','TotalTime ',  &  
                                   'SumReactionSolidX_1 ', 'SumReactionSolidY_1 ', &
                                   'SumReactionSolidX_2 ', 'SumReactionSolidY_2 ', &
                                   'SumReactionSolidX_3 ', 'SumReactionSolidY_3'
                      end if    
                    end if !number of phases
                    end if !number of materials
                        
                end if  
            end do

        end subroutine InitialiseSurfaceReactionOutputFiles

        
        subroutine InitialiseLoadPhaseOutputFiles()
        !**********************************************************************
        !
        !  Function : Opens the text output files for the current load phase
        !
        !**********************************************************************

        implicit none

          call InitialiseMLGFile()
          call InitialiseENGFile()
          call InitialiseINFFile()
          if (CalParams%ApplyEmptyElements) call InitialiseTVMFile()


        end subroutine InitialiseLoadPhaseOutputFiles


        subroutine InitialiseMLGFile()
        !**********************************************************************
        !
        !  Function : Close MLG file of previous load step and open new 
        !             problem_XXX.MLG file where XXX stands for the current load step number
        !
        !**********************************************************************
        
        implicit none

          call CloseFile(LOGunit)
          call FileOpen(LOGunit, trim(CalParams%FileNames%ProjectName)//'_'//trim(CalParams%FileNames%LoadStepExt)//'.MLG')
        
        end subroutine InitialiseMLGFile
        

        subroutine InitialiseENGFile()
        !**********************************************************************
        !
        !  Function : Close ENG file of previous load step and open new 
        !             problem_XXX.ENG file where XXX stands for the current load step number
        !
        !**********************************************************************
        
        implicit none

          call CloseFile(ENGunit)
          call FileOpen(ENGunit, trim(CalParams%FileNames%ProjectName)//'_'//trim(CalParams%FileNames%LoadStepExt)//'.ENG')

        end subroutine InitialiseENGFile


        subroutine InitialiseINFFile()
        !**********************************************************************
        !
        !  Function : Close INF file of previous load step and open new 
        !             problem_XXX.INF file where XXX stands for the current load step number
        !
        !**********************************************************************
        
        implicit none

          call CloseFile(INFunit)
          call FileOpen(INFunit, trim(CalParams%FileNames%ProjectName)//'_'//trim(CalParams%FileNames%LoadStepExt)//'.INF')

        end subroutine InitialiseINFFile

        subroutine InitialiseTVMFile()
        !**********************************************************************
        !
        !  Function : Close TVM file of previous load step and open new 
        !             problem_XXX.TVM file where XXX stands for the current load step number
        !
        !**********************************************************************
        implicit none

          call CloseFile(TVMunit)
          call FileOpen(TVMunit, trim(CalParams%FileNames%ProjectName)//'_'//trim(CalParams%FileNames%LoadStepExt)//'.TVM')

          ! write header
          write(TVMunit, '(26A20)') 'Step ', 'TotalMass ', 'AddedMass ', 'NParticles ', 'NAddedParticles ', 'NRemovedParticles '

        end subroutine InitialiseTVMFile




        subroutine MaterialPointOutput()
        !**********************************************************************
        !
        !  Function : Output of material point data to a text file
        !
        !**********************************************************************

        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: I, J, ParticleIndex, ParticleID
          
          if ( .not.( mod(CalParams%TimeStep, CalParams%OutputCurvesIntervals)==0 .or. CalParams%TimeStep==1 &
                   .or. IsLastTimeStepOfExplicitCalculation() ) ) RETURN

          if ( (CalParams%OutputNumberParticles>0) .and. (CalParams%OutputNumberParticles<=MAXOUTPUTPARTICLES) ) then
            do I = 1, CalParams%OutputNumberParticles

              if (CalParams%OutputParticles(I)>0) then
               if (CalParams%ApplyExcavation) then                 
                do ParticleIndex= 1, Counters%NParticles
                    ParticleID = IDArray(ParticleIndex)
                    if (ParticleID==CalParams%OutputParticles(I)) then
                        exit
                    end if
                end do
               else
                ParticleIndex = CalParams%OutputParticles(I)
               end if 
                if (ParticleIndex<=Counters%NParticles) then
                    write(PARUnit + I, '(I12, I12, 4G12.4, I12, 20G12.4, I12)') &
                      CalParams%IStep, &
                      CalParams%TimeStep, &
                      CalParams%OverallRealTime, &
                      CalParams%Multipliers%SolidACurrent,  &
                      CalParams%Multipliers%GravityCurrent,  &
                      IDArray(ParticleIndex), &
                   
                      (GlobPosArray(ParticleIndex,J), J=1,NVECTOR), & 
                    
                      (UArray(ParticleIndex,J), J=1,NVECTOR), & 
                      
                      (SigmaEffArray(ParticleIndex,J), J=1,NTENSOR), & 
                      
                      Particles(ParticleIndex)%WaterPressure, &

                      (GetEpsI(Particles(ParticleIndex),J), J=1,NTENSOR), & 
                        
                      Particles(ParticleIndex)%IntegrationWeight, &
                      MaterialIDArray(ParticleIndex)
                end if
              end if
            end do
          end if

        end subroutine MaterialPointOutput

        subroutine EnergyOutput()
        !**********************************************************************
        !
        !  Function : Output of material point energy data to a text file
        !
        !**********************************************************************

        implicit none
  
          if (CalParams%TimeStep==1) then  ! write header
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
              write(ENGUnit, '(27A20)')  &
                'LoadStep', &
                'TimeStep', &
                'KineticEnergy',  &
                'InternalWork', &
                'ExternalWork', &
                'Dissipation', &
                'KineticError', &
                'ForceError', &
                'KineticEnergySoil', &
                'InternalWorkSoil', &
                'ExternalWorkSoil', &
                'DissipationSoil', &
                'KineticErrorSoil', &
                'ForceErrorSoil', &
                'KineticEnergyLiquid', &
                'InternalWorkLiquid', &
                'ExternalWorkLiquid', &
                'DissipationLiquid', &
                'KineticErrorLiquid', &
                'ForceErrorLiquid', &
                'TimeInc', &
                'RealTime'
            else
              write(ENGUnit, '(27A20)')  &
                'LoadStep',  &
                'TimeStep', &
                'KineticEnergy', &
                'InternalWork', &
                'ExternalWork', &
                'Dissipation', & 
                'KineticError', &
                'ForceError', &
                'TimeInc', &
                'RealTime'
            end if
          end if
          
          if((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then  ! write data
            write(ENGUnit, '(I20, I20, 25G20.7)') &
              CalParams%IStep,  &
              CalParams%TimeStep, &
              CalParams%ConvergenceCheck%KineticEnergy, &
              CalParams%ConvergenceCheck%InternalWork, &
              CalParams%ConvergenceCheck%ExternalWork, &
              CalParams%ConvergenceCheck%EnergyDissipation, &
              CalParams%ConvergenceCheck%KineticError,  &
              CalParams%ConvergenceCheck%ForceError, &
              CalParams%ConvergenceCheck%KineticEnergySoil, &
              CalParams%ConvergenceCheck%InternalWorkSoil, &
              CalParams%ConvergenceCheck%ExternalWorkSoil, &
              CalParams%ConvergenceCheck%EnergyDissipationSoil, &
              CalParams%ConvergenceCheck%KineticErrorSoil,  &
              CalParams%ConvergenceCheck%ForceErrorSoil, &
              CalParams%ConvergenceCheck%KineticEnergyWater, &
              CalParams%ConvergenceCheck%InternalWorkWater, &
              CalParams%ConvergenceCheck%ExternalWorkWater, &
              CalParams%ConvergenceCheck%EnergyDissipationWater, &
              CalParams%ConvergenceCheck%KineticErrorWater,  &
              CalParams%ConvergenceCheck%ForceErrorWater, &
              CalParams%TimeIncrement,  &
              CalParams%TotalRealTime
          else
            write(ENGUnit, '(I20, I20, 25G20.7)')  &
              CalParams%IStep,  &
              CalParams%TimeStep, &
              CalParams%ConvergenceCheck%KineticEnergy, &
              CalParams%ConvergenceCheck%InternalWork, &
              CalParams%ConvergenceCheck%ExternalWork,  &
              CalParams%ConvergenceCheck%EnergyDissipation, &
              CalParams%ConvergenceCheck%KineticError,  &
              CalParams%ConvergenceCheck%ForceError, &
              CalParams%TimeIncrement,  &
              CalParams%TotalRealTime
          end if
          
        end subroutine EnergyOutput
        
        
        subroutine EnergyOutput2LayForm()
        !**********************************************************************
        !
        !  Function : Output of material point energy data to a text file for 2 Layer formulation
        !
        !**********************************************************************

        implicit none
  
          if (CalParams%TimeStep==1) then  ! write header
              write(ENGUnit, '(27A20)')  &
                'LoadStep', &
                'TimeStep',  &
                'SOLID_KinEnergy',  &
                'SOLID_IntWork', &
                'SOLID_ExtWork',  &
                'LIQUID_KinEnergy', &
                'LIQUID_IntWork', &
                'LIQUID_ExtWork',  &
                'SOLID_KinError', &
                'SOLID_ForceError', &
                'LIQUID_KinError', &
                'LIQUID_ForceError', &
                'TimeInc', &
                'RealTime'
          end if
          
          write(ENGUnit, '(I20, I20, 25G20.7)') &
              CalParams%IStep,  &
              CalParams%TimeStep, &
              CalParams%ConvergenceCheck%KineticEnergySoil, &
              CalParams%ConvergenceCheck%InternalWorkSoil, &
              CalParams%ConvergenceCheck%ExternalWorkSoil, &
              CalParams%ConvergenceCheck%KineticEnergyWater, &
              CalParams%ConvergenceCheck%InternalWorkWater, &
              CalParams%ConvergenceCheck%ExternalWorkWater,  &
              CalParams%ConvergenceCheck%KineticErrorSoil,  &
              CalParams%ConvergenceCheck%ForceErrorSoil, &
              CalParams%ConvergenceCheck%KineticErrorWater,  &
              CalParams%ConvergenceCheck%ForceErrorWater, &
              CalParams%TimeIncrement,  &
              CalParams%TotalRealTime
          
        end subroutine EnergyOutput2LayForm

        
        subroutine SurfaceReactionOutput()
        !*********************************************************************
        !
        ! Function: compute and write on txt file the sum of reaction forces
        !           on selected surfaces
        !
        !*********************************************************************
         implicit none
         real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials) :: ReactionSolid, ReactionWater, ReactionGas
         real(REAL_TYPE), dimension(NVECTOR, CalParams%NumberOfMaterials, Counters%NReactionSurfaceOutput ) :: SumReactionSolid, SumReactionWater, SumReactionGas
         
         if(Counters%NReactionSurfaceOutput<1) RETURN
         
         call ComputeNodalReactions(ReactionSolid, ReactionWater, ReactionGas)
         call SumNodalReactionsOnSurfaces(ReactionSolid, ReactionWater, ReactionGas, SumReactionSolid, SumReactionWater, SumReactionGas)
         call WriteReactionsOnFile(SumReactionSolid, SumReactionWater, SumReactionGas)
         
        end subroutine SurfaceReactionOutput
        
        subroutine ComputeNodalReactions(ReactionSolid, ReactionWater, ReactionGas)
        !*********************************************************************
        !
        ! Function: compute nodal reactions
        !
        !*********************************************************************
        implicit none
        
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials), intent(inout) :: ReactionSolid
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials), intent(inout) :: ReactionWater
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials), intent(inout) :: ReactionGas
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials) :: GravitySolid
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials) :: GravityWater
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials) :: GravityGas
        
        integer(INTEGER_TYPE) :: I, J, IEl, NElemPart, Int, INtGlo, MaterialIndex, Inode, nn, nix
        real(REAL_TYPE):: VPressure, WTN, WPP, GPP, Det
        real(REAL_TYPE) :: Position
        real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues
        real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
        real(REAL_TYPE), dimension(NTENSOR) :: S
        real(REAL_TYPE), dimension(NVECTOR) :: PartGravity, PartGravityWater, PartGravityGas
        logical :: IsUndrEffectiveStress
        
        if((CalParams%TimeStep==1)) then
          call DetermineReactionElements
        end if
        
        ReactionSolid = 0.0
        ReactionWater = 0.0
        ReactionGas = 0.0
        
        !compute ReactionForce
        do I=1, Counters%NElemReactions !loop over the elements that have a node on the reaction surface
          IEl = ConsideredElemReaction(I)
          
          if (.not.IsActiveElement(IEl)) CYCLE
            
          if (IsParticleIntegration(IEl) ) then ! True - material point based integration, false - Gauss point based integration
             NElemPart = NPartEle(IEl)  ! Number of material points in element
          else
             NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
          end if
             
          call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element
        
          do Int = 1, NElemPart ! Loop over number of integration points per element IEl
              
            ! Determine global ID of integration point 
            IntGlo = GetParticleIndex(Int, IEl)
            MaterialIndex = MaterialIDArray(IntGlo)
            WPP = 0.0
            GPP = 0.0
            S = 0.0
            VPressure = 0.0
            ShapeValues = 0.0
            position = 0.0
            
              
            IsUndrEffectiveStress = (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE)
              
            ! Set the integration weight
            if (IsParticleIntegration(IEl)) then 
             if ( ISAXISYMMETRIC ) then
             ! the integration weight of the MP is not corrected due to the axisymmetry
               Position = GlobPosArray(IntGlo, 1) ! INDEX_X = 1
               ShapeValues(:) = ShapeValuesArray(IntGlo, :)
             end if
            ! use the integration weight of material point if it is partially filled, otherwise the weight of gauss point is used
             WTN = Particles(IntGlo)%IntegrationWeight
            else
              if ( ISAXISYMMETRIC ) then
                Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                ShapeValues(:) = GPShapeFunction(Int, :)
                WTN = WTN * Position ! the volume of the element is corrected due to the axisymmetry
              end if
            end if

            ! Determine stress vector for integration point (material point or Gauss point)
            if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) .or.IsUndrEffectiveStress) then
              WPP = Particles(IntGlo)%WaterPressure * WTN 
            end if
            if (CalParams%NumberOfPhases==3) then
              GPP = Particles(IntGlo)%GasPressure * WTN 
            end if
            
           if ((NFORMULATION == 1).or. &
               ((NFORMULATION == 2).and.(MaterialPointTypeArray(IntGlo).ne.MaterialPointTypeLiquid))) then !for liquid material points SigmaEff=WPP, do not sum twice!
            do J = 1, NTENSOR
              S(J) = SigmaEffArray(IntGlo, J) * WTN
            end do
           end if 

            if (CalParams%ApplyBulkViscosityDamping) then
              VPressure = Particles(IntGlo)%DBulkViscousPressure * WTN
            end if
                            
            do iNode=1,ELEMENTNODES ! loop over element nodes
              nn = ElementConnectivities(iNode,iel) ! get global node number
              nix = reducedDof(nn)
                
              if (IsReactionNode(nn)) then
               if ( NVECTOR == 2 ) then ! 2D
                 ! nodal x-load
                 ReactionSolid(nix+1, MaterialIndex) = ReactionSolid(nix+1, MaterialIndex) + B(1,iNode)*S(1) + B(2,iNode)*S(4) 
                 ! nodal y-load 
                 ReactionSolid(nix+2, MaterialIndex) = ReactionSolid(nix+2, MaterialIndex) + B(1,iNode)*S(4) + B(2,iNode)*S(2) 
                 ! nodal r-load (x-load)
                 if ( ISAXISYMMETRIC ) then
                   ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912 
                   ReactionSolid(nix+1, MaterialIndex) = ReactionSolid(nix+1, MaterialIndex) + S(3) * ShapeValues(iNode) / Position
                 end if
                elseif (NVECTOR == 3) then ! 3D case  
                 ReactionSolid(nix + 1, MaterialIndex) = ReactionSolid(nix + 1, MaterialIndex) + B(1, INode) * S(1) + B(2, INode) * S(4) + B(3, INode) * S(6) + B(1, INode) * VPressure
                 ReactionSolid(nix + 2, MaterialIndex) = ReactionSolid(nix + 2, MaterialIndex) + B(1, INode) * S(4) + B(2, INode) * S(2) + B(3, INode) * S(5) + B(2, INode) * VPressure
                 ReactionSolid(nix + 3, MaterialIndex) = ReactionSolid(nix + 3, MaterialIndex) + B(1, INode) * S(6) + B(2, INode) * S(5) + B(3, INode) * S(3) + B(3, INode) * VPressure
                else
                 call GiveError('Dimension is not correct. It must be 2 or 3.')
                end if 
                  
                if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) .or.IsUndrEffectiveStress) then
                  do J = 1, NVECTOR  
                    ReactionWater(nix + J, MaterialIndex) = ReactionWater(nix + J, MaterialIndex) + B(J, INode) * WPP
                  end do 
                end if
                     
                if (CalParams%NumberOfPhases==3) then
                  do J = 1, NVECTOR    
                    ReactionGas(nix + J, MaterialIndex) = ReactionGas(nix + J, MaterialIndex) + B(J, INode) * GPP
                  end do  
                end if
               end if !IsReactionNode
              end do !loop over element nodes
            end do !loop over particles
        end do !loop over elements
        
        !compute gravity load
        GravitySolid = 0.0
        GravityWater = 0.0
        GravityGas = 0.0
        
        do I=1, Counters%NElemReactions !loop over the elements that have a node on the reaction surface
          IEl = ConsideredElemReaction(I)
          NElemPart = NPartEle(IEl)
          PartGravity = 0.0
          PartGravityWater = 0.0
          PartGravityGas = 0.0
          do Int = 1, NElemPart ! Loop over number of integration points per element IEl
            ! Determine global ID of integration point 
            IntGlo = GetParticleIndex(Int, IEl)
            MaterialIndex = MaterialIDArray(IntGlo)
            if (CalParams%ApplySubmergedCalculation) then
             PartGravity = Particles(IntGlo)%FBodyMixed * CalParams%Multipliers%GravityCurrent
            else
             PartGravity = Particles(IntGlo)%FBody * CalParams%Multipliers%GravityCurrent
            end if
            
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3).or.IsUndrEffectiveStress) then
             PartGravityWater = Particles(IntGlo)%FBodyWater * CalParams%Multipliers%GravityCurrent
            end if
            if (CalParams%NumberOfPhases==3) then
             PartGravityGas = Particles(IntGlo)%FBodyGas * CalParams%Multipliers%GravityCurrent
            end if
            
            do iNode=1,ELEMENTNODES ! loop over element nodes
              nn=ElementConnectivities(iNode,iel) ! get global node number
              nix = reducedDof(nn)
              
              if (IsReactionNode(nn)) then
              PartGravity = ShapeValuesArray(IntGlo,INode) * PartGravity
    
              do J = 1, NVECTOR
                GravitySolid(nix + J, MaterialIndex) = GravitySolid(nix + J, MaterialIndex) + PartGravity(J)
              end do  
              
              if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
                PartGravityWater = ShapeValuesArray(IntGlo,INode) * PartGravityWater
                do J = 1, NVECTOR  
                  GravityWater(nix + J, MaterialIndex) = GravityWater(nix + J, MaterialIndex) + PartGravityWater(J)
                end do  
              end if
              if (CalParams%NumberOfPhases==3) then
                PartGravityGas = ShapeValuesArray(IntGlo,INode) * PartGravityGas
                do J = 1, NVECTOR  
                  GravityGas(nix + J, MaterialIndex) = GravityGas(nix + J, MaterialIndex) + PartGravityGas(J)
                end do  
              end if
             end if !is reaction node
            end do !loope over element nodes
           end do !loop over particles
        end do !loop over elements 
        
        ReactionSolid = ReactionSolid - GravitySolid

        if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
          ReactionWater = ReactionWater - GravityWater
        end if
        if (CalParams%NumberOfPhases==3) then
          ReactionGas = ReactionGas - GravityGas
        end if  
    
        end subroutine ComputeNodalReactions
        
        subroutine SumNodalReactionsOnSurfaces(ReactionSolid, ReactionWater, ReactionGas, SumReactionSolid, SumReactionWater, SumReactionGas)
        !**************************************************************************************************
        !
        ! Function: Sum up the reaction forces of the nodes belonging to a predefined surface
        !
        !**************************************************************************************************
        implicit none        
        
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials), intent(in) :: ReactionSolid
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials), intent(in) :: ReactionWater
        real(REAL_TYPE), dimension(Counters%N, CalParams%NumberOfMaterials), intent(in) :: ReactionGas

        real(REAL_TYPE), dimension(NVECTOR, CalParams%NumberOfMaterials, Counters%NReactionSurfaceOutput ), &
            intent(out) :: SumReactionSolid, SumReactionWater, SumReactionGas
        
        integer(INTEGER_TYPE) :: ISurface, MaterialID, INode, dof, IDim
        logical :: IsUndrEffectiveStress
        
        SumReactionSolid = 0.0
        SumReactionWater = 0.0
        SumReactionGas = 0.0
        
        do ISurface = 1, Counters%NReactionSurfaceOutput
         do MaterialID = 1, CalParams%NumberOfMaterials
             IsUndrEffectiveStress = &
            !code version 2016 and previous
            ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialID)%MaterialType)=='2-phase')) .or. &
            !code version 2017.1 and following
            (trim(MatParams(MaterialID)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
           
             do INode = 1, Counters%NodTot
              If (IsReactionNodeSurface(INode,ISurface)) then
               dof = ReducedDof(INode)
               do IDim = 1, NVECTOR
                 SumReactionSolid(IDim,MaterialID,ISurface) = SumReactionSolid(IDim,MaterialID,ISurface) + ReactionSolid(dof+IDim, MaterialID)
               
                 if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) .or.IsUndrEffectiveStress) then
                   SumReactionWater(IDim,MaterialID,ISurface) = SumReactionWater(IDim,MaterialID,ISurface) + ReactionWater(dof+IDim, MaterialID)
                 end if
                 
                 if (CalParams%NumberOfPhases==3) then
                   SumReactionGas(IDim,MaterialID,ISurface) = SumReactionGas(IDim,MaterialID,ISurface) + ReactionGas(dof+IDim, MaterialID)
                 end if
                 
                end do
               end if !is surface node
             end do !node loop
         end do !Materials
         end do !Surfaces                    
                   
        
        
    end subroutine SumNodalreactionsOnSurfaces
    
    
    subroutine WriteReactionsOnFile(SumReactionSolid, SumReactionWater,SumReactionGas)
    !*****************************************************************************************
    !
    !Function: write the sum of reaction forces on a file
    !
    !*****************************************************************************************
    implicit none
    real(REAL_TYPE), dimension(NVECTOR, CalParams%NumberOfMaterials, Counters%NReactionSurfaceOutput):: &
       SumReactionSolid, SumReactionWater, SumReactionGas
    integer(INTEGER_TYPE) :: ISurface, MaterialId, IDim
    logical :: IsUndrEffectiveStress = .false.
    
    do MaterialID=1,CalParams%NumberOfMaterials
       if (.not.IsUndrEffectiveStress) then 
        IsUndrEffectiveStress = &
            !code version 2016 and previous
            ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialID)%MaterialType)=='2-phase')) .or. &
            !code version 2017.1 and following
            (trim(MatParams(MaterialID)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
        end if
    end do
    
    do ISurface = 1, Counters%NReactionSurfaceOutput
        if (CalParams%NumberOfPhases==3) then
          write(SurfReacUnit+ISurface, '(I21, I21, 200G21.4)') CalParams%IStep, CalParams%TimeStep, &
              CalParams%OverallRealTime,  &
            ((SumReactionSolid(IDim, MaterialID, ISurface), &
             IDim=1, NVECTOR), MaterialID=1,CalParams%NumberOfMaterials ), &
            ((SumReactionWater(IDim, MaterialID, ISurface), &
              IDim=1, NVECTOR), MaterialID=1,CalParams%NumberOfMaterials), &
            ((SumReactionGas(IDim, MaterialID, ISurface), &
              IDim=1, NVECTOR), MaterialID=1,CalParams%NumberOfMaterials)
        else if((CalParams%NumberOfPhases==2).or.IsUndrEffectiveStress) then
          write(SurfReacUnit+ISurface, '(I21, I21, 200G21.4)') CalParams%IStep, CalParams%TimeStep, &
              CalParams%OverallRealTime,  &
            ((SumReactionSolid(IDim, MaterialID, ISurface), &
              IDim=1, NVECTOR), MaterialID=1,CalParams%NumberOfMaterials), &
            ((SumReactionWater(IDim, MaterialID, ISurface), &
              IDim=1, NVECTOR), MaterialID=1,CalParams%NumberOfMaterials)
        else
          write(SurfReacUnit+ISurface, '(I21, I21, 200G21.4)') CalParams%IStep, CalParams%TimeStep, &
              CalParams%OverallRealTime,  &
            ((SumReactionSolid(IDim, MaterialID, ISurface), &
              IDim=1, NVECTOR), MaterialID=1,CalParams%NumberOfMaterials)
        end If
    end do   

            
    end subroutine WriteReactionsOnFile
    


        subroutine WriteLoadStepInformation()
        !**********************************************************************
        !
        !  Function : Output of calculation information for each load step to a text file project_XXX.INF
        !
        !**********************************************************************

        implicit none

          character(len=30), dimension(:), allocatable :: DataLabel
          character(len=50), dimension(:), allocatable :: DataValue
          integer(INTEGER_TYPE) :: I, Size

          Size=50
          allocate(DataLabel(Size),DataValue(Size))
          DataLabel=''
          DataValue=''

          DataLabel(1) = 'Project name' 
          DataValue(1) = Trim(CalParams%FileNames%ProjectName)
          DataLabel(2) = 'Start date of calculation'
          DataValue(2) = CalParams%StartDate
          DataLabel(3) = 'Start time of calculation'
          DataValue(3) = CalParams%StartTime
          DataLabel(4) = 'End date of calculation'
          DataValue(4) = CalParams%EndDate
          DataLabel(5) = 'End time of calculation'
          DataValue(5) = CalParams%EndTime
          DataLabel(6) = 'Current load step' 
          DataValue(6) = IntToStr(CalParams%IStep)
          DataLabel(7) = 'Number of elements' 
          DataValue(7) = IntToStr(Counters%NEl)
          DataLabel(8) = 'Number of nodes' 
          DataValue(8) = IntToStr(Counters%NodTot)
          DataLabel(9) = 'Number of degrees of freedom' 
          DataValue(9) = IntToStr(Counters%N)
          DataLabel(10) = 'Number of entities' 
          DataValue(10) = IntToStr(Counters%NEntity)
          DataLabel(11) = 'Number of active elements' 
          DataValue(11) = IntToStr(Counters%NAEl)
          DataLabel(12) = 'Number of material points' 
          DataValue(12) = IntToStr(Counters%NParticles)
          DataLabel(13) = 'Number of solid material points' 
          DataValue(13) = IntToStr(Counters%SolidMaterialPoints)
          DataLabel(14) = 'Number of liquid material points' 
          DataValue(14) = IntToStr(Counters%LiquidMaterialPoints)
          DataLabel(15) = 'Number of phases' 
          DataValue(15) = IntToStr(CalParams%NumberOfPhases)
          
          
          do I = 1,Size
            write(INFunit,'(A30,A50)')DataLabel(i),DataValue(i)
          end do

          write(INFunit, *)'---end-of-file-----------------------------------------------------------------'
          
        end subroutine WriteLoadStepInformation

        subroutine WriteBenchmarkResults()
        !**********************************************************************
        !
        !  Function : Output of calculation results for benchmarking 
        !             for each load step to a text file project_XXX.BMR
        !
        !**********************************************************************

        implicit none

          character(len=40), dimension(:), allocatable :: DataLabel
          character(len=50), dimension(:), allocatable :: DataValue
          integer(INTEGER_TYPE) :: I, J, Size, ParticleIndex, ParticleID

          if (CalParams%IStep==1) then  ! write file header for first load step
            Size=9
            allocate(DataLabel(Size),DataValue(Size))
            DataLabel = ''
            DataValue = ''
            DataLabel(1) = 'Number of load steps' 
            DataValue(1) = IntToStr(CalParams%NLoadSteps)
            DataLabel(2) = 'Number of elements' 
            DataValue(2) = IntToStr(Counters%NEl)
            DataLabel(3) = 'Number of nodes' 
            DataValue(3) = IntToStr(Counters%NodTot)
            DataLabel(4) = 'Number of degrees of freedom' 
            DataValue(4) = IntToStr(Counters%N)
            DataLabel(5) = 'Number of entities' 
            DataValue(5) = IntToStr(Counters%NEntity)
            DataLabel(6) = 'Number of material points' 
            DataValue(6) = IntToStr(Counters%NParticles)
            DataLabel(7) = 'Number of solid material points' 
            DataValue(7) = IntToStr(Counters%SolidMaterialPoints)
            DataLabel(8) = 'Number of liquid material points' 
            DataValue(8) = IntToStr(Counters%LiquidMaterialPoints)
            DataLabel(9) = 'Number of phases' 
            DataValue(9) = IntToStr(CalParams%NumberOfPhases)
            write(BMRunit, *)'---general project data--------------------------------------------------------'
            do I = 1,Size
              write(BMRunit,'(A40,A50)')DataLabel(i),DataValue(i)
            end do
          end if
          
          ! FEM output data
          if (CalParams%OutputNumberNodes>0) then 
            write(BMRunit, *)'---nodal data---------------------------------------------------------'              
            write(BMRunit, *)'LoadStep,TimeStep,Node_ID,U_X,U_Y,U_Z'
            do I = 1, CalParams%OutputNumberNodes
              write(BMRunit,'(I5, I10, I10, 3G20.10E3)')  &
                CalParams%IStep,  &
                CalParams%TimeStep,  &
                CalParams%OutputNodes(I), &
                (TotalDisplacementSoil(ReducedDof(CalParams%OutputNodes(I)) + J), J=1,NVECTOR)
            end do
          end if
 
          ! MPM output data
          if (CalParams%OutputNumberParticles>0) then 
            write(BMRunit, *)'---material point data---------------------------------------------------------'              
            write(BMRunit, *)'LoadStep,TimeStep,MP_ID,X,Y,Z'
            do I = 1, CalParams%OutputNumberParticles
            
                if (CalParams%ApplyExcavation) then
                    do ParticleIndex= 1, Counters%NParticles
                        ParticleID = IDArray(ParticleIndex)
                        if (ParticleID==CalParams%OutputParticles(I)) then
                            exit
                        end if
                    end do
                else
                    ParticleIndex = CalParams%OutputParticles(I)
                end if 
           
                if (ParticleIndex<=Counters%NParticles) then
                    write(BMRunit,'(I5, I10, I10, 3G20.10)') &
                        CalParams%IStep, &
                        CalParams%TimeStep,  &
                        IDArray(ParticleIndex), &
                        (GlobPosArray(ParticleIndex,J), J=1,NVECTOR)
                end if
            end do
          end if


          if (CalParams%IStep==CalParams%NLoadSteps) then  ! write footer for last load step
            write(BMRunit, *)'---end-of-file-----------------------------------------------------------------'
          end if

        end subroutine WriteBenchmarkResults
        
        

        subroutine WriteBenchmarkResultsStressStrain()
        !**********************************************************************
        !
        !  Function : Output of calculation results for benchmarking stress and strains
        !             for each load step to a text file project_XXX.BMS
        !
        !**********************************************************************
        
        implicit none
        
          character(len=40), dimension(:), allocatable :: DataLabel
          character(len=50), dimension(:), allocatable :: DataValue
          integer(INTEGER_TYPE) :: I, J, Size, ParticleIndex, ParticleID
        
          if (CalParams%IStep==1) then  ! write file header for first load step
            Size=9
            allocate(DataLabel(Size),DataValue(Size))
            DataLabel = ''
            DataValue = ''
            DataLabel(1) = 'Number of load steps' 
            DataValue(1) = IntToStr(CalParams%NLoadSteps)
            DataLabel(2) = 'Number of elements' 
            DataValue(2) = IntToStr(Counters%NEl)
            DataLabel(3) = 'Number of nodes' 
            DataValue(3) = IntToStr(Counters%NodTot)
            DataLabel(4) = 'Number of degrees of freedom' 
            DataValue(4) = IntToStr(Counters%N)
            DataLabel(5) = 'Number of entities' 
            DataValue(5) = IntToStr(Counters%NEntity)
            DataLabel(6) = 'Number of material points' 
            DataValue(6) = IntToStr(Counters%NParticles)
            DataLabel(7) = 'Number of solid material points' 
            DataValue(7) = IntToStr(Counters%SolidMaterialPoints)
            DataLabel(8) = 'Number of liquid material points' 
            DataValue(8) = IntToStr(Counters%LiquidMaterialPoints)
            DataLabel(9) = 'Number of phases' 
            DataValue(9) = IntToStr(CalParams%NumberOfPhases)
            write(BMSunit, *)'---general project data--------------------------------------------------------'
            do I = 1,Size
              write(BMSunit,'(A40,A50)')DataLabel(i),DataValue(i)
            end do
          end if
          
          if (CalParams%OutputNumberParticles>0) then 
            write(BMSunit, *)'---material point data---------------------------------------------------------' 
            select case(NDIM)
                case(3)
                    write(BMSunit, '(28A27)')  &
                            'LoadStep ', 'TimeStep ', 'Time', 'PA(1) ','PA(2) ', 'PGravity ', 'ID_MP ', 'X ', 'Y ', 'Z ', &
                            'Ux ', 'Uy ', 'Uz ', 'SigmaXX ', 'SigmaYY ', 'SigmaZZ ', 'SigmaXY ', 'SigmaYZ ', 'SigmaZX ', &
                            'WPressure', 'EpsilonXX ', 'EpsilonYY ', 'EpsilonZZ ', 'GammaXY ', 'GammaYZ ', 'GammaZX ', &
                            'IntWeight ', 'MatID '
                    
                    do I = 1, CalParams%OutputNumberParticles
                        if (CalParams%ApplyExcavation) then
                            do ParticleIndex= 1, Counters%NParticles
                                ParticleID = IDArray(ParticleIndex)
                                if (ParticleID==CalParams%OutputParticles(I)) then
                                    exit
                                end if
                            end do
                        else
                            ParticleIndex = CalParams%OutputParticles(I)
                        end if 
                        if (ParticleIndex<=Counters%NParticles) then
                            write(BMSunit,'(I12, I12, 4G12.4, I12, 20G12.4, I12)') &
                              CalParams%IStep, &
                              CalParams%TimeStep, &
                              CalParams%OverallRealTime, &
                              CalParams%Multipliers%SolidACurrent,  &
                              CalParams%Multipliers%GravityCurrent,  &
                              IDArray(ParticleIndex), &
                              (GlobPosArray(ParticleIndex,J), J=1,NVECTOR), & 
                              (UArray(ParticleIndex,J), J=1,NVECTOR), & 
                              (SigmaEffArray(ParticleIndex,J), J=1,NTENSOR), & 
                              Particles(ParticleIndex)%WaterPressure, &
                              (GetEpsI(Particles(ParticleIndex),J), J=1,NTENSOR), & 
                              Particles(ParticleIndex)%IntegrationWeight, &
                              MaterialIDArray(ParticleIndex)
                        end if
                    end do
                    
                case(2)
                     write(BMSunit,  *)  &
                            'LoadStep ', 'TimeStep ', 'Time', 'PA(1) ','PA(2) ', 'PGravity ', 'ID_MP ', 'X ', 'Y ', &
                            'Ux ', 'Uy ', 'SigmaXX ', 'SigmaYY ', 'SigmaZZ ','SigmaXY ', &
                            'WPressure', 'EpsilonXX ', 'EpsilonYY ', 'EpsilonZZ ','GammaXY ', &
                            'IntWeight ', 'MatID '
                     
                    do I = 1, CalParams%OutputNumberParticles
                        if (CalParams%ApplyExcavation) then
                            do ParticleIndex= 1, Counters%NParticles
                                ParticleID = IDArray(ParticleIndex)
                                if (ParticleID==CalParams%OutputParticles(I)) then
                                    exit
                                end if
                            end do
                        else
                            ParticleIndex = CalParams%OutputParticles(I)
                        end if 
                        if (ParticleIndex<=Counters%NParticles) then
                            write(BMSunit,'(I12, I12, 4G12.4, I12, 14G12.4, I12)') & 
                              CalParams%IStep, &
                              CalParams%TimeStep, &
                              CalParams%OverallRealTime, &
                              CalParams%Multipliers%SolidACurrent,  &
                              CalParams%Multipliers%GravityCurrent,  &
                              IDArray(ParticleIndex), &
                              (GlobPosArray(ParticleIndex,J), J=1,NVECTOR), & 
                              (UArray(ParticleIndex,J), J=1,NVECTOR), & 
                              (SigmaEffArray(ParticleIndex,J), J=1,NTENSOR), & 
                              Particles(ParticleIndex)%WaterPressure, &
                              (GetEpsI(Particles(ParticleIndex),J), J=1,NTENSOR), & 
                              Particles(ParticleIndex)%IntegrationWeight, &
                              MaterialIDArray(ParticleIndex)
                        end if
                    end do
            end select
            
          end if
        
        
          if (CalParams%IStep==CalParams%NLoadSteps) then  ! write footer for last load step
            write(BMSunit, *)'---end-of-file-----------------------------------------------------------------'
          end if
        
        end subroutine WriteBenchmarkResultsStressStrain
        

        function IntToStr(Name)
        !**********************************************************************
        !
        !  Function : Convert integer(INTEGER_TYPE) into string
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: Name
          character(len=50) :: IntToStr, Str
          
          write(Str,*)Name
          IntToStr = Str

        end function IntToStr


        subroutine CloseTextOutputFiles()
        !**********************************************************************
        !
        !    Function:  Close the text files used for output of test data
        !               OUT, TST, MLG, PAR, CTSSum, RX, TVM, ENG, MOM, TIM, TRC, INF
        !
        !**********************************************************************
        
        implicit none
          ! Local variables
          integer(INTEGER_TYPE) :: I

          call CloseFile(OUTunit)
          call CloseFile(TSTunit)
          call CloseFile(BMRunit)
          !call CloseFile(BMSunit)
          call CloseFile(LOGunit)

          if ( (CalParams%OutputNumberParticles>0) .and. (CalParams%OutputNumberParticles<=MAXOUTPUTPARTICLES) ) then
            do I = 1, CalParams%OutputNumberParticles
              if (CalParams%OutputParticles(I)>0) then
                call CloseFile(PARUnit + I)
              end if
            end do
          end if
          do I = 1, Counters%NReactionSurfaceOutput
            call CloseFile(SURFReacUnit + I)
          end do

          call CloseFile(RXunit)
          call CloseFile(TVMunit)
          call CloseFile(ENGUnit)
          call CloseFile(INFunit)

        end subroutine CloseTextOutputFiles

        subroutine CloseFile(IONumber)
        !**********************************************************************
        !
        !  Function : Closes a file with unit number of IONumber
        !
        !**********************************************************************
        
        implicit none
        integer(INTEGER_TYPE), intent(in):: IONumber
        logical :: isOpened

        inquire(UNIT = IONumber, opened = isOpened)
        if (isOpened) then
          endfile(IONumber)
          close(IONumber)
        endif

        end subroutine CloseFile


      end module ModWriteTestData
