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
    !	Copyright (C) 2025  Members of the Anura3D MPM Research Community 
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
	  
	  
	  module ModReadOutputParameters
      !**********************************************************************
      !
      ! Function: Contains routines for reading output data from OPD file
      !
      !  Note: Version archived subroutines are contained at the end of the module
      !
      !     $Revision: 10167 $
      !     $Date: 2024-06-25 15:15:53 +0200 (Tue, 25 Jun 2024) $
      !
      !**********************************************************************
      use ModGlobalConstants
      use ModCounters
      use ModString
      use ModGlobalData
      use ModFeedback
      use ModFileIO
      use ModReadCalculationData
      use ModMPMData
      use ModAdjustParticleDiscretisation

      implicit none
        
        type OutputFrequency ! information required to write output every some time steps
        
          integer(INTEGER_TYPE):: CalculationSteps = 1
        
        end type OutputFrequency
        
        type OutputVariables
            
          logical :: MassMixture = .true.,&
            WeightMixture = .true.,&
            Porosity = .true.,&
            MeanEffectiveStress = .true.,&
            DeviatoricStress = .true.,&
            DeviatoricStrain = .true.,&
            DegreeofSaturation = .true.,&
            StateVariables = .true.,&
            Damping = .true.,&
            MassSolid = .true.,&
            WeightSolid = .true.,&
            VolumetricStrainSolid = .true.,&
            MassLiquid = .true.,&
            WeightLiquid = .true.,&
            VolumetricStrainLiquid = .true.,&
            PressureLiquid = .true.,&
            MassGas = .true.,&
            WeightGas = .true.,&
            PressureGas = .true.,&
            VolumetricStrainGas = .true.,&
            TotalStress = .true.,&
            EffectiveStress = .true.,&
            Strain = .true.,&
            GlobalPosition = .true.,&
            LocalPosition = .true.,&
            BodyForceMixture = .true.,&
            BodyForceSolid = .true.,&
            ExternalForceSolid = .true.,&
            AccelerationSolid = .true.,&
            VelocitySolid = .true.,&
            DisplacementSolid = .true.,&
            BodyForceLiquid = .true.,&
            ExternalForceLiquid = .true.,&
            VelocityLiquid = .true.,&
            DisplacementLiquid = .true.,&  
            BodyForceGas = .true.,&
            ExternalForceGas = .true.,&
            VelocityGas = .true.,&
            DisplacementGas = .true.
        
        end type OutputVariables
        
        type(OutputVariables), allocatable, public, save :: Variables ! Stores variables that manipulate the Dynamic MPM calculation process
        
        
        type OutputMP
            
          integer(INTEGER_TYPE) :: NumberElementsForMPOutput, NumberOfOutputMPPerElement
          integer(INTEGER_TYPE),allocatable :: ElementsForMPOutput(:)
          real (REAL_TYPE) :: TimeBetweenMPOutput

        end type OutputMP
        
        type OutputParams(vsize, dsize)
            
          integer(INTEGER_TYPE), LEN :: vsize, dsize
          type(OutputVariables) :: Variables
          type(OutputMP) :: MP
          type(OutputFrequency) :: Frequency
          
        end type OutputParams
        
        type(OutputParams(:, :)), allocatable, public, save :: OutParams !
        
        
     
    contains ! Routines of this module

        
       subroutine InitialiseOutputParameters()
       !**********************************************************************
       !
       ! Function: Initialises output parameters stored in OpdParams%
       !
       !**********************************************************************
       
       implicit none
       
       ! allocate main OutParams structure
       allocate(OutputParams(NVECTOR, NDIM)::OutParams)
       
       ! Set Default Parameters
       
        ! OutputFrequencyType
       
       OutParams%Frequency%CalculationSteps = 1
        
        ! OutputVariablesType

        OutParams%Variables%MassMixture = .true.
        OutParams%Variables%WeightMixture = .true.
        OutParams%Variables%Porosity = .true.
        OutParams%Variables%MeanEffectiveStress = .true.
        OutParams%Variables%DeviatoricStress = .true.
        OutParams%Variables%DeviatoricStrain = .true.
        OutParams%Variables%DegreeofSaturation = .true.
        OutParams%Variables%StateVariables = .true.
        OutParams%Variables%Damping = .true.
        OutParams%Variables%MassSolid = .true.
        OutParams%Variables%WeightSolid = .true.
        OutParams%Variables%VolumetricStrainSolid = .true.
        OutParams%Variables%MassLiquid = .true.
        OutParams%Variables%WeightLiquid = .true.
        OutParams%Variables%VolumetricStrainLiquid = .true.
        OutParams%Variables%PressureLiquid = .true.
        OutParams%Variables%MassGas = .true.
        OutParams%Variables%WeightGas = .true.
        OutParams%Variables%PressureGas = .true.
        OutParams%Variables%VolumetricStrainGas = .true.
        OutParams%Variables%TotalStress = .true.
        OutParams%Variables%EffectiveStress = .true.
        OutParams%Variables%Strain = .true.
        OutParams%Variables%GlobalPosition = .true.
        OutParams%Variables%LocalPosition = .true.
        OutParams%Variables%BodyForceMixture = .true.
        OutParams%Variables%BodyForceSolid = .true.
        OutParams%Variables%ExternalForceSolid = .true.
        OutParams%Variables%AccelerationSolid = .true.
        OutParams%Variables%VelocitySolid = .true.
        OutParams%Variables%DisplacementSolid = .true.
        OutParams%Variables%BodyForceLiquid = .true.
        OutParams%Variables%ExternalForceLiquid = .true.
        OutParams%Variables%VelocityLiquid = .true.
        OutParams%Variables%DisplacementLiquid = .true.
        OutParams%Variables%BodyForceGas = .true.
        OutParams%Variables%ExternalForceGas = .true.
        OutParams%Variables%VelocityGas = .true.
        OutParams%Variables%DisplacementGas = .true.
        
        ! OutputMP
        
        OutParams%MP%NumberElementsForMPOutput = 0
        OutParams%MP%ElementsForMPOutput = 0
        OutParams%MP%NumberOfOutputMPPerElement = 0
        OutParams%MP%TimeBetweenMPOutput = 0.0
   
       end subroutine InitialiseOutputParameters
              
        
        subroutine ReadOutputParameters()
        !**********************************************************************
        !
        ! Function : Reads output variables from the OPD file         
        !
        !**********************************************************************
        
        implicit none

          character(len=255) :: FileName, FileVersion
          integer(INTEGER_TYPE) :: FileUnit

          ! local variables
          integer(INTEGER_TYPE) :: I, ParticleIndex, NPart, NPart2, IEl, NPartcount
          integer(INTEGER_TYPE) :: ios ! used for error control
          real(REAL_TYPE) :: DumI(1)
          character(len=255) :: DumS
          real(REAL_TYPE) :: DumR(4)
          character(len=255) :: DumFileName
          character(len=21) :: messageIOS = 'OPD file: Can''t read '
          character(len=255) :: TName, BName
          character(len=3) :: Stage
          
          ! OPD file only works with Anura3D_v2025

          FileVersion = CalParams%CPSversion
          FileUnit = TMP_UNIT
          
          if (FileVersion == Anura3D_v2025) then
            write(Stage,'(I0)') CalParams%StageNumber
            FileName = trim(CalParams%FileNames%ProjectName)//OPD_STAGE_FILE_EXTENSION//Stage
          else              
            RETURN
          end if
          
          ! open OPD file
          call FileOpen(FileUnit, trim(FileName))
          
          call GiveMessage('Reading OPD file: ' // trim(FileName) )
          
          do
              
            read(FileUnit, '(A)', iostat=ios) TName
            call Assert( ios == 0, 'OPD file: Can''t read flag from OPD file.' )
            BName = Upcase(TName)
            call WriteInLogFile('Reading flag from OPD file: ' // trim(BName), FEEDBACK_LEVEL_ALWAYS)
            
        ! OUTPUT FREQUENCY
            
            if (trim(BName)=='$$CALCULATION_STEPS_BETWEEN_OUTPUTS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              OutParams%Frequency%CalculationSteps = DumI(1)
            
        ! OUTPUT VARIABLES
        ! Scalar    
            else if (trim(BName)=='$$MASS_MIXTURE') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%MassMixture = .false.
              
            else if (trim(BName)=='$$WEIGHT_MIXTURE') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%WeightMixture = .false.  
              
            else if (trim(BName)=='$$POROSITY') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%Porosity = .false.
              
            else if (trim(BName)=='$$MEAN_EFFECTIVE_STRESS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%MeanEffectiveStress = .false.
              
            else if (trim(BName)=='$$DEVIATORIC_STRESS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%DeviatoricStress = .false.
              
            else if (trim(BName)=='$$DEVIATORIC_STRAIN') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%DeviatoricStrain = .false.
              
            else if (trim(BName)=='$$DEGREE_OF_SATURATION') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%DegreeofSaturation = .false.
              
            else if (trim(BName)=='$$STATE_VARIABLES') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%StateVariables = .false.
              
            else if (trim(BName)=='$$DAMPING') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%Damping = .false.
              
            else if (trim(BName)=='$$MASS_SOLID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%MassSolid = .false.
              
            else if (trim(BName)=='$$WEIGHT_SOLID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%WeightSolid = .false.
              
            else if (trim(BName)=='$$VOLUMETRIC_STRAIN_SOLID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%VolumetricStrainSolid = .false.
              
            else if (trim(BName)=='$$MASS_LIQUID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%MassLiquid = .false.
              
            else if (trim(BName)=='$$WEIGHT_LIQUID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%WeightLiquid = .false.
              
            else if (trim(BName)=='$$VOLUMETRIC_STRAIN_LIQUID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%VolumetricStrainLiquid = .false.
              
            else if (trim(BName)=='$$PRESSURE_LIQUID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%PressureLiquid = .false.
              
            else if (trim(BName)=='$$MASS_GAS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%MassGas = .false.
              
            else if (trim(BName)=='$$WEIGHT_GAS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%WeightGas = .false.
              
            else if (trim(BName)=='$$PRESSURE_GAS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%PressureGas = .false.
              
            else if (trim(BName)=='$$VOLUMETRIC_STRAIN_GAS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%VolumetricStrainGas = .false.        
            
            ! Vector   
            else if (trim(BName)=='$$GLOBAL_POSITION') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%GlobalPosition = .false.
              
            else if (trim(BName)=='$$LOCAL_POSITION') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%LocalPosition = .false.
              
            else if (trim(BName)=='$$BODY_FORCE_MIXTURE') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%BodyForceMixture = .false.
              
            else if (trim(BName)=='$$BODY_FORCE_SOLID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%BodyForceSolid = .false.
              
            else if (trim(BName)=='$$EXTERNAL_FORCE_SOLID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%ExternalForceSolid = .false.
              
            else if (trim(BName)=='$$ACCELERATION_SOLID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%AccelerationSolid = .false.
              
            else if (trim(BName)=='$$VELOCITY_SOLID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%VelocitySolid = .false.
              
            else if (trim(BName)=='$$DISPLACEMENT_SOLID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%DisplacementSolid = .false.
              
            else if (trim(BName)=='$$BODY_FORCE_LIQUID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%BodyForceLiquid = .false.
              
            else if (trim(BName)=='$$EXTERNAL_FORCE_LIQUID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%ExternalForceLiquid = .false.
              
            else if (trim(BName)=='$$VELOCITY_LIQUID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%VelocityLiquid = .false.
              
            else if (trim(BName)=='$$DISPLACEMENT_LIQUID') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%DisplacementLiquid = .false.
              
            else if (trim(BName)=='$$BODY_FORCE_GAS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%BodyForceGas = .false.
              
            else if (trim(BName)=='$$EXTERNAL_FORCE_GAS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%ExternalForceGas = .false.
              
            else if (trim(BName)=='$$VELOCITY_GAS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%VelocityGas = .false.
              
            else if (trim(BName)=='$$DISPLACEMENT_GAS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%DisplacementGas = .false.
              
            ! Tensor
            else if (trim(BName)=='$$TOTAL_STRESS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%TotalStress = .false.              
              
            else if (trim(BName)=='$$EFFECTIVE_STRESS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%EffectiveStress = .false.
              
            else if (trim(BName)=='$$STRAIN') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'OPD file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) OutParams%Variables%Strain = .false.

        ! OUTPUT MP              
            else if (trim(BName)=='$$NUMBER_OF_ELEMENTS_FOR_MP_OUTPUT') then 
              read(FileUnit, *, iostat=ios) OutParams%MP%NumberElementsForMPOutput
              call Assert( ios == 0, messageIOS//trim(BName) )
              
            else if (trim(BName)=='$$NUMBER_OF_OUTPUT_MP_PER_ELEMENT') then 
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              OutParams%MP%NumberOfOutputMPPerElement = DumI(1)
                
            else if (trim(BName)=='$$ELEMENTS_FOR_MP_OUTPUT') then 
              if ((OutParams%MP%NumberElementsForMPOutput>0).and.(.not.IsFollowUpPhase())) then
                allocate(OutParams%MP%ElementsForMPOutput(OutParams%MP%NumberElementsForMPOutput)) 
                do I = 1, OutParams%MP%NumberElementsForMPOutput
                  read(FileUnit, *, iostat=ios) OutParams%MP%ElementsForMPOutput(I)
                  call Assert( ios == 0, messageIOS//trim(BName) )
                end do
                  NPartcount = 1
                  do I = 1, OutParams%MP%NumberElementsForMPOutput
                      IEl=OutParams%MP%ElementsForMPOutput(I)
                      do NPart2 = 1, OutParams%MP%NumberOfOutputMPPerElement ! loop over particles
                          ParticleIndex = GetParticleIndex(NPart2, OutParams%MP%ElementsForMPOutput(I))   
                          CalParams%OutputParticles(NPartCount) = ParticleIndex      ! store particles index in the CalParamsArray for versions previous to v2025
                          NPartCount = NPartcount+1
                      end do
                  end do
                  CalParams%OutputNumberParticles = NPartCount-1
              else
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) ) 
              end if
              
            else if (trim(BName)=='$$TIME_BETWEEN_MP_OUTPUTS') then 
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              OutParams%MP%TimeBetweenMPOutput = DumI(1)
              
        ! END OF OPD-FILE
              
            else if (trim(BName)=='$$END') then
              EXIT
            else
              if (ContainLetters(trim(BName))) then
                call GiveWarning('Unknown flag in OPD file: ' // trim(BName))
              endif
            end if

          end do

        
        end subroutine ReadOutputParameters
        
        end module ModReadOutputParameters




