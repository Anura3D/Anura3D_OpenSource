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
	  !	  Copyright (C) 2023  Members of the Anura3D MPM Research Community 
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
	  
	  
	  module ModAdjustParticleDiscretisation
      !**********************************************************************
      !
      !    Function:  To add and remove particles from the discretisation and updating 
      !               particle housekeeping
      !
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 9808 $
      !     $Date: 2022-10-13 16:48:47 +0200 (do, 13 okt 2022) $
      !
      !**********************************************************************

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

        
        subroutine CorrectParticleMassAndBodyForceAxisymmetric(ParticleIndex)
        !**********************************************************************
        !
        !   Function:   Correct material points mass and body forces for axisymmetric analysis.
        !
        !**********************************************************************
        implicit none
        
          integer, intent(in) :: ParticleIndex
        
          ! local variables
          real(REAL_TYPE) :: r
        
          if ( .not. ISAXISYMMETRIC ) RETURN
        
          r = GlobPosArray(ParticleIndex, 1) ! index 1: radial direction
          MassArray(ParticleIndex)              = MassArray(ParticleIndex)            * r
          MassWaterArray(ParticleIndex)         = MassWaterArray(ParticleIndex)       * r
          Particles(ParticleIndex)%MassMixed    = Particles(ParticleIndex)%MassMixed  * r
          Particles(ParticleIndex)%FBody        = Particles(ParticleIndex)%FBody      * r
          Particles(ParticleIndex)%FBodyWater   = Particles(ParticleIndex)%FBodyWater * r
          Particles(ParticleIndex)%FBodyGas     = Particles(ParticleIndex)%FBodyGas   * r
          Particles(ParticleIndex)%FBodyMixed   = Particles(ParticleIndex)%FBodyMixed * r

        end subroutine CorrectParticleMassAndBodyForceAxisymmetric

        
        subroutine CorrectParticleIntegrationWeightAxisymmetric(ParticleIndex)
        !**********************************************************************
        !
        !   Function:   Correct material points integration weight for axisymmetric analysis.
        !
        !**********************************************************************
        implicit none
        
        integer, intent(in) :: ParticleIndex
        
          ! local variables
          real(REAL_TYPE) :: r
        
          if ( .not. ISAXISYMMETRIC ) RETURN
        
          r = GlobPosArray(ParticleIndex,1) ! index 1: radial direction
          Particles(ParticleIndex)%IntegrationWeight = Particles(ParticleIndex)%IntegrationWeight * r
        
        end subroutine CorrectParticleIntegrationWeightAxisymmetric
        
        
        integer(INTEGER_TYPE) function CheckRemovalParticles()
        !**********************************************************************
        !
        !    Function:  Checks which particles can be removed from the
        !               discretisation.
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, IParticle, &
                     NVirtualParticles, &
                     ParticleIndex
          logical :: ContainsBoundaryParticle
        
          CheckRemovalParticles = 0
        
          do IElement = 1, Counters%NEl
            
            NVirtualParticles = NVirtualParticlesEle(IElement)
            ContainsBoundaryParticle = HasBoundaryParticle(IElement)
            
            if ((NVirtualParticles > 0) .and. &
                ((.not.IsActiveElement(IElement)).or.ContainsBoundaryParticle)) then
              CheckRemovalParticles = CheckRemovalParticles + NVirtualParticles
     
              ! Before removal of virtual particles their stress data should be mapped
              ! to the remaining particles in order to take into consideration their information
              ! about the stress field - averaging of particle data (done also for particles to be removed)

              call AverageElementParticleStateParam(IElement)

              ! Mark virtual particles to be removed
              do IParticle = 1, NPartEle(IElement)
                ParticleIndex = GetParticleIndexFunction(IParticle, IElement)
                if (IsVirtualParticle(ParticleIndex) ) then
                  Particles(ParticleIndex)%Kind = REMOVEDPARTICLE
                end if
              end do
            end if
          end do
        
        end function CheckRemovalParticles
        
        subroutine AdjustParticlesArrays(NUpdateParticles)
        !**********************************************************************
        !
        !    Function:  Adjusts the array Particles and updates the counter
        !               Counters%NParticles to the new number of particles.
        !               Reset house keeping data structure which contains 
        !               incorrect and redundant information.
        !
        !     NUpdateParticles : Change of number of particles
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
     
        implicit none
     
          integer(INTEGER_TYPE), intent(in) :: NUpdateParticles
          ! Local variables
          type(ParticleType(:,:,:)), dimension(:), allocatable :: TempParticles
          integer(INTEGER_TYPE), allocatable::TempElementIDArray(:), TempIDArray(:), TempEntityIDArray(:), TempMaterialIDArray(:),TempMaterialPointTypeArray(:)
          real(REAL_TYPE), allocatable :: TempUArray(:,:), TempUStepArray(:,:), TempShapeValuesArray(:,:), TempUPhaseArray(:,:)
          real(REAL_TYPE), allocatable :: TempGlobPosArray(:,:), TempVelocityArray(:,:)
          real(REAL_TYPE), allocatable :: TempSigmaEffArray(:,:), TempSigmaEff0Array(:,:), TempDShapeValuesArray(:,:,:),TempVelocityWaterArray(:,:) 
          real(REAL_TYPE), allocatable :: TempVelocityGasArray(:,:), TempAccelerationArray(:,:), TempMassArray(:), TempESMstatevArray(:,:), TempMassWaterArray(:)
          !real(REAL_TYPE), allocatable :: TempESMpropsArray(:,:)
 
          integer(INTEGER_TYPE) :: IError, I, NewParticleIndex, IErrorTemp
		  integer(INTEGER_TYPE), dimension(:) , allocatable   :: MPtrackindex

          ! Allocate temporary copy of particle data
          allocate(ParticleType(NTENSOR, NVECTOR, MAX_LOAD_SYSTEMS) :: TempParticles(Counters%NParticles), stat = IErrorTemp)
          
          do i = 1, Counters%NParticles 
              call SetParticleStructureDefault(TempParticles(i))
          end do
          
          allocate(TempElementIDArray(Counters%NParticles), stat = IErrorTemp)
          allocate(TempIDArray(Counters%NParticles), stat = IErrorTemp)
          allocate(TempEntityIDArray(Counters%NParticles), stat = IErrorTemp)
          allocate(TempMaterialIDArray(Counters%NParticles), stat = IErrorTemp)
          allocate(TempUArray(Counters%NParticles,NVECTOR), stat = IErrorTemp)
          allocate(TempMaterialPointTypeArray(Counters%NParticles), stat = IErrorTemp)
          allocate(TempUStepArray(Counters%NParticles,NVECTOR), stat = IErrorTemp)
          allocate(TempShapeValuesArray(Counters%NParticles,ELEMENTNODES), stat = IErrorTemp)
          allocate(TempUPhaseArray(Counters%NParticles,NVECTOR), stat = IErrorTemp)
          allocate(TempGlobPosArray(Counters%NParticles,NVECTOR), stat = IErrorTemp)
          allocate(TempVelocityArray(Counters%NParticles,NVECTOR), stat = IErrorTemp)
          allocate(TempSigmaEffArray(Counters%NParticles,NTENSOR), stat = IErrorTemp)
          allocate(TempSigmaEff0Array(Counters%NParticles,NTENSOR), stat = IErrorTemp)
          allocate(TempDShapeValuesArray(Counters%NParticles,ELEMENTNODES, NVECTOR), stat = IErrorTemp)
          allocate(TempVelocityWaterArray(Counters%NParticles,NVECTOR), stat = IErrorTemp)
          allocate(TempVelocityGasArray(Counters%NParticles,NVECTOR), stat = IErrorTemp)
          allocate(TempAccelerationArray(Counters%NParticles,NVECTOR), stat = IErrorTemp)
          allocate(TempMassArray(Counters%NParticles), stat = IErrorTemp)
          allocate(TempESMstatevArray(Counters%NParticles, NSTATEVAR), stat = IErrorTemp)
          !allocate(TempESMpropsArray(Counters%NParticles, NPROPERTIES), stat = IErrorTemp)
          allocate(TempMassWaterArray(Counters%NParticles), stat = IErrorTemp)

          ! Copy particle data to temporary array
          
          ! TempParticles = Particles
          do i = 1, Counters%NParticles 
              TempParticles(i) = CopyParticle(Particles(i), i)
          end do
          
          TempElementIDArray = ElementIDArray
          TempIDArray = IDArray
          TempEntityIDArray = EntityIDArray
          TempMaterialIDArray = MaterialIDArray
          TempUArray = UArray
          TempMaterialPointTypeArray = MaterialPointTypeArray
          TempUStepArray = UStepArray
          TempShapeValuesArray = ShapeValuesArray
          TempUPhaseArray = UPhaseArray
          TempGlobPosArray = GlobPosArray
          TempVelocityArray = VelocityArray
          TempSigmaEffArray = SigmaEffArray
          TempSigmaEff0Array = SigmaEff0Array
          TempDShapeValuesArray = DShapeValuesArray
          TempVelocityWaterArray = VelocityWaterArray
          TempVelocityGasArray = VelocityGasArray
          TempAccelerationArray = AccelerationArray
          TempMassArray = MassArray
          TempESMstatevArray = ESMstatevArray
          !TempESMpropsArray = ESMpropsArray
          TempMassWaterArray = MassWaterArray

          ! Destroy old particle data
          deallocate(Particles, stat = IError)
          deallocate(ElementIDArray, stat = IError)
          deallocate(IDArray, stat = IError)
          deallocate(EntityIDArray, stat = IError)
          deallocate(MaterialIDArray, stat = IError)
          deallocate(UArray, stat = IError)
          deallocate(MaterialPointTypeArray, stat = IError)
          deallocate(UStepArray, stat = IError)
          deallocate(ShapeValuesArray, stat = IError)
          deallocate(UPhaseArray, stat = IError)
          deallocate(GlobPosArray, stat = IError)
          deallocate(VelocityArray, stat = IError)
          deallocate(SigmaEffArray, stat = IError)
          deallocate(SigmaEff0Array, stat = IError)
          deallocate(DShapeValuesArray, stat = IError)
          deallocate(VelocityWaterArray, stat = IError)
          deallocate(VelocityGasArray, stat = IError)
          deallocate(AccelerationArray, stat = IError)
          deallocate(MassArray, stat = IError)
          deallocate(ESMstatevArray, stat = IError)
          !deallocate(ESMpropsArray, stat = IError)
          deallocate(MassWaterArray, stat = IError)
 
          ! Update number of particles
          Counters%NParticles = Counters%NParticles + NUpdateParticles

          ! Create new particle data array
          allocate(ParticleType(NTENSOR, NVECTOR, MAX_LOAD_SYSTEMS) :: Particles(Counters%NParticles), stat = IError)
          do i = 1, Counters%NParticles 
            call SetParticleStructureDefault(Particles(i))
          end do
                    
          allocate(ElementIDArray(Counters%NParticles), stat = IError)
          allocate(IDArray(Counters%NParticles), stat = IError)
          allocate(EntityIDArray(Counters%NParticles), stat = IError)
          allocate(MaterialIDArray(Counters%NParticles), stat = IError)
          allocate(UArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(MaterialPointTypeArray(Counters%NParticles), stat = IError)
          allocate(UStepArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(ShapeValuesArray(Counters%NParticles,ELEMENTNODES), stat = IError)
          allocate(UPhaseArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(GlobPosArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(VelocityArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(SigmaEffArray(Counters%NParticles,NTENSOR), stat = IError)
          allocate(SigmaEff0Array(Counters%NParticles,NTENSOR), stat = IError)
          allocate(DShapeValuesArray(Counters%NParticles,ELEMENTNODES, NVECTOR), stat = IError)
          allocate(VelocityWaterArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(VelocityGasArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(AccelerationArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(MassArray(Counters%NParticles), stat = IError)
          allocate(ESMstatevArray(Counters%NParticles, NSTATEVAR), stat = IError)
!          allocate(ESMpropsArray(Counters%NParticles, NPROPERTIES), stat = IError)
          allocate(MassWaterArray(Counters%NParticles), stat = IError)
          
          ! Copy old particle array to new particle array, remove deleted particles
          NewParticleIndex = 1
          do I = 1, Counters%NParticles - NUpdateParticles ! Loop over original number of particles
                if ((ANY(CalParams%OutputParticles == I)) .AND. (NUpdateParticles<0)) then
                    MPtrackindex = 0
                    MPtrackindex = MINLOC(abs(CalParams%OutputParticles-I))
                    CalParams%OutputParticles(MPtrackindex) = NewParticleIndex
                end if    
                          
            if (TempParticles(I)%Kind/=REMOVEDPARTICLE) then
              Particles(NewParticleIndex) = CopyParticle(TempParticles(I), TempIDArray(I))
              ElementIDArray(NewParticleIndex) = TempElementIDArray(I) 
              IDArray(NewParticleIndex) = TempIDArray(I)
              EntityIDArray(NewParticleIndex) = TempEntityIDArray(I)
              MaterialIDArray(NewParticleIndex) = TempMaterialIDArray(I)
              UArray(NewParticleIndex,:) = TempUArray(I,:)
              MaterialPointTypeArray(NewParticleIndex) = TempMaterialPointTypeArray(I)
              UStepArray(NewParticleIndex,:) = TempUStepArray(I,:)
              ShapeValuesArray(NewParticleIndex,:) = TempShapeValuesArray(I,:) 
              UPhaseArray(NewParticleIndex,:) = TempUPhaseArray(I,:) 
              GlobPosArray(NewParticleIndex,:) = TempGlobPosArray(I,:)
              VelocityArray(NewParticleIndex,:) = TempVelocityArray(I,:)
              SigmaEffArray(NewParticleIndex,:) = TempSigmaEffArray(I,:)
              SigmaEff0Array(NewParticleIndex,:) = TempSigmaEff0Array(I,:)
              DShapeValuesArray(NewParticleIndex,:,:) = TempDShapeValuesArray(I,:,:)
              VelocityWaterArray(NewParticleIndex,:) = TempVelocityWaterArray(I,:)
              VelocityGasArray(NewParticleIndex,:) = TempVelocityGasArray(I,:)
              AccelerationArray(NewParticleIndex,:) = TempAccelerationArray(I,:)
              MassArray(NewParticleIndex) = TempMassArray(I)
              ESMstatevArray(NewParticleIndex,:) = TempESMstatevArray(I,:)
              !ESMpropsArray(NewParticleIndex,:) = TempESMpropsArray(I,:)
              MassWaterArray(NewParticleIndex) = TempMassWaterArray(I)
              NewParticleIndex = NewParticleIndex + 1
            end if
          end do
          
          ! Destroy temporary copies
          if (allocated(TempParticles)) then
            deallocate(TempParticles, stat = IError)
          end if

          if (allocated(TempElementIDArray)) then
             deallocate(TempElementIDArray, stat = IError)
          end if

          if (allocated(TempIDArray)) then
             deallocate(TempIDArray, stat = IError)
          end if

          if (allocated(TempEntityIDArray)) then
             deallocate(TempEntityIDArray, stat = IError)
          end if

          if (allocated(TempMaterialIDArray)) then
             deallocate(TempMaterialIDArray, stat = IError)
          end if

          if (allocated(TempUArray)) then
             deallocate(TempUArray, stat = IError)
          end if

          if (allocated(TempMaterialPointTypeArray)) then
             deallocate(TempMaterialPointTypeArray, stat = IError)
          end if

          if (allocated(TempUStepArray)) then
             deallocate(TempUStepArray, stat = IError)
          end if

          if (allocated(TempShapeValuesArray)) then
             deallocate(TempShapeValuesArray, stat = IError)
          end if

          if (allocated(TempUPhaseArray)) then
             deallocate(TempUPhaseArray, stat = IError)
          end if

          if (allocated(TempGlobPosArray)) then
             deallocate(TempGlobPosArray, stat = IError)
          end if

          if (allocated(TempVelocityArray)) then
             deallocate(TempVelocityArray, stat = IError)
          end if

          if (allocated(TempSigmaEffArray)) then
             deallocate(TempSigmaEffArray, stat = IError)
          end if

          if (allocated(TempSigmaEff0Array)) then
             deallocate(TempSigmaEff0Array, stat = IError)
          end if

          if (allocated(TempDShapeValuesArray)) then
             deallocate(TempDShapeValuesArray, stat = IError)
          end if

          if (allocated(TempVelocityWaterArray)) then
             deallocate(TempVelocityWaterArray, stat = IError)
          end if

          if (allocated(TempVelocityGasArray)) then
             deallocate(TempVelocityGasArray, stat = IError)
          end if

          if (allocated(TempAccelerationArray)) then
             deallocate(TempAccelerationArray, stat = IError)
          end if

          if (allocated(TempMassArray)) then
             deallocate(TempMassArray, stat = IError)
          end if

          if (allocated(TempESMstatevArray)) then
             deallocate(TempESMstatevArray, stat = IError)
          end if
          
          !if (allocated(TempESMpropsArray)) then
          !   deallocate(TempESMpropsArray, stat = IError)
          !end if

          if (allocated(TempMassWaterArray)) then
             deallocate(TempMassWaterArray, stat = IError)
          end if
          !-----

          ! Change size of array EleParticles and initialise it to -1
          deallocate(EleParticles, stat = IError)
          allocate(EleParticles(Counters%NParticles), stat = IError)
          EleParticles = -1
          
          ! Reset content of array EleParticlesHelp
          EleParticlesHelp = -1
          
          ! Reset content of array NPartEle
          NPartEle = 0
          
          ! Adjust PartFactor to new number of particles
          call SetPartFactor(DecimalFactor(Counters%NParticles))
     
        end subroutine AdjustParticlesArrays
        
        subroutine FillParticlesArray(NUpdateParticles, NRemovedParticles)
        !**********************************************************************
        !
        !    Function:  Places in each activated element void of particles new particles.
        !               Initialises the particle data.
        !
        !     NUpdateParticles : Change of number of particles
        !     NRemovedParticles : Number of removed particles
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NUpdateParticles, NRemovedParticles
          ! Local variables
          integer(INTEGER_TYPE) :: ElementID, ParticleIndex, ICounter, &
                     IParticle
        
          ICounter = 1
          do ElementID = 1, Counters%NEl
            if (IsEmptyElement(ElementID)) then ! Loop over all elements belonging to holes

              do IParticle = 1, CalParams%NVirtualParticles
                ParticleIndex = Counters%NParticles -  &
                                NUpdateParticles - &
                                NRemovedParticles +  &
                                ICounter

                ICounter = ICounter + 1

                call CreateParticleData(ElementID, &
                                        ParticleIndex, &
                                        NEWVIRTUALPARTICLE)
                if (CalParams%OutputDebugData) then
                  call WriteInLogFile('Place virtual particle ' // &
                                      trim(String(IDArray(ParticleIndex))) // &
                                      ' in element '// trim(String(ElementID)))
                end if
              end do
            end if
          end do

          if ((ICounter - 1)/= &
              NUpdateParticles + NRemovedParticles) then
            call GiveError('Not all particles correctly initialised.')
          end if

        end subroutine FillParticlesArray


        subroutine CreateParticleData(ElementID, &
                                      ParticleIndex, &
                                      Kind)
        !**********************************************************************
        !
        !  Function:  Generates the basic data of the added material point inside
        !             the Particle array.
        !
        !     ElementID : ID of the considered element
        !     ParticleIndex : Global index of the particle
        !     Kind : Material or virtual particle
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ElementID, ParticleIndex, Kind
          
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR) :: LocPos
          real(REAL_TYPE) :: DampingFactor=0.0
          logical :: ISSubmerged
          integer(INTEGER_TYPE) :: MaterialID, IDMat, I
          !integer(INTEGER_TYPE) :: IDCounter=0
          real(REAL_TYPE) :: VirtualParticleMass
          integer(INTEGER_TYPE) :: MaterialPointType
        
          ! TODO Consider CalParams%NVirtualParticles particles
          do I = 1, NVECTOR
              LocPos(I) = 0.25
          end do
         
          if (GeoParams%ApplyLocalDampingElement) then
            DampingFactor = GeoParams%LocalDampingFactorElement(ElementID)
          end if    
          
          IsSubmerged = CalParams%ApplySubmergedCalculation
                      
          VirtualParticleMass = 0.0 ! TODO: Take small fraction of mass from surrounding particles
          MaterialID = 1 !By default select first material 
          do IDMat = 1, Counters%NLayers
            if (MaterialElements(IDMat, ElementID)==1) then
              MaterialID = IDMat
              EXIT
            end if
          end do
          
          ! determine type of material point
          if (NFORMULATION==1) then ! 1-layer formulation, material point is of type MIXTURE
            MaterialPointType = MaterialPointTypeMixture
          else ! 2-layer formulation, material point is of either type SOLID or LIQUID
            ! empty element algorithm not yet implemented for 2-layer formulation
            call GiveError('[subroutine CreateParticleData] Empty element algorithm not available for 2-layer formulation.')
          end if
          
          Particles(ParticleIndex) = &
            InitParticle(MaterialID, &
                         MaterialPointType, &
                         0.d0, &
                         LocPos, &
                         Kind, &
                         VirtualParticleMass, &
                         CalParams%GravityData%GAccel, &
                         CalParams%GravityData%GravityVector, &
                         DampingFactor, &
                         IsSubmerged, &
                         MassArray(ParticleIndex),&
                         MassWaterArray(ParticleIndex))
          
         ElementIDArray(ParticleIndex) = ElementID
         IDArray(ParticleIndex) = GetIDCounter()
         EntityIDArray(ParticleIndex) = SOFT_ENTITY
         MaterialIDArray(ParticleIndex) = MaterialID
         MaterialPointTypeArray(ParticleIndex) = MaterialPointType
          
        
          ! Update ShapeValues, DShapeValues
          call SetParticleShapeFunctionData(Particles(ParticleIndex), ParticleIndex)

          ! Update GlobPos
          call DetermineGlobalFromLocalCoord(ParticleIndex, ElementID)
          
          ! Update Mass and integration weight
          call CorrectParticleMassAndBodyForceAxisymmetric(ParticleIndex)
          call CorrectParticleIntegrationWeightAxisymmetric(ParticleIndex)
        
        end subroutine CreateParticleData
        
        subroutine DetermineGlobalFromLocalCoord(ParticleIndex,  &
                                                 IElement)
        !**********************************************************************
        !
        !    Function:  Determines the global coordinates associated with the local
        !               coordinates of particle ParticleIndex inside IElement.
        !
        !     ParticleIndex : ID of the considered particle
        !     IElement : ID of the considered element
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ParticleIndex, IElement
          ! Local variables
          integer(INTEGER_TYPE) :: I, NodeID, IDim
        
          GlobPosArray(ParticleIndex,:) = 0.0 ! Reset coordinates to zero
          do I = 1, ELEMENTNODES ! Loop over nodes of IElement
            NodeID = iabs(ElementConnectivities(I, IElement))
            do IDim = 1, NVECTOR ! Loop over dimensions of IElement
              GlobPosArray(ParticleIndex,IDim) = &
                GlobPosArray(ParticleIndex,IDim) + &
                NodalCoordinatesUpd(NodeID, IDim) *  &
                ShapeValuesArray(ParticleIndex,I)
            end do
          end do
        
        end subroutine DetermineGlobalFromLocalCoord
        
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
             ElementID = ElementIDArray(ParticleIndex)
             
            EleParticles(ParticleIndex) =  &
              ElementID * GetPartFactor() + ParticleIndex
          end do

          ! Sort EleParticles array, update EleParticlesHelp and NPartEle arrays
          call UpdateParticleHouseKeeping()

        end subroutine AdjustParticleHouseKeeping
        
        subroutine CompleteAddedParticlesData()
        !**********************************************************************
        !
        !    Function:  Initialises particle data that involves data from
        !               surrounding elements.
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Local variables
          integer(INTEGER_TYPE), dimension(:), allocatable :: ConsideredParticles
          integer(INTEGER_TYPE) :: NConsideredParticles, IError, ElementID, &
                     ParticleIndex, IParticle

          do ElementID = 1, Counters%NEl
            if (IsEmptyElement(ElementID)) then ! Loop over elements belonging to holes
              do IParticle = 1, NPartEle(ElementID)
                ParticleIndex = GetParticleIndexFunction(IParticle, ElementID)
                  Particles(ParticleIndex)%IsBoundaryParticle = .false.
              end do

              ! Determine state parameters for newly added particles
              if (allocated(ConsideredParticles)) then
                deallocate(ConsideredParticles, stat = IError)
              end if

              call FindConsideredParticlesSurroundingElements(ElementID, &
                                                 ConsideredParticles, &
                                                 NConsideredParticles)

              call InterpolateParticleStateParameters(ElementID, &
                                                 ConsideredParticles, &
                                                 NConsideredParticles)
              call InterpolateParticleMass(ElementID, &
                                           ConsideredParticles, &
                                           NConsideredParticles)
            end if
          end do

          if (allocated(ConsideredParticles) ) then
            deallocate(ConsideredParticles, stat = IError)
          end if

        end subroutine CompleteAddedParticlesData

        subroutine FindConsideredParticlesSurroundingElements(ElementID, &
                                                   ConsideredParticles, &
                                                   NConsideredParticles)
        !**********************************************************************
        !
        !    Function: Determines a list of particles which can be used to
        !              interpolate the required state variable data of newly
        !              inserted particles in ElementID from elements surrounding
        !              and including ElementID.
        !           
        !     ElementID : ID of the considered empty element
        !      
        ! O   ConsideredParticles : Array to store the IDs of the considered particles
        ! O   NConsideredParticles : Number of particles considered in solving the stress polynomial
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ElementID
          integer(INTEGER_TYPE), dimension(:), allocatable,  &
            intent(out) :: ConsideredParticles
          integer(INTEGER_TYPE), intent(out) :: NConsideredParticles
          ! Local variables
          integer(INTEGER_TYPE) :: IError, IDMat, &
                     NElemPart, Counter, ParticleIndex, &
                     NInnerSurroundingElements, MaterialID, I, J, &
                     AdjacentElement, NSphereParticles
          real(REAL_TYPE) :: SphereRadius, ParticleDistance
          logical :: FoundMaterialParticleInsideSphere, FoundMaterialParticle

          NConsideredParticles = 0
          NSphereParticles = 0

          do IDMat = 1, Counters%NLayers
            if (MaterialElements(IDMat, ElementID)==1) then
              MaterialID = IDMat
              EXIT
            end if
          end do

          NInnerSurroundingElements = GetNElmOfElm(ElementID)

          ! Determine sphere radius
          SphereRadius = DetermineMaxSphereRadius(ElementID)

          FoundMaterialParticle = .false.
          FoundMaterialParticleInsideSphere = .false.
          do I = 1, NInnerSurroundingElements
            AdjacentElement = GetElmIOfElm(ElementID, I)
            NElemPart = NPartEle(AdjacentElement)
            do J = 1, NElemPart
              ParticleIndex = GetParticleIndexFunction(J, AdjacentElement)
              if (MaterialID==MaterialIDArray(ParticleIndex)) then
                if (.not.IsNewVirtualParticle(ParticleIndex)) then
                  ParticleDistance = Length(GlobPosArray(ParticleIndex,:) - ElementCentrePoints(ElementID, 1:NVECTOR), NVECTOR)
                  NConsideredParticles = NConsideredParticles + 1
                  if (ParticleDistance<SphereRadius) then
                    NSphereParticles = NSphereParticles + 1
                    if (IsMaterialParticle(ParticleIndex)) then
                      FoundMaterialParticleInsideSphere = .true.
                    end if
                  end if
                  if (IsMaterialParticle(ParticleIndex)) then
                    FoundMaterialParticle = .true.
                  end if
                end if
              end if
            end do
          end do

          if ((NSphereParticles>0).and.FoundMaterialParticleInsideSphere) then
            NConsideredParticles = NSphereParticles
          else
            call WriteInLogFile('No particles found for stress ' // &
                                'initialisation of virtual particles ' // &
                                'inside sphere around element ' // &
                                trim(String(ElementID)))
          end if
          if (NConsideredParticles==0) then
            call GiveError('No particles found for stress ' // &
                           'initialisation of virtual particles ' // &
                           'in element ' //  trim(String(ElementID)))
          end if
          if (.not.FoundMaterialParticle) then
            call WriteInLogFile('No material particle found for stress ' // &
                                'initialisation of virtual particles ' // &
                                'in element '// trim(String(ElementID)))
          end if
          allocate(ConsideredParticles(NConsideredParticles), stat = IError)
            
          Counter = 0
          NInnerSurroundingElements = GetNElmOfElm(ElementID)
          do I = 1, NInnerSurroundingElements
            AdjacentElement = GetElmIOfElm(ElementID, I)
            NElemPart = NPartEle(AdjacentElement)
            do J = 1, NElemPart
              ParticleIndex = GetParticleIndexFunction(J, AdjacentElement)
              if (MaterialID== &
                  MaterialIDArray(ParticleIndex)) then
                if (.not.IsNewVirtualParticle(ParticleIndex)) then
                  ParticleDistance =  &
                    Length(GlobPosArray(ParticleIndex,:) - ElementCentrePoints(ElementID, 1:NVECTOR), NVECTOR)
                  if ((NSphereParticles==0).or.(.not.FoundMaterialParticleInsideSphere)) then
                    Counter = Counter + 1
                    ConsideredParticles(Counter) = ParticleIndex
                  else if (ParticleDistance<SphereRadius) then
                    Counter = Counter + 1
                    ConsideredParticles(Counter) = ParticleIndex
                  end if
                end if
              end if
            end do
          end do
     
        end subroutine FindConsideredParticlesSurroundingElements

        real(REAL_TYPE) function DetermineMaxSphereRadius(ElementID)
        !**********************************************************************
        !
        !    Function: Determine radius of sphere that extends to the average
        !              distance of the nodes of surrounding elements.
        !           
        !     ElementID : ID of the considered element
        !      
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
    
          integer(INTEGER_TYPE), intent(in) :: ElementID
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, K, NodeID, AdjacentElement, &
                     NInnerSurroundingElements
          real(REAL_TYPE) :: NConsideredNodes
          real(REAL_TYPE) :: Distance, AverageNodeDistance, MaxNodeDistance
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: LocalNodes
          logical :: IsInLocalNodes
    
          DetermineMaxSphereRadius = 0.0
          
          LocalNodes = ElementConnectivities(1:ELEMENTNODES, ElementID)
          
          NInnerSurroundingElements = GetNElmOfElm(ElementID)
          AverageNodeDistance = 0.0
          MaxNodeDistance = 0.0
          NConsideredNodes = 0.0
          do I = 1, NInnerSurroundingElements
            AdjacentElement = GetElmIOfElm(ElementID, I)
            do J = 1, ELEMENTNODES
              if (IsCornerNode(J, ELEMENTNODES)) then ! Loop over corner nodes of element
                NodeID = ElementConnectivities(J, AdjacentElement)
                IsInLocalNodes = .false.
                do K = 1, ELEMENTNODES
                  IsInLocalNodes = (NodeID==LocalNodes(K))
                  if (IsInLocalNodes) then
                    EXIT
                  end if
                end do
                if (.not.IsInLocalNodes) then
                  Distance = Length(NodalCoordinatesUpd(NodeID, 1:NVECTOR) - ElementCentrePoints(ElementID, 1:NVECTOR), NVECTOR)
                  if (Distance>MaxNodeDistance) then
                    MaxNodeDistance = Distance
                  end if
                  AverageNodeDistance = AverageNodeDistance + Distance
                  NConsideredNodes = NConsideredNodes + 1.0
                end if
              end if
            end do
          end do          
          
          AverageNodeDistance = AverageNodeDistance / NConsideredNodes
                  
          DetermineMaxSphereRadius = (AverageNodeDistance + MaxNodeDistance) / 2.0
                  
        end function DetermineMaxSphereRadius       
        
        subroutine InterpolateParticleStateParameters(ElementID, &
                                                     ConsideredParticles, &
                                                      NConsideredParticles)
        !**********************************************************************
        !
        !    Function: Determine state parameters of particles inside ElementID from
        !              ConsideredParticles by averaging of state parameters.
        !           
        !     ElementID : ID of the considered element
        !     ConsideredParticles : Array to store the IDs of the considered particles
        !     NConsideredParticles : Number of particles considered in solving the stress polynomial
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: ElementID, NConsideredParticles
          integer(INTEGER_TYPE), dimension(NConsideredParticles),  &
            intent(in) :: ConsideredParticles
          ! Local variables
          integer(INTEGER_TYPE) :: IPart, ParticleIndex, I
          real(REAL_TYPE), dimension(NTENSOR) :: StressAverage
          real(REAL_TYPE) :: WPAverage, GPAverage, ParticleStress, DPBVDAverage
          real(REAL_TYPE) :: ParticleVol, ParticleVols, AverageWeight
          real(REAL_TYPE) :: ParticleStateVariable, &
                              AverageCohesionStSoft, &
                              AveragePhiStSoft, &
                              AveragePsiStSoft, NNonZeroIntegrationWeight
          real(REAL_TYPE) :: AverageHPStateVariables (2), &
                              AverageModifiedHPStateVariables(2), &
                              AverageHPIGStateVariables (7), &
                              AverageEpsP(NTENSOR), &
                              AverageSigmaPrin(NTENSOR), &
                              AveragePP
          
          real(REAL_TYPE), dimension(NSTATEVAR) :: AverageESMstatev

          AverageWeight = 0.0
          NNonZeroIntegrationWeight = 0.0
          do IPart = 1, NConsideredParticles
            ParticleIndex = ConsideredParticles(IPart) ! Get the particle ID
            if (IsMaterialParticle(ParticleIndex)) then
              AverageWeight = AverageWeight +  &
                              Particles(ParticleIndex)%IntegrationWeight
              NNonZeroIntegrationWeight = NNonZeroIntegrationWeight + 1.0
            end if
          end do

          if (NNonZeroIntegrationWeight>0) then
            AverageWeight = AverageWeight / NNonZeroIntegrationWeight
          else 
            AverageWeight = 0.0
          end if

          if (AverageWeight==0.0) then ! If all particles are virtual particles
            AverageWeight = 1.0
          end if
          
          ParticleVols = 0.0
          ParticleVol = 0.0

          StressAverage = 0.0 ! Initialise average stress of element IEl
          WPAverage = 0.0 ! Initialise average water pressure of element IEl
          GPAverage = 0.0 ! Initialise average gas pressure of element IEl
          DPBVDAverage = 0.0
          AverageHPStateVariables   = 0.0
          AverageHPIGStateVariables = 0.0
          AverageModifiedHPStateVariables = 0.0
          AverageEpsP = 0.0
          AverageSigmaPrin = 0.0
          AverageCohesionStSoft = 0.0
          AveragePhiStSoft = 0.0
          AveragePsiStSoft = 0.0
          AveragePP = 0.0
          ! UMAT
          AverageESMstatev = 0.0

          do IPart = 1, NConsideredParticles ! Loop over particles
            ParticleIndex = ConsideredParticles(IPart)
            call Assert(ParticleIndex >= 0, 'ParticleIndex must be greater than or equal to zero')

            if (.not.IsMaterialParticle(ParticleIndex)) then
              ParticleVol = AverageWeight ! Use average weight
            else  
              ParticleVol = Particles(ParticleIndex)%IntegrationWeight
            end if
                
            ! Sum up stress components, water pressure
            ! Stresses all particles
            do I = 1, NTENSOR
              ParticleStress =  SigmaEffArray(ParticleIndex, I)
              StressAverage(I) = StressAverage(I) +  &
                                 ParticleStress * ParticleVol
            end do
                
            WPAverage = WPAverage + &
                        Particles(ParticleIndex)%WaterPressure * &
                        ParticleVol
            GPAverage = GPAverage + &
                        Particles(ParticleIndex)%GasPressure * &
                        ParticleVol
            DPBVDAverage = DPBVDAverage + &
              Particles(ParticleIndex)%DBulkViscousPressure * ParticleVol
            
            ! Sum up state variables


            ! HP model
            do I = 1, 2
              ParticleStateVariable = GetHPStateVariablesI(Particles(ParticleIndex), I)
              AverageHPStateVariables(I) = AverageHPStateVariables(I) +  ParticleStateVariable * ParticleVol
            end do

            ! HPIG model
            do I = 1, 7
              ParticleStateVariable = GetHPIGStateVariablesI(Particles(ParticleIndex), I) 
              AverageHPIGStateVariables(I) = AverageHPIGStateVariables(I) + ParticleStateVariable * ParticleVol
            end do

            ! modified HP model
            do I = 1, size(AverageModifiedHPStateVariables)
              ParticleStateVariable = GetModifiedHPStateVariablesI(Particles(ParticleIndex), I)
              AverageModifiedHPStateVariables(I) = AverageModifiedHPStateVariables(I) +  ParticleStateVariable * ParticleVol
            end do

            do I = 1, NTENSOR ! Strain softening model
              ParticleStateVariable = &
                GetEpsPI(Particles(ParticleIndex), I)
              AverageEpsP(I) = AverageEpsP(I) + &
                ParticleStateVariable * ParticleVol
              ParticleStateVariable = &
                GetSigmaPrinI(Particles(ParticleIndex), I)
              AverageSigmaPrin(I) = AverageSigmaPrin(I) + &
                ParticleStateVariable * ParticleVol
            end do
            AverageCohesionStSoft = AverageCohesionStSoft + &
              Particles(ParticleIndex)%CohesionStSoft * ParticleVol
            AveragePhiStSoft = AveragePhiStSoft + &
              Particles(ParticleIndex)%PhiStSoft * ParticleVol
            AveragePsiStSoft = AveragePsiStSoft + &
              Particles(ParticleIndex)%PsiStSoft * ParticleVol
            ! Average preconsolidation pressure for MCC
            AveragePP = AveragePP +  &
              Particles(ParticleIndex)%PP * ParticleVol
            
            do I=1, NSTATEVAR !user defined model state variables
                ParticleStateVariable = ESMstatevArray(ParticleIndex,I)
                AverageESMstatev(I) = AverageESMstatev(I) + ParticleStateVariable
            end do
             
            ParticleVols = ParticleVols + ParticleVol
                
          end do
                   
          ! Take the average of each stress/state parameter component, the water pressure
          do I = 1, NTENSOR
            StressAverage(I) = StressAverage(I) /  &
                               ParticleVols
          end do
          WPAverage = WPAverage / ParticleVols
          GPAverage = GPAverage / ParticleVols
          DPBVDAverage = DPBVDAverage / ParticleVols




          do I = 1, 2 ! HP model
            AverageHPStateVariables(I) =  AverageHPStateVariables(I) / ParticleVols
          end do

          do I = 1, 7 ! HPIG model
            AverageHPIGStateVariables(I) = AverageHPIGStateVariables(I) / ParticleVols
          end do

          do I = 1, size(AverageModifiedHPStateVariables) ! modified HP model
            AverageModifiedHPStateVariables(I) =  AverageModifiedHPStateVariables(I) / ParticleVols
          end do

          do I = 1, NTENSOR ! Strain softening model
            AverageEpsP(I) = &
              AverageEpsP(I) / ParticleVols
            AverageSigmaPrin(I) = &
              AverageSigmaPrin(I) / ParticleVols
          end do
          AverageCohesionStSoft = AverageCohesionStSoft / ParticleVols
          AveragePhiStSoft = AveragePhiStSoft / ParticleVols
          AveragePsiStSoft = AveragePsiStSoft / ParticleVols
          ! MCC model
          AveragePP = AveragePP / ParticleVols
          do I=1, NSTATEVAR !user defined model state variables
             AverageESMstatev(I) = AverageESMstatev(I) / ParticleVols
          end do

          ! Assign the average stress, water pressure to particles of ElementID
          do IPart = 1, NPartEle(ElementID)
            ParticleIndex = GetParticleIndexFunction(IPart, ElementID) ! Get the particle ID
            
            SigmaEffArray(ParticleIndex,:) = StressAverage(:)
            
            Particles(ParticleIndex)%WaterPressure = WPAverage
            Particles(ParticleIndex)%DBulkViscousPressure = DPBVDAverage



            call SetHPStateVariables(Particles(ParticleIndex),  AverageHPStateVariables)
            call SetHPIGStateVariables(Particles(ParticleIndex), AverageHPIGStateVariables)
            call SetModifiedHPStateVariables(Particles(ParticleIndex),  AverageModifiedHPStateVariables)

            call SetEpsP(Particles(ParticleIndex), AverageEpsP)
            call SetSigmaPrin(Particles(ParticleIndex), AverageSigmaPrin)
            Particles(ParticleIndex)%CohesionStSoft = AverageCohesionStSoft
            Particles(ParticleIndex)%PhiStSoft = AveragePhiStSoft
            Particles(ParticleIndex)%PsiStSoft = AveragePsiStSoft
            Particles(ParticleIndex)%PP = AveragePP
            !User defined material model
            ESMstatevArray(ParticleIndex,:) = AverageESMstatev
          end do

        end subroutine InterpolateParticleStateParameters

        subroutine InterpolateParticleMass(ElementID, &
                                           ConsideredParticles, &
                                           NConsideredParticles)
        !**********************************************************************
        !
        !    Function: Determine mass of particles inside ElementID from
        !              ConsideredParticles by averaging of state parameters.
        !           
        !     ElementID : ID of the considered element
        !     ConsideredParticles : Array to store the IDs of the considered particles
        !     NConsideredParticles : Number of particles considered
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: ElementID, NConsideredParticles
          integer(INTEGER_TYPE), dimension(NConsideredParticles), &
            intent(in) :: ConsideredParticles
          ! Local variables
          integer(INTEGER_TYPE) :: IPart, ParticleIndex
          real(REAL_TYPE) :: MassAverage, MassAverageW

          MassAverage = 0.0 ! Initialise average mass of element IEl
          MassAverageW = 0.0
          do IPart = 1, NConsideredParticles ! Loop over particles
            ParticleIndex = ConsideredParticles(IPart)
            call Assert(ParticleIndex >= 0, 'ParticleIndex must be greater than or equal to zero')

            ! Sum up mass
            MassAverage = MassAverage + &
                          MassArray(ParticleIndex)
            MassAverageW = MassAverageW + &
                          MassWaterArray(ParticleIndex)

          end do

          ! Take the average of the mass
          MassAverage = CalParams%VirtualParticleMassFactor *  &
            MassAverage / dble(NConsideredParticles)
          MassAverageW = CalParams%VirtualParticleMassFactor *  &
            MassAverageW / dble(NConsideredParticles)
                     
          ! Assign the average mass to virtual particles of ElementID
          do IPart = 1, NPartEle(ElementID)
            ParticleIndex = GetParticleIndexFunction(IPart, ElementID) ! Get the particle ID
            MassArray(ParticleIndex) = MassAverage
            MassWaterArray(ParticleIndex) = MassAverageW
          end do

        end subroutine InterpolateParticleMass
        
        subroutine AverageElementParticleStateParam(ElementID)
        !**********************************************************************
        !
        !    Function: Update state parameters of particles inside ElementID
        !              by averaging of state parameters.
        !           
        !     ElementID : ID of the considered element
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ElementID
          ! Local variables
          integer(INTEGER_TYPE) :: NConsideredParticles, IError, I
          integer(INTEGER_TYPE), dimension(:), allocatable :: ConsideredParticles

          NConsideredParticles = NPartEle(ElementID)
          allocate(ConsideredParticles(NConsideredParticles),  &
                   stat = IError)
          do I = 1, NConsideredParticles
            ConsideredParticles(I) = GetParticleIndexFunction(I, ElementID)
          end do

          call InterpolateParticleStateParameters(ElementID, &
                                                  ConsideredParticles, &
                                                  NConsideredParticles)

          deallocate(ConsideredParticles, stat = IError)

        end subroutine AverageElementParticleStateParam

      end module ModAdjustParticleDiscretisation
