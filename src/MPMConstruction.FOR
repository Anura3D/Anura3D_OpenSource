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
    !	Copyright (C) 2024  Members of the Anura3D MPM Research Community 
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


      module ModMPMConstruction
      !**********************************************************************
      !
      ! Function: This module contains all routines which are used for the construction tool
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 07:56:02 -0500 (Thu, 14 Apr 2022) $
      !
      !**********************************************************************
      use ModCounters
      use ModMPMInit
      use ModReadCalculationData
      use ModReadGeometryData
      use ModMPMData
      use ModParticle
      use ModMeshInfo
      use ModAdjustParticleDiscretisation
      use ModGlobalConstants

      contains ! routines of this module


        subroutine ApplyConstruction()
        !**********************************************************************
        !
        ! Function: Applies construction stages
        !
        !**********************************************************************    
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: I, J, k, IElement, IAEl, IEl, NElemPart, IPart ! counters
          integer(INTEGER_TYPE) :: ParticleIndex, MaterialID ! counters
          integer(INTEGER_TYPE) :: ConstVolume, IniConstIStep, EndConstIStep, Vol

          if(.not.CalParams%ApplyConstruction) RETURN        
        
          do I = 1, CalParams%NumberConstructionVolumes
            ConstVolume = CalParams%ConstructionVolumeID(I)
            IniConstIStep = CalParams%ConstructionStages(I,1) 
            EndConstIStep = CalParams%ConstructionStages(I,2)   

            ! Remove all the solids to be constructed
            if (CalParams%IStep < IniConstIStep) then
              do IElement = 1, Counters%NEl ! loop over all elements
                Vol = GeoParams%ConstructionElements(IElement) ! this matrix comes from GOM file
                if (Vol == ConstVolume) then
                  call RemoveParticlesFromNonConstructedElement(IElement)
                end if
              end do

            ! Check if previous constructed elements are still filled
            ! If refill of empty elements is activated 
            elseif ( (CalParams%IStep >= IniConstIStep) .and. (CalParams%IStep <= EndConstIStep).and.(CalParams%ApplyFillEmptyElements) ) then
              do IElement = 1, Counters%NEl ! loop over all elements
                Vol = GeoParams%ConstructionElements(IElement) ! this matrix comes from GOM file
                if ((Vol == ConstVolume).and.(.not.IsActiveElement(IElement))) then ! Check if previous elements are active
                  call AddParticlesToConstructedElement(IElement)
                end if
              end do
            ! If refill of empty elements is NOT activated
            elseif ( (CalParams%IStep == IniConstIStep) .and.(.not.(CalParams%ApplyFillEmptyElements)) ) then
              do IElement = 1, Counters%NEl ! loop over all elements
                Vol = GeoParams%ConstructionElements(IElement) ! this matrix comes from GOM file
                if ((Vol == ConstVolume).and.(.not.IsActiveElement(IElement))) then ! Check if previous elements are active
                  call AddParticlesToConstructedElement(IElement)
                end if
              end do                                       
               
            end if            
          end do
          
          ! Set the MaterialID of elements
          do IAEl = 1,Counters%NAEL                             ! loop through active elements
            IEl = ActiveElement(IAEl)
            NElemPart = NPartEle(IEl)
            do IPart= 1, NElemPart                              ! loop through particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl)      ! Get the particle ID
              MaterialID = MaterialIDArray(ParticleIndex)       ! Get the material ID
              MaterialElements(MaterialID,IEl) = 1
            end do
           end do
                    
          call UpdateParticleHouseKeepingMock()
          call SetActiveElementMock()  
          call SetParticleIndex()
          call ReinitialiseMaterialPointOutputFiles()
          
        
        end subroutine ApplyConstruction
        
       
        subroutine RemoveParticlesFromNonConstructedElement(IElement)
        !****************************************************************************************
        !
        ! Function: Removes material points from the non-constructed elements and update housekeeping
        ! 
        ! IElement: ID of the element that is not to be constructed in the current load step
        !
        !*****************************************************************************************
        
        implicit none 
        
          integer(INTEGER_TYPE), intent (in) :: IElement
          ! local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex ! counters
          integer(INTEGER_TYPE) :: NRemovedParticles ! number of material points to be removed 
          integer(INTEGER_TYPE) :: NegativeNRemovedParticles ! NegativeNRemovedParticles = 0 - NRemovedParticles
        
          NRemovedParticles = 0
        
          if (NPartEle(IElement) <= 0) RETURN ! there is no MP
  
          NRemovedParticles = NPartEle(IElement) 
        
          do IParticle = 1, NPartEle(IElement)
            ParticleIndex = GetParticleIndexFunction(IParticle, IElement)
            Particles(ParticleIndex)%Kind = ADDEDPARTICLE
          end do
        
          NPartEle(IElement) = 0 ! partices inside have to be removed, NPartEle(IElement) is set to 0 
          IsActiveElement(IElement) = .false. ! deactivate the element 
        
          NegativeNRemovedParticles = 0 - NRemovedParticles
        
         ! remove the information of particles to be removed from array Particles, and reset incorrect and redundant particle house keeping information
         call AdjustParticlesArrays(NegativeNRemovedParticles)
        
         ! adjust the remaining particle housekeeping data structure to the new discretisation
         call AdjustParticleHouseKeeping()
        
        end subroutine RemoveParticlesFromNonConstructedElement
    
        subroutine AddParticlesToConstructedElement(IElement)
        !****************************************************************************************
        !
        ! Function: Add elements and particles of the construction volume and update housekeeping
        ! 
        ! IElement: ID of the element to be constructed in the current load step
        !
        !*****************************************************************************************
        
        implicit none 
        
          integer(INTEGER_TYPE), intent (in) :: IElement
          ! local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex ! counters
          integer(INTEGER_TYPE) :: NAddParticles ! number of material points to be add to the constructed elements 
          integer(INTEGER_TYPE) :: PositiveNAddedParticles, maxpartIndex ! 
        
          !GetParticleIndex = InitialGetParticleIndex
          NAddParticles = InitialNPartEle(IElement)
          maxpartIndex = maxval(GetParticleIndex)
          
          PositiveNAddedParticles = NAddParticles ! Number of Added Particles
          
          ! Add the information of particles to be added to array Particles, and update particle house keeping information
          NPartEle(IElement) = NAddParticles    ! update number of particles associated with the IElement 
          IsActiveElement(IElement) = .true.    ! activate the element
          Counters%NAEl = Counters%NAEl + 1     ! increase the counter of active elements
          call SetActiveElement()
          call AdjustParticlesArraysConstruction(PositiveNAddedParticles,IElement)
        
         ! adjust the remaining particle housekeeping data structure to the new discretisation
         call AdjustParticleHouseKeeping()
        
        end subroutine AddParticlestoConstructedElement

        subroutine ReinitialiseMaterialPointOutputFiles()
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
              if ((CalParams%OutputParticles(I)>0).and.(CalParams%OutputParticles(I) <= size(IDArray))) then
                call GetStepExt(IDArray(CalParams%OutputParticles(I)), XXX)
                CompleteFileName = trim(CalParams%FileNames%ProjectName)// '.PAR_' // XXX 
                if ( .not.FExist(trim(CompleteFileName)) ) then

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
                            'LoadStep ', 'TimeStep ', 'Time', 'PA1 ', 'PA2 ', 'PGravity ', 'ID_MP ', 'X ', 'Y ', &
                            'Ux ', 'Uy ', 'SigmaXX ', 'SigmaYY ', 'SigmaZZ ','SigmaXY ', &
                            'WPressure', 'EpsilonXX ', 'EpsilonYY ', 'EpsilonZZ ','GammaXY ', &
                            'IntWeight ', 'MatID '
                  end select

                end if  
              end if
            end do
          end if

        end subroutine ReinitialiseMaterialPointOutputFiles
        
      end module ModMPMConstruction
        

        
        