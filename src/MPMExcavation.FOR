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


      module ModMPMExcavation
      !**********************************************************************
      !
      ! Function: This module contains all routines which are used for the excavation tool
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      use ModCounters
      use ModReadCalculationData
      use ModReadGeometryData
      use ModMPMData
      use ModParticle
      use ModMeshInfo
      use ModAdjustParticleDiscretisation
      use ModGlobalConstants

      contains ! routines of this module


        subroutine ApplyExcavation()
        !**********************************************************************
        !
        ! Function: Applies excavation
        !
        !**********************************************************************    
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: I, IElement ! counters
          integer(INTEGER_TYPE) :: ExcavVolume, IniExcavIStep, EndExcavIStep, Vol
        
          if(.not.CalParams%ApplyExcavation) RETURN
                
          do I = 1, CalParams%NumberExcavatedVolumes
            ExcavVolume = CalParams%ExcavatedVolumeID(I)
            IniExcavIStep = CalParams%ExcavationStages(I,1)
            EndExcavIStep = CalParams%ExcavationStages(I,2)            
            
            if ( (CalParams%IStep >= IniExcavIStep) .and. (CalParams%IStep <= EndExcavIStep) ) then
              do IElement = 1, Counters%NEl ! loop over all elements
                Vol = GeoParams%ExcavatedElements(IElement) ! this matrix comes from GOM file
                if (Vol == ExcavVolume) then
                  call RemoveParticlesFromExcavationElement(IElement)
                end if
              end do
            end if
            
          end do
        
          call UpdateParticleHouseKeepingMock()
          call SetActiveElementMock()  
          call SetParticleIndex()
        
        end subroutine ApplyExcavation
        
       
        subroutine RemoveParticlesFromExcavationElement(IElement)
        !****************************************************************************************
        !
        ! Function: Removes material points from the excavated elements and update housekeeping
        ! 
        ! IElement: Element ID of the element to be excavated
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
            Particles(ParticleIndex)%Kind = REMOVEDPARTICLE
          end do
        
          NPartEle(IElement) = 0 ! partices inside are to be removed, NPartEle(IElement) is set to be 0 
          IsActiveElement(IElement) = .false. ! deactivate the element 
        
          NegativeNRemovedParticles = 0 - NRemovedParticles
        
         ! remove the information of particles to be removed from array Particles, and reset incorrect and redundant particle house keeping information
         call AdjustParticlesArrays(NegativeNRemovedParticles)
        
         ! adjust the remaining particle housekeeping data structure to the new discretisation
         call AdjustParticleHouseKeeping()
        
        end subroutine RemoveParticlesFromExcavationElement
    
          
      end module ModMPMExcavation
        

        
        