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
    !	Copyright (C) 2022  Members of the Anura3D MPM Research Community 
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
	
	
	  module ModMPMData
      !**********************************************************************
      !
      !    Function:  Contains the housekeeping data structure required for the MPM project.
      !
      !               In order to keep the size of this source file reasonably small,
      !               this module only contains routines that are directly related to
      !               the basic manipulation of contained data (maintenance of data structures).
      !
      !               First, routines are provided to access the material point data, followed by routines
      !               for generating the material point house-keeping data structure.
      !
      !               Routines for initialising or updating material point information or 
      !               routines for output of material point data are located in separate modules that access the
      !               data stored in this module.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr. 2022) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluation
      use ModMeshAdjacencies
      use ModParticle
      use ModMeshInfo
      !use NURBS
      
      implicit none

        !> Decimal factor used for defining particle - element connectivities with the array EleParticles
        integer(kind = 8), private, save  :: PartFactor = 10
        
        !> NPartEle(IElement) = Number of particles in element IElement, size = NEl
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: NPartEle !(rank 2 for multipatch)

        !> NSolidEle(IElement) = Number of solid material points in element IElement, size = NEl
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: NSolidEle !(rank 2 for multipatch)

        !> NLiquidEle(IElement) = Number of liquid material points in element IElement, size = NEl
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: NLiquidEle !(rank 2 for multipatch)

        !> List which contains the particle connectivities (particle - element assigments), size = NParticles
        integer(kind = 8), dimension(:), allocatable :: EleParticles

        !> Helper list for EleParticles, index of first particle of element, size = NEl
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: EleParticlesHelp ! (rank 2 for multipatch)

        !> True, if particle based integration is used inside element IElement, size = NEl
        logical, dimension(:, :), allocatable :: IsParticleIntegration ! rank 2 for multipatch generalization
        
        !> True, if liquid free surface Material Point is inside element IElement, size = NEl
        logical, dimension(:), allocatable :: IsElemWithLiquidFreeSurfMP
        
        type(ParticleType(:, :, :)), dimension(:), allocatable :: Particles ! List of particle objects
        
        !> The following arrays correspond to the stripped members of the Particle datatype to improve efficency (spatial locality)  
        integer(INTEGER_TYPE), allocatable          :: IDArray(:)                   !> Unique ID of a particle. This variable should only be set once through the routine InitParticle.
        integer(INTEGER_TYPE), allocatable          :: ElementIDArray(:)            !> Unique ID of the element to whom the particle belongs. A value of -1 indicates no assignment.
        integer(INTEGER_TYPE), allocatable          :: PatchIDArray(:)              !> Patch ID of the element to whon the particle belongs. A value of -1 indicates no assignment. 
        integer(INTEGER_TYPE), allocatable          :: EntityIDArray(:)             !> Unique ID indicating to which entity the particle belongs
        integer(INTEGER_TYPE), allocatable          :: MaterialIDArray(:)           !> Material ID of the particle referring to the material data read from input files
        
        integer(INTEGER_TYPE), allocatable          :: MaterialPointTypeArray(:)    !> Material properties (properties related to the matter represented by particles)
        real(REAL_TYPE), allocatable :: MassArray(:)                 !> Particle mass [kg] (corresponds to the dry density of soil) (fixed during simulation - ensures conservation of mass)
        real(REAL_TYPE), allocatable :: MassWaterArray(:)            !> Particle mass [kg] (corresponds to the water density)(fixed during simulation - ensures conservation of mass)
        
        ! Displacements (solid)
        real(REAL_TYPE), allocatable :: UArray(:,:)                  !> Total displacements in x, y, z direction of the particle
        real(REAL_TYPE), allocatable :: UPhaseArray(:,:)             !> Displacements in x, y, z direction of the particle during a phase
        real(REAL_TYPE), allocatable :: UStepArray(:,:)              !> Displacements in x, y, z direction of the particle during a load step
       
        real(REAL_TYPE), allocatable :: ShapeValuesArray(:,:)        !> Shape function values for the particle, these values are updated when new local coordinates are calculated.
        !Shape function values for NURBS for debugging purposes I am plotting the Xi and Eta shape functions 
        ! -----------------------------------------------------------------------
        
        real(REAL_TYPE), allocatable :: ShapeValuesArray_Xi(:,:)
        real(REAL_TYPE), allocatable :: ShapeValuesArray_Eta(:,:)
        real(REAL_TYPE), allocatable :: ShapeValuesArray_Zeta(:,:)

        
        real(REAL_TYPE), allocatable :: DShapeValuesArray_Xi(:,:)
        real(REAL_TYPE), allocatable :: DShapeValuesArray_Eta(:,:)
        real(REAL_TYPE), allocatable :: DShapeValuesArray_Zeta(:,:)

        ! we honestly do not need this here... weights are not multidimensional
        real(REAL_TYPE), allocatable :: WeightArray_Xi(:) 
        real(REAL_TYPE), allocatable :: WeightArray_Eta(:)
        real(REAL_TYPE), allocatable :: WeightArray_Zeta(:)
        
        ! -----------------------------------------------------------------------
        
        !real(REAL_TYPE), allocatable :: NN_IncludesZeroValues_Print(:)
        !real(REAL_TYPE), allocatable :: dN_dxi_IncludesZeroValues_Print(:)
        !
        !real(REAL_TYPE), allocatable :: MM_IncludesZeroValues_Print(:)
        !real(REAL_TYPE), allocatable :: dM_deta_IncludesZeroValues_Print(:)
        
        ! -----------------------------------------------------------------------

        
        
        
        
        real(REAL_TYPE), allocatable :: GlobPosArray(:,:)            !> Global particle position x = GlobPos(1), y = GlobPos(2), z = GlobPos(3)
        
        real(REAL_TYPE), allocatable :: VelocityArray(:,:)           ! Velocity of the particle (soil)
        real(REAL_TYPE), allocatable :: VelocityWaterArray(:,:)      ! Velocity of the particle (water)
        real(REAL_TYPE), allocatable :: VelocityGasArray(:,:)        ! Velocity of the particle (gas)
        
        ! Stresses
        real(REAL_TYPE), allocatable :: SigmaEffArray(:,:)           !> Effective stresses (Sxx, Syy, Szz, Sxy, Syz, Sxz) of the particle 
        real(REAL_TYPE), allocatable :: SigmaEff0Array(:,:)          !> Initial Effective stresses (Sxx, Syy, Szz, Sxy, Syz, Sxz) of the particle from the last load step
        
        real(REAL_TYPE), allocatable :: DShapeValuesArray(:,:,:)     !> Shape function derivatives for the particle, these values are updated when new local coordinates are calculated.
        real(REAL_TYPE), allocatable :: AccelerationArray(:,:)       !> Acceleration of the particle (soil)
        
        ! User defined model
        real(REAL_TYPE), allocatable :: ESMstatevArray(:,:)         !> array of state variables of the user defined material model

      contains ! Routines of this module


        logical function IsMaterialParticle(MaterialPointIndex)
        !**********************************************************************
        !
        !    Function:  returns .true. if the material point with ID MaterialPointIndex is
        !               a material point
        !
        !    MaterialPointIndex : ID of the considered material point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialPointIndex
          
          IsMaterialParticle =  Particles(MaterialPointIndex)%Kind==MATERIALPARTICLE
          
        end function IsMaterialParticle


        logical function IsVirtualParticle(MaterialPointIndex)
        !**********************************************************************
        !
        !    Function:  Returns .true. if the material point with ID MaterialPointIndex is
        !               a virtual material point
        !
        !    MaterialPointIndex : ID of the considered material point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialPointIndex
          
          IsVirtualParticle = Particles(MaterialPointIndex)%Kind==VIRTUALPARTICLE
          
        end function IsVirtualParticle


        logical function IsRemovedParticle(MaterialPointIndex)
        !**********************************************************************
        !
        !    Function:  Returns .true. if the material point with ID MaterialPointIndex is
        !               a material point to be removed
        !
        !    MaterialPointIndex : ID of the considered material point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialPointIndex
          
          IsRemovedParticle = Particles(MaterialPointIndex)%Kind==REMOVEDPARTICLE
          
        end function IsRemovedParticle


        logical function IsNewVirtualParticle(MaterialPointIndex)
        !**********************************************************************
        !
        !    Function:  Returns .true. if the material point with ID MaterialPointIndex is
        !               a newly initialised virtual material point
        !
        !    MaterialPointIndex : ID of the considered material point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialPointIndex
          
          IsNewVirtualParticle = Particles(MaterialPointIndex)%Kind==NEWVIRTUALPARTICLE
          
        end function IsNewVirtualParticle


        integer(INTEGER_TYPE) function NMaterialParticlesEle(IElement)
        !**********************************************************************
        !
        !    Function:  Returns the number of material points in IElement
        !
        !    IElement : Considered element
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: I, MaterialPointIndex
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary

          NMaterialParticlesEle = 0
          do I = 1, NPartEle(IElement,IPatch_Temporary)
            MaterialPointIndex = GetParticleIndex(I, IElement, IPatch_Temporary)
            
            if (IsMaterialParticle(MaterialPointIndex)) then
              NMaterialParticlesEle = NMaterialParticlesEle + 1
            end if

          end do
        
        end function NMaterialParticlesEle


        integer(INTEGER_TYPE) function NVirtualParticlesEle(IElement)
        !**********************************************************************
        !
        !    Function:  Returns the number of virtual material point in IElement
        !
        !    IElement : Considered element
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: I, ParticleIndex
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary

          NVirtualParticlesEle = 0
          do I = 1, NPartEle(IElement,IPatch_Temporary)
            ParticleIndex = GetParticleIndexFunction(I, IElement, IPatch_Temporary)
            
            if (IsVirtualParticle(ParticleIndex) ) then
              NVirtualParticlesEle = NVirtualParticlesEle + 1
            end if
          end do

        end function NVirtualParticlesEle

        integer(INTEGER_TYPE) function NRemovedParticlesEle(IElement)
        !**********************************************************************
        !
        !    Function:  Returns the number of material points to be removed in IElement
        !
        !    IElement : Considered element
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: I, ParticleIndex
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1

          NRemovedParticlesEle = 0
          do I = 1, NPartEle(IElement,IPatch_Temporary)
            ParticleIndex = GetParticleIndex(I, IElement, IPatch_Temporary)
            
            if (IsRemovedParticle(ParticleIndex) ) then
              NRemovedParticlesEle = NRemovedParticlesEle + 1
            end if
          end do

        end function NRemovedParticlesEle


        logical function HasBoundaryParticle(IElement)
        !**********************************************************************
        !
        !    Function:  Returns .true. if IElement contains a material point marked as
        !               boundary material point
        !
        !     IElement : Considered element
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        
          HasBoundaryParticle = .false.
          do IParticle = 1, NPartEle(IElement,IPatch_Temporary) ! Loop over particles in IElement
            ParticleIndex = GetParticleIndexFunction(IParticle, IElement, IPatch_Temporary)
      
            if (Particles(ParticleIndex)%IsBoundaryParticle) then
              HasBoundaryParticle = .true.
              EXIT
            end if
          end do
        
        end function HasBoundaryParticle


        integer(INTEGER_TYPE) function GetParticleIndexFunction(IParticle, IElement, IPatch) result(res)
        !**********************************************************************
        !
        !    Function:  Returns the index of the material point with the local number
        !               IParticle inside element IElement. The returned index
        !               is the location of the material point in the array 'Particles'.
        !               Note that the index of the material point might not be identical
        !               to the particle ID (if material points might be deleted ... ).
        !               If less than IParticle material points are inside IElement, -1 is
        !               returned.
        !
        !     IParticle : Local number of the material point inside IElement
        !     IElement : ID of the element, that the considered material point is located in
        !
        ! O   GetParticleIndex : Index of the considered material point in the array 'Particles'
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: EleParticleIndex
          integer(kind = 8) :: LocalIElement
          ! Multipatch variables 
          integer(INTEGER_TYPE), intent(in) :: IPatch!_Temporary = 1
          
          if ( (allocated(EleParticlesHelp) ).and.(allocated(EleParticles) ).and. &
               (IParticle>=1).and.(IParticle<=NPartEle(IElement,IPatch))            ) then
            EleParticleIndex = EleParticlesHelp(IElement, IPatch) + IParticle - 1
            LocalIElement = IElement
            res = EleParticles(EleParticleIndex) - LocalIElement * PartFactor
          else
            res = -1
          end if

        end function GetParticleIndexFunction

        subroutine SetParticleIndex()
        !**********************************************************************
        !
        !    Function:  Sets/updates GetParticleIndex.
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        integer(INTEGER_TYPE) :: maxParticle, iError, IAEl, iEl, int, maxEl

        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        integer(INTEGER_TYPE), allocatable, dimension(:,:,:) :: GetParticleIndex_Temporary
        !integer(INTEGER_TYPE) :: maxParticle_previous = 64
        integer(INTEGER_TYPE) :: counter = 1
        integer(INTEGER_TYPE) :: StartingValue
        
        maxParticle = 1
        maxEl       = 1
        do IPatch = 1, NumberOfPatches ! Loop over patches
                do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
        !do IAEl = 1, Counters%NAEl
          iEl = ActiveElement(IAEl, IPatch)
          maxEl = max(maxEl, iEl)
          maxParticle = max(maxParticle, NPartEle(iEl,IPatch))
                end do ! loop over active elements
        end do ! loop over patches
        
        iError = 0

        !!!!!!!!!!!!!!!!!!!!
        !Counters%NParticles
        !allocate(GetParticleIndex_Temporary(maxParticle_previous, maxEl, NumberOfPatches), stat = IError)
        !allocate(GetParticleIndex_Temporary(Counters%NParticles, maxEl, NumberOfPatches), stat = IError)
        allocate(GetParticleIndex_Temporary(Counters%NParticles, Maximum_nel_NURBS, NumberOfPatches), stat = IError)

        GetParticleIndex_Temporary = -huge(GetParticleIndex_Temporary)
        
        if (allocated(GetParticleIndex)) then
        do IPatch = 1, NumberOfPatches ! Loop over patches
            do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
            
          
                iEl = ActiveElement(IAEl, IPatch)
                do int = 1, Counters%NParticles!maxParticle_previous!NPartEle(iEl,IPatch)
                    if (GetParticleIndex(int, iEl, IPatch)>0) then 
                    GetParticleIndex_Temporary(counter, iEl, IPatch) = GetParticleIndex(int, iEl, IPatch)
                    counter = counter + 1
                    end if
                enddo
               counter = 1
            enddo
            counter = 1
        end do!loop over patches
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (allocated(GetParticleIndex)) then
          deallocate(GetParticleIndex, stat = IError)
          call DeAllocationError(iError, 'GetParticleIndex', 'SetParticleIndex')
        endif

        !allocate(GetParticleIndex(maxParticle, maxEl), stat = IError)
        !allocate(GetParticleIndex(maxParticle, maxEl, NumberOfPatches), stat = IError)
        !allocate(GetParticleIndex(Counters%NParticles, maxEl, NumberOfPatches), stat = IError)
        allocate(GetParticleIndex(Counters%NParticles, Maximum_nel_NURBS, NumberOfPatches), stat = IError)

        !Counters%NParticles
        call AllocationError(iError, 'GetParticleIndex', 'SetParticleIndex')
        GetParticleIndex = -huge(GetParticleIndex)

        !do IAEl = 1, Counters%NAEl
        
        StartingValue = 0
        if ((CalParams%IStep == 1) .and. (CalParams%TimeStep == 1) .and. (FirstTimeGetParticleIndex == .false.)) then 
        !if (.not.allocated(GetParticleIndex)) then
        do IPatch = 1, NumberOfPatches ! Loop over patches
                do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
            
          iEl = ActiveElement(IAEl, IPatch)
          do int = 1, NPartEle(iEl,IPatch)
            GetParticleIndex(int, iEl, IPatch) = GetParticleIndexFunction(int, iEl, IPatch) !+ StartingValue
          enddo
                enddo
                
                !StartingValue = GetParticleIndex(int-1, iEl, IPatch)
                
        end do!loop over patches
        
        FirstTimeGetParticleIndex = .true.
        
        else 
            
          !  do IPatch = 1, NumberOfPatches ! Loop over patches
          !      do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
          !  
          !iEl = ActiveElement(IAEl, IPatch)
          !do int = 1, NPartEle(iEl,IPatch)
              
              
              !GetParticleIndex(int, iEl, IPatch) = GetParticleIndex_Temporary(int, iEl, IPatch)
              GetParticleIndex = GetParticleIndex_Temporary

              
              

          !enddo
          !
          !      enddo
                
        !end do!loop over patches
            
            
            
        end if 
        
        !maxParticle_previous = maxParticle
        
        
        end subroutine SetParticleIndex


        integer(INTEGER_TYPE) function GetElementIDFromList(IParticle)
        !**********************************************************************
        !
        !    Function:  Returns the Element ID of the particle stored at index
        !               IParticle in array EleParticles.
        !               
        !     IParticle : Index of list EleParticles
        !
        ! O   GetElementIDFromList : ID of the element that particle at index IParticle is located in
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle

          GetElementIDFromList = int(EleParticles(IParticle) /  &
                                     PartFactor)
        
        end function GetElementIDFromList

        integer(INTEGER_TYPE) function GetParticleIndexFromList(IParticle)
        !**********************************************************************
        !
        !    Function:  Returns the index of the particle for the location of the 
        !               particle inside the list EleParticles. The returned index
        !               is the location of the particle in the array 'Particles'.
        !               Note: This routine is used for updating the particles, where
        !                     a loop over the items of EleParticles is performed.
        !               
        !     IParticle : Index of list EleParticles
        !
        ! O   GetParticleIndexFromList : Index of the considered particle in the array 'Particles'
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          
          if ( (allocated(EleParticles) ).and. &
               (IParticle>=1) ) then
            GetParticleIndexFromList = mod(EleParticles(IParticle), &
                                           PartFactor)
          else
            GetParticleIndexFromList = -1
          end if
        
        end function GetParticleIndexFromList


        integer(kind = 8) function GetPartFactor()
        !**********************************************************************
        !
        !    Function:  Returns PartFactor, which is declared private as it should only
        !               be modified inside this module.
        !
        ! O   GetPartFactor : Decimal factor used internally for storing particle-element connectivities
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
          
          GetPartFactor = PartFactor
        
        end function GetPartFactor


        subroutine SetPartFactor(NewPartFactor)
        !**********************************************************************
        !
        !    Function:  Set PartFactor, which is declared private as it should only
        !               be modified inside this module. PartFactor stores the decimal factor 
        !               used internally for storing particle-element connectivities.
        !
        !     NewPartFactor : Decimal factor read from the BRF file of the previous load phase
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(kind = 8), intent(in) :: NewPartFactor
        
          PartFactor = NewPartFactor
        
        end subroutine SetPartFactor


        integer(INTEGER_TYPE) function NumberOfIntegrationPoints(IElement, IPatch)
        !**********************************************************************
        !
        !    Function: Returns the number of particles in element IElement if the element is 
        !              partially filled and the number of Gauss points per element if the 
        !              element is fully filled.
        !
        !     IElement : ID of the considered element
        !
        ! O   NumberOfIntegrationPoints : Number of Gauss points or particles
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
     
        implicit none
        
          integer(INTEGER_TYPE) :: IElement
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE), intent(in) :: IPatch
        
          if (IsParticleIntegration(IElement, IPatch) ) then ! Partially filled element
            NumberOfIntegrationPoints = NPartEle(IElement,IPatch) ! Number of particles in IElement
          else ! Fully filled element
            NumberOfIntegrationPoints = ELEMENTGAUSSPOINTS ! Number of Gauss points in element
          end if
        
        end function NumberOfIntegrationPoints


        subroutine SetEleParticles(IElement, IParticle, ParticleIndex, IPatch)
        !**********************************************************************
        !
        !    Function:  Sets EleParticles(ParticleIndex) and EleParticlesHelp(IElement)
        !               if ParticleIndex is larger than 0. Else EleParticlesHelp(IElement)
        !               is set to -1 to indicate that IElement is inactive and therefor
        !               contains no particles. The particle-element connectivity stored
        !               in EleParticles consists of the element ID multiplied by PartFactor
        !               and ParticleIndex. EleParticlesHelp contains the ParticleIndex of
        !               the first particle located in IElement. The following NPartEle(IElement)
        !               particles of EleParticles also belong to IElement.
        !
        !     IElement : ID of the considered element
        !     IParticle : Local ID of the considered particle inside IElement
        !     ParticleIndex : Index of the considered particle in the Particles array
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
          implicit none
          
            integer(INTEGER_TYPE), intent(in) :: IElement
            integer(INTEGER_TYPE), intent(in) :: IParticle
            integer(INTEGER_TYPE), intent(in) :: ParticleIndex
            ! Local variables
            integer(kind = 8) :: LocalIElement
            integer(kind = 8) :: LocalParticleIndex
            
            ! Multipatch variables 
            integer(INTEGER_TYPE) :: IPatch_Temporary = 1
            integer(INTEGER_TYPE), intent(in) :: IPatch

            
            if (ParticleIndex>0) then
              LocalIElement = IElement
              LocalParticleIndex = ParticleIndex

              EleParticles(ParticleIndex) = LocalIElement * PartFactor + &
                                            LocalParticleIndex
     
          
              if (IParticle==1) then ! First particle of element
                EleParticlesHelp(IElement, IPatch) = ParticleIndex
              end if
            else ! Inactive element
              EleParticlesHelp(IElement, IPatch) = -1
            end if
        
        end subroutine SetEleParticles


        subroutine SetParticleElementID(IParticle, ParticleIndex, NewElementID)
        !**********************************************************************
        !
        !    Function:  Sets EleParticles(IParticle) and Particles(ParticleIndex)%ElementID.
        !               NOTE: Only this routine should be used to update these two values in order
        !               to guarantee that they are identical!
        !
        !     IParticle : Index of the considered particle in the EleParticles array
        !     ParticleIndex : Index of the considered particle in the Particles array
        !     NewElementID : ID of the new element of the considered particle
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          integer(INTEGER_TYPE), intent(in) :: ParticleIndex
          integer(INTEGER_TYPE), intent(in) :: NewElementID
          ! Local variables
          integer(kind = 8) :: LocalNewElementID
          integer(kind = 8) :: LocalParticleIndex
        
          LocalNewElementID = NewElementID
          LocalParticleIndex = ParticleIndex

          ElementIDArray(ParticleIndex)= NewElementID
          !PatchIDArray(ParticleIndex)= NewPatchID
          EleParticles(IParticle) = LocalNewElementID * PartFactor +  &
                                    LocalParticleIndex
        end subroutine SetParticleElementID


        subroutine InitialiseHouseKeepingArrays()
        !**********************************************************************
        !
        !    Function:  Allocates the arrays storing house-keeping data.
        !               Sets their initial values.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
          
          ! local variables
          integer(INTEGER_TYPE) :: IError
          integer(INTEGER_TYPE) :: I
          
          ! Multipatch variable
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1

          !allocate(EleParticlesHelp(Counters%NEl), stat = IError) ! I need to make this patch dependent         
          allocate(EleParticlesHelp(Maximum_nel_NURBS,NumberOfPatches), stat = IError) ! I need to make this patch dependent
          EleParticlesHelp = -1
          
          allocate(EleParticles(Counters%NParticles), stat = IError)
           
          EleParticles = -1
          
          allocate(ParticleType(NTENSOR, NVECTOR, MAX_LOAD_SYSTEMS)::Particles(Counters%NParticles), stat = IError)
          do i = 1, Counters%NParticles 
              call SetParticleStructureDefault(Particles(i))
          end do
          
          allocate(ElementIDArray(Counters%NParticles), stat = IError)
          
          ! Multipatch generalization 
          allocate(PatchIDArray(Counters%NParticles), stat = IError)
          
          allocate(IDArray(Counters%NParticles), stat = IError)
          allocate(EntityIDArray(Counters%NParticles), stat = IError)
          allocate(MaterialIDArray(Counters%NParticles), stat = IError)
          allocate(UArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(MaterialPointTypeArray(Counters%NParticles), stat = IError)
          allocate(UStepArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(ShapeValuesArray(Counters%NParticles,ELEMENTNODES), stat = IError)
          
          
          ! For NURBS debugging can create these functions here 
          ! ------------------------------------------------------------------------------------

          ! Xi
          allocate(ShapeValuesArray_Xi(NXiKnotOrder(IPatch_Temporary)+1, Counters%NParticles), stat = IError) ! this needs to be a maximum
          !Counters%NParticles,ELEMENTNODES), stat = IError)
          allocate(DShapeValuesArray_Xi(NXiKnotOrder(IPatch_Temporary)+1, Counters%NParticles), stat = IError) ! this needs to be a maximum
          !Counters%NParticles,ELEMENTNODES), stat = IError)
          
          
          ! Eta
          allocate(ShapeValuesArray_Eta(NEtaKnotOrder(IPatch_Temporary)+1, Counters%NParticles), stat = IError) ! this needs to be a maximum
          !Counters%NParticles,ELEMENTNODES), stat = IError)
          allocate(DShapeValuesArray_Eta(NEtaKnotOrder(IPatch_Temporary)+1, Counters%NParticles), stat = IError) ! this needs to be a maximum
          !Counters%NParticles,ELEMENTNODES), stat = IError)
          
          
          ! Weights 
          allocate(WeightArray_Xi(Counters%NParticles), stat = IError) !(NXiKnotOrder+1, Counters%NParticles), stat = IError)
          !Counters%NParticles,ELEMENTNODES), stat = IError)
          allocate(WeightArray_Eta(Counters%NParticles), stat = IError) !(NEtaKnotOrder+1, Counters%NParticles), stat = IError)
          !Counters%NParticles,ELEMENTNODES), stat = IError)
          
          
          ! ------------------------------------------------------------------------------------
          

          !
          !allocate(NN_IncludesZeroValues_Print(NXiKnotEntries-1), stat = IError)
          !allocate(dN_dxi_IncludesZeroValues_Print(NXiKnotEntries-1), stat = IError)
          !
          !
          !allocate(MM_IncludesZeroValues_Print(NEtaKnotEntries-1), stat = IError)
          !allocate(dM_deta_IncludesZeroValues_Print(NEtaKnotEntries-1), stat = IError)
          
          
          
          ! ------------------------------------------------------------------------------------

          
          
          allocate(UPhaseArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(GlobPosArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(VelocityArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(SigmaEffArray(Counters%NParticles,NTENSOR), stat = IError)
          allocate(SigmaEff0Array(Counters%NParticles,NTENSOR), stat = IError)
          allocate(DShapeValuesArray(Counters%NParticles,ELEMENTNODES,NVECTOR), stat = IError)
          allocate(VelocityWaterArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(VelocityGasArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(AccelerationArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(MassArray(Counters%NParticles), stat = IError)
          allocate(ESMstatevArray(Counters%NParticles, NSTATEVAR), stat = IError)
          allocate(MassWaterArray(Counters%NParticles), stat = IError)
          
          !allocate(IsParticleIntegration(Counters%NEl), stat = IError)
          
          !allocate(IsParticleIntegration(Counters%NEl), stat = IError)         
          allocate(IsParticleIntegration(Maximum_nel_NURBS, NumberOfPatches), stat = IError) 

          ! above this should the number of control point elements which is only three 
          
          ! Initialize
          ElementIDArray = -1
          
          ! Multipatch generalization 
          PatchIDArray = -1 
          
          IDArray = - 1
          EntityIDArray = -1
          MaterialIDArray = -1
          MaterialPointTypeArray = MaterialPointTypeUndefined
          UStepArray = 0.0
          UArray = 0.0
          ShapeValuesArray = 0.0
          UPhaseArray = 0.0
          GlobPosArray =  0.0
          VelocityArray = 0.0
          SigmaEffArray = 0.0
          SigmaEff0Array = 0.0
          DShapeValuesArray = 0.0
          VelocityWaterArray = 0.0
          VelocityGasArray = 0.0
          AccelerationArray = 0.0
          MassArray = -1.0
          MassWaterArray = -1.0
          
          if (IsMPMWithMPIntegration()) then
            ! Initially (and always) all elements are partially filled, always use mass point integration
            IsParticleIntegration = .true.
          else
            ! At least initially all elements are fully filled
            IsParticleIntegration = .false.
          end if

        end subroutine InitialiseHouseKeepingArrays


        subroutine DestroyHouseKeeping()
        !**********************************************************************
        !
        !    Function:  Frees the arrays belonging to this module.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError
                
          if (allocated(NPartEle) ) then
            deallocate(NPartEle, stat = IError)
          end if
        
          if (allocated(EleParticlesHelp) ) then
            deallocate(EleParticlesHelp, stat = IError)
          end if
          
          if (allocated(EleParticles) ) then
            deallocate(EleParticles, stat = IError)
          end if
          
          if (allocated(Particles) ) then
            deallocate(Particles, stat = IError)
          end if
         
         !Fileds stripped from Particle datatype
         if (allocated(IDArray) ) then
            deallocate(IDArray, stat = IError)
         end if
                    
          if (allocated(ElementIDArray) ) then
            deallocate(ElementIDArray, stat = IError)
          end if
          
          if (allocated(PatchIDArray) ) then
            deallocate(PatchIDArray, stat = IError)
          end if

         if (allocated(EntityIDArray) ) then
            deallocate(EntityIDArray, stat = IError)
         end if
                    
         if (allocated(MaterialIDArray) ) then
            deallocate(MaterialIDArray, stat = IError)
         end if  
         
         if (allocated(UArray) ) then
            deallocate(UArray, stat = IError)
         end if  
         
        if (allocated(MaterialPointTypeArray) ) then
            deallocate(MaterialPointTypeArray, stat = IError)
        end if
        
        if (allocated(UStepArray) ) then
            deallocate(UStepArray, stat = IError)
        end if 
        
        if (allocated(ShapeValuesArray) ) then
            deallocate(ShapeValuesArray, stat = IError)
        end if
        
         if (allocated(UPhaseArray) ) then
            deallocate(UPhaseArray, stat = IError)
         end if
                 
         if (allocated(GlobPosArray) ) then
            deallocate(GlobPosArray, stat = IError)
         end if
         
         if (allocated(VelocityArray) ) then
            deallocate(VelocityArray, stat = IError)
         end if
         
         if (allocated(SigmaEffArray) ) then
            deallocate(SigmaEffArray, stat = IError)
         end if
         
         if (allocated(SigmaEff0Array) ) then
            deallocate(SigmaEff0Array, stat = IError)
         end if
         
         if (allocated(DShapeValuesArray) ) then
            deallocate(DShapeValuesArray, stat = IError)
         end if
         
         if (allocated(VelocityWaterArray) ) then
            deallocate(VelocityWaterArray, stat = IError)
         end if
         
         if (allocated(VelocityGasArray) ) then
            deallocate(VelocityGasArray, stat = IError)
         end if
         
         if (allocated(AccelerationArray) ) then
            deallocate(AccelerationArray, stat = IError)
         end if
         
         if (allocated(MassArray)) then
            deallocate(MassArray, stat = IError)
         end if
         
         if (allocated(ESMstatevArray)) then
            deallocate(ESMstatevArray, stat = IError)
         end if
         
         if (allocated(MassWaterArray)) then
            deallocate(MassWaterArray, stat = IError)
         end if
         
          if (allocated(IsParticleIntegration) ) then
            deallocate(IsParticleIntegration, stat = IError)
          end if
          
        end subroutine DestroyHouseKeeping


        subroutine DetermineDecimalFactor()
        !**********************************************************************
        !
        !    Function:  Determines the decimal factor PartFactor which is used
        !               to set up the particle house-keeping arrays. Requires
        !               Counters%NParticles. Writes PartFactor.
        !               Determines the next highest decimal number.
        !               For example: Value = 345 return 1000 etc
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          PartFactor = DecimalFactor(Counters%NParticles)
        
        end subroutine DetermineDecimalFactor


        integer(kind = 8) function DecimalFactor(Number)
        !**********************************************************************
        !
        !    Function:  Determines the next highest decimal number of Number
        !               For example: Value = 345 return 1000 etc
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: Number

          DecimalFactor = 10
          do while (DecimalFactor<=Number)
            DecimalFactor = 10 * DecimalFactor
          end do
        
        end function DecimalFactor


        subroutine CoordLocalToGlobal(IElement, LNodalCoordinates, IPatch)
        !**********************************************************************
        !
        !    Function:  Determines the global coordinates associated with the local
        !               coordinates of particles inside a considered element.
        !               Updates the global coordinates of particles inside IElement.
        !
        !     IElement : ID of the considered element
        !     LNodalCoordinates : Nodal coordinates used for updating the integration point locations
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        use ModCounters

        implicit none 

          integer(INTEGER_TYPE), intent(in) :: IElement
          !real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR) :: LNodalCoordinates
          real(REAL_TYPE), dimension(NControlPoints(IPatch), NVECTOR) :: LNodalCoordinates

          ! Local variables
          integer(INTEGER_TYPE) :: NodeID, ParticleIndex
          integer(INTEGER_TYPE) :: I, IParticle, IDim
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
          do IParticle = 1, NPartEle(IElement,IPatch) ! Loop over particles of IElement
            ParticleIndex = GetParticleIndexFunction(IParticle, IElement, IPatch)

            GlobPosArray(ParticleIndex,:) = 0.0 ! Reset coordinates to zero
            do I = 1, ELEMENTNODES ! Loop over nodes of IElement
              !NodeID = iabs(Multipatch_Connecting_Local_To_Global_ControlPoints(ElementConnectivities(I, IElement, IPatch), IPatch ) )
              NodeID = iabs(ElementConnectivities(I, IElement, IPatch) )

              do IDim = 1, NVECTOR ! Loop over dimensions of IElement
                GlobPosArray(ParticleIndex,IDim) = GlobPosArray(ParticleIndex,IDim) + LNodalCoordinates(NodeID, IDim) * ShapeValuesArray(ParticleIndex,I)  ! hardcoded
              end do
            end do
            
          end do

        end subroutine CoordLocalToGlobal

        subroutine SetParticleShapeFunctionData(Particle, ParticleIndex, ElementID, ICon, Co, PatchID)
        !**********************************************************************
        !
        !     Function:  Calculate the shape function values and derivatives
        !               for the new local position of Particle.
        !
        !     Particle : Considered particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        
        implicit none

          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)) :: Particle
          integer(INTEGER_TYPE), intent(in) :: ParticleIndex, ElementID
          ! Local variables
          real(REAL_TYPE), allocatable, dimension(:) :: ShapeValues
          real(REAL_TYPE), allocatable, dimension(:, :) :: DShapeValues
          
          real(REAL_TYPE) :: Xi_ParametricDomain, Eta_ParametricDomain, Zeta_ParametricDomain
          
          ! NURBS Xi output for debugging purposes 
          !real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Xi 
          !real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Xi
          !                       
          !! NURBS Eta output for debugging purposes 
          !real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Eta 
          !real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Eta
          !
          !real(REAL_TYPE), dimension(:), allocatable :: ShapeValuesArray_Temp_Zeta
          !real(REAL_TYPE), dimension(:), allocatable :: DShapeValuesArray_Temp_Zeta
          
          ! Jacobian variables 
          
          !integer(INTEGER_TYPE), intent(in) :: NEl, NodTot, IDim, IElement
          !real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(:, :), intent(in) :: Co
          real(REAL_TYPE), dimension(NDIM, NDIM) :: RJac
          real(REAL_TYPE), dimension(NDIM, NDIM) :: dx_dxi
          real(REAL_TYPE), dimension(NDIM, NDIM) :: dxi_dx
          real(REAL_TYPE), dimension(NDIM, NDIM) :: InvRJac
          real(REAL_TYPE) :: DetJac
          
          !real(REAL_TYPE), dimension(nen_NURBS, NDIM) :: dR_dx
          
          integer(INTEGER_TYPE) :: I, J, INode, NodeID, NNodes
          !real(REAL_TYPE), dimension(:), allocatable :: HS
          !real(REAL_TYPE), dimension(:, :), allocatable :: dHS
          real(REAL_TYPE) :: Det1
          
          ! Hardcoded index to 1 
          real(REAL_TYPE), dimension(nen_NURBS(1),NDIM) :: dR_dx ! IPatch_Temporary
          real(REAL_TYPE), dimension(NDIM,NDIM) :: dxi_dtildexi
          real(REAL_TYPE), dimension(NDIM,NDIM) :: J_mat
          integer(INTEGER_TYPE) :: aa, bb, cc, loc_num, ii, jj
          integer(INTEGER_TYPE) :: ni, nj, nk
                    
          real(REAL_TYPE), dimension(NDIM,NDIM) :: dx_dtildexi
          
          real(REAL_TYPE) :: Det_dx_dxi
          real(REAL_TYPE) :: Det_dxi_dx
          
          ! Multipatch variable 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: PatchID
          !integer(INTEGER_TYPE) :: aa
          
          !ShapeValues_Xi, DShapeValues_Xi, Weights_Xi, &
          !
          !    ShapeValues_Eta, DShapeValues_Eta, Weights_Eta
          
          
          
          !! NURBS
          !real(REAL_TYPE), allocatable, dimension(:,:) :: HS
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS 
          !real(REAL_TYPE), allocatable, dimension(:) :: Wt !classic inout parameters
          !                                  
          !real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Xi 
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Xi
          !real(REAL_TYPE), allocatable, dimension(:) :: Wt_Xi
          !
          !real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Eta 
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Eta
          !real(REAL_TYPE), allocatable, dimension(:) :: Wt_Eta          
          !
          !
          !allocate( Wt(ELEMENTGAUSSPOINTS) )
          !allocate( HS(ELEMENTGAUSSPOINTS, ELEMENTNODES) )
          !allocate( dHS(ELEMENTGAUSSPOINTS, ELEMENTNODES, NDOFL) ) 
          !
          !! Xi direction 
          !allocate( Wt_Xi(ELEMENTBOUNDARYGAUSSPOINTS) )
          !allocate( HS_Xi(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES) )
          !allocate( dHS_Xi(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NDOFL-1) )
          !
          !! Eta direction 
          !allocate( Wt_Eta(ELEMENTBOUNDARYGAUSSPOINTS) )
          !allocate( HS_Eta(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES) )
          !allocate( dHS_Eta(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NDOFL-1) )
        
          
          !Xi_ParametricDomain = Particle%LocPos(1)
          !Eta_ParametricDomain = Particle%LocPos(2)
          !
          !call InitialiseShapeFunctionsQUAD4_NURBS(HS, dHS, Wt, & !classic inout parameters
          !                                          HS_Xi, dHS_Xi, Wt_Xi, &
          !                                          HS_Eta, dHS_Eta, Wt_Eta, &
          !                                          XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
          !                                          EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
          !                                          ni, nj) !NURBS related inputs in the eta direction 
          

          call ShapeFunctionData(Particle%LocPos, ELEMENTNODES, ShapeValues, DShapeValues, ElementID, PatchID) ! ELEMENTNODES,
          
          
          
          
          
          ! finding the NURBS coordinates for the element

          ni = INN(IEN(1,ElementID,PatchID),1,PatchID)    
          nj = INN(IEN(1,ElementID,PatchID),2,PatchID)
          
          	
          if (NDIM == 3) then 
              nk = INN(IEN(1,ElementID,PatchID),3,PatchID)
          end if 
          
          
          
          
          
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating dx/dxi
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          
          call Evaluate_dx_dxi(DShapeValues, Co, ICon, dx_dxi, ElementID)
          !NNodes = nen_NURBS
          !! Determine the Jacobian matrix RJac
          !RJac = 0.0
          !dx_dxi = 0
          !do INode = 1, NNodes ! loop nodes of each element
          !  
          !    !ni = INN(IEN(1,IElement),1)    
          !    !nj = INN(IEN(2,IElement),2)
          !    
          !    !XiKnotEntries(ni) 
          !    !note that the coordinates of the control points are the ones we need here in 'Co'
          !    NodeID = ICon(INode, ElementID)
          !    ! This is gradient mapping from the parameter space to the physical space 
          !  do I = 1, NDIM
          !    do J = 1, NDIM
          !      dx_dxi(I, J) = dx_dxi(I, J) + DShapeValues(INode, I) * Co(NodeID, J) !calculation of the jacobian matrix from the coordinates of one element ! -> this is checked and now it is working
          !    end do
          !  end do
          !end do
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating dxi/dx
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! debugger
      if (dx_dxi(3,3) == 0.0) then 
          aa=1
      end if

          
          
          
          ! Calculate inverse and determinant of Jacobian matrix
          call RJacInv(NDIM, dx_dxi, dxi_dx, Det_dx_dxi, Det_dxi_dx)

          
          
          ! InvRJac -> dxi_dx 
          !
          !
          ! EVERYTHING BELOW HERE SHOULD GO IN THE PARENT SUBROUTINE 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating dR/dx
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !Calculate derivative of basis functions with respect to physical coordinates
          call Get_dR_dx_MappingShapeFunctionDerivativeWRTPhysicalCoordinates(dxi_dx, DShapeValues, dR_dx, PatchID) ! dR_dx = dxi_dx * dR_dxi 
          
          
          
          ! Note this Jacobian needs to be updated but since we are not outputting it then there is no need to add more calcs
          if ( (Det_dx_dxi < 0.0) .or. (Det_dxi_dx < 0.0) ) then
            call WriteInLogFile('Negative determinant '// trim(String(DetJac)) // trim(String(Det1)) // ' element ' // trim(String(ElementID)))
          end if

          
          
          
          
          
          
          
          
          
          
          
          
          
          
          ShapeValuesArray(ParticleIndex, :) = ShapeValues
          DShapeValuesArray(ParticleIndex, :, : ) = dR_dx

        end subroutine SetParticleShapeFunctionData
        
        
        
        
        
        !----------------------------------------------------------------------------------
        ! SINGLE PARTICLE 
        
        
        
        !subroutine SetParticleShapeFunctionData_SingleParticle(Particle, ParticleIndex, ElementID)
        !!**********************************************************************
        !!
        !!     Function:  Calculate the shape function values and derivatives
        !!               for the new local position of Particle.
        !!
        !!     Particle : Considered particle
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)) :: Particle
        !  integer(INTEGER_TYPE), intent(in) :: ParticleIndex, ElementID
        !  ! Local variables
        !  real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues
        !  real(REAL_TYPE), dimension(ELEMENTNODES, NVECTOR) :: DShapeValues
        !  
        !  real(REAL_TYPE) :: Xi_ParametricDomain, Eta_ParametricDomain
        !  
        !  ! NURBS Xi output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Xi 
        !  real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Xi
        !                         
        !  ! NURBS Eta output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Eta 
        !  real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Eta
        !  
        !  real(REAL_TYPE), dimension(:), allocatable :: ShapeValuesArray_Temp_Zeta
        !  real(REAL_TYPE), dimension(:), allocatable :: DShapeValuesArray_Temp_Zeta
        !  
        !  !ShapeValues_Xi, DShapeValues_Xi, Weights_Xi, &
        !  !
        !  !    ShapeValues_Eta, DShapeValues_Eta, Weights_Eta
        !  
        !  
        !  
        !  !! NURBS
        !  !real(REAL_TYPE), allocatable, dimension(:,:) :: HS
        !  !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS 
        !  !real(REAL_TYPE), allocatable, dimension(:) :: Wt !classic inout parameters
        !  !                                  
        !  !real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Xi 
        !  !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Xi
        !  !real(REAL_TYPE), allocatable, dimension(:) :: Wt_Xi
        !  !
        !  !real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Eta 
        !  !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Eta
        !  !real(REAL_TYPE), allocatable, dimension(:) :: Wt_Eta          
        !  !
        !  !
        !  !allocate( Wt(ELEMENTGAUSSPOINTS) )
        !  !allocate( HS(ELEMENTGAUSSPOINTS, ELEMENTNODES) )
        !  !allocate( dHS(ELEMENTGAUSSPOINTS, ELEMENTNODES, NDOFL) ) 
        !  !
        !  !! Xi direction 
        !  !allocate( Wt_Xi(ELEMENTBOUNDARYGAUSSPOINTS) )
        !  !allocate( HS_Xi(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES) )
        !  !allocate( dHS_Xi(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NDOFL-1) )
        !  !
        !  !! Eta direction 
        !  !allocate( Wt_Eta(ELEMENTBOUNDARYGAUSSPOINTS) )
        !  !allocate( HS_Eta(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES) )
        !  !allocate( dHS_Eta(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NDOFL-1) )
        !
        !  
        !  !Xi_ParametricDomain = Particle%LocPos(1)
        !  !Eta_ParametricDomain = Particle%LocPos(2)
        !  !
        !  !call InitialiseShapeFunctionsQUAD4_NURBS(HS, dHS, Wt, & !classic inout parameters
        !  !                                          HS_Xi, dHS_Xi, Wt_Xi, &
        !  !                                          HS_Eta, dHS_Eta, Wt_Eta, &
        !  !                                          XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
        !  !                                          EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
        !  !                                          ni, nj) !NURBS related inputs in the eta direction 
        !  
        !
        !  call ShapeFunctionData_SingleParticle(Particle%LocPos, ELEMENTNODES, ShapeValues, DShapeValues, ElementID, &
        !                         ShapeValuesArray_Temp_Xi, DShapeValuesArray_Temp_Xi, &
        !                         ShapeValuesArray_Temp_Eta, DShapeValuesArray_Temp_Eta, &
        !                         ShapeValuesArray_Temp_Zeta, DShapeValuesArray_Temp_Zeta)
        !  
        !  ! This is where we should save the values of ShapeValuesArray_Xi, DShapeValuesArray_Xi, WeightsArray_Xi,  ShapeValuesArray_Eta, DShapeValuesArray_Eta, WeightsArray_Eta
        !  ! for debugging purposes 
        !  
        !  !ShapeValuesArray_Xi(:,  IElement) = HS_Xi
        !  !DShapeValuesArray_Xi(:, IElement) = dHS_Xi
        !  !WeightArray_Xi(:,  IElement) = Wt_Xi
        !  
        !  !ShapeValues_Xi, DShapeValues_Xi, Weights_Xi, &
        !  !ShapeValues_Eta, DShapeValues_Eta, Weights_Eta
        !  
        !  ! Xi 
        !  ShapeValuesArray_Xi(:,  ParticleIndex) = ShapeValuesArray_Temp_Xi
        !  DShapeValuesArray_Xi(:, ParticleIndex) = DShapeValuesArray_Temp_Xi
        !  !WeightArray_Xi(:,  ParticleIndex) = Weights_Xi    
        !  
        !  
        !  ! Eta
        !  ShapeValuesArray_Eta(:,  ParticleIndex) = ShapeValuesArray_Temp_Eta
        !  DShapeValuesArray_Eta(:, ParticleIndex) = DShapeValuesArray_Temp_Eta
        !  !WeightArray_Eta(:,  ParticleIndex) = Weights_Eta
        !  
        !  
        !  ! Zeta
        !  ! ... need to extend to 3D
        !  
        !  ShapeValuesArray(ParticleIndex, : ) = ShapeValues ! no of rows = no of particles, no of columns = no of nodes per element
        !
        !  DShapeValuesArray(ParticleIndex, :, : ) = DShapeValues
        !
        !end subroutine SetParticleShapeFunctionData_SingleParticle
        
        
        
        ! SINGLE PARTICLE
        !----------------------------------------------------------------------------------------------
        
        
        
        
        
        
        
        
        
        
        
        
        
        


        subroutine DetermineAdjacentParticles(ISide, NParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles from 1 to NParticles
        !               of the element lie next to ISide and how many.
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage! (Which particles of
        !                     an element lie adjacent to any side of the element? )
        !
        !     ISide : Local number of considered element side
        !     NParticles : Number of particles inside the considered element
        !
        ! O   ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ISide, NParticles
          logical, dimension(NParticles), intent(inout) :: ParticleStatus

            select case(ELEMENTTYPE)
                case(TRI3)
                    select case (NParticles)  
                        case(1) ! 1 material point per elmeent / FEM calculation
                          call DetermineAdjacentParticlesTRI_MP1(NParticles, ParticleStatus)
                        case(3) ! 3 material points per element
                          call DetermineAdjacentParticlesTRI_MP3(ISide, NParticles, ParticleStatus)
                        case(6) ! 6 material points per element
                            call DetermineAdjacentParticlesTRI_MP6(ISide, NParticles, ParticleStatus)
                        case(12) ! 12 material points per element
                            call DetermineAdjacentParticlesTRI_MP12(ISide, NParticles, ParticleStatus)
                        case(16) ! 16 material points per element
                            call DetermineAdjacentParticlesTRI_MP16(ISide, NParticles, ParticleStatus)
                        case(25) ! 25 material points per element
                            call DetermineAdjacentParticlesTRI_MP25(ISide, NParticles, ParticleStatus)
                        case(46) ! 46 material points per element
                            call DetermineAdjacentParticlesTRI_MP46(ISide, NParticles, ParticleStatus)
                        case(88) ! 88 material points per element
                            call DetermineAdjacentParticlesTRI_MP88(ISide, NParticles, ParticleStatus)
                        case default
                          call GiveError("Undefined number of material points per element in [subroutine DetermineAdjacentParticles()], 2D case.")
                    end select
                case(TETRAOLD)
                    select case (NParticles)
                        case(1) ! 1 material point per element / FEM calculation
                            call DetermineAdjacentParticlesTETRA_MP1(NParticles, ParticleStatus)
                        case(4) ! 4 material points per element
                            call DetermineAdjacentParticlesTETRA_MP4(ISide, NParticles, ParticleStatus)
                        case(7) ! 7 material points per element
                            call DetermineAdjacentParticlesTETRA_MP7(ISide, NParticles, ParticleStatus)
                        case(8) ! 8 material points per element
                            call DetermineAdjacentParticlesTETRA_MP8(ISide, NParticles, ParticleStatus)
                        case(10) ! 10 material points per element
                            call DetermineAdjacentParticlesTETRA_MP10(ISide, NParticles, ParticleStatus)
                        case(13) ! 13 material points per element
                            call DetermineAdjacentParticlesTETRA_MP13(ISide, NParticles, ParticleStatus)
                        case(20) ! 20 material points per element
                            call DetermineAdjacentParticlesTETRA_MP20(ISide, NParticles, ParticleStatus)
                        case default
                          call GiveError("Undefined number of material points per element in [subroutine DetermineAdjacentParticles()], 3D case.")
                        end select
                        
                        
                        
                        
                        
                        
                        
                        
                case(QUAD4_NURBS)
                    select case (NParticles)
                        case(1) ! 1 material point per element / FEM calculation
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP1(NParticles, ParticleStatus)
                        case(4) ! 4 material points per element
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP4(NParticles, ParticleStatus)
                        case(9) ! 9 material points per element 
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP9(NParticles, ParticleStatus)
                        case(16) ! 1 material point per element / FEM calculation
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP16(NParticles, ParticleStatus)
                        case(25) ! 4 material points per element
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP25(NParticles, ParticleStatus)
                        case(36) ! 9 material points per element 
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP36(NParticles, ParticleStatus)
                        case(49) ! 1 material point per element / FEM calculation
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP49(NParticles, ParticleStatus)
                        case(64) ! 4 material points per element
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP64(NParticles, ParticleStatus)
                        case(81) ! 9 material points per element 
                            call DetermineAdjacentParticlesQUAD4_NURBS_MP81(NParticles, ParticleStatus)
                            
                            !call DetermineAdjacentParticlesQUAD4_NURBS_MP4(ISide, NParticles, ParticleStatus)
                        !case(7) ! 7 material points per element
                        !    call DetermineAdjacentParticlesQUAD4_NURBS_MP7(ISide, NParticles, ParticleStatus)
                        !case(8) ! 8 material points per element
                        !    call DetermineAdjacentParticlesQUAD4_NURBS_MP8(ISide, NParticles, ParticleStatus)
                        !case(10) ! 10 material points per element
                        !    call DetermineAdjacentParticlesQUAD4_NURBS_MP10(ISide, NParticles, ParticleStatus)
                        !case(13) ! 13 material points per element
                        !    call DetermineAdjacentParticlesQUAD4_NURBS_MP13(ISide, NParticles, ParticleStatus)
                        !case(20) ! 20 material points per element
                        !    call DetermineAdjacentParticlesQUAD4_NURBS_MP20(ISide, NParticles, ParticleStatus)
                        !case default
                        !  call GiveError("Undefined number of material points per element in [subroutine DetermineAdjacentParticles()], 3D case.")
                        end select
                
                        
                        
                        
                        
                end select

        end subroutine DetermineAdjacentParticles

        
        integer(INTEGER_TYPE) function DetermineNAdjacentParticles(ParticleStatus, NParticles)
        !**********************************************************************
        !
        !    Function:  Determines how many of the particles from 1 to Flags%NMaterialPoints
        !               of the element lie next to a side by evaluating ParticleStatus.
        !
        !     ParticleStatus : Array storing for each particle whether it was found
        !                      to lie next to a checked element side
        !
        ! O   DetermineNAdjacentParticles : Number of particles adjacent to side(s) of 
        !                                   the element
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
         integer(INTEGER_TYPE), intent (in):: NParticles
          logical, dimension(NParticles), intent(in) :: ParticleStatus
          ! Local variables
          integer(INTEGER_TYPE) :: I

          DetermineNAdjacentParticles = 0
          
          do I = 1, NParticles
            if (ParticleStatus(I) ) then
              DetermineNAdjacentParticles = DetermineNAdjacentParticles + 1
            end if
          end do

        end function DetermineNAdjacentParticles

        subroutine DefineVariableData()
        !**********************************************************************
        !
        !    Function:  Sets Particles(I)%VariableData(1 .. 3) fields. This fields can be
        !               used for debugging purposes to display some particle data
        !               with the Output program.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, IElement
        
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1

          
          do IPatch = 1, NumberOfPatches ! loop over patches
          do IElement = 1, nel_NURBS(IPatch)!Counters%NEl ! loop over elements
        
            do IParticle = 1, NPartEle(IElement,IPatch)

              ParticleIndex = GetParticleIndex(IParticle, IElement, IPatch)

              ! Indicate partially and fully filled elements -> VariableData(1)
              if (IsParticleIntegration(IElement, IPatch) ) then
                call SetVariableDataI(Particles(ParticleIndex), 1, 1.d0)
              else
                call SetVariableDataI(Particles(ParticleIndex), 1, 0.d0)
              end if

              call SetVariableDataI(Particles(ParticleIndex), 2, Particles(ParticleIndex)%IntegrationWeight)

              call SetVariableDataI(Particles(ParticleIndex), 3, (Particles(ParticleIndex)%Density))

            end do
          
          end do ! loop over elements
          end do ! loop over patches
        
        end subroutine DefineVariableData        

      end module ModMPMData
