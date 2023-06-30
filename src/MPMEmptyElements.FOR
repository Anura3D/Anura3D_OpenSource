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
	  
	  
	  module ModEmptyElements
      !**********************************************************************
      !
      !    Function:  Contains routines for detecting gaps inside a
      !               body assumed as fully filled with material.
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 9801 $
      !     $Date: 2022-10-11 15:00:20 +0200 (di, 11 okt 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModMPMData
      use ModMeshAdjacencies
      use ModWriteTestData
      use ModMPMMeshAdjustment
      use ModGlobalConstants
      use ModFileIO

      
      implicit none

        logical, dimension(:), allocatable :: IsDanglingElement ! True, if IElement is a dangling element, size = NEl
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: DanglingBodySizes
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: CheckedElements ! Array containing the group ID's of the found element groups
        integer(INTEGER_TYPE), save, private :: CheckCounter ! Number of found elements connected to IElement
        integer(INTEGER_TYPE), save, private :: NewGroupID ! -1 if no connected element group is found, else the ID of the found connected group
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: GroupStorage ! Temporarily store element IDs of current group
        integer(INTEGER_TYPE) :: NSkipEleActivation ! Number of elements not to be considered part of a hole
        integer(INTEGER_TYPE), dimension(:), allocatable :: SkipEleActivation ! IDs of elements not to be considered part of a hole
        integer(INTEGER_TYPE), dimension(:), allocatable :: Holes ! Array containing for each element the group ID's of the found holes
        integer(INTEGER_TYPE), dimension(:), allocatable :: HoleSizes ! Array containing for each group the number of elements
        integer(INTEGER_TYPE), save, private :: MaxHoleID, MaxDanglingBodyID, ActivationStatus ! Maximum ID of hole found for certain load step
        logical, dimension(:), allocatable :: IsEmptyElement ! True, if IElement is an element with no particles inside, which should be considered fully filled
        integer(INTEGER_TYPE), save :: NEmptyElements ! Number of empty elements

      contains ! Routines of this module

        subroutine InitialiseModEmptyElementsData()
        !**********************************************************************
        !
        !    Function:  Initialisation of the data of this modules.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          allocate(IsEmptyElement(Counters%NEl), stat = IError)
          IsEmptyElement = .false. ! Initially all fully filled elements should be considered filled
          NEmptyElements = 0
           
          allocate(IsDanglingElement(Counters%NEl), stat = IError)
          IsDanglingElement = .false.

          NSkipEleActivation = 0
          NewGroupID = 0
          CheckCounter = 0
                
        end subroutine InitialiseModEmptyElementsData

        subroutine DestroyModEmptyElementsArrays()
        !**********************************************************************
        !
        !    Function:  Deallocates globally accessed allocatable arrays of this module.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(SkipEleActivation) ) then
            deallocate(SkipEleActivation, stat = IError)
          end if

          if (allocated(IsEmptyElement) ) then
            deallocate(IsEmptyElement, stat = IError)
          end if

          if (allocated(IsDanglingElement) ) then
            deallocate(IsDanglingElement, stat = IError)
          end if

        end subroutine DestroyModEmptyElementsArrays


        logical function SkipElementActivation(IElement)
        !**********************************************************************
        !
        !    Function:  Returns .true. if the element is specified as an element
        !               not to be activated if empty.
        !
        !     IElement : ID of considered element
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          !arguments
          integer(INTEGER_TYPE), intent(in) :: IElement
          
          ! Local variables
          integer(INTEGER_TYPE) :: I
          
          SkipElementActivation = .false.
          
          if (ContactSurfaceSoilElements(IElement)) then
            SkipElementActivation = ContactSurfaceSoilElements(IElement)
          else
            do I = 1, NSkipEleActivation
              if (SkipEleActivation(I)==IElement) then
                SkipElementActivation = .true.
                EXIT
              end if
            end do
          end if

        end function SkipElementActivation


        subroutine CheckEmptyElements()
        !**********************************************************************
        !
        !    Function:  Checks the newly determined discretisation of activated/
        !               deactivated elements.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IError

          allocate(Holes(Counters%NEl), stat = IError)
          Holes = 0
          
          allocate(CheckedElements(Counters%NEl), stat = IError)
          CheckedElements = 0 ! No elements checked, no groups found
          allocate(GroupStorage(CalParams%RecursionThreshold), stat = IError)
          GroupStorage = 0 

          MaxHoleID = 0 ! Initial ID of considered group of elements

          ! Search groups of connected deactivated elements in the whole mesh
          ActivationStatus = 0
          call CheckElementGroups(MaxHoleID)

          ! Sort by size groups of deactivated elements
          allocate(HoleSizes(MaxHoleID), stat = IError)
          HoleSizes = 0
          
          call SortElementGroups(MaxHoleID, HoleSizes)
          
          ! Check smaller groups for elements forming a hole
          if (CalParams%ApplyEmptyElements) then
            call CheckHoles()
            call OutputMaximumHoleSizes()      
            call AssignMaterialToHoles()
          end if
          
          ! ... searching for empty elements could also be put inside a loop if not all are found ...
          
          ! Search groups of connected activated elements
          MaxDanglingBodyID = 0
          
          ActivationStatus = 1
          call CheckElementGroups(MaxDanglingBodyID)
          
          ! Sort by size groups of activated elements
          allocate(DanglingBodySizes(MaxDanglingBodyID), stat = IError)
          DanglingBodySizes = 0
          
          call SortElementGroups(MaxDanglingBodyID, DanglingBodySizes)
          
          ! Check smaller groups for elements forming dangling elements
          call CheckDanglingElements()
          
          ! All elements should have been checked, else something went wrong
          call Assert(IsAllElementsChecked(), 'Something went wrong in CheckEmptyElements!, All elements have not been checked!')

          deallocate(DanglingBodySizes, stat = IError)
          deallocate(CheckedElements, stat = IError)
          deallocate(GroupStorage, stat = IError)
          deallocate(HoleSizes, stat = IError)
          deallocate(Holes, stat = IError)

          MaxHoleID = 0
          MaxDanglingBodyID = 0


        end subroutine CheckEmptyElements
 
        subroutine CheckElementGroups(MaxGroupID)
        !**********************************************************************
        !
        !    Function:  Checks the newly determined discretisation of activated and
        !               deactivated elements.
        !               Elements forming holes will be marked as such.
        !
        ! O   MaxGroupID : ID of the currently considered group
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(inout) :: MaxGroupID
          ! Local variables
          integer(INTEGER_TYPE) :: IElement

          do IElement = 1, Counters%NEl
            if ((IsActiveElement(IElement).and.(ActivationStatus==1)).or. &
                (.not.IsActiveElement(IElement).and.(ActivationStatus==0))) then ! Loop over de-/activated elements
              if (CheckedElements(IElement)==0) then ! The element has not yet been checked

                CheckCounter = 0 ! No connected elements found yet for element IElement
                MaxGroupID = MaxGroupID + 1 ! ID of the currently considered group of elements
                NewGroupID = -1 ! ID of found group connected to the considered group of elements
                GroupStorage = 0 ! Reset temporary storage of group element ID's

                ! Recursive check for elements connected to IElement
                call CheckSideElement(IElement, MaxGroupID)
                
                ! Check whether groups should be merged (if NewGroupID is not -1)
                call CheckGroupMerging(MaxGroupID)
                
              end if
            end if
          end do
        
        end subroutine CheckElementGroups

        recursive subroutine CheckSideElement(IElement, GroupID)
        !**********************************************************************
        !
        !    Function:  Recursive method for checking the number of elements
        !               connected to a considered element.
        !
        !     IElement : Element to be considered in current recursion
        !     GroupID : ID of the currently considered group
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement
          integer(INTEGER_TYPE), intent(in) :: GroupID
          ! Local variables
          integer(INTEGER_TYPE) :: NewAdjacentElementID, ISide

          if (CheckCounter<CalParams%RecursionThreshold) then ! Group of elements is not large enough
	      if (IElement>0) then ! There is at least an element, but it might be deactivated
	        if ((IsActiveElement(IElement).and.(ActivationStatus==1)).or. &
                  (.not.IsActiveElement(IElement).and.(ActivationStatus==0))) then ! There is an active/inactive element next to the element, so it is no dangling element
                  if (CheckedElements(IElement)==0) then ! Element not checked yet
                    CheckCounter = CheckCounter + 1
                    CheckedElements(IElement) = GroupID
                    GroupStorage(CheckCounter) = IElement
                    do ISide = 1, ELEMENTSIDES ! Loop over all sides of the considered element
                      ! Retrieve the ID of the adjacent element to be checked next
                      NewAdjacentElementID =  GetAdjacentElement(IElement, ISide)
                      call CheckSideElement(NewAdjacentElementID, GroupID)
                    end do
                  else ! Element has already been checked
                  
                    if (CheckedElements(IElement)/=GroupID) then ! Found element of another group connected to the currently considered one
                      ! Store ID of connected group for later merging, stop recursion at this side
                      NewGroupID = CheckedElements(IElement)
                    end if
                  end if
              end if
            end if
          end if

        end subroutine CheckSideElement

        subroutine CheckGroupMerging(GroupID)
        !**********************************************************************
        !
        !    Function:  Checks whether GroupID should be merged with NewGroupID.
        !               Reduces GroupID if a merging takes place.
        !
        ! O   GroupID : Number of found groups
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(inout) :: GroupID
          ! Local variables
          integer(INTEGER_TYPE) :: I

          if (NewGroupID/=-1) then ! Current group MaxGroupID is connected to NewGroupID
            ! Merge groups
            do I = 1, CalParams%RecursionThreshold
              if (GroupStorage(I)>0) then
                CheckedElements(GroupStorage(I) ) = NewGroupID
              end if
            end do
            
            ! Reduce MaxGroupID in order to keep number of groups small
            GroupID = GroupID - 1
          end if
        
        end subroutine CheckGroupMerging

        subroutine SortElementGroups(MaxGroupID, GroupSizes)
        !**********************************************************************
        !
        !    Function:  Determines the number of found groups and sorts the groups
        !               according to their size.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: MaxGroupID
          integer(INTEGER_TYPE), dimension(MaxGroupID), intent(out) :: GroupSizes
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, GroupCounter

          ! Determine size of groups and number of groups
          GroupSizes = 0
          GroupCounter = 0
          do IElement = 1, Counters%NEl
            if ( ((IsActiveElement(IElement).and.(ActivationStatus==1)).or. &
                  (.not.IsActiveElement(IElement).and.(ActivationStatus==0))).and. &
                 (CheckedElements(IElement)/=0) ) then
              if (GroupSizes(CheckedElements(IElement) )==0) then
                GroupCounter = GroupCounter + 1
              end if
              GroupSizes(CheckedElements(IElement) ) = GroupSizes(CheckedElements(IElement) ) + 1
            end if
          end do
          
          ! Set initial number of void and filled groups, output
          if (ActivationStatus==0) then ! Deactivated elements
            call WriteInLogFile(trim(String(GroupCounter)) // ' void volumes.')
          else ! Activated elements
            call WriteInLogFile(trim(String(GroupCounter)) // ' filled volumes.')
          end if
          
        end subroutine SortElementGroups

        subroutine CheckHoles()
        !**********************************************************************
        !
        !    Function:  Checks which elements belong to holes inside the solid body by
        !               determining which groups have less than CalParams%GroupThreshold elements
        !               and activates them.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
     
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IElement, IGroup, GroupSize, IError, I, J
          integer(INTEGER_TYPE) :: NActivatedElements, AElCounter
          real(REAL_TYPE), dimension(:, :), allocatable :: HoleCoord
     
          Holes = 0 ! Reset holes array
         
          allocate(HoleCoord(MaxHoleID, NVECTOR), stat = IError)
          HoleCoord = 0.0
         
          IsEmptyElement = .false.
          NEmptyElements = 0
          NActivatedElements = 0
          do IElement = 1, Counters%NEl
            if (.not.IsActiveElement(IElement)) then ! Check for deactivated elements
              IGroup = CheckedElements(IElement)
              if (IGroup>0) then ! Element assigned to group
                GroupSize = HoleSizes(IGroup)
                if ( (GroupSize>0).and. &
                     (.not.SkipElementActivation(IElement)).and. &
                     (GroupSize<CalParams%GroupThreshold).and. &
                     (.not.(IsNearBoundaryParticles(IElement))) ) then
                  ! Found a hole, activate elements
                  CheckedElements(IElement) = 0 ! Set to not checked
                  IsActiveElement(IElement) = .true. ! Activate element
                  NActivatedElements = NActivatedElements + 1
                  if (NVirtualParticlesEle(IElement)==0) then ! Deactivated element containing no virtual particles
                    Holes(IElement) = IGroup ! Store group ID of found hole for later use
                    IsEmptyElement(IElement) = .true. ! Mark as empty activated element
                    
                    HoleCoord(IGroup,:) = HoleCoord(IGroup,:) + ElementCentrePoints(IElement,:)
                    
                    NEmptyElements = NEmptyElements + 1
                    if (CalParams%OutputDebugData) then
                      call WriteInLogFile(' Found empty element ' //  trim(String(IElement))) ! Output
                    end if
                  else ! Deactivated element of hole with virtual particles - remove it from the hole
                    HoleSizes(IGroup) = HoleSizes(IGroup) - 1
                  end if ! Else, the element is just activated again, nothing else happens with it
                else
                  if ((GroupSize>0).and.(GroupSize<CalParams%GroupThreshold).and. &
                      IsNearBoundaryParticles(IElement)) then
                    HoleSizes(IGroup) = HoleSizes(IGroup) - 1
                  end if
                end if
                
                if (SkipElementActivation(IElement) ) then
                  HoleSizes(IGroup) = HoleSizes(IGroup) - 1
                end if
                
              end if
            end if
          end do
     
          if (NActivatedElements>0) then
            if (allocated(ActiveElement)) then
              deallocate(ActiveElement, stat = IError)
            end if
            Counters%NAEl = Counters%NAEl + NActivatedElements
            allocate(ActiveElement(Counters%NAEl), stat = IError)
            AElCounter = 0
            do IElement = 1, Counters%NEl
              if (IsActiveElement(IElement)) then
                AElCounter = AElCounter + 1
                ActiveElement(AElCounter) = IElement
              end if
            end do
            
            do I = 1, MaxHoleID
              if (HoleSizes(I)>0) then
                do J = 1, NVECTOR
                  HoleCoord(I, J) = HoleCoord(I, J) / HoleSizes(I)
                end do
              end if
            end do
            
          end if
          if ((CalParams%OutputDebugData).and.(NEmptyElements>0)) then
            call WriteInLogFile('  Found ' //  trim(String(NEmptyElements)) // ' empty elements.')
          end if
          
          deallocate(HoleCoord, stat = IError)

        end subroutine CheckHoles

        logical function IsNearBoundaryParticles(ElementID)
        !**********************************************************************
        !
        !    Function: Identify if element is isnearboundaryParticles
        !
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ElementID
          ! Local variables
          integer(INTEGER_TYPE) :: ISide, AdjacentElementID, NElementsWithBoundaryParticles

          NElementsWithBoundaryParticles = 0
          do ISide = 1, ELEMENTSIDES ! Loop over all sides of the considered element
            AdjacentElementID = GetAdjacentElement(ElementID, ISide)
            if (AdjacentElementID/=0) then
              if (HasBoundaryParticle(AdjacentElementID)) then
                NElementsWithBoundaryParticles = NElementsWithBoundaryParticles + 1
              end if
            end if
          end do

          IsNearBoundaryParticles = (NElementsWithBoundaryParticles>=1)

        end function IsNearBoundaryParticles

        subroutine AssignMaterialToHoles()
        !**********************************************************************
        !
        !    Function:  Assigns a material ID to all elements of the found holes.
        !               The material ID is determined from the elements surrounding the
        !               hole. The material ID that occurs most often in surrounding elements
        !               is used.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IElement, IGroup, IError, CheckElement, ElementCounter
          integer(INTEGER_TYPE), dimension(:), allocatable :: CheckElements
          integer(INTEGER_TYPE), dimension(:), allocatable :: GroupElements

          allocate(CheckElements(Counters%NEl), stat = IError)
          CheckElements = Holes ! Copy Holes array for noting which elements have been checked

          do IElement = 1, Counters%NEl
            if ( (Holes(IElement)>0) .and. (CheckElements(IElement)>0) ) then ! Loop over holes
              
              ! Set up array containing element ID's of the hole
              IGroup = Holes(IElement)
              allocate(GroupElements(HoleSizes(IGroup) ), stat = IError)
              
              ! Determine elements of considered group and mark elements belonging to it as checked
              ElementCounter = 1
              do CheckElement = IElement, Counters%NEl
                if (Holes(CheckElement)==IGroup) then
                  CheckElements(CheckElement) = 0 ! Set to checked
                  GroupElements(ElementCounter) = CheckElement
                  ElementCounter = ElementCounter + 1
                end if
              end do

              ! Assign material to elements of the hole
              call AssignMaterialToHole(HoleSizes(IGroup), GroupElements)
              
              deallocate(GroupElements, stat = IError)
              
            end if
          end do

          deallocate(CheckElements, stat = IError)

        end subroutine AssignMaterialToHoles

        subroutine CheckDanglingElements()
        !**********************************************************************
        !
        !    Function:  Checks which elements are activated outside the solid body by
        !               determining which groups have less than Flags%GroupThreshold elements
        !
        !     GroupSizes : Array containing the size of the groups in SortedElements
        !     MaxGroupID : Maximum number of found groups (most likely less groups were found)
        !
        ! O   IsElm : Element switches
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
     
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IElement, IGroup, GroupSize
     
          IsDanglingElement = .false.
          do IElement = 1, Counters%NEl
            if (IsActiveElement(IElement)) then ! Check for activated elements
              IGroup = CheckedElements(IElement)
              if (IGroup>0) then ! Element assigned to group
                GroupSize = DanglingBodySizes(IGroup)
                if ( (GroupSize>0) .and. (GroupSize<CalParams%GroupThreshold) ) then ! Found a hole, activate elements
                  IsDanglingElement(IElement) = .true. ! Mark as dangling activated element
                  if (CalParams%OutputDebugData) then
                    call WriteInLogFile(' Found dangling element ' // trim(String(IElement))) ! Output
                  end if
                end if
              end if
            end if
          end do
     
        end subroutine CheckDanglingElements

        logical function IsAllElementsChecked()
        !**********************************************************************
        !
        !    Function:  Check whether all elements have been assigned to a group
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IElement

          IsAllElementsChecked = .true.

          do IElement = 1, Counters%NEl
            if (CheckedElements(IElement)==0) then
              IsAllElementsChecked = .false.
              if (CalParams%OutputDebugData) then
                call WriteInLogFile('Element ' //  trim(String(IElement)) // ' not checked ' // 'for empty or dangling elements.')
              end if
            end if
          end do
        
        end function IsAllElementsChecked

        subroutine AssignMaterialToHole(NGroupElements, GroupElements)
        !**********************************************************************
        !
        !    Function:  Assigns a material ID to all elements of the hole whose
        !               elements are specified in GroupElements.
        !               The material ID is determined from the elements surrounding the
        !               hole. The material ID that occurs most often in surrounding elements
        !               is used.
        !
        !     NGroupElements : Number of elements belonging to the hole
        !     GroupElements : ID's of elements forming the hole
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NGroupElements
          integer(INTEGER_TYPE), dimension(NGroupElements), intent(in) :: GroupElements
          ! Local variables
          integer(INTEGER_TYPE) :: NSurroundingElements, IError
          integer(INTEGER_TYPE), dimension(:), allocatable :: SurroundingElements

            ! Search for all activated elements surrounding / connected to the group
            call SearchSurroundingElements(Counters%NEl, NGroupElements, GroupElements, NSurroundingElements, SurroundingElements)

            ! Check material of surrounding elements and assign material to elements of hole
            call SetMaterialHoleElements(NGroupElements, GroupElements, NSurroundingElements, SurroundingElements)

            if (allocated(SurroundingElements) ) then
              deallocate(SurroundingElements, stat = IError)
            end if
        
        end subroutine AssignMaterialToHole

        subroutine SetMaterialHoleElements(NGroupElements, GroupElements, NSurroundingElements, SurroundingElements)
        !**********************************************************************
        !
        !    Function:  Checks which material ID's particles of surrounding elements
        !               carry and assigns the material ID that occurs most often to
        !               elements of the hole.
        !
        !     NGroupElements : Nnumber of elements forming the hole
        !     GroupElements : ID's of the elements forming the hole
        !     NSurroundingElements : Number of elements surrounding the hole
        !     SurroundingElements : ID's of the surrounding elements
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NGroupElements
          integer(INTEGER_TYPE), dimension(NGroupElements), intent(in) :: GroupElements
          integer(INTEGER_TYPE), intent(in) :: NSurroundingElements
          integer(INTEGER_TYPE), dimension(NSurroundingElements), intent(in) :: SurroundingElements
          ! Local variables
          integer(INTEGER_TYPE), dimension(MAXNMATERIALS) :: OccurrenceMat
          integer(INTEGER_TYPE) :: ParticleIndex, IElement, AdjacentElement, IParticle, NMaxMSet, MaxMSet, I, EntityID, MaterialID
        
          ! Determine material that occurs most often in elements surrounding the hole
          OccurrenceMat = 0
          ParticleIndex = -1
          do IElement = 1, NSurroundingElements
            AdjacentElement = SurroundingElements(IElement)
        
            do IParticle = 1, NPartEle(AdjacentElement) ! Loop over particles inside surrounding elements
              ParticleIndex = GetParticleIndex(IParticle, AdjacentElement)
              OccurrenceMat(MaterialIDArray(ParticleIndex)) = OccurrenceMat(MaterialIDArray(ParticleIndex)) + 1
            end do
        
          end do

          NMaxMSet = 0
          MaxMSet = 0
          do I = 1, MAXNMATERIALS          
            if ( (.not.CalParams%ApplyMeshSmoothing) .or. (I/=CalParams%MovingMesh%MovingMaterialID) ) then ! In case of mesh adjustment the structure material should not be considered
              if (OccurrenceMat(I)>NMaxMSet) then
                NMaxMSet = OccurrenceMat(I)
                MaxMSet = I
              end if
            end if
          end do
          
          if (ParticleIndex<0) then
            call WriteInLogFile(' Error determining particle adjacent to hole ' // trim(String(IDArray(ParticleIndex))))
          end if
          
          ! Assign material data that occurs most often to the empty elements 
          do I = 1, NGroupElements
            IElement = GroupElements(I)
            
            do MaterialID = 1, Counters%NLayers
              if (MaterialID==MaxMSet) then
                MaterialElements(MaterialID, IElement) = 1
                If(CalParams%ApplyContactAlgorithm) then
                  if (MaterialID == CalParams%MovingMesh%StructureMaterialID) then
                    EntityID = HARD_ENTITY
                  else
                    EntityID = SOFT_ENTITY
                  end if
                else
                    EntityID = 1
                end if
                
                EntityElements(EntityID, IElement) = 1
              else
                MaterialElements(MaterialID, IElement) = 0
              end if                
            end do 
            
          end do      
        
        end subroutine SetMaterialHoleElements

        subroutine OutputMaximumHoleSizes()
        !**********************************************************************
        !
        !    Function:  Writes sizes of the largest holes to text file.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: I, NLargeHoles, Counter, IError
          integer(INTEGER_TYPE), dimension(:), allocatable :: LargeHoles

          NLargeHoles = 0
          do I = 1, MaxHoleID
            if (HoleSizes(I)>CalParams%GroupThreshold) then
              NLargeHoles = NLargeHoles + 1
            end if
          end do
          
          allocate(LargeHoles(NLargeHoles), stat = IError)
          LargeHoles = 0
          
          Counter = 0
          do I = 1, MaxHoleID
            if (HoleSizes(I)>CalParams%GroupThreshold) then
              Counter = Counter + 1
              LargeHoles(Counter) = HoleSizes(I)
            end if
          end do
          
          call WriteInLogFile('Found hole sizes above threshold: ' // trim(String(LargeHoles,1, NLargeHoles)))
          
          deallocate(LargeHoles, stat = IError)
            
        end subroutine OutputMaximumHoleSizes

        subroutine ReadSHE()
        !**********************************************************************
        !
        !    Function:  Reads ID's of elements who should not be considered as
        !               part of a hole - who will not be activated. 
        !               If this file exists .true. is returned, else .false. .
        !
        !     FileName : Name of the SHE file to open
        !
        ! O   ReadSHE : Returns .true. if the SHE file exists for the considered project
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          character(len = 255) :: CompleteFileName
          integer(INTEGER_TYPE) :: IError, I

          if (.not.CalParams%ApplyEmptyElements) RETURN

          NSkipEleActivation = 0
          
          CompleteFileName = trim(CalParams%FileNames%ProjectName)//'.SHE'

          if (FExist(CompleteFileName) ) then
            call WriteInLogFile('  Reading SHE file')
            call FileOpen(1, CompleteFileName)
            
            read(1, *) NSkipEleActivation
            
            if (allocated(SkipEleActivation) ) then
              deallocate(SkipEleActivation, stat = IError)
            end if
            allocate(SkipEleActivation(NSkipEleActivation), stat = IError)
            
            do I = 1, NSkipEleActivation
              read(1, *) SkipEleActivation(I)
            end do
            
            close(1)
          end if
        
        end subroutine ReadSHE     
       
      end module ModEmptyElements
