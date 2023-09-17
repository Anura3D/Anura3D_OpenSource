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


	  module ModMeshAdjacencies
      !**********************************************************************
      !
      !    Function:  Contains routines for determining element-element and node-element
      !               adjacencies for Finite Element meshes. Furthermore, the centrepoints and
      !               volumes of all elements and a bounding box of the mesh are stored in this module.
      !
      !               The first routines provide access to the adjacency information (GetXXX),
      !               followed by routines for generating the adjacency information from
      !               element connectivity information.
      !
      !               The routine forming the starting point for generating the adjacency
      !               information is: DetermineAdjacencies.
      !               Currently, this modules supports 4- and 10-noded tetrahedral elements.
      !
      !     Implemented in the frame of the MPM project.
      !
      !     $Revision: 8842 $
      !     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModElementEvaluation
      use ModCounters
      use ModMeshInfo
      use ModElementEvaluation
      use ModIsort
      use ModMatrixMath
      
      implicit none
      
        integer(INTEGER_TYPE), parameter :: NRECURSIONS = 50
      
        real(REAL_TYPE), dimension(:, :), allocatable :: ElementCentrePoints ! Centrepoints of elements
        real(REAL_TYPE), dimension(:), allocatable :: ElementSpace ! Surface(2D) or Volume(3D) of the elements
        real(REAL_TYPE), dimension(:), allocatable :: ElementLMin ! minmum altitude of elements
        real(REAL_TYPE), dimension(:,:), allocatable :: ElementStrain ! Element strain
        real(REAL_TYPE), dimension(:, :), allocatable :: StrainSmooth
        real(REAL_TYPE), dimension(:, :), allocatable :: LiquidPressureSmooth 
        real(REAL_TYPE), dimension(:, :), allocatable :: ElementVector 
        real(REAL_TYPE), dimension(:, :), allocatable :: ElementMatrix

        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElementAdjacencies ! Element adjacencies 1 to ELEMENTSIDES for each element
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: EntityElements      !  For each entity, gives the active elements
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: MaterialElements  ! For each material ID, gives the active elements

        integer(INTEGER_TYPE), dimension(:), allocatable :: NodeElementAdjacencies ! Node-element adjacencies
        ! Offsets assembled from the number of elements attached to each node, needed to access NodeElementAdjacencies
        integer(INTEGER_TYPE), dimension(:), allocatable :: NodeElementOffsets
        ! Number of elements attached to each node, needed to access NodeElementAdjacencies
        integer(INTEGER_TYPE), dimension(:), allocatable :: ElementsPerNode

        integer(INTEGER_TYPE), dimension(:), allocatable :: ElementElementAdjacencies ! ID of elements adjacent to each element
        ! Offsets assembled from the number of elements attached to each element, needed to access ElementElementAdjacencies
        integer(INTEGER_TYPE), dimension(:), allocatable :: ElementElementOffsets
        integer(INTEGER_TYPE), dimension(:), allocatable :: ElementsPerElement ! Number of elements adjacent to each element

        integer(INTEGER_TYPE), dimension(:), allocatable :: NodeNodeAdjacencies ! ID of corner (!) nodes adjacent to each corner node
        ! Offsets assembled from the number of corner (!) nodes attached to each corner node
        integer(INTEGER_TYPE), dimension(:), allocatable :: NodeNodeOffsets
        integer(INTEGER_TYPE), dimension(:), allocatable :: NodesPerNode ! Number of corner (!) nodes adjacent to each corner node
        
        logical, dimension(:), allocatable :: IsElementCornerNode ! Stores whether a node is a corner node
        
        real(REAL_TYPE), dimension(:, :), allocatable :: NodeNormals ! Normed normal vectors for each node on the mesh boundary
        
        real(REAL_TYPE), dimension(:), allocatable :: BoundingBoxMin, &
                                                      BoundingBoxMax ! Bounding box of the mesh
        real(REAL_TYPE) :: BoundingBoxDiameter = 0.0

      
    contains ! Routines of this module



        logical function IsInBoundingBox(GlobPos)
                !*************************************************************************************   
        !    FUNCTION:     IsInBoundingBox
        ! 
        !    DESCRIPTION: Returns .true. if GlobPos lies inside or on the bounding box (specified outside)
        !
        !>   @param[in] GlobPos : Point in mesh
        !
        !>   @return IsInBoundingBox : Returns .true. if GlobPos lies inside or on the bounding box
        !
        !*************************************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
          ! Local variables
          integer :: I

          IsInBoundingBox = .true.
          do I = 1, size(GlobPos)
            IsInBoundingBox = IsInBoundingBox .and. (GlobPos(I) >= BoundingBoxMin(I)) .and. (GlobPos(I) <= BoundingBoxMax(I))
          end do          
          
        end function IsInBoundingBox

        
        subroutine SearchSurroundingElements(NEl, NGroupElements, GroupElements, NSurroundingElements, SurroundingElements)
        !**********************************************************************
        !
        !    Function:  Returns a list of elements that enclose the set defined
        !               by GroupElements.
        !
        ! I   NEl : Total number of elements
        ! I   NGroupElements : Number of elements belonging to the set
        ! I   GroupElements : ID's of elements forming the hole
        !
        ! O   NSurroundingElements : Number of elements surrounding the given set of elements
        ! O   SurroundingElements : ID's of the surrounding elements
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NEl, NGroupElements
          integer(INTEGER_TYPE), dimension(NGroupElements), intent(in) :: GroupElements
          integer(INTEGER_TYPE), intent(out) :: NSurroundingElements
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(out) :: SurroundingElements
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, IElement, IError, AdjacentElement, NElementSurroundingElements
          integer(INTEGER_TYPE), dimension(NEl) :: AllElements
        
          if (allocated(SurroundingElements) ) then
            deallocate(SurroundingElements, stat = IError)
          end if
          
          AllElements = 0
        
          ! Single out elements surrounding the set of elements
          do I = 1, NGroupElements
            IElement = GroupElements(I)
            
            NElementSurroundingElements = GetNElmOfElm(IElement)
            do J = 1, NElementSurroundingElements
              AdjacentElement = GetElmIOfElm(IElement, J)
              AllElements(AdjacentElement) =  AllElements(AdjacentElement) + 1
            end do
          end do
        
          ! AllElements now marks all elements belonging to the hole or surrounding it
          do I = 1, NGroupElements
            IElement = GroupElements(I)
            AllElements(IElement) = 0
          end do
        
          ! Determine number of surrounding elements
          NSurroundingElements = 0
          do IElement = 1, NEl
            if (AllElements(IElement)>0) then
              NSurroundingElements = NSurroundingElements + 1
            end if
          end do
        
          ! Write ID's of surrounding elements to array
          allocate(SurroundingElements(NSurroundingElements), stat = IError)
          I = 1
          do IElement = 1, NEl
            if (AllElements(IElement)>0) then
              SurroundingElements(I) = IElement
              I = I + 1
            end if
          end do
        
        end subroutine SearchSurroundingElements

        
        subroutine DetermineElementsAdjacentSurface(IElTyp, NEl, NSurfaceNodes, IConSurface, ICon, IsElm, NElements, IElements, ISides)
        !**********************************************************************
        !
        !    Function:  Returns the number of activated elements adjacent to
        !               the surface with connectivities IConSurface.
        !               For each element, the ID is returned and the ID of the
        !               element side on which the specified surface lies.
        !
        ! I   IElTyp : Number of node connectivities of IElement
        ! I   NEl : Number of elements
        ! I   NSurfaceNodes : Number of surface connectivities
        ! I   IConSurface : Surface connectivities
        ! I   ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        ! I   IsElm : Element switches
        !
        ! O   NElements : Number of elements attached to NodeID
        ! O   IElements : ID's of the found elements
        ! O   ISides : ID's of the sides to which the surface is adjacent
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NSurfaceNodes
          integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(in) :: IConSurface
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          logical, dimension(NEl), intent(in) :: IsElm
          integer(INTEGER_TYPE), intent(out) :: NElements
          integer(INTEGER_TYPE), dimension(2), intent(out) :: IElements ! dimension=2 for 2D and 3D
          integer(INTEGER_TYPE), dimension(2), intent(out) :: ISides ! dimension=2 for 2D and 3D
          ! Local variables
          integer(INTEGER_TYPE) :: IEl, INode, SNode, Count
          
          NElements = 0
          IElements = 0
          ISides = 0
          
          do IEl = 1, NEl
            if (IsElm(IEl)) then ! Loop over all active elements
              
              Count = 0
              do INode = 1, IElTyp ! Loop over all element connectivities
                do SNode = 1, NSurfaceNodes ! Loop over all surface connectivities
                  if (ICon(INode, IEl)==IConSurface(SNode) ) then
                    Count = Count + 1
                  end if
                end do
              end do
              
              if (Count==NSurfaceNodes) then
                NElements = NElements + 1
                IElements(NElements) = IEl
                ISides(NElements) = DetermineElementSideSurface(IEl, NEl, IElTyp, NSurfaceNodes, ICon, IConSurface)
              end if
              
            end if
          end do
     
        end subroutine DetermineElementsAdjacentSurface
        
        
        integer(INTEGER_TYPE) function DetermineElementSideSurface(IEl, NEl, IElTyp, NSurfaceNodes, ICon, IConSurface)
                !*************************************************************************************   
        !    FUNCTION:     DetermineElementSideSurface
        ! 
        !    DESCRIPTION:        
        !>   Returns the ID of the element side which is adjacent to the surface specified by IConSurface
        !
        !>   @param[in] IEl : ID of the considered element
        !>   @param[in] NEl : Number of elements
        !>   @param[in] IElTyp: Number of node connectivities of IElement
        !>   @param[in] NSurfaceNodes : Number of surface connectivities
        !>   @param[in] ICon: Element connectivities ICon(I,J): global node number of local node I in element J
        !>   @param[in] IConSurface : surface connectivities
        !
        !>   @return DetermineElementSideSurface : ID of the element side adjacent to the checked surface
        !
        !*************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IEl, NEl, IElTyp, NSurfaceNodes
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(in) :: IConSurface
          ! Local variables
          integer(INTEGER_TYPE) :: ISide, SNode, Count, i, INeighbour, INodeLoad, INode
          integer(INTEGER_TYPE) :: NNodes
          integer(INTEGER_TYPE), dimension(ELEMENTSIDES) :: ElementSideNeighbours
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodesOfNeighbour
          integer(INTEGER_TYPE), dimension(:,:), allocatable :: CheckEdgeNodes ! (NNodes, ELEMENTSIDES)
          integer(INTEGER_TYPE), dimension(:), allocatable :: CheckNodes ! (NNodes)
          logical :: found

          if ( NDIM == 2 ) then
              ElementSideNeighbours = ElementAdjacencies(IEl, :) !(IElement, ISide)
              DetermineElementSideSurface = 0
              do ISide = 1, ELEMENTSIDES
                INeighbour = ElementSideNeighbours(ISide)  
                if ( INeighbour /= 0 ) then
                   NodesOfNeighbour = ICon(:, INeighbour)  
                   Count = 0
               
                   do INode = 1, ELEMENTNODES
                     found = .false.  
                     do INodeLoad = 1, ELEMENTBOUNDARYNODES  
                       if ( NodesOfNeighbour(INode) == IConSurface(INodeLoad) ) then
                         found = .true.  
                       end if
                     end do  
                     if ( found ) then
                       Count = Count + 1
                     end if  
                   end do
               
                   if ( Count == ELEMENTBOUNDARYNODES ) then
                     DetermineElementSideSurface = ISide
                     RETURN
                   end if 
               
                end if  
              end do    
          end if
          
          if ( NDIM == 3 ) then
              NNodes = 2
              allocate(CheckEdgeNodes(NNodes, ELEMENTSIDES), CheckNodes(NNodes))
          
              call DetermineCheckEdgeNodes(CheckEdgeNodes)

              DetermineElementSideSurface = 0
              do ISide = 1, ELEMENTSIDES
                CheckNodes(:) = ICon(CheckEdgeNodes(:, ISide), IEl)

                Count = 0
                do SNode = 1, NSurfaceNodes
                  found = .false.
                  do i=1,NNodes
                    if (CheckNodes(i) == IConSurface(SNode)) found = .true.
                  enddo
                  if (found) then
                    Count = Count + 1
                  end if
                end do

                if (Count == NNodes) then
                  DetermineElementSideSurface = ISide
                  RETURN
                end if

              end do
        
              deallocate(CheckEdgeNodes, CheckNodes)
           end if   
          
        end function DetermineElementSideSurface


        integer(INTEGER_TYPE) function GetNNodesOfNode(NodeID)
                !*************************************************************************************   
        !    FUNCTION:     GetNNodesOfNode
        ! 
        !    DESCRIPTION:        
        !>   Returns the number of nodes connected to the node with NodeID. If NodeID does not lie in the range of the
        !>   total number of nodes, 0 is returned.
        !
        !>   @param[in] NodeID : ID of the node, whose number of connected nodes is returned
        !
        !>   @return GetNNodesOfNode : Number of nodes attached to NodeID
        !
        !*************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodeID
          
          if ( (allocated(NodesPerNode) ) .and. (NodeID >= 1).and. (NodeID <= ubound(NodesPerNode, 1) ) ) then
            GetNNodesOfNode = NodesPerNode(NodeID)
          else
            GetNNodesOfNode = 0
          end if
        
        end function GetNNodesOfNode


        integer(INTEGER_TYPE) function GetNodeIOfNode(NodeID, INode)
                !*************************************************************************************   
        !    FUNCTION:     GetNodeIOfNode
        ! 
        !    DESCRIPTION:        
        !>   Returns the ID of node INode connected to node NodeID. INode is a local counter ranging from 1 to the 
        !>   number of nodes connected to the node which is determined from GetNNodesOfNode.  If INode exceeds this 
        !>   range, -1 is returned. INode equal 1 returns the lowest node ID.
        !
        !>   @param[in] NodeID : ID of the node, whose attached node INode is returned
        !>   @param[in] INode : Local number of the node connected to the node
        !
        !>   @return GetNodeIOfNode : Node ID of node INode attached to node NodeID
        !
        !*************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodeID, INode
                    
          if ( (allocated(NodesPerNode)) .and. (allocated(NodeNodeOffsets)) .and. (allocated(NodeNodeAdjacencies)) .and. &
               (NodeID >= 1) .and. (NodeID <= ubound(NodesPerNode, 1)) .and. (INode >= 1) .and. (INode <= NodesPerNode(NodeID)) ) then
            GetNodeIOfNode = NodeNodeAdjacencies(NodeNodeOffsets(NodeID) + INode)
          else
            GetNodeIOfNode = -1
          end if
        
        end function GetNodeIOfNode
 

        integer(INTEGER_TYPE) function GetNElmOfNode(NodeID)
                !*************************************************************************************   
        !    FUNCTION:     GetNElmOfNode
        ! 
        !    DESCRIPTION:        
        !>   Returns the number of elements attached to the node with NodeID. If NodeID does not lie in the range of 
        !>   the total number of nodes, 0 is returned.
        !
        !>   @param[in] NodeID : ID of the node, whose number of attached elements is returned
        !
        !>   @return GetNElmOfNode : Number of elements attached to NodeID.
        !
        !*************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodeID
          
          if ( (allocated(ElementsPerNode)) .and. (NodeID >= 1) .and. (NodeID<=ubound(ElementsPerNode, 1)) ) then
            GetNElmOfNode = ElementsPerNode(NodeID)
          else
            GetNElmOfNode = 0
          end if
        
        end function GetNElmOfNode
        

        integer(INTEGER_TYPE) function GetElmIOfNode(NodeID, IElement)
                !*************************************************************************************   
        !    FUNCTION:     GetElmIOfNode
        ! 
        !    DESCRIPTION:        
        !>   Returns the ID of element IElement attached to node NodeID. IElement is a local counter ranging from 1 to
        !>   the number of elements attached to the node which is determined from GetNElmOfNode. If IElement exceeds
        !>   this range, -1 is returned. IElement equal 1 returns the lowest element ID.
        !
        !>   @param[in] NodeID : ID of the node, whose attached element IElement is returned
        !>   @param[in] IElement : Local number of the element attached to the node
        !
        !>   @return GetElmIOfNode : Element ID of element IElement attached to node NodeID
        !
        !*************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodeID, IElement
          
          if ( (allocated(ElementsPerNode)) .and. (allocated(NodeElementOffsets)) .and. (allocated(NodeElementAdjacencies)) .and. &
               (NodeID >= 1) .and. (NodeID <= ubound(ElementsPerNode, 1)) .and. (IElement >= 1) .and. (IElement <= ElementsPerNode(NodeID) ) ) then
            GetElmIOfNode = NodeElementAdjacencies(NodeElementOffsets(NodeID) + IElement)
          else
            GetElmIOfNode = -1
          end if
        
        end function GetElmIOfNode
        

        integer(INTEGER_TYPE) function GetNElmOfElm(ElementID)
                !*************************************************************************************   
        !    FUNCTION:     GetNElmOfElm
        ! 
        !    DESCRIPTION:        
        !>   Returns the number of elements attached to the element with ElementID. If ElementID does not lie in the
        !>   range of the total number of elements, 0 is returned.
        !
        !>   @param[in] ElementID : ID of the element, whose number of attached elements is returned
        !
        !>   @return GetNElmOfElm : Number of elements attached to ElementID
        !
        !*************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ElementID
          
          if ( (allocated(ElementsPerElement)) .and. (ElementID >= 1) .and. (ElementID <= ubound(ElementsPerElement, 1)) ) then
            GetNElmOfElm = ElementsPerElement(ElementID)
          else
            GetNElmOfElm = 0
          end if
        
        end function GetNElmOfElm
        

        integer(INTEGER_TYPE) function GetElmIOfElm(ElementID, IElement)
                !*************************************************************************************   
        !    FUNCTION:     GetElmIOfElm
        ! 
        !    DESCRIPTION:        
        !>   Returns the ID of element IElement attached to element ElementID. IElement is a local counter ranging from
        !>   1 to the number of elements attached to the element which is determined from GetNElmOfElm. If IElement 
        !>   exceeds this range, -1 is returned. IElement equal 1 returns the lowest element ID.
        !
        !>   @param[in] ElementID : ID of the element, whose attached element IElement is returned
        !>   @param[in] IElement : Local number of the element attached to the element
        !
        !>   @return GetElmIOfElm : ElementID of element IElement attached to element ElementID
        !
        !*************************************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: ElementID, IElement
          
          if ( (allocated(ElementsPerElement)) .and. (allocated(ElementElementOffsets)) .and. (allocated(ElementElementAdjacencies)) .and. &
               (ElementID >= 1) .and. (ElementID <= ubound(ElementsPerElement, 1)) .and. (IElement >= 1) .and. (IElement <= ElementsPerElement(ElementID)) ) then
            GetElmIOfElm = ElementElementAdjacencies(ElementElementOffsets(ElementID) +  IElement)
          else
            GetElmIOfElm = -1
          end if

        end function GetElmIOfElm
        
        

        integer(INTEGER_TYPE) function GetAdjacentElement(IElement, ISide)
                !*************************************************************************************   
        !    FUNCTION:     GetAdjacentElement
        ! 
        !    DESCRIPTION:        
        !>   Returns the ID of the element adjacent to ISide of IElement. Returns 0, if no such element exists
        !>   (boundary of the mesh). Returns -999 if something went wrong.
        !
        !>   @param[in] IElement : ID of the considered element
        !>   @param[in] ISide : Local number of the considered side (1 .. ELEMENTSIDES)
        !
        !>   @return GetAdjacentElement : Element ID of the adjacent element (if it exists)
        !
        !*************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IElement, ISide
        
          if ( (allocated(ElementAdjacencies)) .and. (ISide >= 1) .and. (ISide <= ELEMENTSIDES) ) then
            GetAdjacentElement = ElementAdjacencies(IElement, ISide)
          else
            GetAdjacentElement = -999
            call WriteInLogFile('Trying to access non-existing information: GetAdjacentElement')
          end if
        
        end function GetAdjacentElement
        

        integer(INTEGER_TYPE) function BoundaryElementSurface(IElement, ISide, IsElm, NEl)
                !*************************************************************************************   
        !    FUNCTION:     BoundaryElementSurface
        ! 
        !    DESCRIPTION:        
        !>   Returns 0 if side ISide of IElement lies inside a group of activated elements. \n Returns -1 if side ISide
        !>   lies on the boundary of the mesh. \n Returns 1 if the element adjacent to side ISide of IElement is 
        !>   deactivated. \n Returns -999 if something went wrong.
        !
        !>   @param[in] IElement : Global ID of the considered element
        !>   @param[in] ISide : ID of the considered side (1 .. ELEMENTSIDES)
        !>   @param[in] IsElm : Element switches
        !>   @param[in] NEl : Total number of elements
        !
        !>   @return BoundaryElementSurface : Returns the status as described above of the considered element side
        !
        !*************************************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement, ISide, NEl
          logical, dimension(NEl), intent(in) :: IsElm
          ! Local variables
          integer(INTEGER_TYPE) :: AdjacentElement
          
          AdjacentElement = GetAdjacentElement(IElement, ISide)
          
          if (AdjacentElement == 0) then ! No adjacent element exists, boundary of the mesh
            BoundaryElementSurface = -1
          else if (AdjacentElement  >0) then ! Adjacent element exists, check for element switch
            if (.not.IsElm(AdjacentElement)) then
              ! Adjacent element is switched off, ISide lies on the boundary of activated elements
              BoundaryElementSurface = 1
            else
              ! Adjacent element is switched on, ISide lies inside activated elements
              BoundaryElementSurface = 0
            end if
          else ! Something went wrong
            BoundaryElementSurface = -999
            call WriteInLogFile('Error in calling routine: BoundaryElementSurface')
          end if
        
        end function BoundaryElementSurface

        subroutine UpdateMeshAdjacencyInformation(InitialNodalCoord)
        !**********************************************************************
        !
        !    Function:  Updates the element centre points data list, the element volumes data list and
        !               the mesh bounding box data (after nodal coordinates changed).
        !
        ! I   InitialNodalCoord : Nodal coordinates
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR),  intent(in) :: InitialNodalCoord

          call DetermineElementCentrePoints(InitialNodalCoord)
          call DetermineElementSpace(InitialNodalCoord, ElementConnectivities)
          call MeshBoundingBox(InitialNodalCoord, BoundingBoxMin, BoundingBoxMax, BoundingBoxDiameter)
        
        end subroutine UpdateMeshAdjacencyInformation


        subroutine DetermineAdjacencies()
        !**********************************************************************
        !
        !    Function:  Determine the element-element adjacencies and
        !               node-element adjacencies as well as the bounding box
        !               of the mesh as well as the centre point and space
        !               of each element.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), dimension(:,:), allocatable :: ICon
        
          if (ELEMENTTYPE == TETRAOLD) then
            allocate(ICon(N_NODES_HOE, Counters%NEl))  
            ICon = ElementConnectivities10Node ! element connectivities
          else
            allocate(ICon(ELEMENTNODES, Counters%NEl))  
            ICon = ElementConnectivities ! element connectivities  
          end if
          
          call InitialiseAdjacencyArrays()

          call MeshBoundingBox(NodalCoordinates, BoundingBoxMin, BoundingBoxMax, BoundingBoxDiameter)
     
          call DetermineElementCentrePoints(NodalCoordinates)
          
          call InitialiseGPGlobalPositionArrays() ! only required for axisymmetric analysis

          call DetermineElementSpace(NodalCoordinates, ICon)
          
          call DetermineElementMatrix(NodalCoordinates, ElementMatrix, ElementVector)
          
          call DetermineNodeElementAdjacencies()
        
          call DetermineSurroundingElements()

          call DetermineElementSideAdjacencies()

          call DetermineCornerNodes()

          call DetermineNodeNodeConnections()

        end subroutine DetermineAdjacencies


        subroutine InitialiseAdjacencyArrays()
        !**********************************************************************
        !    SUBROUTINE: InitialiseAdjacencyArrays
        ! 
        !    DESCRIPTION: 
        !>   Allocates the arrays needed for storing the adjacency information.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError, NLayers
          
          if (CalParams%ApplyContactMeshBoundary) then
            NLayers = 2
          else
            NLayers = Counters%NLayers ! number of material sets
          end if
          
          call DestroyAdjacencyArrays()
          
          allocate( ElementCentrePoints(Counters%NEl, NVECTOR), stat = IError )
          allocate( ElementSpace(Counters%NEl), stat = IError )
          allocate( ElementLMin(Counters%NEl), stat = IError ) 
          allocate( ElementStrain(Counters%NEl, NTENSOR), stat = IError )
          allocate( EntityElements(Counters%NEntity, Counters%NEl), stat = IError )
          allocate( MaterialElements(NLayers, Counters%NEl), stat = IError )
          allocate( StrainSmooth(Counters%NodTot, 2), stat = IError )
          allocate( LiquidPressureSmooth(Counters%NodTot, 2), stat = IError )
          allocate( ElementVector(Counters%NEl, NVECTOR), stat = IError )
          allocate( ElementMatrix(Counters%NEl, NVECTOR * NVECTOR), stat = IError )
          allocate( ElementAdjacencies(Counters%NEl, ELEMENTSIDES), stat = IError )
          allocate( NodeElementOffsets(Counters%NodTot + 1), stat = IError )
          allocate( ElementsPerNode(Counters%NodTot), stat = IError )
          allocate( ElementElementOffsets(Counters%NEl + 1), stat = IError )
          allocate( ElementsPerElement(Counters%NEl), stat = IError )
          allocate( NodeNodeOffsets(Counters%NodTot + 1), stat = IError )
          allocate( NodesPerNode(Counters%NodTot), stat = IError )
          allocate( IsElementCornerNode(Counters%NodTot), stat = IError )
          allocate( NodeNormals(Counters%NodTot, NVECTOR), stat = IError )          
          allocate( BoundingBoxMin(NVECTOR), BoundingBoxMax(NVECTOR), stat = IError )
          BoundingBoxMin = 0.0
          BoundingBoxMax = 0.0
          
        end subroutine InitialiseAdjacencyArrays
        

        subroutine DestroyAdjacencyArrays()
        !**********************************************************************
        !
        !    Function:  Deallocates the arrays used in this module.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(ElementCentrePoints) ) then
            deallocate(ElementCentrePoints, stat = IError)
          end if

          if (allocated(ElementSpace) ) then
            deallocate(ElementSpace, stat = IError)
          end if
          
           if (allocated(ElementLMin) ) then  
            deallocate(ElementLMin, stat = IError)
           end if
          
          if (allocated(ElementStrain) ) then
            deallocate(ElementStrain, stat = IError)
          end if

           if (allocated(EntityElements) ) then
            deallocate(EntityElements, stat = IError)
          end if
          

          if (allocated(StrainSmooth) ) then
            deallocate(StrainSmooth, stat = IError)
          end if
          
          if (allocated(LiquidPressureSmooth) ) then
            deallocate(LiquidPressureSmooth, stat = IError) 
          end if

          if (allocated(ElementVector) ) then
            deallocate(ElementVector, stat = IError)
          end if

          if (allocated(ElementMatrix) ) then
            deallocate(ElementMatrix, stat = IError)
          end if

          if (allocated(ElementAdjacencies) ) then
            deallocate(ElementAdjacencies, stat = IError)
          end if

          if (allocated(NodeElementAdjacencies) ) then
            deallocate(NodeElementAdjacencies, stat = IError)
          end if

          if (allocated(NodeElementOffsets) ) then
            deallocate(NodeElementOffsets, stat = IError)
          end if

          if (allocated(ElementsPerNode) ) then
            deallocate(ElementsPerNode, stat = IError)
          end if

          if (allocated(ElementElementAdjacencies) ) then
            deallocate(ElementElementAdjacencies, stat = IError)
          end if

          if (allocated(ElementElementOffsets) ) then
            deallocate(ElementElementOffsets, stat = IError)
          end if

          if (allocated(ElementsPerElement) ) then
            deallocate(ElementsPerElement, stat = IError)
          end if

          if (allocated(NodeNodeAdjacencies) ) then
            deallocate(ElementsPerElement, stat = IError)
          end if

          if (allocated(NodeNodeOffsets) ) then
            deallocate(NodeNodeOffsets, stat = IError)
          end if

          if (allocated(NodesPerNode) ) then
            deallocate(NodesPerNode, stat = IError)
          end if
        
          if (allocated(NodeNormals) ) then
            deallocate(NodeNormals, stat = IError)
          end if
        
          if (allocated(IsElementCornerNode)) then
            deallocate(IsElementCornerNode, stat = IError)
          end if
          
        end subroutine DestroyAdjacencyArrays
        

        subroutine DetermineElementSideAdjacencies()
        !**********************************************************************
        !
        !    Function:  Creates a 2D array 1 .. Counters%NEl which stores all ELEMENTSIDES adjacent elements
        !               of an element IElement.
        !               Requires the data about node-element adjacencies.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, ISide, IError

          ! Initialise resulting array
          if (.not.allocated(ElementAdjacencies) ) then
            allocate(ElementAdjacencies(Counters%NEl, ELEMENTSIDES), stat = IError)
          end if
          ElementAdjacencies = -1
          
          do IElement = 1, Counters%NEl ! loop over all elements
            do ISide = 1, ELEMENTSIDES ! loop over all sides of the element
              if (ELEMENTTYPE == TETRAOLD) then ! only 3D, in sprint#2 add 2D functionality
                ElementAdjacencies(IElement, ISide) = CheckEdgeNodeAdjacencies3D(IElement, ISide, N_NODES_HOE, Counters%NEl, ElementConnectivities10Node)
              else
                ElementAdjacencies(IElement, ISide) = CheckEdgeNodeAdjacencies(IElement, ISide) ! only for 3-noded triangle
              end if
            end do
          end do
          
        end subroutine DetermineElementSideAdjacencies

        

        integer(INTEGER_TYPE) function CheckEdgeNodeAdjacencies(IElement, ISide)
                !*************************************************************************************   
        !    FUNCTION:     CheckEdgeNodeAdjacencies
        ! 
        !    DESCRIPTION:        
        !>   Looks for an element which is adjacent to two VERTEX nodes of ISide of IElement and must therefore be 
        !>   adjacent to ISide. The returned value is either: 0 if no adjacent element has been found, ISide lies on 
        !>   the boundary of the mesh, or the element ID of the adjacent element.
        !
        !> @note : for 3-noded triangle and 4-noded quadrilateral
        !
        !>   @param[in] IElement : Considered element
        !>   @param[in] ISide : Considered side of IElement
        !
        !>   @return CheckEdgeNodeAdjacencies
        !
        !*************************************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IElement, ISide
          ! Local variables
          integer(INTEGER_TYPE), parameter :: NNodes = 2 ! number of nodes of each element side, which is a line for 2D
          integer(INTEGER_TYPE), dimension(NNodes, ELEMENTSIDES) :: CheckVertexNodes ! vertex nodes of each element side
          integer(INTEGER_TYPE) :: I, J
          integer(INTEGER_TYPE) :: NodeID1, NodeID2
          integer(INTEGER_TYPE) :: NAdjacentElementsNode1, NAdjacentElementsNode2
          integer(INTEGER_TYPE) :: CheckElementID1, CheckElementID2

          CheckVertexNodes = reshape((/ 1, 2, & ! vertex nodes to check for side 1 spanned by nodes 1-2
                                3, 1, & ! vertex nodes to check for side 2 spanned by nodes 3-4
                                2, 3 /), & ! vertex nodes to check for side 3 spanned by nodes 2-3
                                (/ 2, ELEMENTSIDES /) )

          ! get number of adjacent elements of local vertex node 1, i.e. CheckVertexNodes(1, ISide)
          NodeID1 = ElementConnectivities(CheckVertexNodes(1, ISide), IElement)
          NAdjacentElementsNode1 = GetNElmOfNode(NodeID1)
          ! get number of adjacent elements of local vertex node 2, i.e. CheckVertexNodes(2, ISide)
          NodeID2 = ElementConnectivities(CheckVertexNodes(2, ISide), IElement)
          NAdjacentElementsNode2 = GetNElmOfNode(NodeID2)

          CheckEdgeNodeAdjacencies = 0 ! No adjacent element found

          do I = 1, NAdjacentElementsNode1
            do J = 1, NAdjacentElementsNode2
              CheckElementID1 = GetElmIOfNode(NodeID1, I)
              CheckElementID2 = GetElmIOfNode(NodeID2, J)
              if ( (CheckElementID1 /= IElement) .and. (CheckElementID1==CheckElementID2) ) then ! adjacent element found
                CheckEdgeNodeAdjacencies = CheckElementID1
                EXIT
              end if
            end do
          end do

        end function CheckEdgeNodeAdjacencies


        integer(INTEGER_TYPE) function CheckEdgeNodeAdjacencies3D(IElement, ISide, IElTyp, NEl, ICon)
                
        !*************************************************************************************   
        !    FUNCTION:     CheckEdgeNodeAdjacencies3D
        ! 
        !    DESCRIPTION:        
        !>   Looks for an element which is adjacent to two EDGE nodes of ISide of IElement and must therefore be 
        !>   adjacent to ISide. The returned value is either: 0 if no adjacent element has been found, ISide lies on 
        !>   the boundary of the mesh, or the element ID of the adjacent element.
        !
        !>   @note: 3D function
        !
        !>   @param[in] IElement : Considered element
        !>   @param[in] ISide : Considered side of IElement
        !>   @param[in] IElTyp : Number of nodes per element
        !>   @param[in] NEl : Total number of elements
        !>   @param[in] ICon : Element connectivities
        !
        !>   @return CheckEdgeNodeAdjacencies3D
        !
        !*************************************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IElement, ISide
          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          ! Local variables
          integer(INTEGER_TYPE), parameter :: NNodes = 2 ! for 3D
          integer(INTEGER_TYPE), dimension(NNodes, ELEMENTSIDES) :: CheckEdgeNodes
          integer(INTEGER_TYPE) :: NodeID1, NodeID2, I, J
          integer(INTEGER_TYPE) :: NAdjacentElementsNode1, NAdjacentElementsNode2
          integer(INTEGER_TYPE) :: CheckElementID1, CheckElementID2

          call DetermineCheckEdgeNodes(CheckEdgeNodes)

          ! Get number of adjacent elements of local node CheckEdgeNodes(1, ISide)
          NodeID1 = ICon(CheckEdgeNodes(1, ISide), IElement)
          NAdjacentElementsNode1 = GetNElmOfNode(NodeID1)
          ! Get number of adjacent elements of local node CheckEdgeNodes(2, ISide)
          NodeID2 = ICon(CheckEdgeNodes(2, ISide), IElement)
          NAdjacentElementsNode2 = GetNElmOfNode(NodeID2)

          CheckEdgeNodeAdjacencies3D = 0 ! No adjacent element found

          do I = 1, NAdjacentElementsNode1
            do J = 1, NAdjacentElementsNode2
              CheckElementID1 = GetElmIOfNode(NodeID1, I)
              CheckElementID2 = GetElmIOfNode(NodeID2, J)
              
              if ( (CheckElementID1/=IElement).and. &
                   (CheckElementID1==CheckElementID2) ) then ! Adjacent element found
                CheckEdgeNodeAdjacencies3D = CheckElementID1
                EXIT
              end if
            end do
          end do

        end function CheckEdgeNodeAdjacencies3D


        subroutine DetermineCheckEdgeNodes(CheckEdgeNodes)
        !**********************************************************************
        !
        !    Function:  Determines which edge nodes of each side to check in order
        !               to detect an adjacent element to each side.
        !
        ! O   CheckEdgeNodes : Array containing for each side of the element, 2 local number of edge nodes
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), dimension(:, :), intent(inout) :: CheckEdgeNodes ! size(NNodes, ELEMENTSIDES)
          
          select case(ELEMENTTYPE) ! only for high-order elements for which an edge node exists
            case(TRI6)
              call DetermineCheckEdgeNodesTRI6(CheckEdgeNodes)
            case(QUAD8)
              !call DetermineCheckEdgeNodesQUAD8(CheckEdgeNodes) --> to be added  
            case(TETRA10)
              call DetermineCheckEdgeNodesTETRA10(CheckEdgeNodes)
            case(HEXA20)
              !call DetermineCheckEdgeNodesHEXA20(CheckEdgeNodes) --> to be added  
            case(TETRAOLD)
              call DetermineCheckEdgeNodesTETRA10(CheckEdgeNodes)
          end select
        
        end subroutine DetermineCheckEdgeNodes
        

        subroutine DetermineNodeElementAdjacencies()
        !**********************************************************************
        !
        !    Function:  Determines all element adjacencies of a node.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: NOccurences
          
          call DetermineNodeOccurences(NOccurences)
     
          call DetermineOffsets()
     
          call DetermineNodeElementAdjacenciesFromLists(NOccurences)
     
          call SortNodeElementAdjacencies()
          
        end subroutine DetermineNodeElementAdjacencies
        

        subroutine DetermineCornerNodes()
        !**********************************************************************
        !
        !    Function:  Determines corner nodes of elements
        !
        !**********************************************************************
        
        implicit none
 
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, INode
          
          IsElementCornerNode = .false.
          do IElement = 1, Counters%NEl
            do INode = 1, ELEMENTNODES
              IsElementCornerNode(ElementConnectivities(INode, IElement)) = .true.
            end do
          end do
          
        end subroutine DetermineCornerNodes
        
        
        subroutine DetermineNodeNodeConnections()
        !**********************************************************************
        !
        !    Function:  Determines node-node connections
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: NOccurences
          
          call DetermineNodeNodeOccurences(NOccurences)
          
          call DetermineNodeNodeOffsets()
          
          call DetermineNodeNodeConnectionsFromLists(NOccurences)
          
          call SortNodeNodeConnections()
          
        end subroutine DetermineNodeNodeConnections      
        

        subroutine DetermineNodeNodeOccurences(NOccurences)
        !**********************************************************************
        !
        !    Function:  Determines node-node occurences (only corner nodes)
        !
        ! O    NOccurences : Number of occurences
        !
        !**********************************************************************
 
        implicit none

          integer(INTEGER_TYPE), intent(out) :: NOccurences
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, INode, NodeID, ElmNodeID, NElements, ElementID
          integer(INTEGER_TYPE) :: LowestNodeID, HighestNodeID
          integer(INTEGER_TYPE), dimension(Counters%NodTot) :: TempNodes
        
          NodesPerNode = 0
          do NodeID = 1, Counters%NodTot
            if (IsElementCornerNode(NodeID)) then
              TempNodes = 0
              LowestNodeID = Counters%NodTot
              HighestNodeID = 0
              NElements = GetNElmOfNode(NodeID)
              do IElement = 1, NElements
                ElementID = GetElmIOfNode(NodeID, IElement)
                do INode = 1, ELEMENTNODES
                  ElmNodeID = ElementConnectivities(INode, ElementID)
                  TempNodes(ElmNodeID) = 1
                  if (ElmNodeID<LowestNodeID) then
                    LowestNodeID = ElmNodeID
                  end if
                  if (ElmNodeID>HighestNodeID) then
                    HighestNodeID = ElmNodeID
                  end if
                end do
              end do
              do ElmNodeID = LowestNodeID, HighestNodeID
                if (TempNodes(ElmNodeID)>0) then
                  NodesPerNode(NodeID) = NodesPerNode(NodeID) + 1
                end if
              end do
            end if
          end do
          
          NOccurences = 0
          do INode = 1, Counters%NodTot ! Loop over all entries of NodesPerNode
            NOccurences = NOccurences + NodesPerNode(INode)
          end do
        
        end subroutine DetermineNodeNodeOccurences
        

        subroutine DetermineNodeNodeOffsets()
        !**********************************************************************
        !
        !    Function:  Determines node-node offsets
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: I
          
          NodeNodeOffsets = 0
          do I = 1, Counters%NodTot
            NodeNodeOffsets(I+1) = NodeNodeOffsets(I) + NodesPerNode(I)
          end do
        
        end subroutine DetermineNodeNodeOffsets

        
        subroutine DetermineNodeNodeConnectionsFromLists(NOccurences)
        !**********************************************************************
        !
        !    Function:  Determines node-node connections from lists
        !
        ! I    NOccurences : Number of occurences
        !
        !**********************************************************************
          
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NOccurences
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, INode, NodeID, ElmNodeID, NElements, ElementID, IError
          integer(INTEGER_TYPE) :: LowestNodeID, HighestNodeID, Index
          integer(INTEGER_TYPE), dimension(Counters%NodTot) :: TempNodes

          if (allocated(NodeNodeAdjacencies)) then
            deallocate(NodeNodeAdjacencies, stat = IError)
          end if
          allocate(NodeNodeAdjacencies(NOccurences), stat = IError)
          
          NodesPerNode = 0
          do NodeID = 1, Counters%NodTot
            if (IsElementCornerNode(NodeID)) then
              TempNodes = 0
              LowestNodeID = Counters%NodTot
              HighestNodeID = 0
              NElements = GetNElmOfNode(NodeID)
              do IElement = 1, NElements
                ElementID = GetElmIOfNode(NodeID, IElement)
                do INode = 1, ELEMENTNODES
                  ElmNodeID = ElementConnectivities(INode, ElementID)
                  TempNodes(ElmNodeID) = 1
                  if (ElmNodeID<LowestNodeID) then
                    LowestNodeID = ElmNodeID
                  end if
                  if (ElmNodeID>HighestNodeID) then
                    HighestNodeID = ElmNodeID
                  end if
                end do
              end do
              do ElmNodeID = LowestNodeID, HighestNodeID
                if (TempNodes(ElmNodeID)>0) then
                  NodesPerNode(NodeID) = NodesPerNode(NodeID) + 1
                  Index = NodeNodeOffsets(NodeID) + NodesPerNode(NodeID)
                  NodeNodeAdjacencies(Index) = ElmNodeID
                end if
              end do
            end if
          end do
          
        end subroutine DetermineNodeNodeConnectionsFromLists

        
        subroutine SortNodeNodeConnections()
        !**********************************************************************
        !
        !    Function:  sorts node-node connections
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: INode
        
          do INode = 1, Counters%NodTot
            if (NodesPerNode(INode)/=0) then
              call IHpSort(NodeNodeAdjacencies(NodeNodeOffsets(INode) + 1), NodesPerNode(INode))
            end if
          end do
        
        end subroutine SortNodeNodeConnections
        

        subroutine DetermineNodeOccurences(NOccurences)
        !**********************************************************************
        !
        !    Function:  Determine the number each node occurs in the connectivities
        !               of the elements, the number of elements adjacent to each node.
        !
        ! O   NOccurences : Total number of occurences of nodes in element connectivities
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(out) :: NOccurences
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, INode, NodeID
        
          ElementsPerNode = 0
          do IElement = 1, Counters%NEl ! loop over all elements
            if (ELEMENTTYPE == TETRAOLD) then
            do INode = 1, N_NODES_HOE ! loop over all nodes of the element
              NodeID = ElementConnectivities10Node(INode, IElement)
              ElementsPerNode(NodeID) = ElementsPerNode(NodeID) + 1
            end do
            else
            do INode = 1, ELEMENTNODES ! loop over all nodes of the element
              NodeID = ElementConnectivities(INode, IElement)
              ElementsPerNode(NodeID) = ElementsPerNode(NodeID) + 1
            end do
            end if
          end do
          
          NOccurences = 0
          do INode = 1, Counters%NodTot ! loop over all element nodes of the system
            NOccurences = NOccurences + ElementsPerNode(INode)
          end do
        
        end subroutine DetermineNodeOccurences
        
        
        subroutine DetermineOffsets()
        !**********************************************************************
        !
        !    Function:  Determine the offsets needed to access the information of
        !               NodeElementAdjacencies.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: I

          NodeElementOffsets = 0
          do I = 1, Counters%NodTot ! loop over all element nodes of the system
            NodeElementOffsets(I+1) = NodeElementOffsets(I) + ElementsPerNode(I)
          end do

        end subroutine DetermineOffsets
        

        subroutine DetermineNodeElementAdjacenciesFromLists(NOccurences)
        !**********************************************************************
        !
        !    Function:  Determine the content of the array NodeElementAdjacencies,
        !               which stores at offsets provided by NodeElementOffsets of node
        !               i the IDs of the adjacent elements up to the index provided
        !               by NodeElementOffsets of node i + 1.
        !
        ! I    NOccurences : Total number of elements adjacent to all nodes
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NOccurences
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, INode, NodeID, IError, Index
          
          if ( allocated(NodeElementAdjacencies) ) then
            deallocate(NodeElementAdjacencies, stat = IError)
          end if
          allocate(NodeElementAdjacencies(NOccurences), stat = IError)
          
          ElementsPerNode = 0
          do IElement = 1, Counters%NEl ! loop over all elements
            if (ELEMENTTYPE == TETRAOLD) then
            do INode = 1, N_NODES_HOE ! loop over all nodes of the element
              NodeID = ElementConnectivities10Node(INode, IElement)
              ElementsPerNode(NodeID) = ElementsPerNode(NodeID) + 1
              Index = NodeElementOffsets(NodeID) + ElementsPerNode(NodeID)
              NodeElementAdjacencies(Index) = IElement
            end do
            else    
            do INode = 1, ELEMENTNODES ! loop over all nodes of the element
              NodeID = ElementConnectivities(INode, IElement)
              ElementsPerNode(NodeID) = ElementsPerNode(NodeID) + 1
              Index = NodeElementOffsets(NodeID) + ElementsPerNode(NodeID)
              NodeElementAdjacencies(Index) = IElement
            end do
            end if
          end do

        end subroutine DetermineNodeElementAdjacenciesFromLists

        
        subroutine SortNodeElementAdjacencies()
        !**********************************************************************
        !
        !    Function:  Sort the node element adjacencies. The list of elements
        !               is sorted per node.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: INode
        
          do INode = 1, Counters%NodTot ! loop over all element nodes of the system
            call IHpSort(NodeElementAdjacencies( NodeElementOffsets(INode) + 1),  ElementsPerNode(INode) )
          end do
        
        end subroutine SortNodeElementAdjacencies

        
        subroutine DetermineSurroundingElements()
        !**********************************************************************
        !
        !    Function:  Determines the number of elements surrounding each element and
        !               the ID's of these elements.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError, TotElementCounter
     
          call InitialiseSurroundingElementLists(TotElementCounter)
          
          ! Allocate list for storing the elements adjacent to each element
          if ( allocated(ElementElementAdjacencies) ) then
            deallocate(ElementElementAdjacencies, stat = IError)
          end if
          allocate(ElementElementAdjacencies(TotElementCounter), stat = IError)
        
          call FillElementElementAdjacencies()
          
          call SortElementElementAdjacencies()
        
        end subroutine DetermineSurroundingElements

        
        subroutine InitialiseSurroundingElementLists(TotElementCounter)
        !**********************************************************************
        !
        !    Function:  Determine total number of adjacent elements for all elements.
        !               Determine number of adjacent elements per element (stored in ElementsPerElement list).
        !               Determine offsets for accessing ElementElementAdjacencies (ElementElementOffsets).
        !
        ! O   TotElementCounter : Total number of element adjacencies of the mesh
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(out) :: TotElementCounter
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, ElementCounter, INode, NodeID
          integer(INTEGER_TYPE) :: IAdjacent, AdjacentElementID
          integer(INTEGER_TYPE), dimension(Counters%NEl) :: CheckElements

          TotElementCounter = 0
          ElementElementOffsets = 0
          
          do IElement = 1, Counters%NEl ! loop over all elements
          
            CheckElements = 0
            ElementCounter = 0

            ! Determine which elements are attached to IElement -> CheckElements
            ! Determine the number of elements attached to IElement -> add to ElementCounter
            if (ELEMENTTYPE == TETRAOLD) then
            do INode = 1, N_NODES_HOE
                
              NodeID = ElementConnectivities10Node(INode, IElement)
              if ( IsCornerNode(INode, N_NODES_HOE) ) then ! check adjacent elements only of corner nodes
                do IAdjacent = 1, ElementsPerNode(NodeID)
                  AdjacentElementID = GetElmIOfNode(NodeID, IAdjacent)
                  if ( (AdjacentElementID/=IElement) .and. (CheckElements(AdjacentElementID)==0) ) then ! found new element ID not yet added to ElementCounter
                    CheckElements(AdjacentElementID) = 1
                    ElementCounter = ElementCounter + 1
                  end if
                end do
              end if
              
            end do
            else
            do INode = 1, ELEMENTNODES
                
              NodeID = ElementConnectivities(INode, IElement)
              if ( IsCornerNode(INode, ELEMENTNODES) ) then ! check adjacent elements only of corner nodes
                do IAdjacent = 1, ElementsPerNode(NodeID)
                  AdjacentElementID = GetElmIOfNode(NodeID, IAdjacent)
                  if ( (AdjacentElementID/=IElement) .and. (CheckElements(AdjacentElementID)==0) ) then ! found new element ID not yet added to ElementCounter
                    CheckElements(AdjacentElementID) = 1
                    ElementCounter = ElementCounter + 1
                  end if
                end do
              end if
              
            end do
                
            end if
            

            ElementsPerElement(IElement) = ElementCounter
            ElementElementOffsets(IElement + 1) = ElementElementOffsets(IElement) + ElementsPerElement(IElement)
            TotElementCounter = TotElementCounter + ElementCounter
            
          end do

        end subroutine InitialiseSurroundingElementLists

        
        subroutine FillElementElementAdjacencies()
        !**********************************************************************
        !
        !    Function:  Fill the list of elements adjacent to each element (ElementElementAdjacencies).
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IElement, ElementCounter, INode, NodeID
          integer(INTEGER_TYPE) :: IAdjacent, AdjacentElementID, Index
          integer(INTEGER_TYPE), dimension(Counters%NEl) :: CheckElements

          do IElement = 1, Counters%NEl ! Loop over all elements
          
            CheckElements = 0
            ElementCounter = 0

            ! Add adjacent elements to ElementElementAdjacencies list
            if (ELEMENTTYPE == TETRAOLD) then
            do INode = 1, N_NODES_HOE
                
              NodeID = ElementConnectivities10Node(INode, IElement)
              if ( IsCornerNode(INode, N_NODES_HOE) ) then ! Check adjacent elements only of corner nodes
                do IAdjacent = 1, ElementsPerNode(NodeID)
                  AdjacentElementID = GetElmIOfNode(NodeID, IAdjacent)
                  if ( (AdjacentElementID/=IElement) .and. (CheckElements(AdjacentElementID)==0) ) then ! Found new element ID not yet added to ElementCounter
                    CheckElements(AdjacentElementID) = 1
                    ElementCounter = ElementCounter + 1
                    Index = ElementElementOffsets(IElement) + ElementCounter
                    ElementElementAdjacencies(Index) = AdjacentElementID
                  end if
                end do
              end if
              
            end do
            else
            do INode = 1, ELEMENTNODES
                
              NodeID = ElementConnectivities(INode, IElement)
              if ( IsCornerNode(INode, ELEMENTNODES) ) then ! Check adjacent elements only of corner nodes
                do IAdjacent = 1, ElementsPerNode(NodeID)
                  AdjacentElementID = GetElmIOfNode(NodeID, IAdjacent)
                  if ( (AdjacentElementID/=IElement) .and. (CheckElements(AdjacentElementID)==0) ) then ! Found new element ID not yet added to ElementCounter
                    CheckElements(AdjacentElementID) = 1
                    ElementCounter = ElementCounter + 1
                    Index = ElementElementOffsets(IElement) + ElementCounter
                    ElementElementAdjacencies(Index) = AdjacentElementID
                  end if
                end do
              end if
              
            end do
            end if
            
          end do

        end subroutine FillElementElementAdjacencies

        
        subroutine SortElementElementAdjacencies()
        !**********************************************************************
        !
        !    Function:  Sort the element element adjacencies (ElementElementAdjacencies).
        !               The list of adjacent elements is sorted per element.
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IElement
          
          do IElement = 1, Counters%NEl
            call IHpSort( ElementElementAdjacencies(ElementElementOffsets(IElement) + 1), ElementsPerElement(IElement) )
          end do
        
        end subroutine SortElementElementAdjacencies

        

        subroutine CheckMinMax(Value, Minimum, Maximum)
                !**********************************************************************
        !    SUBROUTINE : CheckMinMax
        !
        !    DESCRIPTION :
        !>   Checks Value(:) against the bounding box, Minimum(:) and Maximum(:),  
        !>   and adjusts the bounding box if the values lie outside
        !     
        !>   @param[in] Value(:) : Values to be checked
        !
        !>   @param[inout] Minimum(:) : Bounding box minimum coordinates
        !>   @param[inout] Maximum(:) : Bounding box maximum coordinates
        !
        !**********************************************************************
     
          real(REAL_TYPE), dimension(:), intent(in) :: Value
          real(REAL_TYPE), dimension(:), intent(inout) :: Minimum, Maximum
     
          ! Local variables
          integer :: I
          
          do I = 1, NVECTOR
              
            if ( Value(I) > Maximum(I) ) then
              Maximum(I) = Value(I)
            end if
            
            if ( Value(I) < Minimum(I) ) then
              Minimum(I) = Value(I)
            end if
            
          end do
     
        end subroutine CheckMinMax
        
        

        function GetBoundingBoxPoint(PointID, Minimum, Maximum) ! 3D function
                !*************************************************************************************   
        !    FUNCTION:     GetBoundingBoxPoint
        ! 
        !    DESCRIPTION:        
        !>   Returns the coordinates of a corner point PointID of the bounding box defined by MinX, ... , MaxZ.
        !
        !>   @note : 3D function
        !
        !>   @param[in] PointID : 1 to 8
        !>   @param[in] Minimum : Bounding box minimum value
        !>   @param[in] Maximum : Bounding box maximum value
        !
        !>   @return CheckBoundingBoxPoint : Returns the corner points of the bounding box
        !
        !*************************************************************************************

        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! 3D function only
          integer(INTEGER_TYPE), intent(in) :: PointID
          real(REAL_TYPE), dimension(IDim), intent(in) :: Minimum, Maximum
          real(REAL_TYPE), dimension(IDim) :: GetBoundingBoxPoint
          
          select case(PointID)
          
            case(1)
              GetBoundingBoxPoint(1) = Minimum(1)
              GetBoundingBoxPoint(2) = Minimum(2)
              GetBoundingBoxPoint(3) = Minimum(3)
            case(2)
              GetBoundingBoxPoint(1) = Maximum(1)
              GetBoundingBoxPoint(2) = Minimum(2)
              GetBoundingBoxPoint(3) = Minimum(3)
            case(3)
              GetBoundingBoxPoint(1) = Maximum(1)
              GetBoundingBoxPoint(2) = Minimum(2)
              GetBoundingBoxPoint(3) = Maximum(3)
            case(4)
              GetBoundingBoxPoint(1) = Minimum(1)
              GetBoundingBoxPoint(2) = Minimum(2)
              GetBoundingBoxPoint(3) = Maximum(3)
            case(5)
              GetBoundingBoxPoint(1) = Minimum(1)
              GetBoundingBoxPoint(2) = Maximum(2)
              GetBoundingBoxPoint(3) = Minimum(3)
            case(6)
              GetBoundingBoxPoint(1) = Maximum(1)
              GetBoundingBoxPoint(2) = Maximum(2)
              GetBoundingBoxPoint(3) = Minimum(3)
            case(7)
              GetBoundingBoxPoint(1) = Maximum(1)
              GetBoundingBoxPoint(2) = Maximum(2)
              GetBoundingBoxPoint(3) = Maximum(3)
            case(8)
              GetBoundingBoxPoint(1) = Minimum(1)
              GetBoundingBoxPoint(2) = Maximum(2)
              GetBoundingBoxPoint(3) = Maximum(3)

          end select

        end function GetBoundingBoxPoint

        
        subroutine ActiveElementsBoundingBox(NodTot, IDim, NodeCoord, NEl, IElTyp, ICon, IsElm, Minimum, Maximum)
        !**********************************************************************
        !
        !    Function:  Determines the bounding box of the activated elements.
        !
        ! I   NodTot : Total number of nodes
        ! I   IDim : Dimension of the mesh
        ! I   NodeCoord : Global nodal coordinates
        ! I   NEl : Number of elements
        ! I   IElTyp : Number of nodes per element
        ! I   ICon : Element node connectivities
        ! I   IsElm : Element switches
        !
        ! O   Minimum : Minimum coordinates
        ! O   Maximum : Maximum coordinates
        !
        !**********************************************************************
     
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, NEl, IElTyp
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          logical, dimension(NEl), intent(in) :: IsElm
          real(REAL_TYPE), dimension(IDim), intent(out) :: Minimum, Maximum
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, NodeID
          
          Minimum = 1E30
          Maximum = -1E30
          
          do I = 1, NEl
            if (IsElm(I)) then ! Active element
              do J = 1, IElTyp
                NodeID = ICon(J, I)
                call CheckMinMax(NodeCoord(NodeID, :), Minimum, Maximum)
              end do
            end if
          end do
     
        end subroutine ActiveElementsBoundingBox

        
        subroutine ElementBoundingBox(ElementID, NodTot, IDim, NodeCoord, NEl, IElTyp, ICon, Minimum, Maximum)
        !**********************************************************************
        !
        !    Function:  Determines the bounding box of the element with ElementID.
        !
        ! I   ElementID : ID of the considered element
        ! I   NodTot : Total number of nodes
        ! I   IDim : Dimension of the mesh
        ! I   NodeCoord : Global nodal coordinates
        ! I   NEl : Number of elements
        ! I   IElTyp : Number of nodes per element
        ! I   ICon : Element connectivities
        !
        ! O   Minimum : Minimum coordinates
        ! O   Maximum : Maximum coordinates
        !
        !**********************************************************************
     
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, NEl, IElTyp, ElementID
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(IDim), intent(out) :: Minimum, Maximum
          ! Local variables
          integer(INTEGER_TYPE) :: I, NodeID
          
          Minimum = 1E30
          Maximum = -1E30
          
          do I = 1, IElTyp ! Loop over element nodes
            NodeID = ICon(I, ElementID)
            call CheckMinMax(NodeCoord(NodeID, :), Minimum, Maximum)
          end do
     
        end subroutine ElementBoundingBox
        

        subroutine CylinderBoundingBox(RotationOrigin, RotationAxis, StartRotationVector, EndRotationVector, &
                                       Minimum, Maximum, MinAxis, MinRadius, StartAngle, MaxAxis, MaxRadius, EndAngle) ! 3D function
        !**********************************************************************
        !
        !    Function:  Translates a bounding box given by MinX, ... , MaxZ into
        !               a bounding box in cylindrical format (MinAxis, ... , EndAngle).
        !
        !               Note: The bounding box is assumed to lie inside the positive sector
        !                     of the local coordinate system defined by RotationOrigin,
        !                     RotationAxis, StartRotationVector and EndRotationVector.
        !                     The bounding box is limited by StartRotationVector and
        !                     EndRotationVector.
        !     Note: 3D Function
        !
        ! I   RotationOrigin : Origin of the cylinder
        ! I   RotationAxis : Axis of the cylinder
        ! I   StartRotationVector,
        ! I   EndRotationVector : Vectors forming a plane perpendicular to RotationAxis
        ! I   Minimum, Maximum : Bounding box
        !
        ! O   MinAxis : Minimum axial coordinate of the bounding box
        ! O   MinRadius : Minimum radius of the bounding box (distance from the axis)
        ! O   StartAngle : Angle between the bounding box and StartRotationVector
        ! O   MaxAxis : Maximum axial coordinate of the bounding box
        ! O   MaxRadius : Maximum radius of the bounding box (distance from the axis)
        ! O   EndAngle : Angle between the bounding box and EndRotationVector
        !
        !**********************************************************************

        implicit none
          
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! 3D function only        
          real(REAL_TYPE), dimension(IDim), intent(in) :: RotationOrigin, RotationAxis, StartRotationVector, EndRotationVector
          real(REAL_TYPE), dimension(IDim), intent(in) :: Minimum, Maximum
          real(REAL_TYPE), intent(out) :: MinAxis, MinRadius,  StartAngle, MaxAxis, MaxRadius, EndAngle
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE), dimension(IDim) :: BoundingBoxPoint, OriginVector, ProjectionVector, Origin
          real(REAL_TYPE) :: ProjectionAxis, ProjectionRadius, StartProjectionAngle, EndProjectionAngle, RotationAngle
        
          MinAxis = 0.0
          MaxAxis = 0.0
          MinRadius = 0.0
          MaxRadius = 0.0
          StartAngle = 0.0
          EndAngle = 0.0
          Origin = 0.0

          do I = 1, 8 ! Loop over 8 points of bounding box
            BoundingBoxPoint =  GetBoundingBoxPoint(I, Minimum, Maximum)
            OriginVector = BoundingBoxPoint - RotationOrigin
          
            ! Maximum value of projection of OriginVector onto RotationAxis (cylinder axis)
            ProjectionAxis = DotProduct(OriginVector, RotationAxis, IDim)
            if (ProjectionAxis<MinAxis) then
              MinAxis = ProjectionAxis
            end if
            if (ProjectionAxis>MaxAxis) then
              MaxAxis = ProjectionAxis
            end if
            
            ! Maximum value of projection of OriginVector onto RotationPlane (cylinder radius)
            ProjectionVector =  VectorPlaneProjection(StartRotationVector, EndRotationVector, OriginVector)
            ProjectionRadius = Length(ProjectionVector, IDim)
            if (ProjectionRadius<MinRadius) then
              MinRadius = ProjectionRadius
            end if
            if (ProjectionRadius>MaxRadius) then
              MaxRadius = ProjectionRadius
            end if
            
            ! Minimum angle between ProjectionVector and StartRotationVector and between ProjectionVector and EndRotationVector
            RotationAngle = dabs(VectorAngle(Origin, StartRotationVector, EndRotationVector, 3)) ! generalise dimension
            StartProjectionAngle = dabs(VectorAngle(Origin, ProjectionVector, StartRotationVector, 3)) ! generalise dimension
            EndProjectionAngle = dabs(VectorAngle(Origin, ProjectionVector, EndRotationVector, 3)) ! generalise dimension

            if ( (StartProjectionAngle<RotationAngle) .and. (EndProjectionAngle<RotationAngle) ) then ! Otherwise the point lies outside the rotation area
              if (StartProjectionAngle<StartAngle) then
                StartAngle = StartProjectionAngle
              end if
            
              if (EndProjectionAngle<EndAngle) then
                EndAngle = EndProjectionAngle
              end if
            end if
          end do
        
        end subroutine CylinderBoundingBox

        !*************************************************************************************   
        !    FUNCTION:     TranslateCylinderToCartesian
        ! 
        !    DESCRIPTION:        
        !>   Translates (RadialFactor, AxialFactor, ArcLengthFactor) into global coordinates (X, Y, Z)
        !
        !>   @note: 3D function
        !
        !>   @param[in] RotationOrigin : Reference point for the rotational sweep
        !>   @param[in] RotationAxis : Rotation axis
        !>   @param[in] StartRotationVector : Vector of the rotation plane setting the start direction
        !>   @param[in] EndRotationVector : Vector of the rotation plane setting the end direction
        !>   @param[in] RadialFactor : Distance of the point from RotationAxis in radial direction
        !>   @param[in] AxialFactor : Distance of the point from RotationOrigin in axial direction
        !>   @param[in] ArcLengthFactor : Arc length distance of the point from StartRotationVector.
        !
        !>   @return TranslateCylinderToCartesian : Enter return variable description here
        !
        !*************************************************************************************
        function TranslateCylinderToCartesian(RotationOrigin, RotationAxis, StartRotationVector, EndRotationVector, &
                                              RadialFactor, AxialFactor, ArcLengthFactor)
     
        implicit none
     
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! 3D function only
          real(REAL_TYPE), dimension(IDim), intent(in) :: RotationOrigin, RotationAxis, StartRotationVector, EndRotationVector
          real(REAL_TYPE), intent(in) :: RadialFactor, AxialFactor, ArcLengthFactor
          real(REAL_TYPE), dimension(IDim) :: TranslateCylinderToCartesian
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StartAngle, EndAngle, RotationAngle
          real(REAL_TYPE), dimension(IDim) :: PlaneVector, Origin
     
          Origin = 0.0
          RotationAngle = dabs(VectorAngle(Origin, StartRotationVector, EndRotationVector, IDim))
          StartAngle = ArcLengthFactor / RadialFactor
          EndAngle = RotationAngle - StartAngle
            
          do I = 1, IDim
            PlaneVector(I) = StartAngle * StartRotationVector(I) + EndAngle * EndRotationVector(I)
          end do
          PlaneVector = VectorNorm(PlaneVector, IDim)
            
          do I = 1, IDim
            TranslateCylinderToCartesian(I) =  RotationOrigin(I) + AxialFactor * RotationAxis(I) +  RadialFactor * PlaneVector(I)
          end do
     
        end function TranslateCylinderToCartesian


        subroutine MeshBoundingBox(NodeCoord, Minimum, Maximum, Diameter)
        !**********************************************************************
        !
        !    Function:  Determines the bounding box of the mesh from NodeCoord, which
        !               might not necessarily correspond to the mesh volume depending on
        !               the shape of the mesh.
        !
        ! I   NodeCoord : Global nodal coordinates
        !
        ! O   Minimum : Minimum coordinates of the mesh
        ! O   Maximum : Maximum coordinates of the mesh
        ! O   Diameter : Diameter of bounding box
        !
        !**********************************************************************
     
        implicit none
        
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          real(REAL_TYPE), dimension(:), intent(inout) :: Minimum, Maximum 

          real(REAL_TYPE), intent(out) :: Diameter
          ! Local variables
          integer(INTEGER_TYPE) :: I
        
          Minimum = 1.E30
          Maximum = -1.E30
          do I = 1, Counters%NodTot ! loop nodes of all elements
            call CheckMinMax( NodeCoord(I, :), Minimum, Maximum )
          end do

          Diameter = 0.0
          do I = 1, size(Minimum)
            Diameter = Diameter + ( Minimum(I) - Maximum(I) ) * ( Minimum(I) - Maximum(I) )
          end do
          Diameter = sqrt(Diameter)     
          
        end subroutine MeshBoundingBox
        

        subroutine DetermineElementCentrePoints(NodeCoord)
        !**********************************************************************
        !
        !    Function:  Determines the centre points for all elements.
        !
        ! I    NodeCoord : Global nodal coordinates
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, INode, NodeID
          real(REAL_TYPE), dimension(ELEMENTVERTICES, NVECTOR) :: Vertices
          
          do IElement = 1, Counters%NEl ! loop all elements
          
            do INode = 1, ELEMENTVERTICES ! loop vertice nodes of each element
                
              NodeID = ElementConnectivities(INode, IElement)
              Vertices(INode, :) = NodeCoord(NodeID, :)
              
            end do
          
            ElementCentrePoints(IElement, :) = PolyhedronCentrePoint(Vertices)
          
          end do

        end subroutine DetermineElementCentrePoints

        
        subroutine DetermineElementSpace(NodeCoord, ICon)
        !**********************************************************************
        !
        !    Function:  Determines the space spanned by each element
        !
        !>  @note : in 2D surface, in 3D volumen
        !
        ! I    NodeCoord : Global nodal coordinates
        ! I    ICon : Element connectivities ICon(I, J)
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, IGaussPoint
          real(REAL_TYPE), dimension(NVECTOR) :: PosGP
          real(REAL_TYPE) :: WeiGP
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: RJac
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: RJacInv
          real(REAL_TYPE) :: DetJac
                  
          ElementSpace = 0.0
          do IElement = 1, Counters%NEl ! loop nodes of all elements
            
            if (ELEMENTTYPE == TETRAOLD) then

              do IGaussPoint = 1, 4 ! loop gauss points of each element (previously was N_GAUSSPOINTS_HOE=4)
                ! determine the location and integration weight of IGaussPoint
                call GaussTETRA_Q4(IGaussPoint, PosGP, WeiGP)
                ! calculate the determinante of the Jacobian matrix
                call DetJacob(PosGP, Counters%NEl, Counters%NodTot, NVECTOR, IElement, ICon, NodeCoord, RJac, RJacInv, DetJac)
                ! determine element volume
                ElementSpace(IElement) = ElementSpace(IElement) + WeiGP * DetJac
              end do
            
            else
              
              do IGaussPoint = 1, ELEMENTGAUSSPOINTS ! loop gauss points of each element
                ! determine the location and integration weight of IGaussPoint
                call GaussPointLocalCoordinates(IGaussPoint, WeiGP, PosGP)
                ! calculate the determinante of the Jacobian matrix
                call DetJacob(PosGP, Counters%NEl, Counters%NodTot, NVECTOR, IElement, ElementConnectivities, NodeCoord, RJac, RJacInv, DetJac)
                
                if ( ISAXISYMMETRIC ) then ! determine element volume for axisymmetric
                  WeiGP = WeiGP * GPGlobalPositionElement(1, IGaussPoint, IElement) ! GPGlobalPositionElement is global variable, index 1 is r-direction
                end if
                
                ! determine element volume
                ElementSpace(IElement) = ElementSpace(IElement) + WeiGP * DetJac
              end do
              
            end if  
          end do
        
        end subroutine DetermineElementSpace
        

!---------------------------------------------------------------------------------
        subroutine DetermineElementMatrix(NodeCoord, EleMat, EleVec)
        !**********************************************************************
        !
        !    Function:  Determines matrix and vector for each element
        !
        ! I    NodeCoord : Global nodal coordinates
        !
        !**********************************************************************
!---------------------------------------------------------------------------------
          implicit none

          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          real(REAL_TYPE), dimension(:, :), intent(inout) :: EleMat
          real(REAL_TYPE), dimension(:, :), intent(inout) :: EleVec
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: M, CofTM, MInv
          real(REAL_TYPE), dimension(NVECTOR) :: X1, MIX1
          real(REAL_TYPE) :: DetV
          integer(INTEGER_TYPE), dimension (ELEMENTNODES) :: NodeID
          integer(INTEGER_TYPE) :: IEl, INode, I, J, K
          
          EleVec = 0.0
          EleMat = 0.0

          do IEl = 1, Counters%NEl   ! Loop over all elements
              
            do INode = 1, ELEMENTNODES
              NodeID(INode) = ElementConnectivities(INode, IEl)    ! Global node ID
            end do
            
            M = 0.0
            X1 = 0.0
            MInv = 0.0
            MIX1 = 0.0
            
            do I = 1, NVECTOR
              X1(I) = NodeCoord(NodeID(1), I)
              do J = 1, NVECTOR
                M(I, J) = NodeCoord(NodeID(J+1), I) - X1(I)
              end do
            end do
            
            call CoFactor(M, NVECTOR, CofTM)
            
            DetV = CalculateDeterminant(M, NVECTOR)
            
            ! calculate inverse matrix
            do I = 1, NVECTOR
              do J = 1, NVECTOR
                MInv(I, J) = CofTM(J, I) / DetV
              end do
            end do
            
            do I = 1, NVECTOR
              do J = 1, NVECTOR
                MIX1(I) = MIX1(I) + MInv(I, J) * X1(J)
              end do
            end do

            K=0
            do I = 1, NVECTOR
            EleVec(IEl, I) = MIX1(I)
              do J = 1, NVECTOR
                K = K + 1
                EleMat(IEl, k) = MInv(I, J)
              end do
            end do
            
          end do
          
        end subroutine DetermineElementMatrix


        subroutine DetermineNodeLocalCS(NodTot, IDim, IElTyp, NEl, NSideNodes, NodeCoord, ICon, ConsideredBoundarySides, NodeNormalN)
        !**********************************************************************
        !
        !    Function:  Determines the normals of nodes lying on the boundary of the
        !               mesh with unit length.
        !               ConsideredBoundarySides can be used to specify which element sides
        !               should be considered for determining the node normals if the normals
        !               are not needed for all nodes - f.e. only for a specific surface.
        !               Afterwards, the tangent vector are determined from the normal vector.
        !
        !               Note: The routines also work with 4-noded low-order tetrahedral elements.
        !                     In this case, NSideNodes is set to 3, so that all loops consider only
        !                     corner nodes.
        !
        !               TODO: Currently for each element only one side can be specified,
        !                     it should work for any assemblage of element sides.
        !
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NSideNodes : Number of nodes per side (4-noded or 10-noded tetrahedral element)
        !     NodeCoord : Global nodal coordinates
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !     ConsideredBoundarySides : Element sides considered for determining the normals
        !
        ! O   NodeNormalN : Normal N for all considered nodes
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, IElTyp, NEl, NSideNodes
          integer(INTEGER_TYPE), dimension(NEl), intent(in) :: ConsideredBoundarySides
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(IDim, NodTot), intent(out) :: NodeNormalN
          ! Local variables
          integer(INTEGER_TYPE) :: NBoundarySides, IError, I
          logical :: DoConsiderAllSides
          integer(INTEGER_TYPE), dimension(:), allocatable :: SideNodeID
          real(REAL_TYPE), dimension(:, :), allocatable :: SideNodeNormals
     
          ! Check whether element sides are specified for determining the node normals
          DoConsiderAllSides = .true.
          do I = 1, NEl
            if (ConsideredBoundarySides(I)>0) then ! Do not consider all element sides
              DoConsiderAllSides = .false.
              EXIT
            end if
          end do
          
          ! Determine number of element sides lying on the boundary of the mesh
          NBoundarySides = DetermineNBoundarySides(NEl, DoConsiderAllSides, ConsideredBoundarySides)
          
          ! Assign to each node of each side a global node ID
          allocate(SideNodeID(NBoundarySides * NSideNodes), stat = IError)
          call DetermineSideNodeID(NEl, IElTyp, ICon, NSideNodes, NBoundarySides * NSideNodes, &
                                   DoConsiderAllSides, ConsideredBoundarySides, SideNodeID)
          
          ! Determine for each node for each side a normal
          allocate(SideNodeNormals(IDim, NBoundarySides * NSideNodes), stat = IError)
          call DetermineSideNodeNormals(NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, NSideNodes, &
                                        NBoundarySides * NSideNodes, DoConsiderAllSides, &
                                        ConsideredBoundarySides, SideNodeNormals)
          
          ! Average the normals and set to unit length
          NodeNormalN = 0.0
          call AverageNodeNormalN(NodTot, IDim, NBoundarySides * NSideNodes, SideNodeID, &
                                  SideNodeNormals, NodeNormalN)
          
          deallocate(SideNodeNormals, stat = IError)
          deallocate(SideNodeID, stat = IError)
     
        end subroutine DetermineNodeLocalCS


        integer(INTEGER_TYPE) function DetermineNBoundarySides(NEl, DoConsiderAllSides, ConsideredBoundarySides)
                !*************************************************************************************   
        !    FUNCTION:     DetermineNBoundarySides
        ! 
        !    DESCRIPTION:        
        !>   Determines the number of element sides lying on the boundary of the mesh. Requires data from array 
        !>   'ElementAdjacencies'.
        !
        !>   @param[in] NEl : Number of elements
        !>   @param[in] DoConsiderAllSides : True, if all element sides should be considered.
        !>   @param[in] ConsideredBoundarySides : Specified which element sides should be considered for determining
        !>   the node normals
        !
        !>   @return DetermineNBoundarySides : Number of element sides lying on the boundary of the mesh
        !
        !*************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NEl
          logical, intent(in) :: DoConsiderAllSides
          integer(INTEGER_TYPE), dimension(NEl), intent(in) :: ConsideredBoundarySides
          ! Local variables
          integer(INTEGER_TYPE) :: I, J

          DetermineNBoundarySides = 0
          do I = 1, NEl
            do J = 1, ELEMENTSIDES
              if ( (DoConsiderAllSides .and. ElementAdjacencies(I, J)==0) .or. (ConsideredBoundarySides(I)==J) ) then
                DetermineNBoundarySides = DetermineNBoundarySides + 1
              end if
            end do
          end do
        
        end function DetermineNBoundarySides


        subroutine DetermineSideNodeID(NEl, IElTyp, ICon, NSideNodes, NBoundarySideNodes, DoConsiderAllSides, ConsideredBoundarySides, SideNodeID)
        !**********************************************************************
        !
        !    Function:  Determines the ID of each node forming the element sides that lie on
        !               the boundary of the mesh.
        !               Requires data of array 'ElementAdjacencies'.
        !
        !     NEl : Number of elements
        !     IElTyp : Number of node connectivities of IElement
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !     NSideNodes : Number of nodes per side
        !     NBoundarySideNodes : Number of nodes of each side lying on the boundary of the mesh
        !     DoConsiderAllSides : True, if all element sides should be considered
        !     ConsideredBoundarySides : Specifies which element sides should be considered for determining the node normals
        !
        ! O   SideNodeID : Global node IDs spanning each side lying on the boundary of the mesh
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NEl, IElTyp, NBoundarySideNodes, NSideNodes
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          logical, intent(in) :: DoConsiderAllSides
          integer(INTEGER_TYPE), dimension(NEl), intent(in) :: ConsideredBoundarySides
          integer(INTEGER_TYPE), dimension(NBoundarySideNodes), intent(out) :: SideNodeID
          ! Local variables
          integer(INTEGER_TYPE) :: IBoundaryNode, I, J, K, LocalNodeID
        
          IBoundaryNode = 0
          if (NDIM == 3) then ! 3D functionality, in sprint#2 add 2D functionality
          do I = 1, NEl
            do J = 1, ELEMENTSIDES
              if ( (DoConsiderAllSides .and. ElementAdjacencies(I, J)==0) .or. (ConsideredBoundarySides(I)==J) ) then
                do K = 1, NSideNodes
                  IBoundaryNode = IBoundaryNode + 1
                  LocalNodeID = DetermineSideNodesTetrahedronHOE(J, K)
                  SideNodeID(IBoundaryNode) = ICon(LocalNodeID, I)
                end do
              end if
            end do
          end do
          end if
        
        end subroutine DetermineSideNodeID

          
        subroutine AverageNodeNormalN(NodTot, IDim, NBoundarySideNodes, SideNodeID, SideNodeNormals, NodeNormalN)
        !**********************************************************************
        !
        !    Function:  Determines the ID of each node forming the element sides that lie on
        !               the boundary of the mesh.
        !
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the vectors
        !     NBoundarySideNodes : Number of nodes of each side lying on the boundary of the mesh
        !     SideNodeID : Global node IDs spanning each side lying on the boundary of the mesh
        !     SideNodeNormals : Normal for each side for each node
        !
        ! O   NodeNormalN : Returns the normal vectors of the considered nodes
        !
        !**********************************************************************
        
          implicit none
          
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, NBoundarySideNodes
          integer(INTEGER_TYPE), dimension(NBoundarySideNodes), intent(in) :: SideNodeID
          real(REAL_TYPE), dimension(IDim, NBoundarySideNodes), intent(in) :: SideNodeNormals
          real(REAL_TYPE), dimension(IDim, NodTot), intent(out) :: NodeNormalN
          ! Local variables
          integer(INTEGER_TYPE) :: I, NodeID
        
          do I = 1, NBoundarySideNodes
            NodeID = SideNodeID(I)
            NodeNormalN(:, NodeID) = NodeNormalN(:, NodeID) + SideNodeNormals(:, I)
          end do

          do I = 1, NodTot
            NodeNormalN(:, I) = -1.0 * VectorNorm(NodeNormalN(:, I), IDim)
          end do
        
        end subroutine AverageNodeNormalN

        
        subroutine DetermineSideNodeNormals(NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, NSideNodes, NBoundarySideNodes, &
                                            DoConsiderAllSides, ConsideredBoundarySides, SideNodeNormals)
        !**********************************************************************
        !
        !    Function:  Determines the normals of nodes for each side lying on the
        !               boundary of the mesh.
        !               Requires data of array 'ElementAdjacencies'.
        !
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NodeCoord : Global nodal coordinates
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !     NSideNodes : Number of nodes per side
        !     NBoundarySideNodes : Number of nodes of each side lying on the boundary of the mesh
        !     DoConsiderAllSides : True, if all element sides should be considered
        !     ConsideredBoundarySides : Specifies which element sides should be considered for determining the node normals
        !
        ! O   SideNodeNormals : Normal for each side for each node
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, IElTyp, NEl,  NBoundarySideNodes, NSideNodes
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          logical, intent(in) :: DoConsiderAllSides
          integer(INTEGER_TYPE), dimension(NEl), intent(in) :: ConsideredBoundarySides
          real(REAL_TYPE), dimension(IDim, NBoundarySideNodes), intent(out) :: SideNodeNormals
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, K, IBoundaryNode

          IBoundaryNode = 0
          if (IDim == 3) then ! 3D functionality, in sprint#2 add 2D functionality
          do I = 1, NEl
            do J = 1, ELEMENTSIDES
              if ( (DoConsiderAllSides .and. ElementAdjacencies(I, J)==0) .or. (ConsideredBoundarySides(I)==J) ) then
                do K = 1, NSideNodes
                  IBoundaryNode = IBoundaryNode + 1

                  call DetermineSideNodeNormalTetrahedron(NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, NSideNodes, &
                                   I, J, K, SideNodeNormals(:, IBoundaryNode) )
                end do
              end if
            end do
          end do
          end if
        
        end subroutine DetermineSideNodeNormals


        subroutine DetermineSideNodeNormalTetrahedron(NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, NSideNodes, &
                                           IElement, ISide, INode, SideNodeNormal) ! 3D function
        !**********************************************************************
        !
        !    Function:  Determines the normal for node INode of side ISide of element
        !               IElement, writes it to SideNodeNormal.
        !    Note: 3D function
        !
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NodeCoord : Global nodal coordinates
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !     NSideNodes : Number of nodes per side
        !     IElement : ID of the considered element
        !     ISide : ID of the considered side
        !     INode : Local number of the considered node
        !
        ! O   SideNodeNormal : Normal at the considered node on the considered side
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, IElTyp, NEl, NSideNodes, IElement, ISide, INode
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(IDim), intent(out) :: SideNodeNormal
          ! Local variables
          integer(INTEGER_TYPE) :: LocalNodeID, I
          real(REAL_TYPE) :: Xi, Eta
          real(REAL_TYPE), dimension(NSideNodes) :: HS
          real(REAL_TYPE), dimension(NSideNodes, 2) :: DHS
          real(REAL_TYPE), dimension(IDim) :: S, T, SideNode

          ! Determine local coordinates of INode
          call DetermineSideNodeLocPosTetrahedronHOE(INode, Xi, Eta)
          
          ! Determine shape function derivatives evaluated at INode
          call ShapeXiEtaT(NSideNodes, Xi, Eta, HS, DHS)

          ! Determine tangential vectors S and T from DN/DXi => S and DN/DEta => T
          S = 0.0
          T = 0.0
          do I = 1, NSideNodes
            LocalNodeID = DetermineSideNodesTetrahedronHOE(ISide, I)
            SideNode(:) = NodeCoord(ICon(LocalNodeID, IElement), :)
            S(:) = S(:) + DHS(I, 1) * SideNode(:)
            T(:) = T(:) + DHS(I, 2) * SideNode(:)
          end do

          ! Cross product S x T (all element side normals point inward, therefore S x T)
          SideNodeNormal = CrossProduct(S, T)

        end subroutine DetermineSideNodeNormalTetrahedron


        subroutine DetermineNodeTangentAxiSymPenetration(NodeNormal, NodeS, NodeT) ! 3D function
        !**********************************************************************
        !
        !    Function:  Determines the tangents NodeS and NodeT from
        !               NodeNormal.
        !    Note: 3D function
        !
        ! I    NodeNormal : Normal vector
        ! O    NodeS : Tangential vector parallel to x-z-plane
        ! O    NodeT : Tangential vector perpendicular to N-S vectors
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed as only 3D functionality
          real(REAL_TYPE), dimension(IDim), intent(in) :: NodeNormal
          real(REAL_TYPE), dimension(IDim), intent(out) :: NodeS, NodeT
          
          NodeS(1) = -1.0 * NodeNormal(3)
          NodeS(2) = 0.0
          NodeS(3) = NodeNormal(1)
          
          NodeS = VectorNorm(NodeS, IDim)
          NodeT = CrossProduct(NodeS, NodeNormal)
          NodeT = VectorNorm(NodeT, IDim)
        
        end subroutine DetermineNodeTangentAxiSymPenetration

        
        subroutine DetermineNodeTangent(NodeNormal, NodeTangentS, NodeTangentT) ! 3D function
        !**********************************************************************
        !
        !    Function:  Determines the tangents NodeS and NodeT from
        !               NodeNormal.
        !    Note: 3D function
        !
        ! I   NodeNormal : Normal vector
        ! O   NodeS : Tangential vector parallel to x-z-plane
        ! O   NodeT : Tangential vector perpendicular to N-S vectors
        !
        !**********************************************************************
          
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3
          real(REAL_TYPE), dimension(IDim), intent(in) :: NodeNormal
          real(REAL_TYPE), dimension(IDim), intent(out) :: NodeTangentS, NodeTangentT

          if ((abs(NodeNormal(1))==0.0).and.(abs(NodeNormal(2))==0.0).and.(abs(NodeNormal(3))==0.0)) then
            call GiveError("Normal with zero length")
          end if
          
          NodeTangentS(1) = -1.0 * NodeNormal(2)
          NodeTangentS(2) = NodeNormal(1)
          NodeTangentS(3) = NodeNormal(3)
          
          NodeTangentT = CrossProduct(NodeTangentS, NodeNormal)
          NodeTangentT = VectorNorm(NodeTangentT, IDim)
          
        end subroutine DetermineNodeTangent
        

        subroutine DetermineGlobPosElement(GlobPos, ElementID, IElTyp, NEl, NodTot, IDim, NodeCoord,  &
                                           ICon, NewElementID, NewLocPos, Success)
        !**********************************************************************
        !
        !    Function:  Determine in which element GlobPos is located
        !               starting the checking from element with ElementID.
        !
        !               Two recursive algorithms have been implemented:
        !               Version 1 checks whether GlobPos lies inside elements surrounding one of the
        !               corner nodes of NewElementID. If no element contains GlobPos the element
        !               whose centrepoint is closest to GlobPos is checked ... .
        !               Version 2 checks whether GlobPos lies inside NewElementID. If not, the element
        !               adjacent to the side, through which GlobPos has left NewElementID is checked ... .
        !               NOTE: Using Version 2 does not require to check whether GlobPos lies inside the
        !               original element of the considered particle!
        !               Version 2 seems to be faster as less elements are checked. Does it always work?
        !
        ! I   GlobPos : Global coordinates of a point inside the mesh
        ! I   ElementID : ID of the element assumed in the vicinity of GlobPos
        ! I   IElTyp : Number of node connectivities of IElement
        ! I   NEl : Number of elements
        ! I   NodTot : Total number of nodes
        ! I   IDim : Dimension of the mesh
        ! I   NodeCoord : Global nodal coordinates
        ! I   ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !
        ! O   NewElementID : Element ID in which GlobPos is located
        ! O   NewLocPos : Local coordinates of the considered point inside NewElementID
        ! O   Success : True, if the element could be detected
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
          integer(INTEGER_TYPE), intent(in) :: ElementID
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          real(REAL_TYPE), dimension(:), intent(inout) :: NewLocPos
          integer(INTEGER_TYPE), intent(out) :: NewElementID
          logical, intent(out) :: Success
          ! Local variables
          logical :: OutsideElement
          integer(INTEGER_TYPE) :: RecursionDepth
          
          NewElementID = ElementID
          Success = .false.
          RecursionDepth = 1
          NewLocPos = 0.2

          if (ELEMENTTYPE == TETRAOLD) then
            select case (IEltyp)
              case(10) ! 10-noded tetrahedral element (curved boundaries)
                call CheckElementForGlobPosVersion1(GlobPos, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, &
                                                     RecursionDepth, NewElementID, NewLocPos, Success)
              case(4) ! 4-noded tetrahedral element
                call CheckElementForGlobPosVersion2(GlobPos, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, &
                                                    RecursionDepth, NewElementID, Success)
            end select
          else
            call CheckElementForGlobPosVersion2(GlobPos, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, &
                                                RecursionDepth, NewElementID, Success)
          end if 

          if (Success) then ! Try determining the local particle coordinates from the global ones
            call GetLocalCoordinates(GlobPos,  NewElementID, IElTyp, NEl,  NodTot, IDim, NodeCoord, ICon, &
                                     NewLocPos, OutsideElement, Success)
          end if
          
        end subroutine DetermineGlobPosElement


        recursive subroutine CheckElementForGlobPosVersion1(GlobPos, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, &
                                                         RecursionDepth, NewElementID, NewLocPos, Success)
        !**********************************************************************
        !
        !    Function:  Recursively checks whether the surrounding elements of NewElementID
        !               contain GlobPos. A maximum of NRECURSIONS recursions are performed.
        !
        ! I   GlobPos : Global coordinates of a point inside the mesh
        ! I   NodTot : Total number of nodes
        ! I   IDim : Dimension of the mesh
        ! I   IElTyp : Number of node connectivities of IElement
        ! I   NEl : Number of elements
        ! I   NodeCoord : Global nodal coordinates
        ! I   ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !
        ! O   RecursiveDepth : Number of performed recursions
        ! O   NewElementID : Element ID in which GlobPos is located,
        !                    initially it stores the ID of the element whose surrounding elements are checked
        ! O   NewLocPos : New local coordinate for GlobPos inside NewElementID
        ! O   Success : True, if the element could be detected
        !
        !**********************************************************************
     
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, IElTyp, NEl
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          integer(INTEGER_TYPE), intent(inout) :: RecursionDepth
          integer(INTEGER_TYPE), intent(inout) :: NewElementID
          real(REAL_TYPE), dimension(IDim), intent(inout) :: NewLocPos
          logical, intent(out) :: Success
          ! Local variables
          real(REAL_TYPE) :: Distance, MinDistance
          integer(INTEGER_TYPE) :: IElement, NElmOfNode, CheckElement, MinDistanceElement
          logical :: OutsideElement
     
          ! Check surrounding elements of corner node with smallest angle whether GlobPos lies inside
          ! and distances between centrepoints of surrounding elements and GlobPos (in case new recursion is necessary)
          MinDistance = 1E30
          MinDistanceElement = -1
          NElmOfNode = GetNElmOfElm(NewElementID)
          do IElement = 1, NElmOfNode
            CheckElement = GetElmIOfElm(NewElementID, IElement)

            if (CheckElement/=NewElementID) then

              ! Check whether GlobPos lies inside CheckElement (also works for curved elements)
              call GetLocalCoordinates(GlobPos, CheckElement, IElTyp, NEl, NodTot, IDim, NodeCoord, ICon, &
                                       NewLocPos, OutsideElement, Success)
              Success = Success.and.(.not.OutsideElement)
              if (Success) then
                NewElementID = CheckElement
                EXIT
              else
                Distance = Length( ElementCentrePoints(CheckElement, :) - GlobPos, IDim)
                if (Distance<MinDistance) then
                  MinDistance = Distance
                  MinDistanceElement = CheckElement
                end if
              end if
            end if
            
          end do
     
          if ( (.not.Success).and.(RecursionDepth<NRECURSIONS) ) then ! Look for next element to check, start next recursion
            RecursionDepth = RecursionDepth + 1
            
            ! Proceed with surrounding element with smallest distance to GlobPos
            NewElementID = MinDistanceElement
            ! Recursive call
            call CheckElementForGlobPosVersion1(GlobPos, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, &
                                                RecursionDepth, NewElementID, NewLocPos, Success)
          end if ! else the recursion ends
     
        end subroutine CheckElementForGlobPosVersion1


        recursive subroutine CheckElementForGlobPosVersion2(GlobPos, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, &
                                                         RecursionDepth, NewElementID, Success)
        !**********************************************************************
        !
        !    Function:  Recursively checks whether the particle lies inside NewElementID.
        !               If not, the element adjacent to the side through which the particle
        !               left the element NewElementID is checked.
        !               A maximum of NRECURSIONS recursions are performed.
        !
        ! I   GlobPos : Global coordinates of a point inside the mesh
        ! I   NodTot : Total number of nodes
        ! I   IDim : Dimension of the mesh
        ! I   IElTyp : Number of node connectivities of IElement
        ! I   NEl : Number of elements
        ! I   NodeCoord : Global nodal coordinates
        ! I   ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !
        ! O   RecursiveDepth : Number of performed recursions
        ! O   NewElementID : Element ID in which GlobPos is located,
        !                    initially it stores the ID of the element whose surrounding elements are checked
        ! O   Success : True, if the element could be detected
        !
        !**********************************************************************
     
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, IElTyp, NEl
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          integer(INTEGER_TYPE), intent(inout) :: RecursionDepth
          integer(INTEGER_TYPE), intent(inout) :: NewElementID
          logical, intent(out) :: Success
          ! Local variables
          integer(INTEGER_TYPE) :: CrossedSide
          real(REAL_TYPE), dimension(IDim) :: CoordCentrePoint
          
          CrossedSide = -1
          CoordCentrePoint = ElementCentrePoints(NewElementID,:)
          
          call CheckForGlobPosPointer(GlobPos, NewElementID, CoordCentrePoint, NodeCoord, ICon, CrossedSide, Success)
          
          if ( (.not.Success).and.(RecursionDepth<NRECURSIONS) ) then ! Start next recursion with adjacent element
            RecursionDepth = RecursionDepth + 1
            NewElementID = GetAdjacentElement(NewElementID, CrossedSide)
            ! Recursive call
            if (NewElementID>0) then
              call CheckElementForGlobPosVersion2(GlobPos, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon, &
                                                  RecursionDepth, NewElementID, Success)
            end if
          else
            if (RecursionDepth==NRECURSIONS) then
              call GiveError('Maximum number of recursions reached')
            endif
          end if ! else the recursion ends

        end subroutine CheckElementForGlobPosVersion2

        
        subroutine InterpolateEdgeNodeCoordinates(NodTot, IDim, LOETyings, NodeCoord)
        !**********************************************************************
        !
        !    Function:  Determines edge node coordinates from corner node coordinates.
        !
        ! I   NodTot : Total number of nodes
        ! I   IDim : Dimension of the mesh
        ! I   LOETyings : Ties edge nodes to adjacent corner nodes
        !
        ! O   NodeCoord : Global nodal coordinates
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim
          integer(INTEGER_TYPE), dimension(2, NodTot), intent(in) :: LOETyings
          real(REAL_TYPE), dimension(NodTot, IDim), intent(inout) :: NodeCoord
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, I1, I3

          do I = 1, NodTot
            I1 = LOETyings(1, I)
            I3 = LOETyings(2, I)
            do J = 1, IDim
              NodeCoord(I, J) = (NodeCoord(I1, J) + NodeCoord(I3, J) ) / 2
            end do
          end do
        
        end subroutine InterpolateEdgeNodeCoordinates
        

        subroutine DetermineElementLMin()   
        !**********************************************************************
        !
        !    FUNCTION:  DetermineElementLMin
        ! 
        !    DESCRIPTION:
        !>   Determines the minimum altitude of each element
        !
        !**********************************************************************
          implicit none
          
          ! local variables
          integer(INTEGER_TYPE) :: NodTot, NEl, IElement
          real(REAL_TYPE), dimension(:,:), allocatable :: NodeCoord
          real(REAL_TYPE), dimension(:), allocatable :: Lmintry
          real(REAL_TYPE) :: Lmin
          
          NodTot = Counters%NodTot ! total number of nodes
          NEl = Counters%NEl ! number of elements
          
          allocate(NodeCoord(NodTot, NDIM), Lmintry(NEl))
          NodeCoord = NodalCoordinatesUpd ! global nodal coordinates
          Lmin = 0.0
          ElementLMin = 0.0
         
          do IElement = 1, NEl ! loop over elements
            call GetMinAltitude(IElement, Lmin, NodeCoord)
            ElementLMin(IElement) = Lmin
          end do

          Lmintry = ElementLMin
        
        end subroutine DetermineElementLMin

        

        subroutine GetMinAltitude(IElement, Lmin, NodeCoord) 
                !**********************************************************************
        !
        !    FUNCTION:  GetMinAltitude
        ! 
        !    DESCRIPTION:
        !>   Determines the minimum altitude of an element
        !
        !>   @param[in] NodeCoord : Global nodal coordinates
        !    @param[in] IElement : the element ID
        !    @param[out] Lmin : minimum altitude of the element
        !
        !**********************************************************************
     
          implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          real(REAL_TYPE), intent (inout):: Lmin
          
          ! local variables
          integer(INTEGER_TYPE) :: NodeNr(ELEMENTNODES)        

          NodeNr(:) = ElementConnectivities(:, IElement) ! get global node number
          
          call GetMinAltitudePointer(NodeNr, NodeCoord, Lmin)

                    
        end subroutine GetMinAltitude
        

        subroutine GetMinAltitudeTetra(NodeNr, NodeCoord, Lmin) 
                !**********************************************************************
        !
        !    FUNCTION:  GetMinAltitudeTetra
        ! 
        !    DESCRIPTION:
        !>   Determines the minimum altitude of an tetrahedral element
        !
        !>   @param[in] NodeNr : node number
        !    @param[in] NodeCoord : Global nodal coordinates
        !    @param[out] Lmin : minimum altitude of the element
        !
        !**********************************************************************
     
          implicit none

          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          real(REAL_TYPE), intent (inout):: Lmin
          
          ! local variables
          integer(INTEGER_TYPE):: NodeNr(:), I
          real(REAL_TYPE), dimension(3):: Vec1, Vec2, Normal, A
          real(REAL_TYPE):: LengthA, LengthN, ADotN, Min
          
          Min = 0.0
          Lmin = 0.0
          A = 0.0
          
          do I = 1, NDIM
          Vec1(I) = NodeCoord(NodeNr(2), I) - NodeCoord(NodeNr(1), I)
          Vec2(I) = NodeCoord(NodeNr(3), I) - NodeCoord(NodeNr(1), I)
          A(I) = NodeCoord(NodeNr(3), I) - NodeCoord(NodeNr(4), I)
          end do
          
          Normal = CrossProduct(Vec1, Vec2)
          LengthA = Length(A, NDIM) 
          LengthN = Length(Normal, NDIM) 
          ADotN = A(1)* Normal(1)+  A(2)* Normal(2)+  A(3)* Normal(3)
          Lmin = abs(ADotN) /(LengthN)
          Min = Lmin
           
          do I = 1, NDIM
          Vec1(I) = NodeCoord(NodeNr(2), I) - NodeCoord(NodeNr(1), I)
          Vec2(I) = NodeCoord(NodeNr(4), I) - NodeCoord(NodeNr(1), I)
          A(I)= NodeCoord(NodeNr(3), I) - NodeCoord(NodeNr(4), I)
          end do
          
          Normal = CrossProduct(Vec1, Vec2)
          LengthA = Length(A, NDIM) 
          LengthN = Length(Normal, NDIM)
          ADotN = A(1)* Normal(1)+  A(2)* Normal(2)+  A(3)* Normal(3)
          Lmin = abs(ADotN) /(LengthN)
           
          if (Lmin<Min) Min = Lmin
           
          do I = 1, NDIM
          Vec1(I) = NodeCoord(NodeNr(3), I) - NodeCoord(NodeNr(2), I)
          Vec2(I) = NodeCoord(NodeNr(4), I) - NodeCoord(NodeNr(2), I)
          A(I)= NodeCoord(NodeNr(1), I) - NodeCoord(NodeNr(4), I)
          end do
          
          Normal = CrossProduct(Vec1, Vec2)
          LengthA = Length(A, NDIM) 
          LengthN = Length(Normal, NDIM)
          ADotN = A(1)* Normal(1)+  A(2)* Normal(2)+  A(3)* Normal(3)
          Lmin = abs(ADotN) /(LengthN)
           
          if (Lmin<Min) Min = Lmin
           
          do I = 1, NDIM
          Vec1(I) = NodeCoord(NodeNr(3), I) - NodeCoord(NodeNr(1), I)
          Vec2(I) = NodeCoord(NodeNr(4), I) - NodeCoord(NodeNr(1), I)
          A(I)= NodeCoord(NodeNr(2), I) - NodeCoord(NodeNr(4), I)
          end do
          
          Normal = CrossProduct(Vec1, Vec2)
          LengthA = Length(A, NDIM)
          LengthN = Length(Normal, NDIM)
          ADotN = A(1)* Normal(1)+  A(2)* Normal(2)+  A(3)* Normal(3)
          Lmin = abs(ADotN) /(LengthN)
           
          if (Lmin<Min) Min = Lmin
           
          Lmin = Min
          
        end subroutine GetMinAltitudeTetra       

        

        subroutine GetMinAltitudeTri(NodeNr, NodeCoord, Lmin) 
                !**********************************************************************
        !
        !    FUNCTION:  GetMinAltitudeTri
        ! 
        !    DESCRIPTION:
        !>   Determines the minimum altitude of an triangular element
        !
        !>   @param[in] NodeNr : node number
        !    @param[in] NodeCoord : Global nodal coordinates
        !    @param[out] Lmin : minimum altitude of the element
        !
        !**********************************************************************
        
          implicit none

          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          real(REAL_TYPE), intent (inout):: Lmin

          ! local variables
          integer(INTEGER_TYPE):: NodeNr(:), I
          real(INTEGER_TYPE):: Idim
          real(REAL_TYPE), dimension(2):: node1, node2, node3, n0, n1, n2
          real(REAL_TYPE) :: d_min, d(3) !d1, d2 and d3 for a triangle


          do IDim = 1, NDIM 
            node1(IDim) = NodeCoord(NodeNr(1), IDim)
            node2(IDim) = NodeCoord(NodeNr(2), IDim)
            node3(IDim) = NodeCoord(NodeNr(3), IDim)
          end do

          d_min = LARGE

          !the distance of node 0 from line formed by node 1 and node 2
          !d = |(x2-x1)*(y1-y0)-(x1-x0)(y2-y1)|/sqrt[(x2-x1)^2+(y2-y1)^2]
          !see: http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html

          n0 = node1
          n1 = node2
          n2 = node3
          d(1) = abs((n2(1)-n1(1))*(n1(2)-n0(2))-(n1(1)-n0(1))*(n2(2)-n1(2)))/&
                sqrt((n2(1)-n1(1))*(n2(1)-n1(1))+(n2(2)-n1(2))*(n2(2)-n1(2)))

          n0 = node2
          n1 = node1
          n2 = node3
          d(2) = abs((n2(1)-n1(1))*(n1(2)-n0(2))-(n1(1)-n0(1))*(n2(2)-n1(2)))/&
                sqrt((n2(1)-n1(1))*(n2(1)-n1(1))+(n2(2)-n1(2))*(n2(2)-n1(2)))

          n0 = node3
          n1 = node1
          n2 = node2
          d(3) = abs((n2(1)-n1(1))*(n1(2)-n0(2))-(n1(1)-n0(1))*(n2(2)-n1(2)))/&
                sqrt((n2(1)-n1(1))*(n2(1)-n1(1))+(n2(2)-n1(2))*(n2(2)-n1(2)))

          do I=1,3
            d_min = min(d(I), d_min)
          end do

          Lmin = d_min
          
        end subroutine
        
        subroutine DetermineReactionElements()
    !********************************************************************************
    !
    !Function: Determine the elements to be considered for the integration of the reaction forces
    !
    !********************************************************************************
    implicit none
    logical, dimension(Counters%NEl) :: IsReactionElm
    integer(INTEGER_TYPE) :: Inode, NElemNod, J, IEle, NReactionElem, Ierror, cont
    
    IsReactionElm = .false.
    NReactionElem = 0
    
    do INode = 1, Counters%NodTot
      if (IsReactionNode(INode)) then
        NElemNod = GetNElmOfNode(INode)
        do J=1,NElemNod
          IEle = GetElmIOfNode(INode,J)
          if (.not.IsReactionElm(IEle)) then
            IsReactionElm (IEle) = .true.
            NReactionElem = NReactionElem + 1
          end if
        end do
       end if
    end do
      
     Counters%NElemReactions = NReactionElem
     
     allocate(ConsideredElemReaction(NReactionElem), stat = Ierror)
     
     cont = 0
    do IEle=1,Counters%NEl
      if(IsReactionElm(IEle))then
        cont = cont + 1
        ConsideredElemReaction(cont) = IEle
      end if
    end do
       
      end subroutine DetermineReactionElements

      end module ModMeshAdjacencies
