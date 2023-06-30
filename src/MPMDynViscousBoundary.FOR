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
	  
	  
	  module ModDynViscousBoudary
      !**********************************************************************
      !
      !  Function : Contains routines for initialising and updating information
      !             related to absorbing boundaries
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModMPMData
      use ModReadCalculationData
      use ModMeshInfo
      use ModCounters      
      use ModRotBoundCond
      use ModMPMStresses
      use ModMeshAdjacencies
      use ModReadGeometryData
      
      implicit none
      
        type ViscousBoundariesParameterType ! File names
          integer(INTEGER_TYPE) :: NVisElemSld = 0, & ! ...
                     NVisElemWat = 0, & ! ...
                     NVisElemGas = 0  ! ...
        end type ViscousBoundariesParameterType 
 
        type(ViscousBoundariesParameterType), public, save :: ViscousBParams

        real(REAL_TYPE), dimension(:, :), allocatable :: NodalTotDisplacementVB ! ...
        real(REAL_TYPE), dimension(:, :), allocatable :: NodalTotDisplacementVBWater ! ...
        real(REAL_TYPE), dimension(:, :), allocatable :: NodalTotDisplacementVBGas    ! ...
        real(REAL_TYPE), dimension(:, :), allocatable :: IncrementalDisplacementGas   ! ...
        ! Arrays related solid
        real(REAL_TYPE), dimension(:, :), allocatable :: VisDampForceSld ! ...
        real(REAL_TYPE), dimension(:), allocatable :: VisNodSurAraSld ! the surface area of viscous nodes
        real(REAL_TYPE), dimension(:), allocatable :: DashpotSld ! dashpot array ! solid
        real(REAL_TYPE), dimension(:), allocatable :: NormDashpotSld ! dashpot array ! solid
        real(REAL_TYPE), dimension(:), allocatable :: SpringSld ! Spring array ! solid
        real(REAL_TYPE), dimension(:), allocatable :: NormSpringSld ! Normalized Spring array ! solid
        real(REAL_TYPE), dimension(:, :), allocatable :: NodeCndVisSld2 ! Node conditions for v2017_2
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: NodeCndVisSld ! Node conditions
        integer(INTEGER_TYPE), dimension(:), allocatable :: DofCndVisSld ! Dof Conditions ! solid
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: VisElmConSld ! connectivities of the viscous face element (solid phase)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElmCondSld ! connectivities of the viscous face element (solid phase)

        ! Arrays relate water
        real(REAL_TYPE), dimension(:, :), allocatable :: VisDampForceWat ! ...
        real(REAL_TYPE), dimension(:), allocatable :: VisNodSurAraWat ! the surface area of viscous nodes
        real(REAL_TYPE), dimension(:), allocatable :: DashpotWat ! dashpot array ! water
        real(REAL_TYPE), dimension(:), allocatable :: SpringWat ! Spring array ! water
        real(REAL_TYPE), dimension(:, :), allocatable :: NodeCndVisWat2 ! Node conditions for v2017_2 
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: NodeCndVisWat ! Node conditions
        integer(INTEGER_TYPE), dimension(:), allocatable :: DofCndVisWat ! Dof Conditions ! water
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: VisElmConWat ! connectivities of the viscous face element (water phase)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElmCondWat ! connectivities of the viscous face element (water phase)

        ! Arrays relate Gas
        real(REAL_TYPE), dimension(:, :), allocatable :: VisDampForceGas ! ...
        real(REAL_TYPE), dimension(:), allocatable :: VisNodSurAraGas ! the surface area of viscous nodes
        real(REAL_TYPE), dimension(:), allocatable :: DashpotGas ! dashpot array ! gas
        real(REAL_TYPE), dimension(:), allocatable :: SpringGas ! Spring array ! gas
        real(REAL_TYPE), dimension(:, :), allocatable :: NodeCndVisGas2 ! Node conditions for v2017_2 
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: NodeCndVisGas ! Node conditions
        integer(INTEGER_TYPE), dimension(:), allocatable :: DofCndVisGas ! Dof Conditions ! gas
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: VisElmConGas ! connectivities of the viscous face element (gas phase)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElmCondGas ! connectivities of the viscous face element (gas phase)

      contains


        subroutine InitialiseAbsorbingBoundaryData
        !**********************************************************************
        !
        !  Function : Contains code for initialising absorbing boundary data
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError, INode, SizeAB, NVisNodSld, NVisNodWat, NVisNodGas
          real(REAL_TYPE), dimension(:), allocatable :: NodeConditionLocal

          call DestroyViscousBoundaryData()

          SizeAB = 1 + 3 * NDIM ! array size for reading absorbing boundary data, see also definition in subroutine ReadGOM_v...()
          allocate(NodeConditionLocal(SizeAB), stat=IError) 
          
          allocate(NodalTotDisplacementVB(Counters%N, Counters%nEntity), stat=IError)
          NodalTotDisplacementVB = 0.0

          allocate(NodalTotDisplacementVBWater(Counters%N, Counters%nEntity), stat=IError)
          NodalTotDisplacementVBWater = 0.0

          allocate(NodalTotDisplacementVBGas(Counters%N, Counters%nEntity), stat=IError)
          NodalTotDisplacementVBGas = 0.0

          allocate(IncrementalDisplacementGas(Counters%N, Counters%nEntity), stat=IError)
          IncrementalDisplacementGas = 0.0

          ! Arrays related to solid phase
          allocate(VisDampForceSld(Counters%N, Counters%nEntity), stat=IError)
          allocate(VisNodSurAraSld(Counters%NodTot), stat=IError)
          allocate(NodeCndVisSld2(Counters%NodTot, SizeAB), stat=IError)
          allocate(DashpotSld(Counters%N), stat=IError)
          allocate(NormDashpotSld(Counters%N), stat=IError)
          allocate(SpringSld(Counters%N), stat=IError)
          allocate(NormSpringSld(Counters%N), stat=IError)
          allocate(DofCndVisSld(Counters%N), stat=IError)
          VisDampForceSld = 0.0
          VisNodSurAraSld = 0.0
          NodeCndVisSld2 = 0
          DashpotSld = 0.0
          SpringSld = 0.0
          NormDashpotSld = 0.0
          NormSpringSld = 0.0
          DofCndVisSld = 0
          NVisNodSld = GeoParams%NABSurfaceNodesSolid
          ViscousBParams%NVisElemSld = NVisNodSld
          
          ! Arrays related to liquid phase
          allocate(VisDampForceWat(Counters%N,Counters%nEntity), stat=IError)
          allocate(VisNodSurAraWat(Counters%NodTot), stat=IError)
          allocate(NodeCndVisWat2(Counters%NodTot, SizeAB), stat=IError)
          allocate(DashpotWat(Counters%N), stat=IError)
          allocate(SpringWat(Counters%N), stat=IError)
          allocate(DofCndVisWat(Counters%N), stat=IError)
          VisDampForceWat = 0.0
          VisNodSurAraWat = 0.0
          NodeCndVisWat2 = 0.0
          DashpotWat = 0.0
          SpringWat = 0.0
          DofCndVisWat = 0.0
          NVisNodWat = GeoParams%NABSurfaceNodesLiquid
          ViscousBParams%NVisElemWat = NVisNodWat

          ! Arrays related to gas phase
          allocate(VisDampForceGas(Counters%N,Counters%nEntity), stat=IError)
          allocate(VisNodSurAraGas(Counters%NodTot), stat=IError)
          allocate(NodeCndVisGas(Counters%NodTot, SizeAB), stat=IError)
          allocate(NodeCndVisGas2(Counters%NodTot, SizeAB), stat=IError)
          allocate(DashpotGas(Counters%N), stat=IError)
          allocate(SpringGas(Counters%N), stat=IError)
          allocate(DofCndVisGas(Counters%N), stat=IError)
          VisDampForceGas = 0.0
          VisNodSurAraGas = 0.0
          NodeCndVisGas = 0.0
          NodeCndVisGas2 = 0.0
          DashpotGas = 0.0
          SpringGas = 0.0
          DofCndVisGas = 0.0
          NVisNodGas = GeoParams%NABSurfaceNodesGas
          ViscousBParams%NVisElemGas = NVisNodGas
          
          ! surfaces solid
          do INode = 1, GeoParams%NABSurfaceNodesSolid  ! loop over viscous elements 
            NodeConditionLocal = GeoParams%AbsorbingBoundariesSurfacesSolid(INode, :) 
            if (IsElementCornerNode(NodeConditionLocal(1))) then
              NodeCndVisSld2(NodeConditionLocal(1), :) = NodeConditionLocal
            end if  
          end do 
          
          ! lines solid
          do INode = 1, GeoParams%NABLineNodesSolid
                NodeConditionLocal = GeoParams%AbsorbingBoundariesLinesSolid(INode, :)
                  if (IsElementCornerNode(NodeConditionLocal(1))) then
                    NodeCndVisSld2(NodeConditionLocal(1), :) = NodeConditionLocal
                  end if  
          end do
            
          ! points solid
          do INode = 1, GeoParams%NABPointNodesSolid
                NodeConditionLocal = GeoParams%AbsorbingBoundariesPointsSolid(INode, :)
                  if (IsElementCornerNode(NodeConditionLocal(1))) then
                    NodeCndVisSld2(NodeConditionLocal(1), :) = NodeConditionLocal
                  end if  
          end do
          
          ! surfaces liquid
          do INode =1, GeoParams%NABSurfaceNodesLiquid  ! loop over viscous elements
               NodeConditionLocal = GeoParams%AbsorbingBoundariesSurfacesLiquid(INode, :)
               if (IsElementCornerNode(NodeConditionLocal(1))) then
                 NodeCndVisWat2(NodeConditionLocal(1), :) = NodeConditionLocal
               end if  
          end do

          ! lines liquid
          do INode = 1, GeoParams%NABLineNodesLiquid
               NodeConditionLocal = GeoParams%AbsorbingBoundariesLinesLiquid(INode, :)
                 if (IsElementCornerNode(NodeConditionLocal(1))) then
                   NodeCndVisWat2(NodeConditionLocal(1), :) = NodeConditionLocal
                 end if  
          end do

          ! points liquid
          do INode = 1, GeoParams%NABPointNodesLiquid               
               NodeConditionLocal = GeoParams%AbsorbingBoundariesPointsLiquid(INode, :)
               if (IsElementCornerNode(NodeConditionLocal(1))) then
                 NodeCndVisWat2(NodeConditionLocal(1), :) = NodeConditionLocal
               end if  
          end do
          
          ! surfaces gas
          do INode =1, GeoParams%NABSurfaceNodesGas  ! loop over viscous elements
               NodeConditionLocal = GeoParams%AbsorbingBoundariesSurfacesGas(INode, :)
               if (IsElementCornerNode(NodeConditionLocal(1))) then
                 NodeCndVisGas2(NodeConditionLocal(1), :) = NodeConditionLocal
               end if  
          end do

          ! lines gas
          do INode = 1, GeoParams%NABLineNodesGas
               NodeConditionLocal = GeoParams%AbsorbingBoundariesLinesGas(INode, :)
                 if (IsElementCornerNode(NodeConditionLocal(1))) then
                   NodeCndVisGas2(NodeConditionLocal(1), :) = NodeConditionLocal
                 end if  
          end do

          ! points gas
          do INode = 1, GeoParams%NABPointNodesGas               
               NodeConditionLocal = GeoParams%AbsorbingBoundariesPointsGas(INode, :)
               if (IsElementCornerNode(NodeConditionLocal(1))) then
                 NodeCndVisGas2(NodeConditionLocal(1), :) = NodeConditionLocal
               end if  
          end do
              
        end subroutine InitialiseAbsorbingBoundaryData
    
        
        subroutine DestroyViscousBoundaryData()
        !**********************************************************************
        !
        !    Function:  Deallocates the arrays used in this module.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(NodalTotDisplacementVB)) then
            deallocate(NodalTotDisplacementVB, stat = IError)
          end if

          if (allocated(NodalTotDisplacementVBWater)) then
            deallocate(NodalTotDisplacementVBWater, stat = IError)
          end if
          
          if (allocated(NodalTotDisplacementVBGas)) then
            deallocate(NodalTotDisplacementVBGas, stat = IError)
          end if
          
          if (allocated(IncrementalDisplacementGas)) then
            deallocate(IncrementalDisplacementGas, stat = IError)
          end if

          if (allocated(VisDampForceSld)) then
            deallocate(VisDampForceSld, stat = IError)
          end if

          if (allocated(VisNodSurAraSld)) then
            deallocate(VisNodSurAraSld, stat = IError)
          end if

          if (allocated(DashpotSld)) then
            deallocate(DashpotSld, stat = IError)
          end if

          if (allocated(SpringSld)) then
            deallocate(SpringSld, stat = IError)
          end if
        
          if (allocated(NormDashpotSld)) then
            deallocate(NormDashpotSld, stat = IError)
          end if

          if (allocated(NormSpringSld)) then
            deallocate(NormSpringSld, stat = IError)
          end if
     
          if (allocated(NodeCndVisSld)) then
            deallocate(NodeCndVisSld, stat = IError)
          end if

          if (allocated(DofCndVisSld)) then
            deallocate(DofCndVisSld, stat = IError)
          end if

          if (allocated(VisElmConSld)) then
            deallocate(VisElmConSld, stat = IError)  
          end if
          
          if (allocated(ElmCondSld)) then
            deallocate(ElmCondSld, stat = IError)  
          end if
          
          if (allocated(VisDampForceWat)) then
            deallocate(VisDampForceWat, stat = IError)
          end if

          if (allocated(VisNodSurAraWat)) then
            deallocate(VisNodSurAraWat, stat = IError)
          end if

          if (allocated(DashpotWat)) then
            deallocate(DashpotWat, stat = IError)
          end if

          if (allocated(SpringWat)) then
            deallocate(SpringWat, stat = IError)
          end if

          if (allocated(NodeCndVisWat)) then
            deallocate(NodeCndVisWat, stat = IError)
          end if

          if (allocated(DofCndVisWat)) then
            deallocate(DofCndVisWat, stat = IError)
          end if

          if (allocated(VisElmConWat)) then 
            deallocate(VisElmConWat, stat = IError)
          end if
          
          if (allocated(ElmCondWat)) then 
            deallocate(ElmCondWat, stat = IError)
          end if

          if (allocated(VisDampForceGas)) then
            deallocate(VisDampForceGas, stat = IError)
          end if

          if (allocated(VisNodSurAraGas)) then
            deallocate(VisNodSurAraGas, stat = IError)
          end if

          if (allocated(DashpotGas)) then
            deallocate(DashpotGas, stat = IError)
          end if

          if (allocated(SpringGas)) then
            deallocate(SpringGas, stat = IError)
          end if

          if (allocated(NodeCndVisGas)) then
            deallocate(NodeCndVisGas, stat = IError)
          end if

          if (allocated(DofCndVisGas)) then
            deallocate(DofCndVisGas, stat = IError)
          end if

          if (allocated(VisElmConGas)) then 
            deallocate(VisElmConGas, stat = IError)
          end if
          
          if (allocated(ElmCondGas)) then 
            deallocate(ElmCondGas, stat = IError)
          end if

        end subroutine DestroyViscousBoundaryData
        
        
        subroutine InitialiseAbsorbingBoundaryDashpotSpring()
        !**********************************************************************
        !   
        !   Function: Initialise data from dashpot and spring in absorbing boundaries
        !
        !**********************************************************************
        
        implicit none
        
          if (.not.CalParams%ApplyAbsorbingBoundary) RETURN ! only if absorbing boundaries are used
          if (.not.NFORMULATION==1) RETURN ! only for single-point formulation
          

            ! used for constant stiffness absorbing boundaries
            call CalculateNodeSurfaceAreaandSetCondSolid()
            call FormDashpotandSpringSolid()
            call CalculateNodeSurfaceAreaandSetCondLiquid()
            call FormDashpotandSpringLiquid()
            call CalculateNodeSurfaceAreaAndSetCondGas()
            call FormDashpotandSpringGas()
 

        end subroutine InitialiseAbsorbingBoundaryDashpotSpring
        

        subroutine CalculateNodeSurfaceAreaandSetCondSolid()
        !**********************************************************************
        !
        !   Function:  To calculate the surface area associated with each node
        !
        ! OUT: populates VisNodSurAraSld(INode), i.e. the surface area of absorbing boundary nodes  
        ! OUT: re-populates NodeCndVisSld2(INode,J), i.e. if absorbing boundary is in normal (1) or shear (2) direction
        !
        !**********************************************************************

        implicit none
        
          ! local variables 
          real(REAL_TYPE), dimension(Counters%NEl) :: ElmSurfaceArea   
          real(REAL_TYPE), dimension(NVECTOR):: Vec1, Vec2, Normal, ScalarProduct, Connector
          integer(INTEGER_TYPE) :: INode, INode2, INode3, ISide, IDim, I, NumberOfElementsNode
          integer(INTEGER_TYPE) :: IsMeshBoundary, ElementID, LocalNodeID
          integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: DummyNode
          logical :: Belongs2AbsorbingBoundary
                                  
          VisNodSurAraSld = 0.0
          ElmSurfaceArea = 0.0   
                 
            do INode = 1, Counters%NodTot                          ! loop all nodes
              if ( NodeCndVisSld2(INode, 1) > 0 ) then             ! check if node is defined as viscous node in GOM-file
                NumberOfElementsNode = ElementsPerNode(INode)      ! get number of elements around this considered node
                do INode2 = 1, NumberOfElementsNode                ! check all elements sharing this node
                  ElementID = NodeElementAdjacencies(NodeElementOffsets(INode)+INode2) ! get the element ID of each neighbouring element
                  do ISide = 1, ELEMENTSIDES                       ! loop around all element sides                       
                    IsMeshBoundary = BoundaryElementSurface(ElementID, ISide, IsActiveElement, Counters%NEl) ! check if elementside is a mesh boundary
                    if ( IsMeshBoundary == -1 ) then               ! check which side of the element is part of the boundary
                      Belongs2AbsorbingBoundary = .true.
                      do INode3 = 1, ELEMENTBOUNDARYNODES          ! get the node IDs of this side
                        if ( NDIM == 3 ) then
                          if ((ELEMENTTYPE == TETRAOLD).OR.(ELEMENTTYPE == TETRA4).OR.(ELEMENTTYPE == TETRA10)) then  
                            LocalNodeID = DetermineSideNodesTetrahedronHOE(ISide, INode3)
                          else
                            call GiveError('Absorbing boundary not implemented for element type '//trim(String(ELEMENTTYPE))) 
                          end if
                          
                        else if ( NDIM == 2 ) then
                          if ((ELEMENTTYPE == TRI3).OR.(ELEMENTTYPE == TRI6)) then
                            LocalNodeID = DetermineSideNodesTRI6(ISide, INode3)
                          !else if ((ELEMENTTYPE == QUAD4).OR.(ELEMENTTYPE == QUAD8)) then
                          !  LocalNodeID = DetermineSideNodesQUAD8(ISide, INode3)
                          else
                            call GiveError('Absorbing boundary not implemented for element type '//trim(String(ELEMENTTYPE))) 
                          end if
                        end if  
                        
                        DummyNode(INode3) = ElementConnectivities(LocalNodeID, ElementID)
                        if ( NodeCndVisSld2(DummyNode(INode3), 1) == 0 ) Belongs2AbsorbingBoundary = .false. ! check if node lies on defined beoundary
                      end do 
                      if (Belongs2AbsorbingBoundary) then
                          
                        if ( NDIM == 2 ) then
                          ! connecting vector between the two nodes spanning the absorbing boundary line
                          Connector = NodalOriginalCoord(DummyNode(2), :) - NodalOriginalCoord(DummyNode(1), :)  
                          ! calculate length of absorbing boundary line 
                          ElmSurfaceArea(ElementID) = Length(Connector, 2)
                          ! calculate contribution to each adjacent node
                          VisNodSurAraSld(INode) = VisNodSurAraSld(INode) + ( ElmSurfaceArea(ElementID) / 2.0 )
                          ! determine if absorbing boundary is in normal or shear direction
                          do I = 1, NVECTOR
                            if ( NodeCndVisSld2(INode, 2+3*(I-1)) /= 0 ) then ! coordinate direction is activated in .GOM
                              if ( Connector(I) == 0 ) then 
                                NodeCndVisSld2(INode, 2+3*(I-1)) = 1 ! normal AB
                              else 
                                NodeCndVisSld2(INode, 2+3*(I-1)) = 2 ! shear AB
                              end if    
                            end if 
                          end do
                        end if  
                        
                        if ( NDIM == 3 )  then  
                          do IDim = 1, NVECTOR ! create connecting vector between the 3 nodes
                            Vec1(IDim) = NodalOriginalCoord(DummyNode(2), IDim) - NodalOriginalCoord(DummyNode(1), IDim)
                            Vec2(IDim) = NodalOriginalCoord(DummyNode(3), IDim) - NodalOriginalCoord(DummyNode(1), IDim)                    
                          end do
                          Normal = CrossProduct(Vec1, Vec2) ! calculate the normal and determine its direction
                          ScalarProduct(1) = DotProduct((/1.0,0.0,0.0/),Normal,3) 
                          ScalarProduct(2) = DotProduct((/0.0,1.0,0.0/),Normal,3) 
                          ScalarProduct(3) = DotProduct((/0.0,0.0,1.0/),Normal,3)
                          ! calculate area of absorbing boundary surface 
                          ElmSurfaceArea(ElementID) = Length(Normal, 3) / 2.0
                          ! calculate contribution to each adjacent node
                          VisNodSurAraSld(INode) = VisNodSurAraSld(INode) + ( ElmSurfaceArea (ElementID) / 3.0 )
                          do I = 1, NVECTOR
                            if ( NodeCndVisSld2(INode, 2+3*(I-1)) /= 0 ) then ! coordinate direction is activated in .GOM
                              if ( ScalarProduct(I) == 0 ) then
                                NodeCndVisSld2(INode, 2+3*(I-1)) = 2 ! shear AB
                              else
                                NodeCndVisSld2(INode, 2+3*(I-1)) = 1 ! normal AB
                              end if
                            end if 
                          end do
                        end if 
                         
                      end if    
                    end if                            
                  end do    
                end do
              end if
            end do
            
        end subroutine CalculateNodeSurfaceAreaandSetCondSolid


        subroutine CalculateNodeSurfaceAreaandSetCondLiquid()
        !**********************************************************************
        !
        !   Function:  To calculate the surface area associated with each node
        !
        ! OUT: populates VisNodSurAraWat(INode), i.e. the surface area of absorbing boundary nodes  
        ! OUT: re-populates NodeCndVisWat2(INode,J), i.e. if absorbing boundary is in normal (1) or shear (2) direction
        !
        !**********************************************************************

        implicit none
        
          ! local variables 
          real(REAL_TYPE), dimension(Counters%NEl) :: ElmSurfaceArea   
          real(REAL_TYPE), dimension(NVECTOR):: Vec1, Vec2, Normal, ScalarProduct, Connector
          integer(INTEGER_TYPE) :: INode, INode2, INode3, ISide, IDim, I, NumberOfElementsNode
          integer(INTEGER_TYPE) :: IsMeshBoundary, ElementID, LocalNodeID
          integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: DummyNode
          logical :: Belongs2AbsorbingBoundary
                                  
          VisNodSurAraWat = 0.0
          ElmSurfaceArea = 0.0   
                 
            do INode = 1, Counters%Nodtot                                ! loop all nodes 
              if ( NodeCndVisWat2(INode, 1) > 0 ) then                   ! check if node is defined as viscous node in .GOM file
                NumberOfElementsNode = ElementsPerNode(INode)            ! get number of elements around this considered node
                do INode2 = 1, NumberOfElementsNode                      ! check all elements sharing this node
                  ElementID = NodeElementAdjacencies(NodeElementOffsets(INode)+INode2) ! get the element ID of each neighbouring element
                  do ISide = 1, ELEMENTSIDES                             ! loop around all element sides        
                    IsMeshBoundary = BoundaryElementSurface(ElementID, ISide, IsActiveElement, Counters%NEl) ! check if elementside is a mesh boundary
                    if ( IsMeshBoundary == -1 ) then                     ! check which side of the element is part of the boundary
                      Belongs2AbsorbingBoundary = .true.
                      do INode3 = 1, ELEMENTBOUNDARYNODES                ! get the node IDs of this side
                        if ( NDIM == 3 ) LocalNodeID = DetermineSideNodesTetrahedronHOE(ISide, INode3)
                        if ( NDIM == 2 ) LocalNodeID = DetermineSideNodesTRI6(ISide, INode3)  
                        DummyNode(INode3) = ElementConnectivities(LocalNodeID, ElementID)
                        if ( NodeCndVisWat2(DummyNode(INode3), 1) == 0 ) Belongs2AbsorbingBoundary = .false. ! check if node lies on defined boundary
                      end do 
                      
                      if (Belongs2AbsorbingBoundary) then
                        
                        if ( NDIM == 2 ) then
                          ! connecting vector between the two nodes spanning the absorbing boundary line
                          Connector = NodalOriginalCoord(DummyNode(2), :) - NodalOriginalCoord(DummyNode(1), :)  
                          ! calculate length of absorbing boundary line 
                          ElmSurfaceArea(ElementID) = Length(Connector, 2)
                          ! calculate contribution to each adjacent node
                          VisNodSurAraWat(INode) = VisNodSurAraWat(INode) + ( ElmSurfaceArea(ElementID) / 2.0 )
                          ! determine if absorbing boundary is in normal or shear direction
                          do I = 1, 2
                            if ( NodeCndVisWat2(INode, 2+3*(I-1)) /= 0 ) then ! coordinate direction is activated in .GOM
                              if ( Connector(I) == 0 ) then 
                                NodeCndVisWat2(INode, 2+3*(I-1)) = 1 ! normal AB
                              else 
                                NodeCndVisWat2(INode, 2+3*(I-1)) = 2 ! shear AB
                              end if    
                            end if 
                          end do
                        end if    
                          
                        if ( NDIM == 3 ) then    
                          do IDim = 1, NVECTOR ! create connecting vector between the 3 nodes
                            Vec1(IDim) = NodalOriginalCoord(DummyNode(2), IDim) - NodalOriginalCoord(DummyNode(1), IDim)
                            Vec2(IDim) = NodalOriginalCoord(DummyNode(3), IDim) - NodalOriginalCoord(DummyNode(1), IDim)                    
                          end do
                          Normal = CrossProduct(Vec1, Vec2) ! calculate the normal and determine its direction
                          ScalarProduct(1) = DotProduct((/1.0,0.0,0.0/),Normal,3)
                          ScalarProduct(2) = DotProduct((/0.0,1.0,0.0/),Normal,3)
                          ScalarProduct(3) = DotProduct((/0.0,0.0,1.0/),Normal,3)
                          ! calculate area of absorbing boundary surface 
                          ElmSurfaceArea(ElementID) = Length(Normal, 3) / 2.0
                          ! calculate contribution to each adjacent node
                          VisNodSurAraWat(INode) = VisNodSurAraWat(INode) + ( ElmSurfaceArea(ElementID) / 3.0 )
                          do I = 1, NVECTOR
                            if ( NodeCndVisWat2(INode, 2+3*(I-1)) /= 0 ) then ! coordinate direction is activated in .GOM
                              if ( ScalarProduct(I) == 0 ) then
                                NodeCndVisWat2(INode, 2+3*(I-1)) = 2 ! shear AB
                              else
                                NodeCndVisWat2(INode, 2+3*(I-1)) = 1 ! normal AB
                              end if
                            end if 
                          end do
                        end if    
                        
                      end if
                    end if                            
                  end do    
                end do
              end if
            end do
                   
        end subroutine CalculateNodeSurfaceAreaandSetCondLiquid
         

        subroutine CalculateNodeSurfaceAreaandSetCondGas()
        !**********************************************************************
        !
        !   Function:  To calculate the surface area associated with each node
        !
        ! OUT: populates VisNodSurAraGas(INode), i.e. the surface area of absorbing boundary nodes  
        ! OUT: re-populates NodeCndVisGas2(INode,J), i.e. if absorbing boundary is in normal (1) or shear (2) direction
        !
        !**********************************************************************

        implicit none
        
          ! local variables 
          real(REAL_TYPE), dimension(Counters%NEl) :: ElmSurfaceArea   
          real(REAL_TYPE), dimension(NVECTOR):: Vec1, Vec2, Normal, ScalarProduct, Connector
          integer(INTEGER_TYPE) :: INode, INode2, INode3, ISide, IDim, I, NumberOfElementsNode
          integer(INTEGER_TYPE) :: IsMeshBoundary, ElementID, LocalNodeID
          integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: DummyNode
          logical :: Belongs2AbsorbingBoundary
                                  
          VisNodSurAraGas = 0.0
          ElmSurfaceArea = 0.0   
                 
            do INode = 1, Counters%Nodtot                                ! loop all nodes 
              if ( NodeCndVisGas2(INode, 1) > 0 ) then                   ! check if node is defined as viscous node in .GOM file
                NumberOfElementsNode = ElementsPerNode(INode)            ! get number of elements around this considered node
                do INode2 = 1, NumberOfElementsNode                      ! check all elements sharing this node
                  ElementID = NodeElementAdjacencies(NodeElementOffsets(INode)+INode2) ! get the element ID of each neighbouring element
                  do ISide = 1, ELEMENTSIDES                             ! loop around all element sides        
                    IsMeshBoundary = BoundaryElementSurface(ElementID, ISide, IsActiveElement, Counters%NEl) ! check if elementside is a mesh boundary
                    if ( IsMeshBoundary == -1 ) then                     ! check which side of the element is part of the boundary
                      Belongs2AbsorbingBoundary = .true.
                      do INode3 = 1, ELEMENTBOUNDARYNODES                ! get the node IDs of this side
                        if ( NDIM == 3 ) LocalNodeID = DetermineSideNodesTetrahedronHOE(ISide, INode3)
                        if ( NDIM == 2 ) LocalNodeID = DetermineSideNodesTRI6(ISide, INode3)  
                        DummyNode(INode3) = ElementConnectivities(LocalNodeID, ElementID)
                        if ( NodeCndVisGas2(DummyNode(INode3), 1) == 0 ) Belongs2AbsorbingBoundary = .false. ! check if node lies on defined boundary
                      end do 
                      
                      if (Belongs2AbsorbingBoundary) then
                        
                        if ( NDIM == 2 ) then
                          ! connecting vector between the two nodes spanning the absorbing boundary line
                          Connector = NodalOriginalCoord(DummyNode(2), :) - NodalOriginalCoord(DummyNode(1), :)  
                          ! calculate length of absorbing boundary line 
                          ElmSurfaceArea(ElementID) = Length(Connector, 2)
                          ! calculate contribution to each adjacent node
                          VisNodSurAraGas(INode) = VisNodSurAraGas(INode) + ( ElmSurfaceArea(ElementID) / 2.0 )
                          ! determine if absorbing boundary is in normal or shear direction
                          do I = 1, 2
                            if ( NodeCndVisGas2(INode, 2+3*(I-1)) /= 0 ) then ! coordinate direction is activated in .GOM
                              if ( Connector(I) == 0 ) then 
                                NodeCndVisGas2(INode, 2+3*(I-1)) = 1 ! normal AB
                              else 
                                NodeCndVisGas2(INode, 2+3*(I-1)) = 2 ! shear AB
                              end if    
                            end if 
                          end do
                        end if    
                          
                        if ( NDIM == 3 ) then    
                          do IDim = 1, NVECTOR ! create connecting vector between the 3 nodes
                            Vec1(IDim) = NodalOriginalCoord(DummyNode(2), IDim) - NodalOriginalCoord(DummyNode(1), IDim)
                            Vec2(IDim) = NodalOriginalCoord(DummyNode(3), IDim) - NodalOriginalCoord(DummyNode(1), IDim)                    
                          end do
                          Normal = CrossProduct(Vec1, Vec2) ! calculate the normal and determine its direction
                          ScalarProduct(1) = DotProduct((/1.0,0.0,0.0/),Normal,3)
                          ScalarProduct(2) = DotProduct((/0.0,1.0,0.0/),Normal,3)
                          ScalarProduct(3) = DotProduct((/0.0,0.0,1.0/),Normal,3)
                          ! calculate area of absorbing boundary surface 
                          ElmSurfaceArea(ElementID) = Length(Normal, 3) / 2.0
                          ! calculate contribution to each adjacent node
                          VisNodSurAraGas(INode) = VisNodSurAraGas(INode) + ( ElmSurfaceArea(ElementID) / 3.0 )
                          do I = 1, NVECTOR
                            if ( NodeCndVisGas2(INode, 2+3*(I-1)) /= 0 ) then ! coordinate direction is activated in .GOM
                              if ( ScalarProduct(I) == 0 ) then
                                NodeCndVisGas2(INode, 2+3*(I-1)) = 2 ! shear AB
                              else
                                NodeCndVisGas2(INode, 2+3*(I-1)) = 1 ! normal AB
                              end if
                            end if 
                          end do
                        end if    
                        
                      end if
                    end if                            
                  end do    
                end do
              end if
            end do
                   
        end subroutine CalculateNodeSurfaceAreaandSetCondGas
        
        
        subroutine FormDashpotandSpringSolid()
        !**********************************************************************
        !
        !   Function:  To calculate the dashpot terms
        !
        ! OUT: DashpotSld:
        ! OUT: SpringSld:
        ! OUT: DoFCndVisSld:
        !        
        !**********************************************************************

        implicit none
        
          ! local variables 
          integer(INTEGER_TYPE) :: INode, IDim, Nix, ICondition, MSet
          real(REAL_TYPE) :: Vp, Vs, rho, A, Eoed, G, Nu, Kw, n
          logical:: IsUndrEffectiveStress
 
          DashpotSld = 0.0
          SpringSld = 0.0
          DofCndVisSld = 0
           
          MSet = CalParams%AbsorbingBoundaries%VBMaterialSet
          rho = MatParams(MSet)%DryWeight / CalParams%GravityData%GAccel ! Mg 
          G = MatParams(MSet)%ShearModulus ! kPa
          Nu = MatParams(MSet)%PoissonRatio
          Kw = MatParams(MSet)%BulkModulusLiquid
          n = MatParams(MSet)%InitialPorosity
          Eoed = ( 2*G*(1-Nu) ) / (1-2*Nu) ! kPa
          IsUndrEffectiveStress = &
            ! code version 2016 and previous
            ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MSet)%MaterialType)=='2-phase')) .or. &
            ! code version 2017.1 and following
            (trim(MatParams(MSet)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
            
          if ( IsUndrEffectivestress .and. ( n > 0.0 ) ) then
            Eoed = Eoed + Kw / n
          end if
            
          Vp = sqrt( Eoed / rho ) ! m/s       
          Vs = sqrt( G / rho ) ! m/s
            
          do INode = 1, Counters%NodTot ! over nodes
              
            A = VisNodSurAraSld(INode) ! surface area
            Nix = ReducedDof(INode) + 1 ! global dof x

            do IDim = 1, NVECTOR ! over coordinates
              ICondition = NodeCndVisSld2(INode, 2+3*(IDim-1))
              DofCndVisSld(Nix+IDim-1) = ICondition
                
              if ( ICondition > 0) then ! if it is 0 then no absorbing boundaries are added

                if ( ICondition == 1) then ! primary
                  DashpotSld(Nix+IDim-1) = rho*Vp*A * NodeCndVisSld2(INode, 3+3*(IDim-1))
                  SpringSld(Nix+IDim-1) = rho*Vp*Vp*A / NodeCndVisSld2(INode, 4+3*(IDim-1))
                else if ( ICondition == 2 ) then ! secondary
                  DashpotSld(Nix+IDim-1) = rho*Vs*A * NodeCndVisSld2(INode, 3+3*(IDim-1)) 
                  SpringSld(Nix+IDim-1) = rho*Vs*Vs*A / NodeCndVisSld2(INode, 4+3*(IDim-1))
                end if
                  
              end if
                
            end do
              
          end do ! over nodes

        end subroutine FormDashpotandSpringSolid


        subroutine FormDashpotandSpringLiquid()
        !**********************************************************************
        !
        !   Function:  To calculate the dashpot terms
        !        
        ! OUT: DashpotWat:
        ! OUT: SpringWat:
        ! OUT: DoFCndVisWat:
        !
        !**********************************************************************

        implicit none
        
          ! local variables 
          integer(INTEGER_TYPE) :: INode, IDim, Nix, ICondition, MSet
          real(REAL_TYPE) :: Vp, rho, A, Kw
                                  
          DashpotWat = 0.0
          SpringWat = 0.0
          DofCndVisWat = 0
           
          MSet = CalParams%AbsorbingBoundaries%VBMaterialSet
          rho = MatParams(MSet)%WeightLiquid / CalParams%GravityData%GAccel ! Mg  
          Kw = MatParams(MSet)%BulkModulusLiquid ! kPa
          if ( rho == 0.0 ) then ! when the object is dry for example pile
            Vp = 0.0
          else
            Vp = sqrt( Kw / rho ) ! m/s
          end if  
            
          do INode = 1, Counters%NodTot ! over nodes
            A = VisNodSurAraWat(INode) ! surface area
            Nix = ReducedDof(INode) + 1  ! global dof x

            do IDim = 1, NVECTOR ! over coordinates
              ICondition = NodeCndVisWat2(INode,2+3*(IDim-1))
              DofCndVisWat(Nix+IDim-1) = ICondition
                
              if ( ICondition > 0 ) then ! if it is 0 then no viscous boundaries are added
                if ( ICondition == 1 ) then ! primary
                  DashpotWat(Nix+IDim-1) = rho*Vp*A * NodeCndVisWat2(INode, 3+3*(IDim-1))
                  SpringWat(Nix+IDim-1) = rho*Vp*Vp*A / NodeCndVisWat2(INode, 4+3*(IDim-1))
                else if ( ICondition == 2 ) then ! secondary
                  DashpotWat(Nix+IDim-1) = rho*Vp*A * NodeCndVisWat2(INode, 3+3*(IDim-1))
                  SpringWat(Nix+IDim-1) = rho*Vp*Vp*A / NodeCndVisWat2(INode, 4+3*(IDim-1))
                end if
              end if
                
            end do
            
          end do ! over nodes
          
        end subroutine FormDashpotandSpringLiquid

        
        subroutine FormDashpotandSpringGas()
        !**********************************************************************
        !
        !   Function:  To calculate the dashpot terms
        !        
        ! OUT: DashpotGas:
        ! OUT: SpringGas:
        ! OUT: DoFCndVisGas:
        !
        !**********************************************************************

        implicit none
        
          ! local variables 
          integer(INTEGER_TYPE) :: INode, IDim, Nix, ICondition, MSet
          real(REAL_TYPE) :: Vp, rho, A, Kg
                                  
          DashpotGas = 0.0
          SpringGas = 0.0
          DofCndVisGas = 0
           
          MSet = CalParams%AbsorbingBoundaries%VBMaterialSet
          rho = MatParams(MSet)%WeightGas / CalParams%GravityData%GAccel ! Mg  
          Kg = MatParams(MSet)%BulkModulusGas ! kPa
          if ( rho == 0.0 ) then ! when the object is non-porous for example pile
            Vp = 0.0
          else
            Vp = sqrt( Kg / rho ) ! m/s
          end if  
            
          do INode = 1, Counters%NodTot ! over nodes
            A = VisNodSurAraGas(INode) ! surface area
            Nix = ReducedDof(INode) + 1  ! global dof x

            do IDim = 1, NVECTOR ! over coordinates
              ICondition = NodeCndVisGas2(INode,2+3*(IDim-1))
              DofCndVisGas(Nix+IDim-1) = ICondition
                
              if ( ICondition > 0 ) then ! if it is 0 then no viscous boundaries are added
                if ( ICondition == 1 ) then ! primary
                  DashpotGas(Nix+IDim-1) = rho*Vp*A * NodeCndVisGas2(INode, 3+3*(IDim-1))
                  SpringGas(Nix+IDim-1) = rho*Vp*Vp*A / NodeCndVisGas2(INode, 4+3*(IDim-1))
                else if ( ICondition == 2 ) then ! secondary
                  DashpotGas(Nix+IDim-1) = rho*Vp*A * NodeCndVisGas2(INode, 3+3*(IDim-1))
                  SpringGas(Nix+IDim-1) = rho*Vp*Vp*A / NodeCndVisGas2(INode, 4+3*(IDim-1))
                end if
              end if
                
            end do
            
          end do ! over nodes
          
        end subroutine FormDashpotandSpringGas
        
        
                

        subroutine GetNodalAccumulatedDisplacementsSolid()
        !**********************************************************************
        !
        !    Function:  To map displacement from particles to the grid points (solid)
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE), dimension(Counters%N) :: ActiveDof
          integer(INTEGER_TYPE) :: IDof
          
          call SetActiveDof(ActiveDof)

          do IDof = 1, Counters%N
            if (DashpotSld(IDof)/=0) then  ! there is VB here
              if (ActiveDof(IDof)==1) then
                NodalTotDisplacementVB(IDof, 1) = NodalTotDisplacementVB(IDof, 1) + IncrementalDisplacementSoil(IDof, 1)
              else
                NodalTotDisplacementVB(IDof, 1) = 0.0
              end if
            end if
          end do
         
        end subroutine GetNodalAccumulatedDisplacementsSolid 
        
         
        subroutine GetNodalAccumulatedDisplacementsWater(NodalIncDisplacementW)
        !**********************************************************************
        !
        !    Function:  To map displacement from particles to the grid points (water)
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent (in) :: NodalIncDisplacementW
        
          ! Local variables
          integer(INTEGER_TYPE), dimension(Counters%N) :: ActiveDof
          integer(INTEGER_TYPE) :: IDof
          
          call SetActiveDof(ActiveDof)

          do IDof = 1, Counters%N
            if (DashpotWat(IDof)/=0) then  ! there is VB here
              if (ActiveDof(IDof)==1) then
                NodalTotDisplacementVBWater(IDof, 1) = NodalTotDisplacementVBWater(IDof, 1) + NodalIncDisplacementW(IDof, 1)
              else
                NodalTotDisplacementVBWater(IDof, 1) = 0.0
              end if
            end if
          end do
         
        end subroutine GetNodalAccumulatedDisplacementsWater 
        
  
        subroutine GetNodalAccumulatedDisplacementsGas(NodalIncDisplacementG)
        !**********************************************************************
        !
        !    Function:  To map displacement from particles to the grid points (gas)
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent (in) :: NodalIncDisplacementG
        
          ! Local variables
          integer(INTEGER_TYPE), dimension(Counters%N) :: ActiveDof
          integer(INTEGER_TYPE) :: IDof
          
          call SetActiveDof(ActiveDof)

          do IDof = 1, Counters%N
            if (DashpotGas(IDof)/=0) then  ! there is VB here
              if (ActiveDof(IDof)==1) then
                NodalTotDisplacementVBGas(IDof, 1) = NodalTotDisplacementVBGas(IDof, 1) + NodalIncDisplacementG(IDof, 1)
              else
                NodalTotDisplacementVBGas(IDof, 1) = 0.0
              end if
            end if
          end do
         
        end subroutine GetNodalAccumulatedDisplacementsGas       
  
        
        subroutine CalcParticleDisplacementsWater()
        !**********************************************************************
        !   
        !   Function: Calculate displacements of the water at the material points
        !
        !**********************************************************************
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity, IParticle, ElementID
          integer(INTEGER_TYPE) :: NodeID, DofID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleWaterDisplacement

          if (CalParams%NumberOfPhases/=2) RETURN
          if (.not.CalParams%ApplyAbsorbingBoundary) RETURN
          
          do IParticle = 1, Counters%NParticles ! Loop over particles

            ElementID = ElementIDArray(IParticle)

            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle) 
            else
              iEntity = 1
            end if

            ParticleWaterDisplacement = 0.0

            do IDof = 1, NVECTOR
              do INode = 1, ELEMENTNODES  
                NodeID = iabs( ElementConnectivities(INode, ElementID) ) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleWaterDisplacement(IDof) = ParticleWaterDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * IncrementalDisplacementWater(DofID,iEntity)
              end do
            end do

            ! Incremental particle displacement
            Particles(IParticle)%UStepWater = ParticleWaterDisplacement
            
          end do

        end subroutine CalcParticleDisplacementsWater
        
         
        subroutine GetIncrementalDisplacementGas(NodalIncDisplacementG)
        !**********************************************************************
        !
        !    Function:  To get the incremental displacements of the Gas
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent (in) :: NodalIncDisplacementG
        
          ! Local variables
          integer(INTEGER_TYPE) :: IDof, IEntity
          
          do IEntity = 1, Counters%NEntity
            do IDof = 1, Counters%N
              IncrementalDisplacementGas(IDof, IEntity) = NodalIncDisplacementG(IDof, IEntity)
            end do
          end do
   
        end subroutine GetIncrementalDisplacementGas        
        

        subroutine SetActiveDof(ActiveDof)
        !**********************************************************************
        !
        !    Function:  set active degrees of freedom (Dof)
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), dimension(Counters%N), intent(inout) :: ActiveDof
          ! Local variables
          integer(INTEGER_TYPE) :: I, IAEl, IEl, INode, NodeID, IDof
          
          ActiveDof = 0
         
          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            do INode = 1, ELEMENTNODES
              NodeID = ElementConnectivities(INode, IEl) ! Global node ID
              IDof = ReducedDof(NodeID)  ! Global storage coordinate of x-val
              do I = 1, NVECTOR
                ActiveDof(IDof+I) = 1
              end do
           end do
         end do

        end subroutine SetActiveDof    

        subroutine GetVisDampingForceSld()
        !**********************************************************************
        !
        !   Function:  To calculate the forces of the viscous boundary
        !        
        !**********************************************************************

        implicit none
       
           !local variables 
           integer(INTEGER_TYPE), dimension(Counters%N) :: ActiveDof
           integer(INTEGER_TYPE) :: IEntity, IDof
           real(REAL_TYPE) :: DeltaT
           real(REAL_TYPE), dimension(Counters%N, Counters%nEntity):: IncrVisDampForceSld
           real(REAL_TYPE), dimension(Counters%N, Counters%nEntity):: CorrectedNodalAcce
           real(REAL_TYPE), dimension(Counters%N, Counters%nEntity):: IncrementalDisplacementSoilVBSld


           if (IsMPMComputation()) then !MPM ...need to get the accumulated displacement at the nodes
             call CorrectIncAccelerationAndIncDisp(CorrectedNodalAcce, IncrementalDisplacementSoilVBSld)
           end if

           call SetActiveDof(ActiveDof)

           DeltaT =  CalParams%TimeIncrement
           IncrVisDampForceSld= 0.0
           
            do IDof = 1, Counters%N ! over dof
              do IEntity = 1, Counters%NEntity
                if (IsMPMComputation()) then !MPM
                  if (ActiveDof(IDof)==1) then
                    if (CalParams%AbsorbingBoundaries%ApplyNodeData) then
                      ! Displacement and velocity at nodes are used,
                      ! boundary conditions are the same as with FEM
                      IncrVisDampForceSld(IDof, IEntity) = SpringSld(IDof) * IncrementalDisplacementSoil(IDof, IEntity) + DashpotSld(IDof)*AccelerationSoil(IDof, IEntity) * DeltaT
                    else
                      ! With MPM displacement and acceleration are corrected
                      IncrVisDampForceSld(IDof, IEntity) = SpringSld(IDof) * IncrementalDisplacementSoilVBSld(IDof, IEntity) + DashpotSld(IDof)*CorrectedNodalAcce(IDof, IEntity) * DeltaT
                    end if
                  else
                    VisDampForceSld(IDof, IEntity) = 0.0
 
                  end if
                else ! FEM
                    IncrVisDampForceSld(IDof, IEntity) = SpringSld(IDof) * IncrementalDisplacementSoil(IDof, IEntity) + DashpotSld(IDof)*AccelerationSoil(IDof, IEntity) * DeltaT
                end if

              end do ! entity
            end do ! over dof
            
             if (IS3DCYLINDRIC) then ! rotation is needed
              do IEntity = 1, Counters%nEntity
                call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, IncrVisDampForceSld(:, IEntity), IncrVisDampForceSld(:, IEntity))
              end do
            end if ! rotation

            VisDampForceSld= VisDampForceSld + IncrVisDampForceSld

        end subroutine GetVisDampingForceSld

         
        subroutine CorrectIncAccelerationAndIncDispWater(CorrectedNodalAcce, IncrementalDisplacementWaterVBSld)
        !**********************************************************************
        !
        !    Function:  Get the acceleration and incremental displacement
        !
        !**********************************************************************

        implicit none
          
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: CorrectedNodalAcce
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: IncrementalDisplacementWaterVBSld
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn
          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: IDof
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration, ParticleDisplacement
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
      
          NoEn = Counters%nEntity
          IncrementalDisplacementWaterVBSld = 0.0
          CorrectedNodalAcce = 0.0

          do IAEl = 1, Counters%NAEl ! Loop over all elements
            IEl = ActiveElement(IAEl)
            do I = 1, NVECTOR
              IDof(I, 1:ELEMENTNODES) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
            end do
            do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              ParticleShape = ShapeValuesArray(ParticleIndex,:)
              ParticleAcceleration = 0.0
              if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex) 
              else
                iEntity = 1
              end if 
              
              do I = 1, NVECTOR 
                ParticleAcceleration(I) = Particles(ParticleIndex)%AccelerationWater(I)
                ParticleDisplacement(I) = Particles(ParticleIndex)%UStepWater(I)
                CorrectedNodalAcce(IDof(I, 1:ELEMENTNODES), iEntity) = CorrectedNodalAcce(IDof(I, 1:ELEMENTNODES), iEntity) + ParticleShape*ParticleAcceleration(I)
                IncrementalDisplacementWaterVBSld(IDof(I, 1:ELEMENTNODES), iEntity) = IncrementalDisplacementWaterVBSld(IDof(I, 1:ELEMENTNODES), iEntity) + ParticleShape*ParticleDisplacement(I)
              end do  
 
            end do ! loop over material points
            
          end do ! loop over elements
          
        end subroutine CorrectIncAccelerationAndIncDispWater

     
        subroutine CorrectIncAccelerationAndIncDisp(CorrectedNodalAcce, IncrementalDisplacementSoilVBSld)
        !**********************************************************************
        !
        !    Function:  Get the acceleration and incremental displacement
        !
        !**********************************************************************

        implicit none
          
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: CorrectedNodalAcce
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: IncrementalDisplacementSoilVBSld
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn
          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: IDof
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration,  ParticleDisplacement
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
      
          NoEn = Counters%nEntity
          IncrementalDisplacementSoilVBSld = 0.0
          CorrectedNodalAcce = 0.0

          do IAEl = 1, Counters%NAEl  ! loop over all active elements
            IEl = ActiveElement(IAEl)
            do I = 1, NVECTOR
              IDof(I, 1:ELEMENTNODES) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
            end do
            
            do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              ParticleShape = ShapeValuesArray(ParticleIndex,:)
              ParticleAcceleration = 0.0
              if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex) 
              else
                iEntity = 1
              end if 
              
              do I = 1, NVECTOR 
                ParticleAcceleration(I) = AccelerationArray(ParticleIndex, I)
              end do  
     
              ParticleDisplacement(:) = UStepArray(ParticleIndex,:)
              
              do I = 1, NVECTOR
                CorrectedNodalAcce(IDof(I,1:ELEMENTNODES), iEntity) = CorrectedNodalAcce(IDof(I,1:ELEMENTNODES), iEntity) + ParticleShape*ParticleAcceleration(I)
                IncrementalDisplacementSoilVBSld(IDof(I,1:ELEMENTNODES), iEntity) = IncrementalDisplacementSoilVBSld(IDof(I,1:ELEMENTNODES), iEntity) + ParticleShape*ParticleDisplacement(I)
              end do 

            end do ! loop material points
            
          end do ! loop elements    
          
        end subroutine CorrectIncAccelerationAndIncDisp

        
        subroutine GetVisDampingForceWat(AccelerationWater)
        !**********************************************************************
        !
        !   Function:  To calculate the forces of the viscous boundary
        !        
        !**********************************************************************

        implicit none
       
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent (in) :: AccelerationWater
          ! local variables 
          integer(INTEGER_TYPE), dimension(Counters%N) :: ActiveDof
          integer(INTEGER_TYPE) :: IEntity, IDof
          real(REAL_TYPE) :: DeltaT
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: IncrVisDampForceWat
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: CorrectedNodalAcce
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: IncrementalDisplacementWaterVBSld
     
          call SetActiveDof(ActiveDof)

          if (IsMPMComputation().and.CalParams%ApplyMeshSmoothing) then !MPM ...need to get the accumulated displacement at the nodes
            call CorrectIncAccelerationAndIncDispWater(CorrectedNodalAcce, IncrementalDisplacementWaterVBSld)
          end if

          DeltaT =  CalParams%TimeIncrement
          IncrVisDampForceWat= 0.0
           
            do IDof = 1, Counters%N ! over dof
              do IEntity = 1, Counters%NEntity
                if (IsMPMComputation()) then !MPM
                  if (ActiveDof(IDof)==1) then
                    if (CalParams%ApplyMeshSmoothing) then
                      IncrVisDampForceWat(IDof, IEntity) = SpringWat(IDof) * IncrementalDisplacementWaterVBSld(IDof, IEntity) + DashpotWat(IDof)*CorrectedNodalAcce(IDof, IEntity) * DeltaT
                    else
                      IncrVisDampForceWat(IDof, IEntity) = SpringWat(IDof) * IncrementalDisplacementWater(IDof, IEntity) + DashpotWat(IDof)*AccelerationWater(IDof, IEntity) * DeltaT
                    end if
                  else
                    VisDampForceWat(IDof, IEntity) = 0.0
 
                  end if
                else ! FEM
                    IncrVisDampForceWat(IDof, IEntity) = SpringWat(IDof) * IncrementalDisplacementWater(IDof, IEntity) + DashpotWat(IDof)*AccelerationWater(IDof, IEntity) * DeltaT
                end if

              end do ! entity
            end do ! over dof

             if (IS3DCYLINDRIC) then ! rotation is needed
              do IEntity = 1, Counters%nEntity
                call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, IncrVisDampForceWat(:, IEntity), IncrVisDampForceWat(:, IEntity))
              end do
            end if ! rotation

            VisDampForceWat= VisDampForceWat + IncrVisDampForceWat

        end subroutine GetVisDampingForceWat
        

        subroutine GetVisDampingForceGas(AccelerationGas)
        !**********************************************************************
        !
        !   Function:  To calculate the forces of the viscous boundary
        !        
        !**********************************************************************

        implicit none
       
           real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent (in) :: AccelerationGas
           !local variables 
           integer(INTEGER_TYPE), dimension(Counters%N) :: ActiveDof
           integer(INTEGER_TYPE) :: IEntity, IDof
           real(REAL_TYPE) :: DeltaT
           real(REAL_TYPE), dimension(Counters%N, Counters%nEntity):: IncrVisDampForceGas
           
           call SetActiveDof(ActiveDof)

           DeltaT = CalParams%TimeIncrement
           IncrVisDampForceGas= 0.0

            do IDof = 1, Counters%N ! over dof
              do IEntity = 1, Counters%NEntity
                if (IsMPMComputation()) then !MPM
                  if (ActiveDof(IDof)==1) then
                    IncrVisDampForceGas(IDof, IEntity) = SpringGas(IDof) * IncrementalDisplacementGas(IDof, IEntity) + DashpotGas(IDof)*AccelerationGas(IDof, IEntity) * DeltaT
                  else
                    VisDampForceGas(IDof, IEntity) = 0.0
                  end if
                else ! FEM
                    IncrVisDampForceGas(IDof, IEntity) = SpringGas(IDof) * IncrementalDisplacementGas(IDof, IEntity) + DashpotGas(IDof)*AccelerationGas(IDof, IEntity) * DeltaT
                end if

              end do ! entity
            end do   ! over dof

             if (IS3DCYLINDRIC) then ! rotation is needed
              do IEntity = 1, Counters%nEntity 
                call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, IncrVisDampForceGas(:, IEntity), IncrVisDampForceGas(:, IEntity))
              end do
            end if ! rotation

            VisDampForceGas= VisDampForceGas + IncrVisDampForceGas

        end subroutine GetVisDampingForceGas
        
        subroutine InitialiseAbsorbingBoundariesForcesAndStiffness()
        !**********************************************************************
        !
        !   Function: Calculate the advective flux of water at the material point
        !             Only if ApplyAbsorbingBoundary
        !   Note : This subroutine works only if NumberOfLayers = 1
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          real(REAL_TYPE), dimension(Counters%NodTot) :: NodalVerticalStressesSoil, NodalPressuresWater
          
          if (.not.CalParams%ApplyAbsorbingBoundary) RETURN
          if (NFORMULATION==2) RETURN 
          
          if ((CalParams%ApplyK0Procedure.and.(.not.IsFollowUpPhase())).or.(CalParams%AbsorbingBoundaries%ApplyVariableStiffAB)) then     
            call ComputeNodalVerticalStressesForK0(NodalOriginalCoord, NodalVerticalStressesSoil, NodalPressuresWater)
          else
            NodalVerticalStressesSoil = 0.0
          end if
          
          if (CalParams%ApplyK0Procedure.and.(.not.IsFollowUpPhase())) then     
            call InitialiseAbsorbingBoundaryForcesForK0(NodalVerticalStressesSoil, NodalPressuresWater)
          end if
         
        end subroutine InitialiseAbsorbingBoundariesForcesAndStiffness
       

        subroutine InitialiseAbsorbingBoundaryForcesForK0(NodalVerticalStressesSoil, NodalPressuresWater)
        !************************************************************************************
        !
        !   Function:  Calculate the dashpot terms --> this is only for 3D, extend for 2D
        !
        !************************************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%Nodtot), intent(in) :: NodalVerticalStressesSoil, NodalPressuresWater
          ! Local variables
          integer(INTEGER_TYPE) :: INod, NiX, IConditionX, IConditionY, IConditionZ, IEntity
          real(REAL_TYPE) :: A, F, FX, FY, FZ, FWP, FWPX, FWPY, FWPZ, K0Value
 
          K0Value = MatParams(CalParams%AbsorbingBoundaries%VBMaterialSet)%K0Value
        
          do INod = 1, Counters%NodTot
          
            A = VisNodSurAraSld(INod) ! Surface area
            if (A>0.0) then ! Get stresses at this node
              F = NodalVerticalStressesSoil(INod) * A
              FWP = NodalPressuresWater(INod) * A
              NiX = ReducedDof(INod) + 1
              IConditionX = NodeCndVisSld(INod, 1)
              IConditionY = NodeCndVisSld(INod, 2)
              IConditionZ = NodeCndVisSld(INod, 3)
              if (IConditionX==1) then ! Normal direction ... give force 
                FX = F * K0Value
                FWPX = FWP
              else
                FX = 0.0
                FWPX = 0.0
              end if
        
              if (IConditionY==1) then ! Normal direction ... give force 
                FY = F
                FWPY = FWP
              else
                FY = 0.0
                FWPY = 0.0
              end if
                 
              if (IConditionZ==1) then ! Normal direction ... give force 
                FZ = F * K0Value
                FWPZ = FWP
              else
                FZ = 0.0
                FWPZ = 0.0
              end if
                  
              do IEntity = 1, Counters%NEntity
                VisDampForceSld(NiX, IEntity) = -FX
                VisDampForceSld(NiX + 1, IEntity) = FY
                VisDampForceSld(NiX + 2, IEntity) = FZ
                VisDampForceWat(NiX, IEntity) = FWPX
                VisDampForceWat(NiX + 1, IEntity) = FWPY
                VisDampForceWat(NiX + 2, IEntity) = FWPZ
              end do
                
            end if
          end do

          if (IS3DCYLINDRIC) then
           do IEntity = 1, Counters%nEntity
             call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, VisDampForceSld(:, IEntity), VisDampForceSld(:, IEntity))
           end do
         end if

        end subroutine InitialiseAbsorbingBoundaryForcesForK0
        

        subroutine ComputeNodalVerticalStressesForK0(NodeCoord, NodalVerticalStressesSoil, NodalPressuresWater)
        !**********************************************************************
        !
        !    Function: Calculate nodal vertical stresses for K0 initialization
        !
        !**********************************************************************   
        implicit none
        
          real(REAL_TYPE), dimension(Counters%NodTot, NDIM), intent(in) :: NodeCoord
          real(REAL_TYPE), dimension(Counters%NodTot), intent(inout) :: NodalVerticalStressesSoil, NodalPressuresWater
          ! Local variables
          integer(INTEGER_TYPE) :: INode, MSet
          real(REAL_TYPE) :: K0Value, UpdatedFluidDensity, SigY, SigWP
            
          MSet = CalParams%AbsorbingBoundaries%VBMaterialSet
          K0Value = MatParams(MSet)%K0Value
          if (K0Value<SKIP_K0_THRESHOLD) RETURN
          
          do INode = 1, Counters%NodTot
            call ComputeVerticalStressForK0(MSet, NodeCoord(INode, 1:2), UpdatedFluidDensity, SigY, SigWP,SoilSurfaceNodeCoordMatrix, PhreaticSurfaceNodeCoordMatrix)
            NodalVerticalStressesSoil(INode) = SigY
            NodalPressuresWater(INode) = SigWP
          end do

        end subroutine ComputeNodalVerticalStressesForK0


      end module ModDynViscousBoudary