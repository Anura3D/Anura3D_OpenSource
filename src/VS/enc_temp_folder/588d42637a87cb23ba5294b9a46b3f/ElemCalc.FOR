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


      module ModElementEvaluation
      !**********************************************************************
      !
      !     ModElementEvaluation:  This module contains routines related to Finite Element
      !                            evaluation independent of the type of element
      !                            (Shape functions, strain interpolation matrix,
      !                            Jacobian matrix).
      !
      !                            Whenever, element specific data is needed it is referred
      !                            to the corresponding source files:
      !                             - ModElementEvaluationTETRA
      !                             - ModElementEvaluationTRI
      !
      !     $Revision: 8842 $
      !     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
      !
      !**********************************************************************
          
      use ModElementEvaluationTETRA
      use ModElementEvaluationTRI
      use ModString
      use ModFeedback
      use ModGlobalConstants
      use ModCounters
      use ModGeometryMath
      
      implicit none
      
        real(REAL_TYPE), dimension(:), allocatable :: GPWeight 
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative

        real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary


      contains ! Routines of this module

      
      subroutine GaussPointLocalCoordinates(IGaussPoint, WeiGP, PosGP)
      !**********************************************************************
      !
      !>   GaussPointLocalCoordinates:  Determines the local coordinates and integration weight assigned to Gauss point
      !                                 IGaussPoint which are returned through WeiGP and PosGP.
      !
      !> IN:
      !> IGaussPoint : Number of the Gauss point
      !
      !> OUT:
      !> WeiGP : Initial weight assigned to Gauss point IGaussPoint
      !> PosGP : Initial local position of Gauss point IGaussPoint
      !
      !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
        
          WeiGP = 0.0
          PosGP = 0.0
    
          call Gauss_Q1Pointer(IGaussPoint, PosGP, WeiGP)     
        
        end subroutine GaussPointLocalCoordinates
        

        integer(INTEGER_TYPE) function GetNSideNodes(IElTyp)
        !**********************************************************************
        !
        !> GetNSideNodes:  Returns the number of nodes on each element side
        !!
        !> Implemented in the frame of the MPM project.
        !!
        !> IN:
        !!
        !>    IElTyp : Number of node connectivities
        !! - = 4 : 4-noded tetrahedral element
        !! - = 10 : 10-noded tetrahedral element
        !!
        !> OUT:
        !!
        !> GetNSideNodes : Number of nodes per element side
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp !< Number of node connectivities

          GetNSideNodes = 0 
          
          select case(IElTyp)
            case(10) ! 10-noded tetrahedral element  ! for 3D only  
              GetNSideNodes = 6
            case(4) ! 4-noded tetrahedral element  ! for 3D only
              GetNSideNodes = 3
          end select
        
        end function GetNSideNodes 

        
        subroutine ShapeXiEtaT(NSideNodes, Xi, Eta, HS, DHS)
        !**********************************************************************
        !
        !    Function:  Determines the nodal shape function values HS and derivatives DHS for
        !               (Xi, Eta) for either a 3-noded or 6-noded tetrahedral element - depending
        !               on NSideNodes (either 3 or 6).
        !
        !     Xi, Eta : Local coordinates
        !
        ! O   HS : Shape function values at Xi, Eta
        ! O   DHS : Shape function derivatives at Xi, Eta
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NSideNodes
          real(REAL_TYPE), intent(in) :: Xi, Eta
          real(REAL_TYPE), dimension(NSideNodes), intent(out) :: HS
          real(REAL_TYPE), dimension(NSideNodes, 2), intent(out) :: DHS
        
          select case(NSideNodes)
            case(6) ! 6-noded triangular element
              call ShapeXiEtaTetrahedronHOE(Xi, Eta, HS, DHS)
            case(3) ! 3-noded triangular element
              call ShapeXiEtaTetrahedronLOE(Xi, Eta, HS, DHS)
          end select
        
        end subroutine ShapeXiEtaT

        
        subroutine ShapeFunctionData(LocPos, IElTyp, ShapeValues, DShapeValues)
        !**********************************************************************
        !
        !    Function:  Evaluates the shape functions and shape function derivatives
        !               at a local coordinate LocPos inside an element.
        !
        !     LocPos : Point inside an element in the local coordinate system
        !     IElTyp : Number of nodes per element
        !     IDim : Number of dimensions
        !
        ! O   ShapeValues : Shape function values at LocPos
        ! O   DShapeValues : Shape function derivatives at LocPos
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          integer(INTEGER_TYPE), intent(in) :: IElTyp
          real(REAL_TYPE), dimension(:), intent(inout) :: ShapeValues
          real(REAL_TYPE), dimension(:, :), intent(inout) :: DShapeValues

        !Temp Solution for TetraOld           
          select case(ELEMENTTYPE)
          case(TETRAOLD)
              select case(IElTyp)
              case(10) ! 10-noded tetrahedral element  ! for 3D only
                ShapeLocPosPointer => ShapeLocPosTETRA10
              case(4) ! 4-noded tetrahedral element  ! for 3D only
                ShapeLocPosPointer => ShapeLocPosTETRA4
              end select
         end select
         call ShapeLocPosPointer(LocPos, ShapeValues, DShapeValues)  

        
        end subroutine ShapeFunctionData

        
        subroutine DetJacob(LocPos, NEl, NodTot, IDim, IElement, ICon, Co, RJac, InvRJac, DetJac)
        !**********************************************************************
        !
        !  Function : Determination of the Jacobian matrix and the determinant of
        !             the Jacobian matrix for the location LocPos inside an element IElement
        !             whose nodal connectivities are defined by ICon and Co depending on
        !             the type of element.
        !
        !  I  LocPos : Local coordinates of the considered point inside an element
        !  I  NEl : Total number of elements
        !  I  NodTot : Total number of nodes
        !  I  IDim : Number of dimensions
        !  I  IElement : ID of the element
        !  I  ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !  I  Co : Global nodal coordinates Co(I, J): j-coordinate of node I
        !
        !  O  RJac : Jacobian matrix
        !  O  InvRJac : Inverse of the Jacobian matrix
        !  O  DetJac : Determinate of the Jacobian matrix
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: NEl, NodTot, IDim, IElement
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(:, :), intent(in) :: Co
          real(REAL_TYPE), dimension(:, :), intent(inout) :: RJac
          real(REAL_TYPE), dimension(:, :), intent(inout) :: InvRJac
          real(REAL_TYPE), intent(inout) :: DetJac
          
          ! local variables
          integer(INTEGER_TYPE) :: I, J, INode, NodeID, NNodes
          real(REAL_TYPE), dimension(:), allocatable :: HS
          real(REAL_TYPE), dimension(:, :), allocatable :: dHS
          real(REAL_TYPE) :: Det1
          
          NNodes = size(ICon,1)
          allocate( HS(NNodes), dHS(NNodes, IDim) )
          
          ! Determine the shape functions HS and shape function derivatives dHS for the element for LocPos.
          call ShapeFunctionData(LocPos, NNodes, HS, dHS)

          ! Determine the Jacobian matrix RJac
          RJac = 0.0
          do INode = 1, NNodes ! loop nodes of each element
            NodeID = ICon(INode, IElement)
            do I = 1, IDim
              do J = 1, IDim
                RJac(I, J) = RJac(I, J) + dHS(INode, I) * Co(NodeID, J)
              end do
            end do
          end do
          
          ! Calculate inverse and determinant of Jacobian matrix
          call RJacInv(IDim, RJac, InvRJac, DetJac, Det1)

          if ( (DetJac < 0.0) .or. (Det1 < 0.0) ) then
            call WriteInLogFile('Negative determinant '// trim(String(DetJac)) // trim(String(Det1)) // ' element ' // trim(String(IElement)))
          end if

        end subroutine DetJacob

        
        subroutine BMatrix(LocPos, IElTyp, NEl, NodTot, IDim, IElement, ICon, NodeCoord, B, DetJac)
        !**********************************************************************
        !
        !    Function:  Determination of the strain interpolation matrix (B matrix)
        !               for the location LocPos inside an element IElement whose nodal
        !               connectivities are defined by ICon and NodeCoord.
        !
        !     LocPos : Position inside an element in local coordinates
        !     IElTyp : Number of nodes per element
        !     NEl : Total number of elements
        !     NodTot : Total number of nodes
        !     IDim : Number of dimensions
        !     IElement : ID of the element
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !     NodeCoord : Global nodal coordinates Co(I, J): J-coordinate of node I
        !
        ! O   B : IDim x IElTyp matrix containing the strain interpolation terms
        ! O   DetJac : Determinante of Jacobian
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
          real(REAL_TYPE), dimension(IDim), intent(in) :: LocPos
          integer(INTEGER_TYPE), intent(in) :: IElement
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          real(REAL_TYPE), dimension(IDim, IElTyp), intent(out) :: B
          real(REAL_TYPE), intent(out) :: DetJac
          
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, K
          real(REAL_TYPE), dimension(IElTyp) :: HS ! Shape functions
          real(REAL_TYPE), dimension(IElTyp, IDim) :: dHS ! Derivatives of shape functions
          real(REAL_TYPE), dimension(IDim, IDim) :: RJac, RJacInv ! Jacobian matrix, inverse of Jacobian matrix

          ! Determine the shape functions HS and shape function derivatives dHS for LocPos.
          call ShapeFunctionData(LocPos, IElTyp, HS, dHS)

          ! Calculate Jacobian matrix RJac and the inverse of the Jacobian matrix RJacInv
          call DetJacob(LocPos, NEl, NodTot, IDim, IElement, ICon, NodeCoord, RJac, RJacInv, DetJac)
        
         ! Assemble B matrix (cartesian derivatives)
         B = 0.0
         do J = 1, IDim
           do I = 1, IElTyp
             do K = 1, IDim
               B(K, I) = B(K,I) + RJacInv(K, J) * dHS(I, J)
             end do
           end do
         end do          
          
        end subroutine BMatrix

        
        subroutine GetLocalCoordinates(GlobPos, IElement, IElTyp, NEl, NodTot, IDim, NodeCoord, ICon, LocPos, OutsideElement, Success)
        !**********************************************************************
        !
        !    Function:  Determination of local coordinates from global coordinates,
        !               assuming that the point with global coordinates lies inside 
        !               the element within IElement.
        !               Success returns .false. if the local position could not be found
        !               inside the element within 10 iterations with sufficient accuracy.
        !
        !     GlobPos : Global coordinates of a point inside IElement
        !     IElement : ID of the element
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     NodeCoord : Global nodal coordinates
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !
        ! O   LocPos : Local coordinates of the considered point inside IElement
        ! O   OutsideElement : True, if the local coordinate lie outside IElement
        ! O   Success : True, if the local position could be found
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: IElement
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(IDim), intent(out) :: LocPos
          logical, intent(out) :: OutsideElement, Success
          
          ! Local variables
          real(REAL_TYPE), dimension(IDim, IDim) :: RJac, RJacInv
          real(REAL_TYPE), dimension(IDim) :: GlobPosIteration
          real(REAL_TYPE), dimension(IDim) :: DeltaGlobPos
          integer(INTEGER_TYPE) :: I, J, Iteration
          real(REAL_TYPE) :: Tolerance, Difference, DetJac
                
          if (IDim==2) then
            Tolerance = 1d-10 * (dabs(GlobPos(1)) + dabs(GlobPos(2)) ) + 1d-10
          else
            Tolerance = 1d-10 * (dabs(GlobPos(1) ) + dabs(GlobPos(2) ) + dabs(GlobPos(3) ) ) + 1d-10  
          end if    
          
          Success = .true.
          OutsideElement = .false.
     
          ! Initial guess at local coordinates
          LocPos = (/0.2, 0.211111, 0.188888/)
          ! Global coordinates for guess
          call GetGlobalCoordinates(LocPos,  &
                                    IElTyp, NEl, NodTot, IDim, &
                                    IElement, ICon, NodeCoord, &
                                    GlobPosIteration)

          Iteration = 1
          do ! Iterate towards equality of local and global coordinates
            ! Calculate Jacobian matrix RJac and the inverse of the Jacobian matrix RJacInv
            call DetJacob(LocPos, NEl, NodTot, IDim, &
                          IElement, ICon, NodeCoord, &
                          RJac, RJacInv, DetJac)

            ! Determine new LocPos
            DeltaGlobPos = GlobPos - GlobPosIteration

            do I = 1, IDim
              do J = 1, IDim
                LocPos(I) = LocPos(I) + RJacInv(J, I) * DeltaGlobPos(J)
              end do
            end do
        
            ! Global coordinates for guess
            call GetGlobalCoordinates(LocPos,  &
                                      IElTyp, NEl, NodTot, IDim, &
                                      IElement, ICon, NodeCoord, &
                                      GlobPosIteration)
        
            ! Check whether loop can be aborted
            Difference = 0.0
            do I = 1, IDim
              Difference = Difference + dabs(GlobPos(I) - GlobPosIteration(I))
            end do
     
            if (Difference < Tolerance) then ! Found local coordinates
              Success = .true.
              EXIT
            else ! Difference greater or equal Tolerance
              if (Iteration >= 10) then ! Too many iterations needed, something went wrong
                Success = .false.
                EXIT
              else ! Continue iteration
                Iteration = Iteration + 1
              end if
            end if
          end do ! Iteration loop

          if (.not.IsInsideElementLocPos(LocPos) ) then
            OutsideElement = .true.
            Success = .true.
          end if

          if (.not.Success) then
            call GiveError('Did not find local coordinates of particle in '// &
                           trim(String(IElement)) // &
                           ' within limit number of iterations.')
          end if

        end subroutine GetLocalCoordinates


        subroutine GetLocalCoordinates3(GlobPos, IDim, IElement, IElTyp, NEl, ICon, LocPos, OutsideElement, MInv, MIX1, CrossedSide)
        !**********************************************************************
        !
        !    Function:  Determination of local coordinates from global coordinates,
        !               assuming that the point with global coordinates lies inside 
        !               the element within IElement. Only for 3D
        !               Success returns .false. if the local position could not be found
        !               inside the element within 10 iterations with sufficient accuracy.
        !
        !     GlobPos : Global coordinates of a point inside IElement
        !     IDim : dimension
        !     IElement : ID of the element
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     NodeCoord : Global nodal coordinates
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !
        ! O   LocPos : Local coordinates of the considered point inside IElement
        ! O   OutsideElement : True, if the local coordinate lie outside IElement
        ! O   Success : True, if the local position could be found
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IDim
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: IElement
          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(IDim), intent(out) :: LocPos
          logical, intent(out) :: OutsideElement
          ! Local variables
          integer(INTEGER_TYPE) , dimension(IElTyp) :: NodeID
          integer(INTEGER_TYPE) :: J, INode
          real(REAL_TYPE), dimension(IDIm, IDim) :: MInv
          real(REAL_TYPE), dimension(IDim) :: MIX1
          integer(INTEGER_TYPE) :: CrossedSide
          
          OutsideElement = .true.

          do INode = 1, IElTyp
            NodeID(INode) = ICon(INode, IElement)
          end do

          
          LocPos = 0.0
          do J = 1, 3
            LocPos(3) = LocPos(3) + MInv(1, J) *  GlobPos(J)
            LocPos(1) = LocPos(1) + MInv(2, J) *  GlobPos(J)
            LocPos(2) = LocPos(2) + MInv(3, J) *  GlobPos(J)
          end do
         
          LocPos(3) = LocPos(3) - MIX1(1)
          LocPos(1) = LocPos(1) - MIX1(2)
          LocPos(2) = LocPos(2) - MIX1(3)

         CrossedSide = -1
        if     (LocPos(1) < 0) then
         CrossedSide = 2
        elseif (LocPos(2) < 0) then
         CrossedSide = 1
        elseif (LocPos(3) < 0) then
         CrossedSide = 3
        elseif ((LocPos(1) + LocPos(2) + LocPos(3)) > 1) then
         CrossedSide = 4
        else
        OutsideElement = .false.
        end if
        
        end subroutine GetLocalCoordinates3

        
        subroutine GetGlobalCoordinates(LocPos, IElTyp, NEl, NodTot, IDim, IElement, ICon, NodeCoord, GlobPos)
        !**********************************************************************
        !
        !    Function:  Determination of global coordinates from global coordinates.
        !
        !     LocPos : Local coordinates of the considered point inside an element
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     IElement : ID of the element that the point is located in
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !     NodeCoord : Global nodal coordinates Co(I, J): j-coordinate of node I
        !
        ! O   GlobPos : Global coordinates of the point inside IElement
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
          real(REAL_TYPE), dimension(IDim), intent(in) :: LocPos
          integer(INTEGER_TYPE), intent(in) :: IElement
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          real(REAL_TYPE), dimension(IDim), intent(out) :: GlobPos
          ! Local variables
          real(REAL_TYPE), dimension(IElTyp) :: ShapeValues ! Shape functions
          real(REAL_TYPE), dimension(IElTyp, IDim) :: DShapeValues ! Derivatives of shape functions
          integer(INTEGER_TYPE) :: I, J, NodeID

          ! Determine the shape functions HS and shape function derivatives dHS for LocPos.
          call ShapeFunctionData(LocPos, IElTyp, ShapeValues, DShapeValues)

          GlobPos = 0.0
          do I = 1, IElTyp
            NodeID = ICon(I, IElement)
            do J = 1, IDim
              GlobPos(J) = GlobPos(J) + ShapeValues(I) * NodeCoord(NodeID, J)
            end do  
          end do

        end subroutine GetGlobalCoordinates

        
        logical function IsInsideElementLocPos(LocPos)
        !**********************************************************************
        !
        !    Function:  Returns .true. if LocPos (local coordinates) lies inside the 
        !               volume of the considered element.
        !
        !     IElTyp : Number of node connectivities of IElement
        !
        ! O   IsInsideElementLocPos : True, if the point lies inside the element.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          
        IsInsideElementLocPos = IsInsideElementLocPosPointer(LocPos)
     
        end function IsInsideElementLocPos

        
        logical function IsInsideElementGlobPos(GlobPos, ElementID, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon)
        !**********************************************************************
        !
        !    Function:  Returns .true. if GlobPos (global coordinates) lies inside the 
        !               volume of the considered element.
        !
        !     GlobPos : Global coordinates of the considered point inside an element
        !     ElementID : ID of the considered element
        !     NodTot : Total number of nodes
        !     IDim : Number of dimensions
        !     IElTyp : Number of nodes per element
        !     NEl : Total number of elements
        !     NodeCoord : Nodal coordinates
        !     ICon : Element connectivities
        !
        ! O   IsInsideElementGlobPos : True, if the point lies inside the element.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, IElTyp, NEl
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: ElementID
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
        
          select case(IElTyp) 
            case(10) ! 10-noded tetrahedral element
              IsInsideElementGlobPos =  IsInsideElementGlobPosTETRA(GlobPos, ElementID, NodTot, IElTyp, NEl, NodeCoord, ICon)
            case(4) ! 4-noded tetrahedral element
              IsInsideElementGlobPos =  IsInsideElementGlobPosTETRA(GlobPos, ElementID, NodTot, IElTyp, NEl, NodeCoord, ICon)
          end select

        end function IsInsideElementGlobPos

        
        logical function IsCornerNode(INode, IElTyp)
        !**********************************************************************
        !
        !    Function:  Returns .true. if INode is a corner node.
        !
        !     INode : Local number of a node (1 .. IElTyp)
        !     IElTyp : Number of node connectivities
        !
        ! O   IsCornerNode : True, if INode is a corner node
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: INode
          integer(INTEGER_TYPE), intent(in) :: IElTyp
          
          IsCornerNode = .false.
          
          select case(ELEMENTTYPE)
              
            case(TRI3) ! 'triangular_3-noded'
              IsCornerNode = .true. ! all nodes are corner nodes
              
            case(TRI6) ! 'triangular_6-noded'
              if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) ) IsCornerNode = .true. ! first three nodes are corner nodes
              
            case(QUAD4) ! 'quadrilateral_4-noded'
              IsCornerNode = .true. ! all nodes are corner nodes
              
            case(QUAD8) ! 'quadrilateral_8-noded'
              if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) .or. (INode == 4) ) IsCornerNode = .true. ! first four nodes are corner nodes 
              
            case(TETRA4) ! 'tetrahedral_4-noded'
              IsCornerNode = .true. ! all nodes are corner nodes
                         
            case(TETRA10) ! 'tetrahedral_10-noded'
              if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) .or. (INode == 4) ) IsCornerNode = .true. ! first four nodes are corner nodes
                          
            case(HEXA8) ! 'hexahedral_8-noded'
              IsCornerNode = .true. ! all nodes are corner nodes
              
            case(HEXA20) ! 'hexahedral_20-noded'
              if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) .or. (INode == 4) .or. &
                   (INode == 5) .or. (INode == 6) .or. (INode == 7) .or. (INode == 8) ) IsCornerNode = .true. ! first eight nodes are corner nodes 
         
            case(TETRAOLD)
              select case(IElTyp) 
                case(10) ! 10-noded tetrahedral element
                  if ( (INode == 1).or. (INode == 2).or. (INode == 3).or. (INode == 4) ) IsCornerNode = .true.
                case(4) ! 4-noded tetrahedral element
                  if ( (INode==1).or. (INode==2).or. (INode==3).or. (INode==4) ) IsCornerNode = .true.
              end select
              
          case default
            call GiveError('Element type not defined in function IsCornerNode().')
            
          end select  
        
        end function IsCornerNode

        
        subroutine RearrangeConnectivitiesTRI3(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 3-noded triangular element.
        !               For this element type no re-ordering is required as only corner nodes exist.
        !
        ! I   IConGlobal : Node connectivities
        ! I   ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO  IConLocal : Rearranged node connectivities
        ! IO  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(2)
          IConLocal(3) = IConGlobal(3)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(2, :)
          ValuesLocal(3, :) = ValuesGlobal(3, :)

        end subroutine RearrangeConnectivitiesTRI3
        
        
        subroutine RearrangeConnectivitiesTRI6(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 6-noded triangular element.
        !               The connectivities are rearranged such that corner nodes
        !               are located in the array first, afterwards the mid-nodes.
        !
        ! I    IConGlobal : Node connectivities
        ! I    ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO   IConLocal : Rearranged node connectivities
        ! IO   ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(3)
          IConLocal(3) = IConGlobal(5)
          IConLocal(4) = IConGlobal(2)
          IConLocal(5) = IConGlobal(4)
          IConLocal(6) = IConGlobal(6)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(3, :)
          ValuesLocal(3, :) = ValuesGlobal(5, :)
          ValuesLocal(4, :) = ValuesGlobal(2, :)
          ValuesLocal(5, :) = ValuesGlobal(4, :)
          ValuesLocal(6, :) = ValuesGlobal(6, :)

        end subroutine RearrangeConnectivitiesTRI6

        
        subroutine RearrangeConnectivitiesLINE2(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 2-noded line element.
        !               For this element type no re-ordering is required as only corner nodes exist.
        !
        ! I   IConGlobal : Node connectivities
        ! I   ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO  IConLocal : Rearranged node connectivities
        ! IO  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(2)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(2, :)

        end subroutine RearrangeConnectivitiesLINE2
        
        
        subroutine RearrangeConnectivitiesLINE3(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 3-noded line element.
        !               The connectivities are rearranged such that corner nodes
        !               are located in the array first, afterwards the mid-nodes.
        !
        ! I   IConGlobal : Node connectivities
        ! I   ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO  IConLocal : Rearranged node connectivities
        ! IO  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(3)
          IConLocal(3) = IConGlobal(2)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(3, :)
          ValuesLocal(3, :) = ValuesGlobal(2, :)

        end subroutine RearrangeConnectivitiesLINE3
        
        
        subroutine RearrangeConnectivitiesQUAD4(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 4-noded quadrilateral element.
        !               For this element type no re-ordering is required as only corner nodes exist.
        !
        ! I   IConGlobal : Node connectivities
        ! I   ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO  IConLocal : Rearranged node connectivities
        ! IO  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(2)
          IConLocal(3) = IConGlobal(3)
          IConLocal(4) = IConGlobal(4)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(2, :)
          ValuesLocal(3, :) = ValuesGlobal(3, :)
          ValuesLocal(4, :) = ValuesGlobal(4, :)

        end subroutine RearrangeConnectivitiesQUAD4 
        
        
        subroutine RearrangeConnectivitiesQUAD8(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 8-noded quadrilateral element.
        !               The connectivities are rearranged such that corner nodes
        !               are located in the array first, afterwards the mid-nodes.
        !
        ! I    IConGlobal : Node connectivities
        ! I    ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO   IConLocal : Rearranged node connectivities
        ! IO   ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(3)
          IConLocal(3) = IConGlobal(5)
          IConLocal(4) = IConGlobal(7)
          IConLocal(5) = IConGlobal(2)
          IConLocal(6) = IConGlobal(4)
          IConLocal(7) = IConGlobal(6)
          IConLocal(7) = IConGlobal(8)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(3, :)
          ValuesLocal(3, :) = ValuesGlobal(5, :)
          ValuesLocal(4, :) = ValuesGlobal(7, :)
          ValuesLocal(5, :) = ValuesGlobal(2, :)
          ValuesLocal(6, :) = ValuesGlobal(4, :)
          ValuesLocal(7, :) = ValuesGlobal(6, :)
          ValuesLocal(7, :) = ValuesGlobal(8, :)

        end subroutine RearrangeConnectivitiesQUAD8
        
        
        function IntegrateVectorSurface(NSurfaceNodes, NGP, NodTot, NodeCoord, IConSurface, NodeValues)
        !**********************************************************************
        !
        !    Function:  Returns the integral of the distributed values defined by NodeValues
        !               for a 6-noded triangular surface (3D) and 2-noded line (2D) by Gauss point integration.
        !
        !     NSurfaceNodes : Number of nodes of the surface (=6 for 6-noded triangular element, =2 for 2-noded linear element)
        !     NGP : Number of Gauss points
        !     NodTot : Total number of nodes
        !     NodeCoord : Nodal coordinates
        !     IConSurface : Surface connectivities
        !     NodeValues : Node values integrated over the surface
        !
        ! O   IntegrateVectorSurface : Integral over surface of NodeValues
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
          implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NSurfaceNodes, NGP, NodTot
          real(REAL_TYPE), dimension(NodTot, 3), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(in) :: IConSurface
          real(REAL_TYPE), dimension(NSurfaceNodes, NDOFL), intent(in) :: NodeValues
          real(REAL_TYPE), dimension(NDOFL) :: IntegrateVectorSurface
          ! Local variables
          real(REAL_TYPE), dimension(NDOFL) :: GPValue, VectorN
          integer(INTEGER_TYPE) :: IGP, INode, IDim, nNode
          real(REAL_TYPE) :: DetJ
          
          IntegrateVectorSurface = 0.0

          do IGP = 1, NGP

            ! Determine determinante of the Jacobian
            select case(NSurfaceNodes)
              case(6) ! 6-noded triangular boundary element, only for 3D
                nNode = 3
                call Normal_T3(IGP, NodeCoord, IConSurface, nNode, GPShapeFunctionDerivativeBoundary, VectorN, DetJ)
				DetJ = DetJ * 0.5
			  case(2) ! 2-noded linear element, for 2D
                nNode = 2
                call NormalOnLine(IGP, NodeCoord, IConSurface, GPShapeFunctionDerivativeBoundary, VectorN, DetJ)
            end select
          
            ! Determine distributed value at Gauss point
            GPValue = 0.0
            do INode = 1, nNode 
              do IDim = 1, NDOFL
                GPValue(IDim) = GPValue(IDim) + GPShapeFunctionBoundary(IGP, INode) * NodeValues(INode, IDim)
              end do
            end do
          
            ! Integrate over the element surface
            do IDim = 1, NDOFL
              IntegrateVectorSurface(IDim) = IntegrateVectorSurface(IDim) + GPValue(IDim) * GPWeightBoundary(IGP) * DetJ
            end do
            
          end do
        
        end function IntegrateVectorSurface

        
        subroutine InitialiseShapeFunctions()
        !*************************************************************************************   
        !    FUNCTION:     Initialise shape functions
        ! 
        !    DESCRIPTION:        
        !>   Initialises the shape functions for all element types.
        !
        !>   @note: 
        !
        !>   @param[in] 
        !
        !>   @return
        !
        !*************************************************************************************
        implicit none
               
        ! allocate shape function variables
        allocate( GPWeightBoundary(ELEMENTBOUNDARYGAUSSPOINTS) )
        allocate( GPShapeFunctionBoundary(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES) )
        allocate( GPShapeFunctionDerivativeBoundary(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NDOFL-1) )
        allocate( GPWeight(ELEMENTGAUSSPOINTS) )
        allocate( GPShapeFunction(ELEMENTGAUSSPOINTS, ELEMENTNODES) )
        allocate( GPShapeFunctionDerivative(ELEMENTGAUSSPOINTS, ELEMENTNODES, NDOFL) )
        
        ! initialise shape function variables
        GPWeightBoundary = 0.0
        GPShapeFunctionBoundary = 0.0
        GPShapeFunctionDerivativeBoundary = 0.0
        GPWeight = 0.0
        GPShapeFunction = 0.0 
        GPShapeFunctionDerivative = 0.0
        
        call InitialiseShapeFunctionsBoundaryPointer(GPShapeFunctionBoundary, GPShapeFunctionDerivativeBoundary, GPWeightBoundary)
        call InitialiseShapeFunctionsPointer(GPShapeFunction, GPShapeFunctionDerivative, GPWeight)
       
        end subroutine InitialiseShapeFunctions
              
      
      subroutine Normal_T3(Int, Co, IConL, IelTyp3, dHS, Vn, Vl) 
!***********************************************************************
!     Determine vector V normal to a plane based on 6-noded flat element
!     Three vertices (T3)
!     first determine A = (dx/dXi , dy/dXi , dz/dXi )
!                 and B = (dx/dEta, dy/dEta, dz/dEta)
!     V1 =  a2*b3 - a3*b2
!     V2 =  a3*b1 - a1*b3
!     V3 =  a1*b2 - a2*b1
!     VL is length of (V1,V2,V3)
!     Finally normalize V
!***********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: int, ieltyp3  
      real(REAL_TYPE), dimension(:, :, :), intent(in) :: dHS 
      real(REAL_TYPE), dimension(:, :), intent(in) :: Co 
      integer(INTEGER_TYPE), dimension(:), intent(in) :: IConL
      real(REAL_TYPE), dimension(:), intent(inout) :: Vn
      real(REAL_TYPE), intent(inout) :: Vl
      
      real(REAL_TYPE), dimension(3) :: A, B
      integer(INTEGER_TYPE) :: k, nn, j
      
      Vn = 0.0
      A = 0.0
      B = 0.0

      Do K=1,IelTyp3
        NN= IConL(K)
        Do J=1,3 
          A(J) = A(J) + dHS(Int,K,1)* Co(NN,J)  ! Sum d()/dXi
          B(J) = B(J) + dHS(Int,K,2)* Co(NN,J)  ! Sum d()/dEta
        End Do
      End Do

      ! Vn is cross product of A and B
      Vn = CrossProduct(A, B)

      ! Vl is vector-length of cross product Vn
      Vl = Length(Vn, 3)

      ! normalise vector Vn
      Vn = VectorNorm(Vn, 3)

      end subroutine Normal_T3 
      
      
      subroutine NormalOnLine(IntegPoint, Co, IConL, dHS, NormalVector, VectorLength)
      !***********************************************************************
      !     Determine size of the vector formed by a line based on 2-noded line element
      !     first determine A = (dx/dXi , dy/dXi , dz/dXi )
      !     VectorLength is length of A
      !***********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IntegPoint
      real(REAL_TYPE), dimension(:, :, :), intent(in) :: dHS
      real(REAL_TYPE), dimension(:, :), intent(in) :: Co
      integer(INTEGER_TYPE), dimension(:), intent(in) :: IConL
      real(REAL_TYPE), dimension(:), intent(inout) :: NormalVector
      real(REAL_TYPE), intent(inout) :: VectorLength
      
      ! local variables
      integer(INTEGER_TYPE) :: NodeNumber, DoF, K
      real(REAL_TYPE), dimension(NVECTOR) :: A

      A = 0.0

      do K = 1, 2 ! 2D shape functions
        NodeNumber = IConL(K)
        do DoF = 1, NVECTOR
          A(DoF) = A(DoF) + dHS(IntegPoint, K, 1) * Co(NodeNumber, DoF)  ! Sum d()/dXi
        end do
      end do

      NormalVector(1) = -A(2)
      NormalVector(2) = A(1)
      
      ! length of normal vector
      VectorLength = sqrt( NormalVector(1)*NormalVector(1) + NormalVector(2)*NormalVector(2) )

      ! normalising the normal vector
      NormalVector = NormalVector / VectorLength

      end subroutine NormalOnLine
      
     subroutine DetermineSideNodes(ISide,SideNodes)
        !**********************************************************************
        !
        !    Function:  Returns the local nodeID of  ISide
        !               works for Tetrahedra and triangular element 
        !    ISide :      ID of the considered side 
        !    SideNodes : ID of the boundary node
        !
        !**********************************************************************

        implicit none

        integer(INTEGER_TYPE), intent(in) :: ISide
        integer(INTEGER_TYPE), dimension(:), intent(out):: SideNodes
        ! Local variables
        integer(INTEGER_TYPE):: J
      
        if ( NDIM == 3 ) then
          if ((ELEMENTTYPE == TETRAOLD).OR.(ELEMENTTYPE == TETRA4).OR.(ELEMENTTYPE == TETRA10)) then
            do J = 1,  ELEMENTBOUNDARYNODES
              SideNodes(J) = DetermineSideNodesTetrahedronHOE(ISide, J)
            end do
            else
               call GiveError('DetermineSideNodes not implemented for element type '//trim(String(ELEMENTTYPE))) 
            end if
        else if ( NDIM == 2 ) then
          if ((ELEMENTTYPE == TRI3).OR.(ELEMENTTYPE == TRI6)) then
            do J = 1,  ELEMENTBOUNDARYNODES  
              SideNodes(J) = DetermineSideNodesTRI6(ISide, J)
            end do            
          else
           call GiveError('determineSideNodes not implemented for element type '//trim(String(ELEMENTTYPE))) 
           end if
        end if
        
        end subroutine DetermineSideNodes
      
      end module ModElementEvaluation
      