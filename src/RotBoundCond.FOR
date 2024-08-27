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
	  
	  
	  module ModRotBoundCond
      !**********************************************************************
      !
      !    Function:  Contains routines for rotating boundary conditions.
      !
      !
      !     $Revision: 8878 $
      !     $Date: 2020-09-17 03:59:14 -0400 (Thu, 17 Sep 2020) $
      !
      !**********************************************************************
      use ModCounters
      use ModReadCalculationData
      use ModGeometryMath
      use ModMatrixMath
      use ModMeshAdjacencies
      use ModMeshInfo

      implicit none
      
        integer(INTEGER_TYPE) :: NRotNodes
        integer(INTEGER_TYPE), dimension(:), allocatable :: IRotation ! For each node, ID of the plane with same degree of rotation
        real(REAL_TYPE), dimension(:, :, :), allocatable :: RotMat ! Rotation matrix for each plane
        real(REAL_TYPE), dimension(:, :, :), allocatable :: IRotMat ! Inverse of the rotation matrix for each plane

        integer(INTEGER_TYPE), private :: NBoundaryPlanes ! Number of planes whose boundary conditions are rotated
        integer(INTEGER_TYPE), private :: NBoundaryNodes ! Number of nodes whose boundary conditions should be rotated around y-axis
        integer(INTEGER_TYPE), private :: NBoundaryRadiusNodes
        integer(INTEGER_TYPE), private :: NRotationMatrices ! Total number of nodes whose boundary conditions should be rotated
        integer(INTEGER_TYPE), private :: NfixY ! Number of nodes whose displacements are fixed in the x-z-plane
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: BoundaryElements ! Element sides whose nodal degrees of freedom will be rotated
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: BoundaryRadiusNodes
        integer(INTEGER_TYPE), dimension(:), allocatable,  private :: BoundaryNodes ! Stores the ID's of nodes whose boundary conditions should be rotated
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: XZFixedNodes ! Stores the ID's of nodes whose displacements should be fixed in the x-z-plane
        integer(INTEGER_TYPE), dimension(:), allocatable :: AxisNodes ! Nodes shared by planes
        integer(INTEGER_TYPE), dimension(:, :), allocatable, private :: BoundaryPlanes ! Stores for each of the planes its 4 edge nodes and the prescribed displacement settings
        real(REAL_TYPE) :: BoundaryNodesRadius
        real(REAL_TYPE), dimension(:, :), allocatable, private :: BoundarySurfacePrescrDisp ! Prescribed displacements applied to boundary surfaces

      contains ! Routines of this module

        subroutine InitialiseRotationBC(NodTot, IElTyp, NEl, ICon, CheckIElTyp) ! 3D function 
        !**********************************************************************
        !
        !    Function:  Initialise data structures for rotating boundary conditions.
        !    Note : only in 3D
        !
        !     NodTot : Total number of nodes
        !     IElTyp : Number of nodes per element (high-order tetrahedral element)
        !     NEl : Total number of elements
        !     ICon : Element connectivities for high-order element
        !     CheckIElTyp : Number of nodes per element (either low- or high-order tetrahedral element)
        !
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodTot, IElTyp, NEl, CheckIElTyp
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          ! Local variables
          integer(INTEGER_TYPE) :: IDim = 3 ! fixed dimension as only 3D functionality
          integer(INTEGER_TYPE) :: IError, NSideNodes
          integer(INTEGER_TYPE), dimension(:), allocatable :: Alpha
          integer(INTEGER_TYPE), dimension(:, :), allocatable :: BoundaryPlaneNodes ! Stores whether a node lies on one of the defined boundary planes
          real(REAL_TYPE), dimension(:, :), allocatable :: NodeNormalN ! Normal vectors for nodes with rotated boundary conditions
        
          NSideNodes = GetNSideNodes(CheckIElTyp)

          allocate(AxisNodes(NodTot), stat = IError)
          AxisNodes = 0

          if (NBoundaryPlanes > 0) then

            allocate(Alpha(NBoundaryPlanes), stat = IError)
            allocate(BoundaryPlaneNodes(NBoundaryPlanes, NodTot), stat = IError)

            call AssembleBoundaryPlaneNodes(NodTot, NodalCoordinates, BoundaryPlaneNodes)
          
            call CheckPlaneAxisNodes(NodTot, BoundaryPlaneNodes)

            ! Rotate boundary plane nodes
            Alpha = 0
            call DeterminePlaneRotationAngles(NodTot, NodalCoordinates, Alpha)
            call YRotationMatrix(Alpha)

          end if

          if (allocated(BoundaryElements).or.allocated(BoundaryRadiusNodes)) then
            allocate(NodeNormalN(3, NodTot), stat = IError)
            NodeNormalN = 0.0
          end if
          
          if (allocated(BoundaryElements)) then
            call DetermineNodeLocalCS(NodTot, IDim, IElTyp, NEl, NSideNodes, NodalCoordinates, ICon, BoundaryElements, NodeNormalN)
            call CheckStructureAxisNodes()
            call ApplyXYZRotationMatrix(NodTot, NodeNormalN)
          end if

          if (NBoundaryPlanes > 0) then
            call SetPlaneBoundaryConditions(NodTot, BoundaryPlaneNodes, Alpha)
          end if
          
          if (allocated(BoundaryRadiusNodes)) then
            call CheckRadiusAxisNodes()
            call DetermineRadiusNodeLocalCS(NodTot, NodeNormalN)
            call ApplyRadiusRotationMatrix(NodTot, NodeNormalN)
          end if
          
          if (allocated(BoundaryElements)) then
            call SetStructureBoundaryConditions()
          end if

          if (NBoundaryPlanes>0) then
            deallocate(Alpha, stat = IError)
            deallocate(BoundaryPlaneNodes, stat = IError)
          end if

          if (allocated(NodeNormalN) ) then
            deallocate(NodeNormalN, stat = IError)
          end if

        end subroutine InitialiseRotationBC

        
        subroutine DetermineRadiusNodeLocalCS(NodTot, NodeNormalN) ! 3D function
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
          integer(INTEGER_TYPE), intent(in) :: NodTot
          real(REAL_TYPE), dimension(IDim, NodTot), intent(inout) :: NodeNormalN
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE), dimension(IDim) :: NormalVector

          do I = 1, NBoundaryRadiusNodes
            NormalVector = NodalCoordinates(BoundaryRadiusNodes(I), :)
            NormalVector(2) = 0.0
            NodeNormalN(:, BoundaryRadiusNodes(I)) = VectorNorm(NormalVector, IDim)
          end do
        
        end subroutine DetermineRadiusNodeLocalCS

        
        subroutine ApplyRadiusRotationMatrix(NodTot, NodeNormalN) ! 3D function
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
          integer(INTEGER_TYPE), intent(in) :: NodTot
          real(REAL_TYPE), dimension(IDim, NodTot), intent(inout) :: NodeNormalN
          ! Local variables
          integer(INTEGER_TYPE) :: I, NodeID, RMatOffset
        
                    
          if (allocated(BoundaryRadiusNodes) ) then
            do I = 1, NBoundaryRadiusNodes
             
              NodeID = BoundaryRadiusNodes(I)
              RMatOffset = I + NBoundaryPlanes + NBoundaryNodes
              AxisNodes(NodeID) = AxisNodes(NodeID) + 1
              
              call DetermineXYZRotationMatrixAxiSymPenetration(.false., NodeID, NodeNormalN(:, NodeID), &
                                                               RotMat(:, :, RMatOffset), IRotMat(:, :, RMatOffset), .false.)
			  IRotation(NodeID) = RMatOffSet
			  NRotNodes = NRotNodes + 1
            end do
          end if
        
        end subroutine ApplyRadiusRotationMatrix

        
        subroutine CheckRadiusAxisNodes() ! 3D function
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: I

          do I = 1, NBoundaryRadiusNodes
            AxisNodes(BoundaryRadiusNodes(I)) = AxisNodes(BoundaryRadiusNodes(I)) + 1
          end do
        
        end subroutine CheckRadiusAxisNodes

        
        subroutine DeterminePlaneRotationAngles(NodTot, NodeCoord, Alpha) ! 3D function
        !**********************************************************************
        !
        !    Function:  Determine rotation angles for the planes whose boundary conditions are rotated
        !               (rotation around y-axis, all planes share common rotation axis).
        !    Note : 3D function
        !
        !     NodTot : Total number of nodes
        !     NodeCoord : Nodal coordinates
        !
        !    Alpha : Rotation angles for the planes
        !
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality 
          integer(INTEGER_TYPE), intent(in) :: NodTot
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(NBoundaryPlanes), intent(inout) :: Alpha
          ! Local variables
          integer(INTEGER_TYPE) :: I, J
          real(REAL_TYPE), dimension(IDim) :: VectorX, AxisPoint, PlanePoint
        
          VectorX = 0.0
          VectorX(1) = 1.0
          
          AxisPoint = 0.0
                  
          do I = 1, NBoundaryPlanes
            do J = 1, 4
              
              PlanePoint = NodeCoord(BoundaryPlanes(I, J), :)
              PlanePoint(2) = 0.0
              
              if (abs(Distance(PlanePoint, AxisPoint, IDim) ) > 1.0E-5) then
                Alpha(I) = nint(VectorAngle(AxisPoint, PlanePoint, VectorX, IDim) * 180.0 / PI)
              end if
              
            end do
          end do
              
          if (NBoundaryPlanes == 3) then
            Alpha(3) = 0
          end if 
                
        end subroutine DeterminePlaneRotationAngles

        
        subroutine YRotationMatrix(Alpha) ! 3D function
        !**********************************************************************
        !
        !    Function:  Sets the rotation matrix and its inverse for each node to be rotated.
        !               (rotation around y-axis)
        !
        !               Alpha is defined counterclockwise positive.
        !               RMat rotates a RHS coordinate system (!) counterclockwise.
        !               So, applying RMat on some vector returns the vector components
        !               in the new coordinate system rotated counterclockwise by Alpha.
        !
        !    Note : 3D function
        !
        !     Alpha : Rotation angles for each node
        !
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), dimension(NBoundaryPlanes), intent(in) :: Alpha
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: C, S
        
          do I = 1, NBoundaryPlanes
            C = dcos(dble(Alpha(I)) * PI / 180.0)
            S = dsin(dble(Alpha(I)) * PI / 180.0)
            
            RotMat(1, 1, I) = C
            RotMat(2, 2, I) = 1
            RotMat(3, 3, I) = C
            RotMat(1, 3, I) = -S
            RotMat(3, 1, I) = S
            
            IRotMat(1, 1, I) = C
            IRotMat(2, 2, I) = 1
            IRotMat(3, 3, I) = C
            IRotMat(1, 3, I) = S
            IRotMat(3, 1, I) = -S
          end do
        
        end subroutine YRotationMatrix

        
        subroutine ApplyXYZRotationMatrix(NodTot, NodeNormalN) ! 3D function
        !**********************************************************************
        !
        !    Function:  Sets the rotation matrix and its inverse for each node to be rotated.
        !               (rotation around arbitrary axes)
        !    Note : 3D function
        !
        !     NodTot : Total number of nodes
        !     NodeCoord : Nodal coordinates
        !
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
          integer(INTEGER_TYPE), intent(in) :: NodTot
          real(REAL_TYPE), dimension(IDim, NodTot), intent(inout) :: NodeNormalN
          ! Local variables
          integer(INTEGER_TYPE) :: I, NodeID, RMatOffset

          if (allocated(BoundaryNodes) ) then
            do I = 1, NBoundaryNodes
             
              NodeID = BoundaryNodes(I)
              RMatOffset = I + NBoundaryPlanes
            
              call DetermineXYZRotationMatrixAxiSymPenetration(.false., NodeID, NodeNormalN(:, NodeID), &
                                                               RotMat(:, :, RMatOffset), IRotMat(:, :, RMatOffset), .false.)
              
            end do
          end if
        
        end subroutine ApplyXYZRotationMatrix

        
        subroutine DetermineXYZRotationMatrixAxiSymPenetration(UseForContact, NodeID, NodeNormal, & ! 3D function
                                                               RotationMatrix, InverseRotationMatrix, IsInterfaceNode)
        !**********************************************************************
        !
        !    Function:  Determines the rotation matrix and its inverse for the specified node.
        !               First rotation around the y-axis, then rotation around the z-axis.
        !    Note : 3D function
        !
        !     NodeCoord : Nodal coordinates
        !     NodeID : ID of the considered node
        !     RMatOffset : Identifies which RMat, RIMat is set
        !
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality 
          integer(INTEGER_TYPE), intent(in) :: NodeID
          logical, intent(in) :: UseForContact, IsInterfaceNode
          real(REAL_TYPE), dimension(IDim), intent(inout) :: NodeNormal
          real(REAL_TYPE), dimension(IDim, IDim), intent(out) :: RotationMatrix, InverseRotationMatrix
          ! Local variables
          real(REAL_TYPE), dimension(IDim) :: Projection, AxisPoint, PlanePoint, VectorX, &
                                              VectorY, VectorZ, NodeTangentS, NodeTangentT
          real(REAL_TYPE) :: Alpha, Beta, Gamma, RotationAngle, C, S, CosY, CosZ, SinY, SinZ
          integer(INTEGER_TYPE) :: AxisNode
          logical :: SwitchSignX
        
          VectorX = 0.0
          VectorX(1) = 1.0

          VectorY = 0.0
          VectorY(2) = 1.0
          
          VectorZ = 0.0
          VectorZ(3) = 1.0

          AxisPoint = 0.0

          if (allocated(AxisNodes)) then
            AxisNode = AxisNodes(NodeID)
          else
            AxisNode = 0
          end if
          
          SwitchSignX = .false.
          if (CalParams%MonoPileSliceContact.and.(NodeNormal(1)<0.0).and.IsInterfaceNode) then
            SwitchSignX = .true.  
            NodeNormal = -1.0 * NodeNormal
          end if         
          
          ! Node lies on the defined boundary planes, improve nodal local coordinate system
          if (CalParams%CorrectNormalsForSlice.and.((AxisNode==2).or.((AxisNode==1).and.UseForContact))) then
          
            Projection = NodeNormal
            Projection(2) = 0.0

            PlanePoint = NodalCoordinates(NodeID, :)
            PlanePoint(2) = 0.0
            
            if (abs(Distance(PlanePoint, AxisPoint, IDim) )>1.0E-5) then
              Alpha = VectorAngle(AxisPoint, Projection, VectorX, IDim) ! Angle between normal vector and x-axis
              Beta = VectorAngle(AxisPoint, PlanePoint, VectorX, IDim) ! Angle between plane point and x-axis
            else
              Alpha = 0.0
              Beta = 0.0
            end if
            
            RotationAngle = Beta - Alpha

            C = dcos(RotationAngle)
            S = dsin(RotationAngle)
            NodeNormal(1) = C * NodeNormal(1) +  S * NodeNormal(3)
            NodeNormal(3) = -S * NodeNormal(1) + C * NodeNormal(3)
          end if

          NodeNormal = VectorNorm(NodeNormal, IDim)
          
          call DetermineNodeTangentAxiSymPenetration(NodeNormal, NodeTangentS, NodeTangentT)

          ! Calculate rotation angles Beta, Gamma, their cosine and sine values
          Projection = NodeNormal
          Projection(2) = 0.0

          Beta = VectorAngle(AxisPoint, Projection, VectorX, IDim) ! Angle between projection of normal vector and x-axis
          Gamma = -1 * VectorAngle(AxisPoint, NodeTangentT, VectorY, IDim) ! Angle between up-vector and y-axis
          
          CosY = dcos(Beta)
          CosZ = dcos(Gamma)
          SinY = dsin(Beta)
          SinZ = dsin(Gamma)
                
          ! Column 1
          RotationMatrix(1, 1) = CosY * CosZ
          RotationMatrix(2, 1) = -1 * CosY * SinZ
          RotationMatrix(3, 1) = SinY
          ! Column 2
          RotationMatrix(1, 2) = SinZ
          RotationMatrix(2, 2) = CosZ
          RotationMatrix(3, 2) = 0
          ! Column 3
          RotationMatrix(1, 3) = -1 * SinY * CosZ
          RotationMatrix(2, 3) = SinY * SinZ
          RotationMatrix(3, 3) = CosY
            
          ! Row 1
          InverseRotationMatrix(1, 1) = RotationMatrix(1, 1)
          InverseRotationMatrix(1, 2) = RotationMatrix(2, 1)
          InverseRotationMatrix(1, 3) = RotationMatrix(3, 1)
          ! Row 2
          InverseRotationMatrix(2, 1) = RotationMatrix(1, 2)
          InverseRotationMatrix(2, 2) = RotationMatrix(2, 2)
          InverseRotationMatrix(2, 3) = RotationMatrix(3, 2)
          ! Row 3
          InverseRotationMatrix(3, 1) = RotationMatrix(1, 3)
          InverseRotationMatrix(3, 2) = RotationMatrix(2, 3)
          InverseRotationMatrix(3, 3) = RotationMatrix(3, 3)

          if (SwitchSignX) then
            NodeNormal = -1.0 * NodeNormal
          end if

        end subroutine DetermineXYZRotationMatrixAxiSymPenetration
 
        
        subroutine DetermineXYZRotationMatrix(NodeNormal, RotationMatrix, InverseRotationMatrix) ! 3D function

        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3  ! fixed dimension as only 3D functionality 
          real(REAL_TYPE), dimension(IDim), intent(inout) :: NodeNormal
          real(REAL_TYPE), dimension(IDim, IDim), intent(out) :: RotationMatrix, InverseRotationMatrix
          ! Local variables
          real(REAL_TYPE), dimension(IDim) :: NodeTangentS, NodeTangentT
        
          call DetermineNodeTangent(NodeNormal, NodeTangentS, NodeTangentT)

          ! Column 1
          RotationMatrix(1, 1) = NodeNormal(1)
          RotationMatrix(2, 1) = NodeNormal(2)
          RotationMatrix(3, 1) = NodeNormal(3)
          ! Column 2
          RotationMatrix(1, 2) = NodeTangentS(1)
          RotationMatrix(2, 2) = NodeTangentS(2)
          RotationMatrix(3, 2) = NodeTangentS(3)
          ! Column 3
          RotationMatrix(1, 3) = NodeTangentT(1)
          RotationMatrix(2, 3) = NodeTangentT(2)
          RotationMatrix(3, 3) = NodeTangentT(3)
            
          ! Row 1
          InverseRotationMatrix(1, 1) = RotationMatrix(1, 1)
          InverseRotationMatrix(1, 2) = RotationMatrix(2, 1)
          InverseRotationMatrix(1, 3) = RotationMatrix(3, 1)
          ! Row 2
          InverseRotationMatrix(2, 1) = RotationMatrix(1, 2)
          InverseRotationMatrix(2, 2) = RotationMatrix(2, 2)
          InverseRotationMatrix(2, 3) = RotationMatrix(3, 2)
          ! Row 3
          InverseRotationMatrix(3, 1) = RotationMatrix(1, 3)
          InverseRotationMatrix(3, 2) = RotationMatrix(2, 3)
          InverseRotationMatrix(3, 3) = RotationMatrix(3, 3)
          
        end subroutine DetermineXYZRotationMatrix
 
        
        logical function IsXZFixedNodes(NodeID) ! 3D function
        !**********************************************************************
        !
        !    Function:  Returns .true. if NodeID is specified as a node fixed in x- and z-direction.
        !    Note : 3Dfunction
        !
        !     NodeID : Checked node
        !
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodeID
          ! Local variables
          integer(INTEGER_TYPE) :: I
          
          IsXZFixedNodes = .false.
          
          do I = 1, NfixY
            if (XZFixedNodes(I) == NodeID) then
              IsXZFixedNodes = .true.
              EXIT
            end if
          end do
        
        end function IsXZFixedNodes

        
        subroutine SetPlaneBoundaryConditions(NodTot, BoundaryPlaneNodes, Alpha) ! 3D function
        !**********************************************************************
        !
        !    Function:  Initialise data structures for rotating boundary conditions of planes.
        !    Note : 3D function
        !
        !     NodTot : Total number of nodes
        !     BoundaryPlaneNodes : Stores for each plane for each node
        !                          that lies on the plane a non-zero value
        !     Alpha : Angle by which the planes should be rotated
        !
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimensiona as only 3D functionality
          integer(INTEGER_TYPE), intent(in) :: NodTot
          integer(INTEGER_TYPE), dimension(NBoundaryPlanes, NodTot), intent(in) :: BoundaryPlaneNodes
          integer(INTEGER_TYPE), dimension(NBoundaryPlanes), intent(in) :: Alpha
          ! Local variables
          integer(INTEGER_TYPE) :: I, NodeID

          do I = 1, NBoundaryPlanes ! Loop over defined planes
            if (Alpha(I)>0) then
              do NodeID = 1, NodTot ! Loop over nodes
                if ( (BoundaryPlaneNodes(I, NodeID)>0).and. &
                     (.not.IsXZFixedNodes(NodeID) ) ) then ! Node NodeID lies on plane I
                
                  if (AxisNodes(NodeID)==1) then ! Node is not fixed in x- and z-direction
                    IRotation(NodeID) = I
                    NRotNodes = NRotNodes + 1           
                  end if                 
                end if
              end do
            end if
          end do
        
        end subroutine SetPlaneBoundaryConditions

        
        subroutine SetStructureBoundaryConditions() ! 3D function
        !**********************************************************************
        !
        !    Function:  Initialise data structures for rotating boundary conditions of structures.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, NodeID

          BoundarySurfacePrescrDisp = 0.0
          if (allocated(BoundaryNodes) ) then
            do I = 1, NBoundaryNodes ! Loop over defined structure nodes
              NodeID = BoundaryNodes(I)

              if (.not.IsXZFixedNodes(NodeID) ) then
                if (AxisNodes(NodeID)/=3) then
                  IRotation(NodeID) = I + NBoundaryPlanes
                  NRotNodes = NRotNodes + 1

                  if (AxisNodes(NodeID)==1) then ! Node does not lie on the mesh boundary planes
                    do J = 1, IDim  
                      NodalPrescibedDisp(NodeID, J) = BoundarySurfacePrescrDisp(1, J) ! 0 = no normal displacement
                    end do
                  else
                    do J = 1, IDim-1
                      NodalPrescibedDisp(NodeID, J) = BoundarySurfacePrescrDisp(1, J) ! 0 = no normal displacement
                    end do
                  end if
                else
                  do J = 1, IDim
                    NodalPrescibedDisp(NodeID, IDim) = 0 ! no displacement
                  end do
                end if
              else
                ! These nodes are required to remain fixed, nodes on the tip
                ! can only move inside the tip surface or they would cause
                ! particles to leave the mesh (above the pile tip surface)
              end if
            
            end do
          end if

        end subroutine SetStructureBoundaryConditions

        
        subroutine AssembleBoundaryPlaneNodes(NodTot, NodeCoord,  BoundaryPlaneNodes) ! 3D function
        !**********************************************************************
        !
        !    Function:  Checks which nodes lie on the defined boundary planes.
        !    Note : 3D function
        !
        !     NodTot : Total number of nodes
        !     NodeCoord : Nodal coordinates
        !     BoundaryPlaneNodes : Stores for each plane for each node
        !                          that lies on the plane a non-zero value 
        !
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
          integer(INTEGER_TYPE), intent(in) :: NodTot
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(NBoundaryPlanes, NodTot), intent(out) :: BoundaryPlaneNodes
          ! Local variables
          integer(INTEGER_TYPE) :: I, J
          real(REAL_TYPE) :: Distance
          real(REAL_TYPE), dimension(IDim) :: Point0, Point2, Point3, Vector1, Vector2, &
                                              VectorA, VectorN, Node, Cross
          real(REAL_TYPE) :: Tolerance = 1E-5
          
          if (BoundingBoxDiameter>5.0) then
            Tolerance = 1E-3
          end if
          
          BoundaryPlaneNodes = 0
          
          do I = 1, NBoundaryPlanes ! Loop over defined planes
            
            Point0 = NodeCoord(BoundaryPlanes(I, 1), :)
            Point2 = NodeCoord(BoundaryPlanes(I, 2), :)
            Point3 = NodeCoord(BoundaryPlanes(I, 3), :)
            
            Vector1 = Point0 - Point2
            Vector2 = Point0 - Point3
            Cross   = CrossProduct(Vector1, Vector2)
            VectorN = VectorNorm(Cross, 3)
            
            do J = 1, NodTot ! Loop over nodes
            
              Node = NodeCoord(J, :)
              VectorA = Node - Point0
              
              Distance = abs(DotProduct(VectorA, VectorN, IDim) )
              if (Distance<Tolerance) then ! Node J lies on plane I
                BoundaryPlaneNodes(I, J) = J
              end if
              
            end do
          end do
        
        end subroutine AssembleBoundaryPlaneNodes

        
        subroutine CheckPlaneAxisNodes(NodTot, BoundaryPlaneNodes) ! 3D function
        !**********************************************************************
        !
        !    Function:  Check which nodes lie on more than one plane (stored in AxisNodes).
        !
        !     NodTot : Total number of nodes
        !     BoundaryPlaneNodes : Stores for each plane for each node
        !                          that lies on the plane a non-zero value
        !
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodTot
          integer(INTEGER_TYPE), dimension(NBoundaryPlanes, NodTot), intent(in) :: BoundaryPlaneNodes
          ! Local variables
          integer(INTEGER_TYPE) :: I, NodeID

          do I = 1, NBoundaryPlanes ! Loop over defined planes
            do NodeID = 1, NodTot ! Loop over nodes
              if (BoundaryPlaneNodes(I, NodeID)>0) then ! Node NodeID lies on plane I
                AxisNodes(NodeID) = AxisNodes(NodeID) + 1
              end if
            end do
          end do
        
        end subroutine CheckPlaneAxisNodes

        
        subroutine CheckStructureAxisNodes() ! 3D function
        !**********************************************************************
        !
        !    Function:  Check which nodes belong to the defined boundary (of a structure). Stored in AxisNodes.
        !
        !
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: NodeID

          if (allocated(BoundaryNodes)) then
            do NodeID = 1, NBoundaryNodes ! Loop over structure nodes
              AxisNodes(BoundaryNodes(NodeID) ) = AxisNodes(BoundaryNodes(NodeID) ) + 1
            end do
          end if
        
        end subroutine CheckStructureAxisNodes

        
	subroutine Initialise3DCylindricalAnalysis() ! 3D function
      !*************************************************************************************
      !    
      ! 
      !    Function: The geometry of the problem is analysed                
      !
      !
      !
      !*************************************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
		  integer(INTEGER_TYPE) :: I, J, K, IError
          real(REAL_TYPE) :: Radius, Teta, Dist
          real(REAL_TYPE) :: Tolerance = 1E-5
		  real(REAL_TYPE), dimension(IDim) :: Point1, Point2
		  real(REAL_TYPE), dimension(2) :: Znodes
		  real(REAL_TYPE), dimension(2, 4, IDim) :: CornerNodes

		  if ( .not.IS3DCYLINDRIC ) RETURN
          if (BoundingBoxDiameter>5.0) then
            Tolerance = 1E-3
          end if 

            ! Assume number of boundary planes is 2
             NBoundaryPlanes=2
              allocate(BoundaryPlanes(NBoundaryPlanes, 4), stat = IError) ! 4 Node IDs defining the plane
             ! fill CornerNodes
			  BoundaryNodesRadius=BoundingBoxMax(1)
			  Znodes(1)=BoundingBoxMin(3)
			  Znodes(2)=BoundingBoxMax(3)
			  do I= 1, NBoundaryPlanes
				  Teta=asin(Znodes(I)/BoundaryNodesRadius)				  
  					  do J= 1, 2						  
                            ! Permuted X and Z coordinates						  
							CornerNodes(I,(2*J-1),1)=BoundingBoxMin(1)*cos(teta)
							CornerNodes(I,(2*J),1)=BoundingBoxMax(1)*cos(teta)
							CornerNodes(I,(2*J-1),3)=BoundingBoxMin(1)*sin(teta)
							CornerNodes(I,(2*J),3)=BoundingBoxMax(1)*sin(teta)
					  end do
					  do J= 1, 2						  
                            ! Permuted y coordinates						  
							CornerNodes(I,J,2)=BoundingBoxMin(2)
							CornerNodes(I,J+2,2)=BoundingBoxMax(2)							
					  end do				  
			  end do
	   
              ! Determine the nodes of each CornerNode
              do I = 1, NBoundaryPlanes
				do J= 1, 4
					Point1=CornerNodes(I, J, :)
					do K=1, Counters%NodTot						
						Point2=NodalCoordinates(K,:)
						Dist=Distance(Point1, Point2, IDim)
						if (abs(Dist) <= Tolerance) then
							BoundaryPlanes(I,J)=K
							Exit
						end if						
					end do
				end do
			  end do   
            
            ! Read all node IDs with a distance of BoundaryNodesRadius from the y-axis
            ! Assemble boundary elements from these nodes
           
            if (BoundaryNodesRadius>0.0) then
              
              NBoundaryRadiusNodes = 0
              do I = 1, Counters%NodTot
                Radius = dsqrt( &
                  NodalCoordinates(I, 1) * NodalCoordinates(I, 1) +  &
                  NodalCoordinates(I, 3) * NodalCoordinates(I, 3))
                if ((abs(Radius - BoundaryNodesRadius)<Tolerance) ) then
                  NBoundaryRadiusNodes = NBoundaryRadiusNodes + 1
                end if 
              end do
              
              allocate(BoundaryRadiusNodes(NBoundaryRadiusNodes), stat = IError)

              NBoundaryRadiusNodes = 0
              do I = 1, Counters%NodTot
                Radius = dsqrt( &
                  NodalCoordinates(I, 1) * NodalCoordinates(I, 1) +  &
                  NodalCoordinates(I, 3) * NodalCoordinates(I, 3))
                if ((abs(Radius - BoundaryNodesRadius)<Tolerance) ) then
                  NBoundaryRadiusNodes = NBoundaryRadiusNodes + 1
                  BoundaryRadiusNodes(NBoundaryRadiusNodes) = I
                end if 
			  end do	  
     
            end if  
            
            NRotationMatrices = NBoundaryPlanes +  &
                                NBoundaryNodes + &
                                NBoundaryRadiusNodes
            
            ! Fill the nodes that are fixed in x-z-plane
            NfixY=0 
              do I = 1, Counters%NodTot
                if (NodalPrescibedDisp(I,2) == 0d0) then
					NfixY=NfixY+1
				end if
			  end do
			  
            if (NfixY>0) then              
              allocate(XZFixedNodes(NfixY), stat = IError)              
              ! Search for the nodes IDs with fixities in the y direction
              XZFixedNodes = 0
			  NfixY=0
              do I = 1, Counters%NodTot
                if (NodalPrescibedDisp(I,2) == 0d0) then
					NfixY=NfixY+1
					XZFixedNodes(NfixY)=I
				end if
              end do              
            end if            
	end subroutine Initialise3DCylindricalAnalysis        
                
        subroutine DestroyRotBoundCond() ! 3D function
        !**********************************************************************
        !
        !    Function:  Deallocates the data structures used for rotating boundary conditions.
        !
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError
                
          if (allocated(BoundaryPlanes) ) then
            deallocate(BoundaryPlanes, stat = IError)
          end if
        
          if (allocated(BoundaryNodes) ) then
            deallocate(BoundaryNodes, stat = IError)
          end if

          if (allocated(BoundaryElements) ) then
            deallocate(BoundaryElements, stat = IError)
          end if

          if (allocated(BoundarySurfacePrescrDisp) ) then
            deallocate(BoundarySurfacePrescrDisp, stat = IError)
          end if
        
          if (allocated(XZFixedNodes) ) then
            deallocate(XZFixedNodes, stat = IError)
          end if

          if (allocated(AxisNodes) ) then
            deallocate(AxisNodes, stat = IError)
          end if

          if (allocated(IRotMat) ) then
            deallocate(IRotMat, stat = IError)
          end if
          
          if (allocated(RotMat) ) then
            deallocate(RotMat, stat = IError)
          end if

          if (allocated(IRotation) ) then
            deallocate(IRotation, stat = IError)
          end if

          if (allocated(BoundaryRadiusNodes) ) then
            deallocate(BoundaryRadiusNodes, stat = IError)
          end if

          NRotNodes = 0
          NBoundaryPlanes = 0
          NBoundaryNodes = 0
          NBoundaryRadiusNodes = 0
          NRotationMatrices = 0
          NfixY = 0
          BoundaryNodesRadius = 0.0

        end subroutine DestroyRotBoundCond

        
        subroutine InitialiseRotationMatrix() ! 3D function
        !**********************************************************************
        !
        !    Function:  To fill the arrays RotMat, IRotMat
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
          integer(INTEGER_TYPE) :: IError
          allocate(IRotation(Counters%NodTot), stat = IError)
          
          IRotation = 0
          NRotNodes = 0

          if (IS3DCYLINDRIC) then
            allocate(RotMat(IDim, IDim, NRotationMatrices), IRotMat(IDim, IDim, NRotationMatrices), stat = IError)
            RotMat = 0.0
            IRotMat = 0.0
            
            call InitialiseRotationBC(Counters%NodTot, N_NODES_HOE, Counters%NEl, &
                                      ElementConnectivities10Node, ELEMENTNODES)
          else
            allocate(RotMat(1,1,1), IRotMat(1,1,1), stat = IError)
            RotMat = 0.0
            IRotMat = 0.0
          end if
                   
        end subroutine InitialiseRotationMatrix

        
      subroutine RotateVector(RMat, VectorIn, VectorOut) ! 3D function
      !****************************************************************************
      !     Function: To rotate a vector with 3 components
      !
      !        RMat: rotation matrix
      !        VectorIn: vector to rotate
      !        VectorOut: rotated output vector
      !
      !****************************************************************************
      implicit none

        integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
        real(REAL_TYPE), dimension(IDim, IDim), intent(in) :: RMat
        real(REAL_TYPE), dimension(IDim), intent(in) :: VectorIn
        real(REAL_TYPE), dimension(IDim), intent(out) :: VectorOut
        ! Local variables
        integer :: I, J
        real(REAL_TYPE) :: Value
        
        do I = 1, IDim
          Value = 0.0
          do J = 1, IDim
            Value = Value + RMat(I, J) * VectorIn(J)
          end do
          VectorOut(I) = Value
        end do
      
      end subroutine RotateVector

      
      end module ModRotBoundCond


      subroutine RotVec(IRot, NRotNodesLoc, RMat, NDofEx, V_In, V_Out) ! 3D function
!****************************************************************************
!     Function:  To compute a (load)vector in adapted coordinates
!     Note : 3D function
!      
!     Input:  : RMat  : Array containing different rotation matrices
!               IRot  : Which rotation matrix
!               V_In  : Vector to rotate
!     Output  : V_Out : Rotated output vector
!
!     When calling it can be called using V_out = V_in !!
!
!****************************************************************************
      use ModCounters
      use ModRotBoundCond
      use ModGlobalConstants
      
      implicit none

      integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
      integer(INTEGER_TYPE) :: NRotNodesLoc
      integer(INTEGER_TYPE), dimension(Counters%NodTot) :: IRot, NDofEx
      real(REAL_TYPE) :: RMat(IDim, IDim, *) ! third dimension can be NRotationMatrices (if ApplyRotBoundCond) or =1 (all other cases)
      real(REAL_TYPE), dimension(Counters%N) :: V_In, V_Out
      ! local variables
      real(REAL_TYPE) :: Rotation(IDim, IDim), F(IDim), FR(IDim)
      integer(INTEGER_TYPE) :: I, J, IDof
      
      ! copy also non-rotated DOFs
       V_Out(1:Counters%N)=V_In(1:Counters%N)
       
      If (NRotNodesLoc == 0) RETURN
      
      Do I = 1, Counters%NodTot
        IDof = NDofEx(I)
        If (IRot(I) /= 0) Then
          Rotation = RMat(:, :, IRot(I))
          Do J = 1, IDim
            F(J) = V_In(IDof+J)
          End Do
          Call RotateVector(Rotation, F, Fr)
          Do J = 1, IDim
            V_Out(IDof+J) = FR(J)
          End Do
        End If
      End Do
 
      end subroutine RotVec