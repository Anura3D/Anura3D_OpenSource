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
	  
	  
	  module ModMPMMeshAdjustment
      !**********************************************************************
      !
      !    Function:  Contains routines for adjusting the finite element discretisation
      !               at the end of each load step in order to implement interface and
      !               structural elements with the quasi-static MPM.
      !
      !               The discretised space is divided into three areas, two 'particle
      !               storage areas' and in between the two a 'fixed mesh area'.
      !               While the particle storage areas are compressed, extended respectively,
      !               the fixed mesh area keeps its size. The areas can also be only
      !               a plane (containing a plane structural element in case of the fixed mesh area.)
      !               The nodal coordinates of all nodes are updated according to the 
      !               displacement of a body of structure elements defined by a set of nodes
      !               that either form a volume or a plane.
      !               The areas can either be defined as wedge- or cube-shaped areas (spanned by
      !               6 or 8 nodes).
      !
      !               In order to keep the size of this source file reasonably small,
      !               this module only contains routines that are directly related to
      !               the manipulation of the finite element grid.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModMPMData
      use ModMeshAdjacencies
      use ModMeshInfo
      use ModReadMPMData
      use ModGlobalConstants
      
      implicit none

        integer(INTEGER_TYPE), parameter :: STORAGE1 = 5 ! Node inside storage area 1
        integer(INTEGER_TYPE), parameter :: STORAGE2 = 4 ! Node inside storage area 2
        integer(INTEGER_TYPE), parameter :: FIXEDMESH = 3 ! Node inside the fixed mesh
        integer(INTEGER_TYPE), parameter :: STORAGE1FIXEDMESH = 8 ! Identifies a node on the boundary between storage area 1 and the fixed mesh
        integer(INTEGER_TYPE), parameter :: STORAGE2FIXEDMESH = 7 ! Identifies a node on the boundary between storage area 2 and the fixed mesh
        integer(INTEGER_TYPE), parameter :: PLANEFIXEDMESH = 12 ! If the fixed mesh consists of a plane, the nodes belonging to the fixed mesh are identified by a value of '12'
        integer(INTEGER_TYPE), parameter :: STORAGE1STORAGE2 = 9 ! Identifies a node on the boundary between storage areas 1 and 2

        integer(INTEGER_TYPE), dimension(:, :), allocatable :: AreaNodeAssignment ! Stores for each node to which area it belongs
        logical, dimension(2) :: IsModifiedMeshAdjustment ! True, if order of areas is: fixed, storage, storage instead of: storage, fixed, storage
        real(REAL_TYPE), dimension(:), allocatable :: MeshMovement
        integer(INTEGER_TYPE), dimension(:), allocatable :: DoFMovingMeshStructure ! Nodes of a structure that is aligned with the background mesh
        real(REAL_TYPE), dimension(:, :), allocatable :: UnrotatedNodalCoordinates

      contains


        subroutine InitialiseMeshAdjustment()
        !**********************************************************************
        !
        !    Function:  Checks whether the PSA file exists. If the file exists, the
        !               mesh adjustment data structure is initialised from the file and
        !               the flag ApplyMeshSmoothing is set to .true.
        !
        !     FileName : Name of the PSA file to open
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE) :: K, IError, I, INode
        
          allocate (MeshMovement(NVECTOR))
          
          call DestroyMeshAdjustment()
                              
          if (CalParams%ApplyMeshSmoothing) then !activate mesh adjustment
            allocate(UnrotatedNodalCoordinates(1, 1), stat = IError)
            UnrotatedNodalCoordinates = 0.0


            if (CalParams%MovingMesh%MovingMaterialID>0) then ! Search all nodes that belong to volume elements
                                               ! with material StructureMaterialID
              call DetermineStructureNodes()
            else if (CalParams%ApplyContactAlgorithm.and.CalParams%ApplyContactMeshBoundary) then
              CalParams%MovingMesh%NStructureNodes = NInterfaceNodes
              allocate(CalParams%MovingMesh%StructureNodes(CalParams%MovingMesh%NStructureNodes), stat = IError)
              I = 1
              do INode = 1, Counters%NodTot
                if (InterfaceNodes(INode)) then
                  CalParams%MovingMesh%StructureNodes(I) = INode
                  I = I + 1
                end if
              end do

            end if

            if (CalParams%MovingMesh%NStructureNodes>0) then            
              ! Determine to which area (particle storage or fixed mesh area) each node belongs
            
              call DetermineAreaNodeAssignment(NodalCoordinates)            
              do K = 1, CalParams%MovingMesh%NMovingMeshDirections
                IsModifiedMeshAdjustment(K) = CheckIsModifiedMeshAdjustment(K)             
              end do
              
            else
              CalParams%ApplyMeshSmoothing = .false.
            end if
            
          else
            allocate(UnrotatedNodalCoordinates(1, 1), stat = IError)
            UnrotatedNodalCoordinates = 0.0
          end if

        end subroutine InitialiseMeshAdjustment


        subroutine DestroyMeshAdjustment()
        !**********************************************************************
        !
        !  Function : Destroy the mesh adjustment data structure
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError
                
          if (allocated(CalParams%MovingMesh%StructureNodes) ) then			  
            deallocate(CalParams%MovingMesh%StructureNodes, stat = IError)
          end if

          if (allocated(AreaNodeAssignment) ) then
            deallocate(AreaNodeAssignment, stat = IError)
          end if

          if (allocated(DoFMovingMeshStructure)) then
            deallocate(DoFMovingMeshStructure, stat = IError)
          end if
          
          if (allocated(UnrotatedNodalCoordinates)) then
            deallocate(UnrotatedNodalCoordinates, stat = IError)
          end if

        end subroutine DestroyMeshAdjustment

        subroutine MeshAdjustment()
        !**********************************************************************
        !
        !    Function:  Adjusts the nodal coordinates from the incremental displacement
        !               of the defined structural element nodes. One particle storage area
        !               is compressed, one extended and the fixed mesh is moved rigidly.
        !               Nodal coordinates of the structure inside the fixed mesh are
        !               updated from nodal displacement increments as with the Updated
        !               Lagrangian FEM.
        !
        !     DUTot : Incremental nodal displacements
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Local variables
          real(REAL_TYPE), dimension(:), allocatable :: InitialStorageAreaLength
          real(REAL_TYPE), dimension(NVECTOR) :: DUVector
          integer(INTEGER_TYPE), dimension(NVECTOR-1) :: ConsideredDirections
          integer(INTEGER_TYPE) :: I, IError
        

            call DetermineDisplacementVector(DUVector, ConsideredDirections)


          do I = 1, CalParams%MovingMesh%NMovingMeshDirections
          
            allocate(InitialStorageAreaLength(CalParams%MovingMesh%NStorageAreas(I)), stat = IError)
          
              call DetermineStorageAreaLength(I, ConsideredDirections(I), NodalCoordinates, InitialStorageAreaLength)
         
			
			! Move the moving mesh with the ground motion. Current implementation is only for 1D movement 
            if (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_FILE) then 
              DUVector(I) = CalParams%Multipliers%VelocitySolidCurrent*CalParams%TimeIncrement 
            end if 

            if (IsModifiedMeshAdjustment(I)) then
              call MoveModifiedMeshAdjustment(DUVector, I, ConsideredDirections(I))
            else
              call MoveFixedMesh(DUVector, I, ConsideredDirections(I))
            end if

            call AdjustParticleStorageAreas(DUVector, I, ConsideredDirections(I), InitialStorageAreaLength)
     
            deallocate(InitialStorageAreaLength, stat = IError)
     
          end do

        end subroutine MeshAdjustment


        subroutine DetermineDoFMovingMeshStructure() 
        !**********************************************************************
        !
        !  Function : Determines the vertical degrees of freedom that belong to a
        !             structure which is aligned to the background mesh -
        !             in case of moving mesh approach.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IEl, IPart, INode, ParticleIndex, NodeID, IError

          if (allocated(DoFMovingMeshStructure)) then
            deallocate(DoFMovingMeshStructure, stat = IError)
          end if

          allocate(DoFMovingMeshStructure(Counters%N), stat = IError)
          DoFMovingMeshStructure = 0
  
          if (CalParams%ApplyContactAlgorithm.and.CalParams%ApplyContactMeshBoundary) then
            do NodeID = 1, Counters%NodTot
              if (InterfaceNodes(NodeID)) then
                DoFMovingMeshStructure(ReducedDof(NodeID) + 2) = 1
              end if
            end do
          else
            do IEl = 1, Counters%NEl
              do IPart = 1, NPartEle(IEl)
                ParticleIndex = GetParticleIndex(IPart, IEl)
                if (EntityIDArray(ParticleIndex)==HARD_ENTITY) then
                  do INode = 1, ELEMENTNODES
                    NodeID = ElementConnectivities(INode, IEl)
                    DoFMovingMeshStructure(ReducedDof(NodeID) + 2) = 1
                  end do
                end if
              end do
            end do
          end if
          
        end subroutine DetermineDoFMovingMeshStructure
                  
        subroutine DetermineDisplacementVector(DUVector, ConsideredDirections)
        !**********************************************************************
        !
        !    Function:  Determines the average incremental displacement of the
        !               nodes belonging to the structure.
        !               Two directions with the maximum incremental displacement are considered.
        !
        ! O   DUVector : Average incremental displacement of the structure nodes in x-, y- or z-direction
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(NVECTOR), intent(out) :: DUVector
          integer(INTEGER_TYPE), dimension(2), intent(out) :: ConsideredDirections
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR) :: SummedStructureDisplacement
          integer(INTEGER_TYPE) :: INode, I, ShiftDirection, SkipDirection
          real(REAL_TYPE) :: DisplacementValue
          
          ! Determine sum of incremental displacements
          SummedStructureDisplacement = 0.0
          do INode = 1, CalParams%MovingMesh%NStructureNodes
            do I = 1, NVECTOR
              SummedStructureDisplacement(I) = SummedStructureDisplacement(I) +  &
                (NodalCoordinatesUpd(CalParams%MovingMesh%StructureNodes(INode), I) -  &
                 NodalCoordinates(CalParams%MovingMesh%StructureNodes(INode), I) )
            end do
          end do
          DUVector = SummedStructureDisplacement / CalParams%MovingMesh%NStructureNodes
          
          if(CalParams%MovingMesh%MovingMeshDirection.gt.0)then !moving mesh direction has been defined
            ShiftDirection = CalParams%MovingMesh%MovingMeshDirection
            DisplacementValue = DUVector(ShiftDirection)
            DUVector = 0.0
            DUVector(ShiftDirection) = DisplacementValue
            ConsideredDirections(1) = ShiftDirection
            ConsideredDirections(2) = 0
              
          else if (CalParams%MovingMesh%NMovingMeshDirections==1) then
            ! Determine maximum component of AverageStructureDisplacement

            if (abs(SummedStructureDisplacement(1) )>abs(SummedStructureDisplacement(2) ) ) then
              ShiftDirection = 1
            else
              ShiftDirection = 2
            end if
            if ((NDIM == 3).and.&
               (abs(SummedStructureDisplacement(3) )>abs(SummedStructureDisplacement(ShiftDirection) ) )) then
              ShiftDirection = 3
            end if
            DisplacementValue = DUVector(ShiftDirection)
            DUVector = 0.0
            DUVector(ShiftDirection) = DisplacementValue
            ConsideredDirections(1) = ShiftDirection
            ConsideredDirections(2) = 0
          else if (NDIM == 3) then! Movement in two directions
            ! Determine minimum component of AverageStructureDisplacement
            if (abs(SummedStructureDisplacement(1) )<abs(SummedStructureDisplacement(2) ) ) then
              SkipDirection = 1
            else
              SkipDirection = 2
            end if
            if (abs(SummedStructureDisplacement(3) )<abs(SummedStructureDisplacement(SkipDirection) ) ) then
              SkipDirection = 3
            end if
            DUVector(SkipDirection) = 0.0
            if (SkipDirection==1) then
              ConsideredDirections(1) = 2
              ConsideredDirections(2) = 3
            else if (SkipDirection==2) then
              ConsideredDirections(1) = 1
              ConsideredDirections(2) = 3
            else if (SkipDirection==3) then
              ConsideredDirections(1) = 1
              ConsideredDirections(2) = 2
            end if
          else
             call GiveError('Moving mesh in more than 1 direction is not implemented for NDIM= '//trim(String(NDIM))) 
          end if
                    
          MeshMovement = DUVector
        
        end subroutine DetermineDisplacementVector
        
        subroutine MoveModifiedMeshAdjustment(DUVector, IMovingMesh, Direction)
        !**********************************************************************
        !
        !    Function:  Moves those nodes which belong to both storage area 1 and 2
        !               by DUVector. Afterwards the storage areas are adjusted.
        !
        !     DUVector : Averaged incremental displacement vector of the structure nodes in x-, y- or z-direction
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: DUVector
          integer(INTEGER_TYPE), intent(in) :: IMovingMesh, Direction
          ! Local variables
          integer(INTEGER_TYPE) :: I
        
          do I = 1, Counters%NodTot
            if (AreaNodeAssignment(IMovingMesh, I)==STORAGE1STORAGE2) then ! Move node by DUVector
              NodalCoordinates(I, Direction) = NodalCoordinates(I, Direction) + DUVector(Direction)
            end if
          end do
        
        end subroutine MoveModifiedMeshAdjustment
        
        subroutine MoveFixedMesh(DUVector, IMovingMesh, Direction)
        !**********************************************************************
        !
        !    Function:  Moves the nodes of the fixed mesh, including those on its
        !               boundary and those belonging to the structure by DUVector.
        !
        !     DUVector : Averaged incremental displacement vector of the structure nodes in x-, y- or z-direction
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: DUVector
          integer(INTEGER_TYPE), intent(in) :: IMovingMesh, Direction
          ! Local variables
          integer(INTEGER_TYPE) :: I
        
          do I = 1, Counters%NodTot
            if ( (AreaNodeAssignment(IMovingMesh, I)==FIXEDMESH).or. &
                 (AreaNodeAssignment(IMovingMesh, I)==STORAGE1FIXEDMESH).or. &
                 (AreaNodeAssignment(IMovingMesh, I)==STORAGE2FIXEDMESH).or. &
                 (AreaNodeAssignment(IMovingMesh, I)==PLANEFIXEDMESH) ) then ! Node belongs to the fixed mesh area, including its boundary
              NodalCoordinates(I, Direction) = NodalCoordinates(I, Direction) + DUVector(Direction)
            end if
          end do
        
        end subroutine MoveFixedMesh
        
        subroutine AdjustStructureNodes(DUTot, DUVector, IMovingMesh, Direction)
        !**********************************************************************
        !
        !    Function:  Nodes belonging to the defined structure are updated
        !               from nodal displacements as with the Updated Lagrangian FEM.
        !               So, these nodes maintain the location they obtained during the
        !               Lagrangian phase.
        !               In case IsModifiedMeshAdjustment is true, the structure 
        !               nodes will be updated by DUTot. If the flag
        !               is false, DUVector will be subtracted first.
        !
        !     DUTot : Incremental nodal displacements
        !     DUVector : Vector by which the fixed mesh has been shifted
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
     
        implicit none
     
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: DUTot
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: DUVector
          integer(INTEGER_TYPE), intent(in) :: IMovingMesh, Direction
     
          ! Local variables
          integer(INTEGER_TYPE) :: I, NodeID, IDoF
     
          do I = 1, CalParams%MovingMesh%NStructureNodes ! Loop over structure nodes
            NodeID = CalParams%MovingMesh%StructureNodes(I)
            
            if (.not.IsModifiedMeshAdjustment(IMovingMesh)) then ! NodeCoord(I) - DUVector
              NodalCoordinates(NodeID, Direction) = NodalCoordinates(NodeID, Direction) - DUVector(Direction)
            end if
            
            ! NodeCoord(I) + DUTot(I)
            IDoF = ReducedDof(NodeID) + Direction
            NodalCoordinates(NodeID, Direction) = NodalCoordinates(NodeID, Direction) + DUTot(IDoF)
          end do
     
        end subroutine AdjustStructureNodes
        
        subroutine AdjustParticleStorageAreas(DUVector, IMovingMesh, ShiftDirection, InitialStorageAreaLength)
        !**********************************************************************
        !
        !    Function:  After moving the fixed mesh, the nodes of the particle storage areas are
        !               adjusted - one area is compressed, one area is elongated in the direction
        !               ShiftDirection.
        !
        !     DUVector : Vector by which the fixed mesh has been shifted
        !     ShiftDirection : Direction into which the mesh is shifted
        !     InitialStorageAreaLength : Initial distances between the two sides of the storage areas
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: DUVector
          integer(INTEGER_TYPE), intent(in) :: IMovingMesh, ShiftDirection
          real(REAL_TYPE), dimension(CalParams%MovingMesh%NStorageAreas(IMovingMesh)), intent(in) :: InitialStorageAreaLength
          ! Local variables
          integer(INTEGER_TYPE) :: AreaID
          
          do AreaID = 1, CalParams%MovingMesh%NStorageAreas(IMovingMesh) ! Loop over storage areas
            if (InitialStorageAreaLength(AreaID)>1E-5) then ! Storage area is not defined as plane
              call ShiftStorageAreaNodes(DUVector, AreaID, IMovingMesh, ShiftDirection, InitialStorageAreaLength(AreaID))
            end if
          end do
        
        end subroutine AdjustParticleStorageAreas

        subroutine DetermineStorageAreaLength(IMovingMesh, ShiftDirection, Coordinates, InitialStorageAreaLength)
        !**********************************************************************
        !
        !    Function:  Determines the distance between the two planes of the rectangular
        !               particle storage areas that lie perpendicular to ShiftDirection.
        !
        !     ShiftDirection : Direction into which the mesh is shifted
        !
        ! O   InitialStorageAreaLength : Distance between the two planes of the considered storage area
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IMovingMesh, ShiftDirection
          real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Coordinates
          real(REAL_TYPE), dimension(CalParams%MovingMesh%NStorageAreas(IMovingMesh)), intent(out) :: InitialStorageAreaLength
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, AreaID
          real(REAL_TYPE) :: Length
          
          InitialStorageAreaLength = -1.0
          
          do AreaID = 1, CalParams%MovingMesh%NStorageAreas(IMovingMesh)
          
            do I = 1, CalParams%MovingMesh%NAreaNodes
              do J = 1, CalParams%MovingMesh%NAreaNodes
                if (I/=J) then
                  Length =  &
                   abs(Coordinates(CalParams%MovingMesh%MeshAreas(IMovingMesh, AreaID, I), ShiftDirection) - &
                       Coordinates(CalParams%MovingMesh%MeshAreas(IMovingMesh, AreaID, J), ShiftDirection) )
                  if (Length>InitialStorageAreaLength(AreaID) ) then
                    InitialStorageAreaLength(AreaID) = Length
                  end if
                end if
              end do
            end do
            
          end do
          
        end subroutine DetermineStorageAreaLength

        subroutine ShiftStorageAreaNodes(DUVector, IArea, IMovingMesh, ShiftDirection, InitialStorageAreaLength)
        !**********************************************************************
        !
        !    Function:  Shifts the nodes of the considered storage area AreaID.
        !
        !     DUVector : Incremental displacement vector by which the fixed mesh has been shifted
        !     IArea : Number of the particle storage area
        !     ShiftDirection : Direction into which the mesh is shifted
        !     InitialStorageAreaLength : Distance between the two planes of the considered storage area
        !                                perpendicular to ShiftDirection in the initial mesh configuration
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          real(REAL_TYPE), dimension(NVECTOR), intent(in) :: DUVector
          integer(INTEGER_TYPE), intent(in) :: IArea
          integer(INTEGER_TYPE), intent(in) :: IMovingMesh, ShiftDirection
          real(REAL_TYPE), intent(in) :: InitialStorageAreaLength
          ! Local variables
          integer(INTEGER_TYPE) :: INode, AreaID 
          real(REAL_TYPE) :: BoundaryNodeDistance
        
          AreaID = GetAreaID(IArea)
        
          do INode = 1, Counters%NodTot

            if (AreaNodeAssignment(IMovingMesh, INode)==AreaID) then

            ! Determine distance between the fixed boundary and the considered node
            BoundaryNodeDistance =  &
               abs(NodalCoordinates(INode, ShiftDirection) -  &
                      NodalCoordinates(CalParams%MovingMesh%MeshAreas(IMovingMesh, IArea, 1), &
                      ShiftDirection) )


              ! Loop over nodes of considered storage area that are not on the fixed boundary
              if (BoundaryNodeDistance>1E-5) then
                
                ! Shift the considered node
                NodalCoordinates(INode, ShiftDirection) = &
                  NodalCoordinates(INode, ShiftDirection) + &
                  DUVector(ShiftDirection) / InitialStorageAreaLength * &
                  BoundaryNodeDistance

              end if
            end if
          end do

        end subroutine ShiftStorageAreaNodes

        subroutine DetermineStructureNodes()
        !**********************************************************************
        !
        !    Function:  Determines the nodes that belong to the volume structure
        !               elements which form the structure that moves through
        !               the particle discretisation. The element material ID that
        !               is compared to StructureMaterialID is taken from the particles
        !               inside the activated elements.
        !
        !               Note: Particles inside structural volume elements should never
        !                     leave these elements.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
      
          ! Local variables
          integer(INTEGER_TYPE) :: IEl, ParticleIndex, INode, NodeID, IError
          integer(INTEGER_TYPE), dimension(:), allocatable :: NodeCheck
      
          ! Determine nodes that belong to structure volume elements
          allocate(NodeCheck(Counters%NodTot), stat = IError)
          NodeCheck = 0
          
          do IEl = 1, Counters%NEl
            if (IsActiveElement(IEl)) then 
              if (NPartEle(IEl)>0) then ! Loop over all activated elements that contain particles
                
                ParticleIndex = GetParticleIndex(1, IEl)
                if (MaterialIDArray(ParticleIndex)== &
                    CalParams%MovingMesh%MovingMaterialID) then ! Found an element that is part of the structure
                  
                  do INode = 1, ELEMENTNODES
                    NodeID = ElementConnectivities(INode, IEl)
                    NodeCheck(NodeID) = 1
                  end do
                
                end if
                
              end if
            end if
          end do
          
          call EvaluateNodeCheck(NodeCheck)
        
          deallocate(NodeCheck, stat = IError)
          
        end subroutine DetermineStructureNodes

        subroutine EvaluateNodeCheck(NodeCheck)
        !**********************************************************************
        !
        !    Function:  Determines the number of structure nodes and the ID's of
        !               nodes belonging to the structure volume elements storing
        !               this data in the globally defined variables NStructureNode and
        !               StructureNodes.
        !
        !     NodeCheck : Array containing a flag for each node belonging to a structure volume 
        !                 element
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), dimension(Counters%NodTot), intent(in) :: NodeCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, INode, IError, count
        
          ! Determine number of structure nodes
          CalParams%MovingMesh%NStructureNodes = 0
          do INode = 1, Counters%NodTot
            
            if (NodeCheck(INode)==1) then
              CalParams%MovingMesh%NStructureNodes = CalParams%MovingMesh%NStructureNodes + 1
            end if
            
          end do
          
          ! Determine array of structure nodes
          count = CalParams%MovingMesh%NStructureNodes
          allocate(CalParams%MovingMesh%StructureNodes(count), stat = IError)
          I = 1
          do INode = 1, Counters%NodTot
          
            if (NodeCheck(INode)==1) then
              CalParams%MovingMesh%StructureNodes(I) = INode
              I = I + 1
            end if
          
          end do
        
        end subroutine EvaluateNodeCheck

        

        subroutine DetermineAreaNodeAssignment(InitialNodeCoordinates)
        !*************************************************************************************
        !    SUBROUTINE: DetermineAreaNodeAssignment
        ! 
        !    DESCRIPTION:
        !>   Determines which node belongs to which area - particle storage or fixed mesh area.
        !
        !>   @param[in] InitialNodeCoordinates : Initial nodal coordinates
        !
        !*************************************************************************************
        implicit none
        
          real(REAL_TYPE), dimension(:, :), intent(in) :: InitialNodeCoordinates
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, K, IDArea, IError
          real(REAL_TYPE), dimension(CalParams%MovingMesh%NAreaNodes, NVECTOR) :: BoxAreaCoordinates

          allocate(AreaNodeAssignment(CalParams%MovingMesh%NMovingMeshDirections, Counters%NodTot), stat = IError)
          
          AreaNodeAssignment = 0
          BoxAreaCoordinates = 0.0
        
          do K = 1, CalParams%MovingMesh%NMovingMeshDirections
          
            ! IDArea = 1: Check particle storage area 1
            ! IDArea = 2: Check particle storage area 2
            do IDArea = 1, CalParams%MovingMesh%NStorageAreas(K)

              do I = 1, size(BoxAreaCoordinates, 1) ! loop corner nodes of the areas
                do J = 1, size(BoxAreaCoordinates, 2) ! loop dimensions
                  BoxAreaCoordinates(I, J) = InitialNodeCoordinates(CalParams%MovingMesh%MeshAreas(K, IDArea, I), J)
                end do  
              end do
            
              call CheckNodeInArea(K, InitialNodeCoordinates, IDArea, BoxAreaCoordinates)
            end do

            ! IDArea = 3: Check fixed mesh area
            IDArea = 3
            do I = 1, size(BoxAreaCoordinates, 1) ! loop corner nodes of the areas
              do J = 1, size(BoxAreaCoordinates, 2) ! loop dimensions
                BoxAreaCoordinates(I, J) = InitialNodeCoordinates(CalParams%MovingMesh%MeshAreas(K, IDArea, I), J)
              end do   
            end do
            
            call CheckNodeInArea(K, InitialNodeCoordinates, IDArea, BoxAreaCoordinates)
        
          end do
          
        end subroutine DetermineAreaNodeAssignment

        

        subroutine CheckNodeInArea(IDMovingMesh, InitialNodeCoordinates, CheckArea, BoxAreaCoordinates)
        !*************************************************************************************
        !    SUBROUTINE: CheckNodeInArea
        ! 
        !    DESCRIPTION:
        !>   Determines which node belongs to which area - particle storage or fixed mesh area.
        !>   Nodes on the boundary between storage area 1 and the fixed mesh area are identified 
        !>   by a value of '4', nodes on the boundary between storage area 2 and the fixed mesh 
        !>   area are identified by a value of '5'. A value of '6' identifies fixed mesh nodes, 
        !>   if they are all located on a plane.
        !
        !>   @note : 
        !
        !>   @param[in] IDMovingMesh : ID of direction of moving mesh
        !>   @param[in] InitialNodeCoordinates : Initial nodal coordinates
        !>   @param[in] CheckArea : Checked area
        !>   @param[in] BoxAreaCoordinates : coordinates of the corner nodes of the mesh areas
        !
        !*************************************************************************************        
        
          implicit none

          integer(INTEGER_TYPE), intent(in) :: IDMovingMesh, CheckArea
          real(REAL_TYPE), dimension(:, :), intent(in) :: InitialNodeCoordinates
          real(REAL_TYPE), dimension(:, :), intent(in) :: BoxAreaCoordinates
          ! Local variables
          integer(INTEGER_TYPE) :: I, IDArea
          real(REAL_TYPE), dimension(NVECTOR) :: NodeCoordinates, BoxNodeCoordinates, BoxMinimum, BoxMaximum

          BoxMinimum = 1.E200
          BoxMaximum = -1.E200

          IDArea = GetAreaID(CheckArea)
          
          ! determine bounding box of the considered area
          do I = 1, (NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh)! loop dimension
            BoxNodeCoordinates(:) = BoxAreaCoordinates(I, :) ! extract single node from the array
            call CheckMinMax(BoxNodeCoordinates, BoxMinimum, BoxMaximum)
          end do

          do I = 1, size(AreaNodeAssignment, 2) ! loop total number of nodes

            NodeCoordinates(:) = InitialNodeCoordinates(I, :)
           
            if ( IsNodeInside(NodeCoordinates, BoxMinimum, BoxMaximum) ) then
              AreaNodeAssignment(IDMovingMesh, I) = AreaNodeAssignment(IDMovingMesh, I) + IDArea
            end if
                        
          end do
        
        end subroutine CheckNodeInArea
        


        logical function CheckIsModifiedMeshAdjustment(IDMovingMesh)
        !*************************************************************************************
        !    FUNCTION: CheckIsModifiedMeshAdjustment
        ! 
        !    DESCRIPTION:
        !>   Returns .true. if the order of the areas is: fixed, storage, storage instead of:
        !>   storage, fixed, storage.
        !
        !>   @param[in] IDMovingMesh : ID of direction of moving mesh 
        !        
        !>   @return CheckIsModifiedMeshAdjustment : 
        !
        !*************************************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IDMovingMesh
          ! local variables
          integer(INTEGER_TYPE) :: I
        
          CheckIsModifiedMeshAdjustment = .false.
        
          do I = 1, size(AreaNodeAssignment, 2) ! loop total number of nodes
            if ( AreaNodeAssignment(IDMovingMesh, I) == STORAGE1STORAGE2 ) then
              CheckIsModifiedMeshAdjustment = .true.
              EXIT
            end if
          end do
        
        end function CheckIsModifiedMeshAdjustment

        

        logical function IsNodeInside(Node, Minimum, Maximum)
        !*************************************************************************************
        !    FUNCTION: IsNodeInside
        ! 
        !    DESCRIPTION:
        !>   Checks whether node lies inside or on the box spanned by <Minimum, Maximum>
        !
        !>   @param[in] Node(:) : coordinates of the node that is checked
        !>   @param[in] Minimum(:) : minimum coordinates of a bounding box
        !>   @param[in] Maximum(:) : maximum coordinates of a bounding box
        !
        !>   @return IsNodeInside : .true. if inside or on boundary, .false. if outside
        !
        !*************************************************************************************
        implicit none

          real(REAL_TYPE), intent(in), dimension(:) :: Node
          real(REAL_TYPE), intent(in), dimension(:) :: Minimum, Maximum
          ! local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: Offset
          
          Offset = 1.e-10
        
          IsNodeInside = .true.
          do I = 1, size(Node)
            IsNodeInside = IsNodeInside .and. ( ( Node(I) >= (Minimum(I) - Offset) ) .and. ( Node(I) <= (Maximum(I) + Offset) ) )                
          end do
          
        end function IsNodeInside
        


        integer(INTEGER_TYPE) function GetAreaID(IDArea)
        !*************************************************************************************
        !    FUNCTION: GetAreaID
        ! 
        !    DESCRIPTION:
        !>   Returns the area ID defined as a parameter above for the defined area.
        !
        !>   @param[in] IDArea : = 1 --> Check particle storage area 1
        !>                       = 2 --> Check particle storage area 2
        !>                       = 3 --> Check fixed mesh area
        !
        !>   @return GetAreaID : ID of the area
        !
        !*************************************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IDArea
          
          select case(IDArea)
            case(1)
              GetAreaID = STORAGE1
            case(2)
              GetAreaID = STORAGE2
            case(3)
              GetAreaID = FIXEDMESH
          end select
        
        end function GetAreaID


      end module ModMPMMeshAdjustment
