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
      !     $Date: 2022-04-14 14:56:02 +0200 (Thu, 14 Apr 2022) $
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

        integer(INTEGER_TYPE), dimension(:, :, :), allocatable :: AreaNodeAssignment ! Stores for each node to which area it belongs --> increased rank for multipatch
        !integer(INTEGER_TYPE), dimension(:, :, :), allocatable :: AreaNodeAssignment ! Stores for each node to which area it belongs --> increased rank for multipatch

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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
        
          allocate(MeshMovement(NVECTOR)) ! this describes the direction of the moving mesh
          
          call DestroyMeshAdjustment() ! deallocates moving mesh variables 
          
          
          !------------------------------------------------------------------------------
          ! Part 1: Find node structures which correspond to the moving mesh
          
          
          if (CalParams%ApplyMeshSmoothing) then ! activate mesh adjustment ! --> we do not necessarily have a flag for this
            
              allocate(UnrotatedNodalCoordinates(1, 1), stat = IError)
              UnrotatedNodalCoordinates = 0.0

              !allocate(AreaNodeAssignment(CalParams%MovingMesh%NMovingMeshDirections, Counters%NodTot), stat = IError)          
              ! NEED TO INCLUDE MULTIPATCH GENERALIZATION --> I moved this here
              allocate(AreaNodeAssignment(CalParams%MovingMesh%NMovingMeshDirections(IPatch_Temporary), NURBS%Maximum_NControlPoints, Counters%NPatches), stat = IError)! , Counters%NPatches ! 
              AreaNodeAssignment = 0
          
             
              do IPatch = 1, Counters%NPatches ! loop over patches

            if (CalParams%MovingMesh%MovingMaterialID(IPatch)>0) then ! Search all nodes that belong to volume elements
                                                              ! with material StructureMaterialID
              
                call DetermineStructureNodes(IPatch) ! determine nodes of the moving mesh
              
            else if (CalParams%ApplyContactAlgorithm.and.CalParams%ApplyContactMeshBoundary) then ! need to make sure this works with the contact algorithm
              CalParams%MovingMesh%NStructureNodes(IPatch) = NInterfaceNodes
              allocate(CalParams%MovingMesh%StructureNodes(CalParams%MovingMesh%NStructureNodes(IPatch), Counters%NPatches), stat = IError)
              I = 1
              !do IPatch = 1, Counters%NPatches ! loop over patches
              do INode = 1, Counters%NodTot(IPatch)!Counters%NodTot
                  
                  GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
                  
                if (InterfaceNodes(INode)) then
                  CalParams%MovingMesh%StructureNodes(I, IPatch) = INode
                  I = I + 1
                end if
              end do
              
              !end do ! loop over patches

            end if ! Search all nodes that belong to volume elements
            
            !------------------------------------------------------------------------------
            ! Part 2: 
            
            !do IPatch = 1, Counters%NPatches ! loop over patches
            
                if (CalParams%MovingMesh%NStructureNodes(IPatch)>0) then ! if moving mesh structure
                    ! Determine to which area (particle storage or fixed mesh area) each node belongs
                    !do IPatch = 1, Counters%NPatches ! I added this loop here (3 Sept 2024)

                    call DetermineAreaNodeAssignment(NodalCoordinates(:,:, IPatch), IPatch) !, IPatch)    

                    !end do 
                
              
                    do K = 1, CalParams%MovingMesh%NMovingMeshDirections(IPatch)
                        IsModifiedMeshAdjustment(K) = CheckIsModifiedMeshAdjustment(K, IPatch)             
                    end do

                else ! if not a moving mesh structure
              
                    ! I disabled this for now 
                    !CalParams%ApplyMeshSmoothing = .false.

                end if ! if moving mesh structure 

            end do ! loop over patches
            

          else ! if moving mesh is not used
              
            allocate(UnrotatedNodalCoordinates(1, 1), stat = IError)
            UnrotatedNodalCoordinates = 0.0
            
            
            
          
          
            
            
          
          end if ! if moving mesh structure
          
          
          
          !end if !activate mesh adjustment
          

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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          
        IPatch_Temporary = 1
        
        
        ! Moving mesh algorithm breakdown: 
        ! 1- determine displacement for the entire patch of the moving mesh (non-deformable)
        ! 2- determine the storage area length so that we can scale the compression or extension (deformable)
        ! 3- move fixed mesh in a normal way by just adding the DUVector
        ! 4- 
        
        
        

         do IPatch = 1, Counters%NPatches   ! loop over patches
             
             ! we move the entire mesh according to the populated DUVector
             call DetermineDisplacementVector(DUVector, ConsideredDirections, IPatch) !--> DUVector needs to be for the moving mesh
          
             do I = 1, CalParams%MovingMesh%NMovingMeshDirections(IPatch)
          
                 ! dynamic length changes 
                 allocate(InitialStorageAreaLength(CalParams%MovingMesh%NStorageAreas(I, IPatch)), stat = IError) !, Counters%NPatches

                 ! this involves calculating the distance between two planes
                 call DetermineStorageAreaLength(I, ConsideredDirections(I), NodalCoordinates(:,:,IPatch), InitialStorageAreaLength, IPatch)
         
			
                 ! Move the moving mesh with the ground motion. Current implementation is only for 1D movement 
                 if (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_FILE) then 
                     DUVector(I) = CalParams%Multipliers%VelocitySolidCurrent*CalParams%TimeIncrement 
                 end if 

            
                 if (IsModifiedMeshAdjustment(I)) then
                     
                     call MoveModifiedMeshAdjustment(DUVector, I, ConsideredDirections(I), IPatch)
                     
                 else
                     
                     call MoveFixedMesh(DUVector, I, ConsideredDirections(I), IPatch)
                     ! --> this is for the moving mesh adjustment
                     
                 end if


                 
                 call AdjustParticleStorageAreas(DUVector, I, ConsideredDirections(I), InitialStorageAreaLength, IPatch) 
                 ! --> this is for the deforming mesh adjustment
     
                 deallocate(InitialStorageAreaLength, stat = IError)
     
             end do
          
          
          end do ! loop over patches

        end subroutine MeshAdjustment
        
        
        
        
        !---------------------------------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------------------------------
        ! I need to breakdown the mesh adjustment to moving mesh and storage areas
        
        subroutine MeshAdjustment_MovingMesh_Multipatch(IPatch_MovingMesh, IPatch_DeformingMesh)
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
          
          ! Multipatch variables 
          !integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch_MovingMesh ! only one moving mesh
          integer(INTEGER_TYPE), intent(in) :: IPatch_DeformingMesh ! 2 correspond to compressing and extending mesh
          
        !IPatch_Temporary = 1
        
        
        ! Moving mesh algorithm breakdown: 
        ! 1- determine displacement for the entire patch of the moving mesh (non-deformable)
        ! 2- determine the storage area length so that we can scale the compression or extension (deformable)
        ! 3- move fixed mesh in a normal way by just adding the DUVector
        ! 4- 
        
        
        

        ! let's make IPatch an input so that we can move the mesh according to the prescribed velocity of the MPs
         !do IPatch = 1, Counters%NPatches   ! loop over patches
             
             ! we move the entire mesh according to the populated DUVector
             call DetermineDisplacementVector(DUVector, ConsideredDirections, IPatch_MovingMesh) !--> DUVector needs to be for the moving mesh
             
             
             do I = 1, CalParams%MovingMesh%NMovingMeshDirections(IPatch_MovingMesh)
          
                 ! dynamic length changes 
                 allocate(InitialStorageAreaLength(CalParams%MovingMesh%NStorageAreas(I, IPatch_DeformingMesh)), stat = IError) !, Counters%NPatches

                 ! this involves calculating the distance between two planes
                 call DetermineStorageAreaLength(I, ConsideredDirections(I), NodalCoordinates(:,:,IPatch_DeformingMesh), InitialStorageAreaLength, IPatch_DeformingMesh)
         
			
                 ! Move the moving mesh with the ground motion. Current implementation is only for 1D movement 
                 if (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_FILE) then 
                     DUVector(I) = CalParams%Multipliers%VelocitySolidCurrent*CalParams%TimeIncrement 
                 end if 

            
                 if (IsModifiedMeshAdjustment(I)) then
                     
                     call MoveModifiedMeshAdjustment(DUVector, I, ConsideredDirections(I), IPatch_MovingMesh)
                     
                 else
                     
                     call MoveFixedMesh(DUVector, I, ConsideredDirections(I), IPatch_MovingMesh)
                     ! --> this is for the moving mesh adjustment
                     
                 end if


                 
                 call AdjustParticleStorageAreas(DUVector, I, ConsideredDirections(I), InitialStorageAreaLength, IPatch_DeformingMesh) 
                 ! --> this is for the deforming mesh adjustment
     
                 deallocate(InitialStorageAreaLength, stat = IError)
     
             end do
          
          
          !end do ! loop over patches

        end subroutine MeshAdjustment_MovingMesh_Multipatch
        
        
        
        
        
        
        
        !---------------------------------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------------------------------


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

          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          
          if (allocated(DoFMovingMeshStructure)) then
            deallocate(DoFMovingMeshStructure, stat = IError)
          end if

          allocate(DoFMovingMeshStructure(Counters%N), stat = IError)
          DoFMovingMeshStructure = 0
  
          if (CalParams%ApplyContactAlgorithm.and.CalParams%ApplyContactMeshBoundary) then
              do IPatch = 1, Counters%NPatches ! loop over patches
            do NodeID = 1, Counters%NodTot(IPatch)!Counters%NodTot
                
                GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(NodeID,IPatch)
                
              if (InterfaceNodes(NodeID)) then
                DoFMovingMeshStructure(ReducedDof(GlobalNodeID) + 2) = 1
              end if
            end do
            
            end do ! loop over patches
          else
              do IPatch = 1, Counters%NPatches ! loop over patches
            do IEl = 1, Counters%NEl(IPatch) !NURBS%nel_NURBS(IPatch)!Counters%NEl ! loop over elements
              do IPart = 1, NPartEle(IEl, IPatch)
                ParticleIndex = GetParticleIndex(IPart, IEl, IPatch)
                if (EntityIDArray(ParticleIndex)==HARD_ENTITY) then
                  do INode = 1, ELEMENTNODES
                    NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(ElementConnectivities(INode, IEl, IPatch), IPatch)
                    DoFMovingMeshStructure(ReducedDof(NodeID) + 2) = 1
                  end do
                end if
              end do
            end do ! loop over elements
            end do ! loop over patches
          end if
          
        end subroutine DetermineDoFMovingMeshStructure
                  
        subroutine DetermineDisplacementVector(DUVector, ConsideredDirections, IPatch)
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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE), intent(in) :: IPatch ! input
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          ! Determine sum of incremental displacements
          SummedStructureDisplacement = 0.0
          !do INode = 1, CalParams%MovingMesh%NStructureNodes
          !  do I = 1, NVECTOR
          !    SummedStructureDisplacement(I) = SummedStructureDisplacement(I) +  &
          !      (NodalCoordinatesUpd(CalParams%MovingMesh%StructureNodes(INode), I) -  &
          !       NodalCoordinates(CalParams%MovingMesh%StructureNodes(INode), I) )
          !  end do
          !end do
          !DUVector = SummedStructureDisplacement / CalParams%MovingMesh%NStructureNodes
          
          !do IPatch = 1, Counters%NPatches ! start loop over patches
          
          do INode = 1, CalParams%MovingMesh%NStructureNodes(IPatch) ! loop over all structured nodes in each patch
            do I = 1, NVECTOR ! loop over each direction
              
                GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode, IPatch)
                
              
                ! this needs to happen to the structure which is carrying the moving mesh structure, but 
                ! we should not include the compressing mesh nodes 
                SummedStructureDisplacement(I) = SummedStructureDisplacement(I) +  &
                    (NodalCoordinatesUpd(CalParams%MovingMesh%StructureNodes(INode, IPatch), I, IPatch) -  &
                    NodalCoordinates(CalParams%MovingMesh%StructureNodes(INode, IPatch), I, IPatch) )
              
            end do ! loop over each direction
          end do ! loop over all structured nodes
          
          !end do ! end loop over patches
          
          
          if (CalParams%MovingMesh%NStructureNodes(IPatch) .le. 0) RETURN 
          DUVector = SummedStructureDisplacement / CalParams%MovingMesh%NStructureNodes(IPatch)
          !end if 
          
          if(CalParams%MovingMesh%MovingMeshDirection(IPatch).gt.0)then !moving mesh direction has been defined
              
            ShiftDirection = CalParams%MovingMesh%MovingMeshDirection(IPatch)
            DisplacementValue = DUVector(ShiftDirection)
            DUVector = 0.0
            DUVector(ShiftDirection) = DisplacementValue
            ConsideredDirections(1) = ShiftDirection
            ConsideredDirections(2) = 0
              
          else if (CalParams%MovingMesh%NMovingMeshDirections(IPatch)==1) then
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
        
        subroutine MoveModifiedMeshAdjustment(DUVector, IMovingMesh, Direction, IPatch)
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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE), intent(in) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
        
          !do IPatch = 1, Counters%NPatches ! loop over patches
          
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            if (AreaNodeAssignment(IMovingMesh, I, IPatch)==STORAGE1STORAGE2) then ! Move node by DUVector !, IPatch
                NodalCoordinates(I, Direction, IPatch) = NodalCoordinates(I, Direction, IPatch) + DUVector(Direction)
            end if
          
          end do
          
          !end do ! loop over patches
        
        end subroutine MoveModifiedMeshAdjustment
        
        subroutine MoveFixedMesh(DUVector, IMovingMesh, Direction, IPatch)
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
          
          ! Multipatch variables 
          !integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          !IPatch_Temporary = 1
        
          !do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over the control points for each patch
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            
              if ( (AreaNodeAssignment(IMovingMesh, I, IPatch)==FIXEDMESH).or. & !, IPatch
                 (AreaNodeAssignment(IMovingMesh, I, IPatch)==STORAGE1FIXEDMESH).or. & !, IPatch
                 (AreaNodeAssignment(IMovingMesh, I, IPatch)==STORAGE2FIXEDMESH).or. & !, IPatch
                 (AreaNodeAssignment(IMovingMesh, I, IPatch)==PLANEFIXEDMESH) ) then ! Node belongs to the fixed mesh area, including its boundary !, IPatch
              
                
                  NodalCoordinates(I, Direction, IPatch) = NodalCoordinates(I, Direction, IPatch) + DUVector(Direction)
            
                 end if
          
          end do ! loop over the control points for each patch
          !end do
        
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

          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

     
          do IPatch = 1, Counters%NPatches ! loop over all patches
              
          do I = 1, CalParams%MovingMesh%NStructureNodes(IPatch) ! Loop over structure nodes
            NodeID = CalParams%MovingMesh%StructureNodes(I, IPatch)
            
            GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(NodeID, IPatch)
            
            if (.not.IsModifiedMeshAdjustment(IMovingMesh)) then ! NodeCoord(I) - DUVector
              NodalCoordinates(NodeID, Direction, IPatch) = NodalCoordinates(NodeID, Direction, IPatch) - DUVector(Direction)
            end if
            
            ! NodeCoord(I) + DUTot(I)
            !IDoF = ReducedDof(NodeID) + Direction
            IDoF = ReducedDof(GlobalNodeID) + Direction
            
            NodalCoordinates(NodeID, Direction, IPatch) = NodalCoordinates(NodeID, Direction, IPatch) + DUTot(IDoF)
            
          end do ! Loop over structure nodes
          
          end do ! loop over all patches
     
        end subroutine AdjustStructureNodes
        
        subroutine AdjustParticleStorageAreas(DUVector, IMovingMesh, ShiftDirection, InitialStorageAreaLength, IPatch)
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
          
          !integer(INTEGER_TYPE) :: IPatch_Temporary
          
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
          real(REAL_TYPE), dimension(CalParams%MovingMesh%NStorageAreas(IMovingMesh,IPatch)), intent(in) :: InitialStorageAreaLength !IPatch_Temporary --> needs generalization 
          
          ! Local variables
          integer(INTEGER_TYPE) :: AreaID
          
          ! temporary initialization 
          !IPatch_Temporary = 1
          
          do AreaID = 1, CalParams%MovingMesh%NStorageAreas(IMovingMesh, IPatch) ! Loop over storage areas !, IPatch
              
            if (InitialStorageAreaLength(AreaID)>1E-5) then ! Storage area is not defined as plane
              call ShiftStorageAreaNodes(DUVector, AreaID, IMovingMesh, ShiftDirection, InitialStorageAreaLength(AreaID), IPatch)
            end if
            
          end do ! Loop over storage areas
        
        end subroutine AdjustParticleStorageAreas

        
        
        
        
        
        subroutine DetermineStorageAreaLength(IMovingMesh, ShiftDirection, Coordinates, InitialStorageAreaLength, IPatch)
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
          
          ! input multipatch generalization 
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
          real(REAL_TYPE), dimension(Counters%NodTot(IPatch), NVECTOR), intent(in) :: Coordinates
          
          
          
          real(REAL_TYPE), dimension(CalParams%MovingMesh%NStorageAreas(IMovingMesh,IPatch)), intent(out) :: InitialStorageAreaLength ! FIXED! --> HARDCODED and needs multipatch generalization
          
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, AreaID
          real(REAL_TYPE) :: Length
          
          integer(INTEGER_TYPE) :: IPatch_Temporary
          
          
          IPatch_Temporary = 1
          
          InitialStorageAreaLength = -1.0
        
          ! This has to do with the moving mesh
          
          do AreaID = 1, CalParams%MovingMesh%NStorageAreas(IMovingMesh,IPatch) 
          
            do I = 1, CalParams%MovingMesh%NAreaNodes(IPatch)
              
                do J = 1, CalParams%MovingMesh%NAreaNodes(IPatch)
                
                    if (I/=J) then
                  
                        Length =  &
                            abs(Coordinates(CalParams%MovingMesh%MeshAreas(IMovingMesh, AreaID, I, IPatch), ShiftDirection) - &
                            Coordinates(CalParams%MovingMesh%MeshAreas(IMovingMesh, AreaID, J, IPatch), ShiftDirection) )
                  
                        if (Length>InitialStorageAreaLength(AreaID) ) then
                    
                            InitialStorageAreaLength(AreaID) = Length
                  
                        end if
                
                    end if
              
                end do
            
            end do
            
          end do
          
        end subroutine DetermineStorageAreaLength
        
        
        

        subroutine ShiftStorageAreaNodes(DUVector, IArea, IMovingMesh, ShiftDirection, InitialStorageAreaLength, IPatch)
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
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE), intent(in) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
        
          AreaID = GetAreaID(IArea)
        
          !do IPatch = 1, Counters%NPatches ! loop over patches
              
          do INode = 1, Counters%NodTot(IPatch)!Counters%NodTot
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)

            if (AreaNodeAssignment(IMovingMesh, INode, IPatch)==AreaID) then !, IPatch

            ! Determine distance between the fixed boundary and the considered node
            BoundaryNodeDistance =  &
               abs(NodalCoordinates(INode, ShiftDirection, IPatch) -  &
                      NodalCoordinates(CalParams%MovingMesh%MeshAreas(IMovingMesh, IArea, 1, IPatch), &
                      ShiftDirection, IPatch) )


              ! Loop over nodes of considered storage area that are not on the fixed boundary
              if (BoundaryNodeDistance>1E-5) then
                
                ! Shift the considered node
                NodalCoordinates(INode, ShiftDirection, IPatch) = &
                  NodalCoordinates(INode, ShiftDirection, IPatch) + &
                  DUVector(ShiftDirection) / InitialStorageAreaLength * &
                  BoundaryNodeDistance 
                ! the (BoundaryNodeDistance/InitialStorageAreaLength) is a ratio to 
                ! take into account the distance of each node from the fixed boundary

              end if
            end if
          end do
          
          !end do ! loop over patches

        end subroutine ShiftStorageAreaNodes

        subroutine DetermineStructureNodes(IPatch)
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
      
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
          ! Determine nodes that belong to structure volume elements
          !allocate(NodeCheck(Counters%NodTot), stat = IError)
          allocate(NodeCheck(Counters%NodTot(IPatch)), stat = IError) !, Counters%NPatches
          NodeCheck = 0
          
          !do IPatch = 1, Counters%NPatches ! loop over patches
              do IEl = 1, Counters%NEl(IPatch) !NURBS%nel_NURBS(IPatch) !Counters%NEl
            
                  if (IsActiveElement(IEl,IPatch)) then ! if active element
                      if (NPartEle(IEl, IPatch)>0) then ! Loop over all activated elements that contain particles
                
                          ParticleIndex = GetParticleIndex(1, IEl, IPatch)
                
                          if (MaterialIDArray(ParticleIndex)== & ! if material ID matches moving mesh material ID
                              CalParams%MovingMesh%MovingMaterialID(IPatch)) then ! Found an element that is part of the structure
                              do INode = 1, ELEMENTNODES ! loop over element nodes
                                  
                                  NodeID = ElementConnectivities(INode, IEl, IPatch) ! --> based on the local node ID in the patch
                                  NodeCheck(NodeID) = 1 !, IPatch
                  
                              end do ! loop over element nodes
  
                              end if ! if material ID matches moving mesh material ID
                
                      end if ! Loop over all activated elements that contain particles
                
                  end if ! if active element
          
              end do ! loop over elements
          
          !end do ! loop over patches
          
          ! this subroutine populates the moving mesh structure global variables generalized for multiple patches
          call EvaluateNodeCheck(NodeCheck, IPatch)
        
          deallocate(NodeCheck, stat = IError)
          
        end subroutine DetermineStructureNodes

        
        
        
        
        subroutine EvaluateNodeCheck(NodeCheck, IPatch)
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

          !integer(INTEGER_TYPE), dimension(Counters%NodTot), intent(in) :: NodeCheck
          integer(INTEGER_TYPE), intent(in) :: IPatch
          integer(INTEGER_TYPE), dimension(Counters%NodTot(IPatch)), intent(in) :: NodeCheck ! this has value of 1 if part of the moving mesh !Maximum_NControlPoints !, Counters%NPatches
          ! Local variables
          integer(INTEGER_TYPE) :: I, INode, IError, count
          
          ! Multipatch variables 
          
          integer(INTEGER_TYPE) :: GlobalNodeID
        
          ! Determine number of structure nodes
          !allocate(CalParams%MovingMesh%NStructureNodes(Counters%NPatches)) ! allocate NStructureNodes according to the number of patches
          CalParams%MovingMesh%NStructureNodes(IPatch) = 0 ! initialize to zero
          
          !do IPatch = 1, Counters%NPatches ! loop over patches
          
              do INode = 1, Counters%NodTot(IPatch) ! loop over nodes in the patch !Counters%NodTot
              
                  GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch) ! obtain global node ID
            
                  if (NodeCheck(INode)==1) then ! if part of the moving mesh portion !, IPatch
               
                      CalParams%MovingMesh%NStructureNodes(IPatch) = CalParams%MovingMesh%NStructureNodes(IPatch) + 1 ! increase count
            
                  end if
          
              end do ! loop over nodes in patch
          
          !end do ! loop over patches
          
          ! I need to find the maximum NStructureNodes
          ! For allocation purposes, we need the maximum number of structure nodes
          count = CalParams%MovingMesh%NStructureNodes(IPatch) ! Alsardi: Multipatch moving mesh  !maxval(
          
          ! Determine array of structure nodes
          !count = CalParams%MovingMesh%NStructureNodes
          
          if (.not. allocated(CalParams%MovingMesh%StructureNodes) ) then !--> this needs generalization 
              allocate( CalParams%MovingMesh%StructureNodes(count, Counters%NPatches), stat = IError) ! allocate structured nodes so that we can store the local ID of each node 
          end if 
          
          CalParams%MovingMesh%StructureNodes(:,IPatch) = 0 ! initialize numbering
          I = 1 ! initialize count
          
          !do IPatch = 1, Counters%NPatches ! loop over patches
         
              do INode = 1, Counters%NodTot(IPatch) ! loop over number of nodes in a patch !Counters%NodTot
              
                  ! obtain global node ID 
                  GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
         
                  if (NodeCheck(INode)==1) then ! check if node belongs to the moving mesh structure !, IPatch
                      CalParams%MovingMesh%StructureNodes(I, IPatch) = INode ! store the local node numbering 
                      I = I + 1 ! increase count
                  end if

              end do ! loop over nodes 
          
              !I = 1
          
          !end do ! loop over patches
          
        
        end subroutine EvaluateNodeCheck

        

        subroutine DetermineAreaNodeAssignment(InitialNodeCoordinates, IPatch)
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
          real(REAL_TYPE), dimension(CalParams%MovingMesh%NAreaNodes(1), NVECTOR) :: BoxAreaCoordinates ! IPatch_Temporary --> HARDCODED needs generalization 

          integer(INTEGER_TYPE) :: IPatch_Temporary
          
          ! Multipatch generalization 
          integer(INTEGER_TYPE), intent(in) :: IPatch 
          
          ! Temporary initialization 
          IPatch_Temporary = 1
          
          
          
          
          !allocate(AreaNodeAssignment(CalParams%MovingMesh%NMovingMeshDirections, Counters%NodTot), stat = IError)          
          ! NEED TO INCLUDE MULTIPATCH GENERALIZATION --> I commented this here
          !allocate(AreaNodeAssignment(CalParams%MovingMesh%NMovingMeshDirections(IPatch), Maximum_NControlPoints, Counters%NPatches), stat = IError)! , Counters%NPatches ! 
          !AreaNodeAssignment = 0
          BoxAreaCoordinates = 0.0
        
          do K = 1, CalParams%MovingMesh%NMovingMeshDirections(IPatch)
          
            ! IDArea = 1: Check particle storage area 1
            ! IDArea = 2: Check particle storage area 2
            do IDArea = 1, CalParams%MovingMesh%NStorageAreas(K, IPatch) ! deforming mesh  !, IPatch

              do I = 1, size(BoxAreaCoordinates, 1) ! loop corner nodes of the areas
                do J = 1, size(BoxAreaCoordinates, 2) ! loop dimensions
                  BoxAreaCoordinates(I, J) = InitialNodeCoordinates(CalParams%MovingMesh%MeshAreas(K, IDArea, I, IPatch), J)
                end do  
              end do
            
              call CheckNodeInArea(K, InitialNodeCoordinates, IDArea, BoxAreaCoordinates, IPatch)
            end do

            ! IDArea = 3: Check fixed mesh area
            IDArea = 3
            do I = 1, size(BoxAreaCoordinates, 1) ! loop corner nodes of the areas
              do J = 1, size(BoxAreaCoordinates, 2) ! loop dimensions
                  ! if patch does not have a fixed area then exit
                  if ( CalParams%MovingMesh%MeshAreas(K, IDArea, I, IPatch) == -1 ) RETURN
                BoxAreaCoordinates(I, J) = InitialNodeCoordinates(CalParams%MovingMesh%MeshAreas(K, IDArea, I, IPatch), J)
              end do   
            end do
            
            call CheckNodeInArea(K, InitialNodeCoordinates, IDArea, BoxAreaCoordinates, IPatch)
        
          end do
          
        end subroutine DetermineAreaNodeAssignment

        

        subroutine CheckNodeInArea(IDMovingMesh, InitialNodeCoordinates, CheckArea, BoxAreaCoordinates, IPatch)
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

          ! Multipatch generalization 
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
          
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
              AreaNodeAssignment(IDMovingMesh, I, IPatch) = AreaNodeAssignment(IDMovingMesh, I, IPatch) + IDArea
            end if
                        
          end do
        
        end subroutine CheckNodeInArea
        


        logical function CheckIsModifiedMeshAdjustment(IDMovingMesh, IPatch)
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
          
          ! Multipatch generalization 
          integer(INTEGER_TYPE), intent(in) :: IPatch
        
          CheckIsModifiedMeshAdjustment = .false.
        
          do I = 1, size(AreaNodeAssignment, 2) ! loop total number of nodes
            if ( AreaNodeAssignment(IDMovingMesh, I, IPatch) == STORAGE1STORAGE2 ) then !, IPatch
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
