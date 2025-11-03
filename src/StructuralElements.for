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
    !   and soil�water�structure interaction using the material point method (MPM)
    !
    !	Copyright (C) 2020  Members of the Anura3D MPM Research Community 
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

    module ModStructuralElements
    
      !**********************************************************************
      !
      !    Function:  Contains definition of variables used for structural elements
      !             Including rigid and deformable motions. 
      !                     
      !*********************************************************************
    
      use ModGlobalConstants
      use ModCounters 
      use ModMeshInfo
      Use ModMPMData
      use ModGlobalData
      use ModMatrixMath
	  use ModGeometryMath
      
      implicit none           
        ! Counters
        integer(INTEGER_TYPE):: NRigidBodies !store how many bodies there are
        integer(INTEGER_TYPE), dimension(:), allocatable ::  QtyRigidBodies !store how many 1D elements are in each body
        character(len=100), dimension(:), allocatable ::  UniqueBodiesNameF !Store the names of each rigid body
        ! Lists

        real(REAL_TYPE), dimension(:,:), allocatable :: NormalsRigidBody
		real(REAL_TYPE):: Kfactor= sqrt(PI*2/9), kappa=0.0
		logical::Stiffness_assigned=.false.
        
        type Material
            real(REAL_TYPE):: ContactStiffness, FrictionCoef
        end type Material
        
        type Elements
            real(REAL_TYPE) :: RigidBody
        end type Elements
        
        type(Elements), dimension(:), allocatable :: ElementsRigidBody
        
        
        type Body
            integer(INTEGER_TYPE), dimension(:), allocatable :: Ids, Nodes
            real(REAL_TYPE), dimension(:,:), allocatable :: MechProp, Normals
            character(len=100) :: Name
            real(REAL_TYPE) :: Thickness, Density, Centroid(2), Inertia, Area, Velocity(2),&
                Acceleration(2), AngularVelocity, AngularAcceleration, Fsum(2), Tsum
            integer :: NumMaterials
            type(Material), dimension(:), allocatable :: Materials
        end type Body
        
        type(Body), dimension(:), allocatable :: DataStructureRigidBody
         
        type RBConditions
            logical :: ApplyPrescribedVelocity = .false., ApplyPrescribedRotation = .false., &
            ApplyPrescribedForce = .false., &
            ApplyPrescribedTorque = .false., &
            IsThereLinearForce = .false.
            integer(INTEGER_TYPE) :: ApplyCenterRotation = 0
            real(REAL_TYPE) :: PrescribedVelocity(2), &
                               VelConstraints(2)=0.0, &
                               PrescribedRotation, &
                               CenterRotation(2), &
                               PrescribedForce(2)=0.0, &
                               ForceConstraints(2)=0.0, &
                               PrescribedTorque=0.0, &
                               LinearInitial, LinearFinal
        end type RBConditions
        
         type(RBConditions), dimension(:), allocatable :: RigidBodyConditions
        
         !**********************************************************************
         contains

         subroutine CreateDataStructureRigidBody()
        !**********************************************************************
        !
        !    Function:  Create the Data Structure of the Rigid Body
        !
        !
        !    Note : DataStructureRigidBody is allocated in this subroutine
        !
        !**********************************************************************
        implicit none 
        !Variables
        integer(INTEGER_TYPE) :: C, D, E, J, I, N, K, L
        integer(INTEGER_TYPE) :: ND(SIZE(ThinRigidCoordinates,1))
        integer(INTEGER_TYPE) :: ElemOffset, PropId, PropIndex, NodeCount
        logical :: alreadyExists 
        allocate(DataStructureRigidBody(NRigidBodies))
        ElemOffset = SIZE(ElementConnectivities,2)

        ! Build per-body material, kinematic and geometric data from GiD input
        do J = 1, NRigidBodies
            ! Locate property row for current body and convert GiD id to local index
            PropId = -1
            PropIndex = 0
            do N = 1, NThinElmProp
                if (UniqueBodiesNameF(J)==ThinRigidElmProp(N,2)) then
                    PropIndex = N
                    read(ThinRigidElmProp(N,1), *) PropId
                    PropId = PropId - ElemOffset !Assumes GiD ids are offset by FE mesh count
                    exit
                end if
            end do
            if (PropIndex == 0) cycle
            if (PropId < 1 .or. PropId > NThinElmProp) then
                PropId = PropIndex
            end if
            
            !Initialize variables:: Note: You might want to change this depending on your needs
            DataStructureRigidBody(J)%Velocity=0.0	
            DataStructureRigidBody(J)%Acceleration=0.0
            DataStructureRigidBody(J)%Fsum=0.0
            DataStructureRigidBody(J)%Tsum=0.0
            DataStructureRigidBody(J)%AngularVelocity=0.0
            DataStructureRigidBody(J)%AngularAcceleration=0.0
			
			
            !This assigns initial imposed velocity
            do N=1, SIZE(ThinRigidElmInitVel,1)
                read(ThinRigidElmInitVel(N,1), *) K
                if (PropId == K - ElemOffset) then
                    read(ThinRigidElmInitVel(N,2), *) L
                    if (L==1) then
                        read(ThinRigidElmInitVel(N,3), *) DataStructureRigidBody(J)%Velocity(1)
                    end if
                    read(ThinRigidElmInitVel(N,4), *) L
                    if (L==1) then
                        read(ThinRigidElmInitVel(N,5), *) DataStructureRigidBody(J)%Velocity(2)
                    end if
                    exit
                end if
			end do
            
			!This assigns initial force acting on elements
            do N=1, SIZE(ThinRigidElmInitFor,1)
                read(ThinRigidElmInitFor(N,1), *) K
                if (PropId == K - ElemOffset) then
                    read(ThinRigidElmInitFor(N,2), *) L
                    if (L==1) then
                        read(ThinRigidElmInitFor(N,3), *) DataStructureRigidBody(J)%Fsum(1)
                    end if
                    read(ThinRigidElmInitFor(N,4), *) L
                    if (L==1) then
                        read(ThinRigidElmInitFor(N,5), *) DataStructureRigidBody(J)%Fsum(2)
                    end if
                    exit
                end if
            end do
            
            if (.not. ISAXISYMMETRIC) then !plane strain
				!This assigns initial rotation
                do N=1, SIZE(ThinRigidElmInitRot,1)
                read(ThinRigidElmInitRot(N,1), *) K
                if (PropId == K - ElemOffset) then
                    read(ThinRigidElmInitRot(N,2), *) L
                    if (L==1) then !If L ==0 then the center of rotation is different, we need to decide
                        read(ThinRigidElmInitRot(N,3), *) DataStructureRigidBody(J)%AngularVelocity
                    end if
                    exit
                end if
				end do
            
			!This assigns initial torque
            do N=1, SIZE(ThinRigidElmInitTor,1)
                read(ThinRigidElmInitTor(N,1), *) K
                if (PropId == K - ElemOffset) then

                    read(ThinRigidElmInitTor(N,2), *) DataStructureRigidBody(J)%Tsum

                    exit
                end if
            end do
               
            
			end if        
            
			!This part fills the thin rigid element data structures
            
            DataStructureRigidBody(J)%Name = ThinRigidElmProp(PropId,2)                  !Name or identifier of rigid body		
            read(ThinRigidElmProp(PropId,3), *) DataStructureRigidBody(J)%Thickness      !Thickness associated to R.B.
            read(ThinRigidElmProp(PropId,4), *) DataStructureRigidBody(J)%Density        !RB material density
            read(ThinRigidElmProp(PropId,5), *) DataStructureRigidBody(J)%NumMaterials   !Number of possible materials in contact (#of soils in contact up to 4 in current version)
            allocate(DataStructureRigidBody(J)%Materials(DataStructureRigidBody(J)%NumMaterials))
            do I = 1, DataStructureRigidBody(J)%NumMaterials !HERE I CAN CALCULATE STIFFNESS DIRECTLY FOR STABILITY
                read(ThinRigidElmProp(PropId,5 + 2* I -1), *) DataStructureRigidBody(J)%Materials(I)%ContactStiffness !Stiffness with respect to material I
                read(ThinRigidElmProp(PropId,5 + 2* I), *) DataStructureRigidBody(J)%Materials(I)%FrictionCoef		 !Friction coefficient with respect to mat I
				!HERE WE COULD ADD ADDITIONAL USER DEFINED PROPERTIES
            end do
            
            allocate(integer(INTEGER_TYPE) :: DataStructureRigidBody(J)%Ids(QtyRigidBodies(J)))
            allocate(real(REAL_TYPE) :: DataStructureRigidBody(J)%MechProp(QtyRigidBodies(J),5))
            
            !assign IDS of each rigid body (MOVING DOWN THIS CODE IS A MESS::NEEDS TO BE IMPROVED)
            C=1
            do N = 1, SIZE(ThinRigidElmProp,1) !These IDs should match ThinElement connectivity index
                if (UniqueBodiesNameF(J)==ThinRigidElmProp(N,2)) then
                    read(ThinRigidElmProp(N,1), *) L
                    L = L - ElemOffset
                    DataStructureRigidBody(J)%Ids(C)= L !PROBLEM: Assumes something about the ordering of EL numbers given by GID
                    C= C+1
                end if
            end do
            
            !assign NODES belongs to rigid body
            ND = 0
            ND(1) = ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(1),1)
            ND(2) = ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(1),2)
            NodeCount = 2
            do I = 2, QtyRigidBodies(J) !loops elements in rigid body (this could be quantityRigidBodies)
                 do E = 1, 2 !loops nodes in element
                     alreadyExists = .false.
                     do D = 1, NodeCount
                         if (ND(D) == ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),E)) then !node already listed
                             alreadyExists = .true.
                             exit
                         end if
                     end do
                     if (.not. alreadyExists) then !node not listed so add it
                         NodeCount = NodeCount + 1
                         ND(NodeCount) = ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),E)
                     end if
                 end do
            end do
            allocate(integer(INTEGER_TYPE) :: DataStructureRigidBody(J)%Nodes(NodeCount))
			!transcribe node IDs to structure
            DataStructureRigidBody(J)%Nodes = ND           
            
				
             DataStructureRigidBody(J)%Area = 0.0
             DataStructureRigidBody(J)%Centroid = 0.0	
			 
             do I = 1, QtyRigidBodies(J) !Loop bar elements in RB
				 !stores area in MechProp 1 for that 1D bar element
                 DataStructureRigidBody(J)%MechProp(I,1)=& 
                     SQRT((ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),2),1)- &
                     ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),1),1))**2+ &
                     (ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),2),2)- &
                     ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),1),2))**2)* &
                     DataStructureRigidBody(J)%Thickness !Area of element
				 
				 !accumulate area
				 DataStructureRigidBody(J)%Area = DataStructureRigidBody(J)%Area + &
                     DataStructureRigidBody(J)%MechProp(I,1)!Total area of body
				 
				if (.not. ISAXISYMMETRIC) then !For axisymmetric calcs. the centroid is not important
                 DataStructureRigidBody(J)%MechProp(I,2)=&
                     (ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),2),1) + &
                     ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),1),1))/2 !x centroid in element
                 DataStructureRigidBody(J)%MechProp(I,3)=&
                     (ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),2),2)+&
                      ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),1),2)   )/2 !y centroid in element
			
			 !Accumulate area moments
				DataStructureRigidBody(J)%Centroid = DataStructureRigidBody(J)%Centroid + &
                     DataStructureRigidBody(J)%MechProp(I,2:3)*   DataStructureRigidBody(J)%MechProp(I,1)
			  end if
			 end do 
			 

			 !Calculates centroid
             DataStructureRigidBody(J)%Centroid =  DataStructureRigidBody(J)%Centroid / DataStructureRigidBody(J)%Area
             !Calculate  Normals and Inertia
             allocate(real(REAL_TYPE) :: DataStructureRigidBody(J)%Normals(QtyRigidBodies(J),2))
             DataStructureRigidBody(J)%Inertia = 0.0
             do I =1,  QtyRigidBodies(J)
			    if (.not. ISAXISYMMETRIC) then !For axisymmetric calcs. the centroid is not important
                    DataStructureRigidBody(J)%MechProp(I,4) = & !Inertia of each element
                        ( (1/12)* DataStructureRigidBody(J)%MechProp(I,1) * DataStructureRigidBody(J)%Density * &
                        ((DataStructureRigidBody(J)%MechProp(I,1)/DataStructureRigidBody(J)%Thickness)**2 + &
                        (DataStructureRigidBody(J)%Thickness)**2) )
                    DataStructureRigidBody(J)%MechProp(I,5) = &!Square of distance between local centroid and global centroid
                        (DataStructureRigidBody(J)%Centroid(1)-DataStructureRigidBody(J)%MechProp(I,2))**2 + & 
                        (DataStructureRigidBody(J)%Centroid(2)-DataStructureRigidBody(J)%MechProp(I,3))**2
                    DataStructureRigidBody(J)%Inertia = DataStructureRigidBody(J)%Inertia + &
                         DataStructureRigidBody(J)%MechProp(I,4) + &
                             DataStructureRigidBody(J)%MechProp(I,1)*DataStructureRigidBody(J)%MechProp(I,5)* &
                             DataStructureRigidBody(J)%Density !accumulating and applying arm rule  
			    end if
				
                DataStructureRigidBody(J)%Normals(I,1) = DataStructureRigidBody(J)%Thickness * &
                    (ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),1),2)-&
                    ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),2),2))/ &
                    DataStructureRigidBody(J)%MechProp(I,1) !normal in x direction
                DataStructureRigidBody(J)%Normals(I,2) = DataStructureRigidBody(J)%Thickness * &
                    (ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),2),1)-&
                    ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I),1),1))/ &
                    DataStructureRigidBody(J)%MechProp(I,1) !normal in y direction
             end do           
            
            
        end do
        allocate(ElementsRigidBody(NThinElements))
        allocate(NormalsRigidBody(NThinElements,2))
        NormalsRigidBody = 0.0
        ! Map element ownership and store corresponding outward normals
        do J = 1, NRigidBodies
            do I = 1, size(DataStructureRigidBody(J)%Ids)
                ElementsRigidBody(DataStructureRigidBody(J)%Ids(I))%RigidBody = J 
                NormalsRigidBody(DataStructureRigidBody(J)%Ids(I),:) = DataStructureRigidBody(J)%Normals(I,:) !duplicated information
            end do
        end do 
        end subroutine  CreateDataStructureRigidBody



        subroutine CreateDataStructureRigidBodyConditions()
        !**********************************************************************
        !
        !    Function:  Create the Data Structure of the Rigid Body Conditions
        !
        !
        !    Note : RigidBodyConditions is allocated in this subroutine
        !
        !**********************************************************************
        implicit none 
        !Variables
        integer(INTEGER_TYPE):: J, I, GlobElmID, M, K, L
        allocate(RigidBodyConditions(NRigidBodies))
        
        do J = 1, NRigidBodies

            !Prescribed velocities and forces are common for all types of analysis

            !Prescribed Velocities (NEEDS TO DECOMPOSED INDIVIDUAL VEL ASSIGNMENTS)
            do I = 1, NThinElmPresVel !loop prescribed vel elements
               read(ThinRigidElmPresVel(I,1), *) GlobElmID !element ID
               GlobElmID = GlobElmID - Counters%NEl !offset for FE elements
               if (GlobElmID == DataStructureRigidBody(J)%Ids(1)) then
                   RigidBodyConditions(J)%ApplyPrescribedVelocity = .true.
                   read(ThinRigidElmPresVel(I,2), *) RigidBodyConditions(J)%VelConstraints(1)
                   read(ThinRigidElmPresVel(I,3), *) RigidBodyConditions(J)%PrescribedVelocity(1)
                   read(ThinRigidElmPresVel(I,4), *) RigidBodyConditions(J)%VelConstraints(2)
                   read(ThinRigidElmPresVel(I,5), *) RigidBodyConditions(J)%PrescribedVelocity(2)
                end if
            end do 

            !Prescribed Forces (NEEDS TO DECOMPOSED INDIVIDUAL FORCE ASSIGNMENTS)
            do I = 1, NThinElmPresFor !loop prescribed force elements
                read(ThinRigidElmPresFor(I,1), *) GlobElmID !element ID
                GlobElmID = GlobElmID - Counters%NEl !offset for FE elements
                read(ThinRigidElmPresFor(I,6), *) L
                if (GlobElmID == DataStructureRigidBody(J)%Ids(1)) then
                     RigidBodyConditions(J)%ApplyPrescribedForce = .true.
                     read(ThinRigidElmPresFor(I,2), *) RigidBodyConditions(J)%ForceConstraints(1)
                     read(ThinRigidElmPresFor(I,3), *) RigidBodyConditions(J)%PrescribedForce(1)
                     read(ThinRigidElmPresFor(I,4), *) RigidBodyConditions(J)%ForceConstraints(2)
                     read(ThinRigidElmPresFor(I,5), *) RigidBodyConditions(J)%PrescribedForce(2)
                     if (L == 1) then
                          RigidBodyConditions(J)%IsThereLinearForce = .true.
                          read(ThinRigidElmPresFor(I,7), *) RigidBodyConditions(J)%LinearInitial
                          read(ThinRigidElmPresFor(I,8), *) RigidBodyConditions(J)%LinearFinal
                     end if
                end if
            end do

            if (.not. ISAXISYMMETRIC) then !plane strain 

            !Prescribed Rotations  
            do I = 1, SIZE(ThinRigidElmPresRot, 1)
               read(ThinRigidElmPresRot(I,1), *) GlobElmID !element ID
               GlobElmID = GlobElmID - SIZE(ElementConnectivities,2)
               if (GlobElmID == DataStructureRigidBody(J)%Ids(1)) then
                   RigidBodyConditions(J)%ApplyPrescribedRotation = .true.
                   read(ThinRigidElmPresRot(I,3), *) RigidBodyConditions(J)%PrescribedRotation
                   read(ThinRigidElmPresRot(I,2), *) K
                   if (K == 0) then
                       RigidBodyConditions(J)%ApplyCenterRotation = 1
                       read(ThinRigidElmPresRot(I,4), *) RigidBodyConditions(J)%CenterRotation(1)
                       read(ThinRigidElmPresRot(I,5), *) RigidBodyConditions(J)%CenterRotation(2)
                   end if
                end if
            end do 

            !Prescribed Torques
            do I = 1, SIZE(ThinRigidElmPresTor, 1)
               read(ThinRigidElmPresTor(I,1), *) GlobElmID !element ID
               GlobElmID = GlobElmID - SIZE(ElementConnectivities,2)
               if (GlobElmID == DataStructureRigidBody(J)%Ids(1)) then
                   RigidBodyConditions(J)%ApplyPrescribedTorque = .true.
                   read(ThinRigidElmPresTor(I,2), *) RigidBodyConditions(J)%PrescribedTorque
                end if
            end do 
            end if            
        end do    
        
        end subroutine CreateDataStructureRigidBodyConditions
    
    
        
        subroutine CreateCoordinatesRenameConnectivitiesIDs()
        !**********************************************************************
        !
        !    Function:  To update RigidElmConnectivities and create ThinRigidCoordinates
        !
        !
        !    Note : ThinRigidCoordinates is allocated in this subroutine
        !
        !**********************************************************************
        implicit none 
        !Variables
        integer(INTEGER_TYPE)::Iel, J, I,M, N, P, Q, C, Ncoord=2, UniqueCoordinates(2*NThinELements)
        logical:: alreadyExists
        !Initializing unique values with the first element

        UniqueCoordinates(1) = ThinRigidElmConnectivities(1,1)
        UniqueCoordinates(2) = ThinRigidElmConnectivities(1,2)

        
        do Iel = 2,  NThinELements !looping all rigid elements
            do J = 1, 2
                do I = 1, Ncoord
                    alreadyExists = .false.
                    if (ThinRigidElmConnectivities(Iel,j) == UniqueCoordinates(I)) then !this means the element is repeated
                        alreadyExists = .true.
                        exit
                    end if
                end do
                ! If the element doesn't exist, increment the counter and store it
                if (.not. alreadyExists) then
                    Ncoord = Ncoord + 1
                    UniqueCoordinates(Ncoord) = ThinRigidElmConnectivities(Iel,j)
                endif
            end do
        end do
        
        !Create ThinRigidCoordinates
        allocate(real(REAL_TYPE) :: ThinRigidCoordinates(Ncoord, 2))
        do J = 1, Ncoord
            ThinRigidCoordinates(J, :) = NodalCoordinates(UniqueCoordinates(J), :)               
        end do
        
        !Update ThinRigidElemConnectivities
        do Iel = 1, NThinELements
            do J = 1, 2
                do I = 1, Ncoord
                    if (ThinRigidElmConnectivities(Iel, J) == UniqueCoordinates(I) ) ThinRigidElmConnectivities(Iel, J) = I                    
                end do
            end do
        end do    
        
      
        end subroutine CreateCoordinatesRenameConnectivitiesIDs
        
        
        
        subroutine NumberOfRigidBodies()
        !**********************************************************************
        !
        !    Function:  To discover the number of rigid bodies and how many
        !               elements are for each rigid body
        !
        !    Note : NRigidBodies is assigned in this subroutine
        !           QtyRigidBodies is allocated in this subroutine
        !           UniqueBodiesNameF is allocated in this subroutine
        !
        !**********************************************************************
        implicit none 
        !Variables
        integer(INTEGER_TYPE):: J, I, Nu=1
        character(len=100):: UniqueBodiesName(NThinELements)
        logical:: alreadyExists

        UniqueBodiesName(1) = ThinRigidElmProp(1,2)
       
        !Counts how many different bodies there are
        do J = 2, NThinElmProp
            do I = 1, Nu
                alreadyExists = .false.
                if (ThinRigidElmProp(J,2) == UniqueBodiesName(I)) then !this means the element is repeated
                    alreadyExists = .true.
                    exit
                end if
            end do
            ! If the element doesn't exist, increment the counter and store it
            if (.not. alreadyExists) then
                Nu = Nu + 1
                UniqueBodiesName(Nu) = ThinRigidElmProp(J,2)
            endif
        end do
        NRigidBodies = Nu
        
        allocate(integer(INTEGER_TYPE) :: QtyRigidBodies(NRigidBodies))
        allocate(character(len=100) :: UniqueBodiesNameF(NRigidBodies))
        
        !Generate a list with the UniqueBodiesNameF
        do I = 1, NRigidBodies
            QtyRigidBodies(I) = 0
            UniqueBodiesNameF(I) = UniqueBodiesName(I)
        end do
        
        !Counts how many elements there are for each rigid body
        do J = 1, NThinElements
            do I = 1, NRigidBodies
                if (ThinRigidElmProp(J,2) == UniqueBodiesName(I)) then
                    QtyRigidBodies(I) =  QtyRigidBodies(I) + 1
                end if
            end do
        end do
        
        end subroutine NumberOfRigidBodies
        
        
        
    subroutine DistanceToRigidBody
    !**********************************************************************
    !
    !    Function:  Determines the euclidean distance of each material point
    !               to the Rigids Bodies
    !    
    !
    ! Implemented in the frame of the MPM project.
    !
    !********************************************************************** 

    implicit none 

    ! Local variables
    integer(INTEGER_TYPE) :: I, J, M
    real(REAL_TYPE) :: Xp(2), X1(2), X2(2), D1, D2, K, V, L_cord, Normal(2), Xpv(2), Aling_element(2), Dot1, Dot2, &
                    X1v(2), Sign

    do I = 1, Counters%NParticles ! Loop over particles ThinRigidElmConnectivities
        Xp = GlobPosArray(I,:) !Coordinate of MP
        X1 = ThinRigidCoordinates(1,:) !Coordinate of node 1 of node 1 (temporary)

        DistanceField(I) = Length(Xp - X1, 2)    !Distance from point to node 1 of element 1 (temporary)
        AffinityArray(I,1) = 1 !1 means node, 0 means element
        AffinityArray(I,2) = 1 !Element ID
        AffinityArray(I,3) = ThinRigidElmConnectivities(1,1) !Node ID

        do J = 1, NThinELements !loop over rigid elements
            X1 = ThinRigidCoordinates(ThinRigidElmConnectivities(J,1),:) !Coordinate of node 1 of element J
            X2 = ThinRigidCoordinates(ThinRigidElmConnectivities(J,2),:) !Coordinate of node 2 of element J

            Aling_element= VectorNorm(X2 - X1, 2) !Vector aligned with the element

            Dot1 = DotProduct(Xp - X1, Aling_element, 2)
            Dot2 = DotProduct(Xp - X2, Aling_element, 2)

            if ((Dot1*Dot2) < 0) then !Perpendicular projection falls inside element
                !Calculate distance to line (norm of (Xp -X1) . normal)
                Normal= NormalsRigidBody(J,:)
                D1= abs(DotProduct(Xp - X1, Normal, 2)) !Distance to line

                if (D1 < DistanceField(I)) then !Current minimum distance
                    DistanceField(I) = D1
                    AffinityArray(I,1) = 0
                    AffinityArray(I,2) = J
                    AffinityArray(I,3) = ThinRigidElmConnectivities(J,1) !for completeness
                end if
            else !Perpendicular projection falls outside element
                D1 = Length(Xp - X1, 2) !Distance to node 1
                D2 = Length(Xp - X2, 2) !Distance to node 2
                if (D1 < DistanceField(I)) then !Distance to node 1 is minimum
                    DistanceField(I) = D1
                    AffinityArray(I,1) = 1
                    AffinityArray(I,2) = J
                    AffinityArray(I,3) = ThinRigidElmConnectivities(J,1)
                end if
                if (D2 < DistanceField(I)) then !Distance to node 2 is minimum
                    DistanceField(I) = D2
                    AffinityArray(I,1) = 1
                    AffinityArray(I,2) = J
                    AffinityArray(I,3) = ThinRigidElmConnectivities(J,2)
                end if   
            end if
        end do
        Particles(I)%ParticleRadius= Kfactor *SQRT(Particles(I)%MASSMIXED/(PI * Particles(I)%DENSITY)) !Calculates particle radius
        DistanceField(I) = DistanceField(I) - Particles(I)%ParticleRadius !Calculate proper distance field
        
        !Now we need to correct particle radius and position in case of undesired initial inter-penetration
        if (DistanceField(I) < 0.0) then !Particle is initially penetrating the rigid body
            !calculate cord length
            L_cord= 2* sqrt(DistanceField(I)**2  - (DistanceField(I)*Particles(I)%ParticleRadius)	) !cord length
            
            if (L_cord > 0.05 * Particles(I)%ParticleRadius) then !This means it is severely crossing the body
                ! Move the particle along the normal
                Normal=AffinityArray(I,1)*(VectorNorm(Xp-ThinRigidCoordinates(AffinityArray(I,3),:), 2))+&
                (1-AffinityArray(I,1))*NormalsRigidBody(AffinityArray(I,2), :) !choses correct normal
                
                Sign=1
            if  (AffinityArray(I,1)==0) then !distance to element, need to check normal direction
        
                !Update Velocity, position and accumulate displacement on MP
                X1v = ThinRigidCoordinates(ThinRigidElmConnectivities(AffinityArray(I,2),1),:) !node 1

                !Dot product to check normal direction
                Dot1= DotProduct(Xp - X1v, Normal, 2)

                if (Dot1 < 0) then !flip normal
                    Sign = -1
                else
                    Sign = 1
                end if
            endif			   
            
            !Kinematic correction				
            !update position
            GlobPosArray(I,:) = GlobPosArray(I,:) + Sign * (abs(DistanceField(I))) * Normal!  
            DistanceField(I)=0.0  
            else !Not so bad, change particle radius
            Particles(I)%ParticleRadius= Particles(I)%ParticleRadius + DistanceField(I)
            DistanceField(I)=0.0
            endif				  
            
        end if

    end do

    end subroutine DistanceToRigidBody
         
        subroutine InitialiseStructuralElements()
        !**********************************************************************
        !
        !    Function:  Determines the euclidean distance of each material point
        !               to the Rigid Bodies
        !    
        !
        ! Implemented in the frame of the MPM project.
        !
		!**********************************************************************

        implicit none 
		if (.not. IsThereRigidBody) RETURN
		
        call CreateCoordinatesRenameConnectivitiesIDs() !creates connectivities and coords list
        call NumberOfRigidBodies() !identifies the number of rigid bodies and ind. elements
        call CreateDataStructureRigidBody() !Creates a structure containing info about the rigid bodies and their properties
        call CreateDataStructureRigidBodyConditions() !Creates a structure containing info about the conditions applied to the rigid bodies
        call DistanceToRigidBody() !Calculates distance to rigid bodies and fixes initial interpenetrations if needed

        end subroutine InitialiseStructuralElements
          
          
		subroutine RigidBodyLagrangianPhase()
        !**********************************************************************
        !
        !    Function: Calculate the new translation and rotation for each rigid body 
        !               
        !    
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none 
        
        ! Local variables
        integer(INTEGER_TYPE) :: J, I, K, NodeID
        real(REAL_TYPE) :: GravityAcceleration(2), Weight(2), THETA, CosTheta, SinTheta
        real(REAL_TYPE) :: PrescribedForce(2), Mass, NodePos(2), RotCenter(2)
        real(REAL_TYPE) :: X1(2), X2(2), DeltaCoord(2)
        logical :: IsProjectErosion
        real(REAL_TYPE), save :: THETATOT = 0.0, SIGN_VAR = -1.0
        real(REAL_TYPE) :: ThetaMax, PI, XC, YC
       
        if (.not. IsThereRigidBody) RETURN
        if (CalParams%TimeStep == 1 ) RETURN !does not calculate for first time step
        
        ! Calculate gravity acceleration (vectorial)
        GravityAcceleration = CalParams%Multipliers%GravityRealised * CalParams%GravityData%GravityVector * &
                              CalParams%GravityData%GAccel
       
        ! Loop over all rigid bodies
        do J = 1, NRigidBodies
            
            ! Calculate body mass and weight (vectorial)
            Mass = DataStructureRigidBody(J)%Density * DataStructureRigidBody(J)%Area
            if (.not. ISAXISYMMETRIC) then
                Weight = Mass * GravityAcceleration / 1000.0 !in KN <== problem if generalized for other units
            else
                Weight = Mass * GravityAcceleration !Problem (Needs to consider rotation: FIX)
            end if
            
            ! Apply linear force ramping if specified (vectorial)
            if (RigidBodyConditions(J)%IsThereLinearForce) then
                PrescribedForce = RigidBodyConditions(J)%PrescribedForce * &
                    (CalParams%TimeStep * CalParams%TimeIncrement) / &
                    (CalParams%TotalTime * CalParams%NLOADSTEPS)
            else 
                PrescribedForce = RigidBodyConditions(J)%PrescribedForce
            end if
            
            ! Apply force constraints (decompose prescribed forces)
            PrescribedForce = PrescribedForce * RigidBodyConditions(J)%ForceConstraints
            
            ! Calculate acceleration (vectorial)
            DataStructureRigidBody(J)%Acceleration = (DataStructureRigidBody(J)%Fsum + Weight + &
                                                       PrescribedForce) / Mass
			! Calculate velocity
			DataStructureRigidBody(J)%Velocity = DataStructureRigidBody(J)%Velocity + &
                                                 CalParams%TimeIncrement * DataStructureRigidBody(J)%Acceleration
            
            ! Update velocity (vectorial with constraints decomposition)
            if (RigidBodyConditions(J)%ApplyPrescribedVelocity) then
                ! Apply velocity constraints (decompose prescribed velocities)
                DataStructureRigidBody(J)%Velocity = RigidBodyConditions(J)%PrescribedVelocity * &
													 RigidBodyConditions(J)%VelConstraints + &
													 DataStructureRigidBody(J)%Velocity * &
													 (1.0 - RigidBodyConditions(J)%VelConstraints)
            end if
            
            ! Update nodal coordinates (vectorial)
            do I = 1, SIZE(DataStructureRigidBody(J)%Nodes, 1)
                NodeID = DataStructureRigidBody(J)%Nodes(I)
                ThinRigidCoordinates(NodeID, :) = ThinRigidCoordinates(NodeID, :) + &
                                                   CalParams%TimeIncrement * DataStructureRigidBody(J)%Velocity
            end do
            
            ! Update centroid position (vectorial)
            DataStructureRigidBody(J)%Centroid = DataStructureRigidBody(J)%Centroid + &
                                                  CalParams%TimeIncrement * DataStructureRigidBody(J)%Velocity

            ! Rotation update (only for plane strain)
            if (.not. ISAXISYMMETRIC) then
                
                ! Calculate angular acceleration
                DataStructureRigidBody(J)%AngularAcceleration = (DataStructureRigidBody(J)%Tsum + &
                    RigidBodyConditions(J)%PrescribedTorque) / (DataStructureRigidBody(J)%Inertia + 0.1)               

				
                if (RigidBodyConditions(J)%ApplyPrescribedRotation) then !fixed rotation rate
                    DataStructureRigidBody(J)%AngularVelocity = RigidBodyConditions(J)%PrescribedRotation
                else !free
                    DataStructureRigidBody(J)%AngularVelocity = DataStructureRigidBody(J)%AngularVelocity + &
                        CalParams%TimeIncrement * DataStructureRigidBody(J)%AngularAcceleration
                end if
                
                ! Calculate rotation angle
                THETA = CalParams%TimeIncrement * DataStructureRigidBody(J)%AngularVelocity
                
                ! PROJECT EROSION (special case - typically disabled)
                IsProjectErosion = .false.
                if (IsProjectErosion) then
                    DataStructureRigidBody(J)%AngularVelocity = 0.05
                    THETA = SIGN_VAR * CalParams%TimeIncrement * DataStructureRigidBody(J)%AngularVelocity
                    THETATOT = THETATOT + THETA
                    XC = 0.0
                    YC = 0.0
                    PI = 22.0 / 7.0
                    ThetaMax = 20.0 * PI / 180.0
                    if (abs(abs(THETATOT) - ThetaMax) < 0.01) then
                        SIGN_VAR = -SIGN_VAR
                        THETATOT = 0.0
                    end if
                    RotCenter = (/XC, YC/)
                else
                    RotCenter = DataStructureRigidBody(J)%Centroid
                end if
                
                ! Apply rotation to all nodes (vectorial)
                CosTheta = COS(THETA)
                SinTheta = SIN(THETA)
                do K = 1, SIZE(DataStructureRigidBody(J)%Nodes, 1)
                    NodeID = DataStructureRigidBody(J)%Nodes(K)
                    NodePos = ThinRigidCoordinates(NodeID, :)
                    DeltaCoord = NodePos - RotCenter
                    ThinRigidCoordinates(NodeID, 1) = RotCenter(1) + DeltaCoord(1) * CosTheta - &
                                                       DeltaCoord(2) * SinTheta
                    ThinRigidCoordinates(NodeID, 2) = RotCenter(2) + DeltaCoord(1) * SinTheta + &
                                                       DeltaCoord(2) * CosTheta
                end do
                
                ! Update normals (vectorial)
                do I = 1, QtyRigidBodies(J)
                    X1 = ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I), 1), :)
                    X2 = ThinRigidCoordinates(ThinRigidElmConnectivities(DataStructureRigidBody(J)%Ids(I), 2), :)
                    DataStructureRigidBody(J)%Normals(I, 1) = DataStructureRigidBody(J)%Thickness * &
                        (X1(2) - X2(2)) / DataStructureRigidBody(J)%MechProp(I, 1)
                    DataStructureRigidBody(J)%Normals(I, 2) = DataStructureRigidBody(J)%Thickness * &
                        (X2(1) - X1(1)) / DataStructureRigidBody(J)%MechProp(I, 1)
                end do
                
                ! Update global normals array
                do I = 1, SIZE(DataStructureRigidBody(J)%Ids)
                    ElementsRigidBody(DataStructureRigidBody(J)%Ids(I))%RigidBody = J 
                    NormalsRigidBody(DataStructureRigidBody(J)%Ids(I), :) = DataStructureRigidBody(J)%Normals(I, :)
                end do
            end if
            
            ! Reset force and torque accumulators (vectorial)
            DataStructureRigidBody(J)%Fsum = 0.0
            DataStructureRigidBody(J)%Tsum = 0.0
            
        end do

	end subroutine RigidBodyLagrangianPhase
		  
	subroutine Get_Contact_Stiffness()
	!**********************************************************************
    !
    !    Function: calculates the maximum contact stiffness for penalty force 
    !               
    !    
    !    Implemented in the frame of the MPM project.
    !
    !**********************************************************************
	implicit none
	integer(INTEGER_TYPE) :: Ipart
	double precision:: kappa_crit
	
	if (.not. IsThereRigidBody) RETURN
	if (Stiffness_assigned) RETURN
	do Ipart=1, Counters%NParticles !loop particles
		kappa_crit= 4.0* Particles(Ipart)%MASSMIXED * PI**2.0 /CalParams%TimeIncrement**2.0  !critical stiffness 			  
		
		if (kappa<kappa_crit * 0.003) kappa=kappa_crit * 0.003
		Stiffness_assigned=.true.
	end do
	
	
	end subroutine Get_Contact_Stiffness
          
          
    subroutine LevelSetContact()
    !**********************************************************************
    !
    !    Function: Calculate the new traslation and rotation for each rigid body 
    !               
    !    
    !
    ! Implemented in the frame of the MPM project.
    !
    !**********************************************************************

    implicit none 

    ! Local variables
    integer(INTEGER_TYPE) :: I, J, MaterialID
    real(REAL_TYPE) :: Xp(2), X1(2), X2(2), D1, D2, Normal(2), Aling_element(2)
    real(REAL_TYPE) :: Dot1, Dot2, Sign, normal_diff(2), RadiusVec(2)
    real(REAL_TYPE) :: RigidBodyVelAtPoint(2), Xcross(2), ProjdX_n(2), Proj1(2)
    real(REAL_TYPE) :: TimeIncrementNew, RelativeVel(2), Tangent(2), mu_actual
    real(REAL_TYPE) :: k_of_d, NormalForce(2), TangentForce(2), Force(2)
    real(REAL_TYPE) :: Arm, Torque, DampingRatio=0.75, X1v(2)

    
    if (.not. IsThereRigidBody) RETURN
    if (CalParams%TimeStep == 1 ) RETURN !do not calculate for first time step
    
    ContactForceArray(:,:)=0.0 !initialize contact force array

    do I = 1, Counters%NParticles ! Loop over particles ThinRigidElmConnectivities
        Xp = GlobPosArray(I,:) !Coordinate of MP
        X1 = ThinRigidCoordinates(1,:) !Coordinate of node 1 of node 1 (temporary)

        DistanceField(I) = Length(Xp - X1, 2)    !Distance from point to node 1 of element 1 (temporary)
        AffinityArray(I,1) = 1 !1 means node, 0 means element
        AffinityArray(I,2) = 1 !Element ID
        AffinityArray(I,3) = ThinRigidElmConnectivities(1,1) !Node ID

        do J = 1, NThinELements !loop over rigid elements
            X1 = ThinRigidCoordinates(ThinRigidElmConnectivities(J,1),:) !Coordinate of node 1 of element J
            X2 = ThinRigidCoordinates(ThinRigidElmConnectivities(J,2),:) !Coordinate of node 2 of element J

            Aling_element= VectorNorm(X2 - X1, 2) !Vector aligned with the element

            Dot1 = DotProduct(Xp - X1, Aling_element, 2)
            Dot2 = DotProduct(Xp - X2, Aling_element, 2)

            if ((Dot1*Dot2) < 0) then !Perpendicular projection falls inside element
                !Calculate distance to line (norm of (Xp -X1) . normal)
                Normal= NormalsRigidBody(J,:)
                D1= abs(DotProduct(Xp - X1, Normal, 2)) !Distance to line

                if (D1 < DistanceField(I)) then !Current minimum distance
                    DistanceField(I) = D1
                    AffinityArray(I,1) = 0
                    AffinityArray(I,2) = J
                    AffinityArray(I,3) = ThinRigidElmConnectivities(J,1) !for completeness
                end if
            else !Perpendicular projection falls outside element
                D1 = Length(Xp - X1, 2) !Distance to node 1
                D2 = Length(Xp - X2, 2) !Distance to node 2
                if (D1 < DistanceField(I)) then !Distance to node 1 is minimum
                    DistanceField(I) = D1
                    AffinityArray(I,1) = 1
                    AffinityArray(I,2) = J
                    AffinityArray(I,3) = ThinRigidElmConnectivities(J,1)
                end if
                if (D2 < DistanceField(I)) then !Distance to node 2 is minimum
                    DistanceField(I) = D2
                    AffinityArray(I,1) = 1
                    AffinityArray(I,2) = J
                    AffinityArray(I,3) = ThinRigidElmConnectivities(J,2)
                end if   
            end if
        end do

    DistanceField(I) = DistanceField(I) - Particles(I)%ParticleRadius !account for radius of influence
        
    if (DistanceField(I)<0.0) then !contact exist between particle and rigid body
        J = ElementsRigidBody(AffinityArray(I,2))%RigidBody !obtain rigid body ID

        Normal=AffinityArray(I,1)*(VectorNorm(Xp-ThinRigidCoordinates(AffinityArray(I,3),:), 2))+&
               (1-AffinityArray(I,1))*NormalsRigidBody(J, :) !choses correct normal

        Sign=1
        if  (AffinityArray(I,1)==0) then !distance to element, need to check normal direction
    
            !Update Velocity, position and accumulate displacement on MP
            X1v = ThinRigidCoordinates(ThinRigidElmConnectivities(AffinityArray(I,2),1),:) !node 1

            !Dot product to check normal direction
            Dot1= DotProduct(Xp - X1v, Normal, 2)

            if (Dot1 < 0) then !flip normal
                Sign = -1
            else
                Sign = 1
            end if
        endif	
        
        ! Calculate rigid body velocity at contact point including rotation
        ! v_point = v_translation + omega x r
        ! In 2D: v_point = v_translation + omega * (-ry, rx)
        if (.not. ISAXISYMMETRIC) then
            normal_diff= (Particles(I)%ParticleRadius + DistanceField(I))  * Sign * Normal !Remainder vector from particle to point of contact
            RadiusVec = (Xp - DataStructureRigidBody(J)%Centroid) - normal_diff !Vector from centroid to contact point
            RigidBodyVelAtPoint(1) = DataStructureRigidBody(J)%Velocity(1) + &
                                        DataStructureRigidBody(J)%AngularVelocity * (RadiusVec(2))
            RigidBodyVelAtPoint(2) = DataStructureRigidBody(J)%Velocity(2) + &
                                        DataStructureRigidBody(J)%AngularVelocity * (- RadiusVec(1))
        else
            ! For axisymmetric, only translational velocity (no rotation in this plane)
            RigidBodyVelAtPoint = DataStructureRigidBody(J)%Velocity
        end if
        
        !Check future crossing for automatic time-stepping control
        Xcross=GlobPosArray(I,:)+(VelocityArray(I,:)-RigidBodyVelAtPoint)*CalParams%TimeIncrement !anticipated position of MP

        ProjdX_n=DotProduct(Xp-ThinRigidCoordinates(AffinityArray(I,3),:), NormalsRigidBody(AffinityArray(I,2), :), 2)&
            *NormalsRigidBody(AffinityArray(I,2), :) !projection of x_rel onto n
        Proj1=DotProduct(Xcross-ThinRigidCoordinates(AffinityArray(I,3),:), NormalsRigidBody(AffinityArray(I,2), :), 2)&
            *NormalsRigidBody(AffinityArray(I,2), :) !projection of x_rel in future onto n
        Dot2=DotProduct(ProjdX_n, Proj1, 2)
        if (Dot2<0.0) then !critical time must be reduced (useful in hypervelocity cases NOT TESTED YET)
            TimeIncrementNew=CalParams%courantNumber* Particles(I)%ParticleRadius/&
                Length(VelocityArray(I,:)-RigidBodyVelAtPoint,2)
            if (CalParams%TimeIncrement > TimeIncrementNew) CalParams%TimeIncrement = TimeIncrementNew 
        endif   
        

        !calculate tangent vector using relative velocity
        RelativeVel = VelocityArray(I,:) - RigidBodyVelAtPoint
        ProjdX_n=DotProduct(RelativeVel, Sign*Normal, 2)*Sign*Normal !projection of v_rel onto n
        Tangent=RelativeVel-ProjdX_n !Use this to compare with friction angle
        mu_actual=Length(Tangent, 2)/Length(ProjdX_n, 2)
        Tangent=Sign*VectorNorm(Tangent, 2) !normal tangent vector

        MaterialID = MaterialIDArray(I)
                    
        k_of_d=-(2.0*log((Distancefield(I)+Particles(I)%ParticleRadius)/(Particles(I)%ParticleRadius))- &
            (Particles(I)%ParticleRadius/(Distancefield(I)+Particles(I)%ParticleRadius))+1.0) !stiffness as function of distance field

        NormalForce=-((DataStructureRigidBody(J)%Materials(MaterialID)%ContactStiffness)*k_of_d*(-1 * Distancefield(I)) - DampingRatio * &
        2.0*(Particles(I)%MASSMIXED*DataStructureRigidBody(J)%Materials(MaterialID)%ContactStiffness)**0.5*(DotProduct(RelativeVel, Sign*Normal, 2)))*Sign*Normal 
        
        if ( ISAXISYMMETRIC) NormalForce= Xp(1) * NormalForce !account for axisymmetric volume             
        
        if (AffinityArray(I,1)==0) then !use friction on wall
            
            if (DotProduct(RelativeVel, Sign*Normal, 2)<0) then !Approaching contact
            if (mu_actual > DataStructureRigidBody(J)%Materials(MaterialID)%FrictionCoef) mu_actual=DataStructureRigidBody(J)%Materials(MaterialID)%FrictionCoef !If exceeds friction limit otherwise is sticky contact

            TangentForce=Length(NormalForce, 2)*mu_actual*Tangent !Tangenttial force on rigid body

            else !Receding contact
                TangentForce=0.0
            end if

            Force=NormalForce+TangentForce  !total force on rigid body
        else !point in contact with node, no friction
            
            Force=NormalForce 
        
        endif

        !update force
        DataStructureRigidBody(J)%Fsum = DataStructureRigidBody(J)%Fsum + Force
        !DataStructureRigidBody(J)%Fsum(1)=0.0 !(Brute force to avoid X direction movement)
        
        !update torque
        if (.not. ISAXISYMMETRIC) then !only for 2D plain strain rotation
            Arm=Length(RadiusVec, 2) !lever arm
            Torque=(RadiusVec(1)*Force(2))-(Force(1)*RadiusVec(2)) !Comes from cross product and follows right hand rule
            DataStructureRigidBody(J)%Tsum = DataStructureRigidBody(J)%Tsum+Torque
        else
            Torque=0.0 
        end if            

        !Add force to material point
        ContactForceArray(I, :)=-Force !Action-reaction
        
    else 
        TangentForceTrack(I,:)=0.0
    end if
            
    end do
    end subroutine LevelSetContact
		   
		   
		   
		subroutine UpdateThinElementImpulse()
        !**********************************************************************
        !
        !    Function: Updates the forces acting on the rigid body 
        !               
        !    
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

     !   use ModCounters

        implicit none 

          ! Local variables
          integer(INTEGER_TYPE) :: I, J, M, NodeID, ElementID, CoUp, CoDown 
          real(REAL_TYPE) :: ProjdX_n(2), Xp(2), X1(2), dot1, dot2, Yp, Y1, X2(2), d, Y2, D1, D2, K, V, DistanceRatio, Normal(2), &
              Xn, Yn, X0, Y0, Nx1, Ny1, Nx2, Ny2, NormalUpID(2,2), NormalDownID(2,2), TempCos, &
              CS, SN, Xf, Yf, Gx, Gy, NormalForce(2), TangentForce(2), Tangent(2), DistanceForce, Force(2),&
              Xg, Yg, L(2), Lu(2), Arm,  MaterialID, Sign, Torque
          logical :: IsDistanceToNode, IsElemConnToNode, IsThereUp, IsThereDown
          
          if (.not. IsThereRigidBody) RETURN
		  if (CalParams%TimeStep == 1 ) RETURN !does not calculate for first time step
		  
		      do I = 1, Counters%NParticles ! Loop over particles ThinRigidElmConnectivities
			  Xp=GlobPosArray(I,:) !particle position			  
			  
              if (DistanceField(I)<0.0) then !contact exist between particle and rigid body
				  Normal=AffinityArray(I,1)*(VectorNorm(Xp-ThinRigidCoordinates(AffinityArray(I,3),:), 2))+&
					  (1-AffinityArray(I,1))*NormalsRigidBody(AffinityArray(I,2), :) !choses correct normal

				Sign=1
                if  (AffinityArray(I,1)==0) then 
				
                  !Update Velocity, position and accumulate displacement on MP
                  X1 = ThinRigidCoordinates(ThinRigidElmConnectivities(AffinityArray(I,2),1),:)
                  !Y1 = ThinRigidCoordinates(ThinRigidElmConnectivities(ElementID,1),2)
                  X2 = ThinRigidCoordinates(ThinRigidElmConnectivities(AffinityArray(I,2),2),:)
                  !Y2 = ThinRigidCoordinates(ThinRigidElmConnectivities(ElementID,2),2)
				  
                   if ((((X1(1)-X2(1))*Normal(2)-(X1(2)-X2(2))*Normal(1)))*&
					   ((X1(1)-X2(1))*(Xp(2)-X1(2))-(X1(2)-X2(2))*(Xp(1)-X1(1)))<0) then !flop normal
                       Sign = -1
                   else
                       Sign = 1
				   end if
				endif
				   
				   
				   !Kinematic correction
				
				   !update position
                   !GlobPosArray(I,:) = GlobPosArray(I,:) + Sign * (abs(DistanceField(I))) * Normal!                   
				   
				   !update velocity array
                   !VelocityArray(I,:) =  VelocityArray(I,:) + Sign * (abs((DistanceField(I)))/(CalParams%TimeIncrement)) * Normal!                   
				   
				   !Update commutative displacement array
                   !UArray(I,:) = UArray(I,:) + Sign * (abs(DistanceField(I))) * Normal!      
				   
				   !Update incremental displacement array
				   !UStepArray(I, :)= Sign * (abs(DistanceField(I))) * Normal!  
				   
				   !Kinetic correction
				   
				   !Calculate Force
                   Force(1) = ((SigmaEffArray(I,1)+Particles(I)%WaterPressure)*Sign*Normal(1)+SigmaEffArray(I,4)*Sign*Normal(2))*2*Kfactor * SQRT(Particles(I)%MASSMIXED/(PI * Particles(I)%DENSITY)) !Important to update for water pressure
                   Force(2) = (SigmaEffArray(I,4)*Sign*Normal(1) + (SigmaEffArray(I,2)+Particles(I)%WaterPressure)*Sign*Normal(2))*2*Kfactor * SQRT(Particles(I)%MASSMIXED/(PI * Particles(I)%DENSITY))
				   !calculate tangent vector
				   ProjdX_n=DotProduct(VelocityArray(I,:), Sign*Normal, 2)*Sign*Normal !projection of v onto n
				   Tangent=VelocityArray(I,:)-ProjdX_n
				   Tangent=VectorNorm(Tangent, 2) !normal tangent vector
				   
				   !Calculate normal force
				   NormalForce=DotProduct(Force, Sign*Normal, 2)*Sign*Normal !projection of F onto n
				   NormalForce=-Length(NormalForce, 2)*Sign*Normal !Ensuring correct direction
				   
				   !calculate tangential force
				   J = ElementsRigidBody(AffinityArray(I,2))%RigidBody
				   if (AffinityArray(I,1)==0) then !use friction on wall
				   
                   MaterialID = MaterialIDArray(I)
				   TangentForce=Length(NormalForce, 2)*DataStructureRigidBody(J)%Materials(MaterialID)%FrictionCoef*Tangent 
				   Force=NormalForce+TangentForce  
				   
				   endif
				   DataStructureRigidBody(J)%Fsum = DataStructureRigidBody(J)%Fsum + Force!update force
				   
				   
				   !update torque
				   L=Xp-DataStructureRigidBody(J)%Centroid
				   Arm=Length(L, 2)-Kfactor * SQRT(Particles(I)%MASSMIXED/(PI * Particles(I)%DENSITY))
				   L=Arm*VectorNorm(L, 2) 
				   Torque=(L(1)*Force(2))-(Force(1)*L(2))
				   DataStructureRigidBody(J)%Tsum = DataStructureRigidBody(J)%Tsum+Torque
				   
				   
				   !Add for to material point
				   ContactForceArray(I, :)=-Force+ContactForceArray(I, :)

             end if
          end do


		end subroutine UpdateThinElementImpulse
		
		
        subroutine DisplacementCorrection()
        !**********************************************************************
        !
        !    Function: Calculate the new traslation and rotation for each rigid body 
        !               
        !    
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

     !   use ModCounters

        implicit none 

          ! Local variables
          integer(INTEGER_TYPE) :: I, J, M, NodeID, ElementID, CoUp, CoDown 
          real(REAL_TYPE) :: ProjdX_n(2), Xp(2), X1(2), dot1, dot2, Yp, Y1, X2(2), d, Y2, D1, D2, K, V, DistanceRatio, Normal(2), &
              Xn, Yn, X0, Y0, Nx1, Ny1, Nx2, Ny2, NormalUpID(2,2), NormalDownID(2,2), TempCos, &
              CS, SN, Xf, Yf, Gx, Gy, NormalForce(2), TangentForce(2), Tangent(2), DistanceForce, Force(2),&
              Xg, Yg, L(2), Lu(2), Arm(2),  MaterialID, Sign, Torque
          logical :: IsDistanceToNode, IsElemConnToNode, IsThereUp, IsThereDown
          
          if (.not. IsThereRigidBody) RETURN
		  if (CalParams%TimeStep == 1 ) RETURN !does not calculate for first time step

          do I = 1, Counters%NParticles ! Loop over particles ThinRigidElmConnectivities
			  Xp=GlobPosArray(I,:) !particle position			  
			  X1=ThinRigidCoordinates(1,:) !local node 1 of rigid element

              DistanceField(I) = Length(X1-Xp, 2) !norm of diff vector
              IsDistanceToNode = .true.
			  AffinityArray(I,1) = 1
              AffinityArray(I,2) = 1
			  AffinityArray(I,3) = ThinRigidElmConnectivities(1,1)

              do J = 1, SIZE(ThinRigidElmConnectivities,1) !loop over rigid elements
				  X1=ThinRigidCoordinates(ThinRigidElmConnectivities(J,1),:) !local node 1 of rigid element
				  X2=ThinRigidCoordinates(ThinRigidElmConnectivities(J,2),:) !local node 2 of rigid element
				  dot1=DotProduct(X2-X1, Xp-X1, 2)
				  dot2=DotProduct(X1-X2, Xp-X2, 2)
                  if (dot1>0 .and. dot2>0) then !point projects to element
					  ProjdX_n=DotProduct(Xp-X1, NormalsRigidBody(J, :), 2)*NormalsRigidBody(J, :) !projection of dx onto n
					  D1=Length(ProjdX_n, 2)!distance to element			
					  
                      if (D1 < DistanceField(I)) then !update distance field
                          DistanceField(I) = D1
						  AffinityArray(I,1) = 0
						  AffinityArray(I,2) = J
						  AffinityArray(I,3) = ThinRigidElmConnectivities(J,1) !for completeness
                       end if
                  else !point projects to node
                      D1 = Length(Xp-X1, 2)!distance to node 1
                      D2 = Length(Xp-X2, 2)!distance to node 2
                      if (D1 < DistanceField(I)) then !point projects to node 1

                         DistanceField(I) = D1
						 AffinityArray(I,1) = 1
						 AffinityArray(I,2) = J
						 AffinityArray(I,3) = ThinRigidElmConnectivities(J,1)

                      end if
                      if (D2 < DistanceField(I)) then

                         DistanceField(I) = D2

						AffinityArray(I,1) = 1
						AffinityArray(I,2) = J
						 AffinityArray(I,3) = ThinRigidElmConnectivities(J,2)

                      end if   
                  end if
              end do
              DistanceField(I) = DistanceField(I) - &
                  Kfactor * SQRT(Particles(I)%MASSMIXED/(PI * Particles(I)%DENSITY)) !account for radius of influence
              if (DistanceField(I)<0.0) then !contact exist between particle and rigid body
				  Normal=AffinityArray(I,1)*(VectorNorm(Xp-ThinRigidCoordinates(AffinityArray(I,3),:), 2))+&
					  (1-AffinityArray(I,1))*NormalsRigidBody(AffinityArray(I,2), :) !choses correct normal

				Sign=1
                if  (AffinityArray(I,1)==0) then 
				
                  !Update Velocity, position and accumulate displacement on MP
                  X1 = ThinRigidCoordinates(ThinRigidElmConnectivities(AffinityArray(I,2),1),:)
                  !Y1 = ThinRigidCoordinates(ThinRigidElmConnectivities(ElementID,1),2)
                  X2 = ThinRigidCoordinates(ThinRigidElmConnectivities(AffinityArray(I,2),2),:)
                  !Y2 = ThinRigidCoordinates(ThinRigidElmConnectivities(ElementID,2),2)
				  
                   if ((((X1(1)-X2(1))*Normal(2)-(X1(2)-X2(2))*Normal(1)))*&
					   ((X1(1)-X2(1))*(Xp(2)-X1(2))-(X1(2)-X2(2))*(Xp(1)-X1(1)))<0) then !flip normal
                       Sign = -1
                   else
                       Sign = 1
				   end if
				endif			   
				   
				   !Kinematic correction
				
				   !update position
                   GlobPosArray(I,:) = GlobPosArray(I,:) + Sign * (abs(DistanceField(I))) * Normal!          
				   

             end if
          end do


		end subroutine DisplacementCorrection 
		
		
		subroutine GetNodalContactForce(Fcontact)
        !**********************************************************************
        !
        !    Function: tracking information of the rigid body for each time step
        !               
        !    
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

     !   use ModCounters
		implicit none
		
		real(REAL_TYPE), intent(out):: Fcontact(Counters%N , 2)
		real(REAL_TYPE):: PartShape(ELEMENTNODES)
		integer(INTEGER_TYPE):: Iel, INode, IPart, ActEl, nn, IDof, NElemPart, IntGlo
		Fcontact=0.0
		!Get number of active elements
		do Iel=1,Counters%NAEl !loop active element
			ActEL=ActiveElement(Iel) !Active element
			NElemPart = NPartEle(ActEL)! number of particles inside element
			
			do Ipart=1,NElemPart !loop particles
			IntGlo = GetParticleIndex(Ipart, ActEL)	!gets global particle index	
			PartShape = ShapeValuesArray(IntGlo,:) !get particle shape fucntions
			!loop nodes of element
			do INode=1,ELEMENTNODES
				nn=ElementConnectivities(iNode,ActEL) !local node connectivities
                IDof = ReducedDof(nn) !global DOF number
				!loop material points
				
				Fcontact(Idof+1, 1)=Fcontact(Idof+1, 1)+ContactForceArray(IntGlo, 1)*PartShape(INode)
				Fcontact(Idof+2, 1)=Fcontact(Idof+2, 1)+ContactForceArray(IntGlo, 2)*PartShape(INode)
			enddo
			enddo
			
	    enddo
		
		end subroutine GetNodalContactForce
		
           
        subroutine WriteRigidBodyData()
        !**********************************************************************
        !
        !    Function: tracking information of the rigid body for each time step
        !               
        !    
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

     !   use ModCounters

        implicit none 

          ! Local variables
          integer(INTEGER_TYPE) :: I, J, M, NodeID, ElementID
          character(len=20) :: filename
          
          write(filename, "('Data/result_', I0, '.txt')") CalParams%TimeStep
          
          open(unit=10, file=filename, status='unknown')
           write(10, *)'Body_Name ', ' ID ', ' Normal X ', ' Normal Y '
          !write normals
          do J = 1, NRigidBodies
              do I = 1, QtyRigidBodies(J)
                  write(10, *) DataStructureRigidBody(J)%Name, DataStructureRigidBody(J)%Ids(I) +  &
                      SIZE(ElementConnectivities,2), DataStructureRigidBody(J)%Normals(I,1), &
                      DataStructureRigidBody(J)%Normals(I,2)
              end do
          end do
          close(10)
          
        end subroutine WriteRigidBodyData
         
    
        end module ModStructuralElements