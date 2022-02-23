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
    !	Copyright (C) 2021  Members of the Anura3D MPM Research Community 
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
	  
	  
	  module ModMPMDynContact
      !**********************************************************************
      !
      !    Function:  This module contains all routines which are used for the Dynamic MPM contact model
      !
      !     $Revision: 9077 $
      !     $Date: 2021-03-05 15:03:34 +0100 (vr, 05 mrt. 2021) $
      !
      !**********************************************************************
      use ModCounters
      use ModReadCalculationData
      use ModMPMData
      use ModParticle
      use ModMeshInfo
      use ModRotBoundCond
      use ModGlobalConstants
      use ModFileIO
      
      implicit none
        
        character(len=MAX_FILENAME_LENGTH) :: FileVersion ! FileVersion for contact data
        character(len=64), dimension(:,:), allocatable :: ContactMaterialName ! array including the Contact Material name
        real(REAL_TYPE), dimension(:), allocatable :: LumpedMassSys ! Lumped mass vector of the entire system
        real(REAL_TYPE), dimension(:), allocatable :: TotalVelocitySys ! Total velocity vector of the entire system
        real(REAL_TYPE), dimension(:), allocatable :: AccelerationSys ! Acceleration vector of the entire system
        real(REAL_TYPE), dimension(:, :), allocatable :: ContactNodeNormals, ContactSurfaces ! Normals at nodes used with the contact algorithm
        logical, dimension(:), allocatable :: InterfaceNodes
        real(REAL_TYPE), dimension(:,:), allocatable :: FrictionVector ! frictioncoef for multiple materials
        integer(INTEGER_TYPE) :: NInterfaceNodes, StructureElement
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: NodeElement ! Node-element relation - for each node a list of the elements containing this global node
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ContactSurfaceElmNodes
        integer(INTEGER_TYPE):: NContactElm, PileEntityMaterial ! entity of the contact volume
        real(REAL_TYPE), dimension(:,:), allocatable :: InterfaceNodesAdhesion ! Adhesion for multiple materials 
        logical, dimension(:), allocatable :: ContactSurfaceSoilElements
        logical :: GlobContElement (:)  ! global copy of element(:)
        allocatable :: GlobContElement
        integer(INTEGER_TYPE) :: NContVolElem         ! global copy MPMInit

      contains ! Routines of this module


        subroutine InitialiseContactData()
        !**********************************************************************
        !
        !    Function:  Contains code for initialising contact data
        !
        !**********************************************************************
        implicit none
        
          call DestroyContactData()
          
          call InitialiseContactArrays()
      
        end subroutine InitialiseContactData


        subroutine DestroyContactData()
        !**********************************************************************
        !
        !    Function:  Deallocates the arrays used in this module
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(LumpedMassSys)) then
            deallocate(LumpedMassSys, stat = IError)
          end if

          if (allocated(TotalVelocitySys)) then
            deallocate(TotalVelocitySys, stat = IError)
          end if

          if (allocated(AccelerationSys)) then
            deallocate(AccelerationSys, stat = IError)
          end if

          if (allocated(ContactNodeNormals) ) then
            deallocate(ContactNodeNormals, stat = IError)
          end if
          
           if (allocated(InterfaceNodes) ) then
            deallocate(InterfaceNodes, stat = IError)
          end if
          
          if (allocated(InterfaceNodesAdhesion) ) then
            deallocate(InterfaceNodesAdhesion, stat = IError)
          end if

          if (allocated(NodeElement)) then
            deallocate(NodeElement, stat = IError)
          end if

          if (allocated(ContactSurfaceSoilElements)) then
            deallocate(ContactSurfaceSoilElements, stat = IError)
          end if

        end subroutine DestroyContactData


        subroutine InitialiseContactArrays()
        !**********************************************************************
        !
        !    Function:  To initialise the arrays relate to contact
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: IError

          if (CalParams%ApplyContactAlgorithm) then
            allocate(LumpedMassSys(Counters%N), stat = IError)
            allocate(TotalVelocitySys(Counters%N), stat = IError)
            allocate(AccelerationSys(Counters%N), stat = IError)
            allocate(ContactNodeNormals(NVECTOR, Counters%NodTot), stat = IError)
            allocate(InterfaceNodes(Counters%NodTot), stat = IError)
            allocate(NodeElement(Counters%NodTot, 100), stat = IError)
            allocate(InterfaceNodesAdhesion(Counters%NodTot,CalParams%NumberOfMaterials), stat = IError)
            allocate(ContactSurfaceSoilElements(Counters%NEl), stat = IError)
          else
            allocate(LumpedMassSys(1), stat = IError)
            allocate(TotalVelocitySys(1), stat = IError)
            allocate(AccelerationSys(1), stat = IError)
            allocate(ContactNodeNormals(1, 1), stat = IError)
            allocate(NodeElement(1, 1), stat = IError)
            allocate(InterfaceNodes(1), stat = IError)
            allocate(InterfaceNodesAdhesion(1,1), stat = IError)    
            allocate(ContactSurfaceSoilElements(1), stat = IError)
          end if

          LumpedMassSys = 0.0
          TotalVelocitySys = 0.0
          AccelerationSys = 0.0
          ContactNodeNormals = 0.0
          InterfaceNodes = .false.
          NodeElement = 0
          InterfaceNodesAdhesion = 0.0
          ContactSurfaceSoilElements = .false.

        end subroutine InitialiseContactArrays


        subroutine ApplyMPMDynamicContact()
        !**********************************************************************
        !
        !    Function:  Do velocity correction to entity boundary nodes
        !     
        ! Note : This function works only if Number of Layers = 1
        !
        !**********************************************************************

        implicit none
        
        ! local variables
        integer(INTEGER_TYPE), dimension(Counters%nEntity,Counters%NodTot) :: EntityNodes      !for each entity a list of the global nodes belonging to elements containing particles belonging to this entity 
        integer(INTEGER_TYPE), dimension(CalParams%NumberOfMaterials,Counters%NodTot) :: MaterialNodes      
        real(REAL_TYPE), dimension(NVECTOR) :: Unormal      !unit normal vector at node for specific entity
        real(REAL_TYPE), dimension(NVECTOR) :: VelDiff      !Difference in entity and system velocity
        real(REAL_TYPE), dimension(3) :: CrossP      !cross product
        real(REAL_TYPE), dimension(3) :: Omega      !
        real(REAL_TYPE) :: AdhesionFactor, FrictionCoef, AdhesionCoef 
        integer(INTEGER_TYPE) :: I, nn, ient  
        integer(INTEGER_TYPE), dimension(NVECTOR) :: IDof
        real(REAL_TYPE) :: DiffDotUn, Traction, vsize, Dmu
        logical :: ApplyContact
        
        if (.not.(NFORMULATION==1)) RETURN
        if (.not.CalParams%ApplyContactAlgorithm) RETURN
          
        ! set up EntityNodes matrix for current particle configuration
        call SetEntityNod(EntityNodes)
        call SetMaterialNod(MaterialNodes)  
        
        ! loop through all nodes in model, in global coordinate system
        DO nn=1,Counters%NodTot
        
          if (InterfaceNodes(nn)) then  ! in case of pile driving the contact surface is defined
            
            !get storage position of node degree-of-freedoms
            do I = 1, NVECTOR
              IDof(I) = ReducedDof(nn) + I
            end do
            
            !loop through all entities
            do ient=1,Counters%nEntity
              !
              if (EntityNodes(ient,nn) == 1) then             !this node (nn) belongs to the entity (ient)

                ApplyContact = .false.
                do I = 1, NVECTOR
                  if (TotalVelocitySoil(IDof(I), ient) /= TotalVelocitySys(IDof(I))) then
                    ApplyContact = .true.
                    exit
                  end if
                end do
                  
                if (ApplyContact) then   !there is contact at this node because the entity velocity is different from system velocity

                !if (CalParams%ApplyMeshSmoothing) then ! the contact surface are predefined and kept defined
                  ! here, I am taking the normals to the pile (entity 2)
                  Unormal(:) = ContactNodeNormals(:, nn)
               
                  if (ient == SOIL_ENTITY) then ! if it is a soil node (entity 1) then change direction
                    Unormal = - Unormal
                  end if
                !end if

                  DiffDotUn = 0.0
                  do I = 1, NVECTOR
                    ! Calculate velocity difference between entity and system
                    VelDiff(I) = TotalVelocitySoil(IDof(I), ient) - TotalVelocitySys(IDof(I))
                    ! VelDiff.dot.Unormal
                    DiffDotUn = DiffDotUn + VelDiff(I) * Unormal(I)
                  end do  
                    
                    !get the traction at the node for this entity if the modified contact model is used
                    if (CalParams%ApplyTractionContact) then               !MODIFIED (TRACTION) CONTACT MODEL IS USED
                      call EntityNodeTraction(nn, ient, Unormal, Traction)
                    else                                                   !STANDARD CONTACT MODEL IS USED
                      Traction = 0.0                                       !set equal to zero so that (DiffDotUn > 0) is only criteria used below
                    end if 
                    !
                    !If the modified contact model is used, correction is done when Traction < 0 (compression)
                    ! The exeption is when it is the first contact and all stresses are still zero: 
                    !    - this will give Traction = 0
                    !    - In this case the standard contact model is used: (TotalVelocitySoil-TotalVelocitySys).dot. Unormal > 0
                    !
                    !If the standard contact model is used, Traction is set to zero (above) and the only condition used...
                    ! ...to check for contact is:  (TotalVelocitySoil-TotalVelocitySys).dot. Unormal > 0
                    !
                    do I = 1, CalParams%NumberOfMaterials-1 
                      if (MaterialNodes(I,nn)==1) then 
                         AdhesionCoef = InterfaceNodesAdhesion(nn,I)
						 exit
					  else
						  AdhesionCoef=1000
                      end if
                    end do                    
                    ! Adhesion
                    !AdhesionCoef = 0.1
                    AdhesionFactor = AdhesionCoef * CalParams%TimeIncrement / (LumpedMassDry(IDof(1), IEnt) * (DiffDotUn + TINY))
                    
                    if (DiffDotUn>0.0) then ! Standard formulation: DiffDotUn greater than 0 means approaching bodies

                      if (NVECTOR == 2) then ! 2D case
                        ! VelDiff x Unormal
                        CrossP(1) = 0.0
                        CrossP(2) = 0.0
                        CrossP(3) = VelDiff(1) * Unormal(2) - VelDiff(2) * Unormal(1)
                      elseif (NVECTOR == 3) then ! 3D case
                        !VelDiff x Unormal
                        CrossP(1) = VelDiff(2) * Unormal(3) - VelDiff(3) * Unormal(2)
                        CrossP(2) = VelDiff(3) * Unormal(1) - VelDiff(1) * Unormal(3)
                        CrossP(3) = VelDiff(1) * Unormal(2) - VelDiff(2) * Unormal(1)
                      else
                        call GiveError('Dimension is not correct. It must be 2 or 3.')
                      end if                        

                      vsize = sqrt( CrossP(1)*CrossP(1) + CrossP(2)*CrossP(2) + CrossP(3)*CrossP(3) )
                      
                      !CrossP/ |CrossP|
                      Omega = 0.0  
                      if (vsize /= 0) then
                        Omega = CrossP/vsize
                      end if
                       
                     !  if (CalParams%ApplyMeshSmoothing .and. CalParams%IStep==1) then ! do the initial stress phase with rough contact
                     !    Dmu = min(1.0, vsize/DiffDotUn)
                     !  else
                     ! Dmu will always be CalParams%fricCoef (because sliding occurs)
                      do I = 1, CalParams%NumberOfMaterials-1 
                        if (MaterialNodes(I,nn)==1) then 
                          FrictionCoef = FrictionVector(nn,I)
						  exit
						else !node belongs to contact but no material is in neighbour element
						  FrictionCoef = 10 ! A high value will enforce that the solution is vsys
                        end if
                      end do

                       Dmu = min(FrictionCoef + AdhesionFactor, vsize/DiffDotUn) ! Slope larger than friction coefficient, then bring back to limit (max tangential force)
                     ! end if

                       !Unormal X Omega - tangential direction normal to n-t-plane
                       if (NVECTOR == 2) then ! 2D case
                        !VelDiff x Unormal
                        CrossP(1) = Unormal(2) * Omega(3)
                        CrossP(2) = - Unormal(1) * Omega(3)
                        CrossP(3) = Unormal(1) * Omega(2) - Unormal(2) * Omega(1)
                      elseif (NVECTOR == 3)  then ! 3D case
                        !VelDiff x Unormal
                        CrossP(1) = Unormal(2) * Omega(3) - Unormal(3) * Omega(2)
                        CrossP(2) = Unormal(3) * Omega(1) - Unormal(1) * Omega(3)
                        CrossP(3) = Unormal(1) * Omega(2) - Unormal(2) * Omega(1)
                      else
                        call GiveError('Dimension is not correct. It must be 2 or 3.')
                      end if
                      
                      do I = 1, NVECTOR
                        ! correct the velocity for soil nodes
                        TotalVelocitySoil(IDof(I),ient) =  TotalVelocitySoil(IDof(I),ient) - DiffDotUn*( Unormal(I) + Dmu*CrossP(I) )
                        ! with the corrected, get the final acceleration
                        AccelerationSoil(IDof(I),ient) = (TotalVelocitySoil(IDof(I),ient) - TotalVelocitySoilPrevious(IDof(I),ient)) / CalParams%TimeIncrement
                      end do
                     
                   else
                   
                     ! Adhesion factor takes into account that in case of separation, the adhesive force might
                     ! still hold back the body. In this case, slip / no-slip are still relevant.
                     if (abs(AdhesionFactor)>=1.0) then
                       ! Bodies stick to each other
                       ! Set velocity of contacting entities to the system velocity
                     end if
                   end if  !(Traction == 0 .and. DiffDotUn > 0 .or......)
                end if     !(TotalVelocitySoil(IDof,ient) /= TotalVelocitySys(IDof)...)
              end if       !(EntityNodes(ient,nn) = 1)
            end do         !ient=1,nEntity
          end if ! interface node
        end do             !nn=1,NodTot
              
        if (IsFEMComputation()) then
          TotalVelocitySoil = TotalVelocitySoilPrevious ! global coordinate system (not changed), not used for MPM
        end if 
          
        end subroutine ApplyMPMDynamicContact
        

        subroutine SetEntityNod(EntityNodes)
        !**********************************************************************
        !
        !    Function:  To set the entity data for the current particle configuration    
        !     
        !
        ! O   EntityNodes   :  For each entity, gives the active nodes
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), dimension(Counters%nEntity,Counters%NodTot), intent(inout) :: EntityNodes
     
          ! local variables
          integer(INTEGER_TYPE) :: iel, iael, iEnt, inode, nn, I
          
          EntityNodes = 0
                          
          !EntityNodes
          !$omp parallel do private(iael, iel, ient, inode, nn)
          do iael=1,Counters%Nael                       !loop through all elements
            iel = ActiveElement(iael)
            do ient=1,Counters%nEntity                    !loop through all entities
              if (EntityElements(ient,iel) == 1) then !entity contains this element
                do inode=1,ELEMENTNODES              !loop through element nodes
                  nn = ElementConnectivities(inode,iel)                  !get global node number
                  EntityNodes(ient,nn) = 1              !set node as active node for entity
                end do
              end if
            end do
          end do       
          !$omp end parallel do          

          if (CalParams%ApplyContactMeshBoundary) then
            do I = 1, Counters%NodTot
              if (InterfaceNodes(I)) then
                EntityNodes(PILE_ENTITY, I) = 1
              end if
            end do
          end if

        end subroutine SetEntityNod
       
       
        subroutine SetMaterialNod(MaterialNodes) 
        !**********************************************************************
        !
        !    Function:  To set the material data for the current particle configuration    
        !     
        !
        ! O   MaterialNodes   :  For each material, gives the active nodes
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), dimension(CalParams%NumberOfMaterials,Counters%NodTot), intent(inout) :: MaterialNodes
     
          ! local variables
          integer(INTEGER_TYPE) :: iel, iael, iEnt, inode, nn
          
          MaterialNodes = 0

          do iael=1,Counters%Nael                       !loop through all elements
            iel = ActiveElement(iael)
            do ient=1,CalParams%NumberOfMaterials                    !loop through all materials
              if (MaterialElements(ient,iel) == 1) then !material contains this element
                do inode=1,ELEMENTNODES              !loop through element nodes
                  nn = ElementConnectivities(inode,iel)                  !get global node number
                  MaterialNodes(ient,nn) = 1              !set node as active node for material
                end do
              end if
            end do
          end do       

        end subroutine SetMaterialNod


        subroutine EntityNodeTraction(inode, ient, Unormal, Traction)
        !**********************************************************************
        !
        ! Function:  For a given node (inode) and entity (ient), calculate (interpolate from material points) the normal traction
        !
        ! I: inode         - global node number
        !   ient          - entity number
        !   Unormal       - the unit normal vector at node (inode) for entity (ient)
        !
        ! O: Traction       - Traction (scalar)
        !
        !**********************************************************************
        
        implicit none
        
        integer(INTEGER_TYPE), intent(in) :: inode, ient
        real(REAL_TYPE), dimension(NVECTOR),intent(in) :: Unormal
        real(REAL_TYPE), intent(out) :: Traction
        
        !local variables
        integer(INTEGER_TYPE) i,iel,nelempart,j,nn,iPart,IntGlo,nix
        real(REAL_TYPE), dimension(NTENSOR) :: Stress
        
      Traction = 0.0   !reset
      
      do i=1,99                                                                 !loop through all elements containing this global node (inode)
        iel = NodeElement(inode,i)                                              !element global number
        if (iel /= 0) then                                                    !NodeElement might contain zero's if node has little elements
          if (EntityElements(ient,iel) == 1) then                             !the enity contains this element
            nelempart = NPartEle(iel)                                           !number of mass points in element
            do j=1,ELEMENTNODES                                              !loop through element nodes
              nn=ElementConnectivities(j,iel)                                                    !get global node number
              if (nn == inode) then                                           !this is the element node corresponding to the global node
                do iPart=1,nelempart                                            !loop through all material points in element
                  IntGlo = GetParticleIndex(iPart,iel)                          !particle global number
                  if (EntityIDArray(IntGlo) == ient ) then               !check if particle belongs to entity
                    Stress = SigmaEffArray(IntGlo,:)

                    if (NTENSOR == 4) then ! 2D case
                      Traction = Traction + (Stress(1) * Unormal(1) * Unormal(1) + &!...the mass is used as interpolation weight
                                             Stress(2) * Unormal(2) * Unormal(2) + &
                                             Stress(4) * Unormal(1) * Unormal(2)) * &
                                             MassArray(IntGlo) * ShapeValuesArray(IntGlo,j)
                    elseif (NTENSOR == 6) then ! 3D case
                      Traction = Traction + (Stress(1)*Unormal(1)*Unormal(1) + & !...the mass is used as interpolation weight
                                             Stress(4)*Unormal(1)*Unormal(2) + &
                                             Stress(6)*Unormal(1)*Unormal(3) + &
                                             Stress(4)*Unormal(2)*Unormal(1) + &
                                             Stress(2)*Unormal(2)*Unormal(2) + &
                                             Stress(5)*Unormal(2)*Unormal(3) + &
                                             Stress(6)*Unormal(3)*Unormal(1) + &
                                             Stress(5)*Unormal(3)*Unormal(2) + &
                                             Stress(3)*Unormal(3)*Unormal(3) )* &
                                             MassArray(IntGlo) * ShapeValuesArray(IntGlo,j)
                    else
                      call GiveError('Dimension is not correct. It must be 2 or 3.')
                    end if                    
                  end if
                end do 
                exit        !element node corresponding to global node already found (can only be one) - go to next element
              end if
            end do
          end if
        end if
      end do
      
      nix=ReducedDof(inode)+1                      !global storage coordinate of x-val dof of the global node
      Traction = Traction/LumpedMassDry(nix,ient)     !divide the Traction by the nodal mass
      
        end subroutine EntityNodeTraction

        
        subroutine CalculateNodeElement()
        !**********************************************************************
        !
        !  Function:  To calculate the NodeElement array
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IEl, INode, NN, I
          
          if (.not.CalParams%ApplyContactAlgorithm) RETURN

          do IEl = 1, Counters%NEl
            do INode=1, ELEMENTNODES
              NN = ElementConnectivities(INode, IEl)
              I = NodeElement(NN, 100) + 1
              NodeElement(NN, I) = IEl
              NodeElement(NN, 100) = I
           end do
         end do
       
      end subroutine CalculateNodeElement
      
        
        subroutine ReadContactData()
        !**********************************************************************
        !
        ! Function: Reads contact data, depending on GOM file version
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables       
          character(len=MAX_FILENAME_LENGTH) :: FileName 
          integer(INTEGER_TYPE) :: FileUnit, NContVolElem
          character(len=255) :: BName, AName
          integer(INTEGER_TYPE) :: ios ! used for error control
          
          FileName = trim(CalParams%FileNames%ProjectName)//GOM_FILE_EXTENSION
          FileUnit = TMP_UNIT
          
          call FileOpen(FileUnit, trim(FileName))
          
          do  
            read(FileUnit,*) AName
              if ((trim(AName)=='$$START_CONTACT_VOLUME').or.(trim(AName)=='$$START_BODY_CONTACT_2D')) then
                read(FileUnit, *) NContVolElem  ! no of element sides forming the contact surface
                if (NContVolElem>0) then
                  read(FileUnit,*) StructureElement
                  CalParams%MovingMesh%StructureMaterialID = (abs(ElementMaterialID(StructureElement)))
                  EXIT
                end if
              else if (trim(AName)=='$$FINISH') then
                EXIT
              end if
          end do
          
          close(FileUnit)
          
          if (.not.CalParams%ApplyContactAlgorithm) RETURN
          
          ! check if GOM file exists in project folder, otherwise give error and stop execution
          if ( FExist(trim(FileName)) ) then
            call GiveMessage('Reading contact data from GOM file: ' // trim(FileName) )  
          else
            call GiveError('GOM file does not exist!' // NEW_LINE('A') // 'required GOM file: ' // trim(FileName) )
          end if

          ! open GOM file
          call FileOpen(FileUnit, trim(FileName))
          
          ! determine current version of GOM file 
          read(FileUnit, '(A)', iostat=ios) BName ! NB: if no version is specified in the header of the GOM file, the default case will be chosen automatically below
          call Assert( ios == 0, 'GOM file: Can''t read flag from GOM file.' )
          FileVersion = trim(BName)
          
          if (NDIM==2) then 
            call ContactSurfaceNodes_2D(FileUnit)     
          else           
            select case (FileVersion) ! read GOM data depending on file version
              case (Anura3D_v2021)              
                call ContactSurfaceNodes_v2021(FileUnit) 
              end select
          end if    

          ! close GOM file    
          close(FileUnit)
      
        end subroutine ReadContactData                                   
        
        
        subroutine ComputeInterfaceNodesAdhesion2D()
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: NodeID
          integer(INTEGER_TYPE) :: IElement, I, J, K, IError     
          integer(INTEGER_TYPE), dimension(NVECTOR) :: SurfaceNodes
          real(REAL_TYPE), dimension(NVECTOR) :: Vec1
          real(REAL_TYPE) :: ElementSurface
          real(REAL_TYPE), dimension(:,:), allocatable :: AdhesionVector
          integer(INTEGER_TYPE) :: CornerNodes
          
          if (NDIM /= 2) RETURN
          
          CornerNodes = NDIM ! this holds only for surface element of tetrahedron and triangle
          
          if (.not.CalParams%ApplyContactAlgorithm) RETURN

          allocate(AdhesionVector(Counters%NodTot, CalParams%NumberOfMaterials), stat=IError) 
          allocate(FrictionVector(Counters%NodTot, CalParams%NumberOfMaterials), stat=IError) 
          
          FrictionVector = 0.0
          AdhesionVector = 0.0
          InterfaceNodesAdhesion = 0.0
        
          do IElement = 1, NContactElm
            do I = 1, CornerNodes
              SurfaceNodes(I) = ContactSurfaceElmNodes(IElement, I + 2)
            end do

          do I = 1, CornerNodes
            Vec1(I) = NodalOriginalCoord(SurfaceNodes(2), I) - NodalOriginalCoord(SurfaceNodes(1), I)
          end do
            ElementSurface = sqrt(Vec1(1)*Vec1(1)+Vec1(2)*Vec1(2)) !2.0 is number of corner nodes of a line
         
            
            if (CalParams%NumberOfMaterials>1) then 
              do J = 1, CalParams%NumberOfMaterials-1 
                do I = 1, CornerNodes
                  InterfaceNodesAdhesion(SurfaceNodes(I),J) = InterfaceNodesAdhesion(SurfaceNodes(I),J) + ElementSurface
                  if ((CalParams%FricCoef==0).and.(allocated(ContactSurfaces))) then
                      FrictionVector(SurfaceNodes(I),J) = ContactSurfaces(IElement,3+J*2)     !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                  else
                      FrictionVector(SurfaceNodes(I),J) = CalParams%FricCoef
                  end if
                  if ((CalParams%Adhesion==0).and.(allocated(ContactSurfaces))) then
                      AdhesionVector(SurfaceNodes(I),J) = ContactSurfaces(IElement,4+J*2)  !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                  else
                      AdhesionVector(SurfaceNodes(I),J) = CalParams%Adhesion
                  end if
                end do
              end do   
            else 
               do I = 1, CornerNodes
                  InterfaceNodesAdhesion(SurfaceNodes(I),1) = InterfaceNodesAdhesion(SurfaceNodes(I),1) + ElementSurface  
                 if ((CalParams%FricCoef==0).and.(allocated(ContactSurfaces))) then
                   FrictionVector(SurfaceNodes(I),1) = ContactSurfaces(IElement,7+1*2)       !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                 else
                   FrictionVector(SurfaceNodes(I),1) = CalParams%FricCoef
                 end if
                 if ((CalParams%Adhesion==0).and.(allocated(ContactSurfaces))) then
                   AdhesionVector(SurfaceNodes(I),1) = ContactSurfaces(IElement,8+1*2)     !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                 else
                   AdhesionVector(SurfaceNodes(I),1) = CalParams%Adhesion
                 end if
               end do
            end if
          end do
                   
          do K = 1, CalParams%NumberOfMaterials 
            do NodeID = 1, Counters%NodTot
            if (InterfaceNodes(NodeID)) then
              InterfaceNodesAdhesion(NodeID,K) = InterfaceNodesAdhesion(NodeID,K) * AdhesionVector(NodeID,K)
            end if
            end do
          end do
        
        end subroutine ComputeInterfaceNodesAdhesion2D
        

        subroutine ComputeInterfaceNodesAdhesion()
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: NodeID
          integer(INTEGER_TYPE) :: IElement, I, J, K, IError     
          integer(INTEGER_TYPE), dimension(NVECTOR) :: SurfaceNodes
          real(REAL_TYPE), dimension(NVECTOR) :: Vec1, Vec2
          real(REAL_TYPE), dimension(NVECTOR) :: Normal
          real(REAL_TYPE) :: ElementSurface
          real(REAL_TYPE), dimension(:,:), allocatable :: AdhesionVector
          integer(INTEGER_TYPE) :: CornerNodes
          
          if (NDIM /= 3) RETURN
          
          CornerNodes = NDIM ! this holds only for surface element of tetrahedron and triangle
          
          if (.not.CalParams%ApplyContactAlgorithm) RETURN

          allocate(AdhesionVector(Counters%NodTot, CalParams%NumberOfMaterials), stat=IError) 
          allocate(FrictionVector(Counters%NodTot, CalParams%NumberOfMaterials), stat=IError) 
          
          FrictionVector = 0.0
          AdhesionVector = 0.0
          InterfaceNodesAdhesion = 0.0
        
          do IElement = 1, NContactElm
            do I = 1, CornerNodes
              SurfaceNodes(I) = ContactSurfaceElmNodes(IElement, I + 2)
            end do

            if (NDIM == 3) then ! 3D case
              do I = 1, CornerNodes
                Vec1(I) = NodalOriginalCoord(SurfaceNodes(2), I) - NodalOriginalCoord(SurfaceNodes(1), I)
                Vec2(I) = NodalOriginalCoord(SurfaceNodes(3), I) - NodalOriginalCoord(SurfaceNodes(1), I)
              end do
              Normal = CrossProduct(Vec1, Vec2)
              ElementSurface = Length(Normal, 3) / 2.0 / 3.0 !3.0 is number of vertices (corner nodes) of a triangle
            else ! 2D case
              do I = 1, CornerNodes
                Vec1(I) = NodalOriginalCoord(SurfaceNodes(2), I) - NodalOriginalCoord(SurfaceNodes(1), I)
              end do
              ElementSurface = sqrt(Vec1(1)*Vec1(1)+Vec1(2)*Vec1(2)) / 2.0 !2.0 is number of corner nodes of a line
            endif            
            
            if (CalParams%NumberOfMaterials>1) then 
              do J = 1, (CalParams%NumberOfMaterials-1)  
                do I = 1, CornerNodes
                  InterfaceNodesAdhesion(SurfaceNodes(I),J) = InterfaceNodesAdhesion(SurfaceNodes(I),J) + ElementSurface
                  if ((CalParams%FricCoef==0).and.(allocated(ContactSurfaces))) then
                      FrictionVector(SurfaceNodes(I),J) = ContactSurfaces(IElement,7+J*2)     !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                  else
                      FrictionVector(SurfaceNodes(I),J) = CalParams%FricCoef
                  end if
                  if ((CalParams%Adhesion==0).and.(allocated(ContactSurfaces))) then
                      AdhesionVector(SurfaceNodes(I),J) = ContactSurfaces(IElement,8+J*2)  !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                  else
                      AdhesionVector(SurfaceNodes(I),J) = CalParams%Adhesion
                  end if
                end do
              end do   
            else 
               do I = 1, CornerNodes
                  InterfaceNodesAdhesion(SurfaceNodes(I),1) = InterfaceNodesAdhesion(SurfaceNodes(I),1) + ElementSurface  
                 if ((CalParams%FricCoef==0).and.(allocated(ContactSurfaces))) then
                   FrictionVector(SurfaceNodes(I),1) = ContactSurfaces(IElement,7+1*2)       !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                 else
                   FrictionVector(SurfaceNodes(I),1) = CalParams%FricCoef
                 end if
                 if ((CalParams%Adhesion==0).and.(allocated(ContactSurfaces))) then
                   AdhesionVector(SurfaceNodes(I),1) = ContactSurfaces(IElement,8+1*2)     !!!!!!!!!!!!!!!!!!!!!!!!!!!check numbers
                 else
                   AdhesionVector(SurfaceNodes(I),1) = CalParams%Adhesion
                 end if
               end do
            end if
          end do
                   
          do K = 1, CalParams%NumberOfMaterials 
            do NodeID = 1, Counters%NodTot
            if (InterfaceNodes(NodeID)) then
              InterfaceNodesAdhesion(NodeID,K) = InterfaceNodesAdhesion(NodeID,K) * AdhesionVector(NodeID,K)
            end if
            end do
          end do
        
        end subroutine ComputeInterfaceNodesAdhesion
        

        subroutine ContactSurfaceNodes_v2021(FileUnit)
        !**********************************************************************
        !
        !    Function:  Contains code for defining contact nodes and computing the node normals
        !
        !  Note: This function is only for 3D case
        !
        !        GOM version 2021 
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! local variables
          character :: BName*255
          integer(INTEGER_TYPE) :: IError, ISurf, J, NContSur, NodID, I, IEl, NPileelements, NFoundNodes, CheckNodeID, &
                    ElementID, SideID, AdjacentElement, IVol, ISide, ContElem, INode, JNode, NSideNodes, &
                    CheckNodeI, II, CountSurf, DummyElem, IMat, JMat, NeighbourElement, MaterialID1, MaterialID2, AI
          real(REAL_TYPE), dimension(:,:), allocatable :: TotSrfaceAndLoad, ContactVolumes, ContactVolumesDummy
          integer(INTEGER_TYPE), dimension(:,:), allocatable :: ContSurfaceNodes, ContVolSurf
          integer(INTEGER_TYPE), dimension(:), allocatable :: ContctSurfaceElmID, BoundaryElements, Sides
          logical, dimension(:), allocatable :: Element, CheckedLoadedSurface, IsContactSurface
          integer(INTEGER_TYPE) :: NFoundSides

          allocate(Element(Counters%NEl), stat=IError)
          allocate(GlobContElement(Counters%NEl), stat=IError) 
        
          if (allocated(ContactSurfaces) ) then
            deallocate(ContactSurfaces, stat = IError)
          end if
          
          if (allocated(ContactSurfaceElmNodes) ) then
            deallocate(ContactSurfaceElmNodes, stat = IError)
          end if        

          do
            read(FileUnit, '(A)') BName
            if ((trim(BName)=='$$START_CONTACT_SURFACE')) then
              read(FileUnit, *) NContSur  ! no of element sides forming the contact surface
              if (NContSur>0) then
                allocate(ContactSurfaces(NContSur, 16), stat=IError)
                ContactSurfaces = 0.0 
                 do ISurf= 1, NContSur
                   read(FileUnit,*)  (ContactSurfaces(ISurf,J), J = 1,16)  ! (ElementID, 6 NodeIDs, no. of slaves, friction and adhesion for each slave)
                 end do
              end if
            else if (trim(BName)=='$$START_CONTACT_VOLUME') then 
                read(FileUnit, *) NContVolElem  
                allocate(ContactVolumes(NContVolElem,10), stat=IError)
                allocate(ContactVolumesDummy(NContVolElem,10), stat=IError)
                allocate(ContactMaterialName(NContVolElem,4), stat=IError)
                ContactVolumes = 0
                ContactVolumesDummy = 0
                do IVol = 1, NContVolElem
                    read(FileUnit,*) ContactVolumes(IVol,1:2), ContactMaterialName(Ivol,1), ContactVolumes(IVol,3:4), &
                    ContactMaterialName(Ivol,2), ContactVolumes(IVol,5:6), ContactMaterialName(Ivol,3), &
                    ContactVolumes(IVol,7:8), ContactMaterialName(Ivol,4), ContactVolumes(IVol,9:10)
                end do 
            else if (trim(BName)=='$$FINISH') then
                  EXIT
             end if 
          end do 
          ! rearrange contact properties due to MaterialIndex
          do IMat = 1, 4 !loop over "sliding" materials
            do JMat = 1, CalParams%NumberOfMaterials  !loop over all materials  
              if (ContactMaterialName(1,IMat)==MatParams(JMat)%MaterialName) then
                ContactVolumesDummy(:, MatParams(JMat)%MaterialIndex * 2 + 1 : MatParams(JMat)%MaterialIndex * 2 + 2) &      
                = ContactVolumes(:,(IMat*2+1):(IMat*2+2))     
              end if
            end do
          end do  
          ContactVolumes(:,3:10) = ContactVolumesDummy(:,3:10)
           
          allocate(ContVolSurf(4*NContVolElem,6), stat=IError) 
          ContVolSurf = 0
          Element = .False.
          GlobContElement = .False.
          NPileelements = 0
          NContSur = 0
          DummyElem = 0
          do IEl = 1, NContVolElem
             NPileelements = NPileelements+1 
             Element(ContactVolumes(IEl,1)) = .True. ! pile elemen      
             GlobContElement(ContactVolumes(IEl,1)) = .True.
               if (abs(ElementMaterialID(ContactVolumes(IEl,1)))/=CalParams%MovingMesh%StructureMaterialID) then
                 call GiveWarning('Body contact can not be assigned to more than one material')
               end if
             PileEntityMaterial = abs(ElementMaterialID(ContactVolumes(IEl,1)))
             if  (IsActiveElement(ContactVolumes(IEl,1))) then ! element belongs to the contact volume and is active
               DummyElem = ContactVolumes(IEl,1)    
               do ISide = 1, 4
                 DummyElem = ContactVolumes(IEl,1)  
                 ContElem = BoundaryElementSurface(DummyElem, ISide, &
                                                IsActiveElement, Counters%NEl)
        !Returns 1 if the element adjacent to side ISide of IElement is deactivated.
                 if (ContElem==1) then 
                   do INode = 1, 6 ! Loop over nodes of side I
                     CheckNodeI = ElementConnectivities10Node(DetermineSideNodesTetrahedronHOE(ISide, INode), ContactVolumes(IEl,1))
                     ContVolSurf((IEl+(ISide-1)+3*(IEl-1)), INode) = CheckNodeI
                   end do
                   NContSur = NContSur + 1
        !Returns 0 if side ISide of IElement lies inside a group of activated elements.  
                 else if (ContElem==0) then
                   NeighbourElement = GetAdjacentElement(DummyElem,ISide)  
                   if ((NeighbourElement/=0).or.(NeighbourElement/=-999)) then
                     MaterialID1 = ElementMaterialID(NeighbourElement)
                     MaterialID2 = ElementMaterialID(DummyElem)! get ParticleID from first particle in DummyElem
                   end if  
                      if (MaterialID1/=MaterialID2) then
                        do JNode = 1, 6  
                          CheckNodeI = ElementConnectivities10Node(DetermineSideNodesTetrahedronHOE(ISide, JNode), ContactVolumes(IEl,1))
                          ContVolSurf((IEl+(ISide-1)+3*(IEl-1)),JNode) = CheckNodeI  
                        end do
                        NContSur = NContSur + 1
                      end if
                 end if   
               end do
             end if
          end do ! loop over elements
          
          allocate(ContactSurfaces(NContSur, 16), stat=IError)
          ContactSurfaces = 0
          CountSurf = 0
          
          do II = 1, 4*NContVolElem
              if (ContVolSurf(II,1)/=0) then
                  AI = INT(II/4)+1
                  CountSurf = CountSurf + 1
                  ContactSurfaces(CountSurf,3:8) = ContVolSurf(II,1:6)
                  ContactSurfaces(CountSurf,9:16) = ContactVolumes(AI,3:10)
              end if
          end do        
             
          allocate(IsContactSurface(NContSur), stat=IError)
          IsContactSurface = .true.
          
          allocate(ContSurfaceNodes(NContSur, 6), stat=IError)
          ContSurfaceNodes= 0.0
          
          do I = 1, NContSur
              do J = 1, 6
                  ContSurfaceNodes(I,J) = ContactSurfaces(I,J+2)
              end do
          end do
          
          
          allocate(ContctSurfaceElmID(NContSur), stat=IError)
          ContctSurfaceElmID= 0
          
          allocate(Sides(NContSur), stat=IError)
          Sides= 0         
                    
          
         if (CalParams%ApplyContactMeshBoundary) then

            allocate(CheckedLoadedSurface(NContSur), stat = IError)
            CheckedLoadedSurface = .false.
            NFoundSides = 0
           do IEl = 1, Counters%NEl
             do I = 1, 4
               AdjacentElement = GetAdjacentElement(IEl, I) 
               if (AdjacentElement==0) then  
                 do ISurf= 1, NContSur
                   if (.not.CheckedLoadedSurface(ISurf)) then
                     NFoundNodes = 0
                     do J = 1, 6 ! Loop over nodes of side I
                       CheckNodeID = ElementConnectivities10Node(DetermineSideNodesTetrahedronHOE(I, J), IEl)
                       if ( (CheckNodeID== ContSurfaceNodes (ISurf, 1)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 2)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 3)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 4)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 5)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 6)))then
                         NFoundNodes = NFoundNodes + 1
                       end if
                     end do
                     if (NFoundNodes==6) then
                       ContctSurfaceElmID(ISurf) = IEl
                       CheckedLoadedSurface(ISurf) = .true.
                       Sides(ISurf) = I
                       NFoundSides = NFoundSides + 1
                       EXIT
                     end if
                   end if
                 end do
               end if
               if (NFoundSides==NContSur) then
                 EXIT
               end if
             end do
             if (NFoundSides==NContSur) then
               EXIT
             end if
           end do
           
           deallocate(CheckedLoadedSurface, stat = IError)
       else
         do ISurf= 1, NContSur
           do IEl = 1, Counters%NEl
             if (Element(IEl)) then 
               do I = 1, 4 ! Loop over sides
                 NFoundNodes = 0
                  do J = 1, 6 ! Loop over nodes of side I
                    CheckNodeID = ElementConnectivities10Node &
                                  (DetermineSideNodesTetrahedronHOE(I, J), IEl)
                if ( (CheckNodeID== ContSurfaceNodes (ISurf, 1)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 2)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 3)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 4)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 5)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 6)))then
                   NFoundNodes = NFoundNodes + 1
                 end if
                 end do
                 if (NFoundNodes==6) then
                  ContctSurfaceElmID(ISurf) = IEl
                  Sides(ISurf) = I
                  EXIT
                 end if
               end do ! sides
              end if ! pile
            end do ! elements
          end do ! contact surfaces
       end if   
       
       allocate(ContactSurfaceElmNodes(NContSur, 8), stat=IError)
       ContactSurfaceElmNodes= 0
        
       NContactElm = NContSur
       
       do ISurf = 1, NContSur
          ContactSurfaceElmNodes (ISurf, 1)= &
                                    INT(ContctSurfaceElmID(ISurf)) ! element ID
            ContactSurfaceElmNodes (ISurf, 2)= Sides (ISurf)
          do J = 1, 6
            ContactSurfaceElmNodes (ISurf, J+2)= &
                                              ContSurfaceNodes(ISurf, J)
          end do
       end do
       
       
          if (allocated(TotSrfaceAndLoad) ) then
            deallocate(TotSrfaceAndLoad, stat = IError)
          end if
          
          if (allocated(Sides) ) then
            deallocate(Sides, stat = IError)
          end if
          
          if (allocated(IsContactSurface) ) then
            deallocate(IsContactSurface, stat = IError)
          end if
          
          if (allocated(ContSurfaceNodes) ) then
            deallocate(ContSurfaceNodes, stat = IError)
          end if
          
          if (allocated(Element) ) then
            deallocate(Element, stat = IError)
          end if
          
          if (allocated(ContctSurfaceElmID) ) then
            deallocate(ContctSurfaceElmID, stat = IError)
          end if

            if (allocated(ContactNodeNormals)) then
              deallocate(ContactNodeNormals, stat = IError)
            end if
            if (allocated(InterfaceNodes)) then
              deallocate(InterfaceNodes, stat = IError)
            end if

            if (NContactElm==0) RETURN

            allocate(ContactNodeNormals(3, Counters%NodTot), stat = IError)
            ContactNodeNormals = 0.0
            allocate(InterfaceNodes(Counters%NodTot), stat = IError)
            InterfaceNodes = .false.
       
               do I = 1, NContactElm
                do J = 1, 3 ! corner nodes
                  NodID = ContactSurfaceElmNodes (I, J+2)
                  InterfaceNodes(NodID) = .true.
                end do
               end do
            NSideNodes = GetNSideNodes(ELEMENTNODES)

            allocate(BoundaryElements(Counters%NEl), stat = IError)
            BoundaryElements = 0
              do I = 1, NContactElm
               ElementID = ContactSurfaceElmNodes (I, 1)
               SideID = ContactSurfaceElmNodes (I, 2)
               BoundaryElements(ElementID) = SideID
              end do

            call DetermineNodeLocalCS(Counters%NodTot, NDIM, &
                                      ELEMENTNODES, Counters%NEl, &
                                      NSideNodes, &
                                      NodalCoordinates, ElementConnectivities, &
                                      BoundaryElements, &
                                      ContactNodeNormals)

            NInterfaceNodes = 0
            do I = 1, Counters%NodTot
              if (InterfaceNodes(I)) then
                NInterfaceNodes = NInterfaceNodes + 1
                if (CalParams%ApplyContactMeshBoundary) then
                  if ((ContactNodeNormals(1, I)>0.0).and.((NodalCoordinates(I, 2))>1E-5)) then
                    ContactNodeNormals(2, I) = 0.0
                    ContactNodeNormals(1:3, I) = VectorNorm(ContactNodeNormals(1:3, I), 3) ! generalise dimension
                  end if
                
                  ContactNodeNormals(1:3, I) = -1.0 * ContactNodeNormals(1:3, I)
                end if
              end if
            end do

            deallocate(BoundaryElements, stat = IError)               
               
         end subroutine ContactSurfaceNodes_v2021
        
        
        subroutine ContactSurfaceNodes_2D(FileUnit)
        !**********************************************************************
        !
        !    Function:  Contains code for defining contact nodes and computing the node normals
        !
        !  Note: This function is only for 2D case
        !
        !        GOM version 2021
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! local variables
          character :: BName*255
          integer(INTEGER_TYPE) :: IError, ISurf, J, NContSur, NodID, I, IEl, NPileelements, NFoundNodes, CheckNodeID, &
                    ElementID, SideID, AdjacentElement, IVol, ISide, ContElem, INode, JNode, NSideNodes, &
                    CheckNodeI, II, CountSurf, DummyElem, IMat, JMat, NeighbourElement, MaterialID1, MaterialID2, AI, IRef, JRef, KRef, NRef
          real(REAL_TYPE), dimension(:,:), allocatable :: TotSrfaceAndLoad, ContactVolumes, ContactVolumesDummy, NNormals, NNormals2,MatNodesCood, DVectors
          integer(INTEGER_TYPE), dimension(:,:), allocatable :: ContSurfaceNodes, ContVolSurf
          integer(INTEGER_TYPE), dimension(:), allocatable :: ContctSurfaceElmID, BoundaryElements, Sides, NNodes, MatNodes
          logical, dimension(:), allocatable :: Element, CheckedLoadedSurface, IsContactSurface
          integer(INTEGER_TYPE) :: NFoundSides
          real :: dummy ,dummy2, dummy3

          allocate(Element(Counters%NEl), stat=IError)
          allocate(GlobContElement(Counters%NEl), stat=IError) 
        
          if (allocated(ContactSurfaces) ) then
            deallocate(ContactSurfaces, stat = IError)
          end if
          
          if (allocated(ContactSurfaceElmNodes) ) then
            deallocate(ContactSurfaceElmNodes, stat = IError)
          end if        

          do
            read(FileUnit, '(A)') BName          
            if (trim(BName)=='$$START_BODY_CONTACT_2D') then 
                read(FileUnit, *) NContVolElem  
                allocate(ContactVolumes(NContVolElem,10), stat=IError)
                allocate(ContactVolumesDummy(NContVolElem,10), stat=IError)
                allocate(ContactMaterialName(NContVolElem,4), stat=IError)
                ContactVolumes = 0
                ContactVolumesDummy = 0
                do IVol = 1, NContVolElem
                    read(FileUnit,*) ContactVolumes(IVol,1:2), ContactMaterialName(Ivol,1), ContactVolumes(IVol,3:4), &
                    ContactMaterialName(Ivol,2), ContactVolumes(IVol,5:6), ContactMaterialName(Ivol,3), &
                    ContactVolumes(IVol,7:8), ContactMaterialName(Ivol,4), ContactVolumes(IVol,9:10)
                end do 
            else if (trim(BName)=='$$FINISH') then
                  EXIT
             end if 
          end do 
          ! rearrange contact properties due to MaterialIndex
          do IMat = 1, 4 !loop over "sliding" materials
            do JMat = 1, CalParams%NumberOfMaterials  !loop over all materials  
              if (ContactMaterialName(1,IMat)==MatParams(JMat)%MaterialName) then
                ContactVolumesDummy(:, MatParams(JMat)%MaterialIndex * 2 + 1 : MatParams(JMat)%MaterialIndex * 2 + 2) &      
                = ContactVolumes(:,(IMat*2+1):(IMat*2+2))     
              end if
            end do
          end do  
          ContactVolumes(:,3:10) = ContactVolumesDummy(:,3:10)  !until here 2D equals 3D
           
          allocate(ContVolSurf(3*NContVolElem,2), stat=IError) 
          ContVolSurf = 0
          Element = .False.
          GlobContElement = .False.
          NPileelements = 0
          NContSur = 0
          DummyElem = 0
          do IEl = 1, NContVolElem
             NPileelements = NPileelements+1 
             Element(ContactVolumes(IEl,1)) = .True. ! pile elemen      
             GlobContElement(ContactVolumes(IEl,1)) = .True.
               if (abs(ElementMaterialID(ContactVolumes(IEl,1)))/=CalParams%MovingMesh%StructureMaterialID) then
                 call GiveWarning('Body contact can not be assigned to more than one material')
               end if
             PileEntityMaterial = abs(ElementMaterialID(ContactVolumes(IEl,1)))
             if  (IsActiveElement(ContactVolumes(IEl,1))) then ! element belongs to the contact volume and is active
               DummyElem = ContactVolumes(IEl,1)    
               do ISide = 1, 3
                 DummyElem = ContactVolumes(IEl,1)  
                 ContElem = BoundaryElementSurface(DummyElem, ISide, &
                                                IsActiveElement, Counters%NEl)
        !Returns 1 if the element adjacent to side ISide of IElement is deactivated.
                 if (ContElem==1) then 
                   do INode = 1, 2 ! Loop over nodes of side I
                     if (ISide == 1) then 
                        CheckNodeI = ElementConnectivities(INode, ContactVolumes(IEl,1)) !in case that 1-2, 3-1, 2-3 is side order
                     elseif (ISide == 2) then
                        CheckNodeI = ElementConnectivities(3-2*(INode-1), ContactVolumes(IEl,1))  !in case that 1-2-3-1 is side order
                     elseif (ISide == 3) then
                        CheckNodeI = ElementConnectivities(INode+1, ContactVolumes(IEl,1))  !in case that 1-2-3-1 is side order 
                     end if   
                     ContVolSurf((IEl+(ISide-1)+2*(IEl-1)), INode) = CheckNodeI
                   end do
                   NContSur = NContSur + 1
        !Returns 0 if side ISide of IElement lies inside a group of activated elements.  
                 else if (ContElem==0) then
                   NeighbourElement = GetAdjacentElement(DummyElem,ISide)  
                   if ((NeighbourElement/=0).or.(NeighbourElement/=-999)) then
                     MaterialID1 = ElementMaterialID(NeighbourElement)
                     MaterialID2 = ElementMaterialID(DummyElem)! get ParticleID from first particle in DummyElem
                   end if  
                      if (MaterialID1/=MaterialID2) then  ! check if 
                        do JNode = 1, 2  
                          if (ISide == 1) then 
                            CheckNodeI = ElementConnectivities(JNode+(ISide-1), ContactVolumes(IEl,1)) !in case that 1-2-3-1 is side order
                          elseif (ISide == 2) then
                            CheckNodeI = ElementConnectivities(3-2*(JNode-1), ContactVolumes(IEl,1))  !in case that 1-2-3-1 is side order
                          elseif (ISide == 3) then
                            CheckNodeI = ElementConnectivities(JNode+1, ContactVolumes(IEl,1))  !in case that 1-2-3-1 is side order 
                          end if   
                          ContVolSurf((IEl+(ISide-1)+2*(IEl-1)), JNode) = CheckNodeI
                        end do
                        NContSur = NContSur + 1
                      end if
                 end if   
               end do
             end if
          end do ! loop over elements
          
          allocate(ContactSurfaces(NContSur, 12), stat=IError)
          ContactSurfaces = 0
          CountSurf = 0
          
          do II = 1, 3*NContVolElem
              if (ContVolSurf(II,1)/=0) then
                  AI = CEILING(II/3.0)
                  CountSurf = CountSurf + 1
                  ContactSurfaces(CountSurf,3:4) = ContVolSurf(II,1:2)
                  ContactSurfaces(CountSurf,5:12) = ContactVolumes(AI,3:10)
                  ContactSurfaces(CountSurf,1:2) = ContactVolumes(AI,1:2)
              end if
          end do        
             
          allocate(IsContactSurface(NContSur), stat=IError)
          IsContactSurface = .true.
          
          allocate(ContSurfaceNodes(NContSur, 2), stat=IError)
          ContSurfaceNodes= 0.0
          
          do I = 1, NContSur
              do J = 1, 2
                  ContSurfaceNodes(I,J) = ContactSurfaces(I,J+2)
              end do
          end do
          
          
          allocate(ContctSurfaceElmID(NContSur), stat=IError)
          ContctSurfaceElmID= 0
          
          allocate(Sides(NContSur), stat=IError)
          Sides= 0         
                    
          
         if (CalParams%ApplyContactMeshBoundary) then

            allocate(CheckedLoadedSurface(NContSur), stat = IError)
            CheckedLoadedSurface = .false.
            NFoundSides = 0
           do IEl = 1, Counters%NEl
             do I = 1, 4
               AdjacentElement = GetAdjacentElement(IEl, I) 
               if (AdjacentElement==0) then  
                 do ISurf= 1, NContSur
                   if (.not.CheckedLoadedSurface(ISurf)) then
                     NFoundNodes = 0
                     do J = 1, 6 ! Loop over nodes of side I
                       CheckNodeID = ElementConnectivities10Node(DetermineSideNodesTetrahedronHOE(I, J), IEl)
                       if ( (CheckNodeID== ContSurfaceNodes (ISurf, 1)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 2)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 3)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 4)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 5)).or. &
                         (CheckNodeID== ContSurfaceNodes (ISurf, 6)))then
                         NFoundNodes = NFoundNodes + 1
                       end if
                     end do
                     if (NFoundNodes==6) then
                       ContctSurfaceElmID(ISurf) = IEl
                       CheckedLoadedSurface(ISurf) = .true.
                       Sides(ISurf) = I
                       NFoundSides = NFoundSides + 1
                       EXIT
                     end if
                   end if
                 end do
               end if
               if (NFoundSides==NContSur) then
                 EXIT
               end if
             end do
             if (NFoundSides==NContSur) then
               EXIT
             end if
           end do
           
           deallocate(CheckedLoadedSurface, stat = IError)
       else
         do ISurf= 1, NContSur
           do IEl = 1, Counters%NEl
             if (Element(IEl)) then 
               do I = 1, 3 ! Loop over sides
                 NFoundNodes = 0
                  do J = 1, 2 ! Loop over nodes of side I
                    if (I<3) then                                                        ! THIS IS DONE TO DETERMINE THE SIDE index of the CONTACT SIDE
                      CheckNodeID = ElementConnectivities(J+(I-1), IEl) !in case that 1-2-3-1 is side order
                    else 
                      CheckNodeID = ElementConnectivities(3-2*(J-1), IEl)  !in case that 1-2-3-1 is side order
                    end if   
                if ( (CheckNodeID== ContSurfaceNodes (ISurf, 1)).or. &
                     (CheckNodeID== ContSurfaceNodes (ISurf, 2))) then
                   NFoundNodes = NFoundNodes + 1
                 end if
                 end do
                 if (NFoundNodes==2) then
                  ContctSurfaceElmID(ISurf) = IEl
                  Sides(ISurf) = I
                  EXIT
                 end if
               end do ! sides
              end if ! pile
            end do ! elements
          end do ! contact surfaces
       end if   
       
       allocate(ContactSurfaceElmNodes(NContSur, 4), stat=IError)
       ContactSurfaceElmNodes= 0
        
       NContactElm = NContSur
       
       do ISurf = 1, NContSur
          ContactSurfaceElmNodes (ISurf, 1)= &
                                    INT(ContctSurfaceElmID(ISurf)) ! element ID
            ContactSurfaceElmNodes (ISurf, 2)= Sides (ISurf)
          do J = 1, 2
            ContactSurfaceElmNodes (ISurf, J+2)= &
                                              ContSurfaceNodes(ISurf, J)
          end do
       end do
       
       
          if (allocated(TotSrfaceAndLoad) ) then
            deallocate(TotSrfaceAndLoad, stat = IError)
          end if
          
          if (allocated(Sides) ) then
            deallocate(Sides, stat = IError)
          end if
          
          if (allocated(IsContactSurface) ) then
            deallocate(IsContactSurface, stat = IError)
          end if
          
          !if (allocated(ContSurfaceNodes) ) then
            !deallocate(ContSurfaceNodes, stat = IError)
          !end if
          
          if (allocated(Element) ) then
            deallocate(Element, stat = IError)
          end if
          
          if (allocated(ContctSurfaceElmID) ) then
            deallocate(ContctSurfaceElmID, stat = IError)
          end if

            if (allocated(ContactNodeNormals)) then
              deallocate(ContactNodeNormals, stat = IError)
            end if
            if (allocated(InterfaceNodes)) then
              deallocate(InterfaceNodes, stat = IError)
            end if

            if (NContactElm==0) RETURN

            allocate(ContactNodeNormals(2, Counters%NodTot), stat = IError)
            allocate(NNodes(2*NContactElm), stat = IError)
            allocate(DVectors(2,2), stat = IError)
            allocate(NNormals(2*NContactElm,2), stat = IError)
            allocate(NNormals2(2,2), stat = IError)
            allocate(MatNodes(NContactElm), stat = IError)
            allocate(MatNodesCood(NContactElm,2), stat = IError)
            
            ContactNodeNormals = 0.0
            allocate(InterfaceNodes(Counters%NodTot), stat = IError)
            InterfaceNodes = .false.
       
               do I = 1, NContactElm
                do J = 1, 2 ! corner nodes
                  NodID = ContactSurfaceElmNodes (I, J+2)
                  InterfaceNodes(NodID) = .true.
                end do
               end do
            NSideNodes = GetNSideNodes(ELEMENTNODES)

            allocate(BoundaryElements(Counters%NEl), stat = IError)
            BoundaryElements = 0
              do I = 1, NContactElm
               ElementID = ContactSurfaceElmNodes (I, 1)
               SideID = ContactSurfaceElmNodes (I, 2)
               BoundaryElements(ElementID) = SideID
              end do 
            do I = 1, NContactElm
                SideID = ContactSurfaceElmNodes (I, 2)
                
             if (SideID == 1) then
                MatNodes(I)= ElementConnectivities(3, ContactSurfaceElmNodes (I, 1))
                MatNodesCood(I,1:2) = NodalCoordinates(MatNodes(I),1:2)
                do J = 1, 2
                   NNodes(J+2*(I-1)) = ElementConnectivities(J, ContactSurfaceElmNodes (I, 1))
                   NNormals(J+2*(I-1),1:2) = NodalCoordinates(NNodes(J+2*(I-1)),1:2)
                end do
             end if
             
             if (SideID == 2) then
                MatNodes(I)= ElementConnectivities(1, ContactSurfaceElmNodes (I, 1))
                MatNodesCood(I,1:2) = NodalCoordinates(MatNodes(I),1:2)
                do J = 1, 2
                   NNodes(J+2*(I-1)) = ElementConnectivities(J+1, ContactSurfaceElmNodes (I, 1))
                   NNormals(J+2*(I-1),1:2) = NodalCoordinates(NNodes(J+2*(I-1)),1:2)
                end do
             end if
             if (SideID == 3) then
                   MatNodes(I)= ElementConnectivities(2, ContactSurfaceElmNodes (I, 1))
                   MatNodesCood(I,1:2) = NodalCoordinates(MatNodes(I),1:2)
                   
                   NNodes(1+2*(I-1)) = ElementConnectivities(3, ContactSurfaceElmNodes (I, 1))
                   NNormals(1+2*(I-1),1:2) = NodalCoordinates(NNodes(1+2*(I-1)),1:2)
                   NNodes(2+2*(I-1)) = ElementConnectivities(1, ContactSurfaceElmNodes (I, 1))
                   NNormals(2+2*(I-1),1:2) = NodalCoordinates(NNodes(2+2*(I-1)),1:2)
             end if
                 
            end do
            dummy3 = 0
            do I = 1, 2*NContactElm
             dummy2 = 1E5
             JRef = 0
             NRef = 0
             if (NNodes(I) /= dummy3) then
              do J = 1, 2*NContactElm   
                dummy = NORM2(NNormals(I,:)-NNormals(J,:))
                    if (((dummy+1E-10) < dummy2).and.(dummy > 0.0)) then
                       dummy2 = dummy 
                       JRef = J
                       IRef = NNodes(JRef)  
                    else if ((abs(dummy-dummy2)<1E-10).and.(NNodes(J) /= IRef).and.(NNodes(J) /= NNodes(I)))  then    
                       KRef = J
                       NRef = NNodes(KRef) 
                    end if 
              end do 
              DVectors(1,1:2) = NNormals(JRef,1:2) - NNormals(I,1:2) 
              if (NRef /= 0) then
                DVectors(2,1:2) = NNormals(KRef,1:2) - NNormals(I,1:2)
              else 
                DVectors(2,1:2) = 0
              end if
              NNormals2(1,1) = -DVectors(1,2)
              NNormals2(1,2) =  DVectors(1,1)
              NNormals2(2,1) = -DVectors(2,2)
              NNormals2(2,2) =  DVectors(2,1)   
              if (dot_product(MatNodesCood(ceiling(I*0.5),1:2) - NodalCoordinates(NNodes(I),1:2),NNormals2(1,:)) > 0.0) then 
                NNormals2(1,1:2) = -1.0*NNormals2(1,1:2)    
              end if
              if (dot_product(MatNodesCood(ceiling(I*0.5),1:2) - NodalCoordinates(NNodes(I),1:2),NNormals2(2,:)) > 0.0) then 
                NNormals2(2,1:2) = -1.0*NNormals2(2,1:2)      
              end if
              ContactNodeNormals(1:2,NNodes(I)) =  NNormals2(1,1:2) + NNormals2(2,1:2) 
              dummy3 = NNodes(I)
             end if
            end do

            NInterfaceNodes = 0
            do I = 1, Counters%NodTot
              if (InterfaceNodes(I)) then
                NInterfaceNodes = NInterfaceNodes + 1
                if (CalParams%ApplyContactMeshBoundary) then
                  if ((ContactNodeNormals(1, I)>0.0).and.((NodalCoordinates(I, 2))>1E-5)) then
                    ContactNodeNormals(2, I) = 0.0
                    ContactNodeNormals(1:3, I) = VectorNorm(ContactNodeNormals(1:3, I), 3) ! generalise dimension
                  end if
                
                  ContactNodeNormals(1:3, I) = -1.0 * ContactNodeNormals(1:3, I)
                end if
              end if
              ContactNodeNormals(1:2, I) = VectorNorm(ContactNodeNormals(1:2, I), 2)
            end do

            deallocate(BoundaryElements, stat = IError)               
               
         end subroutine ContactSurfaceNodes_2D
    
          
      end module ModMPMDynContact