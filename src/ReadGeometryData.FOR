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
	  
	  
	  module ModReadGeometryData
      !**********************************************************************
      !
      ! Function: Contains routines for reading geometry data from GOM file
      !
      ! Note: Version archived subroutines are contained at the end of the module
      !
      !     $Revision: 10002 $
      !     $Date: 2023-06-19 12:44:04 +0200 (ma, 19 jun 2023) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModReadCalculationData
      use ModFileIO
    
      implicit none

        type GeometryParameterType
          real(REAL_TYPE), dimension(:), allocatable :: &
              ExcavatedElements, & ! Excavation Tool: vector to store the relation between element IDs and excavated volumes
              LocalDampingFactorElement ! local damping (defined per element)
          
          real(REAL_TYPE), dimension(:,:), allocatable :: &
              AbsorbingBoundariesSurfacesSolid, & ! Absorbing Boundaries: array storing absorbing conditions in coordinate-directions for solid
              AbsorbingBoundariesLinesSolid, &
              AbsorbingBoundariesPointsSolid, & 
              AbsorbingBoundariesSurfacesLiquid, & ! Absorbing Boundaries: array storing absorbing conditions in coordinate-directions for liquid
              AbsorbingBoundariesLinesLiquid, &
              AbsorbingBoundariesPointsLiquid, &
              AbsorbingBoundariesSurfacesGas, & ! Absorbing Boundaries: array storing absorbing conditions in coordinate-directions for gas
              AbsorbingBoundariesLinesGas, &
              AbsorbingBoundariesPointsGas, &
              ContactProperties, & !
              PrescribedVeloElValue, & ! (PrescribedVeloNElem,vsize) stores the values of prescribed velocity at elements
		      InitialVelocityonMP ! (Number of elments with initial velocities, vsize+1) store the element ID and the initial velocities of particles insde the element
          
          integer(INTEGER_TYPE), allocatable  :: &
              MovingMeshExtendingCorners(:), & ! dsize * (dsize - 1) + ExtraNodesMovingMesh  | (default = -1)
              MovingMeshCompressingCorners(:), & ! dsize * (dsize - 1) + ExtraNodesMovingMesh | (default = -1)
              MovingMeshMovingCorners(:), & ! dsize * (dsize - 1) + ExtraNodesMovingMesh  | (default = -1)   
              ContactLocation(:,:), &
              PrescribedVeloElemID(:),&!stores the ID of element on which velocity is prescribed
              PrescribedVeloElDirection(:,:), & !0/1 define the direction where velocity is prescribed
              WaterSurfaceMaterialID(:) ! Material Index on which a given water surface is assigned
                
          integer(INTEGER_TYPE) :: &    
              MovingMeshReferenceMaterialID, & ! the material_id of the reference material for the moving mesh
              NABSurfaceNodesSolid, & ! total number of absorbing nodes belonging to surfaces, lines or points for solid
              NABLineNodesSolid, &
              NABPointNodesSolid, &
              NABSurfaceNodesLiquid, & ! total number of absorbing nodes belonging to surfaces, lines or points for liquid
              NABLineNodesLiquid, &
              NABPointNodesLiquid, &
              NABSurfaceNodesGas, & ! total number of absorbing nodes belonging to surfaces, lines or points for gas
              NABLineNodesGas, &
              NABPointNodesGas, &
              AbsorbingBoundaryReferenceMaterialID = 0, &  ! the material_id of the reference material for the absorbing boundaries
              PrescribedVeloNElem = 0, & !Number of element on which velocity is prescribed at MP
              NMovingElements, & !Number of elements with initial velocity condition
		      ExtraNodesMovingMesh = 2, & !Number of extra nodes for moving mesh areas
              NumberWaterSurfaceMaterials

          character(len=255) :: &
              MovingMeshDirection
          
          character(len=255), allocatable :: &
              ContactMaterials(:,:)  , &
              WaterSurfaceFileNumber(:)  ! Number of .PSF file including initial water surface

          logical :: &
              ApplyAbsorbingBoundary, & ! absorbing boundaries
              ApplyLocalDampingElement, & ! local damping (defined per element)
              ApplyNodalVelocityVolume, & !.true. if prescribed velcity is assigned to volumes
              ApplyNodalVelocitySurface, & !.true. if prescribed velcity is assigned to surfaces
              ApplyNodalVelocityLine, & !.true. if prescribed velcity is assigned to lines
              ApplyNodalVelocityPoint, & !.true. if prescribed velcity is assigned to nodes
              ApplyMPVelocityVolume, & !.true. if prescribed velcity is assigned to volumes
              ApplyMPVelocitySurface, & !.true. if prescribed velcity is assigned to surfaces
              ApplyInitialWaterSurfaceFromFile
              
              
        end type GeometryParameterType

        type(GeometryParameterType), public, save :: GeoParams ! stores geometry parameters

      
    contains ! subroutines of this module
     
 
        subroutine InitializeGeometryParameters()
        !**********************************************************************
        !
        ! Function: Initializes geometry parameters (GeoParams%) related to moving mesh, absorving boundaries, local damping
        !
        !**********************************************************************
            implicit none
			if ( IS3DCYLINDRIC ) GeoParams%ExtraNodesMovingMesh = 0 ! For 3D Cylindrical coordinates there are only 6 nodes for moving mesh 
            allocate( &
              GeoParams%MovingMeshExtendingCorners(NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh), &
              GeoParams%MovingMeshCompressingCorners(NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh), &
              GeoParams%MovingMeshMovingCorners(NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh), &
              GeoParams%LocalDampingFactorElement(Counters%NEl) )
              
              GeoParams%MovingMeshExtendingCorners = -1
              GeoParams%MovingMeshCompressingCorners = -1
              GeoParams%MovingMeshMovingCorners = -1
              
              GeoParams%MovingMeshReferenceMaterialID = 0
              GeoParams%MovingMeshDirection = 'undefined'
              GeoParams%ApplyAbsorbingBoundary = .false. ! absorbing boundaries
              GeoParams%ApplyLocalDampingElement = .false. ! local damping (defined per element)
              
              GeoParams%LocalDampingFactorElement = 0.0 ! local damping (defined per element)
              
        end subroutine InitializeGeometryParameters
    
    
        subroutine ReadGeometryParameters()
        !**********************************************************************
        !
        ! Function: Determines GOM file version and calls respective GOM reader
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables       
          character(len=MAX_FILENAME_LENGTH) :: FileName, FileVersion
          integer(INTEGER_TYPE) :: FileUnit
          character(len=255) :: BName
          integer(INTEGER_TYPE) :: ios

          ! Initialize GeoParams
          call InitializeGeometryParameters()
          
          FileName = trim(CalParams%FileNames%ProjectName)//GOM_FILE_EXTENSION
          FileUnit = TMP_UNIT
          
          ! check if GOM file exists in project folder, otherwise give error and stop execution
          if ( FExist(trim(FileName)) ) then
            call GiveMessage('Reading GOM file (Geometry): ' // trim(FileName) )  
          else
            call GiveError('GOM file does not exist!' // NEW_LINE('A') // 'required GOM file: ' // trim(FileName) )
          end if
          
          ! open GOM file
          call FileOpen(FileUnit, trim(FileName))
          
          ! determine current version of GOM file 
          read(FileUnit, '(A)', iostat=ios) BName 
          call Assert( ios == 0, 'GOM file: Can''t read flag from GOM file.' )
          FileVersion = trim(BName)
          
          ! read GOM data              
          select case (FileVersion) ! read GOM data depending on file version
              
            case (Anura3D_v2023)             
              call ReadGOM(FileUnit,FileVersion)
            case (Anura3D_v2022)             
              call ReadGOM(FileUnit,FileVersion)
            case (Anura3D_v2021)             
              call ReadGOM(FileUnit,FileVersion)

            case (Anura3D_v2019_2)             
              call ReadGOM(FileUnit,FileVersion)
            case default
              call GiveError('Wrong version of GOM file!' // NEW_LINE('A') // 'Supported CPS versions: ' &
                // trim(Anura3D_v2019_2) // ', ' &
                // trim(Anura3D_v2021) // ', ' &
                // trim(Anura3D_v2022) // ', ' &
                // trim(Anura3D_v2023)     // '.' )

          end select
            
          ! close GOM file    
          close(FileUnit)
      
        end subroutine ReadGeometryParameters

        
        subroutine InitialiseDimension()
        !**********************************************************************
        !
        !    Function:  Contains code for initialising dimension and element type
        !
        !**********************************************************************
        implicit none
      
          ! Local variables
          character(len = 255) :: FileName, TName, Bname, ProjectName
          integer(INTEGER_TYPE) :: ios ! used for error control
          character(len=255) :: DumS = '' ! reading strings
          character(len=21) :: messageIOS = 'GOM file: Can''t read '
          
          integer(INTEGER_TYPE) :: NDimension = -1
          logical :: Axisymmetric = .false.
          logical :: Cylindric = .false.
          character(len=255) :: ReadElementType = ''
          character(len=255) :: Formulation = ''
          
          ! open GOM file
          call getarg(1, ProjectName)
          FileName=Trim(ProjectName)//'.GOM'
          if (FExist(FileName)) open(GOMunit, FILE=FileName)
          
          do ! read dimension and elementtype from gom-file
              
            read(GOMunit,'(A)') TName
            BName = TName
            
            if (trim(BName) == '$$DIMENSION') then
              read(GOMunit, *, iostat=ios) DumS
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == DIM_2D_PLANESTRAIN .or. DumS == DIM_2D_AXISYMMETRIC .or. DumS == DIM_3D_CARTESIAN .or. DumS == DIM_3D_CYLINDRIC, 'GOM file: ' //trim(BName)// ' must be equal to ' &
                           // DIM_2D_PLANESTRAIN // ' or ' // DIM_2D_AXISYMMETRIC // ' or ' // DIM_3D_CARTESIAN // ' or ' // DIM_3D_CYLINDRIC // '.' )
              if ( DumS == DIM_2D_PLANESTRAIN .or. DumS == DIM_2D_AXISYMMETRIC ) NDimension = 2 ! set dimension for 2D plane strain and 2D axisymmetric analysis
              if ( DumS == DIM_3D_CARTESIAN .or. DumS == DIM_3D_CYLINDRIC ) NDimension = 3 ! set dimension for full 3D analysis (cartesian or cylindric)
              if ( DumS == DIM_2D_AXISYMMETRIC ) Axisymmetric = .true. ! set flag for 2D axisymmetric analysis
              if ( DumS == DIM_3D_CYLINDRIC ) Cylindric = .true. ! set flag for 3D analysis with cylindrical coordinates
              
            else if (trim(BName) == '$$ELEMENTTYPE') then
              read(GOMunit, *, iostat=ios) DumS
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == TRI3 .or. DumS == TRI6 .or. DumS == QUAD4 .or. DumS == QUAD8 .or. &
                           DumS == TETRA4 .or. DumS == TETRA10 .or. DumS == HEXA8 .or. DumS == HEXA20 .or. DumS == TETRAOLD, 'GOM file: ' //trim(BName)// &
                           ' must be equal to ' // TRI3 // ' or ' // TRI6 // ' or ' // QUAD4 // ' or ' // QUAD8 // ' or ' // &
                           TETRA4 // ' or ' // TETRA10 // ' or ' // HEXA8 // ' or ' // HEXA20  // ' or ' // TETRAOLD  // '.' )
              ReadElementType = DumS

            else if (trim(BName) == '$$FORMULATION') then
              read(GOMunit, *, iostat=ios) DumS
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == SINGLE_POINT .or. DumS == DOUBLE_POINT, 'GOM file: ' //trim(BName)// ' must be equal to ' & 
                           // SINGLE_POINT // ' or ' // DOUBLE_POINT //'.' )
              Formulation = DumS
              
            else if (trim(BName) == '$$FINISH') then
              close(GOMunit)  
              EXIT
            end if 
              
          end do

          ! old versions of the GOM file (before v2018.2) do not specify the dimension and elementtype, therefore it is set manually              
          if (NDimension == -1) NDimension = 3 ! this is the only available dimension before v2018.2
          if (ReadElementType == '') ReadElementType = TETRAOLD ! the only elementtype that was used before v2018.2
          if (Formulation == '') Formulation = SINGLE_POINT ! if nothing assumed single-point
          
          call SetDimension(NDimension, Axisymmetric, Cylindric)
          call SetElementType(ReadElementType)
          call SetFormulation(Formulation)
          
        end subroutine InitialiseDimension
        
        
        subroutine ReadGOM(FileUnit,FileVersion)
        !**********************************************************************
        !
        ! Function : Reads input variables from the GOM file
        !
        !            GOM version 2022, 2021, 2019.2
        !           
        !
        !**********************************************************************
        implicit none

          character(len=MAX_FILENAME_LENGTH), intent(in) :: FileVersion
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! local variables
          integer(INTEGER_TYPE) :: I, J, SizeAB
          integer(INTEGER_TYPE) :: ios ! used for error control
          integer(INTEGER_TYPE) :: DumI(2)
          real(REAL_TYPE) :: DumR(3)
          character(len=255) :: DumS
          character(len=21) :: messageIOS = 'GOM file: Can''t read '
          character(len=255) :: BName

          
          SizeAB = 1 + 3 * NDIM ! array size for reading absorbing boundary data
		            
          ! set GOM version number
            select case (FileVersion)
            case (Anura3D_v2023)
                call GiveMessage('Reading... ' // Anura3D_v2023)
                CalParams%GOMversion = Anura3D_v2021
            case (Anura3D_v2022)
                call GiveMessage('Reading... ' // Anura3D_v2022)
                CalParams%GOMversion = Anura3D_v2021
            case (Anura3D_v2021)
                call GiveMessage('Reading... ' // Anura3D_v2021)
                CalParams%GOMversion = Anura3D_v2021
                case (Anura3D_v2019_2)
                call GiveMessage('Reading... ' // Anura3D_v2021)
                CalParams%GOMversion = Anura3D_v2021
            end select

          do 
            read(FileUnit, '(A)') BName
            
              !!! DATA FOR CONTACT ALGORITHM
              if (trim(BName)=='$$START_BODY_CONTACT_2D') then ! 2D body contact data exists  
                read(FileUnit,*,iostat=ios) DumI(1)  ! number of body contact elements
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: ' //trim(BName)// ' must be larger than 0.' )
                call InitialiseContact(DumI(1))  
                do I = 1, DumI(1)
                  read(FileUnit,*) GeoParams%ContactLocation(I,1:2), &
                    GeoParams%ContactMaterials(I,1), GeoParams%ContactProperties(I,1:2), &
                    GeoParams%ContactMaterials(I,2), GeoParams%ContactProperties(I,3:4), &
                    GeoParams%ContactMaterials(I,3), GeoParams%ContactProperties(I,5:6), & 
                    GeoParams%ContactMaterials(I,4), GeoParams%ContactProperties(I,7:8)
                end do 
              
              else if (trim(BName)=='$$START_BOUNDARY_CONTACT_2D') then ! 2D boundary contact data exists  
                read(FileUnit,*,iostat=ios) DumI(1)  ! number of boundary contact elements
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: ' //trim(BName)// ' must be larger than 0.' )
                call InitialiseContact(DumI(1))  
                do I = 1, DumI(1)
                  read(FileUnit,*) GeoParams%ContactLocation(I,1:2), &
                    GeoParams%ContactMaterials(I,1), GeoParams%ContactProperties(I,1:2), &
                    GeoParams%ContactMaterials(I,2), GeoParams%ContactProperties(I,3:4), &
                    GeoParams%ContactMaterials(I,3), GeoParams%ContactProperties(I,5:6), & 
                    GeoParams%ContactMaterials(I,4), GeoParams%ContactProperties(I,7:8)
                end do
            

                
              !!! DATA FOR PRESCRIBED VELOCITY
              else if (trim(BName)=='$$PRESCRIBED_MATERIAL_POINT_VELOCITY_VOLUME') then !prescribed velcoity is assigned
                GeoParams%ApplyMPVelocityVolume=.true.
                call ReadPrescribedVelocityMaterialPoint(FileUnit)
         
              else if (trim(BName)=='$$PRESCRIBED_MATERIAL_POINT_VELOCITY_SURFACE') then !prescribed velcoity is assigned
                GeoParams%ApplyMPVelocitySurface=.true.
                call ReadPrescribedVelocityMaterialPoint(FileUnit)
                
              else if (trim(BName)=='$$PRESCRIBED_NODAL_VELOCITY_VOLUME') then !prescribed velcoity is assigned
                GeoParams%ApplyNodalVelocityVolume=.true.
                call ReadPrescribedVelocityNode(FileUnit)
                
              else if (trim(BName)=='$$PRESCRIBED_NODAL_VELOCITY_SURFACE') then !prescribed velcoity is assigned
                GeoParams%ApplyNodalVelocitySurface=.true.
                call ReadPrescribedVelocityNode(FileUnit)
                
              else if (trim(BName)=='$$PRESCRIBED_NODAL_VELOCITY_LINE') then !prescribed velcoity is assigned
                GeoParams%ApplyNodalVelocityLine=.true.
                call ReadPrescribedVelocityNode(FileUnit)
                
              else if (trim(BName)=='$$PRESCRIBED_NODAL_VELOCITY_POINT') then !prescribed velcoity is assigned
                GeoParams%ApplyNodalVelocityPoint=.true.
                call ReadPrescribedVelocityNode(FileUnit)
				
		      !!! DATA FOR INITIAL VELOCITY ON MP
			  else if (trim(BName)=='$$INITIAL_VELOCITY_MATERIAL_POINT') then ! Initial velocity on MP data
                  read(FileUnit,*,iostat=ios) DumI(1)
				  call Assert( ios == 0, messageIOS // trim(BName) )
				  call Assert( DumI(1) > 0, 'GOM-ERROR: The number of elements to initialise velocity on material points must be greater than 0!' )
				  GeoParams%NMovingElements=DumI(1)
				  allocate(GeoParams%InitialVelocityonMP(DumI(1),NVECTOR+1), stat = ios)
				  call Assert( ios == 0, messageIOS // trim(BName) )				  
				  do I = 1, GeoParams%NMovingElements
					read(FileUnit,*,iostat=ios) (GeoParams%InitialVelocityonMP(I, J), J= 1, NVECTOR+1)
					call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )	
                  end do                 
                
              !!! DATA FOR EXTERNAL INITIAL WATER SURFACE
              else if (trim(BName)=='$$INITIAL_WATER_SURFACE_FROM_FILE') then
                  GeoParams%ApplyInitialWaterSurfaceFromFile = .true.
                  read(FileUnit,*,iostat=ios) DumI(1) ! Number of materials on which an initial water surface is assigned
                  GeoParams%NumberWaterSurfaceMaterials = DumI(1)
                  call InitialiseWaterSurfaceFromFileArrays(GeoParams%NumberWaterSurfaceMaterials)  
                  do I = 1, GeoParams%NumberWaterSurfaceMaterials
                  read(FileUnit,*,iostat=ios) GeoParams%WaterSurfaceMaterialID(I) , GeoParams%WaterSurfaceFileNumber(I)                 
                  end do                  
                                     
              !!! DATA FOR MOVING MESH  
              else if (trim(BName)=='$$EXTENDING_MESH_CORNER_NODES') then ! moving mesh data
                do I = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                  read(FileUnit,*,iostat=ios) DumI(1)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                  call Assert( DumI(1) > 0, 'GOM-ERROR: all node numbers of $$EXTENDING_MESH_CORNER_NODES must be positive!' )
                  GeoParams%MovingMeshExtendingCorners(I) = DumI(1)
                end do   
              else if (trim(BName)=='$$COMPRESSING_MESH_CORNER_NODES') then ! moving mesh data
                do I = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                  read(FileUnit,*,iostat=ios) DumI(1)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                  call Assert( DumI(1) > 0, 'GOM-ERROR: all node numbers of $$COMPRESSING_MESH_CORNER_NODES must be positive!' )
                  GeoParams%MovingMeshCompressingCorners(I) = DumI(1)
                end do   
              else if (trim(BName)=='$$MOVING_MESH_CORNER_NODES') then ! moving mesh data
                do I = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                  read(FileUnit,*,iostat=ios) DumI(1)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                  call Assert( DumI(1) > 0, 'GOM-ERROR: all node numbers of $$MOVING_MESH_CORNER_NODES must be positive!' )
                  GeoParams%MovingMeshMovingCorners(I) = DumI(1)
                end do   
              else if (trim(BName)=='$$MOVING_MESH_DIRECTION') then ! moving mesh data
                read(FileUnit,*,iostat=ios) DumS
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( trim(DumS) == 'x-direction' .or. trim(DumS) == 'y-direction' .or. trim(DumS) == 'z-direction', &
                             'GOM-ERROR: $$MOVING_MESH_DIRECTION must be "x-direction", "y-direction" or "z-direction" (only for 3D)!' )
                GeoParams%MovingMeshDirection = DumS
              else if (trim(BName)=='$$MOVING_MESH_REFERENCE_MATERIAL_INDEX') then 
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%MovingMeshReferenceMaterialID = DumI(1)

              !!! DATA FOR ABSORBING BOUNDARY    
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_SURFACE_SOLID') then
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of surfaces in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABSurfaceNodesSolid = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesSurfacesSolid(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesSurfacesSolid(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do 
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_LINE_SOLID') then
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of lines in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABLineNodesSolid = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesLinesSolid(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesLinesSolid(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do 
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_POINT_SOLID') then   
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of points in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABPointNodesSolid = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesPointsSolid(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesPointsSolid(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do 
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_SURFACE_LIQUID') then 
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of surfaces in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABSurfaceNodesLiquid = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesSurfacesLiquid(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesSurfacesLiquid(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do  
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_LINE_LIQUID') then
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of lines in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABLineNodesLiquid = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesLinesLiquid(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesLinesLiquid(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do 
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_POINT_LIQUID') then 
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of points in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABPointNodesLiquid = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesPointsLiquid(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesPointsLiquid(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do   
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_SURFACE_GAS') then 
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of surfaces in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABSurfaceNodesGas = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesSurfacesGas(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesSurfacesGas(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do  
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_LINE_GAS') then
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of lines in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABLineNodesGas = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesLinesGas(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesLinesGas(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do 
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_POINT_GAS') then 
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: number of points in ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%NABPointNodesGas = DumI(1)
                GeoParams%ApplyAbsorbingBoundary = .true.
                allocate(GeoParams%AbsorbingBoundariesPointsGas(DumI(1),SizeAB), stat = ios)
                call Assert( ios == 0, messageIOS // trim(BName) )
                do I = 1, DumI(1)
                  read(FileUnit,*,iostat=ios) (GeoParams%AbsorbingBoundariesPointsGas(I, J), J = 1, SizeAB)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                end do  
              else if (trim(BName)=='$$ABSORBING_BOUNDARY_REFERENCE_MATERIAL_INDEX') then 
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%AbsorbingBoundaryReferenceMaterialID = DumI(1)
                
              !!! DATA FOR LOCAL DAMPING ELEMENT   
              else if (trim(BName)=='$$STARTDAMPING') then
                do I = 1, Counters%NEl ! loop over elements
                  read(FileUnit,*,iostat=ios) DumR(1)
                  call Assert( ios == 0, messageIOS // trim(BName) )
                  if (DumR(1) > 0.0) GeoParams%ApplyLocalDampingElement = .true. 
                  GeoParams%LocalDampingFactorElement(I) = DumR(1) ! read the damping value for each element
                end do
            
             !!! DATA FOR EXCAVATION
              else if (trim(BName)=='$$START_EXCAVATION_SOLID') then ! excavation data exists  
                call InitialiseExcavationData()  
                do 
                  read(FileUnit,*,iostat=ios) DumI(1), DumI(2)
                  if (ios /= 0) then
                      backspace(FileUnit)
                      exit
                 end if
                  GeoParams%ExcavatedElements(DumI(2)) = DumI(1)
                end do    
              !!! END OF GOM FILE  
              else if (trim(BName)=='$$FINISH') then
                EXIT
                
            end if 
          end do 
          
          ! Check input data for correctness and/or completeness
          call CheckMovingMeshData()
          call CheckAbsorbingBoundaryData()
          call CheckPrescribedVelocityData()
        
        end subroutine ReadGOM
        
        
        subroutine InitialiseExcavationData()
        !**********************************************************************
        !
        ! Function:  Contains code for initialising excavation data
        !
        !**********************************************************************
        implicit none
        
          call DestroyExcavationData()
          
          call InitialiseExcavationArray()
      
        end subroutine InitialiseExcavationData


        subroutine DestroyExcavationData()
        !**********************************************************************
        !
        ! Function:  Deallocates the arrays used in this module
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(GeoParams%ExcavatedElements)) then
            deallocate(GeoParams%ExcavatedElements, stat = IError)
          end if

        end subroutine DestroyExcavationData


        subroutine InitialiseExcavationArray()
        !**********************************************************************
        !
        ! Function:  To initialise the arrays relate to excavation
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (CalParams%ApplyExcavation) then
            allocate(GeoParams%ExcavatedElements(Counters%NEl), stat = IError)
          else
            allocate(GeoParams%ExcavatedElements(1), stat = IError)
          end if

          GeoParams%ExcavatedElements = 0.0

        end subroutine InitialiseExcavationArray
        
  
        subroutine InitialiseContact(ArraySize)
        !**********************************************************************
        !
        ! Function:  Contains code for initialising contact data arrays of size ArraySize
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ArraySize
        
          call DeallocateContactArrays()
          
          call AllocateAndInitialiseContactArrays(ArraySize)
      
        end subroutine InitialiseContact
        

        subroutine DeallocateContactArrays()
        !**********************************************************************
        !
        ! Function:  Deallocates the arrays used for contact
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if(allocated(GeoParams%ContactLocation)) deallocate(GeoParams%ContactLocation, stat = IError)
          if(allocated(GeoParams%ContactProperties)) deallocate(GeoParams%ContactProperties, stat = IError)
          if(allocated(GeoParams%ContactMaterials)) deallocate(GeoParams%ContactMaterials, stat = IError)

        end subroutine DeallocateContactArrays 
        
        
        subroutine AllocateAndInitialiseContactArrays(ArraySize)
        !**********************************************************************
        !
        ! Function:  To allocate and initialise the contact arrays with size ArraySize
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ArraySize      
          ! Local variables
          integer(INTEGER_TYPE) :: IError
          
          allocate(GeoParams%ContactLocation(ArraySize,2), stat=IError)
          allocate(GeoParams%ContactProperties(ArraySize,8), stat=IError)
          allocate(GeoParams%ContactMaterials(ArraySize,4), stat=IError)

          GeoParams%ContactLocation = 0
          GeoParams%ContactProperties = 0.0
          GeoParams%ContactMaterials = 'undefined'
                
        end subroutine AllocateAndInitialiseContactArrays 
        
        subroutine InitialiseWaterSurfaceFromFileArrays(ArraySize)
        !**********************************************************************
        !
        ! Function:  To initialise the arrays relate to Initial Water Surface from File
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ArraySize
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          allocate(GeoParams%WaterSurfaceMaterialID(ArraySize), stat=IError)
          allocate(GeoParams%WaterSurfaceFileNumber(ArraySize), stat=IError)
          
          call AllocateAndInitialiseContactArrays(ArraySize)
      
        end subroutine InitialiseWaterSurfaceFromFileArrays
        
        subroutine DestroyPrescribedMPVeloData()
        !**********************************************************************
        !
        ! Function:  Deallocates the arrays used in this module
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(GeoParams%PrescribedVeloElemID)) then
            deallocate(GeoParams%PrescribedVeloElemID, stat = IError)
          end if
          
         if (allocated(GeoParams%PrescribedVeloElValue)) then
            deallocate(GeoParams%PrescribedVeloElValue, stat = IError)
         end if
                    
         if (allocated(GeoParams%PrescribedVeloElDirection)) then
            deallocate(GeoParams%PrescribedVeloElDirection, stat = IError)
         end if

        end subroutine DestroyPrescribedMPVeloData


        subroutine InitialisePrescribedMPVeloData(Nnode)
        !**********************************************************************
        !
        ! Function:  To initialise the arrays relate to prescribed velocity
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError, Nnode

          allocate( GeoParams%PrescribedVeloElemID(Nnode), stat = IError)
             GeoParams%PrescribedVeloElemID = -1
          allocate( GeoParams%PrescribedVeloElValue(Nnode, NVECTOR), stat = IError)
             GeoParams%PrescribedVeloElValue = 0.0
          allocate( GeoParams%PrescribedVeloElDirection(Nnode, NVECTOR), stat = IError)
             GeoParams%PrescribedVeloElDirection = 1


        end subroutine InitialisePrescribedMPVeloData
        
        
        subroutine DestroyPrescribedNodalVeloData()
        !**********************************************************************
        !
        ! Function:  Deallocates the arrays used in this module
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(CalParams%PrescribedVelo%NodePrescribedVelo)) then
            deallocate(CalParams%PrescribedVelo%NodePrescribedVelo, stat = IError)
          end if
          if (allocated(CalParams%PrescribedVelo%NodalPrescribedVelocityValue)) then
            deallocate(CalParams%PrescribedVelo%NodalPrescribedVelocityValue, stat = IError)
          end if
          if (allocated(CalParams%PrescribedVelo%NodalPrescribedVelocityDirection)) then
            deallocate(CalParams%PrescribedVelo%NodalPrescribedVelocityDirection, stat = IError)
          end if

        end subroutine DestroyPrescribedNodalVeloData


        subroutine InitialisePrescribedNodalVeloData(Nnode)
        !**********************************************************************
        !
        ! Function:  To initialise the arrays relate to prescribed velocity
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError, Nnode

          allocate( CalParams%PrescribedVelo%NodePrescribedVelo(Nnode), stat = IError)
          CalParams%PrescribedVelo%NodePrescribedVelo = -1
          allocate( CalParams%PrescribedVelo%NodalPrescribedVelocityValue(Nnode, NVECTOR), stat = IError)
          CalParams%PrescribedVelo%NodalPrescribedVelocityValue = 0.0
          allocate( CalParams%PrescribedVelo%NodalPrescribedVelocityDirection(Nnode, NVECTOR), stat = IError)
          CalParams%PrescribedVelo%NodalPrescribedVelocityDirection = 1


        end subroutine InitialisePrescribedNodalVeloData
        
        
        subroutine ReadPrescribedVelocityMaterialPoint(FileUnit)
        !**********************************************************************
        !
        ! Function:  To read data for prescribed velocity in material points
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
        
          ! local variables
          integer(INTEGER_TYPE) :: I
          integer(INTEGER_TYPE), dimension(NVECTOR+1) :: DumI
          real(REAL_TYPE), dimension(NVECTOR) :: DumR
        
          call DestroyPrescribedMPVeloData()
                
          read(FileUnit,*) DumI(1)
          GeoParams%PrescribedVeloNElem = DumI(1)
                
          call InitialisePrescribedMPVeloData(GeoParams%PrescribedVeloNElem)
                                
          do I = 1, GeoParams%PrescribedVeloNElem
            read(FileUnit,*) DumI(1), DumI(2:NVECTOR+1), DumR(1:NVECTOR)
            GeoParams%PrescribedVeloElemID(I) = DumI(1)
            GeoParams%PrescribedVeloElDirection(I, 1:NVECTOR) = DumI(2:NVECTOR+1)
            GeoParams%PrescribedVeloElValue(I, 1:NVECTOR) = DumR(1:NVECTOR)
          end do 

        end subroutine ReadPrescribedVelocityMaterialPoint
        
        
        subroutine ReadPrescribedVelocityNode(FileUnit)
        !**********************************************************************
        !
        ! Function:  To read data for prescribed velocity in nodes
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
        
          ! local variables
          integer(INTEGER_TYPE) :: I
          integer(INTEGER_TYPE), dimension(NVECTOR+1) :: DumI
          real(REAL_TYPE), dimension(NVECTOR) :: DumR
                
          call DestroyPrescribedNodalVeloData()   
          
          read(FileUnit,*) DumI(1) 
          CalParams%PrescribedVelo%NNodePrescribedVelo = DumI(1)
                
          call InitialisePrescribedNodalVeloData(CalParams%PrescribedVelo%NNodePrescribedVelo)  
                
          do I = 1, CalParams%PrescribedVelo%NNodePrescribedVelo
            read(FileUnit,*) DumI(1), DumI(2:NVECTOR+1), DumR(1:NVECTOR)
            CalParams%PrescribedVelo%NodePrescribedVelo(I) = DumI(1)
            CalParams%PrescribedVelo%NodalPrescribedVelocityDirection(I, 1:NVECTOR) = DumI(2:NVECTOR+1)
            CalParams%PrescribedVelo%NodalPrescribedVelocityValue(I, 1:NVECTOR) = DumR(1:NVECTOR)
          end do  
          
        end subroutine ReadPrescribedVelocityNode  
                
        
        subroutine CheckMovingMeshData()
        !**********************************************************************
        !
        ! Function:  To check the input data for the moving mesh algorithm
        !
        !**********************************************************************
        
        use ModMeshInfo
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: &
              J, K, & 
              direction, &
              FixedNodeExtending, FixedNodeCompressing, &
              PositionNodeExtending, PositionNodeCompressing
          logical, dimension(NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh) :: IsMovingAndCompressingCorner, IsMovingAndExtendingCorner
        
          ! set flag for applying moving mesh algorithm         
          CalParams%ApplyMeshSmoothing = .false. 
          if ( GeoParams%MovingMeshExtendingCorners(1) /= -1 .or. &
               GeoParams%MovingMeshCompressingCorners(1) /= -1 .or. &
               GeoParams%MovingMeshMovingCorners(1) /= -1) CalParams%ApplyMeshSmoothing = .true.

          if ( .not. CalParams%ApplyMeshSmoothing ) RETURN
          
          ! set basic parameters for moving mesh algorithm
          CalParams%MovingMesh%NMovingMeshDirections = 1
          CalParams%MovingMesh%NAreaNodes = NDIM * (NDIM - 1) +GeoParams%ExtraNodesMovingMesh ! 4 for 2D (square), 8 for 3D (cube), 6 for cylindrical 3D

          ! check if extending and/or compressing mesh areas are defined and determine 'NStorageAreas'
          if ( CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_FILE ) then 
              ! no deforming mesh is defined for prescribed velocity from file
              CalParams%MovingMesh%NStorageAreas = 0
          else if ( GeoParams%MovingMeshExtendingCorners(1) == -1 .and. GeoParams%MovingMeshCompressingCorners(1) == -1) then 
            ! error: no deforming mesh is defined  
            call GiveError('Error: Moving Mesh: Unless you are prescribing velocity from file, at least one deforming mesh area has to be specified, i.e. either extending or compressing mesh.')  
          end if  
          if ( GeoParams%MovingMeshExtendingCorners(1) /= -1 .and. GeoParams%MovingMeshCompressingCorners(1) == -1 ) then
            ! only compressing mesh defined  
            CalParams%MovingMesh%NStorageAreas = 1
          end if  
          if ( GeoParams%MovingMeshExtendingCorners(1) == -1 .and. GeoParams%MovingMeshCompressingCorners(1) /= -1 ) then
            ! only extending mesh defined  
            CalParams%MovingMesh%NStorageAreas = 1
          end if  
          if ( GeoParams%MovingMeshExtendingCorners(1) /= -1 .and. GeoParams%MovingMeshCompressingCorners(1) /= -1 ) then
            ! both, extending and compressing mesh defined  
            CalParams%MovingMesh%NStorageAreas = 2
          end if  
          
          ! check if reference material is defined
          call Assert( GeoParams%MovingMeshReferenceMaterialID > 0 , 'Error: Moving Mesh: The reference material has to be defined.' )
          CalParams%MovingMesh%MovingMaterialID = GeoParams%MovingMeshReferenceMaterialID
          if (CalParams%ApplyContactAlgorithm) then 
			CalParams%MovingMesh%StructureMaterialID=CalParams%MovingMesh%MovingMaterialID !Sets StructureMaterialID required for contact and rigid body Algorithm
		  end if
          ! determine moving mesh direction
          if ( GeoParams%MovingMeshDirection == 'x-direction' ) direction = 1
          if ( GeoParams%MovingMeshDirection == 'y-direction' ) direction = 2
          if ( GeoParams%MovingMeshDirection == 'z-direction' ) direction = 3

          CalParams%MovingMesh%MovingMeshDirection = direction
          
          if ( GeoParams%MovingMeshCompressingCorners(1) /= -1) then
            IsMovingAndCompressingCorner = .false.
            do J = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
              do K = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                if (GeoParams%MovingMeshCompressingCorners(J)==GeoParams%MovingMeshMovingCorners(K)) then
                  IsMovingAndCompressingCorner(J) = .true.
                end if
              end do
            end do
          
            do J = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
              if (.not.IsMovingAndCompressingCorner(J)) then
                FixedNodeCompressing = GeoParams%MovingMeshCompressingCorners(J)
                PositionNodeCompressing = J
                EXIT
              end if
            end do 
          end if
          
          if ( GeoParams%MovingMeshExtendingCorners(1) /= -1) then
            IsMovingAndExtendingCorner = .false.
            do J = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
              do K = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                if (GeoParams%MovingMeshExtendingCorners(J)==GeoParams%MovingMeshMovingCorners(K)) then
                  IsMovingAndExtendingCorner(J) = .true.
                end if
              end do
            end do
          
            do J = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
              if (.not.IsMovingAndExtendingCorner(J)) then
                FixedNodeExtending = GeoParams%MovingMeshExtendingCorners(J)
                PositionNodeExtending = J
                EXIT
              end if
            end do 
          end if
          
          ! assign corner node numbers to array CalParams%MovingMesh%MeshAreas(k,j,i)
          K = CalParams%MovingMesh%NMovingMeshDirections
          
          if (CalParams%MovingMesh%NStorageAreas(K) == 1) then

            ! extending mesh
            if ( GeoParams%MovingMeshExtendingCorners(1) /= -1) then
              CalParams%MovingMesh%MeshAreas(K, 1, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh) = GeoParams%MovingMeshExtendingCorners(1:NDIM * (NDIM - 1) + 2) ! fill whole array
              CalParams%MovingMesh%MeshAreas(K, 1, 1) = GeoParams%MovingMeshExtendingCorners(PositionNodeExtending) ! replace first node with the one on the fixed boundary
              CalParams%MovingMesh%MeshAreas(K, 1, PositionNodeExtending) = GeoParams%MovingMeshExtendingCorners(1) ! put the old first node on the empty position 
            end if
            ! compressing mesh
            if ( GeoParams%MovingMeshCompressingCorners(1) /= -1) then 
              CalParams%MovingMesh%MeshAreas(K, 1, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh) = GeoParams%MovingMeshCompressingCorners(1:NDIM * (NDIM - 1) + 2) ! fill whole array
              CalParams%MovingMesh%MeshAreas(K, 1, 1) = GeoParams%MovingMeshCompressingCorners(PositionNodeCompressing) ! replace first node with the one on the fixed boundary
              CalParams%MovingMesh%MeshAreas(K, 1, PositionNodeCompressing) = GeoParams%MovingMeshCompressingCorners(1) ! put the old first node on the empty position 
            end if

          else if (CalParams%MovingMesh%NStorageAreas(K) == 2) then

            ! extending mesh
            CalParams%MovingMesh%MeshAreas(K, 1, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh) = GeoParams%MovingMeshExtendingCorners(1:NDIM * (NDIM - 1) + 2) ! fill whole array
            CalParams%MovingMesh%MeshAreas(K, 1, 1) = GeoParams%MovingMeshExtendingCorners(PositionNodeExtending) ! replace first node with the one on the fixed boundary
            CalParams%MovingMesh%MeshAreas(K, 1, PositionNodeExtending) = GeoParams%MovingMeshExtendingCorners(1) ! put the old first node on the empty position 
            ! compressing mesh
            CalParams%MovingMesh%MeshAreas(K, 2, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh) = GeoParams%MovingMeshCompressingCorners(1:NDIM * (NDIM - 1) + 2) ! fill whole array
            CalParams%MovingMesh%MeshAreas(K, 2, 1) = GeoParams%MovingMeshCompressingCorners(PositionNodeCompressing) ! replace first node with the one on the fixed boundary
            CalParams%MovingMesh%MeshAreas(K, 2, PositionNodeCompressing) = GeoParams%MovingMeshCompressingCorners(1) ! put the old first node on the empty position 
             
          end if
  
          ! moving mesh 
          CalParams%MovingMesh%MeshAreas(K, 3, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh) = GeoParams%MovingMeshMovingCorners(1:NDIM * (NDIM - 1) + 2) ! fill whole array
            
        end subroutine CheckMovingMeshData

        
        
        subroutine CheckAbsorbingBoundaryData()
        !**********************************************************************
        !
        ! Function:  To check the input data for the absorbing boundary algorithm
        !
        !**********************************************************************
        
        implicit none

          if ( .not. GeoParams%ApplyAbsorbingBoundary ) RETURN        

          CalParams%ApplyAbsorbingBoundary = GeoParams%ApplyAbsorbingBoundary
          CalParams%AbsorbingBoundaries%VBMaterialSet = GeoParams%AbsorbingBoundaryReferenceMaterialID
        
        end subroutine CheckAbsorbingBoundaryData
        
        
        subroutine CheckPrescribedVelocityData()
        !*****************************************************************************************
        !
        ! Function: check if there are conflict in the definition of prescribed velocity input
        !
        !*****************************************************************************************
        implicit none

        if ((GeoParams%ApplyNodalVelocityVolume).and.(GeoParams%ApplyNodalVelocitySurface)) then
          call GiveError('Nodal velocity cannot be assigned to both VOLUMES and SURFACES')
        elseif ((GeoParams%ApplyNodalVelocityVolume).and.(GeoParams%ApplyNodalVelocityLine)) then
          call GiveError('Nodal velocity cannot be assigned to both VOLUMES and LINES')
        elseif ((GeoParams%ApplyNodalVelocityVolume).and.(GeoParams%ApplyNodalVelocityPoint)) then 
          call GiveError('Nodal velocity cannot be assigned to both VOLUMES and POINTS')
        elseif ((GeoParams%ApplyNodalVelocitySurface).and.(GeoParams%ApplyNodalVelocityLine)) then
          call GiveError('Nodal velocity cannot be assigned to both SURFACES and LINES')
        elseif ((GeoParams%ApplyNodalVelocitySurface).and.(GeoParams%ApplyNodalVelocityPoint))then
          call GiveError('Nodal velocity cannot be assigned to both SURFACES and POINTS')
        elseif ((GeoParams%ApplyNodalVelocityLine).and.(GeoParams%ApplyNodalVelocityPoint))then
          call GiveError('Nodal velocity cannot be assigned to both LINES and POINTS')
        elseif ((GeoParams%ApplyMPVelocityVolume).and.(GeoParams%ApplyMPVelocitySurface)) then
          call GiveError('Material Point velocity cannot be assigned to both VOLUMES and SURFACES')
        elseif (((GeoParams%ApplyMPVelocityVolume).and.(.not.CalParams%ApplyMeshSmoothing)).or.&
        ((GeoParams%ApplyMPVelocitySurface).and.(.not.CalParams%ApplyMeshSmoothing))) then
          call GiveWarning('Prescribed material point velocity should be used in combination with moving mesh')
        end if
 
        end subroutine CheckPrescribedVelocityData

        
      end module ModReadGeometryData        