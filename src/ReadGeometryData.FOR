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
    !	Copyright (C) 2022  Members of the Anura3D MPM Research Community 
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
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr. 2022) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModReadCalculationData
      use ModFileIO
      use ModNURBS
    
      implicit none

        type GeometryParameterType
          real(REAL_TYPE), dimension(:), allocatable :: &
              ExcavatedElements!, & ! Excavation Tool: vector to store the relation between element IDs and excavated volumes
          real(REAL_TYPE), dimension(:,:), allocatable :: &
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
		      InitialVelocityonMP ! (Number of elments with initial velocities, vsize+1) store the element ID and the initial velocities of particles insde the element
              !PrescribedVeloElValue, & ! (PrescribedVeloNElem,vsize) stores the values of prescribed velocity at elements  
          
          !I moved this here to add a rank for the multipatch
          real(REAL_TYPE), dimension(:,:, :), allocatable :: PrescribedVeloElValue
          
          integer(INTEGER_TYPE), allocatable  :: &
              MovingMeshExtendingCorners(:, :), & ! dsize * (dsize - 1) + ExtraNodesMovingMesh  | (default = -1) ! multipatch generalization added a rank
              MovingMeshCompressingCorners(:, :), & ! dsize * (dsize - 1) + ExtraNodesMovingMesh | (default = -1) ! multipatch generalization added a rank 
              MovingMeshMovingCorners(:, :), & ! dsize * (dsize - 1) + ExtraNodesMovingMesh  | (default = -1)    ! multipatch generalization added a rank
              ContactLocation(:,:), &
              PrescribedVeloElemID(:, :),&!stores the ID of element on which velocity is prescribed ! multipatch generalization added a rank ,: ! REMOVED RANK ! Added the rank back again
              PrescribedVeloElDirection(:,:, :), & !0/1 define the direction where velocity is prescribed ! multipatch generalization added a rank
              WaterSurfaceMaterialID(:) ! Material Index on which a given water surface is assigned
          
          
          
          
          integer(INTEGER_TYPE) :: &    
              NABSurfaceNodesSolid, & ! total number of absorbing nodes belonging to surfaces, lines or points for solid
              NABLineNodesSolid, &
              NABPointNodesSolid, &
              NABSurfaceNodesLiquid, & ! total number of absorbing nodes belonging to surfaces, lines or points for liquid
              NABLineNodesLiquid, &
              NABPointNodesLiquid, &
              NABSurfaceNodesGas, & ! total number of absorbing nodes belonging to surfaces, lines or points for gas
              NABLineNodesGas, &
              NABPointNodesGas, &
              ExtraNodesMovingMesh = 2, & !Number of extra nodes for moving mesh areas ! this is kind of a dummy variable to generalize the moving mesh to 2D/3D
              AbsorbingBoundaryReferenceMaterialID = 0, &  ! the material_id of the reference material for the absorbing boundaries    
              NumberWaterSurfaceMaterials!, &
              !PrescribedVeloNElem = 0 !Number of element on which velocity is prescribed at MP ! this probably needs multipatch generalization

          
          integer(INTEGER_TYPE), dimension(:), allocatable :: PrescribedVeloNElem ! --> added a rank so that we can allocate this based on the number of patches
          
              !MovingMeshReferenceMaterialID, & ! the material_id of the reference material for the moving mesh ! needs multipatch generalization. Added a rank
              !NMovingElements, & !Number of elements with initial velocity condition ! this probably needs multipatch generalization
		      !
          
          
          ! I took these out from the previous initialization because they need multipatch generalization 
          integer(INTEGER_TYPE), allocatable :: & 
              MovingMeshReferenceMaterialID(:), &
              NMovingElements(:)
              !PrescribedVeloNElem(:), &
          
          character(len=255), allocatable, dimension(:) :: &
              MovingMeshDirection
          
          character(len=255), allocatable :: &
              ContactMaterials(:,:) , &
              WaterSurfaceFileName(:)  ! Name and location of text file including initial water surface

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
            
            ! allocate variables
            allocate( &
              GeoParams%MovingMeshExtendingCorners(NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, Counters%NPatches), & ! I added a rank for the number of patches
              GeoParams%MovingMeshCompressingCorners(NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, Counters%NPatches), & ! I added a rank for the number of patches
              GeoParams%MovingMeshMovingCorners(NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, Counters%NPatches), & ! I added a rank for the number of patches
              GeoParams%LocalDampingFactorElement(NURBS%Maximum_nel_NURBS, Counters%NPatches) )!(Counters%NEl) ) ! I added a rank for the number of patches
            
            allocate( GeoParams%MovingMeshReferenceMaterialID(Counters%NPatches) ) ! integer type variable 
            
            allocate( GeoParams%MovingMeshDirection(Counters%NPatches) ) ! this is a character type variable 
        
            
            
            ! initialize variables 
              GeoParams%MovingMeshExtendingCorners = -1 ! this has an added rank for multipatch
              GeoParams%MovingMeshCompressingCorners = -1 ! this has an added rank for multipatch
              GeoParams%MovingMeshMovingCorners = -1 ! this has an added rank for multipatch
              
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
              

              ! particularly reading prescribed velocity, moving mesh, absorbing boundaries, and contact details
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
                // trim(Anura3D_v2022)     // '.' )

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
          integer(INTEGER_TYPE) :: ios, stat, IError ! used for error control
          character(len=255) :: DumS = '' ! reading strings
          character(len=21) :: messageIOS = 'GOM file: Can''t read '
          
          integer(INTEGER_TYPE) :: NDimension = -1
          logical :: Axisymmetric = .false.
          logical :: Cylindric = .false.
          character(len=255) :: ReadElementType = ''
          character(len=255) :: Formulation = ''
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: ii
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1 !we need to generalize the inout structure within subroutines
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: ILoadSystem
          character(len=255) :: TractionDimension = ''
          
          integer(INTEGER_TYPE) :: ILoadSystem_Temporary = 1
          
          
          ! I moved this from InitialiseMeshData() to InitialiseDimension() (i.e., here) 
          ! Not all of them are used...
          !integer(INTEGER_TYPE) :: I, J, L, IDElement, IDConnectivity, IFxN, &
          !           NElements, NNodes, NFixNod, FixNodID, NRemFixNod, FixRemNodID, &
          !           NLoadedElementSidesSolidNodes, NLoadedElementSidesWaterNodes, NLoadedElementSidesGasNodes, &
          !           NLoadedElementSidesSolidMatPoints, NLoadedElementSidesWaterMatPoints, NLoadedElementSidesGasmatPoints, &
          !           SoilSurfaceNumberofSides, PhreaticSurfaceNumberofSides, &
          !           NLoadedElementSidesSolidNodesB, NLoadedElementSidesWaterNodesB, NLoadedElementSidesGasNodesB, &
          !           NLoadedElementSidesSolidMatPointsB, NLoadedElementSidesWaterMatPointsB, NLoadedElementSidesGasmatPointsB 
          
          
          
          
          
          
          !! NURBS variables 
          !!real(REAL_TYPE), allocatable, dimension(:) :: XiKnotEntriesUNIQUE
          !!integer(INTEGER_TYPE) :: NXiKnotEntriesUNIQUE
          !real(REAL_TYPE), dimension(NDIM) :: LowestPointOnASurfaceDefiningPhysicalDomain
          !real(REAL_TYPE), dimension(NDIM) :: HighestPointOnASurfaceDefiningPhysicalDomain
          
          
          
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
              
              
              call SetDimension(NDimension, Axisymmetric, Cylindric)
              
              
            
              !if () then 
                  

            
              
             ! ! Physical definition of the domain (lowest point) 
             ! else if (trim(BName) == '$$LOWEST_POINT_ON_A_SURFACE_DEFINING_PHYSICAL_DOMAIN') then 
             !   read(GOMunit, *)    LowestPointOnASurfaceDefiningPhysicalDomain(1:NDIM)
             !   
             ! 
             !   
             !   
             !   
             !! Physical definition of the domain (highest point)   
             ! else if (trim(BName) == '$$HIGHEST_POINT_ON_A_SURFACE_DEFINING_PHYSICAL_DOMAIN') then 
             !     read(GOMunit, *)    HighestPointOnASurfaceDefiningPhysicalDomain(1:NDIM)
             ! 
              
              
              ! read the number of patches as part of the geometry ! MULTIPATCH IMPLEMENTATION 
            else if (trim(BName) == '$$NUMBER_OF_PATCHES') then 
                read(GOMunit, *, iostat=ios) Counters%NPatches !NPATCHES is a protected global variable (similar to NVECTOR)
                ! we specify that Counters%NPatches is an integer in the same place we indentify NXiKnotEntries
                
                ! we need to allocate integer variables pertinent to each patch here 
                ! Xi integer variables 
                allocate( NURBS%NXiKnotEntries(Counters%NPatches)  , stat = IError)
                allocate( NURBS%NXiKnotOrder(Counters%NPatches)   , stat = IError)
                NURBS%NXiKnotEntries = 0 ! initialise
                NURBS%NXiKnotOrder = 0 ! initialise 
                ! Eta integer variables
                allocate( NURBS%NEtaKnotEntries(Counters%NPatches)  , stat = IError)
                allocate( NURBS%NEtaKnotOrder(Counters%NPatches)    , stat = IError)
                NURBS%NEtaKnotEntries = 0 ! initialise
                NURBS%NEtaKnotOrder = 0 ! initialise
                ! Zeta integer variables
                allocate( NURBS%NZetaKnotEntries(Counters%NPatches)  , stat = IError)
                allocate( NURBS%NZetaKnotOrder(Counters%NPatches)    , stat = IError)
                NURBS%NZetaKnotEntries = 0 ! initialise
                NURBS%NZetaKnotOrder = 0 ! initialise 
                
                
                ! Volumetric locking smoothing knots 
                ! we need to allocate integer variables pertinent to each patch here 
                ! Xi integer variables 
                allocate( NURBS%NXiKnotEntries_VolLockSmooth(Counters%NPatches)  , stat = IError)
                allocate( NURBS%NXiKnotOrder_VolLockSmooth(Counters%NPatches)   , stat = IError)
                NURBS%NXiKnotEntries_VolLockSmooth = 0 ! initialise
                NURBS%NXiKnotOrder_VolLockSmooth = 0 ! initialise 
                ! Eta integer variables
                allocate( NURBS%NEtaKnotEntries_VolLockSmooth(Counters%NPatches)  , stat = IError)
                allocate( NURBS%NEtaKnotOrder_VolLockSmooth(Counters%NPatches)    , stat = IError)
                NURBS%NEtaKnotEntries_VolLockSmooth = 0 ! initialise
                NURBS%NEtaKnotOrder_VolLockSmooth = 0 ! initialise
                ! Zeta integer variables
                allocate( NURBS%NZetaKnotEntries_VolLockSmooth(Counters%NPatches)  , stat = IError)
                allocate( NURBS%NZetaKnotOrder_VolLockSmooth(Counters%NPatches)    , stat = IError)
                NURBS%NZetaKnotEntries_VolLockSmooth = 0 ! initialise
                NURBS%NZetaKnotOrder_VolLockSmooth = 0 ! initialise 
                
                
                ! number of element items for each patch
                allocate(Counters%NEl(Counters%NPatches)  , stat = IError) !NURBS%nel_NURBS
                allocate(NURBS%nnp_NURBS(Counters%NPatches)  , stat = IError)
                allocate(NURBS%nen_NURBS(Counters%NPatches)  , stat = IError)
                Counters%NEl = 0 !NURBS%nel_NURBS = 0
                NURBS%nnp_NURBS = 0
                NURBS%nen_NURBS = 0
                
                ! number of element items for each patch --> Reduced order volumetric locking smoothing
                allocate(NURBS%nel_NURBS_VolLockSmooth(Counters%NPatches)  , stat = IError)
                allocate(NURBS%nnp_NURBS_VolLockSmooth(Counters%NPatches)  , stat = IError)
                allocate(NURBS%nen_NURBS_VolLockSmooth(Counters%NPatches)  , stat = IError)
                NURBS%nel_NURBS_VolLockSmooth = 0
                NURBS%nnp_NURBS_VolLockSmooth = 0
                NURBS%nen_NURBS_VolLockSmooth = 0
                
                ! number of active element for each patch
                allocate(Counters%NAEl(Counters%NPatches)  , stat = IError)
                Counters%NAEl = 0
                
                ! number of control points for each patch 
                allocate(Counters%NodTot(Counters%NPatches)  , stat = IError)
                Counters%NodTot = 0 !NURBS%NControlPoints = 0
                
                allocate(NURBS%NControlPoints_VolLockSmooth(Counters%NPatches)  , stat = IError)
                NURBS%NControlPoints_VolLockSmooth = 0
                
                
                ! number of active element for each patch
                !allocate(NFixities(Counters%NPatches), stat = IError)
                !NFixities = 0
                
                
                ! number of elements to prescribed material point velocity --> might need 3D generalization
                allocate(NURBS%NElemMPPrescribedVelocity(Counters%NPatches), stat = IError)
                NURBS%NElemMPPrescribedVelocity = 0
                
                
                
                allocate(NURBS%NNodePrescribedVelocity(Counters%NPatches), stat = IError)
                NURBS%NNodePrescribedVelocity = 0
                
                
                allocate( NURBS%nn_NURBS_NumberOfUnivariateXiKnots(Counters%NPatches), stat = IError)
                allocate( NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(Counters%NPatches), stat = IError)
                allocate( NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(Counters%NPatches), stat = IError)
                NURBS%nn_NURBS_NumberOfUnivariateXiKnots = 0
                NURBS%mm_NURBS_NumberOfUnivariateEtaKnots = 0
                NURBS%oo_NURBS_NumberOfUnivariateZetaKnots = 0
                
                ! --> volumetric locking smoothing
                allocate( NURBS%nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(Counters%NPatches), stat = IError)
                allocate( NURBS%mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(Counters%NPatches), stat = IError)
                allocate( NURBS%oo_NURBS_NumberOfUnivariateZetaKnots_VolLockSmooth(Counters%NPatches), stat = IError)
                NURBS%nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth = 0
                NURBS%mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth = 0
                NURBS%oo_NURBS_NumberOfUnivariateZetaKnots_VolLockSmooth = 0
                
                
                allocate( NURBS%NXiKnotEntriesUNIQUE(Counters%NPatches), stat = IError)
                NURBS%NXiKnotEntriesUNIQUE = 0
                allocate( NURBS%NEtaKnotEntriesUNIQUE(Counters%NPatches), stat = IError)
                NURBS%NEtaKnotEntriesUNIQUE = 0
                allocate( NURBS%NZetaKnotEntriesUNIQUE(Counters%NPatches), stat = IError)
                NURBS%NZetaKnotEntriesUNIQUE = 0
                
                
                               
                !allocate( NumberOfGlobalControlPointsUniqueMultipatch(Counters%NPatches), stat = IError)
                NURBS%NumberOfGlobalControlPointsUniqueMultipatch = 0
                
                allocate( NURBS%NumberOfUniqueNodes(Counters%NPatches), stat = IError)
                NURBS%NumberOfUniqueNodes = 0
                
                
                !Number of load systems for each patch 
                allocate( NURBS%NLoadSystems_Traction(Counters%NPatches), stat = IError)
                NURBS%NLoadSystems_Traction = 0
                
                
                
            
                
                
                
                
                
                
                
                !--------------------------------------------------------------------------
            else if (trim(BName) == '$$PATCH_ID') then 
                ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                IPatch = -1
                ! read patch ID from the GOM file
                read(GOMunit, *) IPatch ! --> now we knot patch index so we can read in the following: 
                                        ! 1) Number of control points for each patch
                                        ! 1) Number of knot entries (xi, eta, zeta) 
                                        ! 2) Order of knot (xi, eta, zeta)
                                        ! 3) 
                
                
                
            else if (trim(BName) == '$$NUMBER_OF_LOAD_SYSTEMS') then 
                ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                !ILoadSystem = -1
                ! read patch ID from the GOM file
                
                !alloca
                
                read(GOMunit, *) NURBS%NLoadSystems_Traction(IPatch)
                
                
                
           
            
                
                
                
                !------------------------------------------------------------------------
                
                
                
                
  
            else if (trim(BName) == '$$FORMULATION') then
              read(GOMunit, *, iostat=ios) DumS
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == SINGLE_POINT .or. DumS == DOUBLE_POINT, 'GOM file: ' //trim(BName)// ' must be equal to ' & 
                           // SINGLE_POINT // ' or ' // DOUBLE_POINT //'.' )
              Formulation = DumS
              
              call SetFormulation(Formulation)
            
              
              ! we need to know the maximum number of load systems 
              NURBS%Maximum_NLoadSystems_Traction = maxval(NURBS%NLoadSystems_Traction)
              

              allocate( NURBS%NControlPoints_Traction(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              NURBS%NControlPoints_Traction = 0
              
              allocate( NURBS%NURBS_Traction_Dimension(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              NURBS%NURBS_Traction_Dimension = 0
              
              allocate( NURBS%nel_NURBS_Traction(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              NURBS%nel_NURBS_Traction = 0
              
              allocate( NURBS%nnp_NURBS_Traction(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              NURBS%nnp_NURBS_Traction = 0
              
              allocate( NURBS%nen_NURBS_Traction(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              NURBS%nen_NURBS_Traction = 0  
              

              !------------------------------------------------------------------------
              
              allocate( NURBS%NLoadSystems_Traction_NumberOfElements(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              NURBS%NLoadSystems_Traction_NumberOfElements = 0
              
                !------------------------------------------------------------------------
                
              allocate( Counters%NLoadedElementSidesSolidNodes(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              Counters%NLoadedElementSidesSolidNodes = 0 
              
              allocate( Counters%NLoadedElementSidesWaterNodes(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              Counters%NLoadedElementSidesWaterNodes = 0
              
              allocate( Counters%NLoadedElementSidesSolidMatPoints(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              Counters%NLoadedElementSidesSolidMatPoints = 0
              
              allocate( Counters%NLoadedElementSidesWaterMatPoints(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              Counters%NLoadedElementSidesWaterMatPoints = 0
                
                
              
              
              allocate( Counters%NLoadedElementSidesSolid(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              Counters%NLoadedElementSidesSolid = 0
              
              allocate( Counters%NLoadedElementSidesWater(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError )
              Counters%NLoadedElementSidesWater = 0
              
              
                
                !------------------------------------------------------------------------
                
            else if (trim(BName) == '$$PATCH_ID') then
                ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                IPatch = -1
                ! read patch ID from the GOM file
                read(GOMunit, *) IPatch ! --> now we knot patch index so we can read in the following: 
                                        ! 1) Number of control points for each patch
                                        ! 1) Number of knot entries (xi, eta, zeta) 
                                        ! 2) Order of knot (xi, eta, zeta)
                                        ! 3) 
                
                
                
                
                ! actual solid shape function knot vectors
                    ! Xi vector 
            else if (trim(BName) == '$$XI_NUMBER_OF_KNOTS_AND_ORDER') then 
                read(GOMunit, *) NURBS%NXiKnotEntries(IPatch), NURBS%NXiKnotOrder(IPatch)  ! reading 1) number of knot entries, 2) knot order 
                ! allocating allocatable arrays
                !allocate(XiKnotEntries(NXiKnotEntries(IPatch_Temporary)), stat = IError)
                !
                !! find unique knots items in the knot vector 
                !call FindUniqueEnteriesInVector(XiKnotEntries, NXiKnotEntries, XiKnotEntriesUNIQUE, NXiKnotEntriesUNIQUE)
                
                ! number of univariate basis functions for the Xi direction 
                NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch) = NURBS%NXiKnotEntries(IPatch) - NURBS%NXiKnotOrder(IPatch) - 1
                
                    
                ! Eta vector
            else if (trim(BName) == '$$ETA_NUMBER_OF_KNOTS_AND_ORDER') then 
                read(GOMunit, *) NURBS%NEtaKnotEntries(IPatch), NURBS%NEtaKnotOrder(IPatch)  ! reading 1) number of knot entries, 2) knot order 
                ! allocating allocatable arrays
                !allocate(EtaKnotEntries(NEtaKnotEntries(IPatch_Temporary)), stat = IError)
                ! number of univariate basis functions for the Eta direction 
                NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch) = NURBS%NEtaKnotEntries(IPatch) - NURBS%NEtaKnotOrder(IPatch) - 1
                
                
                ! Zeta vector
            else if (trim(BName) == '$$ZETA_NUMBER_OF_KNOTS_AND_ORDER') then 
                read(GOMunit, *) NURBS%NZetaKnotEntries(IPatch), NURBS%NZetaKnotOrder(IPatch)  ! reading 1) number of knot entries, 2) knot order 
                ! allocating allocatable arrays
                !allocate(ZetaKnotEntries(NZetaKnotEntries(IPatch_Temporary)), stat = IError)
                ! number of univariate basis functions for the Zeta direction 
                NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch) = NURBS%NZetaKnotEntries(IPatch) - NURBS%NZetaKnotOrder(IPatch) - 1
     
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                
                else if (trim(BName) == '$$ELEMENTTYPE') then
              read(GOMunit, *, iostat=ios) DumS
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == TRI3 .or. DumS == TRI6 .or. DumS == QUAD4 .or. DumS == QUAD4_NURBS .or. DumS == HEXA_NURBS .or. DumS == QUAD8 .or. &
                           DumS == TETRA4 .or. DumS == TETRA10 .or. DumS == HEXA8 .or. DumS == HEXA20 .or. DumS == TETRAOLD, 'GOM file: ' //trim(BName)// &
                           ' must be equal to ' // TRI3 // ' or ' // TRI6 // ' or ' // QUAD4 //  ' or ' // QUAD4_NURBS // ' or ' // QUAD8 // ' or ' // &
                           TETRA4 // ' or ' // TETRA10 // ' or ' // HEXA_NURBS // ' or ' // HEXA8 // ' or ' // HEXA20  // ' or ' // TETRAOLD  // '.' )
              ReadElementType = DumS
              
              call SetElementType(ReadElementType, NURBS%NXiKnotOrder(IPatch), NURBS%NEtaKnotOrder(IPatch), NURBS%NZetaKnotOrder(IPatch), &
                  NURBS%NXiKnotOrder_VolLockSmooth(IPatch), NURBS%NEtaKnotOrder_VolLockSmooth(IPatch), NURBS%NZetaKnotOrder_VolLockSmooth(IPatch))
              
                
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Volumetric locking orders
                
            else if (trim(BName) == '$$XI_ORDER_VOLLOCK_SMOOTH') then               
                read(GOMunit, *) NURBS%NXiKnotOrder_VolLockSmooth(IPatch)  ! reading 1) number of knot entries, 2) knot order 
                
                ! removing the two tail ends of the knot vector and adding the two tail ends of the reduced order knot vector
                NURBS%NXiKnotEntries_VolLockSmooth(IPatch) = NURBS%NXiKnotEntries(IPatch) -&
                    (2*(NURBS%NXiKnotOrder(IPatch)+1)) +&
                    (2*(NURBS%NXiKnotOrder_VolLockSmooth(IPatch)+1))
                
                allocate(NURBS%XiKnotEntries_VolLockSmooth(NURBS%NXiKnotEntries_VolLockSmooth(IPatch), Counters%NPatches) ) 
                
                ! number of univariate basis functions for the Xi direction 
                NURBS%nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch) = NURBS%NXiKnotEntries_VolLockSmooth(IPatch) -&
                                                                           NURBS%NXiKnotOrder_VolLockSmooth(IPatch) - 1
                
                
                
                
            else if (trim(BName) == '$$ETA_ORDER_VOLLOCK_SMOOTH') then 
                read(GOMunit, *) NURBS%NEtaKnotOrder_VolLockSmooth(IPatch)  ! reading 1) number of knot entries, 2) knot order 
                
                ! removing the two tail ends of the knot vector and adding the two tail ends of the reduced order knot vector
                NURBS%NEtaKnotEntries_VolLockSmooth(IPatch) = NURBS%NEtaKnotEntries(IPatch) -&
                    (2*(NURBS%NEtaKnotOrder(IPatch)+1)) +&
                    (2*(NURBS%NEtaKnotOrder_VolLockSmooth(IPatch)+1))
                
                allocate(NURBS%EtaKnotEntries_VolLockSmooth(NURBS%NEtaKnotEntries_VolLockSmooth(IPatch), Counters%NPatches) )  
                
                NURBS%mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch) = NURBS%NEtaKnotEntries_VolLockSmooth(IPatch) -&
                                                                                        NURBS%NEtaKnotOrder_VolLockSmooth(IPatch) - 1
            
                
            
            else if (trim(BName) == '$$ZETA_ORDER_VOLLOCK_SMOOTH') then 
                read(GOMunit, *) NURBS%NZetaKnotOrder_VolLockSmooth(IPatch)  ! reading 1) number of knot entries, 2) knot order 
                
                
                ! removing the two tail ends of the knot vector and adding the two tail ends of the reduced order knot vector
                NURBS%NZetaKnotEntries_VolLockSmooth(IPatch) = NURBS%NZetaKnotEntries(IPatch) -&
                    (2*(NURBS%NZetaKnotOrder(IPatch)+1)) +&
                    (2*(NURBS%NZetaKnotOrder_VolLockSmooth(IPatch)+1))
                
                allocate(NURBS%ZetaKnotEntries_VolLockSmooth(NURBS%NZetaKnotEntries_VolLockSmooth(IPatch), Counters%NPatches) )  
                
                NURBS%oo_NURBS_NumberOfUnivariateZetaKnots_VolLockSmooth(IPatch) = NURBS%NZetaKnotEntries_VolLockSmooth(IPatch) -&
                                                               NURBS%NZetaKnotOrder_VolLockSmooth(IPatch) - 1
            
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                      
                ! Xi vector 
            !else if (trim(BName) == '$$XI_NUMBER_OF_KNOTS_AND_ORDER_TRACTION') then 
            !    read(GOMunit, *) NXiKnotEntries_Traction, NXiKnotOrder_Traction  ! reading 1) number of knot entries, 2) knot order 
            !    !                   (IPatch)               (IPatch)
            !    ! allocating allocatable arrays
            !    !allocate(XiKnotEntries(NXiKnotEntries(IPatch_Temporary)), stat = IError)
            !    !
            !    !! find unique knots items in the knot vector 
            !    !call FindUniqueEnteriesInVector(XiKnotEntries, NXiKnotEntries, XiKnotEntriesUNIQUE, NXiKnotEntriesUNIQUE)
            !    
            !    ! number of univariate basis functions for the Xi direction 
            !    nn_NURBS_NumberOfUnivariateXiKnots_Traction = NXiKnotEntries_Traction - NXiKnotOrder_Traction - 1
            !    !           (IPatch)                                        (IPatch)                (IPatch)
                    
                ! Eta vector
            !else if (trim(BName) == '$$ETA_NUMBER_OF_KNOTS_AND_ORDER_TRACTION') then 
            !    read(GOMunit, *) NEtaKnotEntries_Traction, NEtaKnotOrder_Traction  ! reading 1) number of knot entries, 2) knot order 
            !                       !     (IPatch)   (IPatch)
            !    ! allocating allocatable arrays
            !    !allocate(EtaKnotEntries(NEtaKnotEntries(IPatch_Temporary)), stat = IError)
            !    ! number of univariate basis functions for the Eta direction 
            !    mm_NURBS_NumberOfUnivariateEtaKnots_Traction = NEtaKnotEntries_Traction - NEtaKnotOrder_Traction - 1
            !    !           (IPatch)                        (IPatch)        (IPatch)    
            !    
            !    
                
                
                
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                
                
                
            else if (trim(BName) == '$$NUMBER_OF_CONTROL_POINTS') then
                read(GOMunit, *) Counters%NodTot(IPatch) !NURBS%NControlPoints(IPatch)  ! reading 1) NElements -> number of elements (knot spans in the case of NURBS), 
                                                 !         2) NNodes -> number of nodes (control points in the case of NURBS)
                
                
                if (NDIM == 2) then 
                    
                Counters%NEl(IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) * & !NURBS%nel_NURBS(IPatch)
                                    (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NEtaKnotOrder(IPatch)) !* &                
                                        !(oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)-NZetaKnotOrder(IPatch))
        
                NURBS%nnp_NURBS(IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch) &
                                    * NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch) !&
                                    !*oo_NURBS_NumberOfUnivariateZetaKnots(IPatch) !number of global basis functions (global here refers to its global domain within the 'super' element)

                NURBS%nen_NURBS(IPatch) = (NURBS%NXiKnotOrder(IPatch)+1) * (NURBS%NEtaKnotOrder(IPatch)+1) !* (NZetaKnotOrder(IPatch)+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                    
                    
                    

                else if (NDIM == 3) then 
                
                Counters%NEl(IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) * & !NURBS%nel_NURBS(IPatch)
                                    (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NEtaKnotOrder(IPatch)) * &                
                                        (NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)-NURBS%NZetaKnotOrder(IPatch))
        
                NURBS%nnp_NURBS(IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch) &
                                    *NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch) &
                                    *NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch) !number of global basis functions (global here refers to its global domain within the 'super' element)

                NURBS%nen_NURBS(IPatch) = (NURBS%NXiKnotOrder(IPatch)+1) * (NURBS%NEtaKnotOrder(IPatch)+1) * (NURBS%NZetaKnotOrder(IPatch)+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                
                
                
                end if 
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                
            else if (trim(BName) == '$$NUMBER_OF_CONTROL_POINTS_VOLLOCK_SMOOTH') then
                read(GOMunit, *) NURBS%NControlPoints_VolLockSmooth(IPatch)
                
                
                if (NDIM == 2) then 
                    
                NURBS%nel_NURBS_VolLockSmooth(IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch)-NURBS%NXiKnotOrder_VolLockSmooth(IPatch)) * &
                                    (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch)-NURBS%NEtaKnotOrder_VolLockSmooth(IPatch)) !* &                
                                        !(oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)-NZetaKnotOrder(IPatch))
            
                NURBS%nnp_NURBS_VolLockSmooth(IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch) &
                                    * NURBS%mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch) !&
                                    !*oo_NURBS_NumberOfUnivariateZetaKnots(IPatch) !number of global basis functions (global here refers to its global domain within the 'super' element)
            
                NURBS%nen_NURBS_VolLockSmooth(IPatch) = (NURBS%NXiKnotOrder_VolLockSmooth(IPatch)+1) * (NURBS%NEtaKnotOrder_VolLockSmooth(IPatch)+1) !* (NZetaKnotOrder(IPatch)+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                    
                    
                    
            
                else if (NDIM == 3) then 
                
                NURBS%nel_NURBS_VolLockSmooth(IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch)-NURBS%NXiKnotOrder_VolLockSmooth(IPatch)) * &
                                    (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch)-NURBS%NEtaKnotOrder_VolLockSmooth(IPatch)) * &                
                                        (NURBS%oo_NURBS_NumberOfUnivariateZetaKnots_VolLockSmooth(IPatch)-NURBS%NZetaKnotOrder_VolLockSmooth(IPatch))
            
                NURBS%nnp_NURBS_VolLockSmooth(IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch) &
                                    *NURBS%mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch) &
                                    *NURBS%oo_NURBS_NumberOfUnivariateZetaKnots_VolLockSmooth(IPatch) !number of global basis functions (global here refers to its global domain within the 'super' element)
            
                NURBS%nen_NURBS_VolLockSmooth(IPatch) = (NURBS%NXiKnotOrder_VolLockSmooth(IPatch)+1) * (NURBS%NEtaKnotOrder_VolLockSmooth(IPatch)+1) * (NURBS%NZetaKnotOrder_VolLockSmooth(IPatch)+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                
                
                
                end if 
                
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! NURBS traction flags

                else if (trim(BName) == '$$TRACTION_LOAD_SYSTEM') then
                ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                ILoadSystem = -1
                ! read patch ID from the GOM file
                read(GOMunit, *) ILoadSystem ! --> now we knot patch index so we can read in the following: 
                                        ! 1) Number of control points for each patch
                                        ! 1) Number of knot entries (xi, eta, zeta) 
                                        ! 2) Order of knot (xi, eta, zeta)
                                        ! 3) 
                
                
                
                else if (trim(BName) == '$$NUMBER_OF_TRACTION_ELEMENTS') then 
                ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                !ILoadSystem = -1
                ! read patch ID from the GOM file
                read(GOMunit, *) NURBS%NLoadSystems_Traction_NumberOfElements(ILoadSystem, IPatch)
                
                
                
                else if (trim(BName) == '$$NUMBER_OF_CONTROL_POINTS_TRACTION') then
                read(GOMunit, *) NURBS%NControlPoints_Traction(ILoadSystem, IPatch)  ! reading 1) NElements -> number of elements (knot spans in the case of NURBS), 
                                !(IPatch)
                read(GOMunit, *) TractionDimension
                
                if (NURBS%NControlPoints_Traction(ILoadSystem, IPatch)==0) cycle
                
                if (NDIM == 2) then !--> 1D
                    ! if NDIM = 2 then we are looking at a line which is 1D
                    
                    
                    if (TractionDimension == 'XI') then 
                        
                        
                        NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) = 1
                    
                    
                    NURBS%nel_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) !_Traction _Traction
                    !   (IPatch)               (IPatch)                           (IPatch)
              
                    
                    NURBS%nnp_NURBS_Traction(ILoadSystem,IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch) !_Traction
                    !(IPatch)                   (IPatch)
                    !number of global basis functions (global here refers to its global domain within the 'super' element)
                    
                
                    NURBS%nen_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%NXiKnotOrder(IPatch)+1)  
                    !(IPatch)               (IPatch)
                    !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                    
                    
                    elseif (TractionDimension == 'ETA') then
                        
                        NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) = 2
                    
                        
                    NURBS%nel_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) !_Traction _Traction
                    !   (IPatch)               (IPatch)                           (IPatch)
              
                    
                    NURBS%nnp_NURBS_Traction(ILoadSystem,IPatch) = NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch) !nn_NURBS_NumberOfUnivariateXiKnots_Traction
                    !(IPatch)                   (IPatch)
                    !number of global basis functions (global here refers to its global domain within the 'super' element)
                    
                
                    NURBS%nen_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%NXiKnotOrder(IPatch)+1)   !_Traction
                    !(IPatch)               (IPatch)
                    !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                    
                    
                    end if 
                    
                    
                else if (NDIM == 3) then 
                
                    
                    
                    if (TractionDimension == 'XI_ETA') then
                        
                        NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) = 1
                        
                    NURBS%nel_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) * &
                                    (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NEtaKnotOrder(IPatch)) 
        
                
                    NURBS%nnp_NURBS_Traction(ILoadSystem,IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch) &
                                    * NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch) !_Traction  
                    !number of global basis functions (global here refers to its global domain within the 'super' element)

                
                    NURBS%nen_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%NXiKnotOrder(IPatch)+1) * (NURBS%NEtaKnotOrder(IPatch)+1) 
                    !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                    
                    
                    
                    else if (TractionDimension == 'XI_ZETA') then
                        
                        NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) = 2
                    
                        
                    NURBS%nel_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) * &
                                    (NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)-NURBS%NZetaKnotOrder(IPatch)) 
        
                
                    NURBS%nnp_NURBS_Traction(ILoadSystem,IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch) &
                                    * NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)  
                    !number of global basis functions (global here refers to its global domain within the 'super' element)

                
                    NURBS%nen_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%NXiKnotOrder(IPatch)+1) * (NURBS%NZetaKnotOrder(IPatch)+1) 
                    !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                        
                        
                        
                    else if (TractionDimension == 'ETA_ZETA') then
                        
                        NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) = 3
                    
                        
                    NURBS%nel_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NEtaKnotOrder(IPatch)) * &
                                    (NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)-NURBS%NZetaKnotOrder(IPatch)) 
        
                
                    NURBS%nnp_NURBS_Traction(ILoadSystem,IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch) &
                                    * NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)  
                    !number of global basis functions (global here refers to its global domain within the 'super' element)

                
                    NURBS%nen_NURBS_Traction(ILoadSystem,IPatch) = (NURBS%NEtaKnotOrder(IPatch)+1) * (NURBS%NZetaKnotOrder(IPatch)+1) 
                    !number of local basis functions (local here refers to a knot span i.e. accross one single element)
                        
                            
                        
                    end if 
                    
                    
                end if 
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! NURBS traction flags
                
            
                
            !            ! Xi vector 
            !else if (trim(BName) == '$$XI_NUMBER_OF_KNOTS_AND_ORDER_TRACTION') then 
            !    read(GOMunit, *) NXiKnotEntries_Traction, NXiKnotOrder_Traction  ! reading 1) number of knot entries, 2) knot order 
            !    ! (IPatch)                                  !(IPatch)
            !    ! allocating allocatable arrays
            !    !allocate(XiKnotEntries(NXiKnotEntries(IPatch_Temporary)), stat = IError)
            !    !
            !    !! find unique knots items in the knot vector 
            !    !call FindUniqueEnteriesInVector(XiKnotEntries, NXiKnotEntries, XiKnotEntriesUNIQUE, NXiKnotEntriesUNIQUE)
            !    
            !    ! number of univariate basis functions for the Xi direction 
            !    nn_NURBS_NumberOfUnivariateXiKnots_Traction = NXiKnotEntries_Traction - NXiKnotOrder_Traction - 1
            !    !(IPatch)                                       (IPatch)                (IPatch)
            !        
            !    ! Eta vector
            !else if (trim(BName) == '$$ETA_NUMBER_OF_KNOTS_AND_ORDER_TRACTION') then 
            !    read(GOMunit, *) NEtaKnotEntries_Traction, NEtaKnotOrder_Traction  ! reading 1) number of knot entries, 2) knot order 
            !        !                   (IPatch)                (IPatch)
            !    ! allocating allocatable arrays
            !    !allocate(EtaKnotEntries(NEtaKnotEntries(IPatch_Temporary)), stat = IError)
            !    ! number of univariate basis functions for the Eta direction 
            !    mm_NURBS_NumberOfUnivariateEtaKnots_Traction = NEtaKnotEntries_Traction - NEtaKnotOrder_Traction - 1
            !        !(IPatch)                                       (IPatch)                (IPatch)
            !    
                
                

                ! Number of load systems in each patch
            !else if (trim(BName) == '$$NUMBER_OF_LOAD_SYSTEMS') then
            !    read(GOMunit, *) NLoadSystems_Traction(IPatch)
                
                
            
                
            
            
            
            !else if (trim(BName) == '$$LOAD_SYSTEM') then
            !    read(GOMunit, *) ILoadSystems, TractionDimension
            
                
                !TractionDimension(ILoadSystem, IPatch) = 
                
                
                
                

            
                
                
                
                
                
                

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !else if (trim(BName) == '$$NUMBER_OF_START_FIXITY_LINE_SOLID') then 
            !    read(GOMunit, *) NFixities(IPatch)
                
                
                
                
                
                
                
                
                
            else if (trim(BName) == '$$NUMBER_OF_PRESCRIBED_MATERIAL_POINT_VELOCITY_SURFACE') then 
                read(GOMunit, *) NURBS%NElemMPPrescribedVelocity(IPatch) 
                ! --> we need to find the maximum for this so that we can allocate the shape of 
                !     of the prescribed velocity elements
    
                
            else if (trim(BName) == '$$NUMBER_OF_PRESCRIBED_NODAL_VELOCITY_SURFACE') then
                read(GOMunit, *) NURBS%NNodePrescribedVelocity(IPatch) 
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                
                
                
                
                
                
            !else if (trim(BName) == '$$COMPRESSING_MESH_PATCH') then 
            !    CompressingMeshPatchID = IPatch 
            !
            !    
            !    
            !    
            !    
            !    
            !
            !else if (trim(BName) == '$$MOVING_MESH_PATCH') then 
            !    MovingMeshPatchID = IPatch
                
                
                
                
                
                
             
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
            else if (trim(BName) == '$$XI_NUMBER_OF_GAUSS_POINTS') then 
                read(GOMunit, *) NXiGaussPoints
            
                
                
                
            else if (trim(BName) == '$$ETA_NUMBER_OF_GAUSS_POINTS') then 
                read(GOMunit, *) NEtaGaussPoints
                
                
                
                               
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                

            else if (trim(BName) == '$$NUMBER_OF_BOUND_MP_PRESCRIBED_VELOCITY_X') then 
                read(GOMunit, *) NBoundMPPrescVelo_X
                
                ! allocate based on number of number of loci prescribed with 5 of the following details:
                ! - Greater or lower than locus (1: greater/ 0: lower)
                ! - Locus limit (in x-direction)
                ! - Velocity in x
                ! - Velocity in y
                ! - Velocity in z (if 2D, then this doesn't matter, but you have to have a value)
                allocate(BoundPrescVelo_X(NBoundMPPrescVelo_X, 5), stat = IError )
                
                
                
                
            else if (trim(BName) == '$$NUMBER_OF_BOUND_MP_PRESCRIBED_VELOCITY_Y') then 
                read(GOMunit, *) NBoundMPPrescVelo_Y
                
                ! allocate based on number of number of loci prescribed with 5 of the following details:
                ! - Greater or lower than locus (1: greater/ 0: lower)
                ! - Locus limit (in y-direction)
                ! - Velocity in x
                ! - Velocity in y
                ! - Velocity in z (if 2D, then this doesn't matter, but you have to have a value)
                allocate(BoundPrescVelo_Y(NBoundMPPrescVelo_Y, 5), stat = IError)
                
                
            else if (trim(BName) == '$$NUMBER_OF_BOUND_MP_PRESCRIBED_VELOCITY_Z') then 
                read(GOMunit, *) NBoundMPPrescVelo_Z
                
                ! allocate based on number of number of loci prescribed with 5 of the following details:
                ! - Greater or lower than locus (1: greater/ 0: lower)
                ! - Locus limit (in y-direction)
                ! - Velocity in x
                ! - Velocity in y
                ! - Velocity in z (if 2D, then this doesn't matter, but you have to have a value)
                allocate(BoundPrescVelo_Z(NBoundMPPrescVelo_Z, 5), stat = IError)
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                       
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                ! we do not need this as a user input... we want this embedded in the calculation
                ! Xi vector gauss points and their location
                
                ! HARDCODED
            !    
            !else if (trim(BName) == '$$XI_NUMBER_OF_GAUSS_POINTS_AND_LOCATION') then 
            !    read(GOMunit, *) NXiGaussPoints  ! reading 1) number of gauss points, 2) gauss point locations                 
            !    allocate(XiGaussPoints(NXiGaussPoints), stat = IError)
            !    read(GOMunit, *) XiGaussPoints
            !    
            !
            !    ! Eta vector gauss points and their location
            !else if (trim(BName) == '$$ETA_NUMBER_OF_GAUSS_POINTS_AND_LOCATION') then 
            !    read(GOMunit, *) NEtaGaussPoints  ! reading 1) number of gauss points, 2) gauss point locations                 
            !    allocate(EtaGaussPoints(NEtaGaussPoints), stat = IError)
            !    read(GOMunit, *) EtaGaussPoints
            !    
            !    ! Zeta vector gauss points and their location
            !else if (trim(BName) == '$$ZETA_NUMBER_OF_GAUSS_POINTS_AND_LOCATION') then 
            !    read(GOMunit, *) NZetaGaussPoints  ! reading 1) number of gauss points, 2) gauss point locations                 
            !    allocate(ZetaGaussPoints(NZetaGaussPoints), stat = IError)
            !    read(GOMunit, *) ZetaGaussPoints
                
                ! Calculate the total number of gauss points in a 3D problem 
                !ELEMENTGAUSSPOINTS = NXiGaussPoints * NEtaGaussPoints * NZetaGaussPoints
            
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
            else if (trim(BName) == '$$FINISH') then
              close(GOMunit)  
              EXIT
            end if 
              
          end do

          ! old versions of the GOM file (before v2018.2) do not specify the dimension and elementtype, therefore it is set manually              
          if (NDimension == -1) NDimension = 3 ! this is the only available dimension before v2018.2
          if (ReadElementType == '') ReadElementType = TETRAOLD ! the only elementtype that was used before v2018.2
          if (Formulation == '') Formulation = SINGLE_POINT ! if nothing assumed single-point
          
          ! MULTIPATCH generlization aspect
          ! I need to know the maximum number of entries in xi, eta, and zeta, so that I can assign
          ! the size of XiKnotEntries, EtaKnotEntries, ZetaKnotEntries
          NURBS%Maximum_NControlPoints = maxval(Counters%NodTot) !NURBS%NControlPoints) ! <- finding maximum value
          NURBS%Maximum_NControlPoints_Traction = maxval(NURBS%NControlPoints_Traction) ! <- finding maximum value
          
          
          NURBS%Maximum_NXiKnotEntries = maxval(NURBS%NXiKnotEntries) ! <- finding maximum value in this NPATCHES vector
          NURBS%Maximum_NEtaKnotEntries = maxval(NURBS%NEtaKnotEntries) ! <- finding maximum value in this NPATCHES vector
          NURBS%Maximum_NZetaKnotEntries = maxval(NURBS%NZetaKnotEntries) ! <- finding maximum value in this NPATCHES vector
          
          
          allocate(NURBS%XiKnotEntries(NURBS%Maximum_NXiKnotEntries,Counters%NPatches), stat = IError)
          allocate(NURBS%EtaKnotEntries(NURBS%Maximum_NEtaKnotEntries,Counters%NPatches), stat = IError)
          allocate(NURBS%ZetaKnotEntries(NURBS%Maximum_NZetaKnotEntries,Counters%NPatches), stat = IError)
          
          !------traction
          !allocate(XiKnotEntries_Traction(NXiKnotEntries_Traction), stat = IError)
          !allocate(EtaKnotEntries_Traction(NEtaKnotEntries_Traction), stat = IError)
          !allocate(ZetaKnotEntries(NZetaKnotEntries_Traction), stat = IError)
          
          
          
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Traction
          
          

          
          if (NDIM == 2) then 
          do IPatch = 1, Counters%NPatches
              !NURBS%nel_NURBS(IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) * (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NEtaKnotOrder(IPatch)) !number of elements -> note 2D implementation = 2 elements in the example 
              Counters%NEl(IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) * (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NEtaKnotOrder(IPatch)) !number of elements -> note 2D implementation = 2 elements in the example 
              NURBS%nnp_NURBS(IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)*NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch) !number of global basis functions (global here refers to its global domain within the 'super' element)
              NURBS%nen_NURBS(IPatch) = (NURBS%NXiKnotOrder(IPatch)+1) * (NURBS%NEtaKnotOrder(IPatch)+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
          
              !nel_NURBS_VolLockSmooth(IPatch) = (nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch)-NXiKnotOrder_VolLockSmooth(IPatch)) * (mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch)-NEtaKnotOrder_VolLockSmooth(IPatch)) !number of elements -> note 2D implementation = 2 elements in the example 
              !nnp_NURBS_VolLockSmooth(IPatch) = nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch)*mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch) !number of global basis functions (global here refers to its global domain within the 'super' element)
              !nen_NURBS_VolLockSmooth(IPatch) = (NXiKnotOrder_VolLockSmooth(IPatch)+1) * (NEtaKnotOrder_VolLockSmooth(IPatch)+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)  
          
          end do 
          elseif (NDIM == 3) then 
          do IPatch = 1, Counters%NPatches
              ! actual material 
              Counters%NEl(IPatch) = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch)) * & !NURBS%nel_NURBS
                  (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NEtaKnotOrder(IPatch)) * &
                  (NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)-NURBS%NZetaKnotOrder(IPatch))
                  !number of elements -> note 2D implementation = 2 elements in the example 
              
              NURBS%nnp_NURBS(IPatch) = NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)*NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)*NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch) 
              !number of global basis functions (global here refers to its global domain within the 'super' element)
              
              NURBS%nen_NURBS(IPatch) = (NURBS%NXiKnotOrder(IPatch)+1) * (NURBS%NEtaKnotOrder(IPatch)+1) * (NURBS%NZetaKnotOrder(IPatch)+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
              
              
              
              
              
              ! volumetric locking variables 
              !nel_NURBS_VolLockSmooth(IPatch) = (nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch)-NXiKnotOrder_VolLockSmooth(IPatch)) * &
              !    (mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch)-NEtaKnotOrder_VolLockSmooth(IPatch)) * &
              !    (oo_NURBS_NumberOfUnivariateZetaKnots_VolLockSmooth(IPatch)-NZetaKnotOrder_VolLockSmooth(IPatch))
                  !number of elements -> note 2D implementation = 2 elements in the example 
              
              !nnp_NURBS_VolLockSmooth(IPatch) = nn_NURBS_NumberOfUnivariateXiKnots_VolLockSmooth(IPatch)*mm_NURBS_NumberOfUnivariateEtaKnots_VolLockSmooth(IPatch)*oo_NURBS_NumberOfUnivariateZetaKnots_VolLockSmooth(IPatch) 
              !!number of global basis functions (global here refers to its global domain within the 'super' element)
              !
              !nen_NURBS_VolLockSmooth(IPatch) = (NXiKnotOrder_VolLockSmooth(IPatch)+1) * (NEtaKnotOrder_VolLockSmooth(IPatch)+1) * (NZetaKnotOrder_VolLockSmooth(IPatch)+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
              
                
          end do     
          end if 
          
          
          ! --> This is repeated twice
          ! MULTIPATCH generlization aspect
                ! I need to know the maximum number of entries in xi, eta, and zeta, so that I can assign
                ! the size of XiKnotEntries, EtaKnotEntries, ZetaKnotEntries
                NURBS%Maximum_NControlPoints = maxval(Counters%NodTot) !NURBS%NControlPoints) ! <- finding maximum value
    
                NURBS%Maximum_nel_NURBS = maxval(Counters%NEl) !NURBS%nel_NURBS)
                NURBS%Maximum_nnp_NURBS = maxval(NURBS%nnp_NURBS)
                NURBS%Maximum_nen_NURBS = maxval(NURBS%nen_NURBS)
                
                
                NURBS%Maximum_nel_NURBS_Traction = maxval(NURBS%nel_NURBS_Traction)
                NURBS%Maximum_nnp_NURBS_Traction = maxval(NURBS%nnp_NURBS_Traction)
                NURBS%Maximum_nen_NURBS_Traction = maxval(NURBS%nen_NURBS_Traction)
                
                
                !Maximum_NFixities = maxval(NFixities)
                NURBS%Maximum_NElemMPPrescribedVelocity = maxval(NURBS%NElemMPPrescribedVelocity) ! --> we need to use this to allocate the prescribed velocity array
                NURBS%Maximum_NNodePrescribedVelocity = maxval(NURBS%NNodePrescribedVelocity)
                
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          ! these should probably be maximum_nnp_NURBS
          allocate(NURBS%INN(NURBS%Maximum_nnp_NURBS, NVECTOR, Counters%NPatches), stat=IError) ! INN has the size of number of control points(or global basis functions x NDIM )
          allocate(ElementConnectivities(NURBS%Maximum_nen_NURBS, NURBS%Maximum_nel_NURBS, Counters%NPatches), stat=IError)  ! IEN has the size of number of local basis functions x NDIM
          NURBS%INN = 0 !NURBS coordinate array (also called INC)
          ElementConnectivities = 0 !connectivity array

          
          
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          
          ! It is a surface in a 3D simulation so we use NVECTOR-1
          allocate(NURBS%INN_Traction(NURBS%Maximum_nnp_NURBS_Traction, & !(ILoadSystem, IPatch_Temporary) _Temporary
              NVECTOR-1, &
              NURBS%Maximum_NLoadSystems_Traction, &
              Counters%NPatches), stat=IError) ! INN has the size of number of control points(or global basis functions x NDIM )
                                    !(IPatch)              , Counters%NPatches 
          NURBS%INN_Traction = 0 !NURBS coordinate array (also called INC)
          
          allocate(NURBS%ElementConnectivities_Traction(NURBS%Maximum_nen_NURBS_Traction, & !(ILoadSystem_Temporary, IPatch_Temporary)
              NURBS%Maximum_nel_NURBS_Traction, & !(ILoadSystem_Temporary, IPatch_Temporary)
              NURBS%Maximum_NLoadSystems_Traction, &
              Counters%NPatches ), stat=IError)  ! IEN has the size of number of local basis functions x NDIM
                                    !(IPatch)       (IPatch)          , Counters%NPatches  
          NURBS%ElementConnectivities_Traction = 0 !connectivity array
          
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          ! these should probably be maximum_nnp_NURBS 
          allocate(NURBS%INN_VolLockSmooth(NURBS%nnp_NURBS_VolLockSmooth(IPatch_Temporary), NVECTOR, Counters%NPatches), stat=IError) ! INN has the size of number of control points(or global basis functions x NDIM )
          allocate(NURBS%ElementConnectivities_VolLockSmooth(NURBS%nen_NURBS_VolLockSmooth(IPatch_Temporary), NURBS%nel_NURBS_VolLockSmooth(IPatch_Temporary), Counters%NPatches), stat=IError)  ! IEN has the size of number of local basis functions x NDIM
          NURBS%INN_VolLockSmooth = 0 !NURBS coordinate array (also called INC)
          NURBS%ElementConnectivities_VolLockSmooth = 0 !connectivity array
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          
          do IPatch = 1, Counters%NPatches ! loop over patches
              
                  
                if (NDIM ==2) then 
                    
                    do ILoadSystem = 1, NURBS%NLoadSystems_Traction(IPatch) ! loop over load systems in each patch
                    
                        ! NEED SOME IF STATEMENTS HERE
                        
                        
                    !1D --> For LINE traction implementation
                    
                        if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
                        
                        call Build_INC_IEN_Array_1D_Traction_XI(ILoadSystem, IPatch)
                    
                        else if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
                            
                        call Build_INC_IEN_Array_1D_Traction_ETA(ILoadSystem, IPatch)
                        
                        end if 
                        
                    
                    end do 
                    
                    !2D --> For SURFACE implementation
                    call Build_INC_IEN_Array(IPatch) ! IEN and INN arrays for NUBS implementation
                    
                    !2D --> For SURFACE implementation --> reduced order volumetric locking implementation
                    call Build_INC_IEN_Array_VolLockSmooth(IPatch) ! IEN and INN arrays for NUBS implementation
                    
                else if (NDIM ==3) then 
                    
                    do ILoadSystem = 1, NURBS%NLoadSystems_Traction(IPatch) ! loop over load systems in each patch
                    
                        ! NEED SOME IF STATEMENTS HERE
                        
                        !2D --> For SURFACE traction implementation 
                        
                        if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
                        
                    call Build_INC_IEN_Array_2D_Traction_XI_ETA(ILoadSystem, IPatch)
                    !2D --> For SURFACE traction implementation 
                    
                    else if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
                    
                    call Build_INC_IEN_Array_2D_Traction_XI_ZETA(ILoadSystem, IPatch)
                    
                    
                    else if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 3) then 
                    
                    !2D --> For SURFACE traction implementation 
                    call Build_INC_IEN_Array_2D_Traction_ETA_ZETA(ILoadSystem, IPatch)
                    
                    end if 
                    
                    end do 
                    
                    !3D --> For SOLID implemetation
                    call Build_INC_IEN_Array_3D(IPatch)
                    
                    !3D --> For SOLID implemetation --> reduced order volumetric locking implementation
                    call Build_INC_IEN_Array_3D_VolLockSmooth(IPatch)
                    
                end if 
          
                !end do ! loop over load systems
          end do ! loop over patches
          
                
                
                
                
                
                
                
                
                
                
              ! overwrite the number of elements by nel_NURBS --> this should be calculated by the code and not given as an input
              !NElements = nel_NURBS(IPatch_Temporary)
              
              
              !note that if NURBS, we need to back calculate how many elements we have from the knot vectors
              
              ! allocating allocatable arrays
              allocate(NodalCoordinates(NURBS%Maximum_NControlPoints, NVECTOR, Counters%NPatches), stat = IError) !size equal to no of nodes, no of dimensions 
              allocate(NodalOriginalCoord(NURBS%Maximum_NControlPoints, NVECTOR, Counters%NPatches), stat = IError) !size equal to no of nodes, no of dimensions 
              allocate(NodalCoordinatesUpd(NURBS%Maximum_NControlPoints, NVECTOR, Counters%NPatches), stat = IError) !size equal to no of nodes, no of dimensions 
              allocate(ControlPoint_Weights(NURBS%Maximum_NControlPoints, Counters%NPatches), stat = IError) ! size equal to no of nodes
              NodalCoordinates = 0.0
              NodalOriginalCoord = 0.0
              NodalCoordinatesUpd = 0.0
              
              ! Volumetric locking smoothing 
              allocate(NodalCoordinates_VolLockSmooth(NURBS%NControlPoints_VolLockSmooth(IPatch_Temporary), NVECTOR, Counters%NPatches), stat = IError) !size equal to no of nodes, no of dimensions 
              allocate(ControlPoint_Weights_VolLockSmooth(NURBS%NControlPoints_VolLockSmooth(IPatch_Temporary), Counters%NPatches), stat = IError) ! size equal to no of nodes
              NodalCoordinates_VolLockSmooth = 0.0
              ControlPoint_Weights_VolLockSmooth = 0.0
              
              ! NURBS Traction 
              !allocate(NodalCoordinates_Traction(NControlPoints_Traction(ILoadSystem_Temporary, IPatch_Temporary), NVECTOR, ILoadSystem_Temporary, IPatch_Temporary), stat = IError) !size equal to no of nodes, no of dimensions 
              allocate(NodalCoordinates_Traction(NURBS%Maximum_NControlPoints_Traction, NVECTOR, NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError ) !size equal to no of nodes, no of dimensions  ! (ILoadSystem_Temporary, IPatch_Temporary)
              !allocate(ControlPoint_Weights_Traction(NControlPoints_Traction(ILoadSystem_Temporary, IPatch_Temporary), ILoadSystem_Temporary, IPatch_Temporary), stat = IError) ! size equal to no of nodes
              allocate(ControlPoint_Weights_Traction(NURBS%Maximum_NControlPoints_Traction, NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError ) ! size equal to no of nodes
              !, Counters%NPatches , Counters%NPatches
              NodalCoordinates_Traction = 0.0
              ControlPoint_Weights_Traction = 0.0
              
              ! NURBS Traction - system B
              !allocate(NodalCoordinates_Traction_B(NControlPoints_Traction_B, NVECTOR), stat = IError) !size equal to no of nodes, no of dimensions 
              !allocate(ControlPoint_Weights_Traction_B(NControlPoints_Traction_B), stat = IError) ! size equal to no of nodes
              !!, Counters%NPatches , Counters%NPatches
              !NodalCoordinates_Traction_B = 0.0
              !ControlPoint_Weights_Traction_B = 0.0
              
              
              ! we assume ELEMENTNODES is the same for all patches
              !allocate(ElementConnectivities(ELEMENTNODES, NURBS%Maximum_nel_NURBS, Counters%NPatches), stat = IError) !size equal to no of elements, no of dimensions ! -> NURBS implementation
              !ElementConnectivities = 0
              
              if (ELEMENTTYPE == TETRAOLD) then ! elementtype before v2018.2  
                allocate(ElementConnectivities10Node(N_NODES_HOE, NURBS%Maximum_nel_NURBS), stat = IError)
                ElementConnectivities10Node = 0
              end if  
              
              allocate(NodalPrescibedDisp(NURBS%Maximum_NControlPoints, NVECTOR, Counters%NPatches), stat = IError) !size equal to no of nodes, no of dimensions 
              !allocate(NodalPrescibedDisp_VolLockSmooth(NControlPoints_VolLockSmooth(IPatch_Temporary), NVECTOR, Counters%NPatches), stat = IError) !size equal to no of nodes, no of dimensions 
              
              
              allocate(NodalPrescibedDispWater(NURBS%Maximum_NControlPoints, NVECTOR, Counters%NPatches), stat = IError) !size equal to no of nodes, no of dimensions 
              allocate(NodalPrescibedDispGas(NURBS%Maximum_NControlPoints, NVECTOR, Counters%NPatches), stat = IError) !size equal to no of nodes, no of dimensions 
              
              NodalPrescibedDisp = 1.d10
              !NodalPrescibedDisp_VolLockSmooth = 1.d10
              NodalPrescibedDispWater = 1.d10          
              NodalPrescibedDispGas = 1.d10 
              
              
              ! Rank 1 index on ElementMaterialID is local element number
              allocate (ElementMaterialID(NURBS%Maximum_nel_NURBS, Counters%NPatches), stat = IError) !size equal to no of elements
              
              
              ! calculate the total number of elements in the system of all patches
              !NURBS%Sum_nel_NURBS = 0
              Counters%Sum_NEl = 0
              do ii = 1, Counters%NPatches
                  !NURBS%Sum_nel_NURBS = NURBS%Sum_nel_NURBS + Counters%NEl(ii) !NURBS%nel_NURBS(ii)
                  Counters%Sum_NEl = Counters%Sum_NEl + Counters%NEl(ii) 
              end do 

              ! calculate the total number of nodes in the system of all patches
              !NURBS%Sum_NControlPoints = 0
              Counters%Sum_NodTot = 0
              do ii = 1, Counters%NPatches
                  !NURBS%Sum_NControlPoints = NURBS%Sum_NControlPoints + Counters%NodTot(ii) !NURBS%NControlPoints(ii)
                  Counters%Sum_NodTot = Counters%Sum_NodTot + Counters%NodTot(ii)  
              end do 
              
              ! allocating ActiveElement based on the total number of elements
              allocate(IsActiveElement(NURBS%Maximum_nel_NURBS,Counters%NPatches), stat = IError) !size equal to no of elements !Maximum_nel_NURBS
              
              ! Active node is based on the global numbering of control points excluding 
              ! I am moving the allocation of ActiveNode to after we read all the control points and we realise which ones are common
              allocate(ActiveNode(NURBS%Maximum_NControlPoints, Counters%NPatches), stat = IError) !size equal to no of nodes 
              
              allocate(HydraulicHeadLoadedElemID(NURBS%Maximum_nel_NURBS), stat = IError) !size equal to no of elements
              ElementMaterialID = 0
              IsActiveElement = .false.
              ActiveNode = .false.
              HydraulicHeadLoadedElemID = .false.
                
          
              ! this serves as a documentation so that we can easily know what is the 
              ! global control point index from inputting the local control point index and its corresponding patch number
              !allocate(Multipatch_Connecting_Local_To_Global_ControlPoints(Counters%NPatches, Maximum_NControlPoints) )
              allocate(Multipatch_Connecting_Local_To_Global_ControlPoints(NURBS%Maximum_NControlPoints , Counters%NPatches ))

              Multipatch_Connecting_Local_To_Global_ControlPoints = 0
          
              
              ! Finding maximum order 
              NURBS%Maximum_NXiKnotOrder = maxval(NURBS%NXiKnotOrder)
              NURBS%Maximum_NEtaKnotOrder = maxval(NURBS%NEtaKnotOrder)
              NURBS%Maximum_NZetaKnotOrder = maxval(NURBS%NZetaKnotOrder)
          
          
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

          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          
          SizeAB = 1 + 3 * NDIM ! array size for reading absorbing boundary data
		            
          ! set GOM version number
            select case (FileVersion)
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
                
                
            ! Multipatch generalization allocations performed here ------------------------------------------    
                
                
                
                
             ! This needs to go to a different subroutine   
             !Moving mesh 
             !allocate(CalParams%MovingMesh%MeshAreas(2, 3, NDIM * (NDIM - 1) + 2, Counters%NPatches))   
             !allocate(CalParams%MovingMesh%NAreaNodes(Counters%NPatches))
             !allocate(CalParams%MovingMesh%NMovingMeshDirections(Counters%NPatches))
             !allocate(CalParams%MovingMesh%MovingMeshDirection(Counters%NPatches))
             !allocate(CalParams%MovingMesh%NStorageAreas(2,Counters%NPatches))
             !allocate(CalParams%MovingMesh%StructureMaterialID(Counters%NPatches))
             !allocate(CalParams%MovingMesh%MovingMaterialID(Counters%NPatches))
             !allocate(CalParams%MovingMesh%NStructureNodes(Counters%NPatches))
             !
             !
             !CalParams%MovingMesh%MeshAreas = 0
             !CalParams%MovingMesh%NAreaNodes = 0
             !CalParams%MovingMesh%NMovingMeshDirections = 0
             !CalParams%MovingMesh%MovingMeshDirection = 0
             !CalParams%MovingMesh%NStorageAreas = 0
             !CalParams%MovingMesh%StructureMaterialID = 0
             !CalParams%MovingMesh%MovingMaterialID = 0
             !CalParams%MovingMesh%NStructureNodes = 0
             

                
                !Moving mesh initialization 
                if (CalParams%ApplyMovingMesh) then 
                    call InitialiseMovingMeshData()
                end if
             
                
                if (CalParams%PrescribedVelo%ApplyPrescribedVelo) then 
                    !Prescribed velocity initialization ---------> Material point 
                    call InitialisePrescribedMPVeloData() ! --> need to add a switch for this

                    !Prescribed velocity initialization ---------> Node
                    call InitialisePrescribedNodalVeloData() ! --> need to add a switch for this
                
                end if
                
             
            !------------------------------------------------------------------------------------------------

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
              else if (trim(BName) == '$$PATCH_ID') then !--> multipatch generalized 
                  ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                  IPatch = -1
                  ! read patch ID from the GOM file
                  read(FileUnit, *) IPatch
              
              else if (trim(BName)=='$$PRESCRIBED_MATERIAL_POINT_VELOCITY_VOLUME') then !prescribed velcoity is assigned
                if (NURBS%NElemMPPrescribedVelocity(IPatch)>0) then
                  GeoParams%ApplyMPVelocityVolume=.true.
                call ReadPrescribedVelocityMaterialPoint(FileUnit, IPatch)
                end if 
         
              else if (trim(BName)=='$$PRESCRIBED_MATERIAL_POINT_VELOCITY_SURFACE') then !prescribed velcoity is assigned
                if (NURBS%NElemMPPrescribedVelocity(IPatch)>0) then  
                  GeoParams%ApplyMPVelocitySurface=.true. ! turn switch on
                call ReadPrescribedVelocityMaterialPoint(FileUnit, IPatch)
                end if 
                
              else if (trim(BName)=='$$PRESCRIBED_NODAL_VELOCITY_VOLUME') then !prescribed velcoity is assigned
                  if (NURBS%NNodePrescribedVelocity(IPatch)>0) then 
                  GeoParams%ApplyNodalVelocityVolume=.true.
                call ReadPrescribedVelocityNode(FileUnit, IPatch)
                  end if
                
              else if (trim(BName)=='$$PRESCRIBED_NODAL_VELOCITY_SURFACE') then !prescribed velcoity is assigned
                if (NURBS%NNodePrescribedVelocity(IPatch)>0) then 
                  GeoParams%ApplyNodalVelocitySurface=.true.
                call ReadPrescribedVelocityNode(FileUnit, IPatch)
                end if 
                
              else if (trim(BName)=='$$PRESCRIBED_NODAL_VELOCITY_LINE') then !prescribed velcoity is assigned
                  if (NURBS%NNodePrescribedVelocity(IPatch)>0) then 
                GeoParams%ApplyNodalVelocityLine=.true.
                call ReadPrescribedVelocityNode(FileUnit, IPatch)
                  end if 
                  
                
              else if (trim(BName)=='$$PRESCRIBED_NODAL_VELOCITY_POINT') then !prescribed velcoity is assigned
                  if (NURBS%NNodePrescribedVelocity(IPatch)>0) then 
                GeoParams%ApplyNodalVelocityPoint=.true.
                call ReadPrescribedVelocityNode(FileUnit, IPatch)
                end if
				
		      !!! DATA FOR INITIAL VELOCITY ON MP
			  !--> not multipatch generalized
              else if (trim(BName)=='$$INITIAL_VELOCITY_MATERIAL_POINT') then ! Initial velocity on MP data
                  read(FileUnit,*,iostat=ios) DumI(1)
				  call Assert( ios == 0, messageIOS // trim(BName) )
				  call Assert( DumI(1) > 0, 'GOM-ERROR: The number of elements to initialise velocity on material points must be greater than 0!' )
				  GeoParams%NMovingElements=DumI(1)
				  allocate(GeoParams%InitialVelocityonMP(DumI(1),NVECTOR+1), stat = ios)
				  call Assert( ios == 0, messageIOS // trim(BName) )				  
				  do I = 1, GeoParams%NMovingElements(IPatch)
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
                  read(FileUnit,*,iostat=ios) GeoParams%WaterSurfaceMaterialID(I), GeoParams%WaterSurfaceFileName(I)                 
                  end do
                  
              !!! DATA FOR MOVING MESH  
              else if (trim(BName) == '$$PATCH_ID') then 
                  ! before we start collecting data for the moving mesh, we need to define the patch ID so that 
                  ! we store it for the right patch
                  ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                  IPatch = -1
                  ! read patch ID from the GOM file
                  read(FileUnit, *) IPatch 
                  
                  !----------- define corners
              else if (trim(BName)=='$$EXTENDING_MESH_CORNER_NODES') then ! moving mesh data
                ExtendingMeshPatchID = IPatch
               do I = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                  read(FileUnit,*,iostat=ios) DumI(1)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                  call Assert( DumI(1) > 0, 'GOM-ERROR: all node numbers of $$EXTENDING_MESH_CORNER_NODES must be positive!' )
                  GeoParams%MovingMeshExtendingCorners(I, IPatch) = DumI(1)
                end do   
              
              else if (trim(BName)=='$$COMPRESSING_MESH_CORNER_NODES') then ! moving mesh data
                  CompressingMeshPatchID = IPatch
                do I = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                  read(FileUnit,*,iostat=ios) DumI(1)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                  call Assert( DumI(1) > 0, 'GOM-ERROR: all node numbers of $$COMPRESSING_MESH_CORNER_NODES must be positive!' )
                  GeoParams%MovingMeshCompressingCorners(I, IPatch) = DumI(1)
                end do   
              
              else if (trim(BName)=='$$MOVING_MESH_CORNER_NODES') then ! moving mesh data
                  MovingMeshPatchID = IPatch
                  do I = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                  read(FileUnit,*,iostat=ios) DumI(1)
                  call Assert( ios == 0, messageIOS // trim(BName) // '(' // trim(String(I)) // ')' )
                  call Assert( DumI(1) > 0, 'GOM-ERROR: all node numbers of $$MOVING_MESH_CORNER_NODES must be positive!' )
                  GeoParams%MovingMeshMovingCorners(I, IPatch) = DumI(1)
                end do   
                !----------- define corners
              
                !----------- define materials
              else if (trim(BName)=='$$MOVING_MESH_DIRECTION') then ! moving mesh data
                read(FileUnit,*,iostat=ios) DumS
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( trim(DumS) == 'x-direction' .or. trim(DumS) == 'y-direction' .or. trim(DumS) == 'z-direction', &
                             'GOM-ERROR: $$MOVING_MESH_DIRECTION must be "x-direction", "y-direction" or "z-direction" (only for 3D)!' )
                GeoParams%MovingMeshDirection(IPatch) = DumS
                
              else if (trim(BName)=='$$MOVING_MESH_REFERENCE_MATERIAL_INDEX') then 
                read(FileUnit,*,iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS // trim(BName) )
                call Assert( DumI(1) > 0 , 'GOM-ERROR: ' //trim(BName)// ' must be larger than 0.' )
                GeoParams%MovingMeshReferenceMaterialID(IPatch) = DumI(1)
                !----------- define materials

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
              else if (trim(BName) == '$$PATCH_ID') then ! this flag is to generalize the local elements to certain elements in a patch
                  ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                  IPatch = -1
                  ! read patch ID from the GOM file
                  read(FileUnit, *) IPatch 
                  
              else if (trim(BName)=='$$STARTDAMPING') then
                do I = 1, Counters%NEl(IPatch) !NURBS%nel_NURBS(IPatch)!Counters%NEl ! loop over elements !
                  read(FileUnit,*,iostat=ios) DumR(1)
                  call Assert( ios == 0, messageIOS // trim(BName) )
                  if (DumR(1) > 0.0) GeoParams%ApplyLocalDampingElement = .true. 
                  GeoParams%LocalDampingFactorElement(I, IPatch) = DumR(1) ! read the damping value for each element
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
            allocate(GeoParams%ExcavatedElements(Counters%Sum_NEl), stat = IError) !Counters%NEl !NURBS%Sum_nel_NURBS
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
          allocate(GeoParams%WaterSurfaceFileName(ArraySize), stat=IError)
          
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


        subroutine InitialisePrescribedMPVeloData()!Nnode)!, PatchID)
        !**********************************************************************
        !
        ! Function:  To initialise the arrays relate to prescribed velocity
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError, Nnode
          
          allocate( GeoParams%PrescribedVeloNElem(Counters%NPatches) )!, stat = IError )
          GeoParams%PrescribedVeloNElem = 0  
          allocate( GeoParams%PrescribedVeloElemID(NURBS%Maximum_NElemMPPrescribedVelocity, Counters%NPatches) )!, stat = IError) !, Counters%NPatches !maximum_nel_NURBS  
          GeoParams%PrescribedVeloElemID = -1
          allocate( GeoParams%PrescribedVeloElValue(NURBS%Maximum_NElemMPPrescribedVelocity, NVECTOR, Counters%NPatches) )!, stat = IError) !, Counters%NPatches !maximum_nel_NURBS
          GeoParams%PrescribedVeloElValue = 0.0
          allocate( GeoParams%PrescribedVeloElDirection(NURBS%Maximum_NElemMPPrescribedVelocity, NVECTOR, Counters%NPatches) )!, stat = IError) !, Counters%NPatches !maximum_nel_NURBS
          GeoParams%PrescribedVeloElDirection = 1
          
          ! commenting this out for this branch ------------------------------------
          ! Multipatch generalization 
          !integer(INTEGER_TYPE), intent(in) :: PatchID
          !allocate( GeoParams%PrescribedVeloElemID(Nnode), stat = IError) !, Counters%NPatches
          !   GeoParams%PrescribedVeloElemID = -1
          !allocate( GeoParams%PrescribedVeloElValue(Nnode, NVECTOR), stat = IError) !, Counters%NPatches
          !   GeoParams%PrescribedVeloElValue = 0.0
          !allocate( GeoParams%PrescribedVeloElDirection(Nnode, NVECTOR), stat = IError) !, Counters%NPatches
          !   GeoParams%PrescribedVeloElDirection = 1
          ! I changed the Nnode to maximum_nel_NURBS so that we can have a more generalized spots for all the element in the patches
          !allocate( GeoParams%PrescribedVeloElemID(maximum_nel_NURBS, Counters%NPatches), stat = IError) !, Counters%NPatches
          !   GeoParams%PrescribedVeloElemID = -1
          !allocate( GeoParams%PrescribedVeloElValue(maximum_nel_NURBS, NVECTOR, Counters%NPatches), stat = IError) !, Counters%NPatches
          !   GeoParams%PrescribedVeloElValue = 0.0
          !allocate( GeoParams%PrescribedVeloElDirection(maximum_nel_NURBS, NVECTOR, Counters%NPatches), stat = IError) !, Counters%NPatches
          !   GeoParams%PrescribedVeloElDirection = 1
          ! I changed the Nnode to maximum_nel_NURBS so that we can have a more generalized spots for all the element in the patches
          ! I changed this to the maximum number of elements where we are prescribing velocity on in a patch
          ! commenting this out for this branch ------------------------------------

        end subroutine InitialisePrescribedMPVeloData
        
        
        subroutine InitialiseMovingMeshData()!Nnode)!, PatchID)
        !**********************************************************************
        !
        ! Function:  To initialise the arrays relate to prescribed velocity
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          !integer(INTEGER_TYPE) :: IError, Nnode
          
          !Moving mesh 
 
          allocate(CalParams%MovingMesh%MeshAreas(2, 3, NDIM * (NDIM - 1) + 2, Counters%NPatches))   
          allocate(CalParams%MovingMesh%NAreaNodes(Counters%NPatches))
          allocate(CalParams%MovingMesh%NMovingMeshDirections(Counters%NPatches))
          allocate(CalParams%MovingMesh%MovingMeshDirection(Counters%NPatches))
          allocate(CalParams%MovingMesh%NStorageAreas(2,Counters%NPatches))
          allocate(CalParams%MovingMesh%StructureMaterialID(Counters%NPatches))
          allocate(CalParams%MovingMesh%MovingMaterialID(Counters%NPatches))
          allocate(CalParams%MovingMesh%NStructureNodes(Counters%NPatches))
 
          CalParams%MovingMesh%MeshAreas = 0
          CalParams%MovingMesh%NAreaNodes = 0
          CalParams%MovingMesh%NMovingMeshDirections = 0
          CalParams%MovingMesh%MovingMeshDirection = 0
          CalParams%MovingMesh%NStorageAreas = 0
          CalParams%MovingMesh%StructureMaterialID = 0
          CalParams%MovingMesh%MovingMaterialID = 0
          CalParams%MovingMesh%NStructureNodes = 0

        end subroutine InitialiseMovingMeshData
        
        
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


        subroutine InitialisePrescribedNodalVeloData()!Nnode)
        !**********************************************************************
        !
        ! Function:  To initialise the arrays relate to prescribed velocity
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError, Nnode

          allocate( CalParams%PrescribedVelo%NNodePrescribedVelo(Counters%NPatches) )
          CalParams%PrescribedVelo%NNodePrescribedVelo = 0 ! --> need multipatch generalization 
          
          !Prescribed velocity ---------> Node
          allocate( CalParams%PrescribedVelo%NodePrescribedVelo(NURBS%Maximum_NNodePrescribedVelocity, Counters%NPatches), stat = IError)
          CalParams%PrescribedVelo%NodePrescribedVelo = -1
          allocate( CalParams%PrescribedVelo%NodalPrescribedVelocityValue(NURBS%Maximum_NNodePrescribedVelocity, NVECTOR, Counters%NPatches), stat = IError)
          CalParams%PrescribedVelo%NodalPrescribedVelocityValue = 0.0
          allocate( CalParams%PrescribedVelo%NodalPrescribedVelocityDirection(NURBS%Maximum_NNodePrescribedVelocity, NVECTOR, Counters%NPatches), stat = IError)
          CalParams%PrescribedVelo%NodalPrescribedVelocityDirection = 1
          
          
          !allocate( CalParams%PrescribedVelo%NodePrescribedVelo(Nnode), stat = IError)
          !CalParams%PrescribedVelo%NodePrescribedVelo = -1
          !allocate( CalParams%PrescribedVelo%NodalPrescribedVelocityValue(Nnode, NVECTOR), stat = IError)
          !CalParams%PrescribedVelo%NodalPrescribedVelocityValue = 0.0
          !allocate( CalParams%PrescribedVelo%NodalPrescribedVelocityDirection(Nnode, NVECTOR), stat = IError)
          !CalParams%PrescribedVelo%NodalPrescribedVelocityDirection = 1

        end subroutine InitialisePrescribedNodalVeloData
        
        
        subroutine ReadPrescribedVelocityMaterialPoint(FileUnit, IPatch)!PatchID)
        !**********************************************************************
        !
        ! Function:  To read data for prescribed velocity in material points
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          
          ! Multipatch generalization 
          integer(INTEGER_TYPE), intent(in) :: IPatch
          !integer(INTEGER_TYPE) :: maximum_PrescribedVeloNElem  
          
          ! local variables
          integer(INTEGER_TYPE) :: I
          integer(INTEGER_TYPE), dimension(NVECTOR+1) :: DumI
          real(REAL_TYPE), dimension(NVECTOR) :: DumR
          
          
          GeoParams%PrescribedVeloNElem(IPatch) = NURBS%NElemMPPrescribedVelocity(IPatch) !NNodePrescribedVelocity!DumI(1) !(PatchID)
                
          ! this was multipatch generalized for so that we can store the following for each patch
          do I = 1, GeoParams%PrescribedVeloNElem(IPatch) ! loop over number of prescribed velocity elements in a patch
              read(FileUnit,*) DumI(1), DumI(2:NVECTOR+1), DumR(1:NVECTOR)
              GeoParams%PrescribedVeloElemID(I, IPatch) = DumI(1) ! store the element ID
              GeoParams%PrescribedVeloElDirection(I, 1:NVECTOR, IPatch) = DumI(2:NVECTOR+1) ! store the switch vector
              GeoParams%PrescribedVeloElValue(I, 1:NVECTOR, IPatch) = DumR(1:NVECTOR) ! store the prescribed velocity value for each direction
          end do 
          
          
          ! commenting this out for this branch ------------------------------------
          !call DestroyPrescribedMPVeloData() ! --> I commented this for now because we are allocating this at the beginning
                
          !read(FileUnit,*) DumI(1) ! read integer 
          
          ! loop over variables and determine the maximum number of elements 
          !maximum_PrescribedVeloNElem = maxval(GeoParams%PrescribedVeloNElem)
          
          ! this needs generalization because we need to know the allocations 
          ! --> maybe we can have transfer variables for this
          !call InitialisePrescribedMPVeloData(GeoParams%PrescribedVeloNElem) !, PatchID    !--> I MOVED THIS ELSE WHERE OUTSIDE THE LOOP
          
          ! Example: Flag
          !$$PRESCRIBED_MATERIAL_POINT_VELOCITY_SURFACE
          !1 -> removed this
          !1  1   0   0.0001   0.0
          !ElementID, ElementDirection(NDIM) --> this is a switch vector, Value of prescribed velocity--> this is the value of the prescribed velocity   
          ! commenting this out for this branch ------------------------------------
          
        end subroutine ReadPrescribedVelocityMaterialPoint
        
        
        subroutine ReadPrescribedVelocityNode(FileUnit, IPatch)
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
                
          ! multipatch generalization 
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
          !call DestroyPrescribedNodalVeloData()   
          
          !read(FileUnit,*) DumI(1) 
          !CalParams%PrescribedVelo%NNodePrescribedVelo(IPatch) = DumI(1)
          
          CalParams%PrescribedVelo%NNodePrescribedVelo(IPatch) = NURBS%NNodePrescribedVelocity(IPatch)
          
                
          ! commenting this out for this branch ------------------------------------
          !call InitialisePrescribedNodalVeloData(CalParams%PrescribedVelo%NNodePrescribedVelo)  ! --> I commented this out for a second
          ! commenting this out for this branch ------------------------------------
          
          do I = 1, CalParams%PrescribedVelo%NNodePrescribedVelo(IPatch)
            read(FileUnit,*) DumI(1), DumI(2:NVECTOR+1), DumR(1:NVECTOR)
            CalParams%PrescribedVelo%NodePrescribedVelo(I, IPatch) = DumI(1)
            CalParams%PrescribedVelo%NodalPrescribedVelocityDirection(I, 1:NVECTOR, IPatch) = DumI(2:NVECTOR+1)
            CalParams%PrescribedVelo%NodalPrescribedVelocityValue(I, 1:NVECTOR, IPatch) = DumR(1:NVECTOR)
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
          logical, dimension(NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, Counters%NPatches) :: IsMovingAndCompressingCorner, IsMovingAndExtendingCorner
        
          ! Multipatch generalization 
          integer(INTEGER_TYPE) :: IPatch
          
          
          ! Loop over patches -----------------------------------------------------------------------------
          
          if (CalParams%ApplyMovingMesh == .true.) then 
              
          do IPatch = 1, Counters%NPatches ! loop over patches
          
          
          ! set flag for applying moving mesh algorithm         
          CalParams%ApplyMeshSmoothing = .false. !(IPatch)
          
          if ( GeoParams%MovingMeshExtendingCorners(1, IPatch) /= -1 .or. &
               GeoParams%MovingMeshCompressingCorners(1, IPatch) /= -1 .or. &
               GeoParams%MovingMeshMovingCorners(1, IPatch) /= -1) CalParams%ApplyMeshSmoothing = .true. !(IPatch)

          if ( .not. CalParams%ApplyMeshSmoothing ) CYCLE !RETURN !(IPatch)
          
          ! set basic parameters for moving mesh algorithm
          CalParams%MovingMesh%NMovingMeshDirections(IPatch) = 1
          CalParams%MovingMesh%NAreaNodes(IPatch) = NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh ! 4 for 2D (square), 8 for 3D (cube), 6 for cylindrical 3D

          ! check if extending and/or compressing mesh areas are defined and determine 'NStorageAreas'
          if ( CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_FILE) then ! .or. ELEMENTTYPE == QUAD4_NURBS  --> I commented this out
              ! no deforming mesh is defined for prescribed velocity from file
              CalParams%MovingMesh%NStorageAreas(:,IPatch) = 0
          else if ( GeoParams%MovingMeshExtendingCorners(1, IPatch) == -1 .and. GeoParams%MovingMeshCompressingCorners(1, IPatch) == -1 .and. Counters%NPatches==1) then 
            ! error: no deforming mesh is defined  
            call GiveError('Error: Moving Mesh: Unless you are prescribing velocity from file, at least one deforming mesh area has to be specified, i.e. either extending or compressing mesh.')  
          end if  
          if ( GeoParams%MovingMeshExtendingCorners(1, IPatch) /= -1 .and. GeoParams%MovingMeshCompressingCorners(1, IPatch) == -1 ) then
            ! only compressing mesh defined  
            CalParams%MovingMesh%NStorageAreas(:,IPatch) = 1
          end if  
          if ( GeoParams%MovingMeshExtendingCorners(1, IPatch) == -1 .and. GeoParams%MovingMeshCompressingCorners(1, IPatch) /= -1 ) then
            ! only extending mesh defined  
            CalParams%MovingMesh%NStorageAreas(:,IPatch) = 1
          end if  
          if ( GeoParams%MovingMeshExtendingCorners(1, IPatch) /= -1 .and. GeoParams%MovingMeshCompressingCorners(1, IPatch) /= -1 ) then
            ! both, extending and compressing mesh defined  
            CalParams%MovingMesh%NStorageAreas(:,IPatch) = 2
          end if  
          
          ! check if reference material is defined
          call Assert( GeoParams%MovingMeshReferenceMaterialID(IPatch) > 0 , 'Error: Moving Mesh: The reference material has to be defined.' )
          CalParams%MovingMesh%MovingMaterialID(IPatch) = GeoParams%MovingMeshReferenceMaterialID(IPatch)
          if (CalParams%ApplyContactAlgorithm) then 
			CalParams%MovingMesh%StructureMaterialID=CalParams%MovingMesh%MovingMaterialID !Sets StructureMaterialID required for contact and rigid body Algorithm
		  end if
          ! determine moving mesh direction
          if ( GeoParams%MovingMeshDirection(IPatch) == 'x-direction' ) direction = 1
          if ( GeoParams%MovingMeshDirection(IPatch) == 'y-direction' ) direction = 2
          if ( GeoParams%MovingMeshDirection(IPatch) == 'z-direction' ) direction = 3

          CalParams%MovingMesh%MovingMeshDirection(IPatch) = direction
          
          if ( GeoParams%MovingMeshCompressingCorners(1, IPatch) /= -1) then
            IsMovingAndCompressingCorner(:,IPatch) = .false. ! initialize
            do J = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
              do K = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                if (GeoParams%MovingMeshCompressingCorners(J, IPatch)==GeoParams%MovingMeshMovingCorners(K, IPatch)) then
                  IsMovingAndCompressingCorner(J, IPatch) = .true. ! assign compressing mesh 
                end if
              end do
            end do
          
            do J = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
              if (.not.IsMovingAndCompressingCorner(J,IPatch)) then
                FixedNodeCompressing = GeoParams%MovingMeshCompressingCorners(J, IPatch)
                PositionNodeCompressing = J
                EXIT
              end if
            end do 
          end if
          
          if ( GeoParams%MovingMeshExtendingCorners(1, IPatch) /= -1) then
            IsMovingAndExtendingCorner(:,IPatch) = .false.
            do J = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
              do K = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
                if (GeoParams%MovingMeshExtendingCorners(J, IPatch)==GeoParams%MovingMeshMovingCorners(K, IPatch)) then
                  IsMovingAndExtendingCorner(J, IPatch) = .true.
                end if
              end do
            end do
          
            do J = 1, NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh
              if (.not.IsMovingAndExtendingCorner(J, IPatch)) then
                FixedNodeExtending = GeoParams%MovingMeshExtendingCorners(J, IPatch)
                PositionNodeExtending = J
                EXIT
              end if
            end do 
          end if
          
          ! assign corner node numbers to array CalParams%MovingMesh%MeshAreas(k,j,i)
          K = CalParams%MovingMesh%NMovingMeshDirections(IPatch)
          
          if (CalParams%MovingMesh%NStorageAreas(K, IPatch) == 1) then

            ! extending mesh
            if ( GeoParams%MovingMeshExtendingCorners(1, IPatch) /= -1) then
              CalParams%MovingMesh%MeshAreas(K, 1, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, IPatch) = GeoParams%MovingMeshExtendingCorners(1:NDIM * (NDIM - 1) + 2, IPatch) ! fill whole array
              CalParams%MovingMesh%MeshAreas(K, 1, 1, IPatch) = GeoParams%MovingMeshExtendingCorners(PositionNodeExtending, IPatch) ! replace first node with the one on the fixed boundary
              CalParams%MovingMesh%MeshAreas(K, 1, PositionNodeExtending, IPatch) = GeoParams%MovingMeshExtendingCorners(1, IPatch) ! put the old first node on the empty position 
            end if
            ! compressing mesh
            if ( GeoParams%MovingMeshCompressingCorners(1, IPatch) /= -1) then 
              CalParams%MovingMesh%MeshAreas(K, 1, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, IPatch) = GeoParams%MovingMeshCompressingCorners(1:NDIM * (NDIM - 1) + 2, IPatch) ! fill whole array
              CalParams%MovingMesh%MeshAreas(K, 1, 1, IPatch) = GeoParams%MovingMeshCompressingCorners(PositionNodeCompressing, IPatch) ! replace first node with the one on the fixed boundary
              CalParams%MovingMesh%MeshAreas(K, 1, PositionNodeCompressing, IPatch) = GeoParams%MovingMeshCompressingCorners(1, IPatch) ! put the old first node on the empty position 
            end if

          else if (CalParams%MovingMesh%NStorageAreas(K, IPatch) == 2) then

            ! extending mesh
            CalParams%MovingMesh%MeshAreas(K, 1, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, IPatch) = GeoParams%MovingMeshExtendingCorners(1:NDIM * (NDIM - 1) + 2, IPatch) ! fill whole array
            CalParams%MovingMesh%MeshAreas(K, 1, 1, IPatch) = GeoParams%MovingMeshExtendingCorners(PositionNodeExtending, IPatch) ! replace first node with the one on the fixed boundary
            CalParams%MovingMesh%MeshAreas(K, 1, PositionNodeExtending, IPatch) = GeoParams%MovingMeshExtendingCorners(1, IPatch) ! put the old first node on the empty position 
            ! compressing mesh
            CalParams%MovingMesh%MeshAreas(K, 2, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, IPatch) = GeoParams%MovingMeshCompressingCorners(1:NDIM * (NDIM - 1) + 2, IPatch) ! fill whole array
            CalParams%MovingMesh%MeshAreas(K, 2, 1, IPatch) = GeoParams%MovingMeshCompressingCorners(PositionNodeCompressing, IPatch) ! replace first node with the one on the fixed boundary
            CalParams%MovingMesh%MeshAreas(K, 2, PositionNodeCompressing, IPatch) = GeoParams%MovingMeshCompressingCorners(1, IPatch) ! put the old first node on the empty position 
             
          end if
  
          ! moving mesh 
          CalParams%MovingMesh%MeshAreas(K, 3, 1:NDIM * (NDIM - 1) + GeoParams%ExtraNodesMovingMesh, IPatch) = GeoParams%MovingMeshMovingCorners(1:NDIM * (NDIM - 1) + 2, IPatch) ! fill whole array
            
          
          
          end do 
          ! Loop over patches -----------------------------------------------------------------------------

          
          end if ! if we are applying moving mesh
          
          
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