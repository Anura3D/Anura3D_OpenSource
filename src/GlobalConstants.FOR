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


	  module ModGlobalConstants
      !**********************************************************************
      !
      !    Function:  Contains definition of global constants.
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 10167 $
      !     $Date: 2024-06-25 15:15:53 +0200 (Tue, 25 Jun 2024) $
      !
      !**********************************************************************
      implicit none

        ! definition of single, double and quad precision
        integer, parameter :: SINGLE_PRECISION = selected_real_kind(6, 37)
        integer, parameter :: DOUBLE_PRECISION = selected_real_kind(15, 307) ! is used as standard precision
        integer, parameter :: QUAD_PRECISION = selected_real_kind(33, 4931)

        ! definition of variable precision for REAL, INTEGER
        integer, parameter :: REAL_TYPE = DOUBLE_PRECISION
        integer, parameter :: INTEGER_TYPE = 4
       
        ! analysis types
        character(len=15), parameter :: DIM_2D_PLANESTRAIN = '2D-plane_strain'
        character(len=15), parameter :: DIM_2D_AXISYMMETRIC = '2D-axisymmetric'
        character(len=12), parameter :: DIM_3D_CARTESIAN = '3D-cartesian'
        character(len=14), parameter :: DIM_3D_CYLINDRIC = '3D-cylindrical'
        
        ! formulation
        character(len=12), parameter :: SINGLE_POINT = 'single-point'
        character(len=12), parameter :: DOUBLE_POINT = 'double-point'
        
        ! dimension constants, the PROTECTED attribute prevents code outside the module from changing its value
        integer(INTEGER_TYPE), protected :: NDIM = -1 ! number of dimensions
        integer(INTEGER_TYPE), protected :: NDOFL = -1 ! number of degrees of freedom
        integer(INTEGER_TYPE), protected :: NDOFLQS = -1 ! number of degrees of freedom implicit quasi-static
        integer(INTEGER_TYPE), protected :: NTENSOR = -1 ! stress and strain tensor
        integer(INTEGER_TYPE), protected :: NVECTOR = -1 ! velocity, displacement
        integer(INTEGER_TYPE), protected :: NPRINCIPAL = -1 ! principal stress and strain vector
        logical, protected :: ISAXISYMMETRIC = .false. ! specifies if analysis is 2D axisymmetric
        logical, protected :: IS3DCYLINDRIC = .false. ! specifies if analysis is 3D with cylindric coordinates
        integer(INTEGER_TYPE), protected :: NFORMULATION = -1 ! single-point formulation = 1, double-point formulation = 2
        
        ! 'old' 3D constants
        integer(INTEGER_TYPE), parameter :: N_NODES_HOE = 10
        integer(INTEGER_TYPE), parameter :: N_BOUNDARY_NODES_HOE = 6
        
        ! element types
        character(len=18), parameter :: TRI3 = 'triangular_3-noded'
        character(len=18), parameter :: TRI6 = 'triangular_6-noded'
        character(len=21), parameter :: QUAD4 = 'quadrilateral_4-noded'
        character(len=21), parameter :: QUAD8 = 'quadrilateral_8-noded'
        character(len=19), parameter :: TETRA4 = 'tetrahedral_4-noded'
        character(len=20), parameter :: TETRA10 = 'tetrahedral_10-noded'
        character(len=18), parameter :: HEXA8 = 'hexahedral_8-noded'
        character(len=19), parameter :: HEXA20 = 'hexahedral_20-noded'
        character(len=19), parameter :: TETRAOLD = 'tetrahedral_old' ! this is the elementtype that was used before v2018.2
        
        ! General constants
        integer(INTEGER_TYPE), parameter :: MAX_EXTENSION_LENGTH = 20
        integer(INTEGER_TYPE), parameter :: MAX_FILENAME_LENGTH = 255
        integer(INTEGER_TYPE), parameter :: INITIAL_STEP = 1
        integer(INTEGER_TYPE), parameter :: MAXIMUM_STEP = 100000
        integer(INTEGER_TYPE), parameter :: ID_UNDEFINED = -1
        integer(INTEGER_TYPE), parameter :: NSTATEVAR = 50
        integer(INTEGER_TYPE), parameter :: NPROPERTIES = 50

        real(REAL_TYPE), dimension(:), allocatable, protected :: DEFAULT_GRAVITY_DIRECTION 
        real(REAL_TYPE), parameter :: DEFAULT_GRAVITY_ACCELERATION = 9.81
        real(REAL_TYPE), parameter :: SMALL = 1d-15
        real(REAL_TYPE), parameter :: TINY  = 1d-60
        real(REAL_TYPE), parameter :: LARGE = 1D10

        real(REAL_TYPE), parameter :: BULK_MODULUS_WATER = 2.0d6 ! kPa

        ! Mathematical constants
        real(REAL_TYPE), parameter :: PI = 3.141592653589793238462643383279502884197169399D0

        ! Computation-related constants
        integer(INTEGER_TYPE), parameter :: MPM_MIXED_INTEGRATION = 0
        integer(INTEGER_TYPE), parameter :: MPM_MP_INTEGRATION = 1
        integer(INTEGER_TYPE), parameter :: FEM = 2
        integer(INTEGER_TYPE), parameter :: UL_FEM = 3
        character(len = 9), parameter :: MPM_MIXED_INTEGRATION_SPECIFIER = 'MPM-MIXED'
        character(len = 6), parameter :: MPM_MP_INTEGRATION_SPECIFIER = 'MPM-MP'
        character(len = 3), parameter :: FEM_SPECIFIER = 'FEM'
        character(len = 6), parameter :: UL_FEM_SPECIFIER = 'UL-FEM'
        character(len = 6), parameter :: LOAD_LINEAR = 'LINEAR'
        character(len = 4), parameter :: LOAD_STEP = 'STEP'
        character(len = 6), parameter :: LOAD_TYPE_LINEAR = 'linear'
        character(len = 4), parameter :: LOAD_TYPE_STEP = 'step'
        character(len = 4), parameter :: LOAD_TYPE_FILE = 'file'
        character(len = 3), parameter :: LOAD_TYPE_OFF = 'off'
        integer(INTEGER_TYPE), parameter :: MAX_LOAD_SYSTEMS = 2


        ! Convergence-related constants
        integer(INTEGER_TYPE), parameter :: CONVERGENCE_ERROR_NOT_USED = 999999
        real(REAL_TYPE), parameter :: N_PERCENT_TOLERATED_INACCURATE_PLASTIC_POINTS = 0.1
        integer(INTEGER_TYPE), parameter :: N_ADDITIONAL_TOLERATED_INACCURATE_PLASTIC_POINTS = 3
        
        ! Input-related constants
        logical, parameter :: DO_USE_STEP_EXTENSION = .true.
        logical, parameter :: DO_NOT_USE_STEP_EXTENSION = .false.
        character(len = 4), parameter :: GIP_FILE_EXTENSION = '.GIP'
        character(len = 5), parameter :: GIP_STEP_FILE_EXTENSION = '.GIP_'
        character(len = 8), parameter :: END_STEP_SPECIFIER = 'EndStep='
        character(len = 10), parameter :: START_STEP_SPECIFIER = 'StartStep='
        character(len = 5), parameter :: CPS_FILE_EXTENSION = '.CPS_'
        character(len = 4), parameter :: GOM_FILE_EXTENSION = '.GOM'
        character(len = 5), parameter :: BRF_FILE_EXTENSION = '.BRF_'
        character(len = 5), parameter :: HYDRAULICHEAD_FILE_EXTENSION = '.HHBF'
        character(len=5), parameter:: MPMAPFROMEXT_FILE_EXTENSION = '.MAPF'
        character(len = 4), parameter :: PRESCRIBEDVELOCITY_FILE_EXTENSION = '.PVF'
        character(len = 4), parameter :: PRESCRIBEDACCELERATION_FILE_EXTENSION = '.PAF'
        character(len = 5), parameter :: PHREATIC_SURFACE_FILE_EXTENSION = '.PSF_'
        

        ! Particle-related constants
        integer(INTEGER_TYPE), parameter :: MATERIALPARTICLE = 1
        integer(INTEGER_TYPE), parameter :: VIRTUALPARTICLE = 2
        integer(INTEGER_TYPE), parameter :: NEWVIRTUALPARTICLE = 3
        integer(INTEGER_TYPE), parameter :: REMOVEDPARTICLE = -1
        integer(INTEGER_TYPE), parameter :: ADDEDPARTICLE = -1
        integer(INTEGER_TYPE), parameter :: ANYPARTICLE = 0

        integer(INTEGER_TYPE), parameter :: TRANSLATIONALSWEEP = 0
        integer(INTEGER_TYPE), parameter :: ROTATIONALSWEEP = 1

        ! Material and calculation type
        character(len=64), parameter :: DRY_SOIL = 'dry_material'
        character(len=64), parameter :: SATURATED_SOIL_DRAINED = 'saturated_material_drained'
        character(len=64), parameter :: SATURATED_SOIL_UNDRAINED_EFFECTIVE = 'saturated_material_undrained_effective'
        character(len=64), parameter :: SATURATED_SOIL_UNDRAINED_TOTAL = 'saturated_material_undrained_total'
        character(len=64), parameter :: SATURATED_SOIL_COUPLED = 'saturated_material_coupled'
        character(len=64), parameter :: UNSATURATED_SOIL_TWOPHASE = 'unsaturated_material_2phase_suction'
        character(len=64), parameter :: UNSATURATED_SOIL_THREEPHASE = 'unsaturated_material_3phase_coupled'
        character(len=64), parameter :: LIQUID = 'liquid'
        
        ! Soil model related constants
        integer(INTEGER_TYPE), parameter :: LINEAR_ELASTICITY = 101
        integer(INTEGER_TYPE), parameter :: MOHR_COULOMB_STANDARD = 102
        integer(INTEGER_TYPE), parameter :: ISOTROPIC_CONE_HARDENING = 103
        integer(INTEGER_TYPE), parameter :: HYPOPLASTICITY_SAND = 104
        integer(INTEGER_TYPE), parameter :: MOHR_COULOMB_STRAIN_SOFTENING = 105
        integer(INTEGER_TYPE), parameter :: MODIFIED_CAM_CLAY = 106
        integer(INTEGER_TYPE), parameter :: SANISAND = 108
        integer(INTEGER_TYPE), parameter :: MOHR_COULOMB_TEUNISSEN = 109
        integer(INTEGER_TYPE), parameter :: External = 149 ! 149 is the maximum allowable number 

        ! Soil model names
        character(len=64), parameter :: ESM_LINEAR_ELASTICITY = 'linear_elasticity'
        character(len=64), parameter :: ESM_MOHR_COULOMB_STANDARD = 'mohr_coulomb'
        character(len=64), parameter :: ESM_MOHR_COULOMB_TEUNISSEN = 'mohr_coulomb_teunissen'
        character(len=64), parameter :: ESM_MOHR_COULOMB_STRAIN_SOFTENING = 'strain_softening_mohr_coulomb'
        character(len=64), parameter :: ESM_HYPOPLASTICITY_SAND = 'hypoplasticity'
        character(len=64), parameter :: ESM_MODIFIED_CAM_CLAY = 'modified_cam_clay'
        character(len=64), parameter :: ESM_SANISAND = 'sani_sand'
        character(len=64), parameter :: ESM_EXTERNAL_SOIL_MODEL = 'external_soil_model'
        character(len=64), parameter :: ESM_NEWTONIAN_LIQUID = 'newtonian_liquid'
        character(len=64), parameter :: ESM_BINGHAM_LIQUID = 'bingham_liquid'
        character(len=64), parameter :: ESM_FRICTIONAL_LIQUID = 'frictional_liquid'
        character(len=64), parameter :: ESM_RIGID_BODY = 'rigid_body'
		  character(len=64), parameter :: ESM_NON_ASSOC_MOHR_COULOMB = 'non_assoc_mohr_coulomb'
        
       ! User-defined hydraulic model names
        character(len=64), parameter :: SWRC_VANGENUCHTEN = 'van_genuchten'
        character(len=64), parameter :: SWRC_LINEAR = 'linear'
        
        ! User-defined hydraulic conductivity model names
        character(len=64), parameter :: HCC_CONSTANT = 'constant'
        character(len=64), parameter :: HCC_HILLEL = 'hillel'
        character(len=64), parameter :: HCC_MUALEM = 'mualem'
        
        ! External soil model libraries
        character(len=64), parameter :: UMAT_LINEAR_ELASTICITY = 'A3DLinearElasticity.dll'
        character(len=64), parameter :: UMAT_MOHR_COULOMB_STANDARD = 'A3DMohrCoulombStandard.dll'
        character(len=64), parameter :: UMAT_MOHR_COULOMB_TEUNISSEN = 'A3DMohrCoulombTeunissen.dll'
        character(len=64), parameter :: UMAT_MOHR_COULOMB_STRAIN_SOFTENING = 'A3DMohrCoulombStrainSoftening.dll'
        character(len=64), parameter :: UMAT_ISOTROPIC_CONE_HARDENING = 'A3DIsotropicConeHardening.dll'
        character(len=64), parameter :: UMAT_HYPOPLASTICITY_SAND = 'A3DHypoplasticitySand.dll'
        character(len=64), parameter :: UMAT_MODIFIED_CAM_CLAY = 'A3DModifiedCamClay.dll'
        character(len=64), parameter :: UMAT_SANISAND = 'A3DSaniSand.dll'
        character(len=64), parameter :: UMAT_BINGHAM = 'A3DBingham.dll'

        ! Hypoplastic model state variables
        integer(INTEGER_TYPE), parameter :: HP_STVAR_MOB_FRICTION = 1
        integer(INTEGER_TYPE), parameter :: HP_STVAR_VOID         = 2

        integer(INTEGER_TYPE), parameter :: MAT_G = 1
        integer(INTEGER_TYPE), parameter :: MAT_NU_EFF = 2
        integer(INTEGER_TYPE), parameter :: MAT_DRY_SOIL_WEIGHT = 9
        integer(INTEGER_TYPE), parameter :: MAT_MIXED_WEIGHT = 10
        integer(INTEGER_TYPE), parameter :: MAT_K0 = 15
        integer(INTEGER_TYPE), parameter :: MAT_POROSITY = 20
        integer(INTEGER_TYPE), parameter :: MAT_BULK_MODULUS_WATER = 25
        integer(INTEGER_TYPE), parameter :: MAT_NU_UNDRAINED = 28
        integer(INTEGER_TYPE), parameter :: MAT_HP_PT = 47

        integer(INTEGER_TYPE), parameter :: MAX_MATERIALS = 10 ! maximum number of different materials in a project

        real(REAL_TYPE), parameter :: SKIP_K0_THRESHOLD = 0.0

        ! Loading-related constants
        real(REAL_TYPE), parameter :: DEFAULT_INITIAL_MULTIPLIER = 0.0
        real(REAL_TYPE), parameter :: DEFAULT_FINAL_MULTIPLIER = 1.0
        real(REAL_TYPE), parameter :: DEFAULT_TIME_MULTIPLIER = 1.0

        ! Load types for material points
        integer(INTEGER_TYPE), parameter :: LOADTYPE_SOLID  = 1
        integer(INTEGER_TYPE), parameter :: LOADTYPE_LIQUID = 2
        integer(INTEGER_TYPE), parameter :: LOADTYPE_GAS    = 3

        ! Structure-related constants
        integer(INTEGER_TYPE), parameter :: PARABOLIC_TRAJECTORY = 1
        integer(INTEGER_TYPE), parameter :: RECORDED_TRAJECTORY = 2

        ! Output-related constants
        integer(INTEGER_TYPE), parameter :: MAXOUTPUTPARTICLES = 99
        integer(INTEGER_TYPE), parameter :: MAX_OUTPUT_ELEMENTS = 99
        integer(INTEGER_TYPE), parameter :: MAXOUTPUTNODES = 10
        integer(INTEGER_TYPE), parameter :: MAXREACTIONHISTORYNODES = 200
        integer(INTEGER_TYPE), parameter :: MAXNMATERIALS = 30
        integer(INTEGER_TYPE), parameter :: MAXOUTPUTSUBSTEPS = 99
        integer(INTEGER_TYPE), parameter :: MAXOUTPUTTIMESTEPS = 99
        integer(INTEGER_TYPE), parameter :: MAXOUTPUTREALTIME = 99
        integer(INTEGER_TYPE), parameter :: MAXNUMBERSOILLAYERS = 10
        integer(INTEGER_TYPE), parameter :: MAXOUTPUTSURFACES = 50
      
        ! Used file units
        integer(INTEGER_TYPE), parameter :: TMP_UNIT= 21
        integer(INTEGER_TYPE), parameter :: OUTunit = 22
        integer(INTEGER_TYPE), parameter :: TSTunit = 24
        integer(INTEGER_TYPE), parameter :: LOGunit = 25
        integer(INTEGER_TYPE), parameter :: BMRunit = 29
        integer(INTEGER_TYPE), parameter :: BMSunit = 347
        
        integer(INTEGER_TYPE), parameter :: CPSunit = 13
        integer(INTEGER_TYPE), parameter :: BRFunit = 14
        integer(INTEGER_TYPE), parameter :: CPSINITUNIT = 15

        integer(INTEGER_TYPE), parameter :: RXunit = 47
        
        integer(INTEGER_TYPE), parameter :: ENGunit = 41
        integer(INTEGER_TYPE), parameter :: INFunit = 48
        integer(INTEGER_TYPE), parameter :: TVMunit = 50
        integer(INTEGER_TYPE), parameter :: VTKunit = 600

        integer(INTEGER_TYPE), parameter :: DMPunit = 43
        integer(INTEGER_TYPE), parameter :: GIPunit = 45
        integer(INTEGER_TYPE), parameter :: TRJunit = 52

        integer(INTEGER_TYPE), parameter :: TMP1Unit = 60
        integer(INTEGER_TYPE), parameter :: TMP2Unit = 61

        integer(INTEGER_TYPE), parameter :: GOMunit = 100

        integer(INTEGER_TYPE), parameter :: ELEMUnit= 800 ! first unit for element output (max until 899)
        integer(INTEGER_TYPE), parameter :: PARUnit = 900 ! first unit for particle output (max until 999)
        integer(INTEGER_TYPE), parameter :: SURFReacUnit = 1000
        
            ! Prescribed velocity from file 
        integer(INTEGER_TYPE), parameter :: PVFunit = 70
        
        ! Prescribed acceleration from file
        integer(INTEGER_TYPE), parameter :: PAFunit = 71
        
        ! Hydraulic Head from file 
        integer(INTEGER_TYPE), parameter :: HHBFUnit = 72
        
        !Material Point Mapping from external file
        integer(INTEGER_TYPE), parameter :: MAPUnit = 73
        
           ! Initial phreatic surface from file (for K0 procedure) 
        integer(INTEGER_TYPE), parameter :: PSFUnit = 74

        ! Temporarily introduced constant in the frame of code cleanup
        integer(INTEGER_TYPE), parameter :: HARD_ENTITY = 2
        integer(INTEGER_TYPE), parameter :: SOFT_ENTITY = 1

        ! 2 Layer formulation
        real(REAL_TYPE), parameter :: POROSITYTHRES = 0.95

        ! error handling
        integer(INTEGER_TYPE), parameter :: ERROR   = 1
        integer(INTEGER_TYPE), parameter :: NO_ERROR= 0
        integer(INTEGER_TYPE), parameter :: WARNING = -1

        ! Feedback level
        integer(INTEGER_TYPE), parameter :: FEEDBACK_LEVEL_RELEASE = 1
        integer(INTEGER_TYPE), parameter :: FEEDBACK_LEVEL_DEBUG   = 2
        integer(INTEGER_TYPE), parameter :: FEEDBACK_LEVEL_ALWAYS  = huge(FEEDBACK_LEVEL_RELEASE)
        
        ! Anura3D file versions
        character(len=64), parameter :: Anura3D_v2024 = '### Anura3D_2024 ###' ! internal release/debug
        character(len=64), parameter :: Anura3D_v2023 = '### Anura3D_2023 ###' ! internal release/debug
        character(len=64), parameter :: Anura3D_v2022 = '### Anura3D_2022 ###' ! internal release/debug
        character(len=64), parameter :: Anura3D_v2021 = '### Anura3D_2021 ###' ! internal release/debug
        character(len=64), parameter :: Anura3D_v2019_2 = '### Anura3D_2019.2 ###' ! v2019.2 external and internal version release/debug
        
        ! MaterialPointType parameters
        integer(INTEGER_TYPE), parameter :: MaterialPointTypeUndefined  = -1
        integer(INTEGER_TYPE), parameter :: MaterialPointTypeSolid       = 1
        integer(INTEGER_TYPE), parameter :: MaterialPointTypeLiquid      = 2
        integer(INTEGER_TYPE), parameter :: MaterialPointTypeSolidLiquid = 3
        integer(INTEGER_TYPE), parameter :: MaterialPointTypeMixture     = 4
        
        ! ContainedMaterialTypes parameters
        integer(INTEGER_TYPE), parameter :: ContainedMaterialTypeUNDEFINED   = -1
        integer(INTEGER_TYPE), parameter :: ContainedMaterialTypeSOLID       =  1
        integer(INTEGER_TYPE), parameter :: ContainedMaterialTypeLIQUID      =  2
        integer(INTEGER_TYPE), parameter :: ContainedMaterialTypeSOLIDLIQUID =  3
        
        ! PhaseStatus parameters 
        integer(INTEGER_TYPE), parameter :: PhaseStatusUNDEFINED = -1
        integer(INTEGER_TYPE), parameter :: PhaseStatusSOLID = 1
        integer(INTEGER_TYPE), parameter :: PhaseStatusLIQUID = 2
        
        ! elementtype constants, the PROTECTED attribute prevents code outside the module from changing its value
        integer(INTEGER_TYPE), protected :: ELEMENTNODES = -1 ! number of nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTVERTICES = -1 ! number of vertice nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTGAUSSPOINTS = -1 ! number of gauss points per element
        integer(INTEGER_TYPE), protected :: ELEMENTBOUNDARYNODES = -1 ! number of nodes of the element boundary
        integer(INTEGER_TYPE), protected :: ELEMENTBOUNDARYGAUSSPOINTS = -1 ! number of gauss ponts of the element boundary
        integer(INTEGER_TYPE), protected :: ELEMENTSIDES = -1 ! number of sides per element (lines in 2D, surfaces in 3D)
        integer(INTEGER_TYPE), protected :: VTK_CELL = -1 ! cell type for VTK output
        character(len=255), protected :: ELEMENTTYPE = '' ! elementtype
        
        procedure(DummyCheckForGlobPosPointer), pointer:: CheckForGlobPosPointer => NULL()
        procedure(DummyGauss_Q1Pointer), pointer:: Gauss_Q1Pointer => NULL()
        procedure(DummyInitialiseShapeFunctionsPointer), pointer:: InitialiseShapeFunctionsPointer => NULL()
        procedure(DummyInitialiseShapeFunctionsPointer), pointer:: InitialiseShapeFunctionsBoundaryPointer => NULL()
        procedure(DummyIsInsideElementLocPosPointer), pointer:: IsInsideElementLocPosPointer => NULL()
        procedure(DummyGetMinAltitudePointer), pointer:: GetMinAltitudePointer => NULL()
        procedure(DummyInitialLocalMaterialPointCoordinatesPointer), pointer :: InitialLocalMaterialPointCoordinatesPointer => NULL()
        procedure(DummyShapeLocPosPointer), pointer:: ShapeLocPosPointer => NULL() 
        procedure(DummyRearrangeConnectivitiesPointer), pointer:: RearrangeConnectivitiesPointer => NULL() 
        
        ! Output data format 
        character(len=26) :: VTK_VERSION = '# vtk DataFile Version 3.0'
        character(len=6)  :: BINARY_FORMAT = 'BINARY'
        character(len=5)  :: ASCII_FORMAT  = 'ASCII'
        
        ! Pointers for writing VTK output files
        procedure(DummyInitialiseVTKFilePointer), pointer :: InitialiseVTKFilePointer => NULL()  
        procedure(DummyWriteVTKFloatScalarDataPointer), pointer :: WriteVTKFloatScalarDataPointer => NULL()
        procedure(DummyWriteVTKIntegerScalarDataPointer), pointer :: WriteVTKIntegerScalarDataPointer => NULL()  
        procedure(DummyWriteVTKFloatVectorDataPointer), pointer :: WriteVTKFloatVectorDataPointer => NULL()
        procedure(DummyWriteVTKIntegerVectorDataPointer), pointer :: WriteVTKIntegerVectorDataPointer => NULL()
        procedure(DummyWriteVTKFloatTensorDataPointer), pointer :: WriteVTKFloatTensorDataPointer => NULL()
        procedure(DummyWriteVTKIntegerTensorDataPointer), pointer :: WriteVTKIntegerTensorDataPointer => NULL()
        
        interface WriteVTKScalarDataPointer 
            procedure(), pointer :: WriteVTKFloatScalarDataPointer
            procedure(), pointer :: WriteVTKIntegerScalarDataPointer
        end interface WriteVTKScalarDataPointer

        interface WriteVTKVectorDataPointer
            procedure(), pointer :: WriteVTKFloatVectorDataPointer 
            procedure(), pointer :: WriteVTKIntegerVectorDataPointer 
        end interface WriteVTKVectorDataPointer
    
        interface WriteVTKTensorDataPointer
            procedure(), pointer :: WriteVTKFloatTensorDataPointer
            procedure(), pointer :: WriteVTKIntegerTensorDataPointer
        end interface WriteVTKTensorDataPointer 
    
      contains
    
        subroutine DummyShapeLocPosPointer(LocPos, HS, dHS) 
        implicit none
           real(REAL_TYPE), dimension(:), intent(in) :: LocPos
           real(REAL_TYPE), dimension(:), intent(out) :: HS 
           real(REAL_TYPE), dimension(:, :), intent(out) :: dHS 
        end subroutine DummyShapeLocPosPointer
       
        subroutine DummyInitialLocalMaterialPointCoordinatesPointer(IParticle, SolidPointsElement, LiquidPointsElement, WeiGP, PosGP)
           implicit none
           integer(INTEGER_TYPE), intent(in) :: IParticle, SolidPointsElement, LiquidPointsElement 
           real(REAL_TYPE), intent(inout) :: WeiGP 
           real(REAL_TYPE), dimension(:), intent(inout) :: PosGP 
        end subroutine DummyInitialLocalMaterialPointCoordinatesPointer 
    
        subroutine DummyGetMinAltitudePointer(NodeNr, NodeCoord, Lmin)
          implicit none
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord 
          real(REAL_TYPE), intent (inout):: Lmin 
          integer(INTEGER_TYPE):: NodeNr(:) 
        end subroutine DummyGetMinAltitudePointer
    
        subroutine DummyInitialiseShapeFunctionsPointer(HS, dHS, Wt)
           implicit none
           real(REAL_TYPE), dimension(:, :), intent(inout) :: HS 
           real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS 
           real(REAL_TYPE), dimension(:), intent(inout) :: Wt
        end subroutine DummyInitialiseShapeFunctionsPointer
      
        subroutine DummyCheckForGlobPosPointer(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
           implicit none 
           real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
           integer(INTEGER_TYPE), intent(in) :: ElementID
           real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
           real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
           integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
           integer(INTEGER_TYPE), intent(out) :: CrossedSide
           logical, intent(out) :: IsInside                                             
        end subroutine DummyCheckForGlobPosPointer
                     
        subroutine DummyGauss_Q1Pointer(IGaussPoint, PosGP, WeiGP)
           implicit none
           integer(INTEGER_TYPE), intent(in) :: IGaussPoint
           real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
           real(REAL_TYPE), intent(inout) :: WeiGP
        end subroutine DummyGauss_Q1Pointer
        
        logical function DummyIsInsideElementLocPosPointer(LocPos)
           implicit none
           real(REAL_TYPE), dimension(:), intent(in) :: LocPos
        end function DummyIsInsideElementLocPosPointer
        
        subroutine DummyRearrangeConnectivitiesPointer(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
           implicit none
           integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
           real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
           integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
           real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal
        end subroutine DummyRearrangeConnectivitiesPointer

        subroutine DummyInitialiseVTKFilePointer(NumberMaterialPoints, VTKFileName)
           implicit none
           integer(INTEGER_TYPE), intent(in) :: NumberMaterialPoints
           character*(*), intent(in) :: VTKFileName
        end subroutine DummyInitialiseVTKFilePointer
    
        subroutine DummyWriteVTKIntegerScalarDataPointer(DataSetName,DataSet)
           implicit none
           character(len=*), intent(in) :: DataSetName
           integer(INTEGER_TYPE), dimension(:), intent(in) :: DataSet
        end subroutine DummyWriteVTKIntegerScalarDataPointer
    
        subroutine DummyWriteVTKFloatScalarDataPointer(DataSetName,DataSet)
           implicit none
           character(len=*), intent(in) :: DataSetName
           real(REAL_TYPE), dimension(:), intent(in) :: DataSet
        end subroutine DummyWriteVTKFloatScalarDataPointer
    
        subroutine DummyWriteVTKIntegerVectorDataPointer(DataSetName,DataSet)
           implicit none
           character(len=*), intent(in) :: DataSetName
           integer(INTEGER_TYPE), dimension(:,:), intent(in) :: DataSet
        end subroutine DummyWriteVTKIntegerVectorDataPointer
    
        subroutine DummyWriteVTKFloatVectorDataPointer(DataSetName, DataSet)
           implicit none
           character(len=*), intent(in) :: DataSetName
           real(REAL_TYPE), dimension(:, :), intent(in) :: DataSet
        end subroutine DummyWriteVTKFloatVectorDataPointer
    
        subroutine DummyWriteVTKFloatTensorDataPointer(DataSetName,DataSet)
           implicit none
           character(len=*), intent(in) :: DataSetName
           real(REAL_TYPE), dimension(:,:,:), intent(in) :: DataSet
        end subroutine DummyWriteVTKFloatTensorDataPointer
    
        subroutine DummyWriteVTKIntegerTensorDataPointer(DataSetName,DataSet)
           implicit none
           character(len=*), intent(in) :: DataSetName
           integer(INTEGER_TYPE), dimension(:,:,:), intent(in) :: DataSet
        end subroutine DummyWriteVTKIntegerTensorDataPointer
                               
        
        subroutine SetDimension(NDimension, Axisymmetric, Cylindric)
        !**********************************************************************
        !
        !    Function : Sets constants that depend on the analysis dimension
        !
        ! I  NDimension : 2D analysis = 2, 3D analysis = 3
        ! I  Axisymmetric: .true. if 2D axisymmetric analysis
        ! I  Cylindric: : .true. if 3D analysis with cylindrical coordinates
        !
        !**********************************************************************
        implicit none
        
        integer(INTEGER_TYPE), intent(in) :: NDimension
        logical, intent(in) :: Axisymmetric, Cylindric
        
        select case (NDimension)
            
          case(2) ! 2D CASE
            NDIM = 2 
            NDOFL = 2 
            NDOFLQS = 3 ! NDOFLQS = NDOFL+1
            NTENSOR = 4 ! T_xx,T_yy,T_zz,T_xy
            NVECTOR = 2 ! V_x,V_y
            NPRINCIPAL = 3 ! T_1,T_2,T_3
            allocate(DEFAULT_GRAVITY_DIRECTION(2))
            DEFAULT_GRAVITY_DIRECTION = (/ 0.0, -1.0 /) 
            if ( Axisymmetric ) ISAXISYMMETRIC = .true.
            if ( Cylindric ) call GiveError('Cylindrical coordinates can only be used in a 3D analysis. [subroutine SetDimension()].')
            
          case(3) ! 3D CASE
            NDIM = 3 
            NDOFL = 3 
            NDOFLQS = 4 ! NDOFLQS = NDOFL+1
            NTENSOR = 6 ! T_xx,T_yy,T_zz,T_xy,T_yz,T_xz
            NVECTOR = 3 ! V_x,V_y,V_z
            NPRINCIPAL = 3 ! T_1,T_2,T_3
            allocate(DEFAULT_GRAVITY_DIRECTION(3))
            DEFAULT_GRAVITY_DIRECTION = (/ 0.0, -1.0,  0.0/) 
            if ( Axisymmetric ) call GiveError('Axisymmetric calculation can only be performed in a 2D analysis. [subroutine SetDimension()].')
            if ( Cylindric ) IS3DCYLINDRIC = .true.
            
          case default ! not defined
            call GiveError('Dimension not defined. [subroutine SetDimension()].')
          
        end select  
        
        end subroutine SetDimension
        
        
        subroutine SetFormulation(Formulation)
        !**********************************************************************
        !
        !    Function : Sets constant that depends on formulation
        !
        ! I  Formulation: single-point or double-point 
        ! 
        !**********************************************************************
        
        implicit none
        
        character(len=*), intent(in) :: Formulation
        
        select case (Formulation)
            
          case(SINGLE_POINT) 
            NFORMULATION = 1
                        
          case(DOUBLE_POINT) 
            NFORMULATION = 2
            
          case default ! not defined
            call GiveError('Formulation not defined. [subroutine SetFormulation()].')
          
        end select  
        
        end subroutine SetFormulation


        subroutine SetElementType(ReadElementType)
        !**********************************************************************
        !
        !    Function : Sets constants that depend on the element type 
        !               Note: Only the following element types are available, 
        !                     triangular 3-noded (TRI3)
        !                     tetrahedral_old (TETRAOLD)
        !
        ! I  ReadElementType
        ! 
        !**********************************************************************
        
        implicit none
        
        character(len=*), intent(in) :: ReadElementType
                
        select case (ReadElementType)
            
          case(TRI3) ! 'triangular_3-noded'
            ELEMENTNODES = 3
            ELEMENTVERTICES = 3
            ELEMENTGAUSSPOINTS = 1
            ELEMENTBOUNDARYNODES = 2
            ELEMENTBOUNDARYGAUSSPOINTS = 1
            ELEMENTSIDES = 3
            VTK_CELL = 5 ! VTK_TRIANGLE(=5)
            ELEMENTTYPE = TRI3
            
          case(TETRAOLD) ! 'tetrahedral_old' 
            ELEMENTNODES = 4 
            ELEMENTVERTICES = 4
            ELEMENTGAUSSPOINTS = 1
            ELEMENTBOUNDARYNODES = 3
            ELEMENTBOUNDARYGAUSSPOINTS = 1
            ELEMENTSIDES = 4
            VTK_CELL = 10 ! VTK_TETRA(=10)
            ELEMENTTYPE = TETRAOLD
            
          !**************************NOT AVAILABLE***************************
          !case(TRI6) ! 'triangular_6-noded'
          !  ELEMENTNODES = 6
          !  ELEMENTVERTICES = 3
          !  ELEMENTGAUSSPOINTS = 3
          !  ELEMENTBOUNDARYNODES = 3
          !  ELEMENTBOUNDARYGAUSSPOINTS = 2
          !  ELEMENTSIDES = 3
          !  VTK_CELL = 22 ! VTK_QUADRATIC_TRIANGLE(=22)
          !  ELEMENTTYPE = TRI6
          !    
          !case(QUAD4) ! 'quadrilateral_4-noded'
          !  ELEMENTNODES = 4
          !  ELEMENTVERTICES = 4
          !  ELEMENTGAUSSPOINTS = 1
          !  ELEMENTBOUNDARYNODES = 2
          !  ELEMENTBOUNDARYGAUSSPOINTS = 1
          !  ELEMENTSIDES = 4
          !  VTK_CELL = 9 ! VTK_QUAD(=9)
          !  ELEMENTTYPE = QUAD4
          !    
          !case(QUAD8) ! 'quadrilateral_8-noded'
          !  ELEMENTNODES = 8
          !  ELEMENTVERTICES = 4
          !  ELEMENTGAUSSPOINTS = 4
          !  ELEMENTBOUNDARYNODES = 3
          !  ELEMENTBOUNDARYGAUSSPOINTS = 2
          !  ELEMENTSIDES = 4
          !  VTK_CELL = 23 ! VTK_QUADRATIC_QUAD(=23)
          !  ELEMENTTYPE = QUAD8
          !    
          !case(TETRA4) ! 'tetrahedral_4-noded'
          !  ELEMENTNODES = 4 
          !  ELEMENTVERTICES = 4
          !  ELEMENTGAUSSPOINTS = 1
          !  ELEMENTBOUNDARYNODES = 3
          !  ELEMENTBOUNDARYGAUSSPOINTS = 1
          !  ELEMENTSIDES = 4
          !  VTK_CELL = 10 ! VTK_TETRA(=10)
          !  ELEMENTTYPE = TETRA4 
          !    
          !case(TETRA10) ! 'tetrahedral_10-noded'
          !  ELEMENTNODES = 10
          !  ELEMENTVERTICES = 4
          !  ELEMENTGAUSSPOINTS = 4
          !  ELEMENTBOUNDARYNODES = 6 
          !  ELEMENTBOUNDARYGAUSSPOINTS = 3
          !  ELEMENTSIDES = 4
          !  VTK_CELL = 24 ! VTK_QUADRATIC_TETRA(=24)
          !  ELEMENTTYPE = TETRA10
          !    
          !case(HEXA8) ! 'hexahedral_8-noded'
          !  ELEMENTNODES = 8
          !  ELEMENTVERTICES = 8
          !  ELEMENTGAUSSPOINTS = 1
          !  ELEMENTBOUNDARYNODES = 4
          !  ELEMENTBOUNDARYGAUSSPOINTS = 1
          !  ELEMENTSIDES = 6
          !  VTK_CELL = 12 ! VTK_HEXAHEDRON(=12)
          !  ELEMENTTYPE = HEXA8
          !    
          !case(HEXA20) ! 'hexahedral_20-noded'
          !  ELEMENTNODES = 20
          !  ELEMENTVERTICES = 8
          !  ELEMENTGAUSSPOINTS = 4
          !  ELEMENTBOUNDARYNODES = 8
          !  ELEMENTBOUNDARYGAUSSPOINTS = 4
          !  ELEMENTSIDES = 6
          !  VTK_CELL = 25 ! VTK_QUADRATIC_HEXAHEDRON(=25)
          !  ELEMENTTYPE = HEXA20
          !******************************************************************
          
          case default  ! not defined
            call GiveError('Element type not defined. [subroutine SetElementType()].')
            
        end select
 
        end subroutine SetElementType

       
      end module ModGlobalConstants

      !--------------------------------------------------------------------------------------------
      module ModGlobalData
      use ModGlobalConstants
      implicit none
      private
      logical :: IsAnError = .False.
      integer(INTEGER_TYPE), parameter :: MessageSize = 4096
      character(len = MessageSize) :: ErrorMessage = ''
      character(len = 5), parameter :: UninitializedMessage = '     '

      public ::  HasAnErrorHappend, SetAnErrorHappend, MessageSize

      contains

      logical function HasAnErrorHappend()
      !**********************************************************************
      !
      ! Function:  check if an error has just happened
      !
      !**********************************************************************
      implicit none

      HasAnErrorHappend = IsAnError

      end function HasAnErrorHappend

      subroutine SetAnErrorHappend(HasErrorHappend)
      implicit none
      logical, intent(in):: HasErrorHappend

      IsAnError = HasErrorHappend

      end subroutine SetAnErrorHappend

      end module ModGlobalData

