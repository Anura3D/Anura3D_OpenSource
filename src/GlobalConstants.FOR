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


	  module ModGlobalConstants
      !**********************************************************************
      !
      !    Function:  Contains definition of global constants.
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr. 2022) $
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
        character(len=27), parameter :: QUAD4_NURBS = 'quadrilateral_4-noded_NURBS'
        character(len=21), parameter :: QUAD8 = 'quadrilateral_8-noded'
        character(len=19), parameter :: TETRA4 = 'tetrahedral_4-noded'
        character(len=20), parameter :: TETRA10 = 'tetrahedral_10-noded'
        character(len=18), parameter :: HEXA8 = 'hexahedral_8-noded'
        
        character(len=18), parameter :: HEXA_NURBS = 'hexahedral_NURBS' 
        
        character(len=19), parameter :: HEXA20 = 'hexahedral_20-noded'
        character(len=19), parameter :: TETRAOLD = 'tetrahedral_old' ! this is the elementtype that was used before v2018.2
        
        ! General constants
        integer(INTEGER_TYPE), parameter :: MAX_EXTENSION_LENGTH = 20
        integer(INTEGER_TYPE), parameter :: MAX_FILENAME_LENGTH = 200
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

        ! MP initialization from file constants 
        integer(INTEGER_TYPE), parameter :: MP_INIT_GAUSSPT = 0
        integer(INTEGER_TYPE), parameter :: MP_INIT_FILE = 1
        
        character(len = 11), parameter :: MP_INITIALIZATION_GAUSSPOINT_SPECIFIER = 'gauss_point'
        character(len = 4), parameter :: MP_INITIALIZATION_FILE_SPECIFIER = 'file'
        
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
        character(len = 4), parameter :: PRESCRIBEDVELOCITY_FILE_EXTENSION = '.PVF'
        character(len = 4), parameter :: PRESCRIBEDACCELERATION_FILE_EXTENSION = '.PAF'
        
        character(len = 8), parameter :: &
            MPINIT_FILE_EXTENSION = '.MP_INIT' ! contains MP global position and its material assignment
        character(len = 9), parameter :: &
            MATINFO_FILE_EXTENSION = '.MAT_INFO' ! contains material information
        

        
        ! Particle-related constants
        integer(INTEGER_TYPE), parameter :: MATERIALPARTICLE = 1
        integer(INTEGER_TYPE), parameter :: VIRTUALPARTICLE = 2
        integer(INTEGER_TYPE), parameter :: NEWVIRTUALPARTICLE = 3
        integer(INTEGER_TYPE), parameter :: REMOVEDPARTICLE = -1
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
        character(len=64), parameter :: ESM_ARB_Model = 'Arb_Model'
		
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
        integer(INTEGER_TYPE), parameter :: MAXOUTPUTPARTICLES = 999
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
        
        ! Initial Water Surface from file
        integer(INTEGER_TYPE), parameter :: LIQUIDSurfUnit = 10000
        
        ! Prescribed velocity from file 
        integer(INTEGER_TYPE), parameter :: PVFunit = 70
        
        ! Prescribed acceleration from file
        integer(INTEGER_TYPE), parameter :: PAFunit = 71
        
        ! Hydraulic Head from file 
        integer(INTEGER_TYPE), parameter :: HHBFUnit = 72
        
        ! Material point initialization from file 
        integer(INTEGER_TYPE), parameter :: MPInitUnit = 73 

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
        integer(INTEGER_TYPE), protected :: ELEMENTSURFACES = -1 ! number of surfaces per element
        integer(INTEGER_TYPE), protected :: ELEMENTVERTICES = -1 ! number of vertice nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTGAUSSPOINTS = -1 ! number of gauss points per element
        integer(INTEGER_TYPE), protected :: ELEMENTBOUNDARYNODES = -1 ! number of nodes of the element boundary
        integer(INTEGER_TYPE) :: ELEMENTBOUNDARYGAUSSPOINTS = -1 ! number of gauss ponts of the element boundary !, protected
        integer(INTEGER_TYPE), protected :: ELEMENTSIDES = -1 ! number of sides per element (lines in 2D, surfaces in 3D)
        integer(INTEGER_TYPE), protected :: VTK_CELL = -1 ! cell type for VTK output
        integer(INTEGER_TYPE), protected :: VERTICESINPARAMETERSPACE = -1
        character(len=255), protected :: ELEMENTTYPE = '' ! elementtype
        
        integer(INTEGER_TYPE) :: NumberOfControlNetGrid = -1 ! this is read from the GOM file under the $$CONTROL_NET_CONNECTIVITY flag       
        integer(INTEGER_TYPE) :: ControlNetGrid = 4 ! it is a structured square net with 4 nodes. Not sure if this is hardcoded or always the case but I think the latter. 
        
        ! Temporary parameter
        !integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        
        
        ! elementtype constants for Xi direction in case of NURBS
        integer(INTEGER_TYPE), protected :: ELEMENTNODES_XI = 4!NXiKnotOrder(IPatch_Temporary)+1 ! number of nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTVERTICES_XI = 2 ! number of vertice nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTVERTICES_XI_TRACTION = 2 ! number of vertice nodes per element

        integer(INTEGER_TYPE), protected :: ELEMENTGAUSSPOINTS_XI = 1 ! number of gauss points per element
        integer(INTEGER_TYPE), protected :: ELEMENTBOUNDARYNODES_XI = 4!NXiKnotOrder(IPatch_Temporary)+1 ! number of nodes of the element boundary !
        integer(INTEGER_TYPE) :: ELEMENTBOUNDARYGAUSSPOINTS_XI = 1 ! number of gauss ponts of the element boundary !, protected
        integer(INTEGER_TYPE), protected :: ELEMENTSIDES_XI = 1
        
        ! elementtype constants for Eta direction in case of NURBS
        integer(INTEGER_TYPE), protected :: ELEMENTNODES_ETA = 4!NEtaKnotOrder(IPatch_Temporary)+1 ! number of nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTVERTICES_ETA = 2 ! number of vertice nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTVERTICES_ETA_TRACTION = 2 ! number of vertice nodes per element

        integer(INTEGER_TYPE), protected :: ELEMENTGAUSSPOINTS_ETA = 1 ! number of gauss points per element
        integer(INTEGER_TYPE), protected :: ELEMENTBOUNDARYNODES_ETA = 4!NEtaKnotOrder(IPatch_Temporary)+1 ! number of nodes of the element boundary !
        integer(INTEGER_TYPE) :: ELEMENTBOUNDARYGAUSSPOINTS_ETA = 1 ! number of gauss ponts of the element boundary !, protected
        integer(INTEGER_TYPE), protected :: ELEMENTSIDES_ETA = 1 ! number of sides per element (lines in 2D, surfaces in 3D)
        
        ! elementtype constants for Zeta direction in case of NURBS
        integer(INTEGER_TYPE), protected :: ELEMENTNODES_ZETA = 4!NZetaKnotOrder(IPatch_Temporary)+1 ! number of nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTVERTICES_ZETA = 2 ! number of vertice nodes per element
        integer(INTEGER_TYPE), protected :: ELEMENTGAUSSPOINTS_ZETA = 1 ! number of gauss points per element
        integer(INTEGER_TYPE), protected :: ELEMENTBOUNDARYNODES_ZETA = 4 !NZetaKnotOrder(IPatch_Temporary)+1 ! number of nodes of the element boundary
        integer(INTEGER_TYPE), protected :: ELEMENTBOUNDARYGAUSSPOINTS_ZETA = 1 ! number of gauss ponts of the element boundary
        integer(INTEGER_TYPE), protected :: ELEMENTSIDES_ZETA = 1 ! number of sides per element (lines in 2D, surfaces in 3D)
        
        
        
        
        
        
              !-----------------------------------------------------------------------------------------------------------------------------------------------
              real(REAL_TYPE), dimension(:), allocatable :: GPWeight 
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative
        
        
        ! Originaly in Anura3D, GPShapeFunctionDerivative is the same for all elements. 
        ! So we do not need to calculate it for each element. But in NURBS, we have a super element
        ! so we need to know for what element we are doing this for. We can add a dimension with
        ! element number
        real(REAL_TYPE), dimension(:,:), allocatable :: GPWeight_AllElements  
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunction_AllElements
        real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivative_AllElements

        ! NURBS xi direction  
        real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary_XI ! for each boundary (1D) there are a number of weights 
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary_XI ! similarly here for shape function evaluation (1D) 
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_XI ! similarly here for shape function derivative (1D) 
        
        !NURBS eta direction 
        !real(REAL_TYPE), dimension(:), allocatable :: GPWeight_ETA 
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction_ETA
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative_ETA

        real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary_ETA
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary_ETA
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ETA
        
        !NURBS zeta direction 
        !real(REAL_TYPE), dimension(:), allocatable :: GPWeight_ZETA 
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction_ZETA
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative_ZETA

        real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary_ZETA
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary_ZETA
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ZETA
        
        

        ! ------------------ all elements needed for traction implementation -------------------------------------------------

        
        real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_XI_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_XI_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_XI_AllElements
        
        
        real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_ETA_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_ETA_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ETA_AllElements
        
        
        real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_ZETA_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_ZETA_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ZETA_AllElements
              
              
              
              
              !-------------------------------------------------------------------------------------------------------------------------------------------------------
              
              
              
              
              
        
        
        
        
        
        
        
        
        procedure(DummyCheckForGlobPosPointer), pointer:: CheckForGlobPosPointer => NULL()
        procedure(DummyGauss_Q1Pointer), pointer:: Gauss_Q1Pointer => NULL()
        procedure(DummyInitialiseShapeFunctionsPointer), pointer:: InitialiseShapeFunctionsPointer => NULL()
        procedure(DummyInitialiseShapeFunctionsPointer), pointer:: InitialiseShapeFunctionsBoundaryPointer => NULL()
        procedure(DummyIsInsideElementLocPosPointer), pointer:: IsInsideElementLocPosPointer => NULL()
        procedure(DummyGetMinAltitudePointer), pointer:: GetMinAltitudePointer => NULL()
        procedure(DummyInitialLocalMaterialPointCoordinatesPointer), pointer :: InitialLocalMaterialPointCoordinatesPointer => NULL()
        procedure(DummyShapeLocPosPointer), pointer:: ShapeLocPosPointer => NULL() 
        procedure(DummyRearrangeConnectivitiesPointer), pointer:: RearrangeConnectivitiesPointer => NULL() 
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! I moved this from MeshInfo.FOR to here... 
        
        
        ! Nodal coordinates, updated nodal coordinates, The initial coordinates
        ! Note that these are also considered as control points in the context of NURBS
        real(REAL_TYPE), dimension(:, :, :),  & ! I added an extra rank for multipatch implementation
          allocatable :: NodalCoordinates, NodalCoordinatesUpd, NodalOriginalCoord
        
        real(REAL_TYPE), dimension(:, :),  & ! I added an extra rank for multipatch implementation
          allocatable :: NodalCoordinates_Traction!, NodalCoordinatesUpd, NodalOriginalCoord
        
        ! Global nodal coordinates 
        real(REAL_TYPE), dimension(:, :),  &
            allocatable :: GlobalMultipatchNodalCoordinates
        
        ! Input local control point number and patch number
        integer(INTEGER_TYPE), dimension(:, :),  &
          allocatable :: Multipatch_Connecting_Local_To_Global_ControlPoints
        
        ! Connecting the traction surface to solid 3D
        integer(INTEGER_TYPE), dimension(:),  &
          allocatable :: Traction_Connecting_Local_To_Global_ControlPoints !--> needs multipatch generalization
        integer(INTEGER_TYPE), dimension(:),  &
          allocatable :: Traction_Connecting_Local_To_Global_Elements !--> needs multipatch generalization
        integer(INTEGER_TYPE), dimension(:),  &
          allocatable :: Traction_Connecting_Local_To_Global_Sides !--> needs multipatch generalization
        
        !real(REAL_TYPE), dimension(:), &
        !    allocatable :: ControlPoint_Weights ! NURBS weights - each node has a weight
        !real(REAL_TYPE), dimension(:, :),  &
        !  allocatable :: XiKnotEntries, EtaKnotEntries  
        real(REAL_TYPE), dimension(:, :, :),  & ! rank 3 for Counters%NPatches
          allocatable :: NodalPrescibedDisp ! Nodal prescribed displacement in 3 directions (Solid)
        real(REAL_TYPE), dimension(:, :, :),  & ! rank 3 for Counters%NPatches
          allocatable :: NodalPrescibedDispWater ! Nodal prescribed displacement in 3 directions (Water)
        real(REAL_TYPE), dimension(:, :, :),  & ! rank 3 for Counters%NPatches
          allocatable :: NodalPrescibedDispGas ! Nodal prescribed displacement in 3 directions (Gas)
        real(REAL_TYPE), dimension(:, :, :), &
          allocatable :: LoadValuesOnNodesSolid ! Load values (applied on nodes) for each distributed element face (solid)
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: LoadValuesOnNodesSolid_NURBS ! Load values (applied on nodes) for each distributed element face (solid)
        ! In NURBS implementation we need to know which elements below to the traction you are applying 
        ! In the future, this needs to be generalized to let the code deduce automatically which element 
        ! we are in if we are applying this traction.
        integer(INTEGER_TYPE), dimension(:), &
          allocatable :: ElementUponWhichLoadIsApplied_NURBS
        
        
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnNodesWater ! Load values (applied on nodes) for each distributed element face (water)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnNodesGas ! Load values (applied on nodes) for each distributed element face (gas)
        real(REAL_TYPE), dimension(:, :, :), &
          allocatable :: LoadValuesOnNodesSolidB ! Load values (applied on nodes) for each distributed element face (solid, load system B)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnNodesWaterB ! Load values (applied on nodes) for each distributed element face (water)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnNodesGasB ! Load values (applied on nodes) for each distributed element face (gas)
          real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: HydraulicHeadLoad ! Load values (applied on nodes) for each distributed element face (hydraulic head)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnMatPointsSolid ! Load values (applied on material points) for each distributed element face (solid)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnMatPointsWater ! Load values (applied on material points) for each distributed element face (water)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnMatPointsGas ! Load values (applied on material points) for each distributed element face (gas)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnMatPointsSolidB ! Load values (applied on material points) for each distributed element face (solid)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnMatPointsWaterB ! Load values (applied on material points) for each distributed element face (water)
        real(REAL_TYPE), dimension(:, :, :),  &
          allocatable :: LoadValuesOnMatPointsGasB ! Load values (applied on material points) for each distributed element face (gas)
        ! Prescribed boundary condition for the solid phase (1 for free, 0 for fixed)
        real(REAL_TYPE), dimension(:), allocatable :: PBoundary
        ! 4 degrees of freedom
        integer(INTEGER_TYPE), dimension(:), allocatable :: PBoundaryQuasiStatic 
        ! Prescribed boundary condition for the water phase (1 for free, 0 for fixed)
        real(REAL_TYPE), dimension(:), allocatable :: PBoundaryWater
        ! Prescribed boundary condition for the gas phase (1 for free, 0 for fixed)
        real(REAL_TYPE), dimension(:), allocatable :: PBoundaryGas
        ! Nodal load vector of external load (solid)
        real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVector
        ! Nodal load vector of external load (solid)
        real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorB
        ! Nodal load vector of external load (water)
        real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorWater
        ! Nodal load vector of external load (water)
        real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorWaterB
        ! Nodal load vector of external load (gas)
        real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorGas
        ! Nodal load vector of external load (gas)
        real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorGasB
        ! Nodal load vector of external load (HydraulicHead)
        real(REAL_TYPE), dimension(:), allocatable :: HydraulicHeadVector
        ! Element connectivities (4 node tetrahedral element)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElementConnectivities
        
        

        
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ControlNetConnectivity

        
        
        ! Element connectivities for 10-noded tetrahedral element
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElementConnectivities10Node
        ! Ties edge nodes of high-order tetrahedral elements to corner nodes
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: EdgeNodeTyingsHOE
        ! Distributed element face load connectivity (6-noded triangular element)(solid)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesSolid
        ! Distributed element face load connectivity (6-noded triangular element)(solid, load system B)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesSolidB
        ! Distributed element face load connectivity (6-noded triangular element)(water)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesWater
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesWaterB
        ! Distributed element face load connectivity (6-noded triangular element)(gas)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesGas
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesGasB
        !On the soil surface element face connectivity (6-noded triangular element)(soil surface for K0)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: SoilSurfaceNodesConnectivities
        !On the phratic surface element face connectivity (6-noded triangular element)(PHREATIC surface for K0)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: PhreaticSurfaceNodesConnectivities
        ! Distributed element face load connectivity (6-noded triangular element)(solid)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesSolid
         integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesSolidB
        ! Distributed element face load connectivity (6-noded triangular element)(water)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesWater
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesWaterB
        ! Distributed element face load connectivity (6-noded triangular element)(gas)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesGas
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesGasB
        ! Distributed element face load connectivity (6-noded triangular element)(hydraulic head)
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: HydraulicHeadNodesConnectivities
        integer(INTEGER_TYPE), dimension(:,:), allocatable :: ElementMaterialID ! Material ID for the element (rank 3 for multipatches)
        integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDof ! The global storage of degree of freedom in x-direction
        integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDofQuasiStatic ! The global storage of degree of freedom in ux,uy,uz,p
        logical, dimension(:, :), allocatable :: IsActiveElement ! Element switch, 1:Active element and 0:Inactive element (rank 2 for multipatch --> I removed this rank now)
        real(REAL_TYPE), dimension(:), allocatable :: NodalDensity 
        real(REAL_TYPE), dimension(:), allocatable :: ActiveNodeElementVolume 
        logical, dimension(:, :), allocatable :: ActiveNode ! rank 2 for multiple patches
        
        ! Multipatch implementation 
        integer(INTEGER_TYPE) :: NumberOfCommonControlPoints
        
        logical, dimension(:), allocatable :: HydraulicHeadLoadedElemID
        
        ! Variables which are updated within calculation process
        real(REAL_TYPE), dimension(:, :), allocatable :: ExtLoad ! Nodal load array of external load
        ! Nodal load array of external load (Total traction applied on nodes)
        real(REAL_TYPE), dimension(:, :, :), allocatable :: ExtLoadTotal
        real(REAL_TYPE), dimension(:, :), allocatable :: HydraulicHeadLoadTotal
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoad ! Nodal load array of internal load
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadPrevious ! Nodal load array of internal load at previous time step
        real(REAL_TYPE), dimension(:, :), allocatable :: GravityLoad ! Nodal load array of gravity load
        real(REAL_TYPE), dimension(:, :), allocatable :: BulkViscLoad
        real(REAL_TYPE), dimension(:, :), allocatable :: LumpedMassDry ! Lumped mass vector of dry soil
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: TotalVelocitySoil ! Nodal total velocity
        real(REAL_TYPE), dimension(:,:), &
          allocatable :: TotalVelocitySoilPrevious ! Nodal total velocity of the previous time step
        real(REAL_TYPE), dimension(:), &
          allocatable :: TotalDisplacementSoil ! Nodal total displacement
        real(REAL_TYPE), dimension(:), &
          allocatable :: PhaseDisplacementSoil ! Nodal total displacement of the load phase
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: IncrementalDisplacementSoil ! Nodal incremental displacement
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: AccumulatedIncDisplacementSoil ! Nodal incremental displacement
        real(REAL_TYPE), dimension(:),  &
          allocatable :: AccumulatedDisplacementSoil ! Accumulated displacements of a load phase
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: AccelerationSoil  ! Nodal acceleration of soil
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: RateofMomentum  ! Rate of momentum
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ActiveElement ! ID of the active elements (rank 2 for number of patches)
        real(REAL_TYPE), dimension(:), &
          allocatable :: TotalDisplacementWater ! Nodal total displacement of water (FEM only)
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: TotalVelocityWater ! Nodal total velocity for water
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: TotalVelocityWaterPrevious ! Nodal total velocity for water for the previous time step
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: IncrementalDisplacementWater ! Nodal incremental displacement
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: IncrementalDisplacementWaterPrevious ! Nodal incremental displacement at previous time step
        real(REAL_TYPE), dimension(:),  &
          allocatable :: AccumulatedDisplacementWater ! Accumulated displacements of a load phase
        real(REAL_TYPE), dimension(:), &
          allocatable :: PhaseDisplacementWater ! Nodal total displacement of the load phase
        real(REAL_TYPE), dimension(:), &
          allocatable :: TotalDisplacementGas ! Nodal total displacement of gas (FEM only)
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: TotalVelocityGas ! Nodal total velocity for gas
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: NonAdvectiveFluxAirInWater ! Nodal Non advective flux of Air in the Fluid 
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: NonAdvectiveFluxVapourInGas ! Nodal Non advective flux of Vapour in the Gas 
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: AdvectiveFluxDarcyWater ! Nodal advective flux of water 
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: AdvectiveFluxDarcyAir ! Nodal advective flux of air  
        real(REAL_TYPE), dimension(:, :), &
        allocatable :: ThermalConductionFlux ! Nodal Heat Conduction flux
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: FReaction ! Nodal reactions on a user-defined mesh boundary
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: FReactionWater
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: FReactionGas
        real(REAL_TYPE), dimension(:), allocatable :: RateVolStrain
        real(REAL_TYPE), dimension(:, :), &
          allocatable :: NodalUnitMassGradient
        real(REAL_TYPE), dimension(:,:), &
            allocatable :: SoilSurfaceNodeCoordMatrix
        real(REAL_TYPE), dimension(:,:), &
            allocatable :: PhreaticSurfaceNodeCoordMatrix
        real(REAL_TYPE), dimension(:,:), &
            allocatable :: HydrHeadArea
        real(REAL_TYPE), dimension(:,:), &
            allocatable :: SeepageArea
        real(REAL_TYPE), dimension(:,:), &
            allocatable :: InfiltrationArea
        real(REAL_TYPE), dimension(:), allocatable :: InfiltrationRate
    

        !!> stores material points indeces
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: GetParticleIndex

        !> stores determinant of elements
        real(REAL_TYPE), dimension(:, :), allocatable :: ElementDeterminant
        
        ! quasi static additions
        real(REAL_TYPE), dimension(:), allocatable :: ExtFlow ! mass balance steady flow
        real(REAL_TYPE), dimension(:), allocatable :: IntFlow ! correction term
        real(REAL_TYPE), dimension(:), allocatable :: RateOfFlux ! flow residual
        real(REAL_TYPE), dimension(:), allocatable :: IncrementalPressure ! incremental water pressure
        real(REAL_TYPE), dimension(:), allocatable :: TotalPressure ! total water pressure
        real(REAL_TYPE), dimension(:), allocatable :: SubIncrementalDisplacement
        real(REAL_TYPE), dimension(:), allocatable :: SubIncrementalPressure
        
        integer(INTEGER_TYPE), dimension(:), allocatable :: ConsideredElemReaction !vector containing the ElementID of those elements considered in the integration of reaction forces
        logical, dimension(:,:), allocatable :: IsReactionNodeSurface
        character (len=255), dimension(:), allocatable :: OutputSurfaceName
        logical, dimension(:), allocatable :: IsReactionNode!Determine if a node is on the reaction surface
        
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPGlobalPositionElement		
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! NURBS related global variables 
        ! Multipatch variables 
        !integer(INTEGER_TYPE) :: Counters%NPatches                        = -1 ! number of patches in the multipatch generalization
        ! Multipatch --> maximum allocation values
        integer(INTEGER_TYPE) :: Maximum_NControlPoints ! temporarily initialized here
        integer(INTEGER_TYPE) :: Maximum_NXiKnotEntries ! temporarily initialized here
        integer(INTEGER_TYPE) :: Maximum_NEtaKnotEntries ! temporarily initialized here
        integer(INTEGER_TYPE) :: Maximum_NZetaKnotEntries ! temporarily initialized here
          
        
        integer(INTEGER_TYPE) :: Maximum_NXiKnotOrder
        integer(INTEGER_TYPE) :: Maximum_NEtaKnotOrder
        integer(INTEGER_TYPE) :: Maximum_NZetaKnotOrder
        
        ! Number of control points 
        integer(INTEGER_TYPE), allocatable, dimension(:) :: NControlPoints
        integer(INTEGER_TYPE) :: NControlPoints_Traction ! allocatable,
        
        ! Xi direction 
        integer(INTEGER_TYPE), allocatable, dimension(:) :: NXiKnotEntries                           != -1 ! Xi vector number of entries
        integer(INTEGER_TYPE), allocatable, dimension(:) :: NXiKnotOrder                             != -1 ! Xi vector order        
        integer(INTEGER_TYPE), allocatable, dimension(:) :: nn_NURBS_NumberOfUnivariateXiKnots       != -1 ! number of univariate basis functions for eta direction 
        integer(INTEGER_TYPE) :: NXiGaussPoints                           = -1 ! number of gauss points in xi direction 
        ! I added a rank for multipatch generalization
        real(REAL_TYPE), allocatable, dimension(:, :) :: XiKnotEntries            ! Xi knot vector 
        real(REAL_TYPE), allocatable, dimension(:) :: XiGaussPoints            ! Xi gauss points 

        
        ! ----traction
        integer(INTEGER_TYPE):: NXiKnotEntries_Traction                           != -1 ! Xi vector number of entries
        !, allocatable, dimension(:)
        integer(INTEGER_TYPE):: NXiKnotOrder_Traction                             != -1 ! Xi vector order        
        !, allocatable , dimension(:)
        integer(INTEGER_TYPE) :: nn_NURBS_NumberOfUnivariateXiKnots_Traction       != -1 ! number of univariate basis functions for eta direction 
        !, allocatable , dimension(:)
        real(REAL_TYPE), allocatable, dimension(:) :: XiKnotEntries_Traction            ! Xi knot vector 
                                                !, :
        
        
        ! Eta direction 
        integer(INTEGER_TYPE), allocatable, dimension(:) :: NEtaKnotEntries                          != -1 ! Eta vector number of entries
        integer(INTEGER_TYPE), allocatable, dimension(:) :: NEtaKnotOrder                            != -1 ! Eta vector order
        integer(INTEGER_TYPE), allocatable, dimension(:) :: mm_NURBS_NumberOfUnivariateEtaKnots      != -1 ! number of univariate basis functions for eta direction 
        integer(INTEGER_TYPE) :: NEtaGaussPoints                          = -1 ! number of gauss points in eta direction 
        ! I added a rank for multipatch generalization
        real(REAL_TYPE), allocatable, dimension(:, :) :: EtaKnotEntries           ! Eta knot vector 
        real(REAL_TYPE), allocatable, dimension(:) :: EtaGaussPoints           ! Eta gauss points 

        !------traction
        integer(INTEGER_TYPE) :: NEtaKnotEntries_Traction !, allocatable                           != -1 ! Xi vector number of entries
        !, dimension(:)
        integer(INTEGER_TYPE) :: NEtaKnotOrder_Traction !, allocatable                             != -1 ! Xi vector order        
        !, dimension(:)
        integer(INTEGER_TYPE) :: mm_NURBS_NumberOfUnivariateEtaKnots_Traction !, allocatable       != -1 ! number of univariate basis functions for eta direction 
        !, dimension(:)
        real(REAL_TYPE), allocatable, dimension(:) :: EtaKnotEntries_Traction            ! Xi knot vector 
                                                !, :
        
        
        ! Zeta direction 
        integer(INTEGER_TYPE), allocatable, dimension(:) :: NZetaKnotEntries                         != -1 ! Zeta vector number of entries
        integer(INTEGER_TYPE), allocatable, dimension(:) :: NZetaKnotOrder                          != -1 ! Zeta vector order
        integer(INTEGER_TYPE), allocatable, dimension(:) :: oo_NURBS_NumberOfUnivariateZetaKnots     != -1 ! number of univariate basis functions for zeta direction 
        integer(INTEGER_TYPE) :: NZetaGaussPoints                         = -1 ! number of gauss points in zeta direction 
        ! I added a rank for multipatch generalization
        real(REAL_TYPE), allocatable, dimension(:, :) :: ZetaKnotEntries          ! Zeta knot vector           
        real(REAL_TYPE), allocatable, dimension(:) :: ZetaGaussPoints          ! Zeta knot vector           

        ! Knot connectivity array (multipatch rank)
        integer(INTEGER_TYPE), allocatable, dimension(:,:, :) :: INN ! NURBS coordinates telling you (early) knots that define an element 
        integer(INTEGER_TYPE), allocatable, dimension(:,:, :) :: IEN ! Knot connectivity
        
        !------- traction
        integer(INTEGER_TYPE), allocatable, dimension(:,:) :: INN_Traction ! NURBS coordinates telling you (early) knots that define an element 
        integer(INTEGER_TYPE), allocatable, dimension(:,:) :: IEN_Traction ! Knot connectivity
        
        ! I moved this here so that I can access in the shape function subroutine to find the weights
        integer(INTEGER_TYPE), dimension(:,:, :), allocatable :: ElementConnectivities
        
        ! Miscellaneous (multipatch generalization rank for nel_NURBS, nnp_NURBS, nen_NURBS)
        integer(INTEGER_TYPE), allocatable, dimension(:) :: nel_NURBS                        ! number of elements 
        
        integer(INTEGER_TYPE), allocatable, dimension(:) :: nnp_NURBS                        ! number of global basis functions 
                                                                                ! (global here refers to its global domain within the 'super' element)
        integer(INTEGER_TYPE), allocatable, dimension(:) :: nen_NURBS                        ! number of local basis functions
                                                                                ! (local here refers to a knot span i.e. accross one single element)  
        
        !------traction
        integer(INTEGER_TYPE) :: nel_NURBS_Traction                        ! number of elements 
        
        !, allocatable, dimension(:)
        integer(INTEGER_TYPE) :: nnp_NURBS_Traction                        ! number of global basis functions 
                                                                                ! (global here refers to its global domain within the 'super' element)
        !, allocatable, dimension(:) 
        integer(INTEGER_TYPE):: nen_NURBS_Traction                        ! number of local basis functions
                                                                                ! (local here refers to a knot span i.e. accross one single element)  
        !, allocatable, dimension(:) 
        
        ! Number of active elements
        integer(INTEGER_TYPE), allocatable, dimension(:) :: nael_NURBS
        
        integer(INTEGER_TYPE) :: Maximum_nel_NURBS
        integer(INTEGER_TYPE) :: Maximum_nnp_NURBS
        integer(INTEGER_TYPE) :: Maximum_nen_NURBS
        
        integer(INTEGER_TYPE) :: Sum_nel_NURBS
        integer(INTEGER_TYPE) :: Sum_nael_NURBS
        integer(INTEGER_TYPE) :: Sum_NControlPoints
        
        !integer(INTEGER_TYPE), allocatable, dimension(:) :: NumberOfGlobalControlPointsUniqueMultipatch
               
        integer(INTEGER_TYPE) :: NumberOfGlobalControlPointsUniqueMultipatch

        
        integer(INTEGER_TYPE), allocatable, dimension(:) :: NumberOfUniqueNodes
        !Sum_el_NURBS = 0
        
        real(REAL_TYPE), allocatable :: NN_IncludesZeroValues_Print(:,:) !:,  --> I removed the rank associated with the gauss point number... always one..
        real(REAL_TYPE), allocatable :: dN_dxi_IncludesZeroValues_Print(:,:) !:,  --> I removed the rank associated with the gauss point number... always one..
        
        real(REAL_TYPE), allocatable :: MM_IncludesZeroValues_Print(:,:) !:, --> I removed the rank associated with the gauss point number... always one..
        real(REAL_TYPE), allocatable :: dM_deta_IncludesZeroValues_Print(:,:) !:, --> I removed the rank associated with the gauss point number... always one..
        
        real(REAL_TYPE), allocatable :: OO_IncludesZeroValues_Print(:,:) !:, --> I removed the rank associated with the gauss point number... always one..
        real(REAL_TYPE), allocatable :: dO_dzeta_IncludesZeroValues_Print(:,:) !:,  --> I removed the rank associated with the gauss point number... always one..
        
        
        ! Maximum particle previous 
        integer(INTEGER_TYPE) :: maxParticle_previous =0
        integer(INTEGER_TYPE) :: GetParticleIndexCounter        
        logical :: FirstTimeGetParticleIndex 
        ! I added this here from MeshInfo.FOR
        integer(INTEGER_TYPE), dimension(:,:, :), allocatable :: GetParticleIndex ! (rank 3 for multipatch implementation)
        
        
        ! Pointers for NURBS 
        !procedure(DummyCheckForGlobPosPointer), pointer:: CheckForGlobPosPointer => NULL()
        !procedure(DummyGauss_Q1Pointer), pointer:: Gauss_Q1Pointer => NULL()
        procedure(DummyInitialiseShapeFunctionsPointer_NURBS_2D), pointer:: InitialiseShapeFunctionsPointer_NURBS => NULL() ! 2D
        procedure(DummyInitialiseShapeFunctionsPointer_NURBS_3D), pointer:: InitialiseShapeFunctionsPointer_NURBS_3D => NULL() ! 3D
        
        !procedure(DummyInitialiseShapeFunctionsPointer_NURBS), pointer:: InitialiseShapeFunctionsBoundaryPointer_NURBS => NULL() ! 1D
        !procedure(DummyIsInsideElementLocPosPointer), pointer:: IsInsideElementLocPosPointer => NULL()
        !procedure(DummyGetMinAltitudePointer), pointer:: GetMinAltitudePointer => NULL()
        !procedure(DummyInitialLocalMaterialPointCoordinatesPointer), pointer :: InitialLocalMaterialPointCoordinatesPointer_NURBS => NULL()
        !procedure(DummyShapeLocPosPointer), pointer:: ShapeLocPosPointer => NULL() 
        !procedure(DummyRearrangeConnectivitiesPointer), pointer:: RearrangeConnectivitiesPointer => NULL() 
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NURBS
        

        ! ----------------------------------- 1D                                                    

        
        subroutine DummyInitialiseShapeFunctionsPointer_NURBS_2D(HS, dHS, Wt, &
                                                              HS_Xi, dHS_Xi, Wt_Xi, &
                                                              HS_Eta, dHS_Eta, Wt_Eta, &
                                                              XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, &
                                                              EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
                                                              ni, nj, &
                                                              IPatch)
           implicit none
           real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HS 
           real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: dHS 
           real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt ! this should be read with the nodes... size = number of nodes 


           real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HS_Xi 
           real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: dHS_Xi 
           real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt_Xi ! this should be read with the nodes... size = number of nodes
           
           
           real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HS_Eta
           real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: dHS_Eta 
           real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt_Eta ! this should be read with the nodes... size = number of nodes
           
           
           integer(INTEGER_TYPE), intent(inout) :: ni, nj
           
  
           !NURBS related inputs in the xi direction 
           integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
           integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
           real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
           real(REAL_TYPE), intent(in) :: Xi_ParametricDomain !, dimension(NXiGaussPoints)
          
           !NURBS related inputs in the eta direction 
           integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder
           integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntries
           real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: EtaKnotEntries
           real(REAL_TYPE), intent(in) :: Eta_ParametricDomain !, dimension(NXiGaussPoints)
           !
           !InitialiseShapeFunctionsLINE2_NURBS (HS, dHS, Wt, &               !classic inout parameters
           !                                             XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder)
           
           ! Multipatch variable
           integer(INTEGER_TYPE), intent(in) :: IPatch
           
           
           ! NURBS specific parameters 
           !integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
           !real(REAL_TYPE), intent(in), dimension(NXiGaussPoints) :: Xi_ParametricDomain 
           !integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder 
           !real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
           
           ! what other parameters do I need here 
           !xi 
           !order 
           !control point weights
           
                                                              
                                                              
        
                                                              end subroutine DummyInitialiseShapeFunctionsPointer_NURBS_2D
                                                              
           
                                                              
                                                              
                                                              
                                                              
                                                              
                                                              
          ! ----------------------------------- 2D          
                                                              
          ! ----------------------------------- 3D
                                                              
                                                              
          subroutine DummyInitialiseShapeFunctionsPointer_NURBS_3D(HS, dHS, Wt, &
                                                              HS_Xi, dHS_Xi, Wt_Xi, &
                                                              HS_Eta, dHS_Eta, Wt_Eta, &
                                                              HS_Zeta, dHS_Zeta, Wt_Zeta, &
                                                              XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, &
                                                              EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
                                                              ZetaKnotEntries, NZetaKnotEntries, Zeta_ParametricDomain, NZetaKnotOrder, &
                                                              NumberOfGaussPoints_X, NumberOfGaussPoints_Y, NumberOfGaussPoints_Z, &
                                                              ni, nj, nk, &
                                                              IElement, &
                                                              IPatch)
           implicit none
           real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HS 
           real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: dHS 
           real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt ! this should be read with the nodes... size = number of nodes 


           real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HS_Xi 
           real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: dHS_Xi 
           real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt_Xi ! this should be read with the nodes... size = number of nodes
           
           
           real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HS_Eta
           real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: dHS_Eta 
           real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt_Eta ! this should be read with the nodes... size = number of nodes
           
           
           real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: HS_Zeta
           real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: dHS_Zeta 
           real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt_Zeta ! this should be read with the nodes... size = number of nodes
           
           
           integer(INTEGER_TYPE), intent(inout) :: ni, nj, nk
           
           
           integer(INTEGER_TYPE), intent(in) :: IElement
           
           integer(INTEGER_TYPE), intent(in) :: IPatch
  
           !NURBS related inputs in the xi direction 
           integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
           integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
           real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
           real(REAL_TYPE), intent(in):: Xi_ParametricDomain !, dimension(NXiGaussPoints) 
          
           !NURBS related inputs in the eta direction 
           integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder
           integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntries
           real(REAL_TYPE), dimension(NEtaKnotEntries), intent(in) :: EtaKnotEntries
           real(REAL_TYPE), intent(in) :: Eta_ParametricDomain !, dimension(NXiGaussPoints)
           
           !NURBS related inputs in the Zeta direction 
           integer(INTEGER_TYPE), intent(in) :: NZetaKnotOrder
           integer(INTEGER_TYPE), intent(in) :: NZetaKnotEntries
           real(REAL_TYPE), dimension(NZetaKnotEntries), intent(in) :: ZetaKnotEntries
           real(REAL_TYPE), intent(in) :: Zeta_ParametricDomain !, dimension(NXiGaussPoints) 
           
           integer(INTEGER_TYPE), intent(in) :: NumberOfGaussPoints_X
           integer(INTEGER_TYPE), intent(in) :: NumberOfGaussPoints_Y
           integer(INTEGER_TYPE), intent(in) :: NumberOfGaussPoints_Z
           
           !
           !InitialiseShapeFunctionsLINE2_NURBS (HS, dHS, Wt, &               !classic inout parameters
           !                                             XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder)
           
           
           
           ! NURBS specific parameters 
           !integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
           !real(REAL_TYPE), intent(in), dimension(NXiGaussPoints) :: Xi_ParametricDomain 
           !integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder 
           !real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
           
           ! what other parameters do I need here 
           !xi 
           !order 
           !control point weights
           
                                                              
          end subroutine DummyInitialiseShapeFunctionsPointer_NURBS_3D     
                                                              
                                                              
        !                                                      
        !  subroutine DummyInitialiseShapeFunctionsPointer_NURBS(HS, dHS, Wt, &
        !                                                      XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, &
        !                                                      EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder)
        !   implicit none
        !   real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: HS 
        !   real(REAL_TYPE), allocatable, dimension(:, :), intent(inout) :: dHS 
        !   real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt !  ! this should be read with the nodes... size = number of nodes 
        !   
        !   
        !   !
        !   !InitialiseShapeFunctionsLINE2_NURBS (HS, dHS, Wt, &               !classic inout parameters
        !   !                                             XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder)
        !   
        !   
        !   
        !   ! NURBS specific parameters 
        !   ! xi
        !   integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
        !   real(REAL_TYPE), intent(in), dimension(NXiGaussPoints) :: Xi_ParametricDomain 
        !   integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder 
        !   real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
        !   ! eta
        !   integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntries
        !   real(REAL_TYPE), intent(in), dimension(NEtaGaussPoints) :: Eta_ParametricDomain 
        !   integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder 
        !   real(REAL_TYPE), dimension(NEtaKnotEntries), intent(in) :: EtaKnotEntries
        !   
        !   ! what other parameters do I need here 
        !   !xi 
        !   !order 
        !   !control point weights
        !   
        !end subroutine DummyInitialiseShapeFunctionsPointer_NURBS                                                  
                                                              
      
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NURBS
        
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
        
        ! Multipatch variable
        integer(INTEGER_TYPE) :: IPatch_Temporary = 1
                
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
            
          case(QUAD4) ! 'quadrilateral_4-noded'
              ELEMENTNODES = 4
              ELEMENTVERTICES = 4
              ELEMENTGAUSSPOINTS = 1 !I should be able to have as many gauss points as I want... 
                                     !why do we have this difference between gauss point and material points
              ELEMENTBOUNDARYNODES = 2
              ELEMENTBOUNDARYGAUSSPOINTS = 1
              ELEMENTSIDES = 4
              VTK_CELL = 9 ! VTK_QUAD(=9)
              ELEMENTTYPE = QUAD4
              
              
          ! Currently working on this for the implementation of NURBS   
          case(QUAD4_NURBS) ! 'quadrilateral_4-noded_NURBS'
              
              ! consider generalizing this to multiple gauss points -> why is there a difference between ELEMENTGAUSSPOINTS and ELEMENTBOUNDARYGAUSSPOINTS
              ELEMENTBOUNDARYGAUSSPOINTS_XI = NXiGaussPoints 
              ELEMENTBOUNDARYGAUSSPOINTS_ETA = NEtaGaussPoints
              !ELEMENTBOUNDARYGAUSSPOINTS_ZETA = 4!NZetaGaussPoints
              ELEMENTBOUNDARYGAUSSPOINTS = ELEMENTBOUNDARYGAUSSPOINTS_XI
              ELEMENTGAUSSPOINTS = ELEMENTBOUNDARYGAUSSPOINTS_XI *ELEMENTBOUNDARYGAUSSPOINTS_ETA
              
              ELEMENTBOUNDARYNODES_XI = NXiKnotOrder(IPatch_Temporary)+1!NXiKnotOrder_Traction+1 ! Xi  !(IPatch_Temporary)
              ELEMENTBOUNDARYNODES_ETA = NEtaKnotOrder(IPatch_Temporary)+1!NEtaKnotOrder_Traction+1 ! Eta !(IPatch_Temporary)
              !ELEMENTBOUNDARYNODES_ZETA = NZetaKnotOrder(IPatch_Temporary)+1 ! Zeta
              
              ELEMENTBOUNDARYNODES = (NXiKnotOrder_Traction+1)
              
              ELEMENTNODES_XI = (NDIM * NXiKnotOrder(IPatch_Temporary)) 
              ELEMENTNODES_ETA = (NDIM * NEtaKnotOrder(IPatch_Temporary))
              
              ELEMENTNODES = ELEMENTBOUNDARYNODES_XI * ELEMENTBOUNDARYNODES_ETA ! if 2GP by 2GP then 4GPs per element
              
              
              ELEMENTVERTICES_XI = (NDIM * NXiKnotOrder(IPatch_Temporary)) 
              ELEMENTVERTICES_ETA = (NDIM * NEtaKnotOrder(IPatch_Temporary))
              ELEMENTVERTICES = (ELEMENTBOUNDARYNODES_XI-1) * (ELEMENTBOUNDARYNODES_ETA-1)
              
              ! consider generalizing this to multiple gauss points
              !ELEMENTGAUSSPOINTS_XI = NXiGaussPoints ! Xi - this should be read somewhere  
              !ELEMENTGAUSSPOINTS_ETA = NEtaGaussPoints ! Eta
              !ELEMENTGAUSSPOINTS_ZETA = NZetaGaussPoints ! Zeta              
              !ELEMENTGAUSSPOINTS = ELEMENTGAUSSPOINTS_XI * ELEMENTGAUSSPOINTS_ETA ! 2D implementation 
              
              


              !I should be able to have as many gauss points as I want... why do we have this difference between gauss point and material points
              
              
              
              
              
              
              
              ELEMENTSIDES = 4 ! -> not sure if this needs to change because of the NURBS implementation
              VTK_CELL = 9 ! VTK_QUAD(=9)
              ELEMENTTYPE = QUAD4_NURBS
              
              
              
              
              NumberOfControlNetGrid = -1 ! this is read from the GOM file under the $$CONTROL_NET_CONNECTIVITY flag 
              ControlNetGrid = 4 ! it is a structured square net with 4 nodes. Not sure if this is hardcoded or always the case but I think the latter. 
            
              
              VERTICESINPARAMETERSPACE = 4
              
              
              
              
              
              
              
              
              
          case(HEXA_NURBS) ! 'hexahedral_NURBS'
              ! note that the below need to be calculated based on the order inputs from the GOM file 
              !ELEMENTNODES = 8
              !ELEMENTVERTICES = 8
              !ELEMENTGAUSSPOINTS = 1
              
              
              ELEMENTBOUNDARYGAUSSPOINTS_XI = NXiGaussPoints !4!4!NXiGaussPoints 
              ELEMENTBOUNDARYGAUSSPOINTS_ETA = NEtaGaussPoints!4!1!4!NEtaGaussPoints
              ELEMENTBOUNDARYGAUSSPOINTS_ZETA = NEtaGaussPoints!4!1!4!NZetaGaussPoints !---> this should be zeta??
              
              
              
              
              ELEMENTBOUNDARYNODES = (NXiKnotOrder_Traction+1)*(NEtaKnotOrder_Traction+1) !(IPatch_Temporary)*!4
              ELEMENTBOUNDARYGAUSSPOINTS = ELEMENTBOUNDARYGAUSSPOINTS_XI*ELEMENTBOUNDARYGAUSSPOINTS_ETA !ELEMENTBOUNDARYGAUSSPOINTS_XI*ELEMENTBOUNDARYGAUSSPOINTS_ETA!4*4 !1 --> 2D surface
              ELEMENTGAUSSPOINTS = ELEMENTBOUNDARYGAUSSPOINTS_XI*ELEMENTBOUNDARYGAUSSPOINTS_ETA*ELEMENTBOUNDARYGAUSSPOINTS_ZETA !4*4*4
              !ELEMENTSIDES = 6
              !VTK_CELL = 12 ! VTK_HEXAHEDRON(=12)
              !ELEMENTTYPE = HEXA_NURBS
              
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! need to remove above 
              
              !ELEMENTNODES_XI = NXiKnotOrder_Traction+1!(NDIM * NXiKnotOrder(IPatch_Temporary)) 
              !ELEMENTNODES_ETA = NEtaKnotOrder_Traction+1!(NDIM * NEtaKnotOrder(IPatch_Temporary))
              !ELEMENTNODES_ZETA = NZetaKnotOrder_Traction+1!(NDIM * NZetaKnotOrder(IPatch_Temporary))
              
              ELEMENTVERTICES_XI = NXiKnotOrder(IPatch_Temporary)+1!(NDIM * NXiKnotOrder(IPatch_Temporary)) 
              ELEMENTVERTICES_ETA = NEtaKnotOrder(IPatch_Temporary)+1!(NDIM * NEtaKnotOrder(IPatch_Temporary))
              !ELEMENTVERTICES_ZETA = NZetaKnotOrder_Traction+1!(NDIM * NZetaKnotOrder(IPatch_Temporary))
              ! consider generalizing this to multiple gauss points
              
              !ELEMENTGAUSSPOINTS = 4*4*4!ELEMENTBOUNDARYGAUSSPOINTS_XI*ELEMENTBOUNDARYGAUSSPOINTS_ETA*ELEMENTBOUNDARYGAUSSPOINTS_ZETA!4*4*4 --> 3D volume
              
              !ELEMENTGAUSSPOINTS_XI = 4!NXiGaussPoints ! Xi - this should be read somewhere  
              !ELEMENTGAUSSPOINTS_ETA = 4!NEtaGaussPoints ! Eta
              !ELEMENTGAUSSPOINTS_ZETA = 4!NZetaGaussPoints ! Zeta              
              
              !ELEMENTGAUSSPOINTS = 343!1!64 !hardcoded for now as I am now looking at the accuracy in determining the element space calculation for NURBS irregular geometries
              !8 ! 2*2*2 --> ELEMENTGAUSSPOINTS_XI * ELEMENTGAUSSPOINTS_ETA * ELEMENTGAUSSPOINTS_ZETA ! 2D implementation 
              
              


              !I should be able to have as many gauss points as I want... why do we have this difference between gauss point and material points
              
              ELEMENTBOUNDARYNODES_XI = NXiKnotOrder(IPatch_Temporary)+1!NXiKnotOrder(IPatch_Temporary)+1 ! Xi  !_Traction
              ELEMENTBOUNDARYNODES_ETA = NEtaKnotOrder(IPatch_Temporary)+1!(IPatch_Temporary)+1 ! Eta !_Traction
              ELEMENTBOUNDARYNODES_ZETA = NZetaKnotOrder(IPatch_Temporary)+1
              !ELEMENTBOUNDARYNODES_ZETA = NZetaKnotOrder_Traction+1!(IPatch_Temporary)+1 ! Zeta
              
              ! consider generalizing this to multiple gauss points -> why is there a difference between ELEMENTGAUSSPOINTS and ELEMENTBOUNDARYGAUSSPOINTS
              ! HARDCODED: I COMMENTED THIS FOR NOW 
              
              
              ELEMENTNODES = ELEMENTBOUNDARYNODES_XI * ELEMENTBOUNDARYNODES_ETA * ELEMENTBOUNDARYNODES_ZETA ! if 2GP by 2GP then 4GPs per element
              ELEMENTVERTICES = (ELEMENTBOUNDARYNODES_XI-1) * (ELEMENTBOUNDARYNODES_ETA-1) * (ELEMENTBOUNDARYNODES_ZETA-1)
              
              ELEMENTSIDES = 6!6 ! -> not sure if this needs to change because of the NURBS implementation
              ELEMENTSURFACES = 6
              VTK_CELL = 12 ! VTK_HEXA(=2)
              ELEMENTTYPE = HEXA_NURBS
              
              !ELEMENTBOUNDARYNODES = 8
              
              
              !LowestPointOnASurfaceDefiningPhysicalDomain
              !
              !HighestPointOnASurfaceDefiningPhysicalDomain
              
              
              NumberOfControlNetGrid = -1 ! this is read from the GOM file under the $$CONTROL_NET_CONNECTIVITY flag 
              ControlNetGrid = 8 ! it is a structured square net with 4 nodes. Not sure if this is hardcoded or always the case but I think the latter. 
              
              VERTICESINPARAMETERSPACE = 8!4
        !      
        !      !-----------------------------------------------------------------------------------------------------------------------------------------------
        !      real(REAL_TYPE), dimension(:), allocatable :: GPWeight 
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative
        !
        !
        !! Originaly in Anura3D, GPShapeFunctionDerivative is the same for all elements. 
        !! So we do not need to calculate it for each element. But in NURBS, we have a super element
        !! so we need to know for what element we are doing this for. We can add a dimension with
        !! element number
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPWeight_AllElements  
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunction_AllElements
        !real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivative_AllElements
        !
        !! NURBS xi direction  
        !real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary ! for each boundary (1D) there are a number of weights 
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary ! similarly here for shape function evaluation (1D) 
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary ! similarly here for shape function derivative (1D) 
        !
        !!NURBS eta direction 
        !!real(REAL_TYPE), dimension(:), allocatable :: GPWeight_ETA 
        !!real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction_ETA
        !!real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative_ETA
        !
        !real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary_ETA
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary_ETA
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ETA
        !
        !!NURBS zeta direction 
        !!real(REAL_TYPE), dimension(:), allocatable :: GPWeight_ZETA 
        !!real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction_ZETA
        !!real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative_ZETA
        !
        !real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary_ZETA
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary_ZETA
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ZETA
        !
        !
        !
        !! ------------------ all elements needed for traction implementation -------------------------------------------------
        !
        !
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_XI_AllElements
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_XI_AllElements
        !real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_XI_AllElements
        !
        !
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_ETA_AllElements
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_ETA_AllElements
        !real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ETA_AllElements
        !
        !
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_ZETA_AllElements
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_ZETA_AllElements
        !real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ZETA_AllElements
        !      
        !      
        !      
        !      
        !      !-------------------------------------------------------------------------------------------------------------------------------------------------------
        !      
        !      
        !      
        !      
        !      
              
              
              
              
              
              
              
              
              
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

