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
	  
	  
	  module ModMeshInfo
      !**********************************************************************
      !
      !    Function:  Contains routines for initialising mesh data by allocating arrays
      !
      !    Implemented in the frame of the MPM project.
      !
      !    Should this be a module 
      !
      !     $Revision: 8842 $
      !     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
      !
      !**********************************************************************
    
      use ModFileIO
      use ModCounters
      use ModReadCalculationData
      use ModReadMaterialData
      use ModElementEvaluation
      use ModGlobalConstants
      !use ModFileIO
      use ModNURBS
      !use ModMPMInit
      !use ModMPMData
      
      !External :: FExist
      
      implicit none

      
        !The below variables were moved to GlobalConstant.FOR (at the beginning of the module). 
        !This is because they are allocated there because that's where we read the information and 
        !loop accross patched in the GOM file. 
        
        !! Nodal coordinates, updated nodal coordinates, The initial coordinates
        !! Note that these are also considered as control points in the context of NURBS
        !real(REAL_TYPE), dimension(:, :),  &
        !  allocatable :: NodalCoordinates, NodalCoordinatesUpd, NodalOriginalCoord
        !!real(REAL_TYPE), dimension(:), &
        !!    allocatable :: ControlPoint_Weights ! NURBS weights - each node has a weight
        !!real(REAL_TYPE), dimension(:, :),  &
        !!  allocatable :: XiKnotEntries, EtaKnotEntries  
        !real(REAL_TYPE), dimension(:, :),  &
        !  allocatable :: NodalPrescibedDisp ! Nodal prescribed displacement in 3 directions (Solid)
        !real(REAL_TYPE), dimension(:, :),  &
        !  allocatable :: NodalPrescibedDispWater ! Nodal prescribed displacement in 3 directions (Water)
        !real(REAL_TYPE), dimension(:, :),  &
        !  allocatable :: NodalPrescibedDispGas ! Nodal prescribed displacement in 3 directions (Gas)
        !real(REAL_TYPE), dimension(:, :, :), &
        !  allocatable :: LoadValuesOnNodesSolid ! Load values (applied on nodes) for each distributed element face (solid)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnNodesWater ! Load values (applied on nodes) for each distributed element face (water)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnNodesGas ! Load values (applied on nodes) for each distributed element face (gas)
        !real(REAL_TYPE), dimension(:, :, :), &
        !  allocatable :: LoadValuesOnNodesSolidB ! Load values (applied on nodes) for each distributed element face (solid, load system B)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnNodesWaterB ! Load values (applied on nodes) for each distributed element face (water)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnNodesGasB ! Load values (applied on nodes) for each distributed element face (gas)
        !  real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: HydraulicHeadLoad ! Load values (applied on nodes) for each distributed element face (hydraulic head)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnMatPointsSolid ! Load values (applied on material points) for each distributed element face (solid)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnMatPointsWater ! Load values (applied on material points) for each distributed element face (water)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnMatPointsGas ! Load values (applied on material points) for each distributed element face (gas)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnMatPointsSolidB ! Load values (applied on material points) for each distributed element face (solid)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnMatPointsWaterB ! Load values (applied on material points) for each distributed element face (water)
        !real(REAL_TYPE), dimension(:, :, :),  &
        !  allocatable :: LoadValuesOnMatPointsGasB ! Load values (applied on material points) for each distributed element face (gas)
        !! Prescribed boundary condition for the solid phase (1 for free, 0 for fixed)
        !real(REAL_TYPE), dimension(:), allocatable :: PBoundary
        !! 4 degrees of freedom
        !integer(INTEGER_TYPE), dimension(:), allocatable :: PBoundaryQuasiStatic 
        !! Prescribed boundary condition for the water phase (1 for free, 0 for fixed)
        !real(REAL_TYPE), dimension(:), allocatable :: PBoundaryWater
        !! Prescribed boundary condition for the gas phase (1 for free, 0 for fixed)
        !real(REAL_TYPE), dimension(:), allocatable :: PBoundaryGas
        !! Nodal load vector of external load (solid)
        !real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVector
        !! Nodal load vector of external load (solid)
        !real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorB
        !! Nodal load vector of external load (water)
        !real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorWater
        !! Nodal load vector of external load (water)
        !real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorWaterB
        !! Nodal load vector of external load (gas)
        !real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorGas
        !! Nodal load vector of external load (gas)
        !real(REAL_TYPE), dimension(:), allocatable :: ExtLoadVectorGasB
        !! Nodal load vector of external load (HydraulicHead)
        !real(REAL_TYPE), dimension(:), allocatable :: HydraulicHeadVector
        !! Element connectivities (4 node tetrahedral element)
        !!integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElementConnectivities
        !
        !
        !
        !
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: ControlNetConnectivity
        !
        !
        !
        !! Element connectivities for 10-noded tetrahedral element
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElementConnectivities10Node
        !! Ties edge nodes of high-order tetrahedral elements to corner nodes
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: EdgeNodeTyingsHOE
        !! Distributed element face load connectivity (6-noded triangular element)(solid)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesSolid
        !! Distributed element face load connectivity (6-noded triangular element)(solid, load system B)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesSolidB
        !! Distributed element face load connectivity (6-noded triangular element)(water)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesWater
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesWaterB
        !! Distributed element face load connectivity (6-noded triangular element)(gas)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesGas
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnNodesConnectivitiesGasB
        !!On the soil surface element face connectivity (6-noded triangular element)(soil surface for K0)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: SoilSurfaceNodesConnectivities
        !!On the phratic surface element face connectivity (6-noded triangular element)(PHREATIC surface for K0)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: PhreaticSurfaceNodesConnectivities
        !! Distributed element face load connectivity (6-noded triangular element)(solid)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesSolid
        ! integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesSolidB
        !! Distributed element face load connectivity (6-noded triangular element)(water)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesWater
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesWaterB
        !! Distributed element face load connectivity (6-noded triangular element)(gas)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesGas
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: LoadOnMatPointsConnectivitiesGasB
        !! Distributed element face load connectivity (6-noded triangular element)(hydraulic head)
        !integer(INTEGER_TYPE), dimension(:, :), allocatable :: HydraulicHeadNodesConnectivities
        !integer(INTEGER_TYPE), dimension(:), allocatable :: ElementMaterialID ! Material ID for the element
        !integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDof ! The global storage of degree of freedom in x-direction
        !integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDofQuasiStatic ! The global storage of degree of freedom in ux,uy,uz,p
        !logical, dimension(:), allocatable :: IsActiveElement ! Element switch, 1:Active element and 0:Inactive element
        !real(REAL_TYPE), dimension(:), allocatable :: NodalDensity 
        !real(REAL_TYPE), dimension(:), allocatable :: ActiveNodeElementVolume 
        !logical, dimension(:), allocatable :: ActiveNode
        !logical, dimension(:), allocatable :: HydraulicHeadLoadedElemID
        !
        !! Variables which are updated within calculation process
        !real(REAL_TYPE), dimension(:, :), allocatable :: ExtLoad ! Nodal load array of external load
        !! Nodal load array of external load (Total traction applied on nodes)
        !real(REAL_TYPE), dimension(:, :, :), allocatable :: ExtLoadTotal
        !real(REAL_TYPE), dimension(:, :), allocatable :: HydraulicHeadLoadTotal
        !real(REAL_TYPE), dimension(:, :), allocatable :: IntLoad ! Nodal load array of internal load
        !real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadPrevious ! Nodal load array of internal load at previous time step
        !real(REAL_TYPE), dimension(:, :), allocatable :: GravityLoad ! Nodal load array of gravity load
        !real(REAL_TYPE), dimension(:, :), allocatable :: BulkViscLoad
        !real(REAL_TYPE), dimension(:, :), allocatable :: LumpedMassDry ! Lumped mass vector of dry soil
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: TotalVelocitySoil ! Nodal total velocity
        !real(REAL_TYPE), dimension(:,:), &
        !  allocatable :: TotalVelocitySoilPrevious ! Nodal total velocity of the previous time step
        !real(REAL_TYPE), dimension(:), &
        !  allocatable :: TotalDisplacementSoil ! Nodal total displacement
        !real(REAL_TYPE), dimension(:), &
        !  allocatable :: PhaseDisplacementSoil ! Nodal total displacement of the load phase
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: IncrementalDisplacementSoil ! Nodal incremental displacement
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: AccumulatedIncDisplacementSoil ! Nodal incremental displacement
        !real(REAL_TYPE), dimension(:),  &
        !  allocatable :: AccumulatedDisplacementSoil ! Accumulated displacements of a load phase
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: AccelerationSoil  ! Nodal acceleration of soil
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: RateofMomentum  ! Rate of momentum
        !integer(INTEGER_TYPE), dimension(:), allocatable :: ActiveElement ! ID of the active elements
        !real(REAL_TYPE), dimension(:), &
        !  allocatable :: TotalDisplacementWater ! Nodal total displacement of water (FEM only)
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: TotalVelocityWater ! Nodal total velocity for water
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: TotalVelocityWaterPrevious ! Nodal total velocity for water for the previous time step
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: IncrementalDisplacementWater ! Nodal incremental displacement
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: IncrementalDisplacementWaterPrevious ! Nodal incremental displacement at previous time step
        !real(REAL_TYPE), dimension(:),  &
        !  allocatable :: AccumulatedDisplacementWater ! Accumulated displacements of a load phase
        !real(REAL_TYPE), dimension(:), &
        !  allocatable :: PhaseDisplacementWater ! Nodal total displacement of the load phase
        !real(REAL_TYPE), dimension(:), &
        !  allocatable :: TotalDisplacementGas ! Nodal total displacement of gas (FEM only)
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: TotalVelocityGas ! Nodal total velocity for gas
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: NonAdvectiveFluxAirInWater ! Nodal Non advective flux of Air in the Fluid 
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: NonAdvectiveFluxVapourInGas ! Nodal Non advective flux of Vapour in the Gas 
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: AdvectiveFluxDarcyWater ! Nodal advective flux of water 
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: AdvectiveFluxDarcyAir ! Nodal advective flux of air  
        !real(REAL_TYPE), dimension(:, :), &
        !allocatable :: ThermalConductionFlux ! Nodal Heat Conduction flux
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: FReaction ! Nodal reactions on a user-defined mesh boundary
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: FReactionWater
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: FReactionGas
        !real(REAL_TYPE), dimension(:), allocatable :: RateVolStrain
        !real(REAL_TYPE), dimension(:, :), &
        !  allocatable :: NodalUnitMassGradient
        !real(REAL_TYPE), dimension(:,:), &
        !    allocatable :: SoilSurfaceNodeCoordMatrix
        !real(REAL_TYPE), dimension(:,:), &
        !    allocatable :: PhreaticSurfaceNodeCoordMatrix
        !real(REAL_TYPE), dimension(:,:), &
        !    allocatable :: HydrHeadArea
        !real(REAL_TYPE), dimension(:,:), &
        !    allocatable :: SeepageArea
        !real(REAL_TYPE), dimension(:,:), &
        !    allocatable :: InfiltrationArea
        !real(REAL_TYPE), dimension(:), allocatable :: InfiltrationRate
        !
        !
        !> stores material points indeces
        !!!!integer(INTEGER_TYPE), dimension(:,:, :), allocatable :: GetParticleIndex ! (rank 3 for multipatch implementation) ! I COMMENTED THIS 5/26/2024
        !
        !!> stores determinant of elements
        !real(REAL_TYPE), dimension(:), allocatable :: ElementDeterminant
        !
        !! quasi static additions
        !real(REAL_TYPE), dimension(:), allocatable :: ExtFlow ! mass balance steady flow
        !real(REAL_TYPE), dimension(:), allocatable :: IntFlow ! correction term
        !real(REAL_TYPE), dimension(:), allocatable :: RateOfFlux ! flow residual
        !real(REAL_TYPE), dimension(:), allocatable :: IncrementalPressure ! incremental water pressure
        !real(REAL_TYPE), dimension(:), allocatable :: TotalPressure ! total water pressure
        !real(REAL_TYPE), dimension(:), allocatable :: SubIncrementalDisplacement
        !real(REAL_TYPE), dimension(:), allocatable :: SubIncrementalPressure
        !
        !integer(INTEGER_TYPE), dimension(:), allocatable :: ConsideredElemReaction !vector containing the ElementID of those elements considered in the integration of reaction forces
        !logical, dimension(:,:), allocatable :: IsReactionNodeSurface
        !character (len=255), dimension(:), allocatable :: OutputSurfaceName
        !logical, dimension(:), allocatable :: IsReactionNode!Determine if a node is on the reaction surface
        !
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPGlobalPositionElement		
        !
        !
        
        !! NURBS variables 
        !real(REAL_TYPE), allocatable, dimension(:) :: XiKnotEntriesUNIQUE
        !integer(INTEGER_TYPE) :: NXiKnotEntriesUNIQUE
        !
        !real(REAL_TYPE), allocatable, dimension(:) :: EtaKnotEntriesUNIQUE
        !integer(INTEGER_TYPE) :: NEtaKnotEntriesUNIQUE
        !
        !
        !real(REAL_TYPE), allocatable, dimension(:,:) :: KnotBezierMesh 
        !integer(INTEGER_TYPE), allocatable, dimension(:, :) :: KnotBezierMeshConnectivity ! number of elements x number of nodes in element = 4 x 4

        
        
		
      contains ! Routines of this module

        
        logical function isLinearElastic(MaterialID) result(res)
        implicit none
        integer(INTEGER_TYPE), intent(in) :: MaterialID
        res = GetConstitutiveModel(MaterialID) == ESM_LINEAR_ELASTICITY
        end function isLinearElastic

        logical function isModifiedCamClay(MaterialID) result(res)
        implicit none
        integer(INTEGER_TYPE), intent(in) :: MaterialID
        res = GetConstitutiveModel(MaterialID) == ESM_MODIFIED_CAM_CLAY
        end function isModifiedCamClay

        logical function isStrainSofteningMohrCoulomb(MaterialID) result(res)
        implicit none
        integer(INTEGER_TYPE), intent(in) :: MaterialID
        res = GetConstitutiveModel(MaterialID) == ESM_MOHR_COULOMB_STRAIN_SOFTENING
        end function isStrainSofteningMohrCoulomb

        logical function isMohrCoulombStandard(MaterialID) result(res)
        implicit none
        integer(INTEGER_TYPE), intent(in) :: MaterialID
        res = GetConstitutiveModel(MaterialID) == ESM_MOHR_COULOMB_STANDARD
        end function isMohrCoulombStandard

        character(len=32) function GetConstitutiveModel(MaterialID) result(SoilModel)
        implicit none
        integer(INTEGER_TYPE), intent(in) :: MaterialID
        ! name of constitutive model as specified in GOM-file
        SoilModel = MatParams(MaterialID)%MaterialModel
        end function GetConstitutiveModel

        
        subroutine InitialiseMeshData()
        !**********************************************************************
        !
        !    Function:  Contains code for initialising mesh data by allocating arrays
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          character(len = 255) :: FileName, TName, Bname
          
          integer(INTEGER_TYPE) :: BoundCounter_X = 0
          integer(INTEGER_TYPE) :: BoundCounter_Y = 0
          integer(INTEGER_TYPE) :: BoundCounter_Z = 0

                    
          
          integer(INTEGER_TYPE) :: I, J, L, IDElement, IDConnectivity, IFxN, &
                     NElements, NNodes, NFixNod, FixNodID, NRemFixNod, FixRemNodID, &
                     NLoadedElementSidesGasNodes, & !NLoadedElementSidesSolidNodes, NLoadedElementSidesWaterNodes, 
                     NLoadedElementSidesGasmatPoints, & !NLoadedElementSidesSolidMatPoints, NLoadedElementSidesWaterMatPoints, 
                     SoilSurfaceNumberofSides, PhreaticSurfaceNumberofSides, &
                     NLoadedElementSidesSolidNodesB, NLoadedElementSidesWaterNodesB, NLoadedElementSidesGasNodesB, &
                     NLoadedElementSidesSolidMatPointsB, NLoadedElementSidesWaterMatPointsB, NLoadedElementSidesGasmatPointsB 
          
          integer(INTEGER_TYPE), dimension(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches) :: NLoadedElementSidesSolidNodes
          
          integer(INTEGER_TYPE), dimension(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches) :: NLoadedElementSidesWaterNodes
          
          integer(INTEGER_TYPE), dimension(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches) :: NLoadedElementSidesSolidMatPoints
          
          integer(INTEGER_TYPE), dimension(NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches) :: NLoadedElementSidesWaterMatPoints
          
          integer(INTEGER_TYPE) :: IDNControlGrid
          
          
          integer(INTEGER_TYPE) :: IError ! used for error control
          integer(INTEGER_TYPE), dimension(NVECTOR) :: NodFixBC, NodRemBC
          real(REAL_TYPE) :: DumR(4)
          
          
          ! Multipatch variables 
          !integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: IPatch_Other ! comparative integer index
          real(REAL_TYPE), dimension(NVECTOR) :: NodeCoordinatesSelected
          integer(INTEGER_TYPE) :: INode
          integer(INTEGER_TYPE) :: INode_Other
          integer(INTEGER_TYPE) :: GlobalMultipatchNodalIndex
          integer(INTEGER_TYPE) :: ControlPointIndex
          integer(INTEGER_TYPE) :: GlobalNodeID
          integer(INTEGER_TYPE) :: NoOfItems
          logical :: logic_test
          logical :: repetition_test
          integer(INTEGER_TYPE) :: ii
          
          !call DestroyMeshData() ! --> I commented this
          
          ! Multipatch variable
          integer(INTEGER_TYPE), dimension(NURBS%Maximum_NControlPoints, Counters%NPatches, Counters%NPatches) :: NumberOfRepetitions
          logical :: NonUnique
          
          logical :: IsItUnique
          logical :: IsItWritten
          integer(INTEGER_TYPE) :: IPatch_Loop_1
          integer(INTEGER_TYPE) :: IPatch_Loop_2
          integer(INTEGER_TYPE) :: INode_Loop_1
          integer(INTEGER_TYPE) :: INode_Loop_2
          
          
          ! TRACTION VARIABLES
          integer(INTEGER_TYPE) :: IControlPoint_Traction
          integer(INTEGER_TYPE) :: IControlPoint_ActualMaterial
          
          ! TRACTION ON MPs VARIABLES
          integer(INTEGER_TYPE) :: IElement_Traction
          integer(INTEGER_TYPE) :: IElement_ActualMaterial
          integer(INTEGER_TYPE) :: IElementCon_Traction
          integer(INTEGER_TYPE) :: IElementCon_ActualMaterial
          integer(INTEGER_TYPE) :: ICount
          

          real(REAL_TYPE), allocatable, dimension(:,:) :: NodalCoordinatesTractionSide
          real(REAL_TYPE), allocatable, dimension(:,:) :: NodalCoordinatesGlobalSide
          
          !real(REAL_TYPE), dimension(:,:) :: NodalCoordinatesTractionSide

          integer(INTEGER_TYPE), allocatable, dimension(:, :) :: VerticesID !(2, 4)

          integer(INTEGER_TYPE) :: GlobalElementForTraction
          integer(INTEGER_TYPE) :: ISide
          
          logical :: element_success
          logical, allocatable, dimension(:,:) :: side_test
          
          
          logical:: ParticleStatus (:)
          integer(INTEGER_TYPE) :: NParticles
          allocatable :: ParticleStatus 
              
          
          integer(INTEGER_TYPE) :: NumberOfElements_Xi
          integer(INTEGER_TYPE) :: NumberOfElements_Eta
              
          
          
          ! Volumetric locking variables 
          integer(INTEGER_TYPE) :: OrderDifference
          integer(INTEGER_TYPE) :: I_LowerOrder
              
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1 

          ! Load system generalization 
          integer(INTEGER_TYPE) :: ILoadSystem_Temporary = 1
          integer(INTEGER_TYPE) :: ILoadSystem
              
          !External :: FExist    

          !integer(INTEGER_TYPE) :: IPatch_Temporary
              
          IPatch_Temporary = 1

          !allocate( ElementUponWhichLoadIsApplied_NURBS(NLoadedElementSidesSolidNodes), stat = IError)      
          allocate( ElementUponWhichLoadIsApplied_NURBS(NURBS%Maximum_nel_NURBS_Traction, NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches), stat = IError)      
          ElementUponWhichLoadIsApplied_NURBS = 0.0
          
          
          
          
          NumberOfRepetitions = 0
          repetition_test = .false.
          IsItUnique = .true.
          !side_test = .true.
          !IPatch_Loop_1 = 
          
          
          ! initialization
          NLoadedElementSidesSolidNodes = 0
          NLoadedElementSidesWaterNodes = 0
          NLoadedElementSidesGasNodes = 0
          NLoadedElementSidesSolidMatPoints = 0
          NLoadedElementSidesWaterMatPoints = 0
          NLoadedElementSidesGasMatPoints = 0
          SoilSurfaceNumberofSides = 0
          PhreaticSurfaceNumberofSides = 0
          
          NLoadedElementSidesSolidNodesB = 0
          NLoadedElementSidesWaterNodesB = 0
          NLoadedElementSidesGasNodesB = 0
          NLoadedElementSidesSolidMatPointsB = 0
          NLoadedElementSidesWaterMatPointsB = 0
          NLoadedElementSidesGasMatPointsB = 0
          
          !ControlNetGrid = 0
          
          ! open GOM file
          FileName=Trim(CalParams%FileNames%ProjectName)//'.GOM'
          if (FExist(FileName)) open(GOMunit, FILE=FileName)

          do
            read(GOMunit,'(A)') TName !reading from GOM file 
            BName = TName
            
            
            if (trim(BName) == '$$PATCH_ID') then 
                ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                IPatch = -1
                ! read patch ID from the GOM file
                read(GOMunit, *) IPatch ! --> now we knot patch index so we can read in the following: 
                                        ! 1) Number of control points for each patch
                                        ! 1) Number of knot entries (xi, eta, zeta) 
                                        ! 2) Order of knot (xi, eta, zeta)
              
                
                
            else if (trim(BName) == '$$NURBS_CONTROL_POINTS') then
              do I = 1, Counters%NodTot(IPatch)!NNodes ! loop over all nodes
                  
                  !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
                  
                read (GOMunit, *) (NodalCoordinates(I, J, IPatch), J = 1, NVECTOR) ! nodal coordinates
                if ( ISAXISYMMETRIC ) call Assert(NodalCoordinates(I, 1, IPatch) >= 0.0, 'GOM file: For axysimmetric analysis the x-coordinate (radial direction) has to be zero or positive.' )  
				if ( IS3DCYLINDRIC ) call Assert(NodalCoordinates(I, 1, IPatch) >= 0.0, 'GOM file: For 3D Cylindrical analysis the x-coordinate (radial direction) has to be zero or positive.' )
              end do ! loop over all nodes
              NodalCoordinatesUpd = NodalCoordinates
              NodalOriginalCoord = NodalCoordinates
              
              
            else if (trim(BName) == '$$NURBS_CONTROL_POINT_WEIGHTS') then
                do I = 1, Counters%NodTot(IPatch)!NNodes ! loop over all nodes
                    
                    !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
                    
                    read (GOMunit, *) ControlPoint_Weights(I,IPatch) ! for every control point there must be a weight
                end do 
                
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Volumetric locking
                
            else if (trim(BName) == '$$NURBS_CONTROL_POINTS_VOLLOCK_SMOOTH') then
        
                do I = 1, NURBS%NControlPoints_VolLockSmooth(IPatch)!NNodes ! loop over all nodes
                  !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
                    read (GOMunit, *) (NodalCoordinates_VolLockSmooth(I, J, IPatch), J = 1, NVECTOR) ! nodal coordinates
                    if ( ISAXISYMMETRIC ) call Assert(NodalCoordinates_VolLockSmooth(I, 1, IPatch) >= 0.0, 'GOM file: For axysimmetric analysis the x-coordinate (radial direction) has to be zero or positive.' )  
				    if ( IS3DCYLINDRIC ) call Assert(NodalCoordinates_VolLockSmooth(I, 1, IPatch) >= 0.0, 'GOM file: For 3D Cylindrical analysis the x-coordinate (radial direction) has to be zero or positive.' )
                end do ! loop over all nodes
              
              
            else if (trim(BName) == '$$NURBS_CONTROL_POINT_WEIGHTS_VOLLOCK_SMOOTH') then
                do I = 1, NURBS%NControlPoints_VolLockSmooth(IPatch)!NNodes ! loop over all nodes
                    read (GOMunit, *) ControlPoint_Weights_VolLockSmooth(I,IPatch) ! for every control point there must be a weight
                end do 
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                
                
                
                
                
                
            else if (trim(BName) == '$$TRACTION_LOAD_SYSTEM') then 
                ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                ILoadSystem = -1
                ! read patch ID from the GOM file
                read(GOMunit, *) ILoadSystem ! --> now we knot patch index so we can read in the following: 
                                             ! 1) Number of control points for each patch
                                             ! 1) Number of knot entries (xi, eta, zeta) 
                                             ! 2) Order of knot (xi, eta, zeta)
                
                
                
                
              
                
            else if (trim(BName) == '$$NURBS_CONTROL_POINTS_TRACTION') then
                
                !element_success = .false.
                if (NURBS%NControlPoints_Traction(ILoadSystem, IPatch)==0) cycle
                
              do I = 1, NURBS%NControlPoints_Traction(ILoadSystem, IPatch) !(IPatch)!NNodes ! loop over all nodes
                  
                  !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
                  
                read (GOMunit, *) (NodalCoordinates_Traction(I, J, ILoadSystem, IPatch), J = 1, NVECTOR) ! nodal coordinates !, IPatch
                if ( ISAXISYMMETRIC ) call Assert(NodalCoordinates_Traction(I, 1, ILoadSystem, IPatch) >= 0.0, 'GOM file: For axysimmetric analysis the x-coordinate (radial direction) has to be zero or positive.' )  !, IPatch
				if ( IS3DCYLINDRIC ) call Assert(NodalCoordinates_Traction(I, 1, ILoadSystem, IPatch) >= 0.0, 'GOM file: For 3D Cylindrical analysis the x-coordinate (radial direction) has to be zero or positive.' ) !, IPatch
              
                !element_success = .true.
              end do ! loop over all nodes
              !NodalCoordinatesUpd = NodalCoordinates
              !NodalOriginalCoord = NodalCoordinates
              

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              ! this basically gives us the global number given that we input the number of the local control point for that surface
              !allocate( Traction_Connecting_Local_To_Global_ControlPoints(NControlPoints_Traction(ILoadSystem_Temporary, IPatch_Temporary)),  stat = IError)
              !Traction_Connecting_Local_To_Global_ControlPoints = 0
              allocate( Traction_Connecting_Local_To_Global_ControlPoints(NURBS%Maximum_NControlPoints_Traction, NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches),  stat = IError)
              Traction_Connecting_Local_To_Global_ControlPoints = 0
              
              !allocate( Traction_Connecting_Local_To_Global_Elements(nel_NURBS_Traction(ILoadSystem_Temporary, IPatch_Temporary)),  stat = IError)
              !Traction_Connecting_Local_To_Global_Elements = 0
              allocate( Traction_Connecting_Local_To_Global_Elements(NURBS%Maximum_nel_NURBS_Traction, NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches),  stat = IError)
              Traction_Connecting_Local_To_Global_Elements = 0
              
              !allocate( Traction_Connecting_Local_To_Global_Sides(nel_NURBS_Traction(ILoadSystem_Temporary, IPatch_Temporary)),  stat = IError)
              !Traction_Connecting_Local_To_Global_Sides = 0
              allocate( Traction_Connecting_Local_To_Global_Sides(NURBS%Maximum_nel_NURBS_Traction, NURBS%Maximum_NLoadSystems_Traction, Counters%NPatches),  stat = IError)
              Traction_Connecting_Local_To_Global_Sides = 0
              
              !if (element_success) then 
              !    element_success = .false.  
              !! THIS IS FOR APPLYING TRACTION ON THE NODES 
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !do IPatch = 1, Counters%NPatches
              !    do ILoadSystem = 1, NLoadSystems_Traction(IPatch)
              !
              !      
              !        call MatchControlPointsNodes(ILoadSystem, IPatch) ! traction on the control points 
              !
              !
              !        !! THIS IS FOR APPLYING TRACTION ON THE MATERIAL POINTS 
              !        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !        call MatchControlPointsElements(ILoadSystem, IPatch) ! traction on the elements
              !        call MatchTractionElementToActualMaterialSide(ILoadSystem, IPatch) ! compare traction side to all sides of actual material. 
              !                                                ! We need to know side our traction is on. 
              !!end if
              !
              !    end do 
              !end do 
              
              
            !else if (trim(BName) == '$$NURBS_CONTROL_POINT_WEIGHTS_TRACTION') then
            !    do I = 1, NControlPoints_Traction(ILoadSystem, IPatch) !(IPatch)!NNodes ! loop over all nodes
            !        
            !        !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
            !        
            !        read (GOMunit, *) ControlPoint_Weights_Traction(I, ILoadSystem, IPatch) ! for every control point there must be a weight
            !    end do 
                
                
                

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                
                
                
                
                ! APPLY TRACTION ON THE NODES
            elseif (trim(BName) == '$$TRACTION_APPLY_LOCAL_ELEMENTS') then 
                
                !read(GOMunit,*) NLoadedElementSidesSolidNodes(ILoadSystem, IPatch) ! Loaded faces
                
                NLoadedElementSidesSolidNodes = NURBS%NLoadSystems_Traction_NumberOfElements(ILoadSystem, IPatch)
                
                if (NLoadedElementSidesSolidNodes(ILoadSystem, IPatch) > 0) then 
                    
                    do L = 1, NLoadedElementSidesSolidNodes(ILoadSystem, IPatch)
                        
                        read(GOMunit,*) ElementUponWhichLoadIsApplied_NURBS(L, ILoadSystem, IPatch)
                        
                    end do 
                    
                end if
                
                
                
            elseif (trim(BName) == '$$TRACTION_APPLY_LOCAL_ELEMENTS_WATER') then 
                
                !read(GOMunit,*) NLoadedElementSidesWaterNodes ! Loaded faces
                
                NLoadedElementSidesWaterNodes = NURBS%NLoadSystems_Traction_NumberOfElements(ILoadSystem, IPatch)
                
                !allocate( ElementUponWhichLoadIsApplied_NURBS(NLoadedElementSidesWaterNodes), stat = IError)
                !ElementUponWhichLoadIsApplied_NURBS = 0.0
                
                if (NLoadedElementSidesWaterNodes(ILoadSystem, IPatch) > 0) then 
                    
                    do L = 1, NLoadedElementSidesWaterNodes(ILoadSystem, IPatch)
                        
                        read(GOMunit,*) ElementUponWhichLoadIsApplied_NURBS(L, ILoadSystem, IPatch)
                        
                    end do 
                    
                end if
                    
                
                
                ! APPLY TRACTION ON THE MATERIAL POINTS 
            elseif (trim(BName) == '$$TRACTION_APPLY_LOCAL_ELEMENTS_ON_MPs') then 
                
                !read(GOMunit,*) NLoadedElementSidesSolidMatPoints() ! Loaded faces
                
                NLoadedElementSidesSolidMatPoints = NURBS%NLoadSystems_Traction_NumberOfElements
                
                !allocate( ElementUponWhichLoadIsApplied_NURBS(NLoadedElementSidesSolidMatPoints), stat = IError)
                !ElementUponWhichLoadIsApplied_NURBS = 0.0
                
                if (NLoadedElementSidesSolidMatPoints(ILoadSystem, IPatch) > 0) then 
                    
                    do L = 1, NLoadedElementSidesSolidMatPoints(ILoadSystem, IPatch)
                        
                        read(GOMunit,*) ElementUponWhichLoadIsApplied_NURBS(L, ILoadSystem, IPatch)
                        
                    end do 
                    
                end if
            
                
                
                
                
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            else if (trim(BName) == '$$XI_KNOT') then
                do I = 1, NURBS%NXiKnotEntries(IPatch) ! loop over all nodes
                    read (GOMunit, *) NURBS%XiKnotEntries(I,IPatch) ! Xi knot vector     
                end do ! loop over all nodes
                
                call FindNumberOfUniqueEnteriesInVector(NURBS%XiKnotEntries(:,IPatch), NURBS%NXiKnotEntries(IPatch), NURBS%NXiKnotEntriesUNIQUE(IPatch))
                
                ! find the maximum of the unique entries - note that this will update with every patch
                NURBS%Maximum_NXiKnotEntriesUNIQUE = maxval(NURBS%NXiKnotEntriesUNIQUE)
                
                if (CalParams%ApplyFBar) then 
                    NURBS%XiKnotEntries_VolLockSmooth = 0.0
                    
                    OrderDifference = NURBS%NXiKnotOrder(IPatch) - NURBS%NXiKnotOrder_VolLockSmooth(IPatch)
                    
                    do I_LowerOrder = NURBS%NXiKnotOrder_VolLockSmooth(IPatch)+2, NURBS%NXiKnotEntries_VolLockSmooth(IPatch)
                    
                        NURBS%XiKnotEntries_VolLockSmooth(I_LowerOrder, IPatch) = NURBS%XiKnotEntries(I_LowerOrder+OrderDifference, IPatch) 
                        
                    end do 
                    
                end if 
                
                
            else if (trim(BName) == '$$ETA_KNOT') then
                do I = 1, NURBS%NEtaKnotEntries(IPatch) ! loop over all nodes
                    read (GOMunit, *) NURBS%EtaKnotEntries(I,IPatch) ! Eta knot vector 
                end do ! loop over all nodes
                
                ! find unique knots items in the knot vector                
                !call FindUniqueEnteriesInVector(EtaKnotEntries(:,IPatch), NEtaKnotEntries(IPatch), EtaKnotEntriesUNIQUE, NEtaKnotEntriesUNIQUE)
                call FindNumberOfUniqueEnteriesInVector(NURBS%EtaKnotEntries(:,IPatch), NURBS%NEtaKnotEntries(IPatch), NURBS%NEtaKnotEntriesUNIQUE(IPatch))
                
                ! find the maximum of the unique entries - note that this will update with every patch
                NURBS%Maximum_NEtaKnotEntriesUNIQUE = maxval(NURBS%NEtaKnotEntriesUNIQUE)
                
                
                if (CalParams%ApplyFBar) then 
                    NURBS%EtaKnotEntries_VolLockSmooth = 0.0
                    
                    OrderDifference = NURBS%NEtaKnotOrder(IPatch) - NURBS%NEtaKnotOrder_VolLockSmooth(IPatch)
                    
                    do I_LowerOrder = NURBS%NEtaKnotOrder_VolLockSmooth(IPatch)+2, NURBS%NEtaKnotEntries_VolLockSmooth(IPatch)
                    
                        NURBS%EtaKnotEntries_VolLockSmooth(I_LowerOrder, IPatch) = NURBS%EtaKnotEntries(I_LowerOrder+OrderDifference, IPatch) 
                        
                    end do 
                    
                end if
                
            
            else if (trim(BName) == '$$ZETA_KNOT') then
                do I = 1, NURBS%NZetaKnotEntries(IPatch) ! loop over all nodes
                    read (GOMunit, *) NURBS%ZetaKnotEntries(I,IPatch) ! Eta knot vector
                end do ! loop over all nodes

                call FindNumberOfUniqueEnteriesInVector(NURBS%ZetaKnotEntries(:,IPatch), NURBS%NZetaKnotEntries(IPatch), NURBS%NZetaKnotEntriesUNIQUE(IPatch))

                ! find the maximum of the unique entries - note that this will update with every patch
                NURBS%Maximum_NZetaKnotEntriesUNIQUE = maxval(NURBS%NZetaKnotEntriesUNIQUE)
                
                if (CalParams%ApplyFBar) then 
                    NURBS%ZetaKnotEntries_VolLockSmooth = 0.0
                    
                    OrderDifference = NURBS%NZetaKnotOrder(IPatch) - NURBS%NZetaKnotOrder_VolLockSmooth(IPatch)
                    
                    do I_LowerOrder = NURBS%NZetaKnotOrder_VolLockSmooth(IPatch)+2, NURBS%NZetaKnotEntries_VolLockSmooth(IPatch)
                    
                        NURBS%ZetaKnotEntries_VolLockSmooth(I_LowerOrder, IPatch) = NURBS%ZetaKnotEntries(I_LowerOrder+OrderDifference, IPatch) 
                        
                    end do 
                    
                end if
                
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                !else if (trim(BName) == '$$XI_KNOT_TRACTION') then
                !do I = 1, NXiKnotEntries_Traction ! loop over all nodes
                !            !(IPatch)
                !    read (GOMunit, *) XiKnotEntries_Traction(I) ! Xi knot vector
                !                        !(I,IPatch)
                !    
                !end do ! loop over all nodes
                
                
            !else if (trim(BName) == '$$ETA_KNOT_TRACTION') then
            !    do I = 1, NEtaKnotEntries_Traction ! loop over all nodes
            !        read (GOMunit, *) EtaKnotEntries_Traction(I) ! Eta knot vector 
            !                            !,IPatch
            !
            !    end do ! loop over all nodes
              
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                
                
                
                
            ! note that for NURBS, we need to have both the IEN array (which really tells you the connectivity in the parametric space) 
            ! and the control point coordinates      
            else if (trim(BName) == '$$STARTELEMCON') then
              if (ELEMENTTYPE == TETRAOLD) then ! elementtype before v2018.2  
                
                  !do IPatch = 1, Counters%NPatches
                  do IDElement = 1, Counters%NEl(IPatch) !NURBS%nel_NURBS(IPatch)!NElements ! loop over elements
                      ! in the GOM-file always 10 nodes per element are written
                      read (GOMunit,*) (ElementConnectivities10Node(IDConnectivity, IDElement), IDConnectivity = 1, N_NODES_HOE) 
                  do J = 1, ELEMENTNODES ! loop over element nodes
                    ElementConnectivities(J, IDElement, IPatch) = ElementConnectivities10Node(J, IDElement)
                  end do ! loop over element nodes
                  end do  ! loop over elements
                  !end do ! loop over patches
              
                
              else ! elementtypes since v2018.2
                
                  !do IPatch = 1, Counters%NPatches ! loop over patches
                  do IDElement = 1, Counters%NEl(IPatch)!NURBS%nel_NURBS(IPatch)!NElements ! loop over elements
                  read (GOMunit,*) (ElementConnectivities(IDConnectivity, IDElement, IPatch), IDConnectivity = 1, ELEMENTNODES)
                  end do ! loop over elements 
                  !end do ! loop over patches
              end if  

            
              
              
              
              
              
              
            else if (trim(BName) == '$$CONTROL_NET_CONNECTIVITY') then
                
                ControlNetGrid = 4
                
                read(GOMunit, *) NumberOfControlNetGrid
                
                allocate(ControlNetConnectivity(ControlNetGrid, NumberOfControlNetGrid), stat = IError) !size equal to no of elements, no of dimensions ! -> NURBS implementation
                
                do IDNControlGrid = 1, NumberOfControlNetGrid ! loop over control net grid
                    read (GOMunit,*) (ControlNetConnectivity(IDConnectivity, IDNControlGrid), IDConnectivity = 1, ControlNetGrid) ! harcoded control net tessalation
                end do
                
                
            elseif (trim(BName) == '$$START_FIXITY_SURFACE_SOLID') then
              read(GOMunit,*) NFixNod		
              do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDisp(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if  
				end do
              end do ! loop over fixed nodes  

            elseif (trim(BName) == '$$START_FIXITY_LINE_SOLID') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDisp(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if
                end do  
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_FIXITY_POINT_SOLID') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDisp(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if
                end do
              end do ! loop over fixed nodes 
              
              !----------------------------------------------------------------
              
    !          elseif (trim(BName) == '$$START_FIXITY_SURFACE_SOLID_VOLLOCK_SMOOTH') then
    !          read(GOMunit,*) NFixNod		
    !          do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
    !            read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
    !            do I = 1, NVECTOR ! loop over vector size
    !              if (NodFixBC(I) == 1) then 
    !                NodalPrescibedDisp_VolLockSmooth(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
    !              end if  
				!end do
    !          end do ! loop over fixed nodes  
    !
    !        elseif (trim(BName) == '$$START_FIXITY_LINE_SOLID_VOLLOCK_SMOOTH') then
    !          read(GOMunit,*) NFixNod
    !          do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
    !            read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
    !            do I = 1, NVECTOR ! loop over vector size
    !              if (NodFixBC(I) == 1) then 
    !                NodalPrescibedDisp_VolLockSmooth(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
    !              end if
    !            end do  
    !          end do ! loop over fixed nodes 
    !
    !        elseif (trim(BName) == '$$START_FIXITY_POINT_SOLID_VOLLOCK_SMOOTH') then
    !          read(GOMunit,*) NFixNod
    !          do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
    !            read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
    !            do I = 1, NVECTOR ! loop over vector size
    !              if (NodFixBC(I) == 1) then 
    !                NodalPrescibedDisp_VolLockSmooth(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
    !              end if
    !            end do
    !          end do ! loop over fixed nodes 
              
              
              !----------------------------------------------------------------

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_SURFACE_SOLID') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (solid)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveSolidFixities) then
                  do I = 1, NVECTOR ! loop over vector size
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDisp(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do
                end if
              end do ! loop over removed fixed nodes

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_LINE_SOLID') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (solid)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveSolidFixities) then
                  do I = 1, NVECTOR ! loop over vector size
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDisp(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if
              end do ! loop over removed fixed nodes

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_POINT_SOLID') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (solid)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveSolidFixities) then
                  do I = 1, NVECTOR ! loop over degrees of freedom
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDisp(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if
              end do ! loop over removed fixed nodes

            elseif (trim(BName) == '$$START_FIXITY_SURFACE_LIQUID') then 
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (liquid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispWater(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if
                end do
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_FIXITY_LINE_LIQUID') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (liquid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over degrees of freedom
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispWater(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if
                end do
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_FIXITY_POINT_LIQUID') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (liquid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispWater(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if
                end do
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_SURFACE_LIQUID') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (liquid)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveLiquidFixities) then
                  do I = 1, NVECTOR ! loop over vector size
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDispWater(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if  
              end do ! loop over removed fixed nodes

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_LINE_LIQUID') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (liquid)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveLiquidFixities) then
                  do I = 1, NVECTOR ! loop over vector size
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDispWater(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if
              end do ! loop over removed fixed nodes

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_POINT_LIQUID') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (liquid)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveLiquidFixities) then
                  do I = 1, NVECTOR ! loop over vector size
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDispWater(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if
              end do ! loop over removed fixed nodes

            elseif (trim(BName) == '$$START_FIXITY_SURFACE_GAS') then 
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (gas)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispGas(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if
                end do  
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_FIXITY_LINE_GAS') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (gas)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over degrees of freedom
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispGas(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if
                end do  
              end do ! loop over fixed nodes

            elseif (trim(BName) == '$$START_FIXITY_POINT_GAS') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (gas)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispGas(FixNodID,I,IPatch) = 0d0 ! zero prescribed displacement
                  end if
                end do  
              end do ! loop over fixed nodes

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_SURFACE_GAS') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (gas)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveGasFixities) then
                  do I = 1, NVECTOR ! loop over degrees of freedom
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDispGas(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if
              end do ! loop over removed fixed nodes

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_LINE_GAS') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (gas)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveGasFixities) then
                  do I = 1, NVECTOR ! loop over vector size
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDispGas(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if
              end do ! loop over removed fixed nodes

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_POINT_GAS') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (gas)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveGasFixities) then
                  do I = 1, NVECTOR ! loop over vector size
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDispGas(FixRemNodID,I,IPatch) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if
              end do ! loop over removed fixed nodes

            !else if (trim(BName)=='$$START_LOAD_ON_NODES_SOLID') then
            !  read(GOMunit,*) NLoadedElementSidesSolidNodes ! Loaded faces
            !  if (NLoadedElementSidesSolidNodes > 0) then
            !    if ( NDIM == 3 ) then
            !        
            !        ! not multipatch generalized and goes for the first patch only
            !        !---------------------------------------------
            !        ! hardcoded for xi direction traction and needs generalization
            !        allocate(LoadOnNodesConnectivitiesSolid(NXiKnotOrder(IPatch_Temporary)+1, NLoadedElementSidesSolidNodes), stat = IError) 
            !        allocate(LoadValuesOnNodesSolid_NURBS(NDOFL, NLoadedElementSidesSolidNodes), stat = IError)
            !        ! I am considering this as an input from the GOM file but this needs to be generalized to a 
            !        ! more generic search which element belongs to these control points
            !        allocate(ElementUponWhichLoadIsApplied_NURBS(NLoadedElementSidesSolidNodes) )
            !        
            !        
            !        
            !        do L = 1, NLoadedElementSidesSolidNodes ! loop over loaded faces
            !        read(GOMunit,*) LoadOnNodesConnectivitiesSolid(1,L), LoadOnNodesConnectivitiesSolid(2,L),  &
            !                        LoadOnNodesConnectivitiesSolid(3,L), LoadOnNodesConnectivitiesSolid(4,L), &
            !                        ElementUponWhichLoadIsApplied_NURBS(L), &
            !                        LoadValuesOnNodesSolid_NURBS(1, L), &
            !                        LoadValuesOnNodesSolid_NURBS(2, L), &
            !                        LoadValuesOnNodesSolid_NURBS(3, L)
            !                        !LoadValuesOnNodesSolid(1,1,L), LoadValuesOnNodesSolid(1,2,L), LoadValuesOnNodesSolid(1,3,L), &
            !                        !LoadValuesOnNodesSolid(3,1,L), LoadValuesOnNodesSolid(3,2,L), LoadValuesOnNodesSolid(3,3,L), &
            !                        !LoadValuesOnNodesSolid(5,1,L), LoadValuesOnNodesSolid(5,2,L), LoadValuesOnNodesSolid(5,3,L), &
            !                        !LoadValuesOnNodesSolid(2,1,L), LoadValuesOnNodesSolid(2,2,L), LoadValuesOnNodesSolid(2,3,L), &
            !                        !LoadValuesOnNodesSolid(4,1,L), LoadValuesOnNodesSolid(4,2,L), LoadValuesOnNodesSolid(4,3,L), &
            !                        !LoadValuesOnNodesSolid(6,1,L), LoadValuesOnNodesSolid(6,2,L), LoadValuesOnNodesSolid(6,3,L)
            !      end do ! loop over loaded faces
            !        
            !        
            !        
            !        
            !        !---------------------------------------------
            !        
            !        
            !      !allocate(LoadOnNodesConnectivitiesSolid(N_BOUNDARY_NODES_HOE, NLoadedElementSidesSolidNodes), stat = IError)
            !      !allocate(LoadValuesOnNodesSolid(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesSolidNodes), stat = IError)
            !      !do L = 1, NLoadedElementSidesSolidNodes ! loop over loaded faces
            !      !  read(GOMunit,*) LoadOnNodesConnectivitiesSolid(1,L), LoadOnNodesConnectivitiesSolid(3,L),  &
            !      !                  LoadOnNodesConnectivitiesSolid(5,L), LoadOnNodesConnectivitiesSolid(2,L), &
            !      !                  LoadOnNodesConnectivitiesSolid(4,L), LoadOnNodesConnectivitiesSolid(6,L), &
            !      !                  LoadValuesOnNodesSolid(1,1,L), LoadValuesOnNodesSolid(1,2,L), LoadValuesOnNodesSolid(1,3,L), &
            !      !                  LoadValuesOnNodesSolid(3,1,L), LoadValuesOnNodesSolid(3,2,L), LoadValuesOnNodesSolid(3,3,L), &
            !      !                  LoadValuesOnNodesSolid(5,1,L), LoadValuesOnNodesSolid(5,2,L), LoadValuesOnNodesSolid(5,3,L), &
            !      !                  LoadValuesOnNodesSolid(2,1,L), LoadValuesOnNodesSolid(2,2,L), LoadValuesOnNodesSolid(2,3,L), &
            !      !                  LoadValuesOnNodesSolid(4,1,L), LoadValuesOnNodesSolid(4,2,L), LoadValuesOnNodesSolid(4,3,L), &
            !      !                  LoadValuesOnNodesSolid(6,1,L), LoadValuesOnNodesSolid(6,2,L), LoadValuesOnNodesSolid(6,3,L)
            !      !end do ! loop over loaded faces
            !    else if ( NDIM == 2 ) then
            !
            !        ! hardcoded for xi direction traction and needs generalization
            !      allocate(LoadOnNodesConnectivitiesSolid(ELEMENTBOUNDARYNODES, NLoadedElementSidesSolidNodes), stat = IError) 
            !      
            !      allocate(LoadValuesOnNodesSolid(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesSolidNodes), stat = IError)
            !
            !      do L = 1, NLoadedElementSidesSolidNodes ! loop over loaded element sides
            !
            !        read(GOMunit,*) LoadOnNodesConnectivitiesSolid(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnNodesSolid(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
            !
            !      end do ! loop over loaded element sides  
            !    end if
            !  end if
              
              
              
            else if ((trim(BName)=='$$START_LOAD_ON_NODES_SOLID_B')) then
              read(GOMunit,*) NLoadedElementSidesSolidNodesB ! Loaded faces
              if (NLoadedElementSidesSolidNodesB > 0) then
                  Counters%NSoilLoadSystems = 2
                if ( NDIM == 3 ) then
                  allocate(LoadOnNodesConnectivitiesSolidB(N_BOUNDARY_NODES_HOE, NLoadedElementSidesSolidNodesB), stat = IError)
                  allocate(LoadValuesOnNodesSolidB(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesSolidNodesB), stat = IError)
                  do L = 1, NLoadedElementSidesSolidNodesB ! loop over loaded faces
                    read(GOMunit,*) LoadOnNodesConnectivitiesSolidB(1,L), LoadOnNodesConnectivitiesSolidB(3,L),  &
                                    LoadOnNodesConnectivitiesSolidB(5,L), LoadOnNodesConnectivitiesSolidB(2,L), &
                                    LoadOnNodesConnectivitiesSolidB(4,L), LoadOnNodesConnectivitiesSolidB(6,L), &
                                    LoadValuesOnNodesSolidB(1,1,L), LoadValuesOnNodesSolidB(1,2,L), LoadValuesOnNodesSolidB(1,3,L), &
                                    LoadValuesOnNodesSolidB(3,1,L), LoadValuesOnNodesSolidB(3,2,L), LoadValuesOnNodesSolidB(3,3,L), &
                                    LoadValuesOnNodesSolidB(5,1,L), LoadValuesOnNodesSolidB(5,2,L), LoadValuesOnNodesSolidB(5,3,L), &
                                    LoadValuesOnNodesSolidB(2,1,L), LoadValuesOnNodesSolidB(2,2,L), LoadValuesOnNodesSolidB(2,3,L), &
                                    LoadValuesOnNodesSolidB(4,1,L), LoadValuesOnNodesSolidB(4,2,L), LoadValuesOnNodesSolidB(4,3,L), &
                                    LoadValuesOnNodesSolidB(6,1,L), LoadValuesOnNodesSolidB(6,2,L), LoadValuesOnNodesSolidB(6,3,L)
                  end do ! loop over loaded faces
                else if ( NDIM == 2 ) then
                  allocate(LoadOnNodesConnectivitiesSolidB(ELEMENTBOUNDARYNODES, NLoadedElementSidesSolidNodesB), stat = IError)
                  allocate(LoadValuesOnNodesSolidB(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesSolidNodesB), stat = IError)
                  do L = 1, NLoadedElementSidesSolidNodesB ! loop over loaded element sides
                    read(GOMunit,*) LoadOnNodesConnectivitiesSolidB(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnNodesSolidB(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
                  end do ! loop over loaded element sides  
                end if
              end if  

    !        else if (trim(BName)=='$$START_LOAD_ON_MATERIAL_POINTS_SOLID') then ! Traction boundary condition in GOM file starts here
    !          read(GOMunit,*) NLoadedElementSidesSolidMatPoints ! Loaded faces
    !          if (NLoadedElementSidesSolidMatPoints > 0) then
    !            if (.not.IsMPMComputation()) then
    !              call GiveWarning('A load on material points (solid) is specified while not using MPM.')
    !            end if
				!if ( NDIM == 3 ) then
    !            	allocate(LoadOnMatPointsConnectivitiesSolid(N_BOUNDARY_NODES_HOE, NLoadedElementSidesSolidMatPoints), stat = IError)
    !            	allocate(LoadValuesOnMatPointsSolid(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesSolidMatPoints), stat = IError)
    !            	do L = 1, NLoadedElementSidesSolidMatPoints ! loop over loaded faces
    !              	read(GOMunit,*) LoadOnMatPointsConnectivitiesSolid(1,L), LoadOnMatPointsConnectivitiesSolid(3,L),  &
    !                              	LoadOnMatPointsConnectivitiesSolid(5,L), LoadOnMatPointsConnectivitiesSolid(2,L), &
    !                              	LoadOnMatPointsConnectivitiesSolid(4,L), LoadOnMatPointsConnectivitiesSolid(6,L), &
    !                        		LoadValuesOnMatPointsSolid(1,1,L), LoadValuesOnMatPointsSolid(1,2,L), LoadValuesOnMatPointsSolid(1,3,L), &
    !                        		LoadValuesOnMatPointsSolid(3,1,L), LoadValuesOnMatPointsSolid(3,2,L), LoadValuesOnMatPointsSolid(3,3,L), &
    !                        		LoadValuesOnMatPointsSolid(5,1,L), LoadValuesOnMatPointsSolid(5,2,L), LoadValuesOnMatPointsSolid(5,3,L), &
    !                        		LoadValuesOnMatPointsSolid(2,1,L), LoadValuesOnMatPointsSolid(2,2,L), LoadValuesOnMatPointsSolid(2,3,L), &
    !                        		LoadValuesOnMatPointsSolid(4,1,L), LoadValuesOnMatPointsSolid(4,2,L), LoadValuesOnMatPointsSolid(4,3,L), &
    !                        		LoadValuesOnMatPointsSolid(6,1,L), LoadValuesOnMatPointsSolid(6,2,L), LoadValuesOnMatPointsSolid(6,3,L)
    !            	end do ! loop over loaded faces
				!else if ( NDIM == 2 ) then
    !              allocate(LoadOnMatPointsConnectivitiesSolid(ELEMENTBOUNDARYNODES, NLoadedElementSidesSolidMatPoints), stat = IError)
    !              allocate(LoadValuesOnMatPointsSolid(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesSolidMatPoints), stat = IError)
    !              do L = 1, NLoadedElementSidesSolidMatPoints ! loop over loaded element sides
    !                read(GOMunit,*) LoadOnMatPointsConnectivitiesSolid(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnMatPointsSolid(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
    !              end do ! loop over loaded element sides
    !            end if
    !          end if
              
             else if (trim(BName)=='$$START_LOAD_ON_MATERIAL_POINTS_SOLID_B') then
              read(GOMunit,*) NLoadedElementSidesSolidMatPointsB ! Loaded faces
              if (NLoadedElementSidesSolidMatPointsB > 0) then
                  Counters%NSoilLoadSystems = 2
                if (.not.IsMPMComputation()) then
                  call GiveWarning('A load on material points (solid) is specified while not using MPM.')
                end if
				if ( NDIM == 3 ) then
                	allocate(LoadOnMatPointsConnectivitiesSolidB(N_BOUNDARY_NODES_HOE, NLoadedElementSidesSolidMatPointsB), stat = IError)
                	allocate(LoadValuesOnMatPointsSolidB(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesSolidMatPointsB), stat = IError)
                	do L = 1, NLoadedElementSidesSolidMatPointsB ! loop over loaded faces
                  	read(GOMunit,*) LoadOnMatPointsConnectivitiesSolidB(1,L), LoadOnMatPointsConnectivitiesSolidB(3,L),  &
                                  	LoadOnMatPointsConnectivitiesSolidB(5,L), LoadOnMatPointsConnectivitiesSolidB(2,L), &
                                  	LoadOnMatPointsConnectivitiesSolidB(4,L), LoadOnMatPointsConnectivitiesSolidB(6,L), &
                            		LoadValuesOnMatPointsSolidB(1,1,L), LoadValuesOnMatPointsSolidB(1,2,L), LoadValuesOnMatPointsSolidB(1,3,L), &
                            		LoadValuesOnMatPointsSolidB(3,1,L), LoadValuesOnMatPointsSolidB(3,2,L), LoadValuesOnMatPointsSolidB(3,3,L), &
                            		LoadValuesOnMatPointsSolidB(5,1,L), LoadValuesOnMatPointsSolidB(5,2,L), LoadValuesOnMatPointsSolidB(5,3,L), &
                            		LoadValuesOnMatPointsSolidB(2,1,L), LoadValuesOnMatPointsSolidB(2,2,L), LoadValuesOnMatPointsSolidB(2,3,L), &
                            		LoadValuesOnMatPointsSolidB(4,1,L), LoadValuesOnMatPointsSolidB(4,2,L), LoadValuesOnMatPointsSolidB(4,3,L), &
                            		LoadValuesOnMatPointsSolidB(6,1,L), LoadValuesOnMatPointsSolidB(6,2,L), LoadValuesOnMatPointsSolidB(6,3,L)
                	end do ! loop over loaded faces
				else if ( NDIM == 2 ) then
                  allocate(LoadOnMatPointsConnectivitiesSolidB(ELEMENTBOUNDARYNODES, NLoadedElementSidesSolidMatPointsB), stat = IError)
                  allocate(LoadValuesOnMatPointsSolidB(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesSolidMatPointsB), stat = IError)
                  do L = 1, NLoadedElementSidesSolidMatPointsB ! loop over loaded element sides
                    read(GOMunit,*) LoadOnMatPointsConnectivitiesSolidB(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnMatPointsSolidB(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
                  end do ! loop over loaded element sides
                end if
              end if 

            !else if ((trim(BName)=='$$START_LOAD_ON_NODES_LIQUID')) then
            !  read(GOMunit,*) NLoadedElementSidesWaterNodes ! Loaded faces
            !  if (NLoadedElementSidesWaterNodes > 0) then
            !      if (NDIM == 3) then 
            !        allocate(LoadOnNodesConnectivitiesWater(N_BOUNDARY_NODES_HOE, NLoadedElementSidesWaterNodes), stat = IError)
            !        allocate(LoadValuesOnNodesWater(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesWaterNodes), stat = IError)
            !        do L = 1, NLoadedElementSidesWaterNodes ! loop over loaded faces
            !            read(GOMunit,*) LoadOnNodesConnectivitiesWater(1,L), LoadOnNodesConnectivitiesWater(3,L),  &
            !                    LoadOnNodesConnectivitiesWater(5,L), LoadOnNodesConnectivitiesWater(2,L), &
            !                    LoadOnNodesConnectivitiesWater(4,L), LoadOnNodesConnectivitiesWater(6,L), &
            !                    LoadValuesOnNodesWater(1,1,L), LoadValuesOnNodesWater(1,2,L), LoadValuesOnNodesWater(1,3,L), &
            !                    LoadValuesOnNodesWater(3,1,L), LoadValuesOnNodesWater(3,2,L), LoadValuesOnNodesWater(3,3,L), &
            !                    LoadValuesOnNodesWater(5,1,L), LoadValuesOnNodesWater(5,2,L), LoadValuesOnNodesWater(5,3,L), &
            !                    LoadValuesOnNodesWater(2,1,L), LoadValuesOnNodesWater(2,2,L), LoadValuesOnNodesWater(2,3,L), &
            !                    LoadValuesOnNodesWater(4,1,L), LoadValuesOnNodesWater(4,2,L), LoadValuesOnNodesWater(4,3,L), &
            !                    LoadValuesOnNodesWater(6,1,L), LoadValuesOnNodesWater(6,2,L), LoadValuesOnNodesWater(6,3,L)
            !        end do ! loop over loaded faces
            !     else if ( NDIM == 2 ) then
            !        allocate(LoadOnNodesConnectivitiesWater(ELEMENTBOUNDARYNODES, NLoadedElementSidesWaterNodes), stat = IError)
            !        allocate(LoadValuesOnNodesWater(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesWaterNodes), stat = IError)
            !        do L = 1, NLoadedElementSidesWaterNodes ! loop over loaded element sides
            !          read(GOMunit,*) LoadOnNodesConnectivitiesWater(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnNodesWater(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
            !        end do ! loop over loaded element sides  
            !    end if
            !  end if
              else if ((trim(BName)=='$$START_LOAD_ON_NODES_LIQUID_B')) then
              read(GOMunit,*) NLoadedElementSidesWaterNodesB ! Loaded faces
              if (NLoadedElementSidesWaterNodesB > 0) then
                  Counters%NWaterLoadSystems = 2
                  if (NDIM == 3) then 
                    allocate(LoadOnNodesConnectivitiesWaterB(N_BOUNDARY_NODES_HOE, NLoadedElementSidesWaterNodesB), stat = IError)
                    allocate(LoadValuesOnNodesWaterB(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesWaterNodesB), stat = IError)
                    do L = 1, NLoadedElementSidesWaterNodesB ! loop over loaded faces
                        read(GOMunit,*) LoadOnNodesConnectivitiesWaterB(1,L), LoadOnNodesConnectivitiesWaterB(3,L),  &
                                LoadOnNodesConnectivitiesWaterB(5,L), LoadOnNodesConnectivitiesWaterB(2,L), &
                                LoadOnNodesConnectivitiesWaterB(4,L), LoadOnNodesConnectivitiesWaterB(6,L), &
                                LoadValuesOnNodesWaterB(1,1,L), LoadValuesOnNodesWaterB(1,2,L), LoadValuesOnNodesWaterB(1,3,L), &
                                LoadValuesOnNodesWaterB(3,1,L), LoadValuesOnNodesWaterB(3,2,L), LoadValuesOnNodesWaterB(3,3,L), &
                                LoadValuesOnNodesWaterB(5,1,L), LoadValuesOnNodesWaterB(5,2,L), LoadValuesOnNodesWaterB(5,3,L), &
                                LoadValuesOnNodesWaterB(2,1,L), LoadValuesOnNodesWaterB(2,2,L), LoadValuesOnNodesWaterB(2,3,L), &
                                LoadValuesOnNodesWaterB(4,1,L), LoadValuesOnNodesWaterB(4,2,L), LoadValuesOnNodesWaterB(4,3,L), &
                                LoadValuesOnNodesWaterB(6,1,L), LoadValuesOnNodesWaterB(6,2,L), LoadValuesOnNodesWaterB(6,3,L)
                    end do ! loop over loaded faces
                 else if ( NDIM == 2 ) then
                    allocate(LoadOnNodesConnectivitiesWaterB(ELEMENTBOUNDARYNODES, NLoadedElementSidesWaterNodesB), stat = IError)
                    allocate(LoadValuesOnNodesWaterB(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesWaterNodesB), stat = IError)
                    do L = 1, NLoadedElementSidesWaterNodesB ! loop over loaded element sides
                      read(GOMunit,*) LoadOnNodesConnectivitiesWaterB(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnNodesWaterB(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
                    end do ! loop over loaded element sides  
                end if
              end if

            !else if ((trim(BName)=='$$START_LOAD_ON_MATERIAL_POINTS_LIQUID')) then
            !  read(GOMunit,*) NLoadedElementSidesWaterMatPoints ! Loaded faces
            !  if (NLoadedElementSidesWaterMatPoints > 0) then
            !      
            !    if (.not.IsMPMComputation()) then
            !      call GiveWarning('A load on material points (liquid) is specified while not using MPM')
            !    end if
            !  if (NDIM == 3) then 
            !    allocate(LoadOnMatPointsConnectivitiesWater(N_BOUNDARY_NODES_HOE, NLoadedElementSidesWaterMatPoints), stat = IError)
            !    allocate(LoadValuesOnMatPointsWater(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesWaterMatPoints), stat = IError)
            !    do L = 1, NLoadedElementSidesWaterMatPoints ! loop over loaded faces
            !      read(GOMunit,*) LoadOnMatPointsConnectivitiesWater(1,L), LoadOnMatPointsConnectivitiesWater(3,L),  &
            !                      LoadOnMatPointsConnectivitiesWater(5,L), LoadOnMatPointsConnectivitiesWater(2,L), &
            !                      LoadOnMatPointsConnectivitiesWater(4,L), LoadOnMatPointsConnectivitiesWater(6,L), &
            !                LoadValuesOnMatPointsWater(1,1,L), LoadValuesOnMatPointsWater(1,2,L), LoadValuesOnMatPointsWater(1,3,L), &
            !                LoadValuesOnMatPointsWater(3,1,L), LoadValuesOnMatPointsWater(3,2,L), LoadValuesOnMatPointsWater(3,3,L), &
            !                LoadValuesOnMatPointsWater(5,1,L), LoadValuesOnMatPointsWater(5,2,L), LoadValuesOnMatPointsWater(5,3,L), &
            !                LoadValuesOnMatPointsWater(2,1,L), LoadValuesOnMatPointsWater(2,2,L), LoadValuesOnMatPointsWater(2,3,L), &
            !                LoadValuesOnMatPointsWater(4,1,L), LoadValuesOnMatPointsWater(4,2,L), LoadValuesOnMatPointsWater(4,3,L), &
            !                LoadValuesOnMatPointsWater(6,1,L), LoadValuesOnMatPointsWater(6,2,L), LoadValuesOnMatPointsWater(6,3,L)
            !    end do ! loop over loaded faces
            !  else if ( NDIM == 2 ) then
            !    allocate(LoadOnMatPointsConnectivitiesWater(ELEMENTBOUNDARYNODES, NLoadedElementSidesWaterMatPoints), stat = IError)
            !    allocate(LoadValuesOnMatPointsWater(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesWaterMatPoints), stat = IError)
            !      do L = 1, NLoadedElementSidesWaterMatPoints ! loop over loaded element sides
            !        read(GOMunit,*) LoadOnMatPointsConnectivitiesWater(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnMatPointsWater(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
            !      end do ! loop over loaded element sides
            !  end if
            !  end if
              
              else if ((trim(BName)=='$$START_LOAD_ON_MATERIAL_POINTS_LIQUID_B')) then
              read(GOMunit,*) NLoadedElementSidesWaterMatPointsB ! Loaded faces
              if (NLoadedElementSidesWaterMatPointsB > 0) then
                  Counters%NWaterLoadSystems = 2
                if (.not.IsMPMComputation()) then
                  call GiveWarning('A load on material points (liquid) is specified while not using MPM')
                end if
              if (NDIM == 3) then 
                allocate(LoadOnMatPointsConnectivitiesWaterB(N_BOUNDARY_NODES_HOE, NLoadedElementSidesWaterMatPointsB), stat = IError)
                allocate(LoadValuesOnMatPointsWaterB(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesWaterMatPointsB), stat = IError)
                do L = 1, NLoadedElementSidesWaterMatPointsB ! loop over loaded faces
                  read(GOMunit,*) LoadOnMatPointsConnectivitiesWaterB(1,L), LoadOnMatPointsConnectivitiesWaterB(3,L),  &
                                  LoadOnMatPointsConnectivitiesWaterB(5,L), LoadOnMatPointsConnectivitiesWaterB(2,L), &
                                  LoadOnMatPointsConnectivitiesWaterB(4,L), LoadOnMatPointsConnectivitiesWaterB(6,L), &
                            LoadValuesOnMatPointsWaterB(1,1,L), LoadValuesOnMatPointsWaterB(1,2,L), LoadValuesOnMatPointsWaterB(1,3,L), &
                            LoadValuesOnMatPointsWaterB(3,1,L), LoadValuesOnMatPointsWaterB(3,2,L), LoadValuesOnMatPointsWaterB(3,3,L), &
                            LoadValuesOnMatPointsWaterB(5,1,L), LoadValuesOnMatPointsWaterB(5,2,L), LoadValuesOnMatPointsWaterB(5,3,L), &
                            LoadValuesOnMatPointsWaterB(2,1,L), LoadValuesOnMatPointsWaterB(2,2,L), LoadValuesOnMatPointsWaterB(2,3,L), &
                            LoadValuesOnMatPointsWaterB(4,1,L), LoadValuesOnMatPointsWaterB(4,2,L), LoadValuesOnMatPointsWaterB(4,3,L), &
                            LoadValuesOnMatPointsWaterB(6,1,L), LoadValuesOnMatPointsWaterB(6,2,L), LoadValuesOnMatPointsWaterB(6,3,L)
                end do ! loop over loaded faces
              else if ( NDIM == 2 ) then
                allocate(LoadOnMatPointsConnectivitiesWaterB(ELEMENTBOUNDARYNODES, NLoadedElementSidesWaterMatPointsB), stat = IError)
                allocate(LoadValuesOnMatPointsWaterB(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesWaterMatPointsB), stat = IError)
                  do L = 1, NLoadedElementSidesWaterMatPointsB ! loop over loaded element sides
                    read(GOMunit,*) LoadOnMatPointsConnectivitiesWaterB(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnMatPointsWaterB(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
                  end do ! loop over loaded element sides
              end if
            end if
              
            else if (trim(BName)=='$$START_LOAD_ON_NODES_GAS') then
              read(GOMunit,*) NLoadedElementSidesGasNodes ! Loaded faces
              if (NLoadedElementSidesGasNodes > 0) then
              allocate(LoadOnNodesConnectivitiesGas(N_BOUNDARY_NODES_HOE, NLoadedElementSidesGasNodes), stat = IError)
              allocate(LoadValuesOnNodesGas(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesGasNodes), stat = IError)
              do L = 1, NLoadedElementSidesGasNodes ! loop over loaded faces
                read(GOMunit,*) LoadOnNodesConnectivitiesGas(1,L), LoadOnNodesConnectivitiesGas(3,L),  &
                                LoadOnNodesConnectivitiesGas(5,L), LoadOnNodesConnectivitiesGas(2,L), &
                                LoadOnNodesConnectivitiesGas(4,L), LoadOnNodesConnectivitiesGas(6,L), &
                                LoadValuesOnNodesGas(1,1,L), LoadValuesOnNodesGas(1,2,L), LoadValuesOnNodesGas(1,3,L), &
                                LoadValuesOnNodesGas(3,1,L), LoadValuesOnNodesGas(3,2,L), LoadValuesOnNodesGas(3,3,L), &
                                LoadValuesOnNodesGas(5,1,L), LoadValuesOnNodesGas(5,2,L), LoadValuesOnNodesGas(5,3,L), &
                                LoadValuesOnNodesGas(2,1,L), LoadValuesOnNodesGas(2,2,L), LoadValuesOnNodesGas(2,3,L), &
                                LoadValuesOnNodesGas(4,1,L), LoadValuesOnNodesGas(4,2,L), LoadValuesOnNodesGas(4,3,L), &
                                LoadValuesOnNodesGas(6,1,L), LoadValuesOnNodesGas(6,2,L), LoadValuesOnNodesGas(6,3,L)
              end do ! loop over loaded faces
              end if
              
              else if (trim(BName)=='$$START_LOAD_ON_NODES_GAS_B') then
              read(GOMunit,*) NLoadedElementSidesGasNodesB ! Loaded faces
              if (NLoadedElementSidesGasNodesB > 0) then
                  Counters%NGasLoadSystems = 2
              allocate(LoadOnNodesConnectivitiesGasB(N_BOUNDARY_NODES_HOE, NLoadedElementSidesGasNodes), stat = IError)
              allocate(LoadValuesOnNodesGasB(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesGasNodes), stat = IError)
              do L = 1, NLoadedElementSidesGasNodesB ! loop over loaded faces
                read(GOMunit,*) LoadOnNodesConnectivitiesGasB(1,L), LoadOnNodesConnectivitiesGasB(3,L),  &
                                LoadOnNodesConnectivitiesGasB(5,L), LoadOnNodesConnectivitiesGasB(2,L), &
                                LoadOnNodesConnectivitiesGasB(4,L), LoadOnNodesConnectivitiesGasB(6,L), &
                                LoadValuesOnNodesGasB(1,1,L), LoadValuesOnNodesGasB(1,2,L), LoadValuesOnNodesGasB(1,3,L), &
                                LoadValuesOnNodesGasB(3,1,L), LoadValuesOnNodesGasB(3,2,L), LoadValuesOnNodesGasB(3,3,L), &
                                LoadValuesOnNodesGasB(5,1,L), LoadValuesOnNodesGasB(5,2,L), LoadValuesOnNodesGasB(5,3,L), &
                                LoadValuesOnNodesGasB(2,1,L), LoadValuesOnNodesGasB(2,2,L), LoadValuesOnNodesGasB(2,3,L), &
                                LoadValuesOnNodesGasB(4,1,L), LoadValuesOnNodesGasB(4,2,L), LoadValuesOnNodesGasB(4,3,L), &
                                LoadValuesOnNodesGasB(6,1,L), LoadValuesOnNodesGasB(6,2,L), LoadValuesOnNodesGasB(6,3,L)
              end do ! loop over loaded faces
              end if
              
              
            else if (trim(BName)=='$$START_LOAD_ON_MATERIAL_POINTS_GAS') then
              read(GOMunit,*) NLoadedElementSidesGasMatPoints ! Loaded faces
              if (NLoadedElementSidesGasMatPoints > 0) then
                if (.not.IsMPMComputation()) then
                  call GiveWarning('A load on material points (gas) is specified while not using MPM')
                end if
                allocate(LoadOnMatPointsConnectivitiesGas(N_BOUNDARY_NODES_HOE, NLoadedElementSidesGasMatPoints), stat = IError)
                allocate(LoadValuesOnMatPointsGas(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesGasMatPoints), stat = IError)
                do L = 1, NLoadedElementSidesGasMatPoints ! loop over loaded faces
                  read(GOMunit,*) LoadOnMatPointsConnectivitiesGas(1,L), LoadOnMatPointsConnectivitiesGas(3,L), &
                                  LoadOnMatPointsConnectivitiesGas(5,L), LoadOnMatPointsConnectivitiesGas(2,L), &
                                  LoadOnMatPointsConnectivitiesGas(4,L), LoadOnMatPointsConnectivitiesGas(6,L), &
                                  LoadValuesOnMatPointsGas(1,1,L), LoadValuesOnMatPointsGas(1,2,L), LoadValuesOnMatPointsGas(1,3,L), &
                                  LoadValuesOnMatPointsGas(3,1,L), LoadValuesOnMatPointsGas(3,2,L), LoadValuesOnMatPointsGas(3,3,L), &
                                  LoadValuesOnMatPointsGas(5,1,L), LoadValuesOnMatPointsGas(5,2,L), LoadValuesOnMatPointsGas(5,3,L), &
                                  LoadValuesOnMatPointsGas(2,1,L), LoadValuesOnMatPointsGas(2,2,L), LoadValuesOnMatPointsGas(2,3,L), &
                                  LoadValuesOnMatPointsGas(4,1,L), LoadValuesOnMatPointsGas(4,2,L), LoadValuesOnMatPointsGas(4,3,L), &
                                  LoadValuesOnMatPointsGas(6,1,L), LoadValuesOnMatPointsGas(6,2,L), LoadValuesOnMatPointsGas(6,3,L)
                end do ! loop over loaded faces
              end if
              
             else if (trim(BName)=='$$START_LOAD_ON_MATERIAL_POINTS_GAS_B') then
              read(GOMunit,*) NLoadedElementSidesGasMatPointsB ! Loaded faces
              if (NLoadedElementSidesGasMatPointsB > 0) then
                  Counters%NGasLoadSystems = 2
                if (.not.IsMPMComputation()) then
                  call GiveWarning('A load on material points (gas) is specified while not using MPM')
                end if
                allocate(LoadOnMatPointsConnectivitiesGasB(N_BOUNDARY_NODES_HOE, NLoadedElementSidesGasMatPointsB), stat = IError)
                allocate(LoadValuesOnMatPointsGasB(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesGasMatPointsB), stat = IError)
                do L = 1, NLoadedElementSidesGasMatPointsB ! loop over loaded faces
                  read(GOMunit,*) LoadOnMatPointsConnectivitiesGasB(1,L), LoadOnMatPointsConnectivitiesGasB(3,L), &
                                  LoadOnMatPointsConnectivitiesGasB(5,L), LoadOnMatPointsConnectivitiesGasB(2,L), &
                                  LoadOnMatPointsConnectivitiesGasB(4,L), LoadOnMatPointsConnectivitiesGasB(6,L), &
                                  LoadValuesOnMatPointsGasB(1,1,L), LoadValuesOnMatPointsGasB(1,2,L), LoadValuesOnMatPointsGasB(1,3,L), &
                                  LoadValuesOnMatPointsGasB(3,1,L), LoadValuesOnMatPointsGasB(3,2,L), LoadValuesOnMatPointsGasB(3,3,L), &
                                  LoadValuesOnMatPointsGasB(5,1,L), LoadValuesOnMatPointsGasB(5,2,L), LoadValuesOnMatPointsGasB(5,3,L), &
                                  LoadValuesOnMatPointsGasB(2,1,L), LoadValuesOnMatPointsGasB(2,2,L), LoadValuesOnMatPointsGasB(2,3,L), &
                                  LoadValuesOnMatPointsGasB(4,1,L), LoadValuesOnMatPointsGasB(4,2,L), LoadValuesOnMatPointsGasB(4,3,L), &
                                  LoadValuesOnMatPointsGasB(6,1,L), LoadValuesOnMatPointsGasB(6,2,L), LoadValuesOnMatPointsGasB(6,3,L)
                end do ! loop over loaded faces
              end if
              

              
            else if (trim(BName)=='$$START_SOIL_SURFACE_NODES') then
                read(GOMunit,*) SoilSurfaceNumberofSides
                if (SoilSurfaceNumberofSides > 0) then
                    if (NDIM == 3) then
                        allocate(SoilSurfaceNodesConnectivities(N_BOUNDARY_NODES_HOE, SoilSurfaceNumberofSides), stat = IError)
                        do L = 1, SoilSurfaceNumberofSides ! loop over soil surface faces
                            read(GOMunit,*) SoilSurfaceNodesConnectivities(1,L), SoilSurfaceNodesConnectivities(3,L),  &
                                SoilSurfaceNodesConnectivities(5,L), SoilSurfaceNodesConnectivities(2,L), &
                                SoilSurfaceNodesConnectivities(4,L), SoilSurfaceNodesConnectivities(6,L)
                        end do ! loop over soil surface faces
                    else if ( NDIM == 2 ) then
                        allocate(SoilSurfaceNodesConnectivities(ELEMENTBOUNDARYNODES, SoilSurfaceNumberofSides), stat = IError)
                        do L = 1, SoilSurfaceNumberofSides
                            read(GOMunit,*) SoilSurfaceNodesConnectivities(1:ELEMENTBOUNDARYNODES, L)
                        end do
                    end if
                end if
                
                  else if (trim(BName)=='$$START_PHREATIC_SURFACE_NODES') then
                read(GOMunit,*) PhreaticSurfaceNumberofSides
                if (PhreaticSurfaceNumberofSides > 0) then
                    if (NDIM == 3) then
                        allocate(PhreaticSurfaceNodesConnectivities(N_BOUNDARY_NODES_HOE, PhreaticSurfaceNumberofSides), stat = IError)
                        do L = 1, PhreaticSurfaceNumberofSides ! loop over phreatic surface faces
                            read(GOMunit,*) PhreaticSurfaceNodesConnectivities(1,L), PhreaticSurfaceNodesConnectivities(3,L),  &
                                PhreaticSurfaceNodesConnectivities(5,L), PhreaticSurfaceNodesConnectivities(2,L), &
                                PhreaticSurfaceNodesConnectivities(4,L), PhreaticSurfaceNodesConnectivities(6,L)
                        end do ! loop over phreatic surface faces
                    else if ( NDIM == 2 ) then
                        allocate(PhreaticSurfaceNodesConnectivities(ELEMENTBOUNDARYNODES, PhreaticSurfaceNumberofSides), stat = IError)
                        do L = 1, PhreaticSurfaceNumberofSides
                            read(GOMunit,*) PhreaticSurfaceNodesConnectivities(1:ELEMENTBOUNDARYNODES, L)
                        end do
                    end if
                end if      
                
                              
                  else if (trim(BName)=='$$BOUNDARY_HYDRAULIC_HEAD_AREA') then
                      if (NDIM == 3) then
                          allocate(HydrHeadArea(NDIM,2), stat = IError)
                          do L=1, 2! loop over x,y,z bounding coordinate
                              read(GOMunit,*) HydrHeadArea(1,L), HydrHeadArea(2,L), HydrHeadArea(3,L)
                          end do
                      else if (NDIM == 2) then
                          allocate(HydrHeadArea(NDIM,2), stat = IError)
                          do L=1, 2! loop over x,y,z bounding coordinate
                              read(GOMunit,*) HydrHeadArea(1,L), HydrHeadArea(2,L)
                          end do
                      end if
           
                    
          else if (trim(BName)=='$$BOUNDARY_SEEPAGE_AREA') then
                  if (NDIM == 3) then
                          allocate(SeepageArea(NDIM,2), stat = IError)
                          do L=1, 2! loop over x,y,z bounding coordinate
                              read(GOMunit,*) SeepageArea(1,L), SeepageArea(2,L), SeepageArea(3,L)
                          end do
                      else if (NDIM == 2) then
                          allocate(SeepageArea(NDIM,2), stat = IError)
                          do L=1, 2! loop over x,y,z bounding coordinate
                              read(GOMunit,*) SeepageArea(1,L), SeepageArea(2,L)
                          end do
                      end if
            
        
      
      
        else if (trim(BName)=='$$BOUNDARY_INFILTRATION_AREA') then
         if (NDIM == 3) then
                          allocate(InfiltrationArea(NDIM,2), stat = IError)
                          do L=1, 2! loop over x,y,z bounding coordinate
                              read(GOMunit,*) InfiltrationArea(1,L), InfiltrationArea(2,L), InfiltrationArea(3,L)
                          end do
                      else if (NDIM == 2) then
                          allocate(InfiltrationArea(NDIM,2), stat = IError)
                          do L=1, 2! loop over x,y,z bounding coordinate
                              read(GOMunit,*) InfiltrationArea(1,L), InfiltrationArea(2,L)
                          end do
                      end if
               
        
      else if (trim(BName)=='$$INFILTRATION_RATE') then
          allocate(InfiltrationRate(NDIM), stat = IError)
          !do L=1, NDIM
              read(GOMunit,*) InfiltrationRate(1:NDIM)
           !   CalParams%BoundaryConditions%InfiltrationRate = InfiltrationRate
          !end do

        
                       
      else if (trim(BName)=='$$STARTELMMAT') then 
          
          ! I need to warn the user here if MP is from an external file then we 
          ! do not need this flag here
          
          !do IPatch = 1, Counters%NPatches ! loop over patches
          do IDElement = 1, Counters%NEl(IPatch)!NURBS%nel_NURBS(IPatch)!NElements ! loop over elements
                read (GOMunit,*) ElementMaterialID(IDElement,IPatch) ! get the material ID of the element
                if (ElementMaterialID(IDElement,IPatch)==0) then
                  ElementMaterialID(IDElement,IPatch) = -1 ! set the initially deactivated elements
                end if
          end do ! loop over elements  
          !end do ! loop over patches
          
          
      else if (trim(BName)=='$$START_BOUND_MP_PRESCRIBED_VELOCITY_X') then 
          BoundCounter_X = BoundCounter_X + 1
          read (GOMunit,*) BoundPrescVelo_X(BoundCounter_X,1), BoundPrescVelo_X(BoundCounter_X,2), &
                            BoundPrescVelo_X(BoundCounter_X,3), BoundPrescVelo_X(BoundCounter_X,4), BoundPrescVelo_X(BoundCounter_X,5)
          !read (GOMunit,*) BoundHigherOrLower, CoordinateLimit, Velocity_X, Velocity_Y, Velocity_Z
          ! 1- BoundHigherOrLower --> 1 then higher than the bound
          !                    --> 0 then lower than the bound 
          ! 2- Locus in respective direction 
          ! 3- Velocity X
          ! 4- Velocity Y
          ! 5- Velocity Z
          

        
          
          
          
      else if (trim(BName)=='$$START_BOUND_MP_PRESCRIBED_VELOCITY_Y') then     
          BoundCounter_Y = BoundCounter_Y + 1
          read (GOMunit,*) BoundPrescVelo_Y(BoundCounter_Y,1), BoundPrescVelo_Y(BoundCounter_Y,2), &
                            BoundPrescVelo_Y(BoundCounter_Y,3), BoundPrescVelo_Y(BoundCounter_Y,4), BoundPrescVelo_Y(BoundCounter_Y,5)
          !read (GOMunit,*) BoundHigherOrLower, CoordinateLimit, Velocity_X, Velocity_Y, Velocity_Z
          ! 1- BoundHigherOrLower --> 1 then higher than the bound
          !                    --> 0 then lower than the bound 
          ! 2- Locus in respective direction 
          ! 3- Velocity X
          ! 4- Velocity Y
          ! 5- Velocity Z
          
          
          
          
          
      else if (trim(BName)=='$$START_BOUND_MP_PRESCRIBED_VELOCITY_Z') then 
          BoundCounter_Z = BoundCounter_Z + 1
          read (GOMunit,*) BoundPrescVelo_Z(BoundCounter_Z,1), BoundPrescVelo_Z(BoundCounter_Z,2), &
                            BoundPrescVelo_Z(BoundCounter_Z,3), BoundPrescVelo_Z(BoundCounter_Z,4), BoundPrescVelo_Z(BoundCounter_Z,5)
          !read (GOMunit,*) BoundHigherOrLower, CoordinateLimit, Velocity_X, Velocity_Y, Velocity_Z
          ! 1- BoundHigherOrLower --> 1 then higher than the bound
          !                    --> 0 then lower than the bound 
          ! 2- Locus in respective direction 
          ! 3- Velocity X
          ! 4- Velocity Y
          ! 5- Velocity Z
          
          
          
            else if (trim(BName) == '$$FINISH') then
                    EXIT
            end if 
            
          end do 
          
          
          
          ! THIS IS FOR APPLYING TRACTION ON THE NODES 
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do IPatch = 1, Counters%NPatches
                  do ILoadSystem = 1, NURBS%NLoadSystems_Traction(IPatch)
              
                    
                      call MatchControlPointsNodes(ILoadSystem, IPatch) ! traction on the control points 
              

                      !! THIS IS FOR APPLYING TRACTION ON THE MATERIAL POINTS 
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      call MatchControlPointsElements(ILoadSystem, IPatch) ! traction on the elements
                      call MatchTractionElementToActualMaterialSide(ILoadSystem, IPatch) ! compare traction side to all sides of actual material. 
                                                              ! We need to know side our traction is on. 
              !end if
              
                  end do 
              end do 
          
          
          
          
          
          
          
          close(GOMunit)

          ! we need to transfer these variabels into counters type
         ! NLoadedElementSidesSolidNodes = NURBS%NLoadSystems_Traction_NumberOfElements 
          
          call InitCountersFromFile(CalParams%ApplyContactAlgorithm, &
                                    Counters%Sum_NEl, Counters%Sum_NodTot, &!NElements, NNodes, & !NURBS%sum_nel_NURBS, NURBS%sum_NControlPoints
                                    CalParams%NumberOfMaterials, &
                                    NLoadedElementSidesSolidNodes,  &
                                    NLoadedElementSidesSolidNodesB,  &
                                    NLoadedElementSidesWaterNodes, &
                                    NLoadedElementSidesWaterNodesB, &
                                    NLoadedElementSidesGasNodes, &
                                    NLoadedElementSidesGasNodesB, &
                                    NLoadedElementSidesSolidMatPoints,  &
                                    NLoadedElementSidesSolidMatPointsB, &
                                    NLoadedElementSidesWaterMatPoints, &
                                    NLoadedElementSidesWaterMatPointsB, &
                                    NLoadedElementSidesGasMatPoints, &
                                    NLoadedElementSidesGasMatPointsB, &
                                    SoilSurfaceNumberofSides, &
                                    PhreaticSurfaceNumberofSides, &
                                    CalParams%NumberOfMaterials, &
                                    NURBS%Maximum_NLoadSystems_Traction)
          
          ! overwriting these two counters based on multipatch NURBS implementation
          !Counters%NEl = NURBS%sum_nel_NURBS
          !Counters%NodTot = NURBS%sum_NControlPoints
        
          if (CalParams%MPInitialization%MPInitFormat == MP_INIT_GAUSSPT) then 
              call InitialiseNumberOfActiveElements()
          end if 
          
          
          !call InitialiseEdgeNodeTyings()
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! I need to uncomment the below subroutines 
          ! need to know the bounds of the knots from the unique set of knot vectors 
          !
          
          ! Allocating the unique entries vectors based on the maximum read unique vectors for multiple patches and the number of patches
          allocate(NURBS%XiKnotEntriesUNIQUE(NURBS%Maximum_NXiKnotEntriesUNIQUE, Counters%NPatches) )
          allocate(NURBS%EtaKnotEntriesUNIQUE(NURBS%Maximum_NEtaKnotEntriesUNIQUE, Counters%NPatches) )
          allocate(NURBS%ZetaKnotEntriesUNIQUE(NURBS%Maximum_NZetaKnotEntriesUNIQUE, Counters%NPatches) )
          
          ! Populate XiKnotEntriesUNIQUE, EtaKnotEntriesUNIQUE, ZetaKnotEntriesUNIQUE for each patch 
          do IPatch = 1, Counters%NPatches ! loop over patches
              
              call FindUniqueEnteriesInVector(NURBS%XiKnotEntries(:,IPatch), NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntriesUNIQUE(:,IPatch))
              call FindUniqueEnteriesInVector(NURBS%EtaKnotEntries(:,IPatch), NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntriesUNIQUE(:,IPatch))
              call FindUniqueEnteriesInVector(NURBS%ZetaKnotEntries(:,IPatch), NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntriesUNIQUE(:,IPatch))
              
          end do ! loop over patches
          
          if (NDIM ==  2) then 
        
              NoOfItems = NURBS%Maximum_NXiKnotEntriesUNIQUE * &
                  NURBS%Maximum_NEtaKnotEntriesUNIQUE
              
          else if (NDIM == 3) then 
              
          NoOfItems = NURBS%Maximum_NXiKnotEntriesUNIQUE * &
                  NURBS%Maximum_NEtaKnotEntriesUNIQUE * &
                  NURBS%Maximum_NZetaKnotEntriesUNIQUE
        
          end if 
          
          ! initialize
          allocate(KnotBezierMesh(NoOfItems, NDIM, Counters%NPatches), stat = IError)
          KnotBezierMesh = 0.0
          
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          ! I need to generalize below bezier mesh for multiple patches!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (NDIM == 2) then 
              call FormBezierMeshUsingUniqueKnotVectors(NURBS%XiKnotEntriesUNIQUE(:,IPatch), NURBS%NXiKnotEntriesUNIQUE(IPatch), & 
                                                    NURBS%EtaKnotEntriesUNIQUE(:,IPatch), NURBS%NEtaKnotEntriesUNIQUE(IPatch), &
                                                    KnotBezierMesh(:,:,IPatch)) 
              
          else if (NDIM ==3) then 
              
              
              
              
              call FormBezierMeshUsingUniqueKnotVectors_3D(NURBS%XiKnotEntriesUNIQUE(:,IPatch), NURBS%NXiKnotEntriesUNIQUE(IPatch), & 
                                                    NURBS%EtaKnotEntriesUNIQUE(:,IPatch), NURBS%NEtaKnotEntriesUNIQUE(IPatch), &
                                                    NURBS%ZetaKnotEntriesUNIQUE(:,IPatch), NURBS%NZetaKnotEntriesUNIQUE(IPatch), &
                                                    IPatch)!KnotBezierMesh(:,:,IPatch)) 
              
          end if
          end do ! loop over patches
          
          
        
          
          if (NDIM == 2) then 
              
              ! I changed this so that it can take into account repeated knots in the middle of the knot vectors
              NumberOfElements_Xi = NURBS%NXiKnotEntries(IPatch_Temporary) - (2*NURBS%NXiKnotOrder(IPatch_Temporary) + 1)!NXiKnotEntriesUNIQUE(IPatch_Temporary)-1
              NumberOfElements_Eta = NURBS%NEtaKnotEntries(IPatch_Temporary) - (2*NURBS%NEtaKnotOrder(IPatch_Temporary) + 1)!NEtaKnotEntriesUNIQUE(IPatch_Temporary)-1
        

            !  allocate(KnotBezierMeshConnectivity( &
            !NumberOfElements_Xi*&
            !NumberOfElements_Eta, 4, &
            !Counters%NPatches) )
              
              allocate(KnotBezierMeshConnectivity( &
            NURBS%Maximum_nel_NURBS, 4, &
            Counters%NPatches) )

              
          else if (NDIM == 3) then 
              
          !allocate(KnotBezierMeshConnectivity( &
          !  (Maximum_NXiKnotEntriesUNIQUE-1)*&
          !  (Maximum_NEtaKnotEntriesUNIQUE-1)*(Maximum_NZetaKnotEntriesUNIQUE-1) , 8, &
          !  Counters%NPatches) )
          
          allocate(KnotBezierMeshConnectivity( &
            NURBS%Maximum_nel_NURBS , 8, &
            Counters%NPatches) )
          
          
          end if 
          
          KnotBezierMeshConnectivity = 0
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          if (NDIM == 2) then
              
              ! For 2D, we need to generalize for multiple patches
              call FormBezierMeshKnotConnectivity(KnotBezierMeshConnectivity, IPatch)!(:,:,IPatch))
              
          else if (NDIM == 3) then 
              
              
              call FormBezierMeshKnotConnectivity_3D(NURBS%NXiKnotEntriesUNIQUE(IPatch), &
                  NURBS%NEtaKnotEntriesUNIQUE(IPatch), &
                  NURBS%NZetaKnotEntriesUNIQUE(IPatch), &
                   IPatch)  
                  !KnotBezierMeshConnectivity(:,:,IPatch), &
                   
          
          end if 
          end do ! loop over patches
          
          
          !do IPatch = 1, Counters%NPatches ! loop over patches
          !if (NDIM == 2) then 
          !    
          !call Generate_STARTELEMCON_ElementConnectivities(NURBS%NXiKnotEntriesUNIQUE(IPatch), NURBS%NEtaKnotEntriesUNIQUE(IPatch), &
          !                                                      NURBS%NXiKnotOrder(IPatch), NURBS%NEtaKnotOrder(IPatch), & 
          !                                                      NURBS%NControlPoints(IPatch), IPatch)!NNodes)
          !
          !elseif (NDIM == 3) then 
          !    
          !call Generate_STARTELEMCON_ElementConnectivities_3D(NURBS%NXiKnotEntriesUNIQUE(IPatch), NURBS%NEtaKnotEntriesUNIQUE(IPatch), NURBS%NZetaKnotEntriesUNIQUE(IPatch), &
          !                                                      NURBS%NXiKnotOrder(IPatch), NURBS%NEtaKnotOrder(IPatch), NURBS%NZetaKnotOrder(IPatch), & 
          !                                                      NURBS%NControlPoints(IPatch))!NNodes) 
          !end if 
          !
          !end do ! loop over patches
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! this is my attempt to generalize this global multipatch collection
          ! First, we need to know the count of the nodes to allocate the global multipatch connectivity
          ! Step 1: the counter should take into consideration of all the control points in patch 1
          !NumberOfUniqueNodes = NumberOfNodesInPatch1 ! treat them all as one patch 
          !IPatch = 1
          !! loop over nodes of patch 1 
          !do ii = 1, NumerOfNodesInPatch1
          !
          !    ! loop over patches
          !    do IPatch_Other = 1, Counters%NPatches
          !        
          !        
          !        ! if they are equal then cycle
          !        if (IPatch == IPatch_Other) then 
          !                ! if equal, cycle to the next patch
          !            cycle ! move onto the other value of IPatch_Other
          !        else
          !            !... move on to compare them    
          !        end if
          !
          !        
          !            
          !        do jj = 1, NumberOfNodesInPatch2
          !
          !            
          !            NodalCoor
          !            
          !            
          !        
          !    
          !        end do 
          !    
          !
          !
          !end do 
          
          !!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!! we need to allocate ActiveNode here after knowing how many unique control points there is
          !!!! we need to do the following: 
          !!!! - loop over the nodes of each patch and see if any of it is in common
          !!!! - we need to count 
          !!!! locally, each node as a number. globally, we still do not have that numbering
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! looping to find nodal sets for multipatch analysis 
          ! - loop over the first nodes and count them 
          ! - move into patch 2 check if fits in any from patch 1, if not, then count it as a new set 
          ! - move into patch 3 and check with patch 1 and patch 2... Does it fit in any there? if not, count it as a new set
          ! - move into patch 4 and check with patch 1, 2, and 3, does it fit in any there?? 
          
          ! so we loop 2 times 
          ! 1) loop accorss patches 2 to 4--> nn 
          ! 2) for each patch, loop accross all previous patches (i.e., 1 for patch 2, 1 and 2 for patch 3)
          
          !logical :: IsItUnique
          !integer(INTEGER_TYPE) :: IPatch_Loop_1
          !integer(INTEGER_TYPE) :: IPatch_Loop_2
          !integer(INTEGER_TYPE) :: INode_Loop_1
          !integer(INTEGER_TYPE) :: INode_Loop_2
          
          !do IPatch = 1, Counters%NPatches
              ! hardcoded
          
          
          call FindGlobalControlPointsInMultipatchSystem(NumberOfGlobalControlPointsUniqueMultipatch)
          
          
          
          
          call MatchGlobalToLocalControlPointsInMultipatchSystem(NumberOfGlobalControlPointsUniqueMultipatch)
          
          
              
          ! I moved this here because it needs to populate the multipatch data
          call InitialiseEdgeNodeTyings()
          
      
        end subroutine InitialiseMeshData
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine MatchGlobalToLocalControlPointsInMultipatchSystem(NumberOfGlobalControlPointsUniqueMultipatch)
              
              implicit none 
              
              integer(INTEGER_TYPE), intent(in) :: NumberOfGlobalControlPointsUniqueMultipatch
              integer(INTEGER_TYPE) :: ControlPointIndex
              
              logical :: IsItWritten
              logical :: logic_test
              logical :: repetition_test
              
              integer(INTEGER_TYPE) :: IPatch, IPatch_Other, INode, INode_Other, stat, IError
              
              integer(INTEGER_TYPE) :: ii
              
              integer(INTEGER_TYPE) :: GlobalMultipatchNodalIndex
              
              real(REAL_TYPE), dimension(NVECTOR) :: NodeCoordinatesSelected
              
              
              
              
          ! Allocating global nodal coordinates to the number of unique control points
          !allocate(GlobalMultipatchNodalCoordinates(NumberOfGlobalControlPointsUniqueMultipatch, NVECTOR),  stat = IError)
          allocate(GlobalMultipatchNodalCoordinates(NumberOfGlobalControlPointsUniqueMultipatch, NVECTOR),  stat = IError)
          GlobalMultipatchNodalCoordinates = 0
          
                 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          IPatch = 1 
          ! Arranging global nodal coordinates in the GlobalMultipatchNodalCoordinates that we just allocated     
          GlobalMultipatchNodalCoordinates(1:Counters%NodTot(IPatch), 1:NVECTOR) = NodalCoordinates(:, 1:NVECTOR, IPatch)
          
          ! Since we stored everything from the first patch
          GlobalMultipatchNodalIndex = Counters%NodTot(IPatch)+1 ! starting from the second patch
          
          
          ! this is for the first patch
          do ControlPointIndex = 1, Counters%NodTot(1)
              Multipatch_Connecting_Local_To_Global_ControlPoints(ControlPointIndex,1) = ControlPointIndex
              
          end do 
          
          !Multipatch_Connecting_Local_To_Global_ControlPointsIndex = 
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          !counter_Multipatch_Connecting_Local_Global = Counters%NodTot(1) + 1
              
          logic_test = .true.
          repetition_test = .false.
          ! start from the second patch
          ! finding the number of nodes that are in common, and the number of nodes that are unique  
          do IPatch = 2, Counters%NPatches !loop over patches first time 
              do INode = 1, Counters%NodTot(IPatch) !loop over NodTot in first patch
                  
                  
                  IsItWritten = .false. 
                  ! we pick up trial nodal coordinates
                  ! nodal coordinates correspond to IPatch it correspnding dof NodTot_IPatch
                  NodeCoordinatesSelected = NodalCoordinates(INode,:,IPatch)  !IPatch 1
                  
                  do IPatch_Other = 1, (IPatch-1) !loop over patches second time
                      
                      ! we do not want to loop over the same patch
                      ! if the first patch we are looping against is equal to the second patch
                      !if (IPatch == IPatch_Other) then 
                      !    ! if equal, cycle to the next patch
                      !    cycle
                      !else 
                      !    !... move on to compare them
                      !end if
                      
                      logic_test=.true.
                      ! loop over the NodTot of the second patch 
                      do INode_Other = 1, Counters%NodTot(IPatch_Other) !loop over control points in other patch
                          
                          if (NDIM == 2) then 
                              
                              
                              ! if Trial_NodalCoordinates equals any NodalCoordinates from 
                          if ( ( NodeCoordinatesSelected(1) .eq. NodalCoordinates(INode_Other,1,IPatch_Other) ) .and. &
                              ( NodeCoordinatesSelected(2) .eq. NodalCoordinates(INode_Other,2,IPatch_Other) ) ) then ! if equal !NodTot_IPatch_Previous_SamePatch
                              ! if you enter here, it means that this node is repeated and we have already encountered it and
                              ! we do not need to renumber this
                              
                              !Counters%CommonNodes = Counters%CommonNodes + 1
                              ! record the "unique" node coordinate
                              
                              ! record the index of the old repeated control point
                              ! Multipatch_Connecting_Local_Global connects the local control point numbering of a patch
                              !     to the global control point numbering within ALL patches
                              
                              !if (NumberOfRepetition(INode_Other, IPatch, IPatch_Other)>1) then 
                                  
                              !end if 
                              
                              !Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch) = INode_Other
                              Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch) = Multipatch_Connecting_Local_To_Global_ControlPoints(INode_Other,IPatch_Other)
                              IsItWritten = .true.
                              !RepetitionRealised = RepetitionRealised + 1
                              !
                              !
                              !if (RepetitionRealised>1) then 
                              !    do ii = 2, RepetitionRealised
                              !    cycle 
                              !    end do 
                              !end if 
                              
                              
                              do ii = 1, (INode-1) ! loop over patches 
                                 
                                  if (Multipatch_Connecting_Local_To_Global_ControlPoints(ii,IPatch) == INode_Other ) then 
                                      
                                      repetition_test = .true. 
                                      
                                      exit
                                      
                                  else 
                                      
                                      repetition_test = .false.
                                      !logic_test = .false.
                                  end if 
                                  
                                  
                              end do 
                              
                              
                              
                              
                              !cycle
                              if (repetition_test == .true.) then 
                                  repetition_test = .false. 
                                  cycle 
                              else 
                              logic_test = .false. 
                              exit
                              end if 
                              
                              !Multipatch_Connecting_Global_To_Local_ControlPoints(INode_Other,1) = IPatch
                              !Multipatch_Connecting_Global_To_Local_ControlPoints(INode_Other,2) = INode
                              
                              
                              
                              end if 
                              
                              
                              
                              
                          else if (NDIM == 3) then 
                              
                          
                          
                          ! if Trial_NodalCoordinates equals any NodalCoordinates from 
                          if ( ( NodeCoordinatesSelected(1) .eq. NodalCoordinates(INode_Other,1,IPatch_Other) ) .and. &
                              ( NodeCoordinatesSelected(2) .eq. NodalCoordinates(INode_Other,2,IPatch_Other) ) .and. &
                              ( NodeCoordinatesSelected(3) .eq. NodalCoordinates(INode_Other,3,IPatch_Other) )) then ! if equal !NodTot_IPatch_Previous_SamePatch
                              ! if you enter here, it means that this node is repeated and we have already encountered it and
                              ! we do not need to renumber this
                              
                              !Counters%CommonNodes = Counters%CommonNodes + 1
                              ! record the "unique" node coordinate
                              
                              ! record the index of the old repeated control point
                              ! Multipatch_Connecting_Local_Global connects the local control point numbering of a patch
                              !     to the global control point numbering within ALL patches
                              
                              !if (NumberOfRepetition(INode_Other, IPatch, IPatch_Other)>1) then 
                                  
                              !end if 
                              
                              !Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch) = INode_Other
                              Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch) = Multipatch_Connecting_Local_To_Global_ControlPoints(INode_Other,IPatch_Other)
                              IsItWritten = .true.
                              !RepetitionRealised = RepetitionRealised + 1
                              !
                              !
                              !if (RepetitionRealised>1) then 
                              !    do ii = 2, RepetitionRealised
                              !    cycle 
                              !    end do 
                              !end if 
                              
                              
                              do ii = 1, (INode-1) ! loop over patches 
                                 
                                  if (Multipatch_Connecting_Local_To_Global_ControlPoints(ii,IPatch) == INode_Other ) then 
                                      
                                      repetition_test = .true. 
                                      
                                      exit
                                      
                                  else 
                                      
                                      repetition_test = .false.
                                      !logic_test = .false.
                                  end if 
                                  
                                  
                              end do 
                              
                              
                              
                              
                              !cycle
                              if (repetition_test == .true.) then 
                                  repetition_test = .false. 
                                  cycle 
                              else 
                              logic_test = .false. 
                              exit
                              end if 
                              
                              
                              end if 
                          !else ! otherwise not equal 
                          end if 
                          
                      end do 
                      
                          !end if ! if equal                     
                      !end do ! loop over the NodTot of the second patch 
                  end do !loop over patches second time 
                  
                  
                  
                  
                  if ((logic_test==.true.).and.(IsItWritten==.false.)) then
                      
                      GlobalMultipatchNodalCoordinates(GlobalMultipatchNodalIndex,1:NVECTOR) = &
                                  NodeCoordinatesSelected
                      
                      Multipatch_Connecting_Local_To_Global_ControlPoints(INode, IPatch) = GlobalMultipatchNodalIndex
                      
                      GlobalMultipatchNodalIndex = GlobalMultipatchNodalIndex + 1
                      
                      end if 
                  
                  
              end do !loop over dofs in a patch
          end do !loop over patches

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          end subroutine MatchGlobalToLocalControlPointsInMultipatchSystem
        
        
        
        
        
        
        
        
        subroutine FindGlobalControlPointsInMultipatchSystem(NumberOfGlobalControlPointsUniqueMultipatch) 
          !----------------------------------------------------------------------------------------------
          ! Purpose: to find the total number of control points in the multipatch system. 
          ! This can then be used to allocate variables subsequently. 
          !----------------------------------------------------------------------------------------------
          
          implicit none 
          
          integer(INTEGER_TYPE), intent(out) :: NumberOfGlobalControlPointsUniqueMultipatch
          real(REAL_TYPE), dimension(NVECTOR) :: NodeCoordinatesSelected
          
          integer(INTEGER_TYPE) :: IPatch_Loop_1
          integer(INTEGER_TYPE) :: INode_Loop_1
          integer(INTEGER_TYPE) :: IPatch_Loop_2
          integer(INTEGER_TYPE) :: INode_Loop_2
          
          logical :: IsItUnique
          
          ! select the first patch as a base patch 
          NumberOfGlobalControlPointsUniqueMultipatch=Counters%NodTot(1) !(IPatch) 
          !end do 
          
          ! Take the first patch 
          do IPatch_Loop_1 = 2, Counters%NPatches  ! IPatch_Loop_1     
              ! loop accross the control points of the the Patch_Loop_1
              do INode_Loop_1 = 1, Counters%NodTot(IPatch_Loop_1) ! INode_Loop_2
              
                  ! find the control point coordinates 
                  NodeCoordinatesSelected = NodalCoordinates(INode_Loop_1,:,IPatch_Loop_1)
              
                  ! Loop accross all the other patches to detect matches
                  do IPatch_Loop_2 = 1,(IPatch_Loop_1-1) ! IPatch_Loop_2
                      ! loop accross the control points of the the Patch_Loop_2
                      do INode_Loop_2 = 1, Counters%NodTot(IPatch_Loop_2) ! INode_Loop_2
                          
                          ! Two dimensional
                          if (NDIM == 2) then 
                              
                          if ( ( NodeCoordinatesSelected(1) == NodalCoordinates(INode_Loop_2,1,IPatch_Loop_2) ) .and. &
                              ( NodeCoordinatesSelected(2) == NodalCoordinates(INode_Loop_2,2,IPatch_Loop_2) ) ) then 
                              
                              IsItUnique = .false.
                              exit
                              else
                                  cycle
                              end if 
                              
                          ! Three dimensional     
                          else if (NDIM == 3) then 
                              
                          if ( ( NodeCoordinatesSelected(1) == NodalCoordinates(INode_Loop_2,1,IPatch_Loop_2) ) .and. &
                              ( NodeCoordinatesSelected(2) == NodalCoordinates(INode_Loop_2,2,IPatch_Loop_2) ) .and. &
                              ( NodeCoordinatesSelected(3) == NodalCoordinates(INode_Loop_2,3,IPatch_Loop_2) ) ) then 
                              
                              
                              
                              IsItUnique = .false.
                              exit
                              
                              else
                                  
                                  cycle
                                  
                              end if 
                              
                          end if 
                          
                      end do ! INode_Loop_1
                      
              
                  end do ! IPatch_Loop_2
              
                  if (IsItUnique == .true.) then 
                      NumberOfGlobalControlPointsUniqueMultipatch=NumberOfGlobalControlPointsUniqueMultipatch+1    
                      cycle
                      !exit
                      end if 
                  
                      IsItUnique = .true.
                  
              end do ! INode_Loop_2 
              
          end do ! IPatch_Loop_1 
    
          
        end subroutine FindGlobalControlPointsInMultipatchSystem
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine MatchControlPointsNodes(ILoadSystem, IPatch) 
        !-------------------------------------------------------------
        ! Purpose: This subroutine is used to match the control points in the 
        ! traction to the actual material. The main variable that is populated is 
        ! Traction_Connecting_Local_To_Global_ControlPoints which takes IControlPoint_Traction,
        ! traction control index, as an input and returns IControlPoint_ActualMaterial, the
        ! actual material control point index.  
        ! 
        ! Note that this subroutine needs multipatch generalization. 
        ! 
              implicit none 
              
              ! local variables 
              ! variable used to detect whether we found the node
              logical :: Node_success 
              ! array to prevent assigning multiple traction control points 
              ! to the same actual material control point 
              logical, allocatable, dimension(:) :: Traction_Matching_Status 
              ! index local variables 
              integer(INTEGER_TYPE) :: IPatch_Temporary, IControlPoint_Traction, IControlPoint_ActualMaterial
              integer(INTEGER_TYPE), intent(in) :: IPatch
              
              ! Load system generalization 
              integer(INTEGER_TYPE) :: ILoadSystem_Temporary = 1
              integer(INTEGER_TYPE), intent(in) :: ILoadSystem
              
              ! initialize 
              IPatch_Temporary = 1 ! work with patch 1
              
              
              
              allocate(Traction_Matching_Status(Counters%NodTot(IPatch)) ) ! allocate the variable based 
                                                                                    ! on no of control point in a patch
              Traction_Matching_Status = .false. ! initialize as false
              Node_success = .false. ! initialize as false
              
              do IControlPoint_Traction = 1, NURBS%NControlPoints_Traction(ILoadSystem, IPatch) ! loop over traction control points
              
                  
                  Node_success = .false.

                  do IControlPoint_ActualMaterial = 1, Counters%NodTot(IPatch) ! loop over actual material control points

                      
                      
                      if (Node_success) then ! have we detected this traction node and assigned it an actual material node?
                          Node_success = .false. 
                          exit ! if so, then exit and move on to the next control point
                      end if 
                      

                      ! Two dimensional 
                      if (NDIM == 2) then 
                          
                          ! check if we have spotted the matching actual node to the traction node
                          if ( (NodalCoordinates_Traction(IControlPoint_Traction,1, ILoadSystem, IPatch)==NodalCoordinates(IControlPoint_ActualMaterial,1, IPatch)) &
                          .and. (NodalCoordinates_Traction(IControlPoint_Traction,2, ILoadSystem, IPatch)==NodalCoordinates(IControlPoint_ActualMaterial,2, IPatch)) &
                          .and. (Traction_Matching_Status(IControlPoint_ActualMaterial)==.false.)) then 
                          
                              ! record the matching actual node 
                              Traction_Connecting_Local_To_Global_ControlPoints(IControlPoint_Traction, ILoadSystem, IPatch) = IControlPoint_ActualMaterial
                              Node_success = .true. ! set to true so that we can exit the loop for this traction node
                              Traction_Matching_Status(IControlPoint_ActualMaterial) = .true. ! set to true which means that we do not need to 
                                                                                              ! assign more traction nodes to this
                              
                              ControlPoint_Weights_Traction(IControlPoint_Traction, ILoadSystem, IPatch) = ControlPoint_Weights(IControlPoint_ActualMaterial, IPatch)
                              
                          else
                              
                              !... need an error message here
                          
                          end if 
                          
                          
                      ! Three dimensional 
                      elseif (NDIM == 3) then  
                          
                          ! check if we have spotted the matching actual node to the traction node
                          if ((NodalCoordinates_Traction(IControlPoint_Traction,1, ILoadSystem, IPatch)==NodalCoordinates(IControlPoint_ActualMaterial,1, IPatch)) &
                          .and. (NodalCoordinates_Traction(IControlPoint_Traction,2, ILoadSystem, IPatch)==NodalCoordinates(IControlPoint_ActualMaterial,2, IPatch)) &
                          .and. (NodalCoordinates_Traction(IControlPoint_Traction,3, ILoadSystem, IPatch)==NodalCoordinates(IControlPoint_ActualMaterial,3, IPatch)) & 
                          .and. (Traction_Matching_Status(IControlPoint_ActualMaterial)==.false.) ) then 
                          
                              ! record the matching actual node 
                              Traction_Connecting_Local_To_Global_ControlPoints(IControlPoint_Traction, ILoadSystem, IPatch) = IControlPoint_ActualMaterial
                              Node_success = .true. ! set to true so that we can exit the loop for this traction node
                              Traction_Matching_Status(IControlPoint_ActualMaterial) = .true. ! set to true which means that we do not need to 
                                                                                              ! assign more traction nodes to this
                          
                              ControlPoint_Weights_Traction(IControlPoint_Traction, ILoadSystem, IPatch) = ControlPoint_Weights(IControlPoint_ActualMaterial, IPatch)
                              
                          else
                              
                              !... need an error message here
                          
                          end if ! matching if statement
                          
                          
                      end if ! NDIM if statement
                      
                      
                  end do ! loop over actual material control points
                  
                  
              end do ! loop over traction control points
              
              
        
        end subroutine MatchControlPointsNodes
        
        
        
        
        
        
        subroutine MatchControlPointsElements(ILoadSystem, IPatch) ! traction on the material points 
              !---------------------------------------------------------
              !
              !
              !
              !
              !---------------------------------------------------------
              
              integer(INTEGER_TYPE) :: ICount
              integer(INTEGER_TYPE) :: IElement_Traction
              integer(INTEGER_TYPE) :: IElement_ActualMaterial
              
              integer(INTEGER_TYPE) :: IElementCon_Traction
              integer(INTEGER_TYPE) :: IElementCon_ActualMaterial
              logical :: element_success
              
              ! Multipatch generalization 
              integer(INTEGER_TYPE) :: IPatch_Temporary = 1
              integer(INTEGER_TYPE), intent(in) :: IPatch
              
              ! Load system generalization 
              integer(INTEGER_TYPE) :: ILoadSystem_Temporary = 1
              integer(INTEGER_TYPE), intent(in) :: ILoadSystem
              
              ! initialize integers
              IPatch_Temporary = 1
              ICount = 0
              
              if (NDIM ==2) then 
                  
                  
              ! looping accross traction element
              do IElement_Traction = 1, NURBS%nel_NURBS_Traction(ILoadSystem, IPatch)
                  
                  
                  element_success = .false.
                  
                      ! looping accross actual element
                      do IElement_ActualMaterial = 1, Counters%NEl(IPatch)!NURBS%nel_NURBS(IPatch)
                          
                          if (element_success == .true.) then 
                                  exit 
                          end if 
                          
                          do IElementCon_Traction = 1, NURBS%nen_NURBS_Traction(ILoadSystem, IPatch)
                              
                              if (element_success == .true.) then 
                                  exit 
                              end if 
                              
                          
                          do IElementCon_ActualMaterial = 1, NURBS%nen_NURBS(IPatch)                              
                              if ( (NodalCoordinates_Traction(NURBS%ElementConnectivities_Traction( IElementCon_Traction, IElement_Traction, ILoadSystem, IPatch), 1, ILoadSystem, IPatch)==NodalCoordinates(ElementConnectivities( IElementCon_ActualMaterial, IElement_ActualMaterial, IPatch ),1, IPatch ))  .and. &
                                  (NodalCoordinates_Traction(NURBS%ElementConnectivities_Traction( IElementCon_Traction, IElement_Traction, ILoadSystem, IPatch), 2, ILoadSystem, IPatch)==NodalCoordinates(ElementConnectivities( IElementCon_ActualMaterial, IElement_ActualMaterial, IPatch ),2, IPatch ))  ) then 
                
                                    ICount = ICount + 1
                              
                                    if (ICount == NURBS%nen_NURBS_Traction(ILoadSystem, IPatch)) then 
                              
                                        Traction_Connecting_Local_To_Global_Elements(IElement_Traction, ILoadSystem, IPatch) = IElement_ActualMaterial
                                        ICount = 0
                                        element_success = .true. 
                                        exit
                                    end if 
                                    
                                end if 
                          end do    
                          
                          
                          end do 
                          
                          
                          ICount = 0
                  end do 
              end do 
                  
                  
                  
                  
                  
              else if (NDIM ==3) then 
                  
              
              ! looping accross traction element
              do IElement_Traction = 1, NURBS%nel_NURBS_Traction(ILoadSystem, IPatch)
                  
                  
                  element_success = .false.
                  
                      ! looping accross actual element
                      do IElement_ActualMaterial = 1, Counters%NEl(IPatch)!NURBS%nel_NURBS(IPatch)
                          
                          if (element_success == .true.) then 
                                  exit 
                          end if 
                          
                          do IElementCon_Traction = 1, NURBS%nen_NURBS_Traction(ILoadSystem, IPatch)
                              
                              if (element_success == .true.) then 
                                  exit 
                              end if 
                              
                          
                          do IElementCon_ActualMaterial = 1, NURBS%nen_NURBS(IPatch)                              
                              if ( (NodalCoordinates_Traction(NURBS%ElementConnectivities_Traction( IElementCon_Traction, IElement_Traction, ILoadSystem, IPatch), 1, ILoadSystem, IPatch)==NodalCoordinates(ElementConnectivities( IElementCon_ActualMaterial, IElement_ActualMaterial, IPatch ),1, IPatch))  .and. &
                                  (NodalCoordinates_Traction(NURBS%ElementConnectivities_Traction( IElementCon_Traction, IElement_Traction, ILoadSystem, IPatch), 2, ILoadSystem, IPatch)==NodalCoordinates(ElementConnectivities( IElementCon_ActualMaterial, IElement_ActualMaterial, IPatch ),2, IPatch)) .and. &
                                  (NodalCoordinates_Traction(NURBS%ElementConnectivities_Traction( IElementCon_Traction, IElement_Traction, ILoadSystem, IPatch), 3, ILoadSystem, IPatch)==NodalCoordinates(ElementConnectivities( IElementCon_ActualMaterial, IElement_ActualMaterial, IPatch ),3, IPatch)) ) then 
                
                                    ICount = ICount + 1
                              
                                    if (ICount == NURBS%nen_NURBS_Traction(ILoadSystem, IPatch)) then 
                              
                                        Traction_Connecting_Local_To_Global_Elements(IElement_Traction, ILoadSystem, IPatch) = IElement_ActualMaterial
                                        ICount = 0
                                        element_success = .true. 
                                        exit
                                    end if 
                                    
                                end if 
                          end do    
                          
                          
                          end do 
                          
                          ICount = 0
                  end do 
              end do 
                          
              end if 
              
        end subroutine MatchControlPointsElements
        
        
        
        
        
        
        
        
        
        subroutine MatchTractionElementToActualMaterialSide(ILoadSystem, IPatch) 
        !--------------------------------------------------------
        ! The actual material is a hexahedral which has 6 sides. 
        ! We need to identify where our quadrilateral sits on this 
        ! hexahedral. This is particularly important 

        integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES, ELEMENTSIDES) :: VerticesID

        !integer(INTEGER_TYPE) :: IPatch_Temporary
        integer(INTEGER_TYPE) :: ICount
        
        logical, dimension(ELEMENTSIDES, NURBS%nel_NURBS_Traction(1, 1)) :: side_test ! hardcoded 1 but should be generalized based on number of patches
        integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES, NDIM) :: NodalCoordinatesTractionSide
        !integer(INTEGER_TYPE), dimension(ELEMENTSIDES, NDIM) :: NodalCoordinatesGlobalSide

        integer(INTEGER_TYPE), dimension(ELEMENTNODES, NDIM) :: NodalCoordinatesGlobalSide

        
        integer(INTEGER_TYPE) :: GlobalElementForTraction
        integer(INTEGER_TYPE) :: IElement_Traction
        integer(INTEGER_TYPE) :: ISide, INode
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch_Temporary = 1 
        integer(INTEGER_TYPE), intent(in) :: IPatch
        
        ! Load system generalization 
        integer(INTEGER_TYPE) :: ILoadSystem_Temporary = 1
        integer(INTEGER_TYPE), intent(in) :: ILoadSystem
        
        !IPatch_Temporary = 1
              
             if (NDIM == 2) then 
                 
              
              
              if (NURBS%NXiKnotOrder(IPatch_Temporary) == 1) then  
              ! linear 
              VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 2, 4,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 4, 3,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 3, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
              
              else if (NURBS%NXiKnotOrder(IPatch_Temporary) == 2) then 

              ! quadratic 
              VerticesID = reshape( (/ 1, 2, 3,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 3, 4, 5,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 5, 6, 7,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 7, 8, 9/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 3, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.

              
              else if (NURBS%NXiKnotOrder(IPatch_Temporary) == 3) then 

              ! cubic 
              VerticesID = reshape( (/ 1, 2, 3, 4,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 4, 8, 12, 16,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 16, 15, 14, 13,  & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 13, 9, 5, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 4, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
              
              
              else if (NURBS%NXiKnotOrder(IPatch_Temporary) == 4) then 

              
              ! quartic 
              VerticesID = reshape( (/ 1, 2, 3, 4, 5,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 5, 10, 15, 20, 25,  & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 25, 24, 23, 22, 21,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 21, 16, 11, 6, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 5, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.

              end if 
              
             else if (NDIM == 3) then 
                 
                 
              
                 
                 if (NURBS%NXiKnotOrder(IPatch_Temporary) == 1) then  
              ! linear 
              VerticesID = reshape( (/ 1, 2, 6, 5,  & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 1, 5, 7, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 1, 2, 4, 3,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 6, 2, 4, 8, &
                                 5, 6, 8, 7, &
                                 3, 4, 8, 7/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 4, 6 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
              
                                  
                 else if (NURBS%NXiKnotOrder(IPatch_Temporary) == 2) then 
                ! quadratic 
              VerticesID = reshape( (/ 1,2,3, 10,11,12, 19,20,21,  & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 1,10,19, 4,13,22, 7,16,25,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                  1,2,3, 4,5,6, 7,8,9,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 3,12,21, 6,15,24, 9,18,27, &
                                 19,20,21, 22,23,24, 25,26,27, &
                                 7,8,9, 16,17,18, 25,26,27/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 9, 6 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.     
                                  
                     
                     
                 else if (NURBS%NXiKnotOrder(IPatch_Temporary) == 3) then 
                ! cubic
              VerticesID = reshape( (/ 1,2,3,4, 17,18,19,20, 33,34,35,36, 49,50,51,52,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 1,17,33,49, 5,21,37,53, 9,25,41,57, 13,29,45,61,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16,  & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                49,50,51,52, 53,54,55,56, 57,58,59,60, 61,62,63,64, &
                                 52,36,20,4, 56,40,24,8, 60,44,28,12, 64,48,32,16, & 
                                 13,14,15,16,  29,30,31,32, 45,46,47,48, 61,62,63,64/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 16, 6 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.     
                                  
                 
              
              else if (NURBS%NXiKnotOrder(IPatch_Temporary) == 4) then 
                ! quartic
              VerticesID = reshape( (/ 1,2,3,4,5, 26,27,28,29,30, 51,52,53,54,55, 76,77,78,79,80, 101,102,103,104,105,  & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 1,26,51,76,101, 6,31,56,81,106, 11,36,61,86,111, 16,41,66,91,116, 21,46,71,96,121,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20, 21,22,23,24,25,  & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 5,30,55,80,105, 10,35,60,85,110, 15,40,65,90,115, 20,45,70,95,120, 25,50,75,100,125, &
                                 101,102,103,104,105, 106,107,108,109,110, 111,112,113,114,115, 116,117,118,119,120, 121,122,123,124,125, &
                                 21,46,74,96,121, 22,47,72,97,122, 23,48,73,98,123, 24,49,74,99,124, 25,50,75,100,125/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 25, 6 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.     
                     
                     
                                  
                 end if 
                 
                 
             end if 
             
             
             ! matching 
             !call MatchControlPointsNodes()
              
              
              !allocate( NodalCoordinatesTractionSide (ELEMENTBOUNDARYNODES, NDIM) )
              !allocate( NodalCoordinatesGlobalSide (ELEMENTSIDES, NDIM) )
              !allocate( side_test(ELEMENTSIDES, nel_NURBS_Traction) )
              
              ICount = 0
              side_test = .false.
              
              
              if (NDIM == 2) then 
              
              
              ! looping accross traction element
              do IElement_Traction = 1, NURBS%nel_NURBS_Traction(ILoadSystem_Temporary, IPatch_Temporary) ! traction elements
                  
                  NodalCoordinatesTractionSide = NodalCoordinates_Traction( NURBS%ElementConnectivities_Traction(:,IElement_Traction, ILoadSystem_Temporary, IPatch_Temporary), :, ILoadSystem_Temporary, IPatch_Temporary)
                  
                  
                  GlobalElementForTraction = Traction_Connecting_Local_To_Global_Elements(IElement_Traction, ILoadSystem_Temporary, IPatch_Temporary)
                  
                  NodalCoordinatesGlobalSide = NodalCoordinates( ElementConnectivities(:,GlobalElementForTraction, IPatch_Temporary), :, IPatch_Temporary)
                  
                  do ISide = 1, ELEMENTSIDES ! sides
                  
                      do INode = 1, ELEMENTBOUNDARYNODES
                      
                           
                          
                          if ( ( NodalCoordinatesTractionSide(INode,1) == NodalCoordinatesGlobalSide( VerticesID(INode, ISide), 1)) .and. &
                              ( NodalCoordinatesTractionSide(INode,2) == NodalCoordinatesGlobalSide( VerticesID(INode, ISide), 2)) ) then 
                          
                              side_test(ISide , IElement_Traction) = .true.
                              
                              Traction_Connecting_Local_To_Global_Sides(IElement_Traction, ILoadSystem_Temporary, IPatch_Temporary) = ISide
                              
                            end if 
                              
                          
                      end do ! node 
                  
                  end do ! sides
                  !element_success = .false.
                  
                  
              end do ! traction elements
              
              else if (NDIM == 3) then 
              
                  do IElement_Traction = 1, NURBS%nel_NURBS_Traction(ILoadSystem_Temporary, IPatch_Temporary) ! traction elements
                  
                  
                  side_test = .false.
                  
                  NodalCoordinatesTractionSide = NodalCoordinates_Traction( NURBS%ElementConnectivities_Traction(:,IElement_Traction, ILoadSystem_Temporary, IPatch_Temporary), :, ILoadSystem_Temporary, IPatch_Temporary)
                  
                  
                  GlobalElementForTraction = Traction_Connecting_Local_To_Global_Elements(IElement_Traction, ILoadSystem_Temporary, IPatch_Temporary)
                  
                  NodalCoordinatesGlobalSide = NodalCoordinates( ElementConnectivities(:,GlobalElementForTraction, IPatch_Temporary), :, IPatch_Temporary)
                  
                  do ISide = 1, ELEMENTSIDES ! sides
                  
                      do INode = 1, ELEMENTBOUNDARYNODES
                      
                           
                          
                          if ( ( NodalCoordinatesTractionSide(INode,1) == NodalCoordinatesGlobalSide( VerticesID(INode, ISide), 1)) .and. &
                              ( NodalCoordinatesTractionSide(INode,2) == NodalCoordinatesGlobalSide( VerticesID(INode, ISide), 2)) .and. &
                          ( NodalCoordinatesTractionSide(INode,3) == NodalCoordinatesGlobalSide( VerticesID(INode, ISide), 3))) then 
                          
                              side_test(ISide , IElement_Traction) = .true.
                              
                              Traction_Connecting_Local_To_Global_Sides(IElement_Traction, ILoadSystem_Temporary, IPatch_Temporary) = ISide
                              
                              
                              !return
                              
                            end if 
                              
                          if (side_test(ISide , IElement_Traction) == .true.) then 
                              exit 
                              
                          end if 
                          
                          
                      end do ! node 
                      
                      if (side_test(ISide , IElement_Traction) == .true.) then 
                              exit 
                      end if 
                          
                  
                  end do ! sides
                  
                  
              end do ! traction elements
                  
              end if 
              
              end subroutine MatchTractionElementToActualMaterialSide
        
        
        
        
        subroutine FormBezierMeshKnotConnectivity &
            (KnotBezierMeshConnectivity, IPatch)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Description: This subroutine forms the connectivity matrix for forming the Bezier connectivity mesh. 
        !              Note that this works in tandem with FormBezierMeshUsingUniqueKnotVectors to plot the 
        !              control net in the paraview file. 
        !
        !               Need to write more to describe this implementation.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
        !KnotBezierMeshConnectivity
        implicit none 
        
        
        ! local 
        
        
        
        ! input 
        
        
        
        
        ! output
        !real(REAL_TYPE), allocatable, dimension(:,:), intent(out) :: KnotBezierMeshConnectivity
        
        integer(INTEGER_TYPE) :: ii, jj, kk, counter, &
            NumberOfElements_Xi, NumberOfElements_Eta
        integer(INTEGER_TYPE) :: sum_variable
        
        integer(INTEGER_TYPE), allocatable, dimension(:,:,:), &
            intent(inout) :: KnotBezierMeshConnectivity ! number of elements x number of nodes in element = 4 x 4
        
        integer(INTEGER_TYPE), dimension(1,4) :: UnitBezierVector
        integer(INTEGER_TYPE) :: UniqueHorizontalIndex, &
            UniqueHorizontalDummyIndex
        
        ! Multipatch variable 
        integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        integer(INTEGER_TYPE), intent(in) :: IPatch

        
        !KnotBezierMeshConnectivity = 0
        NumberOfElements_Xi = NURBS%NXiKnotEntries(IPatch) - (2*NURBS%NXiKnotOrder(IPatch) + 1)!NXiKnotEntriesUNIQUE(IPatch_Temporary)-1
        NumberOfElements_Eta = NURBS%NEtaKnotEntries(IPatch) - (2*NURBS%NEtaKnotOrder(IPatch) + 1)!NEtaKnotEntriesUNIQUE(IPatch_Temporary)-1
        
        ! 1) For each element start with number of element ee 
        ! 2) use 2ee-1 as a starting element index 
        ! 3) fill the connectivity matrix by +0, +1, +2, +3 --> looping accross the element sides 
        ! 4) 
        ! 5) 
        
        sum_variable = 1
        
        !! Creating the unit vector 
        UnitBezierVector(1, 1) = 1!UniqueHorizontalDummyIndex
        UnitBezierVector(1, 2) = 2!UniqueHorizontalDummyIndex+1
        UnitBezierVector(1, 3) = NumberOfElements_Xi+2! 3 UniqueHorizontalIndex NumberOfElements_Eta
        UnitBezierVector(1, 4) = NumberOfElements_Xi+3! 4 UniqueHorizontalIndex+1
  
        
        KnotBezierMeshConnectivity(1,1:4, IPatch) = UnitBezierVector(1,1:4) ! 4 points form an element 
        
        counter = 1
                    
        

        do ii = 1, NumberOfElements_Eta
            
            do jj = 1, NumberOfElements_Xi

                counter = counter + 1
                if (counter > Counters%NEl(IPatch)) then !NURBS%nel_NURBS(IPatch) 
                    exit 
                end if 
                do kk = 1, 4 ! 4 points in one element
            
                    
                    
                    KnotBezierMeshConnectivity(counter, kk, IPatch) = &
                    UnitBezierVector(1, kk) + sum_variable  
     
                end do  

                if (jj == NumberOfElements_Xi) then 
                    KnotBezierMeshConnectivity(counter,1:4, IPatch) = &
                        KnotBezierMeshConnectivity(counter,1:4, IPatch) + 1
                
                end if
                
                UnitBezierVector(1,1:4) = &
                    KnotBezierMeshConnectivity(counter,1:4, IPatch)

            end do                
        end do 
        
            
        
        end subroutine FormBezierMeshKnotConnectivity
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine FormBezierMeshKnotConnectivity_3D&
            (NXiKnotEntriesUNIQUE, &
            NEtaKnotEntriesUNIQUE, &
            NZetaKnotEntriesUNIQUE, &
            IPatch) ! KnotBezierMeshConnectivity, &
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Description: This subroutine forms the connectivity matrix for forming the Bezier connectivity mesh. 
        !              Note that this works in tandem with FormBezierMeshUsingUniqueKnotVectors to plot the 
        !              control net in the paraview file. 
        !
        !               Need to write more to describe this implementation.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        implicit none 
        
        
        ! local 
        
        
        
        ! input 
        integer(INTEGER_TYPE), intent(in) :: NXiKnotEntriesUNIQUE
        integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntriesUNIQUE
        integer(INTEGER_TYPE), intent(in) :: NZetaKnotEntriesUNIQUE
        
        integer(INTEGER_TYPE), intent(in) :: IPatch
        
        integer(INTEGER_TYPE) :: ii, jj, kk, ll, counter, &
        NumberOfElements_Xi, NumberOfElements_Eta, NumberOfElements_Zeta
        integer(INTEGER_TYPE) :: sum_variable
        
        !integer(INTEGER_TYPE), dimension(:,:), &
        !    intent(inout) :: KnotBezierMeshConnectivity ! number of elements x number of nodes in element = 4 x 4
        
        integer(INTEGER_TYPE), dimension(1,8) :: UnitBezierVector
        integer(INTEGER_TYPE) :: UniqueHorizontalIndex, &
            UniqueHorizontalDummyIndex
        
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        !integer(INTEGER_TYPE), intent(in) :: IPatch

        
        
        !KnotBezierMeshConnectivity = 0
        NumberOfElements_Xi = (NURBS%nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NURBS%NXiKnotOrder(IPatch))!NXiKnotEntriesUNIQUE-1
        NumberOfElements_Eta = (NURBS%mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NURBS%NEtaKnotOrder(IPatch))!NEtaKnotEntriesUNIQUE-1
        NumberOfElements_Zeta = (NURBS%oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)-NURBS%NZetaKnotOrder(IPatch))!NZetaKnotEntriesUNIQUE-1

        
        
        KnotBezierMeshConnectivity(1, 1, IPatch) = 1!UniqueHorizontalDummyIndex
        KnotBezierMeshConnectivity(1, 2, IPatch) = 2!UniqueHorizontalDummyIndex+1
        KnotBezierMeshConnectivity(1, 3, IPatch) = KnotBezierMeshConnectivity(1, 1, IPatch) + NumberOfElements_Xi + 1!3!NumberOfElements_Xi+2! 3 UniqueHorizontalIndex NumberOfElements_Eta
        KnotBezierMeshConnectivity(1, 4, IPatch) = KnotBezierMeshConnectivity(1, 2, IPatch) + NumberOfElements_Xi + 1!4!NumberOfElements_Xi+3! 4 UniqueHorizontalIndex+1
        KnotBezierMeshConnectivity(1, 5, IPatch) = ((NumberOfElements_Xi+1)*(NumberOfElements_Eta+1)) + KnotBezierMeshConnectivity(1, 1, IPatch)
        KnotBezierMeshConnectivity(1, 6, IPatch) = ((NumberOfElements_Xi+1)*(NumberOfElements_Eta+1)) + KnotBezierMeshConnectivity(1, 2, IPatch)
        KnotBezierMeshConnectivity(1, 7, IPatch) = ((NumberOfElements_Xi+1)*(NumberOfElements_Eta+1)) + KnotBezierMeshConnectivity(1, 3, IPatch)
        KnotBezierMeshConnectivity(1, 8, IPatch) = ((NumberOfElements_Xi+1)*(NumberOfElements_Eta+1)) + KnotBezierMeshConnectivity(1, 4, IPatch)
        
        
        do ii = 2, NumberOfElements_Xi
            do jj = 1, 8
                
                KnotBezierMeshConnectivity(ii,jj, IPatch) = KnotBezierMeshConnectivity(ii-1,jj, IPatch)+1
            
            end do
        end do 
        
        do ii = 2, NumberOfElements_Eta
                
                KnotBezierMeshConnectivity(((ii-1)*NumberOfElements_Xi)+1:ii*NumberOfElements_Xi, :, IPatch) = KnotBezierMeshConnectivity(((ii-2)*NumberOfElements_Xi)+1:(ii-1)*NumberOfElements_Xi, :, IPatch) + NumberOfElements_Xi +1
               
        end do
        
        
        do ii = 2, NumberOfElements_Zeta
                
                KnotBezierMeshConnectivity(((ii-1)*NumberOfElements_Xi*NumberOfElements_Eta)+1:ii*NumberOfElements_Xi*NumberOfElements_Eta, :, IPatch) = KnotBezierMeshConnectivity(((ii-2)*NumberOfElements_Xi*NumberOfElements_Eta)+1:(ii-1)*NumberOfElements_Xi*NumberOfElements_Eta, :, IPatch) + ((NumberOfElements_Xi+1) * (NumberOfElements_Eta+1))
               
        end do
        
        
        end subroutine 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine FormBezierMeshUsingUniqueKnotVectors(XiKnotEntriesUNIQUE, NXiKnotEntriesUNIQUE, & 
                                                        EtaKnotEntriesUNIQUE, NEtaKnotEntriesUNIQUE, &
                                                        KnotBezierMesh) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Description: It is necessary to know the boundaries of the knot spans as they do not necessarily 
        !             conform with the control point locations. The KnotBezierMesh tells you the coordinates of
        !             the knot spans in 2D. 
        !             -This subroutine needs to be generalized to 3D. 
        !             -First, you find the unique xi knots 
        !             -Second, you find the unique eta knots 
        !             -Third, you find combos between the xi and eta knots
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none 
        ! Abdelrahman Alsardi: this subroutine needs multipatch generalization 
        ! locals 
        integer(INTEGER_TYPE) :: stat 
        integer(INTEGER_TYPE) :: IError
        integer(INTEGER_TYPE) :: ii, jj, kk
        integer(INTEGER_TYPE) :: NoOfItems
        integer(INTEGER_TYPE) :: counter
        integer(INTEGER_TYPE) :: Xi_loop
        integer(INTEGER_TYPE) :: Eta_loop
        
        ! inputs 
        integer(INTEGER_TYPE), intent(in) :: NXiKnotEntriesUNIQUE
        integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntriesUNIQUE
        
        real(REAL_TYPE), intent(in), dimension(:) :: XiKnotEntriesUNIQUE       
        real(REAL_TYPE), intent(in), dimension(:) :: EtaKnotEntriesUNIQUE

        ! outputs
        !real(REAL_TYPE), allocatable, intent(out), dimension(:,:) :: KnotBezierMesh         
        real(REAL_TYPE), intent(out), dimension(:,:) :: KnotBezierMesh 

        
        NoOfItems = NXiKnotEntriesUNIQUE * NEtaKnotEntriesUNIQUE
        
        
        ! initialize
        !allocate(KnotBezierMesh(NoOfItems, NDIM), stat = IError) !--> we do not need to allocate this here... it can be outside
        KnotBezierMesh = 0.0
        
        counter = 1
        
        Xi_loop = NoOfItems/NXiKnotEntriesUNIQUE !10/2 = 5 = NEtaKnotEntriesUNIQUE
        Eta_loop = NoOfItems/NEtaKnotEntriesUNIQUE !10/5 = 2 = NXiKnotEntriesUNIQUE
        
        ! uknot
        !do ii = 1, NDIM 
            do  kk = 1, NEtaKnotEntriesUNIQUE
                do jj = 1, NXiKnotEntriesUNIQUE
                   KnotBezierMesh(counter, 1) = XiKnotEntriesUNIQUE(jj)
                   counter = counter + 1
                end do 
            end do    
        !end do 
        
        counter = 1
        
        ! vknot 
        !do ii = 1, NDIM 
            !do  jj = 1, Eta_loop
                do kk = 1, NEtaKnotEntriesUNIQUE
                    do jj = 1, NXiKnotEntriesUNIQUE
                        KnotBezierMesh(counter, 2) = EtaKnotEntriesUNIQUE(kk)
                        counter = counter + 1
                   !KnotBezierMesh(counter, 2) = EtaKnotEntriesUNIQUE(kk)
                   !counter = counter + 1
                    end do 
                end do 
                
            !end do    
        !end do 
        
        ! wknot
        ! z-direction
        
        !KnotBezierMesh
        
                                                        end subroutine FormBezierMeshUsingUniqueKnotVectors
        
        

                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                              
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                 
                                                        
                                                        
                                                        
        subroutine FormBezierMeshUsingUniqueKnotVectors_3D(XiKnotEntriesUNIQUE, NXiKnotEntriesUNIQUE, & 
                                                        EtaKnotEntriesUNIQUE, NEtaKnotEntriesUNIQUE, &
                                                        ZetaKnotEntriesUNIQUE, NZetaKnotEntriesUNIQUE, &
                                                        IPatch) !KnotBezierMesh) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Description: It is necessary to know the boundaries of the knot spans as they do not necessarily 
        !             conform with the control point locations. The KnotBezierMesh tells you the coordinates of
        !             the knot spans in 2D. 
        !             -This subroutine needs to be generalized to 3D. 
        !             -First, you find the unique xi knots 
        !             -Second, you find the unique eta knots 
        !             -Third, you find combos between the xi and eta knots
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none 
        
        ! locals 
        integer(INTEGER_TYPE) :: stat 
        integer(INTEGER_TYPE) :: IError
        integer(INTEGER_TYPE) :: ii, jj, kk
        integer(INTEGER_TYPE) :: NoOfItems
        integer(INTEGER_TYPE) :: counter
        integer(INTEGER_TYPE) :: Xi_loop
        integer(INTEGER_TYPE) :: Eta_loop
        integer(INTEGER_TYPE) :: Zeta_loop
        
        ! inputs 
        integer(INTEGER_TYPE), intent(in) :: NXiKnotEntriesUNIQUE
        integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntriesUNIQUE
        integer(INTEGER_TYPE), intent(in) :: NZetaKnotEntriesUNIQUE
        
        integer(INTEGER_TYPE), intent(in) :: IPatch
        
        real(REAL_TYPE), intent(in), dimension(:) :: XiKnotEntriesUNIQUE       
        real(REAL_TYPE), intent(in), dimension(:) :: EtaKnotEntriesUNIQUE
        real(REAL_TYPE), intent(in), dimension(:) :: ZetaKnotEntriesUNIQUE

        ! outputs
        !real(REAL_TYPE), intent(out), dimension(:,:) :: KnotBezierMesh  !allocatable,
        
        
        !This subroutine worked for 1 by 1 by 1 for 3D NURBS
        
        
        NoOfItems = NXiKnotEntriesUNIQUE * NEtaKnotEntriesUNIQUE * NZetaKnotEntriesUNIQUE
        
        
        !! initialize
        !allocate(KnotBezierMesh(NoOfItems, NDIM), stat = IError)
        !KnotBezierMesh = 0.0
        
        counter = 1
        
        Xi_loop = NoOfItems/NXiKnotEntriesUNIQUE !10/2 = 5 = NEtaKnotEntriesUNIQUE
        Eta_loop = NoOfItems/NEtaKnotEntriesUNIQUE !10/5 = 2 = NXiKnotEntriesUNIQUE
        Zeta_loop = NoOfItems/NZetaKnotEntriesUNIQUE !
        ! uknot
        !do ii = 1, NDIM 
        do ii = 1, NZetaKnotEntriesUNIQUE
            do  kk = 1, NEtaKnotEntriesUNIQUE
                do jj = 1, NXiKnotEntriesUNIQUE
                   KnotBezierMesh(counter, 1, IPatch) = XiKnotEntriesUNIQUE(jj) !* HighestPointOnASurfaceDefiningPhysicalDomain(1)
                   counter = counter + 1
                end do 
            end do    
        end do 
        
        counter = 1
        
        ! vknot 
        !do ii = 1, NDIM 
            !do  jj = 1, Eta_loop
       
        do jj = 1, NZetaKnotEntriesUNIQUE
        do  kk = 1, NEtaKnotEntriesUNIQUE

            do ii = 1, NXiKnotEntriesUNIQUE
            
                
                    
                        
                        KnotBezierMesh(counter, 2, IPatch) = EtaKnotEntriesUNIQUE(kk) !* HighestPointOnASurfaceDefiningPhysicalDomain(2)
                        counter = counter + 1
                   !KnotBezierMesh(counter, 2) = EtaKnotEntriesUNIQUE(kk)
                   !counter = counter + 1
                    end do 
                end do
        end do 
        
        
        counter = 1

        
        
        
                
        do jj = 1, NZetaKnotEntriesUNIQUE
            do ii = 1, NXiKnotEntriesUNIQUE
            do  kk = 1, NEtaKnotEntriesUNIQUE
                
                        KnotBezierMesh(counter, 3, IPatch) = ZetaKnotEntriesUNIQUE(jj) !* HighestPointOnASurfaceDefiningPhysicalDomain(3)
                        counter = counter + 1
                   !KnotBezierMesh(counter, 2) = EtaKnotEntriesUNIQUE(kk)
                   !counter = counter + 1
                    end do 
                end do
        end do 
        
        
        
        
        
                
            !end do    
        !end do 
        
        ! wknot
        ! z-direction
        
        !KnotBezierMesh
        
         end subroutine FormBezierMeshUsingUniqueKnotVectors_3D
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        

        subroutine InitialiseEdgeNodeTyings()
        !**********************************************************************
        !
        !    Function:  Ties edge nodes of an element to its corner nodes.
        !> @note :      This is used for the 'old' 10-noded tetrahedral elements
        !               for the elementtype before v2018.2
        !
        !**********************************************************************
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError, I
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

          allocate(EdgeNodeTyingsHOE(2, Counters%Sum_NodTot), stat = IError)
          EdgeNodeTyingsHOE = 0
       
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over control points 
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
              
              EdgeNodeTyingsHOE(1, GlobalNodeID) = GlobalNodeID
              EdgeNodeTyingsHOE(2, GlobalNodeID) = GlobalNodeID
              
            !EdgeNodeTyingsHOE(1, I) = I
            !EdgeNodeTyingsHOE(2, I) = I
          end do ! loop over control points
          end do ! loop over patches
          
          if (ELEMENTTYPE == TETRAOLD) then ! dimension 10 as only for 10-noded tetrahedral element
            do I = 1, Counters%Sum_NEl !NURBS%Sum_nel_NURBS!Counters%NEl!(IPatch_Temporary)
              call TieHOENode(ElementConnectivities10Node(1:10, I), 1, 5, 2)
              call TieHOENode(ElementConnectivities10Node(1:10, I), 2, 6, 3)
              call TieHOENode(ElementConnectivities10Node(1:10, I), 3, 7, 1)
              call TieHOENode(ElementConnectivities10Node(1:10, I), 1, 8, 4)
              call TieHOENode(ElementConnectivities10Node(1:10, I), 2, 9, 4)
              call TieHOENode(ElementConnectivities10Node(1:10, I), 3, 10, 4)
            end do
          end if 
          
        end subroutine InitialiseEdgeNodeTyings
        
        
        subroutine TieHOENode(Connectivities, CornerNode1, EdgeNode, CornerNode2) ! 3D function
        !**********************************************************************
        !
        !    Function:  Ties EdgeNode to CornerNode1 and CornerNode2.
        !>   @note : only used for TETRAOLD elementtype which is 3D olny.
        !
        !     Connectivities : Node IDs of the considered element.
        !     CornerNode1 : ID of a corner node.
        !     EdgeNode : ID of a node between the corner nodes.
        !     CornerNode2 : ID of a corner node.
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), dimension(10), intent(in) :: Connectivities ! dimension 10 as only required for 3D
          integer(INTEGER_TYPE), intent(in) :: CornerNode1, EdgeNode, CornerNode2
          ! Local variables
          integer(INTEGER_TYPE) :: NN

          NN = Connectivities(EdgeNode)
          EdgeNodeTyingsHOE(1, NN) = Connectivities(CornerNode1)
          EdgeNodeTyingsHOE(2, NN) = Connectivities(CornerNode2)

        end subroutine TieHOENode

        
        subroutine DestroyMeshData()
        !**********************************************************************
        !
        !    Function:  Deallocates the arrays used in this module.
        !
        !**********************************************************************
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(NodalCoordinates)) then
            deallocate(NodalCoordinates, stat = IError)
          end if
          
          if (allocated(NodalOriginalCoord)) then 
            deallocate(NodalOriginalCoord, stat = IError)
          end if

          if (allocated(NodalCoordinatesUpd)) then
            deallocate(NodalCoordinatesUpd, stat = IError)
          end if

          if (allocated(ElementConnectivities)) then
            deallocate(ElementConnectivities, stat = IError)
          end if
          
          if (allocated(ElementConnectivities10Node)) then
            deallocate(ElementConnectivities10Node, stat = IError)
          end if

          if (allocated(EdgeNodeTyingsHOE) ) then
            deallocate(EdgeNodeTyingsHOE, stat = IError)
          end if

          if (allocated(LoadOnNodesConnectivitiesSolid)) then
            deallocate(LoadOnNodesConnectivitiesSolid, stat = IError)
          end if

          if (allocated(LoadOnNodesConnectivitiesWater)) then
            deallocate(LoadOnNodesConnectivitiesWater, stat = IError)
          end if

          if (allocated(LoadOnNodesConnectivitiesGas)) then
            deallocate(LoadOnNodesConnectivitiesGas, stat = IError)
          end if
                   
         if (allocated(LoadOnNodesConnectivitiesSolidB)) then
            deallocate(LoadOnNodesConnectivitiesSolidB, stat = IError)
         end if
         
          if (allocated(LoadOnNodesConnectivitiesWaterB)) then
            deallocate(LoadOnNodesConnectivitiesWaterB, stat = IError)
          end if

          if (allocated(LoadOnNodesConnectivitiesGasB)) then
            deallocate(LoadOnNodesConnectivitiesGasB, stat = IError)
          end if

          if (allocated(LoadOnMatPointsConnectivitiesSolid)) then
            deallocate(LoadOnMatPointsConnectivitiesSolid, stat = IError)
          end if

          if (allocated(LoadOnMatPointsConnectivitiesWater)) then
            deallocate(LoadOnMatPointsConnectivitiesWater, stat = IError)
          end if

          if (allocated(LoadOnMatPointsConnectivitiesGas)) then
            deallocate(LoadOnMatPointsConnectivitiesGas, stat = IError)
          end if
          
          if (allocated(LoadOnMatPointsConnectivitiesSolidB)) then
            deallocate(LoadOnMatPointsConnectivitiesSolidB, stat = IError)
          end if

          if (allocated(LoadOnMatPointsConnectivitiesWaterB)) then
            deallocate(LoadOnMatPointsConnectivitiesWaterB, stat = IError)
          end if

          if (allocated(LoadOnMatPointsConnectivitiesGasB)) then
            deallocate(LoadOnMatPointsConnectivitiesGasB, stat = IError)
          end if

          if (allocated(ElementMaterialID)) then
            deallocate(ElementMaterialID, stat = IError)
          end if

          if (allocated(LoadValuesOnNodesSolid)) then
            deallocate(LoadValuesOnNodesSolid, stat = IError)
          end if

          if (allocated(LoadValuesOnNodesWater)) then
            deallocate(LoadValuesOnNodesWater, stat = IError)
          end if

          if (allocated(LoadValuesOnNodesGas)) then
            deallocate(LoadValuesOnNodesGas, stat = IError)
          end if
                  
          if (allocated(LoadValuesOnNodesSolidB)) then
            deallocate(LoadValuesOnNodesSolidB, stat = IError)
          end if
          
          if (allocated(LoadValuesOnNodesWaterB)) then
            deallocate(LoadValuesOnNodesWaterB, stat = IError)
          end if

          if (allocated(LoadValuesOnNodesGasB)) then
            deallocate(LoadValuesOnNodesGasB, stat = IError)
          end if

          if (allocated(LoadValuesOnMatPointsSolid)) then
            deallocate(LoadValuesOnMatPointsSolid, stat = IError)
          end if

          if (allocated(LoadValuesOnMatPointsWater)) then
            deallocate(LoadValuesOnMatPointsWater, stat = IError)
          end if

          if (allocated(LoadValuesOnMatPointsGas)) then
            deallocate(LoadValuesOnMatPointsGas, stat = IError)
          end if
          
          if (allocated(LoadValuesOnMatPointsSolidB)) then
            deallocate(LoadValuesOnMatPointsSolidB, stat = IError)
          end if

          if (allocated(LoadValuesOnMatPointsWaterB)) then
            deallocate(LoadValuesOnMatPointsWaterB, stat = IError)
          end if

          if (allocated(LoadValuesOnMatPointsGasB)) then
            deallocate(LoadValuesOnMatPointsGasB, stat = IError)
          end if

          if (allocated(NodalPrescibedDisp)) then
            deallocate(NodalPrescibedDisp, stat = IError)
          end if

          !if (allocated(NodalPrescibedDisp_VolLockSmooth)) then
          !  deallocate(NodalPrescibedDisp_VolLockSmooth, stat = IError)
          !end if
          
          if (allocated(NodalPrescibedDispWater)) then
            deallocate(NodalPrescibedDispWater, stat = IError)
          end if
          
          if (allocated(NodalPrescibedDispGas)) then
            deallocate(NodalPrescibedDispGas, stat = IError)
          end if
          
          if (allocated(IsActiveElement)) then
            deallocate(IsActiveElement, stat = IError)
          end if
          
          if (allocated(NodalDensity)) then
            deallocate(NodalDensity, stat = IError)
          end if
         
          if (allocated(ActiveNodeElementVolume)) then
            deallocate(ActiveNodeElementVolume, stat = IError)
          end if
           
          if (allocated(ActiveNode)) then
            deallocate(ActiveNode, stat = IError)
          end if
          
          if (allocated(PBoundary)) then
            deallocate(PBoundary, stat = IError)
          end if

          if (allocated(PBoundaryQuasiStatic)) then
            deallocate(PBoundaryQuasiStatic, stat = IError)
          end if
          
           if (allocated(PBoundaryWater)) then
            deallocate(PBoundaryWater, stat = IError)
          end if

           if (allocated(PBoundaryGas)) then
            deallocate(PBoundaryGas, stat = IError)
          end if

          if (allocated(ReducedDof)) then
            deallocate(ReducedDof, stat = IError)
          end if
          
          if (allocated(ReducedDofQuasiStatic)) then
            deallocate(ReducedDofQuasiStatic, stat = IError)
          end if
          
          if (allocated(ExtLoadVector)) then
            deallocate(ExtLoadVector, stat = IError)
          end if
          
          if (allocated(ExtLoadVectorB)) then
            deallocate(ExtLoadVectorB, stat = IError)
          end if
          
          if (allocated(ExtLoadVectorWater)) then
            deallocate(ExtLoadVectorWater, stat = IError)
          end if
          
          if (allocated(ExtLoadVectorWaterB)) then
            deallocate(ExtLoadVectorWaterB, stat = IError)
          end if

          if (allocated(ExtLoadVectorGas)) then
            deallocate(ExtLoadVectorGas, stat = IError)
          end if
          
         if (allocated(ExtLoadVectorGasB)) then
            deallocate(ExtLoadVectorGasB, stat = IError)
          end if

          if (allocated(ExtLoad)) then
            deallocate(ExtLoad, stat = IError)
          end if

          if (allocated(ExtLoadTotal)) then
            deallocate(ExtLoadTotal, stat = IError)
          end if

          if (allocated(IntLoad)) then
            deallocate(IntLoad, stat = IError)
          end if

          if (allocated(IntLoadPrevious)) then
            deallocate(IntLoadPrevious, stat = IError)
          end if
          
          if (allocated(GravityLoad)) then
            deallocate(GravityLoad, stat = IError)
          end if

          if (allocated(BulkViscLoad)) then
            deallocate(BulkViscLoad, stat = IError)
          end if

          if (allocated(LumpedMassDry)) then
            deallocate(LumpedMassDry, stat = IError)
          end if
          
          if (allocated(LumpedMassDry_VolLockSmooth)) then
            deallocate(LumpedMassDry_VolLockSmooth, stat = IError)
          end if
          
          !if (allocated(NodalFJacDet_VolLockSmooth)) then
          !  deallocate(NodalFJacDet_VolLockSmooth, stat = IError)
          !end if

          if (allocated(TotalDisplacementSoil)) then
            deallocate(TotalDisplacementSoil, stat = IError)
          end if

          if (allocated(PhaseDisplacementSoil)) then
            deallocate(PhaseDisplacementSoil, stat = IError)
          end if

          if (allocated(IncrementalDisplacementSoil)) then
            deallocate(IncrementalDisplacementSoil, stat = IError)
          end if

          if (allocated(AccumulatedIncDisplacementSoil)) then
            deallocate(AccumulatedIncDisplacementSoil, stat = IError)
          end if

          if (allocated(AccumulatedDisplacementSoil) ) then
            deallocate(AccumulatedDisplacementSoil, stat = IError)
          end if

          if (allocated(AccelerationSoil)) then
            deallocate(AccelerationSoil, stat = IError)
          end if

          if (allocated(TotalVelocitySoil)) then
            deallocate(TotalVelocitySoil, stat = IError)
          end if
          
          if (allocated(TotalVelocitySoilPrevious)) then
            deallocate(TotalVelocitySoilPrevious, stat = IError)
          end if

          if (allocated(RateofMomentum)) then
            deallocate(RateofMomentum, stat = IError)
          end if

          if (allocated(ActiveElement)) then
            deallocate(ActiveElement, stat = IError)
          end if

          if (allocated(ElementDeterminant)) then
            deallocate(ElementDeterminant, stat = IError)
          end if

          if (allocated(TotalDisplacementWater)) then
            deallocate(TotalDisplacementWater, stat = IError)
          end if
                
          if (allocated(PhaseDisplacementWater)) then
            deallocate(PhaseDisplacementWater, stat = IError)
          end if

          if (allocated(TotalVelocityWater)) then
            deallocate(TotalVelocityWater, stat = IError)
          end if
          
          if (allocated(TotalVelocityWaterPrevious)) then
            deallocate(TotalVelocityWaterPrevious, stat = IError)
          end if
          
          if (allocated(IncrementalDisplacementWater)) then
            deallocate(IncrementalDisplacementWater, stat = IError)
          end if
          
          if (allocated(IncrementalDisplacementWaterPrevious)) then
            deallocate(IncrementalDisplacementWaterPrevious, stat = IError)
          end if
          
          if (allocated(AccumulatedDisplacementWater) ) then
            deallocate(AccumulatedDisplacementWater, stat = IError)
          end if

          if (allocated(TotalDisplacementGas)) then
            deallocate(TotalDisplacementGas, stat = IError)
          end if

          if (allocated(TotalVelocityGas)) then
            deallocate(TotalVelocityGas, stat = IError)
          end if

          if (allocated(NonAdvectiveFluxAirInWater)) then
            deallocate(NonAdvectiveFluxAirInWater, stat = IError)
          end if

          if (allocated(NonAdvectiveFluxVapourInGas)) then
            deallocate(NonAdvectiveFluxVapourInGas, stat = IError)
          end if

          if (allocated(AdvectiveFluxDarcyWater)) then
            deallocate(AdvectiveFluxDarcyWater, stat = IError)
          end if

          if (allocated(AdvectiveFluxDarcyAir)) then
            deallocate(AdvectiveFluxDarcyAir, stat = IError)
          end if
          
          if (allocated(ThermalConductionFlux)) then
            deallocate(ThermalConductionFlux, stat = IError)
          end if

          if (allocated(FReaction)) then
            deallocate(FReaction, stat = IError)
          end if

          if (allocated(FReactionWater)) then
            deallocate(FReactionWater, stat = IError)
          end if

          if (allocated(FReactionGas)) then
            deallocate(FReactionGas, stat = IError)
          end if

          if (allocated(RateVolStrain)) then
            deallocate(RateVolStrain, stat = IError)
          end if
          
          if (allocated(IsReactionNodeSurface)) then
            deallocate(IsReactionNodeSurface, stat = IError)
          end if
          if (allocated(IsReactionNode)) then
            deallocate(IsReactionNode, stat = IError)
          end if
          
          if (allocated(ConsideredElemReaction)) then
            deallocate(ConsideredElemReaction, stat = IError)
          end if
          
          if(allocated(OutputSurfaceName)) then
            deallocate (OutputSurfaceName, stat = IError)
          end if
          
                    if (allocated(IntFlow)) then
            deallocate(IntFlow, stat = IError)
          end if
          
          if (allocated(ExtFlow)) then
            deallocate(ExtFlow, stat = IError)
          end if
                    
          if (allocated(RateOfFlux)) then
            deallocate(RateOfFlux, stat = IError)
          end if
                              
          if (allocated(IncrementalPressure)) then
            deallocate(IncrementalPressure, stat = IError)
          end if
          
          if (allocated(TotalPressure)) then
            deallocate(TotalPressure, stat = IError)
          end if

          if (allocated(SubIncrementalDisplacement)) then
             deallocate(SubIncrementalDisplacement, stat = IError)
          end if
      
          if (allocated(SubIncrementalPressure)) then
             deallocate(SubIncrementalPressure, stat = IError)
          end if
          
        end subroutine DestroyMeshData

        !**********************************************************************
        !    SUBROUTINE: InitialiseDerivedMeshData
        ! 
        !    DESCRIPTION: 
        !>   Initialises mesh data that requires further processing
        !    of mesh data read from the project files.
        !
        !>   @note:  Initialises Counters%N
        !
        !**********************************************************************        
        subroutine InitialiseDerivedMeshData()       
        implicit none

          ! Local variables
          integer(INTEGER_TYPE), dimension(NDOFL*Counters%Sum_NodTot) :: NTreat
          integer(INTEGER_TYPE), dimension(NDOFL*Counters%Sum_NodTot) :: NTreatWater
          integer(INTEGER_TYPE), dimension(NDOFL*Counters%Sum_NodTot) :: NTreatGas
          integer(INTEGER_TYPE) :: IError
          
          NTreat = 1
          NTreatWater = 1
          NTreatGas = 1
          
          call InitialiseReduceDOF(NTreat) ! Initialises Counters%N
          call InitialiseReduceDOF_VolLockSmooth() ! Initialises Counters%N_VolLockSmooth
          
          call InitialiseReduceDofQuasiStatic()
          
          call FillPBoundary(NTreat)
          !call FillPBoundary_VolLockSmooth()
          
          call FillNTreatWater(NTreatWater)
          call FillNPBoundaryWater(NTreatWater)
          !call FillPBoundaryQuasiStatic()

          if (CalParams%ApplyAbsorbingBoundary) then
            PBoundaryWater = PBoundary
          end if
        
          !call FillNTreatGas(NTreatGas)
          !call FillNPBoundaryGas(NTreatGas)
      
          allocate(RateVolStrain(NURBS%Maximum_nel_NURBS), stat = IError)!Counters%NEl
          RateVolStrain = 0.0

          !nel_NURBS(IPatch)
          
          
        end subroutine InitialiseDerivedMeshData

        subroutine DetermineActiveNodes()
        !**********************************************************************
        !
        !    Function:  Determine which nodes belong to an active element
        !
        !**********************************************************************
        
        ! Abdelrahman Alsardi: Stamp of approval for multipatch generlization
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IActiveElement, ElementID, INode, NodeID
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          
          ActiveNode = .false.
          
          !!Duan hardcoded
          !if (Counters%NAEl==39) then
          !ActiveElement(39) = 39
          !end if
          
          do IPatch = 1, Counters%NPatches ! loop over all patches
          do IActiveElement = 1, Counters%NAEl(IPatch) !Counters%NAEl ! loop over all elements
          !do IActiveElement = 1, Counters%NAEl
            ElementID = ActiveElement(IActiveElement, IPatch)!,IPatch_Temporary)
            do INode = 1, ELEMENTNODES!(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)!ELEMENTNODES ! loop on element nodes
              NodeID = ElementConnectivities(INode, ElementID, IPatch)!, IPatch) Multipatch_Connecting_Local_To_Global_ControlPoints(
              ActiveNode(NodeID, IPatch) = .true.
            end do
          end do
          end do 
        
        end subroutine DetermineActiveNodes

        
        subroutine InitialiseNumberOfActiveElements()
        !**********************************************************************
        !
        ! Function:  Calculate the number of active elements
        !
        ! O   NEl : Total number of elements
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IEl, ISet
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch  
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do IEl = 1, Counters%NEl(IPatch)!NURBS%nel_NURBS(IPatch) ! this is local element numbering
            ISet = ElementMaterialID(IEl,IPatch) !IEl = local element number within a patch
            ! Initialise activation status of elements
            if (ISet < 0) then
              IsActiveElement(IEl,IPatch) = .false.
            else
              IsActiveElement(IEl, IPatch) = .true.
              if (.not.IsFollowUpPhase()) then ! Else, this counter is read from file
                
                  !Counters%NAEl = Counters%NAEl + 1 ! for single patch
                
                Counters%NAEl(IPatch) = Counters%NAEl(IPatch) + 1 !NURBS%NAEl_NURBS
              
              end if
            end if
          end do ! loop over elements
          end do ! loop over patches
          
          ! bookmark
          ! calculate the sum of the acitve elements
          !NURBS%sum_nael_NURBS = 0
          Counters%Sum_NAEl = 0
          do IPatch = 1, Counters%NPatches
          !NURBS%sum_nael_NURBS = NURBS%sum_nael_NURBS + Counters%NAEl(IPatch)
          Counters%Sum_NAEl = Counters%Sum_NAEl + Counters%NAEl(IPatch)
          end do 

        end subroutine InitialiseNumberOfActiveElements

        subroutine setElementDeterminant(DShapeValuesArray)
        !**********************************************************************
        !
        !    Function: sets determinant of each element
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none
        integer(INTEGER_TYPE) :: iErr,iEl
        real(REAL_TYPE), external :: GetElementDeterminant
        real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray

        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch

        if (.not.CalParams%AutomaticSkipConvection) RETURN

        if (.not.allocated(ElementDeterminant)) then
          iErr = 0
          !sum_nael_NURBS
          !allocate(ElementDeterminant(Counters%NAEl), stat = iErr) ! Abdelrahman Alsardi: this might need generalization
          allocate(ElementDeterminant(Counters%Sum_NAEl, Counters%NPatches), stat = iErr) ! Abdelrahman Alsardi: this might need generalization !NURBS%sum_nael_NURBS
          call AllocationError(IErr, 'ElementDeterminant', 'setElementDeterminant')
        else
          if (size(ElementDeterminant) /= Counters%NAEl(IPatch)) then !Counters%NAEl
            iErr = 0
            deallocate(ElementDeterminant, stat = iErr)
            call DeAllocationError(IErr, 'ElementDeterminant', 'setElementDeterminant')
            !allocate(ElementDeterminant(Counters%NAEl), stat = iErr)
            allocate(ElementDeterminant(Counters%Sum_NAEl, Counters%NPatches), stat = iErr) !NURBS%sum_nael_NURBS
            call AllocationError(IErr, 'ElementDeterminant', 'setElementDeterminant')
          endif
        endif

        !do iEl=1,Counters%NAEl ! loop active elements
       do IPatch = 1, Counters%NPatches ! loop over all patches
          do iEl = 1, Counters%NAEl(IPatch) !Counters%NAEl ! loop over all elements

            ElementDeterminant(iEl, IPatch) = GetElementDeterminant(iEl, ElementConnectivities(:,:,IPatch), NodalCoordinatesUpd(:,:,IPatch), DShapeValuesArray)
          enddo
          end do

        end subroutine setElementDeterminant

        logical function HasDistortedElement(DShapeValuesArray) result(res)
        !**********************************************************************
        !
        !    Function:  returns true when an element is distorted
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none
        integer(INTEGER_TYPE) :: iEl
        real(REAL_TYPE) :: MinDeterminantRatio, DeterminantRatio
        real(REAL_TYPE), external :: GetElementDeterminant

        real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray

        
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary
        integer(INTEGER_TYPE) :: IPatch
        
        res = .false.
        if (.not.allocated(ElementDeterminant)) then
          call setElementDeterminant(DShapeValuesArray)
          RETURN
        endif

        MinDeterminantRatio = huge(MinDeterminantRatio)
        !do iEl=1,Counters%NAEl
            
        do IPatch = 1, Counters%NPatches ! loop over all patches
          do iEl = 1, Counters%NAEl(IPatch) !Counters%NAEl ! loop over all elements !_NURBS

          if (ElementDeterminant(iEl, IPatch) > 0) then
            DeterminantRatio = GetElementDeterminant(iEl, ElementConnectivities(:,:,IPatch), NodalCoordinatesUpd(:,:,IPatch), DShapeValuesArray) / ElementDeterminant(iEl, IPatch)
            MinDeterminantRatio = min(MinDeterminantRatio, DeterminantRatio)
          endif
          enddo ! loop over all elements
        end do ! loop over all patches

        MinimumDeterminantRatioReached = MinDeterminantRatio
        res = MinDeterminantRatio < CalParams%MinimumDeterminantRatio

        if (res) then
          call WriteInLogFile('Distortion detected at time step: ' // trim(String(CalParams%TimeStep)),  FEEDBACK_LEVEL_ALWAYS)
          call WriteInLogFile('Minimum determinant ratio reached: ' // trim(String(MinDeterminantRatio)), FEEDBACK_LEVEL_ALWAYS)
        endif


        end function HasDistortedElement

        subroutine SetIsDistorted(IsDistortedInitialised, DShapeValuesArray)
        !**********************************************************************
        !
        !    Function:  Sets a global variable to true when an element is distorted
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none
        logical, intent(in), optional :: IsDistortedInitialised
        real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray


        if (.not.CalParams%AutomaticSkipConvection) RETURN
        if (present(IsDistortedInitialised)) then
          IsDistorted = IsDistortedInitialised
        else
          IsDistorted = HasDistortedElement(DShapeValuesArray)
        endif

        end subroutine SetIsDistorted

        !**********************************************************************
        !    SUBROUTINE: FillPBoundary
        ! 
        !    DESCRIPTION: 
        !>   Contains code for preparing boundary conditions (solid)
        !
        !>   @param[inout] NTreat : Fixed degrees of freedom set to zero, else 1
        !
        !**********************************************************************
        subroutine FillPBoundary(NTreat)
        implicit none
        
          integer(INTEGER_TYPE), dimension(:), intent(in) :: NTreat ! dimension is total number of degrees of freedom
          ! Local variables
          real(REAL_TYPE) :: UMax
          integer(INTEGER_TYPE) :: IError, I, J, K, PduFree
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

          allocate(PBoundary(Counters%N), stat = IError)
          PBoundary = 1.0 ! Set to 1
          UMax = 0
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch) !Counters%NodTot ! loop all nodes !NURBS%NControlPoints
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            do J = 1, NDOFL ! loop on element dof's
              UMax = max(UMax, NodalPrescibedDisp(i,j,IPatch))
            end do
          end do ! loop over all nodes 
          end do ! loop over patches

          if (UMax > 9.9d4) then
            PduFree = 9.9d4
          else
            PduFree = 19.9
          end if

          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop all nodes
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDisp(I,J,IPatch) < PduFree) then ! Prescribed displacement detected
                K = ReducedDof(GlobalNodeID) + J!ReducedDof(I) + J
                PBoundary(K) = NodalPrescibedDisp(I,J,IPatch)!NodalPrescibedDisp(I,J,IPatch)
              end if
            end do
          end do ! loop all nodes 
          end do ! loop over patches

          !do I = 1, Counters%N ! loop all dof's
          !  if (NTreat(I) == 0) then
          !    PBoundary(I) = 0.0
          !  end if
          !end do

        end subroutine FillPBoundary
        
        
        !**********************************************************************
        !    SUBROUTINE: FillPBoundary
        ! 
        !    DESCRIPTION: 
        !>   Contains code for preparing boundary conditions (solid)
        !
        !>   @param[inout] NTreat : Fixed degrees of freedom set to zero, else 1
        !
        !**********************************************************************
        !subroutine FillPBoundary_VolLockSmooth()!(NTreat)
        !implicit none
        !
        !  !integer(INTEGER_TYPE), dimension(:), intent(in) :: NTreat ! dimension is total number of degrees of freedom
        !  ! Local variables
        !  real(REAL_TYPE) :: UMax
        !  integer(INTEGER_TYPE) :: IError, I, J, K, PduFree
        !  
        !  ! Multipatch variables 
        !  integer(INTEGER_TYPE) :: IPatch_Temporary
        !  integer(INTEGER_TYPE) :: IPatch
        !  integer(INTEGER_TYPE) :: GlobalNodeID
        !
        !  allocate(PBoundary_VolLockSmooth(Counters%N_VolLockSmooth), stat = IError)
        !  PBoundary_VolLockSmooth = 1.0 ! Set to 1
        !  UMax = 0
        !  do IPatch = 1, Counters%NPatches ! loop over patches
        !  do I = 1, NControlPoints_VolLockSmooth(IPatch)!Counters%NodTot ! loop all nodes
        !      
        !      !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
        !      
        !    do J = 1, NDOFL ! loop on element dof's
        !      UMax = max(UMax, NodalPrescibedDisp_VolLockSmooth(i,j,IPatch))
        !    end do
        !  end do ! loop over all nodes 
        !  end do ! loop over patches
        !
        !  if (UMax > 9.9d4) then
        !    PduFree = 9.9d4
        !  else
        !    PduFree = 19.9
        !  end if
        !
        !  do IPatch = 1, Counters%NPatches ! loop over patches
        !  do I = 1, NControlPoints_VolLockSmooth(IPatch)!Counters%NodTot ! loop all nodes
        !      
        !      !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
        !      
        !    do J = 1, NDOFL ! loop on element dof's
        !      if (NodalPrescibedDisp_VolLockSmooth(I,J,IPatch) < PduFree) then ! Prescribed displacement detected
        !        K = ReducedDof_VolLockSmooth(I) + J!ReducedDof_VolLockSmooth(GlobalNodeID) + J!ReducedDof(I) + J
        !        PBoundary_VolLockSmooth(K) = NodalPrescibedDisp_VolLockSmooth(I,J,IPatch) !NodalPrescibedDisp(I,J,IPatch)
        !      end if
        !    end do
        !  end do ! loop all nodes 
        !  end do ! loop over patches
        !
        !  !do I = 1, Counters%N ! loop all dof's
        !  !  if (NTreat(I) == 0) then
        !  !    PBoundary(I) = 0.0
        !  !  end if
        !  !end do
        !
        !end subroutine FillPBoundary_VolLockSmooth
        
        !**********************************************************************
        !    SUBROUTINE: FillPBoundaryQuasiStatic
        ! 
        !    DESCRIPTION: 
        !>   Contains code for preparing boundary conditions (water)
        !
        !**********************************************************************
        subroutine FillPBoundaryQuasiStatic()
          implicit none

          integer(INTEGER_TYPE) :: I, NodeID, IDof, IDofQuasiStatic, IError
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

          !allocate(PBoundaryQuasiStatic(Counters%NodTot * (NDOFL + 1)), stat = IError)          
          allocate(PBoundaryQuasiStatic(Counters%N), stat = IError) !Counters%NodTot

        
          PBoundaryQuasiStatic = 0
          do IPatch = 1, Counters%NPatches ! loop over patches
          do NodeID = 1, Counters%NodTot(IPatch)!Counters%NodTot
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(NodeID,IPatch)
              
            !if (IsElementCornerNode(NodeID)) then
            do I = 1, NDIM
                IDoF = ReducedDoF(GlobalNodeID) + I!ReducedDoF(NodeID) + I
                IDoFQuasiStatic = (GlobalNodeID-1) * (NDOFL + 1) + I!(NodeID-1) * (NDOFL + 1) + I
                PBoundaryQuasiStatic(IDofQuasiStatic) = PBoundary(IDoF)
            end do
            IDoF = ReducedDoF(GlobalNodeID) + 1!ReducedDoF(NodeID) + 1
            IDoFQuasiStatic = GlobalNodeID * (NDOFL + 1)!NodeID * (NDOFL + 1)
            if (PBoundaryWater(IDoF) /= 0) then
                PBoundaryQuasiStatic(IDofQuasiStatic) = 1 ! 1: natural
            else
                PBoundaryQuasiStatic(IDofQuasiStatic) = 0 ! 0: essential
            end if
            !endif
          end do ! loop over control points 
          end do ! loop over patches

        end subroutine FillPBoundaryQuasiStatic
        
        !**********************************************************************
        !    SUBROUTINE: FillNPBoundaryWater
        ! 
        !    DESCRIPTION: 
        !>   Contains code for preparing boundary conditions (water)
        !
        !>   @param[in] NTreatWater : Fixed degrees of freedom set to zero, else 1
        !
        !**********************************************************************
        subroutine FillNPBoundaryWater(NTreatWater)
          implicit none
        
          integer(INTEGER_TYPE), dimension(:), intent(in) :: NTreatWater ! dimension is total number of degrees of freedom
          ! Local variables
          real(REAL_TYPE) :: UMax
          integer(INTEGER_TYPE) :: IError, I, J, K, PduFree
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

          allocate(PBoundaryWater(Counters%N), stat = IError)
         
          PBoundaryWater = 1.0 ! Set to 1
          UMax = 0
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop all nodes
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
              
            do J = 1, NDOFL ! loop on element dof's
              UMax = max(UMax, NodalPrescibedDispWater(i,j,IPatch))
            end do
          end do ! loop all nodes 
          end do ! loop over patches

          if (UMax > 9.9d4) then
            PduFree = 9.9d4
          else
            PduFree = 19.9
          end if

          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop all nodes
              
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
              
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDispWater(I,J,IPatch) < PduFree) then ! Prescribed displacement detected
                !K = ReducedDof(I) + J
                !PBoundaryWater(K) = NodalPrescibedDispWater(I,J,IPatch)
                K = ReducedDof(GlobalNodeID) + J
                PBoundaryWater(K) = NodalPrescibedDispWater(I,J,IPatch)
              end if
            end do
          end do ! loop all nodes (i.e., control points)
          end do ! loop over patches
 
          do I = 1, Counters%N ! loop all dof's
            if (NTreatWater(I) == 0) then
              PBoundaryWater(I) = 0.0
            end if
          end do
          
        end subroutine FillNPBoundaryWater

        !**********************************************************************
        !    SUBROUTINE: FillNPBoundaryGas
        ! 
        !    DESCRIPTION: 
        !>   Contains code for preparing boundary conditions (water)
        !
        !>   @param[out] NTreatGas : Fixed degrees of freedom set to zero, else 1
        !
        !**********************************************************************
        subroutine FillNPBoundaryGas(NTreatGas)
         implicit none
        
          integer(INTEGER_TYPE), dimension(:), intent(in) :: NTreatGas ! dimension is total number of degrees of freedom
          ! Local variables
          real(REAL_TYPE) :: UMax
          integer(INTEGER_TYPE) :: IError, I, J, K, PduFree
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

          allocate(PBoundaryGas(Counters%N), stat = IError)
         
          PBoundaryGas = 1.0 ! Set to 1
          UMax = 0
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop all nodes
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            do J = 1, NDOFL ! loop on element dof's
              UMax = max(UMax, NodalPrescibedDispGas(i,j,IPatch))
            end do
          end do ! loop all nodes 
          end do ! loop over patches

          if (UMax > 9.9d4) then
            PduFree = 9.9d4
          else
            PduFree = 19.9
          end if

          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop all nodes
              
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
              
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDispGas(I,J,IPatch) < PduFree) then ! Prescribed displacement detected
                !K = ReducedDof(I) + J
                !PBoundaryGas(K) = NodalPrescibedDispGas(I,J,IPatch)
                K = ReducedDof(GlobalNodeID) + J
                PBoundaryGas(K) = NodalPrescibedDispGas(I,J,IPatch)
              end if
            end do
          end do ! loop all nodes
          end do ! loop over patches
 
          do I = 1, Counters%N ! loop all dof's
            if (NTreatGas(I) == 0) then
              PBoundaryGas(I) = 0.0
            end if
          end do
          
        end subroutine FillNPBoundaryGas

        !**********************************************************************
        !    SUBROUTINE: InitialiseReduceDOF
        ! 
        !    DESCRIPTION: 
        !>   To fill the vectors NodAct, ReducedDof
        !
        !>   @note: ReducedDof : Vector containing the offset of dof's per node
        !                        Global DOF-dof of node i = ReducedDof(i)+DOF
        !
        !>   @param[inout] NTreat : Fixed degrees of freedom set to zero, else 1
        !
        !**********************************************************************
        subroutine InitialiseReduceDOF(NTreat)
        
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(inout) :: NTreat
          integer(INTEGER_TYPE) I, J, IError, IDof
          integer(INTEGER_TYPE), dimension(:), allocatable :: NodAct
          real(REAL_TYPE) :: aux
          
          ! Multipatch varibales 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          !allocate(ReducedDof(Counters%NodTot + 1), stat = IError)
          !allocate(NodAct(Counters%NodTot), stat = IError)
          
          !allocate(ReducedDof(NumberOfGlobalControlPointsUniqueMultipatch(IPatch_Temporary) + 1), stat = IError)
          allocate(ReducedDof(NumberOfGlobalControlPointsUniqueMultipatch + 1), stat = IError)

          !allocate(NodAct(NumberOfGlobalControlPointsUniqueMultipatch(IPatch_Temporary)), stat = IError)
          allocate(NodAct(NumberOfGlobalControlPointsUniqueMultipatch), stat = IError)


          NodAct = 1
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot
              
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
              
            aux = 0.0
            do J = 1, NDOFL
              aux = aux + dabs(NodalPrescibedDisp(I, J, IPatch))
            end do
            if (aux == 0) NodAct(I) = 0
          end do ! loop all nodes 
          end do ! loop over patches
          
          ReducedDof =-1
          J = 0
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
              if (ReducedDof(GlobalNodeID) == -1) then !if (ReducedDof(I) == -1) then
              ReducedDof(GlobalNodeID) = J
              J = J + NDOFL
              !else
              !    ReducedDof(GlobalNodeID) = ReducedDof(Multipatch_Connecting_Local_To_Global_ControlPoints(I-1,IPatch))
            end if
          end do ! loop over control points
          end do ! loop over patches
          
          !ReducedDof(NumberOfGlobalControlPointsUniqueMultipatch(IPatch_Temporary) + 1) = J !Counters%NodTot   
          !ReducedDof(GlobalNodeID + 1) = J !Counters%NodTot

          ReducedDof(NumberOfGlobalControlPointsUniqueMultipatch + 1) = NumberOfGlobalControlPointsUniqueMultipatch*NVector
          
          !Counters%N = ReducedDof(NumberOfGlobalControlPointsUniqueMultipatch(IPatch_Temporary) + 1) ! this is the total number of dof's
          Counters%N = ReducedDof(NumberOfGlobalControlPointsUniqueMultipatch + 1)!ReducedDof(GlobalNodeID + 1) ! this is the total number of dof's

          
          
          
          
          NTreat = 1
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            do J = 1, NDOFL
              if (NodalPrescibedDisp(I, J,IPatch) == 0) then
                IDof = ReducedDof(GlobalNodeID) + J!ReducedDof(I) + J
                NTreat(IDof) = 0
              end if
            end do
          end do ! loop over control points
          end do ! loop over patches
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch) !Counters%NodTot
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            if (NodAct(I) == 0) then
              do J = 1, NDOFL
                NTreat(ReducedDof(GlobalNodeID) + J) = 0
              end do
            end if
          end do ! loop over control points 
          end do ! loop over patches
          
          deallocate(NodAct, stat = IError)
      
        end subroutine InitialiseReduceDOF
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !    SUBROUTINE: InitialiseReduceDOF
        ! 
        !    DESCRIPTION: 
        !>   To fill the vectors NodAct, ReducedDof
        !
        !>   @note: ReducedDof : Vector containing the offset of dof's per node
        !                        Global DOF-dof of node i = ReducedDof(i)+DOF
        !
        !>   @param[inout] NTreat : Fixed degrees of freedom set to zero, else 1
        !
        !**********************************************************************
        subroutine InitialiseReduceDOF_VolLockSmooth()!(NTreat)
        
        implicit none
        
          integer(INTEGER_TYPE) I, J, IError, IDof
          real(REAL_TYPE) :: aux
          
          ! Multipatch varibales 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
                    
          allocate(ReducedDof_VolLockSmooth(NURBS%NControlPoints_VolLockSmooth(IPatch_Temporary) + 1), stat = IError)
          
          ReducedDof_VolLockSmooth =-1
          J = 0
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, NURBS%NControlPoints_VolLockSmooth(IPatch)!Counters%NodTot
              
              !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
              if (ReducedDof_VolLockSmooth(I) == -1) then !if (ReducedDof(I) == -1) then
              ReducedDof_VolLockSmooth(I) = J
              J = J + NDOFL
              !else
              !    ReducedDof(GlobalNodeID) = ReducedDof(Multipatch_Connecting_Local_To_Global_ControlPoints(I-1,IPatch))
            end if
          end do ! loop over control points
          end do ! loop over patches
          
          ! Write up the last reduced dof number for the very last control point
          ReducedDof_VolLockSmooth(NURBS%NControlPoints_VolLockSmooth + 1) = NURBS%NControlPoints_VolLockSmooth(IPatch_Temporary)*NVector
          
          !Counters%N = ReducedDof(NumberOfGlobalControlPointsUniqueMultipatch(IPatch_Temporary) + 1) ! this is the total number of dof's
          ! Total number of degrees of freedom for the lower order degree of freedoms
          Counters%N_VolLockSmooth = ReducedDof_VolLockSmooth(NURBS%NControlPoints_VolLockSmooth(IPatch_Temporary) + 1) !ReducedDof(GlobalNodeID + 1) ! this is the total number of dof's
        
          
        
        end subroutine InitialiseReduceDOF_VolLockSmooth
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !    SUBROUTINE: InitialiseReduceDofQuasiStatic
        ! 
        !    DESCRIPTION: 
        !>   To fill the vectors To fill the vectors NodAct, ReducedDofQuasiStatic
        !
        !    @note: ReducedDof : Vector containing the offset of dof's per node
        !                        Global DOF-dof of node i = ReducedDofQuasiStatic(i)+DOF
        !                        Global p-dof of node i = ReducedDofQuasiStatic(i)+DOF+1
        !
        !**********************************************************************
        subroutine InitialiseReduceDofQuasiStatic() 
        implicit none

        integer(INTEGER_TYPE) :: I, J, IError
        
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch
        integer(INTEGER_TYPE) :: GlobalNodeID

        allocate(ReducedDofQuasiStatic(Counters%Sum_NodTot + 1), stat = IError)
        
        ReducedDofQuasiStatic = -1
        J = 0
        do IPatch = 1, Counters%NPatches ! loop over patches
        do I = 1, Counters%NodTot(IPatch)!Counters%NodTot
            
            GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
            
          ReducedDofQuasiStatic(I) = J
          J = J + NDOFL + 1
        end do ! loop over control points 
        end do ! loop over patches
        ReducedDofQuasiStatic(Counters%NodTot + 1) = J

        end subroutine InitialiseReduceDofQuasiStatic

        
        !**********************************************************************
        !    SUBROUTINE: FillNTreatWater
        ! 
        !    DESCRIPTION: 
        !>   To fill the vector NTreatWater
        !
        !>   @param[inout] NTreatWater : Fixed degrees of freedom set to zero, else 1, water phase
        !
        !**********************************************************************        
        subroutine FillNTreatWater(NTreatWater)
          implicit none
        
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: NTreatWater ! dimension is total number of dof's
          integer(INTEGER_TYPE) I, J, IError
          integer(INTEGER_TYPE), dimension(:), allocatable :: NodAct
          real(REAL_TYPE) :: aux
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          allocate(NodAct(Counters%Sum_NodTot), stat = IError)
 
          NodAct = 1
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over control points 
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            aux = 0.0
            do J = 1, NDOFL
              aux = aux + dabs(NodalPrescibedDispWater(I, J,IPatch_Temporary))
            end do
            if (aux == 0) NodAct(I) = 0
          end do ! loop over control points 
          end do ! loop over patches
                   
          NTreatWater = 1
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over control points
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDispWater(I, J,IPatch_Temporary) == 0) then
                NTreatWater(ReducedDof(GlobalNodeID) + J) = 0
              end If
            end do
          end do ! loop over control points 
          end do ! loop over patches
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over control points
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            if (NodAct(I) == 0) then
              do J = 1, NDOFL ! loop on element dof's
                NTreatWater(ReducedDof(GlobalNodeID) + J) = 0
              end do
            end if
          end do ! loop over control points 
          end do ! loop over patches
          
          deallocate(NodAct, stat = IError)
      
        end subroutine FillNTreatWater

        
        !**********************************************************************
        !    SUBROUTINE: FillNTreatGas
        ! 
        !    DESCRIPTION: 
        !>   To fill the vector NTreatWater
        !
        !>   @param[inout] NTreatWater : Fixed degrees of freedom set to zero, else 1, Gas phase
        !
        !**********************************************************************           
        subroutine FillNTreatGas(NTreatGas)
          implicit none

          integer(INTEGER_TYPE), dimension(:), intent(inout) :: NTreatGas
          
          integer(INTEGER_TYPE) I, J, IError
          integer(INTEGER_TYPE), dimension(:), allocatable :: NodAct
          real(REAL_TYPE) :: aux
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          allocate(NodAct(Counters%Sum_NodTot), stat = IError)

          NodAct = 1
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over control points 
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            aux = 0.0
            do J = 1, NDOFL
              aux = aux + dabs(NodalPrescibedDispGas(I, J,IPatch))
            end do
            if (aux == 0) NodAct(I) = 0
          end do  ! loop over control points 
          end do ! loop over patches

          NTreatGas = 1
          do IPatch = 1, Counters%NPatches ! loop over patches 
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over control points
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDispGas(I, J,IPatch) == 0) then
                NTreatGas(ReducedDof(GlobalNodeID) + J) = 0
              end If
            end do
          end do ! loop over control points 
          end do ! loop over patches
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over control points 
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            if (NodAct(I) == 0) then
              do J = 1, NDOFL ! loop on element dof's
                NTreatGas(ReducedDof(GlobalNodeID) + J) = 0
              end do
            end if
          end do ! loop over control points 
          end do ! loop over patches
          
          deallocate(NodAct, stat = IError)
      
        end subroutine FillNTreatGas


        subroutine InitialiseNodalArrays()
        !**********************************************************************
        !
        !  Function : To fill the vectors NodAct, ReducedDof
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: IError

          allocate(ExtLoadVector(Counters%N), stat = IError)
          ExtLoadVector = 0.0
          
          ! hardcode fix
          allocate(RateOfMomentumWNet(Counters%N), stat = IError)
          RateOfMomentumWNet = 0.0
          
          allocate(ExtLoadVectorB(Counters%N), stat = IError)
          ExtLoadVectorB = 0.0
          
          allocate(ExtLoadVectorWater(Counters%N), stat = IError)
          ExtLoadVectorWater = 0.0
          
          allocate(ExtLoadVectorWaterB(Counters%N), stat = IError)
          ExtLoadVectorWaterB = 0.0

          allocate(ExtLoadVectorGas(Counters%N), stat = IError)
          ExtLoadVectorGas = 0.0
          
          allocate(ExtLoadVectorGasB(Counters%N), stat = IError)
          ExtLoadVectorGasB = 0.0
          
          allocate(HydraulicHeadVector(Counters%N), stat = IError)
          HydraulicHeadVector = 0.0

          allocate(ExtLoad(Counters%N, Counters%NEntity), stat = IError)
          ExtLoad = 0.0

          allocate(ExtLoadTotal(Counters%N, Counters%NEntity,Counters%NSoilLoadSystems), stat = IError)
          ExtLoadTotal = 0.0
          
          allocate(HydraulicHeadLoadTotal(Counters%N, Counters%NEntity), stat = IError)
          HydraulicHeadLoadTotal= 0.0

          allocate(IntLoad(Counters%N, Counters%NEntity), stat = IError)
          IntLoad = 0.0
          
          allocate(IntLoadPrevious(Counters%N, Counters%NEntity), stat = IError)
          IntLoadPrevious = 0.0
          
          allocate(GravityLoad(Counters%N, Counters%NEntity), stat = IError)
          GravityLoad = 0.0

          allocate(BulkViscLoad(Counters%N, Counters%NEntity), stat = IError)
          BulkViscLoad = 0.0
          
          allocate(LumpedMassDry(Counters%N, Counters%NEntity), stat = IError)
          LumpedMassDry = 0.0
          
          !allocate(NodalStress_VolLockSmooth(Counters%N, Counters%NEntity), stat = IError) ! nodal stress mapped to the actual mesh
          !NodalStress_VolLockSmooth = 0.0
          
          allocate(LumpedMassDry_VolLockSmooth(Counters%N_VolLockSmooth, Counters%NEntity), stat = IError)
          LumpedMassDry_VolLockSmooth = 0.0
          
          !allocate(WeightedFJacDet_Numerator(Counters%N_VolLockSmooth, Counters%NEntity), stat = IError)
          !WeightedFJacDet_Numerator = 0.0
          
          !allocate(NodalFJacDet_VolLockSmooth(Counters%N_VolLockSmooth, Counters%NEntity), stat = IError)
          !NodalFJacDet_VolLockSmooth = 0.0
          
          if (CalParams%ApplyFBarLowerOrder) then 
 
          allocate(WeightedStrainIncrementTensor_Numerator(NTENSOR, NURBS%NControlPoints_VolLockSmooth(1), Counters%NEntity), stat = IError)
          WeightedStrainIncrementTensor_Numerator = 0.0
          
          allocate(NodalStrainIncrement_VolLockSmooth(NTENSOR, NURBS%NControlPoints_VolLockSmooth(1), Counters%NEntity), stat = IError)
          NodalStrainIncrement_VolLockSmooth = 0.0
          
          
          allocate(WeightedStressIncrementTensor_Numerator(NTENSOR, NURBS%NControlPoints_VolLockSmooth(1), Counters%NEntity), stat = IError)
          WeightedStressIncrementTensor_Numerator = 0.0
          
          allocate(NodalStressIncrement_VolLockSmooth(NTENSOR, NURBS%NControlPoints_VolLockSmooth(1), Counters%NEntity), stat = IError)
          NodalStressIncrement_VolLockSmooth = 0.0
          
          
          ! if 2-phase 
          allocate(WeightedWP_Numerator(NURBS%NControlPoints_VolLockSmooth(1), Counters%NEntity), stat = IError) !NTENSOR, 
          WeightedWP_Numerator = 0.0
          
          allocate(NodalWP_VolLockSmooth(NURBS%NControlPoints_VolLockSmooth(1), Counters%NEntity), stat = IError) !NTENSOR, 
          NodalWP_VolLockSmooth = 0.0
          
          
          
          
          else
              
          allocate(WeightedStrainIncrementTensor_Numerator(NTENSOR, Counters%NodTot(1), Counters%NEntity), stat = IError)
          WeightedStrainIncrementTensor_Numerator = 0.0
          
          allocate(NodalStrainIncrement_VolLockSmooth(NTENSOR, Counters%NodTot(1), Counters%NEntity), stat = IError)
          NodalStrainIncrement_VolLockSmooth = 0.0    
          
          allocate(WeightedStressIncrementTensor_Numerator(NTENSOR, Counters%NodTot(1), Counters%NEntity), stat = IError)
          WeightedStressIncrementTensor_Numerator = 0.0
          
          allocate(NodalStressIncrement_VolLockSmooth(NTENSOR, Counters%NodTot(1), Counters%NEntity), stat = IError)
          NodalStressIncrement_VolLockSmooth = 0.0   
          
          ! if 2-phase 
          allocate(WeightedWP_Numerator(Counters%NodTot(1), Counters%NEntity), stat = IError) !NTENSOR, 
          WeightedWP_Numerator = 0.0
          
          allocate(NodalWP_VolLockSmooth(Counters%NodTot(1), Counters%NEntity), stat = IError) !NTENSOR, 
          NodalWP_VolLockSmooth = 0.0
              
          end if
          
          
          
          
          !allocate(WeightedStressTensor_Numerator(NTENSOR, NControlPoints_VolLockSmooth(1), Counters%NEntity), stat = IError)
          !WeightedStressTensor_Numerator = 0.0
          
          
          !allocate(NodalStress_VolLockSmooth(NTENSOR, NControlPoints_VolLockSmooth(1), Counters%NEntity), stat = IError)
          !NodalStress_VolLockSmooth = 0.0
          !allocate(WeightedNodalVolumetricStressGradient_Numerator(Counters%N_VolLockSmooth, Counters%NEntity), stat = IError)
          !WeightedNodalVolumetricStressGradient_Numerator = 0.0
          
          !allocate(NodalVolumetricStrain_VolLockSmooth(Counters%N_VolLockSmooth, Counters%NEntity), stat = IError)
          !NodalVolumetricStrain_VolLockSmooth = 0.0
          
          !allocate(NodalVolumetricStressGradient_VolLockSmooth(Counters%N_VolLockSmooth, Counters%NEntity), stat = IError)
          !NodalVolumetricStressGradient_VolLockSmooth = 0.0
          
          ! for F-bar technique
          !allocate(NodalStressArray(Counters%N), stat = IError) !--> only volumetric needed so Counters%N is okay
          !NodalStressArray = 0.0 !, Counters%NEntity

          allocate(TotalDisplacementSoil(Counters%N), stat = IError)
          TotalDisplacementSoil = 0.0

          allocate(PhaseDisplacementSoil(Counters%N), stat = IError)
          PhaseDisplacementSoil = 0.0

          allocate(IncrementalDisplacementSoil(Counters%N, Counters%NEntity), stat = IError)
          IncrementalDisplacementSoil = 0.0
          
          !allocate(IncrementalDisplacementSoil_VolLockSmooth(Counters%N_VolLockSmooth, Counters%NEntity), stat = IError)
          !IncrementalDisplacementSoil_VolLockSmooth = 0.0

          allocate(AccumulatedIncDisplacementSoil(Counters%N, Counters%NEntity), stat = IError)
          AccumulatedIncDisplacementSoil = 0.0

          allocate(AccumulatedDisplacementSoil(Counters%N), stat = IError)
          AccumulatedDisplacementSoil = 0.0  
                  
          allocate(AccelerationSoil(Counters%N, Counters%NEntity), stat = IError)
          AccelerationSoil = 0.0
          
          
          
          
          
          
          
          
          
          allocate(AccelerationSoilGeneralizedAlpha(Counters%N, Counters%NEntity), stat = IError)
          AccelerationSoilGeneralizedAlpha = 0.0 !This vector is populated from the acceleration of the material points 

          
          
          
          
          
          
          allocate(TotalVelocitySoil(Counters%N, Counters%NEntity), stat = IError)
          TotalVelocitySoil = 0.0
          
          !allocate(TotalVelocitySoil_VolLockSmooth(Counters%N_VolLockSmooth, Counters%NEntity), stat = IError)
          !TotalVelocitySoil_VolLockSmooth = 0.0
          
          allocate(TotalVelocitySoilPrevious(Counters%N, Counters%NEntity), stat = IError)
          TotalVelocitySoilPrevious = 0.0

          allocate(RateofMomentum(Counters%N, Counters%NEntity), stat = IError)
          RateofMomentum = 0.0

          allocate(TotalDisplacementWater(Counters%N), stat = IError)
          TotalDisplacementWater = 0.0
          
          allocate(PhaseDisplacementWater(Counters%N), stat = IError)
          PhaseDisplacementWater = 0.0

          allocate(TotalVelocityWater(Counters%N, Counters%NEntity), stat = IError)
          TotalVelocityWater = 0.0
          
          
          
          
          
          
          allocate(AccelerationWaterGeneralizedAlpha(Counters%N, Counters%NEntity), stat = IError)
          AccelerationWaterGeneralizedAlpha = 0.0 !This vector is populated from the acceleration of the material points 
          
          
          
          
          
          
          allocate(TotalVelocityWaterPrevious(Counters%N, Counters%NEntity), stat = IError)
          TotalVelocityWaterPrevious = 0.0
          
          allocate(IncrementalDisplacementWater(Counters%N, Counters%NEntity), stat = IError)
          IncrementalDisplacementWater = 0.0
          
          allocate(IncrementalDisplacementWaterPrevious(Counters%N, Counters%NEntity), stat = IError)
          IncrementalDisplacementWaterPrevious = 0.0
          
          allocate(AccumulatedDisplacementWater(Counters%N), stat = IError)
          AccumulatedDisplacementWater = 0.0  
          
          allocate(TotalDisplacementGas(Counters%N), stat = IError)
          TotalDisplacementGas = 0.0

          allocate(TotalVelocityGas(Counters%N, Counters%NEntity), stat = IError)
          TotalVelocityGas = 0.0
          
          allocate(TotalVelocityGasGeneralizedAlpha(Counters%N, Counters%NEntity), stat = IError)
          TotalVelocityGas = 0.0
          
          

          
          allocate(NonAdvectiveFluxAirInWater(Counters%N, Counters%NEntity), stat = IError)
          NonAdvectiveFluxAirInWater = 0.0

          allocate(NonAdvectiveFluxVapourInGas(Counters%N, Counters%NEntity), stat = IError)
          NonAdvectiveFluxVapourInGas = 0.0

          allocate(AdvectiveFluxDarcyWater(Counters%N, Counters%NEntity), stat = IError)
          AdvectiveFluxDarcyWater = 0.0

          allocate(AdvectiveFluxDarcyAir(Counters%N, Counters%NEntity), stat = IError)
          AdvectiveFluxDarcyAir = 0.0
          
          allocate(ThermalConductionFlux(Counters%N, Counters%NEntity), stat = IError)
          ThermalConductionFlux = 0.0

          allocate(FReaction(Counters%N, Counters%NEntity), stat = IError)
          FReaction = 0.0

          allocate(FReactionWater(Counters%N, Counters%NEntity), stat = IError)
          FReactionWater = 0.0

          allocate(FReactionGas(Counters%N, Counters%NEntity), stat = IError)
          FReactionGAs = 0.0
          
          allocate(NodalUnitMassGradient(Counters%Sum_NodTot, NDIM), stat = IError)
          NodalUnitMassGradient = 0.0
          
          allocate(ExtFlow(Counters%Sum_NodTot), stat = IError)
          ExtFlow = 0.0

          allocate(RateOfFlux(Counters%Sum_NodTot), stat = IError)
          RateOfFlux = 0.0
          
          allocate(IntFlow(Counters%Sum_NodTot), stat = IError)
          IntFlow = 0.0
          
          allocate(IncrementalPressure(Counters%Sum_NodTot), stat = IError)
          IncrementalPressure = 0.0

          allocate(TotalPressure(Counters%Sum_NodTot), stat = IError)
          TotalPressure = 0.0
          
          allocate(SubIncrementalDisplacement(Counters%N), stat = IError)
          SubIncrementalDisplacement = 0.0
      
          allocate(SubIncrementalPressure(Counters%Sum_NodTot), stat = IError)
          SubIncrementalPressure = 0.0
          
        end subroutine InitialiseNodalArrays

      
        subroutine InitialiseSurfaceReaction()
        !*********************************************************************************
        !
        !Function: read GOM file and determine the surface reaction force output
        !
        !**********************************************************************************
        implicit none
          ! local variables       
          character(len=MAX_FILENAME_LENGTH) :: FileName
          integer(INTEGER_TYPE) :: I, J, SurfaceID, IError, ios !NSurfaces,
          
          integer(INTEGER_TYPE) :: NSurfaces !, dimension(Counters%NPatches) 
          
          integer(INTEGER_TYPE) :: DumE, DumI, NumSideNodes !--> note DumE is unused
          
          !integer(INTEGER_TYPE), dimension(:), allocatable :: DumN ! --> might need to comment this
          real(REAL_TYPE), dimension(:), allocatable :: DumR ! --> real number for the control points
          
          
          character(len=255) :: DumS
          character(len=255) :: BName
          character(len=255), dimension(MAXOUTPUTSURFACES) :: DumSurfaceName !, Counters%NPatches
          
          ! Patch ID
          integer(INTEGER_TYPE) :: PatchID
          
          ! multipatch variable 
          real(REAL_TYPE), dimension(NVECTOR) :: ControlPointSurface ! control point coordinates
          integer(INTEGER_TYPE) :: IControlPoint_Surface
          integer(INTEGER_TYPE) :: IControlPoint_Full
          integer(INTEGER_TYPE) :: counter
          integer(INTEGER_TYPE) :: FirstIndex, LastIndex
          
          !External :: FExist
          
          ! initialize
          NSurfaces = 0
          
          !allocate(IsReactionNode(Counters%Sum_NodTot), stat = Ierror) ! maybe this needs to be multipatch unique global control points
          allocate(IsReactionNode(NURBS%Maximum_NControlPoints, Counters%NPatches), stat = Ierror) !NumberOfGlobalControlPointsUniqueMultipatch), stat = Ierror) 
          !allocate(IsReactionNode(NumberOfGlobalControlPointsUniqueMultipatch), stat = Ierror) !NumberOfGlobalControlPointsUniqueMultipatch), stat = Ierror) 
          IsReactionNode = .false. ! initialise
            
          
          if (ELEMENTTYPE==TETRAOLD) then
            NumSideNodes = N_BOUNDARY_NODES_HOE
          else
            NumSideNodes = ELEMENTBOUNDARYNODES
          end if
          
          ! This is for control point number 
          !allocate(DumN(NumSideNodes), Stat = Ierror) ! NumSideNodes --> ELEMENTBOUNDARYNODES
          !DumN = -1
          
          ! This is for control point coordinate
          allocate(DumR(NumSideNodes*NVECTOR), Stat = Ierror) ! NumSideNodes --> ELEMENTBOUNDARYNODES
          DumR = -1
          
          
          ! read GOM file
          FileName=Trim(CalParams%FileNames%ProjectName)//'.GOM'
          
          if (FExist(FileName) ) then
            open(GOMunit,FILE=FileName)
          end if
          DumSurfaceName = '/'

          
          !!If the number of surface is not know, then the GOM file has to be read twice
          !*************************************************
          !First read: determine NSurfaces and assign name 
          NSurfaces = 0 
          do
          read(GOMunit, '(A)', iostat=ios) BName 
          !if (trim(BName) == '$$PATCH_ID') then
          !    read(GOMUnit, *) PatchID
          !    
          !else 
              
        
             
          if (trim(BName) == '$$START_OUTPUT_REACTION_FORCES') then
           read(GOMUnit, *) DumI
           if (DumI < 1) EXIT
           do J=1,DumI                   ! I changed this to real number coordinates
            read(GOMunit, *) DumS, DumE, DumR(1:(NumSideNodes*NVECTOR)) !DumN(1:NumSideNodes)
            !DumSurfaceName(J) = DumS
            !do while (trim(DumS) /= '$$END_OUTPUT_REACTION_FORCES')
            do I=1,MAXOUTPUTSURFACES
              if (DumSurfaceName(I) == DumS) then !the surface has already been specified before !, PatchID
                 EXIT
              else if (DumSurfaceName(I) == '/') then !first non-assigned position of the array !, PatchID
                DumSurfaceName(I) = DumS !Assign the surface name to this position !, PatchID
                NSurfaces = NSurfaces + 1 !increment the number of surfaces !(PatchID) !(PatchID)
                EXIT
              end if
            end do
            
            
            counter = 0
            
            ! loop over all control points in the surface
            do IControlPoint_Surface = 1, NumSideNodes
                
                counter = counter + 1
                
                FirstIndex = 1 + ((counter - 1)*2) ! integer for 1, 3, 5, 7, ...
                LastIndex = 2 + ((counter - 1)*2) ! integer for 2, 4, 6, ...
                
                ControlPointSurface = DumR(FirstIndex:LastIndex)
                
                do PatchID = 1, Counters%NPatches
                
                ! loop over all control points in the patch
                do IControlPoint_Full = 1, Counters%NodTot(PatchID) 
                    
                    if (ControlPointSurface(1) == NodalCoordinates(IControlPoint_Full, 1, PatchID) &
                        .and. ControlPointSurface(2) == NodalCoordinates(IControlPoint_Full, 2, PatchID)) then 
               
                        !IsReactionNode(Multipatch_Connecting_Local_To_Global_ControlPoints(IControlPoint_Full, PatchID)) = .true. 
                        IsReactionNode(IControlPoint_Full, PatchID) = .true. 
                        ! turn reaction on 
                        
                    end if 
                    
                    
                end do ! loop over all control points in the patch
                
                end do ! loop over patches
                
            end do ! loop over all control points in the surface
            
            
            
            
            !read(GOMunit, *) DumE, DumN(1:6)
              !IsReactionNode(DumN(1:NumSideNodes), PatchID) = .true.
              !read(GOMUnit,*) DumS, DumE, DumN(1:6)
           end do
           EXIT
          else if (trim(BName) == '$$FINISH') then
           EXIT 
          end if
          end do
          if (NSurfaces < 1) RETURN !(PatchID)
          rewind(GOMUnit)
          Counters%NReactionSurfaceOutput = NSurfaces!maxval()
          
          !allocate vectors
          allocate(IsReactionNodeSurface(NURBS%Maximum_NControlPoints, Counters%NReactionSurfaceOutput, Counters%NPatches), stat = Ierror)
          !allocate(IsReactionNodeSurface(NumberOfGlobalControlPointsUniqueMultipatch, Counters%NReactionSurfaceOutput, Counters%NPatches), stat = Ierror)
          
          allocate(OutputSurfaceName(Counters%NReactionSurfaceOutput), stat = IError) !, Counters%NPatches
          
          OutputSurfaceName(1:NSurfaces) = DumSurfaceName(1:NSurfaces) !, PatchID !, PatchID !(PatchID) !(PatchID)
          IsReactionNodeSurface = .false.
          
          !**************************************************************
          !Second read: Assign the nodes to the surfaces
          do
            read(GOMunit, '(A)', iostat=ios) BName 
          !if (trim(BName) == '$$PATCH_ID') then
          !    read(GOMUnit, *) PatchID
          !  
          !else 
              
            if (trim(BName) == '$$START_OUTPUT_REACTION_FORCES') then 
              read(GOMUnit, *) DumI
              do J=1,DumI
                read(GOMunit, *) DumS, DumE, DumR(1:(NumSideNodes*NVECTOR)) !DumN(1:NumSideNodes)
                do I=1,NSurfaces!(PatchID)
                 if (OutputSurfaceName(I) == DumS) then !the surface has already been specified before !, PatchID
                  SurfaceID = I !get surfaceID
                  
                  
                  
                  
                  
                  counter = 0
            
            ! loop over all control points in the surface
            do IControlPoint_Surface = 1, NumSideNodes
                
                counter = counter + 1
                
                FirstIndex = 1 + ((counter - 1)*2) ! integer for 1, 3, 5, 7, ...
                LastIndex = 2 + ((counter - 1)*2) ! integer for 2, 4, 6, ...
                
                ControlPointSurface = DumR(FirstIndex:LastIndex)
                
                do PatchID = 1, Counters%NPatches
                
                ! loop over all control points in the patch
                do IControlPoint_Full = 1, Counters%NodTot(PatchID) 
                    
                    if (ControlPointSurface(1) == NodalCoordinates(IControlPoint_Full, 1, PatchID) &
                        .and. ControlPointSurface(2) == NodalCoordinates(IControlPoint_Full, 2, PatchID)) then 
               
                        !IsReactionNode(IControlPoint_Surface, PatchID) = .true. 
                        ! turn reaction on 
                        
                        !IsReactionNodeSurface(DumN(1:NumSideNodes),SurfaceID, PatchID) = .true.
                        !IsReactionNodeSurface(Multipatch_Connecting_Local_To_Global_ControlPoints(IControlPoint_Full,PatchID),SurfaceID, PatchID) = &
                        !    .true.
                        IsReactionNodeSurface(IControlPoint_Full,SurfaceID, PatchID) = &
                            .true.
                    end if 
                    
                    
                end do ! loop over all control points in the patch
                
                end do ! loop over patches
                
            end do ! loop over all control points in the surface
            
            
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  !IsReactionNodeSurface(DumN(1:NumSideNodes),SurfaceID, PatchID) = .true.
                  
                  
                  
                  
                  EXIT
                  
                 end if
                end do
              end do
              EXIT
              
            else if (trim(BName) == '$$FINISH') then
              EXIT 
            end if
          end do
          !***************************************************************
          !if (allocated(DumN)) then
          !  deallocate (DumN, stat = IError)
          !end if
          
          if (allocated(DumR)) then
            deallocate (DumR, stat = IError)
          end if
          
          
          !call DetermineReactionElements()        
          
        end subroutine InitialiseSurfaceReaction
    
        
        
        
        
    !    subroutine DetermineReactionElements()
    !!********************************************************************************
    !!
    !!Function: Determine the elements to be considered for the integration of the reaction forces
    !!
    !!********************************************************************************
    !implicit none
    !logical, dimension(Counters%Sum_nel, Counters%NPatches) :: IsReactionElm !Counters%NEl !Counters%Sum_nel
    !integer(INTEGER_TYPE) :: Inode, NElemNod, J, IEle, Ierror, cont !NReactionElem,
    !
    !integer(INTEGER_TYPE), dimension(Counters%NPatches) :: NReactionElem
    !
    !integer(INTEGER_TYPE) :: IPatch, Max_NReactionElem, GlobalNodeID
    !
    !! allocate to the number of patches
    !allocate( Counters%NElemReactions(Counters%NPatches) )
    !
    !! initialise
    !Counters%NElemReactions = 0
    !IsReactionElm = .false.
    !NReactionElem = 0
    !
    !do IPatch = 1, Counters%NPatches ! loop over patches ! --> local numbering 
    !do INode = 1, Counters%NodTot(IPatch) ! loop over control points ! --> local numbering
    !    
    !    !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode, IPatch)
    !    
    !  if (IsReactionNode(INode, IPatch)) then ! --> local numbering 
    !    NElemNod = GetNElmOfNode(INode, IPatch) ! number of elements associated with a node
    !    do J=1,NElemNod ! loop over number of elements associated with the node
    !      IEle = GetElmIOfNode(INode,J, IPatch)
    !      if (.not.IsReactionElm(IEle, IPatch)) then
    !        IsReactionElm (IEle, IPatch) = .true. ! we are turning on reaction element
    !        NReactionElem(IPatch) = NReactionElem(IPatch) + 1 ! counting how many elements that are linked to reaction
    !      end if
    !    end do
    !   end if
    !end do ! loop over control points
    !end do ! loop over patches
    !  
    !! finding the maximum number of reaction elements
    ! Counters%NElemReactions = NReactionElem!maxval()
    ! Max_NReactionElem = MAXVAL(NReactionElem)
    ! 
    ! !allocate(ConsideredElemReaction(NReactionElem), stat = Ierror)
    ! allocate(ConsideredElemReaction(Max_NReactionElem, Counters%NPatches), stat = Ierror)
    ! 
    ! ! initialize
    ! ConsideredElemReaction = 0
    ! cont = 0
    ! 
    ! do IPatch = 1, Counters%NPatches ! Loop over patches
    !do IEle=1, Counters%NEl(IPatch)!NURBS%nel_NURBS(IPatch)!Counters%NEl
    !  if(IsReactionElm(IEle, IPatch))then
    !    cont = cont + 1
    !    ConsideredElemReaction(cont, IPatch) = IEle ! we are setting the element numbers for each patch
    !    ! --> this is local numbering
    !  end if
    !end do ! Loop over elements
    !
    !cont = 0 ! reset counter before moving onto the next patch 
    !
    ! end do ! Loop over patches
    !   
    !    end subroutine DetermineReactionElements
        
        
      !subroutine InitialiseGPGlobalPositionArrays()
      !!*********************************************************************************
      !!
      !! Function: Determine global position of Gauss points for axisymmetric.
      !!
      !!**********************************************************************************
      !implicit none
      !
      !    integer(INTEGER_TYPE) :: iError
      !
      !    if ( .not. ISAXISYMMETRIC ) RETURN
      !    
      !    if ( allocated(GPGlobalPositionElement) ) then
      !      deallocate( GPGlobalPositionElement, stat=iError )
      !      call DeAllocationError(iError, 'GPGlobalPositionElement', 'initialiseGPGlobalPositionArrays')
      !    end if
      !
      !    allocate( GPGlobalPositionElement(NDIM, ELEMENTGAUSSPOINTS, Counters%NEl), stat=iError )
      !    call AllocationError(iError, 'GPGlobalPositionElement', 'initialiseGPGlobalPositionArrays')
      !    call setGPGlobalPositionElement(GPGlobalPositionElement)
      !
      !end subroutine initialiseGPGlobalPositionArrays
      !
      !
      !subroutine setGPGlobalPositionElement(GPGlobalPositionElement)
      !
      !    implicit none
      !
      !    real(REAL_TYPE), dimension(NDIM, ELEMENTGAUSSPOINTS, Counters%NEl), intent(out):: GPGlobalPositionElement
      !
      !    ! local variables
      !    integer(INTEGER_TYPE) :: IElement, IGaussPoint
      !    real(REAL_TYPE), dimension(NDIM, ELEMENTGAUSSPOINTS) :: PosGP
      !    real(REAL_TYPE), dimension(ELEMENTGAUSSPOINTS) :: WeiGP
      !    real(REAL_TYPE), dimension(ELEMENTNODES) :: GaussPointShapeValues
      !    real(REAL_TYPE), dimension(ELEMENTNODES, NDIM, ELEMENTGAUSSPOINTS) :: DShapeValues
      !    real(REAL_TYPE), dimension(NDIM) :: LocalPosition
      !
      !    GPGlobalPositionElement = 0.0
      !    LocalPosition = 0.0
      !
      !    do IGaussPoint = 1, ELEMENTGAUSSPOINTS
      !      ! Detemine local poisiton of the GaussPoint
      !      call GaussPointLocalCoordinates(IGaussPoint, WeiGP(IGaussPoint), PosGP(:,IGaussPoint))
      !      ! Detemine the shape function of the GaussPoint
      !      call ShapeFunctionData(PosGP(:,IGaussPoint), ELEMENTNODES, GaussPointShapeValues, DShapeValues(:,:,IGaussPoint))
      !    end do
      !
      !    do IElement = 1, Counters%NEl
      !      do IGaussPoint = 1, ELEMENTGAUSSPOINTS
      !        ! Determine the global location of IGaussPoint
      !        call CoordGaussPointLocalToGlobal(IElement, NodalCoordinatesUpd, GaussPointShapeValues, LocalPosition) 
      !        GPGlobalPositionElement(:, IGaussPoint,IElement) = LocalPosition
      !      end do
      !    end do
      !
      !end subroutine setGPGlobalPositionElement
      
      
      subroutine CoordGaussPointLocalToGlobal(IElement, LNodalCoordinates, GaussPointShapeValues, GaussPointGlobalCoord, IPatch)
      !**********************************************************************
      !
      !    Function:  Determines the global coordinates associated with the local
      !               coordinates of GaussPoint inside a considered element.
      !               
      !
      !     IElement : ID of the considered element
      !     LNodalCoordinates : Nodal coordinates used for updating the integration point locations
      !
      ! Implemented in the frame of the MPM project.
      !
      !**********************************************************************

      implicit none

        integer(INTEGER_TYPE), intent(in) :: IElement
        real(REAL_TYPE), dimension(:,:), intent(in) :: LNodalCoordinates
        real(REAL_TYPE), dimension(:), intent(in) :: GaussPointShapeValues
        real(REAL_TYPE), dimension(:), intent(inout) :: GaussPointGlobalCoord
        integer(INTEGER_TYPE), intent(in) :: IPatch

        ! Local variables
        integer(INTEGER_TYPE) :: NodeID
        integer(INTEGER_TYPE) :: I, IDim
        
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary
        
        ! initialize global coordinate of Gauss Point
        GaussPointGlobalCoord = 0.0

        do I = 1, ELEMENTNODES ! Loop over nodes of IElement
          NodeID = iabs(ElementConnectivities(I, IElement, IPatch) )
          do IDim = 1, NVECTOR ! Loop over dimensions of IElement
            GaussPointGlobalCoord(IDim) = GaussPointGlobalCoord(IDim) + LNodalCoordinates(NodeID, IDim) * GaussPointShapeValues(I)
          end do
        end do

      end subroutine CoordGaussPointLocalToGlobal
      
       subroutine SurfaceNodesCoordinates(CoordSurfNodes, Connectiv, SidesCount)
        !**********************************************************************
        !
        ! Function:  Store in a matrix the coordinates of the nodes belonging to a user-defined surface (soil or phreatic)
        ! which can be used in the pressure or stress initialization procedure
        ! Currently not completely implemented for 3D problems
        !
        !**********************************************************************
             
        
        integer(INTEGER_TYPE), dimension(:,:), intent(in) :: Connectiv
        integer(INTEGER_TYPE),intent(in) :: SidesCount
        real(REAL_TYPE), dimension(:,:), allocatable, intent(out) :: CoordSurfNodes
        real(REAL_TYPE), dimension(:), allocatable :: Dummy_row
        integer(INTEGER_TYPE) :: I, J, IError, NodeID, k, n, m
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch_Temporary
                          
        if (NDIM == 3) then
            allocate(CoordSurfNodes((SidesCount+2), N_BOUNDARY_NODES_HOE), stat = IError)
            allocate(Dummy_row(N_BOUNDARY_NODES_HOE), stat = IError)
            n = SidesCount+2
            m = N_BOUNDARY_NODES_HOE
        elseif (NDIM == 2) then
            allocate(CoordSurfNodes((SidesCount+1), ELEMENTBOUNDARYNODES), stat = IError)
            allocate(Dummy_row(ELEMENTBOUNDARYNODES), stat = IError)
            n =  SidesCount+1
            m = ELEMENTBOUNDARYNODES
        end if
       
               if (NDIM == 2) then
            do J = 1, NDIM
                do I = 1, SidesCount
                 NodeID = Connectiv(J,I) !
                    if (J == 1) then
                        CoordSurfNodes(I, :) = (NodalCoordinates(NodeID,:,IPatch_Temporary))
                    else
                        if (Connectiv(J,I) /= Connectiv(J-1,I)) then
                            CoordSurfNodes(SidesCount+1,:) = (NodalCoordinates(NodeID,:,IPatch_Temporary))
                        end if
                    end if
                end do
            end do
        else if (NDIM == 3) then
            do J = 1, NDIM
                do I = 1, SidesCount
                NodeID = Connectiv(J,I) !
                    if (J == 1) then
                        CoordSurfNodes(I, :) = (NodalCoordinates(NodeID,:,IPatch_Temporary))
                    else if (J == 2) then
                        if (Connectiv(J,I) /= Connectiv(J-1,I)) then
                            CoordSurfNodes(SidesCount+1,:) = (NodalCoordinates(NodeID,:,IPatch_Temporary))
                        end if
                    else
                        if ((Connectiv(J,I) /= Connectiv(J-1,I)) .and. (Connectiv(J,I) /= Connectiv(J-2,I))) then
                            CoordSurfNodes(SidesCount+2,:) = (NodalCoordinates(NodeID,:,IPatch_Temporary))
                        end if
                    end if
                end do
            end do
        end if
        
        do i = 1,n
            do j = i+1,n
                if ( CoordSurfNodes(j,1) < CoordSurfNodes(i,1) ) then
                    do k=1, m
                        Dummy_row(k) = CoordSurfNodes(i,k)
                        CoordSurfNodes(i,k) =CoordSurfNodes(j,k)
                        CoordSurfNodes(j,k) = Dummy_row(k)
                    end do
                end if
            end do
        end do
        
       end subroutine SurfaceNodesCoordinates
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
        
        
        
        
!        subroutine Generate_STARTELEMCON_ElementConnectivities(NXiKnotEntriesUNIQUE, NEtaKnotEntriesUNIQUE, &
!                                                                NXiKnotOrder, NEtaKnotOrder, &
!                                                                NNodes, &
!                                                                IPatch)
!          ! this subroutine generates the element connectivity for IGA
!          
!          implicit none 
!          
!          ! initiate variables 
!          integer(INTEGER_TYPE) :: ii, jj, kk
!          integer(INTEGER_TYPE) :: counter_Xknot, counter_JumpWithinRows, counter_JumpWithinColumns
!          integer(INTEGER_TYPE) :: KnotSpan_Y_Order 
!          integer(INTEGER_TYPE) :: KnotSpan_X_Order
!          integer(INTEGER_TYPE) :: counter
!          integer(INTEGER_TYPE) :: Size_KnotSpan_X_Unique
!          integer(INTEGER_TYPE) :: Size_KnotSpan_Y_Unique
!          integer(INTEGER_TYPE) :: rows
!          
!          integer(INTEGER_TYPE), intent(in) :: NXiKnotEntriesUNIQUE
!          integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntriesUNIQUE
!          integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
!          integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder
!          integer(INTEGER_TYPE), intent(in) :: NNodes
!          
!          ! Multipatch variables 
!          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
!          integer(INTEGER_TYPE), intent(in) :: IPatch
!          
!          !NXiKnotEntriesUNIQUE
!          !NEtaKnotEntriesUNIQUE
!          
!          !ElementConnectivities = 0
!!          
!!          KnotSpan_X_Order = NXiKnotOrder
!!          KnotSpan_Y_Order = NEtaKnotOrder
!!          
!!          Size_KnotSpan_X_Unique = NXiKnotEntriesUNIQUE
!!          Size_KnotSpan_Y_Unique = NEtaKnotEntriesUNIQUE
!!          
!!          
!!          ii = 1
!!          counter_Xknot = 0 
!!          counter_JumpWithinRows = (Size_KnotSpan_X_Unique-1)+KnotSpan_X_Order 
!!          counter_JumpWithinColumns = KnotSpan_Y_Order+1 
!!          
!!          counter = 1
!!          rows = 0
!!          
!!          !% for ii = 1: (Size_KnotSpan_X_Unique-1)
!!          do jj = 1,KnotSpan_Y_Order+1
!!              do kk = 1,KnotSpan_X_Order+1 
!!            
!! 
!!                  !ElementConnectivities(((KnotSpan_X_Order+1)*(KnotSpan_Y_Order+1))-counter ,  nel_NURBS-rows)  &                       
!!                ElementConnectivities(counter ,  nel_NURBS(IPatch)-rows, IPatch)  &
!!                  = NNodes - counter_Xknot 
!!                  counter = counter + 1 
!!                  counter_Xknot = counter_Xknot + 1 
!!                  
!!              end do
!!              !rows = rows+1
!!              counter_Xknot = counter_Xknot + (Size_KnotSpan_X_Unique-1)-1 
!!              !counter_Xknot = counter_Xknot + Size_KnotSpan_X_Unique
!!          end do       
!!    
!!    
!!    
!!    
!!          counter = 1 
!!          rows = rows+1 !2
!!!% for jj = 1: (Size_KnotSpan_Y_Unique-1)-1 % knot spans - 1 
!!  
!!        
!!          do ii = 1,(Size_KnotSpan_X_Unique-1)-1 !% knot spans - 1
!!              ElementConnectivities(:, nel_NURBS(IPatch)-rows,IPatch) = ElementConnectivities(:, nel_NURBS(IPatch), IPatch) - counter
!!              rows = rows +1 
!!              counter = counter + 1 
!!          end do
!!    
!!          counter = counter + (KnotSpan_Y_Order)
!!    
!!    
!!!% end    
!!    
!!    
!!    
!!    
!!!% counter = 1; 
!!!% rows = 2
!!
!!          do jj = 1, (Size_KnotSpan_Y_Unique-1)-1 !% knot spans - 1 
!!              do ii = 1, (Size_KnotSpan_X_Unique-1) !% knot spans - 1
!!                  ElementConnectivities(:, nel_NURBS(IPatch)-rows,IPatch) = ElementConnectivities(:, nel_NURBS(IPatch), IPatch) - counter 
!!                  rows = rows +1 
!!                  counter = counter + 1     
!!              end do
!!
!!              counter = counter + (KnotSpan_Y_Order)
!!          end do
!!
!!    
!!!%     rows = rows +1; 
!!!%     counter = 0; 
!!!% end        
!!
!!        
!!!%         end  
!!!%     end
!!% end 
!
!          
!          ! I need to comment all above this
!          !ElementConnectivities = IEN
!          
!                                                                end subroutine 
!       
!       
!       
!                                                                
!                                                                
!                                                                
!                                                                
!        subroutine Generate_STARTELEMCON_ElementConnectivities_3D(NXiKnotEntriesUNIQUE, NEtaKnotEntriesUNIQUE, NZetaKnotEntriesUNIQUE, &
!                                                                NXiKnotOrder, NEtaKnotOrder, NZetaKnotOrder, &
!                                                                NNodes)
!          ! this subroutine generates the element connectivity for IGA
!          
!          implicit none 
!          
!          
!          
!          
!          integer(INTEGER_TYPE) :: NumberOfControlPointsInXDirection
!          
!          integer(INTEGER_TYPE) :: NumberOfControlPointsInYDirection
!          
!          integer(INTEGER_TYPE) :: NumberOfControlPointsInZDirection
!          
!          integer(INTEGER_TYPE) :: NumberOfControlPointsInXYPlane
!          
!          integer(INTEGER_TYPE) :: NumberOfControlPointsInXZPlane
!          
!          integer(INTEGER_TYPE) :: NumberOfControlPointsForEachKnotSpan
!          
!          ! initiate variables 
!          integer(INTEGER_TYPE) :: ii, jj, kk, ll
!          integer(INTEGER_TYPE) :: counter_Xknot, counter_JumpWithinRows, counter_JumpWithinColumns
!          integer(INTEGER_TYPE) :: KnotSpan_Y_Order 
!          integer(INTEGER_TYPE) :: KnotSpan_X_Order
!          integer(INTEGER_TYPE) :: KnotSpan_Z_Order
!          integer(INTEGER_TYPE) :: counter
!          integer(INTEGER_TYPE) :: Size_KnotSpan_X_Unique
!          integer(INTEGER_TYPE) :: Size_KnotSpan_Y_Unique
!          integer(INTEGER_TYPE) :: Size_KnotSpan_Z_Unique
!          integer(INTEGER_TYPE) :: rows
!          
!          
!          integer(INTEGER_TYPE), intent(in) :: NXiKnotEntriesUNIQUE
!          integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntriesUNIQUE
!          integer(INTEGER_TYPE), intent(in) :: NZetaKnotEntriesUNIQUE
!          integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
!          integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder
!          integer(INTEGER_TYPE), intent(in) :: NZetaKnotOrder
!          integer(INTEGER_TYPE), intent(in) :: NNodes
!          
!          !NXiKnotEntriesUNIQUE
!          !NEtaKnotEntriesUNIQUE
!          
!          !ElementConnectivities = IEN
!          
!!          KnotSpan_X_Order = NXiKnotOrder
!!          KnotSpan_Y_Order = NEtaKnotOrder
!!          KnotSpan_Z_Order = NZetaKnotOrder
!!          
!!          
!!          Size_KnotSpan_X_Unique = NXiKnotEntriesUNIQUE
!!          Size_KnotSpan_Y_Unique = NEtaKnotEntriesUNIQUE
!!          Size_KnotSpan_Z_Unique = NZetaKnotEntriesUNIQUE
!!          
!!          ii = 1
!!          counter_Xknot = 0 
!!          counter_JumpWithinRows = (Size_KnotSpan_X_Unique-1)+KnotSpan_X_Order 
!!          counter_JumpWithinColumns = KnotSpan_Y_Order+1 
!!          
!!          
!!          
!!          
!!          
!!          
!!          ! Abdelrahman: we are trying to calculate how many control points are in each direction
!!          NumberOfControlPointsInXDirection = NXiKnotEntriesUNIQUE+NXiKnotOrder-1 !4 for the cubic column exmaple 
!!          NumberOfControlPointsInYDirection = NEtaKnotEntriesUNIQUE+NEtaKnotOrder-1 !7 for the cubic column example 
!!          NumberOfControlPointsInZDirection = NZetaKnotEntriesUNIQUE+NZetaKnotOrder-1 !4 for the cubic column example 
!!          
!!          
!!          ! let's look at the x-y plane 
!!          NumberOfControlPointsInXYPlane = NumberOfControlPointsInXDirection*NumberOfControlPointsInYDirection
!!          ! let's look at the x-z plane 
!!          NumberOfControlPointsInXZPlane = NumberOfControlPointsInXDirection*NumberOfControlPointsInZDirection
!!          
!!          ! let's calculate how many control points for one element 
!!          NumberOfControlPointsForEachKnotSpan = (NXiKnotOrder+1)*(NEtaKnotOrder+1)*(NZetaKnotOrder+1)
!!          
!!          counter = 0 
!!
!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          !     1 Element  
!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          ! looping accross the number of control points in the x-z plane 
!!          !do ii = 0, NumberOfControlPointsInXZPlane-1
!!          !    
!!          !    ElementConnectivities(NumberOfControlPointsForEachKnotSpan-counter,1) = ii+1
!!          !    counter = counter + 1
!!          !end do 
!!          
!!          do jj =0,NumberOfControlPointsInZDirection-1
!!              do ii = 0, NumberOfControlPointsInXZPlane-1
!!              
!!                  ElementConnectivities(NumberOfControlPointsForEachKnotSpan-counter,1) = ii+1+(NumberOfControlPointsInXYPlane*jj)
!!                  counter = counter + 1
!!                  
!!                  
!!                  !counter = counter + 1
!!                  !ElementConnectivities(counter,1) = ii+1+(NumberOfControlPointsInXYPlane*jj)
!!                  
!!              end do 
!!          end do 
!!          
!!          
!!          ! loop accross the number of number of control points for each knot span 
!!          !do jj = 1,NumberOfControlPointsInZDirection
!!          !
!!          !    do ii = 0, NumberOfControlPointsInXZPlane-1
!!          !        ElementConnectivities(NumberOfControlPointsForEachKnotSpan-counter,1) = ii+1  
!!          !        counter = counter + 1
!!          !        
!!          !    end do 
!!          !    
!!          !end do 
!!          
!!          
!!          
!!          !ElementConnectivities(NumberOfControlPointsForEachKnotSpan-counter,1)
!!          
!!          
!!          
!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          ! loop over elements so that we can establish the connectivity 
!!          counter = 0
!!          do ii = 2, nel_NURBS
!!              
!!              
!!              ElementConnectivities(:,ii) = ElementConnectivities(:,1) + ( (ii-1)*NumberOfControlPointsInXDirection)
!!              counter = counter + 1
!!              
!!              
!!          end do 
!!          
!!    !      16    12     8     4
!!    !15    11     7     3
!!    !14    10     6     2
!!    !13     9     5     1
!!          
!!    !      ElementConnectivities= 0
!!    !      ElementConnectivities(1,1) = 1
!!    !      ElementConnectivities(2,1) = 5 
!!    !      ElementConnectivities(3,1) = 9
!!    !      ElementConnectivities(4,1) = 13 
!!    !      ElementConnectivities(5,1) = 2
!!    !      ElementConnectivities(6,1) = 6
!!    !      ElementConnectivities(7,1) = 10
!!    !      ElementConnectivities(8,1) = 14
!!    !      ElementConnectivities(9,1) = 3
!!    !      ElementConnectivities(10,1) = 7
!!    !      ElementConnectivities(11,1) = 11
!!    !      ElementConnectivities(12,1) = 15
!!    !      ElementConnectivities(13,1) = 4
!!    !      ElementConnectivities(14,1) = 8
!!    !      ElementConnectivities(15,1) = 12
!!    !      ElementConnectivities(16,1) = 16
!!    !      
!!    !      
!!    !!      32    28    24    20
!!    !!31    27    23    19
!!    !!30    26    22    18
!!    !!29    25    21    17
!!    !      
!!    !      ElementConnectivities(17,1) = 17
!!    !      ElementConnectivities(18,1) = 21
!!    !      ElementConnectivities(19,1) = 25
!!    !      ElementConnectivities(20,1) = 29
!!    !      ElementConnectivities(21,1) = 18
!!    !      ElementConnectivities(22,1) = 22
!!    !      ElementConnectivities(23,1) = 26
!!    !      ElementConnectivities(24,1) = 30
!!    !      ElementConnectivities(25,1) = 19
!!    !      ElementConnectivities(26,1) = 23
!!    !      ElementConnectivities(27,1) = 27
!!    !      ElementConnectivities(28,1) = 31
!!    !      ElementConnectivities(29,1) = 20
!!    !      ElementConnectivities(30,1) = 24
!!    !      ElementConnectivities(31,1) = 28
!!    !      ElementConnectivities(32,1) = 32
!!    !      
!!    !      
!!    !!       48    44    40    36
!!    !!47    43    39    35
!!    !!46    42    38    34
!!    !!45    41    37    33
!!    !      
!!    !      ElementConnectivities(33,1) = 33
!!    !      ElementConnectivities(34,1) = 37 
!!    !      ElementConnectivities(35,1) = 41
!!    !      ElementConnectivities(36,1) = 45
!!    !      ElementConnectivities(37,1) = 34
!!    !      ElementConnectivities(38,1) = 38
!!    !      ElementConnectivities(39,1) = 42
!!    !      ElementConnectivities(40,1) = 46
!!    !      ElementConnectivities(41,1) = 35
!!    !      ElementConnectivities(42,1) = 39
!!    !      ElementConnectivities(43,1) = 43
!!    !      ElementConnectivities(44,1) = 47
!!    !      ElementConnectivities(45,1) = 36
!!    !      ElementConnectivities(46,1) = 40
!!    !      ElementConnectivities(47,1) = 44
!!    !      ElementConnectivities(48,1) = 48
!!    !      
!!    !!      64    60    56    52
!!    !!63    59    55    51
!!    !!62    58    54    50
!!    !!61    57    53    49
!!    !      
!!    !      ElementConnectivities(49,1) = 49
!!    !      ElementConnectivities(50,1) = 53
!!    !      ElementConnectivities(51,1) = 57
!!    !      ElementConnectivities(52,1) = 61
!!    !      ElementConnectivities(53,1) = 50
!!    !      ElementConnectivities(54,1) = 54
!!    !      ElementConnectivities(55,1) = 58
!!    !      ElementConnectivities(56,1) = 62
!!    !      ElementConnectivities(57,1) = 51
!!    !      ElementConnectivities(58,1) = 55
!!    !      ElementConnectivities(59,1) = 59
!!    !      ElementConnectivities(60,1) = 63
!!    !      ElementConnectivities(61,1) = 52
!!    !      ElementConnectivities(62,1) = 56
!!    !      ElementConnectivities(63,1) = 60
!!    !      ElementConnectivities(64,1) = 64
!!          
!!          
!!          
!!    !      ElementConnectivities= 0
!!    !      ElementConnectivities(64,1) = 4
!!    !      ElementConnectivities(63,1) = 8 
!!    !      ElementConnectivities(62,1) = 12
!!    !      ElementConnectivities(61,1) = 16 
!!    !      ElementConnectivities(60,1) = 3
!!    !      ElementConnectivities(59,1) = 7
!!    !      ElementConnectivities(58,1) = 11
!!    !      ElementConnectivities(57,1) = 15
!!    !      ElementConnectivities(56,1) = 2
!!    !      ElementConnectivities(55,1) = 6
!!    !      ElementConnectivities(54,1) = 10
!!    !      ElementConnectivities(53,1) = 14
!!    !      ElementConnectivities(52,1) = 1
!!    !      ElementConnectivities(51,1) = 5
!!    !      ElementConnectivities(50,1) = 9
!!    !      ElementConnectivities(49,1) = 13
!!    !      
!!    !      
!!    !!      32    28    24    20
!!    !!31    27    23    19
!!    !!30    26    22    18
!!    !!29    25    21    17
!!    !      
!!    !      ElementConnectivities(48,1) = 20
!!    !      ElementConnectivities(47,1) = 24
!!    !      ElementConnectivities(46,1) = 28
!!    !      ElementConnectivities(45,1) = 32
!!    !      ElementConnectivities(44,1) = 19
!!    !      ElementConnectivities(43,1) = 23
!!    !      ElementConnectivities(42,1) = 27
!!    !      ElementConnectivities(41,1) = 31
!!    !      ElementConnectivities(40,1) = 18
!!    !      ElementConnectivities(39,1) = 22
!!    !      ElementConnectivities(38,1) = 26
!!    !      ElementConnectivities(37,1) = 30
!!    !      ElementConnectivities(36,1) = 17
!!    !      ElementConnectivities(35,1) = 21
!!    !      ElementConnectivities(34,1) = 25
!!    !      ElementConnectivities(33,1) = 29
!!    !      
!!    !      
!!    !!       48    44    40    36
!!    !!47    43    39    35
!!    !!46    42    38    34
!!    !!45    41    37    33
!!    !      
!!    !      ElementConnectivities(32,1) = 36
!!    !      ElementConnectivities(31,1) = 40 
!!    !      ElementConnectivities(30,1) = 44
!!    !      ElementConnectivities(29,1) = 48
!!    !      ElementConnectivities(28,1) = 35
!!    !      ElementConnectivities(27,1) = 39
!!    !      ElementConnectivities(26,1) = 43
!!    !      ElementConnectivities(25,1) = 47
!!    !      ElementConnectivities(24,1) = 34
!!    !      ElementConnectivities(23,1) = 38
!!    !      ElementConnectivities(22,1) = 42
!!    !      ElementConnectivities(21,1) = 46
!!    !      ElementConnectivities(20,1) = 33
!!    !      ElementConnectivities(19,1) = 37
!!    !      ElementConnectivities(18,1) = 41
!!    !      ElementConnectivities(17,1) = 45
!!    !      
!!    !!      64    60    56    52
!!    !!63    59    55    51
!!    !!62    58    54    50
!!    !!61    57    53    49
!!    !      
!!    !      ElementConnectivities(16,1) = 52
!!    !      ElementConnectivities(15,1) = 56
!!    !      ElementConnectivities(14,1) = 60
!!    !      ElementConnectivities(13,1) = 64
!!    !      ElementConnectivities(12,1) = 51
!!    !      ElementConnectivities(11,1) = 55
!!    !      ElementConnectivities(10,1) = 59
!!    !      ElementConnectivities(9,1) = 63
!!    !      ElementConnectivities(8,1) = 50
!!    !      ElementConnectivities(7,1) = 54
!!    !      ElementConnectivities(6,1) = 58
!!    !      ElementConnectivities(5,1) = 62
!!    !      ElementConnectivities(4,1) = 49
!!    !      ElementConnectivities(3,1) = 53
!!    !      ElementConnectivities(2,1) = 57
!!    !      ElementConnectivities(1,1) = 61
!!          
!!          
!!!          ElementConnectivities(
!!!          ElementConnectivities
!!!          
!!!          
!!!          32    28    24    20
!!!    31    27    23    19
!!!    30    26    22    18
!!!29        25    21    17
!!          
!!          
!!          
!!          
!!          
!!          
!!          
!!    !            ElementConnectivities= 0
!!    !      ElementConnectivities(64,1) = 4
!!    !      ElementConnectivities(63,1) = 3 
!!    !      ElementConnectivities(62,1) = 2
!!    !      ElementConnectivities(61,1) = 1 
!!    !      ElementConnectivities(60,1) = 8
!!    !      ElementConnectivities(59,1) = 7
!!    !      ElementConnectivities(58,1) = 6
!!    !      ElementConnectivities(57,1) = 5
!!    !      ElementConnectivities(56,1) = 12
!!    !      ElementConnectivities(55,1) = 11
!!    !      ElementConnectivities(54,1) = 10
!!    !      ElementConnectivities(53,1) = 9
!!    !      ElementConnectivities(52,1) = 16
!!    !      ElementConnectivities(51,1) = 15
!!    !      ElementConnectivities(50,1) = 14
!!    !      ElementConnectivities(49,1) = 13
!!    !      
!!    !      
!!    !!      32    28    24    20
!!    !!31    27    23    19
!!    !!30    26    22    18
!!    !!29    25    21    17
!!    !      
!!    !      ElementConnectivities(48,1) = 17
!!    !      ElementConnectivities(47,1) = 18
!!    !      ElementConnectivities(46,1) = 19
!!    !      ElementConnectivities(45,1) = 20
!!    !      ElementConnectivities(44,1) = 21
!!    !      ElementConnectivities(43,1) = 22
!!    !      ElementConnectivities(42,1) = 23
!!    !      ElementConnectivities(41,1) = 24
!!    !      ElementConnectivities(40,1) = 25
!!    !      ElementConnectivities(39,1) = 26
!!    !      ElementConnectivities(38,1) = 27
!!    !      ElementConnectivities(37,1) = 28
!!    !      ElementConnectivities(36,1) = 29
!!    !      ElementConnectivities(35,1) = 30
!!    !      ElementConnectivities(34,1) = 31
!!    !      ElementConnectivities(33,1) = 32
!!    !      
!!    !      
!!    !!       48    44    40    36
!!    !!47    43    39    35
!!    !!46    42    38    34
!!    !!45    41    37    33
!!    !      
!!    !      ElementConnectivities(32,1) = 33
!!    !      ElementConnectivities(31,1) = 34 
!!    !      ElementConnectivities(30,1) = 35
!!    !      ElementConnectivities(29,1) = 36
!!    !      ElementConnectivities(28,1) = 37
!!    !      ElementConnectivities(27,1) = 38
!!    !      ElementConnectivities(26,1) = 39
!!    !      ElementConnectivities(25,1) = 40
!!    !      ElementConnectivities(24,1) = 41
!!    !      ElementConnectivities(23,1) = 42
!!    !      ElementConnectivities(22,1) = 43
!!    !      ElementConnectivities(21,1) = 44
!!    !      ElementConnectivities(20,1) = 45
!!    !      ElementConnectivities(19,1) = 46
!!    !      ElementConnectivities(18,1) = 47
!!    !      ElementConnectivities(17,1) = 48
!!    !      
!!    !!      64    60    56    52
!!    !!63    59    55    51
!!    !!62    58    54    50
!!    !!61    57    53    49
!!    !      
!!    !      ElementConnectivities(16,1) = 49
!!    !      ElementConnectivities(15,1) = 50
!!    !      ElementConnectivities(14,1) = 51
!!    !      ElementConnectivities(13,1) = 52
!!    !      ElementConnectivities(12,1) = 53
!!    !      ElementConnectivities(11,1) = 54
!!    !      ElementConnectivities(10,1) = 55
!!    !      ElementConnectivities(9,1) = 56
!!    !      ElementConnectivities(8,1) = 57
!!    !      ElementConnectivities(7,1) = 58
!!    !      ElementConnectivities(6,1) = 59
!!    !      ElementConnectivities(5,1) = 60
!!    !      ElementConnectivities(4,1) = 61
!!    !      ElementConnectivities(3,1) = 62
!!    !      ElementConnectivities(2,1) = 63
!!    !      ElementConnectivities(1,1) = 64
!!          
!!!          counter = 1
!!!          rows = 0
!!!          
!!!          !% for ii = 1: (Size_KnotSpan_X_Unique-1)
!!!          do ll = 1,KnotSpan_Z_Order+1
!!!              do jj = 1,KnotSpan_Y_Order+1
!!!                  do kk = 1,KnotSpan_X_Order+1 
!!!            
!!! 
!!!                  !ElementConnectivities(((KnotSpan_X_Order+1)*(KnotSpan_Y_Order+1))-counter ,  nel_NURBS-rows)  &                       
!!!                ElementConnectivities(counter ,  nel_NURBS-rows)  &
!!!                 !ElementConnectivities(((KnotSpan_X_Order+1)*(KnotSpan_Y_Order+1)*(KnotSpan_Z_Order+1))-counter ,  nel_NURBS-rows) &
!!!                    = NNodes - counter_Xknot 
!!!                  counter = counter + 1 
!!!                  counter_Xknot = counter_Xknot + 1 
!!!                  
!!!              end do
!!!              !rows = rows+1
!!!              counter_Xknot = counter_Xknot + (Size_KnotSpan_X_Unique-1)-1 
!!!              !counter_Xknot = counter_Xknot + Size_KnotSpan_X_Unique
!!!              end do       
!!!          end do 
!!!          
!!!    
!!!    
!!!    
!!!    
!!!          counter = 1 
!!!          rows = rows+1 !2
!!!!% for jj = 1: (Size_KnotSpan_Y_Unique-1)-1 % knot spans - 1 
!!!  
!!!        
!!!          do ii = 1,(Size_KnotSpan_X_Unique-1)-1 !% knot spans - 1
!!!              ElementConnectivities(:, nel_NURBS-rows) = ElementConnectivities(:, nel_NURBS) - counter
!!!              rows = rows +1 
!!!              counter = counter + 1 
!!!          end do
!!!    
!!!          counter = counter + (KnotSpan_Y_Order)
!!!    
!!!    
!!!!% end    
!!!    
!!!    
!!!    
!!!    
!!!!% counter = 1; 
!!!!% rows = 2
!!!
!!!          do jj = 1, (Size_KnotSpan_Y_Unique-1)-1 !% knot spans - 1 
!!!              do ii = 1, (Size_KnotSpan_X_Unique-1) !% knot spans - 1
!!!                  ElementConnectivities(:, nel_NURBS-rows) = ElementConnectivities(:, nel_NURBS) - counter 
!!!                  rows = rows +1 
!!!                  counter = counter + 1     
!!!              end do
!!!
!!!              counter = counter + (KnotSpan_Y_Order)
!!!          end do
!!
!!          
!!          
!!          
!!          
!!          
!!          
!!          
!!    
!!!%     rows = rows +1; 
!!!%     counter = 0; 
!!!% end        
!!
!!        
!!!%         end  
!!!%     end
!!!% end 
!
!          
!          
!          
!          
!                                                                end subroutine 
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                          
 ! MATLAB version of the ElementConnectivitity matrix           
!                                                                % This part is related to NURBS connectivity 
!% building up
!
!KnotSpan_X = [0
!              0
!              0
!              0
!              1
!              2
!              3
!              4
!              5
!              6
!              7
!              8
!              8
!              8
!              8];
!
!KnotSpan_Y = [0
!              0
!              0
!              0
!              1
!              2
!              3
!              4
!              5
!              6
!              7
!              8
!              8
!              8
!              8];
!
!NNodes = 121; % number of control points 
!
!KnotSpan_X_Order = 3; 
!KnotSpan_Y_Order = 3; 
!
!KnotSpan_X_Unique = unique(KnotSpan_X); 
!KnotSpan_Y_Unique = unique(KnotSpan_Y); 
!
!
!Size_KnotSpan_X_Unique = size(KnotSpan_X_Unique,1);
!Size_KnotSpan_Y_Unique = size(KnotSpan_Y_Unique,1);
!
!ElementConnectivities = zeros(  (Size_KnotSpan_X_Unique-1)*(Size_KnotSpan_Y_Unique-1), ... % x-direction size 
!                                (KnotSpan_X_Order+1)*(KnotSpan_Y_Order+1));  % y-direction size 
!
!
!ii = 1;
!counter_Xknot = 0; 
!counter_JumpWithinRows = (Size_KnotSpan_X_Unique-1)+KnotSpan_X_Order; 
!counter_JumpWithinColumns = KnotSpan_Y_Order+1; 
!% need to loop accross 
!% for ii = 1:(Size_KnotSpan_X_Unique-1)*(Size_KnotSpan_Y_Unique-1)
!%     for jj = 1:Size_KnotSpan_Y_Unique-1
!%         for kk = 1: Size_KnotSpan_X_Unique-1 
!%         
!
!
!counter = 0 ;
!rows = 1;
!% for ii = 1: (Size_KnotSpan_X_Unique-1)
!    for jj = 1:KnotSpan_Y_Order+1
!        for kk = 1:KnotSpan_X_Order+1 
!            
!            ElementConnectivities(rows, ((KnotSpan_X_Order+1)*(KnotSpan_Y_Order+1))-counter )  ...
!            = NNodes - counter_Xknot; 
!            counter = counter + 1; 
!            counter_Xknot = counter_Xknot + 1; 
!        end 
!        counter_Xknot = counter_Xknot + (Size_KnotSpan_X_Unique-1)-1; 
!    end        
!    
!    
!    
!counter = 1; 
!rows = 2
!% for jj = 1: (Size_KnotSpan_Y_Unique-1)-1 % knot spans - 1 
!    for ii = 1:(Size_KnotSpan_X_Unique-1)-1 % knot spans - 1
!        
!        ElementConnectivities(rows, :) = ElementConnectivities(1, :) - counter; 
!
!        rows = rows +1; 
!        counter = counter + 1; 
!    end
!    
!    counter = counter + (KnotSpan_Y_Order);
!    
!    
!% end    
!    
!    
!    
!    
!% counter = 1; 
!% rows = 2
!for jj = 1: (Size_KnotSpan_Y_Unique-1)-1 % knot spans - 1 
!    for ii = 1:(Size_KnotSpan_X_Unique-1) % knot spans - 1
!        
!        ElementConnectivities(rows, :) = ElementConnectivities(1, :) - counter; 
!
!        rows = rows +1; 
!        counter = counter + 1; 
!    end
!    
!    counter = counter + (KnotSpan_Y_Order);
!    
!    
!end    
                                                                
                                                                
                                                                
                                                                
                                                                
       
       
       
        
      end module ModMeshInfo
