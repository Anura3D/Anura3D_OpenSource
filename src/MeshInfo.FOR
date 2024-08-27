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
    
      use ModCounters
      use ModReadCalculationData
      use ModReadMaterialData
      use ModElementEvaluation
      use ModGlobalConstants
      use ModFileIO
      
      implicit none

        ! Nodal coordinates, updated nodal coordinates, The initial coordinates
        real(REAL_TYPE), dimension(:, :),  &
          allocatable :: NodalCoordinates, NodalCoordinatesUpd, NodalOriginalCoord
        real(REAL_TYPE), dimension(:, :),  &
          allocatable :: NodalPrescibedDisp ! Nodal prescribed displacement in 3 directions (Solid)
        real(REAL_TYPE), dimension(:, :),  &
          allocatable :: NodalPrescibedDispWater ! Nodal prescribed displacement in 3 directions (Water)
        real(REAL_TYPE), dimension(:, :),  &
          allocatable :: NodalPrescibedDispGas ! Nodal prescribed displacement in 3 directions (Gas)
        real(REAL_TYPE), dimension(:, :, :), &
          allocatable :: LoadValuesOnNodesSolid ! Load values (applied on nodes) for each distributed element face (solid)
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
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElementConnectivities
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
        integer(INTEGER_TYPE), dimension(:), allocatable :: ElementMaterialID ! Material ID for the element
        integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDof ! The global storage of degree of freedom in x-direction
        integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDofQuasiStatic ! The global storage of degree of freedom in ux,uy,uz,p
        logical, dimension(:), allocatable :: IsActiveElement ! Element switch, 1:Active element and 0:Inactive element
        real(REAL_TYPE), dimension(:), allocatable :: NodalDensity 
        real(REAL_TYPE), dimension(:), allocatable :: ActiveNodeElementVolume 
        logical, dimension(:), allocatable :: ActiveNode
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
        integer(INTEGER_TYPE), dimension(:), allocatable :: ActiveElement ! ID of the active elements
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
    

        !> stores material points indeces
        integer(INTEGER_TYPE), dimension(:, :), allocatable :: GetParticleIndex

        !> stores determinant of elements
        real(REAL_TYPE), dimension(:), allocatable :: ElementDeterminant
        
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
          
          integer(INTEGER_TYPE) :: I, J, L, IDElement, IDConnectivity, IFxN, &
                     NElements, NNodes, NFixNod, FixNodID, NRemFixNod, FixRemNodID, &
                     NLoadedElementSidesSolidNodes, NLoadedElementSidesWaterNodes, NLoadedElementSidesGasNodes, &
                     NLoadedElementSidesSolidMatPoints, NLoadedElementSidesWaterMatPoints, NLoadedElementSidesGasmatPoints, &
                     SoilSurfaceNumberofSides, PhreaticSurfaceNumberofSides, &
                     NLoadedElementSidesSolidNodesB, NLoadedElementSidesWaterNodesB, NLoadedElementSidesGasNodesB, &
                    NLoadedElementSidesSolidMatPointsB, NLoadedElementSidesWaterMatPointsB, NLoadedElementSidesGasmatPointsB 

          integer(INTEGER_TYPE) :: IError ! used for error control
          integer(INTEGER_TYPE), dimension(NVECTOR) :: NodFixBC, NodRemBC
          real(REAL_TYPE) :: DumR(4)
          
          call DestroyMeshData()
          
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
          
          ! open GOM file
          FileName=Trim(CalParams%FileNames%ProjectName)//'.GOM'
          if (FExist(FileName)) open(GOMunit, FILE=FileName)

          do
            read(GOMunit,'(A)') TName
            BName = TName
              
            if (trim(BName) == '$$STARTCOUNTERS') then
              read(GOMunit, *) NElements, NNodes  ! number of elements, number of nodes

              ! allocating allocatable arrays
              allocate(NodalCoordinates(NNodes, NVECTOR), stat = IError)
              allocate(NodalOriginalCoord(NNodes, NVECTOR), stat = IError)
              allocate(NodalCoordinatesUpd(NNodes, NVECTOR), stat = IError)
              NodalCoordinates = 0.0
              NodalOriginalCoord = 0.0
              NodalCoordinatesUpd = 0.0

              allocate(ElementConnectivities(ELEMENTNODES, NElements), stat = IError)
              ElementConnectivities = 0
              if (ELEMENTTYPE == TETRAOLD) then ! elementtype before v2018.2  
                allocate(ElementConnectivities10Node(N_NODES_HOE, NElements), stat = IError)
                ElementConnectivities10Node = 0
              end if  

              allocate(NodalPrescibedDisp(NNodes, NVECTOR), stat = IError)
              allocate(NodalPrescibedDispWater(NNodes, NVECTOR), stat = IError)
              allocate(NodalPrescibedDispGas(NNodes, NVECTOR), stat = IError)
              NodalPrescibedDisp = 1.d10
              NodalPrescibedDispWater = 1.d10          
              NodalPrescibedDispGas = 1.d10 

              allocate (ElementMaterialID(NElements), stat = IError)
              allocate(IsActiveElement(NElements), stat = IError)
              allocate(ActiveNode(NNodes), stat = IError)
              allocate(HydraulicHeadLoadedElemID(NElements), stat = IError)
              ElementMaterialID = 0
              IsActiveElement = .false.
              ActiveNode = .false.
              HydraulicHeadLoadedElemID = .false.
              

            else if (trim(BName) == '$$STARTNODES') then
              do I = 1, NNodes ! loop over all nodes
                read (GOMunit, *) (NodalCoordinates(I, J), J = 1, NVECTOR) ! nodal coordinates
                if ( ISAXISYMMETRIC ) call Assert(NodalCoordinates(I, 1) >= 0.0, 'GOM file: For axysimmetric analysis the x-coordinate (radial direction) has to be zero or positive.' )  
				if ( IS3DCYLINDRIC ) call Assert(NodalCoordinates(I, 1) >= 0.0, 'GOM file: For 3D Cylindrical analysis the x-coordinate (radial direction) has to be zero or positive.' )
              end do ! loop over all nodes
              NodalCoordinatesUpd = NodalCoordinates
              NodalOriginalCoord = NodalCoordinates

            else if (trim(BName) == '$$STARTELEMCON') then
              if (ELEMENTTYPE == TETRAOLD) then ! elementtype before v2018.2  
                do IDElement = 1, NElements ! loop over elements
                  read (GOMunit,*) (ElementConnectivities10Node(IDConnectivity, IDElement), IDConnectivity = 1, N_NODES_HOE) ! in the GOM-file always 10 nodes per element are written
                  do J = 1, ELEMENTNODES ! loop over element nodes
                    ElementConnectivities(J, IDElement) = ElementConnectivities10Node(J, IDElement)
                  end do ! loop over element nodes
                end do  ! loop over elements
              else ! elementtypes since v2018.2
                do IDElement = 1, NElements ! loop over elements
                  read (GOMunit,*) (ElementConnectivities(IDConnectivity, IDElement), IDConnectivity = 1, ELEMENTNODES)
                end do ! loop over elements 
              end if  

            elseif (trim(BName) == '$$START_FIXITY_SURFACE_SOLID') then
              read(GOMunit,*) NFixNod		
              do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDisp(FixNodID,I) = 0d0 ! zero prescribed displacement
                  end if  
				end do
              end do ! loop over fixed nodes  

            elseif (trim(BName) == '$$START_FIXITY_LINE_SOLID') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDisp(FixNodID,I) = 0d0 ! zero prescribed displacement
                  end if
                end do  
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_FIXITY_POINT_SOLID') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (solid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDisp(FixNodID,I) = 0d0 ! zero prescribed displacement
                  end if
                end do
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_REMOVE_FIXITY_SURFACE_SOLID') then
              read(GOMunit,*) NRemFixNod
              do IFxN = 1, NRemFixNod  ! loop over removed fixed nodes (solid)
                read(GOMunit,*) FixRemNodID, (NodRemBC(J), J = 1, NVECTOR)
                if(CalParams%ApplyRemoveSolidFixities) then
                  do I = 1, NVECTOR ! loop over vector size
                    if (NodRemBC(I) == -1) then
                      NodalPrescibedDisp(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
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
                      NodalPrescibedDisp(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
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
                      NodalPrescibedDisp(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
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
                    NodalPrescibedDispWater(FixNodID,I) = 0d0 ! zero prescribed displacement
                  end if
                end do
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_FIXITY_LINE_LIQUID') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (liquid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over degrees of freedom
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispWater(FixNodID,I) = 0d0 ! zero prescribed displacement
                  end if
                end do
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_FIXITY_POINT_LIQUID') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (liquid)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispWater(FixNodID,I) = 0d0 ! zero prescribed displacement
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
                      NodalPrescibedDispWater(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
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
                      NodalPrescibedDispWater(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
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
                      NodalPrescibedDispWater(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
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
                    NodalPrescibedDispGas(FixNodID,I) = 0d0 ! zero prescribed displacement
                  end if
                end do  
              end do ! loop over fixed nodes 

            elseif (trim(BName) == '$$START_FIXITY_LINE_GAS') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (gas)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over degrees of freedom
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispGas(FixNodID,I) = 0d0 ! zero prescribed displacement
                  end if
                end do  
              end do ! loop over fixed nodes

            elseif (trim(BName) == '$$START_FIXITY_POINT_GAS') then
              read(GOMunit,*) NFixNod
              do IFxN = 1, NFixNod  ! loop over fixed nodes (gas)
                read(GOMunit,*) FixNodID, (NodFixBC(J), J = 1, NVECTOR)
                do I = 1, NVECTOR ! loop over vector size
                  if (NodFixBC(I) == 1) then 
                    NodalPrescibedDispGas(FixNodID,I) = 0d0 ! zero prescribed displacement
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
                      NodalPrescibedDispGas(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
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
                      NodalPrescibedDispGas(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
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
                      NodalPrescibedDispGas(FixRemNodID,I) = 1d10 ! remove the fixities in x direction
                    end if
                  end do  
                end if
              end do ! loop over removed fixed nodes

            else if (trim(BName)=='$$START_LOAD_ON_NODES_SOLID') then
              read(GOMunit,*) NLoadedElementSidesSolidNodes ! Loaded faces
              if (NLoadedElementSidesSolidNodes > 0) then
                if ( NDIM == 3 ) then
                  allocate(LoadOnNodesConnectivitiesSolid(N_BOUNDARY_NODES_HOE, NLoadedElementSidesSolidNodes), stat = IError)
                  allocate(LoadValuesOnNodesSolid(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesSolidNodes), stat = IError)
                  do L = 1, NLoadedElementSidesSolidNodes ! loop over loaded faces
                    read(GOMunit,*) LoadOnNodesConnectivitiesSolid(1,L), LoadOnNodesConnectivitiesSolid(3,L),  &
                                    LoadOnNodesConnectivitiesSolid(5,L), LoadOnNodesConnectivitiesSolid(2,L), &
                                    LoadOnNodesConnectivitiesSolid(4,L), LoadOnNodesConnectivitiesSolid(6,L), &
                                    LoadValuesOnNodesSolid(1,1,L), LoadValuesOnNodesSolid(1,2,L), LoadValuesOnNodesSolid(1,3,L), &
                                    LoadValuesOnNodesSolid(3,1,L), LoadValuesOnNodesSolid(3,2,L), LoadValuesOnNodesSolid(3,3,L), &
                                    LoadValuesOnNodesSolid(5,1,L), LoadValuesOnNodesSolid(5,2,L), LoadValuesOnNodesSolid(5,3,L), &
                                    LoadValuesOnNodesSolid(2,1,L), LoadValuesOnNodesSolid(2,2,L), LoadValuesOnNodesSolid(2,3,L), &
                                    LoadValuesOnNodesSolid(4,1,L), LoadValuesOnNodesSolid(4,2,L), LoadValuesOnNodesSolid(4,3,L), &
                                    LoadValuesOnNodesSolid(6,1,L), LoadValuesOnNodesSolid(6,2,L), LoadValuesOnNodesSolid(6,3,L)
                  end do ! loop over loaded faces
                else if ( NDIM == 2 ) then
                  allocate(LoadOnNodesConnectivitiesSolid(ELEMENTBOUNDARYNODES, NLoadedElementSidesSolidNodes), stat = IError)
                  allocate(LoadValuesOnNodesSolid(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesSolidNodes), stat = IError)
                  do L = 1, NLoadedElementSidesSolidNodes ! loop over loaded element sides
                    read(GOMunit,*) LoadOnNodesConnectivitiesSolid(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnNodesSolid(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
                  end do ! loop over loaded element sides  
                end if
              end if
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

            else if (trim(BName)=='$$START_LOAD_ON_MATERIAL_POINTS_SOLID') then
              read(GOMunit,*) NLoadedElementSidesSolidMatPoints ! Loaded faces
              if (NLoadedElementSidesSolidMatPoints > 0) then
                if (.not.IsMPMComputation()) then
                  call GiveWarning('A load on material points (solid) is specified while not using MPM.')
                end if
				if ( NDIM == 3 ) then
                	allocate(LoadOnMatPointsConnectivitiesSolid(N_BOUNDARY_NODES_HOE, NLoadedElementSidesSolidMatPoints), stat = IError)
                	allocate(LoadValuesOnMatPointsSolid(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesSolidMatPoints), stat = IError)
                	do L = 1, NLoadedElementSidesSolidMatPoints ! loop over loaded faces
                  	read(GOMunit,*) LoadOnMatPointsConnectivitiesSolid(1,L), LoadOnMatPointsConnectivitiesSolid(3,L),  &
                                  	LoadOnMatPointsConnectivitiesSolid(5,L), LoadOnMatPointsConnectivitiesSolid(2,L), &
                                  	LoadOnMatPointsConnectivitiesSolid(4,L), LoadOnMatPointsConnectivitiesSolid(6,L), &
                            		LoadValuesOnMatPointsSolid(1,1,L), LoadValuesOnMatPointsSolid(1,2,L), LoadValuesOnMatPointsSolid(1,3,L), &
                            		LoadValuesOnMatPointsSolid(3,1,L), LoadValuesOnMatPointsSolid(3,2,L), LoadValuesOnMatPointsSolid(3,3,L), &
                            		LoadValuesOnMatPointsSolid(5,1,L), LoadValuesOnMatPointsSolid(5,2,L), LoadValuesOnMatPointsSolid(5,3,L), &
                            		LoadValuesOnMatPointsSolid(2,1,L), LoadValuesOnMatPointsSolid(2,2,L), LoadValuesOnMatPointsSolid(2,3,L), &
                            		LoadValuesOnMatPointsSolid(4,1,L), LoadValuesOnMatPointsSolid(4,2,L), LoadValuesOnMatPointsSolid(4,3,L), &
                            		LoadValuesOnMatPointsSolid(6,1,L), LoadValuesOnMatPointsSolid(6,2,L), LoadValuesOnMatPointsSolid(6,3,L)
                	end do ! loop over loaded faces
				else if ( NDIM == 2 ) then
                  allocate(LoadOnMatPointsConnectivitiesSolid(ELEMENTBOUNDARYNODES, NLoadedElementSidesSolidMatPoints), stat = IError)
                  allocate(LoadValuesOnMatPointsSolid(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesSolidMatPoints), stat = IError)
                  do L = 1, NLoadedElementSidesSolidMatPoints ! loop over loaded element sides
                    read(GOMunit,*) LoadOnMatPointsConnectivitiesSolid(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnMatPointsSolid(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
                  end do ! loop over loaded element sides
                end if
              end if
              
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

            else if ((trim(BName)=='$$START_LOAD_ON_NODES_LIQUID')) then
              read(GOMunit,*) NLoadedElementSidesWaterNodes ! Loaded faces
              if (NLoadedElementSidesWaterNodes > 0) then
                  if (NDIM == 3) then 
                    allocate(LoadOnNodesConnectivitiesWater(N_BOUNDARY_NODES_HOE, NLoadedElementSidesWaterNodes), stat = IError)
                    allocate(LoadValuesOnNodesWater(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesWaterNodes), stat = IError)
                    do L = 1, NLoadedElementSidesWaterNodes ! loop over loaded faces
                        read(GOMunit,*) LoadOnNodesConnectivitiesWater(1,L), LoadOnNodesConnectivitiesWater(3,L),  &
                                LoadOnNodesConnectivitiesWater(5,L), LoadOnNodesConnectivitiesWater(2,L), &
                                LoadOnNodesConnectivitiesWater(4,L), LoadOnNodesConnectivitiesWater(6,L), &
                                LoadValuesOnNodesWater(1,1,L), LoadValuesOnNodesWater(1,2,L), LoadValuesOnNodesWater(1,3,L), &
                                LoadValuesOnNodesWater(3,1,L), LoadValuesOnNodesWater(3,2,L), LoadValuesOnNodesWater(3,3,L), &
                                LoadValuesOnNodesWater(5,1,L), LoadValuesOnNodesWater(5,2,L), LoadValuesOnNodesWater(5,3,L), &
                                LoadValuesOnNodesWater(2,1,L), LoadValuesOnNodesWater(2,2,L), LoadValuesOnNodesWater(2,3,L), &
                                LoadValuesOnNodesWater(4,1,L), LoadValuesOnNodesWater(4,2,L), LoadValuesOnNodesWater(4,3,L), &
                                LoadValuesOnNodesWater(6,1,L), LoadValuesOnNodesWater(6,2,L), LoadValuesOnNodesWater(6,3,L)
                    end do ! loop over loaded faces
                 else if ( NDIM == 2 ) then
                    allocate(LoadOnNodesConnectivitiesWater(ELEMENTBOUNDARYNODES, NLoadedElementSidesWaterNodes), stat = IError)
                    allocate(LoadValuesOnNodesWater(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesWaterNodes), stat = IError)
                    do L = 1, NLoadedElementSidesWaterNodes ! loop over loaded element sides
                      read(GOMunit,*) LoadOnNodesConnectivitiesWater(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnNodesWater(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
                    end do ! loop over loaded element sides  
                end if
              end if
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

            else if ((trim(BName)=='$$START_LOAD_ON_MATERIAL_POINTS_LIQUID')) then
              read(GOMunit,*) NLoadedElementSidesWaterMatPoints ! Loaded faces
              if (NLoadedElementSidesWaterMatPoints > 0) then
                  
                if (.not.IsMPMComputation()) then
                  call GiveWarning('A load on material points (liquid) is specified while not using MPM')
                end if
              if (NDIM == 3) then 
                allocate(LoadOnMatPointsConnectivitiesWater(N_BOUNDARY_NODES_HOE, NLoadedElementSidesWaterMatPoints), stat = IError)
                allocate(LoadValuesOnMatPointsWater(N_BOUNDARY_NODES_HOE, NDOFL, NLoadedElementSidesWaterMatPoints), stat = IError)
                do L = 1, NLoadedElementSidesWaterMatPoints ! loop over loaded faces
                  read(GOMunit,*) LoadOnMatPointsConnectivitiesWater(1,L), LoadOnMatPointsConnectivitiesWater(3,L),  &
                                  LoadOnMatPointsConnectivitiesWater(5,L), LoadOnMatPointsConnectivitiesWater(2,L), &
                                  LoadOnMatPointsConnectivitiesWater(4,L), LoadOnMatPointsConnectivitiesWater(6,L), &
                            LoadValuesOnMatPointsWater(1,1,L), LoadValuesOnMatPointsWater(1,2,L), LoadValuesOnMatPointsWater(1,3,L), &
                            LoadValuesOnMatPointsWater(3,1,L), LoadValuesOnMatPointsWater(3,2,L), LoadValuesOnMatPointsWater(3,3,L), &
                            LoadValuesOnMatPointsWater(5,1,L), LoadValuesOnMatPointsWater(5,2,L), LoadValuesOnMatPointsWater(5,3,L), &
                            LoadValuesOnMatPointsWater(2,1,L), LoadValuesOnMatPointsWater(2,2,L), LoadValuesOnMatPointsWater(2,3,L), &
                            LoadValuesOnMatPointsWater(4,1,L), LoadValuesOnMatPointsWater(4,2,L), LoadValuesOnMatPointsWater(4,3,L), &
                            LoadValuesOnMatPointsWater(6,1,L), LoadValuesOnMatPointsWater(6,2,L), LoadValuesOnMatPointsWater(6,3,L)
                end do ! loop over loaded faces
              else if ( NDIM == 2 ) then
                allocate(LoadOnMatPointsConnectivitiesWater(ELEMENTBOUNDARYNODES, NLoadedElementSidesWaterMatPoints), stat = IError)
                allocate(LoadValuesOnMatPointsWater(ELEMENTBOUNDARYNODES, NDOFL, NLoadedElementSidesWaterMatPoints), stat = IError)
                  do L = 1, NLoadedElementSidesWaterMatPoints ! loop over loaded element sides
                    read(GOMunit,*) LoadOnMatPointsConnectivitiesWater(1:ELEMENTBOUNDARYNODES, L), ( LoadValuesOnMatPointsWater(I, 1:NDOFL, L), I = 1, ELEMENTBOUNDARYNODES )
                  end do ! loop over loaded element sides
              end if
              end if
              
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
              do IDElement = 1, NElements ! loop over elements
                read (GOMunit,*) ElementMaterialID(IDElement) ! get the material ID of the element
                if (ElementMaterialID(IDElement)==0) then
                  ElementMaterialID(IDElement) = -1 ! set the initially deactivated elements
                end if
              end do ! loop over elements  
              
            else if (trim(BName) == '$$FINISH') then
                    EXIT
            end if 
            
          end do 
          
          close(GOMunit)

          call InitCountersFromFile(CalParams%ApplyContactAlgorithm, &
                                    NElements, NNodes, &
                                    CalParams%NumberOfMaterials, &
                                    NLoadedElementSidesSolidNodes,  &
                                    NLoadedElementSidesSolidNodesB,  &
                                    NLoadedElementSidesWaterNodes, &
                                    NLoadedElementSidesWaterNodesB, &
                                    NLoadedElementSidesGasNodes, &
                                    NLoadedElementSidesGasNodesB, &
                                    NLoadedElementSidesSolidMatPoints,  &
                                    NLoadedElementSidesSolidMatPointsB,  &
                                    NLoadedElementSidesWaterMatPoints, &
                                    NLoadedElementSidesWaterMatPointsB, &
                                    NLoadedElementSidesGasMatPoints, &
                                    NLoadedElementSidesGasMatPointsB, &
                                    SoilSurfaceNumberofSides, &
                                    PhreaticSurfaceNumberofSides, &
                                    CalParams%NumberOfMaterials)
        
          call InitialiseNumberOfActiveElements()
          
          call InitialiseEdgeNodeTyings()
      
        end subroutine InitialiseMeshData


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

          allocate(EdgeNodeTyingsHOE(2, Counters%NodTot), stat = IError)
          EdgeNodeTyingsHOE = 0
       
          do I = 1, Counters%NodTot
            EdgeNodeTyingsHOE(1, I) = I
            EdgeNodeTyingsHOE(2, I) = I
          end do
          
          if (ELEMENTTYPE == TETRAOLD) then ! dimension 10 as only for 10-noded tetrahedral element
            do I = 1, Counters%NEl
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
          integer(INTEGER_TYPE), dimension(NDOFL*Counters%NodTot) :: NTreat
          integer(INTEGER_TYPE), dimension(NDOFL*Counters%NodTot) :: NTreatWater
          integer(INTEGER_TYPE), dimension(NDOFL*Counters%NodTot) :: NTreatGas
          integer(INTEGER_TYPE) :: IError
          
          NTreat = 1
          NTreatWater = 1
          NTreatGas = 1
          
          call InitialiseReduceDOF(NTreat) ! Initialises Counters%N
          call InitialiseReduceDofQuasiStatic()
          
          call FillPBoundary(NTreat)
          call FillNTreatWater(NTreatWater)
          call FillNPBoundaryWater(NTreatWater)
          call FillPBoundaryQuasiStatic()

          if (CalParams%ApplyAbsorbingBoundary) then
            PBoundaryWater = PBoundary
          end if
        
          call FillNTreatGas(NTreatGas)
          call FillNPBoundaryGas(NTreatGas)
      
          allocate(RateVolStrain(Counters%NEl), stat = IError)
          RateVolStrain = 0.0
          
        end subroutine InitialiseDerivedMeshData

        subroutine DetermineActiveNodes()
        !**********************************************************************
        !
        !    Function:  Determine which nodes belong to an active element
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IActiveElement, ElementID, INode, NodeID
          
          ActiveNode = .false.
          
          do IActiveElement = 1, Counters%NAEl
            ElementID = ActiveElement(IActiveElement)
            do INode = 1, ELEMENTNODES ! loop on element nodes
              NodeID = ElementConnectivities(INode, ElementID)
              ActiveNode(NodeID) = .true.
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

          do IEl = 1, Counters%NEl
            ISet = ElementMaterialID(IEl)
            ! Initialise activation status of elements
            if (ISet < 0) then
              IsActiveElement(IEl) = .false.
            else
              IsActiveElement(IEl) = .true.
              if (.not.IsFollowUpPhase()) then ! Else, this counter is read from file
                Counters%NAEl = Counters%NAEl + 1
              end if
            end if
          end do

        end subroutine InitialiseNumberOfActiveElements

        subroutine setElementDeterminant()
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

        if (.not.CalParams%AutomaticSkipConvection) RETURN

        if (.not.allocated(ElementDeterminant)) then
          iErr = 0
          allocate(ElementDeterminant(Counters%NAEl), stat = iErr)
          call AllocationError(IErr, 'ElementDeterminant', 'setElementDeterminant')
        else
          if (size(ElementDeterminant) /= Counters%NAEl) then
            iErr = 0
            deallocate(ElementDeterminant, stat = iErr)
            call DeAllocationError(IErr, 'ElementDeterminant', 'setElementDeterminant')
            allocate(ElementDeterminant(Counters%NAEl), stat = iErr)
            call AllocationError(IErr, 'ElementDeterminant', 'setElementDeterminant')
          endif
        endif

        do iEl=1,Counters%NAEl ! loop active elements
          ElementDeterminant(iEl) = GetElementDeterminant(iEl, ElementConnectivities, NodalCoordinatesUpd)
        enddo

        end subroutine setElementDeterminant

        logical function HasDistortedElement() result(res)
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

        res = .false.
        if (.not.allocated(ElementDeterminant)) then
          call setElementDeterminant()
          RETURN
        endif

        MinDeterminantRatio = huge(MinDeterminantRatio)
        do iEl=1,Counters%NAEl
          if (ElementDeterminant(iEl) > 0) then
            DeterminantRatio = GetElementDeterminant(iEl, ElementConnectivities, NodalCoordinatesUpd) / ElementDeterminant(iEl)
            MinDeterminantRatio = min(MinDeterminantRatio, DeterminantRatio)
          endif
        enddo

        MinimumDeterminantRatioReached = MinDeterminantRatio
        res = MinDeterminantRatio < CalParams%MinimumDeterminantRatio

        if (res) then
          call WriteInLogFile('Distortion detected at time step: ' // trim(String(CalParams%TimeStep)),  FEEDBACK_LEVEL_ALWAYS)
          call WriteInLogFile('Minimum determinant ratio reached: ' // trim(String(MinDeterminantRatio)), FEEDBACK_LEVEL_ALWAYS)
        endif


        end function HasDistortedElement

        subroutine SetIsDistorted(IsDistortedInitialised)
        !**********************************************************************
        !
        !    Function:  Sets a global variable to true when an element is distorted
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none
        logical, intent(in), optional :: IsDistortedInitialised

        if (.not.CalParams%AutomaticSkipConvection) RETURN
        if (present(IsDistortedInitialised)) then
          IsDistorted = IsDistortedInitialised
        else
          IsDistorted = HasDistortedElement()
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

          allocate(PBoundary(Counters%N), stat = IError)
          PBoundary = 1.0 ! Set to 1
          UMax = 0
          do I = 1, Counters%NodTot ! loop all nodes
            do J = 1, NDOFL ! loop on element dof's
              UMax = max(UMax, NodalPrescibedDisp(i,j))
            end do
          end do

          if (UMax > 9.9d4) then
            PduFree = 9.9d4
          else
            PduFree = 19.9
          end if

          do I = 1, Counters%NodTot ! loop all nodes
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDisp(I,J) < PduFree) then ! Prescribed displacement detected
                K = ReducedDof(I) + J
                PBoundary(K) = NodalPrescibedDisp(I,J)
              end if
            end do
          end do

          do I = 1, Counters%N ! loop all dof's
            if (NTreat(I) == 0) then
              PBoundary(I) = 0.0
            end if
          end do

        end subroutine FillPBoundary
        
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

          allocate(PBoundaryQuasiStatic(Counters%NodTot * (NDOFL + 1)), stat = IError)
        
          PBoundaryQuasiStatic = 0
          do NodeID = 1, Counters%NodTot
            !if (IsElementCornerNode(NodeID)) then
            do I = 1, NDIM
                IDoF = ReducedDoF(NodeID) + I
                IDoFQuasiStatic = (NodeID-1) * (NDOFL + 1) + I
                PBoundaryQuasiStatic(IDofQuasiStatic) = PBoundary(IDoF)
            end do
            IDoF = ReducedDoF(NodeID) + 1
            IDoFQuasiStatic = NodeID * (NDOFL + 1)
            if (PBoundaryWater(IDoF) /= 0) then
                PBoundaryQuasiStatic(IDofQuasiStatic) = 1 ! 1: natural
            else
                PBoundaryQuasiStatic(IDofQuasiStatic) = 0 ! 0: essential
            end if
            !endif
          end do

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

          allocate(PBoundaryWater(Counters%N), stat = IError)
         
          PBoundaryWater = 1.0 ! Set to 1
          UMax = 0
          do I = 1, Counters%NodTot ! loop all nodes
            do J = 1, NDOFL ! loop on element dof's
              UMax = max(UMax, NodalPrescibedDispWater(i,j))
            end do
          end do

          if (UMax > 9.9d4) then
            PduFree = 9.9d4
          else
            PduFree = 19.9
          end if

          do I = 1, Counters%NodTot ! loop all nodes
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDispWater(I,J) < PduFree) then ! Prescribed displacement detected
                K = ReducedDof(I) + J
                PBoundaryWater(K) = NodalPrescibedDispWater(I,J)
              end if
            end do
          end do
 
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

          allocate(PBoundaryGas(Counters%N), stat = IError)
         
          PBoundaryGas = 1.0 ! Set to 1
          UMax = 0
          do I = 1, Counters%NodTot ! loop all nodes
            do J = 1, NDOFL ! loop on element dof's
              UMax = max(UMax, NodalPrescibedDispGas(i,j))
            end do
          end do

          if (UMax > 9.9d4) then
            PduFree = 9.9d4
          else
            PduFree = 19.9
          end if

          do I = 1, Counters%NodTot ! loop all nodes
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDispGas(I,J) < PduFree) then ! Prescribed displacement detected
                K = ReducedDof(I) + J
                PBoundaryGas(K) = NodalPrescibedDispGas(I,J)
              end if
            end do
          end do
 
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
          
          allocate(ReducedDof(Counters%NodTot + 1), stat = IError)
          allocate(NodAct(Counters%NodTot), stat = IError)

          NodAct = 1
          do I = 1, Counters%NodTot
            aux = 0.0
            do J = 1, NDOFL
              aux = aux + dabs(NodalPrescibedDisp(I, J))
            end do
            if (aux == 0) NodAct(I) = 0
          end do
          
          ReducedDof =-1
          J = 0
          Do I = 1, Counters%NodTot
            if (ReducedDof(I) == -1) then
              ReducedDof(I) = J
              J = J + NDOFL
            end If
          end do
          ReducedDof(Counters%NodTot + 1) = J
          
          Counters%N = ReducedDof(Counters%NodTot + 1) ! this is the total number of dof's
          
          NTreat = 1
          do I = 1, Counters%NodTot
            do J = 1, NDOFL
              if (NodalPrescibedDisp(I, J) == 0) then
                IDof = ReducedDof(I) + J
                NTreat(IDof) = 0
              end If
            end do
          end do
          
          do I = 1, Counters%NodTot
            if (NodAct(I) == 0) then
              do J = 1, NDOFL
                NTreat(ReducedDof(I) + J) = 0
              end do
            end if
          end do
          
          deallocate(NodAct, stat = IError)
      
        end subroutine InitialiseReduceDOF
        
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

        allocate(ReducedDofQuasiStatic(Counters%NodTot + 1), stat = IError)
        
        ReducedDofQuasiStatic = -1
        J = 0
        do I = 1, Counters%NodTot
          ReducedDofQuasiStatic(I) = J
          J = J + NDOFL + 1
        end do
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
          
          allocate(NodAct(Counters%NodTot), stat = IError)
 
          NodAct = 1
          do I = 1, Counters%NodTot
            aux = 0.0
            do J = 1, NDOFL
              aux = aux + dabs(NodalPrescibedDispWater(I, J))
            end do
            if (aux == 0) NodAct(I) = 0
          end do  
                   
          NTreatWater = 1
          do I = 1, Counters%NodTot
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDispWater(I, J) == 0) then
                NTreatWater(ReducedDof(I) + J) = 0
              end If
            end do
          end do
          
          do I = 1, Counters%NodTot
            if (NodAct(I) == 0) then
              do J = 1, NDOFL ! loop on element dof's
                NTreatWater(ReducedDof(I) + J) = 0
              end do
            end if
          end do
          
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
          
          allocate(NodAct(Counters%NodTot), stat = IError)

          NodAct = 1
          do I = 1, Counters%NodTot
            aux = 0.0
            do J = 1, NDOFL
              aux = aux + dabs(NodalPrescibedDispGas(I, J))
            end do
            if (aux == 0) NodAct(I) = 0
          end do  

          NTreatGas = 1
          do I = 1, Counters%NodTot
            do J = 1, NDOFL ! loop on element dof's
              if (NodalPrescibedDispGas(I, J) == 0) then
                NTreatGas(ReducedDof(I) + J) = 0
              end If
            end do
          end do
          
          do I = 1, Counters%NodTot
            if (NodAct(I) == 0) then
              do J = 1, NDOFL ! loop on element dof's
                NTreatGas(ReducedDof(I) + J) = 0
              end do
            end if
          end do
          
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

          allocate(TotalDisplacementSoil(Counters%N), stat = IError)
          TotalDisplacementSoil = 0.0

          allocate(PhaseDisplacementSoil(Counters%N), stat = IError)
          PhaseDisplacementSoil = 0.0

          allocate(IncrementalDisplacementSoil(Counters%N, Counters%NEntity), stat = IError)
          IncrementalDisplacementSoil = 0.0

          allocate(AccumulatedIncDisplacementSoil(Counters%N, Counters%NEntity), stat = IError)
          AccumulatedIncDisplacementSoil = 0.0

          allocate(AccumulatedDisplacementSoil(Counters%N), stat = IError)
          AccumulatedDisplacementSoil = 0.0  
                  
          allocate(AccelerationSoil(Counters%N, Counters%NEntity), stat = IError)
          AccelerationSoil = 0.0

          allocate(TotalVelocitySoil(Counters%N, Counters%NEntity), stat = IError)
          TotalVelocitySoil = 0.0
          
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
          
          allocate(NodalUnitMassGradient(Counters%NodTot, NDIM), stat = IError)
          NodalUnitMassGradient = 0.0
          allocate(ExtFlow(Counters%NodTot), stat = IError)
          ExtFlow = 0.0

          allocate(RateOfFlux(Counters%NodTot), stat = IError)
          RateOfFlux = 0.0
          
          allocate(IntFlow(Counters%NodTot), stat = IError)
          IntFlow = 0.0
          
          allocate(IncrementalPressure(Counters%NodTot), stat = IError)
          IncrementalPressure = 0.0

          allocate(TotalPressure(Counters%NodTot), stat = IError)
          TotalPressure = 0.0
          
          allocate(SubIncrementalDisplacement(Counters%N), stat = IError)
          SubIncrementalDisplacement = 0.0
      
          allocate(SubIncrementalPressure(Counters%NodTot), stat = IError)
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
          integer(INTEGER_TYPE) :: I, J, SurfaceID, NSurfaces, IError, ios
          integer(INTEGER_TYPE) :: DumE, DUmI, NumSideNodes
          integer(INTEGER_TYPE), dimension(:), allocatable :: DumN
          character(len=255) :: DumS
          character(len=255) :: BName
          character(len=255), dimension(MAXOUTPUTSURFACES) :: DumSurfaceName
          
          allocate(IsReactionNode(Counters%NodTot), stat = Ierror)
          IsReactionNode = .false.
        
          
          if (ELEMENTTYPE==TETRAOLD) then
            NumSideNodes = N_BOUNDARY_NODES_HOE
          else
            NumSideNodes = ELEMENTBOUNDARYNODES
          end if
          allocate(DumN(NumSideNodes), Stat = Ierror)
          DumN = -1
          
          
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
          if (trim(BName) == '$$START_OUTPUT_REACTION_FORCES') then
           read(GOMUnit, *) DumI
           if (DumI < 1) EXIT
           do J=1,DumI
            read(GOMunit, *) DumS, DumE, DumN(1:NumSideNodes)
            !DumSurfaceName(J) = DumS
            !do while (trim(DumS) /= '$$END_OUTPUT_REACTION_FORCES')
            do I=1,MAXOUTPUTSURFACES
              if (DumSurfaceName(I) == DumS) then !the surface has already been specified before
                 EXIT
              else if (DumSurfaceName(I) == '/') then !first non-assigned position of the array
                DumSurfaceName(I) = DumS !Assign the surface name to this position 
                NSurfaces = NSurfaces + 1 !increment the number of surfaces
                EXIT
              end if
             end do
              !read(GOMunit, *) DumE, DumN(1:6)
              IsReactionNode(DumN(1:NumSideNodes)) = .true.
              !read(GOMUnit,*) DumS, DumE, DumN(1:6)
           end do
           EXIT
          else if (trim(BName) == '$$FINISH') then
           EXIT 
          end if
          end do
          if (NSurfaces < 1) RETURN
          rewind(GOMUnit)
          Counters%NReactionSurfaceOutput = NSurfaces
          
          !allocate vectors
          allocate(IsReactionNodeSurface(Counters%NodTot,NSurfaces), stat = Ierror)
          allocate(OutputSurfaceName(NSurfaces), stat = IError)
          
          OutputSurfaceName(1:NSurfaces) = DumSurfaceName(1:NSurfaces)
          IsReactionNodeSurface = .false.
          
          !**************************************************************
          !Second read: Assign the nodes to the surfaces
          do
            read(GOMunit, '(A)', iostat=ios) BName 
            if (trim(BName) == '$$START_OUTPUT_REACTION_FORCES') then 
              read(GOMUnit, *) DumI
              do J=1,DumI
                read(GOMunit, *) DumS, DumE, DumN(1:NumSideNodes)
                do I=1,NSurfaces
                 if (OutputSurfaceName(I) == DumS) then !the surface has already been specified before
                  SurfaceID = I !get surfaceID
                  IsReactionNodeSurface(DumN(1:NumSideNodes),SurfaceID) = .true.
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
          if (allocated(DumN)) then
            deallocate (DumN, stat = IError)
          end if
          
          
          !call DetermineReactionElements()        
          
      end subroutine InitialiseSurfaceReaction
    
        
      subroutine InitialiseGPGlobalPositionArrays()
      !*********************************************************************************
      !
      ! Function: Determine global position of Gauss points for axisymmetric.
      !
      !**********************************************************************************
      implicit none
      
          integer(INTEGER_TYPE) :: iError

          if ( .not. ISAXISYMMETRIC ) RETURN
          
          if ( allocated(GPGlobalPositionElement) ) then
            deallocate( GPGlobalPositionElement, stat=iError )
            call DeAllocationError(iError, 'GPGlobalPositionElement', 'initialiseGPGlobalPositionArrays')
          end if

          allocate( GPGlobalPositionElement(NDIM, ELEMENTGAUSSPOINTS, Counters%NEl), stat=iError )
          call AllocationError(iError, 'GPGlobalPositionElement', 'initialiseGPGlobalPositionArrays')
          call setGPGlobalPositionElement(GPGlobalPositionElement)

      end subroutine initialiseGPGlobalPositionArrays

      
      subroutine setGPGlobalPositionElement(GPGlobalPositionElement)
      
          implicit none
      
          real(REAL_TYPE), dimension(NDIM, ELEMENTGAUSSPOINTS, Counters%NEl), intent(out):: GPGlobalPositionElement

          ! local variables
          integer(INTEGER_TYPE) :: IElement, IGaussPoint
          real(REAL_TYPE), dimension(NDIM, ELEMENTGAUSSPOINTS) :: PosGP
          real(REAL_TYPE), dimension(ELEMENTGAUSSPOINTS) :: WeiGP
          real(REAL_TYPE), dimension(ELEMENTNODES) :: GaussPointShapeValues
          real(REAL_TYPE), dimension(ELEMENTNODES, NDIM, ELEMENTGAUSSPOINTS) :: DShapeValues
          real(REAL_TYPE), dimension(NDIM) :: LocalPosition
      
          GPGlobalPositionElement = 0.0
          LocalPosition = 0.0

          do IGaussPoint = 1, ELEMENTGAUSSPOINTS
            ! Detemine local poisiton of the GaussPoint
            call GaussPointLocalCoordinates(IGaussPoint, WeiGP(IGaussPoint), PosGP(:,IGaussPoint))
            ! Detemine the shape function of the GaussPoint
            call ShapeFunctionData(PosGP(:,IGaussPoint), ELEMENTNODES, GaussPointShapeValues, DShapeValues(:,:,IGaussPoint))
          end do

          do IElement = 1, Counters%NEl
            do IGaussPoint = 1, ELEMENTGAUSSPOINTS
              ! Determine the global location of IGaussPoint
              call CoordGaussPointLocalToGlobal(IElement, NodalCoordinatesUpd, GaussPointShapeValues, LocalPosition) 
              GPGlobalPositionElement(:, IGaussPoint,IElement) = LocalPosition
            end do
          end do

      end subroutine setGPGlobalPositionElement
      
      
      subroutine CoordGaussPointLocalToGlobal(IElement, LNodalCoordinates, GaussPointShapeValues, GaussPointGlobalCoord)
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

        ! Local variables
        integer(INTEGER_TYPE) :: NodeID
        integer(INTEGER_TYPE) :: I, IDim
          
        ! initialize global coordinate of Gauss Point
        GaussPointGlobalCoord = 0.0

        do I = 1, ELEMENTNODES ! Loop over nodes of IElement
          NodeID = iabs(ElementConnectivities(I, IElement) )
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
                        CoordSurfNodes(I, :) = (NodalCoordinates(NodeID,:))
                    else
                        if (Connectiv(J,I) /= Connectiv(J-1,I)) then
                            CoordSurfNodes(SidesCount+1,:) = (NodalCoordinates(NodeID,:))
                        end if
                    end if
                end do
            end do
        else if (NDIM == 3) then
            do J = 1, NDIM
                do I = 1, SidesCount
                NodeID = Connectiv(J,I) !
                    if (J == 1) then
                        CoordSurfNodes(I, :) = (NodalCoordinates(NodeID,:))
                    else if (J == 2) then
                        if (Connectiv(J,I) /= Connectiv(J-1,I)) then
                            CoordSurfNodes(SidesCount+1,:) = (NodalCoordinates(NodeID,:))
                        end if
                    else
                        if ((Connectiv(J,I) /= Connectiv(J-1,I)) .and. (Connectiv(J,I) /= Connectiv(J-2,I))) then
                            CoordSurfNodes(SidesCount+2,:) = (NodalCoordinates(NodeID,:))
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
        
      end module ModMeshInfo
