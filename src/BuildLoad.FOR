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
	  
	  
	  ! Module BuildLoad
      !**********************************************************************
      !
      !     $Revision: 8842 $
      !     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
      !
      !**********************************************************************
      
    
      subroutine InitialiseTractionLoad()
      !**********************************************************************
      !
      !>    Function:  Initialisation of the traction load for all phases.
      !
      !**********************************************************************

      use ModRotBoundCond
      use ModCounters
      use ModGlobalConstants

      
      implicit none

        if (Counters%NLoadedElementSidesSolid > 0) then ! solid
          call InitialiseTractionLoadVector(LOADTYPE_SOLID)
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVector, ExtLoadVector) ! 3D function, only if IS3DCYLINDRIC=.true.
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorB, ExtLoadVectorB) ! 3D function, only if IS3DCYLINDRIC=.true.
        end if
        
        if (((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) .or. (Counters%NLoadedElementSidesWater > 0)) then ! liquid
        !just temporary: need to find a different condition
          call InitialiseTractionLoadVector(LOADTYPE_LIQUID)  
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorwater, ExtLoadVectorWater) ! 3D function, only if IS3DCYLINDRIC=.true.
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorwaterB, ExtLoadVectorWaterB) ! 3D function, only if IS3DCYLINDRIC=.true.
          ! note: the rotation of the BC for the liquid has to be fixed
        end if
        
        if (Counters%NLoadedElementSidesGas > 0) then ! gas
          call InitialiseTractionLoadVector(LOADTYPE_GAS)  
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorGas, ExtLoadVectorGas) ! 3D function, only if IS3DCYLINDRIC=.true.
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorGasB, ExtLoadVectorGasB) ! 3D function, only if IS3DCYLINDRIC=.true.
          ! note: the rotation of the BC for the gas has to be fixed
        end if

      end subroutine InitialiseTractionLoad
      
      
      subroutine InitialiseTractionLoadVector(LoadType)
      !**********************************************************************
      !
      !>    Function:  Initialisation of the traction load vector, splitting  
      !!               the load on nodes and the load on material points.
      !
      !**********************************************************************

      use ModReadCalculationData
      use ModMPMInit
      use ModMeshInfo
      use ModGlobalConstants
      use ModDYNConvectivePhase      
      
      
      implicit none
      
        ! arguments
      integer(INTEGER_TYPE), intent(in) :: LoadType

      ! Local variables
      integer(INTEGER_TYPE) :: I, LoadedSides

      if ( LoadType == LOADTYPE_SOLID ) then ! solid

          ! initialise ExtLoadVector with distributed load on nodes
          if (NDIM == 3) then ! for 3D
            do I = 1, Counters%NLoadedElementSidesSolidNodes
              call Load3D(ExtLoadVector, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType)
            end do
            do I = 1, Counters%NLoadedElementSidesSolidNodesB
              call Load3D(ExtLoadVectorB, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesSolidB(1, I), LoadValuesOnNodesSolidB(1, 1, I), LoadType)
            end do
          elseif (NDIM == 2) then ! for 2D
            do I = 1, Counters%NLoadedElementSidesSolidNodes
              call Load2D(ExtLoadVector, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType)              
            end do
            do I = 1, Counters%NLoadedElementSidesSolidNodesB
              call Load2D(ExtLoadVectorB, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesSolidB(1, I), LoadValuesOnNodesSolidB(1, 1, I), LoadType)              
            end do
          end if  
        
          ! initialise Particle%ExtLoadVector with distributed load on material points (only in first load step)
          if (.not.IsFollowUpPhase()) then
            if (NDIM == 3) then ! for 3D only, sprint #2: add 2D functionality
              do I = 1, Counters%NLoadedElementSidesSolidMatPoints
                call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesSolid(1, I), LoadValuesOnMatPointsSolid(1, 1, I), 1, NodalCoordinates,1)
              end do
              do I = 1, Counters%NLoadedElementSidesSolidMatPointsB
                call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesSolidB(1, I), LoadValuesOnMatPointsSolidB(1, 1, I), 1, NodalCoordinates,2)
			  end do
            elseif (NDIM == 2) then ! for 2D
              do I = 1, Counters%NLoadedElementSidesSolidMatPoints
                call TransferExternalLoadsToParticles(LoadType, ELEMENTBOUNDARYNODES, ELEMENTNODES, LoadOnMatPointsConnectivitiesSolid(1, I), LoadValuesOnMatPointsSolid(1, 1, I), 1, NodalCoordinates,1)
              end do
              do I = 1, Counters%NLoadedElementSidesSolidMatPointsB
                call TransferExternalLoadsToParticles(LoadType, ELEMENTBOUNDARYNODES, ELEMENTNODES, LoadOnMatPointsConnectivitiesSolidB(1, I), LoadValuesOnMatPointsSolidB(1, 1, I), 1, NodalCoordinates,2)
              end do
            end if
          end if
              
      else if ( LoadType == LOADTYPE_LIQUID ) then ! liquid
          ! initialise ExtLoadVectorWater with distributed load on nodes
       ! Hydraulic head on NODES
          
          
          if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
              
              call GetAndUpdateHydraulicHeadLoad(HydraulicHeadLoadedElemID, LoadedSides, HydraulicHeadNodesConnectivities, HydraulicHeadLoad)
              Counters%HydraulicHeadSides = LoadedSides
              
              if (NDIM == 3) then
               do I = 1, Counters%HydraulicHeadSides
                call Load3D(HydraulicHeadVector, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, HydraulicHeadNodesConnectivities(1, I), HydraulicHeadLoad(1, 1, I), LoadType)
            end do
        elseif (NDIM == 2) then
            do I = 1, Counters%HydraulicHeadSides
                call Load2D(HydraulicHeadVector, ReducedDof, NodalCoordinates, 1, HydraulicHeadNodesConnectivities(1, I), HydraulicHeadLoad(1, 1, I), LoadType)
            end do
        end if
          end if
          ! general pressure load on NODES
          
          ! initialise ExtLoadVectorWater with distributed load on nodes
          if (NDIM == 3) then ! for 3D
              do I = 1, Counters%NLoadedElementSidesWaterNodes
                  call Load3D(ExtLoadVectorWater, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesWater(1, I), LoadValuesOnNodesWater(1, 1, I), LoadType)
              end do
              do I = 1, Counters%NLoadedElementSidesWaterNodesB
                  call Load3D(ExtLoadVectorWaterB, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesWaterB(1, I), LoadValuesOnNodesWaterB(1, 1, I), LoadType)
              end do
          else if (NDIM == 2) then ! for2D
              do I = 1, Counters%NLoadedElementSidesWaterNodes
                  call Load2D(ExtLoadVectorWater, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesWater(1, I), LoadValuesOnNodesWater(1, 1, I), LoadType)
              end do 
              do I = 1, Counters%NLoadedElementSidesWaterNodesB
                  call Load2D(ExtLoadVectorWaterB, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesWaterB(1, I), LoadValuesOnNodesWaterB(1, 1, I), LoadType)
              end do
          end if 
        
          ! initialise Particle%ExtLoadVectorWater with distributed load on material points (only in first load step)
          if (.not.IsFollowUpPhase()) then
            if (NDIM == 3) then ! for 3D
              do I = 1, Counters%NLoadedElementSidesWaterMatPoints
                call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesWater(1, I), LoadValuesOnMatPointsWater(1, 1, I), 1, NodalCoordinates,1)
              end do
              do I = 1, Counters%NLoadedElementSidesWaterMatPointsB
                call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesWaterB(1, I), LoadValuesOnMatPointsWaterB(1, 1, I), 1, NodalCoordinates,2)
              end do
            elseif (NDIM == 2) then ! for 2D
              do I = 1, Counters%NLoadedElementSidesWaterMatPoints
                call TransferExternalLoadsToParticles(LoadType, ELEMENTBOUNDARYNODES, ELEMENTNODES, LoadOnMatPointsConnectivitiesWater(1, I), LoadValuesOnMatPointsWater(1, 1, I), 1, NodalCoordinates,1)  
              end do
              do I = 1, Counters%NLoadedElementSidesWaterMatPointsB
                call TransferExternalLoadsToParticles(LoadType, ELEMENTBOUNDARYNODES, ELEMENTNODES, LoadOnMatPointsConnectivitiesWaterB(1, I), LoadValuesOnMatPointsWaterB(1, 1, I), 1, NodalCoordinates,2)  
              end do
            end if
          end if  

        else if ( LoadType == LOADTYPE_GAS ) then ! gas

          ! initialise ExtLoadVectorGas with distributed load on nodes
          if (NDIM == 3) then ! for 3D only, sprint #2: add 2D functionality
            do I = 1, Counters%NLoadedElementSidesGasNodes
              call Load3D(ExtLoadVectorGas, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesGas(1, I), LoadValuesOnNodesGas(1, 1, I), LoadType)
            end do
            do I = 1, Counters%NLoadedElementSidesGasNodesB
              call Load3D(ExtLoadVectorGasB, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesGasB(1, I), LoadValuesOnNodesGasB(1, 1, I), LoadType)
            end do
           
          else if (NDIM == 2) then ! for2D
              do I = 1, Counters%NLoadedElementSidesGasNodes
                  call Load2D(ExtLoadVectorGas, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesGas(1, I), LoadValuesOnNodesGas(1, 1, I), LoadType)
              end do 
              do I = 1, Counters%NLoadedElementSidesWaterNodesB
                  call Load2D(ExtLoadVectorGasB, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesGasB(1, I), LoadValuesOnNodesGasB(1, 1, I), LoadType)
              end do
          end if
        
          ! initialise Particle%ExtLoadVectorGas with distributed load on material points (only in first load step)
          if (.not.IsFollowUpPhase()) then
            if (NDIM == 3) then ! for 3D only, sprint #2: add 2D functionality
              do I = 1, Counters%NLoadedElementSidesGasMatPoints
                call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesGas(1, I), LoadValuesOnMatPointsGas(1, 1, I), 1, NodalCoordinates,1)
              end do
              do I = 1, Counters%NLoadedElementSidesGasMatPointsB
                call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesGasB(1, I), LoadValuesOnMatPointsGasB(1, 1, I), 1, NodalCoordinates,2)
              end do
            end if
          end if  
       
       end if ! LoadType
       
      end subroutine InitialiseTractionLoadVector
      
      
      subroutine Load2D(RLoad, NDof, Coord, NumberOfIntegrationPoints, ILoadCon, LoadValue, LoadType)
      !**********************************************************************
      !
      !>    Function:  Calculate external load vector for a 2-noded (low-order) 
      !>               or 3-noded (high-order) line element
      !
      !     Note : 2D function
      !
      !**********************************************************************
        use ModCounters
        use ModElementEvaluation
        use ModGlobalConstants
        use ModMeshAdjacencies
        use ModMeshInfo
      
        implicit none
      
          ! arguments
          integer(INTEGER_TYPE), intent(in) :: NumberOfIntegrationPoints, LoadType
          integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES), intent(in) :: ILoadCon ! size (ELEMENTBOUNDARYNODES)
          integer(INTEGER_TYPE), dimension(Counters%NodTot), intent(in) :: NDof ! size (Counters%NodTot)
          real(REAL_TYPE), dimension(ELEMENTBOUNDARYNODES, NVECTOR), intent(in) ::  LoadValue ! size (ELEMENTBOUNDARYNODES, NVECTOR)
          real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Coord ! size (Counters%NodTot, NVECTOR)      ! Coordinates vector
          real(REAL_TYPE), dimension(Counters%NodTot*NVECTOR), intent(inout) :: RLoad ! size(Counters%NodTot*NVECTOR)     ! Resulting load vector
      
          ! local variable
          integer(INTEGER_TYPE) :: I, J, K, NodeNumber, ND, NumberOfCornerNodes
          integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: ILoadConLocal
          real(REAL_TYPE) :: VectorLength, Temp, Factor, Radius, xLoad
          real(REAL_TYPE), dimension(NVECTOR) :: NormalVector, ProductTractionNodesNormalVector, Vector
          real(REAL_TYPE), dimension(NumberOfIntegrationPoints, NVECTOR) :: TractionGaussPoint
          real(REAL_TYPE), dimension(ELEMENTBOUNDARYNODES, NVECTOR) :: TractionNodes, SE
          real(REAL_TYPE), dimension(ELEMENTBOUNDARYNODES) :: AppliedPressure
          logical :: found = .false.
          
          integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, IElNode1, IElNode2, ElementOfNode1, ElementOfNode2, ISharedElement, INode
          real(REAL_TYPE) :: ScalarProduct 

          NumberOfCornerNodes = 2 ! number of corner nodes of a line element

          if (ELEMENTBOUNDARYNODES == 2) then ! 2-noded (low-order) line element
            ! no re-ordering needed
            ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
            ILoadConLocal(2) = ILoadCon(2); TractionNodes(2,:) = LoadValue(2,:)
          else if (ELEMENTBOUNDARYNODES == 3) then ! 3-noded (high-order) line element
            ! re-order nodes: corner nodes first, mid-nodes last
            ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
            ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
            ILoadConLocal(3) = ILoadCon(2); TractionNodes(3,:) = LoadValue(2,:)
          end if      

          ! set arrays to zero
          TractionGaussPoint = 0.0
          SE = 0.0

          do I = 1, NumberOfIntegrationPoints
            ! determine vector normal to a line (for 2- and 3-noded line), note that only the two corner nodes are passed
            call NormalOnLine(I, Coord, ILoadConLocal(1:2), GPShapeFunctionDerivativeBoundary, NormalVector, VectorLength)

            Factor = 1.0 ! Factor for length of lines

            if (LoadType == LOADTYPE_LIQUID .or. LoadType == LOADTYPE_GAS) then ! Liquid or Gas, Loading is applied only normal to the surface
                AppliedPressure(:) = LoadValue(:,1) !Pressure is hydrostatic, all components are equal LoadValue(:,1)=LoadValue(:,2)
                !*** check the normal pointing outward of the body **** 
                ! Note that this implementation works only for TRI elements. There is room for vectorization to extend to more types of element. 
                !Correct normal if necessary. LoadValue<0 means compression
                NElemNode1 = GetNElmOfNode(ILoadConLocal(1))
                NElemNode2 = GetNElmOfNode(ILoadConLocal(2))
                do IElNode1 = 1,NElemNode1
                    ElementOfNode1 = GetElmIOfNode(ILoadConLocal(1),IElNode1)
                    do IElNode2 = 1,NElemNode2
                        ElementOfNode2 = GetElmIOfNode(ILoadConLocal(2),IElNode2)
                        if (ElementOfNode1==ElementOfNode2) then !found element connected to line
                            ISharedElement = ElementOfNode2 !also = ElementOfNode1
                            if (IsActiveElement(ISharedElement)) then 
                             do K=1,ELEMENTNODES
                                 INode =ElementConnectivities(K,ISharedElement)
                                 if ((INode/=ILoadConLocal(1)).and.(INode/=ILoadConLocal(2))) then !this is the node not belonging to the line
                                     Vector(:) = Coord(INode,:) - Coord(ILoadConLocal(1),:)
                                     ScalarProduct = DotProduct(Vector,NormalVector,NVECTOR)
                                     if (ScalarProduct > 0) then ! The normal is pointing inside element, reverse
                                        NormalVector = -NormalVector !the normal is pointing outside active element
                                        EXIT !no need to check the other nodes of the element
                                     elseif (ScalarProduct < 0) then ! The normal is pointing outside element, keep direction
                                        EXIT !no need to check the other nodes of the element
                                     end if
                                 end if
                             end do
                             EXIT !No need to check the other elements because already found one which is active
                           end if
                        end if
                    end do
                end do
                
                
                if (ELEMENTBOUNDARYNODES == 3) then 
                    TractionNodes(3,1:NVECTOR) = 0.0 ! remove middle point data
                end if 
                
              ProductTractionNodesNormalVector = 0.0 ! initialize

              do K = 1, NumberOfCornerNodes ! for each corner node

                !do J = 1, NVECTOR ! for each direction
                !  ProductTractionNodesNormalVector(K) = ProductTractionNodesNormalVector(K) + TractionNodes(K,J) * NormalVector(J)
                !end do
                !
                !  do J = 1, NVECTOR
                !    TractionNodes(K, J) = ProductTractionNodesNormalVector(K) * NormalVector(J)
                !  end do
                !  ScalarProduct = DotProduct(TractionNodes(K,:),NormalVector,NVECTOR) !>0 if traction, <0 if compression, normal pointing outward
                !  if (sign(1.d0,ScalarProduct) /= sign(1.d0,LoadValue(K,1))) then !LoadValue(K,1)=LoadValue(K,2)>0 for traction
                !      TractionNodes(K, :) = - TractionNodes(K, :)
                !  end if
                  do J = 1, NVECTOR
                     TractionNodes(K, J) = AppliedPressure(K) * NormalVector(J)
                  end do   
    
              end do
            end if
            
            
            ! calculate traction at Gauss Point
            do J = 1, NVECTOR
              do K = 1, NumberOfCornerNodes ! Corner nodes
                if ( ISAXISYMMETRIC ) then
                  NodeNumber = ILoadConLocal(K)
                  Radius = Coord(NodeNumber, 1) ! index 1 is r-direction
                else
                  Radius = 1.0
                end if
                TractionGaussPoint(I,J) = TractionGaussPoint(I,J) + GPShapeFunctionBoundary(I,K) * TractionNodes(K,J) * Radius
              end do
            end do

            ! calculate traction at Nodes
            do J = 1, NVECTOR
              
              found = .false.
              do K = 1, NumberOfCornerNodes ! Corner nodes
                ! factor = 0.5 for triangle, 1.0 for ractangle
                Temp = GPWeightBoundary(I) * GPShapeFunctionBoundary(I, K) * VectorLength * Factor

                SE(K, J) = SE(K, J) + Temp * TractionGaussPoint(I, J)
                if ( ISAXISYMMETRIC ) then
                  NodeNumber = ILoadConLocal(K)
                  Radius = Coord(NodeNumber, 1) ! index 1 is r-direction
                  if ( Radius < TINY ) then
                    xLoad = SE(K, J)
                    found = .true.
                  end if
                end if
              end do
              
              if ( ISAXISYMMETRIC ) then
                if (found) then
                  xLoad = (xLoad / real((NumberOfCornerNodes -1), REAL_TYPE)) * (7.0 / 8.0)
                  do K = 1, NumberOfCornerNodes ! Corner nodes
                    NodeNumber = ILoadConLocal(k)
                    Radius = Coord(NodeNumber, 1) ! index 1 is r-direction
                    if ( Radius > TINY ) then
                      SE(K, J) = SE(K, J) + xLoad
                    else
                      SE(K, J) = SE(K, J) - xLoad
                    end if
                  end do
                end if
              end if
              
            end do
            
          end do

          ! add to global load vector
          do I = 1, ELEMENTBOUNDARYNODES
            ND = NDof(ILoadConLocal(I))
            do J = 1, NVECTOR
              RLoad(ND + J) = RLoad(ND + J) + SE(I,J)
            end do
          end do

      end subroutine Load2D
      
      
      subroutine Load3D(RLoad, NDof, Coord, NInt, NNod, ILoadCon, LoadValue, LoadType)
      !**********************************************************************
      !
      !>    Function:  Calculate external load vector for a 6-noded triangle 
      !!               or 8-noded rectangle
      !     Note : 3D function
      !
      !**********************************************************************
      
      use ModCounters
      use ModElementEvaluation
      use ModGlobalConstants
      use ModMeshAdjacencies
        use ModMeshInfo
      
      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: NInt, NNod, LoadType 
        integer(INTEGER_TYPE), dimension(NNod) :: ILoadCon
        integer(INTEGER_TYPE), dimension(Counters%NodTot) :: NDof
        real(REAL_TYPE), dimension(NNod, NVECTOR), intent(in) ::  LoadValue
        real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Coord       ! Coordinates vector
        real(REAL_TYPE), dimension(Counters%NodTot*NDOFL), intent(inout) :: RLoad      ! Resulting load vector
      
        ! local variable
        integer(INTEGER_TYPE) :: I, J, K, NN, ND, NNodLOE
        integer(INTEGER_TYPE), dimension(NNod) :: ILoadConLocal
        real(REAL_TYPE) :: VectorLength, Temp
        real(REAL_TYPE), dimension(NVECTOR) :: NormalVector, ProductTractionNodesNormalVector, Vector
        real(REAL_TYPE), dimension(NInt, NVECTOR) :: TractionGaussPoint
        real(REAL_TYPE), dimension(NNod, NVECTOR) :: TractionNodes, SE
        real(REAL_TYPE), dimension(NNod) :: AppliedPressure
        
         integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, NElemNode3, IElNode1, IElNode2, IElNode3, ElementOfNode1, ElementOfNode2, ElementOfNode3, ISharedElement, INode
         real(REAL_TYPE) :: ScalarProduct
      
        NNodLOE = NNod/2 ! number of nodes of a lower order element
      
        if (NNod == 6) Then
          ! corners first, mid-nodes later
          ! 6-noded triangle
          ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
          ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
          ILoadConLocal(3) = ILoadCon(5); TractionNodes(3,:) = LoadValue(5,:)
          ILoadConLocal(4) = ILoadCon(2); TractionNodes(4,:) = LoadValue(2,:)
          ILoadConLocal(5) = ILoadCon(4); TractionNodes(5,:) = LoadValue(4,:)
          ILoadConLocal(6) = ILoadCon(6); TractionNodes(6,:) = LoadValue(6,:)
        else
          ! corners first, mid-nodes later
          ! 8-noded rectangle
          ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
          ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
          ILoadConLocal(3) = ILoadCon(5); TractionNodes(3,:) = LoadValue(5,:)
          ILoadConLocal(4) = ILoadCon(7); TractionNodes(4,:) = LoadValue(7,:)
          ILoadConLocal(5) = ILoadCon(2); TractionNodes(5,:) = LoadValue(2,:)
          ILoadConLocal(6) = ILoadCon(4); TractionNodes(6,:) = LoadValue(4,:)
          ILoadConLocal(7) = ILoadCon(6); TractionNodes(7,:) = LoadValue(6,:)
          ILoadConLocal(8) = ILoadCon(8); TractionNodes(8,:) = LoadValue(8,:)
          !'no 8-noded loads yet'
        end if
        
        ! set arrays to zero
        TractionGaussPoint = 0d0
        SE = 0d0
      
        do I = 1, NInt
          ! determine vector normal to a plane (for 3-noded triangles)
          ! vector length = det(Jacobian_S)
          call Normal_T3( I, Coord, ILoadConLocal, NNodLOE, &
                          GPShapeFunctionDerivativeBoundary, NormalVector, VectorLength)
            
          if((LoadType == LOADTYPE_LIQUID) .or. (LoadType == LOADTYPE_GAS)) then ! Liquid or Gas, Loading is applied only normal to the surface
              AppliedPressure(:) = LoadValue(:,1) !Liquid pressure is hydrostatic, LoadValue(:,1)=LoadValue(:,2)...
           !*** se the normal pointing outward of the body **** 
              !IMPROVE IMPLEMENTATION!!! ONLY WORKS FOR TETRA
              
              NElemNode1 = GetNElmOfNode(ILoadConLocal(1))
              NElemNode2 = GetNElmOfNode(ILoadConLocal(2))
              NElemNode3 = GetNElmOfNode(ILoadConLocal(3))
               
              do IElNode1 = 1,NElemNode1
                ElementOfNode1 = GetElmIOfNode(ILoadConLocal(1),IElNode1)
                do IElNode2 = 1,NElemNode2
                  ElementOfNode2 = GetElmIOfNode(ILoadConLocal(2),IElNode2)
                    if (ElementOfNode1==ElementOfNode2) then !found element connected to line
                      do IElNode3 = 1,NElemNode3
                        ElementOfNode3 = GetElmIOfNode(ILoadConLocal(3),IElNode3)
                        if (ElementOfNode1==ElementOfNode3) then !found element connected to line
                          ISharedElement = ElementOfNode3 !also = ElementOfNode1 =ElementOfNode2
                          if (IsActiveElement(ISharedElement)) then 
                             do K=1,ELEMENTNODES
                                 INode =ElementConnectivities(K,ISharedElement)
                                 if ((INode/=ILoadConLocal(1)).and.(INode/=ILoadConLocal(2)).and.(INode/=ILoadConLocal(3))) then !this is the node not belonging to the line
                                     Vector(:) = Coord(INode,:) - Coord(ILoadConLocal(1),:)
                                     ScalarProduct = DotProduct(Vector,NormalVector,NVECTOR)
                                     if (ScalarProduct > 0) then ! The normal is pointing inside element, reverse
                                        NormalVector = -NormalVector !the normal is pointing outside active element
                                        EXIT !no need to check the other nodes of the element
                                     elseif (ScalarProduct < 0) then ! The normal is pointing outside element, keep direction
                                        EXIT !no need to check the other nodes of the element   
                                     end if
                                 end if
                             end do
                             EXIT !No need to check the other elements because already found one which is active
                          end if
                            
                        end if
                      end do 
                    end if
                  end do
                end do

            TractionNodes(4:6,1:3) = 0d0 ! remove middle point data
            ProductTractionNodesNormalVector = 0d0 ! initialize

            do K = 1, NNodLOE ! for each corner node
              do J=1,NVECTOR  
                TractionNodes(K,J) = AppliedPressure(K) * NormalVector (J)
              end do
              !do J = 1, NVECTOR ! for each direction
              !   ProductTractionNodesNormalVector(K) = ProductTractionNodesNormalVector(K) + TractionNodes(K,j) * NormalVector (J)
              !end do
              !
              !  do J = 1, NVECTOR
              !     TractionNodes(K, J) = ProductTractionNodesNormalVector(K) * NormalVector(J)
              !  end do
              !  ScalarProduct = DotProduct(TractionNodes(K,:),NormalVector,NVECTOR) !>0 if traction, <0 if compression, normal pointing outward
              !  if (sign(1.d0,ScalarProduct) /= sign(1.d0,LoadValue(K,1))) then !LoadValue(K,1)=LoadValue(K,2)>0 for traction
              !    TractionNodes(K, :) = - TractionNodes(K, :)
              !  end if
            end do
          end if
                    
          ! calculate traction at Gauss Point
          do J = 1, 3
            do K = 1, NNodLOE ! Corner nodes
              NN = ILoadConLocal(K)
              TractionGaussPoint(I,J) = TractionGaussPoint(I,J) + GPShapeFunctionBoundary(I,K) * TractionNodes(K,J)
            end do
          end do
          
          ! calculate traction at Nodes
          do J = 1, NVECTOR
            do K = 1, NNodLOE ! Corner nodes
              Temp = GPWeightBoundary(I) * GPShapeFunctionBoundary(I, K) * VectorLength
              SE(K, J) = SE(K, J) + Temp * TractionGaussPoint(I, J)
            end do 
          end do
        end do

        ! add to global load vector
        do I = 1, NNod
          ND = NDof(ILoadConLocal(I))
          do J = 1, NVECTOR
            RLoad(ND + J) = RLoad(ND + J) + SE(I,J)    
          end do
        end do

      end subroutine Load3D
      
