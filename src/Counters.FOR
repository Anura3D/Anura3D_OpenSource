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


      module ModCounters
      !**********************************************************************
      !
      !    Function:  Stores general mesh information.
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      
      use ModGlobalConstants
      
      implicit none

        type CounterType
          integer(INTEGER_TYPE) :: &
                     NEl = 0, & ! Number of elements
                     NAEl = 0, & ! Number of active elements
                     NodTot = 0, & ! Number of nodes
                     N = 0, & ! Number of degrees of freedom
                     NEntity = 0, & ! Dynamic contact - set to 1 if no contact model is used, when contact is used = number of entities  CC
                     NLayers = 0, & ! Number of soil layers/materials
                     NLoadedElementSidesSolid = 0, & ! Number of distributed loads on element sides (solid)
                     NLoadedElementSidesWater = 0, & ! Number of distributed loads on element sides (water)
                     NLoadedElementSidesGas = 0, & ! Number of distributed loads on element sides (gas)
                     NLoadedElementSidesSolidNodes = 0, & ! Number of distributed loads on element sides (solid, on nodes)
                     NLoadedElementSidesWaterNodes = 0, & ! Number of distributed loads on element sides (water, on nodes)
                     NLoadedElementSidesGasNodes = 0, & ! Number of distributed loads on element sides (gas, on nodes)
                     NLoadedElementSidesSolidNodesB = 0, & ! Number of distributed loads on element sides (solid, on nodes, LOAD SYSTEM B)
                     NLoadedElementSidesWaterNodesB = 0, & ! Number of distributed loads on element sides (water, on nodes, LOAD SYSTEM B)
                     NLoadedElementSidesGasNodesB = 0, & ! Number of distributed loads on element sides (gas, on nodes, LOAD SYSTEM B)
                     NLoadedElementSidesSolidMatPoints = 0, & ! Number of distributed loads on element sides (solid, on material points)
                     NLoadedElementSidesWaterMatPoints = 0, & ! Number of distributed loads on element sides (water, on material points)
                     NLoadedElementSidesGasMatPoints = 0, & ! Number of distributed loads on element sides (gas, on material points)
                     NLoadedElementSidesSolidMatPointsB = 0, & ! Number of distributed loads on element sides (solid, on material points)
                     NLoadedElementSidesWaterMatPointsB = 0, & ! Number of distributed loads on element sides (water, on material points)
                     NLoadedElementSidesGasMatPointsB = 0, & ! Number of distributed loads on element sides (gas, on material points)
                     SoilSurfaceNumberofSides = 0, & ! Number of element side on the soil surface ( on nodes)
                     PhreaticSurfaceNumberofSides = 0, & ! Number of element side on the phreatic surface ( on nodes)
                     NParticles = 0, & ! number of initialised material points
                     SolidMaterialPoints = 0, & ! number of solid material points
                     LiquidMaterialPoints = 0, & ! number of liquid material points
                     NReactionSurfaceOutput = 0, & !Number of surfaces for output of reaction forces
                     NElemReactions = 0, & !number of elements to be considered for the output of reaction surface forces
                     HydraulicHeadSides = 0,&
                     NSoilLoadSystems = 1,& !number of load systems (solid)
                     NWaterLoadSystems = 1,& !number of load systems (water)
                     NGasLoadSystems = 1 !number of load systems (gas)
                     
        end type CounterType
        
        type(CounterType), public, save :: Counters ! Stores general counters relevant to the house-keeping

      contains ! Routines of this module


        subroutine InitCountersFromFile(ApplyContactAlgorithm, NEl, NodTot, MaxLayer, &
                                        NLoadedElementSidesSolidNodes, &
                                        NLoadedElementSidesSolidNodesB, &
                                        NLoadedElementSidesWaterNodes, &
                                        NLoadedElementSidesWaterNodesB, &
                                        NLoadedElementSidesGasNodes, &
                                        NLoadedElementSidesGasNodesB, &
                                        NLoadedElementSidesSolidMatPoints, &
                                        NLoadedElementSidesSolidMatPointsB, &
                                        NLoadedElementSidesWaterMatPoints, &
                                        NLoadedElementSidesWaterMatPointsB, &
                                        NLoadedElementSidesGasMatPoints, &
                                        NLoadedElementSidesGasMatPointsB, &
                                        SoilSurfaceNumberofSides, &
                                        PhreaticSurfaceNumberofSides, &
                                        NumberOfMaterials)                                        
        !**********************************************************************
        !
        !    Function:  Assigns the provided values to the data of this module.
        !
        !     ApplyContactAlgorithm : 1 for usage of contact
        !     NEl    : Number of element
        !     NodTot : Number of total nodes
        !     MaxLayer : Number of layers
        !     NLoadedElementSidesSolidNodes : Number of distributed loads on element sides (solid, on nodes)
        !     NLoadedElementSidesWaterNodes : Number of distributed loads on element sides (water, on nodes)
        !     NLoadedElementSidesGasNodes : Number of distributed loads on element sides (gas, on nodes)
        !     NLoadedElementSidesSolidMatPoints : Number of distributed loads on element sides (solid, on material points)
        !     NLoadedElementSidesWaterMatPoints : Number of distributed loads on element sides (water, on material points)
        !     NLoadedElementSidesGasMatPoints : Number of distributed loads on element sides (gas, on material points)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
     
        implicit none
        
          ! Routine parameters
          integer(INTEGER_TYPE), intent(in) :: &
                                 NEl, NodTot, MaxLayer, &
                                 NLoadedElementSidesSolidNodes, &
                                 NLoadedElementSidesSolidNodesB, &
                                 NLoadedElementSidesWaterNodes, &
                                 NLoadedElementSidesWaterNodesB, &
                                 NLoadedElementSidesGasNodes, &
                                 NLoadedElementSidesGasNodesB, &
                                 NLoadedElementSidesSolidMatPoints, &
                                 NLoadedElementSidesSolidMatPointsB, &
                                 NLoadedElementSidesWaterMatPoints, &
                                 NLoadedElementSidesWaterMatPointsB, &
                                 NLoadedElementSidesGasMatPoints, &
                                 NLoadedElementSidesGasMatPointsB, &
                                 NumberOfMaterials, &
                                 SoilSurfaceNumberofSides, &
                                 PhreaticSurfaceNumberofSides
          logical, intent(in) :: ApplyContactAlgorithm
     
          Counters%NEl = NEl
          Counters%NodTot = NodTot
          Counters%NLayers = MaxLayer
          Counters%NLoadedElementSidesSolidNodes = NLoadedElementSidesSolidNodes
          Counters%NLoadedElementSidesWaterNodes = NLoadedElementSidesWaterNodes
          Counters%NLoadedElementSidesGasNodes = NLoadedElementSidesGasNodes
          
          Counters%NLoadedElementSidesSolidNodesB = NLoadedElementSidesSolidNodesB
          Counters%NLoadedElementSidesWaterNodesB = NLoadedElementSidesWaterNodesB
          Counters%NLoadedElementSidesGasNodesB = NLoadedElementSidesGasNodesB
          
          Counters%NLoadedElementSidesSolidMatPoints = NLoadedElementSidesSolidMatPoints
          Counters%NLoadedElementSidesWaterMatPoints = NLoadedElementSidesWaterMatPoints
          Counters%NLoadedElementSidesGasMatPoints = NLoadedElementSidesGasMatPoints
          
          Counters%NLoadedElementSidesSolidMatPointsB = NLoadedElementSidesSolidMatPointsB
          Counters%NLoadedElementSidesWaterMatPointsB = NLoadedElementSidesWaterMatPointsB
          Counters%NLoadedElementSidesGasMatPointsB = NLoadedElementSidesGasMatPointsB
          
          Counters%NLoadedElementSidesSolid = NLoadedElementSidesSolidNodes + NLoadedElementSidesSolidMatPoints &
              + NLoadedElementSidesSolidNodesB + NLoadedElementSidesSolidMatPointsB
          Counters%NLoadedElementSidesWater = NLoadedElementSidesWaterNodes + NLoadedElementSidesWaterMatPoints &
              + NLoadedElementSidesWaterNodesB + NLoadedElementSidesWaterMatPointsB
          Counters%NLoadedElementSidesGas = NLoadedElementSidesGasNodes + NLoadedElementSidesGasMatPoints &
              + NLoadedElementSidesGasNodesB + NLoadedElementSidesGasMatPointsB
          Counters%SoilSurfaceNumberofSides = SoilSurfaceNumberofSides
          Counters%PhreaticSurfaceNumberofSides = PhreaticSurfaceNumberofSides
          if (ApplyContactAlgorithm) then
              Counters%NEntity = 2
          else
              Counters%NEntity = 1
          end if
        
        end subroutine InitCountersFromFile
      
      end module ModCounters