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


      module ModCounters
      !**********************************************************************
      !
      !    Function:  Stores general mesh information.
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 8842 $
      !     $Date: 2020-07-30 13:58:40 +0200 (do, 30 jul. 2020) $
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
                     NLoadedElementSidesSolidMatPoints = 0, & ! Number of distributed loads on element sides (solid, on material points)
                     NLoadedElementSidesWaterMatPoints = 0, & ! Number of distributed loads on element sides (water, on material points)
                     NLoadedElementSidesGasMatPoints = 0, & ! Number of distributed loads on element sides (gas, on material points)
                     NParticles = 0, & ! number of initialised material points
                     SolidMaterialPoints = 0, & ! number of solid material points
                     LiquidMaterialPoints = 0, & ! number of liquid material points
                     NReactionSurfaceOutput = 0, & !Number of surfaces for output of reaction forces
                     NElemReactions = 0 !number of elements to be considered for the output of reaction surface forces

        end type CounterType
        
        type(CounterType), public, save :: Counters ! Stores general counters relevant to the house-keeping

      contains ! Routines of this module


        subroutine InitCountersFromFile(ApplyContactAlgorithm, NEl, NodTot, MaxLayer, &
                                        NLoadedElementSidesSolidNodes, &
                                        NLoadedElementSidesWaterNodes, &
                                        NLoadedElementSidesGasNodes, &
                                        NLoadedElementSidesSolidMatPoints, &
                                        NLoadedElementSidesWaterMatPoints, &
                                        NLoadedElementSidesGasMatPoints, &
                                        NumberOfMaterials)
        !**********************************************************************
        !
        !    Function:  Assigns the provided values to the data of this module.
        !
        !     ApplyContactAlgorithm : 1 for usage of contact
        !     NEl    : "
        !     NodTot : "
        !     MaxLayer : "
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
                                 NLoadedElementSidesWaterNodes, &
                                 NLoadedElementSidesGasNodes, &
                                 NLoadedElementSidesSolidMatPoints, &
                                 NLoadedElementSidesWaterMatPoints, &
                                 NLoadedElementSidesGasMatPoints, &
                                 NumberOfMaterials
          logical, intent(in) :: ApplyContactAlgorithm
     
          Counters%NEl = NEl
          Counters%NodTot = NodTot
          Counters%NLayers = MaxLayer
          Counters%NLoadedElementSidesSolidNodes = NLoadedElementSidesSolidNodes
          Counters%NLoadedElementSidesWaterNodes = NLoadedElementSidesWaterNodes
          Counters%NLoadedElementSidesGasNodes = NLoadedElementSidesGasNodes
          Counters%NLoadedElementSidesSolidMatPoints = NLoadedElementSidesSolidMatPoints
          Counters%NLoadedElementSidesWaterMatPoints = NLoadedElementSidesWaterMatPoints
          Counters%NLoadedElementSidesGasMatPoints = NLoadedElementSidesGasMatPoints
          Counters%NLoadedElementSidesSolid = NLoadedElementSidesSolidNodes + NLoadedElementSidesSolidMatPoints
          Counters%NLoadedElementSidesWater = NLoadedElementSidesWaterNodes + NLoadedElementSidesWaterMatPoints
          Counters%NLoadedElementSidesGas = NLoadedElementSidesGasNodes + NLoadedElementSidesGasMatPoints
          
          if (ApplyContactAlgorithm) then
              Counters%NEntity = 2
          else
              Counters%NEntity = NumberOfMaterials
          end if
        
        end subroutine InitCountersFromFile
      
      end module ModCounters