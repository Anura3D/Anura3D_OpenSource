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
      !use ModMPMData

      
      implicit none

      
      integer(INTEGER_TYPE) :: IPatch, ILoadSystem
      
      
      do IPatch = 1, Counters%NPatches ! loop over patches       
          do ILoadSystem = 1, NURBS%NLoadSystems_Traction(IPatch) ! loop over load systems
      
        if ( Counters%NLoadedElementSidesSolid(ILoadSystem,IPatch)  > 0) then ! solid !Counters%NLoadedElementSidesSolid    !Counters%NLoadedElementSidesSolidNodes
          call InitialiseTractionLoadVector(LOADTYPE_SOLID, ILoadSystem, IPatch)
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVector, ExtLoadVector) ! 3D function, only if IS3DCYLINDRIC=.true.
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorB, ExtLoadVectorB) ! 3D function, only if IS3DCYLINDRIC=.true.
        end if
        
          end do ! loop over load systems
      end do  ! loop over patches
      
      
      do IPatch = 1, Counters%NPatches ! loop over patches       
          do ILoadSystem = 1, NURBS%NLoadSystems_Traction(IPatch) ! loop over load systems
        
        if (((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) .or. (Counters%NLoadedElementSidesWater(ILoadSystem,IPatch) > 0)) then ! liquid
        !just temporary: need to find a different condition
          call InitialiseTractionLoadVector(LOADTYPE_LIQUID, ILoadSystem, IPatch)  
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorwater, ExtLoadVectorWater) ! 3D function, only if IS3DCYLINDRIC=.true.
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorwaterB, ExtLoadVectorWaterB) ! 3D function, only if IS3DCYLINDRIC=.true.
          ! note: the rotation of the BC for the liquid has to be fixed
        end if
      
          end do ! loop over load systems
      end do ! loop over patches
      
        if (Counters%NLoadedElementSidesGas > 0) then ! gas
          call InitialiseTractionLoadVector(LOADTYPE_GAS, ILoadSystem, IPatch)  
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorGas, ExtLoadVectorGas) ! 3D function, only if IS3DCYLINDRIC=.true.
          if (NDIM == 3) call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadVectorGasB, ExtLoadVectorGasB) ! 3D function, only if IS3DCYLINDRIC=.true.
          ! note: the rotation of the BC for the gas has to be fixed
        end if

      end subroutine InitialiseTractionLoad
      
      
      subroutine InitialiseTractionLoadVector(LoadType, ILoadSystem, IPatch)
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
      !use ModMPMData
      
      implicit none
      
        ! arguments
      integer(INTEGER_TYPE), intent(in) :: LoadType

      ! Local variables
      integer(INTEGER_TYPE) :: I, LoadedSides, NNodes, ee
      
      ! Multipatch variables 
      integer(INTEGER_TYPE) :: IPatch_Temporary
      integer(INTEGER_TYPE), intent(in) :: IPatch
      
      
      integer(INTEGER_TYPE) :: ILoadSystem_Temporary
      integer(INTEGER_TYPE), intent(in) :: ILoadSystem
      
      !real(REAL_TYPE), dimension(Counters%N) :: ExtLoadVector_MapToPoints
      
      IPatch_Temporary = 1
      
      !ExtLoadVector_MapToPoints = 0.0
      
      NNodes = NURBS%NXiKnotOrder(IPatch)+1

      if ( LoadType == LOADTYPE_SOLID ) then ! solid

          ! initialise ExtLoadVector with distributed load on nodes
          if (NDIM == 3) then ! for 3D
              
            !do IPatch = 1, Counters%NPatches
            !  
            !    
            !    do ILoadSystem = 1, NLoadSystems_Traction(IPatch)
              
                    
                    
                    
                    
                    
            do I = 1, Counters%NLoadedElementSidesSolidNodes(ILoadSystem, IPatch)
                ! IEN_Traction --> subscripts are hardcoded --> we need to loop over the list of elements
                ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
              
                
                
                ! we need to check the directions of the traction load system and input the knot vectors 
                ! XI_ETA
                if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
                
                call Load3D(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                    ILoadSystem, IPatch, &
                    NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                    NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch) )
            
                ! XI_ZETA
                elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
                    
                call Load3D(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                    ILoadSystem, IPatch, &
                    NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                    NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch) )
                
                
                ! ETA_ZETA
                elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 3) then 
                    
                call Load3D(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                    ILoadSystem, IPatch, &
                    NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch), &
                    NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch) )
                
                    
                end if 
                
            
            end do
            
            
            
            !    end do 
            !    
            !    
            !    
            !    
            !end do 
            
            ! I COMMENTED THE LOAD SYSTEM B FOR NOW
            !do I = 1, Counters%NLoadedElementSidesSolidNodesB
            !    ee = ElementUponWhichLoadIsApplied_NURBS_B(I) 
            !  !call Load3D(ExtLoadVectorB, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesSolidB(1, I), LoadValuesOnNodesSolidB(1, 1, I), LoadType)
            !  !call Load3D(ExtLoadVectorB, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesSolidB(1, I), LoadValuesOnNodesSolidB(1, 1, I), LoadType)
            !  call Load3D(ExtLoadVectorB, ReducedDof, NodalCoordinates_Traction_B(:,:), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, IEN_Traction_B(:,ee), ElementUponWhichLoadIsApplied_NURBS_B(I))  
            !end do
          elseif (NDIM == 2) then ! for 2D
              
              !do IPatch = 1, Counters%NPatches ! loop over patches
              !
              !    do ILoadSystem = 1, NLoadSystems_Traction(IPatch) ! loop over load systems
                      
            do I = 1, Counters%NLoadedElementSidesSolidNodes(ILoadSystem, IPatch) ! loop over elements
              !call Load2D(ExtLoadVector, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code   
              
                ! IEN_Traction --> subscripts are hardcoded --> we need to loop over the list of elements
                ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
                
                 !we need to check the directions of the traction load system and input the knot vectors 
                ! XI
                if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
                
                call Load2D_NURBS(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code
              !call Load2D(ExtLoadVector, ReducedDof, NodalCoordinates, 2, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! let's try two here 

                
                ! ETA
                else if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
                    
                    
                    call Load2D_NURBS(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch) )
                    
                    
                end if 
                
            end do ! loop over elements
            !do I = 1, Counters%NLoadedElementSidesSolidNodesB
            !  call Load2D(ExtLoadVectorB, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, LoadOnNodesConnectivitiesSolidB(1, I), LoadValuesOnNodesSolidB(1, 1, I), LoadType)              
            !end do
            
            !      end do ! loop over load systems
            !
            !end do ! loop over patches
          
          end if  
          
      !end if
      
        
          ! initialise Particle%ExtLoadVector with distributed load on material points (only in first load step)
          if (.not.IsFollowUpPhase()) then
            if (NDIM == 3) then ! for 3D only, sprint #2: add 2D functionality
              
                
                !do IPatch = 1, Counters%NPatches ! loop over patches
                !
                !  do ILoadSystem = 1, NLoadSystems_Traction(IPatch) ! loop over load systems
                      
                
                do I = 1, Counters%NLoadedElementSidesSolidMatPoints(ILoadSystem, IPatch)
                  ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch)                   
                  !call Load3D_NURBS_OnMPs(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES)!, IEN_Traction)!, ElementUponWhichLoadIsApplied_NURBS(I))
                  !call Load3D(ExtLoadVector_MapToPoints, ReducedDof, NodalCoordinates_Traction(:,:), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, IEN_Traction(:,ee), ElementUponWhichLoadIsApplied_NURBS(I))          
                 
                  ! XI_ETA
                  if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
                  
                  call Load3D_NURBS_OnMPs(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch) )
                  
                  
                  ! XI_ZETA
                  elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
                      
                      
                  call Load3D_NURBS_OnMPs(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch) )
                      
                      
                  ! ETA_ZETA
                  elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 3) then 
                      
                      
                      
                  call Load3D_NURBS_OnMPs(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch) )
                      
                      
                      
                      
                  end if 
                  
                  
                  
                end do
                
                
                
                
                !  end do ! loop over load systems
                !  
                !end do ! loop over patches
                
              
              
              !do I = 1, Counters%NLoadedElementSidesSolidMatPoints
              !    ee = ElementUponWhichLoadIsApplied_NURBS(I)               
              !    call TransferLoadFromNodesToMP(ExtLoadVector_MapToPoints, ReducedDof, ee)
              !end do 
              
                  !I need to update the inputs to this subroutine
                  !call TransferExternalLoadsToParticles_NURBS(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, 1, 1)
                
                !call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesSolid(1, I), LoadValuesOnMatPointsSolid(1, 1, I), 1, NodalCoordinates(:,:,IPatch_Temporary),1)

              !end do
              do I = 1, Counters%NLoadedElementSidesSolidMatPointsB
                call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesSolidB(1, I), LoadValuesOnMatPointsSolidB(1, 1, I), 1, NodalCoordinates(:,:,IPatch_Temporary),2)
			  end do
            elseif (NDIM == 2) then ! for 2D
                
                
                
                !do IPatch = 1, Counters%NPatches
                !
                !do ILoadSystem = 1, NLoadSystems_Traction(IPatch)
                
              do I = 1, Counters%NLoadedElementSidesSolidMatPoints(ILoadSystem, IPatch)
                  
                  ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
                  
                  !I need to update the inputs to this subroutine
                  !call TransferExternalLoadsToParticles_NURBS(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, 1, 1)
               
                  !do I = 1, Counters%NLoadedElementSidesSolidNodes
                  !call Load2D(ExtLoadVector, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code   
                  ! IEN_Traction --> subscripts are hardcoded --> we need to loop over the list of elements
                  !Traction_Connecting_Local_To_Global_Elements(I)!

                  ! XI
                  if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then
                  
                  call Load2D_NURBS_onMPs(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code   
                  
                  ! ETA
                  else if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then
                  
                  call Load2D_NURBS_onMPs(ExtLoadVector, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code   
                  
                  
                  end if 
                  
                  
                  

                  !call Load2D(ExtLoadVector, ReducedDof, NodalCoordinates, 2, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! let's try two here 
                  !end do

                  !call TransferExternalLoadsToParticles(LoadType, ELEMENTBOUNDARYNODES, ELEMENTNODES, LoadOnMatPointsConnectivitiesSolid(1, I), LoadValuesOnMatPointsSolid(1, 1, I), 2, NodalCoordinates(:,:,IPatch_Temporary),1)  !hardcoded 2 gauss points

                  ! I need to map ExtLoadVector to the material points 
              
              end do
              
                !end do 
                !
                !end do 
                
              do I = 1, Counters%NLoadedElementSidesSolidMatPointsB
                call TransferExternalLoadsToParticles(LoadType, ELEMENTBOUNDARYNODES, ELEMENTNODES, LoadOnMatPointsConnectivitiesSolidB(1, I), LoadValuesOnMatPointsSolidB(1, 1, I), 1, NodalCoordinates(:,:,IPatch_Temporary),2)
              end do
            end if
          end if
              
      else if ( LoadType == LOADTYPE_LIQUID ) then ! liquid
          ! initialise ExtLoadVectorWater with distributed load on nodes
       ! Hydraulic head on NODES
          
          
        !  if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
        !      
        !      call GetAndUpdateHydraulicHeadLoad(HydraulicHeadLoadedElemID, LoadedSides, HydraulicHeadNodesConnectivities, HydraulicHeadLoad)
        !      Counters%HydraulicHeadSides = LoadedSides
        !      
        !      if (NDIM == 3) then
        !       do I = 1, Counters%HydraulicHeadSides
        !        call Load3D(HydraulicHeadVector, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, HydraulicHeadNodesConnectivities(1, I), HydraulicHeadLoad(1, 1, I), LoadType)
        !    end do
        !elseif (NDIM == 2) then
        !    do I = 1, Counters%HydraulicHeadSides
        !        call Load2D(HydraulicHeadVector, ReducedDof, NodalCoordinates, 1, HydraulicHeadNodesConnectivities(1, I), HydraulicHeadLoad(1, 1, I), LoadType)
        !    end do
        !end if
        !  end if
          ! general pressure load on NODES
          
          ! initialise ExtLoadVectorWater with distributed load on nodes
          if (NDIM == 3) then ! for 3D
              
              !do IPatch = 1, Counters%NPatches
              !    
              !    do ILoadSystem = 1, NloadSystems_Traction(IPatch)
              
              do I = 1, Counters%NLoadedElementSidesWaterNodes(ILoadSystem, IPatch)
                  ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
                  !call Load3D(ExtLoadVectorWater, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesWater(1, I), LoadValuesOnNodesWater(1, 1, I), LoadType)
                  
                  
                  ! we need to check the directions of the traction load system and input the knot vectors 
                  ! XI_ETA
                  if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
                  
                  call Load3D(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch))
                  
                  ! XI_ZETA
                  elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
                  
                  call Load3D(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch))
                  
                  ! ETA_ZETA
                  elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 3) then 
                  
                  call Load3D(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch), &
                      NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch))
                  
                  
                  end if 
                  
                  
                  
                  
                  
              end do
              
              
              !    end do 
              !    
              !end do
              
              !do I = 1, Counters%NLoadedElementSidesWaterNodesB
              !    call Load3D(ExtLoadVectorWaterB, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesWaterB(1, I), LoadValuesOnNodesWaterB(1, 1, I), LoadType)
              !end do
          else if (NDIM == 2) then ! for2D
              
              !do IPatch = 1, Counters%NPatches
              !    
              !    do ILoadSystem = 1, NLoadSystems_Traction(IPatch)
              
              
              do I = 1, Counters%NLoadedElementSidesWaterNodes(ILoadSystem, IPatch)                  
                  ! IEN_Traction --> subscripts are hardcoded --> we need to loop over the list of elements
                ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
                
                
                
                
                ! XI
                if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then
                
                call Load2D_NURBS(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code
                
                ! ETA
                else if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then
                
                call Load2D_NURBS(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code
                
                end if 
                
                
                
                
                
                
                
                
                
                  !call Load2D(ExtLoadVectorWater, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesWater(1, I), LoadValuesOnNodesWater(1, 1, I), LoadType)
              end do 
              !do I = 1, Counters%NLoadedElementSidesWaterNodesB
              !    call Load2D(ExtLoadVectorWaterB, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesWaterB(1, I), LoadValuesOnNodesWaterB(1, 1, I), LoadType)
              !end do
              
              !    end do 
              !
              !end do 
              
          end if 
        
          ! initialise Particle%ExtLoadVectorWater with distributed load on material points (only in first load step)
          if (.not.IsFollowUpPhase()) then
            if (NDIM == 3) then ! for 3D
                
                !do IPatch = 1, Counters%NPatches
                !    
                !    
                !    do ILoadSystem = 1, NLoadSystems_Traction(IPatch)
                
                
              do I = 1, Counters%NLoadedElementSidesWaterMatPoints(ILoadSystem, IPatch) 
                !call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesWater(1, I), LoadValuesOnMatPointsWater(1, 1, I), 1, NodalCoordinates,1)
                                  
                 
                  ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
                  
                  
                  
                  
                  
                  ! we need to check the directions of the traction load system and input the knot vectors 
    ! XI_ETA
    if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
    
        call Load3D_NURBS_OnMPs(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch))
    

    ! XI_ZETA
    elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
        
        call Load3D_NURBS_OnMPs(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch))
    
    
    ! ETA_ZETA
    elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 3) then 
        
        call Load3D_NURBS_OnMPs(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch), &
                      NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch))
    
        
    end if 
                  
                  
                  
                  
                  
                  
                  
                  
                  !call Load3D_NURBS_OnMPs(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, IEN_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                  !    ILoadSystem, IPatch, &
                  !    NXiKnotEntries(IPatch), XiKnotEntries(1:NXiKnotEntries(IPatch), IPatch), &
                  !    NEtaKnotEntries(IPatch), EtaKnotEntries(1:NEtaKnotEntries(IPatch), IPatch))

                !call Load2D_NURBS_onMPs(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, IEN_Traction(:,ee), ElementUponWhichLoadIsApplied_NURBS(I))! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code   

                
              end do
              
                !    end do 
                !    
                !end do 
                
              
              
              do I = 1, Counters%NLoadedElementSidesWaterMatPointsB
                call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesWaterB(1, I), LoadValuesOnMatPointsWaterB(1, 1, I), 1, NodalCoordinates,2)
              end do
            elseif (NDIM == 2) then ! for 2D
              
                !do IPatch = 1, Counters%NPatches ! loop over patches
                !
                !do ILoadSystem = 1, NLoadSystems_Traction(IPatch) ! loop over load systems
                
                
                do I = 1, Counters%NLoadedElementSidesWaterMatPoints(ILoadSystem, IPatch) ! loop over 
                
                  !call TransferExternalLoadsToParticles(LoadType, ELEMENTBOUNDARYNODES, ELEMENTNODES, LoadOnMatPointsConnectivitiesWater(1, I), LoadValuesOnMatPointsWater(1, 1, I), 1, NodalCoordinates,1)  
                    
                  ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
                  
                  ! XI  
                  if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
                  
                call Load2D_NURBS_onMPs(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code   

                ! ETA  
                  else if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
                      
                call Load2D_NURBS_onMPs(ExtLoadVectorWater, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code   
      
                  end if 
                  
                
                !, &
                !    NEtaKnotEntries(IPatch), EtaKnotEntries(1:NEtaKnotEntries(IPatch), IPatch)
                  
                end do
                
                
                !end do 
                !
                !
                !end do
                
              !do I = 1, Counters%NLoadedElementSidesWaterMatPointsB
              !  call TransferExternalLoadsToParticles(LoadType, ELEMENTBOUNDARYNODES, ELEMENTNODES, LoadOnMatPointsConnectivitiesWaterB(1, I), LoadValuesOnMatPointsWaterB(1, 1, I), 1, NodalCoordinates,2)  
              !end do
            end if
          end if  
      
        else if ( LoadType == LOADTYPE_GAS ) then ! gas
      
          ! initialise ExtLoadVectorGas with distributed load on nodes
          if (NDIM == 3) then ! for 3D only, sprint #2: add 2D functionality
            do I = 1, Counters%NLoadedElementSidesGasNodes
              
                !call Load3D(ExtLoadVectorGas, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesGas(1, I), LoadValuesOnNodesGas(1, 1, I), LoadType)
                  ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  ! we need to check the directions of the traction load system and input the knot vectors 
    ! XI_ETA
    if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
    
        call Load3D(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                    NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch) )
    

    ! XI_ZETA
    elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
        
        call Load3D(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                    NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch) )
    
    
    ! ETA_ZETA
    elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 3) then 
        
        call Load3D(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch), &
                    NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch) )
    
    end if 
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                !call Load3D(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, IEN_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                !    NXiKnotEntries(IPatch), XiKnotEntries(1:NXiKnotEntries(IPatch), IPatch), &
                !    NEtaKnotEntries(IPatch), EtaKnotEntries(1:NEtaKnotEntries(IPatch), IPatch) )

            end do
            
            !do I = 1, Counters%NLoadedElementSidesGasNodesB
            !  call Load3D(ExtLoadVectorGasB, ReducedDof, NodalCoordinates, 1, N_BOUNDARY_NODES_HOE, LoadOnNodesConnectivitiesGasB(1, I), LoadValuesOnNodesGasB(1, 1, I), LoadType)
            !end do
           
          else if (NDIM == 2) then ! for2D
              
              do I = 1, Counters%NLoadedElementSidesGasNodes
                  !call Load2D(ExtLoadVectorGas, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesGas(1, I), LoadValuesOnNodesGas(1, 1, I), LoadType)
              
                  ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch) 
                  
                  ! XI
                  if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then
                  
                call Load2D_NURBS(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code
                  
                ! ETA  
                  else if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
                
                call Load2D_NURBS(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS,  ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), ILoadSystem, IPatch, &
                    NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch) )! LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid(1, 1, I), LoadType) ! hardcoded as one in the code
                
                end if 
                  
                  
                  
              end do 
              
              !do I = 1, Counters%NLoadedElementSidesWaterNodesB
              !    call Load2D(ExtLoadVectorGasB, ReducedDof, NodalCoordinates, 1, LoadOnNodesConnectivitiesGasB(1, I), LoadValuesOnNodesGasB(1, 1, I), LoadType)
              !end do
          end if
        
          ! initialise Particle%ExtLoadVectorGas with distributed load on material points (only in first load step)
          if (.not.IsFollowUpPhase()) then
            if (NDIM == 3) then ! for 3D only, sprint #2: add 2D functionality
              do I = 1, Counters%NLoadedElementSidesGasMatPoints
                
                  ee = ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem_Temporary, IPatch_Temporary) 
                  
                  
                  
                  ! we need to check the directions of the traction load system and input the knot vectors 
    ! XI_ETA
    if (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 1) then 
    
        call Load3D_NURBS_OnMPs(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch))
    

    ! XI_ZETA
    elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 2) then 
        

        call Load3D_NURBS_OnMPs(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NXiKnotEntries(IPatch), NURBS%XiKnotEntries(1:NURBS%NXiKnotEntries(IPatch), IPatch), &
                      NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch))

    ! ETA_ZETA
    elseif (NURBS%NURBS_Traction_Dimension(ILoadSystem, IPatch) == 3) then 
        
    
        call Load3D_NURBS_OnMPs(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NURBS%ElementConnectivities_Traction(:,ee, ILoadSystem, IPatch), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem, IPatch), &
                      ILoadSystem, IPatch, &
                      NURBS%NEtaKnotEntries(IPatch), NURBS%EtaKnotEntries(1:NURBS%NEtaKnotEntries(IPatch), IPatch), &
                      NURBS%NZetaKnotEntries(IPatch), NURBS%ZetaKnotEntries(1:NURBS%NZetaKnotEntries(IPatch), IPatch))

        
    end if 
                  
                  
                  
                  
                  
                  !call Load3D_NURBS_OnMPs(ExtLoadVectorGas, ReducedDof, NodalCoordinates_Traction(:,:, ILoadSystem_Temporary, IPatch_Temporary), ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, IEN_Traction(:,ee, ILoadSystem_Temporary, IPatch_Temporary), ElementUponWhichLoadIsApplied_NURBS(I, ILoadSystem_Temporary, IPatch_Temporary), &
                  !    ILoadSystem, IPatch, &
                  !    NXiKnotEntries(IPatch), XiKnotEntries(1:NXiKnotEntries(IPatch), IPatch), &
                  !    NEtaKnotEntries(IPatch), EtaKnotEntries(1:NEtaKnotEntries(IPatch), IPatch))

                  !call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesGas(1, I), LoadValuesOnMatPointsGas(1, 1, I), 1, NodalCoordinates,1)
              end do
              !do I = 1, Counters%NLoadedElementSidesGasMatPointsB
              !  call TransferExternalLoadsToParticles(LoadType, N_BOUNDARY_NODES_HOE, N_NODES_HOE, LoadOnMatPointsConnectivitiesGasB(1, I), LoadValuesOnMatPointsGasB(1, 1, I), 1, NodalCoordinates,2)
              !end do
            end if
          end if  
       
       end if ! LoadType
       
    end subroutine InitialiseTractionLoadVector
    
    
    
    
    
    !subroutine TransferLoadFromNodesToMP(RLoad, NDOF, ElementUponWhichLoadIsApplied)
    !
    !use ModReadCalculationData
    !  use ModMPMInit
    !  use ModMeshInfo
    !  use ModGlobalConstants
    !  use ModDYNConvectivePhase  
    !  
    !implicit none
    !! - loop over elements 
    !! - pick out the dofs associated with each element based on the dofs
    !real(REAL_TYPE), intent(in), dimension(Counters%N) :: RLoad
    !integer(INTEGER_TYPE), intent(in) :: ElementUponWhichLoadIsApplied
    !
    !integer(INTEGER_TYPE), dimension(Counters%NodTot), intent(in) :: NDof
    !
    !real(REAL_TYPE), dimension(ELEMENTBOUNDARYNODES, NVECTOR) :: SE
    !integer(INTEGER_TYPE) :: I, J, ND
    !
    !
    !  !
    !  !use ModCounters
    !  !use ModElementEvaluation
    !  !use ModGlobalConstants
    !  !use ModMeshAdjacencies
    !  !  use ModMeshInfo
    !  !
    !  !implicit none
    !  
    !    ! arguments
    !    !integer(INTEGER_TYPE), intent(in) :: NInt, NNod!, LoadType 
    !    integer(INTEGER_TYPE), dimension(NNod) :: ILoadCon
    !    !integer(INTEGER_TYPE), dimension(Counters%NodTot) :: NDof
    !    !real(REAL_TYPE), dimension(NNod, NVECTOR), intent(in) ::  LoadValue
    !    !real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Coord       ! Coordinates vector
    !    !real(REAL_TYPE), dimension(Counters%NodTot*NDOFL), intent(inout) :: RLoad      ! Resulting load vector
    !    !integer(INTEGER_TYPE), intent(in) :: ElementUponWhichLoadIsApplied
    !  
    !    ! local variable
    !    integer(INTEGER_TYPE) :: K, NN, NNodLOE
    !    integer(INTEGER_TYPE), dimension(NNod) :: ILoadConLocal
    !    real(REAL_TYPE) :: VectorLength, Temp
    !    real(REAL_TYPE), dimension(NVECTOR) :: NormalVector, ProductTractionNodesNormalVector, Vector
    !    real(REAL_TYPE), dimension(NInt, NVECTOR) :: TractionGaussPoint
    !    real(REAL_TYPE), dimension(NNod, NVECTOR) :: TractionNodes
    !    real(REAL_TYPE), dimension(NNod) :: AppliedPressure
    !    
    !    real(REAL_TYPE) :: RLoad_DEBUG
    !    real(REAL_TYPE) :: SE_DEBUG
    !    
    !    integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, NElemNode3, IElNode1, IElNode2, IElNode3, ElementOfNode1, ElementOfNode2, ElementOfNode3, ISharedElement, INode
    !    real(REAL_TYPE) :: ScalarProduct
    !     
    !    
    !    ! arguments
    !    !!integer(INTEGER_TYPE), intent(in) :: NInt, NNod!, LoadType 
    !    !integer(INTEGER_TYPE), dimension(NNod) :: ILoadCon
    !    !integer(INTEGER_TYPE), dimension(Counters%NodTot) :: NDof
    !    !!real(REAL_TYPE), dimension(NNod, NVECTOR), intent(in) ::  LoadValue
    !    !real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Coord       ! Coordinates vector
    !    !real(REAL_TYPE), dimension(Counters%NodTot*NDOFL), intent(inout) :: RLoad      ! Resulting load vector
    !    !integer(INTEGER_TYPE), intent(in) :: ElementUponWhichLoadIsApplied
    !  
    !    
    !    ! local variable
    !    !integer(INTEGER_TYPE) :: I, J, K, NN, ND, NNodLOE
    !    !integer(INTEGER_TYPE), dimension(NNod) :: ILoadConLocal
    !    !real(REAL_TYPE) :: VectorLength, Temp
    !    !real(REAL_TYPE), dimension(NVECTOR) :: NormalVector, ProductTractionNodesNormalVector, Vector
    !    !real(REAL_TYPE), dimension(NInt, NVECTOR) :: TractionGaussPoint
    !    !real(REAL_TYPE), dimension(NNod, NVECTOR) :: TractionNodes, SE
    !    !real(REAL_TYPE), dimension(NNod) :: AppliedPressure
    !    
    !    real(REAL_TYPE) :: RLoad_DEBUG
    !    
    !     integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, NElemNode3, IElNode1, IElNode2, IElNode3, ElementOfNode1, ElementOfNode2, ElementOfNode3, ISharedElement, INode
    !     real(REAL_TYPE) :: ScalarProduct
    !     
    !     !integer(INTEGER_TYPE) :: ELEMENTUPONWHICHLOADISAPPLIED
    !  
    !     ! Multipatch variables 
    !     integer(INTEGER_TYPE) :: IPatch_Temporary = 1
    !     
    !     
    !     
    !     
    !     integer(INTEGER_TYPE) :: LoadedParticleCounter
    !     integer(INTEGER_TYPE) :: IDim
    !     real(REAL_TYPE) :: Factor
    !     integer(INTEGER_TYPE) :: NParticles, ParticleID, GlobalElementForTraction, IError
    !     integer(INTEGER_TYPE) :: NLoadedParticles, IParticle, ParticleIndex
    !              
    !     logical, allocatable, dimension(:) :: ParticleStatus
    !
    !     
    !     real(REAL_TYPE), dimension(NVECTOR) :: LoadPerParticle
    !     
    !     integer(INTEGER_TYPE), dimension(NNod) :: LocalNodeConnectionArray
    !
    !SE = 0.0
    !
    !
    !! add to global load vector
    !do I = 1, ELEMENTBOUNDARYNODES!NNod  ! --> ELEMENTBOUNDARYNODES is the number of nodes in the element (eg. 2*2 if traction on a linear square element)
    !    ND = NDOF(Traction_Connecting_Local_To_Global_ControlPoints(IEN_TRACTION(I,ElementUponWhichLoadIsApplied)))!NDof(ILoadConLocal(I))
    !    do J = 1, NVECTOR
    !        SE(I,J) = RLoad(ND + J)
    !        
    !        !RLoad(ND + J) = RLoad(ND + J) + SE(I,J)    
    !        !RLoad_DEBUG = RLoad_DEBUG + RLoad(ND + J)
    !    end do    
    !end do
    !
    !
    !    !---------------------------------------------------------------------------------------------------------------------
    !    ! TRACTION ON THE MPs
    !    
    !    ! This approach obviously wouldn't work because there is a problem with how this is being set up. 
    !    ! You are mapping these forces from the nodes to  
    !    
    !    
    !    
    !    
    !    !Factor = real(ELEMENTBOUNDARYNODES) / real(ELEMENTBOUNDARYGAUSSPOINTS)
    !    
    !    
    !          NParticles = 0
    !          
    !          
    !          !do IElement_Traction = 1, nel_NURBS_Traction ! traction elements
    !              
    !              ParticleID = 0
    !          LoadedParticleCounter = 0
    !              
    !              GlobalElementForTraction = Traction_Connecting_Local_To_Global_Elements(ElementUponWhichLoadIsApplied)
    !              
    !          NParticles =  NPartEle(GlobalElementForTraction, IPatch_Temporary) ! HARDCODED
    !          
    !          allocate(ParticleStatus(NParticles), stat = IError)
    !          
    !          ParticleStatus = .false.
    !          
    !          !call DetermineAdjacentParticles(Traction_Connecting_Local_To_Global_Sides(ElementUponWhichLoadIsApplied), NParticles, ParticleStatus)
    !                        
    !          call DetermineAdjacentParticles(4, NParticles, ParticleStatus)
    !    
    !          
    !          ! side number hardcoded              
    !          !call DetermineAdjacentParticles(1, NParticles, ParticleStatus)
    !    
    !          ! hardcoded 
    !            !LocalNodeConnectionArray(1) = 16
    !            !LocalNodeConnectionArray(2) = 15            
    !            !LocalNodeConnectionArray(3) = 14            
    !            !LocalNodeConnectionArray(4) = 13            
    !            !LocalNodeConnectionArray(5) = 12            
    !            !LocalNodeConnectionArray(6) = 11            
    !            !LocalNodeConnectionArray(7) = 10            
    !            !LocalNodeConnectionArray(8) = 9             
    !            !LocalNodeConnectionArray(9) = 8             
    !            !LocalNodeConnectionArray(10) = 7             
    !            !LocalNodeConnectionArray(11) = 6             
    !            !LocalNodeConnectionArray(12) = 5             
    !            !LocalNodeConnectionArray(13) = 4             
    !            !LocalNodeConnectionArray(14) = 4             
    !            !LocalNodeConnectionArray(15) = 2             
    !            !LocalNodeConnectionArray(16) = 1             
    !          
    !          
    !          NLoadedParticles = DetermineNAdjacentParticles(ParticleStatus, NParticles)
    !          
    !          !call MarkBoundaryParticles(GlobalElementForTraction, ParticleStatus, NParticles, IPatch_Temporary)
    !          do IParticle = 1, NPartEle(GlobalElementForTraction, IPatch_Temporary) ! loop all MPs of the element (solid and liquid MPs in case of double-point formulation)
    !              ParticleIndex = GetParticleIndex(IParticle, GlobalElementForTraction, IPatch_Temporary)
    !          
    !              !if ( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then
    !                  ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
    !          
    !                  ParticleID = ParticleID + 1
    !          
    !                  if(ParticleStatus(ParticleID)) then ! this MP is loaded
    !          
    !                      LoadedParticleCounter = LoadedParticleCounter + 1
    !                      
    !                      do INode = 1, (ElementBoundaryNodes_Xi*ElementBoundaryNodes_Eta)
    !                      
    !                          !NodeID = IEN_Traction(INode, ElementUponWhichLoadIsApplied)
    !                          
    !                          do IDim = 1, NDIM
    !                          
    !                          !LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1  * Factor  ) !(2/3)
    !                          !LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( SE(INode,IDim)  / real(ELEMENTBOUNDARYGAUSSPOINTS) ) !( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1  * Factor  ) !(2/3)
    !                          
    !                          LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( SE(INode,IDim) * GPShapeFunction_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) ) ! TRACTION BOOKMARK
    !                          !LocalNodeConnectionArray(
    !                          !ElementUponWhichLoadIsApplied
    !                          
    !                          ! HARDCODED HALF * 0.5
    !                          ! HARDCODED MINUS ONE FOR DOWNWARDS
    !                          
    !                          end do 
    !                          
    !                          
    !                          
    !                          
    !                          
    !                      end do 
    !                                      
    !                      call SetFExt(Particles(ParticleIndex), LoadPerParticle)
    !                      !call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle, 1) !ILoadSystem
    !          
    !                      LoadPerParticle = 0
    !                    !LoadPerParticle = DetermineLoadPerParticle(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal, Counters%NodTot, NodeCoord, Loadtype)
    !                    !
    !                    !!if ( LoadType == 1 ) then ! assign to solid phase / solid material point
    !                    !  call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
    !          
    !                    !else if ( LoadType == 2 ) then ! assign to water phase / liquid material point
    !                    !  call IncreaseFExtWater(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
    !                    !
    !                    !else if ( LoadType == 3 ) then ! assign to gas phase
    !                    !  call IncreaseFExtGas(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
    !                    !
    !                    !  if ( NFORMULATION == 2 ) then
    !                    !    call GiveError('Erroneous definition of applied external load: For double-point formulation only load on solid or liquid phase is allowed.')
    !                    !  end if
    !                    !
    !                    !end if ! load type
    !                  end if ! loaded y/n
    !              !end if ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
    !          
    !            end do ! loop over MPs
    !          deallocate(ParticleStatus, stat = IError)
    !          
    !end subroutine TransferLoadFromNodesToMP
    
      
      
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
        use ModNURBS
        !use ModMPMData
      
        implicit none
      
          ! arguments
          integer(INTEGER_TYPE), intent(in) :: NumberOfIntegrationPoints, LoadType
          integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES), intent(in) :: ILoadCon ! size (ELEMENTBOUNDARYNODES)
          !integer(INTEGER_TYPE), dimension(Counters%NodTot), intent(in) :: NDof ! size (Counters%NodTot)
          integer(INTEGER_TYPE), dimension(Counters%Sum_NodTot), intent(in) :: NDof ! size (Counters%NodTot)
    
          real(REAL_TYPE), dimension(ELEMENTBOUNDARYNODES, NVECTOR), intent(in) ::  LoadValue ! size (ELEMENTBOUNDARYNODES, NVECTOR)
          real(REAL_TYPE), dimension(Counters%Sum_NodTot, NVECTOR), intent(in) :: Coord ! size (Counters%NodTot, NVECTOR)      ! Coordinates vector
          real(REAL_TYPE), dimension(Counters%Sum_NodTot*NVECTOR), intent(inout) :: RLoad ! size(Counters%NodTot*NVECTOR)     ! Resulting load vector
      
          ! local variable
          integer(INTEGER_TYPE) :: I, J, K, NodeNumber, ND, NumberOfCornerNodes
          integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: ILoadConLocal
          
          integer(INTEGER_TYPE), allocatable, dimension(:,:) :: IConLElement
          
          real(REAL_TYPE) :: VectorLength, Temp, Factor, Radius, xLoad
          real(REAL_TYPE), dimension(NVECTOR) :: NormalVector, ProductTractionNodesNormalVector, Vector
    
          !real(REAL_TYPE), dimension(NumberOfIntegrationPoints, NVECTOR) :: TractionGaussPoint
    
          real(REAL_TYPE), allocatable, dimension(:, :, :) :: TractionGaussPoint
          real(REAL_TYPE), dimension(ELEMENTBOUNDARYNODES_XI, NVECTOR) :: TractionNodes
          real(REAL_TYPE), allocatable, dimension(:, :, :) :: SE
    
          !real(REAL_TYPE), dimension(ELEMENTBOUNDARYNODES_XI, NVECTOR) :: SE
    
          real(REAL_TYPE), dimension(ELEMENTBOUNDARYNODES_XI) :: AppliedPressure
          logical :: found = .false.
          
          integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, IElNode1, IElNode2, ElementOfNode1, ElementOfNode2, ISharedElement, INode
          real(REAL_TYPE) :: ScalarProduct 
    
          ! NURBS
          integer(INTEGER_TYPE) :: NumberOfControlPointsForNURBSTraction 
          integer(INTEGER_TYPE), allocatable, dimension(:) :: NURBSTractionNodes
          real(REAL_TYPE), allocatable, dimension(:,:) :: NURBSTractionLoadValue
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: TractionNodesElement
          integer(INTEGER_TYPE) :: NumberOfTractionElements, ee, ii, jj, kk
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          
          ! Assuming boundary condition is being applied in the xi direction 
          NumberOfCornerNodes = NURBS%NXiKnotOrder(IPatch_Temporary) + 1
          
          !do ii = 1, ELEMENTBOUNDARYNODES
              
          
          !if (NXiKnotOrder == 3) then 
          !    NumberOfCornerNodes = 4 ! hardcoded for quadratic NURBS with 3 nodes... number of corner nodes of a line element
          !elseif (NXiKnotOrder == 2) then 
          !    NumberOfCornerNodes = 3 ! hardcoded for quadratic NURBS with 3 nodes... number of corner nodes of a line element
          !elseif (NXiKnotOrder == 1) then 
          !    NumberOfCornerNodes = 2
          !end if 
          
          !ELEMENTBOUNDARYNODES = 2 ! hardcoded...
    
          ! We do not know how many control points for the traction
          
          
          !if (NURBS) then 
          !    
          !  !1) fill out the NURBSTractionNodes through a subroutine  
          !  !2) we need to do the same for the traction if each is to have its own value. For now, just use the same for all.              
          !    
          !    
          !    
          !end if 
          !NURBSTractionNodes = 4 ! for cubic
          !call FindControlPointNumbersForTractionApplication( NumberOfControlPointsForNURBSTraction, NURBSTractionNodes, ILoadCon)!, NumberOfTractionElements )
          
          allocate(IConLElement (NumberOfTractionElements, NURBS%NXiKnotOrder(IPatch_Temporary)+1) )
          allocate(TractionNodesElement (NumberOfTractionElements, NURBS%NXiKnotOrder(IPatch_Temporary)+1, NDIM) )
          allocate(SE (NumberOfTractionElements, ELEMENTBOUNDARYNODES_XI, NVECTOR) )
    
          allocate(TractionGaussPoint (NumberOfTractionElements, NumberOfIntegrationPoints, NDIM) )
            
          IConLElement = 0
          
          !if (ELEMENTBOUNDARYNODES_XI == 2) then ! 2-noded (low-order) line element
            ! no re-ordering needed
            ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
            ILoadConLocal(2) = ILoadCon(2); TractionNodes(2,:) = LoadValue(2,:)
          !else if (ELEMENTBOUNDARYNODES_XI == 3) then ! 3-noded (high-order) line element
          !  ! re-order nodes: corner nodes first, mid-nodes last
          !    ! do not re-order for NURBS implementation
          !  ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
          !  ILoadConLocal(2) = ILoadCon(2); TractionNodes(2,:) = LoadValue(2,:)
          !  ILoadConLocal(3) = ILoadCon(3); TractionNodes(3,:) = LoadValue(3,:)
            !ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
            !ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
            !ILoadConLocal(3) = ILoadCon(2); TractionNodes(3,:) = LoadValue(2,:)
          !end if      
    
          ! set arrays to zero
          TractionGaussPoint = 0.0
          SE = 0.0
    
          do ii = 1, NumberOfTractionElements
              do jj = 1, NURBS%NXiKnotOrder(IPatch_Temporary)+1
                  do kk = 1, NDIM
                      TractionNodesElement(ii, jj, kk) = TractionNodes(1,kk)
                  end do 
              end do
          end do
          
          
          
          do I = 1, NumberOfIntegrationPoints
            ! determine vector normal to a line (for 2- and 3-noded line), note that only the two corner nodes are passed
            !call NormalOnLine(I, Coord, ILoadConLocal(1:2), GPShapeFunctionDerivativeBoundary, NormalVector, VectorLength) !note that 3 is hardcoded and should depend on the order of NURBS shape function
             ! Hardcoded : IElements should be removed here
              !call NormalOnLine(I, Coord, ILoadConLocal(1:2), GPShapeFunctionDerivativeBoundary_XI_AllElements, NormalVector, VectorLength, NumberOfTractionElements, NURBSTractionNodes, IElements)!, IConLElement) !note that 3 is hardcoded and should depend on the order of NURBS shape function
            ! hardcoded: need to uncomment the NormalOnLine right above this
            
            
            ! need to change above to AllElements and select specific element
    
            Factor = 1.0 ! Factor for length of lines
    
            ! ----- MULTIPHASE -------------------------------------------------------------------------------------------------------------------
            if (LoadType == LOADTYPE_LIQUID .or. LoadType == LOADTYPE_GAS) then ! Liquid or Gas, Loading is applied only normal to the surface
                AppliedPressure(:) = LoadValue(:,1) !Pressure is hydrostatic, all components are equal LoadValue(:,1)=LoadValue(:,2)
                !*** check the normal pointing outward of the body **** 
                ! Note that this implementation works only for TRI elements. There is room for vectorization to extend to more types of element. 
                !Correct normal if necessary. LoadValue<0 means compression
                NElemNode1 = GetNElmOfNode(ILoadConLocal(1), IPatch_Temporary)
                NElemNode2 = GetNElmOfNode(ILoadConLocal(2), IPatch_Temporary)
                do IElNode1 = 1,NElemNode1
                    ElementOfNode1 = GetElmIOfNode(ILoadConLocal(1),IElNode1, IPatch_Temporary)
                    do IElNode2 = 1,NElemNode2
                        ElementOfNode2 = GetElmIOfNode(ILoadConLocal(2),IElNode2, IPatch_Temporary)
                        if (ElementOfNode1==ElementOfNode2) then !found element connected to line
                            ISharedElement = ElementOfNode2 !also = ElementOfNode1
                            if (IsActiveElement(ISharedElement, IPatch_Temporary)) then 
                             do K=1,ELEMENTNODES
                                 INode =ElementConnectivities(K,ISharedElement,IPatch_Temporary)
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
            
                        
            ! ----- MULTIPHASE -------------------------------------------------------------------------------------------------------------------
            
            do ee = 1, NumberOfTractionElements ! need to loop accross elements 
            ! calculate traction at Gauss Point
            do J = 1, NVECTOR
              do K = 1, NumberOfCornerNodes ! Corner nodes
                if ( ISAXISYMMETRIC ) then
                  NodeNumber = IConLElement(ee,K)
                  Radius = Coord(NodeNumber, 1) ! index 1 is r-direction
                else
                  Radius = 1.0
                end if
                TractionGaussPoint(ee, I,J) = TractionGaussPoint(ee, I,J) + GPShapeFunctionBoundary_XI_AllElements(ee,I,K) * TractionNodesElement(ee,K,J) * Radius
    
                !TractionGaussPoint(I,J) = TractionGaussPoint(I,J) + GPShapeFunctionBoundary_XI_AllElements(ee,I,K) * TractionNodes(K,J) * Radius
    
              end do
            end do
    
            !SE = 0.0
            
            ! calculate traction at Nodes
            
            do K = 1, NumberOfCornerNodes ! Corner nodes
            
            do J = 1, NVECTOR
              
              found = .false.
              
                ! factor = 0.5 for triangle, 1.0 for ractangle
                Temp = GPWeightBoundary_XI(I) * GPShapeFunctionBoundary_XI_AllElements(ee,I,K) * VectorLength * Factor
    
                SE(ee, K, J) = SE(ee, K, J) + Temp * TractionGaussPoint(ee, I, J)
            !    if ( ISAXISYMMETRIC ) then
            !      NodeNumber = IConLElement(ee,K)
            !      Radius = Coord(NodeNumber, 1) ! index 1 is r-direction
            !      if ( Radius < TINY ) then
            !        xLoad = SE(K, J)
            !        found = .true.
            !      end if
            !    end if
            !  end do
            !  
            !  if ( ISAXISYMMETRIC ) then
            !    if (found) then
            !      xLoad = (xLoad / real((NumberOfCornerNodes -1), REAL_TYPE)) * (7.0 / 8.0)
            !      do K = 1, NumberOfCornerNodes ! Corner nodes
            !        NodeNumber = IConLElement(ee,k)
            !        Radius = Coord(NodeNumber, 1) ! index 1 is r-direction
            !        if ( Radius > TINY ) then
            !          SE(K, J) = SE(K, J) + xLoad
            !        else
            !          SE(K, J) = SE(K, J) - xLoad
            !        end if
            !      end do
            !    end if
            !  end if
            !  
            end do
            end do
            
            
            
    
          ! add to global load vector
          
          
            
          
          
            end do
            
          end do
          do ee = 1, NumberOfTractionElements ! need to loop accross elements 
          do II = 1, ELEMENTBOUNDARYNODES_XI
            ND = NDof(IConLElement(ee,II))
            
            do J = 1, NVECTOR
              RLoad(ND + J) = RLoad(ND + J) + SE(ee, II,J)
            
            end do
          
          
          end do
          end do
          
    
    end subroutine Load2D
    
    
    
    
    !subroutine FindControlPointNumbersForTractionApplication( NumberOfControlPointsForNURBSTraction, NURBSTractionNodes, ILoadCon )
    !!call FindControlPointNumbersForTractionApplication( NumberOfControlPointsForNURBSTraction, NURBSTractionNodes, ILoadCon )
    !!use ModCounters
    !!  use ModElementEvaluation
    !!  use ModGlobalConstants
    !!  use ModMeshAdjacencies
    !!    use ModMeshInfo
    !
    !
    !implicit none
    !
    !! initialize variables 
    !! input 
    !integer(INTEGER_TYPE), intent(inout) :: NumberOfControlPointsForNURBSTraction
    !integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: NURBSTractionNodes
    !integer(INTEGER_TYPE), dimension(:), intent(in) :: ILoadCon
    !! local 
    !integer(INTEGER_TYPE) :: MinimumControlPointNumber, MaximumControlPointNumber
    !integer(INTEGER_TYPE) :: ii
    !
    !
    !MinimumControlPointNumber = ILoadCon(1) ! 1
    !MaximumControlPointNumber = ILoadCon(2) ! 3
    !
    !NumberOfControlPointsForNURBSTraction = MaximumControlPointNumber - MinimumControlPointNumber + 1
    !
    !allocate(NURBSTractionNodes ( NumberOfControlPointsForNURBSTraction ) )
    !
    !!NURBSTractionNodes(1) = MinimumControlPointNumber
    !
    !do ii = 1, NumberOfControlPointsForNURBSTraction ! start from 1 
    !    
    !    NURBSTractionNodes(ii) = MinimumControlPointNumber
    !    
    !    MinimumControlPointNumber = MinimumControlPointNumber + 1
    !
    !end do
    !
    !    
    !    
    !!    NumberOfControlPointsForNURBSTraction =  NumberOfControlPointsForNURBSTraction + 1
    !!   
    !!   if (NumberOfControlPointsForNURBSTraction>=MaximumControlPointNumber) then 
    !!       
    !!       exit 
    !!       
    !!       end if
    !!
    !!!end do 
    !
    !
    !
    !
    !
    !end subroutine FindControlPointNumbersForTractionApplication  
    
    
      
      
      subroutine Load3D(RLoad, NDof, Coord, NInt, NNod, ILoadCon, ElementUponWhichLoadIsApplied, ILoadSystem, IPatch, &
          NXiKnotEntries_SideSpecific, XiKnotEntries_SideSpecific, &
          NEtaKnotEntries_SideSpecific, EtaKnotEntries_SideSpecific)
      !call Load3D(ExtLoadVector, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, NNodes, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid_NURBS(:, I), LoadType)
      !**********************************************************************
      !
      !>    Function:  Calculate external load vector for NURBS element. 
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
        integer(INTEGER_TYPE), intent(in) :: NInt, NNod!, LoadType 
        integer(INTEGER_TYPE), dimension(NNod) :: ILoadCon
        integer(INTEGER_TYPE), dimension(Counters%Sum_NodTot) :: NDof
        !real(REAL_TYPE), dimension(NNod, NVECTOR), intent(in) ::  LoadValue
        real(REAL_TYPE), dimension(Counters%Sum_NodTot, NVECTOR), intent(in) :: Coord       ! Coordinates vector !Counters%Sum_NodTot !NURBS_NControlPoints_Traction(ILoadSystem,IPatch)
        real(REAL_TYPE), dimension(Counters%Sum_NodTot*NDOFL), intent(inout) :: RLoad      ! Resulting load vector
        integer(INTEGER_TYPE), intent(in) :: ElementUponWhichLoadIsApplied
      
        ! local variable
        integer(INTEGER_TYPE) :: I, J, K, NN, ND, NNodLOE
        integer(INTEGER_TYPE), dimension(NNod) :: ILoadConLocal
        real(REAL_TYPE) :: VectorLength, Temp
        real(REAL_TYPE), dimension(NVECTOR) :: NormalVector, ProductTractionNodesNormalVector, Vector
        real(REAL_TYPE), dimension(NInt, NVECTOR) :: TractionGaussPoint
        real(REAL_TYPE), dimension(NNod, NVECTOR) :: TractionNodes, SE
        real(REAL_TYPE), dimension(NNod) :: AppliedPressure
        
        real(REAL_TYPE) :: RLoad_DEBUG
        real(REAL_TYPE) :: SE_DEBUG
        
        integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, NElemNode3, IElNode1, IElNode2, IElNode3, ElementOfNode1, ElementOfNode2, ElementOfNode3, ISharedElement, INode
        real(REAL_TYPE) :: ScalarProduct
        
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary
        integer(INTEGER_TYPE) :: ILoadSystem_Temporary
        integer(INTEGER_TYPE), intent(in) :: IPatch
        integer(INTEGER_TYPE), intent(in) :: ILoadSystem
         
        integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries_SideSpecific
        real(REAL_TYPE), dimension(NXiKnotEntries_SideSpecific), intent(in) :: XiKnotEntries_SideSpecific
        
        
        integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntries_SideSpecific
        real(REAL_TYPE), dimension(NEtaKnotEntries_SideSpecific), intent(in) :: EtaKnotEntries_SideSpecific
        
        
        
        RLoad_DEBUG = 0.0
        SE_DEBUG = 0.0
        
        ! set arrays to zero
        TractionGaussPoint = 0d0
        SE = 0d0
      
        do I = 1, ELEMENTBOUNDARYGAUSSPOINTS
            
            call NormalOnQUAD_NURBS(NodalCoordinates_Traction(:,:, ILoadSystem, IPatch), GPShapeFunctionDerivative_AllElements(ElementUponWhichLoadIsApplied,I,:,:, ILoadSystem, IPatch), &
                NormalVector, VectorLength, ElementUponWhichLoadIsApplied, &
                ILoadSystem, IPatch, &
                NXiKnotEntries_SideSpecific, XiKnotEntries_SideSpecific, &
                NEtaKnotEntries_SideSpecific, EtaKnotEntries_SideSpecific)
          
          ! calculate traction at Nodes
          do K = 1, ELEMENTBOUNDARYNODES!NNodLOE ! Corner nodes
            do J = 1, NVECTOR  
                
                Temp = GPWeight_AllElements(ElementUponWhichLoadIsApplied,I, ILoadSystem, IPatch) * GPShapeFunction_AllElements(ElementUponWhichLoadIsApplied, I, K, ILoadSystem, IPatch) * VectorLength * NormalVector(J) & 
                    * -1 
                ! hardcoded negative sign to make traction go inwards
                SE(K,J) = SE(K,J) + ( Temp )!* TractionGaussPoint(K, J) )
          
                ! Normal outside factor --> -1 if it is pointing inwards 
                ! Normal outside factor -->  1 if it is pointing outwards
              
            end do 
          end do
          
          
        end do

        ! add to global load vector
        do I = 1, ELEMENTBOUNDARYNODES!NNod  ! --> ELEMENTBOUNDARYNODES is the number of nodes in the element (eg. 2*2 if traction on a linear square element)
          ND = NDOF(Traction_Connecting_Local_To_Global_ControlPoints(NURBS%ElementConnectivities_TRACTION(I,ElementUponWhichLoadIsApplied, ILoadSystem, IPatch), ILoadSystem, IPatch))!NDof(ILoadConLocal(I))
          do J = 1, NVECTOR
            RLoad(ND + J) = RLoad(ND + J) + SE(I,J)    
            RLoad_DEBUG = RLoad_DEBUG + RLoad(ND + J)
            SE_DEBUG = SE_DEBUG + (SE(I,J))**2
          end do
        end do

    end subroutine Load3D
      

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    