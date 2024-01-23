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
	  
	  
	  module ModMPMInit
      !**********************************************************************
      !
      !  Function : Contains the initialisation routines for the housekeeping data structure
      !               
      !             This module should only contain routines that are directly related to
      !             the initialisation of material point data.
      !
      !     $Revision: 9728 $
      !     $Date: 2022-04-18 22:49:38 +0200 (ma, 18 apr. 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluation
      use ModMeshAdjacencies
      use ModMPMData
      use ModConvectivePhase
      use ModReadMPMData
      use ModWriteMPMData
      use ModMeshInfo
      use ModDynViscousBoudary
      use ModEmptyElements
      use ModAdjustParticleDiscretisation
      use ModMPMDynContact
      use ModGlobalConstants
      use ModReadGeometryData
      
      implicit none

      contains


        subroutine InitialiseMaterialPointHousekeeping()
        !**********************************************************************
        !
        !  Function : Initialises the list of material points
        !
        !             Material point information is initialised in case of the calculation of an initial phase.
        !             If a 'follow-up' phase is calculated, material point information is reinitialised
        !             from data from the BRF file of the previous load phase. 
        !
        !**********************************************************************

        implicit none

          call DestroyHouseKeeping()

          if (.not.IsFollowUpPhase()) then ! Initial phase: initialise material point data

            call InitialiseMaterialPoints() ! Initialise material point data structure

          else ! Follow-up phase: reinitialise material point data from the BRF file of the previous load phase

            call ReinitialiseMaterialPoints() ! Reinitialise material point data structure
            CalParams%ApplyPrescribedDisplacements = ReadPPD() ! Initialise prescribed displacements if the PPD file exists

          end if

          call SetActiveElement() ! allocate and assign array ActiveElement !hardcoded

          FirstTimeGetParticleIndex = .false.
          call SetParticleIndex()

        end subroutine InitialiseMaterialPointHousekeeping
        

        subroutine InitialiseMaterialPoints()
        !**********************************************************************
        !
        !  Function : Initialise list of material points and housekeeping arrays.
        !             This is where we initialize the shape functions for the first time. 
        !
        !**********************************************************************

        implicit none

          ! local variables
          integer(INTEGER_TYPE) :: IElement, ParticleIndex, IParticle, Index, EntityID, IDCounter
          integer(INTEGER_TYPE), dimension(maximum_nel_NURBS, NumberOfPatches) :: InitialSolidMaterialPoints, InitialLiquidMaterialPoints !Counters%NEl
          real(REAL_TYPE) :: IntegrationWeight, DetJac, LocPos(NVECTOR), DampingFactor=0.0
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: RJac, InvRJac
          integer(INTEGER_TYPE) :: MaterialPointType
          
          real(REAL_TYPE), allocatable, dimension(:) :: HS
          real(REAL_TYPE), allocatable, dimension(:,:) :: dHS 
          
          !real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: dxi_dx
          !real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: dtildexi_dxi
          
          ! NURBS jacobian parameters
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: dxi_dx
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: dtildexi_dxi
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: dx_dxi
          real(REAL_TYPE) :: Det_dx_dtildexi
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: dx_dtildexi
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: dtildexi_dx
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          
          InitialSolidMaterialPoints = 0
          InitialLiquidMaterialPoints = 0
          IDCounter = 0


            if((CalParams%ComputationMethod==FEM).or.(CalParams%ComputationMethod==UL_FEM)) then !for FEM calculation
              CalParams%NMaterialPoints = 1
              InitialSolidMaterialPoints = 1
              InitialLiquidMaterialPoints = 0
            else ! for MPM calculation, as specified in GOM file
              call ReadNumberOfMaterialPointsPerElement(InitialSolidMaterialPoints, InitialLiquidMaterialPoints)
              call CheckNumberOfMaterialPointsPerElement(InitialSolidMaterialPoints, InitialLiquidMaterialPoints)
            end if
            
            ! determine initial number of material points: fill NPartEle, NSolidEle, NLiquidEle, determine Counters%NParticles
            call InitialiseNumberMaterialPoints(InitialSolidMaterialPoints, InitialLiquidMaterialPoints)

          ! initialise material point housekeeping arrays:
          ! EleParticlesHelp, EleParticles, Particles, IsParticleIntegration (global variables in ModMPMData)
          call InitialiseHouseKeepingArrays()

          if (CalParams%ApplyEmptyElements.or.CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix) then
            call InitialiseModEmptyElementsData()
          end if

          ! determine decimal factor for material point housekeeping arrays: PartFactor
          call DetermineDecimalFactor()

          ParticleIndex = 0

          do IPatch = 1, NumberOfPatches ! loop over patches
          do IElement = 1, nel_NURBS(IPatch)!Counters%NEl ! loop over all elements
            if (IsActiveElement(IElement,IPatch)) then ! only if active element

              ! loop over total number of material points (solid+liquid) per element 
              ! -> we also calculate their shape functions here in this loop 
              do IParticle = 1, NSolidEle(IElement,IPatch) + NLiquidEle(IElement,IPatch)
                ParticleIndex = ParticleIndex + 1 ! counting from 1 to Counters%NParticles

              
                call InitialLocalMaterialPointCoordinates(IParticle, NSolidEle(IElement,IPatch), NLiquidEle(IElement,IPatch), IntegrationWeight, LocPos)


                !evaluate shape functions here and input into the jacobian evaluation   
                call ShapeFunctionData(LocPos, ELEMENTNODES, HS, dHS, IElement, IPatch)!NNodes, 
              
              
                ! determine the determinant of the Jacobian
                !call DetJacob(HS, dHS, &
                !LocPos, Counters%NEl, Counters%NodTot, NDIM, &
                !    IElement, ElementConnectivities, NodalCoordinates, RJac, InvRJac, DetJac, &
                !  dxi_dx, dtildexi_dxi)
                !call DetJacob(HS, dHS, &
                !LocPos, Counters%NEl, Counters%NodTot, NDIM, &
                !    IElement, ElementConnectivities(:,:, IPatch), NodalCoordinates(:,:,IPatch), dx_dxi, dx_dtildexi, dtildexi_dx, Det_dx_dtildexi, &
                !  dxi_dx, dtildexi_dxi, &
                !    IPatch)
                
                call DetJacob(HS, dHS, &
                    LocPos, nel_NURBS(IPatch), NControlPoints(IPatch), NDIM, &
                    IElement, ElementConnectivities(:,1:nel_NURBS(IPatch), IPatch), NodalCoordinates(1:NControlPoints(IPatch),:,IPatch), dx_dxi, dx_dtildexi, dtildexi_dx, Det_dx_dtildexi, &
                    dxi_dx, dtildexi_dxi, &
                    IPatch)

                ! determine global integration weight of material point
                !IntegrationWeight = IntegrationWeight * DetJac
                IntegrationWeight = IntegrationWeight * Det_dx_dtildexi

                
                if (GeoParams%ApplyLocalDampingElement) then
                  DampingFactor = GeoParams%LocalDampingFactorElement(IElement) ! only for absorbing boundaries
                end if  
                
                ! determine type of material point
                if (NFORMULATION==1) then ! 1-layer formulation, material point is of type MIXTURE
                  if (MatParams(ElementMaterialID(IElement,IPatch))%MaterialType=='1-phase-liquid'.or.MatParams(ElementMaterialID(IElement,IPatch))%MaterialPhases=='1-phase-liquid') then
                    MaterialPointType = MaterialPointTypeLiquid
                  else
                    MaterialPointType = MaterialPointTypeMixture
                  end if
                else ! 2-layer formulation, material point is of either type SOLID or LIQUID
                  if ((IParticle<=NSolidEle(IElement,IPatch)).and.(NSolidEle(IElement,IPatch)>0)) then ! material point is of type SOLID
                    MaterialPointType = MaterialPointTypeSolid
                    IntegrationWeight = ElementSpace(IElement,IPatch) / NSolidEle(IElement,IPatch) !Temporary!!!
                  else ! material point is of type LIQUID
                    MaterialPointType = MaterialPointTypeLiquid
                    IntegrationWeight = ElementSpace(IElement,IPatch) / NLiquidEle(IElement,IPatch) !Temporary!!!
                  end if
                end if

                ! set IEntityPile = HARD_ENTITY otherwise SOFT_ENTITY
                
                if (allocated (GlobContElement)) then                 
                  if (GlobContElement(IElement)) then
                    EntityID = HARD_ENTITY
                  else 
                    EntityID = SOFT_ENTITY 
                  end if
                else 
                  EntityID = SOFT_ENTITY
                end if

                ! assign (basic) data to the material point
                Particles(ParticleIndex) = InitParticle(abs(ElementMaterialID(IElement,IPatch)), & ! Material ID
                                           MaterialPointType, & ! MaterialPointType: mixture, solid, liquid, 
                                           IntegrationWeight, & ! integration weight
                                           LocPos, & ! local coordinate
                                           MATERIALPARTICLE, & ! material or virtual material point
                                           -1.0, & !  mass assigned to material point in case of virtual material point
                                           CalParams%GravityData%GAccel, & ! gravity acceleration
                                           CalParams%GravityData%GravityVector, & ! direction-vector of gravitational force
                                           DampingFactor, & ! element based damping factor
                                           CalParams%ApplySubmergedCalculation,&
                                           MassArray(ParticleIndex), &
                                           MassWaterArray(ParticleIndex)) ! submerged calculation
                !New ElementIDArray
                ElementIDArray(ParticleIndex) = IElement
                
                PatchIDArray(ParticleIndex) =  IPatch
                
                IDCounter = IDCounter + 1
                IDArray(ParticleIndex) = IDCounter
                EntityIDArray(ParticleIndex) = EntityID
                MaterialIDArray(ParticleIndex) = abs(ElementMaterialID(IElement,IPatch))
                MaterialPointTypeArray(ParticleIndex) = MaterialPointType

                ! initialise the state parameters for the ICH and HP models in case they are used
                call AllocateMaterialModelsStateVariables(ParticleIndex)

                ! set material point shape function values and derivatives
                call SetParticleShapeFunctionData(Particles(ParticleIndex),ParticleIndex, IElement, ElementConnectivities(:,:, IPatch), NodalCoordinates(:,:,IPatch), IPatch) ! -> this is the first time we call the shape functions in the code

                ! set material point housekeeping arrays EleParticle and EleParticleHelp
                call SetEleParticles(IElement, IParticle, ParticleIndex, IPatch)
                
                
                !call CoordLocalToGlobal(IElement, NodalCoordinates(:,:,IPatch), IPatch)

              end do ! loop over material points

              ! determine global position of material points from local material point coordinates
              ! Abdelrahman Alsardi: there was a bug here which I fixed using "1:NControlPoints(IPatch)" 
              call CoordLocalToGlobal(IElement, NodalCoordinates(1:NControlPoints(IPatch),:,IPatch), IPatch)
              !call CoordLocalToGlobal(IElement, GlobalMultipatchNodalCoordinates, IPatch)


              if ( ISAXISYMMETRIC ) then ! update mass and integration weight
                do IParticle = 1, NPartEle(IElement,IPatch)
                  Index = ParticleIndex - NPartEle(IElement,IPatch) + IParticle
                  call correctParticleMassAndBodyForceAxisymmetric(Index)
                  call CorrectParticleIntegrationWeightAxisymmetric(Index)
                end do
              end if  
              
            else ! element is inactive, no material points in inactive element

              ! set material point housekeeping arrays EleParticle and EleParticleHelp
              call SetEleParticles(IElement, IParticle, -1, IPatch)
              
            end if
            
          end do ! loop over all elements
          end do ! loop over patches

          ! initialise Particles(I)%IsBoundaryParticle flag
          call DetermineBoundaryParticles()
                    
        end subroutine InitialiseMaterialPoints


        subroutine ReinitialiseMaterialPoints()
        !**********************************************************************
        !
        !    Function:  Reinitialises the list of particles and the house-keeping arrays
        !               from particle data from the BRF file of the previous load phase.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IError
          character(len = MAX_FILENAME_LENGTH) :: FileName
        
          FileName = trim(CalParams%FileNames%ProjectName)//BRF_FILE_EXTENSION//trim(CalParams%FileNames%PreviousStepExt)

          ! Read number of particles from the BRF file of the previous load phase
          call ReadHouseKeepingCountersFromFile(FileName)

          ! Initialise particle house-keeping arrays with the counters from the BRF file of the previous load phase
          call InitialiseHouseKeepingArrays()
          if (CalParams%ApplyEmptyElements.or.CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix) then
            call InitialiseModEmptyElementsData()
          end if
          
          allocate(NPartEle(Maximum_nel_NURBS, NumberOfPatches), stat = IError) !Counters%NEl
          NPartEle = 0 ! Initialise to zero

          ! Fill particle house-keeping arrays with data from the BRF file of the previous load phase
          call ReadHouseKeepingDataFromFile(FileName)
          if (IsMPMWithMPIntegration()) then !Correction is required if change from MIXED- to MP-integration is considered in latest CPS
            IsParticleIntegration = .true.
          end if      
          
          ! Fill particle list with particle data from the BRF file of the previous load phase
          call ReadParticleDataFromFile(FileName)

          ! Fill particle list with state parameters from the BRF file of the previous load phase
          call ReadParticleStateParametersFromFile(FileName)

          ! Fill particle list with dynamic related information from the BRF file of the previous load phase
          call ReadParticleDynamicDataFromFile(FileName)

          ! Fill the particle information related to consolidation from the BRF file of the previous load phase
          call ReadConsolidationDataFromFile(FileName)

          if(CalParams%NumberOfPhases==3) then
            ! Fill the particle information related to unsaturated material from the BRF file of the previous load phase
            call ReadConsolidationDataGasFromFile(FileName)
          end if

          ! rewrite the particle material parameters from new GOM file 
          ! if it is not the first phase and ApplyMaterialUpdate is true
          if (CalParams%ApplyMaterialUpdate) then
            call RewriteParticleMaterialParameters()
          end if
          
          

          !ActiveElement
          
          ! Adapt IsActiveElement to list of elements which contain particles
          call AdaptIsElemToParticleData()

          ! Set the particle IDCounter to the ID of the last read particle
          call SetIDCounter(MAXVAL(IDArray))

          if (.not.allocated(ActiveElement)) then
             call SetActiveElement() !hardcoded
          endif
          
          
          ! Set particle shape function values and derivatives
          call InitialiseParticleShapeFunctionData()
          
          !First

          call SetParticleIndex() ! BOOKMARK
          
          if (CalParams%ApplyEmptyElements.or.CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix) then
            ! Check whether elements are deactivated that should be considered fully filled
            call CheckEmptyElements()
          end if
          if (CalParams%ApplyEmptyElements) then
            call AdjustParticleDiscretisation()
          end if
          
          if (CalParams%ApplyContactAlgorithm) then
            call ReAssignEntityID()
          end if
        
        end subroutine ReinitialiseMaterialPoints
        
        
        subroutine ReAssignEntityID()        
        !**********************************************************************
        !
        !    Function:  Reset EntityID
        !
        !**********************************************************************
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: ParticleIndex, IElement    
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          
          do ParticleIndex = 1, Counters%NParticles
            
            IElement = ElementIDArray(ParticleIndex) 
            
            IPatch = PatchIDArray(ParticleIndex)
            
            if (allocated (GlobContElement)) then                 
              if (GlobContElement(IElement)) then ! set IEntityPile = HARD_ENTITY otherwise SOFT_ENTITY
                EntityIDArray(ParticleIndex) = HARD_ENTITY
              else 
                EntityIDArray(ParticleIndex) = SOFT_ENTITY
              end if
            else 
              EntityIDArray(ParticleIndex)= getMaterialEntity(abs(ElementMaterialID(IElement,IPatch_Temporary)))
            end if          
          end do 
          
        end subroutine ReAssignEntityID
        
        
        subroutine ResetMaterialPointDisplacements()       
        !**********************************************************************
        !
        !    Function:  If CalParams%ApplyResetDisplacements = .true. then, 
        !               displacements accumulated at the material points are set to zero
        !
        !**********************************************************************
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: ParticleIndex
          
          if(.not.CalParams%ApplyResetDisplacements) RETURN     
          UArray(:,:)=0d0
          UStepArray(:,:) = 0d0
          UPhaseArray(:,:) =0d0
          do ParticleIndex = 1, Counters%NParticles
            Particles(ParticleIndex)%UW=0
            Particles(ParticleIndex)%UG=0
            Particles(ParticleIndex)%UStepWater=0
            
          end do         

        end subroutine ResetMaterialPointDisplacements
        
        
        subroutine ReinitialiseUpdatedNodes()
        !**********************************************************************
        !
        !    Function:  If the current phase is a follow-up phase, the updated initial
        !               nodal coordinates of the previous phase are read from the BRF 
        !               file of the previous load phase if mesh adjustment was applied 
        !               in the previous phase.
        !
        !**********************************************************************

        implicit none

          ! Local variables
          character(len = MAX_FILENAME_LENGTH) :: FileName
          integer(INTEGER_TYPE) :: I, IDoF, J
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

          if (.not.IsFollowUpPhase()) RETURN
          
          if (IsULFEMComputation()) then
            NodalCoordinatesUpd = NodalCoordinates

            do IPatch = 1, NumberOfPatches ! loop over patches
            
            do I = 1, NControlPoints(IPatch)!Counters%NodTot ! loop over control points
                
                GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
                
              IDoF = ReducedDoF(GlobalNodeID)!I)
              do J = 1, NVECTOR
                NodalCoordinatesUpd(I, J, IPatch) = NodalCoordinatesUpd(I, J, IPatch) + TotalDisplacementSoil(IDoF + J)
              end do
            end do
            
            end do ! loop over patches

          else
            FileName = trim(CalParams%FileNames%ProjectName)//BRF_FILE_EXTENSION//trim(CalParams%FileNames%PreviousStepExt)
        
            do IPatch = 1, NumberOfPatches ! loop over patches
            !call ReadNodalCoordinates(FileName, Counters%NodTot, NDIM, NodalCoordinates(:,:,IPatch))
            call ReadNodalCoordinates(FileName, NControlPoints(IPatch), NDIM, NodalCoordinates(:,:,IPatch))
        
            NodalCoordinatesUpd = NodalCoordinates
        
            call UpdateMeshAdjacencyInformation(NodalCoordinates(:,:,IPatch))
            end do 
          end if
          
        end subroutine ReinitialiseUpdatedNodes


        subroutine InitialiseParticleShapeFunctionData()
        !**********************************************************************
        !
        !    Function:  Initialises the particle shape function values and shape
        !               function derivatives.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
        integer(INTEGER_TYPE) :: ElementID
        
          ! Local variables
          integer(INTEGER_TYPE) :: ParticleIndex
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: PatchID
        
          do ParticleIndex = 1, Counters%NParticles
              ! NURBS
              ElementID = ElementIDArray(ParticleIndex)  
              PatchID = PatchIDArray(ParticleIndex)
            call SetParticleShapeFunctionData(Particles(ParticleIndex),ParticleIndex,ElementID, ElementConnectivities(:,:, PatchID), NodalCoordinates(:,:,PatchID), PatchID) ! This needs to be done properly here!!!
          end do
        
        end subroutine InitialiseParticleShapeFunctionData

        
        subroutine AdaptIsElemToParticleData()
        !**********************************************************************
        !
        !    Function:  Set IsActiveElement(IElement) to .true. or .false. (active or inactive)
        !               depending on whether particles are located inside element IElement.
        !               Used for reinitialising the mesh configuration from particle
        !               data read from the BRF file of the previous load phase.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IElement
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch

          do IPatch = 1, NumberOfPatches ! loop over patches
          do IElement = 1, nel_NURBS(IPatch)!Counters%NEl
            if (NPartEle(IElement,IPatch)>0) then ! There is a particle inside the element
              IsActiveElement(IElement,IPatch) = .true.
              
              NAEl_NURBS(IPatch) = NAEl_NURBS(IPatch) + 1 
              
            else
              IsActiveElement(IElement,IPatch) = .false.
            end if
          end do ! loop over elements
          end do ! loop over patches
      
        end subroutine AdaptIsElemToParticleData

        
        subroutine InitialiseNumberMaterialPoints(InitialSolidMaterialPoints, InitialLiquidMaterialPoints)
        !**********************************************************************
        !
        !  Function : Determines the initial number of material points for each element
        !             stored in NPartEle(IDElement), NSolidEle(IDElement), NLiquidEle(IDElement)  
        !             and number of solid and liquid material points 
        !             stored in Counters%SolidMaterialPoints, Counters%LiquidMaterialPoints
        !  
        !  I  InitialSolidMaterialPoints : Initial number of solid material points per element
        !  I  InitialLiquidMaterialPoints : Initial number of liquid material points per element
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: InitialSolidMaterialPoints ! rank 2 for multipatch implementation
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: InitialLiquidMaterialPoints
 
          ! local variables
          integer(INTEGER_TYPE) :: ElementID, IError
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch

          !allocate(NPartEle(Counters%NEl), stat = IError)
          allocate(NPartEle(Maximum_nel_NURBS, NumberOfPatches), stat = IError)
          !allocate(NSolidEle(Counters%NEl), stat = IError)
          allocate(NSolidEle(Maximum_nel_NURBS, NumberOfPatches), stat = IError)
          !allocate(NLiquidEle(Counters%NEl), stat = IError)
          allocate(NLiquidEle(Maximum_nel_NURBS, NumberOfPatches), stat = IError)

          
          NPartEle = 0
          NSolidEle = 0
          NLiquidEle = 0
          Counters%NParticles = 0
          Counters%SolidMaterialPoints = 0
          Counters%LiquidMaterialPoints = 0
          
          do IPatch = 1, NumberOfPatches ! loop over patches
          do ElementID = 1, nel_NURBS(IPatch)!Counters%NEl ! loop over elements 
            if (IsActiveElement(ElementID,IPatch)) then ! only if active element
              NSolidEle(ElementID, IPatch) = InitialSolidMaterialPoints(ElementID, IPatch)
              NLiquidEle(ElementID, IPatch) = InitialLiquidMaterialPoints(ElementID, IPatch)
              NPartEle(ElementID, IPatch) =  NSolidEle(ElementID, IPatch) + NLiquidEle(ElementID, IPatch)
              ! add initial number of material points of ElementID
              Counters%SolidMaterialPoints = Counters%SolidMaterialPoints + NSolidEle(ElementID, IPatch)
              Counters%LiquidMaterialPoints = Counters%LiquidMaterialPoints + NLiquidEle(ElementID, IPatch)
              Counters%NParticles = Counters%SolidMaterialPoints + Counters%LiquidMaterialPoints
            end if
          end do !elements
          end do ! patches

        end subroutine InitialiseNumberMaterialPoints

        
        subroutine DetermineBoundaryParticles()
        !**********************************************************************
        !
        !  Function:  Determines which material points of the initial discretisation lie
        !             on the boundary of the body of activated elements. The material point
        !             status is stored with the flag IsBoundaryParticle.
        !             Material points next to the boundary of the mesh are not considered as
        !             boundary material points, only those next to deactivated elements.
        !
        !**********************************************************************
        
        implicit none 
        
          ! local variables
          logical:: ParticleStatus (:)
          allocatable :: ParticleStatus 
          integer(INTEGER_TYPE) :: IElement, ISide, BoundarySurface, NParticles, IError

          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          
          NParticles = 0

          do IPatch = 1, NumberOfPatches ! loop over patches
          do IElement = 1, nel_NURBS(IPatch) !Counters%NEl ! loop all elements
            if (IsActiveElement(IElement,IPatch)) then ! only initially active element

              if ( NFORMULATION == 1 ) then
                NParticles =  NPartEle(IElement, IPatch)
              else if ( NFORMULATION == 2 ) then
                NParticles =  NSolidEle(IElement, IPatch) ! only solid material points are considered
                if ( NParticles == 0 ) CYCLE ! an active element which is only filled with liquid material points is skipped
              else
                call GiveError("Formulation has to be single-point or double-point in [subroutine DetermineBoundaryParticles()].")
              end if
             
              allocate(ParticleStatus(NParticles), stat = IError)
            
              ParticleStatus = .false. ! No particles of element marked as boundary particle
            
              do ISide = 1, ELEMENTSIDES ! Loop over all sides of active elements

                ! Determine whether ISide lies on the boundary of the active elements
                BoundarySurface = BoundaryElementSurface(IElement, ISide, IsActiveElement(:,IPatch), Counters%NEl)
              
                if (BoundarySurface==1) then
                  ! Surface between activated and deactivated elements:
                  ! (0 ... inside, -1 ... boundary of mesh, -999 ... something went wrong)
                  ! Determine which particles lie next to the detected boundary side
                  call DetermineAdjacentParticles(ISide, NParticles, ParticleStatus)
                end if
              end do
              
              ! Set flag IsBoundaryParticle for particles of IElement
              call MarkBoundaryParticles(IElement, ParticleStatus, NParticles, IPatch)
              
              deallocate(ParticleStatus, stat = IError)
              
            end if ! active element
          end do ! loop all elements
          end do ! patches
        
        end subroutine DetermineBoundaryParticles


        subroutine MarkBoundaryParticles(IElement, ParticleStatus, NParticles, IPatch)
        !**********************************************************************
        !
        !    Function:  Transfer status from ParticleStatus to Particle(I)%IsBoundaryParticle.
        !
        !     IElement : ID of the considered element
        !     ParticleStatus : For each particle of the considered element in its initial configuration
        !                      .true. if the particle is a boundary particle, else .false.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IElement, NParticles
          logical, dimension(NParticles), intent(in) :: ParticleStatus

          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex
          
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
          do IParticle = 1, NParticles ! Loop over all particles of element (initial configuration)
            ParticleIndex = GetParticleIndexFunction(IParticle, IElement, IPatch)
            Particles(ParticleIndex)%IsBoundaryParticle = ParticleStatus(IParticle)
          end do
        
        end subroutine MarkBoundaryParticles
        
        
        
        
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
        
        
        
        subroutine TransferExternalLoadsToParticles_NURBS(LoadType, NSurfaceNodes, IElTyp,   NGaussPoints, ILoadSystem)
        !(LoadType, NSurfaceNodes, IElTyp, IConGlobal, LoadsGlobal, NGaussPoints, NodeCoord,ILoadSystem)
        !**********************************************************************
        !
        !    Function:  Transfer external distributed loads from element surfaces to
        !               particles closest to the considered surface - all particles get
        !               the same share of the distributed load.
        !               So far, only 6-noded triangular element sides are considered with
        !               constant distributed loads.
        !               If the loaded surface lies between two activated elements,
        !               particles of both elements get an equal share of the load.
        !
        !     LoadType : 1 solid 2 water 3 gas
        !     NSurfaceNodes : Number of nodes of the surface element
        !     IElTyp : Number of nodes per element
        !     IConGlobal : Connectivities of the surface element
        !     LoadsGlobal : Nodal load values of the surface element
        !     NGaussPoints : Number of Gauss points of the surface element
        !     NodeCoord : Nodal coordinates
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NSurfaceNodes, NGaussPoints, IElTyp, LoadType, ILoadSystem !
          integer(INTEGER_TYPE), dimension(NSurfaceNodes) :: IConGlobal !, intent(in)
          real(REAL_TYPE), dimension(NSurfaceNodes, NVECTOR) :: LoadsGlobal !, intent(in)
          real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR) :: NodeCoord !, intent(in)
        
          ! local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, IElement, NElements, NLoadedParticles, NParticles, IError, ParticleID
          integer(INTEGER_TYPE) :: ElementNumber, SideNumber
          integer(INTEGER_TYPE), dimension(2) :: IElements, ISides
          integer(INTEGER_TYPE), dimension(NSurfaceNodes) :: IConLocal
          integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElCon
          real(REAL_TYPE), dimension(NSurfaceNodes, NVECTOR) :: LoadsLocal
          logical, dimension(:), allocatable :: ParticleStatus
          real(REAL_TYPE), dimension(NVECTOR) ::  LoadPerParticle !TotalSurfaceLoad,
        
          real(REAL_TYPE), dimension(NVECTOR, 4) :: TotalSurfaceLoad
          
          integer(INTEGER_TYPE) :: NumberOfTractionElements
          
          integer(INTEGER_TYPE), dimension(:), allocatable :: NURBSTractionNodes
          
          integer(INTEGER_TYPE) :: ii
          
          ! if CornerControlPoints is 1 then it is a corner node 
          integer(INTEGER_TYPE), dimension(4) :: CornerControlPoints ! this should be cubic order + 1 = 3 + 1 = 4
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          
          
          
          
          
          
          
          
          !!CornerControlPoints = 0
          !!CornerControlPoints(1) = 1 ! 1 is a corner node 
          !!CornerControlPoints(4) = 1 ! 4 is a corner node with nothing to the 
          
          ! copy load connectivities to local arrays
          !select case(NSurfaceNodes)
          !  case(6) ! 6-noded triangular element, 3D case
          !    call RearrangeConnectivitiesTRI6(IConGlobal, LoadsGlobal, IConLocal, LoadsLocal)
          !  case(2) ! 2-noded line element, 2D case
          !    call RearrangeConnectivitiesLINE2(IConGlobal, LoadsGlobal, IConLocal, LoadsLocal)
          !  case(3) ! 3-noded line element, 2D case
          !    call RearrangeConnectivitiesLINE3(IConGlobal, LoadsGlobal, IConLocal, LoadsLocal)
          !  end select
            
            IConLocal = IConGlobal 
            LoadsLocal = LoadsGlobal 
            !For NURBS, we probably do not need to worry about the rearrangement part. 
        
            !we have applied the traction load on control points 1-4
            !What do we need from each of these 
            ! Node 1 - control point X coordinates = 0.025 m (corner CP)
            ! Node 2 - control point X coordinate = 0.0166 m (middle CP) 
            ! Node 3 - control point X coordinate = 0.00833 m (middle CP)
            ! Node 4 - control point X coordinate = 0.0 m (corner CP)
            
          !if (IsLoadApplied(LoadsLocal, NSurfaceNodes) ) then ! load values are not all zero
        
            !if ( ELEMENTTYPE == TETRAOLD ) then
            !  allocate( ElCon(N_NODES_HOE, Counters%NEl) )
            !  ElCon = ElementConnectivities10Node
            !else
              ! for NURBS we can apply this too
              allocate( ElCon(ELEMENTNODES, Counters%NEl) )
              ElCon = ElementConnectivities(:,:, IPatch_Temporary) ! we are just rewriting the connectivites
                
            !end if
         
            ! no need for change here for NURBS
            ! determine which elements lie next to the loaded surface and which side of the elements are adjacent to the loaded surface
            call DetermineElementsAdjacentSurface(IElTyp, Counters%NEl, NSurfaceNodes, IConLocal, ElCon, IsActiveElement(:,IPatch_Temporary), NElements, IElements, ISides)
            ! IElements(1) has the element number where the traction is applied. That is the only thing we are getting from above subroutine. 
            ! IElements(1) = 40
            
            allocate(NURBSTractionNodes(IConLocal(2))) ! IConLocal(2) is equal to 4
            NURBSTractionNodes(1) = 1
            !NURBSTractionNodes(2) = 2
            !NURBSTractionNodes(3) = 3
            NumberOfTractionElements = 1
          
            !IConSurface(1) -> start 
            !IConSurface(2) -> end 
          
          do ii = 2,IConLocal(2)
              NURBSTractionNodes(ii) = ii ! --> [1,2,3,4]    
          end do 
          
          !NURBSTRactionNodes = [1, 2, 3, 4]
          
           ! We know that the NURBS Traction Nodes go from 1 to 4. 
            
          ! NElements = 1 
            if ( NElements > 0 ) then ! found at least one ACTIVE element, maybe 2
        
              ! determine total load in x, y and z direction acting on the surface
              TotalSurfaceLoad = IntegrateVectorSurface(NSurfaceNodes, NGaussPoints, Counters%NodTot, NodeCoord, IConLocal, LoadsLocal, IElements, NumberOfTractionElements, NURBSTractionNodes)
        
              do IElement = 1, NElements ! 1 or 2
              
                ElementNumber = IElements(IElement) ! retrieve element number
                SideNumber = ISides(IElement) ! retrieve elment side number
                
                if ( NFORMULATION == 1 ) then
                  NParticles =  NPartEle(ElementNumber, IPatch_Temporary)
                else if ( NFORMULATION == 2 ) then
                  NParticles =  NSolidEle(ElementNumber, IPatch_Temporary) ! load is applied only on solid MPs
                else
                  call GiveError("Formulation has to be single-point or double-point in [subroutine TransferExternalLoadsToParticles()].")
                end if
        
                ! Determine the MPs next to the considered element side
                if ( allocated(ParticleStatus) )deallocate(ParticleStatus, stat = IError)
                allocate( ParticleStatus(NParticles), stat = IError )
                ParticleStatus = .false.
        
                call DetermineAdjacentParticles(SideNumber, NParticles, ParticleStatus)
        
                ! Determine the number of MPs next to the considered element side
                NLoadedParticles = DetermineNAdjacentParticles(ParticleStatus, NParticles)
        
                ParticleID = 0
                do IParticle = 1, NPartEle(ElementNumber, IPatch_Temporary) ! loop all MPs of the element (solid and liquid MPs in case of double-point formulation)
                  ParticleIndex = GetParticleIndex(IParticle, ElementNumber, IPatch_Temporary)
        
                  if ( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then
                      ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
        
                      ParticleID = ParticleID + 1
        
                      if(ParticleStatus(ParticleID)) then ! this MP is loaded
        
                        LoadPerParticle = DetermineLoadPerParticle(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal, Counters%NodTot, NodeCoord, Loadtype, IElements, NumberOfTractionElements, NURBSTractionNodes)
        
                        if ( LoadType == 1 ) then ! assign to solid phase / solid material point
                          call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        
                        else if ( LoadType == 2 ) then ! assign to water phase / liquid material point
                          call IncreaseFExtWater(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        
                        else if ( LoadType == 3 ) then ! assign to gas phase
                          call IncreaseFExtGas(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        
                          if ( NFORMULATION == 2 ) then
                            call GiveError('Erroneous definition of applied external load: For double-point formulation only load on solid or liquid phase is allowed.')
                          end if
        
                        end if ! load type
                      end if ! loaded y/n
                  end if ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
        
                end do ! loop over MPs
        
              end do ! loop over elements
            end if ! active element
          !end if ! load applied
        
        end subroutine TransferExternalLoadsToParticles_NURBS
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !------------------------------------------------------------------------------------------------------------------------
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        
        subroutine TransferExternalLoadsToParticles(LoadType, NSurfaceNodes, IElTyp, IConGlobal, LoadsGlobal, NGaussPoints, NodeCoord,ILoadSystem)
        !**********************************************************************
        !
        !    Function:  Transfer external distributed loads from element surfaces to
        !               particles closest to the considered surface - all particles get
        !               the same share of the distributed load.
        !               So far, only 6-noded triangular element sides are considered with
        !               constant distributed loads.
        !               If the loaded surface lies between two activated elements,
        !               particles of both elements get an equal share of the load.
        !
        !     LoadType : 1 solid 2 water 3 gas
        !     NSurfaceNodes : Number of nodes of the surface element
        !     IElTyp : Number of nodes per element
        !     IConGlobal : Connectivities of the surface element
        !     LoadsGlobal : Nodal load values of the surface element
        !     NGaussPoints : Number of Gauss points of the surface element
        !     NodeCoord : Nodal coordinates
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NSurfaceNodes, NGaussPoints, IElTyp, LoadType, ILoadSystem
          integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(NSurfaceNodes, NVECTOR), intent(in) :: LoadsGlobal
          real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: NodeCoord
        
          ! local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, IElement, NElements, NLoadedParticles, NParticles, IError, ParticleID
          integer(INTEGER_TYPE) :: ElementNumber, SideNumber
          integer(INTEGER_TYPE), dimension(2) :: IElements, ISides
          integer(INTEGER_TYPE), dimension(NSurfaceNodes) :: IConLocal
          integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElCon
          real(REAL_TYPE), dimension(NSurfaceNodes, NVECTOR) :: LoadsLocal
          logical, dimension(:), allocatable :: ParticleStatus
          real(REAL_TYPE), dimension(NVECTOR) ::  LoadPerParticle !TotalSurfaceLoad,
        
          real(REAL_TYPE), dimension(NVECTOR, 4) :: TotalSurfaceLoad
          
          integer(INTEGER_TYPE) :: NumberOfTractionElements
          
          integer(INTEGER_TYPE), dimension(:), allocatable :: NURBSTractionNodes
          
          integer(INTEGER_TYPE) :: ii
          
          ! if CornerControlPoints is 1 then it is a corner node 
          !integer(INTEGER_TYPE), dimension(4) :: CornerControlPoints ! this should be cubic order + 1 = 3 + 1 = 4
          integer(INTEGER_TYPE), dimension(NXiKnotOrder_Traction+1) :: CornerControlPoints ! this should be cubic order + 1 = 3 + 1 = 4

          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          
          !CornerControlPoints = 0
          !CornerControlPoints(1) = 1 ! 1 is a corner node 
          !CornerControlPoints(4) = 1 ! 4 is a corner node with nothing to the 
          !
          ! copy load connectivities to local arrays
          !select case(NSurfaceNodes)
          !  case(6) ! 6-noded triangular element, 3D case
          !    call RearrangeConnectivitiesTRI6(IConGlobal, LoadsGlobal, IConLocal, LoadsLocal)
          !  case(2) ! 2-noded line element, 2D case
          !    call RearrangeConnectivitiesLINE2(IConGlobal, LoadsGlobal, IConLocal, LoadsLocal)
          !  case(3) ! 3-noded line element, 2D case
          !    call RearrangeConnectivitiesLINE3(IConGlobal, LoadsGlobal, IConLocal, LoadsLocal)
          !  end select
            !
            !IConLocal = IConGlobal 
            !LoadsLocal = LoadsGlobal 
            !For NURBS, we probably do not need to worry about the rearrangement part. 
        
            !we have applied the traction load on control points 1-4
            !What do we need from each of these 
            ! Node 1 - control point X coordinates = 0.025 m (corner CP)
            ! Node 2 - control point X coordinate = 0.0166 m (middle CP) 
            ! Node 3 - control point X coordinate = 0.00833 m (middle CP)
            ! Node 4 - control point X coordinate = 0.0 m (corner CP)
            
          ! I commented this because I am not needing to check if they are zero or not
          !if (IsLoadApplied(LoadsLocal, NSurfaceNodes) ) then ! load values are not all zero
        
            !if ( ELEMENTTYPE == TETRAOLD ) then
            !  allocate( ElCon(N_NODES_HOE, Counters%NEl) )
            !  ElCon = ElementConnectivities10Node
            !else
            !  ! for NURBS we can apply this too
            !  allocate( ElCon(ELEMENTNODES, Counters%NEl) )
            !  ElCon = ElementConnectivities(:,:, IPatch_Temporary) ! we are just rewriting the connectivites
            !    
            !end if
         
            ! no need for change here for NURBS
            ! determine which elements lie next to the loaded surface and which side of the elements are adjacent to the loaded surface
            call DetermineElementsAdjacentSurface(IElTyp, Counters%NEl, NSurfaceNodes, IConLocal, ElCon, IsActiveElement(:,IPatch_Temporary), NElements, IElements, ISides)
            
            ! this needs fixing here 
            ! IElements(1) has the element number where the traction is applied. That is the only thing we are getting from above subroutine. 
            ! IElements(1) = 40
            
            allocate(NURBSTractionNodes(IConLocal(2))) ! IConLocal(2) is equal to 4
            NURBSTractionNodes(1) = 1
            !NURBSTractionNodes(2) = 2
            !NURBSTractionNodes(3) = 3
            NumberOfTractionElements = 1
          
            !IConSurface(1) -> start 
            !IConSurface(2) -> end 
          
          do ii = 2,IConLocal(2)
              NURBSTractionNodes(ii) = ii ! --> [1,2,3,4]    
          end do 
          
          !NURBSTRactionNodes = [1, 2, 3, 4]
          
           ! We know that the NURBS Traction Nodes go from 1 to 4. 
            
          ! NElements = 1 
            if ( NElements > 0 ) then ! found at least one ACTIVE element, maybe 2
        
              ! determine total load in x, y and z direction acting on the surface
              TotalSurfaceLoad = IntegrateVectorSurface(NSurfaceNodes, NGaussPoints, Counters%NodTot, NodeCoord, IConLocal, LoadsLocal, IElements, NumberOfTractionElements, NURBSTractionNodes)
        
              do IElement = 1, NElements ! 1 or 2
              
                ElementNumber = IElements(IElement) ! retrieve element number
                SideNumber = ISides(IElement) ! retrieve elment side number
                
                if ( NFORMULATION == 1 ) then
                  NParticles =  NPartEle(ElementNumber, IPatch_Temporary)
                else if ( NFORMULATION == 2 ) then
                  NParticles =  NSolidEle(ElementNumber, IPatch_Temporary) ! load is applied only on solid MPs
                else
                  call GiveError("Formulation has to be single-point or double-point in [subroutine TransferExternalLoadsToParticles()].")
                end if
        
                ! Determine the MPs next to the considered element side
                if ( allocated(ParticleStatus) )deallocate(ParticleStatus, stat = IError)
                allocate( ParticleStatus(NParticles), stat = IError )
                ParticleStatus = .false.
        
                call DetermineAdjacentParticles(SideNumber, NParticles, ParticleStatus)
        
                ! Determine the number of MPs next to the considered element side
                NLoadedParticles = DetermineNAdjacentParticles(ParticleStatus, NParticles)
        
                ParticleID = 0
                do IParticle = 1, NPartEle(ElementNumber, IPatch_Temporary) ! loop all MPs of the element (solid and liquid MPs in case of double-point formulation)
                  ParticleIndex = GetParticleIndex(IParticle, ElementNumber, IPatch_Temporary)
        
                  if ( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then
                      ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
        
                      ParticleID = ParticleID + 1
        
                      if(ParticleStatus(ParticleID)) then ! this MP is loaded
        
                        LoadPerParticle = DetermineLoadPerParticle(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal, Counters%NodTot, NodeCoord, Loadtype, IElements, NumberOfTractionElements, NURBSTractionNodes)
        
                        if ( LoadType == 1 ) then ! assign to solid phase / solid material point
                          call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        
                        else if ( LoadType == 2 ) then ! assign to water phase / liquid material point
                          call IncreaseFExtWater(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        
                        else if ( LoadType == 3 ) then ! assign to gas phase
                          call IncreaseFExtGas(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        
                          if ( NFORMULATION == 2 ) then
                            call GiveError('Erroneous definition of applied external load: For double-point formulation only load on solid or liquid phase is allowed.')
                          end if
        
                        end if ! load type
                      end if ! loaded y/n
                  end if ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
        
                end do ! loop over MPs
        
              end do ! loop over elements
            end if ! active element
          !end if ! load applied
        
        end subroutine TransferExternalLoadsToParticles


        logical function IsLoadApplied(LoadsLocal, NSurfaceNodes)
        !**********************************************************************
        !
        !    Function:  Check whether any of the load values are larger than zero.
        !
        !     LoadsLocal : Values of the specified distributed load
        !     NSurfaceNodes : Number of nodes of the element surface
        !
        ! O   IsLoadApplied : .true. if a load value is non-zero
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          real(REAL_TYPE), dimension(NSurfaceNodes, NVECTOR), intent(in) :: LoadsLocal
          integer(INTEGER_TYPE), intent(in) :: NSurfaceNodes
          ! Local variables
          integer(INTEGER_TYPE) :: I, J

          IsLoadApplied = .false.
          
          do I = 1, NSurfaceNodes
            do J = 1, NVECTOR
              if (abs(LoadsLocal(I, J) )>1.0E-4) then
                IsLoadApplied = .true.
                EXIT
              end if
            end do
            if (IsLoadApplied) then
              EXIT
            end if
          end do
        
        end function IsLoadApplied


        function DetermineLoadPerParticle(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal,NodTot,NodeCoord,LoadType, IElements, NumberOfTractionElements, NURBSTractionElements)
        !**********************************************************************
        !
        !    Function:  Returns the load to be applied on each particle.
        !               The total load acting on the surface is divided by
        !               the number of particles adjacent to the surface and,
        !               eventually, the number of adjacent elements.
        !
        !     I : Local particle number adjacent to the loaded surface (from 1 to NLoadedParticles)
        !     NLoadedParticles : Number of particles adjacent to the loaded surface
        !     IConLocal : Rearranged surface node connectivities (corner first and then middle point)
        !     LoadsLocal : Values of the specified distributed load
        !     NodTot : Total number of nodes
        !     NodeCoord : Nodal Coordinates
        !
        ! O   DetermineLoadPerParticle : Load to be assigned to each particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none
         
          integer(INTEGER_TYPE) :: NLoadedParticles, NSurfaceNodes, nNod3el, LoadType
          real(REAL_TYPE), dimension(NSurfaceNodes, NVECTOR) :: LoadsLocal
          integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(in) :: IConLocal
          integer(INTEGER_TYPE), intent(in) :: NodTot
          real(REAL_TYPE), dimension(NodTot, NVECTOR), intent(in) :: NodeCoord
          real(REAL_TYPE), dimension(NVECTOR) :: DetermineLoadPerParticle, PartDistrLoading
          ! Local variables
          integer(INTEGER_TYPE) :: j, k, nn
          real(REAL_TYPE) ::  Factor, Radius !DetJ,
          real(REAL_TYPE), dimension(NVECTOR) :: VectorN      
          integer(INTEGER_TYPE), dimension(:), intent(in) :: IElements
          integer(INTEGER_TYPE), intent(in) :: NumberOfTractionElements
          integer(INTEGER_TYPE), dimension(:), intent(in) :: NURBSTractionElements
          real(REAL_TYPE), dimension(4) :: DetJ
          ! hardcoded for cubic order. We are looking at each control point seperately. 
          ! we are also assuming that the effective length 
          !NSurfaceNodes = 3
          integer(INTEGER_TYPE) :: NGP, IGP, IDim
          integer(INTEGER_TYPE) :: INode, nNode
          NGP = 2
          nNode = 4
		  select case(NSurfaceNodes)
    !      	case(6) ! 6-noded triangular element
				!nNod3el = NSurfaceNodes/2
    !      		call Normal_T3(1, NodeCoord, IConLocal, nNod3el, GPShapeFunctionDerivativeBoundary, VectorN, DetJ)
				!Factor = 0.5d0 !Factor for triangles (the area is half of the rectangle one)
			case(2) ! 2- noded linear element
				nNod3el = NSurfaceNodes
                !call NormalOnLine(1, NodeCoord, IConLocal, GPShapeFunctionDerivativeBoundary_XI_AllElements, VectorN, DetJ, NumberOfTractionElements, NURBSTractionElements, IElements)
				!call NormalOnLine(IGP, NodeCoord, IConSurface, GPShapeFunctionDerivativeBoundary_XI_AllElements, VectorN, DetJ, NumberOfTractionElements, NURBSTractionNodes, IElements)!, IConLElement)
                Factor = 2.0d0 !the lenght of linear element is 2*DetJ
            !case(3) 
            !    nNod3el = NSurfaceNodes
            end select
          ! determine vector Vn normal to a plane (for 3-noded triangles)
          ! detJ = det(Jacobian_S)
          
          ! Compute the distributed loading at the center of the surface
          PartDistrLoading = 0.0
            !Do j = 1, NVECTOR !direction in space
            !  Do k = 1, nNod3el ! over the surface element
            !    if ( ISAXISYMMETRIC ) then
            !        nn = IConLocal(k)
            !        Radius = NodeCoord(nn, 1) ! index 1 is r-direction
            !    else
            !        Radius = 1.0
            !    end if
            !    if (LoadType == 2 .or. LoadType == 3) then ! Liquid or Gas, Loading is applied only normal to the surface
            !      PartDistrLoading(j) = PartDistrLoading(j) + GPShapeFunctionBoundary(1,k) * LoadsLocal(k,j) * VectorN(j) * Radius
            !    else
            !      PartDistrLoading(j) = PartDistrLoading(j) + GPShapeFunctionBoundary(1,k) * LoadsLocal(k,j) * Radius
            !    end if
            !  end do
            !end do
            
            do IGP = 1, NGP
                do INode = 1, nNode
                    do IDim = 1, NVECTOR
                        PartDistrLoading(IDim) = PartDistrLoading(IDim) + GPShapeFunctionBoundary_XI(IGP,INode) * LoadsLocal(1,IDim) * DetJ(INode)
                    end do 
                end do
            end do 
            
           
          ! Compute the particle loading
          if (LoadType == LOADTYPE_LIQUID .or. LoadType == LOADTYPE_GAS) then ! Liquid or Gas, Loading is applied only normal to the surface
            Do J = 1, NVECTOR !direction in space
              DetermineLoadPerParticle(J) = PartDistrLoading(J) * DetJ(J) * factor/NLoadedParticles ! Hardcoded: J index on DetJ
            End do  
          else    
              !!!! NEED TO FIX HERE DETJ
            DetermineLoadPerParticle = PartDistrLoading * factor/NLoadedParticles !* 0.5 !HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This is hardcoded here... DetJ needs to be 0.5. We need to include the NURBS space jacobian before when we evaluate DetJ 
          end if
        end function DetermineLoadPerParticle
        
        
        
        
        
        
        
   !     function DetermineLoadPerParticle_NURBS(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal,NodTot,NodeCoord,LoadType, IElements, NumberOfTractionElements, NURBSTractionElements)
   !     !**********************************************************************
   !     !
   !     !    Function:  Returns the load to be applied on each particle.
   !     !               The total load acting on the surface is divided by
   !     !               the number of particles adjacent to the surface and,
   !     !               eventually, the number of adjacent elements.
   !     !
   !     !     I : Local particle number adjacent to the loaded surface (from 1 to NLoadedParticles)
   !     !     NLoadedParticles : Number of particles adjacent to the loaded surface
   !     !     IConLocal : Rearranged surface node connectivities (corner first and then middle point)
   !     !     LoadsLocal : Values of the specified distributed load
   !     !     NodTot : Total number of nodes
   !     !     NodeCoord : Nodal Coordinates
   !     !
   !     ! O   DetermineLoadPerParticle : Load to be assigned to each particle
   !     !
   !     ! Implemented in the frame of the MPM project.
   !     !
   !     !**********************************************************************
   !     implicit none
   !      
   !       integer(INTEGER_TYPE) :: NLoadedParticles, NSurfaceNodes, nNod3el, LoadType
   !       real(REAL_TYPE), dimension(NSurfaceNodes, NVECTOR) :: LoadsLocal
   !       integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(in) :: IConLocal
   !       integer(INTEGER_TYPE), intent(in) :: NodTot
   !       real(REAL_TYPE), dimension(NodTot, NVECTOR), intent(in) :: NodeCoord
   !       real(REAL_TYPE), dimension(NVECTOR) :: DetermineLoadPerParticle, PartDistrLoading
   !       ! Local variables
   !       integer(INTEGER_TYPE) :: j, k, nn
   !       real(REAL_TYPE) ::  Factor, Radius !DetJ,
   !       real(REAL_TYPE), dimension(NVECTOR) :: VectorN      
   !       integer(INTEGER_TYPE), dimension(:), intent(in) :: IElements
   !       integer(INTEGER_TYPE), intent(in) :: NumberOfTractionElements
   !       integer(INTEGER_TYPE), dimension(:), intent(in) :: NURBSTractionElements
   !       real(REAL_TYPE), dimension(4) :: DetJ
   !       ! hardcoded for cubic order. We are looking at each control point seperately. 
   !       ! we are also assuming that the effective length 
   !       !NSurfaceNodes = 3
   !       integer(INTEGER_TYPE) :: NGP, IGP, IDim
   !       integer(INTEGER_TYPE) :: INode, nNode
   !       NGP = 2
   !       nNode = 4
		 ! select case(NSurfaceNodes)
   ! !      	case(6) ! 6-noded triangular element
			!	!nNod3el = NSurfaceNodes/2
   ! !      		call Normal_T3(1, NodeCoord, IConLocal, nNod3el, GPShapeFunctionDerivativeBoundary, VectorN, DetJ)
			!	!Factor = 0.5d0 !Factor for triangles (the area is half of the rectangle one)
			!case(2) ! 2- noded linear element
			!	nNod3el = NSurfaceNodes
   !             !call NormalOnLine(1, NodeCoord, IConLocal, GPShapeFunctionDerivativeBoundary_XI_AllElements, VectorN, DetJ, NumberOfTractionElements, NURBSTractionElements, IElements)
			!	!call NormalOnLine(IGP, NodeCoord, IConSurface, GPShapeFunctionDerivativeBoundary_XI_AllElements, VectorN, DetJ, NumberOfTractionElements, NURBSTractionNodes, IElements)!, IConLElement)
   !             Factor = 2.0d0 !the lenght of linear element is 2*DetJ
   !         !case(3) 
   !         !    nNod3el = NSurfaceNodes
   !         end select
   !       ! determine vector Vn normal to a plane (for 3-noded triangles)
   !       ! detJ = det(Jacobian_S)
   !       
   !       ! Compute the distributed loading at the center of the surface
   !       PartDistrLoading = 0.0
   !         !Do j = 1, NVECTOR !direction in space
   !         !  Do k = 1, nNod3el ! over the surface element
   !         !    if ( ISAXISYMMETRIC ) then
   !         !        nn = IConLocal(k)
   !         !        Radius = NodeCoord(nn, 1) ! index 1 is r-direction
   !         !    else
   !         !        Radius = 1.0
   !         !    end if
   !         !    if (LoadType == 2 .or. LoadType == 3) then ! Liquid or Gas, Loading is applied only normal to the surface
   !         !      PartDistrLoading(j) = PartDistrLoading(j) + GPShapeFunctionBoundary(1,k) * LoadsLocal(k,j) * VectorN(j) * Radius
   !         !    else
   !         !      PartDistrLoading(j) = PartDistrLoading(j) + GPShapeFunctionBoundary(1,k) * LoadsLocal(k,j) * Radius
   !         !    end if
   !         !  end do
   !         !end do
   !         
   !         do IGP = 1, 4 ! 4 gauss points
   !             do INode = 1, nNode
   !                 do IDim = 1, NVECTOR
   !                     PartDistrLoading(IDim) = PartDistrLoading(IDim) + GPShapeFunctionBoundary_XI(IGP,INode) * LoadsLocal(1,IDim) * DetJ(INode)
   !                 end do 
   !             end do
   !         end do 
   !         
   !        
   !       ! Compute the particle loading
   !       !if (LoadType == LOADTYPE_LIQUID .or. LoadType == LOADTYPE_GAS) then ! Liquid or Gas, Loading is applied only normal to the surface
   !       !  Do J = 1, NVECTOR !direction in space
   !       !    DetermineLoadPerParticle(J) = PartDistrLoading(J) * DetJ(J) * factor/NLoadedParticles ! Hardcoded: J index on DetJ
   !       !  End do  
   !       !else    
   !           !!!! NEED TO FIX HERE DETJ
   !         !DetermineLoadPerParticle = PartDistrLoading * factor/NLoadedParticles !* 0.5 !HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !         ! This is hardcoded here... DetJ needs to be 0.5. We need to include the NURBS space jacobian before when we evaluate DetJ 
   !       !end if
   !     end function DetermineLoadPerParticle_NURBS
        
        
        
        
        
        
        
        
        
        
        
        
        

        
        subroutine AllocateMaterialModelsStateVariables(ParticleIndex)
        !**********************************************************************
        !
        !  Function:  Allocation of the state parameters for ICH and HP models
        !
        !     ParticleIndex : the index of the considered particle
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE) :: ParticleIndex
     
          ! Local variables
          integer(INTEGER_TYPE) :: I, MaterialIndex
          character(len=64) :: SoilModel ! name of the constitutive model
          
          MaterialIndex = MaterialIDArray(ParticleIndex) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
          SoilModel = MatParams(MaterialIndex)%MaterialModel ! name of constitutive model as specified in GOM-file
          

          if (SoilModel == ESM_HYPOPLASTICITY_SAND) then ! HP model

            ! initialise the state parameter PhiCS (it is initialised inside the HP routine)
            call SetHPStateVariablesI(Particles(ParticleIndex), 1, 0.d0)

            ! initialise the state parameter Void ratio (it is initialised inside the HP routine)
            call SetHPStateVariablesI(Particles(ParticleIndex), 2, 0.d0)

            if (MatParams(MaterialIndex)%IGSmr/=(0.d0).and.MatParams(MaterialIndex)%IGSmt/=(0.d0)) then
              ! if (mr and mt) > 0.0 then intergranular strain
              do I = 1, 7  ! no. of state variables
                ! initialise the state parameters (they are initialised inside the HP routine)
                call SetHPIGStateVariablesI(Particles(ParticleIndex), I, 0.d0)
              end do
            end if

          end if



          if (SoilModel==ESM_MOHR_COULOMB_STRAIN_SOFTENING) then ! Strain Softening MC model
            do I = 1, NTENSOR
              call SetEpsPI(Particles(ParticleIndex),I,0.d0)
            end do
          end if

        end subroutine AllocateMaterialModelsStateVariables

        
        subroutine RewriteParticleMaterialParameters()
        !**********************************************************************
        !
        !    Function:  update the material of particles using the ratio of the new properties
        !
        !               G, mass, materialweight, mixedweight, cohesion, phi, 
        !               porosity, degree of saturation, conductivities, bulkmodulus..
        !
        !**********************************************************************
          implicit none
          
          integer(INTEGER_TYPE) I, J, K
          real(REAL_TYPE) sinphi, sinpsi, cohesion
          real(REAL_TYPE) oldMaterialWeight, newMaterialWeight, rateMaterialWeight

          real(REAL_TYPE) oldMixedWeight, newMixedWeight, rateMixedWeight
          
          do I = 1, Counters%NParticles
            J = MaterialIDArray(I)
            
            if (IsMaterialParticle(I)) then
              oldMaterialWeight = Particles(I)%MaterialWeight
              newMaterialWeight = MatParams(J)%DryWeight
              rateMaterialWeight = newMaterialWeight/oldMaterialWeight

              oldMixedWeight = Particles(I)%MixedWeight
              newMixedWeight = MatParams(J)%WeightMixture
              rateMixedWeight = newMixedWeight/oldMixedWeight
            else
              rateMaterialWeight = 1.0
              rateMixedWeight = 1.0
            end if
            
            Particles(I)%ShearModulus = MatParams(J)%ShearModulus
            
            sinphi = SIN(MatParams(J)%FrictionAngle*(Pi/180.0))
            cohesion = MatParams(J)%Cohesion
            sinpsi = SIN(MatParams(J)%DilatancyAngle*(Pi/180.0))

            Particles(I)%CohesionCosPhi = cohesion * dsqrt(1D0 - SinPhi * SinPhi)
     
            if (sinphi>1D-10) then
              Particles(I)%SFail = Particles(I)%CohesionCosPhi / sinphi
            else
              Particles(I)%SFail = 1D10
            end if  

            Particles(I)%Density = MatParams(J)%DryWeight / CalParams%GravityData%GAccel
            Particles(I)%MaterialWeight = Particles(I)%MaterialWeight * rateMaterialWeight
            
          ! Determine the mass of the particle
            MassArray(I) = MassArray(I) * rateMaterialWeight

          ! Determine dry weight
            do K = 1, NVECTOR
              Particles(I)%FBody(K) = Particles(I)%FBody(K) * rateMaterialWeight
            end do  
                                                               
            Particles(I)%Conductivity = MatParams(J)%HydraulicConductivityLiquid
            Particles(I)%ConductivityGas = MatParams(J)%HydraulicConductivityGas
            !Particles(I)%Porosity = MatParams(J)%InitialPorosity
            Particles(I)%InitialPorosity = MatParams(J)%InitialPorosity
           ! Particles(I)%DegreeSaturation = MatParams(J)%InitialDegreeOfSaturation
            Particles(I)%BulkWater = MatParams(J)%BulkModulusLiquid
            Particles(I)%BulkGas = MatParams(J)%BulkModulusGas
              
             if (CalParams%ApplyPartialSaturation) then
           !  Hydraulic Conductivity at MP depends on Sr, Material Hydraulic Conductivity has been updated but it is not mapped at particles yet
               call RewriteParticleHydraulicConductivity(I,J)
           end if
              
           ! Determine mixed (saturated) weight
            do K = 1, NVECTOR
              Particles(I)%FBodyMixed(K) = Particles(I)%FBodyMixed(K) * rateMixedWeight
            end do  
     
          end do

        end subroutine RewriteParticleMaterialParameters
        
        
        subroutine ReadNumberOfMaterialPointsPerElement(NumberOfSolidMaterialPoints, NumberOfLiquidMaterialPoints)
        !**************************************************************************************************************************
        !
        !  Function : Read number of solid and liquid material points per element as specified in GOM file
        !
        !  I/O  NumberOfSolidMaterialPoints : number of solid material points per element read from GOM-file
        !  I/O  NumberOfLiquidMaterialPoints : number of liquid material points per element read from GOM-file
        !
        !**************************************************************************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), dimension(:, :), intent(inout) :: NumberOfSolidMaterialPoints, NumberOfLiquidMaterialPoints
          
          ! Local variables
          character(len=1023) :: FileName
          character(len=255) :: BName
          integer(INTEGER_TYPE) :: ElementID, ios, DumI1, DumI2
          character(len=21) :: messageIOS = 'GOM file: Can''t read '
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
        
          FileName=Trim(CalParams%FileNames%ProjectName)//'.GOM'
         
          if (FExist(FileName)) open(GOMunit, FILE = FileName)
          
          do
            read(GOMunit, '(A)') BName
            if (trim(BName) == '$$PATCH_ID') then 
                ! assign IPatch = -1 incase it does not read appropriately we know it went wrong here
                IPatch = -1
                ! read patch ID from the GOM file
                read(GOMunit, *) IPatch ! --> now we knot patch index so we can read in the following: 
                                        ! 1) Number of control points for each patch
                                        ! 1) Number of knot entries (xi, eta, zeta) 
                                        ! 2) Order of knot (xi, eta, zeta)
                                        ! 3) 
            
            
            
            
            
            else if (trim(BName)=='$$START_NUMBER_OF_MATERIAL_POINTS') then
              !do IPatch = 1, NumberOfPatches ! loop over patches
                do ElementID = 1, nel_NURBS(IPatch)!Counters%NEl ! loop over elements
                read (GOMunit,*,iostat=ios) DumI1, DumI2
                call Assert( ios == 0, messageIOS // trim(BName) )   
                select case(ELEMENTTYPE)
                  case(TETRAOLD)  
                    call Assert( DumI1==0 .or. DumI1==1 .or. DumI1==4 .or. DumI1==7 .or. DumI1==8 .or. DumI1==10 .or. DumI1==13 .or. DumI1==20, 'GOM-ERROR: $$START_NUMBER_OF_MATERIAL_POINTS (solid) must be 0,1,4,7,8,10,13,20 for tetrahedral elements!' )
                    call Assert( DumI2==0 .or. DumI2==1 .or. DumI2==4 .or. DumI2==7 .or. DumI2==8 .or. DumI2==10 .or. DumI2==13 .or. DumI2==20, 'GOM-ERROR: $$START_NUMBER_OF_MATERIAL_POINTS (liquid) must be 0,1,4,7,8,10,13,20 for tetrahedral elements!' )
                  case(TRI3)
                    call Assert( DumI1==0 .or. DumI1==1 .or. DumI1==3 .or. DumI1==6 .or. DumI1==12 .or. DumI1==16 .or. DumI1==25 .or. DumI1==46 .or. DumI1==88, 'GOM-ERROR: $$START_NUMBER_OF_MATERIAL_POINTS (solid) must be 0,1,3,6,12,16,25,46,88 for triangular elements!' )
                    call Assert( DumI1==0 .or. DumI1==1 .or. DumI1==3 .or. DumI1==6 .or. DumI1==12 .or. DumI1==16 .or. DumI1==25 .or. DumI1==46 .or. DumI1==88, 'GOM-ERROR: $$START_NUMBER_OF_MATERIAL_POINTS (liquid) must be 0,1,3,6,12,16,25,46,88 for triangular elements!' )
                end select      
                NumberOfSolidMaterialPoints(ElementID, IPatch) = DumI1 ! number of solid material points for the element
                NumberOfLiquidMaterialPoints(ElementID, IPatch) = DumI2 ! number of liquid material points for the element
                !ELEMENTBOUNDARYGAUSSPOINTS = NumberOfSolidMaterialPoints(ElementID, IPatch) ** 0.5
                !ELEMENTBOUNDARYGAUSSPOINTS_XI = NumberOfSolidMaterialPoints(ElementID, IPatch) ** 0.5
                !ELEMENTBOUNDARYGAUSSPOINTS_ETA = NumberOfSolidMaterialPoints(ElementID, IPatch) ** 0.5
                end do ! loop over element
                !end do ! loop over patches
            else if (trim(BName)=='$$FINISH') then
              EXIT
            end if 
          end do      
          
          close(GOMunit)

        end subroutine ReadNumberOfMaterialPointsPerElement

        
        subroutine CheckNumberOfMaterialPointsPerElement(NumberOfSolidMaterialPoints, NumberOfLiquidMaterialPoints)
        !--------------------------------------------------------------------
        !
        !  Function:  Automatic check of number of solid and liquid material points per element with material type.
        !             Calculation is aborted in case of mismatch.
        !
        !   First, distinction between 1- and 2-layer formulation is made
        !   Then, each element is determined to be (in)active.
        !   When inactive, the number of solid and liquid material points must be equal to zero.
        !   When active, the number of solid and liquid material points is checked according to table:
        !
        !                                  solid mp's     liquid mp's
        !   1-layer all material types     >0             =0
        !   2-layer 1-phase-solid          >0             =0
        !   2-layer 1-phase-liquid         =0             >0
        !   2-layer 2-phase-solid          >0             >0
        !
        !--------------------------------------------------------------------
            integer(INTEGER_TYPE), dimension(:, :), intent(inout) :: NumberOfSolidMaterialPoints, NumberOfLiquidMaterialPoints
            ! local variables 
            integer(INTEGER_TYPE) :: ElementID, I, NSolid, NLiquid
            
            ! Multipatch variables 
            integer(INTEGER_TYPE) :: IPatch_Temporary = 1
            integer(INTEGER_TYPE) :: IPatch
            
            
            do IPatch = 1, NumberOfPatches ! loop over patches
                
                ! determine if single- or double-point formulation
            NSolid = 0
            NLiquid = 0
                
            do ElementID = 1, nel_NURBS(IPatch)!Counters%NEl ! loop over elements
              NSolid = NSolid + NumberOfSolidMaterialPoints(ElementID, IPatch)
              NLiquid = NLiquid + NumberOfLiquidMaterialPoints(ElementID, IPatch)
            end do  ! loop over elements
            !end do ! loop over patches
            
            ! NSolid = 0 AND NLiquid = 0 --> error
            ! NSolid = 0 AND NLiquid > 0 --> single-point formulation --> liquid only
            ! NSolid > 0 AND NLiquid = 0 --> single-point formulation --> can be 1-/2-/3-phase
            ! NSolid > 0 AND NLiquid > 0 --> double-point formulation
            
            call SetFormulation(SINGLE_POINT)
            !if ( NSolid == 0 .AND. NLiquid == 0 ) call GiveError('No material points are specified.') 
            !--> I commented this because we have no MPs in a patch
            if ( NSolid > 0 .AND. NLiquid > 0 ) call SetFormulation(DOUBLE_POINT) ! double-point formulation
            
            if (NFORMULATION == 1) then ! single-point formulation   
              

                !do IPatch = 1, NumberOfPatches ! loop over patches
                do ElementID = 1, nel_NURBS(IPatch)!Counters%NEl ! loop over elements
        
                if (ElementMaterialID(ElementID,IPatch) <= 0) then ! inactive element
                  
                  if((NumberOfSolidMaterialPoints(ElementID, IPatch) > 0) .or. (NumberOfLiquidMaterialPoints(ElementID, IPatch) > 0)) then
                    call GiveError('Material point specification of element '// trim(String(ElementID)) // ' does not coincide with material type.')
                  end if
         
                else ! active element
                  
                  if((NumberOfSolidMaterialPoints(ElementID, IPatch) <= 0) .or. (NumberOfLiquidMaterialPoints(ElementID, IPatch) > 0)) then
                    call GiveError('Material point specification of element '// trim(String(ElementID)) // ' does not coincide with material type.')
                  end if
              
                end if
                end do ! end loop over elements
                !end do ! end loop over patches

            else ! double-point formulation  
            
                !do IPatch = 1, NumberOfPatches ! loop over patches
              do ElementID = 1, nel_NURBS(IPatch)!Counters%NEl ! loop over elements
                if (ElementMaterialID(ElementID,IPatch) <= 0) then ! inactive element
                  
                  if((NumberOfSolidMaterialPoints(ElementID, IPatch) > 0) .or. (NumberOfLiquidMaterialPoints(ElementID, IPatch) > 0)) then
                    call GiveError('Material point specification of element '// trim(String(ElementID)) // ' does not coincide with material type.')
                  end if
         
                else ! active element
              
                  do I = 1, CalParams%NumberOfMaterials ! loop over number of materials
                    if (ElementMaterialID(ElementID,IPatch) == MatParams(I)%MaterialIndex) then
                      if ((MatParams(I)%MaterialType=='1-phase-solid') .or. (MatParams(I)%MaterialType==DRY_SOIL)) then
                        if((NumberOfSolidMaterialPoints(ElementID, IPatch) <= 0).or.(NumberOfLiquidMaterialPoints(ElementID, IPatch) > 0)) then
                          call GiveError('Material point specification of element '// trim(String(ElementID)) // ' does not coincide with material type.')
                        else
                          CYCLE
                        end if
                      else if ((MatParams(I)%MaterialType=='1-phase-liquid') .or. (MatParams(I)%MaterialType==LIQUID)) then
                        if((NumberOfSolidMaterialPoints(ElementID,IPatch) > 0).or.(NumberOfLiquidMaterialPoints(ElementID,IPatch) <= 0)) then
                          call GiveError('Material point specification of element '// trim(String(ElementID)) // ' does not coincide with material type.')
                        else
                          CYCLE
                        end if
                      else if ((MatParams(I)%MaterialType=='2-phase') .or. (MatParams(I)%MaterialType==SATURATED_SOIL_COUPLED)) then
                        if((NumberOfSolidMaterialPoints(ElementID,IPatch) <= 0)) then
                          call GiveError('Material point specification of element '// trim(String(ElementID)) // ' does not coincide with material type.')
                        else
                          CYCLE
                        end if
                      else
                        call GiveError('Wrong material type of element '// trim(String(ElementID))) 
                      end if
                    else
                      CYCLE
                    end if
                  end do ! end loop over number of materials
                
                end if
              
              end do ! end loop over elements  
              !end do ! end loop over patches
            
            end if
            end do ! end loop over patches
            
        end subroutine CheckNumberOfMaterialPointsPerElement
        
                
        integer(INTEGER_TYPE) function getMaterialEntity(MaterialSetID) result(IEntity)
        implicit none
        integer(INTEGER_TYPE), intent(in):: MaterialSetID
        
        ! set IEntityPile = HARD_ENTITY otherwise SOFT_ENTITY
        if (MaterialSetID == CalParams%MovingMesh%StructureMaterialID) then
          IEntity = HARD_ENTITY
        else
          IEntity = SOFT_ENTITY
        end if
      
      end function getMaterialEntity

      
      subroutine DetermineParticlesInLoadedElement(IConGlobal, LoadsGlobal, LoadedParticle)
        !**********************************************************************
        !
        !    Function: determine particles in loaded elements
        !
        !**********************************************************************
        implicit none

          ! arguments
          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal ! size (NSurfaceNodes)
          real(REAL_TYPE), dimension(:, :), intent(in) :: LoadsGlobal ! size (NSurfaceNodes, NVECTOR)
          logical, dimension(:), optional, intent(inout) :: LoadedParticle ! size (Counters%NParticles)
          
          ! Local variables
          integer(INTEGER_TYPE) :: NSurfaceNodes, IElTyp
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, IElement, NElements
          integer(INTEGER_TYPE), dimension(2) :: IElements, ISides ! 2 for 2D and 3D
          integer(INTEGER_TYPE), dimension(:), allocatable :: IConLocal
          real(REAL_TYPE), dimension(:, :), allocatable :: LoadsLocal
          integer(INTEGER_TYPE), dimension(:, :), allocatable :: ElCon
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          
          NSurfaceNodes = size( IConGlobal )  ! =ELEMENTBOUNDARYNODES
          allocate( IConLocal(NSurfaceNodes), LoadsLocal(NSurfaceNodes, NVECTOR) )

          if ( ELEMENTTYPE == TETRAOLD ) then 
            allocate( ElCon( N_NODES_HOE, Counters%NEl) )
            ElCon = ElementConnectivities10Node
            IElTyp = N_NODES_HOE
          else
            allocate( ElCon( ELEMENTNODES, Counters%NEl) )
            ElCon = ElementConnectivities(:,:, IPatch_Temporary)
            IElTyp = ELEMENTNODES
          end if

          ! copy load connectivities to local arrays depending on element type
          call RearrangeConnectivitiesPointer(IConGlobal, LoadsGlobal, IConLocal, LoadsLocal)
          
          if ( IsLoadApplied(LoadsLocal, NSurfaceNodes) ) then ! load values are not all zero

            ! determine which elements are located next to the loaded surface and which elementsides are adjacent to the loaded surface
            call DetermineElementsAdjacentSurface(IElTyp, Counters%NEl, NSurfaceNodes, IConLocal, ElCon, IsActiveElement(:,IPatch_Temporary), NElements, IElements, ISides)

            if ( NElements > 0 ) then ! found at least one ACTIVE element, maybe 2

              do IElement = 1, NElements
                do IParticle = 1, NPartEle(IElements(IElement), IPatch_Temporary)
                  
				  ParticleIndex = GetParticleIndex(IParticle, IElements(IElement), IPatch_Temporary)                 
                  LoadedParticle(ParticleIndex) = .true.  
                 
                end do ! loop over particles
              end do
            end if ! active element
          end if ! load applied

      end subroutine DetermineParticlesInLoadedElement  
       
      
        subroutine InitialLocalMaterialPointCoordinates(IParticle, SolidPointsElement, LiquidPointsElement, WeiGP, PosGP)
        !**********************************************************************
        !
        !  Function : Determines the local coordinates and integration weight assigned to material point with ID
        !             IParticle which are returned through WeiGP and PosGP.
        !             Currently, all material points are placed at the same local positions.
        !
        !  I  IParticle : Number of the material point
        !  I  SolidPointsElement : Number of solid material points per element
        !  I  LiquidPointsElement : Number of liquid material points per element
        !  O  WeiGP : Initial weight assigned to material point IParticle
        !  O  PosGP : Initial local position of material point IParticle
        !
        !**********************************************************************
        implicit none
        
        integer(INTEGER_TYPE), intent(in) :: IParticle, SolidPointsElement, LiquidPointsElement
        real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
        real(REAL_TYPE), intent(inout) :: WeiGP

        call InitialLocalMaterialPointCoordinatesPointer(IParticle, SolidPointsElement, LiquidPointsElement, WeiGP, PosGP)

        end subroutine InitialLocalMaterialPointCoordinates
        
        subroutine InitialiseMaterialPointPrescribedVelocity()
        !**********************************************************************
        !
        ! Function:  assign prescribed velocity to MP
        !
        !**********************************************************************
        implicit none
        integer(INTEGER_TYPE) :: J, ElID, ParticleIndex, IPart, Dir
        real(REAL_TYPE) :: velocity(NVECTOR)
        logical :: PrescribedDirection(NVECTOR)
        
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary
        
        PrescribedDirection(:) = .false.
        do J = 1, GeoParams%PrescribedVeloNElem
          ElID = GeoParams%PrescribedVeloElemID(J)
          If (ElementMaterialID(ElID,IPatch_Temporary) /= CalParams%MovingMesh%MovingMaterialID) then
            call GiveWarning('The material on which a prescribed velocity is assigned and the reference material for moving mesh should be the same')
          end if   
          
          velocity (1:NVECTOR) = GeoParams%PrescribedVeloElValue(J,1:NVECTOR)
          
          do Dir = 1, NVECTOR
            if (GeoParams%PrescribedVeloElDirection(J,Dir)==1) then  
              PrescribedDirection(Dir) = .true.
            end if 
          end do  
          do IPart = 1, NPartEle(ElID, IPatch_Temporary)
            ParticleIndex = GetParticleIndex(IPart,ElID, IPatch_Temporary)
            Particles(ParticleIndex)%MPPrescribedVeloDir(:) = PrescribedDirection(:)
            Particles(ParticleIndex)%PrescrVelo(:) = velocity(:)
          end do
        end do        
        
		end subroutine InitialiseMaterialPointPrescribedVelocity
		
		subroutine InitialiseVelocityonMP()
		!*************************************************************************************
        !    SUBROUTINE: InitialiseVelocityonMP
        ! 
        !    DESCRIPTION:
        !>   Initial velocity is set to the material point level
        !
        !>   @note : After A.Yerro
        !
        !
        !************************************************************************************
		implicit none
		integer(INTEGER_TYPE) :: IElm, ElID, IPart, ParticleIndex
		real(REAL_TYPE) :: velocity(NVECTOR)
        ! Multipatch variables   
        integer(INTEGER_TYPE) :: IPatch_Temporary
        integer(INTEGER_TYPE) :: IPatch = 1
		
		if ( .not.CalParams%ApplyInitialVelocityonMP ) RETURN
		if ( CalParams%ApplyQuasiStatic) then
			call GiveError('Quasi-Static convergence analisys can not be performed with initial velocity on material points')
		else
			do IElm= 1, GeoParams%NMovingElements ! loop over number of elements with initial v on MP
				ElID= int(GeoParams%InitialVelocityonMP(IElm, 1))
				do IPart= 1, NPartEle(ElID, IPatch_Temporary) ! loop over particles within element
					ParticleIndex=GetParticleIndex(IPart, ElID, IPatch_Temporary)
					velocity=GeoParams%InitialVelocityonMP(IElm, 2:NVECTOR+1)
					VelocityArray(ParticleIndex, :)= velocity !sets the velocity				
				end do			
			end do
			if ( CalParams%RigidBody%IsRigidBody ) then
				CalParams%RigidBody%Velocity= GeoParams%InitialVelocityonMP(GeoParams%NMovingElements, 2:NVECTOR+1)
			end if
			deallocate(GeoParams%InitialVelocityonMP)
		end if	
		
		end subroutine InitialiseVelocityonMP
        
        subroutine RewriteParticleHydraulicConductivity(ParticleIndex,ISet)
        !**********************************************************************
        !
        !    Function:  Updates hydraulic conductivity (K) of the particle 
        !    as a function of Degree of Saturation (Sr)
        !               
        !
        !**********************************************************************
        implicit none
        ! Local variables
        integer(INTEGER_TYPE), intent(in) :: ParticleIndex, ISet
        
          
            
            if (MatParams(ISet)%HydraulicConductivityCurve==HCC_CONSTANT) then
            ! do nothing, Kw doesn't vary
              else if (MatParams(ISet)%HydraulicConductivityCurve==HCC_HILLEL) then
                  call RewriteParticleHydraulicConductivityHillel(ParticleIndex,ISet)
             else if (MatParams(ISet)%HydraulicConductivityCurve==HCC_MUALEM) then
                  call RewriteParticleHydraulicConductivityMualem(ParticleIndex,ISet)  
              end if
          

        
        end subroutine RewriteParticleHydraulicConductivity
        
        subroutine RewriteParticleHydraulicConductivityHillel(ParticleIndex,ISet)
        !**********************************************************************
        !
        !    Function:  Updates hydraulic conductivity (K) of the particle 
        !    as a function of Degree of Saturation (Sr)
        !   The expression of the Degree of Saturation is the RETENTION CURVE, 
        !   while the expression of K is based on   Hillel (1971)
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 

        implicit none
        ! Local variables
        integer(INTEGER_TYPE), intent (in):: ParticleIndex, ISet
        real(REAL_TYPE) :: Sr, r, ksat
        
        Sr = Particles(ParticleIndex)%DegreeSaturation
        if (Sr >=1.0)  return
        !        !
        r = MatParams(ISet)%rexponentHillel
        ksat = MatParams(ISet)%HydraulicConductivityLiquid !% this is the saturated permeability
        
        Particles(ParticleIndex)%Conductivity = ksat * (Sr**(r))
        
        end subroutine RewriteParticleHydraulicConductivityHillel
        
        subroutine RewriteParticleHydraulicConductivityMualem(ParticleIndex,ISet)
        !**********************************************************************
        !
        !    Function:  Updates hydraulic conductivity (K) of the particle 
        !    as a function of Degree of Saturation (Sr)
        !   The expression of the Degree of Saturation is the RETENTION CURVE, 
        !   while the expression of K is based on   Mualem (1976)
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 

        implicit none
        ! Local variables
        integer(INTEGER_TYPE), intent (in):: ParticleIndex, ISet
        real(REAL_TYPE) :: Sr, ksat, krel
        real(REAL_TYPE) :: Smin, Smax, L
        real(REAL_TYPE) :: n00, w1, w2, w3, w4, w5, w6
        
        Sr = Particles(ParticleIndex)%DegreeSaturation
         if (Sr >=1.0)  return
        !
        ksat = MatParams(ISet)%HydraulicConductivityLiquid
        Smin = MatParams(ISet)%Smin
        Smax = MatParams(ISet)%Smax
        L = MatParams(ISet)%L
        
        n00 = 1.0d0
        w1 = sqrt(Sr)
        w3 = n00/L
        w2 = Sr**(w3)
        w4 = (1-w2)**(L)
        w5 = 1-w4
        w6 = w5*w5
        
        krel = w1 * w6
        Particles(ParticleIndex)%Conductivity = ksat * krel
        
        end subroutine RewriteParticleHydraulicConductivityMualem
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    subroutine Load2D_NURBS(RLoad, NDof, Coord, NInt, NNod, ILoadCon, ElementUponWhichLoadIsApplied)
      !LoadType,
      ! LoadValue, 
      !call Load3D(ExtLoadVector, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, NNodes, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid_NURBS(:, I), LoadType)
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
        integer(INTEGER_TYPE), intent(in) :: NInt, NNod!, LoadType 
        integer(INTEGER_TYPE), dimension(NNod) :: ILoadCon
        integer(INTEGER_TYPE), dimension(Counters%NodTot) :: NDof
        !real(REAL_TYPE), dimension(NNod, NVECTOR), intent(in) ::  LoadValue
        real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Coord       ! Coordinates vector
        real(REAL_TYPE), dimension(Counters%NodTot*NDOFL), intent(inout) :: RLoad      ! Resulting load vector
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
        
         integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, NElemNode3, IElNode1, IElNode2, IElNode3, ElementOfNode1, ElementOfNode2, ElementOfNode3, ISharedElement, INode
         real(REAL_TYPE) :: ScalarProduct
         
         !integer(INTEGER_TYPE) :: ELEMENTUPONWHICHLOADISAPPLIED
      
         ! Multipatch variables 
         integer(INTEGER_TYPE) :: IPatch_Temporary = 1
         
         ! MP Traction 
         integer(INTEGER_TYPE) :: NParticles
         integer(INTEGER_TYPE) :: IElement_Traction
         integer(INTEGER_TYPE) :: GlobalElementForTraction
         logical, allocatable, dimension(:) :: ParticleStatus
         integer(INTEGER_TYPE) :: IError
         integer(INTEGER_TYPE) :: NLoadedParticles
         
         integer(INTEGER_TYPE) :: IParticle
		 integer(INTEGER_TYPE) :: ParticleIndex
		 integer(INTEGER_TYPE) :: ParticleID
         real(REAL_TYPE), dimension(NVECTOR) :: TotalSurfaceLoad, LoadPerParticle
         
         !integer(INTEGER_TYPE) :: INode
         integer(INTEGER_TYPE) :: LoadedParticleCounter
         integer(INTEGER_TYPE) :: IDim
         
         RLoad_DEBUG = 0.0
         !LoadOnNodesConnectivitiesSolid
         !ElementUponWhichLoadIsApplied_NURBS
         !LoadValuesOnNodesSolid_NURBS
         
         
         !$$START_LOAD_ON_NODES_SOLID
         !8
         !2 7 12 17	1	0         1            0
         !27 32 37 42	1	0	  1	       0 		
         !52 57 62 67	1	0	  1	       0
         !77 82 87 91	1	0	  1	       0	
         !7 12 17 22	2	0	  1	       0
         !32 37 42 47	2	0	  1	       0
         !57 62 67 72	2	0	  1	       0
         !82 87 92 97	2	0	  1	       0
         
         
         !GPShapeFunctionBoundary
         !GPShapeFunctionDerivativeBoundary
         
         
        !NNodLOE = NXiKnotOrder(IPatch_Temporary)+1!NNod!/2 ! number of nodes of a lower order element
        
        
        ! ILoadConLocal needs to have 4 components 
        !ILoadConLocal(1) = ILoadCon(1)
        !ILoadConLocal(2) = ILoadCon(2)
        !ILoadConLocal(3) = ILoadCon(3)
        !ILoadConLocal(4) = ILoadCon(4)
        
        
       ! TractionNodes = LoadValue
        
        
        !if (NNod == 6) Then
        !  ! corners first, mid-nodes later
        !  ! 6-noded triangle
        !  ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
        !  ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
        !  ILoadConLocal(3) = ILoadCon(5); TractionNodes(3,:) = LoadValue(5,:)
        !  ILoadConLocal(4) = ILoadCon(2); TractionNodes(4,:) = LoadValue(2,:)
        !  ILoadConLocal(5) = ILoadCon(4); TractionNodes(5,:) = LoadValue(4,:)
        !  ILoadConLocal(6) = ILoadCon(6); TractionNodes(6,:) = LoadValue(6,:)
        !
        !else
        !  ! corners first, mid-nodes later
        !  ! 8-noded rectangle
        !  ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
        !  ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
        !  ILoadConLocal(3) = ILoadCon(5); TractionNodes(3,:) = LoadValue(5,:)
        !  ILoadConLocal(4) = ILoadCon(7); TractionNodes(4,:) = LoadValue(7,:)
        !  ILoadConLocal(5) = ILoadCon(2); TractionNodes(5,:) = LoadValue(2,:)
        !  ILoadConLocal(6) = ILoadCon(4); TractionNodes(6,:) = LoadValue(4,:)
        !  ILoadConLocal(7) = ILoadCon(6); TractionNodes(7,:) = LoadValue(6,:)
        !  ILoadConLocal(8) = ILoadCon(8); TractionNodes(8,:) = LoadValue(8,:)
        !  !'no 8-noded loads yet'
        !end if
        
        
        
        
        
        
        ! set arrays to zero
        TractionGaussPoint = 0d0
        SE = 0d0
        
        
!1) we need local GaussPointLocalCoordinates in 1D
!2) we have the shape functions evaluated there
!3)  
        
        
        
        
      
        do I = 1, ELEMENTBOUNDARYGAUSSPOINTS
            
            call NormalOnLINE_NURBS(NodalCoordinates_Traction, GPShapeFunctionDerivative_AllElements(ElementUponWhichLoadIsApplied,I,:,:), &
                NormalVector, VectorLength, ElementUponWhichLoadIsApplied)
            
            !call Normal_T3( I, Coord, ILoadConLocal, NNodLOE, &  
            !                GPShapeFunctionDerivative_AllElements(ElementUponWhichLoadIsApplied,:,:,:), &
            !                NormalVector, VectorLength, ElementUponWhichLoadIsApplied)
          
            !GPShapeFunctionDerivativeBoundary_ZETA_AllElements(ElementUponWhichLoadIsApplied,:,:,:), &
            
          ! determine vector normal to a plane (for 3-noded triangles)
          ! vector length = det(Jacobian_S)
            
          !(ElementUponWhichLoadIsApplied,:,:)
          
          !if((LoadType == LOADTYPE_LIQUID) .or. (LoadType == LOADTYPE_GAS)) then ! Liquid or Gas, Loading is applied only normal to the surface
          !    AppliedPressure(:) = LoadValue(:,1) !Liquid pressure is hydrostatic, LoadValue(:,1)=LoadValue(:,2)...
          ! !*** se the normal pointing outward of the body **** 
          !    !IMPROVE IMPLEMENTATION!!! ONLY WORKS FOR TETRA
          !    
          !    NElemNode1 = GetNElmOfNode(ILoadConLocal(1), IPatch_Temporary)
          !    NElemNode2 = GetNElmOfNode(ILoadConLocal(2), IPatch_Temporary)
          !    NElemNode3 = GetNElmOfNode(ILoadConLocal(3), IPatch_Temporary)
          !     
          !    do IElNode1 = 1,NElemNode1
          !      ElementOfNode1 = GetElmIOfNode(ILoadConLocal(1),IElNode1, IPatch_Temporary)
          !      do IElNode2 = 1,NElemNode2
          !        ElementOfNode2 = GetElmIOfNode(ILoadConLocal(2),IElNode2, IPatch_Temporary)
          !          if (ElementOfNode1==ElementOfNode2) then !found element connected to line
          !            do IElNode3 = 1,NElemNode3
          !              ElementOfNode3 = GetElmIOfNode(ILoadConLocal(3),IElNode3, IPatch_Temporary)
          !              if (ElementOfNode1==ElementOfNode3) then !found element connected to line
          !                ISharedElement = ElementOfNode3 !also = ElementOfNode1 =ElementOfNode2
          !                if (IsActiveElement(ISharedElement, IPatch_Temporary)) then 
          !                   do K=1,ELEMENTNODES
          !                       INode =ElementConnectivities(K,ISharedElement,IPatch_Temporary)
          !                       if ((INode/=ILoadConLocal(1)).and.(INode/=ILoadConLocal(2)).and.(INode/=ILoadConLocal(3))) then !this is the node not belonging to the line
          !                           Vector(:) = Coord(INode,:) - Coord(ILoadConLocal(1),:)
          !                           ScalarProduct = DotProduct(Vector,NormalVector,NVECTOR)
          !                           if (ScalarProduct > 0) then ! The normal is pointing inside element, reverse
          !                              NormalVector = -NormalVector !the normal is pointing outside active element
          !                              EXIT !no need to check the other nodes of the element
          !                           elseif (ScalarProduct < 0) then ! The normal is pointing outside element, keep direction
          !                              EXIT !no need to check the other nodes of the element   
          !                           end if
          !                       end if
          !                   end do
          !                   EXIT !No need to check the other elements because already found one which is active
          !                end if
          !                  
          !              end if
          !            end do 
          !          end if
          !        end do
          !      end do
          !
          !  TractionNodes(4:6,1:3) = 0d0 ! remove middle point data
          !  ProductTractionNodesNormalVector = 0d0 ! initialize
          !
          !  do K = 1, NNodLOE ! for each corner node
          !    do J=1,NVECTOR  
          !      TractionNodes(K,J) = AppliedPressure(K) * NormalVector (J)
          !    end do
          !    !do J = 1, NVECTOR ! for each direction
          !    !   ProductTractionNodesNormalVector(K) = ProductTractionNodesNormalVector(K) + TractionNodes(K,j) * NormalVector (J)
          !    !end do
          !    !
          !    !  do J = 1, NVECTOR
          !    !     TractionNodes(K, J) = ProductTractionNodesNormalVector(K) * NormalVector(J)
          !    !  end do
          !    !  ScalarProduct = DotProduct(TractionNodes(K,:),NormalVector,NVECTOR) !>0 if traction, <0 if compression, normal pointing outward
          !    !  if (sign(1.d0,ScalarProduct) /= sign(1.d0,LoadValue(K,1))) then !LoadValue(K,1)=LoadValue(K,2)>0 for traction
          !    !    TractionNodes(K, :) = - TractionNodes(K, :)
          !    !  end if
          !  end do
          !end if
                    
          ! calculate traction at Gauss Point
          !do K = 1, ELEMENTBOUNDARYNODES !NNodLOE ! Corner nodes
          !  do J = 1, 3!ELEMENTBOUNDARYNODES!3
          !    NN = IEN_TRACTION(K,ElementUponWhichLoadIsApplied) !ILoadConLocal(K)
          !    TractionGaussPoint(K,J) = TractionGaussPoint(K,J) + GPShapeFunction_AllElements(ElementUponWhichLoadIsApplied,I,K) * NormalVector(J)!TractionNodes(J,1)
          !  end do
          !end do
          
          ! calculate traction at Nodes
          do K = 1, ELEMENTBOUNDARYNODES!NNodLOE ! Corner nodes
            do J = 1, NVECTOR
            
                !ND = IEN_TRACTION(K,ElementUponWhichLoadIsApplied)!ILoadConLocal(K)
                ! HARDCODED
                !if ( (J==1) .or. (J==2)) then 
                
                Temp = GPWeight_AllElements(ElementUponWhichLoadIsApplied,I) * GPShapeFunction_AllElements(ElementUponWhichLoadIsApplied, I, K) * VectorLength * NormalVector(J)
                SE(K,J) = SE(K,J) + ( Temp )!* TractionGaussPoint(K, J) )
                
                
              !Temp =   VectorLength * GPShapeFunctionBoundary_XI_AllElements(ElementUponWhichLoadIsApplied,I,K)!GPShapeFunctionBoundary(I,K) !TractionNodes(J,1)!GPShapeFunctionBoundary(I, K) * !GPWeightBoundary(I) *
              !SE(K, J) = SE(K, J) + (Temp* (NormalVector(J)/ sqrt(NormalVector(1)**2 +NormalVector(2)**2 +NormalVector(3)**2 ))) !* TractionGaussPoint(I, J)
              
                !elseif (J==3) then 
                !    SE(K, J) = 0 
                !    
                !end if 
                
                    
                    
              !(NodalCoordinates(ND,J,IPatch_Temporary)/sqrt(NodalCoordinates(ND,1,IPatch_Temporary)**2 + NodalCoordinates(ND,2,IPatch_Temporary)**2 ))
              !
              
            end do 
          end do
          
          
        end do

        
        
        ! THIS IS NEEDED FOR NODAL FORCES TRACTION 
        ! add to global load vector
        do I = 1, ELEMENTBOUNDARYNODES!NNod  ! --> ELEMENTBOUNDARYNODES is the number of nodes in the element (eg. 2*2 if traction on a linear square element)
          ND = NDOF(Traction_Connecting_Local_To_Global_ControlPoints(IEN_TRACTION(I,ElementUponWhichLoadIsApplied)))!NDof(ILoadConLocal(I))
          do J = 1, NVECTOR
            RLoad(ND + J) = RLoad(ND + J) + (SE(I,J)   * -1) ! HARDCODED NEGATIVE 1 
            RLoad_DEBUG = RLoad_DEBUG + RLoad(ND + J)
          end do
        end do
        
        
        
        
        ! THIS IS NEEDED FOR MP TRACTION 
        
        ! I need to map the SE to the material points 
        
        !GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, MPNumber, NodeNumber) 
        !
        ! I need to move this to elsewhere after Load3D
              !NParticles = 0
              !
              !
              !!do IElement_Traction = 1, nel_NURBS_Traction ! traction elements
              !    
              !    ParticleID = 0
              !LoadedParticleCounter = 0
              !    
              !    GlobalElementForTraction = Traction_Connecting_Local_To_Global_Elements(ElementUponWhichLoadIsApplied)
              !    
              !NParticles =  NPartEle(GlobalElementForTraction, IPatch_Temporary) ! HARDCODED
              !
              !allocate(ParticleStatus(NParticles), stat = IError)
              !
              !ParticleStatus = .false.
              !
              !call DetermineAdjacentParticles(Traction_Connecting_Local_To_Global_Sides(ElementUponWhichLoadIsApplied), NParticles, ParticleStatus)
              !
              !NLoadedParticles = DetermineNAdjacentParticles(ParticleStatus, NParticles)
              !
              !!call MarkBoundaryParticles(GlobalElementForTraction, ParticleStatus, NParticles, IPatch_Temporary)
              !do IParticle = 1, NPartEle(GlobalElementForTraction, IPatch_Temporary) ! loop all MPs of the element (solid and liquid MPs in case of double-point formulation)
              !    ParticleIndex = GetParticleIndex(IParticle, GlobalElementForTraction, IPatch_Temporary)
              !
              !    !if ( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then
              !        ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
              !
              !        ParticleID = ParticleID + 1
              !
              !        if(ParticleStatus(ParticleID)) then ! this MP is loaded
              !
              !            LoadedParticleCounter = LoadedParticleCounter + 1
              !            
              !            do INode = 1, ElementBoundaryNodes_Xi
              !            
              !                do IDim = 1, NDIM
              !                
              !                LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1 * 0.5 )
              !                ! HARDCODED HALF
              !                ! HARDCODED MINUS ONE FOR DOWNWARDS
              !                
              !                end do 
              !                
              !                
              !                
              !                
              !                
              !            end do 
              !                            
              !            call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle, 1) !ILoadSystem
              !
              !            LoadPerParticle = 0
              !          !LoadPerParticle = DetermineLoadPerParticle(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal, Counters%NodTot, NodeCoord, Loadtype)
              !          !
              !          !!if ( LoadType == 1 ) then ! assign to solid phase / solid material point
              !          !  call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
              !
              !          !else if ( LoadType == 2 ) then ! assign to water phase / liquid material point
              !          !  call IncreaseFExtWater(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
              !          !
              !          !else if ( LoadType == 3 ) then ! assign to gas phase
              !          !  call IncreaseFExtGas(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
              !          !
              !          !  if ( NFORMULATION == 2 ) then
              !          !    call GiveError('Erroneous definition of applied external load: For double-point formulation only load on solid or liquid phase is allowed.')
              !          !  end if
              !          !
              !          !end if ! load type
              !        end if ! loaded y/n
              !    !end if ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
              !
              !  end do ! loop over MPs
              !deallocate(ParticleStatus, stat = IError)
              
              !end do 
        
        
        
        
        

    end subroutine Load2D_NURBS
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    subroutine Load2D_NURBS_OnMPs(RLoad, NDof, Coord, NInt, NNod, ILoadCon, ElementUponWhichLoadIsApplied)
      !LoadType,
      ! LoadValue, 
      !call Load3D(ExtLoadVector, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, NNodes, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid_NURBS(:, I), LoadType)
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
        integer(INTEGER_TYPE), intent(in) :: NInt, NNod!, LoadType 
        integer(INTEGER_TYPE), dimension(NNod) :: ILoadCon
        integer(INTEGER_TYPE), dimension(Counters%NodTot) :: NDof
        !real(REAL_TYPE), dimension(NNod, NVECTOR), intent(in) ::  LoadValue
        real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Coord       ! Coordinates vector
        real(REAL_TYPE), dimension(Counters%NodTot*NDOFL), intent(inout) :: RLoad      ! Resulting load vector
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
        
         integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, NElemNode3, IElNode1, IElNode2, IElNode3, ElementOfNode1, ElementOfNode2, ElementOfNode3, ISharedElement, INode
         real(REAL_TYPE) :: ScalarProduct
         
         !integer(INTEGER_TYPE) :: ELEMENTUPONWHICHLOADISAPPLIED
      
         ! Multipatch variables 
         integer(INTEGER_TYPE) :: IPatch_Temporary = 1
         
         ! MP Traction 
         integer(INTEGER_TYPE) :: NParticles
         integer(INTEGER_TYPE) :: IElement_Traction
         integer(INTEGER_TYPE) :: GlobalElementForTraction
         logical, allocatable, dimension(:) :: ParticleStatus
         integer(INTEGER_TYPE) :: IError
         integer(INTEGER_TYPE) :: NLoadedParticles
         
         integer(INTEGER_TYPE) :: IParticle
		 integer(INTEGER_TYPE) :: ParticleIndex
		 integer(INTEGER_TYPE) :: ParticleID
         real(REAL_TYPE), dimension(NVECTOR) :: TotalSurfaceLoad, LoadPerParticle
         
         !integer(INTEGER_TYPE) :: INode
         integer(INTEGER_TYPE) :: LoadedParticleCounter
         integer(INTEGER_TYPE) :: IDim
         real(REAL_TYPE) :: Factor
         
         RLoad_DEBUG = 0.0
         !LoadOnNodesConnectivitiesSolid
         !ElementUponWhichLoadIsApplied_NURBS
         !LoadValuesOnNodesSolid_NURBS
         
         
         !$$START_LOAD_ON_NODES_SOLID
         !8
         !2 7 12 17	1	0         1            0
         !27 32 37 42	1	0	  1	       0 		
         !52 57 62 67	1	0	  1	       0
         !77 82 87 91	1	0	  1	       0	
         !7 12 17 22	2	0	  1	       0
         !32 37 42 47	2	0	  1	       0
         !57 62 67 72	2	0	  1	       0
         !82 87 92 97	2	0	  1	       0
         
         
         !GPShapeFunctionBoundary
         !GPShapeFunctionDerivativeBoundary
         
         
        !NNodLOE = NXiKnotOrder(IPatch_Temporary)+1!NNod!/2 ! number of nodes of a lower order element
        
        
        ! ILoadConLocal needs to have 4 components 
        !ILoadConLocal(1) = ILoadCon(1)
        !ILoadConLocal(2) = ILoadCon(2)
        !ILoadConLocal(3) = ILoadCon(3)
        !ILoadConLocal(4) = ILoadCon(4)
        
        
       ! TractionNodes = LoadValue
        
        
        !if (NNod == 6) Then
        !  ! corners first, mid-nodes later
        !  ! 6-noded triangle
        !  ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
        !  ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
        !  ILoadConLocal(3) = ILoadCon(5); TractionNodes(3,:) = LoadValue(5,:)
        !  ILoadConLocal(4) = ILoadCon(2); TractionNodes(4,:) = LoadValue(2,:)
        !  ILoadConLocal(5) = ILoadCon(4); TractionNodes(5,:) = LoadValue(4,:)
        !  ILoadConLocal(6) = ILoadCon(6); TractionNodes(6,:) = LoadValue(6,:)
        !
        !else
        !  ! corners first, mid-nodes later
        !  ! 8-noded rectangle
        !  ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
        !  ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
        !  ILoadConLocal(3) = ILoadCon(5); TractionNodes(3,:) = LoadValue(5,:)
        !  ILoadConLocal(4) = ILoadCon(7); TractionNodes(4,:) = LoadValue(7,:)
        !  ILoadConLocal(5) = ILoadCon(2); TractionNodes(5,:) = LoadValue(2,:)
        !  ILoadConLocal(6) = ILoadCon(4); TractionNodes(6,:) = LoadValue(4,:)
        !  ILoadConLocal(7) = ILoadCon(6); TractionNodes(7,:) = LoadValue(6,:)
        !  ILoadConLocal(8) = ILoadCon(8); TractionNodes(8,:) = LoadValue(8,:)
        !  !'no 8-noded loads yet'
        !end if
        
        
        
        
        
        
        ! set arrays to zero
        TractionGaussPoint = 0d0
        SE = 0d0
        
        
!1) we need local GaussPointLocalCoordinates in 1D
!2) we have the shape functions evaluated there
!3)  
        
        
        
        
      
        do I = 1, ELEMENTBOUNDARYGAUSSPOINTS
            
            call NormalOnLINE_NURBS(NodalCoordinates_Traction, GPShapeFunctionDerivative_AllElements(ElementUponWhichLoadIsApplied,I,:,:), &
                NormalVector, VectorLength, ElementUponWhichLoadIsApplied)
            
            !call Normal_T3( I, Coord, ILoadConLocal, NNodLOE, &  
            !                GPShapeFunctionDerivative_AllElements(ElementUponWhichLoadIsApplied,:,:,:), &
            !                NormalVector, VectorLength, ElementUponWhichLoadIsApplied)
          
            !GPShapeFunctionDerivativeBoundary_ZETA_AllElements(ElementUponWhichLoadIsApplied,:,:,:), &
            
          ! determine vector normal to a plane (for 3-noded triangles)
          ! vector length = det(Jacobian_S)
            
          !(ElementUponWhichLoadIsApplied,:,:)
          
          !if((LoadType == LOADTYPE_LIQUID) .or. (LoadType == LOADTYPE_GAS)) then ! Liquid or Gas, Loading is applied only normal to the surface
          !    AppliedPressure(:) = LoadValue(:,1) !Liquid pressure is hydrostatic, LoadValue(:,1)=LoadValue(:,2)...
          ! !*** se the normal pointing outward of the body **** 
          !    !IMPROVE IMPLEMENTATION!!! ONLY WORKS FOR TETRA
          !    
          !    NElemNode1 = GetNElmOfNode(ILoadConLocal(1), IPatch_Temporary)
          !    NElemNode2 = GetNElmOfNode(ILoadConLocal(2), IPatch_Temporary)
          !    NElemNode3 = GetNElmOfNode(ILoadConLocal(3), IPatch_Temporary)
          !     
          !    do IElNode1 = 1,NElemNode1
          !      ElementOfNode1 = GetElmIOfNode(ILoadConLocal(1),IElNode1, IPatch_Temporary)
          !      do IElNode2 = 1,NElemNode2
          !        ElementOfNode2 = GetElmIOfNode(ILoadConLocal(2),IElNode2, IPatch_Temporary)
          !          if (ElementOfNode1==ElementOfNode2) then !found element connected to line
          !            do IElNode3 = 1,NElemNode3
          !              ElementOfNode3 = GetElmIOfNode(ILoadConLocal(3),IElNode3, IPatch_Temporary)
          !              if (ElementOfNode1==ElementOfNode3) then !found element connected to line
          !                ISharedElement = ElementOfNode3 !also = ElementOfNode1 =ElementOfNode2
          !                if (IsActiveElement(ISharedElement, IPatch_Temporary)) then 
          !                   do K=1,ELEMENTNODES
          !                       INode =ElementConnectivities(K,ISharedElement,IPatch_Temporary)
          !                       if ((INode/=ILoadConLocal(1)).and.(INode/=ILoadConLocal(2)).and.(INode/=ILoadConLocal(3))) then !this is the node not belonging to the line
          !                           Vector(:) = Coord(INode,:) - Coord(ILoadConLocal(1),:)
          !                           ScalarProduct = DotProduct(Vector,NormalVector,NVECTOR)
          !                           if (ScalarProduct > 0) then ! The normal is pointing inside element, reverse
          !                              NormalVector = -NormalVector !the normal is pointing outside active element
          !                              EXIT !no need to check the other nodes of the element
          !                           elseif (ScalarProduct < 0) then ! The normal is pointing outside element, keep direction
          !                              EXIT !no need to check the other nodes of the element   
          !                           end if
          !                       end if
          !                   end do
          !                   EXIT !No need to check the other elements because already found one which is active
          !                end if
          !                  
          !              end if
          !            end do 
          !          end if
          !        end do
          !      end do
          !
          !  TractionNodes(4:6,1:3) = 0d0 ! remove middle point data
          !  ProductTractionNodesNormalVector = 0d0 ! initialize
          !
          !  do K = 1, NNodLOE ! for each corner node
          !    do J=1,NVECTOR  
          !      TractionNodes(K,J) = AppliedPressure(K) * NormalVector (J)
          !    end do
          !    !do J = 1, NVECTOR ! for each direction
          !    !   ProductTractionNodesNormalVector(K) = ProductTractionNodesNormalVector(K) + TractionNodes(K,j) * NormalVector (J)
          !    !end do
          !    !
          !    !  do J = 1, NVECTOR
          !    !     TractionNodes(K, J) = ProductTractionNodesNormalVector(K) * NormalVector(J)
          !    !  end do
          !    !  ScalarProduct = DotProduct(TractionNodes(K,:),NormalVector,NVECTOR) !>0 if traction, <0 if compression, normal pointing outward
          !    !  if (sign(1.d0,ScalarProduct) /= sign(1.d0,LoadValue(K,1))) then !LoadValue(K,1)=LoadValue(K,2)>0 for traction
          !    !    TractionNodes(K, :) = - TractionNodes(K, :)
          !    !  end if
          !  end do
          !end if
                    
          ! calculate traction at Gauss Point
          !do K = 1, ELEMENTBOUNDARYNODES !NNodLOE ! Corner nodes
          !  do J = 1, 3!ELEMENTBOUNDARYNODES!3
          !    NN = IEN_TRACTION(K,ElementUponWhichLoadIsApplied) !ILoadConLocal(K)
          !    TractionGaussPoint(K,J) = TractionGaussPoint(K,J) + GPShapeFunction_AllElements(ElementUponWhichLoadIsApplied,I,K) * NormalVector(J)!TractionNodes(J,1)
          !  end do
          !end do
          
          ! calculate traction at Nodes
          do K = 1, ELEMENTBOUNDARYNODES!NNodLOE ! Corner nodes
            do J = 1, NVECTOR
            
                !ND = IEN_TRACTION(K,ElementUponWhichLoadIsApplied)!ILoadConLocal(K)
                ! HARDCODED
                !if ( (J==1) .or. (J==2)) then 
                
                Temp = GPWeight_AllElements(ElementUponWhichLoadIsApplied,I) * GPShapeFunction_AllElements(ElementUponWhichLoadIsApplied, I, K) * VectorLength * NormalVector(J)
                SE(K,J) = SE(K,J) + ( Temp )!* TractionGaussPoint(K, J) )
                
                
              !Temp =   VectorLength * GPShapeFunctionBoundary_XI_AllElements(ElementUponWhichLoadIsApplied,I,K)!GPShapeFunctionBoundary(I,K) !TractionNodes(J,1)!GPShapeFunctionBoundary(I, K) * !GPWeightBoundary(I) *
              !SE(K, J) = SE(K, J) + (Temp* (NormalVector(J)/ sqrt(NormalVector(1)**2 +NormalVector(2)**2 +NormalVector(3)**2 ))) !* TractionGaussPoint(I, J)
              
                !elseif (J==3) then 
                !    SE(K, J) = 0 
                !    
                !end if 
                
                    
                    
              !(NodalCoordinates(ND,J,IPatch_Temporary)/sqrt(NodalCoordinates(ND,1,IPatch_Temporary)**2 + NodalCoordinates(ND,2,IPatch_Temporary)**2 ))
              !
              
            end do 
          end do
          
          
        end do

        
        
        ! THIS IS NEEDED FOR NODAL FORCES TRACTION 
        ! add to global load vector
        do I = 1, ELEMENTBOUNDARYNODES!NNod  ! --> ELEMENTBOUNDARYNODES is the number of nodes in the element (eg. 2*2 if traction on a linear square element)
          ND = NDOF(Traction_Connecting_Local_To_Global_ControlPoints(IEN_TRACTION(I,ElementUponWhichLoadIsApplied)))!NDof(ILoadConLocal(I))
          do J = 1, NVECTOR
            !RLoad(ND + J) = RLoad(ND + J) + (SE(I,J)   * -1) ! HARDCODED NEGATIVE 1 
            RLoad_DEBUG = RLoad_DEBUG + (SE(I,J)   * -1) !RLoad(ND + J)
          end do
        end do
        
        
        
        
        ! THIS IS NEEDED FOR MP TRACTION 
        
        ! I need to map the SE to the material points 
        
        !GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, MPNumber, NodeNumber) 
        !
        ! I need to move this to elsewhere after Load3D
        
        !Traction_Connecting_Local_To_Global_Elements(1) = 91
        !Traction_Connecting_Local_To_Global_Elements(2) = 92
        !Traction_Connecting_Local_To_Global_Elements(3) = 93 
        !Traction_Connecting_Local_To_Global_Elements(4) = 94
        !Traction_Connecting_Local_To_Global_Elements(5) = 95
        !Traction_Connecting_Local_To_Global_Elements(6) = 96
        !Traction_Connecting_Local_To_Global_Elements(7) = 97
        !Traction_Connecting_Local_To_Global_Elements(8) = 98
        !Traction_Connecting_Local_To_Global_Elements(9) = 99
        !Traction_Connecting_Local_To_Global_Elements(10) = 100
        
        
        
        
        
        
        
        Factor = real(ELEMENTBOUNDARYNODES) / real(ELEMENTBOUNDARYGAUSSPOINTS)
        
        
              NParticles = 0
              
              
              !do IElement_Traction = 1, nel_NURBS_Traction ! traction elements
                  
                  ParticleID = 0
              LoadedParticleCounter = 0
                  
                  GlobalElementForTraction = Traction_Connecting_Local_To_Global_Elements(ElementUponWhichLoadIsApplied)
                  
              NParticles =  NPartEle(GlobalElementForTraction, IPatch_Temporary) ! HARDCODED
              
              allocate(ParticleStatus(NParticles), stat = IError)
              
              ParticleStatus = .false.
              
              call DetermineAdjacentParticles(Traction_Connecting_Local_To_Global_Sides(ElementUponWhichLoadIsApplied), NParticles, ParticleStatus)
              
              ! side number hardcoded              
              !call DetermineAdjacentParticles(1, NParticles, ParticleStatus)

              
              NLoadedParticles = DetermineNAdjacentParticles(ParticleStatus, NParticles)
              
              !call MarkBoundaryParticles(GlobalElementForTraction, ParticleStatus, NParticles, IPatch_Temporary)
              do IParticle = 1, NPartEle(GlobalElementForTraction, IPatch_Temporary) ! loop all MPs of the element (solid and liquid MPs in case of double-point formulation)
                  ParticleIndex = GetParticleIndex(IParticle, GlobalElementForTraction, IPatch_Temporary)
              
                  !if ( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then
                      ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
              
                      ParticleID = ParticleID + 1
              
                      if(ParticleStatus(ParticleID)) then ! this MP is loaded
              
                          LoadedParticleCounter = LoadedParticleCounter + 1
                          
                          do INode = 1, ElementBoundaryNodes_Xi
                          
                              do IDim = 1, NDIM
                              
                              !LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1  * Factor  ) !(2/3)
                              LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( SE(INode,IDim) * -1 / real(ELEMENTBOUNDARYGAUSSPOINTS) ) !( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1  * Factor  ) !(2/3)
                              
                              ! HARDCODED HALF * 0.5
                              ! HARDCODED MINUS ONE FOR DOWNWARDS
                              
                              end do 
                              
                              
                              
                              
                              
                          end do 
                                          
                          call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle, 1) !ILoadSystem
              
                          LoadPerParticle = 0
                        !LoadPerParticle = DetermineLoadPerParticle(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal, Counters%NodTot, NodeCoord, Loadtype)
                        !
                        !!if ( LoadType == 1 ) then ! assign to solid phase / solid material point
                        !  call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
              
                        !else if ( LoadType == 2 ) then ! assign to water phase / liquid material point
                        !  call IncreaseFExtWater(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
                        !
                        !else if ( LoadType == 3 ) then ! assign to gas phase
                        !  call IncreaseFExtGas(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
                        !
                        !  if ( NFORMULATION == 2 ) then
                        !    call GiveError('Erroneous definition of applied external load: For double-point formulation only load on solid or liquid phase is allowed.')
                        !  end if
                        !
                        !end if ! load type
                      end if ! loaded y/n
                  !end if ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
              
                end do ! loop over MPs
              deallocate(ParticleStatus, stat = IError)
              
              !end do 
        
        
        
        
        

    end subroutine Load2D_NURBS_OnMPs
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !------------------------------------------------------------------------------------------------
    ! TRACTION ON THE MPs 
    
    subroutine Load3D_NURBS_OnMPs(RLoad, NDof, Coord, NInt, NNod, ILoadCon, ElementUponWhichLoadIsApplied)
      !LoadType,
      ! LoadValue, 
      !call Load3D(ExtLoadVector, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, NNodes, LoadOnNodesConnectivitiesSolid(1, I), LoadValuesOnNodesSolid_NURBS(:, I), LoadType)
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
        integer(INTEGER_TYPE), intent(in) :: NInt, NNod!, LoadType 
        integer(INTEGER_TYPE), dimension(NNod) :: ILoadCon
        integer(INTEGER_TYPE), dimension(Counters%NodTot) :: NDof
        !real(REAL_TYPE), dimension(NNod, NVECTOR), intent(in) ::  LoadValue
        real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Coord       ! Coordinates vector
        real(REAL_TYPE), dimension(Counters%NodTot*NDOFL), intent(inout) :: RLoad      ! Resulting load vector
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
        
         integer(INTEGER_TYPE) :: NElemNode1, NElemNode2, NElemNode3, IElNode1, IElNode2, IElNode3, ElementOfNode1, ElementOfNode2, ElementOfNode3, ISharedElement, INode
         real(REAL_TYPE) :: ScalarProduct
         
         !integer(INTEGER_TYPE) :: ELEMENTUPONWHICHLOADISAPPLIED
      
         ! Multipatch variables 
         integer(INTEGER_TYPE) :: IPatch_Temporary = 1
         
         
         
         
         integer(INTEGER_TYPE) :: LoadedParticleCounter
         integer(INTEGER_TYPE) :: IDim
         real(REAL_TYPE) :: Factor
         integer(INTEGER_TYPE) :: NParticles, ParticleID, GlobalElementForTraction, IError
         integer(INTEGER_TYPE) :: NLoadedParticles, IParticle, ParticleIndex
                  
         logical, allocatable, dimension(:) :: ParticleStatus

         
         real(REAL_TYPE), dimension(NVECTOR) :: LoadPerParticle
         
         RLoad_DEBUG = 0.0
         !LoadOnNodesConnectivitiesSolid
         !ElementUponWhichLoadIsApplied_NURBS
         !LoadValuesOnNodesSolid_NURBS
         
         
         !$$START_LOAD_ON_NODES_SOLID
         !8
         !2 7 12 17	1	0         1            0
         !27 32 37 42	1	0	  1	       0 		
         !52 57 62 67	1	0	  1	       0
         !77 82 87 91	1	0	  1	       0	
         !7 12 17 22	2	0	  1	       0
         !32 37 42 47	2	0	  1	       0
         !57 62 67 72	2	0	  1	       0
         !82 87 92 97	2	0	  1	       0
         
         
         !GPShapeFunctionBoundary
         !GPShapeFunctionDerivativeBoundary
         
         
        !NNodLOE = NXiKnotOrder(IPatch_Temporary)+1!NNod!/2 ! number of nodes of a lower order element
        
        
        ! ILoadConLocal needs to have 4 components 
        !ILoadConLocal(1) = ILoadCon(1)
        !ILoadConLocal(2) = ILoadCon(2)
        !ILoadConLocal(3) = ILoadCon(3)
        !ILoadConLocal(4) = ILoadCon(4)
        
        
       ! TractionNodes = LoadValue
        
        
        !if (NNod == 6) Then
        !  ! corners first, mid-nodes later
        !  ! 6-noded triangle
        !  ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
        !  ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
        !  ILoadConLocal(3) = ILoadCon(5); TractionNodes(3,:) = LoadValue(5,:)
        !  ILoadConLocal(4) = ILoadCon(2); TractionNodes(4,:) = LoadValue(2,:)
        !  ILoadConLocal(5) = ILoadCon(4); TractionNodes(5,:) = LoadValue(4,:)
        !  ILoadConLocal(6) = ILoadCon(6); TractionNodes(6,:) = LoadValue(6,:)
        !
        !else
        !  ! corners first, mid-nodes later
        !  ! 8-noded rectangle
        !  ILoadConLocal(1) = ILoadCon(1); TractionNodes(1,:) = LoadValue(1,:)
        !  ILoadConLocal(2) = ILoadCon(3); TractionNodes(2,:) = LoadValue(3,:)
        !  ILoadConLocal(3) = ILoadCon(5); TractionNodes(3,:) = LoadValue(5,:)
        !  ILoadConLocal(4) = ILoadCon(7); TractionNodes(4,:) = LoadValue(7,:)
        !  ILoadConLocal(5) = ILoadCon(2); TractionNodes(5,:) = LoadValue(2,:)
        !  ILoadConLocal(6) = ILoadCon(4); TractionNodes(6,:) = LoadValue(4,:)
        !  ILoadConLocal(7) = ILoadCon(6); TractionNodes(7,:) = LoadValue(6,:)
        !  ILoadConLocal(8) = ILoadCon(8); TractionNodes(8,:) = LoadValue(8,:)
        !  !'no 8-noded loads yet'
        !end if
        
        
        
        
        
        
        ! set arrays to zero
        TractionGaussPoint = 0d0
        SE = 0d0
        
        
!1) we need local GaussPointLocalCoordinates in 1D
!2) we have the shape functions evaluated there
!3)  
        
        
        
        
      
        do I = 1, ELEMENTBOUNDARYGAUSSPOINTS
            
            call NormalOnQUAD_NURBS(NodalCoordinates_Traction, GPShapeFunctionDerivative_AllElements(ElementUponWhichLoadIsApplied,I,:,:), &
                NormalVector, VectorLength, ElementUponWhichLoadIsApplied)
            
            !call Normal_T3( I, Coord, ILoadConLocal, NNodLOE, &  
            !                GPShapeFunctionDerivative_AllElements(ElementUponWhichLoadIsApplied,:,:,:), &
            !                NormalVector, VectorLength, ElementUponWhichLoadIsApplied)
          
            !GPShapeFunctionDerivativeBoundary_ZETA_AllElements(ElementUponWhichLoadIsApplied,:,:,:), &
            
          ! determine vector normal to a plane (for 3-noded triangles)
          ! vector length = det(Jacobian_S)
            
          !(ElementUponWhichLoadIsApplied,:,:)
          
          !if((LoadType == LOADTYPE_LIQUID) .or. (LoadType == LOADTYPE_GAS)) then ! Liquid or Gas, Loading is applied only normal to the surface
          !    AppliedPressure(:) = LoadValue(:,1) !Liquid pressure is hydrostatic, LoadValue(:,1)=LoadValue(:,2)...
          ! !*** se the normal pointing outward of the body **** 
          !    !IMPROVE IMPLEMENTATION!!! ONLY WORKS FOR TETRA
          !    
          !    NElemNode1 = GetNElmOfNode(ILoadConLocal(1), IPatch_Temporary)
          !    NElemNode2 = GetNElmOfNode(ILoadConLocal(2), IPatch_Temporary)
          !    NElemNode3 = GetNElmOfNode(ILoadConLocal(3), IPatch_Temporary)
          !     
          !    do IElNode1 = 1,NElemNode1
          !      ElementOfNode1 = GetElmIOfNode(ILoadConLocal(1),IElNode1, IPatch_Temporary)
          !      do IElNode2 = 1,NElemNode2
          !        ElementOfNode2 = GetElmIOfNode(ILoadConLocal(2),IElNode2, IPatch_Temporary)
          !          if (ElementOfNode1==ElementOfNode2) then !found element connected to line
          !            do IElNode3 = 1,NElemNode3
          !              ElementOfNode3 = GetElmIOfNode(ILoadConLocal(3),IElNode3, IPatch_Temporary)
          !              if (ElementOfNode1==ElementOfNode3) then !found element connected to line
          !                ISharedElement = ElementOfNode3 !also = ElementOfNode1 =ElementOfNode2
          !                if (IsActiveElement(ISharedElement, IPatch_Temporary)) then 
          !                   do K=1,ELEMENTNODES
          !                       INode =ElementConnectivities(K,ISharedElement,IPatch_Temporary)
          !                       if ((INode/=ILoadConLocal(1)).and.(INode/=ILoadConLocal(2)).and.(INode/=ILoadConLocal(3))) then !this is the node not belonging to the line
          !                           Vector(:) = Coord(INode,:) - Coord(ILoadConLocal(1),:)
          !                           ScalarProduct = DotProduct(Vector,NormalVector,NVECTOR)
          !                           if (ScalarProduct > 0) then ! The normal is pointing inside element, reverse
          !                              NormalVector = -NormalVector !the normal is pointing outside active element
          !                              EXIT !no need to check the other nodes of the element
          !                           elseif (ScalarProduct < 0) then ! The normal is pointing outside element, keep direction
          !                              EXIT !no need to check the other nodes of the element   
          !                           end if
          !                       end if
          !                   end do
          !                   EXIT !No need to check the other elements because already found one which is active
          !                end if
          !                  
          !              end if
          !            end do 
          !          end if
          !        end do
          !      end do
          !
          !  TractionNodes(4:6,1:3) = 0d0 ! remove middle point data
          !  ProductTractionNodesNormalVector = 0d0 ! initialize
          !
          !  do K = 1, NNodLOE ! for each corner node
          !    do J=1,NVECTOR  
          !      TractionNodes(K,J) = AppliedPressure(K) * NormalVector (J)
          !    end do
          !    !do J = 1, NVECTOR ! for each direction
          !    !   ProductTractionNodesNormalVector(K) = ProductTractionNodesNormalVector(K) + TractionNodes(K,j) * NormalVector (J)
          !    !end do
          !    !
          !    !  do J = 1, NVECTOR
          !    !     TractionNodes(K, J) = ProductTractionNodesNormalVector(K) * NormalVector(J)
          !    !  end do
          !    !  ScalarProduct = DotProduct(TractionNodes(K,:),NormalVector,NVECTOR) !>0 if traction, <0 if compression, normal pointing outward
          !    !  if (sign(1.d0,ScalarProduct) /= sign(1.d0,LoadValue(K,1))) then !LoadValue(K,1)=LoadValue(K,2)>0 for traction
          !    !    TractionNodes(K, :) = - TractionNodes(K, :)
          !    !  end if
          !  end do
          !end if
                    
          ! calculate traction at Gauss Point
          !do K = 1, ELEMENTBOUNDARYNODES !NNodLOE ! Corner nodes
          !  do J = 1, 3!ELEMENTBOUNDARYNODES!3
          !    NN = IEN_TRACTION(K,ElementUponWhichLoadIsApplied) !ILoadConLocal(K)
          !    TractionGaussPoint(K,J) = TractionGaussPoint(K,J) + GPShapeFunction_AllElements(ElementUponWhichLoadIsApplied,I,K) * NormalVector(J)!TractionNodes(J,1)
          !  end do
          !end do
          
          ! calculate traction at Nodes
          do K = 1, ELEMENTBOUNDARYNODES!NNodLOE ! Corner nodes
            do J = 1, NVECTOR
            
                !ND = IEN_TRACTION(K,ElementUponWhichLoadIsApplied)!ILoadConLocal(K)
                ! HARDCODED
                !if ( (J==1) .or. (J==2)) then 
                
                Temp = GPWeight_AllElements(ElementUponWhichLoadIsApplied,I) * GPShapeFunction_AllElements(ElementUponWhichLoadIsApplied, I, K) * VectorLength * NormalVector(J) * -1 ! hardcoded negative sign to make traction go inwards
                SE(K,J) = SE(K,J) + ( Temp )!* TractionGaussPoint(K, J) )
          
                
                ! Normal outside factor --> -1 if it is pointing inwards 
                ! Normal outside factor -->  1 if it is pointing outwards
                
                
              !Temp =   VectorLength * GPShapeFunctionBoundary_XI_AllElements(ElementUponWhichLoadIsApplied,I,K)!GPShapeFunctionBoundary(I,K) !TractionNodes(J,1)!GPShapeFunctionBoundary(I, K) * !GPWeightBoundary(I) *
              !SE(K, J) = SE(K, J) + (Temp* (NormalVector(J)/ sqrt(NormalVector(1)**2 +NormalVector(2)**2 +NormalVector(3)**2 ))) !* TractionGaussPoint(I, J)
              
                !elseif (J==3) then 
                !    SE(K, J) = 0 
                !    
                !end if 
                
                    
                    
              !(NodalCoordinates(ND,J,IPatch_Temporary)/sqrt(NodalCoordinates(ND,1,IPatch_Temporary)**2 + NodalCoordinates(ND,2,IPatch_Temporary)**2 ))
              !
              
            end do 
          end do
          
          
        end do

        ! add to global load vector
        do I = 1, ELEMENTBOUNDARYNODES!NNod  ! --> ELEMENTBOUNDARYNODES is the number of nodes in the element (eg. 2*2 if traction on a linear square element)
          ND = NDOF(Traction_Connecting_Local_To_Global_ControlPoints(IEN_TRACTION(I,ElementUponWhichLoadIsApplied)))!NDof(ILoadConLocal(I))
          do J = 1, NVECTOR
            !RLoad(ND + J) = RLoad(ND + J) + SE(I,J)    
            RLoad_DEBUG = RLoad_DEBUG + RLoad(ND + J)
          end do
        end do
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------------------------
        
        
        !Factor = real(ELEMENTBOUNDARYNODES) / real(ELEMENTBOUNDARYGAUSSPOINTS)
        !
        !
        !      NParticles = 0
        !      
        !      
        !      !do IElement_Traction = 1, nel_NURBS_Traction ! traction elements
        !          
        !          ParticleID = 0
        !      LoadedParticleCounter = 0
        !          
        !          GlobalElementForTraction = Traction_Connecting_Local_To_Global_Elements(ElementUponWhichLoadIsApplied)
        !          
        !      NParticles =  NPartEle(GlobalElementForTraction, IPatch_Temporary) ! HARDCODED
        !      
        !      allocate(ParticleStatus(NParticles), stat = IError)
        !      
        !      ParticleStatus = .false.
        !      
        !      call DetermineAdjacentParticles(Traction_Connecting_Local_To_Global_Sides(ElementUponWhichLoadIsApplied), NParticles, ParticleStatus)
        !      
        !      ! side number hardcoded              
        !      !call DetermineAdjacentParticles(1, NParticles, ParticleStatus)
        !
        !      
        !      NLoadedParticles = DetermineNAdjacentParticles(ParticleStatus, NParticles)
        !      
        !      !call MarkBoundaryParticles(GlobalElementForTraction, ParticleStatus, NParticles, IPatch_Temporary)
        !      do IParticle = 1, NPartEle(GlobalElementForTraction, IPatch_Temporary) ! loop all MPs of the element (solid and liquid MPs in case of double-point formulation)
        !          ParticleIndex = GetParticleIndex(IParticle, GlobalElementForTraction, IPatch_Temporary)
        !      
        !          !if ( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then
        !              ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
        !      
        !              ParticleID = ParticleID + 1
        !      
        !              if(ParticleStatus(ParticleID)) then ! this MP is loaded
        !      
        !                  LoadedParticleCounter = LoadedParticleCounter + 1
        !                  
        !                  do INode = 1, ElementBoundaryNodes_Xi
        !                  
        !                      do IDim = 1, NDIM
        !                      
        !                      !LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1  * Factor  ) !(2/3)
        !                      LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( SE(INode,IDim) * -1 / real(ELEMENTBOUNDARYGAUSSPOINTS) ) !( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1  * Factor  ) !(2/3)
        !                      
        !                      ! HARDCODED HALF * 0.5
        !                      ! HARDCODED MINUS ONE FOR DOWNWARDS
        !                      
        !                      end do 
        !                      
        !                      
        !                      
        !                      
        !                      
        !                  end do 
        !                                  
        !                  call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle, 1) !ILoadSystem
        !      
        !                  LoadPerParticle = 0
        !                !LoadPerParticle = DetermineLoadPerParticle(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal, Counters%NodTot, NodeCoord, Loadtype)
        !                !
        !                !!if ( LoadType == 1 ) then ! assign to solid phase / solid material point
        !                !  call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        !      
        !                !else if ( LoadType == 2 ) then ! assign to water phase / liquid material point
        !                !  call IncreaseFExtWater(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        !                !
        !                !else if ( LoadType == 3 ) then ! assign to gas phase
        !                !  call IncreaseFExtGas(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
        !                !
        !                !  if ( NFORMULATION == 2 ) then
        !                !    call GiveError('Erroneous definition of applied external load: For double-point formulation only load on solid or liquid phase is allowed.')
        !                !  end if
        !                !
        !                !end if ! load type
        !              end if ! loaded y/n
        !          !end if ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
        !      
        !        end do ! loop over MPs
        !      deallocate(ParticleStatus, stat = IError)
        !      
              !end do 
        
        
        !---------------------------------------------------------------------------------------------------------------------
        ! TRACTION ON THE MPs
        
        
        
        
        
        
        Factor = real(ELEMENTBOUNDARYNODES) / real(ELEMENTBOUNDARYGAUSSPOINTS)
        
        
              NParticles = 0
              
              
              !do IElement_Traction = 1, nel_NURBS_Traction ! traction elements
                  
                  ParticleID = 0
              LoadedParticleCounter = 0
                  
                  GlobalElementForTraction = Traction_Connecting_Local_To_Global_Elements(ElementUponWhichLoadIsApplied)
                  
              NParticles =  NPartEle(GlobalElementForTraction, IPatch_Temporary) ! HARDCODED
              
              allocate(ParticleStatus(NParticles), stat = IError)
              
              ParticleStatus = .false.
              
              call DetermineAdjacentParticles(Traction_Connecting_Local_To_Global_Sides(ElementUponWhichLoadIsApplied), NParticles, ParticleStatus)
              
              ! side number hardcoded              
              !call DetermineAdjacentParticles(1, NParticles, ParticleStatus)

              
              NLoadedParticles = DetermineNAdjacentParticles(ParticleStatus, NParticles)
              
              !call MarkBoundaryParticles(GlobalElementForTraction, ParticleStatus, NParticles, IPatch_Temporary)
              do IParticle = 1, NPartEle(GlobalElementForTraction, IPatch_Temporary) ! loop all MPs of the element (solid and liquid MPs in case of double-point formulation)
                  ParticleIndex = GetParticleIndex(IParticle, GlobalElementForTraction, IPatch_Temporary)
              
                  !if ( (NFORMULATION==1) .or. (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then
                      ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
              
                      ParticleID = ParticleID + 1
              
                      if(ParticleStatus(ParticleID)) then ! this MP is loaded
              
                          LoadedParticleCounter = LoadedParticleCounter + 1
                          
                          do INode = 1, (ElementBoundaryNodes_Xi*ElementBoundaryNodes_Eta)
                          
                              do IDim = 1, NDIM
                              
                              !LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1  * Factor  ) !(2/3)
                              LoadPerParticle(IDim) = LoadPerParticle(IDim) + ( SE(INode,IDim) * 1 / real(ELEMENTBOUNDARYGAUSSPOINTS) ) !( GPShapeFunctionBoundary_Xi_AllElements(ElementUponWhichLoadIsApplied, LoadedParticleCounter, INode) * SE(INode, IDim) * -1  * Factor  ) !(2/3)
                              
                              ! HARDCODED HALF * 0.5
                              ! HARDCODED MINUS ONE FOR DOWNWARDS
                              
                              end do 
                              
                              
                              
                              
                              
                          end do 
                                          
                          call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle, 1) !ILoadSystem
              
                          LoadPerParticle = 0
                        !LoadPerParticle = DetermineLoadPerParticle(NLoadedParticles, NSurfaceNodes, IConLocal, LoadsLocal, Counters%NodTot, NodeCoord, Loadtype)
                        !
                        !!if ( LoadType == 1 ) then ! assign to solid phase / solid material point
                        !  call IncreaseFExt(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
              
                        !else if ( LoadType == 2 ) then ! assign to water phase / liquid material point
                        !  call IncreaseFExtWater(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
                        !
                        !else if ( LoadType == 3 ) then ! assign to gas phase
                        !  call IncreaseFExtGas(Particles(ParticleIndex), LoadPerParticle,ILoadSystem)
                        !
                        !  if ( NFORMULATION == 2 ) then
                        !    call GiveError('Erroneous definition of applied external load: For double-point formulation only load on solid or liquid phase is allowed.')
                        !  end if
                        !
                        !end if ! load type
                      end if ! loaded y/n
                  !end if ! all MPs in 1-point formulation or only solid MPs in 2-point formulation
              
                end do ! loop over MPs
              deallocate(ParticleStatus, stat = IError)
        
        
        
        
        
        
        
        !-------------------------------------------------------------------------------------------------------
        

    end subroutine Load3D_NURBS_OnMPs
    
    
    
    

        
      end module ModMPMInit
