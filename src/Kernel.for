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


      subroutine Kernel()
      !************************************************************************
      !
      !  Function : Main program file of the Kernel:
      !             1. Kernel initialisation
      !             2. Initalisation of mesh related data
      !             3. Initialisation of material points related data
      !             4. a) Load phase loops 
      !                b) Time step / iteration loops
      !                c) Result writing
      !             5. Kernel shut down
      !
      !     $Revision: 6798 $
      !     $Date: 2018-03-29 14:40:25 +0200 (Thu, 29 Mar 2018) $
      !
      !**********************************************************************

#ifdef __INTEL_COMPILER
      use IFCore ! Intel Fortran library required for checking key event
      use IFPort ! Intel Fortran library required for setting stack size
#endif      
      use ModInitialiseKernel
      use ModInitialiseElementType
      use ModReadCalculationData
      use ModReadGeometryData
      use ModReadMaterialData
      use ModWriteTestData
      use ModMPMData
      use ModMeshInfo
      use ModRotBoundCond
      use ModMPMMeshAdjustment
      use ModReadMPMData
      use ModElementEvaluation
      use ModDYNConvectivePhase
      use ModMPMDynContact
      use ModMPMExcavation
      use ModCounters
      use ModEmptyElements
      use ModLiquid
      use ModWriteVTKOutput
      use ModParticle
      use ModWriteResultData
      use ModTwoLayerFormulation
      use ModDynamicExplicit
      use ModQuasiStaticImplicit
      use ModString
      use ModTiming

      implicit none

      logical(4) :: UserPressedKey
      integer :: IDTimerInitialisation, IDTimerLoadStep, IDTimerWriteResults

      UserPressedKey = .false.

      !********** 1 - kernel initialisation ******************************
      call InitialiseCalculationParameters() ! initialises the calculation paramters (CalParams)
      call InitialiseElementType() ! initialises the element type in global variables
      call OpenTextOutputFiles() ! open TextOutputFiles
      call InitialiseCalculationParameters() ! initialises the calculation paramters (CalParams)
      call ShowDisclaimer() ! shows disclaimer on screen
      call ShowKernelInformation() ! shows kernel information on screen
      call ReadCommandLineParameters() ! read project name CalParams%FileNames%ProjectName
      call DetermineLoadStep() ! determine load step number CalParams%IStep
      call ReadCalculationParameters() ! read CPS-file and assign data into CalParams%...
      call SetVTKPointers() ! set pointers for Output Files

      call startTimer('Main', IDTimerMain)
      call startTimer('Initialisation', IDTimerInitialisation)

      call ReadMaterialParameters() ! read material data from GOM-file and assign data into MatParams(I)%...
      call InitialiseTextOutputFiles() ! create OUT, TST, MLG, CTSSum, RX files for writing test/debug output data
 
#ifdef USE_OPENMP
      call omp_set_num_threads(8) ! temporary solution of hardcoding the maximum number of threads
      call kmp_set_blocktime(0)
      call kmp_set_stacksize(120777216)
#endif      

      !********** 2 - mesh data initialisation ******************************
      call InitialiseShapeFunctions() ! initialise shape functions
      call InitialiseMeshData() ! allocate and assign mesh related arrays by reading GOM file
      call ReadGeometryParameters() ! read geometry data from GOM-file and assign data into GeoParams%...
      call DetermineAdjacencies() ! determine mesh and element properties
      call ReadSHE() ! only if ApplyEmptyElements
      call Initialise3DCylindricalAnalysis() ! only for 3D Cylindrical Analysis
      call InitialiseRotationMatrix() ! allocate (zero) matrices for rotational boundaries only if ApplyRotBoundCond .TRUE. -> 3D edge calculation
      call InitialiseDerivedMeshData() ! initialise Counters%N (number of DoF) and nodal fixities at boundaries
      call InitialiseConvectivePhaseData() ! allocate (zero) nodal array "TemporaryMappingVector" 
      call InitialiseContactData() ! allocate (zero) nodal arrays used in contact algorithm
      call ReadContactData() ! define contact nodes and node normals
      call DetermineContactSurfaceSoilElements() ! Determine elements on the contact surface for Contact Algorithm
      call InitialiseTwoPhaseData() ! allocate (zero) nodal arrays for two phase calculation (liquid and mixture)
      call InitialiseThreePhaseData() ! allocate (zero) additional nodal arrays for three phase calculation (gas)
      call InitialiseLiquidData() ! allocate (zero) free liquid related arrays 
      call InitialiseAbsorbingBoundaryData() ! allocate and assign arrays for absorbing boundaries
      call InitialiseNodalArrays() ! allocate (zero) nodal arrays (load, displacment, velocity, acceleration, momentum, etc.)
      call ReadNodalDataFromFile() ! assign data from previous load step (only if IsFollowUpPhase)
      call ReinitialiseUpdatedNodes() ! for updated mesh, only if IsFollowUpPhase
      call DetermineElementLMin() ! calulate minimum element altitude

      ! ********** 3 - material point data initialisation ******************************
      call InitialiseMaterialPointHousekeeping() ! initialise material points and their housekeeping arrays, fill Particles(ID)%...
      call InitialiseMaterialPointPrescribedVelocity() ! only with Moving Mesh
      call TwoLayerData%Initialise() !For Double Point formulation
      call ResetMaterialPointDisplacements() ! only if .CalParams%ApplyResetDisplacements
      call InitialiseMaterialPointOutputFiles() ! create PAR_XXX files for data output 
      call ComputeInterfaceNodesAdhesion() ! only if ApplyContactAlgorithm: read the normals for contact algorithm
      call InitialiseMeshAdjustment() ! only if ApplyMeshSmoothing: for moving mesh algorithm
      call DetermineDoFMovingMeshStructure() ! only if ApplyMeshSmoothing: for moving mesh algorithm
      call InitialiseTractionLoad() ! if traction load is applied (only if NLoadedElementSides>0)
      call AssignTractionToEntity() ! distribute traction load to entities
      call CalculateNodeElement() ! only if ApplyContactAlgorithm
      call SetUpEntityElements() ! create lists storing which material points and elements related to different entities
      call SetUpMaterialElements() !create lists storing which material points and elements related to different materials
      call InitialiseAbsorbingBoundaryDashpotSpring() ! only if ApplyAbsorbingBoundary
      call MapDataFromNodesToParticles() ! only if ApplyFEMtoMPM: map velocity and displacement to particles
      call InitialiseMaterialPointsForK0Stresses() ! only if ApplyK0Procedure and .not.IsFollowUpPhase
      call InitialiseAbsorbingBoundariesForcesAndStiffness() ! only if ApplyAbsorbingBoundary
      call TwoLayerData%DetermineConcentrationRatios() !For Double Point formulation
      call TwoLayerData%DetermineTwoLayerStatus() ! assign a Liquid or Solid status to the MP
      call InitialiseQuasiStaticImplicit() ! contain calls to subroutine use in Quasi-Static procedure
	  call InitialiseVelocityonMP() ! only if ApplyInitialVelocityonMP
      call InitialiseRigidBody() ! only if IsRigidBody
      call InitialiseSurfaceReaction() !read GOM file and determine surface reactions
      call InitialiseSurfaceReactionOutputFiles() ! create RSurf_XXX files for output of reaction surfaces

      !********** 4a - LOAD PHASE LOOP ******************************
      do while(NotFinishedComputation().and.(.not.CalParams%ConvergenceCheck%DoesDiverge))
        call startTimer('LoadStep', IDTimerLoadStep)
        call GiveMessage('Calculation of load step ' // trim(String(CalParams%IStep)))
        call DATE_AND_TIME(CalParams%StartDate,CalParams%StartTime)
        CalParams%ConvergenceCheck%DoesConverge = .false.
        call UpdateLoadPhaseCounters()
        if (.not.CalParams%ApplyImplicitQuasiStatic) then
          call UpdateMultipliersForLoadStep()
        end if
        call InitialiseLoadPhaseOutputFiles()

        AccumulatedDisplacementSoil = 0.0
        AccumulatedDisplacementWater = 0.0
        call ResetDisplacements() ! only if ApplyResetDisplacements=true  

        if (CalParams%IStep>1) then
          call CalculateCriticalTimeStep()
        end if

        call ApplyExcavation()

        !********** 4b - TIME STEP / ITERATION LOOP ******************************
        if (CalParams%ApplyImplicitQuasiStatic) then ! Iteration loop quasi-static MPM
          call RunImplicitQuasiStaticLoadStep()
        else ! Time step loop dynamic MPM
          call RunExplicitDynamicLoadStep()
        end if

#ifdef __INTEL_COMPILER        
        UserPressedKey = PeekCharQQ()
        if (UserPressedKey) then
          ! Ask whether the user really wants to abort the calculation
        end if
#endif      
        !********** 4c - RESULT WRITING ******************************
         call startTimer('WriteResults', IDTimerWriteResults)
         call WriteLoadPhaseResults() ! only if OutputEachLoadStep
         call WriteGIPFile(DO_USE_STEP_EXTENSION)
         call DATE_AND_TIME(CalParams%EndDate,CalParams%EndTime)
         call WriteLoadStepInformation()
         call WriteBenchmarkResults()
         call WriteBenchmarkResultsStressStrain()
        CalParams%IStep = CalParams%IStep + 1

        call GetStepExt(CalParams%IStep, CalParams%FileNames%LoadStepExt)
        call WriteCalculationParameters() !always output this
        call finishTimer(IDTimerWriteResults)
        call finishTimer(IDTimerLoadStep)

      end do ! ----- end LOAD PHASE LOOP -----

      call finishTimer(IDTimerMain)
      call giveTimerTable()

      ! ********** 5 - shut down kernel ******************************
      call DestroyHouseKeeping()
      call DestroyMeshData()
      call DestroyTwoPhaseData()
      call DestroyThreePhaseData()
      call DestroyLiquidData()
      call DestroyContactData()
      call DestroyExcavationData()
      call DestroyConvectivePhaseData()
      call DestroyModEmptyElementsArrays()
      call DestroyRotBoundCond()
      call DestroyMeshAdjustment()
      call DestroyMaterialParameters()
      call TwoLayerData%Destroy()
      call DestroyQuasiStaticImplicit()
      call GiveMessage('Calculation finished.')
      call CloseTextOutputFiles()
      call DestroyPrescribedNodalVeloData()
      call DestroyPrescribedMPVeloData()

      end subroutine Kernel