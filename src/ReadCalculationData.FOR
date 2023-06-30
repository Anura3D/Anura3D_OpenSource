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
	  
	  
	  module ModReadCalculationData
      !**********************************************************************
      !
      ! Function: Contains routines for reading calculation data from CPS file
      !
      !  Note: Version archived subroutines are contained at the end of the module
      !
      !     $Revision: 10002 $
      !     $Date: 2023-06-19 12:44:04 +0200 (ma, 19 jun 2023) $
      !
      !**********************************************************************
      use ModGlobalConstants
      use ModCounters
      use ModString
      use ModGlobalData
      use ModFeedback
      use ModFileIO

      implicit none

        logical, public :: IsDistorted = .false.
        real(REAL_TYPE), public :: MinimumDeterminantRatioReached = 1.0d0
        
        ! Initiating matrices to read hydraulic head from file
        real(REAL_TYPE), dimension(:,:), allocatable :: HydraulicHeadFileMatrix
        ! Initiating matrices to read prescribed velocity/acceleration from file
        ! These matrices are only allocated if prescribed velocity/acceleration from file feature is used
        real(REAL_TYPE), dimension(:,:), allocatable :: PrescribedVelocityFileMatrix
        
        type LoadMultiplierType ! Load multipliers
        
          ! Note: "Realised" refers to the multiplier at the end of a previous set of load steps (for example, steps 0 to 15)
          !       "Current" refers to the multiplier during a set of load steps (for example, at step 18 coming from step 15)
          !       "Final" refers to the multiplier at the end of a set of load steps (for example, at step 20 coming from step 15)
          !       "Increment" refers to the change of multiplier between two load steps 
          !        (in case of time-dependency, the change of amplitude)
          !       "TimeDependency" refers to the change of the multiplier inside a load step as a function of time
          real(REAL_TYPE), dimension(MAX_LOAD_SYSTEMS) :: SolidARealised, & ! Realised multiplier external load A due to soil traction
                             SolidAPrevious, & ! Previous load step multiplier external load A due to soil traction
                             SolidACurrent, & ! Current multiplier external load A due to soil traction
                             SolidALoadStepIncrement, & ! Multiplier increment from previous load step multiplier for the current load step
                             SolidAFinal, & ! Final multiplier external load A due to soil traction
                             SolidAIncrement, & ! Multiplier increment from released multiplier external load A due to soil traction
                              WaterARealised, & ! Realised multiplier external water pressure A
                              WaterAPrevious, & ! Previous load step multiplier
                              WaterACurrent, & ! Current multiplier external water pressure A
                              WaterALoadStepIncrement, & ! Multiplier increment from previous load step multiplier for the current load step
                              WaterAFinal, & ! Final multiplier external water pressure A
                              WaterAIncrement, & ! Multiplier increment external water pressure A
                              WaterATimeDependency, & ! Multiplier for time dependency external water pressure A
                              WaterBRealised, & ! Realised multiplier external water pressure B (MP)
                              WaterBPrevious, & ! Previous load step multiplier
                              WaterBCurrent, & ! Current multiplier external water pressure B
                              WaterBLoadStepIncrement, & ! Multiplier increment from previous load step multiplier for the current load step
                              WaterBFinal, & ! Final multiplier external water pressure B
                              WaterBIncrement, & ! Multiplier increment external water pressure B
                              WaterBTimeDependency, & ! Multiplier for time dependency external water pressure B
                              GasARealised, & ! Realised multiplier external gas pressure A
                              GasAPrevious, & ! Previous load step multiplier
                              GasAFinal, & ! Final multiplier external gas pressure A
                              GasAIncrement, & ! Multiplier increment external gas pressure A
                              GasATimeDependency, & ! Multiplier for time dependency external gas pressure A
                              GasACurrent, & ! Current multiplier external gas pressure A
                              GasALoadStepIncrement
          
          real(REAL_TYPE) ::  HydraulicHeadPrevious, & ! Previous load step multiplier
                              HydraulicHeadRealised, & ! Realised multiplier external HydraulicHead (NODE)
                              HydraulicHeadCurrent, & ! Current multiplier external water pressure A
                              !HydraulicHeadFinal, & ! Final multiplier external HydraulicHead
                              GravityRealised, & ! Realised multiplier gravity
                              GravityPrevious, & ! Previous load step multiplier
                              GravityCurrent, & ! Current multiplier gravity
                              GravityLoadStepIncrement, & ! Multiplier increment from previous load step multiplier for the current load step
                              GravityFinal, & ! Final multiplier gravity
                              GravityIncrement, & ! Multiplier increment gravity
                              GravityTimeDependent, &
                              VelocitySolidRealised, &
                              VelocitySolidPrevious, & ! Previous load step multiplier
                              VelocitySolidCurrent, & ! Current multiplier 
                              VelocitySolidLoadStepIncrement, & ! Multiplier increment from previous load step multiplier for the current load step
                              VelocitySolidFinal, & ! Final multiplier 
                              VelocitySolidIncrement,&
                              AccelerationSolid, & !used to derive acceleration from prescribed velocity variable in time
                              WaterASpaceDependency !used to build the vector of water external load in ModMPMDYNConsolidation

          character(len = 6), dimension(MAX_LOAD_SYSTEMS)  :: SolidALoadType, & ! LOAD_TYPE_OFF, LOAD_TYPE_STEP, LOAD_TYPE_LINEAR
                                WaterALoadType, &
                                WaterBLoadType, &
                                GasALoadType
                                
     
          character(len = 6) :: GravityLoadType, & ! LOAD_TYPE_OFF, LOAD_TYPE_STEP, LOAD_TYPE_LINEAR
                                VelocitySolidLoadType, &
                                HydraulicHeadType
        end type LoadMultiplierType
        
        type MovingMeshType(dsize)
          integer(INTEGER_TYPE), LEN :: dsize  
          integer(INTEGER_TYPE), dimension(2, 3, dsize * (dsize - 1) + 2) :: MeshAreas ! Node ID's of corner nodes ! defining the two particle storage and fixed mesh areas
          integer(INTEGER_TYPE) :: NAreaNodes ! Number of nodes that define an area (6 for triangular shape or 8 for cube)
          integer(INTEGER_TYPE) :: NMovingMeshDirections ! Number of directions the moving mesh can move (i.e., = 1 means that it can only move parallel to X, Y or Z axis) 
          integer(INTEGER_TYPE) :: MovingMeshDirection ! Default direction of movement
          integer(INTEGER_TYPE), dimension(2) :: NStorageAreas ! Number of storage areas
          integer(INTEGER_TYPE) :: StructureMaterialID ! Material ID of the structure material (-1 if not defined)
          integer(INTEGER_TYPE) :: MovingMaterialID ! Material ID to follow with the moving mesh
          integer(INTEGER_TYPE) :: NStructureNodes ! Number of nodes that define the structure
          integer(INTEGER_TYPE), dimension(:), allocatable :: StructureNodes ! Set of nodes defining the structure
        
        end type MovingMeshType
        
        type PrescribedVeloType
          logical :: ApplyPrescribedVelo = .false.
          integer(INTEGER_TYPE) :: NNodePrescribedVelo = 0 ! Number of nodes on which the velocity is prescribed
          integer(INTEGER_TYPE), dimension(:), allocatable :: NodePrescribedVelo ! Set of nodes on which the velocity is prescribed
          real(REAL_TYPE), dimension(:, :), allocatable :: NodalPrescribedVelocityValue ! Prescribed velocity at the nodes
          real(REAL_TYPE), dimension(:, :), allocatable :: NodalPrescribedVelocityDirection ! 0/1 prescribed direction
          real(REAL_TYPE) :: FileMultiplier = 0
          character(len=255) :: FileNamePrescribedVelocity = ''
          integer(INTEGER_TYPE) :: PrescribedVelocityFileNumberOfLines = 0 
          character(len=255) :: FilePrescribedVelocitySwitch = 'off' ! on/off file switch 
          real(REAL_TYPE), allocatable, dimension(:) :: TransientBoundaryConditionStaringTime ! Start time of the transient boundary condition
        
        
        end type PrescribedVeloType
        
         type HydraulicHeadType
          logical :: HydraulicHead  = .false.
          !real(REAL_TYPE) :: FileMultiplier = 0
          character(len=255) :: FileNameHydraulicHead = 'X'
          integer(INTEGER_TYPE) :: HydraulicHeadFileNumberOfLines = 0
          
        end type HydraulicHeadType
            
        type ParabolicTrajectoryType(vsize)
          integer(INTEGER_TYPE), len :: vsize 
          real(REAL_TYPE) :: InitialAngle
          real(REAL_TYPE) :: DeltaY
          real(REAL_TYPE), dimension(vsize) :: RotationPoint
          real(REAL_TYPE) :: InitialVelocity
          logical :: ApplyRotation 
          real(REAL_TYPE), dimension(vsize) :: TrajectoryPrescribedVelocity
          logical :: Restart 

        end type ParabolicTrajectoryType

        type FileNameType ! File names
          
          character(len = MAX_FILENAME_LENGTH) :: ProjectName ! Name of project file
          character(len = 80) :: ResultFileHeader ! Header of the result file
          character(len = MAX_EXTENSION_LENGTH) :: PreviousStepExt ! File extension for previous step in case of restart
          character(len = MAX_EXTENSION_LENGTH) :: LoadStepExt ! Extension of current load phase result file
          character(len = MAX_EXTENSION_LENGTH) :: TimeStepExt ! Extension of current time step result file
        
        end type FileNameType 
        
        type ConvergenceCheckType ! Check convergence
        
          logical :: DoesConverge, & ! True - end of time stepping
                      ApplyDivergenceCheck, & ! Flag for divergence check (default = true)
                     DoesDiverge ! True - end of calculation
          real(REAL_TYPE) :: KineticEnergy, & ! System kinetic energy
                              KineticEnergy0, & ! System initial kinetic energy     
                              KineticError, & ! Tolerance error kinetic energy
                              ExternalWork, & ! System external work
                              InternalWork, & ! System internal work
                              EnergyDissipation, & ! System (unphysical) energy dissipation
                              ForceError, & ! Out-of-balance forces error
                              KineticEnergySoil, & ! System kinetic energy solid phase
                              KineticEnergySoil0, & ! System initial kinetic energy solid phase 
                              KineticErrorSoil, & ! Tolerance error kinetic energy solid phase
                              ExternalWorkSoil, & ! System external work solid phase
                              InternalWorkSoil, & ! System internal work solid phase
                              EnergyDissipationSoil, & ! System (unphysical) energy dissipation solid phase
                              ForceErrorSoil, & ! Out-of-balance forces error solid phase
                              KineticEnergyWater, & ! System kinetic energy water phase
                              KineticEnergyWater0, & ! System initial kinetic energy water phase 
                              KineticErrorWater, & ! Tolerance error kinetic energy water phase
                              ExternalWorkWater, & ! System external work water phase
                              InternalWorkWater, & ! System internal work water phase
                              EnergyDissipationWater, & ! System (unphysical) energy dissipation water phase
                              ForceErrorWater, & ! Out-of-balance forces error water phase
                              SumLocalError, &
                              SumIntegrationPointWeights
          integer(INTEGER_TYPE) :: NInaccuratePlasticPoints
        
        end type ConvergenceCheckType

        type GravityDataType(dsize) ! Gravity data
          integer(INTEGER_TYPE), LEN :: dsize    
          real(REAL_TYPE) :: GAccel, & ! Gravitational accleration
                             GravityVector(dsize) ! Gravity by default in negative y direction
        
        end type GravityDataType

        type IntegrationPointDataType ! Integration point data
        
          integer(INTEGER_TYPE) :: NPlasticPoints, &
                     NNegativePlasticPoints, &
                     NApexPoints, &
                     NTensionCutOffPoints
        
        end type IntegrationPointDataType

        type AbsorbingBoundariesType ! Absorbing Boundaries
        
          real(REAL_TYPE) :: ANormalSolid, & ! (an) coefficient for the normal dashpot in the solid phase
                              AShearSolid, & ! (as) coefficient for the shear dashpot in the solid phase 
                              AWater, & ! (aw) coefficient for the dashpot in the water phase
                              AGas, & ! (ag) coefficient for the dashpot in the gas phase
                              DeltaNormalSolid, & ! (Deltan) the virtual thickness of the normal spring in the solid phase
                              DeltaShearSolid, & ! (Deltas) the virtual thickness of the shear spring in the solid phase
                              DeltaWater, & ! (Deltaw) the virtual thickness of the spring in the water phase 
                              DeltaGas, & ! (DeltaG) the virtual thickness of the spring in the gas phase 
                              VBMaterialSet ! The set of the boundary material --> should be integer?
          logical :: ApplyVariableStiffAB, &
                     UpdateBoundaryStiff, &
                     ApplyNodeData
     
        end type AbsorbingBoundariesType

        type RigidBodyType(vsize, dsize) ! Rigid Body

          integer(INTEGER_TYPE), LEN :: vsize, dsize
          logical :: IsRigidBody ! True - consider rigid body

          real(REAL_TYPE) :: InternalForce(vsize), & ! internal forces assigned to the representative point of the rigid body
                              Mass, & ! mass assigned to the representative point of the rigid body
                              TractionForce(dsize), & ! traction force assigned to the representative point of the rigid body
                              GravityForce(dsize), & ! gravity force assigned to the representative point of the rigid body
                              Acceleration(dsize), & ! acceleration assigned to the representative point of the rigid body
                              Velocity(dsize), & ! velocity assigned to the representative point of the rigid body
                              CurrentTraction(dsize), & ! traction after the time multiplier
                              CurrentGravity(dsize), &
                             InAcceleration(dsize)
          integer(INTEGER_TYPE) :: RigidEntity, & ! the entity that should be treated as a rigid body
								   Constrains(vsize) !Constrained displacements

        end type RigidBodyType

        type ImplicitIntegrationType
          
          integer(INTEGER_TYPE) :: Iteration
          integer(INTEGER_TYPE) :: MinDesiredIterations 
          integer(INTEGER_TYPE) :: MaxDesiredIterations
          real(REAL_TYPE) :: UpScaleFactor
          real(REAL_TYPE) :: DownScaleFactor
          integer(INTEGER_TYPE) :: MaxIterations
          logical :: DoUseZLS
          logical :: IsZeroLoadStep
          logical :: DoUseArcLengthControl
          logical :: DoUseEnhancedStiffnessMatrix
          real(REAL_TYPE) :: StiffnessIncreaseFactor
          real(REAL_TYPE) :: DanglingElementFactor
          logical :: DoUseStepExtrapolation
          logical :: DoCheckChangeJacobian
          real(REAL_TYPE) :: InitialMultiplier
          logical :: DoUseAutomaticLoadStepping
          real(REAL_TYPE) :: ScaleFactor
          real(REAL_TYPE) :: MaximumJacobianRatio
          real(REAL_TYPE) :: LimitJacobianChange
          
        end type ImplicitIntegrationType

        
        type BoundaryConditionsType !(vsize,dsize)
            !integer(INTEGER_TYPE), LEN :: vsize, dsize
            
            !real(REAL_TYPE):: !InfiltrationArea(2,dsize),& !Boundary of infiltration area (min:max,dim)
                              !InfiltrationRate(vsize),& !Value of infiltration rate [m/s]
                              !SeepageArea(2,dsize),& !Boundary of free seepage area (min:max,dim)
                              !HydrHeadArea(2,dsize) !Boundary of Hydraulic Head loaded elements area (min:max,dim)
            logical :: ApplyInfiltrationRate, ApplySeepageFace, UseSeepageFace, UseInfiltration
            
        end type BoundaryConditionsType  

        type CalParamType(vsize, dsize)
        
          integer(INTEGER_TYPE), LEN :: vsize, dsize
          type(LoadMultiplierType) :: Multipliers ! Multipliers
          type(FileNameType) :: FileNames ! File names
          type(ConvergenceCheckType) :: ConvergenceCheck ! Convergence check
          type(GravityDataType(:)), allocatable :: GravityData ! Gravity data
          type(IntegrationPointDataType) :: IntegrationPointData ! Integration point data
          type(AbsorbingBoundariesType) :: AbsorbingBoundaries ! Absorbing boundaries
          type(RigidBodyType(:,:)), allocatable :: RigidBody ! rigid body
          type(MovingMeshType(:)), allocatable :: MovingMesh ! Moving mesh
          type(ParabolicTrajectoryType(:)), allocatable :: ParabolicTrajectory
          type(ImplicitIntegrationType) :: ImplicitIntegration
          type(PrescribedVeloType):: PrescribedVelo
          type(HydraulicHeadType):: PrescribedHead
          type(BoundaryConditionsType) :: BoundaryConditions !(:,:)), allocatable :: BoundaryConditions
     
          real(REAL_TYPE) :: TimeIncrement, & ! The size of time increment
                              TotalTime, & ! The total time of simulation
                              OverallRealTime, &
                              DampingFactor, & ! damping Factor
                              ScalingMassFactor, & ! scaling mass Factor
                              ToleratedErrorEnergy, & ! Tolerated error for energy convergence criterion for the mixture
                              ToleratedErrorForce, & ! Tolerated error for out off balance force for the mixture
                              ToleratedErrorEnergyWater, & ! Tolerated error for energy convergence criterion for the water phase
                              ToleratedErrorForceWater, & ! Tolerated error for out off balance force for the water phase
                              ToleratedDivergence, & ! Tolerated divergence for divergence criterion (in J)
                              TotalRealTime, & ! the real time
                              FricCoef, & ! friction coefficient used on contact model
                              Adhesion, & ! adhesion used in contact formulation
                              RequiredDegreeOfFilling, & ! filling factor deciding whether to use Gauss Point or Material Point integration:
                                                         ! 0.0 - always Gauss Point integration, 
                                                         ! 1.0 - always Material Point integration
                              FacStiffnessIncrease, & ! Factor for setting the percentage of the summed diagonal elastic element
                                                      ! stiffness terms added to the element stiffness matrix
                              SectorAngle, & ! Regarding pile driving
                              XMin, & ! Regarding pile driving
                              XMax, & ! Regarding pile driving
                              YMin, & ! Regarding pile driving
                              ZMin, & ! Regarding pile driving
                              EMax, & ! maximum oedometer stiffness
                              CourantNumber, & ! Fraction of critical time step
                              Fix(2), & ! Fraction of critical time step
                              SoilSurfacePoint(2), & ! used for K0-procedure, an arbitrary point on the soil surface
                              LiquidSurfacePoint(2), & ! used for K0-procedure, an arbitrary point on the liquid surface
                              WaveAlpha, & ! used for K0-procedure (non-horizontal surface, sine-wave)
                              WaveBeta, & ! used for K0-procedure (non-horizontal surface, sine-wave)
                              InitialVerticalLoadK0, & ! External load applied on a horizontal surface considered with K0 procedure
                              K0MaxSurfaceSuction, & !max suction at soil surface with K0 and 2phase+suction formulation
                              PrescribedVelocity(dsize) , &
                              SurfacePrescribedVelocity(vsize), &
                              SurfacePrescribedVelocityCoordMin(vsize), &
                              SurfacePrescribedVelocityCoordMax(vsize), &
                              LiquidPressureCavitationThreshold, & ! Liquid pressure threshold for cavitation (positive value)
                              FreeSurfaceFactor, & ! Factor used for detection of free surface (should be positive and < 1)
                              VirtualParticleMassFactor, &
                              OutputRealTimeID(MAXOUTPUTREALTIME), & ! array of real time values for additional output-files 
                              RigidCoordinateMin(vsize), & ! minimum coordinates of area in which the rigid loading surface is
                              RigidCoordinateMax(vsize), & ! maximum coordinates of area in which the rigid loading surface is
                              LimitPorosity, & ! 2-layer formulation
                              LiquidThresholdDensity, & ! 2-layer formulation
                              BinghamYieldStress, & ! Bingham fluid, yield shear stress
                              BinghamViscosity, & ! Bingham fluid, viscosity
                              BinghamYoungModulus, & ! Bingham fluid, elastic young modulus
                              BinghamPoissonRatio, & ! Bingham fluid, elastic poisson ratio
                              FrictionalFluidFrictionAngle, & ! friction angle of frictional fluid (to define yield stress)
                              TwoLayerErgunLawDiameterPartic, & ! Average soil particle diameter for Ergun Law
                              TwoLayerErgunLawDiameterPartic2, & ! 2nd Average soil partilce diameter for Ergu Law 
                              IntrinsicPermeability, & ! Intrinsic permeability for the double point formulation
                              ERGUNCONSTANTA, & ! Constant A for Ergun Law
                              ERGUNCONSTANTB, & ! Constant B for Ergun Law
                              ThicknessSoilLayer(MAXNUMBERSOILLAYERS), &
							  BoreHoleData(MAXNUMBERSOILLAYERS,8), &
                              InitialWaterPressure, &
                              BulkViscosityDamping1, & ! bulk viscosity damping factors beta1 and beta2
                              BulkViscosityDamping2, &
                              MinimumDeterminantRatio, &
                              ShrinkageMateriaPointPositionFactor, &
                              GammaWater, &
							  ContactScalingFactor, &
                              ManuallyDefinedVectors(MAXNUMBERSOILLAYERS,1+vsize) !manually modified normal vectors


          integer(INTEGER_TYPE) :: OutputNumberParticles, & ! number of material points for which output data is exported
                     OutputNumberElements, & ! number of elements for which output data is exported
                     OutputParticles(MAXOUTPUTPARTICLES), & ! array of material point IDs for which output data is exported
                     OutputElements(MAX_OUTPUT_ELEMENTS), & ! array of elements IDs for which output data is exported
                     OutputNumberNodes, & ! number of nodes for which output data is exported
                     OutputNodes(MAXOUTPUTNODES), & ! array of node IDs for which output data is exported
                     NReactionHistoryNodes, & ! number of nodes for which reaction forces are written to the NRF file
                     ReactionHistoryNodes(MAXREACTIONHISTORYNODES), & ! array of node IDs for which reaction forces are written
                     OutputNumberOfSubsteps, & ! number of substeps for which additional output-files are generated
                     OutputNumberOfTimeSteps, & ! number of timesteps for which additional output-files are generated
                     OutputTimeStepID(MAXOUTPUTTIMESTEPS), & ! array of time step numbers for additional output-files 
                     OutputNumberOfRealTime, & ! number of timesteps for which additional output-files are generated
                     OutputNumberOfLoadSteps, &
                     TimeStep, & ! Current time step
                     MaxTimeSteps, & ! The maximum time step per load step in quasi-static calculation
                     NThreads, & ! Number of parallelized threads
                     NMaterialPoints, & ! initial number of all (solid+liquid) material points per element (default = 4)     
                     NLoadSteps, & ! Initial number of load phases
                     IStep, & ! Number of the current load phase
                     PreviouslyRealisedLoadStep, & ! Number of last previously computed step
                     PKernel32, & ! ...
                     CalculationType, & ! 1 - quasi-static, 2 - dynamic     
                     RefParticleID, &
                     MaxFileLength, & ! Maximum number of rows for an output text file
                     RowCounter, & ! Current number of rows written to the output file for CPT data
                     FileCounter, & ! Current number of the file for storing CPT data
                     ParticleFileCounter, & ! Current number of the particle file for storing particle data
                     ElementFileCounter, & ! Current number of the element file for storing element data
                     GroupThreshold, & ! Limit number of elements below which a new group is defined
                     RecursionThreshold, & ! Number of elements that are checked for belonging to a group
                     NVirtualParticles, & ! Number of virtual particles added for each empty element
                     TimeStepOutputInterval, &
                     ComputationMethod, &
                     TrajectoryType, &
                     NumberOfPhases, & ! Number of phases considered in the calculation, 1: 1-phase (solid), 
                                       !                                                 2: 2-phases (solid + liquid), 
                                       !                                                 3: 3-phase (solid + liquid + gas) 
                     !NumberOfLayers, & ! Number of layers of material points, 1: one layer of material points carrying solid and liquid and gas together, 
                                       !                                      2: one layer of material points carrying solid + one layer of material points carrying liquid + gas, 
                                       !                                      3: not implemented 
                     RigidDirection(vsize), & ! direction in which a fixity is applied in case of rigid loading
                     NumberOfMaterials, & ! number of materials
                     FirstSolidMaterialIndex, & ! Two layer formulation- the material index of the material firstly assignd with a diameter in CPS
                     SecondSolidMaterialIndex, & ! Two layer formulation- the material index of the material secondly assignd with a diameter in CPS
                     NumberSubmergedCalculation, & ! number of load steps for which a submerged calculation is performed
                     NumberExcavatedVolumes, & ! Number of volumes to excavate
                     ExcavatedVolumeID(20), & ! ID of excavated volume (an integer is associated to each volume)
                     ExcavationStages(20,2), & ! List of Loadsteps in which each volume must be excavated 
                     NumberSoilLayers, & ! Number of soil layers. Needed for K0 procedure
                     OutputCurvesIntervals, & ! number of time steps for which curve results are written in output test files
                     NumberSkipConvection, &
                     DumpElement, &
                     FeedbackLevel, &
					 NCorrectedNormals !Number of normals to be modified in the contact

          logical :: ApplyContactMeshBoundary, &
                     ApplyQuasiStatic, & ! True - apply quasti-static convergence criteria
                     TwoLayerApplyUpdatePermeability, & ! True - update permeability according to porosity,  False - use initial intrinsic permeability
                     TwoLayerApplyErgunLaw, & !True - use Ergun Flow law,  False - use Darcy law
                     TwoLayerApplyNoTensStressLiqMPwLiqStatus, & ! True - set mean stress to zero in LiquidMP with Liquid Phase Status,  False - allow tensile stress
                     ApplyStrainSmoothing, & ! True - apply strain smoothing (usually used for solids), False - off
                     ApplyStrainSmoothingLiquidTwoLayer, & ! True - apply strain smoothing , False - off
                     ApplyDetectLiquidFreeSurface, & ! True - apply detection of free surface, False - off
                     ApplyContactAlgorithm, & ! True - use the dynamic contact model !CC 
                     ApplyTractionContact, & ! True - use the normal traction to judge if contact occures, otherwise based on normal velocities
                     AssembleSoilLoadVectorFromMP, & ! Assemble external soil load vector from material point loads, else take initial nodal load vector
                     AssembleWaterLoadVectorFromMP, & ! Assemble external water load vector from material point loads, else take initial nodal load vector
                     AssembleGasLoadVectorFromMP, & ! Assemble external gas load vector from material point loads, else take initial nodal load vector
                     ApplyPrescribedDisplacements, & ! TRUE - PPD file exists, else (by default) FALSE
                     ApplySurfacePrescribedVelocity(vsize), & ! TRUE - consider prescribed velocity applied on surface
                     ApplyMeshSmoothing, & ! TRUE - if moving mesh is used
                     ApplyUpdateWeights, & ! Updating of particle weights: FALSE - off, TRUE - on
                     ApplyAbsorbingBoundary, & ! True - consider the absorbing boundary: False-apply rigid boundary
                     ApplyCPSDamping, & ! if true the damping as defined in CPS file is used for all elements; if false the damping from GOM file is used
                     ApplyFEMtoMPM, & ! if true  the state of the solution is updated from nodes to particles from a previous FEM calculation
                     ApplyResetDisplacements, & ! if ture then the displacement will be set to zero
                     ApplySubmergedCalculation, & ! if ture then the material weight will use submerged weight; if false, use saturated weight
                     ApplyWaveLoad, &  ! if ture use inpulseload to generate the wave load
                     ApplyMaterialUpdate, &  ! if ture then read material parameters from GOM
                     ApplyObjectiveStress, & ! if true then objective stresses are included
                     ApplyRemoveSolidFixities, & ! id true if some boundary conditions of solid have to be read from GOMfile and removed
                     ApplyRemoveLiquidFixities, & ! id true if some boundary conditions of liquid have to be read from GOMfile and removed
                     ApplyRemoveGasFixities, & ! id true if some boundary conditions of gas have to be read from GOMfile and removed
                     ApplyK0Procedure, &   ! Switch for applying K0    
                     ApplyFixParticlesK0, &   ! Switch for applying fixing particles in K0-procedure
                     ApplyEmptyElements, & ! Fill holes inside a soil body with virtual particles if .true.
                     ApplyEffectiveStressAnalysis, & ! Switch to apply effective stress analysis 
                     OutputEachLoadStep = .true., & ! TRUE - write standard output each load step: FALSE - write only the final load step.
                     OutputDebugData, &  ! if true then debug data is written to the corresponding output files
                     ApplyQuickCheckOutput = .true., & ! Write results for first time step to file for checking
                     ApplyMassScaling, &     
                     ApplyFixMeshBottom, &
                     ApplyDampingForSoilStructure, &
                     ApplyPrescribedVelocity(vsize), &
                     CorrectNormalsForSlice, &				     
                     ApplyStructureTrajectory, &
                     ApplyPorosityUpdate, &
                     ApplyRigidLoading, &  ! TRUE - apply a rigid load, i.e. fixities in direction which is perpendicular to loading direction
                     ApplyBinghamFluid, & ! treat the material as a Bingham fluid 
                     ApplyFrictionalFluid, & ! like the Bingham fluid but with stress dipendent yield stress     
                     ApplyImplicitQuasiStatic, &
                     ApplyExcavation, & ! TRUE - remove material points inside geometry
                     ApplyCloseBoundaryforDike, &
                     ExternalSoilModel, & ! TRUE - soil model implemented in the NGGS framework is used 
                     ApplyFixedSolidSkeleton, & ! TRUE - solid momentum = 0
                     SkipConvection, &
                     AutomaticSkipConvection, &
                     ConsiderTimingProcedures, &
                     ApplyBulkViscosityDamping, & ! TRUE - apply bulk viscosity damping
                     OutputBasicData, &
                     SkipWaterBC, &
                     ConsiderPTInLoadK0, &
                     DoMapStresses, &
                     MonoPileSliceContact, &
                     IsVTKBinary, &
                     ApplyOneLayerCavitationThreshold, &
                     UseUniformMaterialPointWeight, &
					 ApplyInitialVelocityonMP = .false.,& !True if initial velocity is set in GOM file, false otherwise
                     ApplyPartialSaturation = .false.,& !True if partial saturation is taken into account in 2-phase single-point formulation
                     ApplySmootheningLiquidPressureIncrement = .false.,& !True if liquid pressure incement is smoothened, useful in unsaturated formulation
			         ApplyContactVelocityScaling= .false.,& !True for scaling wrong velocities at the contact
					 ApplyContactNormalCorrection= .false. !Allows to modify some normal vectors in the contact algorithm


          character(len=64) :: CPSversion = 'not defined'
          character(len=64) :: GOMversion = 'not defined'
          character(len=8) :: StartDate = 'N/A' ! starting date of calculation of current load step
          character(len=8) :: EndDate = 'N/A' ! ending date of calculation of current load step
          character(len=10) :: StartTime = 'N/A' ! starting time of calculation of current load step
          character(len=10) :: EndTime = 'N/A' ! ending time of calculation of current load step
          character(len=24) :: Visualization = 'N/A' ! type of visualization
     
        end type CalParamType

        type(CalParamType(:, :)), allocatable, public, save :: CalParams ! Stores variables that manipulate the Dynamic MPM calculation process
        
        character(len=255), dimension(999), public, save :: ActiveParams = 'not defined' ! stores the labels of calculation parameters that are active (dimension = number of lines in CPS file, len = number of columns in CPS file)
           
    contains ! Routines of this module

        
       subroutine InitialiseCalculationParameters()
       !**********************************************************************
       !
       ! Function: Initialises calculation parameters stored in CalParams%
       !
       !**********************************************************************
       
       implicit none
       
       ! allocate main CalParams structure
       allocate(CalParamType(NVECTOR, NDIM)::CalParams)
       
       ! Set Default Parameters
       
        ! LoadMultiplierType
        CalParams%Multipliers%SolidARealised = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%SolidACurrent = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%SolidAFinal = DEFAULT_FINAL_MULTIPLIER
        CalParams%Multipliers%SolidAIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%SolidALoadStepIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%SolidAPrevious = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterARealised = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterACurrent = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterAFinal = DEFAULT_FINAL_MULTIPLIER
        CalParams%Multipliers%WaterAIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterATimeDependency = DEFAULT_TIME_MULTIPLIER
        CalParams%Multipliers%WaterALoadStepIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterAPrevious = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterBRealised = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterBCurrent = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterBFinal = DEFAULT_FINAL_MULTIPLIER
        CalParams%Multipliers%WaterBIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterBTimeDependency = DEFAULT_TIME_MULTIPLIER
        CalParams%Multipliers%WaterBLoadStepIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%WaterBPrevious = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%HydraulicHeadRealised = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%HydraulicHeadCurrent = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%HydraulicHeadPrevious = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GasARealised = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GasAFinal = DEFAULT_FINAL_MULTIPLIER
        CalParams%Multipliers%GasAIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GasATimeDependency = DEFAULT_TIME_MULTIPLIER
        CalParams%Multipliers%GasACurrent = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GasALoadStepIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GasAPrevious = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GravityRealised = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GravityCurrent = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GravityFinal = DEFAULT_FINAL_MULTIPLIER
        CalParams%Multipliers%GravityIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GravityLoadStepIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%GravityTimeDependent = DEFAULT_FINAL_MULTIPLIER
        CalParams%Multipliers%GravityPrevious = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%VelocitySolidRealised = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%VelocitySolidCurrent = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%VelocitySolidFinal = DEFAULT_FINAL_MULTIPLIER
        CalParams%Multipliers%VelocitySolidIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%VelocitySolidPrevious = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%VelocitySolidLoadStepIncrement = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%AccelerationSolid = DEFAULT_INITIAL_MULTIPLIER
        CalParams%Multipliers%SolidALoadType = LOAD_TYPE_STEP
        CalParams%Multipliers%WaterALoadType = LOAD_TYPE_STEP
        CalParams%Multipliers%WaterBLoadType = LOAD_TYPE_STEP
        CalParams%Multipliers%HydraulicHeadType = LOAD_TYPE_OFF
        CalParams%Multipliers%GasALoadType = LOAD_TYPE_STEP
        CalParams%Multipliers%GravityLoadType = LOAD_TYPE_STEP
        CalParams%Multipliers%VelocitySolidLoadType = LOAD_TYPE_STEP
        CalParams%Multipliers%WaterASpaceDependency = 1.0
        
        ! MovingMeshType
        ALLOCATE(MovingMeshType(NDIM)::CalParams%MovingMesh)
        CalParams%MovingMesh%MeshAreas = 1
        CalParams%MovingMesh%NAreaNodes = NDIM * (NDIM - 1) + 2
        CalParams%MovingMesh%NMovingMeshDirections = 1
        CalParams%MovingMesh%MovingMeshDirection = -1
        CalParams%MovingMesh%NStorageAreas = 2
        CalParams%MovingMesh%StructureMaterialID = -1
        CalParams%MovingMesh%MovingMaterialID = -1
        CalParams%MovingMesh%NStructureNodes = 0
        
        ! ParabolicTrajectoryType
        ALLOCATE(ParabolicTrajectoryType(NVECTOR)::CalParams%ParabolicTrajectory)
        CalParams%ParabolicTrajectory%ApplyRotation = .false.
        CalParams%ParabolicTrajectory%TrajectoryPrescribedVelocity = 0.0
        CalParams%ParabolicTrajectory%Restart = .false.
        
              
        ! FileNameType
        CalParams%FileNames%ProjectName = ' '
        CalParams%FileNames%ResultFileHeader = ' '
        CalParams%FileNames%PreviousStepExt = ''
        CalParams%FileNames%LoadStepExt = ''
        CalParams%FileNames%TimeStepExt = ''
        
        ! ConvergenceCheckType
        CalParams%ConvergenceCheck%DoesConverge = .false.
        CalParams%ConvergenceCheck%ApplyDivergenceCheck = .false.
        CalParams%ConvergenceCheck%DoesDiverge = .false.
        CalParams%ConvergenceCheck%KineticEnergy = 0.0
        CalParams%ConvergenceCheck%KineticEnergy0 = -1.0
        CalParams%ConvergenceCheck%KineticError = 0.0
        CalParams%ConvergenceCheck%ExternalWork = 0.0
        CalParams%ConvergenceCheck%InternalWork = 0.0
        CalParams%ConvergenceCheck%EnergyDissipation = 0.0
        CalParams%ConvergenceCheck%ForceError = 0.0
        CalParams%ConvergenceCheck%KineticEnergySoil = 0.0
        CalParams%ConvergenceCheck%KineticEnergySoil0 = -1.0
        CalParams%ConvergenceCheck%KineticErrorSoil = 0.0
        CalParams%ConvergenceCheck%ExternalWorkSoil = 0.0
        CalParams%ConvergenceCheck%InternalWorkSoil = 0.0
        CalParams%ConvergenceCheck%EnergyDissipationSoil = 0.0
        CalParams%ConvergenceCheck%ForceErrorSoil = 0.0
        CalParams%ConvergenceCheck%KineticEnergyWater = 0.0
        CalParams%ConvergenceCheck%KineticEnergyWater0 = -1.0
        CalParams%ConvergenceCheck%KineticErrorWater = 0.0
        CalParams%ConvergenceCheck%ExternalWorkWater = 0.0
        CalParams%ConvergenceCheck%InternalWorkWater = 0.0
        CalParams%ConvergenceCheck%EnergyDissipationWater = 0.0
        CalParams%ConvergenceCheck%ForceErrorWater = 0.0
        CalParams%ConvergenceCheck%SumLocalError = 0.0
        CalParams%ConvergenceCheck%SumIntegrationPointWeights = 0.0
        CalParams%ConvergenceCheck%NInaccuratePlasticPoints = 0
        
        ! GravityDataType
        ALLOCATE(GravityDataType(NDIM)::CalParams%GravityData)
        CalParams%GravityData%GAccel = DEFAULT_GRAVITY_ACCELERATION
        CalParams%GravityData%GravityVector = DEFAULT_GRAVITY_DIRECTION

        ! IntegrationPointDataType
        CalParams%IntegrationPointData%NPlasticPoints = 0
        CalParams%IntegrationPointData%NNegativePlasticPoints = 0
        CalParams%IntegrationPointData%NApexPoints = 0
        CalParams%IntegrationPointData%NTensionCutOffPoints = 0
        
        ! AbsorbingBoundariesType
        CalParams%AbsorbingBoundaries%ANormalSolid = 1.0
        CalParams%AbsorbingBoundaries%AShearSolid = 1.0
        CalParams%AbsorbingBoundaries%AWater = 1.0
        CalParams%AbsorbingBoundaries%AGas = 1.0
        CalParams%AbsorbingBoundaries%DeltaNormalSolid = 0.5
        CalParams%AbsorbingBoundaries%DeltaShearSolid = 0.5
        CalParams%AbsorbingBoundaries%DeltaWater = 0.5
        CalParams%AbsorbingBoundaries%DeltaGas = 0.5
        CalParams%AbsorbingBoundaries%VBMaterialSet = 0.5

        CalParams%AbsorbingBoundaries%ApplyNodeData = .false.

        ! RigidBodyType(vsize, dsize)
        ALLOCATE(RigidBodyType(NVECTOR, NDIM)::CalParams%RigidBody)
        CalParams%RigidBody%IsRigidBody = .false.
        CalParams%RigidBody%InternalForce = 0.0
        CalParams%RigidBody%Mass = 0.0
        CalParams%RigidBody%TractionForce = 0.0
        CalParams%RigidBody%GravityForce = 0.0
        CalParams%RigidBody%Acceleration = 0.0
        CalParams%RigidBody%Velocity = 0.0
        CalParams%RigidBody%CurrentTraction = 0.0
        CalParams%RigidBody%CurrentGravity = 0.0
        CalParams%RigidBody%InAcceleration = 0.0
        CalParams%RigidBody%RigidEntity = ID_UNDEFINED
        
        ! BoundaryConditionsType(vsize, dsize)
        CalParams%BoundaryConditions%ApplyInfiltrationRate = .false.
        CalParams%BoundaryConditions%ApplySeepageFace = .false.
        CalParams%BoundaryConditions%UseInfiltration = .false.
        CalParams%BoundaryConditions%UseSeepageFace = .false.

        ! ImplicitIntegrationType
        CalParams%ImplicitIntegration%Iteration = 0
        CalParams%ImplicitIntegration%MinDesiredIterations = 6
        CalParams%ImplicitIntegration%MaxDesiredIterations = 15
        CalParams%ImplicitIntegration%UpScaleFactor = 2.0
        CalParams%ImplicitIntegration%DownScaleFactor = 0.5
        CalParams%ImplicitIntegration%MaxIterations = 60
        CalParams%ImplicitIntegration%DoUseZLS = .false.
        CalParams%ImplicitIntegration%IsZeroLoadStep = .false.
        CalParams%ImplicitIntegration%DoUseArcLengthControl = .false.
        CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix = .false.
        CalParams%ImplicitIntegration%StiffnessIncreaseFactor = 0.0
        CalParams%ImplicitIntegration%DanglingElementFactor = 0.0
        CalParams%ImplicitIntegration%DoUseStepExtrapolation = .false.
        CalParams%ImplicitIntegration%DoCheckChangeJacobian = .false.
        CalParams%ImplicitIntegration%InitialMultiplier = 0.0
        CalParams%ImplicitIntegration%DoUseAutomaticLoadStepping = .false.
        CalParams%ImplicitIntegration%ScaleFactor = 1.0
        CalParams%ImplicitIntegration%MaximumJacobianRatio = 0.0
        CalParams%ImplicitIntegration%LimitJacobianChange = 0.0

        ! reals
        CalParams%TimeIncrement = 1.0e-15
        CalParams%TotalTime = 0.5
        CalParams%OverallRealTime = 0.0
        CalParams%DampingFactor = 0.0
        CalParams%ScalingMassFactor = 1.0
        CalParams%ToleratedErrorEnergy = 0.01
        CalParams%ToleratedErrorForce = 0.01
        CalParams%ToleratedErrorEnergyWater = 0.01
        CalParams%ToleratedErrorForceWater = 0.01
        CalParams%ToleratedDivergence = 1.0
        CalParams%TotalRealTime = 0.0
        CalParams%FricCoef = 0.0
        CalParams%Adhesion = 0.0
        CalParams%RequiredDegreeOfFilling = 0.9
        CalParams%FacStiffnessIncrease = 0.0
        CalParams%SectorAngle = 0.0
        CalParams%XMin = 0.0
        CalParams%XMax = 0.0
        CalParams%YMin = 0.0
        CalParams%ZMin = 0.0
        CalParams%CourantNumber = 0.98
        CalParams%Fix(2) = 0.0
        CalParams%K0MaxSurfaceSuction = 1d10
        CalParams%SoilSurfacePoint(2) = 0.0
        CalParams%LiquidSurfacePoint(2) = 0.0
        CalParams%InitialVerticalLoadK0 = 0.0
        CalParams%PrescribedVelocity = 0.0
        CalParams%SurfacePrescribedVelocity = 0.0
        CalParams%SurfacePrescribedVelocityCoordMin = 0.0
        CalParams%SurfacePrescribedVelocityCoordMax = 0.0
        CalParams%LiquidPressureCavitationThreshold = 0.0
        CalParams%FreeSurfaceFactor = 0.0
        CalParams%VirtualParticleMassFactor = 1.0
        CalParams%OutputRealTimeID = -1
        CalParams%LimitPorosity = 0.0
        CalParams%LiquidThresholdDensity = 0.0
        CalParams%BinghamYieldStress = 0.0
        CalParams%BinghamViscosity = 0.0
        CalParams%BinghamYoungModulus = 0.0
        CalParams%BinghamPoissonRatio = 0.0
        CalParams%FrictionalFluidFrictionAngle = 0.0
        CalParams%TwoLayerErgunLawDiameterPartic = 0.0
        CalParams%TwoLayerErgunLawDiameterPartic2 = 0.0
        CalParams%IntrinsicPermeability = 0.0
        CalParams%ERGUNCONSTANTA = 0.0
        CalParams%ERGUNCONSTANTB = 0.0
        CalParams%ThicknessSoilLayer(MAXNUMBERSOILLAYERS) = -1.0
        CalParams%BoreHoleData= 0.0
        CalParams%InitialWaterPressure = 0.0
        CalParams%BulkViscosityDamping1 = 0.0
        CalParams%BulkViscosityDamping2 = 0.0
        CalParams%MinimumDeterminantRatio = 0.75
        CalParams%ShrinkageMateriaPointPositionFactor = 0.0
        CalParams%GammaWater = 0.0
        
        ! integers
        CalParams%OutputNumberParticles = -1
        CalParams%OutputParticles(MAXOUTPUTPARTICLES) = -1
        CalParams%OutputElements(MAX_OUTPUT_ELEMENTS) = -1
        CalParams%OutputNumberNodes = -1
        CalParams%OutputNodes(1:MAXOUTPUTNODES) = -1
        CalParams%ReactionHistoryNodes(MAXREACTIONHISTORYNODES) = -1
        CalParams%OutputNumberOfTimeSteps = -1
        CalParams%OutputTimeStepID(MAXOUTPUTTIMESTEPS) = -1
        CalParams%OutputNumberOfLoadSteps = -1
        CalParams%TimeStep = 1
        CalParams%MaxTimeSteps = 10000
        CalParams%NThreads = 1
        CalParams%NMaterialPoints = 4
        CalParams%NLoadSteps = 1
        CalParams%IStep = 1
        CalParams%PreviouslyRealisedLoadStep = 0
        CalParams%PKernel32 = 0
        CalParams%CalculationType = 2
        CalParams%RowCounter = 0
        CalParams%FileCounter = 0
        CalParams%ParticleFileCounter = 0
        CalParams%GroupThreshold = 40
        CalParams%RecursionThreshold = 60
        CalParams%NVirtualParticles = 1
        CalParams%ComputationMethod = 0        
        CalParams%TrajectoryType = -1
        CalParams%NumberOfPhases = 1
        CalParams%NumberOfMaterials = -1
        CalParams%FirstSolidMaterialIndex = -1
        CalParams%SecondSolidMaterialIndex = -1
        CalParams%NumberSubmergedCalculation = 1
        CalParams%NumberExcavatedVolumes = -1
        CalParams%ExcavatedVolumeID(20) = -1
        CalParams%ExcavationStages(20,2) = -1
        CalParams%NumberSoilLayers = 1
        CalParams%OutputCurvesIntervals = 1
        CalParams%NumberSkipConvection = 1
#ifdef Debug
        CalParams%FeedbackLevel = FEEDBACK_LEVEL_DEBUG
#else
        CalParams%FeedbackLevel = -1
#endif       

        ! logicals
        CalParams%ApplyContactMeshBoundary = .false.
        CalParams%ApplyQuasiStatic = .false.
        CalParams%TwoLayerApplyUpdatePermeability = .false.
        CalParams%TwoLayerApplyErgunLaw = .false.
        CalParams%TwoLayerApplyNoTensStressLiqMPwLiqStatus = .false.
        CalParams%ApplyStrainSmoothing = .false.
        CalParams%ApplyStrainSmoothingLiquidTwoLayer = .false.
        CalParams%ApplyDetectLiquidFreeSurface = .false.
        CalParams%ApplyContactAlgorithm = .false.
        CalParams%ApplyTractionContact = .false.
        CalParams%AssembleSoilLoadVectorFromMP = .false.
        CalParams%AssembleWaterLoadVectorFromMP = .false.
        CalParams%AssembleGasLoadVectorFromMP = .false.
        CalParams%ApplyPrescribedDisplacements = .false.
        CalParams%ApplySurfacePrescribedVelocity = .false.
        CalParams%ApplyMeshSmoothing = .false.
        CalParams%ApplyUpdateWeights = .false.
        CalParams%ApplyAbsorbingBoundary = .false.
        CalParams%ApplyCPSDamping = .false.
        CalParams%ApplyFEMtoMPM = .false.
        CalParams%ApplyResetDisplacements = .false.
        CalParams%ApplySubmergedCalculation = .false.
        CalParams%ApplyMaterialUpdate = .false.
        CalParams%ApplyObjectiveStress = .false.
        CalParams%ApplyRemoveSolidFixities = .false.
        CalParams%ApplyRemoveLiquidFixities = .false.
        CalParams%ApplyRemoveGasFixities = .false.
        CalParams%ApplyK0Procedure = .false.
        CalParams%ApplyFixParticlesK0 = .false.
        CalParams%ApplyEmptyElements = .false.
        CalParams%ApplyEffectiveStressAnalysis = .false.
        CalParams%OutputEachLoadStep = .true.
        CalParams%OutputDebugData = .false.
        CalParams%ApplyQuickCheckOutput = .true.
        CalParams%ApplyMassScaling = .false.
        CalParams%ApplyPrescribedVelocity = .false.
        CalParams%CorrectNormalsForSlice = .false.
        CalParams%ApplyStructureTrajectory = .false.
        CalParams%ApplyPorosityUpdate = .false.
        CalParams%ApplyBinghamFluid = .false.
        CalParams%ApplyFrictionalFluid = .false.
        CalParams%ApplyImplicitQuasiStatic = .false.
        CalParams%ApplyExcavation = .false.
        CalParams%ExternalSoilModel = .false.
        CalParams%ApplyFixedSolidSkeleton = .false.
        CalParams%SkipConvection  = .false.
        CalParams%AutomaticSkipConvection = .false.
        CalParams%ApplyBulkViscosityDamping = .false.
        CalParams%OutputBasicData = .false.
        CalParams%DoMapStresses = .false.
        CalParams%MonoPileSliceContact = .false.
        CalParams%UseUniformMaterialPointWeight = .false.
        CalParams%ApplyInitialVelocityonMP = .false.
        CalParams%ApplyPartialSaturation = .false. 
        CalParams%IsVTKBinary = .true.
        CalParams%ApplySmootheningLiquidPressureIncrement = .false.
         
        ! strings
        CalParams%CPSversion = 'not defined'
        CalParams%GOMversion = 'not defined'
        CalParams%StartDate = 'N/A'
        CalParams%EndDate = 'N/A'
        CalParams%StartTime = 'N/A'
        CalParams%EndTime = 'N/A'
        CalParams%Visualization = 'not defined'
       end subroutine InitialiseCalculationParameters
    
    
       logical function IsMPMComputation()
          IsMPMComputation =  &
            (CalParams%ComputationMethod==MPM_MIXED_INTEGRATION).or. &
            (CalParams%ComputationMethod==MPM_MP_INTEGRATION)
        end function IsMPMComputation

        logical function IsSmallDeformationFEMComputation()
          IsSmallDeformationFEMComputation = (CalParams%ComputationMethod==FEM)
        end function IsSmallDeformationFEMComputation
        
        logical function IsFEMComputation()
          IsFEMComputation = (CalParams%ComputationMethod==FEM).or.(CalParams%ComputationMethod==UL_FEM)
        end function IsFEMComputation

        logical function IsULFEMComputation()
          IsULFEMComputation = CalParams%ComputationMethod==UL_FEM
        end function IsULFEMComputation
        
        logical function IsMPMWithMPIntegration()
          IsMPMWithMPIntegration = CalParams%ComputationMethod==MPM_MP_INTEGRATION
        end function IsMPMWithMPIntegration

        logical function IsMPMWithMixedIntegration()
          IsMPMWithMixedIntegration = CalParams%ComputationMethod==MPM_MIXED_INTEGRATION
        end function IsMPMWithMixedIntegration

        logical function IsFollowUpPhase()
          IsFollowUpPhase = (CalParams%IStep/=1)
        end function IsFollowUpPhase
        
        subroutine UpdateLoadPhaseCounters()
        !**********************************************************************
        !
        !   Function : Increases the load phase counters
        !
        !**********************************************************************

        implicit none

          CalParams%TimeStep = 0
          CalParams%TotalRealTime = 0.0

        end subroutine UpdateLoadPhaseCounters

        
        subroutine ReadCalculationParameters()
        !**********************************************************************
        !
        ! Function: Determines CPS file version and calls respective CPS reader
        !
        !**********************************************************************
        
        implicit none
        
          !  local variables       
          character(len=MAX_FILENAME_LENGTH) :: FileName, FileVersion
          integer(INTEGER_TYPE) :: FileUnit
          character(len=255) :: BName
          integer(INTEGER_TYPE) :: ios

          FileName = trim(CalParams%FileNames%ProjectName)//CPS_FILE_EXTENSION//CalParams%FileNames%LoadStepExt
          FileUnit = TMP_UNIT
          
          ! check if CPS file exists in project folder, otherwise give error and stop execution
          if ( FExist(trim(FileName)) ) then
            call GiveMessage('Reading CPS file: ' // trim(FileName) )  
          else
            call GiveError('CPS file does not exist!' // NEW_LINE('A') // 'required CPS file: ' // trim(FileName) )
          end if
          
          ! open CPS file
          call FileOpen(FileUnit, trim(FileName))
          
          ! determine current version of CPS file 
          read(FileUnit, '(A)', iostat=ios) BName ! NB: if no version is specified in the header of the CPS file, the default case will be chosen automatically below
          call Assert( ios == 0, 'CPS file: Can''t read flag from CPS file.' )
          FileVersion = trim(BName)
          
          ! read CPS data
          select case (FileVersion) ! read CPS data depending on file version
              
            case (Anura3D_v2023)
              call ReadCPS(FileUnit,FileVersion)
            case (Anura3D_v2022)
              call ReadCPS(FileUnit,FileVersion)
            case (Anura3D_v2021)
              call ReadCPS(FileUnit,FileVersion)   
              
            case (Anura3D_v2019_2)
              call ReadCPS(FileUnit,FileVersion)                       
              
            case default  
              call GiveError('Wrong version of CPS file!' // NEW_LINE('A') // 'Supported CPS versions: ' &
                // trim(Anura3D_v2019_2) // ', ' &
                // trim(Anura3D_v2021) // ', ' &
                // trim(Anura3D_v2022) // ', ' &
                // trim(Anura3D_v2023)     // '.' )
              
          end select

          ! close CPS file    
          close(FileUnit)
      
        end subroutine ReadCalculationParameters
              
        
        subroutine ReadCPS(FileUnit,FileVersion)
        !**********************************************************************
        !
        ! Function : Reads input variables from the CPS file
        !
        !            CPS version 2022, 2021, 2019.2
        !            
        !
        !**********************************************************************
        
        implicit none

          character(len=MAX_FILENAME_LENGTH), intent(in) :: FileVersion
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! local variables
          integer(INTEGER_TYPE) :: I
          integer(INTEGER_TYPE) :: ios ! used for error control
          integer(INTEGER_TYPE) :: DumI(3)
          character(len=255) :: DumS
          real(REAL_TYPE) :: DumR(4)
          character(len=255) :: DumFileName
          character(len=21) :: messageIOS = 'CPS file: Can''t read '
          character(len=255) :: TName, BName

          ! set CPS version number
            select case (FileVersion)
            case (Anura3D_v2023) 
                call GiveMessage('Reading... ' // Anura3D_v2023)
                CalParams%CPSversion = Anura3D_v2021
            case (Anura3D_v2022) 
                call GiveMessage('Reading... ' // Anura3D_v2022)
                CalParams%CPSversion = Anura3D_v2021 
            case (Anura3D_v2021) 
                call GiveMessage('Reading... ' // Anura3D_v2021)
                CalParams%CPSversion = Anura3D_v2021 
             case (Anura3D_v2019_2) 
                call GiveMessage('Reading... ' // Anura3D_v2019_2)
                CalParams%CPSversion = Anura3D_v2021 ! external version
            end select
                    
          do
              
            read(FileUnit, '(A)', iostat=ios) TName
            call Assert( ios == 0, 'CPS file: Can''t read flag from CPS file.' )
            BName = Upcase(TName)
            call WriteInLogFile('Reading flag from CPS file: ' // trim(BName), FEEDBACK_LEVEL_ALWAYS)

            ! TIME DATA
            if (trim(BName)=='$$NUMBER_OF_LOADSTEPS') then ! I [ 0 < I < 1000 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) > 0, 'CPS file: ' //trim(BName)// ' must be larger than 0.' )
              call AssertWarning( DumI(1) <= 100, 'CPS file: ' //trim(BName)// ' is larger than 100.' )
              CalParams%NLoadSteps = DumI(1)
              
            else if (trim(BName)=='$$TIME_PER_LOADSTEP') then ! R [ 0.0 < R ]
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) > 0.0, 'CPS file: ' //trim(BName)// ' must be larger than 0.0.' )
              call AssertWarning( DumR(1) <= 10.0, 'CPS file: ' //trim(BName)// ' is larger than 10.0.' )
              CalParams%TotalTime = DumR(1)
              
            else if (trim(BName)=='$$TOTAL_TIME') then ! R [ 0.0 <= R ] ---RELEASE----> this shouldn't be a parameter in CPS but taken from ENG file if necessary.
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) >= 0.0, 'CPS file: ' //trim(BName)// ' must be larger than or equal to 0.0.' )
              CalParams%OverallRealTime = DumR(1)
              
            else if (trim(BName)=='$$COURANT_NUMBER') then ! R [ 0.0 < R <= 1.0]
              read(FileUnit, *, iostat=ios) DumR(1) 
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) > 0.0, 'CPS file: ' //trim(BName)// ' must be larger than 0.0.' )
              CalParams%CourantNumber = DumR(1) 

            ! CALCULATION DATA  
            else if (trim(BName)=='$$COMPUTATION_METHOD') then ! S [ MPM-MIXED, MPM-MP, FEM, UL-FEM ]
              read(FileUnit, *, iostat=ios) DumS
              call Assert( ios == 0, messageIOS//trim(BName) )
              if (trim(Upcase(DumS))==MPM_MIXED_INTEGRATION_SPECIFIER) then
                CalParams%ComputationMethod = MPM_MIXED_INTEGRATION
              else if (trim(Upcase(DumS))==MPM_MP_INTEGRATION_SPECIFIER) then
                CalParams%ComputationMethod = MPM_MP_INTEGRATION
              else if (trim(Upcase(DumS))==FEM_SPECIFIER) then
                CalParams%ComputationMethod = FEM
              else if (trim(Upcase(DumS))==UL_FEM_SPECIFIER) then
                CalParams%ComputationMethod = UL_FEM
              else
                call GiveError( 'CPS file: ' //trim(BName)// ' must be ' //MPM_MIXED_INTEGRATION_SPECIFIER// ' or ' &
                                 //MPM_MP_INTEGRATION_SPECIFIER// ' or ' //FEM_SPECIFIER// ' or ' //UL_FEM_SPECIFIER// '.' )
              end if
              if (IsFollowUpPhase()) call CheckChangeFEMtoMPM()
              
            ! GRAVITY DATA
            else if (trim(BName)=='$$GRAVITY_ACCELERATION') then ! R [ 0.0 <= R ]
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) >= 0.0, 'CPS file: ' //trim(BName)// ' must be larger than or equal to 0.0.' )
              call AssertWarning( DumR(1) == DEFAULT_GRAVITY_ACCELERATION, 'CPS file: Default value of ' //trim(BName)//' changed.')
              CalParams%GravityData%GAccel = DumR(1)
              
            else if (trim(BName)=='$$GRAVITY_VECTOR') then ! R(1:NDIM) [ -1.0 <= R(i) <= +1.0 ]
              DumR = 0.0
              read(FileUnit, *, iostat=ios) DumR(1:NDIM)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) >= -1.0, 'CPS file: ' //trim(BName)// ' must be in the range -1.0 to +1.0.' )
              call Assert( DumR(1) <= +1.0, 'CPS file: ' //trim(BName)// ' must be in the range -1.0 to +1.0.' )
              call Assert( DumR(2) >= -1.0, 'CPS file: ' //trim(BName)// ' must be in the range -1.0 to +1.0.' )
              call Assert( DumR(2) <= +1.0, 'CPS file: ' //trim(BName)// ' must be in the range -1.0 to +1.0.' )
              call Assert( DumR(3) >= -1.0, 'CPS file: ' //trim(BName)// ' must be in the range -1.0 to +1.0.' ) ! DumR(3) = 0 for 2D case
              call Assert( DumR(3) <= +1.0, 'CPS file: ' //trim(BName)// ' must be in the range -1.0 to +1.0.' ) ! DumR(3) = 0 for 2D case
              DumR(4) = sqrt( DumR(1)**2 + DumR(2)**2 + DumR(3)**2 ) ! determine length of vector
              call Assert( DumR(4) <= 1.01, 'CPS file: ' //trim(BName)// ' must have unit length.' )
              CalParams%GravityData%GravityVector(1:NDIM) = DumR(1:NDIM)
              
            ! LOADING DATA  

            else if (trim(BName)=='$$GRAVITY_LOAD') then ! R(1:2) [ 0.0 <= R(i) ]
              read(FileUnit, *, iostat=ios) DumS, DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == LOAD_TYPE_LINEAR .or. DumS == LOAD_TYPE_STEP .or. DumS == LOAD_TYPE_OFF, 'CPS file: Load type in ' //trim(BName)// ' must be "' & 
                  //LOAD_TYPE_LINEAR// '" or "' //LOAD_TYPE_STEP// '" or "' //LOAD_TYPE_OFF// '".' )
              CalParams%Multipliers%GravityLoadType = DumS
              if (CalParams%Multipliers%GravityLoadType == LOAD_TYPE_OFF) then
                CalParams%Multipliers%GravityRealised = 0.0
                CalParams%Multipliers%GravityFinal = 0.0
              else
                CalParams%Multipliers%GravityRealised = DumR(1)
                CalParams%Multipliers%GravityFinal = DumR(2)
              end if
              

            else if (trim(BName)=='$$SOLID_TRACTION') then ! R(1:2) [ 0.0 <= R(i) ]
              read(FileUnit, *, iostat=ios) DumS, DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == LOAD_TYPE_LINEAR .or. DumS == LOAD_TYPE_STEP .or. DumS == LOAD_TYPE_OFF, 'CPS file: Load type in ' //trim(BName)// ' must be "' & 
                  //LOAD_TYPE_LINEAR// '" or "' //LOAD_TYPE_STEP// '" or "' //LOAD_TYPE_OFF// '".' )
              CalParams%Multipliers%SolidALoadType(1) = DumS
              if (CalParams%Multipliers%SolidALoadType(1) == LOAD_TYPE_OFF) then
                CalParams%Multipliers%SolidARealised(1) = 0.0
                CalParams%Multipliers%SolidAFinal(1) = 0.0
              else    
                CalParams%Multipliers%SolidARealised(1) = DumR(1)
                CalParams%Multipliers%SolidAFinal(1) = DumR(2)
              end if
              
           else if (trim(BName)=='$$SOLID_TRACTION_B') then ! R(1:2) [ 0.0 <= R(i) ]
              read(FileUnit, *, iostat=ios) DumS, DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == LOAD_TYPE_LINEAR .or. DumS == LOAD_TYPE_STEP .or. DumS == LOAD_TYPE_OFF, 'CPS file: Load type in ' //trim(BName)// ' must be "' & 
                  //LOAD_TYPE_LINEAR// '" or "' //LOAD_TYPE_STEP// '" or "' //LOAD_TYPE_OFF// '".' )
              CalParams%Multipliers%SolidALoadType(2) = DumS
              if (CalParams%Multipliers%SolidALoadType(2) == LOAD_TYPE_OFF) then
                CalParams%Multipliers%SolidARealised(2) = 0.0
                CalParams%Multipliers%SolidAFinal(2) = 0.0
              else    
                CalParams%Multipliers%SolidARealised(2) = DumR(1)
                CalParams%Multipliers%SolidAFinal(2) = DumR(2)
              end if
              

            else if (trim(BName)=='$$LIQUID_PRESSURE') then ! R(1:2) [ 0.0 <= R(i) ]
              read(FileUnit, *, iostat=ios) DumS, DumR(1:2)
      
              CalParams%Multipliers%WaterALoadType(1) = DumS
              if (CalParams%Multipliers%WaterALoadType(1) == LOAD_TYPE_OFF) then
                CalParams%Multipliers%WaterARealised(1) = 0.0
                CalParams%Multipliers%WaterAFinal(1) = 0.0
              else  
                CalParams%Multipliers%WaterARealised(1) = DumR(1)
                CalParams%Multipliers%WaterAFinal(1) = DumR(2)
              end if
              
              else if (trim(BName)=='$$LIQUID_PRESSURE_B') then ! R(1:2) [ 0.0 <= R(i) ]
                  read(FileUnit, *, iostat=ios) DumS, DumR(1:2)
                  call Assert( ios == 0, messageIOS//trim(BName) )
                  call Assert( DumS == LOAD_TYPE_LINEAR .or. DumS == LOAD_TYPE_STEP .or. DumS == LOAD_TYPE_OFF, 'CPS file: Load type in ' //trim(BName)// ' must be "' &
                      //LOAD_TYPE_LINEAR// '" or "' //LOAD_TYPE_STEP// '" or "' //LOAD_TYPE_OFF// '".' )
                  !call Assert( DumR(1) >= 0.0, 'CPS file: Multipliers in ' //trim(BName)// ' must be larger than or equal to 0.0.' )
                  !call Assert( DumR(2) >= 0.0, 'CPS file: Multipliers in ' //trim(BName)// ' must be larger than or equal to 0.0.' )
                  CalParams%Multipliers%WaterALoadType(2) = DumS
                  if (CalParams%Multipliers%WaterALoadType(2) == LOAD_TYPE_OFF) then
                      CalParams%Multipliers%WaterARealised(2) = 0.0
                      CalParams%Multipliers%WaterAFinal(2) = 0.0
                  else
                      CalParams%Multipliers%WaterARealised(2) = DumR(1)
                      CalParams%Multipliers%WaterAFinal(2) = DumR(2)
                  end if
              
              
            else if (trim(BName)=='$$GAS_PRESSURE') then ! R(1:2) [ 0.0 <= R(i) ]
              read(FileUnit, *, iostat=ios) DumS, DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == LOAD_TYPE_LINEAR .or. DumS == LOAD_TYPE_STEP .or. DumS == LOAD_TYPE_OFF, 'CPS file: Load type in ' //trim(BName)// ' must be "' & 
                  //LOAD_TYPE_LINEAR// '" or "' //LOAD_TYPE_STEP// '" or "' //LOAD_TYPE_OFF// '".' )
              CalParams%Multipliers%GasALoadType(1) = DumS
              if (CalParams%Multipliers%GasALoadType(1) == LOAD_TYPE_OFF) then
                CalParams%Multipliers%GasARealised(1) = 0.0
                CalParams%Multipliers%GasAFinal(1) = 0.0
              else  
                CalParams%Multipliers%GasARealised(1) = DumR(1)
                CalParams%Multipliers%GasAFinal(1) = DumR(2)
              end if
              
                          else if (trim(BName)=='$$GAS_PRESSURE_B') then ! R(1:2) [ 0.0 <= R(i) ]
              read(FileUnit, *, iostat=ios) DumS, DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumS == LOAD_TYPE_LINEAR .or. DumS == LOAD_TYPE_STEP .or. DumS == LOAD_TYPE_OFF, 'CPS file: Load type in ' //trim(BName)// ' must be "' & 
                  //LOAD_TYPE_LINEAR// '" or "' //LOAD_TYPE_STEP// '" or "' //LOAD_TYPE_OFF// '".' )
              CalParams%Multipliers%GasALoadType(2) = DumS
              if (CalParams%Multipliers%GasALoadType(2) == LOAD_TYPE_OFF) then
                CalParams%Multipliers%GasARealised(2) = 0.0
                CalParams%Multipliers%GasAFinal(2) = 0.0
              else  
                CalParams%Multipliers%GasARealised(2) = DumR(1)
                CalParams%Multipliers%GasAFinal(2) = DumR(2)
              end if
              
               else if (trim(BName)=='$$HYDRAULIC_HEAD') then ! R(1:2) [ 0.0 <= R(i) ]
                 read(FileUnit, *, iostat=ios) DumS, DumR(1) !:2) ! [  file, off ]
                 if ( (trim(DumS) == LOAD_TYPE_FILE) ) then

                     CalParams%PrescribedHead%HydraulicHead = .true.
                     CalParams%Multipliers%HydraulicHeadType  = DumS
                  else if (trim(DumS) == LOAD_TYPE_OFF) then
                     CalParams%PrescribedHead%HydraulicHead = .false.
                 else
                    call GiveError( 'CPS file: ' //trim(BName)// ' must be equal to ' &
                         // LOAD_TYPE_FILE // ' or ' // LOAD_TYPE_OFF // '.' )
                 end if
                 
              
                 
                 if (CalParams%Multipliers%HydraulicHeadType  == LOAD_TYPE_OFF) then
                     CalParams%Multipliers%HydraulicHeadRealised = 0.0
                 else
                     CalParams%Multipliers%HydraulicHeadRealised = DumR(1)
                 end if
              
               else if (trim(BName)=='$$APPLY_SEEPAGE_FACE') then ! I [ 0, 1 ]
                   read(FileUnit, *, iostat=ios) DumI(1)
                   call Assert( ios == 0, messageIOS//trim(BName) )
                   call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
                   if ( DumI(1) == 1 ) CalParams%BoundaryConditions%UseSeepageFace = .true.

             
              
          else if (trim(BName)=='$$APPLY_INFILTRATION') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%BoundaryConditions%UseInfiltration = .true.

        
                  
                              
            elseif (trim(BName)=='$$PRESCRIBED_VELOCITY') then 
              read(FileUnit, *, iostat=ios) DumS, DumR(1:2) !, DumFileName ! [ step, linear, file, off ] 
              if ( (trim(DumS) == LOAD_TYPE_LINEAR).or.(trim(DumS) == LOAD_TYPE_STEP).or.(trim(DumS) == LOAD_TYPE_FILE) ) then
                CalParams%PrescribedVelo%ApplyPrescribedVelo = .true.
              else if (trim(DumS) == LOAD_TYPE_OFF) then
                CalParams%PrescribedVelo%ApplyPrescribedVelo = .false.
              else
                call GiveError( 'CPS file: ' //trim(BName)// ' must be equal to ' &
                           // LOAD_TYPE_LINEAR // ' or ' // LOAD_TYPE_STEP // ' or ' // LOAD_TYPE_FILE // ' or ' // LOAD_TYPE_OFF // '.' )  
              end if    
              CalParams%Multipliers%VelocitySolidLoadType = DumS
              if (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_OFF) then
                CalParams%Multipliers%VelocitySolidRealised = 0.0
                CalParams%Multipliers%VelocitySolidRealised = 0.0

              else  
                CalParams%Multipliers%VelocitySolidRealised = DumR(1)
                CalParams%Multipliers%VelocitySolidFinal = DumR(2)
              end if
              
                          
            ! QUASI-STATIC CALCULATION
            else if (trim(BName)=='$$QUASISTATIC_CONVERGENCE') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyQuasiStatic = .true.
              
            else if (trim(BName)=='$$TOLERATED_ERROR_SOLID') then ! R(1:2) [ 0.0 < R(i) ]
              read(FileUnit, *, iostat=ios) DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) > 0.0, 'CPS file: ' //trim(BName)// ' must be larger than 0.0.' )
              call Assert( DumR(2) > 0.0, 'CPS file: ' //trim(BName)// ' must be larger than 0.0.' )
              call AssertWarning( DumR(1) < 1.0, 'CPS file: ' //trim(BName)// ' is larger than 1.0.' )
              call AssertWarning( DumR(2) < 1.0, 'CPS file: ' //trim(BName)// ' is larger than 1.0.' )
              CalParams%ToleratedErrorEnergy = DumR(1)
              CalParams%ToleratedErrorForce = DumR(2)
              
            else if (trim(BName)=='$$TOLERATED_ERROR_LIQUID') then ! R(1:2) [ 0.0 < R(i) ]
              read(FileUnit, *, iostat=ios) DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) > 0.0, 'CPS file: ' //trim(BName)// ' must be larger than 0.0.' )
              call Assert( DumR(2) > 0.0, 'CPS file: ' //trim(BName)// ' must be larger than 0.0.' )
              call AssertWarning( DumR(1) < 1.0, 'CPS file: ' //trim(BName)// ' is larger than 1.0.' )
              call AssertWarning( DumR(2) < 1.0, 'CPS file: ' //trim(BName)// ' is larger than 1.0.' )
              CalParams%ToleratedErrorEnergyWater = DumR(1)
              CalParams%ToleratedErrorForceWater = DumR(2)
              
            else if (trim(BName)=='$$MAXIMUM_TIME_STEPS') then ! I [ 0 < I ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) > 0, 'CPS file: ' //trim(BName)// ' must be larger than 0.' )
              call AssertWarning( DumI(1) <= 1e5, 'CPS file: ' //trim(BName)// ' is larger than 100 000.' )
              CalParams%MaxTimeSteps = DumI(1)
                      
            ! QUASI-STATIC IMPLICIT INTEGRATION
            else if (trim(BName)=='$$APPLY_IMPLICIT_QUASISTATIC') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyImplicitQuasiStatic = .true.
                
            else if (trim(BName)=='$$QS_MIN_ITER') then 
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ImplicitIntegration%MinDesiredIterations = DumI(1)
                
            else if (trim(BName)=='$$QS_MAX_ITER') then 
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ImplicitIntegration%MaxDesiredIterations = DumI(1)
                
            else if (trim(BName)=='$$QS_UPSCALE_FACTOR') then 
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ImplicitIntegration%UpScaleFactor = DumR(1)
                
            else if (trim(BName)=='$$QS_DOWNSCALE_FACTOR') then 
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ImplicitIntegration%DownScaleFactor = DumR(1)
                
            else if (trim(BName)=='$$QS_MAX_ITERATIONS') then 
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ImplicitIntegration%MaxIterations = DumI(1)
                
            else if (trim(BName)=='$$APPLY_ZLS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ImplicitIntegration%DoUseZLS = .true.
                
            else if (trim(BName)=='$$APPLY_ARC_LENGTH_CONTROL') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ImplicitIntegration%DoUseArcLengthControl = .true.
                
            else if (trim(BName)=='$$APPLY_ENHANCED_STIFFNESS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix = .true.
                
            else if (trim(BName)=='$$QS_STIFFNESS_INCREASE_FACTOR') then 
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ImplicitIntegration%StiffnessIncreaseFactor = DumR(1)
                
            else if (trim(BName)=='$$QS_DANGLING_ELEMENT_FACTOR') then 
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ImplicitIntegration%DanglingElementFactor = DumR(1)
                
            else if (trim(BName)=='$$APPLY_STEP_EXTRAPOLATION') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ImplicitIntegration%DoUseStepExtrapolation = .true.
                
            else if (trim(BName)=='$$APPLY_JACOBIAN_CHECK') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ImplicitIntegration%DoCheckChangeJacobian = .true.
                
            else if (trim(BName)=='$$LIMIT_JACOBIAN_CHANGE') then 
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ImplicitIntegration%LimitJacobianChange = DumR(1)
                
            else if (trim(BName)=='$$APPLY_AUTOMATIC_LOADSTEPPING') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ImplicitIntegration%DoUseAutomaticLoadStepping = .true.
                
            else if (trim(BName)=='$$INCR_MULTIPLIER_SOILA') then
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%Multipliers%SolidAIncrement = DumR(1)
                               
            else if (trim(BName)=='$$INCR_MULTIPLIER_GRAVITY') then
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%Multipliers%GravityIncrement = DumR(1)
                
            else if (trim(BName)=='$$INCR_MULTIPLIER_WATERA') then
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%Multipliers%WaterAIncrement = DumR(1)
                
            else if (trim(BName)=='$$INCR_MULTIPLIER_GASA') then
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%Multipliers%GasAIncrement = DumR(1)
                
            ! MASS SCALING
            else if (trim(BName)=='$$MASS_SCALING') then ! I [ 0, 1 ] // R [ 1.0 <= R <= 10 000 ]
              read(FileUnit, *, iostat=ios) DumI(1), DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              !call Assert( DumR(1) >= 1.0, 'CPS file: ' //trim(BName)// ' must be larger than or equal to 1.0.' )
              !call Assert( DumR(1) <= 10000.0, 'CPS file: ' //trim(BName)// ' must be smaller than or equal to 10 000.' )

              if ( DumI(1) == 1 ) then
                  CalParams%ApplyMassScaling = .true.
                  CalParams%ScalingMassFactor = DumR(1)
              else
                  CalParams%ApplyMassScaling = .false.
                  CalParams%ScalingMassFactor = 1.0
              end if
              
            ! DAMPING
            else if (trim(BName)=='$$HOMOGENEOUS_LOCAL_DAMPING') then ! I [ 0, 1 ] // R [ 0.0 <= R < 1.0 ]
              read(FileUnit, *, iostat=ios) DumI(1), DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              call Assert( DumR(1) >= 0.0, 'CPS file: ' //trim(BName)// ' must be larger than or equal to 0.0.' )
              call Assert( DumR(1) < 1.0, 'CPS file: ' //trim(BName)// ' must be smaller than 1.0.' )
              if ( DumI(1) == 1 ) CalParams%ApplyCPSDamping = .true.
              CalParams%DampingFactor = DumR(1)

            ! BULK VISCOSITY DAMPING
            else if (trim(BName)=='$$BULK_VISCOSITY_DAMPING') then  ! I [ 0, 1 ] // R(1:2) [ 0.0 < R(i) ]
              read(FileUnit, *, iostat=ios) DumI(1), DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) then
                CalParams%ApplyBulkViscosityDamping = .true.
                call Assert( DumR(1) > 0.0, 'CPS file: ' //trim(BName)// ' must be larger than 0.0.' )
                call Assert( DumR(2) > 0.0, 'CPS file: ' //trim(BName)// ' must be larger than 0.0.' )
                CalParams%BulkViscosityDamping1 = DumR(1)
                CalParams%BulkViscosityDamping2 = DumR(2)
              end if  
              
            ! DIVERGENCE CHECK
            else if (trim(BName)=='$$APPLY_DIVERGENCE_CHECK') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ConvergenceCheck%ApplyDivergenceCheck = .true.

            ! OPTIMIZATION: Speeding up
            else if (trim(BName)=='$$SKIP_CONVECTION_TIME_STEPS') then
              read (FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) >= 1, 'CPS file: ' //trim(BName)// ' must be larger than or equal to 1.' )
              call GiveWarning('"Automatic Skip Convection" Functionality currently under development!')
              
              if ( DumI(1) > 1 ) CalParams%SkipConvection = .true.
              if ( DumI(1) > 1 ) CalParams%NumberSkipConvection =  DumI(1)

            else if (trim(BName)=='$$AUTOMATIC_SKIP_CONVECTION') then
              read (FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%AutomaticSkipConvection = DumI(1) == 1
              call GiveWarning('"Automatic Skip Convection" Functionality currently under development!')
              
            else if (trim(BName)=='$$MINIMUM_DETERMINANT_RATIO') then
              read (FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%MinimumDeterminantRatio = DumR(1)
             
            ! PARALLELIZATION
            else if (trim(BName)=='$$NUMBER_OF_THREADS') then
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%NThreads = DumI(1)

            else if (trim(BName)=='$$APPLY_NODE_DATA_AB' ) then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%AbsorbingBoundaries%ApplyNodeData = .true.

            ! CONTACT ALGORITHM
            else if (trim(BName)=='$$CONTACT_FORMULATION') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyContactAlgorithm = .true.

            else if (trim(BName)=='$$APPLY_CONTACT_MESH_BOUNDARY') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyContactMeshBoundary = .true.

            else if (trim(BName)=='$$APPLY_TRACTION_CONTACT') then ! I [ 0, 1]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyTractionContact = .true.

            else if (trim(BName)=='$$MONOPILE_SLICE_CONTACT') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%MonoPileSliceContact = .true.
			else if (trim(BName)=='$$APPLY_CONTACT_VELOCITY_SCALING' ) then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyContactVelocityScaling = .true.
			  
			else if (trim(BName)=='$$VELOCITY_SCALING_FACTOR'  ) then ! for external release this is specified in GOM file
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ContactScalingFactor = DumR(1) !Note that the scaling factor must be restricted between 1 and 0, I'm allowing to play with it
			
			else if (trim(BName)=='$$APPLY_CONTACT_NORMAL_CORRECTION') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 )  CalParams%ApplyContactNormalCorrection = .true.
			  
			else if (trim(BName)=='$$NUMBER_OF_CORRECTED_NORMALS' ) then
              read(FileUnit, *, iostat=ios)  CalParams%NCorrectedNormals
              call Assert( ios == 0, messageIOS//trim(BName) )

            else if (trim(BName)=='$$NODE_NORMAL_DATA' ) then
              if (CalParams%NCorrectedNormals>0) then
                do I = 1, CalParams%NCorrectedNormals
                  read(FileUnit, *, iostat=ios) DumR(1:(NVECTOR+1))				  
				  CalParams%ManuallyDefinedVectors(I,1:(NVECTOR+1))=DumR(1:(NVECTOR+1))
                  !call Assert( ios == 0, messageIOS//trim(BName) )
                end do
              end if

            ! EXCAVATION            
            else if (trim(BName)=='$$EXCAVATION') then
              read(FileUnit, *, iostat=ios) DumI(1)
              if ( DumI(1) > 0 ) CalParams%ApplyExcavation = .true.
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%NumberExcavatedVolumes = DumI(1) ! maximum number of excavated volumes = 20
              if (DumI(1) > 20) then
                  call GiveError('Number of excavated volumes is more than 20 (maximum value)' )
              end if
              do I = 1, CalParams%NumberExcavatedVolumes
                read(FileUnit, *, iostat=ios) DumI(1:3)
                call Assert( ios == 0, messageIOS//trim(BName) )       
                CalParams%ExcavatedVolumeID(I) = DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
                CalParams%ExcavationStages(I,1:2) = DumI(2:3)
              end do
 
            ! FIXED SOLID SKELETON
            else if (trim(BName)=='$$FIX_SOLID_SKELETON') then ! Input is an integer either [ 0 = false (solid MPs are allowed to move), 1 = true (solid MP position update is skipped)]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyFixedSolidSkeleton = .true.
 
            ! STRAIN SMOOTHING
            else if (trim(BName)=='$$STRAIN_SMOOTHING') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyStrainSmoothing = .true.
              
            ! DOUBLE-POINT FORMULATION  
            ! STRAIN SMOOTHING LIQUID 2 LAYER FORM
            else if (trim(BName)=='$$APPLY_STRAIN_SMOOTHING_LIQUID_TWOLAYERFORM') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyStrainSmoothingLiquidTwoLayer = .true.

            ! LIQUID PRESSURE CAVITATION THRESHOLD
            else if (trim(BName)=='$$LIQUID_CAVITATION') then
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) >= 0.0, 'CPS file: ' //trim(BName)// ' must be greater than or equal to zero.')
              CalParams%LiquidPressureCavitationThreshold = DumR(1)

            ! TWO LAYER FORMULATION :: CHECK IF APPLY UPDATE PERMEABILITY              
            else if (trim(BName)=='$$PERMEABILITY_UPDATE') then  
              read(FileUnit, *, iostat=ios) DumS
              call Assert( ios == 0, messageIOS//trim(BName) )
              if (trim(Upcase(DumS))=="CONSTANT_PERMEABILITY") then
                continue
              else if (trim(Upcase(DumS))=="DARCY_UPDATE") then
                CalParams%TwoLayerApplyUpdatePermeability = .true.
              else if (trim(Upcase(DumS))=="ERGUN_UPDATE") then
                CalParams%TwoLayerApplyUpdatePermeability = .true.
                CalParams%TwoLayerApplyErgunLaw = .true.
              end if

            else if (trim(BName)=='$$INTRINSIC_PERMEABILITY') then  
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%IntrinsicPermeability = DumR(1)
              
            else if (trim(BName)=='$$GRAIN_SIZE_DIAMETER') then  
              read(FileUnit, *, iostat=ios) DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%TwoLayerErgunLawDiameterPartic = DumR(1)
              CalParams%TwoLayerErgunLawDiameterPartic2 = DumR(2)
            
            else if (trim(BName)=='$$ERGUN_CONSTANTS') then  
              read(FileUnit, *, iostat=ios) DumR(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%ERGUNCONSTANTA = DumR(1)
              CalParams%ERGUNCONSTANTB = DumR(2)
              
            ! TWO LAYER FORMULATION :: NO TENSILE STRESS IN LIQUID MATERIAL POINT WITH LIQUID STATUS
            else if (trim(BName)=='$$NO_TENSILE_STRESS_LIQUID_MP_WITH_LIQUID_STATUS') then  ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%TwoLayerApplyNoTensStressLiqMPwLiqStatus = .true.

            ! TWO LAYER FORMULATION :: LIMIT POROSITY
            else if (trim(BName)=='$$MAXIMUM_POROSITY') then 
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) >= 0.0, 'CPS file: ' //trim(BName)// ' must be greater than or equal to 0.0.')
              CalParams%LimitPorosity = DumR(1)

            ! DETECT FREE SURFACE
            else if (trim(BName)=='$$DETECT_FREE_SURFACE') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1), DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyDetectLiquidFreeSurface = .true.
              call Assert( ( DumR(1) < 1.0 ) .and. ( DumR(1) >= 0.0 ), 'CPS file: ' //trim(BName)// ' must be 0.0 <= FreeSurfaceFactor < 1.0.')
              CalParams%FreeSurfaceFactor = DumR(1)

            ! RESET DISPLACEMENTS  
            else if (trim(BName)=='$$RESET_DISPLACEMENTS') then ! I [ 0, 1 ]
              read (FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyResetDisplacements = .true.
                
             
            else if (trim(BName)=='$$REMOVE_FIXITIES') then ! I [ 0, 1 ] 
              read(FileUnit, *, iostat=ios) DumI(1:3) ! solid / liquid / gas
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              call Assert( DumI(2) == 0 .or. DumI(2) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              call Assert( DumI(3) == 0 .or. DumI(3) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyRemoveSolidFixities = .true.
              if ( DumI(2) == 1 ) CalParams%ApplyRemoveLiquidFixities = .true.
              if ( DumI(3) == 1 ) CalParams%ApplyRemoveGasFixities = .true.
              
            ! K0-PROCEDURE
            else if (trim(BName)=='$$K0_PROCEDURE') then ! I [ 0, 1 ] 
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyK0Procedure = .true.
             
              else if (trim(BName)=='$$K0_MAX_SUCTION') then
              DumR(1) = 0.0 
              read(FileUnit, *, iostat=ios) DumR(1) 
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%K0MaxSurfaceSuction = DumR(1) !Value should be specified in kPa        
              
            else if (trim(BName)=='$$SURFACE_ELEVATION') then
              DumR(1) = 0.0 ! currently x-value is always equal to 0.0
              read(FileUnit, *, iostat=ios) DumR(2) ! y-value
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%SoilSurfacePoint(1) = DumR(1) ! x-value
              CalParams%SoilSurfacePoint(2) = DumR(2) ! y-value
              
            else if (trim(BName)=='$$INITIAL_VERTICAL_LOAD_K0' ) then   ! Flag read in CPS_001
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%InitialVerticalLoadK0 = DumR(1)
              
            else if (trim(BName)=='$$CONSIDERED_LOAD_K0' ) then         ! Flag read in the subsequent CPS files
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%InitialVerticalLoadK0 = DumR(1)                      
              
            else if (trim(BName)=='$$NUMBER_SOIL_LAYERS') then
              read(FileUnit, *, iostat=ios)  CalParams%NumberSoilLayers
              call Assert( ios == 0, messageIOS//trim(BName) )
                
            else if (trim(BName)=='$$THICKNESS_SOIL_LAYERS') then
              if (CalParams%NumberSoilLayers>0) then
                do I = 1, CalParams%NumberSoilLayers
                  read(FileUnit, *, iostat=ios) CalParams%ThicknessSoilLayer(I)
                  call Assert( ios == 0, messageIOS//trim(BName) )
                end do
              end if

            else if (trim(BName)=='$$LIQUID_SURFACE') then
              DumR(1) = 0 ! currently x-value is always equal to 0.0
              read(FileUnit, *, iostat=ios) DumR(2) ! y-value
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%LiquidSurfacePoint(1) = DumR(1) ! x-value
              CalParams%LiquidSurfacePoint(2) = DumR(2) ! y-value
                
            ! INITIAL WATER PRESSURE
            else if (trim(BName)=='$$INITIAL_WATER_PRESSURE') then
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%InitialWaterPressure = DumR(1)

            ! MATERIAL DATA
            else if (trim(BName)=='$$APPLY_MATERIAL_UPDATE') then
              read (FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyMaterialUpdate = .true.

            else if (trim(BName)=='$$APPLY_POROSITY_UPDATE') then
              read (FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyPorosityUpdate = .true.

        
            else if (trim(BName)=='$$SUBMERGED_CALCULATION') then
              read (FileUnit, *, iostat=ios) DumI(1:2)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplySubmergedCalculation = .true.
              CalParams%NumberSubmergedCalculation = DumI(2)
			  
			  !INITIAL VELOCITY ON MP
			else if (trim(BName)=='$$INITIAL_VELOCITY') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyInitialVelocityonMP = .true.
			  
            ! MPM SPECIFIC DATA
            else if (trim(BName)=='$$APPLY_OBJECTIVE_STRESS') then ! I [ 0, 1 ]
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyObjectiveStress = .true.

            else if (trim(BName)=='$$DEGREE_OF_FILLING') then ! R [ 0.0 <= R <= 1.0 ]
              read(FileUnit, *, iostat=ios) DumR(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) >= 0.0, 'CPS file: ' //trim(BName)// ' must be greater than or equal to 0.0.')
              call Assert( DumR(1) <= 1.0, 'CPS file: ' //trim(BName)// ' must be less than or equal to 1.0.')
              CalParams%RequiredDegreeOfFilling = DumR(1)
              
            else if (trim(BName)=='$$NUMBER_OF_ACTIVE_ELEMENTS') then ! I [ 0 <= I ] ! ---RELEASE-----> this should be stored elsewhere
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumR(1) >= 0.0, 'CPS file: ' //trim(BName)// ' must be greater than or equal to 0.')
              if ( CalParams%IStep==1 ) then 
                Counters%NAEl = 0
              else
                Counters%NAEl = DumI(1)
              end if   

            ! EMPTY ELEMENTS
            else if (trim(BName)=='$$APPLY_EMPTY_ELEMENTS') then
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyEmptyElements = .true.

            else if (trim(BName)=='$$MASS_FACTOR_VIRTUAL_PARTICLES') then
              read(FileUnit, *, iostat=ios) CalParams%VirtualParticleMassFactor
              call Assert( ios == 0, messageIOS//trim(BName) )

            else if (trim(BName)=='$$GROUP_THRESHOLD') then
              read(FileUnit, *, iostat=ios) CalParams%GroupThreshold
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%RecursionThreshold = max(CalParams%RecursionThreshold, int(3 * CalParams%GroupThreshold))
            
            ! VISUALIZATION
            else if (trim(BName)=='$$VISUALIZATION_OPTION') then ! S [ Paraview, GiD-ASCII, GiD-Binary, Paraview-GiD ]
              read(FileUnit, *, iostat=ios) DumS
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%Visualization = DumS
 
            ! OUTPUT
            else if (trim(BName)=='$$OUTPUT_NUMBER_OF_MATERIAL_POINTS') then
              read(FileUnit, *, iostat=ios) CalParams%OutputNumberParticles
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert(CalParams%OutputNumberParticles <= MAXOUTPUTPARTICLES, 'CPS-ERROR: maximum value for $$OUTPUT_NUMBER_OF_MATERIAL_POINTS is '//trim(String(MAXOUTPUTPARTICLES)))

            else if (trim(BName)=='$$OUTPUT_MATERIAL_POINTS') then
              if (CalParams%OutputNumberParticles>0) then
                do I = 1, CalParams%OutputNumberParticles
                  read(FileUnit, *, iostat=ios) CalParams%OutputParticles(I)
                  call Assert( ios == 0, messageIOS//trim(BName) )
                end do
              else
                ! read the dummy value
                read(FileUnit, *, iostat=ios) DumI(1)
                call Assert( ios == 0, messageIOS//trim(BName) )
              end if

            else if (trim(BName)=='$$OUTPUT_BASIC_DATA') then
              read (FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%OutputBasicData = dumI(1) == 1
              
            else if (trim(BName)=='$$OUTPUT_CURVES_INTERVAL') then
              read(FileUnit, *, iostat=ios) CalParams%OutputCurvesIntervals
              call Assert( ios == 0, messageIOS//trim(BName) )
              CalParams%OutputCurvesIntervals = max(CalParams%OutputCurvesIntervals, 1)

            else if (trim(BName)=='$$OUTPUT_NUMBER_OF_NODES') then
              read(FileUnit, *, iostat=ios) CalParams%OutputNumberNodes
              call Assert( ios == 0, messageIOS//trim(BName) )
              
            else if (trim(BName)=='$$OUTPUT_NODES') then
              read(FileUnit, *, iostat=ios) CalParams%OutputNodes(1), &
                                    CalParams%OutputNodes(2), & 
                                    CalParams%OutputNodes(3), &
                                    CalParams%OutputNodes(4), &
                                    CalParams%OutputNodes(5), &
                                    CalParams%OutputNodes(6), &
                                    CalParams%OutputNodes(7), &
                                    CalParams%OutputNodes(8), &
                                    CalParams%OutputNodes(9), &
                                    CalParams%OutputNodes(10)
              call Assert( ios == 0, messageIOS//trim(BName) )
              
            else if (trim(BName)=='$$OUTPUT_DEBUG_DATA') then
              read (FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%OutputDebugData = .true.

            else if (trim(BName)=='$$OUTPUT_EACH_LOAD_STEP') then
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 0 ) CalParams%OutputEachLoadStep = .false.

            else if (trim(BName)=='$$QUICK_CHECK_OUTPUT') then
              read(FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplyQuickCheckOutput = .true.

            else if (trim(BName)=='$$FILE_COUNTER') then
              read(FileUnit, *, iostat=ios) CalParams%FileCounter
              call Assert( ios == 0, messageIOS//trim(BName) )

            else if (trim(BName)=='$$PARTICLE_FILE_COUNTER') then
              read(FileUnit, *, iostat=ios) CalParams%ParticleFileCounter
              call Assert( ios == 0, messageIOS//trim(BName) )

            else if (trim(BName)=='$$FEEDBACK_LEVEL') then
              read(FileUnit, *, iostat=ios) CalParams%FeedbackLevel
              call Assert( ios == 0, messageIOS//trim(BName) )
              
            else if (trim(BName)=='$$APPLY_SMOOTHENING_LIQUID_PRESSURE_INCREMENT') then
              read (FileUnit, *, iostat=ios) DumI(1)
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI(1) == 0 .or. DumI(1) == 1, 'CPS file: ' //trim(BName)// ' must be 0 or 1.' )
              if ( DumI(1) == 1 ) CalParams%ApplySmootheningLiquidPressureIncrement = .true. 
              

            ! END OF CPS-FILE
            else if (trim(BName)=='$$END') then
              EXIT
            else
              if (ContainLetters(trim(BName))) then
                call GiveWarning('Unknown flag in CPS file: ' // trim(BName))
              endif
            end if

          end do
          
          call InitialiseDependentFlags(.true.,.true.,.true.)

        
        end subroutine ReadCPS
  
        subroutine CheckChangeFEMtoMPM()
        !**********************************************************************
        !
        ! Function :  Automatic detection of change of computational method from (UL)FEM to MPM.
        !             CalParams%ApplyFEMtoMPM is set to .true. when this is the case.
        !    
        !   First, the new computational method is checked to be MPM
        !   Second, the previous CPS file is scanned for the previous computational method
        !   Third, the previous computational methods is checked to be (UL)FEM
        !
        !**********************************************************************
            implicit none
            integer(INTEGER_TYPE) :: prev_unit, PreviousComputationalMethod
            character(len = 255) :: FileName, TName, BName, InputString

        !   The new computational method is checked to be MPM
            if((CalParams%ComputationMethod /= MPM_MIXED_INTEGRATION)  &
                              .and. (CalParams%ComputationMethod /= MPM_MIXED_INTEGRATION)) then
                return
            end if
        
        !   The previous CPS file is scanned for the previous computational method
            FileName = trim(CalParams%FileNames%ProjectName)//CPS_FILE_EXTENSION//CalParams%FileNames%PreviousStepExt

            if (.not.FExist(trim(FileName))) then
              call GiveError('Previous CPS file does not exist!' //NEW_LINE('A')// &
                             'File name: ' // trim(FileName))
            end if

            prev_unit = 88
            call FileOpen(prev_unit, trim(FileName))

            PreviousComputationalMethod = CalParams%ComputationMethod
            do
              read(prev_unit, '(A)') TName
              BName = Upcase(TName)
        
              if ((trim(BName)==('$$COMP_METHOD')).or.(trim(BName)==('$$COMPUTATION_METHOD'))) then
                read(prev_unit, *) InputString
                if (trim(Upcase(InputString))==MPM_MIXED_INTEGRATION_SPECIFIER) then
                  PreviousComputationalMethod = MPM_MIXED_INTEGRATION
                else if (trim(Upcase(InputString))==MPM_MP_INTEGRATION_SPECIFIER) then
                  PreviousComputationalMethod = MPM_MP_INTEGRATION
                else if (trim(Upcase(InputString))==FEM_SPECIFIER) then
                  PreviousComputationalMethod = FEM
                else if (trim(Upcase(InputString))==UL_FEM_SPECIFIER) then
                  PreviousComputationalMethod = UL_FEM
                else
                  call GiveError('$$COMP_METHOD not properly specified in previous CPS file.')
                end if
                EXIT
              else if (trim(BName)=='$$END') then
                call GiveError('$$COMP_METHOD not specified in previous CPS file.')
                EXIT
              end if
            end do
            
            close(prev_unit)
         
        !   The previous computational method is checked to be (UL)FEM
            if((PreviousComputationalMethod == FEM) .or. (CalParams%ComputationMethod == UL_FEM)) then
              CalParams%ApplyFEMtoMPM = .true.
            end if
            
        end subroutine CheckChangeFEMtoMPM
        
        
        
        logical function isThreePhaseCalculation() result(res)
        !**********************************************************************
        !
        !    Function: Checks whether the calculation is three phase or not.
        !    
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none

        res = CalParams%NumberOfPhases == 3

        end function isThreePhaseCalculation
        

        
        logical function isTwoPhaseCalculation() result(res)
        !**********************************************************************
        !
        !    Function: Checks whether the calculation is two phase or not.
        !    
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none

        res = CalParams%NumberOfPhases == 2

        end function isTwoPhaseCalculation
        
       
        
        logical function isContactActive() result(res)
        !**********************************************************************
        !
        !    Function: Checks whether the contact algorithm is applied or not.
        !    
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none

        res = CalParams%ApplyContactAlgorithm

        end function isContactActive



        subroutine InitialiseDependentFlags(IsAssemblageSoilLoadVectorSpecified,  &
                                            IsAssemblageWaterLoadVectorSpecified, &
                                            IsAssemblageGasLoadVectorSpecified)
        !**********************************************************************
        !
        !     Function:  Sets default values for assemblage of load vectors
        !                if assemblage of load vectors is not specified
        !
        !**********************************************************************

        implicit none

          ! Arguments
          logical, intent(in) :: IsAssemblageSoilLoadVectorSpecified, &
                                 IsAssemblageWaterLoadVectorSpecified, &
                                 IsAssemblageGasLoadVectorSpecified

          ! For FEM computations and MPM with moving mesh, load vectors are by default not assembled from MP forces

          if (.not.IsAssemblageSoilLoadVectorSpecified) then ! Consider default settings
            CalParams%AssembleSoilLoadVectorFromMP = IsMPMComputation().and.(.not.CalParams%ApplyMeshSmoothing)
          end if

          if (.not.IsAssemblageWaterLoadVectorSpecified) then ! Consider default settings
            CalParams%AssembleWaterLoadVectorFromMP = IsMPMComputation().and.(.not.CalParams%ApplyMeshSmoothing)
          end if

          if (.not.IsAssemblageGasLoadVectorSpecified) then ! Consider default settings
            CalParams%AssembleGasLoadVectorFromMP = IsMPMComputation().and.(.not.CalParams%ApplyMeshSmoothing)
          end if

          if (CalParams%AutomaticSkipConvection) then
            CalParams%SkipConvection = .true.
            if (CalParams%NumberSkipConvection <= 0) then
              CalParams%NumberSkipConvection =  huge(CalParams%NumberSkipConvection)
            else
              if (CalParams%NumberSkipConvection == 1) then
                call GiveWarning('Number of skip convection cannot be equal to 1 when automatic skip convection is used')
                call GiveWarning('The number of skip convection is ignored!')
                CalParams%NumberSkipConvection =  huge(CalParams%NumberSkipConvection)
              endif
            endif
          endif

        end subroutine InitialiseDependentFlags

        
        subroutine CreateResultFileHeader(LoadStep)
        !**********************************************************************
        !
        !    Function:  Determines the header of the BRF files.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none

          ! Arguments
          character(len = *), intent(in) :: LoadStep
          
          ! Local variables
          character(len = 80) :: Header
        
          Header = 'Anura3D'
          Header = trim(Header) // ' ' // trim('Version 1.0')
          Header = trim(Header) // ' ' // trim('Loadstep ' // LoadStep)

          CalParams%FileNames%ResultFileHeader = Header
      
        end subroutine CreateResultFileHeader


        subroutine GetStepExt(StepNumber, XXX)
        !**********************************************************************
        !
        !     Function:  Translates integer StepNumber into string XXX.
        !
        ! I   StepNumber : Integer
        ! O   XXX : File extension string
        !
        !**********************************************************************
        
        implicit none
        
          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: StepNumber
          character(len = *), intent(inout) :: XXX
          
          if (StepNumber<=999) then
            write(XXX, '(I3.3)') StepNumber
          else if (StepNumber<=9999) Then
            write(XXX, '(I4.4)') StepNumber
          else if (StepNumber<=99999) Then
            write(XXX, '(I5.5)') StepNumber
          else if (StepNumber<=999999) then
            write(XXX, '(I6.6)') StepNumber
          else if (StepNumber<=9999999) then
            write(XXX, '(I7.7)') StepNumber
          end if

        end subroutine GetStepExt


        subroutine GetTimeStep(StepNumber, XXX)
        !**********************************************************************
        !
        !  Function : Translates integer TimeStepNumber into string XXX
        !
        !  I   StepNumber : Integer
        !  O   XXX : File extension string
        !
        !**********************************************************************
        
        implicit none
        
          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: StepNumber
          character(len=*), intent(inout) :: XXX
          
          if (StepNumber<=999999) then
            write(XXX,'(I6.6)')StepNumber
          else
            call GiveWarning('Time step larger than allowed for output purposes.')
          end if

        end subroutine GetTimeStep

        
        logical function NotFinishedComputation()
        !**********************************************************************
        !
        !    Function: Checks whether the multipliers for external loads, water pressures,
        !              gas pressures, gravity loading have reached the final value.
        !    
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
          
          if (CalParams%ApplyImplicitQuasiStatic) then
            NotFinishedComputation =  &
              (CalParams%Multipliers%SolidACurrent(1)<CalParams%Multipliers%SolidAFinal(1)).or. &
              (CalParams%Multipliers%WaterACurrent(1)<CalParams%Multipliers%WaterAFinal(1)).or. &
              (CalParams%Multipliers%GravityCurrent<CalParams%Multipliers%GravityFinal)
          else
            NotFinishedComputation = CalParams%IStep<=CalParams%NLoadSteps
          end if

        end function NotFinishedComputation
        
        
        subroutine UpdateMultipliersForLoadStep()
        !**********************************************************************
        !
        !    Function: Updates the multipliers for external loads, water pressures,
        !              gas pressures, gravity loading from the current value of CalParams%IStep.
        !    
        !              The computation of the multiplier increments could be done
        !              outside the load step loop. In anticipation of possibly
        !              variable increments and for easing readability of the code
        !              it is kept in this routine. (No high computation effort.)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          real(REAL_TYPE) :: Multiplier, MultiplierPrevious
          real(REAL_TYPE) :: PreviousMultiplierSolidAIncrement
          real(REAL_TYPE) :: PreviousMultiplierWaterAIncrement
          real(REAL_TYPE) :: PreviousMultiplierWaterBIncrement
          real(REAL_TYPE) :: PreviousMultiplierGasAIncrement 
          real(REAL_TYPE) :: PreviousMultiplierGravityIncrement
          real(REAL_TYPE) :: PreviousMultiplierSoilVelocityIncrement
          integer(INTEGER_TYPE) :: ILoadSystem
          
          Multiplier =  &
            dble(CalParams%IStep - CalParams%PreviouslyRealisedLoadStep) /  &
            dble(CalParams%NLoadSteps - CalParams%PreviouslyRealisedLoadStep)
          MultiplierPrevious =  &
            dble(CalParams%IStep - 1 - CalParams%PreviouslyRealisedLoadStep) /  &
            dble(CalParams%NLoadSteps - CalParams%PreviouslyRealisedLoadStep)
          
          !Loop over load systems
          do ILoadSystem = 1,Counters%NSoilLoadSystems
          !Calculate multipliers solid
          if (CalParams%Multipliers%SolidALoadType(ILoadSystem)==LOAD_TYPE_OFF) then
              CalParams%Multipliers%SolidACurrent(ILoadSystem) = 0.0
              CalParams%Multipliers%SolidARealised(ILoadSystem) = 0.0
          else
            CalParams%Multipliers%SolidAIncrement(ILoadSystem) = &
            Multiplier * (CalParams%Multipliers%SolidAFinal(ILoadSystem) - CalParams%Multipliers%SolidARealised(ILoadSystem))
            CalParams%Multipliers%SolidACurrent(ILoadSystem) = CalParams%Multipliers%SolidARealised(ILoadSystem) + CalParams%Multipliers%SolidAIncrement(ILoadSystem)
            PreviousMultiplierSolidAIncrement = MultiplierPrevious*(CalParams%Multipliers%SolidAFinal(ILoadSystem) - CalParams%Multipliers%SolidARealised(ILoadSystem))
            CalParams%Multipliers%SolidAPrevious(ILoadSystem) = CalParams%Multipliers%SolidARealised(ILoadSystem) + PreviousMultiplierSolidAIncrement
            CalParams%Multipliers%SolidALoadStepIncrement(ILoadSystem) = CalParams%Multipliers%SolidACurrent(ILoadSystem) - CalParams%Multipliers%SolidAPrevious(ILoadSystem)
          end if
          end do
          
          !Loop over load systems
          do ILoadSystem = 1,Counters%NWaterLoadSystems
          !calculate multipliers liquid
          if (CalParams%Multipliers%WaterALoadType(ILoadSystem)==LOAD_TYPE_OFF) then
              CalParams%Multipliers%WaterACurrent(ILoadSystem) = 0.0
              CalParams%Multipliers%WaterARealised(ILoadSystem) = 0.0
          else 

              CalParams%Multipliers%WaterAIncrement(ILoadSystem) =  &
              Multiplier * (CalParams%Multipliers%WaterAFinal(ILoadSystem) - CalParams%Multipliers%WaterARealised(ILoadSystem))
            
            CalParams%Multipliers%WaterACurrent(ILoadSystem) = CalParams%Multipliers%WaterARealised(ILoadSystem) + CalParams%Multipliers%WaterAIncrement(ILoadSystem)
            PreviousMultiplierWaterAIncrement = MultiplierPrevious*(CalParams%Multipliers%WaterAFinal(ILoadSystem) - CalParams%Multipliers%WaterARealised(ILoadSystem))
            
            CalParams%Multipliers%WaterAPrevious(ILoadSystem) = CalParams%Multipliers%WaterARealised(ILoadSystem) + PreviousMultiplierWaterAIncrement
            CalParams%Multipliers%WaterALoadStepIncrement(ILoadSystem) = CalParams%Multipliers%WaterACurrent(ILoadSystem) - CalParams%Multipliers%WaterAPrevious(ILoadSystem)
           
          end if
          end do
        
          do ILoadSystem = 1,Counters%NGasLoadSystems
          if (CalParams%Multipliers%GasALoadType(ILoadSystem)==LOAD_TYPE_OFF) then
            CalParams%Multipliers%GasACurrent(ILoadSystem) = 0.0
            CalParams%Multipliers%GasARealised(ILoadSystem) = 0.0
          else      
            CalParams%Multipliers%GasAIncrement(ILoadSystem) =  &
            Multiplier * (CalParams%Multipliers%GasAFinal(ILoadSystem) - CalParams%Multipliers%GasARealised(ILoadSystem))
            
            CalParams%Multipliers%GasACurrent(ILoadSystem) = CalParams%Multipliers%GasARealised(ILoadSystem) + CalParams%Multipliers%GasAIncrement(ILoadSystem)
            
            PreviousMultiplierGasAIncrement = MultiplierPrevious*(CalParams%Multipliers%GasAFinal(ILoadSystem) - CalParams%Multipliers%GasARealised(ILoadSystem))
            CalParams%Multipliers%GasAPrevious(ILoadSystem) = CalParams%Multipliers%GasARealised(ILoadSystem) + PreviousMultiplierGasAIncrement
            CalParams%Multipliers%GasALoadStepIncrement(ILoadSystem) = CalParams%Multipliers%GasACurrent(ILoadSystem) - CalParams%Multipliers%GasAPrevious(ILoadSystem)
          end if
          end do !loop over load systems
          
          if (CalParams%Multipliers%GravityLoadType==LOAD_TYPE_OFF) then
            CalParams%Multipliers%GravityCurrent = 0.0
            CalParams%Multipliers%GravityRealised = 0.0
          else 
            CalParams%Multipliers%GravityIncrement =  &
            Multiplier * (CalParams%Multipliers%GravityFinal - CalParams%Multipliers%GravityRealised)
            
            CalParams%Multipliers%GravityCurrent = CalParams%Multipliers%GravityRealised + CalParams%Multipliers%GravityIncrement
            
            PreviousMultiplierGravityIncrement = MultiplierPrevious*(CalParams%Multipliers%GravityFinal - CalParams%Multipliers%GravityRealised)
            CalParams%Multipliers%GravityPrevious = CalParams%Multipliers%GravityRealised + PreviousMultiplierGravityIncrement
            CalParams%Multipliers%GravityLoadStepIncrement = CalParams%Multipliers%GravityCurrent - CalParams%Multipliers%GravityPrevious
          end if
          
          if (CalParams%Multipliers%VelocitySolidLoadType==LOAD_TYPE_OFF) then
            CalParams%Multipliers%VelocitySolidCurrent = 0.0
            CalParams%Multipliers%VelocitySolidRealised = 0.0
          else if (CalParams%Multipliers%VelocitySolidLoadType /= LOAD_TYPE_FILE) then 
              ! Prescribing velocity from file is time dependent only and is not dependent on load step
            CalParams%Multipliers%VelocitySolidIncrement =  &
            Multiplier * (CalParams%Multipliers%VelocitySolidFinal - CalParams%Multipliers%VelocitySolidRealised)
            
            CalParams%Multipliers%VelocitySolidCurrent = CalParams%Multipliers%VelocitySolidRealised + CalParams%Multipliers%VelocitySolidIncrement
            PreviousMultiplierSoilVelocityIncrement = MultiplierPrevious*(CalParams%Multipliers%VelocitySolidFinal - CalParams%Multipliers%VelocitySolidRealised)
            
            CalParams%Multipliers%VelocitySolidPrevious = CalParams%Multipliers%VelocitySolidRealised + PreviousMultiplierSoilVelocityIncrement
            CalParams%Multipliers%VelocitySolidLoadStepIncrement = CalParams%Multipliers%VelocitySolidCurrent - CalParams%Multipliers%VelocitySolidPrevious
          end if

          if (CalParams%RigidBody%IsRigidBody) then
            CalParams%RigidBody%CurrentTraction = CalParams%RigidBody%TractionForce * CalParams%Multipliers%SolidACurrent
            CalParams%RigidBody%CurrentGravity = CalParams%RigidBody%GravityForce * CalParams%Multipliers%GravityCurrent
          end if

        end subroutine UpdateMultipliersForLoadStep


        subroutine UpdateMultipliersForTimeDepencency()
        !**********************************************************************
        !
        !    Function: Updates the external load for time-dependency.
        !
        !              The increment CalParams%Multipliers%XXXXIncrement traces the
        !              change of amplitude, not the change of the multiplier with time!
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE) :: Increment
          integer(INTEGER_TYPE) :: ILoadSystem
                
          !Loop over load systems
          do ILoadSystem = 1,Counters%NSoilloadSystems
          if (CalParams%Multipliers%SolidALoadType(ILoadSystem) == LOAD_TYPE_LINEAR) then !Solid
            Increment = (CalParams%Multipliers%SolidALoadStepIncrement(ILoadSystem)) * &
                CalParams%TotalRealTime / CalParams%TotalTime
            
            CalParams%Multipliers%SolidACurrent(ILoadSystem) = CalParams%Multipliers%SolidAPrevious(ILoadSystem) + Increment
          end if  
          end do
          
          do ILoadSystem = 1,Counters%NWaterLoadSystems
          if (CalParams%Multipliers%WaterALoadType(ILoadSystem) == LOAD_TYPE_LINEAR) then !Water (NODES)
            Increment = (CalParams%Multipliers%WaterALoadStepIncrement(ILoadSystem)) * &
                CalParams%TotalRealTime / CalParams%TotalTime
            
            CalParams%Multipliers%WaterACurrent(ILoadSystem) = CalParams%Multipliers%WaterAPrevious(ILoadSystem) + Increment
          end if
          end do

          do ILoadSystem = 1,Counters%NGasLoadSystems
         if (CalParams%Multipliers%GasALoadType(ILoadSystem) == LOAD_TYPE_LINEAR) then !Gas
            Increment = (CalParams%Multipliers%GasALoadStepIncrement(ILoadSystem)) * &
                CalParams%TotalRealTime / CalParams%TotalTime
            
            CalParams%Multipliers%GasACurrent(ILoadSystem) = CalParams%Multipliers%GasAPrevious(ILoadSystem) + Increment
          end if
          end do

        
          if (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE) then !Water
              !reading the hydraulic heads file 
              call ReadHydraulicHeadFromFile
              !finding the multiplier specific using the overall real time
              call GetHydraulicHeadForTimeStep
          end if
          
            
         if (CalParams%Multipliers%GravityLoadType == LOAD_TYPE_LINEAR) then !Gravity
            Increment = (CalParams%Multipliers%GravityLoadStepIncrement) * &
                CalParams%TotalRealTime / CalParams%TotalTime
            
            CalParams%Multipliers%GravityCurrent = CalParams%Multipliers%GravityPrevious + Increment  
          end if
         if (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_LINEAR) then !Velocity
            Increment = (CalParams%Multipliers%VelocitySolidLoadStepIncrement) * &
                CalParams%TotalRealTime / CalParams%TotalTime
            
            CalParams%Multipliers%VelocitySolidCurrent = CalParams%Multipliers%VelocitySolidPrevious + Increment  
            CalParams%Multipliers%AccelerationSolid = CalParams%Multipliers%VelocitySolidLoadStepIncrement/CalParams%TotalTime
          end if
          if (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_STEP) then !velocity
            if(CalParams%TimeStep==1) then !Velocity
              CalParams%Multipliers%AccelerationSolid = &
                    CalParams%Multipliers%VelocitySolidLoadStepIncrement/CalParams%TimeIncrement 
            else
              CalParams%Multipliers%AccelerationSolid = 0.0
            end if
          end if
          if (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_FILE) then !Velocity
              !reading the prescribed velocity file 
              call ReadPrescribedAccelerationOrVelocityFromFile
              !finding the multiplier specific using the overall real time
              call GetPrescribedVelocityOrAccelerationForTimeStep
          end if
          

          
          if (CalParams%RigidBody%IsRigidBody) then
            CalParams%RigidBody%CurrentTraction = CalParams%RigidBody%TractionForce * CalParams%Multipliers%SolidACurrent(1)
          end if
                  
        end subroutine UpdateMultipliersForTimeDepencency




subroutine ReadPrescribedAccelerationOrVelocityFromFile()
        !**********************************************************************
        !
        ! Function:  Read the time history to provide a matrix of the multipliers. 
        !            This function is to be used to read both acceleration and/or velocity time history. 
        !            Function currently written for 1D prescribed velocity or acceleration from file.
        !            
        !
        ! File extension: PVF for prescribed velocity from file 
        !                 PAF for prescribed acceleration from file
        !
        ! File organisation:
        !            Column 1 - time (seconds)
        !            Column 2 - prescribed velocity or acceleration (dimensionless)
        !
        !**********************************************************************        
        implicit none
     
        ! local variables 
        integer(INTEGER_TYPE) :: ios, ii!, PVFunit
        character(len=255) :: PVFileName
        real(REAL_TYPE), dimension(1,2) :: Temp ! temporary variable for reading number of lines
        
        !PVFunit = 70

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! START: READING VELOCITY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        ! Open Prescribed Velocity File 
        if ((CalParams%PrescribedVelo%ApplyPrescribedVelo == .true.) &
            .and. (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_FILE)) then 
            
            ! Prepare file name for opening
            PVFileName = trim(CalParams%PrescribedVelo%FileNamePrescribedVelocity)//PRESCRIBEDVELOCITY_FILE_EXTENSION ! Concatenate prescribed velocity file extention (.PVF)
        
            

            if ((.not. allocated(PrescribedVelocityFileMatrix))) then
                ! Check if PVF file exists in project folder, otherwise give error and stop execution
                if ( FExist(trim(PVFileName)) ) then
                    call GiveMessage('Reading prescribed velocity file: ' // trim(PVFileName) )                     
                    call FileOpen(PVFunit,trim(PVFileName))
                else    
                    ! Give warning if file is not available
                    call GiveError('Prescribed velocity file does not exist!' // NEW_LINE('A') // 'required PVF file: ' // trim(PVFileName) )      
                end if
                
                ! Read Prescribed Velocity File in the matrix line by line
                do
                    read(PVFunit,*,iostat=ios) Temp ! Reading each row
                    if (ios/=0) then 
                        exit ! ios becomes non-zero when end of file was reached
                    end if 
                    CalParams%PrescribedVelo%PrescribedVelocityFileNumberOfLines = &
                        CalParams%PrescribedVelo%PrescribedVelocityFileNumberOfLines + 1 ! Counting the number of lines 
                end do
            
                ! Rewind file for reading in allocated matrix 
                rewind(PVFunit) 
            
                ! Allocate and read the matrix of velocity time history 
                allocate (PrescribedVelocityFileMatrix(CalParams%PrescribedVelo%PrescribedVelocityFileNumberOfLines,2)) ! Allocate the size of the matrix (specifying 2 columns)

                do ii = 1,CalParams%PrescribedVelo%PrescribedVelocityFileNumberOfLines
                    read(PVFunit,*,iostat=ios) (PrescribedVelocityFileMatrix(ii,1),PrescribedVelocityFileMatrix(ii,2))
                end do 
                close(PVFunit)
            end if 
        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END: READING VELOCITY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            


        end subroutine ReadPrescribedAccelerationOrVelocityFromFile
        
        
        subroutine GetPrescribedVelocityOrAccelerationForTimeStep()
        !**********************************************************************
        !
        ! Function:  obtains the prescribed veloctiy/acceleration from file for the current time step (at real time)
        !            using linear interpolation. 
        !
        !**********************************************************************     
         
        implicit none 
         
        !local variables 
        integer(INTEGER_TYPE) :: ii
        real(REAL_TYPE) :: Increment ! 1D which would become a matrix if extended to 2D or 3D file prescribed acceleration or velocity
        real(REAL_TYPE), dimension(2,2) :: CoordinatesForInterpolation ! Hardcoded (needs to be generalized for different dimensions e.g. 3D shaking)
         
        
        !!!!!!!! START: GET PRESCRIBED VELOCITY FROM FILE FOR OVERALL REAL TIME !!!!!!!!
 
        ! Prescribed velocity for the particular timestep
        if ((CalParams%PrescribedVelo%ApplyPrescribedVelo == .true.) &
            .and. (CalParams%Multipliers%VelocitySolidLoadType == LOAD_TYPE_FILE)) then
        
            ! Record the start time of the application of the transient boundary condition
            if( .not. allocated(CalParams%PrescribedVelo%TransientBoundaryConditionStaringTime)) then
                allocate(CalParams%PrescribedVelo%TransientBoundaryConditionStaringTime(1)) ! Allocate the size of the matrix (specifying 2 columns)
                CalParams%PrescribedVelo%TransientBoundaryConditionStaringTime = CalParams%OverallRealTime
                
                do ii=1,CalParams%PrescribedVelo%PrescribedVelocityFileNumberOfLines
                    !Correcting the starting time of the data read from file (hardcoded to 1D prescribed boundary condition) 
                    PrescribedVelocityFileMatrix(ii,1) = &
                        PrescribedVelocityFileMatrix(ii,1) + CalParams%PrescribedVelo%TransientBoundaryConditionStaringTime(1) - CalParams%TimeIncrement
                end do 
                
            end if
        
            ! Looks into the matrix first column (time) and finds the first point higher than the total real time
            do ii=1,CalParams%PrescribedVelo%PrescribedVelocityFileNumberOfLines !loop accross the rows from the presctibed velocity file
            
                if (CalParams%OverallRealTime > &
                    PrescribedVelocityFileMatrix(CalParams%PrescribedVelo%PrescribedVelocityFileNumberOfLines,1)) then
                    CalParams%Multipliers%VelocitySolidCurrent = 0 ! Stop prescribed velocity
                    CalParams%Multipliers%AccelerationSolid = 0
         
                else if (PrescribedVelocityFileMatrix(ii,1)>CalParams%OverallRealTime) then 
                    CoordinatesForInterpolation(1,:) = PrescribedVelocityFileMatrix(ii-1,:) ! Store from file: time (1st column) and velocity multiplier (2nd column) in previous line 
                    CoordinatesForInterpolation(2,:) = PrescribedVelocityFileMatrix(ii,:) ! Store from file: time(1st column) and velocity multiplier (2nd column) in current line
                 
                    ! Linear interpolation 
                    Increment = (((CoordinatesForInterpolation(2,2)-CoordinatesForInterpolation(1,2))&
                        /(CoordinatesForInterpolation(2,1)-CoordinatesForInterpolation(1,1))) * &
                        (CalParams%OverallRealTime - CoordinatesForInterpolation(1,1)))
         
                    ! For the first timestep VelocitySolidPrevious does not exist. 
                    ! Hence, special consideration is made for the first step by using the reading before the timestep occurs
                    if ((CalParams%TimeStep == 1) .and. (CalParams%IStep == 1)) then ! If it is the first timestep in the first load step (i.e. no CalParams%Multipliers%VelocitySolidPrevious)
                        CalParams%Multipliers%VelocitySolidPrevious = CoordinatesForInterpolation(1,2)
                    else ! If it is after the first step 
                        CalParams%Multipliers%VelocitySolidPrevious = CalParams%Multipliers%VelocitySolidCurrent
                    end if     
                                        
                    CalParams%Multipliers%VelocitySolidCurrent = CoordinatesForInterpolation(1,2) + Increment
                    CalParams%Multipliers%AccelerationSolid = (CalParams%Multipliers%VelocitySolidCurrent- &
                                                              CalParams%Multipliers%VelocitySolidPrevious)&
                                                              /CalParams%TimeIncrement 
                    
                    exit 
                end if 
            end do
            
        
    
        else         
                CalParams%Multipliers%VelocitySolidCurrent = 0
                !CalParams%Multipliers%AccelerationSolid = 0
        end if

        !!!!!!!! END: GET PRESCRIBED VELOCITY FROM FILE FOR OVERALL REAL TIME !!!!!!!!    
        end subroutine GetPrescribedVelocityOrAccelerationForTimeStep


        subroutine ReadHydraulicHeadFromFile()
        !**********************************************************************
        !
        ! Function:  Read the time history to provide a matrix of the multipliers. 
        !            This function is to be used to read hydraulic head time history. 
        !            
        !
        ! File extension: HHBF for Hydraulic Head from file 
        !                 
        !
        ! File organisation:
        !            Column 1 - time (seconds)
        !            Column 2 - hydraulic head values (m) !!!!Later necessary check that reference for hydraulic heads is equal to the one for boundary node location
        !
        !**********************************************************************        
        implicit none
        
            
        ! local variables 
        integer(INTEGER_TYPE) :: ios, ii, IError
        character(len=MAX_FILENAME_LENGTH) :: FileName
        real(REAL_TYPE), dimension(1,2) :: Temp ! temporary variable for reading number of lines
        integer(INTEGER_TYPE) :: NumberOfLines
        

        
        if ((CalParams%PrescribedHead%HydraulicHead == .false.) .and. (.not.(CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE))) RETURN
        
        ! Prepare file name for opening
        FileName = trim(CalParams%FileNames%ProjectName)//HYDRAULICHEAD_FILE_EXTENSION ! 
        if (FExist(FileName)) call FileOpen(HHBFUnit,trim(FileName))
        
        if (.not. allocated(HydraulicHeadFileMatrix)) then
    
          ! Read hydraulic head File in the matrix line by line
          NumberOfLines = 0  
          do
            read (HHBFUnit, *, iostat=ios) Temp! Reading each row
            if (ios/=0) exit ! ios becomes non-zero when end of file was reached
            NumberOfLines = NumberOfLines + 1 ! Counting the number of lines 
            
          end do
          CalParams%PrescribedHead%HydraulicHeadFileNumberOfLines = NumberOfLines
        
          ! Rewind file for reading in allocated matrix 
          rewind(HHBFUnit) 
            
          ! Allocate and read the matrix 
          allocate (HydraulicHeadFileMatrix(NumberOfLines,2), stat = IError) !allocate the size of the matrix (specifying 2 columns for 1D shaking)
          HydraulicHeadFileMatrix = 0.0
          do ii = 1, NumberOfLines
             read(HHBFUnit, *, iostat=ios) (HydraulicHeadFileMatrix(ii,1), HydraulicHeadFileMatrix(ii,2))
          end do 
          close(HHBFUnit)
        end if      

  
        end subroutine ReadHydraulicHeadFromFile
        
        
        
        subroutine GetHydraulicHeadForTimeStep()
        !**********************************************************************
        !
        ! Function:  obtains the HYDRAULIC HEAD  from file for the current time step (at real time)
        !            using linear interpolation. 
        !
        !**********************************************************************     
         
        implicit none 
         
        !local variables 
        integer(INTEGER_TYPE) :: ii
        real(REAL_TYPE) :: Increment ! 1D which would become a matrix if extended to 2D or 3D file prescribed acceleration or velocity
        real(REAL_TYPE), dimension(2,2) :: CoordinatesForInterpolation ! Hardcoded (needs to be generalized for different dimensions e.g. 3D shaking)
        real(REAL_TYPE) :: LastRow
         
        
        LastRow = 0.0
        !!!!!!!! START: GET HYDRAULIC HEAD FROM FILE FOR OVERALL REAL TIME !!!!!!!!
        
        ! Hydaulic head for the particular timestep
        if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
             
        
        LastRow = CalParams%PrescribedHead%HydraulicHeadFileNumberOfLines
        
            ! Looks into the matrix first column (time) and finds the first point higher than the total real time
            do ii=1,CalParams%PrescribedHead%HydraulicHeadFileNumberOfLines  ! Loop accross the rows from the hydraulic head file
            
                if (CalParams%OverallRealTime>HydraulicHeadFileMatrix(CalParams%PrescribedHead%HydraulicHeadFileNumberOfLines,1)) then
                CalParams%Multipliers%HydraulicHeadCurrent = HydraulicHeadFileMatrix(LastRow,2)
                                    
                else if (HydraulicHeadFileMatrix(ii,1)>CalParams%OverallRealTime) then 
                    CoordinatesForInterpolation(1,:) = HydraulicHeadFileMatrix(ii-1,:) ! Store from file: time (1st column) and HYDR HEAD (2nd column) in previous line 
                    CoordinatesForInterpolation(2,:) = HydraulicHeadFileMatrix(ii,:) ! Store from file: time (1st column) and HYDR HEAD (2nd column) in current line
                                
                    
                    ! Linear interpolation 
                    Increment = (((CoordinatesForInterpolation(2,2)-CoordinatesForInterpolation(1,2))/(CoordinatesForInterpolation(2,1)-CoordinatesForInterpolation(1,1))) * &
                        (CalParams%OverallRealTime - CoordinatesForInterpolation(1,1)))
                 
                    CalParams%Multipliers%HydraulicHeadCurrent = CoordinatesForInterpolation(1,2) + Increment
                    exit 
                end if 
            end do
        end if 
        end subroutine GetHydraulicHeadForTimeStep
        
        
        subroutine WriteCalculationParameters() 
        !**********************************************************************
        !
        ! Function:  Write the CPS file at the end of each load phase for the next step
        !
        !**********************************************************************
                
        implicit none
      
          ! local variables
          character(len=MAX_FILENAME_LENGTH) :: FileName, FileVersion 
          integer(INTEGER_TYPE) :: FileUnit
 
          FileName = trim( CalParams%FileNames%ProjectName) // CPS_FILE_EXTENSION // trim(CalParams%FileNames%LoadStepExt )
          FileUnit = CPSUnit
          FileVersion = CalParams%CPSversion
          
          call FileOpen(FileUnit, FileName) ! open the CPS file for the next load step

          ! write CPS data depending on file version
          select case (FileVersion)
          
          case(Anura3D_v2023)  
            call WriteCPS(FileUnit,FileVersion)
          case(Anura3D_v2022)  
            call WriteCPS(FileUnit,FileVersion)
          case(Anura3D_v2021)  
            call WriteCPS(FileUnit,FileVersion)
          case(Anura3D_v2019_2)  
            call WriteCPS(FileUnit,FileVersion) 
          case default
            call GiveError('Writing CPS file: Version not properly defined.' // NEW_LINE('A') // 'CPS version: ' // trim(FileVersion) )
          end select  

          close(FileUnit) ! close the CPS file for the next load step
        
        end subroutine WriteCalculationParameters
        
        
        subroutine WriteCPS(FileUnit,FileVersion)
        !**********************************************************************
        !
        ! Function:  Write the CPS file at the end of each load phase for the next step
        !
        !            CPS version 2022, 2021, 2019.2 
        !   
        !**********************************************************************
   
        implicit none
         
          character(len=MAX_FILENAME_LENGTH), intent(in) :: FileVersion
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: I
          
          select case (FileVersion)
              case(Anura3D_v2023)  
                  write(FileUnit, '(A)') Anura3D_v2023
              case(Anura3D_v2022)  
                  write(FileUnit, '(A)') Anura3D_v2022
              case(Anura3D_v2021)  
                  write(FileUnit, '(A)') Anura3D_v2021
              case(Anura3D_v2019_2)  
                  write(FileUnit, '(A)') Anura3D_v2021
          end select
          
          ! TIME AND LOAD STEP DATA
          call WriteInteger(FileUnit, '$$NUMBER_OF_LOADSTEPS', [CalParams%NLoadSteps], 1, .true.)
          call WriteReal(FileUnit, '$$TIME_PER_LOADSTEP', [CalParams%TotalTime], 1, .true.)
          call WriteReal(FileUnit, '$$TOTAL_TIME', [CalParams%OverallRealTime], 1, .true.)
          call WriteReal(FileUnit, '$$COURANT_NUMBER', [CalParams%CourantNumber], 1, .true.)

          ! CALCULATION DATA
          select case(CalParams%ComputationMethod)
            case(MPM_MIXED_INTEGRATION) 
              call WriteString(FileUnit, '$$COMPUTATION_METHOD', MPM_MIXED_INTEGRATION_SPECIFIER, 1, .true.)
            case(MPM_MP_INTEGRATION) 
              call WriteString(FileUnit, '$$COMPUTATION_METHOD', MPM_MP_INTEGRATION_SPECIFIER, 1, .true.)
            case(FEM) 
              call WriteString(FileUnit, '$$COMPUTATION_METHOD', FEM_SPECIFIER, 1, .true.)
            case(UL_FEM) 
              call WriteString(FileUnit, '$$COMPUTATION_METHOD', UL_FEM_SPECIFIER, 1, .true.)
          end select

          ! GRAVITY DATA
          call WriteReal(FileUnit, '$$GRAVITY_ACCELERATION', [CalParams%GravityData%GAccel], 1, .true.)
          call WriteReal(FileUnit, '$$GRAVITY_VECTOR', [CalParams%GravityData%GravityVector(1:NDIM)], NDIM, .true.)

          ! LOAD/VELOCITY MULTIPLIERS AND TPYE
          call WriteStringReal(FileUnit, '$$GRAVITY_LOAD', CalParams%Multipliers%GravityLoadType, [CalParams%Multipliers%GravityCurrent, CalParams%Multipliers%GravityFinal], .true.)
        
          call WriteStringReal(FileUnit, '$$SOLID_TRACTION', CalParams%Multipliers%SolidALoadType(1), [CalParams%Multipliers%SolidACurrent(1), CalParams%Multipliers%SolidAFinal(1)], .true.)
          call WriteStringReal(FileUnit, '$$SOLID_TRACTION_B', CalParams%Multipliers%SolidALoadType(2), [CalParams%Multipliers%SolidACurrent(2), CalParams%Multipliers%SolidAFinal(2)], .true.)
          call WriteStringReal(FileUnit, '$$LIQUID_PRESSURE', CalParams%Multipliers%WaterALoadType(1), [CalParams%Multipliers%WaterACurrent(1), CalParams%Multipliers%WaterAFinal(1)], .true.)
		  call WriteStringReal(FileUnit, '$$LIQUID_PRESSURE_B', CalParams%Multipliers%WaterALoadType(2), [CalParams%Multipliers%WaterACurrent(2), CalParams%Multipliers%WaterAFinal(2)], .true.)
          call WriteStringReal(FileUnit, '$$GAS_PRESSURE', CalParams%Multipliers%GasALoadType(1), [CalParams%Multipliers%GasACurrent(1), CalParams%Multipliers%GasAFinal(1)], .false.)
          call WriteStringReal(FileUnit, '$$GAS_PRESSURE_B', CalParams%Multipliers%GasALoadType(2), [CalParams%Multipliers%GasACurrent(2), CalParams%Multipliers%GasAFinal(2)], .false.)
		   call WriteStringReal(FileUnit, '$$HYDRAULIC_HEAD', CalParams%Multipliers%HydraulicHeadType, [CalParams%Multipliers%HydraulicHeadCurrent], .true.)
          call WriteLogic(FileUnit, '$$APPLY_SEEPAGE_FACE', [CalParams%BoundaryConditions%UseSeepageFace], 1, .true.) 
          call WriteLogic(FileUnit, '$$APPLY_INFILTRATION', [CalParams%BoundaryConditions%UseInfiltration], 1, .true.) 
          call WriteStringReal(FileUnit, '$$PRESCRIBED_VELOCITY', CalParams%Multipliers%VelocitySolidLoadType, [CalParams%Multipliers%VelocitySolidCurrent, CalParams%Multipliers%VelocitySolidFinal], .true.)

          ! QUASI-STATIC CALCULATION
          call WriteLogic(FileUnit, '$$QUASISTATIC_CONVERGENCE', [CalParams%ApplyQuasiStatic], 1, .true.)
          call WriteReal(FileUnit, '$$TOLERATED_ERROR_SOLID', [CalParams%ToleratedErrorEnergy, CalParams%ToleratedErrorForce], 2, .true.)
          call WriteReal(FileUnit, '$$TOLERATED_ERROR_LIQUID', [CalParams%ToleratedErrorEnergyWater, CalParams%ToleratedErrorForceWater], 2, .true.)
          call WriteInteger(FileUnit, '$$MAXIMUM_TIME_STEPS', [CalParams%MaxTimeSteps], 1, .true.)

          ! IMPLICIT INTEGRATION QUASI-STATIC
          call WriteLogic(FileUnit, '$$APPLY_IMPLICIT_QUASISTATIC', [CalParams%ApplyImplicitQuasiStatic], 1, .true.) 
          call WriteInteger(FileUnit, '$$QS_MIN_ITER', [CalParams%ImplicitIntegration%MinDesiredIterations], 1, .true.) 
          call WriteInteger(FileUnit, '$$QS_MAX_ITER', [CalParams%ImplicitIntegration%MaxDesiredIterations], 1, .true.) 
          call WriteReal(FileUnit, '$$QS_UPSCALE_FACTOR', [CalParams%ImplicitIntegration%UpScaleFactor], 1, .true.) 
          call WriteReal(FileUnit, '$$QS_DOWNSCALE_FACTOR', [CalParams%ImplicitIntegration%DownScaleFactor], 1, .true.) 
          call WriteInteger(FileUnit, '$$QS_MAX_ITERATIONS', [CalParams%ImplicitIntegration%MaxIterations], 1, .true.) 
          call WriteLogic(FileUnit, '$$APPLY_ZLS', [CalParams%ImplicitIntegration%DoUseZLS], 1, .true.) 
          call WriteLogic(FileUnit, '$$APPLY_ARC_LENGTH_CONTROL', [CalParams%ImplicitIntegration%DoUseArcLengthControl], 1, .true.) 
          call WriteLogic(FileUnit, '$$APPLY_ENHANCED_STIFFNESS', [CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix], 1, .true.) 
          call WriteReal(FileUnit, '$$QS_STIFFNESS_INCREASE_FACTOR', [CalParams%ImplicitIntegration%StiffnessIncreaseFactor], 1, .true.) 
          call WriteReal(FileUnit, '$$QS_DANGLING_ELEMENT_FACTOR', [CalParams%ImplicitIntegration%DanglingElementFactor], 1, .true.) 
          call WriteLogic(FileUnit, '$$APPLY_STEP_EXTRAPOLATION', [CalParams%ImplicitIntegration%DoUseStepExtrapolation], 1, .true.) 
          call WriteLogic(FileUnit, '$$APPLY_JACOBIAN_CHECK', [CalParams%ImplicitIntegration%DoCheckChangeJacobian], 1, .true.) 
          call WriteReal(FileUnit, '$$LIMIT_JACOBIAN_CHANGE', [CalParams%ImplicitIntegration%LimitJacobianChange], 1, .true.) 
          call WriteLogic(FileUnit, '$$APPLY_AUTOMATIC_LOADSTEPPING', [CalParams%ImplicitIntegration%DoUseAutomaticLoadStepping], 1, .true.) 
          call WriteReal(FileUnit, '$$INCR_MULTIPLIER_SOILA', [CalParams%Multipliers%SolidAIncrement], 1, .true.) 
          call WriteReal(FileUnit, '$$INCR_MULTIPLIER_GRAVITY', [CalParams%Multipliers%GravityIncrement], 1, .true.) 
          call WriteReal(FileUnit, '$$INCR_MULTIPLIER_WATERA', [CalParams%Multipliers%WaterAIncrement], 1, .true.) 
          call WriteReal(FileUnit, '$$INCR_MULTIPLIER_GASA', [CalParams%Multipliers%GasAIncrement], 1, .true.) 


          ! MASS SCALING
          call WriteLogicReal(FileUnit, '$$MASS_SCALING', [CalParams%ApplyMassScaling], [CalParams%ScalingMassFactor], 1,.true.)
 
          ! DAMPING
          call WriteLogicReal(FileUnit, '$$HOMOGENEOUS_LOCAL_DAMPING', [CalParams%ApplyCPSDamping], [CalParams%DampingFactor], 1, .true.)
          
          ! BULK VISCOSITY DAMPING
          call WriteLogicReal(FileUnit, '$$BULK_VISCOSITY_DAMPING', [CalParams%ApplyBulkViscosityDamping], [CalParams%BulkViscosityDamping1, CalParams%BulkViscosityDamping2], 2, .true.)

          ! DIVERGENCE CHECK
          call WriteLogic(FileUnit, '$$APPLY_DIVERGENCE_CHECK', [CalParams%ConvergenceCheck%ApplyDivergenceCheck], 1, .true.)  
          
          ! OPTIMIZATION: Speeding up
          call WriteReal(FileUnit, '$$MINIMUM_DETERMINANT_RATIO', [CalParams%MinimumDeterminantRatio], 1, .true.)

          ! PARALLELIZATION
          call WriteInteger(FileUnit, '$$NUMBER_OF_THREADS', [CalParams%NThreads], 1, .true.)
            
          ! ABSORBING BOUNDARIES

          call WriteLogic(FileUnit, '$$APPLY_NODE_DATA_AB', [CalParams%AbsorbingBoundaries%ApplyNodeData], 1, .true.)
          
          ! CONTACT ALGORITHM
          call WriteLogic(FileUnit, '$$CONTACT_FORMULATION', [CalParams%ApplyContactAlgorithm], 1, .true.)
          
          call WriteLogic(FileUnit, '$$APPLY_CONTACT_MESH_BOUNDARY', [CalParams%ApplyContactMeshBoundary], 1, .true.)
          call WriteLogic(FileUnit, '$$APPLY_TRACTION_CONTACT', [CalParams%ApplyTractionContact], 1, .true.) 
          call WriteLogic(FileUnit, '$$MONOPILE_SLICE_CONTACT', [CalParams%MonoPileSliceContact], 1, .true.)
          call WriteLogic(FileUnit, '$$APPLY_CONTACT_VELOCITY_SCALING', [CalParams%ApplyContactVelocityScaling], 1, .true.)
		  call WriteReal(FileUnit, '$$VELOCITY_SCALING_FACTOR', [CalParams%ContactScalingFactor], 1, .true.)
		  call WriteLogic(FileUnit, '$$APPLY_CONTACT_NORMAL_CORRECTION', [CalParams%ApplyContactNormalCorrection], 1, .true.)
		  call WriteInteger(FileUnit, '$$NUMBER_OF_CORRECTED_NORMALS', [CalParams%NCorrectedNormals], 1, .true.)
		   if ((CalParams%NCorrectedNormals>0)) then
              write(FileUnit, '(A)') '$$NODE_NORMAL_DATA'
              do I = 1, CalParams%NumberSoilLayers
                write(FileUnit, *) CalParams%ManuallyDefinedVectors(I,1:(1+NVECTOR))
              end do
            end if 
          ! APPLY EXCAVATION
          call WriteInteger(FileUnit, '$$EXCAVATION', [CalParams%NumberExcavatedVolumes], 1, .true.) 
          do I = 1, CalParams%NumberExcavatedVolumes
            call WriteInteger(FileUnit, '', [CalParams%ExcavatedVolumeID(I),CalParams%ExcavationStages(I,1:2)], 3, .true.) 
          end do 
          
          ! FIXED SOLID SKELETON
          call WriteLogic(FileUnit, '$$FIX_SOLID_SKELETON', [CalParams%ApplyFixedSolidSkeleton], 1, .true.) 

          ! STRAIN SMOOTHING
          call WriteLogic(FileUnit, '$$STRAIN_SMOOTHING', [CalParams%ApplyStrainSmoothing], 1, .true.)
          
          ! DOUBLE POINT :: STRAIN SMOOTHING LIQUID 
          call WriteLogic(FileUnit, '$$APPLY_STRAIN_SMOOTHING_LIQUID_TWOLAYERFORM', [CalParams%ApplyStrainSmoothingLiquidTwoLayer], 1, .true.)
          
          ! DOUBLE POINT :: LIQUID PRESSURE CAVITATION THRESHOLD
          call WriteReal(FileUnit, '$$LIQUID_CAVITATION', [CalParams%LiquidPressureCavitationThreshold], 1, .true.) 
          
          ! DOUBLE POINT :: APPLY_ERGUN_LAW
          call WriteLogic(FileUnit, '$$UPDATE_PERMEABILITY_DARCY_ERGUN', [CalParams%TwoLayerApplyUpdatePermeability, CalParams%TwoLayerApplyErgunLaw], 2, .true.)
          call WriteReal(FileUnit, '$$GRAIN_SIZE_DIAMETER', [CalParams%TwoLayerErgunLawDiameterPartic, CalParams%TwoLayerErgunLawDiameterPartic2], 2, .true.) 
          call WriteReal(FileUnit, '$$ERGUN_CONSTANTS', [CalParams%ERGUNCONSTANTA, CalParams%ERGUNCONSTANTB], 2, .true.)
            
          ! DOUBLE POINT :: NO TENSILE STRESS IN LIQUID MATERIAL POINT WITH LIQUID STATUS
          call WriteLogic(FileUnit, '$$NO_TENSILE_STRESS_LIQUID_MP_WITH_LIQUID_STATUS', [CalParams%TwoLayerApplyNoTensStressLiqMPwLiqStatus], 1, .true.)

          ! DOUBLE POINT :: MAXIMUM POROSITY
          call WriteReal(FileUnit, '$$MAXIMUM_POROSITY', [CalParams%LimitPorosity], 1, .true.)
            
          ! DOUBLE POINT :: DETECT FREE SURFACE
          call WriteLogicReal(FileUnit, '$$DETECT_FREE_SURFACE', [CalParams%ApplyDetectLiquidFreeSurface], [CalParams%FreeSurfaceFactor], 1, .true.)
          
          call WriteLogic(FileUnit, '$$RESET_DISPLACEMENTS', [.false.], 1, .true.) 
         
          call WriteLogic(FileUnit, '$$REMOVE_FIXITIES', [CalParams%ApplyRemoveSolidFixities, CalParams%ApplyRemoveLiquidFixities, CalParams%ApplyRemoveGasFixities], 3, .true.)
            
          ! K0-PROCEDURE
          call WriteLogic(FileUnit, '$$K0_PROCEDURE', [CalParams%ApplyK0Procedure], 1, .true.)
          call WriteReal(FileUnit, '$$K0_MAX_SUCTION', [CalParams%K0MaxSurfaceSuction], 1, .true.)
          call WriteReal(FileUnit, '$$SURFACE_ELEVATION', [CalParams%SoilSurfacePoint(1:2)], 2, .true.)
          
          call WriteInteger(FileUnit, '$$NUMBER_SOIL_LAYERS', [CalParams%NumberSoilLayers], 1, .true.) 
            if (CalParams%NumberSoilLayers>0) then
              write(FileUnit, '(A)') '$$THICKNESS_SOIL_LAYERS'
              do I = 1, CalParams%NumberSoilLayers
                write(FileUnit, *) CalParams%ThicknessSoilLayer(I)
              end do
            end if
          
          call WriteReal(FileUnit, '$$LIQUID_SURFACE', [CalParams%LiquidSurfacePoint(1:2)], 2, .true.) 
          call WriteReal(FileUnit, '$$CONSIDERED_LOAD_K0', [CalParams%InitialVerticalLoadK0], 1, .true.) 
          call WriteReal(FileUnit, '$$INITIAL_WATER_PRESSURE', [CalParams%InitialWaterPressure], 1, .true.)
          
          ! MATERIAL DATA
          call WriteLogic(FileUnit, '$$APPLY_MATERIAL_UPDATE', [.false.], 1, .true.)
          call WriteLogic(FileUnit, '$$APPLY_POROSITY_UPDATE', [CalParams%ApplyPorosityUpdate], 1, .true.)
                   
          call WriteLogicInteger(FileUnit, '$$SUBMERGED_CALCULATION', [CalParams%ApplySubmergedCalculation], [CalParams%NumberSubmergedCalculation], 1, .true.)
                   
          ! INITIAL VELOCITY ON MATERIAL POINTS
		  call WriteLogic(FileUnit, '$$INITIAL_VELOCITY', [.false.], 1, .true.)
		  
		  call WriteLogic(FileUnit, '$$APPLY_OBJECTIVE_STRESS', [CalParams%ApplyObjectiveStress], 1, .true.)
          
          call WriteReal(FileUnit, '$$DEGREE_OF_FILLING', [CalParams%RequiredDegreeOfFilling], 1, .true.)
            
          call WriteInteger(FileUnit, '$$NUMBER_OF_ACTIVE_ELEMENTS', [Counters%NAEl], 1, .true.)
          
          ! EMPTY ELEMENTS
          call WriteLogic(FileUnit, '$$APPLY_EMPTY_ELEMENTS', [CalParams%ApplyEmptyElements], 1, .true.)
          call WriteReal(FileUnit, '$$MASS_FACTOR_VIRTUAL_PARTICLES', [CalParams%VirtualParticleMassFactor], 1, .true.)
          call WriteInteger(FileUnit, '$$GROUP_THRESHOLD', [CalParams%GroupThreshold], 1, .true.)
          
          ! VISUALIZATION
          call WriteString(FileUnit, '$$VISUALIZATION_OPTION', CalParams%Visualization, 1, .true.)
 
          ! OUTPUT
          call WriteInteger(FileUnit, '$$OUTPUT_NUMBER_OF_MATERIAL_POINTS', [CalParams%OutputNumberParticles], 1, .true.)
          if (CalParams%OutputNumberParticles > 0) then
            call WriteInteger(FileUnit, '$$OUTPUT_MATERIAL_POINTS', [CalParams%OutputParticles(1)], 1, .true.)
            do I = 2, CalParams%OutputNumberParticles
              call WriteInteger(FileUnit, '', [CalParams%OutputParticles(I)], 1, .true.) 
            end do
          else
            call WriteInteger(FileUnit, '$$OUTPUT_MATERIAL_POINTS', [-1], 1, .true.)
          end if

          call WriteLogic(FileUnit, '$$OUTPUT_BASIC_DATA', [CalParams%OutputBasicData], 1, .true.)
          call WriteInteger(FileUnit, '$$OUTPUT_CURVES_INTERVAL', [CalParams%OutputCurvesIntervals], 1, .true.)

          call WriteInteger(FileUnit, '$$OUTPUT_NUMBER_OF_NODES', [CalParams%OutputNumberNodes], 1, .true.) 
          call WriteInteger(FileUnit, '$$OUTPUT_NODES', [CalParams%OutputNodes(1:10)], 10, .true.)

          call WriteLogic(FileUnit, '$$OUTPUT_DEBUG_DATA', [CalParams%OutputDebugData], 1, .true.)
          call WriteLogic(FileUnit, '$$OUTPUT_EACH_LOAD_STEP', [CalParams%OutputEachLoadStep], 1, .true.)
          call WriteLogic(FileUnit, '$$QUICK_CHECK_OUTPUT', [CalParams%ApplyQuickCheckOutput], 1, .true.)

          call WriteInteger(FileUnit, '$$FILE_COUNTER', [CalParams%FileCounter], 1, .true.) 

          call WriteInteger(FileUnit, '$$PARTICLE_FILE_COUNTER', [CalParams%ParticleFileCounter], 1, .true.) 

          ! Feedback
          call WriteInteger(FileUnit, '$$FEEDBACK_LEVEL', [CalParams%FeedbackLevel], 1, .true.) 
            
          call WriteLogic(FileUnit, '$$APPLY_SMOOTHENING_LIQUID_PRESSURE_INCREMENT', [CalParams%ApplySmootheningLiquidPressureIncrement], 1, .true.) ! omit for external version   

          ! END OF CPS FILE
          write(FileUnit, '(A)') '$$END'
          
        end subroutine WriteCPS
        
        
        subroutine DetermineLoadStep()
        !**********************************************************************
        !
        !    Function:  Retrieves the last computed step from the GIP file
        !               by searching in the project directory for the GIP
        !               file with the highest step number, GIP_099 ....
        !               For a new computation the GIP file is generated for
        !               post-processing purposes.
        !
        !**********************************************************************
        
        implicit none

          ! Local variables
          character(len = MAX_FILENAME_LENGTH) :: FileName
          character(len = MAX_EXTENSION_LENGTH) :: StepString
          integer(INTEGER_TYPE) :: CheckedStep
       
          CheckedStep = 1
          CalParams%PreviouslyRealisedLoadStep = 0
          do
            call GetStepExt(CheckedStep, StepString)
            FileName = trim(CalParams%FileNames%ProjectName)//GIP_STEP_FILE_EXTENSION//trim(StepString)
            if (FExist(trim(FileName))) then
              CalParams%PreviouslyRealisedLoadStep = CheckedStep
              CheckedStep = CheckedStep + 1
            else
              EXIT
            end if
          end do

          call GetStepExt(CalParams%PreviouslyRealisedLoadStep, StepString)
          CalParams%FileNames%PreviousStepExt = trim(StepString)
          CalParams%IStep = CalParams%PreviouslyRealisedLoadStep + 1 ! IStep will be 1 for new computation
          call GetStepExt(CalParams%IStep, StepString)
          CalParams%FileNames%LoadStepExt = trim(StepString)

          if (.not.IsFollowUpPhase()) then
            call WriteGIPFile(DO_NOT_USE_STEP_EXTENSION)
          end if

        end subroutine DetermineLoadStep


        subroutine WriteGIPFile(DoUseStepExtension)
        !**********************************************************************
        !
        !    Function:  Writes the GIP file using the provided extension.
        !               If DoUseStepExtension is true, the GIP_XXX extension
        !               is used, otherwise the GIP extension.
        !
        !**********************************************************************
        
        implicit none
          
          ! Arguments
          logical, intent(in) :: DoUseStepExtension
          
          ! Local variables
          character(len = MAX_EXTENSION_LENGTH) :: FileExtension
          integer(INTEGER_TYPE) :: StartStep, EndStep
          
          if (DoUseStepExtension) then
            FileExtension = GIP_STEP_FILE_EXTENSION//CalParams%FileNames%LoadStepExt
          else
            FileExtension = GIP_FILE_EXTENSION
          end if
          
          call FileOpen(TMP_UNIT, trim(CalParams%FileNames%ProjectName)//trim(FileExtension))

          if (DoUseStepExtension) then
            StartStep = CalParams%PreviouslyRealisedLoadStep + 1
            EndStep = CalParams%IStep
          else
            StartStep = INITIAL_STEP
            EndStep = MAXIMUM_STEP
          end if

          write(TMP_UNIT, 101) '[Phase0]' ! Presently not used
          write(TMP_UNIT, 102) START_STEP_SPECIFIER, StartStep
          write(TMP_UNIT, 102) END_STEP_SPECIFIER, EndStep

          close(TMP_UNIT)

  101     format(A)
  102     format(A, I)
        
        end subroutine WriteGIPFile

        
        
        logical function IsLastTimeStepOfExplicitCalculation() result(res)
        !**********************************************************************
        !
        !    Function: Checks whether is the last time step of the explicit calculation or not.
        !    
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none

        res = ((CalParams%TotalTime - CalParams%TotalRealTime) < SMALL) .and. (.not.CalParams%ApplyQuasiStatic)

        end function IsLastTimeStepOfExplicitCalculation


        
        logical function IsMPMSkipConvection() result(res)
        !**********************************************************************
        !
        !    Function: Checks whether is the MPM convection phase can be skipped or not.
        !    
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none

        if (CalParams%SkipConvection) then
          res = .not.IsDistorted .and. .not.(     mod(CalParams%TimeStep, CalParams%NumberSkipConvection) == 0 &
                                             .or. CalParams%TimeStep==1 &
                                             .or. IsLastTimeStepOfExplicitCalculation())
        else
          res = .false.
        endif

        end function IsMPMSkipConvection

        
      
      
        subroutine WriteLogic(unit, name, value, ndim, IsWritten)
        !**********************************************************************
        !
        ! Function:  Write a logical (array) value into the CPS file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(*), intent(in) :: name
          integer(INTEGER_TYPE), intent(in) :: ndim
          logical, dimension(:), intent(in) :: value
          logical, intent(in) :: IsWritten
          
          ! local variables
          integer(INTEGER_TYPE) :: ios

          if (.not.IsWritten) RETURN
          
          if ( trim(name) /= '' ) then
            write(unit, '(A)', iostat=ios) trim(name)
            call Assert( ios == 0, 'CPS file: Can''t write flag ' // trim(name) // ' into CPS file.')
          else
            ! ...do nothing
          end if    

          select case (ndim)
            case(1) 
              if (value(1)) write(unit, '(A2)', iostat=ios) ' 1'
              if (.not.value(1)) write(unit, '(A2)', iostat=ios) ' 0'
            case(2)
              if (value(1).AND.value(2)) write(unit, '(A4)', iostat=ios) ' 1 1'
              if (value(1).AND..not.value(2)) write(unit, '(A4)', iostat=ios) ' 1 0'
              if (.not.value(1).AND.value(2)) write(unit, '(A4)', iostat=ios) ' 0 1'
              if (.not.value(1).AND..not.value(2)) write(unit, '(A4)', iostat=ios) ' 0 0'
            case(3)
              if (value(1).AND.value(2).AND.value(3)) write(unit, '(A6)', iostat=ios) ' 1 1 1'
              if (value(1).AND.value(2).AND..not.value(3)) write(unit, '(A6)', iostat=ios) ' 1 1 0'
              if (value(1).AND..not.value(2).AND.value(3)) write(unit, '(A6)', iostat=ios) ' 1 0 1'
              if (value(1).AND..not.value(2).AND..not.value(3)) write(unit, '(A6)', iostat=ios) ' 1 0 0'
              if (.not.value(1).AND.value(2).AND.value(3)) write(unit, '(A6)', iostat=ios) ' 0 1 1'
              if (.not.value(1).AND.value(2).AND..not.value(3)) write(unit, '(A6)', iostat=ios) ' 0 1 0'
              if (.not.value(1).AND..not.value(2).AND.value(3)) write(unit, '(A6)', iostat=ios) ' 0 0 1'
              if (.not.value(1).AND..not.value(2).AND..not.value(3)) write(unit, '(A6)', iostat=ios) ' 0 0 0'
            case default
              call GiveError('CPS file: Writing logic value not possible. Wrong dimension.')
          end select
            
          call Assert( ios == 0, 'CPS file: Can''t write value of ' // trim(name) // ' into CPS file.')

        end subroutine WriteLogic

        
        subroutine WriteInteger(unit, name, value, ndim, IsWritten)
        !**********************************************************************
        !
        ! Function:  Write an integer (array) value to the CPS file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(*), intent(in) :: name
          integer(INTEGER_TYPE), intent(in) :: ndim
          integer(INTEGER_TYPE), dimension(:), intent(in) :: value
          logical, intent(in) :: IsWritten
          
          ! local variables
          integer(INTEGER_TYPE) :: ios
          
          if (.not.IsWritten) RETURN

          if ( trim(name) /= '' ) then
            write(unit, '(A)', iostat=ios) trim(name)
            call Assert( ios == 0, 'CPS file: Can''t write flag ' // trim(name) // ' into CPS file.')
          else
            ! ...do nothing
          end if    

          select case (ndim)
            case(1)
              write(unit, *, iostat=ios) value(1)
            case(2)
              write(unit, *, iostat=ios) value(1:2)
            case(3)
              write(unit, *, iostat=ios) value(1:3)
            case(10)  
              write(unit, *, iostat=ios) value(1:10)
            case default  
              call GiveError('CPS file: Writing integer value not possible. Wrong dimension.')
          end select
          
          call Assert( ios == 0, 'CPS file: Can''t write value of ' // trim(name) // ' into CPS file.')

        end subroutine WriteInteger      

        
        subroutine WriteReal(unit, name, value, ndim, IsWritten)
        !**********************************************************************
        !
        ! Function:  Write a real (array) value to the CPS file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(*), intent(in) :: name
          integer(INTEGER_TYPE), intent(in) :: ndim
          real(REAL_TYPE), dimension(:), intent(in) :: value
          logical, intent(in) :: IsWritten
          
          ! local variables
          integer(INTEGER_TYPE) :: ios
          
          if (.not.IsWritten) RETURN 
          
          if ( trim(name) /= '' ) then
            write(unit, '(A)', iostat=ios) trim(name)
            call Assert( ios == 0, 'CPS file: Can''t write flag ' // trim(name) // ' into CPS file.')
          else
            ! ...do nothing
          end if    
          
          write(unit, *, iostat=ios) value(1:ndim)
            
          call Assert( ios == 0, 'CPS file: Can''t write value of ' // trim(name) // ' into CPS file.')

        end subroutine WriteReal 

        
        subroutine WriteString(unit, name, value, ndim, IsWritten)
        !**********************************************************************
        !
        ! Function:  Write a real (array) value to the CPS file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(*), intent(in) :: name
          integer(INTEGER_TYPE), intent(in) :: ndim
          character(*), intent(in) :: value
          logical, intent(in) :: IsWritten
          
          ! local variables
          integer(INTEGER_TYPE) :: ios
          
          if (.not.IsWritten) RETURN 
          
          if ( trim(name) /= '' ) then
            write(unit, '(A)', iostat=ios) trim(name)
            call Assert( ios == 0, 'CPS file: Can''t write flag ' // trim(name) // ' into CPS file.')
          else
            ! ...do nothing
          end if    
          
          select case (ndim)
            case(1)
              write(unit, *, iostat=ios) value
            case default  
              call GiveError('CPS file: Writing string value not possible. Wrong dimension.')
          end select
          
          call Assert( ios == 0, 'CPS file: Can''t write value of ' // trim(name) // ' into CPS file.')

        end subroutine WriteString
        
        
        subroutine WriteStringReal(unit, name, valueS, valueR, IsWritten)
        !**********************************************************************
        !
        ! Function:  Write a combination of string(scalar) and real(array) value into the CPS file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(*), intent(in) :: name
          character(*), intent(in) :: valueS
          real(REAL_TYPE), dimension(:), intent(in) :: valueR
          logical, intent(in) :: IsWritten
          
          ! local variables
          integer(INTEGER_TYPE) :: ios

          if (.not.IsWritten) RETURN 
          
          if ( trim(name) /= '' ) then
            write(unit, '(A)', iostat=ios) trim(name)
            call Assert( ios == 0, 'CPS file: Can''t write flag ' // trim(name) // ' into CPS file.')
          else
            ! ...do nothing
          end if    

          write(unit, *, iostat=ios) valueS, valueR
          call Assert( ios == 0, 'CPS file: Can''t write value of ' // trim(name) // ' into CPS file.')

        end subroutine WriteStringReal        
        
        
        subroutine WriteLogicReal(unit, name, valueL, valueR, ndim, IsWritten)
        !**********************************************************************
        !
        ! Function:  Write a combination of logical and real value into the CPS file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(*), intent(in) :: name
          integer(INTEGER_TYPE), intent(in) :: ndim
          logical, dimension(:), intent(in) :: valueL
          real(REAL_TYPE), dimension(:), intent(in) :: valueR
          logical, intent(in) :: IsWritten
          
          ! local variables
          integer(INTEGER_TYPE) :: ios

          if (.not.IsWritten) RETURN 
          
          if ( trim(name) /= '' ) then
            write(unit, '(A)', iostat=ios) trim(name)
            call Assert( ios == 0, 'CPS file: Can''t write flag ' // trim(name) // ' into CPS file.')
          else
            ! ...do nothing
          end if    

          select case (ndim)
            case(1) 
              if (valueL(1)) write(unit, '(A2, F)', iostat=ios) ' 1', valueR(1)
              if (.not.valueL(1)) write(unit, '(A2, F)', iostat=ios) ' 0', valueR(1)
            case(2) 
              if (valueL(1)) write(unit, '(A2, 2F)', iostat=ios) ' 1', valueR(1), valueR(2)
              if (.not.valueL(1)) write(unit, '(A2, 2F)', iostat=ios) ' 0', valueR(1), valueR(2)
            case(3) 
              if (valueL(1)) write(unit, '(A2, 3F)', iostat=ios) ' 1', valueR(1), valueR(2), valueR(3)
              if (.not.valueL(1)) write(unit, '(A2, 3F)', iostat=ios) ' 0', valueR(1), valueR(2), valueR(3)
            case(4) 
              if (valueL(1)) write(unit, '(A2, 4F)', iostat=ios) ' 1', valueR(1), valueR(2), valueR(3), valueR(4)
              if (.not.valueL(1)) write(unit, '(A2, 4F)', iostat=ios) ' 0', valueR(1), valueR(2), valueR(3), valueR(4)
            case default
              call GiveError('CPS file: Writing logic and real value not possible. Wrong dimension.')
          end select
            
          call Assert( ios == 0, 'CPS file: Can''t write value of ' // trim(name) // ' into CPS file.')

        end subroutine WriteLogicReal

        
        subroutine WriteLogicInteger(unit, name, valueL, valueI, ndim, IsWritten)
        !**********************************************************************
        !
        ! Function:  Write a combination of logical and integer value into the CPS file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(*), intent(in) :: name
          integer(INTEGER_TYPE), intent(in) :: ndim
          logical, dimension(:), intent(in) :: valueL
          integer(INTEGER_TYPE), dimension(:), intent(in) :: valueI
          logical, intent(in) :: IsWritten
          
          ! local variables
          integer(INTEGER_TYPE) :: ios

          if (.not.IsWritten) RETURN 
          
          if ( trim(name) /= '' ) then
            write(unit, '(A)', iostat=ios) trim(name)
            call Assert( ios == 0, 'CPS file: Can''t write flag ' // trim(name) // ' into CPS file.')
          else
            ! ...do nothing
          end if    

          select case (ndim)
            case(1) 
              if (valueL(1)) write(unit, '(A2, I)', iostat=ios) ' 1', valueI(1)
              if (.not.valueL(1)) write(unit, '(A2, I)', iostat=ios) ' 0', valueI(1)
            case(2) 
              if (valueL(1)) write(unit, '(A2, 2I)', iostat=ios) ' 1', valueI(1), valueI(2)
              if (.not.valueL(1)) write(unit, '(A2, 2I)', iostat=ios) ' 0', valueI(1), valueI(2)
            case(3) 
              if (valueL(1)) write(unit, '(A2, 3I)', iostat=ios) ' 1', valueI(1), valueI(2), valueI(3)
              if (.not.valueL(1)) write(unit, '(A2, 3I)', iostat=ios) ' 0', valueI(1), valueI(2), valueI(3)
            case default
              call GiveError('CPS file: Writing logic and integer value not possible. Wrong dimension.')
          end select
            
          call Assert( ios == 0, 'CPS file: Can''t write value of ' // trim(name) // ' into CPS file.')

        end subroutine WriteLogicInteger                

              
      end module ModReadCalculationData

        
        
      
      integer(INTEGER_TYPE) function GetFeedbackLevel() result(res)
      !**********************************************************************
      !
      ! Function:  It reads CalParams%FeedbackLevel
      !
      !**********************************************************************
      use ModReadCalculationData, only: CalParams
      use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
      implicit none

      res = CalParams%FeedbackLevel

      end function GetFeedbackLevel

      
        