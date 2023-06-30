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
	  
	  
	  module ModMPMDYN2PhaseSP
      !**********************************************************************
      !
      !  Function : Contains the routine for forming the consolidation matrices
      !
      !     $Revision: 9434 $
      !     $Date: 2022-03-08 10:22:04 +0100 (mar, 08 mar 2022) $
      !
      !**********************************************************************
      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluationTETRA
      use ModElementEvaluationTRI
      use ModElementEvaluation
      use ModMPMData
      use ModMeshInfo
      use ModRotBoundCond
      use ModTwoLayerFormulation
   
      
      implicit none

        !> Lumped mass vector of water
        real(REAL_TYPE), dimension(:, :), allocatable :: LumpedMassWater
        !> Lumped mass vector of water * porosity       (or * porosity * DegreeSaturation, if is Unsaturated Calculation)
        real(REAL_TYPE), dimension(:, :), allocatable :: LumpedMassNWater
        !> Nodal load array of water gravity load
        real(REAL_TYPE), dimension(:, :), allocatable :: GravityLoadWater
        !> Nodal load array of water gravity load taking into accont porosity
        real(REAL_TYPE), dimension(:, :), allocatable :: GravityLoadWaterPorosity
        !> Nodal load array of mixture gravity load
        real(REAL_TYPE), dimension(:, :), allocatable :: GravityLoadMixture
        !> Conductivity matrix
        real(REAL_TYPE), dimension(:, :), allocatable :: ConductivityMatrix
        !> Conductivity matrix * porosity     (or * porosity * DegreeSaturation, if is Unsaturated Calculation)
        real(REAL_TYPE), dimension(:, :), allocatable :: ConductivityMatrixPorosity
        !> Nodal drag force
        real(REAL_TYPE), dimension(:, :), allocatable :: QVW
        !> Nodal drag force * Porosity
        real(REAL_TYPE), dimension(:, :), allocatable :: QVWPorosity
        !> Nodal drag force * Porosity at previous time step
        real(REAL_TYPE), dimension(:, :), allocatable :: QVWPorosityPrevious
        !> Internal load vector correspond to water
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadWater
        !> Internal load vector correspond to water at previous time step
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadWaterPrevious
        !> Internal load vector correspond to water at previous time step, correspond to Bishop load in partially saturated conditions
        real(REAL_TYPE), dimension(:, :), allocatable :: DummyIntLoadWaterPrevious
        !> Internal load vector correspond to water taking into account porosity
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadWaterPorosity
        !> Internal load vector correspond to water taking into account porosity at previous time step
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadWaterPorosityPrevious
        !> Nodal acceleration for water
        real(REAL_TYPE), dimension(:, :), allocatable :: AccelerationWater
        !> Nodal load array of external water load
        real(REAL_TYPE), dimension(:, :), allocatable :: ExtLoadWater
        !> Nodal load array of external water load taking into accont porosity
        real(REAL_TYPE), dimension(:, :), allocatable :: ExtLoadWaterPorosity
        !> Nodal load array of external water load (Total traction applied on nodes)              
        real(REAL_TYPE), dimension(:, :,:), allocatable :: ExtLoadWaterTotal
        real(REAL_TYPE), dimension(:, :), allocatable :: SpaceTimeExTLoadWater
        !real(REAL_TYPE), dimension(:, :), allocatable :: HydraulicHeadLoadTotal

        logical, dimension(:), allocatable :: IsInfiltrationNode !Determine if a node is on the infiltration surface
        !> Nodal load array of porosity*DegreeSaturation             
        real(REAL_TYPE), dimension(:, :), allocatable :: LumpedNodalPorosityDegSat
        logical, dimension(:), allocatable :: IsSeepageFaceNode !Determine if a node is on the seepage surface
        logical, dimension(:), allocatable :: IsHydraulicHeadNode !Determine if a node is on the seepage surface
      contains

        subroutine InitialiseTwoPhaseData()
        !**********************************************************************
        !
        !  Function : Contains code for initialising data for two-phase calculation
        !
        !**********************************************************************
        implicit none
        
          call DestroyTwoPhaseData()
          
          call InitialiseTwoPhaseArrays()
      
        end subroutine InitialiseTwoPhaseData


        subroutine DestroyTwoPhaseData()
        !**********************************************************************
        !
        !  Function : Deallocates the arrays used in this module
        !
        !**********************************************************************       
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(LumpedMassWater)) then
            deallocate(LumpedMassWater, stat = IError)
          end if

          if (allocated(LumpedMassNWater)) then
            deallocate(LumpedMassNWater, stat = IError)
          end if

          if (allocated(GravityLoadWater)) then
            deallocate(GravityLoadWater, stat = IError)
          end if
          
          if (allocated(GravityLoadWaterPorosity)) then
            deallocate(GravityLoadWaterPorosity, stat = IError)
          end if

          if (allocated(GravityLoadMixture)) then
            deallocate(GravityLoadMixture, stat = IError)
          end if

          if (allocated(ConductivityMatrix)) then
            deallocate(ConductivityMatrix, stat = IError)
          end if
          
          if (allocated(ConductivityMatrixPorosity)) then
            deallocate(ConductivityMatrixPorosity, stat = IError)
          end if
          
          if (allocated(QVW)) then
            deallocate(QVW, stat = IError)
          end if
          
          if (allocated(QVWPorosity)) then
            deallocate(QVWPorosity, stat = IError)
          end if
          
          if (allocated(QVWPorosityPrevious)) then
            deallocate(QVWPorosityPrevious, stat = IError)
          end if
          
          if (allocated(IntLoadWater)) then
            deallocate(IntLoadWater, stat = IError)
          end if
          
          if (allocated(IntLoadWaterPrevious)) then
            deallocate(IntLoadWaterPrevious, stat = IError)
          end if
          
         if (allocated(DummyIntLoadWaterPrevious)) then
            deallocate(DummyIntLoadWaterPrevious, stat = IError)
          end if
          
          if (allocated(IntLoadWaterPorosity)) then
            deallocate(IntLoadWaterPorosity, stat = IError)
          end if

          if (allocated(IntLoadWaterPorosityPrevious)) then
            deallocate(IntLoadWaterPorosityPrevious, stat = IError)
          end if

          if (allocated(AccelerationWater)) then
            deallocate(AccelerationWater, stat = IError)
          end if
          
          if (allocated(ExtLoadWater)) then
            deallocate(ExtLoadWater, stat = IError)
          end if
          
          if (allocated(ExtLoadWaterPorosity)) then
            deallocate(ExtLoadWaterPorosity, stat = IError)
          end if

          if (allocated(ExtLoadWaterTotal)) then
            deallocate(ExtLoadWaterTotal, stat = IError)
          end if
          
           if (allocated(HydraulicHeadLoadTotal)) then
            deallocate(HydraulicHeadLoadTotal, stat = IError)
          end if
                    
          if (allocated(SpaceTimeExTLoadWater)) then
            deallocate(SpaceTimeExTLoadWater, stat = IError)
          end if        
          
          if (allocated(IsInfiltrationNode)) then
            deallocate(IsInfiltrationNode, stat = IError)
          end if
          
          if (allocated(LumpedNodalPorosityDegSat)) then
              deallocate(LumpedNodalPorosityDegSat, stat = IError)
          end if
          
          if (allocated(IsSeepageFaceNode)) then
            deallocate(IsSeepageFaceNode, stat = IError)
          end if
          
             if (allocated(IsHydraulicHeadNode)) then
            deallocate(IsHydraulicHeadNode, stat = IError)
          end if
          

        end subroutine DestroyTwoPhaseData


        subroutine InitialiseTwoPhaseArrays()
        !**********************************************************************
        !
        !  Function : Initialise the arrays related to two-phase calculation
        !
        !**********************************************************************

        implicit none

          ! local variables
          integer(INTEGER_TYPE) :: IError
          
          if (NFORMULATION==1) then !mixture
              if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
                allocate(LumpedMassWater(Counters%N, Counters%NEntity), stat = IError)
                allocate(LumpedMassNWater(Counters%N, Counters%NEntity), stat = IError)
                allocate(GravityLoadWater(Counters%N, Counters%NEntity), stat = IError)
                allocate(GravityLoadWaterPorosity(Counters%N, Counters%NEntity), stat = IError)
                allocate(GravityLoadMixture(Counters%N, Counters%NEntity), stat = IError)
                allocate(ConductivityMatrix(Counters%N, Counters%NEntity), stat = IError)
                allocate(ConductivityMatrixPorosity(Counters%N, Counters%NEntity), stat = IError)
                allocate(QVW(Counters%N, Counters%NEntity), stat = IError)
                allocate(QVWPorosity(Counters%N, Counters%NEntity), stat = IError)
                allocate(QVWPorosityPrevious(Counters%N, Counters%NEntity), stat = IError)
                allocate(IntLoadWater(Counters%N, Counters%NEntity), stat = IError)
                allocate(IntLoadWaterPrevious(Counters%N, Counters%NEntity), stat = IError)
                allocate(DummyIntLoadWaterPrevious(Counters%N, Counters%NEntity), stat = IError)
                allocate(IntLoadWaterPorosity(Counters%N, Counters%NEntity), stat = IError)
                allocate(IntLoadWaterPorosityPrevious(Counters%N, Counters%NEntity), stat = IError)
                allocate(AccelerationWater(Counters%N, Counters%NEntity), stat = IError)
                allocate(ExtLoadWater(Counters%N, Counters%NEntity), stat = IError)
                allocate(ExtLoadWaterPorosity(Counters%N, Counters%NEntity), stat = IError)
                allocate(ExtLoadWaterTotal(Counters%N, Counters%NEntity,Counters%NWaterLoadSystems), stat = IError)
                allocate(HydraulicHeadLoadTotal(Counters%N, Counters%NEntity), stat = IError)
 				allocate(SpaceTimeExTLoadWater(Counters%N, Counters%NEntity), stat = IError)
                allocate(IsInfiltrationNode(Counters%NodTot), stat = IError)
                allocate(LumpedNodalPorosityDegSat(Counters%N, Counters%NEntity), stat = IError)
                allocate(IsSeepageFaceNode(Counters%NodTot), stat = IError)
                allocate(IsHydraulicHeadNode(Counters%NodTot), stat = IError)
              else
                allocate(LumpedMassWater(1, 1), stat = IError)
                allocate(LumpedMassNWater(1, 1), stat = IError)
                allocate(GravityLoadWater(1, 1), stat = IError)
                allocate(GravityLoadWaterPorosity(1, 1), stat = IError)
                allocate(GravityLoadMixture(1, 1), stat = IError)
                allocate(ConductivityMatrix(1, 1), stat = IError)
                allocate(ConductivityMatrixPorosity(1, 1), stat = IError)
                allocate(QVW(1, 1), stat = IError)
                allocate(QVWPorosity(1, 1), stat = IError)
                allocate(QVWPorosityPrevious(1, 1), stat = IError)
                allocate(IntLoadWater(1, 1), stat = IError)
                allocate(IntLoadWaterPrevious(1, 1), stat = IError)
                allocate(DummyIntLoadWaterPrevious(1, 1), stat = IError)
                allocate(IntLoadWaterPorosity(1, 1), stat = IError)
                allocate(IntLoadWaterPorosityPrevious(1, 1), stat = IError)
                allocate(AccelerationWater(1, 1), stat = IError)
                allocate(ExtLoadWater(1, 1), stat = IError)
                allocate(ExtLoadWaterPorosity(1, 1), stat = IError)
                allocate(ExtLoadWaterTotal(1, 1, 1), stat = IError)
                allocate(HydraulicHeadLoadTotal(1, 1), stat = IError)
				allocate(SpaceTimeExTLoadWater(1, 1), stat = IError)                
                allocate(IsInfiltrationNode(1), stat = IError)
                allocate(LumpedNodalPorosityDegSat(1, 1), stat = IError)
                allocate(IsSeepageFaceNode(1), stat = IError)
                allocate(IsHydraulicHeadNode(1), stat = IError)
              end if
          else !constituents
              allocate(LumpedMassWater(Counters%N, Counters%NEntity), stat = IError)
              allocate(LumpedMassNWater(Counters%N, Counters%NEntity), stat = IError)
              allocate(GravityLoadWater(Counters%N, Counters%NEntity), stat = IError)
              allocate(GravityLoadWaterPorosity(1, 1), stat = IError)
              allocate(GravityLoadMixture(Counters%N, Counters%NEntity), stat = IError)
              allocate(ConductivityMatrix(Counters%N, Counters%NEntity), stat = IError)
              allocate(ConductivityMatrixPorosity(1, 1), stat = IError)
              allocate(QVW(Counters%N, Counters%NEntity), stat = IError)
              allocate(QVWPorosity(1, 1), stat = IError)
              allocate(QVWPorosityPrevious(1, 1), stat = IError)
              allocate(IntLoadWater(Counters%N, Counters%NEntity), stat = IError)
              allocate(IntLoadWaterPrevious(Counters%N, Counters%NEntity), stat = IError)
              allocate(IntLoadWaterPorosity(1, 1), stat = IError)
              allocate(IntLoadWaterPorosityPrevious(1, 1), stat = IError)
              allocate(AccelerationWater(Counters%N, Counters%NEntity), stat = IError)
              allocate(ExtLoadWater(Counters%N, Counters%NEntity), stat = IError)
              allocate(ExtLoadWaterPorosity(1, 1), stat = IError)
              allocate(ExtLoadWaterTotal(Counters%N, Counters%NEntity, Counters%NWaterLoadSystems), stat = IError) 
              allocate(HydraulicHeadLoadTotal(Counters%N, Counters%NEntity), stat = IError)
			  allocate(SpaceTimeExTLoadWater(Counters%N, Counters%NEntity), stat = IError)
			  allocate(IsInfiltrationNode(Counters%NodTot), stat = IError)
              allocate(LumpedNodalPorosityDegSat(Counters%N, Counters%NEntity), stat = IError)
          end if
              
          LumpedMassWater = 0.0
          LumpedMassNWater = 0.0
          GravityLoadWater = 0.0
          GravityLoadWaterPorosity = 0.0
          GravityLoadMixture = 0.0
          ConductivityMatrix = 0.0
          ConductivityMatrixPorosity = 0.0
          QVW = 0.0
          QVWPorosity = 0.0
          QVWPorosityPrevious = 0.0
          IntLoadWater = 0.0
          IntLoadWaterPrevious = 0.0
          DummyIntLoadWaterPrevious = 0.0
          IntLoadWaterPorosity = 0.0
          IntLoadWaterPorosityPrevious = 0.0
          AccelerationWater = 0.0
          ExtLoadWater = 0.0
          ExtLoadWaterPorosity = 0.0
          ExtLoadWaterTotal = 0.0
          SpaceTimeExTLoadWater = 0.0
          HydraulicHeadLoadTotal = 0.0
          IsInfiltrationNode = .false.
          IsSeepageFaceNode = .false.
          IsHydraulicHeadNode = .false.
          LumpedNodalPorosityDegSat = 0.0

        end subroutine InitialiseTwoPhaseArrays


        subroutine FormConsolidationMatrices(TotVelSoil,TotVelWater) 
        !**********************************************************************
        !
        !  Function :  To extrapolate loads and masses from particles to nodes using 
        !              the shape function values evaluated at the particles local position
        !
        !**********************************************************************
      
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity),  &
                 intent(in):: TotVelSoil, TotVelWater

          integer(INTEGER_TYPE) :: IEntity, ILoadSystem

          call FormMatrices(LumpedMassWater, &
                            LumpedMassNWater, &
                            TotVelSoil, &
                            TotVelWater, &
                            ConductivityMatrix, &
                            ConductivityMatrixPorosity)

          if(CalParams%TimeStep==1) then ! first time step of load step
            call ConsolidationForces(ExtLoadWater, &
                                     GravityLoadWater, &
                                     GravityLoadMixture, &
                                     IntLoadWater)
          else
            call ConsolidationExtForces(ExtLoadWater, &
                                        GravityLoadWater, &
                                        GravityLoadMixture)
          end if

       
          !!!!!!    if pressures on boundary nodes depend on hydraulic head from file  !!!!!!  
          if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
              call UpdateMultipliersWaterForSpaceDependency(SpaceTimeExTLoadWater, ReducedDof) 
          
          end if
          do ILoadSystem =1, Counters%NWaterLoadSystems                      
            ExtLoadWater = ExtLoadWater + ExtLoadWaterTotal(:,:,ILoadSystem) * CalParams%Multipliers%WaterACurrent(ILoadSystem)
          end do
          ExtLoadWater = ExtLoadWater + SpaceTimeExTLoadWater

          ! Rotate vectors from global to local coordinate system
          if (IS3DCYLINDRIC) then
            do IEntity = 1, Counters%nEntity 
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof,  &
                           ExtLoadWater(:, IEntity), ExtLoadWater(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof,  &
                           GravityLoadWater(:, IEntity), GravityLoadWater(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof,  &
                           IntLoadWater(:, IEntity), IntLoadWater(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof,  &
                           GravityLoadMixture(:, IEntity), GravityLoadMixture(:, IEntity))
            end do
          end if ! rotation

        end subroutine FormConsolidationMatrices

        subroutine MapWaterMomentumFromMaterialPointsToNodes(Momentum)
        !**********************************************************************
        !
        !  Function :  To map water momentum from material points to the grid points (nodes)
        !
        !  I/O  Momentum :  Nodal water momentum vector, the output of this subroutine
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum

          ! Local variables
          integer(INTEGER_TYPE) :: I, IAEl, IEl, IPart, INode, Nix, ParticleIndex, NodeID, iEntity
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE) :: ParticleMass
          
          Momentum = 0.0

          do IAEl = 1, Counters%NAEl ! loop over all elements
            IEl = ActiveElement(IAEl)
            do IPart = 1, NPartEle(IEl) ! loop over all particles in element
              ParticleIndex = GetParticleIndex(IPart,IEl) ! get the particle ID
              if((NFORMULATION==1).or. &
                        (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid)) then !Material Point can be MIXTURE or LIQUID
              
                if(NFORMULATION==1) then ! 1 Constituent
                  ParticleVelocity = VelocityWaterArray(ParticleIndex,:) ! get particle water velocity vector
                  ParticleMass     = MassWaterArray(ParticleIndex)
                else ! 2 Constituents
                  ParticleVelocity = VelocityArray(ParticleIndex,:) ! get particle velocity vector Particle%Velocity 
                  ParticleMass     = MassArray(ParticleIndex)
                end if
                        
                if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex) ! entity to which particle belongs
                else
                  iEntity = 1
                end if    
                  
                do INode = 1, ELEMENTNODES ! loop over nodes
                
                  NodeID = ElementConnectivities(INode, IEl) ! global node ID
                  do I = 1, NVECTOR
                    ! global storage coordinate
                    Nix = ReducedDof(NodeID) + I
                    ! update nodal momentum
                    Momentum(Nix,iEntity) = Momentum(Nix,iEntity) + ParticleMass * ShapeValuesArray(ParticleIndex,INode) * ParticleVelocity(I)
                  end do
                  
                end do ! loop over nodes
                
               end if ! NumbOfLayers=1 or LIQUID MatPoint
            
            end do ! loop over particles
          end do ! loop over elements
         
        end subroutine MapWaterMomentumFromMaterialPointsToNodes
        
        subroutine GetNodalWaterVelocityFromNodalWaterMomentum(Momentum)
        !**********************************************************************
        !
        !    Function:  To calculate the nodal water velocities from nodal water mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        !**********************************************************************

        implicit none
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity),  &
                                            intent(in) :: Momentum
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          
               do IDOF = 1, Counters%N                                
            do J = 1, Counters%nEntity !loop through all entities    
                if(LumpedMassWater(IDOF,J)/=0) then
                    TotalVelocityWater(IDOF,J) = ( Momentum(IDOF,J) /  &
                                        LumpedMassWater(IDOF,J) )*  &
                                        PboundaryWater(IDOF)
                else
                    TotalVelocityWater(IDOF,J) = 0.0
                end if
            end do
          end do
          
        end subroutine GetNodalWaterVelocityFromNodalWaterMomentum
        
        
        subroutine GetQVWArray(QVW)
        !**********************************************************************
        !
        !    Function:  To calculate the QVW array
        !
        !    O  QVW : The output of this subroutine
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity),      &
                                intent(inout) :: QVW    
              
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          real(REAL_TYPE) :: VminusW
          
          if(NFORMULATION==2) then
          
            QVW = 0.0
          
            do IDOF = 1, Counters%N       !loop over all degrees of freedom
              do J = 1, Counters%nEntity  !loop over all entities
                VminusW = TotalVelocitySoil(IDOF,J) - &
                          TotalVelocityWater(IDOF,J)
                QVW(IDOF,J) =  ConductivityMatrix(IDOF,J) * VminusW
              end do
            end do
          else
            QVW = 0.0
            QVWPorosity = 0.0
          
            do IDOF = 1, Counters%N       !loop over all degrees of freedom
              do J = 1, Counters%nEntity  !loop over all entities
                VminusW = TotalVelocitySoil(IDOF,J) - &
                          TotalVelocityWater(IDOF,J)
                QVW(IDOF,J) =  ConductivityMatrix(IDOF,J) * VminusW
                QVWPorosity(IDOF,J) =  ConductivityMatrixPorosity(IDOF,J) * VminusW
              end do
            end do
          end if
            
        end subroutine GetQVWArray       
        
        subroutine GetQVWArrayPorosity(QVWPorosity)    
        !**********************************************************************
        !
        !    Function:  To calculate the QVW array
        !
        !    O  QVW : The output of this subroutine
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity),      &
                                intent(inout) :: QVWPorosity    
              
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          real(REAL_TYPE) :: VminusW
          
          QVWPorosity = 0.0
   
          do IDOF = 1, Counters%N       !Loop over all degrees of freedom
            do J = 1, Counters%nEntity  !loop over all entities
              VminusW = TotalVelocitySoil(IDOF,J) - &
                        TotalVelocityWater(IDOF,J)
              QVWPorosity(IDOF,J) =  ConductivityMatrixPorosity(IDOF,J) * VminusW
            end do
          end do
            
        end subroutine GetQVWArrayPorosity
               
        subroutine GetWaterInertiaArray(WaterInertia)
        
        !**********************************************************************
        !
        !    Function:  To calculate the WaterInertia array
        !
        !    O  WaterInertia : The output of this subroutine
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity),  &
                                intent(inout) :: WaterInertia    
                   
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          
          WaterInertia = 0.0

          do IDOF = 1, Counters%N       !Loop over all degrees of freedom
            do J = 1, Counters%nEntity  !loop over all entities
              WaterInertia(IDOF,J) = LumpedMassNWater(IDOF,J) * &
                                      AccelerationWater(IDOF,J)
            end do
          end do
            
        end subroutine GetWaterInertiaArray
        
        
        subroutine CalculateWaterIncrementalNodalAcceleration(RateofMomentum)
                                          
        !**********************************************************************
        !
        !  Function :  To calculate the incremental nodal accelerations (liquid phase)
        !
        !  I  RateofMomentum :  The array which stores the rate of momentum (liquid) 
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: RateofMomentum 

          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          
           do IDOF = 1, Counters%N ! loop over all degrees of freedom
            do J = 1, Counters%nEntity  ! loop over all entities
              if( LumpedMassWater(IDOF,J)/=0 ) then
                AccelerationWater(IDOF,J) = ( RateofMomentum(IDOF,J) / LumpedMassWater(IDOF,J) ) * PboundaryWater(IDOF)
              else
                AccelerationWater(IDOF,J) = 0.0
              end if
            end do
          end do
          
        end subroutine CalculateWaterIncrementalNodalAcceleration
        
        
        subroutine UpdateParticleWaterVelocityAndMapMomentumW(Momentum)                        
        !**********************************************************************
        !
        !    Function:  To update particles total water velocities and accelerations.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum  
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
          integer(INTEGER_TYPE), dimension(ELEMENTNODES, NVECTOR) :: IDof
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleIncrementalVelocity
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity, NVECTOR) :: NodAcc
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) :: Time, ParticleMass
                              
          Momentum = 0.0
           
          Time = CalParams%TimeIncrement

          NoEn = Counters%nEntity

             do IAEl = 1, Counters%NAEl                                      ! Loop over all elements
            IEl = ActiveElement(IAEl)
            
            do I = 1, NVECTOR
              ! Loop over dimensions
              IDof(1:ELEMENTNODES, I) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
              NodAcc(1:ELEMENTNODES, 1:NoEn, I) = AccelerationWater(IDof(1:ELEMENTNODES, I), 1:NoEn)
            end do
            
              do IPart = 1, NPartEle(IEl)                                 ! Loop over all particles in element
                ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                
                if((NFORMULATION==1).or. &
                        (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid)) then 
              
                  if(NFORMULATION==1) then ! 1 Constituent
                    ParticleVelocity = VelocityWaterArray(ParticleIndex,:)  ! get particle water velocity vector
                    ParticleMass     = MassWaterArray(ParticleIndex)
                  else ! 2 Constituents
                    ParticleVelocity = VelocityArray(ParticleIndex,:) ! get particle velocity vector 
                    ParticleMass     = MassArray(ParticleIndex)
                  end if
                  ParticleShape = ShapeValuesArray(ParticleIndex,:)
                  ParticleIncrementalVelocity = 0.0
                  ParticleAcceleration = 0.0
                 if (CalParams%ApplyContactAlgorithm) then
                   iEntity = EntityIDArray(ParticleIndex)
                  else
                   iEntity = 1
                 end if 
                  
                 do I = 1, NVECTOR ! Loop over dimensions
                   do INode = 1, ELEMENTNODES ! loop over element nodes
                     ! Particle acceleration
                     ParticleAcceleration(I) = ParticleAcceleration(I) + ParticleShape(INode) * NodAcc(INode, iEntity, I)
                     ! Particle velocity
                     ParticleIncrementalVelocity(I) = ParticleIncrementalVelocity(I) + Time * ParticleShape(INode) * NodAcc(INode, iEntity, I)
                   end do  !Loop over dimensions
                 end do !Loop over nodes
               
                 ParticleVelocity = ParticleVelocity + ParticleIncrementalVelocity

                 ! nodal momentum
                 do I = 1, NVECTOR ! Loop over dimensions
                   Momentum(IDof(1:ELEMENTNODES, I), iEntity) = Momentum(IDof(1:ELEMENTNODES, I), iEntity) + ParticleMass * ParticleShape * ParticleVelocity(I)
                 end do !Loop over dimensions

     
                 if(NFORMULATION==1) then ! 1 constituent
                   VelocityWaterArray(ParticleIndex,:) = ParticleVelocity
                 else  ! 2 constituents
                   VelocityArray(ParticleIndex,:) = ParticleVelocity
                   AccelerationArray(ParticleIndex,:) =  ParticleAcceleration
                 end if ! end constituent
     
              end if !NumbOfLayers=1 or LIQUID MatPoint
              
            end do !Loop over particles
          end do !elements
          
        end subroutine UpdateParticleWaterVelocityAndMapMomentumW  
        
        
        subroutine UpdateParticleWaterAcceleration()
        !**********************************************************************
        !
        !    Function:  update acceleration of the water in the material points
        !
        !*********************************************************************  
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
          integer(INTEGER_TYPE), dimension(ELEMENTNODES, NVECTOR) :: Nix
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration
          real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity, NVECTOR) :: NodAcc
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          
          if (CalParams%NumberOfPhases/=2) RETURN
          if (.not.CalParams%ApplyAbsorbingBoundary) RETURN 

          NoEn= Counters%NEntity
          
          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            do I = 1, NVECTOR
              ! Loop over dimensions
              Nix(1:ELEMENTNODES, I) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
              NodAcc(1:ELEMENTNODES, 1:NoEn, I) = AccelerationWater(Nix(1:ELEMENTNODES, I), 1:NoEn)
            end do
            
              do IPart = 1, NPartEle(IEl)                
                ParticleIndex = GetParticleIndex(IPart, IEl)
                ParticleShape = ShapeValuesArray(ParticleIndex,:)
                ParticleAcceleration = 0.0
                if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex) 
                else
                  iEntity = 1
                end if 
                  
                do INode = 1, ELEMENTNODES
                  do I = 1, NVECTOR ! Loop over dimensions
                    ! Particle acceleration
                    ParticleAcceleration(I) = ParticleAcceleration(I) + ParticleShape(INode) * NodAcc(INode, iEntity, I)
                  end do
                end do
               
               Particles(ParticleIndex)%AccelerationWater = ParticleAcceleration
            end do
          end do
          
        end subroutine UpdateParticleWaterAcceleration  

        subroutine CalculateWaterVolumetricStrain(DUTotAll)
        !**********************************************************************
        !
        !    Function:  Calculate the volumetric strain (water).
        !
        !     DUTot : Incremental nodal displacements (water0
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
          implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTotAll
          !Local variables
          real(REAL_TYPE), dimension(Counters%N) :: DUTotEnt
          real(REAL_TYPE), dimension(NTENSOR) :: Eps
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac
          integer(INTEGER_TYPE) :: IAElement, IElement, NElemPart, IParticle, ParticleIndex, iEntity, iEntityDefault
     
          !========================================================
          !set the nodal displacement increments to the first entity
          !only need to change if the contact model is used and 
          !the next particle belongs to a different entity
          iEntityDefault = 1    !set default to first entity
          DUTotEnt(1:Counters%N) = DUTotAll(1:Counters%N, iEntityDefault)
          !=======================================================

          do IAElement = 1, Counters%NAEl     !loop over all elements
            IElement = ActiveElement(IAElement)
     
              NElemPart = NPartEle(IElement)    !number of particles in element
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if
              
              do IParticle = 1, NElemPart       !loop over all particles of the element
              
                ParticleIndex = GetParticleIndex(IParticle, IElement)   !particle global ID
                
                if((NFORMULATION==1).or.(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid)) then
                  !NumberOfLayers = 1 or LIQUID MatPoint

                  ! Determine B matrix of deformed element at the particle local position
                  call BMatrix(Particles(ParticleIndex)%LocPos, &
                             ELEMENTNODES, Counters%NEl,  &
                             Counters%NodTot, NVECTOR, &
                             IElement, ElementConnectivities, &
                             NodalCoordinatesUpd,  &
                             BMatrixDeformed, DetJac)
     
                  !get the nodal displacement increments for the specific entity
                  !=================================
                  if (CalParams%ApplyContactAlgorithm) then
                    iEntity = EntityIDArray(ParticleIndex) 
                    if (iEntity /= iEntityDefault) then
                      DUTotEnt(1:Counters%N) = DUTotAll(1:Counters%N,iEntity)
                      ! set the entity ID for the next particle loop
                      iEntityDefault = iEntity
                    end if
                  else
                    DUTotEnt(1:Counters%N) = DUTotAll(1:Counters%N, 1)
                  end if
                  !====================================

                  ! Eps = (Exx, Eyy, Ezz, Gxy, Gyz, Gzx)
                  call Get_Strain(IElement, IParticle, ElementConnectivities, BMatrixDeformed, DUTotEnt, ReducedDof, Eps)

                  Particles(ParticleIndex)%WaterVolumetricStrain = Eps(1) + Eps(2) + Eps(3) ! valid for 2D and 3D
                    
                end if !NumberOfLayers = 1 or LIQUID MatPoint
              
              end do
          end  do

         end subroutine CalculateWaterVolumetricStrain
        
         subroutine GetGravityLoadWaterPorosity(GravityLoadWaterLoc)
        !**********************************************************************
        !
        !    Function:  To extrapolate real water load from particles to nodes using 
        !               the shape function values evaluated at the particles local 
        !               position.
        !
        !**********************************************************************

        implicit none
       
         real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: GravityLoadWaterLoc
         ! local variables
          integer(INTEGER_TYPE) :: I, IAElement, IPart, INode, IElement, ParticleIndex, NodeID, GlobDof, iEntity
          real(REAL_TYPE), dimension(NVECTOR) :: PartGravity

          !must set to zero here, because used in summation
          GravityLoadWaterLoc = 0.0

    
          do IAElement = 1, Counters%NAEl ! Loop over all elements
           IElement = ActiveElement(IAElement)
          
          !reset to zero
          PartGravity = 0.0

          do IPart = 1, NPartEle(IElement) ! Loop over all particles in element
            ParticleIndex = GetParticleIndex(IPart, IElement) ! Get the particle ID
            if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex) 
            else
                iEntity = 1
            endif  
            
            do INode = 1, ELEMENTNODES ! Loop over all nodes in element IEl
              NodeID = iabs(ElementConnectivities(INode, IElement) ) ! Nodal global ID
              GlobDof = ReducedDof(NodeID) 

              if (CalParams%NumberOfPhases==3) then
                PartGravity =                                      & ! real water gravity
                  ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyWater * Particles(ParticleIndex)%Porosity * Particles(ParticleIndex)%DegreeSaturation  
              else 
                PartGravity =                                      & ! real water gravity
                  ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyWater * Particles(ParticleIndex)%Porosity    ! real water gravity
              end if

              do I = 1, NVECTOR
                GravityLoadWaterLoc(GlobDof + I,iEntity) = GravityLoadWaterLoc(GlobDof + I,iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
              end do
                  
            end do ! nodes

          end do ! particles
         
          end do ! elements
          
         if (IS3DCYLINDRIC) then ! rotation is needed
           do IEntity = 1, Counters%nEntity 
             call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, GravityLoadWater(:, IEntity), GravityLoadWater(:, IEntity))
           end do
         end if ! rotation

        end subroutine GetGravityLoadWaterPorosity
         

        subroutine FormMatrices(LumpedMassWaterLoc, LumpedMassNWaterLoc, TotVelSoilLoc, &
                                  TotVelWaterLoc, ConductivityMatrixLoc, ConductivityMatrixPorosityLoc)
        !**********************************************************************
        !
        !    Function:  To Form the Lumped mass matrix (water)
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: &
            LumpedMassWaterLoc, LumpedMassNWaterLoc,ConductivityMatrixLoc, ConductivityMatrixPorosityLoc, TotVelSoilLoc, TotVelWaterLoc

          ! Local variables
          integer(INTEGER_TYPE) :: IAEl, IEl, IPart, INode, IDof, ParticleIndex, NodeID, iEntity, INodeNew, IDofNew, NodeIDNew
          integer(INTEGER_TYPE) :: I, J
          real(REAL_TYPE) :: MassEntry, MassNEntry, ConductivityEntry, ConductivityNEntry
          real(REAL_TYPE) :: LiquidViscosity, IntrinsicPermeability
          real(REAL_TYPE) :: ConductivityEntry_1, ConductivityEntry_2, ErgunConstantF
          real(REAL_TYPE) :: AbsRelativeVelocityNodal, AbsRelativeVelocityMatPoint
          real(REAL_TYPE) :: ConstDensityLiquid, PorosityMP, UpdIntrinsicPermeability, IntrinsicPermeabilityLiquid
          real(REAL_TYPE) :: TwoLayerErgunLawDiameterPartic
          real(REAL_TYPE) :: ConcentraRatioLiquidS
                       
          LumpedMassWaterLoc = 0.0
          LumpedMassNWaterLoc = 0.0
          ConductivityMatrixLoc = 0.0
          if(NFORMULATION==1) then !1 Constituent
            ConductivityMatrixPorosityLoc = 0.0
          end if
          
          !--- Calculate  data  ----
           ConstDensityLiquid = 0.0
           LiquidViscosity = 0.0
                
           do J = 1, Counters%NLayers ! loop over all material points in element
             if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
                 ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                 LiquidViscosity = MatParams(J)%ViscosityLiquid
             else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
                 ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                 LiquidViscosity = MatParams(J)%ViscosityLiquid 
             else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-solid'.or.MatParams(j)%MaterialPhases=='1-phase-solid') then
                 IntrinsicPermeabilityLiquid = MatParams(J)%IntrinsicPermeabilityLiquid
             end if
           end do 
           ! ---------------------------------------------------
           
          do IAEl = 1, Counters%NAEl                                    ! Loop over all elements
            IEl = ActiveElement(IAEl)
            do IPart = 1, NPartEle(IEl)                                 ! Loop over all particles in element
              
              ParticleIndex = GetParticleIndex(IPart, IEl)              ! Get the particle ID      
              !------------------------------------------------------------------------ 
              if(NFORMULATION==1) then !1 Constituent
                
                if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex) ! get entity to which particle belongs
                else
                  iEntity = 1
                end if
              
                do INode = 1, ELEMENTNODES                           ! loop over element nodes
                    NodeID = ElementConnectivities(INode, IEl)              ! Global node ID
                    IDof = ReducedDof(NodeID)                              ! Global storage coordinate of x-val
                    
                    MassEntry =  MassWaterArray(ParticleIndex) * ShapeValuesArray(ParticleIndex,INode)           !calculate mass contribution
                    ! nodal mass 
                    do I = 1, NVECTOR
                      !add mass contribution
                      LumpedMassWaterLoc(IDof+I, iEntity) = LumpedMassWaterLoc(IDof+I, iEntity) + MassEntry
                    end do
                    
                    if ((CalParams%NumberOfPhases==3).or. &
                        ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation)))then                 !calculate mass contribution
                      MassNEntry = MassEntry * Particles(ParticleIndex)%Porosity * Particles(ParticleIndex)%DegreeSaturation
                    else
                      MassNEntry = MassEntry * Particles(ParticleIndex)%Porosity
                    end if
     
                    ! nodal mass
                    do I = 1, NVECTOR
                      ! add mass contribution
                      LumpedMassNWaterLoc(IDof+I, iEntity) = LumpedMassNWaterLoc(IDof+I, iEntity) + MassNEntry
                    end do

                    if (Particles(ParticleIndex)%Conductivity/=0.0) then  
                      if ((CalParams%NumberOfPhases==3).or. &
                        ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation)))then                              !calculate conductivity contribution
                        ConductivityEntry = MassWaterArray(ParticleIndex) * &
                                           Particles(ParticleIndex)%Porosity * & ! In unsaturated calculation the degree of saturation has to be taken in account 
                                           Particles(ParticleIndex)%DegreeSaturation * &
                                           CalParams%GravityData%GAccel * &
                                           (1.0 / Particles(ParticleIndex)%Conductivity) * &
                                           ShapeValuesArray(ParticleIndex,INode)
                      else
                        ConductivityEntry = MassWaterArray(ParticleIndex) * &
                                             Particles(ParticleIndex)%Porosity * &
                                             CalParams%GravityData%GAccel * &
                                             (1.0 / Particles(ParticleIndex)%Conductivity) * &
                                             ShapeValuesArray(ParticleIndex,INode)
                      end if
                    
                      ! nodal conductivity
                      do I = 1, NVECTOR
                        ! add conductivity contribution
                        ConductivityMatrixLoc(IDof+I, iEntity) = ConductivityMatrixLoc(IDof+I, iEntity) + ConductivityEntry
                      end do
      
                      if ((CalParams%NumberOfPhases==3).or. &
                        ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation)))then                                !calculate conductivity contribution
                         ConductivityNEntry =  ConductivityEntry *   &
                                               Particles(ParticleIndex)%Porosity * &
                                               Particles(ParticleIndex)%DegreeSaturation
                      else
                         ConductivityNEntry =  ConductivityEntry *   &
                                               Particles(ParticleIndex)%Porosity
                      end if
     
                      ! nodal conductivity
                      do I = 1, NVECTOR
                        ! add conductivity contribution
                        ConductivityMatrixPorosityLoc(IDof+I, iEntity) = ConductivityMatrixPorosityLoc(IDof+I, iEntity) + ConductivityNEntry
                      end do
 
                  end if 
                end do !Loop over nodes
              end if !1 Constituent
              
              !------------------------------------------------------------------------
              if(NFORMULATION==2) then !2 Constituents
                  
                iEntity = 1 
                  
                if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then !only if LIQUID
              
                  do INode = 1, ELEMENTNODES                         ! loop over element nodes
                    NodeID = ElementConnectivities(INode, IEl)          ! Global node ID
                    IDof = ReducedDof(NodeID)                           ! Global storage coordinate of x-val

                    MassEntry = MassArray(ParticleIndex) * ShapeValuesArray(ParticleIndex,INode)       !calculate mass contribution
                    ! nodal mass 
                    do I = 1, NVECTOR
                      !add mass contribution
                      LumpedMassWaterLoc(IDof+I, iEntity) = LumpedMassWaterLoc(IDof+I, iEntity) + MassEntry
                    end do
                
                    MassNEntry = MassEntry                              !calculate mass contribution
                    ! nodal mass
                    do I = 1, NVECTOR
                      ! add mass contribution
                      LumpedMassNWaterLoc(IDof+I, iEntity) = LumpedMassNWaterLoc(IDof+I, iEntity) + MassNEntry
                    end do
     
                  end do ! loop nodes
                end if !only if LIQUID
                
                if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) then !only if SOLID
                   
                  !! ---------- Define conductivity if Conductivity=0 ------------------
                  if(Particles(ParticleIndex)%Conductivity<= 0.0) then
                    Particles(ParticleIndex)%Conductivity = ConstDensityLiquid * CalParams%GravityData%GAccel * IntrinsicPermeabilityLiquid / LiquidViscosity
                  end if

                  !! ---------- Update permeability ------------------
                  if(CalParams%TwoLayerApplyUpdatePermeability) then    
                     PorosityMP = Particles(ParticleIndex)%EffPorosity
                     if(PorosityMP>POROSITYTHRES) PorosityMP = POROSITYTHRES
                     
                     if (MaterialIDArray(ParticleIndex) == CalParams%FirstSolidMaterialIndex) then 
                         TwoLayerErgunLawDiameterPartic = CalParams%TwoLayerErgunLawDiameterPartic
                     else if (MaterialIDArray(ParticleIndex) == CalParams%SecondSolidMaterialIndex) then
                         TwoLayerErgunLawDiameterPartic = CalParams%TwoLayerErgunLawDiameterPartic2
                     end if 

                     UpdIntrinsicPermeability = (TwoLayerErgunLawDiameterPartic * TwoLayerErgunLawDiameterPartic * &
                            PorosityMP * PorosityMP * PorosityMP ) / ( CalParams%ERGUNCONSTANTA * ( 1.0 - PorosityMP ) * ( 1.0 - PorosityMP ))

                      Particles(ParticleIndex)%Conductivity = ConstDensityLiquid * CalParams%GravityData%GAccel * UpdIntrinsicPermeability / LiquidViscosity
                  end if ! Update permeability
                  
                  !! ---------- Calculate conductivity matrix ------------------
                  if (TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLIDLIQUID) then   
                   do INode = 1, ELEMENTNODES                           ! loop over element nodes
                      NodeID = ElementConnectivities(INode, IEl)           ! Global node ID
                      IDof = ReducedDof(NodeID)                            ! Global storage coordinate of x-val
                      
                      IntrinsicPermeability = Particles(ParticleIndex)%Conductivity * LiquidViscosity / (CalParams%GravityData%GAccel * ConstDensityLiquid)
     
                      ConcentraRatioLiquidS = Particles(ParticleIndex)%ConcentrationRatioLiquidS
                      if(ConcentraRatioLiquidS>Particles(ParticleIndex)%EffPorosity) then
                         ConcentraRatioLiquidS = Particles(ParticleIndex)%EffPorosity
                      end if
                      
                      ConductivityEntry_1 = LiquidViscosity * ConcentraRatioLiquidS * &
                            ConcentraRatioLiquidS / IntrinsicPermeability * Particles(ParticleIndex)%IntegrationWeight * &
                            ShapeValuesArray(ParticleIndex,INode)     
     
                      ConductivityEntry_2 = 0.0
                               
                      if(CalParams%TwoLayerApplyErgunLaw) then !Ergun Flow
                          
                        PorosityMP = Particles(ParticleIndex)%EffPorosity
                        if(PorosityMP>POROSITYTHRES) PorosityMP = POROSITYTHRES
                
                          ! Calculated absolute value of relative velocity per material point
                          AbsRelativeVelocityNodal = 0.0
                          AbsRelativeVelocityMatPoint = 0.0
                          do INodeNew = 1, ELEMENTNODES                           ! loop over element nodes
                             NodeIDNew = ElementConnectivities(INodeNew, IEl)        ! Global node ID
                             IDofNew = ReducedDof(NodeIDNew)+1                        ! Global storage coordinate of x-val
                             AbsRelativeVelocityNodal = 0.0
                             do I = 1, NVECTOR
                               AbsRelativeVelocityNodal = AbsRelativeVelocityNodal + &
                                                          (TotVelSoilLoc(IDofNew + I, iEntity) - TotVelWaterLoc(IDofNew + I, iEntity)) * &
                                                          (TotVelSoilLoc(IDofNew + I, iEntity) - TotVelWaterLoc(IDofNew + I, iEntity))
                             end do
                             
                             AbsRelativeVelocityNodal = sqrt(AbsRelativeVelocityNodal)
                             AbsRelativeVelocityMatPoint = AbsRelativeVelocityMatPoint +  AbsRelativeVelocityNodal * ShapeValuesArray(ParticleIndex,INodeNew)
                             
                          end do
                    
                          ErgunConstantF = CalParams%ERGUNCONSTANTB / sqrt(CalParams%ERGUNCONSTANTA * PorosityMP * PorosityMP * PorosityMP )

                          ConductivityEntry_2 = ConstDensityLiquid * &
                              ErgunConstantF * ConcentraRatioLiquidS * ConcentraRatioLiquidS * ConcentraRatioLiquidS / &
                              sqrt(IntrinsicPermeability) * AbsRelativeVelocityMatPoint * ShapeValuesArray(ParticleIndex,INode) * &
                              Particles(ParticleIndex)%IntegrationWeight

                      end if ! Ergun law
                          
                      ConductivityEntry = ConductivityEntry_1 + ConductivityEntry_2
     
                      do I = 1, NVECTOR
                        ! add conductivity contribution
                        ConductivityMatrixLoc(IDof+I, iEntity) = ConductivityMatrixLoc(IDof+I, iEntity) + ConductivityEntry
                      end do
     
                   end do !Loop over nodes
                end if ! ContainedMaterialTypes==SOIL_FLUID
              end if !only if SOLID
              end if !2 Constituent           
            end do !Loop over particles
          end do !elements

        end subroutine FormMatrices

        
        subroutine ConsolidationForces(ExtLoadWaterLoc, &
                                       GravityLoadWaterLoc, &
                                       GravityLoadMixtureLoc, &
                                       IntLoadWaterLoc)
        !**********************************************************************
        !
        !    Function:  Calculation of the equivalent nodal forces due to
        !               a given stress InternalLD = Integral {BT*Sigma}.
        !
        !    O ExternalLD : External load
        !    O InternalLD : Internal load
        !    O GravityLD : Gravity load
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), &
            intent(inout) :: ExtLoadWaterLoc, IntLoadWaterLoc, GravityLoadWaterLoc, GravityLoadMixtureLoc
          
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
          integer(INTEGER_TYPE) :: IntGlo, IEl, Int, I, J, NN, NElemPart, iEntityID, iNode, IDof
          real(REAL_TYPE) :: WtN, Det, S, TotIntWeightLiqMP
          integer(INTEGER_TYPE) :: IPart, IAEl, ParticleIndex, NodeID, ILoadSystem
          real(REAL_TYPE), dimension(NVECTOR) :: PartLoad, PartGravity, PartGravityMix
          real(REAL_TYPE), dimension(NTENSOR) :: SW
          real(REAL_TYPE) :: Position
          real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)
          logical :: IsLoadOnMP

          B = 0.0
          IntGlo = 0
          IEl = 0
          Int = 0
          J = 0
          NN = 0
          NElemPart = 0
          iEntityID = 0
          iNode = 0
          WtN = 0.0
          Det = 0.0
          S = 0.0
          TotIntWeightLiqMP = 0.0
          IPart = 0
          IAEl = 0
          ParticleIndex = 0
          NodeID = 0
          PartLoad = 0.0
          PartGravity = 0.0
          PartGravityMix = 0.0
          
          SW = 0.0
          ExtLoadWaterLoc = 0.0
          IntLoadWaterLoc = 0.0
          GravityLoadWaterLoc = 0.0
          GravityLoadMixtureLoc = 0.0
          IsLoadOnMP =(Counters%NLoadedElementSidesWaterMatPoints+Counters%NLoadedElementSidesWaterMatPointsB)>0

            do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            
            !========================================================================
            !=====================1 Layer============================================
            !========================================================================
            
            if(NFORMULATION==1) then !1 Constituents

!//////////////////// External force calculation //////////////////////// 
              do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
                ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                PartLoad = 0.0
                PartGravity = 0.0
                PartGravityMix = 0.0
              
                if (CalParams%ApplyContactAlgorithm) then
                  iEntityID = EntityIDArray(ParticleIndex)  !CC entity to which particle belongs
                else
                  iEntityID = 1
                endif  

                do INode = 1, ELEMENTNODES 
                  NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                  IDof = ReducedDof(NodeID)
 
                  if(IsLoadOnMP) then
                    do ILoadSystem = 1, Counters%NWaterLoadSystems
                      PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExtWater(:,ILoadSystem) * CalParams%Multipliers%WaterACurrent(ILoadSystem)
                      do I = 1, NVECTOR
                        ExtLoadWaterLoc(IDof+I, iEntityID) =  ExtLoadWaterLoc(IDof+I, iEntityID) + PartLoad(I)
                      end do
                    end do
                  end if

                  PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyWater

                  do I = 1, NVECTOR
                    GravityLoadWaterLoc(IDof+I, iEntityID) =  GravityLoadWaterLoc(IDof+I, iEntityID) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                  end do

                  PartGravityMix = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyMixed

                  do I = 1, NVECTOR
                    GravityLoadMixtureLoc(IDof+I, iEntityID) =  GravityLoadMixtureLoc(IDof+I, iEntityID) + PartGravityMix(I) * CalParams%Multipliers%GravityCurrent
                  end do

                end do
              end do
            
!//////////////////// Internal force calculation ////////////////////////
              ! Determine number of integration points inside element
              if (IsParticleIntegration(IEl) ) then ! True - particle based integration, false - Gauss point based integration
                NElemPart = NPartEle(IEl)  ! Number of particles in element
              else
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if
              
              call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element
        
              !------------------------------------ INTEGRATION POINT LOOP --------------------------------
              do Int = 1, NElemPart ! Loop over number of integration points per element IEl
                
                ! Determine global ID of integration point 
                IntGlo = GetParticleIndex(Int, IEl)
              
                ! Set the integration weight
                if (IsParticleIntegration(IEl) ) then
                  ! use the integration weight of particle if it is partially filled, otherwise the weight of gauss point is used
                  WTN = Particles(IntGlo)%IntegrationWeight
                  if ( ISAXISYMMETRIC ) then ! the integration weight of the MP is not corrected due to axisymmetry
                    Position = GlobPosArray(IntGlo, 1) ! index 1 is readial direction
                    ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                  end if
                else
                  if ( ISAXISYMMETRIC ) then
                    Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                    ShapeValues(:) = GPShapeFunction(Int, :)
                    WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                  end if  
                end if

                S = Particles(IntGlo)%WaterPressure * WTN

                ! get particle entity
                if (.not.CalParams%ApplyContactAlgorithm) then
                  iEntityID = 1
                else
                  iEntityID = EntityIDArray(IntGlo) 
                end if

                do INode = 1,ELEMENTNODES  ! loop through element nodes
                  
                  nn = ElementConnectivities(iNode,iel) ! get global node number
                  IDof = ReducedDof(nn)
                  ! nodal load
                  
                  do I = 1, NVECTOR
                    IntLoadWaterLoc(IDof+I, iEntityID) = IntLoadWaterLoc(IDof+I, iEntityID) + B(I, INode) * S
                  end do
                  
                  if ( ISAXISYMMETRIC ) then
                    ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                    IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1, iEntityID) + S * ShapeValues(INode) / Position
                  end if
                  
                end do ! loop element nodes
                
              end do ! loop element points
              !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
            end if ! 1 Constituent
            
            !========================================================================
            !====================2 Layers============================================
            !========================================================================
            
            if(NFORMULATION==2) then !2 Constituents
              if(CalParams%NumberOfPhases>1) then ! 2 Phases

!//////////////////// External force calculation //////////////////////// 
                do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
                  ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              
                  if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then ! Only LIQUID material Point
  
                    PartLoad = 0.0
                    PartGravity = 0.0
                    PartGravityMix = 0.0  ! Not needed
              
                    iEntityID = 1
               
                    do INode = 1, ELEMENTNODES ! Element Nodes
                      NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                      IDof = ReducedDof(NodeID)

                      ! Assemble External Load on Liquid Material Point
                      if(IsLoadOnMP) then 
                        do ILoadSystem = 1, Counters%NWaterLoadSystems  
                        !  ILoadSystem = 1
                          PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExtWater(:,ILoadSystem)* CalParams%Multipliers%WaterACurrent(ILoadSystem)
                          do I = 1, NVECTOR
                            ExtLoadWaterLoc(IDof+I, iEntityID) =  ExtLoadWaterLoc(IDof+I, iEntityID) + PartLoad(I) 
                          end do
                        end do
                      end if
 
                      PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBody 

                      do I = 1, NVECTOR
                        GravityLoadWaterLoc(IDof+I, iEntityID) =  GravityLoadWaterLoc(IDof+I, iEntityID) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                      end do

                    end do ! Element Nodes
                  end if ! Only LIQUID material Point
                end do ! Loop MatPoints
                
!//////////////////// Internal force calculation ////////////////////////
                ! Determine number of integration points inside element
                if (IsParticleIntegration(IEl) ) then ! True - particle based integration, false - Gauss point based integration
                  NElemPart = NPartEle(IEl)  ! Number of particles in element
                else
                  NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
                end if
              
                TotIntWeightLiqMP = 0.0
                do J = 1, NPartEle(IEl)
                  ParticleIndex = GetParticleIndex(J, IEl)
                  if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                    TotIntWeightLiqMP = TotIntWeightLiqMP + Particles(ParticleIndex)%IntegrationWeight
                  end if
                end do
                
                call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element
        
                !------------------------------------ INTEGRATION POINT LOOP --------------------------------
                do Int = 1, NElemPart ! Loop over number of integration points per element IEl
              
                  ! Determine global ID of integration point 
                  IntGlo = GetParticleIndex(Int, IEl)
                  
                  if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then ! Only LIQUID material Point
              
                    ! Set the integration weight
                    WTN = Particles(IntGlo)%IntegrationWeight
                    if(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0) then
                      WTN = ElementSpace(IEl) / TotIntWeightLiqMP * Particles(IntGlo)%IntegrationWeight
                      if(TwoLayerData%Elements(IEl)%IsBoundaryOfLiquidElements) then ! Boundary of Liquid domain
                        if((TotIntWeightLiqMP/ElementSpace(IEl))<CalParams%RequiredDegreeOfFilling) then
                          WTN = Particles(IntGlo)%IntegrationWeight 
                        end if
                      end if
                    end if
                    
                    if ( ISAXISYMMETRIC ) then
                      if (IsParticleIntegration(IEl) ) then
                        ! the integration weight of the MP is not corrected due to axisymmetry
                        Position = GlobPosArray(IntGlo, 1) ! index 1 is readial direction
                        ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                      else
                        Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                        ShapeValues(:) = GPShapeFunction(Int, :)
                        WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                      end if
                    end if  
                    
                    do J = 1, NTENSOR
                      SW(J) = SigmaEffArray(IntGlo,J)  * WTN * (1.0 - TwoLayerData%ELEMENTS(IEl)%ConcentrationRatioSolidL)
                    end do
                    
                    !get material point entity
                    iEntityID = 1
                    
                    do iNode = 1,ELEMENTNODES ! loop over element nodes
                      nn = ElementConnectivities(iNode,iel) ! get global node number
                      IDof = ReducedDof(nn) ! global storage coordinate of x-val
               
                      if (NVECTOR == 3) then ! 3D case
                        IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1,iEntityID) + B(1,iNode)*SW(1) + B(2,iNode)*SW(4) + B(3,iNode)*SW(6)
                        IntLoadWaterLoc(IDof+2, iEntityID) = IntLoadWaterLoc(IDof+2,iEntityID) + B(1,iNode)*SW(4) + B(2,iNode)*SW(2) + B(3,iNode)*SW(5)
                        IntLoadWaterLoc(IDof+3, iEntityID) = IntLoadWaterLoc(IDof+3,iEntityID) + B(1,iNode)*SW(6) + B(2,iNode)*SW(5) + B(3,iNode)*SW(3)
                      else if (NVECTOR == 2) then ! 2D case
                        if ( ISAXISYMMETRIC ) then
                          call GiveError("2D axisymmetric not implemented in [subroutine ConsolidationForces()] for double-point formulation.")   
                        end if
                        IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1,iEntityID) + B(1,iNode)*SW(1) + B(2,iNode)*SW(4)
                        IntLoadWaterLoc(IDof+2, iEntityID) = IntLoadWaterLoc(IDof+2,iEntityID) + B(1,iNode)*SW(4) + B(2,iNode)*SW(2)
                        ! axisymmetric to be added
                      end if  
                        
                    end do ! loop over element nodes
 
                  end if ! Only LIQUID material Point
                end do ! Loop over number of integration points per element IEl
!------------------------------------ END INTEGRATION POINT LOOP -----------------------------  
              end if ! 2 Phases
              
              !============================================================================================
              if(CalParams%NumberOfPhases==1) then ! 1 Phase

!//////////////////// External force calculation //////////////////////// 
                do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
                                    
                  ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              
                  if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then ! Only LIQUID material Point
  
                    PartLoad = 0.0
                    PartGravity = 0.0
                    PartGravityMix = 0.0   ! Not needed
              
                    iEntityID = 1
               
                    do INode = 1, ELEMENTNODES ! Element Nodes
                      NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                      IDof = ReducedDof(NodeID) 

                      PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBody
                      do I = 1, NVECTOR
                        GravityLoadWaterLoc(IDof+I, iEntityID) =  GravityLoadWaterLoc(IDof+I, iEntityID) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                      end do

                    end do ! Element Nodes
                  end if ! Only LIQUID material Point
                end do ! Loop MatPoints

!//////////////////// Internal force calculation ////////////////////////
                ! Determine number of integration points inside element
                if (IsParticleIntegration(IEl) ) then ! True - particle based integration, false - Gauss point based integration
                  NElemPart = NPartEle(IEl)  ! Number of particles in element
                else
                  NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
                end if
              
                call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element
        
                !------------------------------------ INTEGRATION POINT LOOP --------------------------------
                do Int = 1, NElemPart ! Loop over number of integration points per element IEl
              
                  ! Determine global ID of integration point 
                  IntGlo = GetParticleIndex(Int, IEl)
                  
                  if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then ! Only LIQUID material Point
              
                    ! Set the integration weight
                    if (IsParticleIntegration(IEl) ) then
                      ! use the integration weight of particle if it is partially filled, otherwise the weight of gauss point is used
                      WTN = Particles(IntGlo)%IntegrationWeight
                      if ( ISAXISYMMETRIC ) then ! the integration weight of the MP is not corrected due to axisymmetry
                        Position = GlobPosArray(IntGlo, 1) ! index 1 is readial directoin
                        ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                      end if
                    else
                      if ( ISAXISYMMETRIC ) then
                        Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                        ShapeValues(:) = GPShapeFunction(Int, :)
                        WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                      end if  
                    end if
            
                    SW(:) = SigmaEffArray(IntGlo,:)  * WTN

                    !get material point entity
                    iEntityID = 1

                    do iNode = 1,ELEMENTNODES ! loop over element nodes
                      nn = ElementConnectivities(iNode,iel) ! get global node number
                      IDof = ReducedDof(nn) ! global storage coordinate of x-val

                      if (NVECTOR == 3) then ! 3D case
                        IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1, iEntityID) + B(1,iNode)*SW(1) + B(2,iNode)*SW(4) + B(3,iNode)*SW(6)
                        IntLoadWaterLoc(IDof+2, iEntityID) = IntLoadWaterLoc(IDof+2, iEntityID) + B(1,iNode)*SW(4) + B(2,iNode)*SW(2) + B(3,iNode)*SW(5)
                        IntLoadWaterLoc(IDof+3, iEntityID) = IntLoadWaterLoc(IDof+3, iEntityID) + B(1,iNode)*SW(6) + B(2,iNode)*SW(5) + B(3,iNode)*SW(3)  
                      else if (NVECTOR == 2) then ! 2D case
                        if ( ISAXISYMMETRIC ) then
                          call GiveError("2D axisymmetric not implemented in [subroutine ConsolidationForces()] for double-point formulation.")   
                        end if
                        IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1, iEntityID) + B(1,iNode)*SW(1) + B(2,iNode)*SW(4)
                        IntLoadWaterLoc(IDof+2, iEntityID) = IntLoadWaterLoc(IDof+2, iEntityID) + B(1,iNode)*SW(4) + B(2,iNode)*SW(2)
                        ! axisymmetric to be added
                      end if
                      
                    end do ! loop over element nodes
                  end if ! Only LIQUID material Point
                end do ! Loop over number of integration points per element IEl
!------------------------------------ END INTEGRATION POINT LOOP -----------------------------  
              end if ! 1 Phase
            end if !2 Constituents
          end do ! Loop over elements

        end subroutine ConsolidationForces  
        
        subroutine ComputeWaterSpecificWeight()

        real(REAL_TYPE) :: Gravity, ConstDensityLiquid
        integer(INTEGER_TYPE) :: I

        do I = 1, Counters%NLayers ! loop over all material points in element
            if(trim(MatParams(I)%MATERIALTYPE)=='2-phase'.or.MatParams(I)%MaterialPhases=='2-phase') then
                ConstDensityLiquid = (MatParams(I)%DensityLiquid/1000)
            else if(trim(MatParams(I)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(I)%MaterialPhases=='1-phase-liquid') then
                ConstDensityLiquid = (MatParams(I)%DensityLiquid/1000)
            end if
        end do

        Gravity = (CalParams%GravityData%GAccel) * (CalParams%Multipliers%GravityCurrent)
        CalParams%GammaWater = (Gravity) * ConstDensityLiquid

        end subroutine ComputeWaterSpecificWeight

        subroutine UpdateMultipliersWaterForSpaceDependency(SpaceTimeExTLoadWater, NDof) 
           
        real(REAL_TYPE), dimension(Counters%N, 1), intent(out) :: SpaceTimeExTLoadWater
        integer(INTEGER_TYPE), dimension(Counters%NodTot), intent(in) :: NDof   
                 
        ! Local variables
        integer(INTEGER_TYPE), dimension(NDIM) :: NodeID
        integer(INTEGER_TYPE) :: I, J, K, ND
        real(REAL_TYPE) :: HydrHead, yNode, ExternalSuction, WaterAuxSpaceDependency !Gravity, ConstDensityLiquid,
        real(REAL_TYPE),dimension(Counters%NodTot*NVECTOR):: WaterAuxVector
        real(REAL_TYPE) :: GammaW    
        
        ! Preparing  auxiliary vectors for final extForce vector assembling
        WaterAuxSpaceDependency = 0.0
        WaterAuxVector = 0.0   
        
        GammaW = CalParams%GammaWater
        if (CalParams%Multipliers%HydraulicHeadCurrent /= 0) then
         HydrHead = CalParams%Multipliers%HydraulicHeadCurrent ! The current value of hydraulic head read from the file &
                                                               ! (corresponds to the actual multiplier)
        else
         HydrHead = 0.000001
        end if
        !!! Boundary nodes subjected to the hydraulic head conditions will have unmodified pressure value & 
        !!!(based on initial water pressure prescription) if they are located above current hydraulic head
        ExternalSuction = CalParams%InitialWaterPressure
               
        do I = 1, Counters%HydraulicHeadSides

            do J = 1, ELEMENTBOUNDARYNODES

                NodeID(J) = HydraulicHeadNodesConnectivities(J, I) !
                ND = NDof(NodeID(J))

                do K = 1, NDIM

                    if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
                        yNode = NodalCoordinates(NodeID(J), 2)  ! y-coordinate of a node belonging to loaded boundary
                        if (yNode <= HydrHead) then  ! check if the node is above/below the water level
                            WaterAuxSpaceDependency = - GammaW * (HydrHead - yNode)/HydrHead
                        else
                           WaterAuxSpaceDependency = 0.0
                       
                        end if
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        WaterAuxVector(ND + K) =  WaterAuxSpaceDependency
                    end if
                end do

            end do
        end do
        
        do I = 1, Counters%N
         SpaceTimeExTLoadWater(I,1) = HydraulicHeadLoadTotal(I,1) * WaterAuxVector(I)* HydrHead   
        end do     

       end subroutine UpdateMultipliersWaterForSpaceDependency  
       
                
        subroutine ConsolidationIntForces(IntLoadWaterLoc)
        !**********************************************************************
        !
        !    Function:  Calculation of the equivalent nodal internal forces 
        !               due to a given stress 
        !               InternalLD = Integral {BT*Sigma}.
        !
        !    O InternalLD : Internal load
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: IntLoadWaterLoc

          ! Local variables
          real(REAL_TYPE), dimension(NDIM, ELEMENTNODES) :: B
          integer(INTEGER_TYPE) :: I, J, IntGlo, IEl, Int, NN, NElemPart, iEntityID, iNode, IDof
          real(REAL_TYPE) :: WtN, Det, S, TotIntWeightLiqMP
          integer(INTEGER_TYPE) :: IAEl, ParticleIndex
          real(REAL_TYPE), dimension(NTENSOR) :: SW
          real(REAL_TYPE) :: Position
          real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)
          
          B = 0.0
          IntGlo = 0
          IEl = 0
          Int = 0
          J = 0
          NN = 0
          NElemPart = 0
          iEntityID = 0
          iNode = 0
          WtN = 0.0
          Det = 0.0
          S = 0.0
          TotIntWeightLiqMP = 0.0
          IAEl = 0
          ParticleIndex = 0
          SW = 0.0
          IntLoadWaterLoc = 0.0

          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)

            !========================================================================
            !==========================1 Layer=======================================
            !========================================================================
            
            if(NFORMULATION==1) then !1 Constituents

              ! Determine number of integration points inside element
              if (IsParticleIntegration(IEl) ) then ! True - particle based integration, false - Gauss point based integration
                NElemPart = NPartEle(IEl)  ! Number of particles in element
              else
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if
              
              call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element
        
              !------------------------------------ INTEGRATION POINT LOOP --------------------------------
              do Int = 1, NElemPart ! Loop over number of integration points per element IEl
              
                ! Determine global ID of integration point 
                IntGlo = GetParticleIndex(Int, IEl)
              
                if (IsParticleIntegration(IEl) ) then
                  ! use the integration weight of particle if it is partially filled, otherwise the weight of gauss point is used
                  WTN = Particles(IntGlo)%IntegrationWeight
                  if ( ISAXISYMMETRIC ) then
                    ! the integration weight of the MP is not corrected due to axisymmetry
                    Position = GlobPosArray(IntGlo, 1) ! index 1 is readial directoin
                    ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                  end if
                else
                  if ( ISAXISYMMETRIC ) then
                    Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                    ShapeValues(:) = GPShapeFunction(Int, :)
                    WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                  end if  
                end if
                
                S = Particles(IntGlo)%WaterPressure * WTN
        
                !get particle entity
                if (.not.CalParams%ApplyContactAlgorithm) then
                  iEntityID = 1
                else
                  iEntityID = EntityIDArray(IntGlo)
                end if

                do INode = 1,ELEMENTNODES  !loop through element nodes

                  nn = ElementConnectivities(iNode,iel) ! get global node number
                  IDof = ReducedDof(nn)   ! global storage coordinate of x-val
            
                  ! nodal x-load 
                  do I = 1, NVECTOR
                    IntLoadWaterLoc(IDof+I, iEntityID) = IntLoadWaterLoc(IDof+I, iEntityID) + B(I, INode) * S
                  end do
                  
                  if ( ISAXISYMMETRIC ) then
                    ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                    IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1, iEntityID) + S * ShapeValues(INode) / Position
                  end if
                                             
                end do !node loop             
              end do ! Loop over 1 .. NElemPart
              !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
            end if ! 1 Constituent
            
            !========================================================================
            !===========================2 Layers=====================================
            !========================================================================
            
            if(NFORMULATION==2) then !2 Constituents
              if(CalParams%NumberOfPhases>1) then ! 2 Phases
                  
                ! Determine number of integration points inside element
                if (IsParticleIntegration(IEl) ) then ! True - particle based integration, false - Gauss point based integration
                  NElemPart = NPartEle(IEl)  ! Number of particles in element
                else
                  NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
                end if
              
                TotIntWeightLiqMP = 0.0
                do J = 1, NPartEle(IEl)
                  ParticleIndex = GetParticleIndex(J, IEl)
                  if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                    TotIntWeightLiqMP = TotIntWeightLiqMP + Particles(ParticleIndex)%IntegrationWeight
                  end if
                end do
              
                call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element
        
                !------------------------------------ INTEGRATION POINT LOOP --------------------------------
                do Int = 1, NElemPart ! Loop over number of integration points per element IEl
              
                  ! Determine global ID of integration point 
                  IntGlo = GetParticleIndex(Int, IEl)
                  
                  if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then ! Only LIQUID material Point
              
                    ! Set the integration weight
                    WTN = Particles(IntGlo)%IntegrationWeight
                    if(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0) then
                      WTN = ElementSpace(IEl) / TotIntWeightLiqMP * Particles(IntGlo)%IntegrationWeight
                      if(TwoLayerData%Elements(IEl)%IsBoundaryOfLiquidElements) then ! Boundary of Liquid domain
                        if((TotIntWeightLiqMP/ElementSpace(IEl))<CalParams%RequiredDegreeOfFilling) then
                          WTN = Particles(IntGlo)%IntegrationWeight 
                        end if
                      end if
                    end if
                    
                    if ( ISAXISYMMETRIC ) then
                      if (IsParticleIntegration(IEl) ) then
                        ! the integration weight of the MP is not corrected due to axisymmetry
                        Position = GlobPosArray(IntGlo, 1) ! index 1 is readial directoin
                        ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                      else
                        Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                        ShapeValues(:) = GPShapeFunction(Int, :)
                        WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                      end if
                    end if

                    SW(:) = SigmaEffArray(IntGlo,:)  * WTN * (1.0 - TwoLayerData%ELEMENTS(IEl)%ConcentrationRatioSolidL)
              
                    !get material point entity
                    iEntityID = 1
                            
                    do iNode = 1,ELEMENTNODES ! loop over element nodes
                      nn = ElementConnectivities(iNode,iel) ! get global node number
                      IDof = ReducedDof(nn) ! global storage coordinate of x-val
               
                      if (NVECTOR == 3) then ! 3D case
                        IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1, iEntityID) + B(1,iNode)*SW(1) + B(2,iNode)*SW(4) + B(3,iNode)*SW(6)
                        IntLoadWaterLoc(IDof+2, iEntityID) = IntLoadWaterLoc(IDof+2, iEntityID) + B(1,iNode)*SW(4) + B(2,iNode)*SW(2) + B(3,iNode)*SW(5)
                        IntLoadWaterLoc(IDof+3, iEntityID) = IntLoadWaterLoc(IDof+3, iEntityID) + B(1,iNode)*SW(6) + B(2,iNode)*SW(5) + B(3,iNode)*SW(3)
                      else if (NVECTOR == 2) then ! 2D case
                        if ( ISAXISYMMETRIC ) then
                          call GiveError("2D axisymmetric not implemented in [subroutine ConsolidationForces()] for double-point formulation.")   
                        end if
                        IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1,iEntityID) + B(1,iNode)*SW(1) + B(2,iNode)*SW(4)
                        IntLoadWaterLoc(IDof+2, iEntityID) = IntLoadWaterLoc(IDof+2,iEntityID) + B(1,iNode)*SW(4) + B(2,iNode)*SW(2)
                      end if
                        
                    end do ! loop over element nodes
  
                  end if ! Only LIQUID material Point
                end do ! Loop over number of integration points per element IEl
                ! ------------------------------------ END INTEGRATION POINT LOOP -----------------------------  
              end if ! 2 Phases
              
              !============================================================================================
              
              if(CalParams%NumberOfPhases==1) then ! 1 Phase

                ! Determine number of integration points inside element
                if (IsParticleIntegration(IEl) ) then ! True - particle based integration, false - Gauss point based integration
                  NElemPart = NPartEle(IEl)  ! Number of particles in element
                else
                  NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
                end if
              
                call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element
        
                !------------------------------------ INTEGRATION POINT LOOP --------------------------------
                do Int = 1, NElemPart ! Loop over number of integration points per element IEl
              
                  ! Determine global ID of integration point 
                  IntGlo = GetParticleIndex(Int, IEl)
                  
                  if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then ! Only LIQUID material Point
              
                    ! Set the integration weight
                    if (IsParticleIntegration(IEl) ) then
                      ! use the integration weight of particle if it is partially filled,
                      ! otherwise the weight of gauss point is used
                      WTN = Particles(IntGlo)%IntegrationWeight
                    end if
                    
                    if ( ISAXISYMMETRIC ) then
                      if (IsParticleIntegration(IEl) ) then
                        ! the integration weight of the MP is not corrected due to axisymmetry
                        Position = GlobPosArray(IntGlo, 1) ! index 1 is readial directoin
                        ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                      else
                        Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                        ShapeValues(:) = GPShapeFunction(Int, :)
                        WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                      end if
                    end if

                    SW(:) = SigmaEffArray(IntGlo, :)  * WTN

                    !get material point entity
                    iEntityID = 1

                    do iNode = 1,ELEMENTNODES ! loop over element nodes
                      nn = ElementConnectivities(iNode,iel) ! get global node number
                      IDof = ReducedDof(nn) ! global storage coordinate of x-val
               
                      if (NVECTOR == 3) then ! 3D case
                        IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1, iEntityID) + B(1,iNode)*SW(1) + B(2,iNode)*SW(4) + B(3,iNode)*SW(6)
                        IntLoadWaterLoc(IDof+2, iEntityID) = IntLoadWaterLoc(IDof+2, iEntityID) + B(1,iNode)*SW(4) + B(2,iNode)*SW(2) + B(3,iNode)*SW(5)
                        IntLoadWaterLoc(IDof+3, iEntityID) = IntLoadWaterLoc(IDof+3, iEntityID) + B(1,iNode)*SW(6) + B(2,iNode)*SW(5) + B(3,iNode)*SW(3)
                      else if (NVECTOR == 2) then ! 2D case
                        if ( ISAXISYMMETRIC ) then
                          call GiveError("2D axisymmetric not implemented in [subroutine ConsolidationForces()] for double-point formulation.")   
                        end if
                        IntLoadWaterLoc(IDof+1, iEntityID) = IntLoadWaterLoc(IDof+1,iEntityID) + B(1,iNode)*SW(1) + B(2,iNode)*SW(4)
                        IntLoadWaterLoc(IDof+2, iEntityID) = IntLoadWaterLoc(IDof+2,iEntityID) + B(1,iNode)*SW(4) + B(2,iNode)*SW(2)
                      end if
                      
                    end do ! loop over element nodes
                    
                  end if ! Only LIQUID material Point
                end do ! Loop over number of integration points per element IEl
                !------------------------------------ END INTEGRATION POINT LOOP -----------------------------  
              end if ! 1 Phase
            end if !2 Constituents
            
          end do ! Loop over elements

      end subroutine ConsolidationIntForces
     
          
      subroutine ConsolidationExtForces(ExtLoadWaterLoc, &
                                        GravityLoadWaterLoc, &
                                        GravityLoadMixtureLoc)
        !**********************************************************************
        !
        !    Function:  Calculation of the equivalent nodal external forces
        !
        !   O ExternalLD : External load
        !   O GravityLD : Gravity load
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: ExtLoadWaterLoc, GravityLoadWaterLoc, GravityLoadMixtureLoc
          
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, iNode, IPart, IAEl, ParticleIndex, NodeID, IDof, iEntity, ILoadSystem
          real(REAL_TYPE), dimension(NVECTOR) :: PartLoad, PartGravity, PartGravityMix
          logical :: IsLoadOnMP

          IEl = 0
          iNode = 0
          IPart = 0
          IAEl = 0
          ParticleIndex = 0
          NodeID = 0
          iEntity = 0
          PartLoad = 0.0
          PartGravity = 0.0
          PartGravityMix = 0.0
          ExtLoadWaterLoc = 0.0
          GravityLoadWaterLoc = 0.0
          GravityLoadMixtureLoc = 0.0
          IsLoadOnMP =(Counters%NLoadedElementSidesWaterMatPoints+Counters%NLoadedElementSidesWaterMatPointsB)>0

                 do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            
            !========================================================================
            !==========================1 Layer=======================================
            !========================================================================
            
            if(NFORMULATION==1) then !1 Constituents

              do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
                
                ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                PartLoad = 0.0
                PartGravity = 0.0
                PartGravityMix = 0.0
                
                if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex) !CC entity to which particle belongs
                else
                  iEntity = 1
                endif  

                do INode = 1, ELEMENTNODES 
                  NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                  IDof = ReducedDof(NodeID)

                  if(IsLoadOnMP) then 
                   do ILoadSystem=1, Counters%NWaterLoadSystems   
                   !   ILoadSystem=1
                    PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExtWater(:,ILoadSystem)* CalParams%Multipliers%WaterACurrent(ILoadSystem)
                    do I = 1, NVECTOR
                     ExtLoadWaterLoc(IDof+I, iEntity) =  ExtLoadWaterLoc(IDof+I, iEntity) + PartLoad(I) 
                    end do
                   end do
                  end if

                  PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyWater

                  do I = 1, NVECTOR
                    GravityLoadWaterLoc(IDof+I, iEntity) =  GravityLoadWaterLoc(IDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                  end do
                  
                  PartGravityMix = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyMixed

                  do I = 1, NVECTOR
                    GravityLoadMixtureLoc(IDof+I, iEntity) =  GravityLoadMixtureLoc(IDof+I, iEntity) + PartGravityMix(I) * CalParams%Multipliers%GravityCurrent
                  end do

                end do
              end do
            end if ! 1 Constituent
            
            !========================================================================
            !==========================2 Layers======================================
            !========================================================================
           
            if(NFORMULATION==2) then !2 Constituents
            
              if(CalParams%NumberOfPhases>1) then ! 2 Phases

                do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
                  ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              
                  if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then ! Only LIQUID material Point
  
                    PartLoad = 0.0
                    PartGravity = 0.0
                    PartGravityMix = 0.0 ! Not needed
              
                    iEntity = 1
               
                    do INode = 1, ELEMENTNODES ! Element Nodes
                      NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                      IDof = ReducedDof(NodeID)

                     ! Assemble External Load on Liquid Material Point
                     if(IsLoadOnMP) then 
                       do ILoadSystem=1, Counters%NWaterLoadSystems  
                       !  ILoadSystem=1
                         PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExtWater(:,ILoadSystem) * CalParams%Multipliers%WaterACurrent(ILoadSystem)
                         do I = 1, NVECTOR
                           ExtLoadWaterLoc(IDof+I, iEntity) =  ExtLoadWaterLoc(IDof+I, iEntity) + PartLoad(I)
                         end do
                       end do
                     end if
                    
                     PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBody 

                     do I = 1, NVECTOR
                       GravityLoadWaterLoc(IDof+I, iEntity) =  GravityLoadWaterLoc(IDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                     end do
 
                    end do ! Element Nodes
                  end if ! Only LIQUID material Point
                end do ! Loop MatPoints
              end if ! 2 Phases
              
              !============================================================================================
              
              if(CalParams%NumberOfPhases==1) then ! 1 Phase

                do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
                  ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              
                  if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then ! Only LIQUID material Point
  
                    PartLoad = 0.0
                    PartGravity = 0.0
                    PartGravityMix = 0.0 ! Not needed
              
                    iEntity = 1
               
                    do INode = 1, ELEMENTNODES ! Element Nodes
                      NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                      IDof = ReducedDof(NodeID)

                      PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBody

                      do I = 1, NVECTOR
                       GravityLoadWaterLoc(IDof+I, iEntity) =  GravityLoadWaterLoc(IDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                      end do
                      
                    end do ! Element Nodes
                  end if ! Only LIQUID material Point
                end do ! Loop MatPoints
              end if ! 1 Phase
            
            end if !2 Constituents
          end do ! Loop over elements

      end subroutine ConsolidationExtForces

      subroutine CorrectWaterForPrescribedVelocity(InterfaceNodes,InterfaceNodeNormals,PrescrVelocity)
        !**********************************************************************
        !
        !    Function:  Corrects the water velocities or accelerations at nodes
        !               shared by a structure and soil.
        !
        !    InterfaceNodes : ID's of nodes shared by structure and soil
        !    GlobalPrescribedValues : Velocities or accelerations prescribed on the structure
        !    O  WaterArray : Water velocities or accelerations at nodes to be corrected
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE):: nn
          logical, dimension(:), intent(in) :: InterfaceNodes
          real(REAL_TYPE), dimension(:, :), intent(in) :: InterfaceNodeNormals
          real(REAL_TYPE), dimension(:,:), intent(in) :: PrescrVelocity
          ! Local variables
          integer(INTEGER_TYPE) :: IDoF, I, J
          real(REAL_TYPE), dimension(NVECTOR) :: Unormal, VelDiff, VelocityValues, PrescribedValues
          real(REAL_TYPE) :: DiffDotUn


          do nn = 1, Counters%NodTot

          if (.not.InterfaceNodes(nn)) CYCLE

           Unormal(:) = -InterfaceNodeNormals (:, nn)

           IDof = ReducedDof(nn)

           do J = 1, NVECTOR
              VelocityValues(J) = TotalVelocityWater(IDoF + J, SOFT_ENTITY)
              PrescribedValues(J) = PrescrVelocity(nn, J)
           end do

           !Calculate velocity difference between entity and system
           VelDiff = VelocityValues - PrescribedValues

           !VelDiff.dot.Unormal
           DiffDotUn = 0
           do I = 1, NVECTOR
             DiffDotUn = DiffDotUn + VelDiff(i) * Unormal(i)
           end do

           !if (DiffDotUn>0.0) then
              VelocityValues = VelocityValues - DiffDotUn*Unormal
           !end if

           do J = 1, NVECTOR
              TotalVelocityWater(IDoF + J, SOFT_ENTITY) = VelocityValues(J)
              AccelerationWater(IDoF + J, SOFT_ENTITY) =  (TotalVelocityWater(IDoF + J, SOFT_ENTITY) - TotalVelocityWaterPrevious(IDoF + J, SOFT_ENTITY)) /  CalParams%TimeIncrement
           end do
           end do

      end subroutine CorrectWaterForPrescribedVelocity


      subroutine ConsolidationForcesPorosity(ExtLoadWaterPorosityLoc, &
                                             GravityLoadWaterPorosityLoc, &
                                             IntLoadWaterPorosityLoc)
        !**********************************************************************
        !
        !    Function:  Calculation of the equivalent nodal forces due to
        !               a given stress InternalLD = Integral {BT*Sigma}.
        !
        !   O ExternalLDPorosity : Real External load
        !   O InternalLDPorosity : Real Internal load
        !   O GravityLDPorosity : Real Gravity load
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: ExtLoadWaterPorosityLoc, IntLoadWaterPorosityLoc, GravityLoadWaterPorosityLoc
          
          ! Local variables
          real(REAL_TYPE), dimension(NDIM, ELEMENTNODES) :: B
          integer(INTEGER_TYPE) :: I, IntGlo, IEl, Int, NN, NElemPart, iEntityID, iNode, IDof, ILoadSystem
          real(REAL_TYPE) :: WtN, Det, S
          real(REAL_TYPE) :: PartDegreeSat, PartPorosity
          integer(INTEGER_TYPE) :: IPart, IAEl, ParticleIndex, NodeID, iEntity
          real(REAL_TYPE), dimension(NVECTOR) :: PartLoad, PartGravity
          real(REAL_TYPE) :: Position
          real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)
          logical :: IsLoadOnMP

          ExtLoadWaterPorosityLoc = 0.0
          IntLoadWaterPorosityLoc = 0.0
          GravityLoadWaterPorosityLoc = 0.0
          IsLoadOnMP =(Counters%NLoadedElementSidesWaterMatPoints+Counters%NLoadedElementSidesWaterMatPointsB)>0

          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            
!//////////////////// External force calculation //////////////////////// 
            do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              PartLoad = 0.0
              PartGravity = 0.0
              
              PartPorosity = Particles(ParticleIndex)%Porosity
               if (CalParams%ApplyPartialSaturation) then
                   PartDegreeSat = Particles(ParticleIndex)%DegreeSaturation
               else    
                   PartDegreeSat = 1.0
               end if    

              if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex)   !CC entity to which particle belongs
              else
                  iEntity = 1
              endif  

              do INode = 1, ELEMENTNODES 
                NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                IDof = ReducedDof(NodeID)

                if(IsLoadOnMP) then 
                  do ILoadSystem=1, Counters%NWaterLoadSystems 
                  !  ILoadSystem=1
                    PartLoad = ShapeValuesArray(ParticleIndex,INode) *PartPorosity * PartDegreeSat* &
                        Particles(ParticleIndex)%FExtWater(:,ILoadSystem) * CalParams%Multipliers%WaterACurrent(ILoadSystem)
                    do I = 1, NVECTOR
                      ExtLoadWaterPorosityLoc(IDof+I, iEntity) = ExtLoadWaterPorosityLoc(IDof+I, iEntity) + PartLoad(I) 
                    end do
                  end do
                end if

                PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyWater * PartPorosity * PartDegreeSat

                do I = 1, NVECTOR
                  GravityLoadWaterPorosityLoc(IDof+I, iEntity) = GravityLoadWaterPorosityLoc(IDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                end do
                
              end do
            end do
            
!//////////////////// Internal force calculation ////////////////////////
            ! Determine number of integration points inside element
            if (IsParticleIntegration(IEl) ) then ! True - particle based integration, false - Gauss point based integration
              NElemPart = NPartEle(IEl)  ! Number of particles in element
            else
              NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
            end if
              
            call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element

            !------------------------------------ INTEGRATION POINT LOOP --------------------------------
            do Int = 1, NElemPart ! Loop over number of integration points per element IEl
              
              ! Determine global ID of integration point 
              IntGlo = GetParticleIndex(Int, IEl)
              
               PartPorosity = Particles(ParticleIndex)%Porosity
               if (CalParams%ApplyPartialSaturation) then
                   PartDegreeSat = Particles(ParticleIndex)%DegreeSaturation
               else    
                   PartDegreeSat = 1.0
               end if
              
              if (IsParticleIntegration(IEl) ) then
                  ! use the integration weight of particle if it is partially filled, otherwise the weight of gauss point is used
                  WTN = Particles(IntGlo)%IntegrationWeight
                  if ( ISAXISYMMETRIC ) then
                    ! the integration weight of the MP is not corrected due to axisymmetry
                    Position = GlobPosArray(IntGlo, 1) ! index 1 is readial directoin
                    ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                  end if
              else
                  if ( ISAXISYMMETRIC ) then
                    Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                    ShapeValues(:) = GPShapeFunction(Int, :)
                    WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                  end if  
              end if

              S = Particles(IntGlo)%WaterPressure * WTN * PartPorosity * PartDegreeSat

              !get particle entity
              if (.not.CalParams%ApplyContactAlgorithm) then
                iEntityID = 1
              else
                iEntityID = EntityIDArray(IntGlo) 
              end if

              do INode = 1,ELEMENTNODES  !loop through element nodes

                nn = ElementConnectivities(iNode,iel) ! get global node number
                IDof = ReducedDof(nn)   ! global storage coordinate of x-val

                do I = 1, NVECTOR
                  IntLoadWaterPorosityLoc(IDof+I, iEntityID) = IntLoadWaterPorosityLoc(IDof+I, iEntityID) + B(I, INode) * S
                end do

                if ( ISAXISYMMETRIC ) then
                  ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                  IntLoadWaterPorosityLoc(IDof+1, iEntityID) = IntLoadWaterPorosityLoc(IDof+1, iEntityID) + S * ShapeValues(INode) / Position
                end if
                
              end do !node loop
            end do ! Loop over 1 .. NElemPart
            !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
          end do ! Loop over elements

           if (IS3DCYLINDRIC) then ! Rotate vectors from global to local coordinate system
             do IEntity = 1, Counters%nEntity 
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadWaterPorosityLoc(:, IEntity), ExtLoadWaterPorosityLoc(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, IntLoadWaterPorosityLoc(:, IEntity), IntLoadWaterPorosityLoc(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, GravityLoadWaterPorosityLoc(:, IEntity), GravityLoadWaterPorosityLoc(:, IEntity))
             end do
           end if ! only for RotBoundCond

      end subroutine ConsolidationForcesPorosity
     
      
      subroutine ConsolidationIntForcesPorosity(IntLoadWaterPorosityLoc)
        !**********************************************************************
        !
        !    Function:  Calculation of the equivalent nodal internal forces 
        !               due to a given stress
        !               InternalLD = Integral {BT*Sigma}.
        !
        !    O InternalLDPorosity : Real Internal load
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: IntLoadWaterPorosityLoc

          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
          integer(INTEGER_TYPE) :: I, IntGlo, IEl, Int, NN, NElemPart, iEntityID, iNode, IDof
          real(REAL_TYPE) :: WtN, Det, S
          real(REAL_TYPE) :: PartDegreeSat, PartPorosity
          integer(INTEGER_TYPE) :: IAEl
          real(REAL_TYPE) :: Position
          real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)

          IntLoadWaterPorosityLoc = 0.0

    
          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)

            ! Determine number of integration points inside element
            if (IsParticleIntegration(IEl) ) then ! True - particle based integration, false - Gauss point based integration
              NElemPart = NPartEle(IEl)  ! Number of particles in element
            else
              NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
            end if

            call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element

            !------------------------------------ INTEGRATION POINT LOOP --------------------------------
            do Int = 1, NElemPart ! Loop over number of integration points per element IEl

              ! Determine global ID of integration point 
              IntGlo = GetParticleIndex(Int, IEl)
              
              PartPorosity = Particles(IntGlo)%Porosity
               if (CalParams%ApplyPartialSaturation) then
                   PartDegreeSat = Particles(IntGlo)%DegreeSaturation
               else    
                   PartDegreeSat = 1.0
               end if

              if (IsParticleIntegration(IEl) ) then
                  ! use the integration weight of particle if it is partially filled, otherwise the weight of gauss point is used
                  WTN = Particles(IntGlo)%IntegrationWeight
                  if ( ISAXISYMMETRIC ) then
                    ! the integration weight of the MP is not corrected due to axisymmetry
                    Position = GlobPosArray(IntGlo, 1) ! index 1 is readial directoin
                    ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                  end if
              else
                  if ( ISAXISYMMETRIC ) then
                    Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                    ShapeValues(:) = GPShapeFunction(Int, :)
                    WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                  end if  
              end if

              S = Particles(IntGlo)%WaterPressure * WTN * PartPorosity * PartDegreeSat

              !get particle entity
              if (.not.CalParams%ApplyContactAlgorithm) then
                iEntityID = 1
              else
                iEntityID = EntityIDArray(IntGlo) 
              end if

              do INode = 1,ELEMENTNODES  !loop through element nodes
                nn = ElementConnectivities(iNode,iel) ! get global node number
                IDof = ReducedDof(nn)    ! global storage coordinate of x-val

                do I = 1, NVECTOR
                  IntLoadWaterPorosityLoc(IDof+I, iEntityID) = IntLoadWaterPorosityLoc(IDof+I, iEntityID) + B(I, INode) * S
                end do
                
                if ( ISAXISYMMETRIC ) then
                    ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                    IntLoadWaterPorosityLoc(IDof+1, iEntityID) = IntLoadWaterPorosityLoc(IDof+1, iEntityID) + S * ShapeValues(INode) / Position
                end if

              end do !node loop
            end do ! Loop over 1 .. NElemPart
            !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
          end do ! Loop over elements
          
      end subroutine ConsolidationIntForcesPorosity

      
       subroutine ConsolidationExtForcesPorosity(ExtLoadWaterPorosityLoc, GravityLoadWaterPorosityLoc)
         !**********************************************************************
         !
         !    Function:  Calculation of the equivalent nodal external forces 
         !
         !   O ExternalLDPorosity : Real External load
         !   O GravityLDPorosity : Real Gravity load
         !
         !**********************************************************************

        implicit none

        real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: ExtLoadWaterPorosityLoc, GravityLoadWaterPorosityLoc

        ! Local variables
        integer(INTEGER_TYPE) :: I, IEl, iNode, IPart, IAEl, ParticleIndex, NodeID, IDof, iEntity, ILoadSystem
        real(REAL_TYPE), dimension(NVECTOR) :: PartLoad, PartGravity
        real(REAL_TYPE) :: PartDegreeSat, PartPorosity
        logical :: IsLoadOnMP

        ExtLoadWaterPorosityLoc = 0.0
        GravityLoadWaterPorosityLoc = 0.0
        IsLoadOnMP =(Counters%NLoadedElementSidesWaterMatPoints+Counters%NLoadedElementSidesWaterMatPointsB)>0

        do IAEl = 1, Counters%NAEl
           IEl = ActiveElement(IAEl)

           do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
             ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
             PartLoad = 0.0
             PartGravity = 0.0

             PartPorosity = Particles(ParticleIndex)%Porosity
             if (CalParams%ApplyPartialSaturation) then
                 PartDegreeSat = Particles(ParticleIndex)%DegreeSaturation
             else    
                 PartDegreeSat = 1.0
             end if
               
               
             if (CalParams%ApplyContactAlgorithm) then
               iEntity = EntityIDArray(ParticleIndex)   !CC entity to which particle belongs
             else
               iEntity = 1
             end if  

             do INode = 1, ELEMENTNODES 
                 NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                 IDof = ReducedDof(NodeID)

                if(IsLoadOnMP) then 
                  do ILoadSystem =1, Counters%NWaterLoadSystems
                    !ILoadSystem =1
                    PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExtWater(:,ILoadSystem) *  &
                             PartPorosity * PartDegreeSat*CalParams%Multipliers%WaterACurrent(ILoadSystem)
                   do I = 1, NVECTOR
                     ExtLoadWaterPorosityLoc(IDof+I, iEntity) = ExtLoadWaterPorosityLoc(IDof+I, iEntity) + PartLoad(I)
                   end do
                 end do
               end if

               PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyWater * PartPorosity * PartDegreeSat* CalParams%Multipliers%GravityCurrent

               do I = 1, NVECTOR
                 GravityLoadWaterPorosityLoc(IDof+I, iEntity) = GravityLoadWaterPorosityLoc(IDof+I, iEntity) + PartGravity(I) 
               end do

             end do  ! Loop over nodes
           end do ! Loop over particles

        end do ! Loop over elements

        if (IS3DCYLINDRIC) then ! Rotate vectors from global to local coordinate system
           do IEntity = 1, Counters%nEntity 
             call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoadWaterPorosityLoc(:, IEntity), ExtLoadWaterPorosityLoc(:, IEntity))
             call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, GravityLoadWaterPorosityLoc(:, IEntity), GravityLoadWaterPorosityLoc(:, IEntity))
           end do
        end if ! only for RotBoundCond

       end subroutine ConsolidationExtForcesPorosity
       
       subroutine DetectInfiltrationNodes(IsInfiltrationNode)   
       !**************************************************************************************
       !
       !    Function:  Determine which nodes are inside the infiltration zone and 
       !               on the free surface used to appli infiltration velocity at these nodes
       !
       !**************************************************************************************
       
       implicit none
       
       logical, dimension(:), intent(inout) :: IsInfiltrationNode
       integer(INTEGER_TYPE) :: IEl, ElemID, ISide, NodeID, J
       !integer(INTEGER_TYPE) :: LocalNode
       integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: LocalNodes
       real(REAL_TYPE) :: IsBoundarySide
       logical :: flag

       flag=.false.
       IsInfiltrationNode = .false.
       do IEl=1,Counters%NAEl !loop over active lement
         ElemID = ActiveElement(IEl)
           
         do ISide =1,ELEMENTSIDES !loop over side nodes
             IsBoundarySide = BoundaryElementSurface(ElemID,ISide,IsActiveElement, Counters%NEl) !Give 1 if ISide is adiacent to inactive element
             
             if (IsBoundarySide==1) then !side is on a free surface
               call DetermineSideNodes(ISide,LocalNodes) !Local ID (1, 2, 3...) of nodes at the boundary of ISide, TRI3,TRI6 = 2 nodes, TETRA = 3 nodes 
               
               do J=1,ELEMENTBOUNDARYNODES 
                 NodeID = ElementConnectivities(LocalNodes(J),ElemID) !Global name of boundary node
                 IsInfiltrationNode(NodeID) = IsInsideArea(NodeID,InfiltrationArea) !.true. if node is inside infiltration area
                 flag = flag.or.IsInfiltrationNode(NodeID)
               end do
             end if            
             
         end do
         
       end do
       CalParams%BoundaryConditions%ApplyInfiltrationRate = flag
       end subroutine DetectInfiltrationNodes
       

        logical function IsInsideArea(NodeID,BoundaryArea)
        !*************************************************************************************   
        !    function:     IsInsideInfiltrationArea        
        !
        !>   @param[in] NodeID : Node of the mesh
        !
        !>   @return IsInsideInfiltrationArea : Returns .true. if NodeID lies inside or on the bounding box defined by CalParams%InfiltrationArea
        !
        !    GlobPos: global coordinate of the nodes
        !    Infiltration area: boundary of the infiltration area Xmin Xmax/ Ymin Ymax / Zmin Zmax
        !*************************************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodeID 
          real(REAL_TYPE), dimension(:,:), intent(in) :: BoundaryArea !(min:max , NDim)
          real(REAL_TYPE), dimension(NDIM) :: GlobPos       
          integer :: I

          IsInsideArea = .TRUE.
          GlobPos = NodalCoordinates(NodeID,:)
          
          do I = 1, size(GlobPos)
            IsInsideArea = IsInsideArea .and. (GlobPos(I) >= BoundaryArea(1,I)) .and. (GlobPos(I) <= BoundaryArea(2,I))
          end do          
          
        end function IsInsideArea
        
       subroutine DetectSeepageFaceNodes(IsSeepageNode)
       !**************************************************************************************
       !
       !    Function: determine which nodes are inside the free seepage zone and on the free surface
       !              used to apply correct liquid velocity at these nodes
       !
       !**************************************************************************************
       
       implicit none
       
       logical, dimension(:), intent(inout) :: IsSeepageNode
       integer(INTEGER_TYPE) :: IEl, ElemID, ISide, NodeID, J
       !integer(INTEGER_TYPE) :: LocalNode
       integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: LocalNodes
       real(REAL_TYPE) :: IsBoundarySide
       real(REAL_TYPE) :: HydrHead, yNode
       logical :: flag
       
       HydrHead = CalParams%Multipliers%HydraulicHeadCurrent
       
       flag=.false.
       IsSeepageNode = .false.
       do IEl=1,Counters%NAEl !loop over active lement
         ElemID = ActiveElement(IEl)
           
         do ISide =1,ELEMENTSIDES !loop over side nodes
             IsBoundarySide = BoundaryElementSurface(ElemID,ISide,IsActiveElement, Counters%NEl) !Give 1 if ISide is adiacent to inactive element
             
             if (IsBoundarySide==1) then !side is on a free surface
                 call DetermineSideNodes(ISide,LocalNodes) !Local ID (1, 2, 3...) of nodes at the boundary of ISide, TRI3,TRI6 = 2 nodes, TETRA = 3 nodes

                 do J=1,ELEMENTBOUNDARYNODES
                
                     NodeID = ElementConnectivities(LocalNodes(J),ElemID) !Global name of boundary node
                     yNode = NodalCoordinates(NodeID, 2)
                     if ((IsHydraulicHeadNode(NodeID)) .and. (yNode > HydrHead)) then ! check if the node is above/below the water level
                         IsSeepageNode(NodeID) = IsInsideArea(NodeID,SeepageArea)
                          flag = flag.or.IsSeepageNode(NodeID)
                     else if (IsHydraulicHeadNode(NodeID) == .false.) then
                         IsSeepageNode(NodeID) = IsInsideArea(NodeID,SeepageArea)
                          flag = flag.or.IsSeepageNode(NodeID)  !.true. if node is inside seepage area
                     end if
                                      
                 end do
             end if
             
         end do
         
       end do
       CalParams%BoundaryConditions%ApplySeepageFace = flag
       end subroutine DetectSeepageFaceNodes
       
       subroutine DetectHydraulicHeadNode(IsHeadNode)
       !**************************************************************************************
       !
       !    Function: determine which nodes are inside the hydraulic head zone and on the free surface
       !
       !**************************************************************************************
       
       implicit none
       
       logical, dimension(:), intent(inout) :: IsHeadNode
       integer(INTEGER_TYPE) :: IEl, ElemID, ISide, NodeID, J
       !integer(INTEGER_TYPE) :: LocalNode
       integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: LocalNodes
       real(REAL_TYPE) :: IsBoundarySide
       logical :: flag
    
       flag=.false.
       IsHeadNode = .false.
       do IEl=1,Counters%NAEl !loop over active lement
         ElemID = ActiveElement(IEl)
           
         do ISide =1,ELEMENTSIDES !loop over side nodes
             IsBoundarySide = BoundaryElementSurface(ElemID,ISide,IsActiveElement, Counters%NEl) !Give 1 if ISide is adiacent to inactive element
             
             if (IsBoundarySide==1) then !side is on a free surface
                 call DetermineSideNodes(ISide,LocalNodes) !Local ID (1, 2, 3...) of nodes at the boundary of ISide, TRI3,TRI6 = 2 nodes, TETRA = 3 nodes

                 do J=1,ELEMENTBOUNDARYNODES

                     NodeID = ElementConnectivities(LocalNodes(J),ElemID) !Global name of boundary node

                     IsHeadNode(NodeID) = IsInsideArea(NodeID, HydrHeadArea)
                     flag = flag.or.IsHeadNode(NodeID)  !.true. if node is inside hydr head area

                 end do
             end if
             
         end do
         
       end do
    
       end subroutine DetectHydraulicHeadNode
       
       subroutine ApplyNodalInfiltrationRate(VelocityLiquid,VelocitySoil,PorosityDegSat,NodalUnitMassGradient)
       !*********************************************************************************************************************
       ! 
       !     Function: correct true liquid velocity and solid velocity considering the Infiltration Rate
       !               correction is applied only if NetInfiltrationRate is negative = water enters the system at 0 pressure
       !
       !*********************************************************************************************************************
                      
       implicit none
       !logical, dimension(:), intent(in) :: IsInfiltrationNode
       real(REAL_TYPE), dimension(:,:), intent(in):: PorosityDegSat, NodalUnitMassGradient
       
       real(REAL_TYPE), dimension(:,:), intent(inout):: VelocityLiquid,VelocitySoil
        
       real(REAL_TYPE), dimension(NVECTOR):: NetInfiltrationVelocity!InfiltrationRate
       real(REAL_TYPE), dimension(NDIM)::  Normal
       real(REAL_TYPE) :: NetInfiltrationRate, MassLiquid, MassSolid, PorosityDegSatNode
       real(REAL_TYPE) :: CorrectionLiquidVelocity, CorrectionSolidVelocity
       integer(INTEGER_TYPE):: INode, IDof, IDim, IEntity
   
        do INode=1,Counters%NodTot
          if (IsInfiltrationNode(INode)) then
            IDof = ReducedDof(INode) 
            do IEntity=1,Counters%NEntity
              NetInfiltrationRate = 0.0 
              do IDim =1,NVECTOR 
                 Normal(IDim) = NodalUnitMassGradient(INode,IDim)  
                 NetInfiltrationVelocity(IDim) = PorosityDegSat(IDof+IDim, IEntity)*&
                       (VelocityLiquid(IDof+IDim,IEntity)- VelocitySoil(IDof+IDim,IEntity))- &
                              InfiltrationRate(IDim)
                 NetInfiltrationRate = NetInfiltrationRate + NetInfiltrationVelocity(IDim)*Normal(IDim)
              end do
              !hardcoding to apply normal flux
              !  NetInfiltrationRate = NetInfiltrationRate-InfiltrationRate(2) 
              if(NetInfiltrationRate<0.0) then !apply correction
                 MassLiquid = LumpedMassNWater(IDof+1,IEntity)
                 MassSolid = LumpedMassDry(IDof+1,IEntity) 
                 PorosityDegSatNode = PorosityDegSat(IDof+1, IEntity)
                 CorrectionLiquidVelocity =   0.0
                 CorrectionSolidVelocity = 0.0
                  
                 if (PorosityDegSatNode/=0) then
                   CorrectionLiquidVelocity = -MassSolid*NetInfiltrationRate/ &
                       (PorosityDegSatNode*(MassLiquid+MassSolid)) !Normal component of the correction
                   CorrectionSolidVelocity = MassLiquid*NetInfiltrationRate/ &
                       (PorosityDegSatNode*(MassLiquid+MassSolid)) !Normal component of the correction
                  ! CorrectionLiquidVelocity = -NetInfiltrationRate/ &
                  !     (PorosityDegSatNode) !Normal component of the correction
                 else
                    CorrectionLiquidVelocity =   0.0
                    CorrectionSolidVelocity = 0.0
                 end if
                   
                 do IDim=1,NVECTOR
                   VelocityLiquid(IDof+IDim,IEntity) = VelocityLiquid(IDof+IDim,IEntity) + CorrectionLiquidVelocity*Normal(IDim)
                   VelocitySoil(IDof+IDim,IEntity) = VelocitySoil(IDof+IDim,IEntity) + CorrectionSolidVelocity*Normal(IDim)
                 end do
              else
                  !no correction is necessary, water flows out of the soil at 0 pressure (ponding or seepage surface)
              end if
                
            end do 
            
            elseif (IsSeepageFaceNode(INode)) then
            IDof = ReducedDof(INode) 
            do IEntity=1,Counters%NEntity
              NetInfiltrationRate = 0.0 
              do IDim =1,NVECTOR 
                 Normal(IDim) = NodalUnitMassGradient(INode,IDim)  
                 NetInfiltrationVelocity(IDim) = PorosityDegSat(IDof+IDim, IEntity)*&
                       (VelocityLiquid(IDof+IDim,IEntity)- VelocitySoil(IDof+IDim,IEntity))
                 NetInfiltrationRate = NetInfiltrationRate + NetInfiltrationVelocity(IDim)*Normal(IDim)
              end do
                
              if(NetInfiltrationRate<0.0) then !apply correction because water is flowing inside soil 
                 MassLiquid = LumpedMassNWater(IDof+1,IEntity)
                 MassSolid = LumpedMassDry(IDof+1,IEntity) 
                 PorosityDegSatNode = PorosityDegSat(IDof+1, IEntity)
                 CorrectionLiquidVelocity =   0.0
                 CorrectionSolidVelocity = 0.0
                  
                 if (PorosityDegSatNode/=0) then
                   CorrectionLiquidVelocity = -MassSolid*NetInfiltrationRate/ &
                       (PorosityDegSatNode*(MassLiquid+MassSolid)) !Normal component of the correction
                   CorrectionSolidVelocity = MassLiquid*NetInfiltrationRate/ &
                       (PorosityDegSatNode*(MassLiquid+MassSolid)) !Normal component of the correction
                 else
                    CorrectionLiquidVelocity =   0.0
                    CorrectionSolidVelocity = 0.0
                 end if
                   
                 do IDim=1,NVECTOR
                   VelocityLiquid(IDof+IDim,IEntity) = VelocityLiquid(IDof+IDim,IEntity) + CorrectionLiquidVelocity*Normal(IDim)
                   VelocitySoil(IDof+IDim,IEntity) = VelocitySoil(IDof+IDim,IEntity) + CorrectionSolidVelocity*Normal(IDim)
                 end do
              else
                  !no correction is necessary, water flows out of the soil at 0 pressure (ponding or seepage surface)
              end if
                
            end do  
          end if
        end do
        if (CalParams%ApplyFixedSolidSkeleton) then ! Solid is fixed, compresibility of solid skeleton is neglected.
             VelocitySoil = 0.0
         end if
        
        end subroutine ApplyNodalInfiltrationRate
        
        subroutine MapMPPorosityDegSatToNodes(LumpedNodalPorosityDegSat)
        !*****************************************************************************************************************************
        !
        !        Function: map MP porosity and degree of saturation to the nodes of the mesh
        !
        !****************************************************************************************************************************
        implicit none
        real(REAL_TYPE), dimension(:,:), intent(inout):: LumpedNodalPorosityDegSat
        real(REAL_TYPE), dimension(Counters%N,Counters%NEntity) :: SumIntegrationWeight
        real(REAL_TYPE) :: Porosity, DegreeSaturation, PorosityDegSatEntry, WTN
        integer(INTEGER_TYPE) :: IAel, IEl, IPart, ParticleIndex, iEntity, INode, IDof, I, NodeID
        
        SumIntegrationWeight= 0.0
        LumpedNodalPorosityDegSat = 0.0
        do IAEl = 1, Counters%NAEl        ! Loop over all elements
            IEl = ActiveElement(IAEl)
            
            do IPart = 1, NPartEle(IEl)  ! Loop over all particles in element
              
              ParticleIndex = GetParticleIndex(IPart, IEl)    ! Get the particle ID      
              Porosity = Particles(ParticleIndex)%Porosity
              DegreeSaturation = Particles(ParticleIndex)%DegreeSaturation              
              WTN = Particles(ParticleIndex)%IntegrationWeight

                
              if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex) ! get entity to which particle belongs
              else
                iEntity = 1
              end if
              
              do INode = 1, ELEMENTNODES                           ! loop over element nodes
                  NodeID = ElementConnectivities(INode, IEl)       ! Global node ID
                  IDof = ReducedDof(NodeID)                        ! Global storage coordinate of x-val
                    
                  PorosityDegSatEntry =  Porosity * DegreeSaturation * WTN           !calculate contribution
                  do I = 1, NVECTOR
                    !add contribution
                    LumpedNodalPorosityDegSat(IDof+I, iEntity) = LumpedNodalPorosityDegSat(IDof+I, iEntity) + PorosityDegSatEntry
                    SumIntegrationWeight(IDof+I, iEntity) = SumIntegrationWeight(IDof+I, iEntity) + WTN
                  end do
              end do
            end do
        end do
        
        do iEntity=1,Counters%NEntity
          do IDof =1, Counters%N
              
            if (SumIntegrationWeight(IDof, iEntity) /= 0.0) then
             LumpedNodalPorosityDegSat(IDof, iEntity) = LumpedNodalPorosityDegSat(IDof, iEntity)/SumIntegrationWeight(IDof, iEntity)
            end if

          end do
        end do
            
        end subroutine MapMPPorosityDegSatToNodes
    
        subroutine DynUpdateParticleDegreeOfSaturationVanGenuchten(ParticleIndex,ISet,Sr)
        !**********************************************************************
        !
        !    Function:  Updates Degree of Saturation (Sr) of the particle 
        !               The expression of the Degree of Saturation is the RETENTION CURVE, based
        !               on the VAN GENUCHTEN MODEL (1980)
        !
        !**********************************************************************

        implicit none
        ! Local variables
        integer(INTEGER_TYPE), intent (in):: ParticleIndex, ISet
        real(REAL_TYPE), intent(inout) ::  Sr
        real(REAL_TYPE) ::  Pw, Pg, N, e
        real(REAL_TYPE) ::  Suc, Smin, Smax
        real(REAL_TYPE) ::  L, P, P0
        real(REAL_TYPE) ::  n00, n1, n2, n3, n4, n5, n6, n7, n8

        Pw = Particles(ParticleIndex)%WaterPressure
        Pg = Particles(ParticleIndex)%GasPressure
        N = Particles(ParticleIndex)%Porosity

        n00 = 1.0d0
        n1 = 0.625d0
        n2 = 0.2358d0
        n3 = 1.256d0
        n4 = 0.0019406d0
        n5 = 0.05d0
        n6 = 360.0d0
        n7 = 374.15d0
        n8 = 647.3d0

        Suc = Pw-Pg  !If Suction is > 0 (here) --> Sr < 1

        !If Smin=0.0 and Smax=1.0, the Effective Degree of Saturation (Se) is
        ! equivalent to the "liquid phase" Degree of Saturation (Sr)
        Smin = MatParams(ISet)%Smin_SWRC

        !If Smin=0.0 and Smax=1.0, the Effective Degree of Saturation (Se) is
        ! equivalent to the "liquid phase" Degree of Saturation (Sr)   
        Smax = MatParams(ISet)%Smax_SWRC

        ! P0=15000 kPa (mudstone model) / P0=18000 kPa (bentonite model) / P0=15 kPa (sand) 
        P0 = MatParams(ISet)%P0_SWRC

        ! L=0.36 (mudstone model) / P0=0.38 (bentonite model)
        L = MatParams(ISet)%Lambda_SWRC

        P = P0

        if (Suc<=0.0) then
            Sr = 1.0d0 !Saturated

        else
            e =  1.0d0/(1.0d0-L)
            Sr = Smin + (Smax-Smin)*((n00 + ((Suc/P)**(e)))**(-L))
        end if
        
         Particles(ParticleIndex)%DegreeSaturation = Sr
         
        end subroutine DynUpdateParticleDegreeOfSaturationVanGenuchten


        subroutine DynUpdateParticleDegreeOfSaturationLinear(ParticleIndex,ISet,Sr)
        !****************************************************************************************
        !
        !    Function:  Updates Degree of Saturation (Sr) of the particle 
        !               The expression of the Degree of Saturation is the RETENTION CURVE, based
        !               on a LINEAR Model
        !
        !****************************************************************************************
        
        implicit none
        ! Local variables
        integer(INTEGER_TYPE), intent (in):: ParticleIndex, ISet
        real(REAL_TYPE) ::  Pw, Pg, Suc, av
        real(REAL_TYPE), intent(inout) ::  Sr

            Pw = Particles(ParticleIndex)%WaterPressure
            Pg = Particles(ParticleIndex)%GasPressure
          
            Suc = Pw-Pg  !If Suction is > 0 (here) --> Sr < 1
            av = MatParams(ISet)%av_SWRC   ! m2/KN
            
            if (Suc<=0.0) then
                Sr=1.0d0 !Saturated
            else
                Sr = 1 - av*Suc
            end if

            Particles(ParticleIndex)%DegreeSaturation = Sr
            
            if (Particles(ParticleIndex)%DegreeSaturation<0.0) then
                Particles(ParticleIndex)%DegreeSaturation = 0.0
            end if
              
        end subroutine DynUpdateParticleDegreeOfSaturationLinear

        
          
      end module ModMPMDYN2PhaseSP