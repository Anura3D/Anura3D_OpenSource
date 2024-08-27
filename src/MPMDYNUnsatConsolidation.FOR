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
	  
	  
	  module ModMPMDYNUnsatConsolidation
      !**********************************************************************
      !
      !    Function:  Contains the routine for forming the unsaturated consolidation matrices
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 9262 $
      !     $Date: 2021-04-29 05:31:48 +0200 (do, 29 apr 2021) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluationTETRA
      use ModElementEvaluationTRI
      use ModMPMData
      use ModMeshInfo
      use ModRotBoundCond
      
      implicit none

        real(REAL_TYPE), dimension(:, :), allocatable :: LumpedMassGas ! Lumped mass vector of gas
        real(REAL_TYPE), dimension(:, :), allocatable :: LumpedMassNGas ! Lumped mass vector of gas * porosity  
        real(REAL_TYPE), dimension(:, :), allocatable :: GravityLoadGas ! Nodal load array of water gravity load
        real(REAL_TYPE), dimension(:, :), allocatable :: GravityLoadGasPorosityDegreeSat ! Nodal load array of gas gravity load taking into accont porosity and the degree of saturation
        real(REAL_TYPE), dimension(:, :), allocatable :: GravityLoadWaterPorosityDegreeSat ! Nodal load array of water gravity load taking into accont porosity and the degree of saturation
        real(REAL_TYPE), dimension(:, :), allocatable :: ConductivityGasMatrix ! Conductivity of gas matrix
        real(REAL_TYPE), dimension(:, :), allocatable :: ConductivityGasMatrixPorosityDegreeSat ! Conductivity of gas matrix * porosity * degree of saturation
        real(REAL_TYPE), dimension(:, :), allocatable :: QVWg ! Nodal drag force
        real(REAL_TYPE), dimension(:, :), allocatable :: QVWgPorosityDegreeSat ! Nodal drag force * Porosity * degree of saturation
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadGas ! Internal load vector correspond to gas
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadGasPorosityDegreeSat ! Internal load vector correspond to gas taking into accont porosity and the degree of saturation
        real(REAL_TYPE), dimension(:, :), allocatable :: IntLoadWaterPorosityDegreeSat ! Internal load vector correspond to water taking into accont porosity and the degree of saturation
        real(REAL_TYPE), dimension(:, :), allocatable :: BishopIntLoad ! 
        real(REAL_TYPE), dimension(:, :), allocatable :: AccelerationGas ! Nodal acceleration for gas
        real(REAL_TYPE), dimension(:, :), allocatable :: ExtLoadGas ! Nodal load array of external gas load
        real(REAL_TYPE), dimension(:, :), allocatable :: ExtLoadGasPorosityDegreeSat ! External load vector correspond to gas taking into accont porosity and the degree of saturation
        real(REAL_TYPE), dimension(:, :), allocatable :: ExtLoadWaterPorosityDegreeSat ! External load vector correspond to water taking into accont porosity and the degree of saturation
        real(REAL_TYPE), dimension(:, :), allocatable :: ExtLoadGasTotal ! Nodal load array of external gas load (Total pressure applied on nodes)
        
      contains ! Routines of this module

      
        subroutine InitialiseThreePhaseData()
        !**********************************************************************
        !
        ! Function: Contains code for initialising data for three-phase calculation and two-phase with suction
        !
        !**********************************************************************
        implicit none
        
          call DestroyThreePhaseData()
          
          call InitialiseThreePhaseArrays()
      
        end subroutine InitialiseThreePhaseData


        subroutine DestroyThreePhaseData()
        !**********************************************************************
        !
        ! Function: Deallocates the arrays used in this module
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError

          if (allocated(LumpedMassGas)) then
            deallocate(LumpedMassGas, stat = IError)
          end if

          if (allocated(LumpedMassNGas)) then
            deallocate(LumpedMassNGas, stat = IError)
          end if

          if (allocated(GravityLoadGas)) then
            deallocate(GravityLoadGas, stat = IError)
          end if

          if (allocated(GravityLoadGasPorosityDegreeSat)) then
            deallocate(GravityLoadGasPorosityDegreeSat, stat = IError)
          end if

          if (allocated(GravityLoadWaterPorosityDegreeSat)) then
            deallocate(GravityLoadWaterPorosityDegreeSat, stat = IError)
          end if

          if (allocated(ConductivityGasMatrix)) then
            deallocate(ConductivityGasMatrix, stat = IError)
          end if
          
          if (allocated(ConductivityGasMatrixPorosityDegreeSat)) then
            deallocate(ConductivityGasMatrixPorosityDegreeSat, stat = IError)
          end if

          if (allocated(QVWg)) then
            deallocate(QVWg, stat = IError)
          end if
          if (allocated(QVWgPorosityDegreeSat)) then
            deallocate(QVWgPorosityDegreeSat, stat = IError)
          end if
          
          if (allocated(IntLoadGas)) then
            deallocate(IntLoadGas, stat = IError)
          end if
          
          if (allocated(IntLoadGasPorosityDegreeSat)) then
            deallocate(IntLoadGasPorosityDegreeSat, stat = IError)
          end if

          if (allocated(IntLoadWaterPorosityDegreeSat)) then
            deallocate(IntLoadWaterPorosityDegreeSat, stat = IError)
          end if

          if (allocated(BishopIntLoad)) then
            deallocate(BishopIntLoad, stat = IError)
          end if

          if (allocated(AccelerationGas)) then
            deallocate(AccelerationGas, stat = IError)
          end if

          if (allocated(ExtLoadGas)) then
            deallocate(ExtLoadGas, stat = IError)
          end if

          if (allocated(ExtLoadGasPorosityDegreeSat)) then
            deallocate(ExtLoadGasPorosityDegreeSat, stat = IError)
          end if

          if (allocated(ExtLoadWaterPorosityDegreeSat)) then
            deallocate(ExtLoadWaterPorosityDegreeSat, stat = IError)
          end if

          if (allocated(ExtLoadGasTotal)) then
            deallocate(ExtLoadGasTotal, stat = IError)
          end if

        end subroutine DestroyThreePhaseData


        subroutine InitialiseThreePhaseArrays()
        !**********************************************************************
        !
        ! Function: Initialise the arrays related to three-phase calculation
        !
        !**********************************************************************

        implicit none

          ! local variables
          integer(INTEGER_TYPE) :: IError

          if (CalParams%NumberOfPhases==3) then
            allocate(LumpedMassGas(Counters%N, Counters%NEntity), stat = IError)
            allocate(LumpedMassNGas(Counters%N, Counters%NEntity), stat = IError)
            allocate(GravityLoadGas(Counters%N, Counters%NEntity), stat = IError)
            allocate(GravityLoadGasPorosityDegreeSat(Counters%N, Counters%NEntity), stat = IError)
            allocate(GravityLoadWaterPorosityDegreeSat(Counters%N, Counters%NEntity), stat = IError)
            allocate(ConductivityGasMatrix(Counters%N, Counters%NEntity), stat = IError)
            allocate(ConductivityGasMatrixPorosityDegreeSat(Counters%N, Counters%NEntity), stat = IError)
            allocate(QVWg(Counters%N, Counters%NEntity), stat = IError)
            allocate(QVWgPorosityDegreeSat(Counters%N, Counters%NEntity), stat = IError)
            allocate(IntLoadGas(Counters%N, Counters%NEntity), stat = IError)
            allocate(IntLoadGasPorosityDegreeSat(Counters%N, Counters%NEntity), stat = IError)
            allocate(IntLoadWaterPorosityDegreeSat(Counters%N, Counters%NEntity), stat = IError)
            allocate(BishopIntLoad(Counters%N, Counters%NEntity), stat = IError)
            allocate(AccelerationGas(Counters%N, Counters%NEntity), stat = IError)
            allocate(ExtLoadGas(Counters%N, Counters%NEntity), stat = IError)
            allocate(ExtLoadGasPorosityDegreeSat(Counters%N, Counters%NEntity), stat = IError)
            allocate(ExtLoadWaterPorosityDegreeSat(Counters%N, Counters%NEntity), stat = IError)
            allocate(ExtLoadGasTotal(Counters%N, Counters%NEntity), stat = IError)
          elseif ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation)) then
            allocate(BishopIntLoad(Counters%N, Counters%NEntity), stat = IError)
          else    
            allocate(LumpedMassGas(1, 1), stat = IError)
            allocate(LumpedMassNGas(1, 1), stat = IError)
            allocate(GravityLoadGas(1, 1), stat = IError)
            allocate(GravityLoadGasPorosityDegreeSat(1, 1), stat = IError)
            allocate(GravityLoadWaterPorosityDegreeSat(1, 1), stat = IError)
            allocate(ConductivityGasMatrix(1, 1), stat = IError)
            allocate(ConductivityGasMatrixPorosityDegreeSat(1, 1), stat = IError)
            allocate(QVWg(1, 1), stat = IError)
            allocate(QVWgPorosityDegreeSat(1, 1), stat = IError)
            allocate(IntLoadGas(1, 1), stat = IError)
            allocate(IntLoadGasPorosityDegreeSat(1, 1), stat = IError)
            allocate(IntLoadWaterPorosityDegreeSat(1, 1), stat = IError)
            allocate(BishopIntLoad(1, 1), stat = IError)
            allocate(AccelerationGas(1, 1), stat = IError)
            allocate(ExtLoadGas(1, 1), stat = IError)
            allocate(ExtLoadGasPorosityDegreeSat(1, 1), stat = IError)
            allocate(ExtLoadWaterPorosityDegreeSat(1, 1), stat = IError)
            allocate(ExtLoadGasTotal(1, 1), stat = IError)
          end if
          
          LumpedMassGas = 0.0
          LumpedMassNGas = 0.0
          GravityLoadGas = 0.0
          GravityLoadGasPorosityDegreeSat = 0.0
          GravityLoadWaterPorosityDegreeSat = 0.0
          ConductivityGasMatrix = 0.0
          ConductivityGasMatrixPorosityDegreeSat = 0.0
          QVWg = 0.0
          QVWgPorosityDegreeSat = 0.0
          IntLoadGas = 0.0
          IntLoadGasPorosityDegreeSat = 0.0
          IntLoadWaterPorosityDegreeSat = 0.0
          BishopIntLoad = 0.0
          AccelerationGas = 0.0
          ExtLoadGas = 0.0
          ExtLoadGasPorosityDegreeSat = 0.0
          ExtLoadWaterPorosityDegreeSat = 0.0
          ExtLoadGasTotal = 0.0

        end subroutine InitialiseThreePhaseArrays


        subroutine FormConsolidationMatricesGas() 
        !**********************************************************************
        !
        !    Function:  To extrapolate loads and masses from particles to nodes using 
        !               the shape function values evaluated at the particles local 
        !               position.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
      implicit none
      
      integer(INTEGER_TYPE) :: IEntity

          call FormMatricesGas(LumpedMassGas, &
                               LumpedMassNGas, &
                               ConductivityGasMatrix, &
                               ConductivityGasMatrixPorosityDegreeSat)

          call ConsolidationForcesGas(ExtLoadGas, &
                                      GravityLoadGas, &
                                      IntLoadGas)
          
          if (Counters%NLoadedElementSidesGasNodes > 0) then
            ExtLoadGas = ExtLoadGas + ExtLoadGasTotal * CalParams%Multipliers%GasACurrent
          end if
     
          if (IS3DCYLINDRIC) then ! rotation is needed
            do IEntity = 1, Counters%nEntity
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, &
                           ExtLoadGas(:, IEntity), ExtLoadGas(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, &
                           GravityLoadGas(:, IEntity), GravityLoadGas(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, &
                           IntLoadGas(:, IEntity), IntLoadGas(:, IEntity))
            end do
          end if ! rotation
      
         end subroutine FormConsolidationMatricesGas

        
         subroutine MapGasMomentumFromParticlesToNodes(Momentum)
        !**********************************************************************
        !
        !    Function:  To map gas momentum from particles to the grid points (nodes)
        !
        !
        ! O  Momentum : Nodal gas momentum vector, the output of this subroutine
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum  
          ! Local variables
          integer(INTEGER_TYPE) :: I, IAEl, IEl, IPart, INode, IDof, ParticleIndex, NodeID, iEntity
          real(REAL_TYPE), dimension (NVECTOR) :: ParticleVelocity
          
          Momentum = 0.0
            
          do IAEl = 1, Counters%NAEl                                      !loop over all elements
            IEl = ActiveElement(IAEl)
            do IPart = 1, NPartEle(IEl)                                 !loop over all particles in element
                ParticleIndex = GetParticleIndex(IPart,IEl)             !get the particle ID
                ParticleVelocity = VelocityGasArray(ParticleIndex,:)    !get particle water velocity vector
                
                if (CalParams%ApplyContactAlgorithm) then
                    iEntity = EntityIDArray(ParticleIndex)  !entity to which particle belongs
                else
                    iEntity = 1
                end if    
                  
                do INode = 1, ELEMENTNODES                           !loop over nodes
                
                      NodeID = ElementConnectivities(INode, IEl)                         !Global node ID
                      IDof = ReducedDof(NodeID)
                    
                      do I = 1, NVECTOR
                        Momentum(IDof+I,iEntity)= Momentum(IDof+I,iEntity) + Particles(ParticleIndex)%MassGas * ShapeValuesArray(ParticleIndex,INode) * ParticleVelocity(I)
                      end do  
               
               end do !Loop over nodes
            end do !Loop over particles
          end do !elements
         
         end subroutine MapGasMomentumFromParticlesToNodes
        
        subroutine GetNodalGasVelocityFromNodalGasMomentum(Momentum)
        !**********************************************************************
        !
        !    Function:  To calculate the nodal gas velocities from nodal gas mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Momentum
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          
          do IDOF = 1, Counters%N                                
            do J = 1, Counters%nEntity !loop through all entities    
                if(LumpedMassGas(IDOF,J)/=0) then
                    TotalVelocityGas(IDOF,J) = ( Momentum(IDOF,J) / LumpedMassGas(IDOF,J) ) * PboundaryGas(IDOF)
                else
                    TotalVelocityGas(IDOF,J) = 0.0
                end if
            end do
          end do
          
        end subroutine GetNodalGasVelocityFromNodalGasMomentum
        
        subroutine GetQVWgArray()
        
        !**********************************************************************
        !
        !    Function:  To calculate the QVWg array
        !
        ! O  QVWg : The output of this subroutine
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          real(REAL_TYPE) :: VminusWg
          
          QVWg = 0.0
          
          do IDOF = 1, Counters%N       !Loop over all degrees of freedom
            do J = 1, Counters%nEntity  !loop over all entities
              VminusWg = TotalVelocitySoil(IDOF,J) - TotalVelocityGas(IDOF,J)
              QVWg(IDOF,J) =  ConductivityGasMatrix(IDOF,J) * VminusWg
              QVWgPorosityDegreeSat(IDOF,J) =  ConductivityGasMatrixPorosityDegreeSat(IDOF,J) * VminusWg
            end do
          end do
            
        end subroutine GetQVWgArray
               
        subroutine GetGasInertiaArray(GasInertia)
        
        !**********************************************************************
        !
        !    Function:  To calculate the GasInertia array
        !
        ! O  GasInertia : The output of this subroutine
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(inout) :: GasInertia    
                   
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          
          GasInertia = 0.0
          
          do IDOF = 1, Counters%N       !Loop over all degrees of freedom
            do J = 1, Counters%nEntity  !loop over all entities
              GasInertia(IDOF,J) = LumpedMassNGas(IDOF,J) * AccelerationGas(IDOF,J)
            end do
          end do
            
        end subroutine GetGasInertiaArray
        
        subroutine CalculateGasIncrementalNodalAcceleration(RateofMomentum)
        !**********************************************************************
        !
        !    Function:  To calculate the incremental nodal accelerations (gas phase)
        !
        !    RateofMomentum: The arry which stores the rate of momentum (gas) 
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
         real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: RateofMomentum           
     
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          
           do IDOF = 1, Counters%N       !Loop over all degrees of freedom
            do J = 1, Counters%nEntity  !loop over all entities
               if(LumpedMassGas(IDOF,J)/=0) then
                 AccelerationGas(IDOF,J)=( RateofMomentum(IDOF,J)/ LumpedMassGas(IDOF,J) ) * PboundaryGas(IDOF)
                else
                  AccelerationGas(IDOF,J) = 0.0
                end if
              end do
            end do
                    
         end subroutine CalculateGasIncrementalNodalAcceleration
        
         subroutine UpdateParticleGasVelocityAndMapMomentumG(Momentum)
                              
        !**********************************************************************
        !
        !    Function:  To update particles total water velocities and Accelerations
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum  
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: IDof
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleIncrementalVelocity
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: NodAcc
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) :: Time, PartilceMass
          
          Momentum = 0.0
           
          Time = CalParams%TimeIncrement

          NoEn = Counters%nEntity

          do IAEl = 1, Counters%NAEl                                      ! Loop over all elements
            IEl = ActiveElement(IAEl)
            
            do I = 1, NVECTOR
              IDof(I, 1:ELEMENTNODES) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
              NodAcc(I, 1:ELEMENTNODES, 1:Counters%nEntity) = AccelerationGas(IDof(I, 1:ELEMENTNODES), 1:NoEn)
            end do
            
              do IPart = 1, NPartEle(IEl)                                 ! Loop over all particles in element
                ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                ParticleVelocity =  VelocityGasArray(ParticleIndex,:)
                PartilceMass = Particles(ParticleIndex)%MassGas
                ParticleShape = ShapeValuesArray(ParticleIndex,:)
                ParticleIncrementalVelocity = 0.0
                ParticleAcceleration = 0.0
               if (CalParams%ApplyContactAlgorithm) then
                 iEntity = EntityIDArray(ParticleIndex) 
                else
                 iEntity = 1
               end if 
                  
               do INode = 1, ELEMENTNODES                         ! loop over element nodes
  
                 do I = 1, NVECTOR  
                   ! Particle accelerations
                   ParticleAcceleration(I) = ParticleAcceleration(I) + ParticleShape(INode) * NodAcc(I, INode, iEntity)
                   ! Particle x-velocity
                   ParticleIncrementalVelocity(I) = ParticleIncrementalVelocity(I) + Time * ParticleShape(INode) * NodAcc(I, INode, iEntity)
                 end do    

               end do !Loop over nodes
               
                 ParticleVelocity = ParticleVelocity +  ParticleIncrementalVelocity

                 do I = 1, NVECTOR ! nodal x-momentum
                   Momentum(IDof(I,1:ELEMENTNODES), iEntity) = Momentum(IDof(I,1:ELEMENTNODES), iEntity) + PartilceMass * ParticleShape * ParticleVelocity(I)
                 end do
                 
                    
                 if (Particles(ParticleIndex)%DegreeSaturation>=1) then
                     ParticleVelocity = 0.0
                     do I = 1, NVECTOR ! nodal x-momentum
                       Momentum(IDof(I,1:ELEMENTNODES), iEntity) = 0.0
                     end do
                 end if
                                      
                  VelocityGasArray(ParticleIndex,:) = ParticleVelocity
         
            end do !Loop over particles
          end do !elements
          
        end subroutine UpdateParticleGasVelocityAndMapMomentumG  
        
        subroutine CalculateGasVolumetricStrain(DUTotAll)
        !**********************************************************************
        !
        !    Function:  Calculate the volumetric strain (gas).
        !
        !     DUTot : Incremental nodal displacements (gas)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
          implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTotAll              
          !Local variables
          real(REAL_TYPE), dimension(Counters%N) :: DUTotEnt        
          real(REAL_TYPE), dimension (NTENSOR) :: Eps
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac
          integer(INTEGER_TYPE) :: IAElement, IElement, NElemPart, IParticle, ParticleIndex, I, iEntity, iEntityDefault   
     
          !========================================================            
          !set the nodal displacement increments to the first entity
          !only need to change if the contact model is used and 
          !the next particle belongs to a different entity
          iEntityDefault = 1    !set default to first entity
          !=======================================================

          do IAElement = 1, Counters%NAEl     !loop over all elements
            IElement = ActiveElement(IAElement)
     
              NElemPart = NPartEle(IElement)    !number of particles in element
              if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
              end if
              
              do IParticle = 1, NElemPart       !loop over all particles of the element
              
                ParticleIndex = GetParticleIndex(IParticle, IElement)   !particle global ID

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
                    do I=1,Counters%N                              
                      DUTotEnt(I) = DUTotAll(I,iEntity)
                    end do
                    iEntityDefault = iEntity                                !set the entity ID for the next particle loop
                  end if
                else
                  do I = 1,Counters%N
                    DUTotEnt(I) = DUTotAll(I, 1)
                  end do
                end if
                !====================================         
     
                ! Eps = (Exx, Eyy, Ezz, Gxy, Gyz, Gzx)
                call Get_Strain(IElement, IParticle, &
                                ElementConnectivities, & !GetStrain.FOR
                                BMatrixDeformed, DUTotEnt, ReducedDof, &
                                Eps)
                 
                    Particles(ParticleIndex)%GasVolumetricStrain = Eps(1) + Eps(2) + Eps(3) ! for 2D and 3D valid
                    
              end do
          end  do

        end subroutine CalculateGasVolumetricStrain
        
        
        subroutine CalculateParticleWaterAdvectiveFlux(IncrementalDisplacementWater)
        !**********************************************************************
        !   
        !   Function: Calculate the advective flux of water at the material point
        !
        !**********************************************************************
       
          implicit none
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: IncrementalDisplacementWater
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac
          real(REAL_TYPE), dimension(NVECTOR) :: Eps, WminusV
          integer(INTEGER_TYPE) :: I, IAElement, IElement, NElemPart, IParticle, ParticleIndex, J, iEntity,NN

          do IAElement = 1, Counters%NAEl     !loop over all elements
            IElement = ActiveElement(IAElement)
     
              NElemPart = NPartEle(IElement)    !number of particles in element
              
              do IParticle = 1, NElemPart       !loop over all particles of the element
              
                ParticleIndex = GetParticleIndex(IParticle, IElement)   !particle global ID
                Particles(ParticleIndex)%WaterAdvectiveFlux = 0

                ! Determine B matrix of deformed element at the particle local position
                call BMatrix(Particles(ParticleIndex)%LocPos, &
                             ELEMENTNODES, Counters%NEl,  &
                             Counters%NodTot, NVECTOR, &
                             IElement, ElementConnectivities, &
                             NodalCoordinatesUpd,  &
                             BMatrixDeformed, DetJac)

                if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex) 
                else
                  ientity = 1
                end if

                Eps = 0
                WminusV = 0

                do J= 1, ELEMENTNODES !loop over all nodes of the element
                  NN = ElementConnectivities(J,IElement)
                 
                  do I = 1, NVECTOR
                    WminusV(I) = IncrementalDisplacementWater(ReducedDof(NN)+I,ientity) - IncrementalDisplacementSoil(ReducedDof(NN)+I,ientity)
                    Eps(I) = BMatrixDeformed(I,J) * AdvectiveFluxDarcyWater(ReducedDof(NN)+I,ientity) * WminusV(I)
                  end do  
                  
                  if (NVECTOR == 3) then ! 3D
                    Particles(ParticleIndex)%WaterAdvectiveFlux = Particles(ParticleIndex)%WaterAdvectiveFlux + Eps(1) + Eps(2) + Eps(3) 
                  elseif (NVECTOR == 2) then ! 2D
                    Particles(ParticleIndex)%WaterAdvectiveFlux = Particles(ParticleIndex)%WaterAdvectiveFlux + Eps(1) + Eps(2)
                  end if  

                end do  !loop over all nodes of the element
              end do !loop over all particles of the element
          end do !loop over all elements

        end subroutine CalculateParticleWaterAdvectiveFlux

        subroutine CalculateParticleAirAdvectiveFlux(NodalIncDisplacementG)
        !**********************************************************************
        !   
        !   Function: Calculate the advective flux of air at the material point
        !
        !**********************************************************************
       
          implicit none
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: NodalIncDisplacementG
          real(REAL_TYPE), dimension(NDIM, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac
          real(REAL_TYPE), dimension (NVECTOR) :: Eps, WgminusV
          integer(INTEGER_TYPE) :: I, IAElement, IElement, NElemPart, IParticle, ParticleIndex, J, iEntity, NN

          do IAElement = 1, Counters%NAEl     !loop over all elements
            IElement = ActiveElement(IAElement)

            NElemPart = NPartEle(IElement)    !number of particles in element
              
            do IParticle = 1, NElemPart       !loop over all particles of the element

              ParticleIndex = GetParticleIndex(IParticle, IElement)   !particle global ID
              Particles(ParticleIndex)%AirAdvectiveFlux = 0

              ! Determine B matrix of deformed element at the particle local position
              call BMatrix(Particles(ParticleIndex)%LocPos, &
                             ELEMENTNODES, Counters%NEl, &
                             Counters%NodTot, NVECTOR, &
                             IElement, ElementConnectivities, &
                             NodalCoordinatesUpd, &
                             BMatrixDeformed, DetJac)


              if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex) 
              else
                ientity = 1
              end if

              Eps = 0
              WgminusV = 0

              do J= 1, ELEMENTNODES !loop over all nodes of the element
                NN = ElementConnectivities(J,IElement)

                do I = 1, NVECTOR
                  WgminusV(I) = NodalIncDisplacementG(ReducedDof(NN)+I, ientity) - IncrementalDisplacementSoil(ReducedDof(NN)+I, ientity)
                  Eps(I) = BMatrixDeformed(I,J) * AdvectiveFluxDarcyAir(ReducedDof(NN)+I, ientity) * WgminusV(I)
                end do  

                if (NVECTOR == 3) then ! 3D
                  Particles(ParticleIndex)%AirAdvectiveFlux = Particles(ParticleIndex)%AirAdvectiveFlux + Eps(1) + Eps(2) + Eps(3)
                elseif (NVECTOR == 2) then
                  Particles(ParticleIndex)%AirAdvectiveFlux = Particles(ParticleIndex)%AirAdvectiveFlux + Eps(1) + Eps(2)
                end if  

              end do  !loop over all nodes of the element
            end do !loop over all particles of the element
          end do !loop over all elements

        end subroutine CalculateParticleAirAdvectiveFlux

        
         subroutine GetGravityLoadGasPorosity(GravityLoadGasLoc)
        !**********************************************************************
        !
        !    Function:  To extrapolate real gas load from particles to nodes using 
        !               the shape function values evaluated at the particles local 
        !               position.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
       
         real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: GravityLoadGasLoc
         
         ! local variables
          integer(INTEGER_TYPE) :: IAElement, IPart, INode, IElement, I, ParticleIndex, NodeID, GlobDof, iEntity
          real(REAL_TYPE), dimension(NVECTOR) :: PartGravity

          !must set to zero here, because used in summation
          GravityLoadGasLoc = 0.0

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
              PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyGas * Particles(ParticleIndex)%Porosity * (1.0 - Particles(ParticleIndex)%DegreeSaturation)

              do I = 1, NVECTOR
                GravityLoadGasLoc(GlobDof+I, iEntity) =  GravityLoadGasLoc(GlobDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
              end do 

              end do ! nodes

          end do ! particles
         
          end do ! elements
          
          if (IS3DCYLINDRIC) then ! rotation is needed
            do IEntity = 1, Counters%nEntity
              call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, GravityLoadGas(:, IEntity), GravityLoadGas(:, IEntity))
            end do
          end if ! rotation

        end subroutine GetGravityLoadGasPorosity

        subroutine FormMatricesGas(LumpedMassGasLoc, LumpedMassNGasLoc, ConductivityGasMatrixLoc, ConductivityGasMatrixPorosityDegreeSatLoc)
        !**********************************************************************
        !
        !    Function:  To Form the Lumped mass matrix (gas)
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: &
                             LumpedMassGasLoc,  &
                             LumpedMassNGasLoc, &
                             ConductivityGasMatrixLoc, &
                             ConductivityGasMatrixPorosityDegreeSatLoc

          ! Local variables
          integer(INTEGER_TYPE) :: I, IAEl, IEl, IPart, INode, IDof, ParticleIndex, NodeID, iEntity
     
          real(REAL_TYPE) :: MassEntry, MassNEntry, ConductivityGasEntry, ConductivityGasEntryPorosityDegreeSatLoc
                       
          LumpedMassGasLoc = 0.0
          LumpedMassNGasLoc = 0.0
          ConductivityGasMatrixLoc = 0.0
          ConductivityGasMatrixPorosityDegreeSatLoc = 0.0
           
          do IAEl = 1, Counters%NAEl                                      ! Loop over all elements
            IEl = ActiveElement(IAEl)
            do IPart = 1, NPartEle(IEl)                                   ! Loop over all particles in element
                
                ParticleIndex = GetParticleIndex(IPart, IEl)              ! Get the particle ID
                if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex)            ! get entity to which particle belongs
                else
                  iEntity = 1
                end if
              
                do INode = 1, ELEMENTNODES                             ! loop over element nodes
                  NodeID = ElementConnectivities(INode, IEl)              ! Global node ID
                  IDof = ReducedDof(NodeID)
                    
                  MassEntry =  Particles(ParticleIndex)%MassGas * ShapeValuesArray(ParticleIndex,INode)           ! calculate mass contribution
                    
                  do I = 1, NVECTOR       ! nodal mass 
                    LumpedMassGasLoc(IDof+I,iEntity) = LumpedMassGasLoc(IDof+I,iEntity) + MassEntry
                  end do  
                  ! calculate mass contribution   
                  MassNEntry = MassEntry * Particles(ParticleIndex)%Porosity * (1 - Particles(ParticleIndex)%DegreeSaturation)       
     
                  do I = 1, NVECTOR       ! nodal mass 
                    LumpedMassNGasLoc(IDof+I,iEntity) = LumpedMassNGasLoc(IDof+I,iEntity) + MassNEntry
                  end do  

                  if (Particles(ParticleIndex)%ConductivityGas/=0.0) then
                    ! calculate conductivity contribution
                      ConductivityGasEntry =   &
                      Particles(ParticleIndex)%MassGas * &
                      Particles(ParticleIndex)%Porosity * &
                      (1.0 - Particles(ParticleIndex)%DegreeSaturation) * &
                      CalParams%GravityData%GAccel *  &
                      (1.0 / Particles(ParticleIndex)%ConductivityGas) * &
                      ShapeValuesArray(ParticleIndex,INode)         
                    
                    do I = 1, NVECTOR ! nodal conductivity
                      ConductivityGasMatrixLoc(IDof+I,iEntity)= ConductivityGasMatrixLoc(IDof+I,iEntity) + ConductivityGasEntry
                    end do
                    
                    ! calculate conductivity contribution
                    ConductivityGasEntryPorosityDegreeSatLoc =   &
                      ConductivityGasEntry * &
                      Particles(ParticleIndex)%Porosity * &
                      (1.0 - Particles(ParticleIndex)%DegreeSaturation)    
                    
                    do I = 1, NVECTOR! nodal conductivity
                      ConductivityGasMatrixPorosityDegreeSatLoc(IDof+I,iEntity)=  ConductivityGasMatrixPorosityDegreeSatLoc(IDof+I,iEntity) + ConductivityGasEntryPorosityDegreeSatLoc
                    end do

                  end if
               
                end do !Loop over nodes
            end do !Loop over particles
          end do !Loop over elements

        end subroutine FormMatricesGas

        subroutine ConsolidationForcesGas(ExtLoadGasLoc, GravityLoadGasLoc, IntLoadGasLoc)
        !**********************************************************************
        !
        !    Function:  Calculation of the equivalent nodal forces due to
        !               a given stress InternalLD = Integral {BT*Sigma}.
        !
        ! O ExternalLD : External load
        ! O InternalLD : Internal load
        ! O GravityLD : Gravity load
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: ExtLoadGasLoc, IntLoadGasLoc, GravityLoadGasLoc
          ! Local variables
          real(REAL_TYPE), &
            dimension(NVECTOR, ELEMENTNODES) :: B
          integer(INTEGER_TYPE) :: IntGlo, IEl, Int, I, NN, NElemPart, iEntityID, iNode, IDof
          real(REAL_TYPE) :: WtN, Det, S
          integer(INTEGER_TYPE) :: IPart, IAEl, ParticleIndex, iEntity
          real(REAL_TYPE), dimension(NVECTOR) :: PartLoad, PartGravity, PartGravityMix
          real(REAL_TYPE) :: Position
          real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)                   

          ExtLoadGasLoc = 0.0
          IntLoadGasLoc = 0.0
          GravityLoadGasLoc = 0.0
          
          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
!//////////////////// External force calculation //////////////////////// 
          do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
            ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
            PartLoad = 0.0
            PartGravity = 0.0
            PartGravityMix = 0.0
            if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex)
            else
                iEntity = 1
            endif  

               do INode = 1, ELEMENTNODES 
                    nn = ElementConnectivities(INode, IEl)  ! Nodal global ID
                    IDof = ReducedDof(nn)

                  if(Counters%NLoadedElementSidesGasMatPoints > 0) then 
                    PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExtGas
                    do I = 1, NVECTOR
                      ExtLoadGasLoc(IDof+I, iEntity) =  ExtLoadGasLoc(IDof+I, iEntity) + PartLoad(I) * CalParams%Multipliers%GasACurrent
                    end do  
                  end if

                    PartGravity =  ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyGas

                do I = 1, NVECTOR
                 GravityLoadGasLoc(IDof+I, iEntity) =  GravityLoadGasLoc(IDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
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
                  if ( ISAXISYMMETRIC ) then
                    ! the integration weight of the MP is not corrected due to axisymmetry
                    Position = GlobPosArray(IntGlo, 1) ! index 1 is readial directoin
                    ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                  end if
                else
                  if ( ISAXISYMMETRIC ) then
                    !Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                    !ShapeValues(:) = GPShapeFunction(Int, :)
                    WTN = WTN * Position ! the volume of the element is corrected due to axisymmetry
                  end if  
                end if
                
                  S = (Particles(IntGlo)%GasPressure) * WTN
        
                  !get particle entity
                  if (.not.CalParams%ApplyContactAlgorithm) then      
                    iEntityID = 1
                  else                                           
                    iEntityID = EntityIDArray(IntGlo)           
                  end if        
                
                  do INode=1,ELEMENTNODES  !loop through element nodes

                    nn=ElementConnectivities(iNode,iel) ! get global node number
                    IDof = ReducedDof(nn)   ! global storage coordinate 
                    
                    do I = 1, NVECTOR ! nodal load 
                      IntLoadGasLoc(IDof+I, iEntityID) =  IntLoadGasLoc(IDof+I, iEntityID) + B(I,INode)*S 
                    end do
                    
                    if ( ISAXISYMMETRIC ) then
                      ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                      IntLoadGasLoc(IDof+1, iEntityID) = IntLoadGasLoc(IDof+1, iEntityID) + S * ShapeValues(INode) / Position
                    end if
                    
                end do !node loop             
              end do ! Loop over 1 .. NElemPart
              !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
          end do ! Loop over elements

        end subroutine ConsolidationForcesGas
        
        subroutine ConsolidationForcesBishop(BishopIntLoadLoc)
        !**********************************************************************
        !
        !    Function:  Calculation of the equivalent nodal forces due to
        !               a given stress NodalForce = Integral {BT*Sigma}.
        !
        ! O BishopIntLoadLoc : 
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: BishopIntLoadLoc

          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
          integer(INTEGER_TYPE) :: I, IntGlo, IEl, Int, NN,  NElemPart, iEntityID, iNode, IDof
          real(REAL_TYPE) :: WtN, Det, S, S1, S2, BishopParameter, Pl, Pg
          integer(INTEGER_TYPE) :: IAEl
          real(REAL_TYPE) :: Position
          real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)
                   
          BishopIntLoadLoc = 0.0


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
                  
                ! Set the integration weight
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
                
                  Pl = Particles(IntGlo)%WaterPressure
                  Pg = Particles(IntGlo)%GasPressure
        
                  call CalculateBishopParameter_SrAlpha(IntGlo, BishopParameter)
                  !call CalculateBishopParameter_Gesto(IntGlo, BishopParameter)
                    
                  S1 = ((1.0d0 - BishopParameter) * Pg + BishopParameter * Pl) * WTN !unsaturated material point --->  TotalSig = EffSig + Pg - BishopParameter*(Pg-Pl) = EffSig + S1
                  S2 = Pl * WTN   !saturated material point -------------------------------------------------------->  TotalSig = EffSig + Pl = EffSig + S2
                  
                  S = min(S1,S2)
                  
                  !get particle entity
                  if (.not.CalParams%ApplyContactAlgorithm) then      
                    iEntityID = 1
                  else                                           
                    iEntityID = EntityIDArray(IntGlo) 
                  end if        
                
                  do INode=1,ELEMENTNODES  !loop through element nodes

                    nn=ElementConnectivities(iNode,iel) ! get global node number
                    IDof = ReducedDof(nn)   ! global storage coordinate
                    
                    ! Calculation of Force 
                    do I = 1, NVECTOR
                      BishopIntLoadLoc(IDof+I, iEntityID) =  BishopIntLoadLoc(IDof+I, iEntityID) + B(I,INode)*S 
                    end do  
                    
                    if ( ISAXISYMMETRIC ) then
                      ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                      BishopIntLoadLoc(IDof+1, iEntityID) = BishopIntLoadLoc(IDof+1, iEntityID) + S * ShapeValues(INode) / Position
                    end if
                                             
                  end do !node loop             
              end do ! Loop over 1 .. NElemPart
              !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
          end do ! Loop over elements

        end subroutine ConsolidationForcesBishop

  
        subroutine CalculateBishopParameter_Gesto(IntGlo,BP)
        !**********************************************************************
        !
        !    Function: Calculates the Bishop Parameter (BP), which is the "Effective" Degree of Saturation
        !               The current used expression is based on Gesto et al. (2011)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        integer(INTEGER_TYPE), intent(in) :: IntGlo
        real(REAL_TYPE), intent(inout) :: BP
        ! Local variables
        real(REAL_TYPE) :: Sr, N
        real(REAL_TYPE) :: Smooth !Degree of smoothing 
        real(REAL_TYPE) :: X !Microstructural State Variable = Sr_min 
        real(REAL_TYPE) :: e, em !Void Ratio and Microstructural Void Ratio
        real(REAL_TYPE) :: v1, v2, v3

            Sr = Particles(IntGlo)%DegreeSaturation
            N = Particles(IntGlo)%Porosity
            
            e = N/(1-N)         !Void ratio
            em = 0.1         !999999 Microstructural Void Ratio
            X = em/e            !Microstructural State Variable = Sr_min
            Smooth = 5.0d0      !Degree of smoothing, can be 3,5,10... (the higher the number-->the lower the smoothing)
            
            
            v1 = Sr-X
            v2 = 1.0d0 - X
            v3 = v1/v3
            
            BP = v3+(1/N)*log(1+exp(-Smooth*v3))
            

        end subroutine CalculateBishopParameter_Gesto
        
        
        
        subroutine CalculateBishopParameter_SrAlpha(IntGlo,BP)
        !**********************************************************************
        !
        !    Function: Calculates the Bishop Parameter (BP), which is the "Effective" Degree of Saturation
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        integer(INTEGER_TYPE), intent(in) :: IntGlo
        real(REAL_TYPE), intent(inout) :: BP
        ! Local variables
        real(REAL_TYPE) :: Sr, Alpha
        integer(INTEGER_TYPE) :: ISet

          ISet = MaterialIDArray(IntGlo)

          Sr = Particles(IntGlo)%DegreeSaturation
          Alpha = MatParams(ISet)%BishopsAlpha
          BP = Sr**Alpha
            
        end subroutine CalculateBishopParameter_SrAlpha
        
        
        
        subroutine ConsolidationForcesPorosityDegreeSat(ExtLoadWaterPorosityLoc,ExtLoadGasPorosityLoc, &
                                       GravityLoadWaterPorosityLoc,GravityLoadGasPorosityLoc, &
                                       IntLoadWaterPorosityLoc,IntLoadGasPorosityLoc)
        !**********************************************************************
        !
        !    Function:  Calculation of the equivalent nodal forces due to
        !               a given stress InternalLD = Integral {BT*Sigma}.
        !
        ! O ExternalLDPorosity : Real External load
        ! O InternalLDPorosity : Real Internal load
        ! O GravityLDPorosity : Real Gravity load
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), &
            intent(inout) :: ExtLoadWaterPorosityLoc, IntLoadWaterPorosityLoc, GravityLoadWaterPorosityLoc, &
                             ExtLoadGasPorosityLoc, IntLoadGasPorosityLoc, GravityLoadGasPorosityLoc
          ! Local variables
          real(REAL_TYPE), dimension(NDIM, ELEMENTNODES) :: B
          integer(INTEGER_TYPE) :: IntGlo, IEl, Int, I, NN, NElemPart, iEntityID, iNode, IDof
          real(REAL_TYPE) :: WtN, Det, SWater, SGas
          integer(INTEGER_TYPE) :: IPart, IAEl, ParticleIndex, NodeID, GlobDof, iEntity
          real(REAL_TYPE), dimension(NVECTOR) :: PartLoad, PartGravity
          real(REAL_TYPE) :: Position
          real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)
                   
          ExtLoadWaterPorosityLoc = 0.0
          IntLoadWaterPorosityLoc = 0.0
          GravityLoadWaterPorosityLoc = 0.0
          ExtLoadGasPorosityLoc = 0.0
          IntLoadGasPorosityLoc = 0.0
          GravityLoadGasPorosityLoc = 0.0

          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
!//////////////////// External force calculation //////////////////////// 
          do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
            ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
            PartLoad = 0.0
            PartGravity = 0.0
            if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex) !
            else
                iEntity = 1
            endif  

               do INode = 1, ELEMENTNODES 
                    NodeID = ElementConnectivities(INode, IEl)  ! Nodal global ID
                    GlobDof = ReducedDof(NodeID)

                  if(Counters%NLoadedElementSidesWaterMatPoints > 0) then 
                    PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExtWater * Particles(ParticleIndex)%Porosity * Particles(ParticleIndex)%DegreeSaturation
                    do I = 1, NVECTOR
                      ExtLoadWaterPorosityLoc(GlobDof+I, iEntity) = ExtLoadWaterPorosityLoc(GlobDof+I, iEntity) + PartLoad(I) * CalParams%Multipliers%WaterACurrent
                    end do
                  end if

                  if(Counters%NLoadedElementSidesGasMatPoints > 0) then 
                    PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExtGas * Particles(ParticleIndex)%Porosity * (1.0 - Particles(ParticleIndex)%DegreeSaturation)
                    do I = 1, NVECTOR
                      ExtLoadGasPorosityLoc(GlobDof+I, iEntity) = ExtLoadGasPorosityLoc(GlobDof+I, iEntity) + PartLoad(I) * CalParams%Multipliers%GasACurrent
                    end do
                  end if

                    PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyWater *  Particles(ParticleIndex)%Porosity * Particles(ParticleIndex)%DegreeSaturation

                  do I = 1, NVECTOR  
                    GravityLoadWaterPorosityLoc(GlobDof+I, iEntity) = GravityLoadWaterPorosityLoc(GlobDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                  end do
                    
                  PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyGas * Particles(ParticleIndex)%Porosity * (1.0 - Particles(ParticleIndex)%DegreeSaturation)

                  do I = 1, NVECTOR
                    GravityLoadGasPorosityLoc(GlobDof+I, iEntity) = GravityLoadGasPorosityLoc(GlobDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
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
                
                  SWater = Particles(IntGlo)%WaterPressure * WTN * Particles(ParticleIndex)%Porosity * Particles(ParticleIndex)%DegreeSaturation
                  SGas = Particles(IntGlo)%GasPressure * WTN * Particles(ParticleIndex)%Porosity * (1.0 - Particles(ParticleIndex)%DegreeSaturation)
        
                  !get particle entity
                  if (.not.CalParams%ApplyContactAlgorithm) then      
                    iEntityID = 1
                  else                                           
                    iEntityID = EntityIDArray(IntGlo)            
                  end if        
                
                  do INode=1,ELEMENTNODES  !loop through element nodes

                    nn=ElementConnectivities(iNode,iel) ! get global node number
                    IDof = ReducedDof(nn)  ! global storage coordinate of x-val
                
                    do I = 1, NVECTOR ! nodal load 
                      IntLoadWaterPorosityLoc(IDof+I, iEntityID) = IntLoadWaterPorosityLoc(IDof+I, iEntityID) + B(I,INode)*SWater 
                      IntLoadGasPorosityLoc(IDof+I, iEntityID) = IntLoadGasPorosityLoc(IDof+I, iEntityID) + B(I,INode)*SGas 
                    end do
                    
                    if ( ISAXISYMMETRIC ) then
                      ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                      IntLoadWaterPorosityLoc(IDof+1, iEntityID) = IntLoadWaterPorosityLoc(IDof+1, iEntityID) + SWater * ShapeValues(INode) / Position
                      IntLoadGasPorosityLoc(IDof+1, iEntityID) = IntLoadGasPorosityLoc(IDof+1, iEntityID) + SGas * ShapeValues(INode) / Position
                    end if
                      
                end do !node loop             
              end do ! Loop over 1 .. NElemPart
              !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
          end do ! Loop over elements

        end subroutine ConsolidationForcesPorosityDegreeSat

      end module ModMPMDYNUnsatConsolidation